"""Diagnose the E6 E2c regression mismatch without regenerating layouts."""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[3]
S = ROOT / "paper/strengthening"
SOURCE = S / "source"
QA = S / "qa/E6"
CGM = ROOT
sys.path[:0] = [str(SOURCE), str(ROOT / "paper"), str(ROOT / "src")]

import numpy as np
import pandas as pd

from e1b_neighbors import neighborhood_metrics, spearman_rho
from e2c_analysis import ACTIVITIES, curate_activity_table, subset_positions, summarize_entries
from chemgridmap.reliability import grid_neighbor_sets, representation_neighbor_sets
from chemgridmap.representations import morgan_fingerprints

K = 10
LAYOUTS = ("projection", "minimum_displacement", "neighbor_refined")
E2C_RESULTS = S / "output/E2c/layout_by_target.csv"


def maximum_difference(observed, expected) -> float:
    return float(np.max(np.abs(
        np.asarray(observed, dtype=float) - np.asarray(expected, dtype=float)
    )))


def diagnose_dataset(payload) -> tuple[tuple[str, str], list[dict]]:
    key, items = payload
    target, endpoint = key
    path = QA / "e2c" / f"{target}_{endpoint}.npz"
    if not path.exists():
        raise FileNotFoundError(path)
    with np.load(path, allow_pickle=False) as cache:
        activity = cache["activity"]
        cached_keys = cache["molecule_identity_key"].astype(str)
    records = curate_activity_table(
        pd.read_csv(
            ACTIVITIES / f"{target}_{endpoint}.csv",
            low_memory=False,
        ),
        target,
        endpoint,
    ).retained_records
    entries = summarize_entries(records)
    positions = subset_positions(entries.molecule_identity_key.to_numpy())
    selected = entries.iloc[positions].reset_index(drop=True)
    if not np.array_equal(cached_keys, selected.molecule_identity_key.to_numpy(str)):
        raise ValueError(f"Identity mismatch: {target}/{endpoint}")
    fingerprints = morgan_fingerprints(selected.canonical_smiles.tolist())
    reference = representation_neighbor_sets(fingerprints, k=K)
    rows = []
    for item in items:
        with np.load(path, allow_pickle=False) as cache:
            coords = cache[item["layout"]]
        metrics = neighborhood_metrics(
            activity, reference, grid_neighbor_sets(coords, k=K)
        )
        row = {
            "target": target,
            "endpoint": endpoint,
            "layout": item["layout"],
            "mean_recall_observed": float(metrics.recall.mean()),
            "mean_e_grid_observed": float(metrics.e_grid.mean()),
            "mean_e_fp_observed": float(metrics.e_fp.mean()),
            "rho_observed": spearman_rho(metrics.recall, metrics.excess_error),
            "mean_recall_frozen": float(item["mean_recall"]),
            "mean_e_grid_frozen": float(item["mean_e_grid"]),
            "mean_e_fp_frozen": float(item["mean_e_fp"]),
            "rho_frozen": float(item["rho"]),
        }
        row["recall_difference"] = row["mean_recall_observed"] - row["mean_recall_frozen"]
        row["e_grid_difference"] = row["mean_e_grid_observed"] - row["mean_e_grid_frozen"]
        row["e_fp_difference"] = row["mean_e_fp_observed"] - row["mean_e_fp_frozen"]
        row["rho_difference"] = row["rho_observed"] - row["rho_frozen"]
        rows.append(row)
    return (target, endpoint), rows


def run() -> pd.DataFrame:
    frozen = pd.read_csv(E2C_RESULTS, float_precision="round_trip")
    included = frozen.loc[
        frozen.included.astype(bool) & frozen.layout.isin(LAYOUTS)
    ].copy()
    payloads = [
        (key, group.to_dict("records"))
        for key, group in included.groupby(["target", "endpoint"], sort=True)
    ]
    rows = []
    with ProcessPoolExecutor(max_workers=4) as executor:
        futures = {
            executor.submit(diagnose_dataset, payload): payload[0]
            for payload in payloads
        }
        for index, future in enumerate(as_completed(futures), 1):
            key, dataset_rows = future.result()
            rows.extend(dataset_rows)
            print(
                f"E6 regression diagnosis {index}/77 {key[0]} {key[1]}",
                flush=True,
            )
    table = pd.DataFrame(rows)
    table["within_tolerance"] = (
        table[["recall_difference", "e_grid_difference", "e_fp_difference", "rho_difference"]]
        .abs()
        .max(axis=1)
        <= 1e-9
    )
    table.to_csv(QA / "e2c_regression_rows.csv", index=False)
    mismatches = table.loc[~table.within_tolerance]
    mismatch_dataset_keys = [
        f"{target}/{endpoint}"
        for target, endpoint in mismatches[["target", "endpoint"]]
        .drop_duplicates()
        .itertuples(index=False, name=None)
    ]
    passed = mismatches.empty
    summary = {
        "rows": int(len(table)),
        "mismatch_rows": int(len(mismatches)),
        "mismatch_datasets": len(mismatch_dataset_keys),
        "mismatch_dataset_keys": mismatch_dataset_keys,
        "mismatches_by_layout": {
            str(key): int(value)
            for key, value in mismatches.layout.value_counts().items()
        },
        "max_absolute_differences": {
            "mean_recall": maximum_difference(
                table.mean_recall_observed, table.mean_recall_frozen
            ),
            "mean_e_grid": maximum_difference(
                table.mean_e_grid_observed, table.mean_e_grid_frozen
            ),
            "mean_e_fp": maximum_difference(
                table.mean_e_fp_observed, table.mean_e_fp_frozen
            ),
            "rho": maximum_difference(table.rho_observed, table.rho_frozen),
        },
        "status": "passed" if passed else "mismatch_recorded",
        "action": (
            "proceed with E2c fixed-k and L1/L2 analyses"
            if passed
            else "review listed mismatch datasets before E2c fixed-k and L1/L2 analyses"
        ),
    }
    (QA / "e2c_regression_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return table


if __name__ == "__main__":
    result = run()
    print(
        f"E6 E2c mismatch rows={int((~result.within_tolerance).sum())}/{len(result)}",
        flush=True,
    )
    print(
        result.loc[
            ~result.within_tolerance,
            ["target", "endpoint", "layout", "recall_difference", "e_grid_difference", "rho_difference"],
        ].to_string(index=False),
        flush=True,
    )
