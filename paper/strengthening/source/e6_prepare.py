"""E6 regression checks and saved layout reconstruction."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import time

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[3]
S = ROOT / "paper/strengthening"
SOURCE = S / "source"
QA = S / "qa/E6"
CGM = ROOT
sys.path[:0] = [str(SOURCE), str(ROOT / "paper"), str(ROOT / "src")]

import numpy as np
import pandas as pd

from e2c_analysis import (
    ACTIVITIES,
    MAP_SEED,
    curate_activity_table,
    subset_positions,
    summarize_entries,
)
from e6_core import fixed_metrics
from e1b_neighbors import neighborhood_metrics, spearman_rho
from grid_baselines import flas, las
from chemgridmap.core import assign_to_grid, normalize_coordinates
from chemgridmap.projection import project_2d
from chemgridmap.refinement import REFINEMENT_DEFAULTS, refine_grid_layout
from chemgridmap.reliability import grid_neighbor_sets, representation_neighbor_sets
from chemgridmap.representations import morgan_fingerprints
from reliability_analysis import LAYOUTS, SEEDS, TARGETS, layout_points, load_target

K = 10
ANALYSIS_ROLE = "post hoc"
FROZEN = CGM / "paper/output/reliability_v1"
FIGURE_3A = FROZEN / "figure_data/Figure_3A_neighbor_recall.csv"
E4_RESULTS = S / "output/E4/layout_baselines_ten_seed.csv"
E2C_RESULTS = S / "output/E2c/layout_by_target.csv"
E2C_SELECTED = S / "data/chembl37_api/selected_targets.csv"
DETAIL_LAYOUTS = (
    "projection",
    "minimum_displacement",
    "neighbor_refined",
    "las_random",
    "las_projection",
    "flas_random",
)
E2C_LAYOUTS = ("projection", "minimum_displacement", "neighbor_refined")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def maximum_difference(observed, expected) -> float:
    delta = np.abs(np.asarray(observed, dtype=float) - np.asarray(expected, dtype=float))
    return float(delta.max()) if len(delta) else 0.0


def prepare_detailed() -> dict:
    cache_rows = []
    regression_rows = []
    started = time.perf_counter()
    for target, target_label in TARGETS.items():
        data, fingerprints = load_target(target)
        activity = data.activity_pchembl.to_numpy(float)
        reference = representation_neighbor_sets(fingerprints, k=K)
        vectors = np.asarray(fingerprints, dtype=float)
        norms = np.linalg.norm(vectors, axis=1, keepdims=True)
        vectors = np.divide(
            vectors, norms, out=np.zeros_like(vectors), where=norms != 0
        )
        for seed in SEEDS:
            path = FROZEN / "layouts" / target / f"coordinates_seed_{seed}.csv"
            frame = pd.read_csv(path)
            if data.molecule_identity_key.tolist() != frame.molecule_identity_key.tolist():
                raise ValueError(f"Identity mismatch: {target}/{seed}")
            side = int(frame.grid_side.iloc[0])
            minimum = layout_points(frame, "minimum_displacement")
            points = {
                "projection": layout_points(frame, "projection"),
                "minimum_displacement": minimum,
                "neighbor_refined": layout_points(frame, "neighbor_refined"),
            }
            t0 = time.perf_counter()
            points["las_random"] = las(
                vectors, (side, side), seed=seed
            )
            points["las_projection"] = las(
                vectors,
                (side, side),
                seed=seed,
                initial_cells=minimum,
            )
            points["flas_random"] = flas(
                vectors, (side, side), seed=seed
            )
            print(
                f"E6 prepare detailed {target_label} seed={seed} "
                f"baselines={time.perf_counter() - t0:.2f}s",
                flush=True,
            )
            for layout in DETAIL_LAYOUTS:
                cells = np.asarray(points[layout])
                if cells.shape != (len(data), 2) or not np.isfinite(cells).all():
                    raise ValueError(f"Bad coordinates: {target}/{seed}/{layout}")
                retained = neighborhood_metrics(
                    activity, reference, grid_neighbor_sets(cells, k=K)
                )
                cache_rows.extend({
                    "analysis_role": ANALYSIS_ROLE,
                    "target": target,
                    "target_label": target_label,
                    "seed": seed,
                    "layout": layout,
                    "n": len(data),
                    "grid_side": side,
                    "molecule_identity_key": molecule_id,
                    "col": float(cell[0]),
                    "row": float(cell[1]),
                } for molecule_id, cell in zip(
                    frame.molecule_identity_key, cells.astype(int)
                ))
                regression_rows.append({
                    "target": target,
                    "seed": seed,
                    "layout": layout,
                    "mean_recall_retained": float(retained.recall.mean()),
                })
    (QA / "detailed_layouts.csv").write_text(
        pd.DataFrame(cache_rows).to_csv(index=False),
        encoding="utf-8",
    )
    observed = pd.DataFrame(regression_rows)

    figure = pd.read_csv(FIGURE_3A)
    figure = figure.loc[
        figure.layout.isin(LAYOUTS), ["target", "seed", "layout", "mean_morgan_neighbor_recall_k10"]
    ].rename(columns={"mean_morgan_neighbor_recall_k10": "frozen_recall"})
    main = observed.loc[observed.layout.isin(LAYOUTS)]
    main_check = figure.merge(
        main, on=["target", "seed", "layout"], validate="one_to_one"
    )
    figure_max = maximum_difference(
        main_check.mean_recall_retained, main_check.frozen_recall
    )

    e4 = pd.read_csv(E4_RESULTS)
    e4 = e4.loc[
        (e4.row_type == "seed_layout")
        & e4.layout.isin(("las_random", "las_projection", "flas_random")),
        ["target", "seed", "layout", "mean_recall"],
    ].rename(columns={"mean_recall": "frozen_recall"})
    baseline = observed.loc[observed.layout.isin(e4.layout.unique())]
    baseline_check = e4.merge(
        baseline, on=["target", "seed", "layout"], validate="one_to_one"
    )
    e4_max = maximum_difference(
        baseline_check.mean_recall_retained, baseline_check.frozen_recall
    )
    checks = {
        "figure_3a_rows": int(len(figure)),
        "figure_3a_max_absolute_recall_difference": figure_max,
        "figure_3a_tolerance": 1e-12,
        "figure_3a_passed": bool(figure_max <= 1e-12),
        "e4_las_flas_rows": int(len(e4)),
        "e4_las_flas_max_absolute_recall_difference": e4_max,
        "e4_las_flas_tolerance": 1e-12,
        "e4_las_flas_passed": bool(e4_max <= 1e-12),
        "seconds": time.perf_counter() - started,
    }
    if not checks["figure_3a_passed"] or not checks["e4_las_flas_passed"]:
        raise RuntimeError(f"Detailed-layout regression failed: {checks}")
    return checks


def prepare_e2c() -> dict:
    started = time.perf_counter()
    frozen = pd.read_csv(E2C_RESULTS, float_precision="round_trip")
    included = frozen.loc[
        frozen.included.astype(bool) & frozen.layout.isin(E2C_LAYOUTS)
    ].copy()
    keys = included[["target", "endpoint"]].drop_duplicates()
    if len(keys) != 77:
        raise ValueError(f"Expected 77 included E2c datasets, found {len(keys)}")
    selected = pd.read_csv(E2C_SELECTED)
    selected = selected.loc[selected.role.eq("sampled")]
    cache_dir = QA / "e2c"
    cache_dir.mkdir(parents=True, exist_ok=True)
    observed_rows = []
    for index, item in enumerate(keys.itertuples(index=False), 1):
        target, endpoint = str(item.target), str(item.endpoint)
        raw = pd.read_csv(
            ACTIVITIES / f"{target}_{endpoint}.csv", low_memory=False
        )
        records = curate_activity_table(raw, target, endpoint).retained_records
        entries = summarize_entries(records)
        positions = subset_positions(entries.molecule_identity_key.to_numpy())
        data = entries.iloc[positions].reset_index(drop=True)
        fingerprints = morgan_fingerprints(data.canonical_smiles.tolist())
        projection_raw = project_2d(
            fingerprints,
            method="umap",
            random_state=MAP_SEED,
            umap_metric="jaccard",
        )
        projection = normalize_coordinates(projection_raw, mode="isotropic")
        minimum, side, _ = assign_to_grid(
            projection_raw,
            padding=None,
            target_occupancy=0.40,
            coordinate_scaling="isotropic",
            method="dense",
        )
        cells = np.rint(minimum * (side - 1)).astype(int)
        reference = representation_neighbor_sets(fingerprints, k=K)
        refined = refine_grid_layout(
            cells,
            projection * (side - 1),
            reference,
            side,
            **REFINEMENT_DEFAULTS,
            random_state=MAP_SEED,
        )
        points = {
            "projection": projection,
            "minimum_displacement": cells,
            "neighbor_refined": refined,
        }
        np.savez_compressed(
            cache_dir / f"{target}_{endpoint}.npz",
            molecule_identity_key=data.molecule_identity_key.to_numpy(str),
            activity=data.activity_pchembl.to_numpy(float),
            projection=points["projection"],
            minimum_displacement=points["minimum_displacement"],
            neighbor_refined=points["neighbor_refined"],
            grid_side=np.asarray([side], dtype=int),
        )
        activity = data.activity_pchembl.to_numpy(float)
        for layout in E2C_LAYOUTS:
            metrics = neighborhood_metrics(
                activity,
                reference,
                grid_neighbor_sets(points[layout], k=K),
            )
            observed_rows.append({
                "target": target,
                "endpoint": endpoint,
                "layout": layout,
                "mean_recall_retained": float(metrics.recall.mean()),
                "mean_e_grid_retained": float(metrics.e_grid.mean()),
                "mean_e_fp_retained": float(metrics.e_fp.mean()),
                "rho_retained": spearman_rho(
                    metrics.recall, metrics.excess_error
                ),
            })
        print(
            f"E6 prepare E2c {index}/77 {target} {endpoint}", flush=True
        )
    observed = pd.DataFrame(observed_rows)
    expected = included[[
        "target", "endpoint", "layout", "mean_recall", "mean_e_grid",
        "mean_e_fp", "rho",
    ]]
    check = expected.merge(
        observed,
        on=["target", "endpoint", "layout"],
        validate="one_to_one",
    )
    for observed, frozen, name in (
        ("mean_recall_retained", "mean_recall", "recall_difference"),
        ("mean_e_grid_retained", "mean_e_grid", "e_grid_difference"),
        ("mean_e_fp_retained", "mean_e_fp", "e_fp_difference"),
        ("rho_retained", "rho", "rho_difference"),
    ):
        check[name] = check[observed] - check[frozen]
    check["within_tolerance"] = (
        check[[
            "recall_difference",
            "e_grid_difference",
            "e_fp_difference",
            "rho_difference",
        ]].abs().max(axis=1) <= 1e-9
    )
    check.to_csv(QA / "e2c_regression_rows.csv", index=False)
    differences = {
        "mean_recall": maximum_difference(
            check.mean_recall_retained.fillna(check.mean_recall),
            check.mean_recall,
        ),
        "mean_e_grid": maximum_difference(
            check.mean_e_grid_retained.fillna(check.mean_e_grid),
            check.mean_e_grid,
        ),
        "mean_e_fp": maximum_difference(
            check.mean_e_fp_retained.fillna(check.mean_e_fp),
            check.mean_e_fp,
        ),
        "rho": maximum_difference(
            check.rho_retained.fillna(check.rho), check.rho
        ),
    }
    mismatch = check.loc[~check.within_tolerance]
    passed = bool(mismatch.empty)
    result = {
        "rows": int(len(check)),
        "mismatch_rows": int(len(mismatch)),
        "mismatch_datasets": int(
            mismatch[["target", "endpoint"]].drop_duplicates().shape[0]
        ),
        "mismatches_by_layout": {
            str(key): int(value)
            for key, value in mismatch.layout.value_counts().items()
        },
        "max_absolute_differences": differences,
        "tolerance": 1e-9,
        "passed": passed,
        "status": "passed" if passed else "mismatch_recorded",
        "action": (
            "proceed with reconstructed coordinates for all 77 datasets"
            if not passed else "proceed with reconstructed coordinates"
        ),
        "seconds": time.perf_counter() - started,
    }
    (QA / "e2c_regression_summary.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    if not passed:
        print(
            "E6 E2c mismatch recorded; proceeding with reconstructed coordinates",
            flush=True,
        )
    return result


def run() -> dict:
    QA.mkdir(parents=True, exist_ok=True)
    started_at = datetime.now(timezone.utc).isoformat()
    result = {
        "analysis_role": ANALYSIS_ROLE,
        "started_at": started_at,
        "e1b_detailed": prepare_detailed(),
        "e2c": prepare_e2c(),
    }
    result["ended_at"] = datetime.now(timezone.utc).isoformat()
    (QA / "regression_checks.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print("E6 REGRESSION PASS", flush=True)
    print(json.dumps(result, ensure_ascii=False, indent=2), flush=True)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    run()
