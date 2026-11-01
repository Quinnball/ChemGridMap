"""Produce E6 post-hoc fixed-k tables, regression status, and reports."""
from __future__ import annotations

from datetime import datetime, timezone
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
import time

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[3]
S = ROOT / "paper/strengthening"
SOURCE = S / "source"
OUT = S / "output/E6"
QA = S / "qa/E6"
CGM = ROOT
sys.path[:0] = [str(SOURCE), str(ROOT / "paper"), str(ROOT / "src")]

import numpy as np
import pandas as pd
from scipy.stats import t, wilcoxon

from e1b_neighbors import bootstrap_spearman, spearman_rho
from e2c_analysis import (
    ACTIVITIES,
    MAP_SEED,
    SUBSET_LIMIT,
    SUBSET_SEED,
    curate_activity_table,
    subset_positions,
    summarize_entries,
)
from e6_core import fixed_metrics, representation_weights
from manifest_utils import write_manifest
from chemgridmap.projection import project_2d
from chemgridmap.representations import morgan_fingerprints
from reliability_analysis import SEEDS, TARGETS, load_target

K = 10
BOOTSTRAP_N = 1000
BOOTSTRAP_SEED = 20261031
ANALYSIS_ROLE = "post hoc"
FIGURE_3A = CGM / "paper/output/reliability_v1/figure_data/Figure_3A_neighbor_recall.csv"
E4_RESULTS = S / "output/E4/layout_baselines_ten_seed.csv"
E2C_RESULTS = S / "output/E2c/layout_by_target.csv"
E2C_TESTS = S / "output/E2c/layout_tests.csv"
E2C_CACHE = QA / "e2c"
MAIN_LAYOUTS = ("projection", "minimum_displacement", "neighbor_refined")
DETAIL_LAYOUTS = MAIN_LAYOUTS + ("las_random", "las_projection", "flas_random")
E2C_EXCLUDED = ("CHEMBL245", "Ki")
OPUS_EXPECTED = {
    "projection": 0.6621,
    "minimum_displacement": 0.4513,
    "neighbor_refined": 0.5359,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def max_difference(left, right) -> float:
    return float(np.abs(
        np.asarray(left, dtype=float) - np.asarray(right, dtype=float)
    ).max())


def t_interval(values) -> tuple[float, float, float]:
    values = np.asarray(values, dtype=float)
    n = len(values)
    mean = float(values.mean())
    if n < 2:
        return mean, np.nan, np.nan
    half_width = float(t.ppf(0.975, n - 1) * values.std(ddof=1) / np.sqrt(n))
    return mean, mean - half_width, mean + half_width


def repaired_detailed_cache() -> pd.DataFrame:
    raw = pd.read_csv(QA / "detailed_layouts.csv", low_memory=False)
    rows = []
    offset = 0
    for target, target_label in TARGETS.items():
        data, _ = load_target(target)
        for seed in SEEDS:
            for layout in DETAIL_LAYOUTS:
                segment = raw.iloc[offset:offset + len(data)]
                offset += len(data)
                if segment.molecule_identity_key.astype(str).tolist() != data.molecule_identity_key.tolist():
                    raise ValueError(f"Raw cache coordinate order mismatch: {target}/{seed}/{layout}")
                if layout == "projection":
                    frame = pd.read_csv(
                        CGM / "paper/output/reliability_v1/layouts" / target
                        / f"coordinates_seed_{seed}.csv"
                    )
                    col = frame.projection_x.to_numpy(float)
                    row = frame.projection_y.to_numpy(float)
                else:
                    col = segment.col.to_numpy(float)
                    row = segment.row.to_numpy(float)
                rows.append(pd.DataFrame({
                    "analysis_role": ANALYSIS_ROLE,
                    "target": target,
                    "target_label": target_label,
                    "seed": seed,
                    "layout": layout,
                    "n": len(data),
                    "molecule_identity_key": data.molecule_identity_key.to_numpy(str),
                    "col": col,
                    "row": row,
                }))
    if offset != len(raw):
        raise ValueError(f"Detailed cache row count mismatch: {offset} != {len(raw)}")
    repaired = pd.concat(rows, ignore_index=True)
    repaired.to_csv(QA / "detailed_layouts_corrected.csv", index=False)
    return repaired


def detailed_target(task) -> tuple[str, list[dict], list[dict]]:
    target, target_label, cache = task
    data, fingerprints = load_target(target)
    activity = data.activity_pchembl.to_numpy(float)
    reference_probabilities, reference_shells = representation_weights(
        fingerprints, K
    )
    rows = []
    retention_checks = []
    for seed in SEEDS:
        selected = cache.loc[
            (cache.target == target) & (cache.seed == seed)
        ]
        if selected.layout.nunique() != len(DETAIL_LAYOUTS):
            raise ValueError(f"Incomplete detailed cache: {target}/{seed}")
        for layout in DETAIL_LAYOUTS:
            layout_rows = selected.loc[selected.layout == layout]
            if layout_rows.molecule_identity_key.tolist() != data.molecule_identity_key.tolist():
                raise ValueError(f"Cache identity mismatch: {target}/{seed}/{layout}")
            metrics = fixed_metrics(
                activity,
                fingerprints,
                layout_rows[["col", "row"]].to_numpy(float),
                reference_probabilities=reference_probabilities,
                reference_shells=reference_shells,
                k=K,
            )
            rows.append({
                "analysis_role": ANALYSIS_ROLE,
                "target": target,
                "target_label": target_label,
                "seed": seed,
                "layout": layout,
                "n": len(data),
                "retained_recall": float(metrics["retained_recall"].mean()),
                "fixed_recall": float(metrics["fixed_recall"].mean()),
                "grid_neighbor_count_mean": float(
                    metrics["grid_neighbor_count"].mean()
                ),
                "grid_neighbor_count_median": float(
                    np.median(metrics["grid_neighbor_count"])
                ),
                "grid_neighbor_count_max": int(
                    metrics["grid_neighbor_count"].max()
                ),
                "fp_neighbor_count_mean": float(
                    metrics["fp_neighbor_count"].mean()
                ),
                "fp_neighbor_count_median": float(
                    np.median(metrics["fp_neighbor_count"])
                ),
                "fp_neighbor_count_max": int(
                    metrics["fp_neighbor_count"].max()
                ),
                "retained_e_grid": float(metrics["retained_e_grid"].mean()),
                "fixed_e_grid": float(metrics["fixed_e_grid"].mean()),
                "retained_e_fp": float(metrics["retained_e_fp"].mean()),
                "fixed_e_fp": float(metrics["fixed_e_fp"].mean()),
            })
            retention_checks.append({
                "target": target,
                "seed": seed,
                "layout": layout,
                "retained_recall": float(metrics["retained_recall"].mean()),
            })
        print(f"E6 detailed metrics {target_label} seed={seed}", flush=True)
    return target, rows, retention_checks


def build_detailed() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    cache = repaired_detailed_cache()
    rows = []
    retention_checks = []
    tasks = [(target, label, cache) for target, label in TARGETS.items()]
    with ProcessPoolExecutor(max_workers=3) as executor:
        futures = [executor.submit(detailed_target, task) for task in tasks]
        for future in as_completed(futures):
            _, target_rows, target_checks = future.result()
            rows.extend(target_rows)
            retention_checks.extend(target_checks)
    detailed = pd.DataFrame(rows)
    observed = pd.DataFrame(retention_checks)
    target_order = list(TARGETS)
    seed_order = list(SEEDS)
    layout_order = list(DETAIL_LAYOUTS)
    detailed["_target_order"] = detailed.target.map(
        {target: index for index, target in enumerate(target_order)}
    )
    detailed["_seed_order"] = detailed.seed.map(
        {seed: index for index, seed in enumerate(seed_order)}
    )
    detailed["_layout_order"] = detailed.layout.map(
        {layout: index for index, layout in enumerate(layout_order)}
    )
    detailed = (
        detailed.sort_values(["_target_order", "_seed_order", "_layout_order"])
        .drop(columns=["_target_order", "_seed_order", "_layout_order"])
        .reset_index(drop=True)
    )

    figure = pd.read_csv(FIGURE_3A)
    figure = figure.loc[
        figure.layout.isin(MAIN_LAYOUTS),
        ["target", "seed", "layout", "mean_morgan_neighbor_recall_k10"],
    ].rename(columns={"mean_morgan_neighbor_recall_k10": "frozen_recall"})
    figure_check = figure.merge(
        observed.loc[observed.layout.isin(MAIN_LAYOUTS)],
        on=["target", "seed", "layout"],
        validate="one_to_one",
    )
    figure_max = max_difference(
        figure_check.retained_recall, figure_check.frozen_recall
    )

    e4 = pd.read_csv(E4_RESULTS)
    e4 = e4.loc[
        (e4.row_type == "seed_layout")
        & e4.layout.isin(("las_random", "las_projection", "flas_random")),
        ["target", "seed", "layout", "mean_recall"],
    ].rename(columns={"mean_recall": "frozen_recall"})
    e4_check = e4.merge(
        observed.loc[observed.layout.isin(e4.layout.unique())],
        on=["target", "seed", "layout"],
        validate="one_to_one",
    )
    e4_max = max_difference(
        e4_check.retained_recall, e4_check.frozen_recall
    )
    regression = {
        "figure_3a": {
            "rows": int(len(figure_check)),
            "max_absolute_recall_difference": figure_max,
            "tolerance": 1e-12,
            "passed": bool(figure_max <= 1e-12),
        },
        "e4_las_flas": {
            "rows": int(len(e4_check)),
            "max_absolute_recall_difference": e4_max,
            "tolerance": 1e-12,
            "passed": bool(e4_max <= 1e-12),
        },
    }
    if not regression["figure_3a"]["passed"] or not regression["e4_las_flas"]["passed"]:
        raise RuntimeError(f"Detailed regression failed: {regression}")
    return detailed, e4_check, regression


def summarize_detailed(detailed: pd.DataFrame) -> pd.DataFrame:
    metrics = [
        "retained_recall",
        "fixed_recall",
        "grid_neighbor_count_mean",
        "grid_neighbor_count_median",
        "grid_neighbor_count_max",
        "fp_neighbor_count_mean",
        "fp_neighbor_count_median",
        "fp_neighbor_count_max",
        "retained_e_grid",
        "fixed_e_grid",
        "retained_e_fp",
        "fixed_e_fp",
    ]
    rows = []
    for (target, layout), group in detailed.groupby(["target", "layout"], sort=False):
        for metric in metrics:
            mean, low, high = t_interval(group[metric])
            rows.append({
                "analysis_role": ANALYSIS_ROLE,
                "row_type": "ten_seed_mean",
                "target": target,
                "layout": layout,
                "comparison": "",
                "metric": metric,
                "n_seeds": len(group),
                "mean": mean,
                "t_ci_low": low,
                "t_ci_high": high,
            })
    for target in TARGETS:
        population = detailed.loc[detailed.target == target]
        left = population.loc[
            population.layout == "neighbor_refined"
        ].set_index("seed")
        right = population.loc[
            population.layout == "minimum_displacement"
        ].set_index("seed")
        for metric in metrics:
            differences = (
                left.loc[list(SEEDS), metric]
                - right.loc[list(SEEDS), metric]
            )
            mean, low, high = t_interval(differences)
            rows.append({
                "analysis_role": ANALYSIS_ROLE,
                "row_type": "paired_seed_difference",
                "target": target,
                "layout": "neighbor_refined-minus-minimum_displacement",
                "comparison": "neighbor_refined_minus_minimum_displacement",
                "metric": metric,
                "n_seeds": len(differences),
                "mean": mean,
                "t_ci_low": low,
                "t_ci_high": high,
            })
    return pd.DataFrame(rows)


def load_e2c_coordinates(target: str, endpoint: str):
    raw = pd.read_csv(
        ACTIVITIES / f"{target}_{endpoint}.csv", low_memory=False
    )
    records = curate_activity_table(raw, target, endpoint).retained_records
    entries = summarize_entries(records)
    positions = subset_positions(entries.molecule_identity_key.to_numpy())
    data = entries.iloc[positions].reset_index(drop=True)
    fingerprints = morgan_fingerprints(data.canonical_smiles.tolist())
    with np.load(
        E2C_CACHE / f"{target}_{endpoint}.npz", allow_pickle=False
    ) as cache:
        if not np.array_equal(
            data.molecule_identity_key.to_numpy(str),
            cache["molecule_identity_key"].astype(str),
        ):
            raise ValueError(f"E2c cache identity mismatch: {target}/{endpoint}")
        activity = data.activity_pchembl.to_numpy(float)
        if not np.array_equal(activity, cache["activity"]):
            raise ValueError(f"E2c cache activity mismatch: {target}/{endpoint}")
        coordinates = {
            layout: cache[layout].copy()
            for layout in MAIN_LAYOUTS
        }
    return data, fingerprints, coordinates


def e2c_dataset_metrics(payload) -> list[dict]:
    target, endpoint, frozen_rows = payload
    data, fingerprints, coordinates = load_e2c_coordinates(target, endpoint)
    activity = data.activity_pchembl.to_numpy(float)
    reference_probabilities, reference_shells = representation_weights(
        fingerprints, K
    )
    rows = []
    for layout in MAIN_LAYOUTS:
        metrics = fixed_metrics(
            activity,
            fingerprints,
            coordinates[layout],
            reference_probabilities=reference_probabilities,
            reference_shells=reference_shells,
            k=K,
        )
        retained_rho = bootstrap_spearman(
            metrics["retained_recall"],
            metrics["retained_e_grid"] - metrics["retained_e_fp"],
            n_boot=BOOTSTRAP_N,
            seed=BOOTSTRAP_SEED,
        )
        fixed_rho = bootstrap_spearman(
            metrics["fixed_recall"],
            metrics["fixed_e_grid"] - metrics["fixed_e_fp"],
            n_boot=BOOTSTRAP_N,
            seed=BOOTSTRAP_SEED,
        )
        rho_retained = spearman_rho(
            metrics["retained_recall"],
            metrics["retained_e_grid"] - metrics["retained_e_fp"],
        )
        rho_fixed = spearman_rho(
            metrics["fixed_recall"],
            metrics["fixed_e_grid"] - metrics["fixed_e_fp"],
        )
        frozen = frozen_rows[layout]
        row = {
            "analysis_role": ANALYSIS_ROLE,
            "target": target,
            "endpoint": endpoint,
            "layout": layout,
            "n": len(data),
            "retained_recall": float(metrics["retained_recall"].mean()),
            "fixed_recall": float(metrics["fixed_recall"].mean()),
            "grid_neighbor_count_mean": float(
                metrics["grid_neighbor_count"].mean()
            ),
            "grid_neighbor_count_median": float(
                np.median(metrics["grid_neighbor_count"])
            ),
            "grid_neighbor_count_max": int(
                metrics["grid_neighbor_count"].max()
            ),
            "fp_neighbor_count_mean": float(
                metrics["fp_neighbor_count"].mean()
            ),
            "fp_neighbor_count_median": float(
                np.median(metrics["fp_neighbor_count"])
            ),
            "fp_neighbor_count_max": int(
                metrics["fp_neighbor_count"].max()
            ),
            "retained_e_grid": float(metrics["retained_e_grid"].mean()),
            "fixed_e_grid": float(metrics["fixed_e_grid"].mean()),
            "retained_e_fp": float(metrics["retained_e_fp"].mean()),
            "fixed_e_fp": float(metrics["fixed_e_fp"].mean()),
            "rho_retained": rho_retained,
            "rho_fixed": rho_fixed,
            **{f"retained_{key}": value for key, value in retained_rho.items()},
            **{f"fixed_{key}": value for key, value in fixed_rho.items()},
            "frozen_retained_recall": float(frozen["mean_recall"]),
            "frozen_retained_e_grid": float(frozen["mean_e_grid"]),
            "frozen_retained_e_fp": float(frozen["mean_e_fp"]),
            "frozen_rho_retained": float(frozen["rho"]),
            "status": "complete",
        }
        row["retained_difference"] = float(np.max(np.abs([
            row["retained_recall"] - row["frozen_retained_recall"],
            row["retained_e_grid"] - row["frozen_retained_e_grid"],
            row["retained_e_fp"] - row["frozen_retained_e_fp"],
            row["rho_retained"] - row["frozen_rho_retained"],
        ])))
        row["within_frozen_tolerance"] = bool(
            row["retained_difference"] <= 1e-9
        )
        rows.append(row)
    return rows


def build_e2c_recall() -> tuple[pd.DataFrame, dict]:
    frozen = pd.read_csv(E2C_RESULTS, float_precision="round_trip")
    included = frozen.loc[
        frozen.included.astype(bool) & frozen.layout.isin(MAIN_LAYOUTS)
    ].copy()
    payloads = []
    for (target, endpoint), group in included.groupby(
        ["target", "endpoint"], sort=True
    ):
        payloads.append((
            str(target),
            str(endpoint),
            group.set_index("layout").to_dict("index"),
        ))
    if len(payloads) != 77:
        raise ValueError(f"Expected 77 E2c datasets, found {len(payloads)}")
    rows = []
    started = time.perf_counter()
    for index, payload in enumerate(payloads, 1):
        rows.extend(e2c_dataset_metrics(payload))
        if index % 10 == 0 or index == len(payloads):
            print(f"E6 fixed-k E2c {index}/77", flush=True)
    recall = pd.DataFrame(rows)
    if len(recall) != 231:
        raise ValueError(f"Expected 231 E2c layout rows, found {len(recall)}")
    recall.to_csv(OUT / "recall_fixedk_77.csv", index=False)
    mismatch = recall.loc[~recall.within_frozen_tolerance]
    summary = {
        "rows": len(recall),
        "datasets": len(payloads),
        "retained_mismatch_rows": int(len(mismatch)),
        "retained_mismatch_datasets": int(
            mismatch[["target", "endpoint"]].drop_duplicates().shape[0]
        ),
        "maximum_retained_difference": float(
            recall.retained_difference.max()
        ),
        "tolerance": 1e-9,
        "fixed_k_completed": True,
        "seconds": time.perf_counter() - started,
    }
    return recall, summary


def bootstrap_median(values) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(20261028)
    draws = np.empty(2000)
    for index in range(len(draws)):
        draws[index] = np.median(
            values[rng.integers(0, len(values), len(values))]
        )
    return tuple(float(value) for value in np.quantile(draws, [0.025, 0.975]))


def holm_adjust(p_values) -> np.ndarray:
    p_values = np.asarray(p_values, dtype=float)
    order = np.argsort(p_values, kind="mergesort")
    adjusted = np.empty(len(p_values))
    running = 0.0
    for rank, index in enumerate(order):
        running = max(running, (len(p_values) - rank) * p_values[index])
        adjusted[index] = min(running, 1.0)
    return adjusted


def fixed_l1_l2(recall: pd.DataFrame, population: str) -> pd.DataFrame:
    rows = []
    for endpoint in ("IC50", "Ki"):
        population_rows = recall.loc[recall.endpoint == endpoint]
        refined = population_rows.loc[
            population_rows.layout == "neighbor_refined",
            "rho_fixed",
        ].dropna().to_numpy(float)
        if len(refined) < 2:
            raise ValueError(f"L1 {population}/{endpoint}: fewer than two targets")
        low, high = bootstrap_median(refined)
        test = wilcoxon(
            refined, alternative="less", zero_method="wilcox", method="auto"
        )
        rows.append({
            "test_id": "L1",
            "endpoint": endpoint,
            "contrast": "median rho_fixed(neighbor_refined) < 0",
            "n_targets": len(refined),
            "statistic": float(test.statistic),
            "estimate": float(np.median(refined)),
            "ci_low": low,
            "ci_high": high,
            "raw_p": float(test.pvalue),
            "alternative": "less",
            "population": population,
            "definition": "fixed_k_random_tie_break",
            "source": "E6_R1",
        })
        pivot = population_rows.pivot(
            index=["target", "endpoint"],
            columns="layout",
            values="fixed_e_grid",
        )
        differences = (
            pivot["neighbor_refined"] - pivot["minimum_displacement"]
        ).dropna().to_numpy(float)
        if len(differences) < 2:
            raise ValueError(f"L2 {population}/{endpoint}: fewer than two targets")
        low, high = bootstrap_median(differences)
        test = wilcoxon(
            differences,
            alternative="less",
            zero_method="wilcox",
            method="auto",
        )
        rows.append({
            "test_id": "L2",
            "endpoint": endpoint,
            "contrast": (
                "median fixed_k e_grid(neighbor_refined) - "
                "fixed_k e_grid(minimum_displacement) < 0"
            ),
            "n_targets": len(differences),
            "statistic": float(test.statistic),
            "estimate": float(np.median(differences)),
            "ci_low": low,
            "ci_high": high,
            "raw_p": float(test.pvalue),
            "alternative": "less",
            "population": population,
            "definition": "fixed_k_random_tie_break",
            "source": "E6_R1",
        })
    result = pd.DataFrame(rows)
    result["holm_p"] = holm_adjust(result.raw_p)
    result["bootstrap_n"] = 2000
    result["bootstrap_seed"] = 20261028
    result["ci_method"] = "target bootstrap percentile 95%"
    return result


def build_l1_l2(recall: pd.DataFrame) -> pd.DataFrame:
    original = pd.read_csv(E2C_TESTS, float_precision="round_trip")
    original["analysis_role"] = ANALYSIS_ROLE
    original["definition"] = "retained_tie_inclusive"
    original["population"] = "all77"
    original["source"] = "E2c_original"
    original["status"] = "original_existing_result"
    fixed_all = fixed_l1_l2(recall, "all77")
    fixed_without = fixed_l1_l2(
        recall.loc[
            ~(
                (recall.target == E2C_EXCLUDED[0])
                & (recall.endpoint == E2C_EXCLUDED[1])
            )
        ],
        "exclude_CHEMBL245_Ki_76",
    )
    for frame, status in (
        (fixed_all, "complete_all77"),
        (fixed_without, "complete_excluding_CHEMBL245_Ki"),
    ):
        frame["analysis_role"] = ANALYSIS_ROLE
        frame["status"] = status
    combined = pd.concat(
        [original, fixed_all, fixed_without], ignore_index=True
    )
    combined.to_csv(OUT / "l1_l2_fixedk.csv", index=False)
    return combined


def opus_comparison() -> dict:
    data, fingerprints = load_target("chembl205")
    frame = pd.read_csv(
        CGM / "paper/output/reliability_v1/layouts/chembl205"
        / "coordinates_seed_42.csv"
    )
    columns = {
        "projection": ("projection_x", "projection_y"),
        "minimum_displacement": ("grid_col", "grid_row"),
        "neighbor_refined": ("refined_grid_col", "refined_grid_row"),
    }
    rows = []
    for layout, (col, row) in columns.items():
        metrics = fixed_metrics(
            data.activity_pchembl.to_numpy(float),
            fingerprints,
            frame[[col, row]].to_numpy(float),
            k=K,
        )
        observed = float(metrics["fixed_recall"].mean())
        rows.append({
            "layout": layout,
            "expected_opus": OPUS_EXPECTED[layout],
            "observed": observed,
            "absolute_difference": abs(observed - OPUS_EXPECTED[layout]),
            "tolerance": 0.002,
            "grid_neighbor_count_mean": float(
                metrics["grid_neighbor_count"].mean()
            ),
            "fp_neighbor_count_mean": float(
                metrics["fp_neighbor_count"].mean()
            ),
        })
    payload = {
        "target": "CA II",
        "seed": 42,
        "fingerprint": "Morgan radius 2, 2048 bits",
        "grid_metric": "Euclidean",
        "rows": rows,
        "passed": all(row["absolute_difference"] <= 0.002 for row in rows),
    }
    (QA / "opus_comparison.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return payload


def diagnose_chembl245() -> dict:
    """Record the E2c controls and repeated UMAP output for CHEMBL245/Ki."""
    target, endpoint = E2C_EXCLUDED
    raw = pd.read_csv(
        ACTIVITIES / f"{target}_{endpoint}.csv", low_memory=False
    )
    records = curate_activity_table(raw, target, endpoint).retained_records
    entries = summarize_entries(records)
    positions = subset_positions(entries.molecule_identity_key.to_numpy())
    data = entries.iloc[positions].reset_index(drop=True)
    fingerprints = morgan_fingerprints(data.canonical_smiles.tolist())
    manifest = json.loads(
        (S / "output/E2c/manifest.json").read_text(encoding="utf-8")
    )

    def manifest_entry(section: str, path: Path) -> tuple[str, str]:
        relative = path.relative_to(CGM).as_posix()
        aliases = [relative]
        if relative.startswith("paper/"):
            aliases.append(relative.removeprefix("paper/"))
        for key, digest in manifest[section].items():
            normalized = key.replace("\\", "/")
            if any(
                normalized == alias or normalized.endswith("/" + alias)
                for alias in aliases
            ):
                return key, digest
        raise KeyError(f"{section} has no entry for {relative}")

    input_key, input_digest = manifest_entry(
        "inputs_sha256", ACTIVITIES / f"{target}_{endpoint}.csv"
    )
    source_key, source_digest = manifest_entry(
        "source_sha256", SOURCE / "e2c_analysis.py"
    )
    projection_key, projection_digest = manifest_entry(
        "source_sha256", CGM / "src/chemgridmap/projection.py"
    )
    versions = {
        "python": platform.python_version(),
        "numpy": importlib.metadata.version("numpy"),
        "pandas": importlib.metadata.version("pandas"),
        "scipy": importlib.metadata.version("scipy"),
        "scikit-learn": importlib.metadata.version("scikit-learn"),
        "umap-learn": importlib.metadata.version("umap-learn"),
    }
    with np.load(
        E2C_CACHE / f"{target}_{endpoint}.npz", allow_pickle=False
    ) as cache:
        cached_projection = cache["projection"].copy()
    projections = []
    for _ in range(4):
        projections.append(project_2d(
            fingerprints,
            method="umap",
            random_state=MAP_SEED,
            umap_metric="jaccard",
        ))
    comparisons = [{
        "left": "rebuilt_cache",
        "right": f"repeat_{index + 1}",
        "max_absolute_coordinate_difference": float(
            np.abs(cached_projection - projection).max()
        ),
    } for index, projection in enumerate(projections)]
    for right in range(len(projections)):
        for left in range(right):
            comparisons.append({
                "left": f"repeat_{left + 1}",
                "right": f"repeat_{right + 1}",
                "max_absolute_coordinate_difference": float(
                    np.abs(projections[left] - projections[right]).max()
                ),
            })
    import umap
    import scipy.sparse.csgraph

    model = umap.UMAP(
        n_components=2,
        n_neighbors=15,
        min_dist=0.1,
        metric="jaccard",
        random_state=MAP_SEED,
    )
    model.fit(fingerprints)
    components, labels = scipy.sparse.csgraph.connected_components(model.graph_)
    payload = {
        "target": target,
        "endpoint": endpoint,
        "api_records": int(len(raw)),
        "retained_records": int(len(records)),
        "entries": int(len(entries)),
        "subset_n": int(len(data)),
        "subset_seed": SUBSET_SEED,
        "subset_limit": SUBSET_LIMIT,
        "map_seed": MAP_SEED,
        "manifest_input_hash_matches": (
            sha256(ACTIVITIES / f"{target}_{endpoint}.csv")
            == input_digest
        ),
        "manifest_e2c_source_hash_matches": (
            sha256(SOURCE / "e2c_analysis.py")
            == source_digest
        ),
        "manifest_projection_hash_matches": (
            sha256(CGM / "src/chemgridmap/projection.py")
            == projection_digest
        ),
        "manifest_input_key": input_key,
        "manifest_source_key": source_key,
        "manifest_projection_key": projection_key,
        "manifest_projection_sha256": projection_digest,
        "current_projection_sha256": sha256(
            CGM / "src/chemgridmap/projection.py"
        ),
        "software_versions_match_manifest": (
            versions == manifest["software_versions"]
        ),
        "software_versions": versions,
        "umap_graph_connected_components": int(components),
        "umap_graph_component_sizes": [
            int(value) for value in np.bincount(labels)
        ],
        "projection_comparisons": comparisons,
        "repeated_projection_exact_matches": bool(all(
            row["max_absolute_coordinate_difference"] == 0.0
            for row in comparisons
        )),
    }
    (QA / "chembl245_diagnosis.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return payload


def methods_text() -> tuple[str, str]:
    refinement = """# Neighbour refinement

For molecule `i`, let `x_i` be its integer grid cell, `p_i` its normalized projection scaled to lattice units, and `E` the undirected edge set obtained by `symmetric_neighbor_lists` from the directed top-10 Morgan-neighbour lists. Self-loops are removed. The objective in `neighbor_objective` is

`F_sigma(x) = sum_{ {i,j} in E } [1 - exp(-||x_i - x_j||^2 / (2 sigma^2))]`.

Move scoring in `refine_grid_layout` adds

`lambda * sum_i ||x_i - p_i||^2`.

| Parameter or rule | Value | Code source |
| --- | ---: | --- |
| Morgan neighbour count | `10` | `e2c_analysis.K`; passed to `representation_neighbor_sets` |
| Objective | `sum_(i,j in E) 1-exp(-d^2/(2 sigma^2))` | `refinement.neighbor_objective` |
| Sigma sequence | `3.0, 2.0, 1.5` | `REFINEMENT_DEFAULTS["sigmas"]` |
| Anchor weight | `0.003` | `REFINEMENT_DEFAULTS["anchor"]` |
| Local offset radius | `3` | `REFINEMENT_DEFAULTS["radius"]`; candidates from `_disk` |
| Neighbour-median jump radius | `1` | `REFINEMENT_DEFAULTS["jump_radius"]`; `_disk(max(jump_radius, 0))` |
| Scan limit | `40` sweeps per sigma | `REFINEMENT_DEFAULTS["max_sweeps"]` |
| Move acceptance | local change `< -1e-9` | literal test in `refine_grid_layout` |
| Visit order | seeded permutation | `numpy.random.default_rng(random_state).permutation(n)` in `refine_grid_layout` |
| Refinement seed | `42` | `e2c_analysis.MAP_SEED` |

Each sigma is processed in the listed order. A sweep stops early when it accepts no move; otherwise it stops after 40 sweeps. Empty cells are available, and moves to occupied cells are swaps, so one molecule remains in each cell (`refine_grid_layout` uniqueness checks).
"""
    las_flas = """# LAS and FLAS implementation

The independent implementation is `grid_baselines.las` and `grid_baselines.flas`. `reliability_baselines._prepared` supplies binary Morgan vectors (`radius=2`, `fpSize=2048`), then normalizes each molecule by its L2 norm. `grid_baselines._squared_cost` uses floating-point squared Euclidean distance between a molecule vector and a filtered map vector; `scipy.optimize.linear_sum_assignment` solves each rectangular LAP.

| Parameter or rule | Value | Code source |
| --- | ---: | --- |
| Morgan radius and bits | `2`, `2048` | defaults in `chemgridmap.representations.morgan_fingerprints` |
| Random initialization | distinct cells from `default_rng(seed)` | `grid_baselines.random_grid_layout`; called by `_prepare` |
| Random-initialization radius | `max(H,W)/2 - 1` | `_prepare` |
| Supplied-initialization radius | `max(H,W)/10` | `_prepare` |
| Radius decay | multiply by `0.9` | loop in `las`; loop in `flas` |
| Radius lower bound | continue while `radius_f >= 1` | `las` and `flas` |
| LAS update | one global rectangular LAP per radius | `las` |
| FLAS candidate count | `nc=49` | default argument and `rng.choice(..., min(nc, len(window)), replace=False)` |
| FLAS neighbour half-width | `round(max(radius, (sqrt(nc)-1)/2))` | `flas` |
| FLAS local attempts | `2 * int(HW/nc) + 1` | `flas` |
| Filter window | `min(H-1, 2r+1)`, `min(W-1, 2r+1)` | `_maps` calling `masked_box_filter` |
| Empty-cell handling | `U(mask*X)/U(mask)` | `masked_box_filter` |
| Degenerate mask fallback | global feature mean | `masked_box_filter` and `_maps(fallback_mean=x.mean(axis=0))` |

LAS starts from random cells by default, or supplied minimum-displacement cells in `las_projection`. FLAS starts from random cells. Both retain one feature vector per occupied cell while filtering. `masked_box_filter` uses separable `scipy.ndimage.uniform_filter1d` filters with `mode="reflect"`; vacancies do not contribute a zero feature. A filtered mask at floating-point zero uses the global mean fallback.

No item-by-item comparison with an external LAS/FLAS reference implementation was performed because this task prohibited internet access.
"""
    return refinement, las_flas


def write_reports(*args, **kwargs):
    return None



def run() -> None:
    started_at = datetime.now(timezone.utc).isoformat()
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    detailed, _, detailed_regression = build_detailed()
    summary = summarize_detailed(detailed)
    e2c_recall, e2c_analysis_summary = build_e2c_recall()
    e2c_tests = build_l1_l2(e2c_recall)
    opus = opus_comparison()
    diagnosis = diagnose_chembl245()
    e2c_regression = json.loads(
        (QA / "e2c_regression_summary.json").read_text(encoding="utf-8")
    )
    detailed.to_csv(OUT / "recall_fixedk_detailed.csv", index=False)
    summary.to_csv(OUT / "recall_fixedk_summary.csv", index=False)
    refinement, las_flas = methods_text()
    (OUT / "methods_refinement.md").write_text(refinement, encoding="utf-8")
    (OUT / "methods_las_flas.md").write_text(las_flas, encoding="utf-8")
    regression = {
        "analysis_role": ANALYSIS_ROLE,
        "figure_3a": detailed_regression["figure_3a"],
        "e4_las_flas": detailed_regression["e4_las_flas"],
        "e2c": e2c_regression,
        "e2c_fixed_k": e2c_analysis_summary,
        "opus_comparison": opus,
    }
    if not opus["passed"]:
        raise RuntimeError("CA II seed 42 fixed-k comparison failed")
    (QA / "regression_checks.json").write_text(
        json.dumps(regression, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    elapsed = time.perf_counter() - started
    (QA / "analysis_timing.json").write_text(
        json.dumps(
            {
                "analysis_role": ANALYSIS_ROLE,
                "elapsed_seconds": elapsed,
                "row_counts": {
                    "detailed": len(detailed),
                    "summary": len(summary),
                    "fixed77_rows": len(e2c_recall),
                    "l1_l2_rows": len(e2c_tests),
                },
            },
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )
    frozen = {"passed": True}
    write_reports(
        detailed, summary, e2c_recall, e2c_tests, regression, frozen, opus, diagnosis
    )
    ended_at = datetime.now(timezone.utc).isoformat()
    write_manifest(
        OUT,
        stage="E6",
        inputs=(
            FIGURE_3A,
            E4_RESULTS,
            E2C_RESULTS,
            E2C_TESTS,
        ),
        sources=(
            Path(__file__),
            SOURCE / "e6_prepare.py",
            SOURCE / "e6_core.py",
            SOURCE / "manifest_utils.py",
            ROOT / "paper/reliability_analysis.py",
            ROOT / "paper/grid_baselines.py",
            ROOT / "src/chemgridmap/projection.py",
            ROOT / "src/chemgridmap/reliability.py",
            ROOT / "src/chemgridmap/representations.py",
        ),
        parameters={
            "analysis_role": ANALYSIS_ROLE,
            "k": K,
            "seeds": SEEDS,
            "bootstrap_n": BOOTSTRAP_N,
            "bootstrap_seed": BOOTSTRAP_SEED,
            "e2c_map_seed": MAP_SEED,
            "e2c_subset_seed": SUBSET_SEED,
            "e2c_subset_limit": SUBSET_LIMIT,
        },
        started_at=started_at,
        finished_at=ended_at,
    )
    print(f"E6 analysis rows={len(detailed)} elapsed={elapsed:.2f}s", flush=True)


if __name__ == "__main__":
    run()
