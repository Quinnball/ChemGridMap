"""Compare pre-0.5 and 0.5 disconnected-layout results."""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import warnings

sys.dont_write_bytecode = True
os.environ.setdefault("NUMBA_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[3]
S = ROOT / "paper/strengthening"
SOURCE = S / "source"
OUT = S / "output/E10"
QA = S / "qa/E10"
AS_EXECUTED = S / "output/as_executed_pre05"
sys.path[:0] = [str(SOURCE), str(ROOT / "paper"), str(ROOT / "src")]

import numpy as np
import pandas as pd
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr

from chemgridmap.core import assign_to_grid, normalize_coordinates
from chemgridmap.projection import project_2d
from chemgridmap.refinement import REFINEMENT_DEFAULTS, refine_grid_layout
from chemgridmap.reliability import representation_neighbor_sets
from chemgridmap.representations import morgan_fingerprints
from e2c_analysis import (
    ACTIVITIES,
    MAP_SEED,
    MIN_ENTRIES,
    curate_activity_table,
    subset_positions,
    summarize_entries,
)
from e6_core import fixed_metrics, representation_weights
from grid_baselines import random_grid_layout
from manifest_utils import write_manifest


ANALYSIS_ROLE = "post hoc"
K = 10
TOLERANCE = 1e-12
MAIN_LAYOUTS = ("projection", "minimum_displacement", "neighbor_refined")
OLD_E2C_RESULTS = AS_EXECUTED / "E2c/layout_by_target.csv"
OLD_E2C_TESTS = AS_EXECUTED / "E2c/layout_tests.csv"
OLD_E6_RESULTS = AS_EXECUTED / "E6/recall_fixedk_77.csv"
OLD_E6_TESTS = AS_EXECUTED / "E6/l1_l2_fixedk.csv"
E2C_RESULTS = S / "output/E2c/layout_by_target.csv"
E2C_TESTS = S / "output/E2c/layout_tests.csv"
E6_RESULTS = S / "output/E6/recall_fixedk_77.csv"
E6_TESTS = S / "output/E6/l1_l2_fixedk.csv"
CONNECTED_CONTROLS = (
    ("CHEMBL1163125", "IC50"),
    ("CHEMBL1833", "Ki"),
    ("CHEMBL1841", "IC50"),
)
E2C_METRIC_COLUMNS = (
    "n",
    "mean_recall",
    "mean_e_grid",
    "mean_e_fp",
    "rho",
    "rho_ci_low",
    "rho_ci_high",
    "bootstrap_valid_count",
    "grid_side",
    "bootstrap_n",
    "lambda",
)
E6_METRIC_COLUMNS = (
    "retained_recall",
    "fixed_recall",
    "retained_e_grid",
    "fixed_e_grid",
    "retained_e_fp",
    "fixed_e_fp",
    "rho_retained",
    "rho_fixed",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def array_sha256(array: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def json_safe(value):
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            json_safe(value),
            indent=2,
            ensure_ascii=False,
            allow_nan=False,
        )
        + "\n",
        encoding="utf-8",
    )


def load_dataset(target: str, endpoint: str):
    raw = pd.read_csv(
        ACTIVITIES / f"{target}_{endpoint}.csv",
        low_memory=False,
    )
    records = curate_activity_table(raw, target, endpoint).retained_records
    entries = summarize_entries(records)
    if len(entries) < MIN_ENTRIES:
        raise ValueError(f"{target}/{endpoint} has fewer than {MIN_ENTRIES} entries")
    positions = subset_positions(entries.molecule_identity_key.to_numpy())
    data = entries.iloc[positions].reset_index(drop=True)
    fingerprints = morgan_fingerprints(data.canonical_smiles.tolist())
    return data, fingerprints


def projection_parameters(n_entries: int) -> dict:
    return {
        "n_components": 2,
        "n_neighbors": min(15, max(2, n_entries - 1)),
        "min_dist": 0.1,
        "metric": "jaccard",
        "random_state": MAP_SEED,
    }


def component_summary(fingerprints: np.ndarray) -> tuple[int, list[int]]:
    """Fit the UMAP graph used by project_2d and return component sizes."""
    import umap

    matrix = np.asarray(fingerprints, dtype=np.float32)
    reducer = umap.UMAP(
        init="random",
        n_epochs=1,
        **projection_parameters(len(matrix)),
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        reducer.fit(matrix)
    count, labels = connected_components(reducer.graph_)
    sizes = sorted(
        (int(value) for value in np.bincount(labels)),
        reverse=True,
    )
    return int(count), sizes


def pair_difference(left, right) -> float:
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    if left.shape != right.shape:
        return float("inf")
    if not len(left):
        return 0.0
    both_nan = np.isnan(left) & np.isnan(right)
    one_nan = np.isnan(left) ^ np.isnan(right)
    delta = np.abs(left - right)
    delta[both_nan] = 0.0
    delta[one_nan] = np.inf
    return float(np.max(delta))


def changed_e2c_datasets(
    old: pd.DataFrame,
    new: pd.DataFrame,
) -> tuple[set[tuple[str, str]], pd.DataFrame]:
    """Return datasets whose E2c metric columns changed beyond tolerance."""
    keys = {
        (str(target), str(endpoint))
        for target, endpoint in new.loc[
            new.included.astype(bool), ["target", "endpoint"]
        ].itertuples(index=False, name=None)
    }
    if len(keys) != 77:
        raise ValueError(f"Expected 77 included E2c datasets, found {len(keys)}")
    rows = []
    changed = set()
    for target, endpoint in sorted(keys):
        old_rows = old.loc[
            old.target.eq(target) & old.endpoint.eq(endpoint)
        ]
        new_rows = new.loc[
            new.target.eq(target) & new.endpoint.eq(endpoint)
        ]
        maximum = 0.0
        for layout in MAIN_LAYOUTS:
            old_layout = old_rows.loc[old_rows.layout.eq(layout)]
            new_layout = new_rows.loc[new_rows.layout.eq(layout)]
            if len(old_layout) != 1 or len(new_layout) != 1:
                maximum = float("inf")
                break
            columns = [
                column
                for column in E2C_METRIC_COLUMNS
                if column in old_layout and column in new_layout
            ]
            maximum = max(
                maximum,
                pair_difference(
                    old_layout.iloc[0][columns],
                    new_layout.iloc[0][columns],
                ),
            )
        is_changed = maximum > TOLERANCE
        if is_changed:
            changed.add((target, endpoint))
        rows.append({
            "target": target,
            "endpoint": endpoint,
            "max_absolute_metric_difference": maximum,
            "changed_beyond_1e-12": bool(is_changed),
        })
    return changed, pd.DataFrame(rows)


def layouts_from_raw(raw_projection: np.ndarray, fingerprints: np.ndarray):
    projection = normalize_coordinates(raw_projection, mode="isotropic")
    minimum, side, _ = assign_to_grid(
        raw_projection,
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
    random = random_grid_layout(len(fingerprints), (side, side), seed=MAP_SEED)
    return {
        "projection": projection,
        "minimum_displacement": cells,
        "neighbor_refined": refined,
        "random": random,
        "grid_side": int(side),
    }


def metric_row(
    target: str,
    endpoint: str,
    layout: str,
    data: pd.DataFrame,
    fingerprints: np.ndarray,
    points: np.ndarray,
    grid_side: int,
    reference_probabilities: np.ndarray,
    reference_shells,
    random_fixed_e_grid: float,
    random_retained_e_grid: float,
) -> dict:
    metrics = fixed_metrics(
        data.activity_pchembl.to_numpy(float),
        fingerprints,
        points,
        reference_probabilities=reference_probabilities,
        reference_shells=reference_shells,
        k=K,
    )
    return {
        "target": target,
        "endpoint": endpoint,
        "layout": layout,
        "n": int(len(data)),
        "grid_side": int(grid_side),
        "retained_recall": float(metrics["retained_recall"].mean()),
        "fixed_recall": float(metrics["fixed_recall"].mean()),
        "retained_e_grid": float(metrics["retained_e_grid"].mean()),
        "fixed_e_grid": float(metrics["fixed_e_grid"].mean()),
        "retained_e_fp": float(metrics["retained_e_fp"].mean()),
        "fixed_e_fp": float(metrics["fixed_e_fp"].mean()),
        "rho_retained": float(spearmanr(
            metrics["retained_recall"],
            metrics["retained_e_grid"] - metrics["retained_e_fp"],
        ).statistic),
        "rho_fixed": float(spearmanr(
            metrics["fixed_recall"],
            metrics["fixed_e_grid"] - metrics["fixed_e_fp"],
        ).statistic),
        "lambda_retained": (
            1.0 - float(metrics["retained_e_fp"].mean()) / random_retained_e_grid
        ),
        "lambda_fixed": (
            1.0 - float(metrics["fixed_e_fp"].mean()) / random_fixed_e_grid
        ),
    }


def explicit_spectral_projection(fingerprints: np.ndarray) -> np.ndarray:
    """Execute the historical spectral path without ChemGridMap fallback."""
    import umap

    matrix = np.asarray(fingerprints, dtype=np.float32)
    reducer = umap.UMAP(
        init="spectral",
        **projection_parameters(len(matrix)),
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.asarray(reducer.fit_transform(matrix), dtype=np.float32)


def connected_controls() -> tuple[pd.DataFrame, dict]:
    e2c = pd.read_csv(E2C_RESULTS, float_precision="round_trip")
    e6 = pd.read_csv(E6_RESULTS, float_precision="round_trip")
    rows = []
    for target, endpoint in CONNECTED_CONTROLS:
        print(f"E10 connected control {target}/{endpoint}", flush=True)
        data, fingerprints = load_dataset(target, endpoint)
        current_raw = np.asarray(
            project_2d(
                fingerprints,
                method="umap",
                random_state=MAP_SEED,
                umap_metric="jaccard",
            ),
            dtype=np.float32,
        )
        spectral_raw = explicit_spectral_projection(fingerprints)
        points = layouts_from_raw(current_raw, fingerprints)
        reference_probabilities, reference_shells = representation_weights(
            fingerprints, K
        )
        random_metrics = fixed_metrics(
            data.activity_pchembl.to_numpy(float),
            fingerprints,
            points["random"],
            reference_probabilities=reference_probabilities,
            reference_shells=reference_shells,
            k=K,
        )
        for layout in MAIN_LAYOUTS:
            row = metric_row(
                target,
                endpoint,
                layout,
                data,
                fingerprints,
                points[layout],
                points["grid_side"],
                reference_probabilities,
                reference_shells,
                float(random_metrics["fixed_e_grid"].mean()),
                float(random_metrics["retained_e_grid"].mean()),
            )
            e6_row = e6.loc[
                e6.target.eq(target)
                & e6.endpoint.eq(endpoint)
                & e6.layout.eq(layout)
            ].iloc[0]
            e2c_row = e2c.loc[
                e2c.target.eq(target)
                & e2c.endpoint.eq(endpoint)
                & e2c.layout.eq(layout)
            ].iloc[0]
            e6_difference = max(
                abs(float(row[column]) - float(e6_row[column]))
                for column in E6_METRIC_COLUMNS
            )
            e2c_difference = max(
                abs(float(row["retained_recall"]) - float(e2c_row.mean_recall)),
                abs(float(row["retained_e_grid"]) - float(e2c_row.mean_e_grid)),
                abs(float(row["retained_e_fp"]) - float(e2c_row.mean_e_fp)),
                abs(float(row["rho_retained"]) - float(e2c_row.rho)),
            )
            row.update({
                "coordinate_exact_vs_legacy_spectral_path": bool(
                    np.array_equal(current_raw, spectral_raw)
                ),
                "coordinate_max_absolute_difference_vs_legacy_spectral_path": float(
                    np.max(np.abs(current_raw - spectral_raw))
                ),
                "e6_metric_max_absolute_difference": e6_difference,
                "e2c_metric_max_absolute_difference": e2c_difference,
            })
            rows.append(row)
    result = pd.DataFrame(rows)
    checks = {
        "datasets": len(CONNECTED_CONTROLS),
        "rows": int(len(result)),
        "coordinates_exact": bool(
            result.coordinate_exact_vs_legacy_spectral_path.all()
        ),
        "e6_metrics_within_1e-12": bool(
            (result.e6_metric_max_absolute_difference <= TOLERANCE).all()
        ),
        "e2c_metrics_within_1e-12": bool(
            (result.e2c_metric_max_absolute_difference <= TOLERANCE).all()
        ),
    }
    if not all((
        checks["coordinates_exact"],
        checks["e6_metrics_within_1e-12"],
        checks["e2c_metrics_within_1e-12"],
    )):
        raise RuntimeError(f"Connected-control regression failed: {checks}")
    return result, checks


def comparison_tables(
    disconnected: set[tuple[str, str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    old = pd.read_csv(OLD_E6_RESULTS, float_precision="round_trip")
    new = pd.read_csv(E6_RESULTS, float_precision="round_trip")
    keys = ["target", "endpoint", "layout"]
    compared = old.loc[
        old[["target", "endpoint"]].apply(tuple, axis=1).isin(disconnected)
        & old.layout.isin(MAIN_LAYOUTS),
        keys + list(E6_METRIC_COLUMNS),
    ].merge(
        new.loc[
            new[["target", "endpoint"]].apply(tuple, axis=1).isin(disconnected)
            & new.layout.isin(MAIN_LAYOUTS),
            keys + list(E6_METRIC_COLUMNS),
        ],
        on=keys,
        suffixes=("_old", "_new"),
        validate="one_to_one",
    )
    expected_rows = len(disconnected) * len(MAIN_LAYOUTS)
    if len(compared) != expected_rows:
        raise ValueError(
            f"Expected {expected_rows} old/new rows, found {len(compared)}"
        )
    for column in E6_METRIC_COLUMNS:
        compared[f"{column}_difference"] = (
            compared[f"{column}_new"] - compared[f"{column}_old"]
        )
    rho_changes = compared.loc[
        compared.layout.eq("neighbor_refined"),
        [
            "target",
            "endpoint",
            "rho_retained_old",
            "rho_retained_new",
            "rho_retained_difference",
            "rho_fixed_old",
            "rho_fixed_new",
            "rho_fixed_difference",
        ],
    ].copy()
    return compared, rho_changes


def test_rows(path: Path, definition: str) -> pd.DataFrame:
    frame = pd.read_csv(path, float_precision="round_trip")
    if "definition" not in frame:
        frame["definition"] = definition
    elif "population" in frame:
        frame = frame.loc[frame.population.eq("all77")].copy()
    return frame.loc[
        frame.test_id.isin(("L1", "L2"))
        & frame.definition.eq(definition),
        [
            "test_id",
            "endpoint",
            "estimate",
            "ci_low",
            "ci_high",
            "holm_p",
        ],
    ]


def l1_l2_changes() -> pd.DataFrame:
    pairs = (
        (
            test_rows(OLD_E2C_TESTS, "retained_tie_inclusive"),
            test_rows(E2C_TESTS, "retained_tie_inclusive"),
            "retained_tie_inclusive",
        ),
        (
            test_rows(OLD_E6_TESTS, "fixed_k_random_tie_break"),
            test_rows(E6_TESTS, "fixed_k_random_tie_break"),
            "fixed_k_random_tie_break",
        ),
    )
    rows = []
    for old, new, definition in pairs:
        merged = old.merge(
            new,
            on=["test_id", "endpoint"],
            suffixes=("_old", "_new"),
            validate="one_to_one",
        )
        for item in merged.itertuples(index=False):
            sign_changed = bool(
                np.sign(item.estimate_old) != np.sign(item.estimate_new)
            )
            significance_changed = bool(
                (item.holm_p_old < 0.05) != (item.holm_p_new < 0.05)
            )
            rows.append({
                "test_id": item.test_id,
                "endpoint": item.endpoint,
                "definition": definition,
                "old_estimate": float(item.estimate_old),
                "new_estimate": float(item.estimate_new),
                "estimate_difference": float(
                    item.estimate_new - item.estimate_old
                ),
                "old_ci_low": float(item.ci_low_old),
                "new_ci_low": float(item.ci_low_new),
                "old_ci_high": float(item.ci_high_old),
                "new_ci_high": float(item.ci_high_new),
                "old_holm_p": float(item.holm_p_old),
                "new_holm_p": float(item.holm_p_new),
                "sign_changed": sign_changed,
                "significance_at_0.05_changed": significance_changed,
                "conclusion_changed": bool(
                    sign_changed or significance_changed
                ),
            })
    result = pd.DataFrame(rows)
    if len(result) != 8:
        raise ValueError(f"Expected eight L1/L2 rows, found {len(result)}")
    return result


def descriptive_changes() -> pd.DataFrame:
    old = pd.read_csv(OLD_E2C_RESULTS, float_precision="round_trip")
    new = pd.read_csv(E2C_RESULTS, float_precision="round_trip")
    rows = []
    for version, frame in (
        ("as_executed_pre05", old),
        ("v0.5", new),
    ):
        refined = frame.loc[
            frame.included.astype(bool)
            & frame.layout.eq("neighbor_refined")
        ].copy()
        for endpoint in ("IC50", "Ki"):
            population = refined.loc[refined.endpoint.eq(endpoint)]
            above = population.loc[
                population.rho_ci_low > 0, "target"
            ].astype(str).tolist()
            rows.append({
                "version": version,
                "endpoint": endpoint,
                "n_targets": int(len(population)),
                "rho_ci_below_zero": int(
                    (population.rho_ci_high < 0).sum()
                ),
                "rho_ci_includes_zero": int(
                    (
                        (population.rho_ci_low <= 0)
                        & (population.rho_ci_high >= 0)
                    ).sum()
                ),
                "rho_ci_above_zero": int(
                    (population.rho_ci_low > 0).sum()
                ),
                "targets_above_zero": ";".join(above),
                "rho_lambda_spearman": float(
                    spearmanr(population.rho, population["lambda"]).statistic
                ),
                "lambda_median": float(population["lambda"].median()),
            })
    return pd.DataFrame(rows)


def run() -> None:
    started_at = datetime.now(timezone.utc).isoformat()
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)

    old_e2c = pd.read_csv(OLD_E2C_RESULTS, float_precision="round_trip")
    new_e2c = pd.read_csv(E2C_RESULTS, float_precision="round_trip")
    keys = (
        new_e2c.loc[
            new_e2c.included.astype(bool), ["target", "endpoint"]
        ]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    if len(keys) != 77:
        raise ValueError(f"Expected 77 E2c datasets, found {len(keys)}")

    graph_rows = []
    for index, item in enumerate(keys.itertuples(index=False), 1):
        target, endpoint = str(item.target), str(item.endpoint)
        _, fingerprints = load_dataset(target, endpoint)
        count, sizes = component_summary(fingerprints)
        graph_rows.append({
            "target": target,
            "endpoint": endpoint,
            "n": int(len(fingerprints)),
            "components": count,
            "component_sizes": json.dumps(sizes),
        })
        print(
            f"E10 graph {index}/77 {target}/{endpoint}: {count} components",
            flush=True,
        )
    graph = pd.DataFrame(graph_rows)
    graph.to_csv(OUT / "graph_components.csv", index=False)
    disconnected = {
        (row.target, row.endpoint)
        for row in graph.loc[graph.components > 1].itertuples(index=False)
    }
    changed, change_rows = changed_e2c_datasets(old_e2c, new_e2c)
    change_rows.to_csv(QA / "e2c_changed_datasets.csv", index=False)
    if changed != disconnected:
        raise AssertionError(
            "E2c changed datasets do not equal disconnected datasets: "
            f"changed_only={sorted(changed - disconnected)}, "
            f"disconnected_only={sorted(disconnected - changed)}"
        )

    component_lookup = {
        (row.target, row.endpoint): row
        for row in graph.itertuples(index=False)
    }
    repeat_rows = []
    for index, (target, endpoint) in enumerate(
        sorted(disconnected), 1
    ):
        print(
            f"E10 repeat {index}/{len(disconnected)} {target}/{endpoint}",
            flush=True,
        )
        _, fingerprints = load_dataset(target, endpoint)
        first = np.asarray(
            project_2d(
                fingerprints,
                method="umap",
                random_state=MAP_SEED,
                umap_metric="jaccard",
            ),
            dtype=np.float32,
        )
        second = np.asarray(
            project_2d(
                fingerprints,
                method="umap",
                random_state=MAP_SEED,
                umap_metric="jaccard",
            ),
            dtype=np.float32,
        )
        component = component_lookup[(target, endpoint)]
        exact = bool(np.array_equal(first, second))
        repeat_rows.append({
            "target": target,
            "endpoint": endpoint,
            "components": int(component.components),
            "component_sizes": component.component_sizes,
            "first_projection_sha256": array_sha256(first),
            "second_projection_sha256": array_sha256(second),
            "exact_equal": exact,
            "max_absolute_difference": float(
                np.max(np.abs(first - second))
            ),
        })
        if not exact:
            raise RuntimeError(
                f"Projection repeatability failed: {target}/{endpoint}"
            )
    reproducibility = pd.DataFrame(repeat_rows)
    reproducibility.to_csv(OUT / "reproducibility.csv", index=False)

    comparison, rho_changes = comparison_tables(disconnected)
    comparison.to_csv(OUT / "old_new_comparison.csv", index=False)
    rho_changes.to_csv(OUT / "rho_refined_changes.csv", index=False)

    controls, control_checks = connected_controls()
    controls.to_csv(OUT / "connected_controls.csv", index=False)

    l1_l2 = l1_l2_changes()
    l1_l2.to_csv(OUT / "l1_l2_change_summary.csv", index=False)
    descriptive = descriptive_changes()
    descriptive.to_csv(OUT / "descriptive_changes.csv", index=False)

    comparison_max = max(
        float(comparison[f"{column}_difference"].abs().max())
        for column in E6_METRIC_COLUMNS
    )
    write_json(QA / "comparison_checks.json", {
        "datasets": int(len(keys)),
        "disconnected_datasets": sorted(
            f"{target}/{endpoint}" for target, endpoint in disconnected
        ),
        "e2c_changed_datasets": sorted(
            f"{target}/{endpoint}" for target, endpoint in changed
        ),
        "e2c_changed_equals_disconnected": changed == disconnected,
        "reprojection_exact_matches": int(
            reproducibility.exact_equal.sum()
        ),
        "reprojection_datasets": int(len(reproducibility)),
        "old_new_rows": int(len(comparison)),
        "old_new_max_absolute_difference": comparison_max,
        "connected_controls": control_checks,
        "l1_l2_rows": int(len(l1_l2)),
        "conclusion_changes": int(l1_l2.conclusion_changed.sum()),
        "graph_parameter_note": (
            "UMAP graph uses the project_2d parameters with init=random and "
            "n_epochs=1; neither setting changes the nearest-neighbor graph."
        ),
    })
    elapsed = time.perf_counter() - started
    write_json(QA / "timing.json", {
        "analysis_role": ANALYSIS_ROLE,
        "elapsed_seconds": elapsed,
        "row_counts": {
            "graph_components": int(len(graph)),
            "reproducibility": int(len(reproducibility)),
            "old_new_comparison": int(len(comparison)),
            "rho_refined_changes": int(len(rho_changes)),
            "connected_controls": int(len(controls)),
            "l1_l2_change_summary": int(len(l1_l2)),
            "descriptive_changes": int(len(descriptive)),
        },
    })
    finished_at = datetime.now(timezone.utc).isoformat()
    manifest = write_manifest(
        OUT,
        stage="E10",
        inputs=(
            OLD_E2C_RESULTS,
            OLD_E2C_TESTS,
            OLD_E6_RESULTS,
            OLD_E6_TESTS,
            E2C_RESULTS,
            E2C_TESTS,
            E6_RESULTS,
            E6_TESTS,
            S / "data/chembl37_api/selected_targets.csv",
        ),
        sources=(
            Path(__file__),
            SOURCE / "test_e10.py",
            SOURCE / "e2c_analysis.py",
            SOURCE / "e6_analysis.py",
            SOURCE / "e6_core.py",
            SOURCE / "manifest_utils.py",
            ROOT / "paper/grid_baselines.py",
            ROOT / "src/chemgridmap/core.py",
            ROOT / "src/chemgridmap/projection.py",
            ROOT / "src/chemgridmap/refinement.py",
            ROOT / "src/chemgridmap/reliability.py",
            ROOT / "src/chemgridmap/representations.py",
        ),
        parameters={
            "analysis_role": ANALYSIS_ROLE,
            "k": K,
            "map_seed": MAP_SEED,
            "subset_seed": 20261102,
            "subset_limit": 1200,
            "minimum_entries": MIN_ENTRIES,
            "projection": projection_parameters(1200),
            "graph_projection": {
                **projection_parameters(1200),
                "init": "random",
                "n_epochs": 1,
            },
            "connected_controls": [
                {"target": target, "endpoint": endpoint}
                for target, endpoint in CONNECTED_CONTROLS
            ],
            "tolerance": TOLERANCE,
        },
        started_at=started_at,
        finished_at=finished_at,
    )
    print(
        f"E10 COMPLETE elapsed={elapsed:.1f}s "
        f"disconnected={len(disconnected)} outputs={len(manifest['outputs_sha256'])}",
        flush=True,
    )


if __name__ == "__main__":
    run()
