"""E4d: occupancy-0.4 layout baselines across ten seeds."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t


ROOT = Path(__file__).resolve().parents[4]
REV = ROOT / "project1/JCIM_reliability_revision_20261024"
TASK_S = REV / "strengthening"
E1B_SOURCE = REV / "strengthening/source/e1b_neighbors.py"
E1B_ERRORS = REV / "strengthening/output/E1B/layout_activity_error.csv"
LAYOUT_ROOT = ROOT / "chemgridmap/paper/output/reliability_v1"
CHECKPOINT_ROOT = LAYOUT_ROOT / "baselines/checkpoints"
OUTPUT = TASK_S / "output/E4/layout_baselines_ten_seed.csv"
VALIDATION = TASK_S / "qa/E4/e4d_validation.json"
K = 10
OCCUPANCY = 0.4
ANALYSIS_ROLE = "post hoc"
BASELINES = {
    "las_random": "LAS (random initialization)",
    "las_projection": "LAS (projection initialization)",
    "flas_random": "FLAS (random initialization)",
}
LAYOUT_LABELS = {
    **BASELINES,
    "neighbor_refined": "Refined grid",
    "minimum_displacement": "Minimum-displacement grid",
}


def _load_e1b():
    spec = importlib.util.spec_from_file_location("e4d_e1b_neighbors", E1B_SOURCE)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import {E1B_SOURCE}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def t_interval(differences: np.ndarray) -> dict[str, float | int]:
    """Two-sided 95% Student t interval for paired seed differences."""
    values = np.asarray(differences, dtype=float)
    if values.ndim != 1 or len(values) < 2 or not np.isfinite(values).all():
        raise ValueError("Paired differences must be a finite vector with n >= 2")
    n = len(values)
    mean = float(values.mean())
    sd = float(values.std(ddof=1))
    critical = float(t.ppf(0.975, n - 1))
    half_width = critical * sd / np.sqrt(n)
    return {
        "n_seeds": n,
        "paired_difference_mean": mean,
        "paired_difference_sd": sd,
        "t_critical": critical,
        "t_ci_low": mean - half_width,
        "t_ci_high": mean + half_width,
        "t_degrees_of_freedom": n - 1,
    }


def _layout_metrics(e1b, activity: np.ndarray, reference, points: np.ndarray) -> dict[str, float]:
    if points.ndim != 2 or points.shape[1] != 2 or not np.isfinite(points).all():
        raise ValueError("Grid coordinates must be a finite two-column array")
    metrics = e1b.neighborhood_metrics(
        activity, reference, e1b.grid_neighbor_sets(points, k=K)
    )
    return {
        "mean_e_grid": float(metrics.e_grid.mean()),
        "mean_recall": float(metrics.recall.mean()),
    }


def _summary_rows(target: str, target_label: str, raw: pd.DataFrame) -> list[dict]:
    contrasts = [
        (baseline, "neighbor_refined", f"{baseline}_minus_neighbor_refined")
        for baseline in BASELINES
    ]
    contrasts.append(
        (
            "neighbor_refined",
            "minimum_displacement",
            "neighbor_refined_minus_minimum_displacement",
        )
    )
    rows = []
    for left, right, comparison in contrasts:
        left_rows = raw.loc[raw.layout == left].set_index("seed")
        right_rows = raw.loc[raw.layout == right].set_index("seed")
        seeds = sorted(set(left_rows.index) & set(right_rows.index))
        if len(seeds) != raw.seed.nunique():
            raise ValueError(f"Incomplete paired seeds: {target}/{comparison}")
        for metric in ("mean_e_grid", "mean_recall"):
            differences = (
                left_rows.loc[seeds, metric].to_numpy(float)
                - right_rows.loc[seeds, metric].to_numpy(float)
            )
            interval = t_interval(differences)
            rows.append(
                {
                    "analysis_role": ANALYSIS_ROLE,
                    "row_type": "paired_t_interval",
                    "target": target,
                    "target_label": target_label,
                    "occupancy": OCCUPANCY,
                    "comparison": comparison,
                    "left_layout": left,
                    "left_layout_label": LAYOUT_LABELS[left],
                    "right_layout": right,
                    "right_layout_label": LAYOUT_LABELS[right],
                    "metric": metric,
                    "left_mean": float(left_rows.loc[seeds, metric].mean()),
                    "right_mean": float(right_rows.loc[seeds, metric].mean()),
                    **interval,
                    "ci_level": 0.95,
                    "ci_method": "two-sided Student t interval across paired seeds",
                }
            )
    return rows


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_table() -> tuple[pd.DataFrame, dict]:
    e1b = _load_e1b()
    raw_rows = []
    identities = []
    input_paths = [E1B_SOURCE, E1B_ERRORS]
    for target, target_label in e1b.TARGETS.items():
        data, fingerprints = e1b.load_target(target)
        activity = data.activity_pchembl.to_numpy(float)
        reference = e1b.representation_neighbor_sets(fingerprints, k=K)
        input_paths.append(e1b.FROZEN / target / "map_molecules.csv")
        for seed in e1b.SEEDS:
            layout_path = LAYOUT_ROOT / "layouts" / target / f"coordinates_seed_{seed}.csv"
            checkpoint_path = (
                CHECKPOINT_ROOT / f"{target}_occupancy{OCCUPANCY}_seed{seed}.npz"
            )
            input_paths.extend([layout_path, checkpoint_path])
            layout_frame = pd.read_csv(layout_path)
            if data.molecule_identity_key.tolist() != layout_frame.molecule_identity_key.tolist():
                raise ValueError(f"Refined layout identity/order mismatch: {target}/{seed}")
            points = {
                "neighbor_refined": layout_frame[
                    ["refined_grid_col", "refined_grid_row"]
                ].to_numpy(float),
                "minimum_displacement": layout_frame[["grid_col", "grid_row"]].to_numpy(float),
            }
            with np.load(checkpoint_path, allow_pickle=False) as checkpoint:
                checkpoint_identity = checkpoint["molecule_identity_key"].astype(str)
                if not np.array_equal(checkpoint_identity, data.molecule_identity_key.to_numpy(str)):
                    raise ValueError(f"Checkpoint identity/order mismatch: {target}/{seed}")
                for baseline in BASELINES:
                    points[baseline] = checkpoint[baseline]
                identities.append(
                    {
                        "target": target,
                        "seed": seed,
                        "n": len(data),
                        "checkpoint_aligned": True,
                        "refined_layout_aligned": True,
                    }
                )
            if points["neighbor_refined"].shape != (len(data), 2):
                raise ValueError(f"Bad refined layout shape: {target}/{seed}")
            for layout, coordinates in points.items():
                raw_rows.append(
                    {
                        "analysis_role": ANALYSIS_ROLE,
                        "row_type": "seed_layout",
                        "target": target,
                        "target_label": target_label,
                        "occupancy": OCCUPANCY,
                        "seed": seed,
                        "layout": layout,
                        "layout_label": LAYOUT_LABELS[layout],
                        "n": len(data),
                        **_layout_metrics(e1b, activity, reference, coordinates),
                    }
                )
        print(f"E4d {target_label}: 10 seeds complete", flush=True)

    raw = pd.DataFrame(raw_rows)
    summary_rows = []
    for target, target_label in e1b.TARGETS.items():
        summary_rows.extend(
            _summary_rows(target, target_label, raw.loc[raw.target == target])
        )
    table = pd.concat([raw, pd.DataFrame(summary_rows)], ignore_index=True)

    # The main-layout values must reproduce the accepted E1-B output exactly.
    frozen = pd.read_csv(E1B_ERRORS, float_precision="round_trip")
    frozen = frozen.loc[
        (frozen.row_type == "seed_layout")
        & frozen.layout.isin(["neighbor_refined", "minimum_displacement"]),
        ["target", "seed", "layout", "mean_e_grid", "mean_recall"],
    ]
    ours = raw.loc[
        raw.layout.isin(["neighbor_refined", "minimum_displacement"]),
        ["target", "seed", "layout", "mean_e_grid", "mean_recall"],
    ]
    key = ["target", "seed", "layout"]
    frozen = frozen.set_index(key).sort_index()
    ours = ours.set_index(key).sort_index()
    if not frozen.index.equals(ours.index):
        raise ValueError("E1-B main-layout comparison keys differ")
    max_e_grid = float(np.abs(ours.mean_e_grid - frozen.mean_e_grid).max())
    max_recall = float(np.abs(ours.mean_recall - frozen.mean_recall).max())
    tolerance = 1e-12
    if max(max_e_grid, max_recall) > tolerance:
        raise ValueError(
            f"E1-B main-layout regression failed: e_grid={max_e_grid}, recall={max_recall}"
        )
    validation = {
        "analysis_role": ANALYSIS_ROLE,
        "seed_count_per_target": len(e1b.SEEDS),
        "seeds": list(e1b.SEEDS),
        "targets": list(e1b.TARGETS),
        "baselines": list(BASELINES),
        "occupancy": OCCUPANCY,
        "k": K,
        "n_rows": len(table),
        "n_seed_rows": len(raw),
        "n_summary_rows": len(summary_rows),
        "identity_checks": len(identities),
        "e1b_regression": {
            "rows": len(ours),
            "max_absolute_e_grid_difference": max_e_grid,
            "max_absolute_recall_difference": max_recall,
            "tolerance": tolerance,
            "passed": True,
        },
        "t_interval": "two-sided 95% Student t interval; paired target-seed differences",
        "input_sha256": {
            str(path.relative_to(ROOT)): _sha256(path) for path in sorted(set(input_paths))
        },
    }
    return table, validation


def run(output: Path = OUTPUT, validation_path: Path = VALIDATION) -> pd.DataFrame:
    table, validation = build_table()
    output.parent.mkdir(parents=True, exist_ok=True)
    validation_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output, index=False)
    validation["output_sha256"] = _sha256(output)
    validation_path.write_text(
        json.dumps(validation, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    )
    print(f"Wrote {len(table)} rows to {output}", flush=True)
    print(json.dumps(validation["e1b_regression"], ensure_ascii=False), flush=True)
    return table


if __name__ == "__main__":
    run()
