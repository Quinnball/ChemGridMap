"""Profile the main stages of the CHEMBL240 large-scale validation."""

from __future__ import annotations

import json
from pathlib import Path
import resource
import sys
import time

import numpy as np
import pandas as pd

from chemgridmap.core import (
    assign_to_grid,
    compute_metrics,
    normalize_coordinates,
    prepare_molecules,
)
from chemgridmap.plotting import render_grid_map
from chemgridmap.projection import project_2d
from chemgridmap.representations import compute_representation


ROOT = Path("validation_data/chembl37_chembl240_large_scale")
INPUT = ROOT / "chembl37_chembl240_ic50_molecule_level.csv"


def peak_rss_mib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / (1024.0 * 1024.0) if sys.platform == "darwin" else value / 1024.0


def record(results, name, started):
    results[name] = {
        "seconds": time.perf_counter() - started,
        "peak_rss_mib": peak_rss_mib(),
    }
    print(name, results[name], flush=True)


def main() -> int:
    results = {}

    started = time.perf_counter()
    source = pd.read_csv(INPUT)
    prepared, _ = prepare_molecules(
        source,
        smiles_col="canonical_smiles",
        value_col="activity_pchembl",
        label_col="activity_class",
    )
    record(results, "read_and_prepare", started)

    started = time.perf_counter()
    representation, _ = compute_representation(prepared, representation="morgan")
    record(results, "morgan_fingerprint", started)

    started = time.perf_counter()
    projection = project_2d(representation, method="pca", random_state=42)
    record(results, "pca_projection", started)

    started = time.perf_counter()
    grid, rows, cols = assign_to_grid(projection, method="sparse")
    record(results, "sparse_grid_assignment", started)

    started = time.perf_counter()
    values = prepared["activity_pchembl"].to_numpy(dtype=float)
    labels = prepared["activity_class"].astype(str).to_numpy()
    metrics = compute_metrics(
        representation,
        projection,
        grid,
        values=values,
        labels=labels,
        k=10,
        trustworthiness_sample_size=3000,
        random_state=42,
    )
    record(results, "quality_metrics", started)

    started = time.perf_counter()
    result_data = prepared.copy()
    projection_01 = normalize_coordinates(projection)
    result_data["projection_x"] = projection_01[:, 0]
    result_data["projection_y"] = projection_01[:, 1]
    result_data["grid_x"] = grid[:, 0]
    result_data["grid_y"] = grid[:, 1]
    result_data["grid_col"] = np.rint(grid[:, 0] * (cols - 1)).astype(int)
    result_data["grid_row"] = np.rint(grid[:, 1] * (rows - 1)).astype(int)
    render_grid_map(
        result_data,
        ROOT / "profile_render_overview",
        detail="overview",
        formats=["png"],
    )
    record(results, "overview_render", started)

    payload = {
        "n_molecules": int(len(prepared)),
        "stages": results,
        "metrics": metrics.iloc[0].to_dict(),
    }
    (ROOT / "large_scale_stage_profile.json").write_text(
        json.dumps(payload, indent=2, default=float), encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
