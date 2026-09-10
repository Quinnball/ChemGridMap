"""Audit adaptive matching against dense controls and a full-size UMAP case."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd

from chemgridmap.core import assign_to_grid, compute_metrics
from chemgridmap.projection import project_2d
from chemgridmap.representations import morgan_fingerprints
from chemgridmap.provenance import environment_versions
from run_revision_validation import coordinates, evaluate, dump, ROOT


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=ROOT / "paper/output/revision_v4")
    parser.add_argument("--large", action="store_true")
    args = parser.parse_args()
    out = args.output
    rows = []
    for path in sorted(out.glob("chembl*/coordinates_seed_*.csv")):
        data = pd.read_csv(path)
        xy = data[["projection_x_raw", "projection_y_raw"]].to_numpy()
        dense = data[["grid_x", "grid_y"]].to_numpy()
        info = {}
        start = time.perf_counter()
        grid, _, _ = assign_to_grid(xy, method="sparse", diagnostics=info)
        elapsed = time.perf_counter() - start
        baseline = evaluate(xy, dense, data)
        metrics = evaluate(xy, grid, data)
        rows.append({"target": path.parent.name.upper(), "seed": int(path.stem.split("_")[-1]),
                     "n": len(data), "assignment_seconds": elapsed,
                     "dense_cost_ratio": metrics["mean_squared_grid_displacement"] /
                     baseline["mean_squared_grid_displacement"],
                     "same_cell_fraction": float(np.all(grid == dense, axis=1).mean()),
                     **metrics, **{k: (json.dumps(v) if isinstance(v, list) else v) for k, v in info.items()}})
        print(path.parent.name, rows[-1]["seed"], info["nearest_candidates"],
              rows[-1]["dense_cost_ratio"], info["converged"], flush=True)
    pd.DataFrame(rows).to_csv(out / "adaptive_validation.csv", index=False)
    if args.large:
        data = pd.read_csv(out / "chembl240/chembl240_molecule_level.csv")
        coordinate_path = out / "chembl240/full_umap_coordinates.csv"
        projection_seconds = None
        if coordinate_path.exists():
            xy = pd.read_csv(coordinate_path)[["projection_x_raw", "projection_y_raw"]].to_numpy()
        else:
            print(f"Full-size parent-normalized hERG UMAP: n={len(data)}", flush=True)
            t0 = time.perf_counter()
            fp = morgan_fingerprints(data.canonical_smiles.tolist())
            xy = project_2d(fp, method="umap", random_state=42, umap_metric="jaccard")
            projection_seconds = time.perf_counter() - t0
            pd.DataFrame(xy, columns=["projection_x_raw", "projection_y_raw"]).to_csv(coordinate_path, index=False)
        scale_rows = []
        for n in [1000, 2500, 5000, len(data)]:
            idx = np.sort(np.random.default_rng(20260908).choice(len(data), n, replace=False))
            for repeat in range(3):
                info = {}
                t0 = time.perf_counter()
                grid, side, _ = assign_to_grid(xy[idx], method="sparse", diagnostics=info)
                elapsed = time.perf_counter() - t0
                scale_rows.append({"n": n, "repeat": repeat, "seconds": elapsed,
                                   **{k: (json.dumps(v) if isinstance(v, list) else v) for k, v in info.items()}})
            print("Scale", n, elapsed, info["nearest_candidates"], info["converged"], flush=True)
            if n == len(data):
                metrics = compute_metrics(xy, xy, grid, trustworthiness_sample_size=3000)
                metrics.to_csv(out / "large_scale_metrics.csv", index=False)
                coordinates(data, xy, grid, side).to_csv(out / "chembl240/full_grid_coordinates.csv", index=False)
        pd.DataFrame(scale_rows).to_csv(out / "large_scale_runtime.csv", index=False)
        dump(out / "large_scale_manifest.json", {"n": len(data), "projection_and_fingerprint_seconds": projection_seconds,
                                                "environment": environment_versions(),
                                                "annotation": "Geometry-only scaling check, no claim of biological domain discovery.",
                                                "seed": 42, "initial_candidates": 128, "cap": 1024})


if __name__ == "__main__":
    main()
