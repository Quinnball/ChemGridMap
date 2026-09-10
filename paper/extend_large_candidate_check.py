"""Check whether a larger candidate cap resolves the recorded full-scale warning."""
import json
import time
from pathlib import Path

import pandas as pd

from chemgridmap.core import assign_to_grid, compute_metrics
from run_revision_validation import coordinates, dump, ROOT


if __name__ == "__main__":
    out = ROOT / "paper/output/revision_v4"
    data = pd.read_csv(out / "chembl240/chembl240_molecule_level.csv")
    xy = pd.read_csv(out / "chembl240/full_umap_coordinates.csv").to_numpy()
    diagnostics = {}
    start = time.perf_counter()
    grid, side, _ = assign_to_grid(xy, method="sparse", sparse_max_neighbors=2048,
                                  diagnostics=diagnostics)
    elapsed = time.perf_counter() - start
    metrics = compute_metrics(xy, xy, grid, trustworthiness_sample_size=3000).iloc[0].to_dict()
    dump(out / "large_scale_extension.json", {"n": len(data), "seconds": elapsed,
                                             "initial_candidates": 128, "cap": 2048,
                                             "assignment": diagnostics, "metrics": metrics})
    coordinates(data, xy, grid, side).to_csv(out / "chembl240/full_grid_coordinates_cap2048.csv", index=False)
    print(json.dumps({"seconds": elapsed, **diagnostics}, indent=2), flush=True)
