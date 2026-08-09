"""Select the sparse grid candidate count on the CHEMBL240 projection."""

from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
from sklearn.manifold import trustworthiness

from chemgridmap.core import (
    assign_to_grid,
    knn_overlap,
    label_purity,
    normalize_coordinates,
    prepare_molecules,
)
from chemgridmap.projection import project_2d
from chemgridmap.representations import compute_representation


ROOT = Path("validation_data/chembl37_chembl240_large_scale")
INPUT = ROOT / "chembl37_chembl240_ic50_molecule_level.csv"


def main() -> int:
    source = pd.read_csv(INPUT)
    prepared, _ = prepare_molecules(
        source,
        smiles_col="canonical_smiles",
        value_col="activity_pchembl",
        label_col="activity_class",
    )
    representation, _ = compute_representation(prepared, representation="morgan")
    projection = project_2d(representation, method="pca", random_state=42)
    projection_01 = normalize_coordinates(projection)
    labels = prepared["activity_class"].astype(str).to_numpy()
    sample = np.sort(
        np.random.RandomState(42).choice(len(prepared), size=3000, replace=False)
    )

    rows = []
    for neighbors in [32, 64, 128, 256, 512]:
        started = time.perf_counter()
        grid, _, _ = assign_to_grid(
            projection,
            method="sparse",
            sparse_neighbors=neighbors,
        )
        elapsed = time.perf_counter() - started
        displacement = np.linalg.norm(projection_01 - grid, axis=1)
        rows.append(
            {
                "sparse_neighbors": neighbors,
                "assignment_seconds": elapsed,
                "mean_displacement": float(np.mean(displacement)),
                "p95_displacement": float(np.quantile(displacement, 0.95)),
                "projection_grid_knn_overlap_k10": knn_overlap(
                    projection_01, grid, k=10
                ),
                "grid_label_purity_k10": label_purity(grid, labels, k=10),
                "trustworthiness_k10_sample3000": float(
                    trustworthiness(
                        projection_01[sample],
                        grid[sample],
                        n_neighbors=10,
                    )
                ),
            }
        )
        print(rows[-1], flush=True)

    frame = pd.DataFrame(rows)
    frame.to_csv(ROOT / "sparse_neighbor_benchmark.csv", index=False)
    (ROOT / "sparse_neighbor_benchmark.json").write_text(
        json.dumps(rows, indent=2), encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
