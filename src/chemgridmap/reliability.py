"""Per-molecule reliability of grid neighbourhoods.

The grid shows each molecule next to a handful of others. This module measures,
for every molecule, how many of its representation-space neighbours (for
example Morgan/Tanimoto neighbours) are also its grid neighbours, so that the
map can display where spatial proximity is trustworthy.
"""
from __future__ import annotations

from typing import List

import numpy as np
from scipy.spatial import cKDTree
from sklearn.neighbors import NearestNeighbors

RELATIVE_TOLERANCE = 1e-9
ABSOLUTE_TOLERANCE = 1e-12


def _shell(row: np.ndarray, k: int) -> np.ndarray:
    radius = np.partition(row, k - 1)[k - 1]
    return np.flatnonzero(row <= radius * (1 + RELATIVE_TOLERANCE) + ABSOLUTE_TOLERANCE)


def representation_neighbor_sets(matrix: np.ndarray, k: int = 10, metric: str = "jaccard",
                                 batch_size: int = 1024) -> List[np.ndarray]:
    """Tie-inclusive k-nearest-neighbour sets, excluding self.

    Jaccard (Tanimoto) distances on binary fingerprints are computed exactly
    from integer bit counts with matrix products, which scales to tens of
    thousands of molecules. Other metrics use Euclidean distances.
    """
    matrix = np.asarray(matrix)
    n = len(matrix)
    if n < 2:
        return [np.empty(0, dtype=int) for _ in range(n)]
    k = min(max(1, int(k)), n - 1)
    neighbors = []
    if metric == "jaccard":
        bits = (matrix != 0).astype(np.float32)
        counts = bits.sum(axis=1).astype(np.float64)
        for start in range(0, n, int(batch_size)):
            stop = min(n, start + int(batch_size))
            shared = np.rint(bits[start:stop] @ bits.T).astype(np.float64)
            union = counts[start:stop, None] + counts[None, :] - shared
            distance = np.where(union > 0, 1.0 - shared / np.where(union > 0, union, 1.0), 0.0)
            for local, row in enumerate(distance):
                row[start + local] = np.inf
                neighbors.append(_shell(row, k))
        return neighbors
    if metric != "euclidean":
        raise ValueError("metric must be 'jaccard' or 'euclidean'.")
    points = np.asarray(matrix, dtype=np.float64)
    model = NearestNeighbors(n_neighbors=k + 1).fit(points)
    distances, indices = model.kneighbors(points)
    kth = np.array([d[idx != i][k - 1] for i, (d, idx) in enumerate(zip(distances, indices))])
    tree = cKDTree(points)
    balls = tree.query_ball_point(points, kth * (1 + RELATIVE_TOLERANCE) + ABSOLUTE_TOLERANCE)
    return [np.asarray(sorted(j for j in row if j != i), dtype=int) for i, row in enumerate(balls)]


def grid_neighbor_sets(cells: np.ndarray, k: int = 10) -> List[np.ndarray]:
    """Tie-inclusive k-nearest-neighbour sets on the lattice (Euclidean cell distance)."""
    return representation_neighbor_sets(np.asarray(cells, dtype=np.float64), k=k, metric="euclidean")


def neighbor_recall(reference: List[np.ndarray], displayed: List[np.ndarray]) -> np.ndarray:
    """Fraction of each molecule's reference neighbours that are also displayed neighbours."""
    return np.array([len(np.intersect1d(r, d, assume_unique=True)) / len(r) if len(r) else np.nan
                     for r, d in zip(reference, displayed)])


def reliability_columns(reference: List[np.ndarray], displayed: List[np.ndarray]) -> dict:
    shared = [len(np.intersect1d(r, d, assume_unique=True)) for r, d in zip(reference, displayed)]
    return {"neighbor_recall": neighbor_recall(reference, displayed),
            "neighbor_shared_count": np.array(shared, dtype=int),
            "representation_neighbor_count": np.array([len(r) for r in reference], dtype=int),
            "grid_neighbor_count": np.array([len(d) for d in displayed], dtype=int)}
