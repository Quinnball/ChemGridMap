"""Per-molecule reliability of grid neighbourhoods.

The grid shows each molecule next to a handful of others. This module measures,
for every molecule, how many of its representation-space neighbours (for
example Morgan/Tanimoto neighbours) are also its grid neighbours, so that the
map can display where spatial proximity is trustworthy.
"""
from __future__ import annotations

from typing import List, Mapping, Sequence

import numpy as np
from scipy.spatial import cKDTree
from sklearn.neighbors import NearestNeighbors

RELATIVE_TOLERANCE = 1e-9
ABSOLUTE_TOLERANCE = 1e-12


def _shell(row: np.ndarray, k: int) -> np.ndarray:
    radius = np.partition(row, k - 1)[k - 1]
    return np.flatnonzero(row <= radius * (1 + RELATIVE_TOLERANCE) + ABSOLUTE_TOLERANCE)


def _weights(row: np.ndarray, shell: np.ndarray, k: int, *, index_space=None) -> dict:
    """Expected membership for random tie-breaking within the kth-distance shell."""
    radius = float(np.partition(row, k - 1)[k - 1])
    margin = RELATIVE_TOLERANCE * abs(radius) + ABSOLUTE_TOLERANCE
    strict = row < radius - margin
    lower = int(np.count_nonzero(strict))
    tied = shell[~strict[shell]]
    if index_space is None:
        strict_indices = np.flatnonzero(strict)
        tied_indices = tied
    else:
        index_space = np.asarray(index_space)
        strict_indices = index_space[np.flatnonzero(strict)]
        tied_indices = index_space[tied]
    weights = {int(index): 1.0 for index in strict_indices}
    tied_weight = (k - lower) / len(tied)
    weights.update({int(index): tied_weight for index in tied_indices})
    if abs(sum(weights.values()) - k) > 1e-9:
        raise ValueError("Fixed-k neighbour weights must sum to k.")
    return weights


def representation_neighbor_weights(matrix: np.ndarray, k: int = 10, metric: str = "jaccard",
                                    batch_size: int = 1024, return_kth: bool = False):
    """Tie-inclusive neighbour sets and random tie-breaking weights.

    Jaccard (Tanimoto) distances on binary fingerprints are computed exactly
    from integer bit counts with matrix products, which scales to tens of
    thousands of molecules. Other metrics use Euclidean distances.
    """
    matrix = np.asarray(matrix)
    n = len(matrix)
    if n < 2:
        result = ([np.empty(0, dtype=int) for _ in range(n)], [{} for _ in range(n)])
        return result + (np.full(n, np.nan),) if return_kth else result
    k = min(max(1, int(k)), n - 1)
    neighbors, weights, kth_values = [], [], []
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
                shell = _shell(row, k)
                neighbors.append(shell)
                weights.append(_weights(row, shell, k))
                kth_values.append(float(np.partition(row, k - 1)[k - 1]))
        result = (neighbors, weights)
        return result + (np.asarray(kth_values),) if return_kth else result
    if metric != "euclidean":
        raise ValueError("metric must be 'jaccard' or 'euclidean'.")
    points = np.asarray(matrix, dtype=np.float64)
    model = NearestNeighbors(n_neighbors=k + 1).fit(points)
    distances, indices = model.kneighbors(points)
    kth = np.array([d[idx != i][k - 1] for i, (d, idx) in enumerate(zip(distances, indices))])
    tree = cKDTree(points)
    balls = tree.query_ball_point(points, kth * (1 + RELATIVE_TOLERANCE) + ABSOLUTE_TOLERANCE)
    for i, row in enumerate(balls):
        shell = np.asarray(sorted(j for j in row if j != i), dtype=int)
        shell_distances = np.linalg.norm(points[shell] - points[i], axis=1)
        neighbors.append(shell)
        weights.append(_weights(
            shell_distances, np.arange(len(shell)), k, index_space=shell
        ))
    result = (neighbors, weights)
    return result + (kth,) if return_kth else result


def representation_neighbor_sets(matrix: np.ndarray, k: int = 10, metric: str = "jaccard",
                                 batch_size: int = 1024) -> List[np.ndarray]:
    """Tie-inclusive k-nearest-neighbour sets, excluding self."""
    neighbors, _ = representation_neighbor_weights(matrix, k=k, metric=metric, batch_size=batch_size)
    return neighbors


def grid_neighbor_sets(cells: np.ndarray, k: int = 10) -> List[np.ndarray]:
    """Tie-inclusive k-nearest-neighbour sets on the lattice (Euclidean cell distance)."""
    return representation_neighbor_sets(np.asarray(cells, dtype=np.float64), k=k, metric="euclidean")


def grid_neighbor_weights(cells: np.ndarray, k: int = 10):
    """Tie-inclusive lattice neighbours and fixed-k tie-breaking weights."""
    return representation_neighbor_weights(np.asarray(cells, dtype=np.float64), k=k, metric="euclidean")


def neighbor_recall(reference: List[np.ndarray], displayed: List[np.ndarray]) -> np.ndarray:
    """Fraction of each molecule's reference neighbours that are also displayed neighbours."""
    return np.array([len(np.intersect1d(r, d, assume_unique=True)) / len(r) if len(r) else np.nan
                     for r, d in zip(reference, displayed)])


def neighbor_recall_fixed_k(reference_weights: Sequence[Mapping[int, float]],
                            displayed_weights: Sequence[Mapping[int, float]],
                            *, k: int = 10) -> np.ndarray:
    """Expected recall after independently breaking each row's distance ties.

    Each input row assigns weights to the complete tie-inclusive kth-distance
    shell. Strictly closer neighbours receive weight 1, tied boundary
    neighbours share the remaining weight, and each row sums to ``k`` (or
    ``n - 1`` for datasets smaller than ``k + 1``).
    """
    if len(reference_weights) != len(displayed_weights):
        raise ValueError("Reference and displayed neighbourhoods must have the same length.")
    recall = []
    for reference, displayed in zip(reference_weights, displayed_weights):
        overlap = sum(weight * displayed.get(index, 0.0) for index, weight in reference.items())
        denominator = min(float(k), sum(reference.values()))
        recall.append(overlap / denominator if denominator else np.nan)
    return np.asarray(recall, dtype=float)


def reliability_columns(reference: List[np.ndarray], displayed: List[np.ndarray],
                        *, reference_weights=None, displayed_weights=None, k: int = 10,
                        kth_values=None, kth_column: str = "neighbor_kth_similarity") -> dict:
    """Counts, tie-inclusive recall and fixed-k expected recall."""
    shared = [len(np.intersect1d(r, d, assume_unique=True)) for r, d in zip(reference, displayed)]
    columns = {"neighbor_recall": neighbor_recall(reference, displayed),
               "neighbor_shared_count": np.array(shared, dtype=int),
               "representation_neighbor_count": np.array([len(r) for r in reference], dtype=int),
               "grid_neighbor_count": np.array([len(d) for d in displayed], dtype=int)}
    if reference_weights is not None and displayed_weights is not None:
        columns["neighbor_recall_fixed_k"] = neighbor_recall_fixed_k(
            reference_weights, displayed_weights, k=k
        )
    if kth_values is not None:
        columns[kth_column] = np.asarray(kth_values, dtype=float)
    return columns
