"""Paper-only LAS, FLAS, and DPQ; package defaults remain unchanged.

Based on Barthel et al., CGF 42(1), 261–276 (2023), doi:10.1111/cgf.14718,
Equations 9–13, Algorithms 3–4, and the authors' LAS_FLAS reference code.
This independent implementation uses floating-point LAP and extends filtering
to vacant cells through mask normalization. Callers prepare features:
L2-normalized Morgan vectors for molecules, raw colors for RGB tests.
"""
from __future__ import annotations

from typing import Literal

import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.optimize import linear_sum_assignment
from scipy.sparse import csr_matrix
from scipy.spatial.distance import cdist


def random_grid_layout(n: int, grid_shape: tuple[int, int], seed: int = 42) -> np.ndarray:
    """Return distinct (col, row) cells in input order for a (height, width) grid."""
    height, width = grid_shape
    selected = np.random.default_rng(seed).choice(height * width, n, replace=False)
    return np.column_stack((selected % width, selected // width))


def masked_box_filter(
    grid_vectors: np.ndarray,
    occupied_mask: np.ndarray,
    size_y: int,
    size_x: int,
    fallback_mean: np.ndarray | None = None,
) -> np.ndarray:
    """Filter as U(mask X)/U(mask) with reflected edges and a global-mean fallback."""
    values = np.asarray(grid_vectors, dtype=float)
    mask = np.asarray(occupied_mask, dtype=float)
    if fallback_mean is None:
        fallback_mean = values[mask.astype(bool)].mean(axis=0)
    numerator = values * mask[..., None]
    denominator = mask
    for axis, size in ((1, size_x), (0, size_y)):
        numerator = uniform_filter1d(numerator, max(1, size), axis=axis, mode="reflect")
        denominator = uniform_filter1d(denominator, max(1, size), axis=axis, mode="reflect")
    # Floating-point filtering can leave rounding residuals near zero.
    populated = denominator > np.finfo(float).eps
    result = np.broadcast_to(fallback_mean, values.shape).copy()
    np.divide(numerator, denominator[..., None], out=result, where=populated[..., None])
    return result


def _prepare(features, grid_shape, seed, initial_cells):
    x = np.asarray(features, dtype=float)
    height, width = map(int, grid_shape)
    if x.ndim != 2 or not len(x) or height < 1 or width < 1 or len(x) > height * width:
        raise ValueError("Features must be a nonempty 2D matrix that fits in the grid.")
    if initial_cells is None:
        cells = random_grid_layout(len(x), (height, width), seed)
        radius = max(height, width) / 2 - 1
    else:
        cells = np.asarray(initial_cells, dtype=int).copy()
        radius = max(height, width) / 10
    if cells.shape != (len(x), 2) or len(np.unique(cells, axis=0)) != len(x):
        raise ValueError("Initial cells must be an N x 2 array of distinct coordinates.")
    if np.any(cells < 0) or np.any(cells[:, 0] >= width) or np.any(cells[:, 1] >= height):
        raise ValueError("Initial cells exceed the grid bounds.")
    placement = np.full(height * width, -1, dtype=int)
    placement[cells[:, 1] * width + cells[:, 0]] = np.arange(len(x))
    # CSR accelerates sparse Morgan dot products without changing the cost.
    sparse_x = csr_matrix(x) if x.shape[1] > 64 and np.count_nonzero(x) < 0.2 * x.size else None
    return x, (height, width), placement, radius, sparse_x


def _maps(x, placement, grid_shape, radius):
    height, width = grid_shape
    occupied = placement >= 0
    vectors = np.zeros((height * width, x.shape[1]), dtype=float)
    vectors[occupied] = x[placement[occupied]]
    return masked_box_filter(
        vectors.reshape(height, width, -1), occupied.reshape(height, width),
        min(height - 1, 2 * radius + 1), min(width - 1, 2 * radius + 1),
        fallback_mean=x.mean(axis=0),
    ).reshape(height * width, -1)


def _squared_cost(x, maps, sparse_x=None):
    if sparse_x is None and x.shape[1] <= 64:
        return cdist(x, maps, metric="sqeuclidean")
    product = (sparse_x if sparse_x is not None else x) @ maps.T
    cost = np.einsum("ij,ij->i", x, x)[:, None] + np.einsum("ij,ij->i", maps, maps)[None, :]
    cost -= 2 * product
    return np.maximum(cost, 0)


def _coordinates(placement, n, width):
    occupied = np.flatnonzero(placement >= 0)
    result = np.empty((n, 2), dtype=int)
    result[placement[occupied], 0] = occupied % width
    result[placement[occupied], 1] = occupied // width
    return result


def las(
    features: np.ndarray,
    grid_shape: tuple[int, int],
    *,
    seed: int = 42,
    initial_cells: np.ndarray | None = None,
) -> np.ndarray:
    """Algorithm 3: solve one global rectangular LAP per radius; retain input order."""
    x, shape, placement, radius_f, sparse_x = _prepare(features, grid_shape, seed, initial_cells)
    while radius_f >= 1:
        maps = _maps(x, placement, shape, round(radius_f))
        rows, columns = linear_sum_assignment(_squared_cost(x, maps, sparse_x))
        placement.fill(-1)
        placement[columns] = rows
        radius_f *= 0.9
    return _coordinates(placement, len(x), shape[1])


def flas(
    features: np.ndarray,
    grid_shape: tuple[int, int],
    *,
    seed: int = 42,
    initial_cells: np.ndarray | None = None,
    nc: int = 49,
) -> np.ndarray:
    """Algorithm 4: fix maps per radius and solve LAPs over random local cells."""
    x, shape, placement, radius_f, sparse_x = _prepare(features, grid_shape, seed, initial_cells)
    height, width = shape
    rng = np.random.default_rng(seed)
    # Match initialization's random stream by consuming the same cell selection.
    if initial_cells is None:
        rng.choice(height * width, len(x), replace=False)
    while radius_f >= 1:
        radius = int(radius_f)
        maps = _maps(x, placement, shape, radius)
        # Fixed maps allow costs to be reused across windows at this radius.
        cost = _squared_cost(x, maps, sparse_x)
        half_width = round(max(radius, (np.sqrt(nc) - 1) / 2))
        window_width = min(width, 2 * half_width + 1)
        window_height = min(height, 2 * half_width + 1)
        for _ in range(2 * int(height * width / nc) + 1):
            selected = int(rng.integers(height * width))
            center_y, center_x = divmod(selected, width)
            left = min(max(center_x - half_width, 0), width - window_width)
            top = min(max(center_y - half_width, 0), height - window_height)
            window = (np.arange(top, top + window_height)[:, None] * width
                      + np.arange(left, left + window_width)[None, :]).ravel()
            candidates = rng.choice(window, min(nc, len(window)), replace=False)
            molecules = placement[candidates]
            molecules = molecules[molecules >= 0]
            if not len(molecules):
                continue
            rows, columns = linear_sum_assignment(cost[np.ix_(molecules, candidates)])
            placement[candidates] = -1
            placement[candidates[columns]] = molecules[rows]
        radius_f *= 0.9
    return _coordinates(placement, len(x), width)


def dpq_scores(
    high_distances: np.ndarray,
    coordinates: np.ndarray,
    powers: tuple[int, ...] = (16, 2),
    *,
    ties: Literal["sorted", "mean"] = "sorted",
) -> dict[int, float]:
    """Equations 9–12 over occupied coordinates, with optimistic or mean ties.

    For equal 2D distances, ``sorted`` orders HD distances (DPQ_p), while
    ``mean`` replaces them with their group mean (pessimistic DPQ^-_p).
    """
    if ties not in ("sorted", "mean"):
        raise ValueError("ties must be 'sorted' or 'mean'.")
    distances = np.asarray(high_distances, dtype=float)
    n = len(distances)
    if distances.shape != (n, n) or len(coordinates) != n or n < 2:
        raise ValueError("DPQ requires at least two items and an N x N HD distance matrix.")
    keep = ~np.eye(n, dtype=bool)
    high = distances[keep].reshape(n, n - 1)
    low = cdist(np.asarray(coordinates, dtype=float), np.asarray(coordinates, dtype=float))
    low = low[keep].reshape(n, n - 1)
    average = float(high.mean())
    if average == 0:
        return {p: 1.0 for p in powers}
    high_sorted = np.sort(high, axis=1)
    low_order = np.lexsort((high, low), axis=1) if ties == "sorted" else np.argsort(low, axis=1)
    low_sorted = np.take_along_axis(high, low_order, axis=1)
    if ties == "mean":
        low_distances_sorted = np.take_along_axis(low, low_order, axis=1)
        for row, spatial_distances in zip(low_sorted, low_distances_sorted):
            # Integer grid coordinates yield exactly equal distances within shells.
            starts = np.r_[0, np.flatnonzero(np.diff(spatial_distances)) + 1]
            sizes = np.diff(np.r_[starts, n - 1])
            row[:] = np.repeat(np.add.reduceat(row, starts) / sizes, sizes)
    counts = np.arange(1, n, dtype=float)
    high_gain = (average - (np.cumsum(high_sorted, axis=1) / counts).mean(axis=0)) / average
    low_gain = np.maximum((average - (np.cumsum(low_sorted, axis=1) / counts).mean(axis=0)) / average, 0)
    return {p: float(np.linalg.norm(low_gain, ord=p) / np.linalg.norm(high_gain, ord=p))
            if np.linalg.norm(high_gain, ord=p) > 0 else 1.0 for p in powers}


def distance_preservation_quality(
    features: np.ndarray,
    coordinates: np.ndarray,
    *,
    p: int = 16,
    metric: str = "euclidean",
    distance_matrix: np.ndarray | None = None,
    ties: Literal["sorted", "mean"] = "sorted",
) -> float:
    """Compute DPQ from features or HD distances; use jaccard explicitly for Morgan."""
    distances = cdist(features, features, metric=metric) if distance_matrix is None else distance_matrix
    return dpq_scores(distances, coordinates, powers=(p,), ties=ties)[p]
