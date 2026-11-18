"""Check color layouts, vacancies, reproducibility, and DPQ for T2a."""
from itertools import groupby

import numpy as np
import pytest
from scipy.ndimage import uniform_filter1d
from scipy.spatial.distance import cdist
from scipy.sparse import csr_matrix

from grid_baselines import (
    distance_preservation_quality, dpq_scores, flas, las, masked_box_filter, random_grid_layout,
)


def adjacent_distance(features, cells, shape):
    height, width = shape
    image = np.empty((height, width, features.shape[1]))
    image[cells[:, 1], cells[:, 0]] = features
    distances = [np.linalg.norm(np.diff(image, axis=axis), axis=2).ravel() for axis in (0, 1)]
    return np.concatenate(distances).mean()


def test_dpq_exact_gradient_and_random_layout():
    yy, xx = np.indices((20, 20))
    coordinates = np.column_stack((xx.ravel(), yy.ravel()))
    colors = coordinates.astype(float) / 19
    for ties in ("sorted", "mean"):
        assert distance_preservation_quality(colors, coordinates, ties=ties) == pytest.approx(1.0, abs=1e-9)
    random_cells = random_grid_layout(len(colors), (20, 20), seed=42)
    optimistic = distance_preservation_quality(colors, random_cells)
    pessimistic = distance_preservation_quality(colors, random_cells, ties="mean")
    assert pessimistic < 0.1, pessimistic
    assert pessimistic < optimistic < 0.6, optimistic


def test_las_rgb_neighbor_quality():
    # Preserve the reference experiment's legacy NumPy RNG and seed=3.
    np.random.seed(3)
    colors = np.random.random((32 * 32, 3))
    random_cells = random_grid_layout(len(colors), (32, 32), seed=3)
    result = las(colors, (32, 32), seed=3)
    ratio = adjacent_distance(colors, result, (32, 32)) / adjacent_distance(colors, random_cells, (32, 32))
    quality = distance_preservation_quality(colors, result)
    assert ratio <= 0.5, ratio
    assert quality >= 0.7, quality


@pytest.mark.parametrize("algorithm", [las, flas])
@pytest.mark.parametrize("projection_initialized", [False, True])
def test_vacancies_unique_and_reproducible(algorithm, projection_initialized):
    shape = (20, 20)
    colors = np.random.default_rng(7).random((160, 3))
    initial = random_grid_layout(len(colors), shape, seed=21) if projection_initialized else None
    one = algorithm(colors, shape, seed=42, initial_cells=initial)
    two = algorithm(colors, shape, seed=42, initial_cells=initial)
    np.testing.assert_array_equal(one, two)
    assert one.shape == (160, 2)
    assert len(np.unique(one, axis=0)) == 160
    assert np.all(one >= 0)
    assert np.all(one[:, 0] < shape[1]) and np.all(one[:, 1] < shape[0])


def test_full_mask_matches_box_filter():
    values = np.random.default_rng(42).random((12, 15, 3))
    ordinary = uniform_filter1d(values, 7, axis=1, mode="reflect")
    ordinary = uniform_filter1d(ordinary, 5, axis=0, mode="reflect")
    masked = masked_box_filter(values, np.ones((12, 15), dtype=bool), 5, 7)
    np.testing.assert_allclose(masked, ordinary, atol=1e-14, rtol=0)


def test_empty_filter_window_uses_global_mean():
    values = np.zeros((10, 10, 2))
    occupied = np.zeros((10, 10), dtype=bool)
    occupied[0, 0] = True
    values[0, 0] = (0.3, 0.7)
    filtered = masked_box_filter(values, occupied, 3, 3)
    np.testing.assert_allclose(filtered, np.broadcast_to([0.3, 0.7], filtered.shape), atol=1e-14)


@pytest.mark.parametrize("ties", ["sorted", "mean"])
def test_dpq_matches_independent_scalar_definition(ties):
    # Python sorting and scalar means are independent of the vectorized code.
    colors = np.array([[0., 0.], [1., 0.], [0., 2.], [3., 1.], [2., 4.]])
    cells = np.array([[0, 0], [1, 1], [0, 1], [1, 0], [2, 0]])
    high = cdist(colors, colors)
    low = cdist(cells, cells)
    n = len(colors)
    dbar = sum(high[i, j] for i in range(n) for j in range(n) if i != j) / (n * (n - 1))
    ranks_high = [sorted((j for j in range(n) if i != j), key=lambda j: high[i, j]) for i in range(n)]
    ranked_low_distances = []
    for i in range(n):
        neighbors = sorted((j for j in range(n) if i != j), key=lambda j: low[i, j])
        values = []
        for _, group in groupby(neighbors, key=lambda j: low[i, j]):
            group_distances = [high[i, j] for j in group]
            if ties == "mean":
                values.extend([sum(group_distances) / len(group_distances)] * len(group_distances))
            else:
                values.extend(sorted(group_distances))
        ranked_low_distances.append(values)
    gains_high, gains_low = [], []
    for k in range(1, n):
        dh = sum(sum(high[i, j] for j in ranks_high[i][:k]) / k for i in range(n)) / n
        dl = sum(sum(values[:k]) / k for values in ranked_low_distances) / n
        gains_high.append((dbar - dh) / dbar)
        gains_low.append(max((dbar - dl) / dbar, 0))
    for p in (2, 16):
        expected = (sum(value ** p for value in gains_low) / sum(value ** p for value in gains_high)) ** (1 / p)
        assert distance_preservation_quality(colors, cells, p=p, ties=ties) == pytest.approx(expected, abs=1e-12)
        assert dpq_scores(high, cells, ties=ties)[p] == pytest.approx(expected, abs=1e-12)


def test_dpq_rejects_unknown_tie_rule():
    with pytest.raises(ValueError, match="ties must be"):
        dpq_scores(np.array([[0., 1.], [1., 0.]]), np.array([[0, 0], [1, 0]]), ties="unknown")


def test_sparse_morgan_cost_matches_direct_euclidean():
    from grid_baselines import _squared_cost
    rng = np.random.default_rng(42)
    features = (rng.random((40, 128)) < .05).astype(float)
    features /= np.linalg.norm(features, axis=1, keepdims=True)
    maps = rng.random((70, 128)) / 20
    cost = _squared_cost(features, maps, csr_matrix(features))
    np.testing.assert_allclose(cost, cdist(features, maps, metric="sqeuclidean"), atol=1e-13)
