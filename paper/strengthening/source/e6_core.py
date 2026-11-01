"""E6 post-hoc fixed-k neighbourhood metrics."""
from __future__ import annotations

from pathlib import Path
import sys
from typing import Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / "paper"), str(ROOT / "src")]

from chemgridmap.reliability import (  # noqa: E402
    ABSOLUTE_TOLERANCE,
    RELATIVE_TOLERANCE,
    neighbor_recall,
)


def _distance_weights(row: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray, float]:
    """Return fixed-k probabilities where distances are sorted ascending."""
    radius = float(np.partition(row, k - 1)[k - 1])
    margin = RELATIVE_TOLERANCE * abs(radius) + ABSOLUTE_TOLERANCE
    cutoff = radius * (1.0 + RELATIVE_TOLERANCE) + ABSOLUTE_TOLERANCE
    shell = np.flatnonzero(row <= cutoff)
    strict = row < radius - margin
    lower = int(np.count_nonzero(strict))
    tied = shell[~strict[shell]]
    weights = np.zeros_like(row, dtype=float)
    weights[strict] = 1.0
    weights[tied] = (k - lower) / len(tied)
    row_sum = float(weights.sum())
    assert abs(row_sum - k) <= 1e-9, (k, lower, len(tied), row_sum)
    return weights, shell, row_sum


def tie_weights(points: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    """Return independent tie-breaking probabilities and the tie-inclusive set."""
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or len(points) < 2:
        raise ValueError("points must be a finite N x 2 array with N >= 2")
    if not np.isfinite(points).all():
        raise ValueError("points must be finite")
    k = min(max(1, int(k)), len(points) - 1)
    probabilities = np.zeros((len(points), len(points)), dtype=float)
    rows = []
    for start in range(0, len(points), 1024):
        stop = min(len(points), start + 1024)
        delta = np.sqrt(
            ((points[start:stop, None, :] - points[None, :, :]) ** 2).sum(axis=2)
        )
        for local, row in enumerate(delta):
            row[start + local] = np.inf
            weights, shell, _ = _distance_weights(row, k)
            probabilities[start + local] = weights
            rows.append(shell)
    assert np.allclose(probabilities.sum(axis=1), k, rtol=0.0, atol=1e-9)
    return probabilities, np.asarray(rows, dtype=object)


def representation_weights(
    fingerprints: np.ndarray, k: int
) -> tuple[np.ndarray, list[np.ndarray]]:
    """Tie-breaking probabilities and tie-inclusive sets in Morgan space."""
    fingerprints = np.asarray(fingerprints)
    bits = (fingerprints != 0).astype(np.float32)
    counts = bits.sum(axis=1).astype(np.float64)
    probabilities = np.zeros((len(bits), len(bits)), dtype=float)
    shells: list[np.ndarray] = []
    k = min(max(1, int(k)), len(bits) - 1)
    for start in range(0, len(bits), 1024):
        stop = min(len(bits), start + 1024)
        shared = np.rint(bits[start:stop] @ bits.T).astype(np.float64)
        union = counts[start:stop, None] + counts[None, :] - shared
        distance = 1.0 - shared / np.where(union > 0, union, 1.0)
        for local, row in enumerate(distance):
            row[start + local] = np.inf
            weights, shell, _ = _distance_weights(row, k)
            probabilities[start + local] = weights
            shells.append(shell)
    assert np.allclose(probabilities.sum(axis=1), k, rtol=0.0, atol=1e-9)
    return probabilities, shells


def grid_weights(cells: np.ndarray, k: int) -> tuple[np.ndarray, list[np.ndarray]]:
    """Tie-breaking probabilities and tie-inclusive sets in lattice space."""
    cells = np.asarray(cells, dtype=float)
    k = min(max(1, int(k)), len(cells) - 1)
    probabilities = np.zeros((len(cells), len(cells)), dtype=float)
    shells: list[np.ndarray] = []
    for start in range(0, len(cells), 1024):
        stop = min(len(cells), start + 1024)
        delta = np.sqrt(
            ((cells[start:stop, None, :] - cells[None, :, :]) ** 2).sum(axis=2)
        )
        for local, row in enumerate(delta):
            row[start + local] = np.inf
            weights, shell, _ = _distance_weights(row, k)
            probabilities[start + local] = weights
            shells.append(shell)
    assert np.allclose(probabilities.sum(axis=1), k, rtol=0.0, atol=1e-9)
    return probabilities, shells


def expected_activity_error(
    activity: np.ndarray, probabilities: np.ndarray, k: int
) -> np.ndarray:
    """Expected mean absolute pChEMBL difference under independent tie-breaking."""
    activity = np.asarray(activity, dtype=float)
    error = np.abs(activity[:, None] - activity[None, :])
    return np.einsum("ij,ij->i", probabilities, error) / k


def fixed_metrics(
    activity: np.ndarray,
    fingerprints: np.ndarray,
    cells: np.ndarray,
    reference_probabilities: np.ndarray | None = None,
    reference_shells: Sequence[np.ndarray] | None = None,
    reference_sets: Sequence[np.ndarray] | None = None,
    k: int = 10,
) -> dict[str, np.ndarray]:
    """Compute retained and fixed-k recall, neighbour sizes, and errors."""
    activity = np.asarray(activity, dtype=float)
    fingerprints = np.asarray(fingerprints)
    cells = np.asarray(cells, dtype=float)
    if reference_probabilities is None or reference_shells is None:
        reference_probabilities, reference_shells = representation_weights(
            fingerprints, k
        )
    reference_shells = list(reference_shells)
    grid_probabilities, grid_shells = grid_weights(cells, k)
    retained_grid_sets = list(grid_shells)
    retained_reference_sets = (
        list(reference_sets)
        if reference_sets is not None
        else list(reference_shells)
    )
    retained_recall = neighbor_recall(
        retained_reference_sets, retained_grid_sets
    )
    fixed_recall = np.einsum(
        "ij,ij->i", reference_probabilities, grid_probabilities
    ) / k
    return {
        "retained_recall": retained_recall,
        "fixed_recall": fixed_recall,
        "retained_e_grid": np.array(
            [np.abs(activity[group] - activity[i]).mean()
             for i, group in enumerate(retained_grid_sets)]
        ),
        "retained_e_fp": np.array(
            [np.abs(activity[group] - activity[i]).mean()
             for i, group in enumerate(retained_reference_sets)]
        ),
        "fixed_e_grid": expected_activity_error(
            activity, grid_probabilities, k
        ),
        "fixed_e_fp": expected_activity_error(
            activity, reference_probabilities, k
        ),
        "grid_neighbor_count": np.array(
            [len(group) for group in retained_grid_sets], dtype=int
        ),
        "fp_neighbor_count": np.array(
            [len(group) for group in retained_reference_sets], dtype=int
        ),
    }
