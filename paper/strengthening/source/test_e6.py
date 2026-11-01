from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from e6_core import (
    fixed_metrics,
    grid_weights,
    representation_weights,
    tie_weights,
)
from reliability_analysis import load_target

ROOT = Path(__file__).resolve().parents[3]
CGM = ROOT


def test_fixed_weights_with_strictly_closer_neighbor_and_tie():
    cells = np.array([[0.0, 0.0], [0.5, 0.0], [1.0, 0.0], [-1.0, 0.0]])
    probabilities, shells = tie_weights(cells, k=2)
    assert shells[0].tolist() == [1, 2, 3]
    assert probabilities[0, 1] == 1.0
    assert probabilities[0, 2] == 0.5
    assert probabilities[0, 3] == 0.5
    assert probabilities[0].sum() == 2.0


def test_each_weight_row_sums_to_k_with_and_without_ties():
    cells = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]])
    fingerprints = np.array([
        [1, 0, 0, 0],
        [0, 1, 1, 0],
        [1, 1, 1, 0],
        [1, 1, 1, 1],
    ])
    for probabilities in (
        tie_weights(cells, k=2)[0],
        grid_weights(cells, k=2)[0],
        representation_weights(fingerprints, k=2)[0],
    ):
        np.testing.assert_allclose(
            probabilities.sum(axis=1), 2.0, rtol=0.0, atol=1e-9
        )


def test_projection_without_ties_has_equal_recall_definitions():
    activity = np.array([0.0, 1.0, 2.0, 3.0])
    fingerprints = np.array([
        [1, 0, 0, 0],
        [0, 1, 1, 0],
        [1, 1, 1, 0],
        [1, 1, 1, 1],
    ])
    cells = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]])
    metrics = fixed_metrics(activity, fingerprints, cells, k=2)
    assert np.allclose(metrics["retained_recall"], metrics["fixed_recall"])
    assert np.allclose(metrics["retained_e_grid"], metrics["fixed_e_grid"])


def test_ca_ii_seed42_matches_opus_independent_fixed_recall():
    layout_path = (
        CGM
        / "paper/output/reliability_v1/layouts/chembl205"
        / "coordinates_seed_42.csv"
    )
    if not layout_path.is_file():
        return pytest.skip("Restore the frozen reliability layouts first.")
    data, fingerprints = load_target("chembl205")
    activity = data.activity_pchembl.to_numpy(float)
    frame = pd.read_csv(layout_path)
    if data.molecule_identity_key.tolist() != frame.molecule_identity_key.tolist():
        raise AssertionError("CA II seed 42 coordinate identities do not match")
    expected = {
        "projection": 0.6621,
        "minimum_displacement": 0.4513,
        "neighbor_refined": 0.5359,
    }
    columns = {
        "projection": ("projection_x", "projection_y"),
        "minimum_displacement": ("grid_col", "grid_row"),
        "neighbor_refined": ("refined_grid_col", "refined_grid_row"),
    }
    for layout, (col, row) in columns.items():
        metrics = fixed_metrics(
            activity, fingerprints, frame[[col, row]].to_numpy(float), k=10
        )
        observed = float(metrics["fixed_recall"].mean())
        assert abs(observed - expected[layout]) <= 0.002, (
            layout,
            observed,
            expected[layout],
        )
