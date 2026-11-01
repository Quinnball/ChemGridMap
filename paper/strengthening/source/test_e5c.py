"""Tests for E5c offset-aware support and leakage controls."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from e5b_models import grid_prediction
from e5c_models import (
    delta_document_map,
    unit_deltas,
    unit_sigmas,
    weighted_gaussian_prediction,
    weighted_grid_prediction,
)


def test_zero_delta_g_f_equals_g():
    values = np.array([6.91, 7.13, 8.00])
    sigma = 0.74
    support, predictive, mean, variance = weighted_gaussian_prediction(
        values, np.full(len(values), sigma**2), sigma
    )
    expected_support = grid_prediction(values, sigma, np.inf, None)[0]
    expected_predictive = grid_prediction(values, sigma, np.inf, None)[1]
    assert np.allclose(support, expected_support, atol=1e-4, rtol=0)
    assert np.allclose(
        predictive, expected_predictive, atol=1e-4, rtol=0
    )
    assert mean == pytest.approx(values.mean(), abs=1e-12)
    assert variance == pytest.approx(sigma**2 / len(values), abs=1e-12)


def test_document_mark_uses_training_entries_only():
    observations = pd.DataFrame(
        {
            "document_a": ["train-a"] * 3 + ["late"] * 3,
            "document_b": ["train-b"] * 3 + ["held-only"] * 3,
            "entry": [
                "train-entry-1",
                "train-entry-2",
                "train-entry-3",
                "held-entry-1",
                "held-entry-2",
                "held-entry-3",
            ],
            "difference": [1.5, 1.6, 1.4, 1.5, 1.6, 1.4],
        }
    )
    marked = delta_document_map(
        observations,
        training_entries={
            "train-entry-1",
            "train-entry-2",
            "train-entry-3",
        },
    )
    assert marked == {"train-a": 1.5, "train-b": 1.5}
    assert "late" not in marked
    assert "held-only" not in marked


def test_weighted_merge_matches_hand_calculation():
    values = np.array([5.0, 8.0])
    sigmas = np.array([1.0, 2.0])
    base_sigma = 0.5
    variances = sigmas**2
    weights = 1.0 / variances
    expected_mean = np.sum(weights * values) / weights.sum()
    expected_variance = 1.0 / weights.sum()
    _, _, mean, variance = weighted_gaussian_prediction(
        values, variances, base_sigma
    )
    assert mean == pytest.approx(expected_mean, abs=1e-12)
    assert variance == pytest.approx(expected_variance, abs=1e-12)


def test_weighted_grid_without_marks_equals_g_np_t():
    values = np.array([6.91, 7.13])
    sigma = 0.74
    prior = np.full(601, 1.0 / 601)
    sigmas = unit_sigmas(np.zeros(len(values)), sigma)
    actual = weighted_grid_prediction(
        values, sigmas, sigma, prior, df=4.0
    )
    expected = grid_prediction(values, sigma, 4.0, prior)
    assert np.allclose(actual[0], expected[0], atol=1e-12, rtol=0)
    assert np.allclose(actual[1], expected[1], atol=1e-12, rtol=0)
    assert actual[2] == pytest.approx(expected[2], abs=1e-12)


def test_unit_delta_uses_largest_matched_document_offset():
    units = pd.DataFrame(
        {
            "documents": ["a;b", "c", "d"],
            "document_value": [5.0, 7.0, 9.0],
        }
    )
    delta, touched = unit_deltas(
        units, {"a": 1.2, "b": 2.1, "c": 1.6}
    )
    assert np.allclose(delta, [2.1, 1.6, 0.0])
    assert touched.tolist() == [True, True, False]
