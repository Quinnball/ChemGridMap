"""Focused tests for E7 analysis primitives."""
from __future__ import annotations

import numpy as np
import pandas as pd

from e7_analysis import (
    connected_components,
    dataset_counts,
    entry_equal_weights,
    pooled_stratified_auc,
    recompute_hypotheses,
    two_stage_weights,
    weighted_auc,
)


def test_two_stage_weights_resample_datasets_and_entries():
    frame = pd.DataFrame({
        "target": ["A"] * 4 + ["B"] * 4,
        "endpoint": ["IC50"] * 8,
        "entry": ["a1", "a1", "a2", "a2", "b1", "b1", "b2", "b2"],
    })
    weights = two_stage_weights(frame, seed=20261103, repetitions=300)
    dataset_totals = np.column_stack([
        weights[:, :4].sum(axis=1),
        weights[:, 4:].sum(axis=1),
    ])
    assert len(np.unique(dataset_totals, axis=0)) > 1
    assert np.any(dataset_totals == 0)
    assert np.array_equal(weights[:, 0], weights[:, 1])
    assert np.array_equal(weights[:, 2], weights[:, 3])
    assert np.array_equal(weights[:, 4], weights[:, 5])
    assert np.array_equal(weights[:, 6], weights[:, 7])


def test_dataset_counts_are_multinomial():
    n_datasets = 2
    repetitions = 100
    counts = dataset_counts(
        n_datasets=n_datasets,
        rng=np.random.default_rng(1),
        repetitions=repetitions,
    )
    assert counts.shape == (repetitions, n_datasets)
    assert np.all(counts.sum(axis=1) == n_datasets)
    assert np.var(counts[:, 0]) > 0
    assert np.any(counts[:, 0] != counts[:, 1])


def test_stratified_auc_equals_plain_auc_for_one_dataset():
    rows = pd.DataFrame({
        "n_agree": [4],
        "n_disagree": [5],
        "auc": [0.75],
        "stratum_weight": [20],
    })
    assert pooled_stratified_auc(rows) == 0.75


def test_weighted_auc_matches_mann_whitney_without_weights():
    correct = np.array([1, 0, 1, 0, 1], dtype=bool)
    score = np.array([0.9, 0.8, 0.7, 0.4, 0.3])
    weights = np.ones((1, 5))
    from e1a_statistics import auc_mann_whitney

    assert np.isclose(
        weighted_auc(correct, score, weights)[0],
        auc_mann_whitney(correct, score),
    )


def test_connected_components():
    components = connected_components(
        ["a", "b", "c", "d"],
        [("a", "b"), ("b", "c")],
    )
    assert components == [["a", "b", "c"], ["d"]]


def test_entry_equal_weights_sum_to_entry_count():
    frame = pd.DataFrame({
        "target": ["A", "A", "A", "A"],
        "endpoint": ["IC50"] * 4,
        "entry": ["a1", "a1", "a2", "a3"],
    })
    assert np.isclose(entry_equal_weights(frame).sum(), frame.entry.nunique())


def test_h7_excludes_unpaired_zero_differences():
    frame = pd.DataFrame({
        "target": ["paired", "unpaired"],
        "endpoint": ["IC50", "IC50"],
        "h3_difference": [0.2, 0.1],
        "h4_difference": [-0.1, -0.1],
        "h5_difference": [0.1, 0.1],
        "h6_difference": [0.1, 0.1],
        "h7_difference": [0.25, 0.0],
        "h8_difference": [0.1, 0.1],
    })

    result = recompute_hypotheses(
        frame,
        "synthetic",
        h7_eligible={"IC50": {"paired"}, "Ki": set()},
    )
    h7 = result.loc[
        (result.hypothesis == "H7") & (result.endpoint == "IC50")
    ].iloc[0]

    assert h7.n_datasets == 1
    assert h7.median_difference == 0.25
