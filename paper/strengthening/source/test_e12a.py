"""Regression and unit tests for the E12a analysis."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import e12a_rule_calibration as e12a


def test_chemgridmap_imports_from_repository():
    assert Path(__file__).resolve().parents[3] in e12a.CHEMGRIDMAP_FILE.parents


def test_rule_states_form_a_partition():
    frame = pd.DataFrame({
        "support": [0.9, 0.7, 0.9, 0.2],
        "flag_touched": [False, False, True, True],
    })
    result = e12a.add_rule_states(frame)
    assert result.g_status.tolist() == [
        "supported", "uncertain", "supported", "uncertain"
    ]
    assert result.p_status.tolist() == [
        "supported",
        "uncertain_low_support",
        "uncertain_touched",
        "uncertain_touched",
    ]
    assert result.p_status.eq("supported").sum() == 1
    assert result.p_status.ne("supported").sum() == 3


def test_two_sided_wilcoxon_and_holm_adjustment():
    values = np.arange(1.0, 9.0)
    _, p_value, n = e12a.two_sided_wilcoxon(values)
    expected = e12a.stats.wilcoxon(
        values, alternative="two-sided", zero_method="wilcox"
    ).pvalue
    assert n == 8
    assert p_value == pytest.approx(expected)
    adjusted = e12a.holm_adjust([0.01, 0.04, 0.03])
    np.testing.assert_allclose(adjusted, [0.03, 0.06, 0.06])


def test_component_values_average_only_eligible_datasets():
    differences = pd.DataFrame({
        "target": ["A", "B", "C"],
        "eligible": [True, True, False],
        "difference": [0.1, 0.3, 0.9],
    })
    membership = pd.DataFrame({
        "endpoint": ["IC50"] * 3,
        "component": ["c1", "c1", "c2"],
        "target": ["A", "B", "C"],
    })
    result = e12a.component_values(differences, membership)
    assert result.component.tolist() == ["c1"]
    assert result.difference.iloc[0] == pytest.approx(0.2)
    assert result.n_eligible_datasets.iloc[0] == 2


def test_rule_test_counts_match_manuscript():
    if not (e12a.E5C / "predictions_e5c.csv").is_file():
        pytest.skip("Run E5c before the rule-count regression.")
    _, designs = e12a.load_designs()
    members = pd.read_csv(e12a.E7 / "shared_document_cluster_members.csv")
    tests = e12a.rule_tests(designs, members)
    e12a.rule_test_regression(tests)


def test_calibration_output_regresses_to_e9():
    output = e12a.OUT / "calibration_all_designs.csv"
    if not output.exists():
        pytest.skip("run e12a_rule_calibration.py before the output regression")
    calibration = pd.read_csv(output)
    result = e12a.calibration_regression(calibration)
    assert result["n_overlapping_rows"] == 200
    assert result["max_absolute_error"] <= 1e-12
