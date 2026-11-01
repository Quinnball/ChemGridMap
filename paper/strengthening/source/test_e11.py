"""Tests for the sealed E11 validation implementation."""
from __future__ import annotations

import numpy as np
import pytest

from e11_analysis import (
    E5C,
    hypothesis_result,
    run_regression,
)
from e7_analysis import holm_adjust


def test_e11_regression_matches_e5c() -> None:
    if not (E5C / "predictions_e5c.csv").is_file():
        pytest.skip("Run E5c before the E11 regression.")
    result = run_regression()
    assert result["passed"]
    assert result["n_matched"] == result["n_actual"] == result["n_expected"]
    assert result["maximum_probability_error"] <= 1e-12
    assert result["agree_identical"]
    assert result["flag_touched_identical"]
    assert sum(
        row["n_flag_touched"] > 0
        for row in result["selected_datasets"]
    ) >= 2


def test_holm_adjustment() -> None:
    adjusted = holm_adjust([0.01, 0.04, 0.03])
    assert np.allclose(adjusted, [0.03, 0.06, 0.06], atol=1e-15)
    with_nan = holm_adjust([0.01, np.nan, 0.03])
    assert np.allclose(with_nan[[0, 2]], [0.02, 0.03], atol=1e-15)
    assert np.isnan(with_nan[1])


def test_fewer_than_eight_datasets_skips_test() -> None:
    result = hypothesis_result(np.arange(7, dtype=float) + 1)
    assert result["n_datasets"] == 7
    assert not result["test_performed"]
    assert np.isnan(result["p_value"])
