"""E5b grid, NPMLE, calibration and t-variance tests."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from chemgridmap.evidence import class_probabilities
from e5b_models import (
    entry_values,
    grid_prediction,
    student_t_scale,
)


def test_gaussian_grid_matches_analytic_g():
    values = np.array([6.91, 7.13])
    sigma = 0.74
    support, predictive, posterior_mean, _ = grid_prediction(
        values, sigma, np.inf, None
    )
    m = len(values)
    ybar = float(values.mean())
    expected_support = class_probabilities(
        [ybar], [sigma / np.sqrt(m)], np.inf, 6.0, 7.0
    )[0]
    expected_predictive = class_probabilities(
        [ybar], [sigma * np.sqrt(1.0 + 1.0 / m)], np.inf, 6.0, 7.0
    )[0]
    assert np.allclose(support, expected_support, atol=1e-4, rtol=0)
    assert np.allclose(
        predictive, expected_predictive, atol=1e-4, rtol=0
    )
    assert posterior_mean == pytest.approx(ybar, abs=1e-4)


def test_npmle_training_helper_excludes_test_entries():
    frame = pd.DataFrame(
        {
            "entry": ["a", "a", "b", "b", "c"],
            "document_value": [5.0, 5.2, 7.5, 7.7, 9.0],
        }
    )
    values = entry_values(frame, excluded_entries={"b"})
    assert len(values) == 2
    assert np.allclose(values[0], [5.0, 5.2])
    assert np.allclose(values[1], [9.0])


def test_student_t_scale_has_requested_variance():
    sigma = 0.83
    for df in (3, 4, 8):
        scale = student_t_scale(sigma, df)
        variance = scale**2 * df / (df - 2.0)
        assert variance == pytest.approx(sigma**2)
        assert stats.t.var(df, scale=scale) == pytest.approx(sigma**2)


def test_calibration_training_excludes_held_dataset():
    from e5b_analysis import (
        calibration_features,
        crossfit_gcal,
        fit_calibrator,
    )

    rows = []
    for target in ("a", "b", "c"):
        for probability, agreement in (
            (0.9, 1),
            (0.8, 0),
            (0.6, 0),
            (0.2, 1),
        ):
            shown = 2 if probability >= 0.7 else 0
            rows.append(
                {
                    "split_design": "entry",
                    "target": target,
                    "endpoint": "IC50",
                    "fold": 0,
                    "entry": f"{target}{probability}{agreement}",
                    "held_out_unit": f"D{target}{probability}{agreement}",
                    "p0": 1.0 if shown == 0 else 0.0,
                    "p1": 0.0,
                    "p2": 1.0 if shown == 2 else 0.0,
                    "p_displayed": probability,
                    "support": probability,
                    "display_status": (
                        "supported"
                        if probability >= 0.8
                        else "uncertain"
                    ),
                    "displayed": shown,
                    "m": 1,
                    "agree": bool(agreement),
                    "method": "G",
                }
            )
    frame = pd.DataFrame(rows)
    held_target = "a"
    train = frame.loc[frame.target != held_target]
    first = fit_calibrator(train)
    second = fit_calibrator(train.assign(agree=~train.agree))
    assert calibration_features(frame.iloc[[0]]).shape[0] == 1
    assert not np.allclose(first, second)
    original, original_fits = crossfit_gcal("entry", frame, frame)
    changed = frame.copy()
    changed.loc[changed.target == held_target, "agree"] = (
        ~changed.loc[changed.target == held_target, "agree"]
    )
    changed.loc[changed.target == held_target, "p_displayed"] = 0.35
    _output, changed_fits = crossfit_gcal("entry", changed, changed)
    original_row = original_fits[0]
    changed_row = changed_fits[0]
    original_coefficients = [
        original_row[key]
        for key in original_row
        if key.startswith("coefficient_")
    ]
    changed_coefficients = [
        changed_row[key]
        for key in changed_row
        if key.startswith("coefficient_")
    ]
    assert original_row["target"] == held_target
    assert original_coefficients == changed_coefficients
    assert len(original) == len(frame)
