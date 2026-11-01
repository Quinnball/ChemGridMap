import numpy as np
import pandas as pd
import pytest

import e7_analysis
from e9_analysis import (
    calibration_bin_index,
    dataset_bias_stats,
    random_effect_details,
    two_stage_weights,
)


def test_binomial_variance_positive_at_zero_and_one():
    for agreement in (0.0, 1.0):
        frame = pd.DataFrame(
            {
                "target": ["A"] * 20,
                "agree": [bool(agreement)] * 20,
                "p_displayed": [0.5] * 20,
            }
        )
        result = random_effect_details(frame, minimum_predictions=20)
        assert result["n_datasets"] == 1
        assert np.isfinite(result["binomial_max_weight_fraction"])
        epsilon = 0.25 / 20 / 20
        assert epsilon > 0


def test_two_stage_weights_is_e7_function_itself():
    assert two_stage_weights is e7_analysis.two_stage_weights


def test_dataset_median_respects_minimum_predictions():
    frame = pd.DataFrame(
        {
            "target": ["A"] * 9 + ["B"] * 10,
            "p_displayed": [0.0] * 9 + [0.5] * 10,
            "agree": [1.0] * 9 + [0.0] * 10,
        }
    )
    result = dataset_bias_stats(frame, minimum_predictions=10)
    assert result["dataset_median_n"] == 1
    assert result["bias_median"] == 0.5


def test_probability_bin_boundaries_include_right_endpoint():
    assert calibration_bin_index([0.0, 0.2, 0.4, 0.8, 1.0]).tolist() == [
        0,
        1,
        2,
        4,
        4,
    ]


def test_saved_old_method_regressions():
    from e9_analysis import OUT

    if not (OUT / "old_dl_diagnostic.csv").is_file():
        pytest.skip("Run E9 before the archived output regression.")
    section1 = pd.read_csv(OUT / "old_dl_diagnostic.csv")
    section4 = pd.read_csv(OUT / "rule2_recomputed.csv")
    time_tests = pd.read_csv(OUT / "time_split_tests.csv")
    assert section1.absolute_error.max() <= 1e-9
    assert section4.old_dl_absolute_error.max() <= 1e-9
    assert time_tests.reference_max_absolute_error.max() <= 1e-9
