import numpy as np
import pandas as pd

import e7_analysis
from e5d_models import add_rule_states, b_straddle, corrected_two_stage_weights


def test_rule_p_only_removes_supported_predictions():
    frame = pd.DataFrame(
        {
            "support": [0.9, 0.9, 0.7, 0.7],
            "flag_touched": [False, True, False, True],
        }
    )
    result = add_rule_states(frame)
    assert result.g_status.tolist() == [
        "supported",
        "supported",
        "uncertain",
        "uncertain",
    ]
    assert result.p_status.tolist() == [
        "supported",
        "uncertain_touched",
        "uncertain_low_support",
        "uncertain_touched",
    ]
    assert not (
        (result.g_status == "uncertain")
        & (result.p_status == "supported")
    ).any()


def test_two_stage_bootstrap_uses_e7_r1_dataset_resampling():
    assert (
        corrected_two_stage_weights.__globals__["two_stage_weights"]
        is e7_analysis.two_stage_weights
    )
    frame = pd.DataFrame(
        {
            "target": ["A"] * 4 + ["B"] * 4,
            "endpoint": ["IC50"] * 8,
            "entry": ["a1", "a1", "a2", "a2", "b1", "b1", "b2", "b2"],
        }
    )
    weights = corrected_two_stage_weights(
        frame, repetitions=300, seed=20261105
    )
    totals = np.column_stack(
        [
            weights[:, :4].sum(axis=1),
            weights[:, 4:].sum(axis=1),
        ]
    )
    assert np.any(totals[:, 0] == 0)
    assert np.any(totals[:, 1] == 0)
    assert np.any(totals[:, 0] != totals[:, 1])
    assert np.any(np.diff(totals, axis=0) != 0)


def test_b_straddle_toy_examples():
    assert b_straddle([5.0, 6.0, 6.5]) == 1
    assert b_straddle([5.0, 7.0]) == 0
    assert b_straddle([6.1, 6.9]) == 1
    assert b_straddle([]) == 1
