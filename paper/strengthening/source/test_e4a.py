"""E4a 的公式单元测试与冻结输出回归检查。"""
from __future__ import annotations

import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True

import numpy as np
import pandas as pd
import pytest

S = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(S / "source")]

from e4a_analysis import (
    BOOTSTRAP_SEED,
    CURVE_MIN_PREDICTIONS,
    distance_score,
    distance_scores,
    displayed_ece,
)


def test_distance_score_by_displayed_class():
    assert distance_score(5.5, 0) == pytest.approx(0.5)
    assert distance_score(6.5, 1) == pytest.approx(0.5)
    assert distance_score(7.5, 2) == pytest.approx(0.5)
    assert distance_score(6.1, 1) == pytest.approx(0.1)
    assert distance_scores([5.5, 6.5, 7.5], [0, 1, 2]) == pytest.approx(
        [0.5, 0.5, 0.5]
    )


def test_displayed_ece_uses_shown_probability_and_agreement():
    score = np.array([0.1, 0.2, 0.8, 1.0])
    agree = np.array([0.0, 0.0, 1.0, 0.0])
    # 每一箱贡献 0.1、0.2、0.2、1.0，按四条命中平均。
    assert displayed_ece(score, agree)[0] == pytest.approx(1.5 / 4)


def test_e4a_outputs_and_frozen_regression():
    out = S / "output/E4"
    qa = S / "qa/E4/e4a_checks.json"
    expected = {
        "comparators_posthoc.csv",
        "calibration_displayed.csv",
        "calibration_curve_displayed.csv",
    }
    if not qa.is_file() or not expected.issubset(
        {path.name for path in out.iterdir()} if out.is_dir() else set()
    ):
        pytest.skip("Run e4a_analysis.py first.")
    assert expected.issubset({path.name for path in out.iterdir()})

    checks = json.loads(qa.read_text())
    assert checks["passed"] is True
    assert checks["included_datasets"] == 52
    assert checks["maximum_metric_auc_error"] <= 1e-12
    assert checks["seeds"] == {
        "bootstrap": BOOTSTRAP_SEED,
        "bootstrap_repetitions": 2000,
    }
    assert checks["forest_g_ci_counts"] == {
        "IC50": {"n_datasets": 29, "g_ci_fully_above_half": 24},
        "Ki": {"n_datasets": 23, "g_ci_fully_above_half": 19},
    }
    assert checks["ki_excluding_original_targets"]["n_datasets"] == 21
    assert checks["ki_excluding_original_targets"]["median_auc"] == pytest.approx(
        0.7468357726964946, abs=1e-12
    )
    assert checks["ki_excluding_original_targets"]["wilcoxon_p"] == pytest.approx(
        4.76837158203125e-07, rel=1e-12
    )
    summary = pd.read_csv(
        S / "output/E4/e4a_comparison_summary.csv",
        float_precision="round_trip",
    )
    b_sd = summary.loc[
        (summary.row_type == "auc_ci_above_half")
        & (summary.comparison == "B-sd")
        & (summary.endpoint == "Ki")
    ].iloc[0]
    assert b_sd.n_datasets == 23
    assert b_sd.n_auc_valid == 22


@pytest.mark.parametrize(
    "name",
    [
        "comparators_posthoc.csv",
        "calibration_displayed.csv",
        "calibration_curve_displayed.csv",
    ],
)
def test_required_outputs_are_post_hoc(name):
    if not (S / "output/E4" / name).is_file():
        pytest.skip("Run e4a_analysis.py first.")
    frame = pd.read_csv(
        S / "output/E4" / name,
        float_precision="round_trip",
    )
    assert not frame.empty
    assert set(frame.analysis_role) == {"post hoc"}


def test_calibration_curve_rule_and_bin_counts():
    if not (S / "output/E4/calibration_curve_displayed.csv").is_file():
        pytest.skip("Run e4a_analysis.py first.")
    curve = pd.read_csv(
        S / "output/E4/calibration_curve_displayed.csv",
        float_precision="round_trip",
    )
    assert set(curve.endpoint) == {"IC50", "Ki"}
    assert set(curve.model) == {"G", "B0"}
    assert curve.groupby(["endpoint", "model"]).bin.count().eq(5).all()
    assert (curve.minimum_predictions == CURVE_MIN_PREDICTIONS).all()
    assert (
        curve.eligible
        == (curve.n_predictions >= CURVE_MIN_PREDICTIONS)
    ).all()
    assert curve.loc[curve.eligible, "mean_probability"].notna().all()
    assert curve.loc[curve.eligible, "observed_agreement"].notna().all()
