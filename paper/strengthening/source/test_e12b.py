"""Small regression and unit tests for E12b."""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "paper")]
import chemgridmap

assert ROOT in Path(chemgridmap.__file__).resolve().parents

from e12b_analysis import (
    build_ece_reference,
    build_random_effect_alternatives,
    build_scores_h4_h7,
    build_time_populations,
    holm_adjust,
    read_predictions,
    score_probabilities,
    two_sided_wilcoxon,
)


def test_score_functions_on_toy_probabilities():
    frame = pd.DataFrame({
        "p0": [0.8, 0.2],
        "p1": [0.1, 0.3],
        "p2": [0.1, 0.5],
        "observed": [0, 2],
    })
    assert np.allclose(score_probabilities(frame, "brier"), [0.06, 0.04 + 0.09 + 0.25])
    expected = -np.log(np.clip([0.8, 0.5], 1e-15, 1.0))
    assert np.allclose(score_probabilities(frame, "log_loss_1e-15"), expected)


def test_two_sided_and_holm():
    values = [0.2, 0.3, 0.4, 0.5, 0.6]
    _, two_p, n = two_sided_wilcoxon(values)
    assert n == 5 and 0 < two_p < 1
    assert np.allclose(holm_adjust([0.01, 0.04, 0.03]), [0.03, 0.06, 0.06])


def test_section_1_regression():
    if not (ROOT / "paper/strengthening/output/E2/predictions.csv").is_file():
        pytest.skip("Run E2 before the section 1 regression.")
    _, regression = build_scores_h4_h7(read_predictions())
    assert regression["passed"]


def test_section_2_regression():
    if not (ROOT / "paper/strengthening/output/E2/hypothesis_tests.csv").is_file():
        pytest.skip("Run E2 before the section 2 regression.")
    _, regression = build_random_effect_alternatives()
    assert regression["passed"]


def test_section_4_regression():
    if not (ROOT / "paper/strengthening/output/E2/predictions.csv").is_file():
        pytest.skip("Run E2 before the section 4 regression.")
    _, regression = build_ece_reference(read_predictions())
    assert regression["passed"]


def test_time_population_regression():
    if not (ROOT / "paper/strengthening/output/E7/time_split_by_dataset.csv").is_file():
        pytest.skip("Run E7 before the time-population regression.")
    _, regression = build_time_populations()
    assert regression["passed"]
