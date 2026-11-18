"""Regression checks for the E12c stereochemistry sensitivity."""
from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd
import pytest

SOURCE = Path(__file__).resolve().parent
sys.path.insert(0, str(SOURCE))
from e12c_core_stereo import ROOT, S, compute_settings

REQUIRED_INPUTS = (
    S / "output/E1B/core_group_members.csv",
    S / "output/T1/core_timing_by_group.csv",
    S / "output/E12b/core_chirality.csv",
)
needs_validation_data = pytest.mark.skipif(
    not all(path.is_file() for path in REQUIRED_INPUTS),
    reason="Validation data are not restored.",
)


@needs_validation_data
def test_first_two_settings_reproduce_e12b():
    grouped, _, checks = compute_settings()
    assert checks["e12b_two_setting_regression"]
    assert len(grouped) == 36


@needs_validation_data
def test_no_stereochemistry_counts():
    _, _, checks = compute_settings()
    assert checks["no_chirality_strict_majority"] == {
        "coverage": 470,
        "most_frequent": 272,
        "fmcs": 480,
    }
    assert checks["changed_group_rule_pairs"] == 9
    assert checks["changed_pairs_by_rule"] == {
        "coverage": 0,
        "most_frequent": 9,
        "fmcs": 0,
    }


def test_repository_package_and_archived_output():
    import chemgridmap

    assert Path(chemgridmap.__file__).resolve().is_relative_to(ROOT / "src")
    output = S / "output/E12c/core_stereo_settings.csv"
    if output.exists():
        frame = pd.read_csv(output)
        assert len(frame) == 36
