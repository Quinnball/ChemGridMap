from __future__ import annotations

import sys
from pathlib import Path

import pytest


SOURCE_DIR = Path(__file__).resolve().parent / "figures/source"
sys.path.insert(0, str(SOURCE_DIR))

from fig6 import _pale_legend_text  # noqa: E402


def test_pale_legend_comes_from_capture_report():
    report = {
        "checks": {
            "legend": (
                "Inactive ≤6 uncertain Medium (6, 7) uncertain "
                "Pale: uncertain (class support below 0.8)."
            )
        }
    }
    assert _pale_legend_text(report) == (
        "Pale: uncertain (class support below 0.8)."
    )


def test_pale_legend_requires_complete_capture_sentence():
    with pytest.raises(AssertionError, match="Pale:"):
        _pale_legend_text(
            {"checks": {"legend": "Inactive ≤6 uncertain Medium (6, 7) uncertain"}}
        )
