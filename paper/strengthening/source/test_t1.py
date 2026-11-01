"""Fast checks for the T1 helpers; the full timing benchmark is not rerun."""

from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import t1_scale

sys.path.insert(0, str(t1_scale.CGM / "src"))


def test_scaling_subsets_are_unique_and_nested():
    herg = [f"H{index:03d}" for index in range(10)]
    extra = [f"X{index:03d}" for index in range(4)]
    subsets = t1_scale.nested_size_selection(herg, extra, (4, 7, 10, 11, 12, 13))
    assert [len(subsets[size]) for size in subsets] == [4, 7, 10, 11, 12, 13]
    ordered = sorted(subsets)
    for smaller, larger in zip(ordered, ordered[1:]):
        assert set(subsets[smaller]).issubset(subsets[larger])
    assert subsets[13][-3:] == extra[:3]


def test_fmcs_timeout_is_kept_in_the_call_record(monkeypatch):
    def timed_out(*_args, **_kwargs):
        return {
            "core_smarts": "[#6]",
            "core_heavy_atoms": 1,
            "coverage_count": 2,
            "coverage_fraction": 1.0,
            "timed_out": True,
            "seconds": 1.0,
        }

    monkeypatch.setattr("chemgridmap.cores.mcs_core", timed_out)
    result = t1_scale.measure_fmcs([], timeout=1)
    assert result["timed_out"] is True
    assert result["coverage_count"] == 2
    assert result["seconds"] >= 0


def test_grouping_limit_error_is_recognized():
    assert t1_scale.expected_limit_error(
        {
            "size_entries": 10001,
            "failure": {
                "message": (
                    "Exact grouping is limited to 10,000 molecules; "
                    "compare an explicit subset."
                )
            },
        }
    )
    assert not t1_scale.expected_limit_error(
        {"size_entries": 10000, "failure": {"message": "another error"}}
    )
