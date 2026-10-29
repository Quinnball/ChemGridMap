"""Focused checks for the E4e release/evaluation alignment."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from e4e_analysis import OUT, RESULT, SOURCE, WORKSPACE


IS_REPOSITORY_COPY = SOURCE.parents[1].name == "paper"
REQUIRED_SOURCE = {
    "f0_build.py",
    "test_f0.py",
    "e4a_analysis.py",
    "e4b_analysis.py",
    "e4c_analysis.py",
    "e4d_analysis.py",
    "e4e_analysis.py",
    "test_e4a.py",
    "test_e4b.py",
    "test_e4c.py",
    "test_e4d.py",
    "test_e4e.py",
}


def require_e4e_output():
    if not RESULT.is_file():
        pytest.skip("Run e4e_analysis.py first.")


@pytest.mark.skipif(
    IS_REPOSITORY_COPY,
    reason="Analysis outputs are delivered separately from the source repository copy.",
)
def test_regression_output_has_complete_scope_and_posthoc_role():
    require_e4e_output()
    result = pd.read_csv(RESULT)
    assert result.analysis_role.eq("post hoc").all()
    assert set(result.check) == {
        "mapped_class_support",
        "sigma_doc_full_data",
        "cli_layout_coordinate",
        "cli_layout_recall",
    }
    assert result.loc[result.check.eq("sigma_doc_full_data"), "n_items"].size == 52
    assert set(
        result.loc[result.check.eq("mapped_class_support"), "scope"]
    ) == {
        "chembl204 (Thrombin)",
        "chembl205 (CA II)",
        "chembl240 (hERG)",
    }


def test_repository_source_copy_is_complete_and_portable():
    target = WORKSPACE / "chemgridmap/paper/strengthening/source"
    if not target.is_dir():
        pytest.skip("Run from the research workspace that contains chemgridmap/paper/strengthening.")
    assert REQUIRED_SOURCE <= {path.name for path in target.glob("*.py")}
    for name in REQUIRED_SOURCE:
        text = (target / name).read_text(encoding="utf-8")
        executable_text = "\n".join(
            line for line in text.splitlines() if "assert " not in line
        )
        assert "/Users/" not in executable_text
        assert "qiuzhenming" not in executable_text


@pytest.mark.skipif(
    IS_REPOSITORY_COPY,
    reason="Analysis outputs are delivered separately from the source repository copy.",
)
def test_e4e_outputs_stay_in_allocated_directories():
    require_e4e_output()
    assert RESULT.parent == OUT
    assert RESULT.is_file()
    assert Path(__file__).resolve().parent == SOURCE
