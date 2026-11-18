"""Portability checks for the complete analysis and figure source set."""
from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
STRENGTHENING = ROOT / "paper/strengthening"
FIGURES = ROOT / "paper/figures"

ANALYSIS_FILES = {
    "e5_analysis.py",
    "e5_models.py",
    "e5b_analysis.py",
    "e5b_models.py",
    "e5c_analysis.py",
    "e5c_models.py",
    "e5d_analysis.py",
    "e5d_models.py",
    "e6_analysis.py",
    "e6_core.py",
    "e6_prepare.py",
    "e6_regression_diag.py",
    "e7_analysis.py",
    "e9_analysis.py",
    "e10_analysis.py",
    "e11_analysis.py",
    "download_chembl37_validation.py",
    "e12a_rule_calibration.py",
    "e12b_analysis.py",
    "t1_scale.py",
    "t1_ui_capture.cjs",
}

FIGURE_FILES = {
    "build_figures.py",
    "common.py",
    "data_final.py",
    "fig2.py",
    "fig3.py",
    "fig4.py",
    "fig5.py",
    "fig6.py",
    "fig6_capture.cjs",
    "molecule.py",
    "qa_figures.py",
    "style.py",
}

TEST_FILES = {
    "test_e5.py",
    "test_e5b.py",
    "test_e5c.py",
    "test_e5d.py",
    "test_e6.py",
    "test_e7.py",
    "test_e9.py",
    "test_e10.py",
    "test_e11.py",
    "test_e12a.py",
    "test_e12b.py",
    "test_t1.py",
}


def test_required_portable_sources_exist():
    source = STRENGTHENING / "source"
    figure_source = FIGURES / "source"
    assert ANALYSIS_FILES <= {path.name for path in source.iterdir()}
    assert TEST_FILES <= {path.name for path in source.iterdir()}
    assert FIGURE_FILES <= {path.name for path in figure_source.iterdir()}
    assert (FIGURES / "data/target_symbols.csv").is_file()


def test_portable_sources_do_not_embed_workspace_paths():
    forbidden = (
        "/" + "Users/",
        "qiu" + "zhenming",
        "codex" + "_tasks",
        "pro" + "ject1",
        "JCIM_" + "reliability_revision",
    )
    for base in (STRENGTHENING / "source", FIGURES / "source"):
        for path in base.iterdir():
            if not path.is_file() or path.suffix not in {".py", ".cjs"}:
                continue
            text = path.read_text(encoding="utf-8")
            executable = "\n".join(
                line for line in text.splitlines() if "assert " not in line
            )
            for token in forbidden:
                assert token not in executable, (path.name, token)


def test_sealed_validation_plan_matches_recorded_hash():
    plan = STRENGTHENING / "plans/E11_validation_plan_sealed.md"
    digest = hashlib.sha256(plan.read_bytes()).hexdigest()
    recorded = (
        STRENGTHENING / "plans/E11_validation_plan_sealed.sha256"
    ).read_text(encoding="utf-8").split()[0]
    assert digest == recorded
