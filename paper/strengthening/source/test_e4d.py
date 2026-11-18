from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


SOURCE = Path(__file__).with_name("e4d_analysis.py")
spec = importlib.util.spec_from_file_location("e4d_analysis", SOURCE)
assert spec is not None and spec.loader is not None
e4d = importlib.util.module_from_spec(spec)
spec.loader.exec_module(e4d)


def test_t_interval_for_paired_seed_differences():
    differences = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], dtype=float)
    result = e4d.t_interval(differences)
    half_width = result["t_critical"] * differences.std(ddof=1) / np.sqrt(10)
    assert result["paired_difference_mean"] == 5.5
    assert np.isclose(result["t_ci_low"], 5.5 - half_width)
    assert np.isclose(result["t_ci_high"], 5.5 + half_width)
    assert result["n_seeds"] == 10
    assert result["t_degrees_of_freedom"] == 9


def test_summary_uses_left_minus_right_and_preserves_ten_seeds():
    raw = pd.DataFrame(
        [
            {"seed": seed, "layout": layout, "mean_e_grid": value, "mean_recall": value / 10}
            for seed, value in zip(range(1, 11), range(1, 11))
            for layout in (
                "las_random",
                "las_projection",
                "flas_random",
                "neighbor_refined",
                "minimum_displacement",
            )
        ]
    )
    rows = e4d._summary_rows("target", "Target", raw)
    e_grid = next(
        row
        for row in rows
        if row["comparison"] == "las_random_minus_neighbor_refined"
        and row["metric"] == "mean_e_grid"
    )
    assert e_grid["n_seeds"] == 10
    assert np.isclose(e_grid["paired_difference_mean"], 0.0)
    assert e_grid["t_ci_low"] <= e_grid["paired_difference_mean"] <= e_grid["t_ci_high"]


def test_generated_output_has_complete_post_hoc_table():
    output = e4d.OUTPUT
    if not output.is_file():
        pytest.skip(f"Run {SOURCE.name} first.")
    table = pd.read_csv(output, float_precision="round_trip")
    assert set(table.analysis_role) == {"post hoc"}
    assert len(table.loc[table.row_type == "seed_layout"]) == 150
    assert len(table.loc[table.row_type == "paired_t_interval"]) == 24
    raw = table.loc[table.row_type == "seed_layout"]
    assert raw.groupby(["target", "layout"]).seed.nunique().eq(10).all()
    assert set(raw.layout) == {
        "las_random",
        "las_projection",
        "flas_random",
        "neighbor_refined",
        "minimum_displacement",
    }
    summaries = table.loc[table.row_type == "paired_t_interval"]
    assert summaries.n_seeds.eq(10).all()
    assert summaries.t_ci_low.le(summaries.paired_difference_mean).all()
    assert summaries.paired_difference_mean.le(summaries.t_ci_high).all()
