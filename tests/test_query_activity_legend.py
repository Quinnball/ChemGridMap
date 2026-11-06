"""Observed-activity query legends describe all classes and the actual cutoffs."""

from pathlib import Path
import xml.etree.ElementTree as ET

import pandas as pd
import pytest

from chemgridmap.pipeline import build_grid_map
from chemgridmap.plotting import DEFAULT_COLORS, render_grid_map


@pytest.mark.parametrize("lower,upper", [(6.0, 7.0), (5.4, 8.2)])
def test_query_activity_legend_uses_actual_thresholds_without_evidence(tmp_path, lower, upper):
    data = pd.DataFrame({"smiles": ["CCC", "CCO", "CCN"],
                         "activity": [lower - 1, (lower + upper) / 2, upper + 1]})
    result = build_grid_map(
        data, value_col="activity", projection="pca", grouping="none",
        grid_refinement="none", lower_threshold=lower, upper_threshold=upper,
        query_data=pd.DataFrame({"smiles": ["CCF"]}),
        output_dir=tmp_path, output_formats=["svg"],
    )
    assert "activity_confidence_fit" not in result.data.attrs
    root = ET.parse(result.output_files["svg"]).getroot()
    legend = next(node for node in root.iter() if node.get("data-color-legend") == "activity")
    labels = [node.text for node in legend if node.tag.endswith("text") and node.text]
    assert labels == [f"Inactive ≤ {lower:g}", f"Medium ({lower:g}, {upper:g})",
                      f"Active ≥ {upper:g}", "Not measured"]
    swatches = [node for node in legend if node.tag.endswith("rect")]
    assert [node.get("fill") for node in swatches] == [
        DEFAULT_COLORS["inactive"], DEFAULT_COLORS["medium"], DEFAULT_COLORS["active"], "#FFFFFF",
    ]
    assert swatches[-1].get("stroke") == "#8A99A6"


def test_measured_only_activity_render_ignores_query_threshold_metadata(tmp_path):
    data = pd.DataFrame({"canonical_smiles": ["CCC", "CCO"],
                         "activity_class": ["inactive", "active"],
                         "grid_row": [0, 0], "grid_col": [0, 1]})
    before = render_grid_map(data, tmp_path / "before", formats=["svg", "png"], detail="overview")
    data.attrs["activity_thresholds"] = {"lower": 5.4, "upper": 8.2}
    after = render_grid_map(data, tmp_path / "after", formats=["svg", "png"], detail="overview")
    for extension in ("svg", "png"):
        assert Path(before[extension]).read_bytes() == Path(after[extension]).read_bytes()
    assert 'data-color-legend="activity"' not in Path(after["svg"]).read_text(encoding="utf-8")
