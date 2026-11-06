"""Actual SVG and PNG rendering of measured references and unmeasured molecules."""

import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import pytest
from PIL import Image

from chemgridmap.plotting import (
    DEFAULT_COLORS,
    NOT_MEASURED_COLOR,
    NOT_MEASURED_EDGE_COLOR,
    RELIABILITY_BINS,
    render_grid_map,
)


@pytest.fixture
def query_grid():
    return pd.DataFrame({
        "molecule_identity_key": ["R1", "Q1", "R2"],
        "canonical_smiles": ["CCN", "CCO", "CCC"],
        "measurement_status": ["measured", "not_measured", "measured"],
        "activity_class": ["active", "not_measured", "inactive"],
        "class_status": ["supported", np.nan, "supported"],
        "class_support": [.9, np.nan, .95],
        "offset_pair_touched": [False, np.nan, False],
        "is_conflicted": [False, np.nan, False],
        "neighbor_recall": [.1, .5, .9],
        "structure_group": ["G001", "G001", "G002"],
        "structure_group_color": ["#2B6CB0", "#2B6CB0", "#D97706"],
        "grid_col": [0, 1, 2],
        "grid_row": [0, 0, 0],
    })


def _hex_pixels(value):
    return np.array([int(value[index:index + 2], 16) for index in (1, 3, 5)])


@pytest.mark.parametrize("color_by", ["activity", "class_support"])
@pytest.mark.parametrize("detail", ["overview", "full"])
def test_unmeasured_activity_tiles_are_white_with_gray_outline(query_grid, tmp_path, color_by, detail):
    files = render_grid_map(
        query_grid, tmp_path / f"{color_by}_{detail}", color_by=color_by,
        detail=detail, formats=["svg", "png"],
    )
    root = ET.parse(files["svg"]).getroot()
    query = next(node for node in root.iter() if node.get("data-molecule-id") == "Q1")
    tile = next(node for node in query if node.tag.endswith("rect"))
    assert tile.get("fill") == NOT_MEASURED_COLOR
    assert tile.get("stroke") == NOT_MEASURED_EDGE_COLOR
    assert tile.get("stroke-width") == "1"
    assert not any(node.get("data-offset-pair-marker") or node.get("data-conflict-marker")
                   for node in query.iter())
    legend = next(node for node in root.iter() if node.get("data-color-legend") == color_by)
    assert any(node.text == "Not measured" for node in legend)
    assert any(node.get("fill") == NOT_MEASURED_COLOR
               and node.get("stroke") == NOT_MEASURED_EDGE_COLOR for node in legend)

    image = np.asarray(Image.open(files["png"]).convert("RGB"))
    # Only inspect the map above the legend. The query lies between these two
    # measured tiles, so background and legend whites cannot satisfy this check.
    upper = image[:int(image.shape[0] * .7)]
    active = np.argwhere(np.all(upper == _hex_pixels(DEFAULT_COLORS["active"]), axis=2))
    inactive = np.argwhere(np.all(upper == _hex_pixels(DEFAULT_COLORS["inactive"]), axis=2))
    assert len(active) and len(inactive)
    left, right = int(active[:, 1].max()), int(inactive[:, 1].min())
    top, bottom = int(active[:, 0].min()), int(active[:, 0].max())
    query_interior = upper[top + 3:bottom - 3, left + 3:right - 3]
    assert np.mean(np.all(query_interior == 255, axis=2)) > .75
    query_with_border = upper[max(0, top - 3):bottom + 3, left - 3:right + 3]
    assert np.any(np.all(query_with_border == _hex_pixels(NOT_MEASURED_EDGE_COLOR), axis=2))


@pytest.mark.parametrize("color_by,expected", [
    ("neighbor_recall", RELIABILITY_BINS[2][1]),
    ("structure_group", "#2B6CB0"),
])
def test_unmeasured_structural_layers_use_existing_colors(query_grid, tmp_path, color_by, expected):
    files = render_grid_map(
        query_grid, tmp_path / color_by, color_by=color_by, detail="overview", formats=["svg", "png"],
    )
    root = ET.parse(files["svg"]).getroot()
    query = next(node for node in root.iter() if node.get("data-molecule-id") == "Q1")
    tile = next(node for node in query if node.tag.endswith("rect"))
    assert tile.get("fill") == expected
    assert tile.get("stroke") == "#FFFFFF"
    assert "Not measured" not in "".join(root.itertext())
    image = np.asarray(Image.open(files["png"]).convert("RGB"))
    assert np.any(np.all(image[:int(image.shape[0] * .7)] == _hex_pixels(expected), axis=2))


@pytest.mark.parametrize("missing", [None, ""])
def test_ungrouped_query_keeps_existing_no_group_color(query_grid, tmp_path, missing):
    query_grid.loc[1, ["structure_group", "structure_group_color"]] = missing
    files = render_grid_map(query_grid, tmp_path / "ungrouped", color_by="structure_group", formats=["svg", "png"])
    root = ET.parse(files["svg"]).getroot()
    query = next(node for node in root.iter() if node.get("data-molecule-id") == "Q1")
    tile = next(node for node in query if node.tag.endswith("rect"))
    assert tile.get("fill") == DEFAULT_COLORS["unlabelled"]
    assert "No group" in "".join(root.itertext())
    image = np.asarray(Image.open(files["png"]).convert("RGB"))
    assert np.any(np.all(image[:int(image.shape[0] * .7)] == _hex_pixels(DEFAULT_COLORS["unlabelled"]), axis=2))
