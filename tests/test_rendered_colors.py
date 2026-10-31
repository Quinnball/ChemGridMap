"""Color-vision checks for the rendered reliability layers."""

from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
from PIL import Image

from chemgridmap.plotting import RELIABILITY_BINS, render_grid_map


DEUTERANOPIA_1_0 = np.array([
    [0.367322, 0.860646, -0.227968],
    [0.280085, 0.672501, 0.047413],
    [-0.011820, 0.042940, 0.968881],
])


def _linear(rgb):
    rgb = np.asarray(rgb, dtype=float)
    return np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)


def _gamma(rgb):
    rgb = np.clip(rgb, 0, 1)
    return np.where(rgb <= 0.0031308, 12.92 * rgb, 1.055 * rgb ** (1 / 2.4) - 0.055)


def _hex_rgb(value):
    return np.array([int(value[index:index + 2], 16) / 255 for index in (1, 3, 5)])


def srgb_to_oklab(color):
    """Convert an sRGB hex string or RGB triplet to OKLab."""
    rgb = _hex_rgb(color) if isinstance(color, str) else np.asarray(color, dtype=float)
    red, green, blue = _linear(rgb)
    l = 0.4122214708 * red + 0.5363325363 * green + 0.0514459929 * blue
    m = 0.2119034982 * red + 0.6806995451 * green + 0.1073969566 * blue
    s = 0.0883024619 * red + 0.2817188376 * green + 0.6299787005 * blue
    l, m, s = np.cbrt([l, m, s])
    return np.array([
        0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s,
        1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s,
        0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s,
    ])


def delta_e_100(first, second):
    """OKLab Euclidean difference multiplied by 100."""
    return float(np.linalg.norm(srgb_to_oklab(first) - srgb_to_oklab(second)) * 100)


def simulate_deuteranopia(color):
    """Machado et al. (2009), severity 1.0 deuteranopia simulation."""
    linear = DEUTERANOPIA_1_0 @ _linear(_hex_rgb(color) if isinstance(color, str) else color)
    return _gamma(linear)


def sampled_png_color(path, expected=None):
    """Return the dominant rendered color above the legend, sampled from its center."""
    image = np.asarray(Image.open(path).convert("RGB"))
    upper = image[: int(image.shape[0] * .9)]
    if expected is not None:
        expected_rgb = (_hex_rgb(expected) * 255).round().astype(int)
        mask = np.max(np.abs(upper.astype(int) - expected_rgb), axis=2) <= 2
        points = np.argwhere(mask)
        assert len(points), f"{expected} was not rendered in {path}"
        # The topmost exact match lies in the cell's colored frame, not the legend.
        center = points[np.argmin(points[:, 0])]
        return tuple(int(value) for value in upper[center[0], center[1]])
    colors, counts = np.unique(upper.reshape(-1, 3), axis=0, return_counts=True)
    white_distance = np.linalg.norm(colors.astype(float) - 255, axis=1)
    chromatic = np.ptp(colors, axis=1) > 15
    candidates = colors[(white_distance > 18) & chromatic]
    candidate_counts = counts[(white_distance > 18) & chromatic]
    chosen = candidates[int(np.argmax(candidate_counts))]
    mask = np.all(upper == chosen, axis=2)
    points = np.argwhere(mask)
    center = np.median(points, axis=0).astype(int)
    return tuple(int(value) for value in upper[center[0], center[1]])


def test_rendered_class_support_and_recall_colors(tmp_path):
    data = pd.DataFrame({
        "canonical_smiles": ["CCO", "CCC", "CCN"],
        "activity_class": ["active", "inactive", "active"],
        "class_status": ["supported", "supported", "uncertain"],
        "neighbor_recall": [.1, .5, .9],
        "grid_col": [0, 1, 2],
        "grid_row": [0, 0, 0],
    })
    data.attrs["activity_confidence_fit"] = {"confidence_threshold": .8, "offset_display": "uncertain"}
    overview = render_grid_map(
        data, tmp_path / "support_overview", color_by="class_support",
        detail="overview", formats=["svg", "png"],
    )
    full = render_grid_map(
        data, tmp_path / "support_full", color_by="class_support",
        detail="full", formats=["svg", "png"],
    )
    active = "#%02x%02x%02x" % sampled_png_color(overview["png"])
    assert delta_e_100(simulate_deuteranopia("#2A9D8F"), simulate_deuteranopia("#C44E52")) >= 8
    assert delta_e_100("#2A9D8F", "#b8e0db") >= 15
    assert delta_e_100("#C44E52", "#ebc0c1") >= 15
    assert active == "#2a9d8f"
    assert sampled_png_color(full["png"], "#2A9D8F") == tuple(round(value * 255) for value in _hex_rgb("#2A9D8F"))

    recalls = [(bin_value, color) for bin_value, color in RELIABILITY_BINS] + [(1.0, "#ffffff")]
    for (_, first), (_, second) in zip(recalls[:4], recalls[1:5]):
        assert delta_e_100(first, second) >= 10
    assert delta_e_100(RELIABILITY_BINS[0][1], "#ffffff") >= 10

    root = ET.parse(overview["svg"]).getroot()
    cells = [node for node in root.iter() if "data-grid-col" in node.attrib]
    assert all(
        next(child for child in cell if child.tag.endswith("rect")).get("fill-opacity") == "1"
        for cell in cells
    )
    legend = next(node for node in root.iter() if node.get("data-color-legend") == "class_support")
    swatches = [node.get("fill") for node in legend.iter() if node.tag.endswith("rect")]
    cell_fills = [
        next(child for child in cell if child.tag.endswith("rect")).get("fill")
        for cell in cells
    ]
    assert set(cell_fills) <= set(swatches)
