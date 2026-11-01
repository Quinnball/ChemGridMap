"""Check size, text, geometry, data, and vectors for Figures 2-5."""

from __future__ import annotations

import itertools
import json
import math
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

import matplotlib.text as mtext
import numpy as np
import pandas as pd
import pymupdf
from PIL import Image
from matplotlib.collections import Collection, LineCollection, PathCollection, PolyCollection
from matplotlib.colors import to_hex, to_rgb
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from chemgridmap.plotting import DEFAULT_COLORS

import common as c


SOURCE_HEIGHT_MM = {2: 125.0, 3: 160.0, 4: 130.0, 5: 150.0, 6: 220.0}
MIN_FONT_PT = 6.0
MIN_LINE_PT = 0.4


def _bbox_intersection(a, b, tolerance: float = 0.0) -> bool:
    width = min(a.x1, b.x1) - max(a.x0, b.x0)
    height = min(a.y1, b.y1) - max(a.y0, b.y0)
    return width > tolerance and height > tolerance


def _audio_texts(fig: Figure) -> list[mtext.Text]:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    result = []
    for text in fig.findobj(mtext.Text):
        if not text.get_visible() or not text.get_text().strip():
            continue
        bbox = text.get_window_extent(renderer)
        if bbox.width > 0 and bbox.height > 0:
            result.append(text)
    return result


def _text_overlap_matplotlib(fig: Figure) -> list[list[object]]:
    renderer = fig.canvas.get_renderer()
    texts = _audio_texts(fig)
    failures = []
    for first, second in itertools.combinations(texts, 2):
        a = first.get_window_extent(renderer)
        b = second.get_window_extent(renderer)
        if _bbox_intersection(a, b, tolerance=0.05):
            failures.append(
                [
                    first.get_text(),
                    second.get_text(),
                    [round(a.x0, 2), round(a.y0, 2), round(a.x1, 2), round(a.y1, 2)],
                    [round(b.x0, 2), round(b.y0, 2), round(b.x1, 2), round(b.y1, 2)],
                ]
            )
    return failures


def _bbox_gap_mm(a, b, dpi: float) -> float:
    dx = max(a.x0 - b.x1, b.x0 - a.x1, 0.0)
    dy = max(a.y0 - b.y1, b.y0 - a.y1, 0.0)
    if dx == 0.0 and dy == 0.0:
        overlap = min(
            min(a.x1, b.x1) - max(a.x0, b.x0),
            min(a.y1, b.y1) - max(a.y0, b.y0),
        )
        return -overlap * 25.4 / dpi
    return float(np.hypot(dx, dy) * 25.4 / dpi)


def _text_spacing_check(fig: Figure) -> dict[str, object]:
    """Text bounding boxes must be at least 0.5 mm apart."""
    renderer = fig.canvas.get_renderer()
    texts = _audio_texts(fig)
    minimum = None
    failures = []
    for first, second in itertools.combinations(texts, 2):
        first_box = first.get_window_extent(renderer)
        second_box = second.get_window_extent(renderer)
        gap_mm = _bbox_gap_mm(first_box, second_box, fig.dpi)
        minimum = gap_mm if minimum is None else min(minimum, gap_mm)
        if gap_mm < 0.5 - 1e-9:
            failures.append(
                {
                    "first": first.get_text(),
                    "second": second.get_text(),
                    "gap_mm": round(gap_mm, 6),
                    "first_bbox_mm": [
                        round(value, 4)
                        for value in (
                            first_box.x0 * 25.4 / fig.dpi,
                            first_box.y0 * 25.4 / fig.dpi,
                            first_box.x1 * 25.4 / fig.dpi,
                            first_box.y1 * 25.4 / fig.dpi,
                        )
                    ],
                    "second_bbox_mm": [
                        round(value, 4)
                        for value in (
                            second_box.x0 * 25.4 / fig.dpi,
                            second_box.y0 * 25.4 / fig.dpi,
                            second_box.x1 * 25.4 / fig.dpi,
                            second_box.y1 * 25.4 / fig.dpi,
                        )
                    ],
                }
            )
    return {
        "texts_checked": len(texts),
        "minimum_required_mm": 0.5,
        "minimum_gap_mm": round(minimum, 6) if minimum is not None else None,
        "failures": failures,
        "passed": not failures,
    }


def _legend_texts(fig: Figure) -> list[mtext.Text]:
    texts = []
    legends = list(fig.legends)
    legends.extend(ax.get_legend() for ax in fig.axes if ax.get_legend() is not None)
    for legend in legends:
        texts.extend(legend.get_texts())
    texts.extend(
        text
        for text in fig.findobj(mtext.Text)
        if getattr(text, "_qa_legend_text", False)
    )
    return list({id(text): text for text in texts}.values())


def _legend_title_spacing_check(fig: Figure) -> dict[str, object]:
    """Legend text and panel titles must be at least 1.5 mm apart."""
    renderer = fig.canvas.get_renderer()
    legend_texts = _legend_texts(fig)
    title_texts = [
        text
        for text in fig.findobj(mtext.Text)
        if getattr(text, "_qa_panel_title", False)
    ]
    failures = []
    minimum = None
    for legend_text in legend_texts:
        legend_box = legend_text.get_window_extent(renderer)
        for title_text in title_texts:
            title_box = title_text.get_window_extent(renderer)
            gap_mm = _bbox_gap_mm(legend_box, title_box, fig.dpi)
            minimum = gap_mm if minimum is None else min(minimum, gap_mm)
            if gap_mm < 1.5 - 1e-9:
                failures.append(
                    {
                        "legend": legend_text.get_text(),
                        "title": title_text.get_text(),
                        "gap_mm": round(gap_mm, 6),
                        "legend_bbox_mm": [
                            round(value, 4)
                            for value in (
                                legend_box.x0 * 25.4 / fig.dpi,
                                legend_box.y0 * 25.4 / fig.dpi,
                                legend_box.x1 * 25.4 / fig.dpi,
                                legend_box.y1 * 25.4 / fig.dpi,
                            )
                        ],
                        "title_bbox_mm": [
                            round(value, 4)
                            for value in (
                                title_box.x0 * 25.4 / fig.dpi,
                                title_box.y0 * 25.4 / fig.dpi,
                                title_box.x1 * 25.4 / fig.dpi,
                                title_box.y1 * 25.4 / fig.dpi,
                            )
                        ],
                    }
                )
    return {
        "legend_texts_checked": len(legend_texts),
        "panel_titles_checked": len(title_texts),
        "minimum_required_mm": 1.5,
        "minimum_gap_mm": round(minimum, 6) if minimum is not None else None,
        "failures": failures,
        "passed": not failures,
    }


def _fig2c_spacing_check(fig: Figure) -> dict[str, object]:
    """The Figure 2c mean label must stay at least 1.5 mm below the colour row."""
    renderer = fig.canvas.get_renderer()
    swatches = [
        patch
        for ax in fig.axes
        for patch in ax.patches
        if getattr(patch, "_qa_fig2c_swatch", False)
    ]
    labels = [
        text
        for text in fig.findobj(mtext.Text)
        if text.get_text().startswith(("Mean 0.", "Fixed k 0."))
    ]
    failures = []
    measurements = []
    for label in labels:
        label_box = label.get_window_extent(renderer)
        nearest = min(
            swatches,
            key=lambda patch: _bbox_gap_mm(
                label_box,
                patch.get_window_extent(renderer),
                fig.dpi,
            ),
        )
        gap_mm = _bbox_gap_mm(
            label_box,
            nearest.get_window_extent(renderer),
            fig.dpi,
        )
        item = {"text": label.get_text(), "gap_mm": round(gap_mm, 6)}
        measurements.append(item)
        if gap_mm < 1.5 - 1e-9:
            failures.append(item)
    return {
        "labels_checked": len(labels),
        "swatches_checked": len(swatches),
        "minimum_required_mm": 1.5,
        "measurements": measurements,
        "failures": failures,
        "passed": len(labels) == 2 and len(swatches) == 5 and not failures,
    }


def _svg_grid_fills(path: Path) -> dict[tuple[int, int], str]:
    root = ET.parse(path).getroot()
    fills = {}
    for element in root.iter():
        if "data-grid-row" not in element.attrib:
            continue
        rect = next(
            (
                child
                for child in element
                if child.tag.rsplit("}", 1)[-1] == "rect"
            ),
            None,
        )
        if rect is None:
            continue
        key = (int(element.attrib["data-grid-row"]), int(element.attrib["data-grid-col"]))
        fills[key] = str(rect.attrib["fill"]).lower()
    return fills


def _fig2_color_check(fig: Figure) -> dict[str, object]:
    """Check panels a-c cell by cell against the software SVG, legend, and colours."""
    svg_dir = (
        c.FIG_DIR / "data_fig2" / "default"
    )
    svg_files = {
        "a": svg_dir / "caii_seed42_activity.svg",
        "b": svg_dir / "caii_seed42_class_support.svg",
        "c": svg_dir / "caii_seed42_neighbor_recall.svg",
    }
    tiles: dict[str, dict[tuple[int, int], object]] = {
        "a": {},
        "b": {},
        "c": {},
    }
    swatches = []
    triangles = []
    legend_triangles = []
    for ax in fig.axes:
        for patch in ax.patches:
            if hasattr(patch, "_qa_fig2_grid_panel"):
                tiles[patch._qa_fig2_grid_panel][
                    (patch._qa_fig2_grid_row, patch._qa_fig2_grid_col)
                ] = patch
            if hasattr(patch, "_qa_fig2_color_panel"):
                swatches.append(patch)
            if hasattr(patch, "_qa_fig2_triangle_leg"):
                triangles.append(patch)
            if getattr(patch, "_qa_fig2_legend_triangle", False):
                legend_triangles.append(patch)

    failures = []
    color_counts = {}
    for panel, path in svg_files.items():
        expected = _svg_grid_fills(path)
        actual = tiles[panel]
        if set(expected) != set(actual):
            failures.append(
                {
                    "panel": panel,
                    "problem": "grid keys differ",
                    "missing": len(set(expected) - set(actual)),
                    "extra": len(set(actual) - set(expected)),
                }
            )
            continue
        mismatches = [
            {
                "row": key[0],
                "col": key[1],
                "svg": expected[key],
                "figure": to_hex(actual[key].get_facecolor(), keep_alpha=False).lower(),
                "expected": actual[key]._qa_fig2_expected_fill,
            }
            for key in expected
            if expected[key]
            != to_hex(actual[key].get_facecolor(), keep_alpha=False).lower()
            or expected[key] != actual[key]._qa_fig2_expected_fill
        ]
        color_counts[panel] = dict(Counter(expected.values()))
        failures.extend(
            [{"panel": panel, **mismatch} for mismatch in mismatches]
        )

    swatch_failures = [
        {
            "panel": patch._qa_fig2_color_panel,
            "expected": patch._qa_fig2_expected_fill,
            "actual": to_hex(patch.get_facecolor(), keep_alpha=False).lower(),
        }
        for patch in swatches
        if to_hex(patch.get_facecolor(), keep_alpha=False).lower()
        != patch._qa_fig2_expected_fill
    ]
    triangle_panels = Counter(patch._qa_fig2_triangle_panel for patch in triangles)
    triangle_legs = sorted(patch._qa_fig2_triangle_leg for patch in triangles)

    def upper_right_angle(patch) -> bool:
        vertices = np.asarray(patch._qa_fig2_triangle_vertices)  # type: ignore[attr-defined]
        right = float(vertices[:, 0].max())
        top = float(vertices[:, 1].max())
        return any(
            math.isclose(x, right, abs_tol=1e-9)
            and math.isclose(y, top, abs_tol=1e-9)
            for x, y in vertices
        )

    triangle_styles_ok = all(
        math.isclose(
            patch._qa_fig2_triangle_edgewidth,  # type: ignore[attr-defined]
            0.3,
            abs_tol=1e-9,
        )
        and to_hex(patch.get_edgecolor(), keep_alpha=False).lower() == "#ffffff"
        and upper_right_angle(patch)
        for patch in [*triangles, *legend_triangles]
    )
    triangle_cells_ok = all(
        all(
            abs(float(x) - col) <= 0.5 + 1e-9
            and abs(float(y) - row) <= 0.5 + 1e-9
            for x, y in patch._qa_fig2_triangle_vertices  # type: ignore[attr-defined]
        )
        for patch in triangles
        for col, row in [patch._qa_fig2_triangle_cell]  # type: ignore[attr-defined]
    )
    triangles_ok = (
        triangle_panels == {"a": 17, "b": 17}
        and triangle_legs
        and all(math.isclose(value, 0.62) for value in triangle_legs)
        and len(legend_triangles) == 1
        and triangle_styles_ok
        and triangle_cells_ok
    )
    return {
        "cells_compared": sum(len(value) for value in tiles.values()),
        "expected_cells": 889 * 3,
        "color_counts": color_counts,
        "grid_mismatches": failures[:80],
        "legend_swatches_checked": len(swatches),
        "legend_failures": swatch_failures,
        "triangles": {
            "counts": dict(triangle_panels),
            "legs_grid_cells": sorted(set(triangle_legs)),
            "legend_count": len(legend_triangles),
            "same_upper_right_corner": triangle_styles_ok,
            "inside_cells": triangle_cells_ok,
        },
        "passed": (
            not failures
            and not swatch_failures
            and len(swatches) == 13
            and triangles_ok
        ),
    }


def _fig4a_check(fig: Figure) -> dict[str, object]:
    """Compare Figure 4a curves point by point with the class calibration table."""
    panels = [
        ax for ax in fig.axes if getattr(ax, "_qa_fig4a", False)
    ]
    n_labels = [
        text.get_text()
        for ax in panels
        for text in ax.texts
        if text.get_text().startswith(("n=", "n\n"))
    ]
    calibration = c.load_csv("4a_calibration_by_class.csv")
    expected_styles = {
        "inactive": (to_hex(DEFAULT_COLORS["inactive"], keep_alpha=False).lower(), "o"),
        "medium": (to_hex(DEFAULT_COLORS["medium"], keep_alpha=False).lower(), "s"),
        "active": (to_hex(DEFAULT_COLORS["active"], keep_alpha=False).lower(), "^"),
    }
    comparison_failures = []
    actual_counts = {}
    for panel in panels:
        endpoint = panel._qa_fig4a_endpoint  # type: ignore[attr-defined]
        actual = sorted(
            (
                item["display_class"],
                item["probability_bin"],
                item["mean_p"],
                item["agreement_rate"],
                item["wilson_ci_low"],
                item["wilson_ci_high"],
                item["color"],
                item["marker"],
            )
            for item in panel._qa_fig4a_records  # type: ignore[attr-defined]
        )
        visible = calibration.loc[
            calibration["endpoint"].eq(endpoint)
            & ~calibration["omit_in_figure"].astype(bool)
        ]
        expected = sorted(
            (
                str(row.display_class),
                int(row.probability_bin),
                float(row.mean_p),
                float(row.agreement_rate),
                float(row.wilson_ci_low),
                float(row.wilson_ci_high),
                expected_styles[str(row.display_class)][0],
                expected_styles[str(row.display_class)][1],
            )
            for row in visible.itertuples(index=False)
        )
        actual_counts[endpoint] = len(actual)
        if len(actual) != len(expected):
            comparison_failures.append(
                {
                    "endpoint": endpoint,
                    "problem": "point count",
                    "actual": len(actual),
                    "expected": len(expected),
                }
            )
            continue
        for actual_row, expected_row in zip(actual, expected, strict=True):
            if (
                actual_row[:2] != expected_row[:2]
                or actual_row[6:] != expected_row[6:]
                or not np.allclose(
                    actual_row[2:6],
                    expected_row[2:6],
                    rtol=0.0,
                    atol=1e-15,
                )
            ):
                comparison_failures.append(
                    {
                        "endpoint": endpoint,
                        "actual": actual_row,
                        "expected": expected_row,
                    }
                )
    limits = {
        panel.get_label(): {
            "xlim": [float(value) for value in panel.get_xlim()],
            "ylim": [float(value) for value in panel.get_ylim()],
        }
        for panel in panels
    }
    limits_ok = all(
        np.allclose(panel.get_xlim(), (-0.03, 1.03))
        and np.allclose(panel.get_ylim(), (-0.03, 1.03))
        for panel in panels
    )
    legend_labels = {
        panel._qa_fig4a_endpoint: (  # type: ignore[attr-defined]
            [text.get_text() for text in panel.get_legend().get_texts()]
            if panel.get_legend() is not None
            else []
        )
        for panel in panels
    }
    return {
        "panels_checked": len(panels),
        "n_labels": n_labels,
        "limits": limits,
        "actual_point_counts": actual_counts,
        "comparison_failures": comparison_failures[:40],
        "legend_labels": legend_labels,
        "passed": (
            len(panels) == 2
            and not n_labels
            and limits_ok
            and actual_counts == {"IC50": 12, "Ki": 11}
            and not comparison_failures
            and legend_labels == {
                "IC50": ["Inactive", "Medium", "Active"],
                "Ki": [],
            }
        ),
    }


def _fig4d_check(fig: Figure, png_path: Path) -> dict[str, object]:
    panels = [ax for ax in fig.axes if getattr(ax, "_qa_fig4d", False)]
    bars = [
        patch
        for ax in panels
        for patch in ax.patches
        if hasattr(patch, "_qa_fig4d_state")
    ]
    omitted = []
    expected_bars = 0
    for panel in panels:
        for row in panel._qa_fig4d_rows:  # type: ignore[attr-defined]
            if int(row["n_predictions"]) < 8:
                omitted.append(
                    {
                        "endpoint": panel._qa_fig4d_endpoint,  # type: ignore[attr-defined]
                        "group": row["earlier_units_group"],
                        "state": row["state"],
                        "n_predictions": int(row["n_predictions"]),
                    }
                )
            else:
                expected_bars += 1
    counts = Counter(patch._qa_fig4d_state for patch in bars)
    expected_counts = {
        "P_supported": 4,
        "uncertain_low_support": 4,
        "touched_all": 2,
    }
    colors = {
        state: sorted({patch._qa_fig4d_color for patch in bars if patch._qa_fig4d_state == state})
        for state in expected_counts
    }
    expected_omissions = [
        {"endpoint": "IC50", "group": "m_eq_1", "state": "touched_all", "n_predictions": 2},
        {"endpoint": "Ki", "group": "m_eq_1", "state": "touched_all", "n_predictions": 3},
    ]
    intervals = [
        collection
        for panel in panels
        for collection in panel.collections
        if getattr(collection, "_qa_fig4d_interval", False)
    ]
    bar_zorder = max((patch.get_zorder() for patch in bars), default=np.inf)
    interval_zorder = min(
        (collection.get_zorder() for collection in intervals),
        default=-np.inf,
    )
    count_labels = sorted(
        [
            tuple(text._qa_fig4d_count)  # type: ignore[attr-defined]
            for panel in panels
            for text in panel.texts
            if hasattr(text, "_qa_fig4d_count")
        ]
    )
    expected_count_labels = sorted(
        [
            (
                panel._qa_fig4d_endpoint,  # type: ignore[attr-defined]
                row["earlier_units_group"],
                row["state"],
                int(row["n_predictions"]),
                int(row["n_datasets"]),
            )
            for panel in panels
            for row in panel._qa_fig4d_rows  # type: ignore[attr-defined]
            if int(row["n_predictions"]) >= 8
        ]
    )
    legend_patches = [
        handle
        for legend in fig.legends
        for handle in (
            getattr(legend, "legend_handles", None)
            or getattr(legend, "legendHandles", [])
        )
        if hasattr(handle, "_qa_fig4d_legend_state")
    ]
    bar_alpha = {
        state: sorted({float(patch.get_alpha()) for patch in bars if patch._qa_fig4d_state == state})
        for state in expected_counts
    }
    legend_alpha = {
        handle._qa_fig4d_legend_state: float(handle.get_alpha())
        for handle in legend_patches
    }
    rgb = np.asarray(Image.open(png_path).convert("RGB"))
    renderer = fig.canvas.get_renderer()
    pixel_comparisons = []
    for handle in legend_patches:
        state = handle._qa_fig4d_legend_state
        bar = next(
            patch
            for patch in bars
            if patch._qa_fig4d_state == state
        )
        legend_box = handle.get_window_extent(renderer)
        legend_x = round((legend_box.x0 + legend_box.x1) / 2)
        legend_y = round(fig.bbox.height - (legend_box.y0 + legend_box.y1) / 2)
        bar_box = bar.get_window_extent(renderer)
        bar_x = round(bar_box.x0 + 0.25 * bar_box.width)
        bar_y = round(fig.bbox.height - (bar_box.y0 + 0.25 * bar_box.height))
        legend_rgb = rgb[legend_y, legend_x].astype(int).tolist()
        bar_rgb = rgb[bar_y, bar_x].astype(int).tolist()
        channel_differences = np.abs(
            np.asarray(legend_rgb) - np.asarray(bar_rgb)
        ).tolist()
        pixel_comparisons.append(
            {
                "state": state,
                "legend_rgb": legend_rgb,
                "bar_rgb": bar_rgb,
                "channel_differences": channel_differences,
                "passed": bool(max(channel_differences) <= 2),
            }
        )
    return {
        "panels_checked": len(panels),
        "bars_expected": expected_bars,
        "bars_plotted": len(bars),
        "state_counts": dict(counts),
        "colors": colors,
        "omitted_groups": omitted,
        "intervals_checked": len(intervals),
        "bar_zorder": bar_zorder,
        "interval_zorder": interval_zorder,
        "count_labels": count_labels,
        "bar_alpha": bar_alpha,
        "legend_alpha": legend_alpha,
        "legend_pixel_comparisons": pixel_comparisons,
        "passed": (
            len(panels) == 2
            and len(bars) == expected_bars == 10
            and counts == expected_counts
            and omitted == expected_omissions
            and len(intervals) == 20
            and interval_zorder > bar_zorder
            and count_labels == expected_count_labels
            and all(len(values) == 1 for values in bar_alpha.values())
            and len(legend_alpha) == 3
            and all(
                len(bar_alpha[state]) == 1
                and np.isclose(
                    bar_alpha[state][0],
                    legend_alpha.get(state, np.nan),
                )
                for state in expected_counts
            )
            and len(pixel_comparisons) == 3
            and all(item["passed"] for item in pixel_comparisons)
        ),
    }


def _text_margin_failures(fig: Figure) -> list[list[object]]:
    renderer = fig.canvas.get_renderer()
    margin_px = 0.5 * fig.dpi / 25.4
    boundary = (0, 0, fig.bbox.width, fig.bbox.height)
    failures = []
    for text in _audio_texts(fig):
        bbox = text.get_window_extent(renderer)
        if (
            bbox.x0 < boundary[0] + margin_px
            or bbox.y0 < boundary[1] + margin_px
            or bbox.x1 > boundary[2] - margin_px
            or bbox.y1 > boundary[3] - margin_px
        ):
            failures.append(
                [
                    text.get_text(),
                    [round(bbox.x0, 2), round(bbox.y0, 2), round(bbox.x1, 2), round(bbox.y1, 2)],
                    round(margin_px, 2),
                ]
            )
    return failures


def _text_owner_axis(ax, text: mtext.Text) -> bool:
    """Test whether text belongs to an axes, excluding cross-axes figure text."""
    candidates = [
        *ax.texts,
        *ax.get_xticklabels(),
        *ax.get_yticklabels(),
        ax.title,
        ax.xaxis.label,
        ax.yaxis.label,
    ]
    if ax.get_legend() is not None:
        candidates.extend(ax.get_legend().get_texts())
    return any(text is candidate for candidate in candidates)


def _text_geometry_check(fig: Figure) -> dict[str, object]:
    """Check that text stays on the canvas, outside other axes, and does not overlap text."""
    renderer = fig.canvas.get_renderer()
    texts = _audio_texts(fig)
    owners = {}
    for text in texts:
        owner = next(
            (ax for ax in fig.axes if _text_owner_axis(ax, text)),
            None,
        )
        if owner is None:
            # Treat figure-level coordinates in one plot area as that axes annotation.
            intersecting_axes = [
                ax
                for ax in fig.axes
                if ax.get_visible()
                and _bbox_intersection(
                    text.get_window_extent(renderer),
                    ax.bbox,
                    tolerance=0.05,
                )
            ]
            if len(intersecting_axes) == 1:
                owner = intersecting_axes[0]
        owners[id(text)] = owner

    def box_mm(box) -> list[float]:
        return [
            round(float(value * 25.4 / fig.dpi), 6)
            for value in (box.x0, box.y0, box.x1, box.y1)
        ]

    canvas = fig.bbox
    canvas_failures = []
    axis_intersections = []
    for text in texts:
        box = text.get_window_extent(renderer)
        if (
            box.x0 < canvas.x0 - 0.05
            or box.y0 < canvas.y0 - 0.05
            or box.x1 > canvas.x1 + 0.05
            or box.y1 > canvas.y1 + 0.05
        ):
            canvas_failures.append(
                {"text": text.get_text(), "bbox_mm": box_mm(box)}
            )
        for ax in fig.axes:
            if not ax.get_visible() or owners[id(text)] is ax:
                continue
            if _bbox_intersection(box, ax.bbox, tolerance=0.05):
                axis_intersections.append(
                    {
                        "text": text.get_text(),
                        "text_bbox_mm": box_mm(box),
                        "other_axis": ax.get_label()
                        or f"axes@{ax.get_position().bounds}",
                        "other_axis_bbox_mm": box_mm(ax.bbox),
                    }
                )

    text_overlaps = []
    same_axis_exceptions = []
    for first, second in itertools.combinations(texts, 2):
        first_box = first.get_window_extent(renderer)
        second_box = second.get_window_extent(renderer)
        if not _bbox_intersection(first_box, second_box, tolerance=0.05):
            continue
        first_owner = owners[id(first)]
        second_owner = owners[id(second)]
        item = {
            "first": first.get_text(),
            "second": second.get_text(),
            "first_bbox_mm": box_mm(first_box),
            "second_bbox_mm": box_mm(second_box),
        }
        if first_owner is not None and first_owner is second_owner:
            item["axis"] = (
                first_owner.get_label()
                or f"axes@{first_owner.get_position().bounds}"
            )
            same_axis_exceptions.append(item)
        else:
            text_overlaps.append(item)

    return {
        "texts_checked": len(texts),
        "canvas_failures": canvas_failures,
        "cross_axis_intersections": axis_intersections,
        "text_overlaps": text_overlaps,
        "same_axis_allowed_exceptions": same_axis_exceptions,
        "passed": (
            not canvas_failures
            and not axis_intersections
            and not text_overlaps
        ),
    }


def _sample_segment(p0, p1, spacing_px: float):
    length = float(np.hypot(p1[0] - p0[0], p1[1] - p0[1]))
    count = max(1, int(np.ceil(length / spacing_px)))
    for fraction in np.linspace(0, 1, count + 1):
        yield (
            p0[0] + (p1[0] - p0[0]) * fraction,
            p0[1] + (p1[1] - p0[1]) * fraction,
        )


def _line_geometry(line, spacing_px: float):
    transform = line.get_transform()
    vertices = transform.transform(line.get_xydata())
    if line.get_linestyle() not in {"None", "none", "", " "}:
        for p0, p1 in zip(vertices[:-1], vertices[1:], strict=False):
            yield from _sample_segment(p0, p1, spacing_px)
    marker = line.get_marker()
    if marker not in {None, "None", "none", "", " "}:
        radius = line.get_markersize() / 2 * line.figure.dpi / 72
        for x, y in vertices:
            yield x, y, radius


def _collection_geometry(collection, spacing_px: float):
    if isinstance(collection, PathCollection):
        transform = collection.get_offset_transform()
        offsets = transform.transform(np.asarray(collection.get_offsets()))
        sizes = collection.get_sizes()
        for index, (x, y) in enumerate(offsets):
            area = float(sizes[index] if len(sizes) > 1 else sizes[0]) if len(sizes) else 9.0
            yield x, y, math.sqrt(area) / 2 * collection.figure.dpi / 72
        return
    transform = collection.get_transform()
    for segment in collection.get_segments():
        points = transform.transform(np.asarray(segment))
        for p0, p1 in zip(points[:-1], points[1:], strict=False):
            yield from _sample_segment(p0, p1, spacing_px)


def _text_data_overlap(fig: Figure) -> list[dict[str, object]]:
    renderer = fig.canvas.get_renderer()
    spacing_px = 0.2 * fig.dpi / 25.4
    failures = []
    for ax in fig.axes:
        if not ax.get_visible():
            continue
        texts = []
        for text in _audio_texts(fig):
            if getattr(text, "_qa_cell_label", False):
                continue
            bbox = text.get_window_extent(renderer)
            center = ((bbox.x0 + bbox.x1) / 2, (bbox.y0 + bbox.y1) / 2)
            if ax.bbox.contains(*center):
                texts.append((text, bbox))
        if not texts:
            continue

        geometry = []
        for line in ax.lines:
            geometry.extend(_line_geometry(line, spacing_px))
        for collection in ax.collections:
            if isinstance(collection, (LineCollection, PathCollection, PolyCollection)):
                geometry.extend(_collection_geometry(collection, spacing_px))
        for patch in ax.patches:
            if not patch.get_visible():
                continue
            patch_bbox = patch.get_window_extent(renderer)
            geometry.append(patch_bbox)

        for text, text_bbox in texts:
            for item in geometry:
                if hasattr(item, "x0"):
                    hit = _bbox_intersection(text_bbox, item, tolerance=0.05)
                elif len(item) == 3:
                    x, y, radius = item
                    hit = (
                        text_bbox.x0 < x + radius
                        and text_bbox.x1 > x - radius
                        and text_bbox.y0 < y + radius
                        and text_bbox.y1 > y - radius
                    )
                else:
                    x, y = item
                    hit = (
                        text_bbox.x0 < x < text_bbox.x1
                        and text_bbox.y0 < y < text_bbox.y1
                    )
                if hit:
                    failures.append(
                        {
                            "axis": ax.get_label() or f"axes@{ax.get_position().bounds}",
                            "text": text.get_text(),
                            "text_bbox": [
                                round(text_bbox.x0, 2),
                                round(text_bbox.y0, 2),
                                round(text_bbox.x1, 2),
                                round(text_bbox.y1, 2),
                            ],
                        }
                    )
                    break
    return failures


def _path_line_geometry(transform, path, spacing_px: float):
    """Convert path line segments to sampled display coordinates."""
    vertices = transform.transform_path(path).vertices
    if len(vertices) < 2:
        return
    for p0, p1 in zip(vertices[:-1], vertices[1:], strict=False):
        yield from _sample_segment(p0, p1, spacing_px)


def _line_text_intersection_check(fig: Figure) -> dict[str, object]:
    """Visible lines, patch edges, and axes borders must not cross text boxes."""
    renderer = fig.canvas.get_renderer()
    spacing_px = 0.05 * fig.dpi / 25.4
    text_boxes = [
        (text, text.get_window_extent(renderer)) for text in _audio_texts(fig)
    ]
    checked_artists = 0
    failures = []

    for artist in fig.findobj():
        if not artist.get_visible():
            continue
        geometry = []
        if isinstance(artist, Line2D):
            vertices = artist.get_transform().transform(artist.get_xydata())
            if artist.get_linestyle() not in {"None", "none", "", " "}:
                for p0, p1 in zip(vertices[:-1], vertices[1:], strict=False):
                    geometry.extend(_sample_segment(p0, p1, spacing_px))
        elif isinstance(artist, Collection):
            if isinstance(artist, PathCollection):
                offsets = artist.get_offset_transform().transform(
                    np.asarray(artist.get_offsets())
                )
                sizes = artist.get_sizes()
                for index, (x, y) in enumerate(offsets):
                    area = (
                        float(sizes[index] if len(sizes) > 1 else sizes[0])
                        if len(sizes)
                        else 9.0
                    )
                    geometry.append(
                        (x, y, math.sqrt(max(area, 0.0)) / 2 * fig.dpi / 72)
                    )
            else:
                transform = artist.get_transform()
                for segment in artist.get_segments():
                    points = transform.transform(np.asarray(segment))
                    for p0, p1 in zip(points[:-1], points[1:], strict=False):
                        geometry.extend(_sample_segment(p0, p1, spacing_px))
        elif isinstance(artist, Patch):
            edge_rgba = to_rgb(artist.get_edgecolor())
            has_edge = (
                artist.get_linewidth() > 0
                and artist.get_linestyle() not in {"None", "none", "", " "}
                and max(edge_rgba) > 0
            )
            if has_edge:
                geometry.extend(_path_line_geometry(
                    artist.get_transform(),
                    artist.get_path(),
                    spacing_px,
                ))

        if not geometry:
            continue
        checked_artists += 1
        for text, box in text_boxes:
            for item in geometry:
                if len(item) == 3:
                    x, y, radius = item
                    hit = (
                        box.x0 < x + radius
                        and box.x1 > x - radius
                        and box.y0 < y + radius
                        and box.y1 > y - radius
                    )
                else:
                    x, y = item
                    hit = box.x0 < x < box.x1 and box.y0 < y < box.y1
                if hit:
                    failures.append(
                        {
                            "artist": type(artist).__name__,
                            "text": text.get_text(),
                            "text_bbox_mm": [
                                round(value, 4)
                                for value in (
                                    box.x0 * 25.4 / fig.dpi,
                                    box.y0 * 25.4 / fig.dpi,
                                    box.x1 * 25.4 / fig.dpi,
                                    box.y1 * 25.4 / fig.dpi,
                                )
                            ],
                        }
                    )
                    break
            else:
                continue
            break

    return {
        "artists_checked": checked_artists,
        "texts_checked": len(text_boxes),
        "failures": failures,
        "passed": not failures,
    }


def _axis_clip_failures(fig: Figure) -> list[dict[str, object]]:
    failures = []
    for ax in fig.axes:
        for x, y in getattr(ax, "_qa_points", []):
            display = ax.transData.transform((x, y))
            if not (
                ax.bbox.x0 - 0.1 <= display[0] <= ax.bbox.x1 + 0.1
                and ax.bbox.y0 - 0.1 <= display[1] <= ax.bbox.y1 + 0.1
            ):
                failures.append(
                    {
                        "axis": f"axes@{ax.get_position().bounds}",
                        "data": [x, y],
                        "display": [round(float(display[0]), 2), round(float(display[1]), 2)],
                        "axes_bbox": [
                            round(ax.bbox.x0, 2),
                            round(ax.bbox.y0, 2),
                            round(ax.bbox.x1, 2),
                            round(ax.bbox.y1, 2),
                        ],
                    }
                )
    return failures


def _rotation_failures(fig: Figure) -> list[list[object]]:
    failures = []
    for ax in fig.axes:
        for text in [*ax.get_xticklabels(), *ax.get_yticklabels()]:
            if not text.get_visible() or not text.get_text().strip():
                continue
            rotation = float(text.get_rotation()) % 360
            if not math.isclose(rotation, 0.0, abs_tol=1e-6):
                failures.append([ax.get_label() or str(ax.get_position().bounds), text.get_text(), rotation])
    return failures


def _swarm_spacing_check(fig: Figure) -> dict[str, object]:
    """Check centre distances on a swarm row in physical millimetres."""
    rows: dict[str, dict[str, object]] = {}
    for ax in fig.axes:
        for swarm in getattr(ax, "_qa_swarms", []):
            key = str(swarm["row_key"])
            display = ax.transData.transform(np.asarray(swarm["points"], dtype=float))
            points_mm = display * 25.4 / fig.dpi
            entry = rows.setdefault(
                key,
                {
                    "diameter_pt": float(swarm["marker_diameter_pt"]),
                    "points": [],
                },
            )
            entry["points"].extend(
                {
                    "value": float(data_point[0]),
                    "y": float(data_point[1]),
                    "xy_mm": (float(mm_point[0]), float(mm_point[1])),
                }
                for data_point, mm_point in zip(swarm["points"], points_mm, strict=True)
            )

    minimum = None
    minimum_required = None
    failures = []
    point_count = 0
    for row_key, entry in sorted(rows.items()):
        points = entry["points"]
        point_count += len(points)
        required = 0.95 * float(entry["diameter_pt"]) * 25.4 / 72
        for first, second in itertools.combinations(points, 2):
            distance = float(np.hypot(
                first["xy_mm"][0] - second["xy_mm"][0],
                first["xy_mm"][1] - second["xy_mm"][1],
            ))
            if minimum is None or distance < minimum:
                minimum = distance
                minimum_required = required
            if distance < required - 1e-9:
                failures.append(
                    {
                        "row": row_key,
                        "positions": [
                            [first["value"], first["y"]],
                            [second["value"], second["y"]],
                        ],
                        "distance_mm": round(distance, 6),
                        "required_mm": round(required, 6),
                    }
                )
    return {
        "rows_checked": len(rows),
        "points_checked": point_count,
        "minimum_distance_mm": round(minimum, 6) if minimum is not None else None,
        "required_distance_mm": round(minimum_required, 6) if minimum_required is not None else None,
        "failures": failures,
        "passed": not failures,
    }


def _text_pdf_bboxes(fig: Figure):
    renderer = fig.canvas.get_renderer()
    page_height = fig.get_figheight() * 72
    boxes = []
    for text in _audio_texts(fig):
        bbox = text.get_window_extent(renderer)
        pdf = pymupdf.Rect(
            bbox.x0 / fig.dpi * 72,
            page_height - bbox.y1 / fig.dpi * 72,
            bbox.x1 / fig.dpi * 72,
            page_height - bbox.y0 / fig.dpi * 72,
        )
        boxes.append(pdf)
    return boxes


def _pdf_text_checks(page, fig: Figure) -> dict[str, object]:
    text_boxes = page.get_text("dict")
    line_boxes = []
    for block in text_boxes["blocks"]:
        if block["type"] != 0:
            continue
        for line in block["lines"]:
            bbox = pymupdf.Rect(line["bbox"])
            shrunk = pymupdf.Rect(
                bbox.x0 + 0.15,
                bbox.y0 + 0.15,
                bbox.x1 - 0.15,
                bbox.y1 - 0.15,
            )
            if shrunk.width > 0 and shrunk.height > 0:
                line_boxes.append((bbox, shrunk))

    source_boxes = _text_pdf_bboxes(fig)
    grouped = []
    for bbox, shrunk in line_boxes:
        center = pymupdf.Point((bbox.x0 + bbox.x1) / 2, (bbox.y0 + bbox.y1) / 2)
        candidates = [
            (index, source_box.get_area())
            for index, source_box in enumerate(source_boxes)
            if pymupdf.Rect(
                source_box.x0 - 0.5,
                source_box.y0 - 0.5,
                source_box.x1 + 0.5,
                source_box.y1 + 0.5,
            ).contains(center)
        ]
        group = f"matplotlib:{min(candidates, key=lambda item: item[1])[0]}" if candidates else f"pdf:{len(grouped)}"
        grouped.append((bbox, shrunk, group))

    overlaps = []
    for first, second in itertools.combinations(grouped, 2):
        if first[2] == second[2]:
            continue
        if _bbox_intersection(first[1], second[1], tolerance=0.0):
            overlaps.append([list(first[0]), list(second[0]), first[2], second[2]])

    spans = []
    for block in text_boxes["blocks"]:
        if block["type"] != 0:
            continue
        for line in block["lines"]:
            spans.extend(line["spans"])
    fonts = sorted({span["font"] for span in spans})
    font_failures = [font for font in fonts if "arial" not in font.lower()]
    minimum_font = min((float(span["size"]) for span in spans), default=0.0)
    return {
        "line_boxes": len(line_boxes),
        "overlaps": overlaps,
        "font_names": fonts,
        "non_arial_fonts": font_failures,
        "minimum_font_pt": minimum_font,
    }


def _pdf_line_widths(page) -> list[float]:
    widths = []
    for drawing in page.get_drawings():
        kind = str(drawing.get("type", ""))
        width = float(drawing.get("width") or 0)
        if "s" in kind and width > 1e-8:
            widths.append(width)
    return widths


def _assertion(name: str, expected, actual, passed: bool | None = None) -> dict[str, object]:
    if passed is None:
        passed = expected == actual
    return {
        "name": name,
        "expected": expected,
        "actual": actual,
        "passed": bool(passed),
    }


def _pdf_text_sequence(page) -> list[str]:
    return [
        span["text"]
        for block in page.get_text("dict")["blocks"]
        if block["type"] == 0
        for line in block["lines"]
        for span in line["spans"]
        if span["text"]
    ]


def _crop_check(fig_no: int, page, audit: dict[str, object]) -> dict[str, object]:
    content_bottom_mm = float(audit["content_bottom_mm"])
    expected_height_mm = math.ceil((content_bottom_mm + 1.5) * 2) / 2
    actual_height_mm = float(audit["height_mm"])
    source_size_mm = list(audit["source_size_mm"])
    text_after = _pdf_text_sequence(page)
    text_sequence_unchanged = text_after == audit["text_before"]
    source_size_ok = (
        abs(source_size_mm[0] - 170.0) <= 0.1
        and abs(source_size_mm[1] - SOURCE_HEIGHT_MM[fig_no]) <= 0.1
    )
    cropped_strip_empty = int(audit["cropped_strip_item_count"]) == 0
    height_formula_ok = math.isclose(actual_height_mm, expected_height_mm, abs_tol=1e-9)
    return {
        "height_mm": round(actual_height_mm, 6),
        "source_size_mm": source_size_mm,
        "source_size_ok": source_size_ok,
        "content_bottom_mm": round(content_bottom_mm, 6),
        "bottom_padding_mm": round(float(audit["bottom_padding_mm"]), 6),
        "height_formula_ok": height_formula_ok,
        "cropped_strip_empty": cropped_strip_empty,
        "cropped_strip_item_count": int(audit["cropped_strip_item_count"]),
        "text_span_count": len(text_after),
        "text_sequence_unchanged": text_sequence_unchanged,
        "passed": (
            source_size_ok
            and cropped_strip_empty
            and height_formula_ok
            and text_sequence_unchanged
        ),
    }


def figure_2_supplementary(
    pdf_path: Path,
    measurements: list[dict[str, float]],
) -> dict[str, object]:
    """Check the six Figure 2 structure slots, shared bond length, and key text."""
    expected_heavy_atoms = [31, 12, 22, 9, 17, 17]
    expected_dummy_atoms = [0, 0, 2, 0, 0, 0]
    expected_highlight_atoms = [0, 0, 0, 0, 9, 0]
    expected_highlight_bonds = [0, 0, 0, 0, 10, 0]
    heavy_atoms = [int(item["heavy_atoms"]) for item in measurements]
    dummy_atoms = [int(item["dummy_atoms"]) for item in measurements]
    highlight_atoms = [
        int(item.get("highlight_atom_count", 0)) for item in measurements
    ]
    highlight_bonds = [
        int(item.get("highlight_bond_count", 0)) for item in measurements
    ]
    nitrogen_hydrogens = [
        int(item["max_nitrogen_hydrogens"]) for item in measurements
    ]
    bond_lengths = [float(item["printed_bond_length_mm"]) for item in measurements]
    deviations = [abs(value - 2.4) / 2.4 for value in bond_lengths]
    centers = [17.0, 47.0, 77.0, 108.5, 133.5, 158.5]
    slots = [
        (2.0, 32.0),
        (32.0, 62.0),
        (62.0, 92.0),
        (96.0, 121.0),
        (121.0, 146.0),
        (146.0, 171.0),
    ]
    structure_layout = [
        {
            "center_x_mm": center,
            "left_mm": center - float(item["width_mm"]) / 2,
            "right_mm": center + float(item["width_mm"]) / 2,
            "height_mm": float(item["height_mm"]),
            "slot_mm": slot,
            "within_slot": (
                slot[0] <= center - float(item["width_mm"]) / 2
                and center + float(item["width_mm"]) / 2 <= slot[1]
                and float(item["height_mm"]) <= 20.0
            ),
        }
        for item, center, slot in zip(measurements, centers, slots, strict=True)
    ]
    required_text = (
        "Activity class",
        "Class support",
        "Neighbor recall (k = 10)",
        "Supported: 506 of 889",
        "Mean 0.585",
        "Fixed k 0.536",
        "Wide record range",
        "wavy bonds: attachment points",
        "Most frequent scaffold",
        "31 heavy atoms",
        "in 5 of 42 members",
        "Coverage-target core",
        "12 heavy atoms",
        "in 42 of 42 members",
        "FMCS core (10-s budget)",
        "22 heavy atoms",
        "in 39 of 42 members",
        "9 heavy atoms",
        "in 27 of 28 members",
        "CHEMBL81688",
        "group center",
        "contains the core",
        "CHEMBL85636",
        "benzothiophene",
        "lacks the core",
    )
    with pymupdf.open(pdf_path) as pdf:
        page = pdf[0]
        page_text = page.get_text("text")
        wavy_boxes = page.search_for("wavy bonds: attachment points")
        page_height_mm = page.rect.height * 25.4 / 72
        cropped_bottom_mm = SOURCE_HEIGHT_MM[2] - page_height_mm
        structure_boxes = [
            pymupdf.Rect(
                (float(item["center_x_mm"]) - float(item["width_mm"]) / 2) * 72 / 25.4,
                (
                    SOURCE_HEIGHT_MM[2]
                    - float(item["center_y_mm"])
                    - float(item["height_mm"]) / 2
                    - cropped_bottom_mm
                ) * 72 / 25.4,
                (float(item["center_x_mm"]) + float(item["width_mm"]) / 2) * 72 / 25.4,
                (
                    SOURCE_HEIGHT_MM[2]
                    - float(item["center_y_mm"])
                    + float(item["height_mm"]) / 2
                    - cropped_bottom_mm
                ) * 72 / 25.4,
            )
            for item in measurements
        ]
    wavy_structure_overlaps = [
        index
        for index, structure_box in enumerate(structure_boxes)
        if any(_bbox_intersection(wavy_box, structure_box) for wavy_box in wavy_boxes)
    ]
    wavy_center_x_mm = (
        (wavy_boxes[0].x0 + wavy_boxes[0].x1) / 2 * 25.4 / 72
        if len(wavy_boxes) == 1
        else None
    )
    wavy_width_mm = (
        wavy_boxes[0].width * 25.4 / 72
        if len(wavy_boxes) == 1
        else 0.0
    )
    wavy_overlap_mm = (
        max(
            0.0,
            min(wavy_boxes[0].x1 * 25.4 / 72, 92.0)
            - max(wavy_boxes[0].x0 * 25.4 / 72, 62.0),
        )
        if len(wavy_boxes) == 1
        else 0.0
    )
    wavy_in_fmcs_column = (
        len(wavy_boxes) == 1
        and wavy_center_x_mm is not None
        and 62.0 <= wavy_center_x_mm <= 92.0
        and wavy_width_mm > 0
        and wavy_overlap_mm / wavy_width_mm >= 0.9
    )
    missing_text = [text for text in required_text if text not in page_text]
    passed = (
        len(measurements) == 6
        and heavy_atoms == expected_heavy_atoms
        and dummy_atoms == expected_dummy_atoms
        and highlight_atoms == expected_highlight_atoms
        and highlight_bonds == expected_highlight_bonds
        and nitrogen_hydrogens[2] == 0
        and max(deviations, default=1.0) <= 0.03
        and all(item["within_slot"] for item in structure_layout)
        and wavy_in_fmcs_column
        and not wavy_structure_overlaps
        and not missing_text
    )
    return {
        "name": "图 2 结构、统一键长与关键文字",
        "expected": {
            "heavy_atoms": expected_heavy_atoms,
            "dummy_atoms": expected_dummy_atoms,
            "highlight_atoms": expected_highlight_atoms,
            "highlight_bonds": expected_highlight_bonds,
            "fmcs_max_nitrogen_hydrogens": 0,
            "printed_bond_length_mm": 2.4,
            "maximum_fractional_deviation": 0.03,
        },
        "actual": {
            "heavy_atoms": heavy_atoms,
            "dummy_atoms": dummy_atoms,
            "highlight_atoms": highlight_atoms,
            "highlight_bonds": highlight_bonds,
            "max_nitrogen_hydrogens": nitrogen_hydrogens,
            "printed_bond_length_mm": bond_lengths,
            "maximum_fractional_deviation": max(deviations, default=None),
            "structure_layout": structure_layout,
            "wavy_box_count": len(wavy_boxes),
            "wavy_in_fmcs_column": wavy_in_fmcs_column,
            "wavy_center_x_mm": wavy_center_x_mm,
            "wavy_column_overlap_fraction": (
                wavy_overlap_mm / wavy_width_mm
                if wavy_width_mm > 0
                else 0.0
            ),
            "wavy_structure_overlaps": wavy_structure_overlaps,
            "missing_text": missing_text,
        },
        "passed": passed,
    }


def _description_layout_check(fig: Figure) -> dict[str, object]:
    """Check the six descriptions, including the fourth FMCS caption line."""
    renderer = fig.canvas.get_renderer()
    texts = [
        text
        for text in fig.findobj(mtext.Text)
        if hasattr(text, "_qa_description_block")
    ]
    blocks: dict[int, dict[str, object]] = {}
    slot_failures = []
    baseline_failures = []
    for text in texts:
        block = int(text._qa_description_block)  # type: ignore[attr-defined]
        line = int(text._qa_description_line)  # type: ignore[attr-defined]
        slot = tuple(text._qa_description_slot_mm)  # type: ignore[attr-defined]
        bbox = text.get_window_extent(renderer)
        bbox_mm = (
            bbox.x0 * 25.4 / fig.dpi,
            bbox.x1 * 25.4 / fig.dpi,
        )
        entry = blocks.setdefault(
            block,
            {
                "bbox": bbox,
                "lines": {},
                "slot_mm": slot,
            },
        )
        entry["bbox"] = pymupdf.Rect(
            min(entry["bbox"].x0, bbox.x0),
            min(entry["bbox"].y0, bbox.y0),
            max(entry["bbox"].x1, bbox.x1),
            max(entry["bbox"].y1, bbox.y1),
        )
        entry["lines"][line] = float(text.get_position()[1])
        if (
            not getattr(text, "_qa_fig2_fmcs_wavy", False)
            and (bbox_mm[0] < slot[0] - 1e-6 or bbox_mm[1] > slot[1] + 1e-6)
        ):
            slot_failures.append(
                {
                    "block": block,
                    "line": line,
                    "text": text.get_text(),
                    "bbox_x_mm": list(bbox_mm),
                    "slot_mm": list(slot),
                }
            )

    overlaps = []
    entries = sorted(blocks.items())
    for (first_id, first), (second_id, second) in itertools.combinations(entries, 2):
        if _bbox_intersection(first["bbox"], second["bbox"], tolerance=0.05):
            overlaps.append([first_id, second_id])
    if len(blocks) == 6:
        for line in range(3):
            baselines = [entry["lines"].get(line) for entry in blocks.values()]
            if any(value is None for value in baselines) or not np.allclose(
                baselines,
                baselines[0],
                atol=1e-12,
            ):
                baseline_failures.append({"line": line, "positions": baselines})
    else:
        baseline_failures.append({"line": None, "positions": len(blocks)})

    wavy_texts = [
        text
        for text in texts
        if getattr(text, "_qa_fig2_fmcs_wavy", False)
    ]
    wavy_layout_ok = False
    if len(wavy_texts) == 1 and len(blocks) == 6:
        wavy = wavy_texts[0]
        fmcs = blocks[2]
        wavy_layout_ok = bool(
            wavy.get_text() == "wavy bonds: attachment points"
            and int(wavy._qa_description_block) == 2  # type: ignore[attr-defined]
            and int(wavy._qa_description_line) == 3  # type: ignore[attr-defined]
            and tuple(wavy._qa_description_slot_mm) == (62.0, 92.0)  # type: ignore[attr-defined]
            and np.isclose(wavy.get_position()[0], 77.0 / 170.0)
            and 2 in fmcs["lines"]
            and np.isclose(
                fmcs["lines"][2] - fmcs["lines"][3],
                3.0 / SOURCE_HEIGHT_MM[2],
            )
        )
    return {
        "blocks_checked": len(blocks),
        "texts_checked": len(texts),
        "slot_failures": slot_failures,
        "overlaps": overlaps,
        "baseline_failures": baseline_failures,
        "fmcs_wavy_text_checked": len(wavy_texts),
        "fmcs_wavy_layout_ok": wavy_layout_ok,
        "passed": (
            len(blocks) == 6
            and len(texts) == 19
            and not slot_failures
            and not overlaps
            and not baseline_failures
            and wavy_layout_ok
        ),
    }


def _letter_centering_check(png_path: Path, fig: Figure) -> dict[str, object]:
    """Use a 600 dpi PNG to check letter-ink bounding-box centres in Figure 5c."""
    labels = [
        text
        for text in fig.findobj(mtext.Text)
        if hasattr(text, "_qa_cell_center_mm")
    ]
    rgb = np.asarray(Image.open(png_path).convert("RGB"))
    height_px, width_px = rgb.shape[:2]
    pixels_per_mm = width_px / (fig.get_figwidth() * 25.4)
    failures = []
    measurements = []
    maximum_horizontal = 0.0
    maximum_vertical = 0.0
    for label in labels:
        marker = str(label._qa_cell_letter)  # type: ignore[attr-defined]
        center_x_mm, center_y_mm = label._qa_cell_center_mm  # type: ignore[attr-defined]
        center_x_px = center_x_mm * pixels_per_mm
        center_y_px = (
            fig.get_figheight() * 25.4 - center_y_mm
        ) * pixels_per_mm
        half_window_px = 0.8 * pixels_per_mm
        x0 = max(0, int(math.floor(center_x_px - half_window_px)))
        x1 = min(width_px, int(math.ceil(center_x_px + half_window_px)) + 1)
        y0 = max(0, int(math.floor(center_y_px - half_window_px)))
        y1 = min(height_px, int(math.ceil(center_y_px + half_window_px)) + 1)
        window = rgb[y0:y1, x0:x1]
        color_rgb = np.asarray(to_rgb(label.get_color())) * 255
        if float(color_rgb.mean()) >= 220:
            mask = np.all(window >= 220, axis=2)
        else:
            mask = np.all(window <= 80, axis=2)
        ys, xs = np.where(mask)
        if len(xs) == 0:
            failures.append({"letter": marker, "reason": "未找到字母颜色像素"})
            continue
        ink_center_x_px = x0 + (float(xs.min()) + float(xs.max())) / 2
        ink_center_y_px = y0 + (float(ys.min()) + float(ys.max())) / 2
        error_x_mm = abs(ink_center_x_px - center_x_px) / pixels_per_mm
        error_y_mm = abs(ink_center_y_px - center_y_px) / pixels_per_mm
        maximum_horizontal = max(maximum_horizontal, error_x_mm)
        if marker != "Q":
            maximum_vertical = max(maximum_vertical, error_y_mm)
        item = {
            "letter": marker,
            "cell_center_mm": [center_x_mm, center_y_mm],
            "horizontal_error_mm": error_x_mm,
            "vertical_error_mm": error_y_mm,
            "vertical_checked": marker != "Q",
        }
        measurements.append(item)
        if error_x_mm > 0.08 or (marker != "Q" and error_y_mm > 0.08):
            failures.append(item)
    return {
        "cells_checked": len(labels),
        "legend_cells_checked": 5,
        "maximum_horizontal_error_mm": maximum_horizontal,
        "maximum_vertical_error_mm": maximum_vertical,
        "tolerance_mm": 0.08,
        "measurements": measurements,
        "failures": failures,
        "passed": len(labels) == 44 and not failures,
    }


def figure_6_supplementary(supplementary_checks: dict[str, object]) -> dict[str, object]:
    """Check the four-panel Figure 6 build and screenshot inputs."""
    figure_6 = supplementary_checks.get("figure_6", {})
    final = figure_6.get("build") if isinstance(figure_6, dict) else None
    final_ok = bool(isinstance(final, dict) and final.get("passed"))
    layout = final.get("layout_mm", {}) if isinstance(final, dict) else {}
    source_size = final.get("source_size_mm", []) if isinstance(final, dict) else []
    layout_ok = (
        len(source_size) == 2
        and np.allclose(source_size, [170.0, 220.0], atol=0.1)
        and set(layout) == {"a", "b", "c", "d"}
    )
    passed = final_ok and layout_ok
    return _assertion(
        "图 6 四面板构建与 data_fig6 截图核对",
        {
            "final_build_passed": True,
            "source_size_mm": [170.0, 220.0],
            "layout_panels": ["a", "b", "c", "d"],
        },
        {
            "final_build_passed": final_ok,
            "source_size_mm": source_size,
            "layout_mm": layout,
        },
        passed,
    )


def _fig3b_label_gap_check(fig: Figure) -> dict[str, object]:
    """The IC50 x-axis labels and Ki title must be at least 2.5 mm apart."""
    renderer = fig.canvas.get_renderer()
    upper = next(
        ax for ax in fig.axes if getattr(ax, "_qa_fig3b_ic50", False)
    )
    title = next(
        text
        for text in fig.findobj(mtext.Text)
        if getattr(text, "_qa_fig3b_ki_title", False)
    )
    title_box = title.get_window_extent(renderer)
    measurements = []
    for label in upper.get_xticklabels():
        if not label.get_visible() or not label.get_text().strip():
            continue
        label_box = label.get_window_extent(renderer)
        measurements.append(
            {
                "tick_label": label.get_text(),
                "gap_mm": round(_bbox_gap_mm(label_box, title_box, fig.dpi), 6),
            }
        )
    minimum = min(item["gap_mm"] for item in measurements)
    return {
        "minimum_required_mm": 2.5,
        "minimum_gap_mm": minimum,
        "measurements": measurements,
        "passed": bool(measurements) and minimum >= 2.5 - 1e-9,
    }


def _fig3b_legend_check(fig: Figure) -> dict[str, object]:
    """Figure 3b shows only G and the post hoc distance score D."""
    matching = [
        legend
        for legend in fig.legends
        if any(
            text.get_text() == "G with 95% CI"
            for text in legend.get_texts()
        )
    ]
    labels = (
        [text.get_text() for text in matching[0].get_texts()]
        if len(matching) == 1
        else []
    )
    return {
        "labels": labels,
        "expected": ["G with 95% CI", "D (post hoc)"],
        "passed": labels == ["G with 95% CI", "D (post hoc)"],
    }


def _fig5e_layout_check(fig: Figure) -> dict[str, object]:
    """Check vertical spacing for the Figure 5e letter, legend, note, and title."""
    renderer = fig.canvas.get_renderer()
    letter = next(
        text
        for text in fig.findobj(mtext.Text)
        if getattr(text, "_qa_panel_letter", None) == "e"
    )
    caption = next(
        text
        for text in fig.findobj(mtext.Text)
        if getattr(text, "_qa_fig5_c_caption", False)
    )
    contrast_title = next(
        text
        for text in fig.findobj(mtext.Text)
        if getattr(text, "_qa_fig5b_axis_title", False)
    )
    legend = [
        text
        for text in fig.findobj(mtext.Text)
        if getattr(text, "_qa_fig5e_legend", False)
    ]
    titles = [
        text
        for text in fig.findobj(mtext.Text)
        if getattr(text, "_qa_fig5e_title", False)
    ]
    letter_box = letter.get_window_extent(renderer)
    caption_box = caption.get_window_extent(renderer)
    contrast_title_box = contrast_title.get_window_extent(renderer)
    legend_boxes = [text.get_window_extent(renderer) for text in legend]
    title_boxes = [text.get_window_extent(renderer) for text in titles]
    caption_gaps = [
        _bbox_gap_mm(box, caption_box, fig.dpi) for box in legend_boxes
    ]
    title_gaps = [
        _bbox_gap_mm(legend_box, title_box, fig.dpi)
        for legend_box in legend_boxes
        for title_box in title_boxes
    ]
    letter_gaps = [
        _bbox_gap_mm(letter_box, legend_box, fig.dpi)
        for legend_box in legend_boxes
    ]
    contrast_title_gap = _bbox_gap_mm(
        letter_box,
        contrast_title_box,
        fig.dpi,
    )
    result = {
        "minimum_required_mm": 2.0,
        "legend_items_checked": len(legend),
        "titles_checked": len(titles),
        "minimum_gap_to_CHEMBL81688_mm": round(min(caption_gaps), 6),
        "minimum_gap_to_panel_titles_mm": round(min(title_gaps), 6),
        "letter_e_to_fig5b_axis_title_mm": round(contrast_title_gap, 6),
        "minimum_gap_from_letter_e_mm": round(min(letter_gaps), 6),
    }
    result["passed"] = (
        len(legend) == 5
        and len(titles) == 3
        and result["minimum_gap_to_CHEMBL81688_mm"] >= 2.0 - 1e-9
        and result["minimum_gap_to_panel_titles_mm"] >= 2.0 - 1e-9
        and result["letter_e_to_fig5b_axis_title_mm"] >= 2.0 - 1e-9
        and result["minimum_gap_from_letter_e_mm"] >= 0.5 - 1e-9
    )
    return result


def _source_specs():
    return {
        "3a_B0_disagreement.csv": c.DATA_DIR / "3a_B0_disagreement.csv",
        "3b_IC50.csv": c.DATA_DIR / "3b_IC50.csv",
        "3b_Ki.csv": c.DATA_DIR / "3b_Ki.csv",
        "3c_IC50.csv": c.DATA_DIR / "3c_IC50.csv",
        "3c_Ki.csv": c.DATA_DIR / "3c_Ki.csv",
        "4a_calibration_by_class.csv": c.DATA_DIR / "4a_calibration_by_class.csv",
        "4b_hERG_patent_pair.csv": c.DATA_DIR / "4b_hERG_patent_pair.csv",
        "4c_IC50.csv": c.DATA_DIR / "4c_IC50.csv",
        "4c_Ki.csv": c.DATA_DIR / "4c_Ki.csv",
        "4d_IC50.csv": c.DATA_DIR / "4d_IC50.csv",
        "4d_Ki.csv": c.DATA_DIR / "4d_Ki.csv",
        "target_symbols.csv": c.DATA_DIR / "target_symbols.csv",
        "5a_recall.csv": c.DATA_DIR / "5a_recall.csv",
        "5a_fixedk_means.csv": c.DATA_DIR / "5a_fixedk_means.csv",
        "5b_layout_contrasts.csv": c.DATA_DIR / "5b_layout_contrasts.csv",
        "5d_layout_correlation.csv": c.DATA_DIR / "5d_layout_correlation.csv",
        "5e_core_coverage_by_group_size.csv": c.DATA_DIR / "5e_core_coverage_by_group_size.csv",
        "Figure_3C_CHEMBL81688_neighbors.csv": c.CGM_DIR / "paper/output/reliability_v1/figure_data/Figure_3C_CHEMBL81688_neighbors.csv",
        "coordinates_seed_42.csv": c.CGM_DIR / "paper/output/reliability_v1/layouts/chembl205/coordinates_seed_42.csv",
    }


def _check_plotted_values() -> dict[str, object]:
    sources = {
        name: pd.read_csv(path, float_precision="round_trip")
        for name, path in _source_specs().items()
    }
    expected: dict[tuple[str, str], list[float]] = {}

    def add(panel: str, variable: str, values):
        expected.setdefault((panel, variable), []).extend(
            float(value) for value in np.asarray(values, dtype=float).ravel() if np.isfinite(value)
        )

    frame = sources["3a_B0_disagreement.csv"]
    add("3a", "IC50", frame["IC50"].dropna())
    add("3a", "Ki", frame["Ki"].dropna())

    for endpoint in ("IC50", "Ki"):
        frame = sources[f"3b_{endpoint}.csv"]
        for column in ("g_ci_low", "g_ci_high", "g_auc", "d_auc"):
            add("3b", column, frame[column])

        frame = sources[f"3c_{endpoint}.csv"]
        target = frame[frame["row_type"] == "target"]
        summary = frame[frame["row_type"] == "summary"].iloc[0]
        for column in (
            "g_minus_record_log_loss",
            "g_minus_global_log_loss",
            "g_minus_fixed_log_loss",
            "g_minus_x_log_loss",
        ):
            add("3c", column, target[column].dropna())
            for prefix in ("median", "ci_low", "ci_high"):
                add("3c", f"{prefix}_{column}", [summary[f"{prefix}_{column}"]])

        frame = sources["4a_calibration_by_class.csv"]
        frame = frame[
            frame["endpoint"].eq(endpoint) & ~frame["omit_in_figure"].astype(bool)
        ]
        for column in ("mean_p", "agreement_rate", "wilson_ci_low", "wilson_ci_high"):
            add("4a", column, frame[column])

        frame = sources[f"4c_{endpoint}.csv"]
        for column in (
            "supported_agreement", "supported_ci_low", "supported_ci_high",
            "uncertain_agreement", "uncertain_ci_low", "uncertain_ci_high",
        ):
            add("4c", column, frame[column])

        frame = sources[f"4d_{endpoint}.csv"]
        frame = frame[frame["n_predictions"] >= 8]
        for column in (
            "agreement_rate",
            "agreement_ci_low",
            "agreement_ci_high",
            "mean_probability",
            "n_predictions",
            "n_datasets",
        ):
            add("4d", column, frame[column])

    frame = sources["4b_hERG_patent_pair.csv"]
    for column in ("us_10131663_b2", "us_10590130_b2", "offset_b_minus_a"):
        add("4b", column, frame[column])

    frame = sources["5a_recall.csv"]
    for column in frame.columns:
        if "__r" in column:
            add("5a", column, frame[column])
    frame = sources["5a_fixedk_means.csv"]
    add(
        "5a",
        "fixed_recall_mean",
        frame.loc[
            frame["layout"].isin(
                ["projection", "minimum_displacement", "neighbor_refined"]
            ),
            "fixed_recall_mean",
        ],
    )

    frame = sources["5b_layout_contrasts.csv"]
    for column in ("other_minus_refined_e_grid", "ci_low", "ci_high"):
        add("5b", column, frame[column])

    frame = sources["5d_layout_correlation.csv"]
    for column in ("rho", "rho_ci_low", "rho_ci_high"):
        add("5d", column, frame[column])

    frame = sources["5e_core_coverage_by_group_size.csv"]
    for column in (
        "coverage_ge_0_8", "coverage_gt_0_5_lt_0_8", "coverage_le_0_5",
        "no_core", "groups", "mean_core_heavy_atoms",
    ):
        add("5e", column, frame[column])

    failures = []
    compared = 0
    for fig_no in (3, 4, 5):
        plotted_path = c.QA_DIR / f"Figure_{fig_no}_plotted_values.csv"
        plotted = pd.read_csv(plotted_path, float_precision="round_trip")
        for (panel, variable), expected_values in sorted(expected.items()):
            if not any(row["panel"] == panel for _, row in plotted.iterrows()):
                continue
            actual_values = plotted[(plotted["panel"] == panel) & (plotted["variable"] == variable)]["value"].to_numpy(dtype=float)
            if len(actual_values) == 0:
                continue
            compared += len(expected_values)
            if sorted(actual_values.tolist()) != sorted(expected_values):
                failures.append(
                    {
                        "figure": fig_no,
                        "panel": panel,
                        "variable": variable,
                        "expected_count": len(expected_values),
                        "actual_count": len(actual_values),
                        "max_abs_difference": float(
                            np.max(np.abs(np.sort(actual_values) - np.sort(expected_values)))
                        )
                        if len(actual_values) == len(expected_values)
                        else None,
                    }
                )
    return {"values_compared": compared, "mismatches": failures, "passed": not failures}


def _source_assertions() -> list[dict[str, object]]:
    results = []
    frame_3a = c.load_csv("3a_B0_disagreement.csv")
    for endpoint, median, q1, q3 in (
        ("IC50", 0.2879, 0.1961, 0.3507),
        ("Ki", 0.1636, 0.1214, 0.2257),
    ):
        values = frame_3a[endpoint].dropna()
        actual = [float(values.median()), float(values.quantile(0.25)), float(values.quantile(0.75))]
        results.append(_assertion(f"3a {endpoint} median/quartiles", [median, q1, q3], actual, all(round(a, 4) == e for a, e in zip(actual, (median, q1, q3), strict=True))))

    symbol_frame = c.load_csv("target_symbols.csv")
    symbol_ids = set(symbol_frame["target_id"])
    plotted_3 = pd.read_csv(
        c.QA_DIR / "Figure_3_plotted_values.csv",
        float_precision="round_trip",
    )
    results.append(
        _assertion(
            "3b no B-range-2 plotted values",
            [],
            sorted(
                plotted_3.loc[
                    plotted_3["panel"].eq("3b")
                    & plotted_3["variable"].eq("b_range_2_auc"),
                    "row_identifier",
                ].tolist()
            ),
        )
    )
    for endpoint, expected_count in (("IC50", 29), ("Ki", 23)):
        frame = c.load_csv(f"3b_{endpoint}.csv")
        results.append(_assertion(f"3b {endpoint} rows", expected_count, len(frame)))
        results.append(_assertion(f"3b {endpoint} rank order", list(range(1, expected_count + 1)), sorted(frame["target_rank"].astype(int).tolist())))
        results.append(
            _assertion(
                f"3b {endpoint} target symbols",
                True,
                set(frame["target_id"]) <= symbol_ids,
            )
        )

    columns_3c = (
        "g_minus_record_log_loss",
        "g_minus_global_log_loss",
        "g_minus_fixed_log_loss",
        "g_minus_x_log_loss",
    )
    for endpoint, expected in (("IC50", [29, 29, 29, 15]), ("Ki", [23, 23, 23, 11])):
        frame = c.load_csv(f"3c_{endpoint}.csv")
        target = frame[frame["row_type"] == "target"]
        counts = [int(target[column].notna().sum()) for column in columns_3c]
        finite = target[list(columns_3c)].to_numpy(dtype=float)
        finite = finite[np.isfinite(finite)]
        results.append(_assertion(f"3c {endpoint} row counts", expected, counts))
        results.append(
            _assertion(
                f"3c {endpoint} source range",
                "within [-0.70, 0.20]",
                [float(finite.min()), float(finite.max())],
                finite.min() >= -0.70 and finite.max() <= 0.20,
            )
        )
        gap_count = int(np.sum((finite > -0.40) & (finite < -0.32)))
        results.append(_assertion(f"3c {endpoint} empty break", 0, gap_count))
        extended = finite[finite <= -0.40]
        results.append(
            _assertion(
                f"3c {endpoint} extended segment",
                "every point within [-0.70, -0.40]",
                [float(extended.min()), float(extended.max())],
                bool(len(extended)) and extended.min() >= -0.70,
            )
        )

    calibration = c.load_csv("4a_calibration_by_class.csv")
    for endpoint, expected in (
        ("IC50", {"active": 4, "inactive": 5, "medium": 3}),
        ("Ki", {"active": 3, "inactive": 4, "medium": 4}),
    ):
        local = calibration[
            calibration["endpoint"].eq(endpoint)
            & ~calibration["omit_in_figure"].astype(bool)
        ]
        actual = local.groupby("display_class").size().to_dict()
        results.append(_assertion(f"4a {endpoint} plotted bins", expected, actual))
    results.append(_assertion("4a omitted bins", 7, int(calibration["omit_in_figure"].sum())))
    si_4a = c.load_csv("si_4a_without_CHEMBL5113.csv")
    results.extend(
        [
            _assertion("SI 4a rows", 20, len(si_4a)),
            _assertion(
                "SI 4a omitted bins",
                2,
                int(si_4a["omit_in_figure"].sum()),
            ),
        ]
    )

    pair = c.load_csv("4b_hERG_patent_pair.csv")
    median_offset = float(pair["offset_b_minus_a"].median())
    spearman = float(pair[["us_10131663_b2", "us_10590130_b2"]].corr(method="spearman").iloc[0, 1])
    results.extend(
        [
            _assertion("4b row count", 30, len(pair)),
            _assertion("4b median offset", 3.47, round(median_offset, 2)),
            _assertion("4b Spearman rho", 0.20, round(spearman, 2)),
        ]
    )

    for endpoint, expected_count in (("IC50", 18), ("Ki", 15)):
        frame = c.load_csv(f"4c_{endpoint}.csv")
        results.append(_assertion(f"4c {endpoint} rows", expected_count, len(frame)))
    for endpoint in ("IC50", "Ki"):
        frame = c.load_csv(f"4d_{endpoint}.csv")
        results.extend(
            [
                _assertion(f"4d {endpoint} rows", 6, len(frame)),
                _assertion(
                    f"4d {endpoint} states",
                    ["P_supported", "touched_all", "uncertain_low_support"],
                    sorted(frame["state"].unique().tolist()),
                ),
                _assertion(
                    f"4d {endpoint} omitted n<8",
                    ["touched_all"] if endpoint == "IC50" else ["touched_all"],
                    frame.loc[frame["n_predictions"] < 8, "state"].unique().tolist(),
                ),
                _assertion(
                    f"4d {endpoint} dataset counts",
                    True,
                    bool(frame["n_datasets"].notna().all()),
                ),
            ]
        )

    frame_5a = c.load_csv("5a_recall.csv")
    expected_means = {
        "CA II": [0.657, 0.504, 0.589],
        "Thrombin": [0.644, 0.490, 0.586],
        "hERG": [0.500, 0.402, 0.488],
    }
    layouts = ["projection", "minimum_displacement", "neighbor_refined"]
    for _, row in frame_5a.iterrows():
        means = [
            round(float(row[[f"{layout}__r{i:02d}" for i in range(1, 11)]].mean()), 3)
            for layout in layouts
        ]
        results.append(_assertion(f"5a {row['Row_title']} means", expected_means[row["Row_title"]], means))
    fixed = c.load_csv("5a_fixedk_means.csv")
    results.append(_assertion("5a fixed-k rows", 18, len(fixed)))
    results.append(_assertion("5a fixed-k seeds", [10], sorted(fixed["n_seeds"].unique().tolist())))

    frame_5b = c.load_csv("5b_layout_contrasts.csv")
    results.append(_assertion("5b row count", 12, len(frame_5b)))
    hERG_las = frame_5b[(frame_5b["target_label"] == "hERG") & (frame_5b["comparison"] == "las_random_minus_neighbor_refined")].iloc[0]["other_minus_refined_e_grid"]
    results.append(_assertion("5b hERG LAS random value", -0.0240002267509859, float(hERG_las)))

    neighbors = c.load_csv("Figure_3C_CHEMBL81688_neighbors.csv", root=c.CGM_DIR / "paper/output/reliability_v1/figure_data")
    coordinates = c.load_csv("coordinates_seed_42.csv", root=c.CGM_DIR / "paper/output/reliability_v1/layouts/chembl205")
    coordinate_index = coordinates.set_index("molecule_identity_key")
    results.append(_assertion("5c coordinate rows", 889, len(coordinates)))
    results.append(_assertion("5c neighbor rows", 40, len(neighbors)))
    for layout, expected in (("minimum_displacement", [3, 9, 8]), ("neighbor_refined", [5, 7, 6])):
        local = neighbors[neighbors["layout"] == layout]
        counts = Counter(local["figure_marker"])
        actual = [counts["S"], counts["G"], counts["F"]]
        results.append(_assertion(f"5c {layout} S/G/F", expected, actual))
        col_prefix = "grid_col" if layout == "minimum_displacement" else "refined_grid_col"
        row_prefix = "grid_row" if layout == "minimum_displacement" else "refined_grid_row"
        query = local[local["figure_marker"] == "Q"].iloc[0]
        query_coordinate = coordinate_index.loc[query["molecule_identity_key"]]
        results.append(
            _assertion(
                f"5c {layout} query center",
                [11, 28],
                [
                    int(query_coordinate[col_prefix]),
                    int(query_coordinate[row_prefix]),
                ],
            )
        )
        matches = all(
            (int(coordinate_index.loc[row["molecule_identity_key"], col_prefix]), int(coordinate_index.loc[row["molecule_identity_key"], row_prefix]))
            == (int(row["grid_col"]), int(row["grid_row"]))
            for _, row in local.iterrows()
        )
        results.append(_assertion(f"5c {layout} coordinates match", True, matches))

    frame_5d = c.load_csv("5d_layout_correlation.csv")
    sampled = frame_5d[frame_5d["row_type"] == "sampled_dataset"]
    results.extend(
        [
            _assertion("5d sampled datasets", 77, len(sampled)),
            _assertion("5d IC50 rows", 48, int((sampled["endpoint"] == "IC50").sum())),
            _assertion("5d Ki rows", 29, int((sampled["endpoint"] == "Ki").sum())),
        ]
    )
    for endpoint, expected in (("IC50", -0.090), ("Ki", -0.074)):
        actual = round(float(sampled[sampled["endpoint"] == endpoint]["rho"].median()), 3)
        results.append(_assertion(f"5d {endpoint} median", expected, actual))

    frame_5e = c.load_csv("5e_core_coverage_by_group_size.csv")
    results.append(_assertion("5e rows", 12, len(frame_5e)))
    for rule, local in frame_5e.groupby("rule"):
        results.append(_assertion(f"5e {rule} groups", 480, int(local["groups"].sum())))
        category_sum = local[["coverage_ge_0_8", "coverage_gt_0_5_lt_0_8", "coverage_le_0_5", "no_core"]].sum(axis=1).to_numpy()
        results.append(_assertion(f"5e {rule} category sums", local["groups"].tolist(), category_sum.tolist()))
    return results


def _output_metadata(pdf_path: Path, png_path: Path, preview_path: Path, height_mm: float) -> dict[str, object]:
    with pymupdf.open(pdf_path) as pdf:
        page = pdf[0]
        pdf_size_mm = [page.rect.width * 25.4 / 72, page.rect.height * 25.4 / 72]
    with Image.open(png_path) as image:
        png_size = list(image.size)
        png_dpi = list(image.info.get("dpi", (0, 0)))
    with Image.open(preview_path) as image:
        preview_size = list(image.size)
        preview_dpi = list(image.info.get("dpi", (0, 0)))
    expected_pixels = [round(170 / 25.4 * 600), round(height_mm / 25.4 * 600)]
    return {
        "pdf_size_mm": pdf_size_mm,
        "pdf_size_ok": abs(pdf_size_mm[0] - 170) <= 0.1 and abs(pdf_size_mm[1] - height_mm) <= 0.1,
        "png_size_px": png_size,
        "png_expected_px": expected_pixels,
        "png_size_ok": all(abs(a - b) <= 2 for a, b in zip(png_size, expected_pixels, strict=True)),
        "png_dpi": png_dpi,
        "png_dpi_ok": all(np.isclose(value, 600, atol=0.5) for value in png_dpi),
        "png_rgb": Image.open(png_path).mode == "RGB",
        "preview_size_px": preview_size,
        "preview_dpi": preview_dpi,
    }


def run_qa(
    figures: dict[int, Figure],
    pdf_paths: dict[int, Path],
    png_paths: dict[int, Path],
    preview_paths: dict[int, Path],
    molecule_measurements: dict[str, float] | None = None,
    crop_audits: dict[int, dict[str, object]] | None = None,
    supplementary_checks: dict[str, object] | None = None,
    allowed_raster_counts: dict[int, int] | None = None,
) -> dict[str, object]:
    if crop_audits is None:
        raise ValueError("Missing PDF crop audit")
    figure_results = {}
    swarm_spacing = {}
    for fig_no, fig in figures.items():
        fig.canvas.draw()
        audit = crop_audits[fig_no]
        metadata = _output_metadata(
            pdf_paths[fig_no],
            png_paths[fig_no],
            preview_paths[fig_no],
            float(audit["height_mm"]),
        )
        with pymupdf.open(pdf_paths[fig_no]) as pdf:
            page = pdf[0]
            pdf_text = _pdf_text_checks(page, fig)
            line_widths = _pdf_line_widths(page)
            image_count = len(page.get_images(full=True))
            pdf_size_mm = [page.rect.width * 25.4 / 72, page.rect.height * 25.4 / 72]
            crop = _crop_check(fig_no, page, audit)
        matplotlib_overlaps = _text_overlap_matplotlib(fig)
        text_spacing = _text_spacing_check(fig)
        legend_title_spacing = _legend_title_spacing_check(fig)
        data_overlaps = _text_data_overlap(fig)
        line_text_intersections = _line_text_intersection_check(fig)
        margins = _text_margin_failures(fig)
        text_geometry = _text_geometry_check(fig)
        clipping = _axis_clip_failures(fig)
        rotations = _rotation_failures(fig)
        if any(getattr(ax, "_qa_swarms", []) for ax in fig.axes):
            swarm_spacing[str(fig_no)] = _swarm_spacing_check(fig)
        expected_raster_count = (
            int(allowed_raster_counts.get(fig_no, 0))
            if allowed_raster_counts
            else 0
        )
        figure_results[fig_no] = {
            "size": {
                **metadata,
                "pdf_actual_mm": pdf_size_mm,
                "passed": metadata["pdf_size_ok"] and metadata["png_size_ok"] and metadata["png_dpi_ok"],
            },
            "fonts": {
                **pdf_text,
                "passed": not pdf_text["non_arial_fonts"] and pdf_text["minimum_font_pt"] >= MIN_FONT_PT,
            },
            "text_overlap_matplotlib": {"passed": not matplotlib_overlaps, "failures": matplotlib_overlaps[:80]},
            "text_spacing": {**text_spacing, "failures": text_spacing["failures"][:80]},
            "legend_title_spacing": {
                **legend_title_spacing,
                "failures": legend_title_spacing["failures"][:80],
            },
            "text_overlap_pdf": {"passed": not pdf_text["overlaps"], "failures": pdf_text["overlaps"][:80]},
            "text_data_overlap": {"passed": not data_overlaps, "failures": data_overlaps[:80]},
            "line_text_intersections": {
                **line_text_intersections,
                "failures": line_text_intersections["failures"][:80],
            },
            "canvas_margin": {"minimum_mm": 0.5, "passed": not margins, "failures": margins[:80]},
            "text_geometry": text_geometry,
            "clipping": {"passed": not clipping, "failures": clipping[:80]},
            "tick_rotation": {"passed": not rotations, "failures": rotations[:80]},
            "line_width": {
                "minimum_pt": min(line_widths) if line_widths else None,
                "required_pt": 0.3 if fig_no == 2 else MIN_LINE_PT,
                "passed": (
                    bool(line_widths)
                    and min(line_widths)
                    >= (0.3 if fig_no == 2 else MIN_LINE_PT) - 1e-9
                ),
            },
            "vector": {
                "raster_image_count": image_count,
                "allowed_raster_image_count": expected_raster_count,
                "passed": image_count == expected_raster_count,
            },
            "crop": crop,
        }
        if fig_no == 2:
            figure_results[fig_no]["description_layout"] = _description_layout_check(fig)
            figure_results[fig_no]["fig2c_spacing"] = _fig2c_spacing_check(fig)
            figure_results[fig_no]["fig2_colors"] = _fig2_color_check(fig)
        if fig_no == 4:
            figure_results[fig_no]["fig4a"] = _fig4a_check(fig)
            figure_results[fig_no]["fig4d"] = _fig4d_check(
                fig,
                png_paths[fig_no],
            )
        if fig_no == 3:
            figure_results[fig_no]["fig3b_label_gap"] = _fig3b_label_gap_check(fig)
            figure_results[fig_no]["fig3b_legend"] = _fig3b_legend_check(fig)
        if fig_no == 5:
            figure_results[fig_no]["letter_centering"] = _letter_centering_check(
                png_paths[fig_no],
                fig,
            )
            figure_results[fig_no]["fig5e_layout"] = _fig5e_layout_check(fig)

    data_check = _check_plotted_values()
    assertions = _source_assertions()
    if supplementary_checks:
        assertions.extend(supplementary_checks.get("assertions", []))
        assertions.append(figure_6_supplementary(supplementary_checks))
    assertion_items = {
        "passed": all(item["passed"] for item in assertions),
        "items": assertions,
    }
    check_names = [
        "size",
        "fonts",
        "text_overlap_matplotlib",
        "text_spacing",
        "legend_title_spacing",
        "text_overlap_pdf",
        "text_data_overlap",
        "line_text_intersections",
        "canvas_margin",
        "text_geometry",
        "clipping",
        "tick_rotation",
        "line_width",
        "vector",
        "crop",
    ]
    all_passed = data_check["passed"] and assertion_items["passed"]
    if supplementary_checks:
        all_passed = all_passed and bool(supplementary_checks.get("passed", True))
    all_passed = (
        all_passed
        and all(result["passed"] for result in swarm_spacing.values())
        and all(
            check["passed"]
            for result in figure_results.values()
            for name, check in result.items()
            if name in {
                *check_names,
                "description_layout",
                "fig2c_spacing",
                "fig2_colors",
                "fig4a",
                "fig4d",
                "fig3b_label_gap",
                "fig3b_legend",
                "letter_centering",
                "fig5e_layout",
            }
        )
    )
    report = {
        "all_passed": bool(all_passed),
        "figures": {str(key): value for key, value in figure_results.items()},
        "data_consistency": data_check,
        "swarm_point_spacing": {
            "passed": all(result["passed"] for result in swarm_spacing.values()),
            "figures": swarm_spacing,
        },
        "assertions": assertion_items,
        "supplementary": supplementary_checks or {},
        "molecule_measurements": molecule_measurements or {},
        "standards": {
            "minimum_font_pt": MIN_FONT_PT,
            "minimum_line_pt": MIN_LINE_PT,
            "png_dpi": 600,
            "preview_dpi": 200,
        },
        "allowed_raster_counts": allowed_raster_counts or {},
    }
    output = c.QA_DIR / "F6_qa.json"
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    return report


if __name__ == "__main__":
    raise SystemExit("Run QA through build_figures.py to access Matplotlib objects.")
