"""Shared data, coordinate, layout, and plotted-value registry helpers for Figures 3-5."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Callable, Iterable, Sequence

import matplotlib.pyplot as plt
import matplotlib.text as mtext
import numpy as np
import pandas as pd
from matplotlib import font_manager, ft2font
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from style import apply_style


SOURCE_DIR = Path(__file__).resolve().parent
FIG_DIR = SOURCE_DIR.parent
ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = FIG_DIR / "data_final"
CGM_DIR = ROOT
QA_DIR = FIG_DIR / "qa"

MM_PER_INCH = 25.4
POINT_FORMATTER = {"format": "png"}

INK = "#222222"
GRAY = "#8C8C8C"
BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
LIGHT_GRAY = "#EDEDED"
MUTED = "#595959"

PLOTTED_VALUES: list[dict[str, object]] = []


def init() -> None:
    apply_style()
    QA_DIR.mkdir(parents=True, exist_ok=True)


def mm_figure(width_mm: float, height_mm: float) -> Figure:
    """Create a white canvas with an exact physical size."""
    return plt.figure(
        figsize=(width_mm / MM_PER_INCH, height_mm / MM_PER_INCH),
        dpi=600,
        facecolor="white",
    )


def cap_height_mm(
    font_size_pt: float,
    *,
    family: str = "Arial",
    weight: str = "bold",
) -> float:
    """Measure Arial Bold cap height from the H glyph and convert it to millimetres."""
    path = font_manager.findfont(
        font_manager.FontProperties(family=family, weight=weight)
    )
    font = ft2font.FT2Font(path)
    glyph = font.load_char(
        ord("H"),
        flags=ft2font.LoadFlags.NO_HINTING,
    )
    cap_height_em = glyph.height / font.units_per_EM
    return float(cap_height_em * font_size_pt * MM_PER_INCH / 72)


def add_axes_mm(
    fig: Figure,
    rect_mm: Sequence[float],
    canvas_mm: tuple[float, float],
) -> Axes:
    """Interpret millimetre coordinates with the origin at the lower-left corner."""
    x, y, width, height = map(float, rect_mm)
    canvas_width, canvas_height = canvas_mm
    ax = fig.add_axes(
        [
            x / canvas_width,
            y / canvas_height,
            width / canvas_width,
            height / canvas_height,
        ]
    )
    ax.set_facecolor("white")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK)
        ax.spines[side].set_linewidth(0.6)
    ax.tick_params(
        axis="both",
        which="major",
        direction="out",
        length=2.5,
        width=0.6,
        colors=INK,
        labelsize=6.5,
        pad=1.5,
    )
    ax.tick_params(
        axis="both",
        which="minor",
        direction="out",
        length=1.7,
        width=0.6,
        colors=INK,
        pad=1.5,
    )
    ax.xaxis.label.set_size(7)
    ax.yaxis.label.set_size(7)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)
    ax.set_axisbelow(True)
    ax._qa_points = []  # type: ignore[attr-defined]
    ax._qa_swarms = []  # type: ignore[attr-defined]
    return ax


def panel_letter(
    fig: Figure,
    letter: str,
    x_mm: float,
    y_mm: float,
    canvas_mm: tuple[float, float],
) -> mtext.Text:
    """Place panel letters from their upper-left anchors."""
    width, height = canvas_mm
    text = fig.text(
        x_mm / width,
        y_mm / height,
        letter,
        ha="left",
        va="top",
        fontsize=9,
        fontweight="bold",
        color=INK,
    )
    text._qa_panel_letter = letter  # type: ignore[attr-defined]
    return text


def mark_panel_title(text: mtext.Text) -> mtext.Text:
    text._qa_panel_title = True  # type: ignore[attr-defined]
    return text


def mark_legend_text(text: mtext.Text) -> mtext.Text:
    text._qa_legend_text = True  # type: ignore[attr-defined]
    return text


def set_bounds(ax: Axes, *, x: tuple[float, float] | None = None,
               y: tuple[float, float] | None = None) -> None:
    if x is not None:
        ax.spines["bottom"].set_bounds(*x)
    if y is not None:
        ax.spines["left"].set_bounds(*y)


def load_csv(name: str, *, root: Path | None = None) -> pd.DataFrame:
    """Read every source CSV with round-trip floating-point precision."""
    path = (root or DATA_DIR) / name
    return pd.read_csv(path, float_precision="round_trip")


def record_value(
    fig_no: int,
    panel: str,
    series: str,
    row_id: str,
    variable: str,
    value: float,
    source_file: str,
) -> None:
    """Record values that enter plotted objects; never disguise a missing value as zero."""
    numeric = float(value)
    if not np.isfinite(numeric):
        raise ValueError(f"Plotted value is not finite: {panel} {series} {row_id} {variable}")
    PLOTTED_VALUES.append(
        {
            "figure": fig_no,
            "panel": panel,
            "series": series,
            "row_identifier": str(row_id),
            "variable": variable,
            "source_file": source_file,
            "value": numeric,
        }
    )


def write_plotted_values(fig_no: int) -> None:
    path = QA_DIR / f"Figure_{fig_no}_plotted_values.csv"
    rows = [row for row in PLOTTED_VALUES if row["figure"] == fig_no]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=[
                "figure",
                "panel",
                "series",
                "row_identifier",
                "variable",
                "source_file",
                "value",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def register_axis_points(ax: Axes, x: Iterable[float], y: Iterable[float]) -> None:
    """Expose points and CI endpoints so QA can check them against the displayed axes."""
    for xv, yv in zip(x, y, strict=True):
        if np.isfinite(xv) and np.isfinite(yv):
            ax._qa_points.append((float(xv), float(yv)))  # type: ignore[attr-defined]


def register_swarm_points(
    ax: Axes,
    row_key: str,
    x: Iterable[float],
    y: Iterable[float],
    marker_diameter_pt: float,
) -> None:
    """Register points on a swarm row so QA can check centre distances in millimetres."""
    points = [
        (float(xv), float(yv))
        for xv, yv in zip(x, y, strict=True)
        if np.isfinite(xv) and np.isfinite(yv)
    ]
    ax._qa_swarms.append(  # type: ignore[attr-defined]
        {
            "row_key": row_key,
            "marker_diameter_pt": float(marker_diameter_pt),
            "points": points,
        }
    )


def shade(color: str, fraction: float = 0.35) -> str:
    """Return an opaque hexadecimal colour made by blending a layout colour with white."""
    color = color.lstrip("#")
    rgb = np.array([int(color[i:i + 2], 16) for i in (0, 2, 4)], dtype=float)
    mixed = fraction * rgb + (1 - fraction) * 255
    return "#" + "".join(f"{round(v):02X}" for v in mixed)


def swarm_offsets(
    values: Sequence[float],
    *,
    x_to_mm: Callable[[float], float],
    y_center: float,
    max_half_height_data: float,
    data_to_y_mm: float,
    marker_diameter_pt: float = 2.8,
) -> list[float]:
    """
    Deterministic horizontal swarm layout.

    Values are unchanged. Candidate offsets start at zero and at tangent positions
    of placed circles, then the smallest absolute offset is selected. Distances are
    checked in millimetres and must be at least 0.95 symbol diameters.
    """
    diameter_mm = marker_diameter_pt * 25.4 / 72
    minimum_mm = 0.95 * diameter_mm
    indexed = sorted(enumerate(values), key=lambda item: (float(item[1]), item[0]))
    result = [0.0] * len(values)
    placed: list[tuple[float, float]] = []
    center_mm = y_center * data_to_y_mm
    for index, value in indexed:
        x_mm = x_to_mm(float(value))
        candidates = [center_mm]
        for old_x_mm, old_y_mm in placed:
            dx = abs(x_mm - old_x_mm)
            if dx < minimum_mm:
                rise = float(np.sqrt(max(minimum_mm**2 - dx**2, 0.0)))
                candidates.extend([old_y_mm + rise, old_y_mm - rise])

        valid = []
        for y_mm in candidates:
            if abs(y_mm - center_mm) > max_half_height_data * data_to_y_mm + 1e-9:
                continue
            if all(
                np.hypot(x_mm - old_x_mm, y_mm - old_y_mm) >= minimum_mm - 1e-9
                for old_x_mm, old_y_mm in placed
            ):
                valid.append(float(y_mm))
        if not valid:
            raise RuntimeError(
                f"swarm layout within {max_half_height_data=} cannot avoid overlap "
                f"at point {index + 1}, value={value}"
            )
        chosen_mm = min(
            valid,
            key=lambda y_mm: (abs(y_mm - center_mm), y_mm < center_mm),
        )
        result[index] = chosen_mm / data_to_y_mm - y_center
        placed.append((x_mm, chosen_mm))
    return result


def ci_handle(
    *,
    marker: str = "o",
    marker_size: float = 3.0,
    color: str = "black",
    hollow: bool = False,
    line: bool = True,
) -> Line2D:
    return Line2D(
        [0],
        [0],
        color=color,
        linewidth=0.7 if line else 0,
        linestyle="-" if line else "none",
        marker=marker,
        markersize=marker_size,
        markerfacecolor="none" if hollow else color,
        markeredgecolor=color,
        markeredgewidth=0.7,
        solid_capstyle="butt",
    )


def add_break_marks(ax: Axes, side: str) -> None:
    """Draw parallel break marks about 1.2 mm long at 60 degrees along the broken axis."""
    axis_width_mm = ax.get_position().width * ax.figure.get_figwidth() * 25.4
    axis_height_mm = ax.get_position().height * ax.figure.get_figheight() * 25.4
    dx = 1.2 * np.cos(np.deg2rad(60)) / axis_width_mm
    dy = 1.2 * np.sin(np.deg2rad(60)) / axis_height_mm
    x = 0 if side == "left" else 1
    ax.plot(
        (x - dx / 2, x + dx / 2),
        (-dy / 2, dy / 2),
        transform=ax.transAxes,
        color=INK,
        lw=0.6,
        clip_on=False,
        solid_capstyle="butt",
    )


def legend_patch(
    color: str,
    edge: str = "none",
    linewidth: float = 0,
    alpha: float = 1.0,
) -> Patch:
    return Patch(
        facecolor=color,
        edgecolor=edge,
        linewidth=linewidth,
        alpha=alpha,
    )
