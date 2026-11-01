"""图 3–5 共用的数据、坐标、排版与绘图值登记工具。"""

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
    """创建物理尺寸严格的白色画布。"""
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
    """按 Arial Bold 的 H 字形读取 cap height，并换算为毫米。"""
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
    """按左下角原点解释毫米坐标。"""
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
    """面板字母按左上角定位。"""
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
    """所有源 CSV 必须以 round-trip 浮点精度读取。"""
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
    """记录真正进入绘图对象的数值；缺失值不伪装成零。"""
    numeric = float(value)
    if not np.isfinite(numeric):
        raise ValueError(f"未绘图数值不是有限值：{panel} {series} {row_id} {variable}")
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
    """供 QA 检查点或 CI 端点是否在坐标轴显示范围内。"""
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
    """登记同一蜂群行，供 QA 在毫米空间检查圆心距。"""
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
    """返回布局色与白色按指定比例混合后的不透明十六进制色。"""
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
    确定性横向蜂群。

    数值本身不改变。从零及已放置圆的相切位置产生连续候选偏移，
    按绝对值最小选择；已在毫米空间中的距离须至少为 0.95 个符号直径。
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
                f"蜂群排布在 {max_half_height_data=} 内无法避免重叠，"
                f"第 {index + 1} 个点 value={value}"
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
    """在断开轴底边各画一条同向“/”，长约 1.2 mm、倾角约 60°。"""
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


def legend_patch(color: str, edge: str = "none", linewidth: float = 0) -> Patch:
    return Patch(facecolor=color, edgecolor=edge, linewidth=linewidth)
