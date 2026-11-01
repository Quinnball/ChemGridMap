"""Rebuild Figure 3: held-out evaluation across 52 datasets."""

from __future__ import annotations

import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

import common as c


W, H = 170.0, 160.0


def _row_axis(ax, frame: pd.DataFrame, label_column: str) -> None:
    ranks = frame[label_column].astype(int).to_numpy()
    count = len(ranks)
    padding = (ax.get_position().height * ax.figure.get_figheight() * 25.4 / 2.54 - (count - 1)) / 2
    ax.set_ylim(count + padding, 1 - padding)
    ax.set_yticks(ranks)
    ax.set_yticklabels(
        [
            (
                f"{symbol}† ({target})"
                if target in {"CHEMBL205", "CHEMBL204"}
                else f"{symbol} ({target})"
            )
            for symbol, target in zip(
                frame["symbol"],
                frame["target_id"],
                strict=True,
            )
        ],
        fontsize=6,
    )
    ax.tick_params(axis="y", length=0, pad=2)
    ax.spines["left"].set_visible(False)


def _draw_forest_panel(
    fig: Figure,
    frame: pd.DataFrame,
    rect: tuple[float, float, float, float],
    endpoint: str,
) -> None:
    ax = c.add_axes_mm(fig, rect, (W, H))
    _row_axis(ax, frame, "target_rank")
    ax.set_xlim(-0.02, 1.02)
    ax.set_xticks(np.arange(0, 1.01, 0.2))
    c.set_bounds(ax, x=(0, 1))

    for _, row in frame.iterrows():
        y = float(row["target_rank"])
        ax.axhline(y, xmin=0, xmax=1, color=c.LIGHT_GRAY, lw=0.4, zorder=0)
    ax.axvline(0.5, color=c.GRAY, lw=0.6, ls=(0, (3, 2)), zorder=1)

    for _, row in frame.iterrows():
        y = float(row["target_rank"])
        ax.hlines(
            y,
            float(row["g_ci_low"]),
            float(row["g_ci_high"]),
            color="black",
            lw=0.7,
            zorder=2,
        )
        points = [
            float(row["g_ci_low"]),
            float(row["g_ci_high"]),
            float(row["g_auc"]),
            float(row["d_auc"]),
        ]
        c.register_axis_points(ax, points, [y] * len(points))
        c.record_value(3, "3b", "G 95% CI", row["target_id"], "g_ci_low", row["g_ci_low"], f"3b_{endpoint}.csv")
        c.record_value(3, "3b", "G 95% CI", row["target_id"], "g_ci_high", row["g_ci_high"], f"3b_{endpoint}.csv")

    ax.plot(
        frame["g_auc"],
        frame["target_rank"],
        linestyle="none",
        marker="o",
        markersize=3.0,
        markerfacecolor="black",
        markeredgecolor="black",
        markeredgewidth=0.5,
        zorder=4,
    )
    ax.plot(
        frame["d_auc"],
        frame["target_rank"],
        linestyle="none",
        marker="o",
        markersize=3.4,
        markerfacecolor="none",
        markeredgecolor=c.BLUE,
        markeredgewidth=0.7,
        zorder=5,
    )

    for _, row in frame.iterrows():
        row_id = row["target_id"]
        for variable, series in [
            ("g_auc", "G"),
            ("d_auc", "D"),
        ]:
            if np.isfinite(row[variable]):
                c.record_value(3, "3b", series, row_id, variable, row[variable], f"3b_{endpoint}.csv")
    return ax


def build() -> Figure:
    c.init()
    fig = c.mm_figure(W, H)
    c.panel_letter(fig, "a", 1, 159, (W, H))
    c.panel_letter(fig, "b", 96, 159, (W, H))
    c.panel_letter(fig, "c", 1, 114, (W, H))

    # 3a: B0 disagreement rate
    frame_3a = c.load_csv("3a_B0_disagreement.csv")
    ax_3a = c.add_axes_mm(fig, (14, 126.5, 76, 24), (W, H))
    ax_3a.set_xlim(-0.01, 0.71)
    ax_3a.set_xticks(np.arange(0, 0.71, 0.1))
    ax_3a.set_ylim(-0.5, 1.5)
    ax_3a.set_yticks([1, 0])
    ax_3a.set_yticklabels(["IC50", "Ki"], fontsize=6.5)
    ax_3a.tick_params(axis="y", length=0, pad=2)
    ax_3a.spines["left"].set_visible(False)
    ax_3a.set_xlabel("B0 disagreement rate")
    c.set_bounds(ax_3a, x=(0, 0.7))

    scale_x = fig.bbox.width / fig.dpi * c.MM_PER_INCH / W
    scale_y = fig.bbox.height / fig.dpi * c.MM_PER_INCH / H
    data_to_y_mm = (24 / 2) * scale_y
    x_to_mm = lambda value: (float(value) + 0.01) / 0.72 * 76 * scale_x
    for center, column, endpoint in ((1.12, "IC50", "IC50"), (0.12, "Ki", "Ki")):
        values = frame_3a[column].dropna()
        offsets = c.swarm_offsets(
            values.to_numpy(),
            x_to_mm=x_to_mm,
            y_center=center,
            max_half_height_data=0.30,
            data_to_y_mm=data_to_y_mm,
            marker_diameter_pt=2.4,
        )
        ys = center + np.asarray(offsets)
        ax_3a.plot(
            values.to_numpy(),
            ys,
            linestyle="none",
            marker="o",
            markersize=2.4,
            markerfacecolor="#6E6E6E",
            markeredgecolor="white",
            markeredgewidth=0.4,
            zorder=2,
        )
        c.register_axis_points(ax_3a, values.to_numpy(), ys)
        c.register_swarm_points(ax_3a, f"3a:{endpoint}", values.to_numpy(), ys, 2.4)
        for row_id, value in zip(frame_3a.loc[values.index, "target_id"], values, strict=True):
            c.record_value(3, "3a", endpoint, row_id, column, value, "3a_B0_disagreement.csv")

        median = float(values.median())
        q1 = float(values.quantile(0.25))
        q3 = float(values.quantile(0.75))
        summary_y = center - 0.46
        ax_3a.hlines(summary_y, q1, q3, color="black", lw=0.9, linestyle="none", zorder=3)
        ax_3a.plot(median, summary_y, marker="o", ms=3.2, color="black", ls="none", zorder=4)
        c.register_axis_points(ax_3a, [q1, median, q3], [summary_y] * 3)
        c.record_value(3, "3a", f"{endpoint} median", "summary", "median", median, "3a_B0_disagreement.csv")
        c.record_value(3, "3a", f"{endpoint} quartiles", "summary", "q1", q1, "3a_B0_disagreement.csv")
        c.record_value(3, "3a", f"{endpoint} quartiles", "summary", "q3", q3, "3a_B0_disagreement.csv")

    # 3b: per-dataset AUC
    frame_ic50 = c.load_csv("3b_IC50.csv")
    frame_ki = c.load_csv("3b_Ki.csv")
    # The long 5b labels need about 25 mm on the left; narrow and shift the plot,
    # keeping clear separation from the 5c plotting area.
    ax_ic50 = _draw_forest_panel(fig, frame_ic50, (120.0, 77.75, 48.5, 72.85), "IC50")
    ax_ic50._qa_fig3b_ic50 = True  # type: ignore[attr-defined]
    ax_ki = _draw_forest_panel(fig, frame_ki, (120.0, 14.3, 48.5, 57.8), "Ki")
    ax_ki.set_xlabel("AUC")

    group_style = dict(ha="left", va="bottom", fontsize=6.5, color=c.INK)
    c.mark_panel_title(fig.text(94.5 / W, 151.3 / H, "IC50 (29 datasets)", **group_style))
    ki_title = c.mark_panel_title(
        fig.text(94.4 / W, 72.7 / H, "Ki (23 datasets)", **group_style)
    )
    ki_title._qa_fig3b_ki_title = True  # type: ignore[attr-defined]
    fig.legend(
        handles=[
            c.ci_handle(marker="o", marker_size=3.0, color="black", line=True),
            c.ci_handle(marker="o", marker_size=3.4, color=c.BLUE, hollow=True, line=False),
        ],
        labels=["G with 95% CI", "D (post hoc)"],
        loc="center left",
        bbox_to_anchor=(100 / W, 157.3 / H),
        frameon=False,
        fontsize=6.5,
        ncol=2,
        handlelength=1.8,
        handletextpad=0.45,
        columnspacing=1.1,
        borderaxespad=0,
    )

    # 3c: alternative-model differences with a broken x-axis
    row_labels = ["G-record", "G-global", "G-fixed", "X"]
    columns = [
        "g_minus_record_log_loss",
        "g_minus_global_log_loss",
        "g_minus_fixed_log_loss",
        "g_minus_x_log_loss",
    ]
    for endpoint, y0 in (("IC50", 68.0), ("Ki", 17.0)):
        frame = c.load_csv(f"3c_{endpoint}.csv")
        target = frame[frame["row_type"] == "target"]
        summary = frame[frame["row_type"] == "summary"].iloc[0]
        all_values = target[columns].to_numpy(dtype=float)
        finite = all_values[np.isfinite(all_values)]
        assert finite.min() >= -0.70 and finite.max() <= 0.20
        assert not np.any((finite > -0.40) & (finite < -0.32))

        ext = c.add_axes_mm(fig, (14, y0, 15, 40), (W, H))
        main = c.add_axes_mm(fig, (31, y0, 59, 40), (W, H))
        for ax in (main, ext):
            ax.set_ylim(-0.5, 3.5)
            ax.set_yticks([0, 1, 2, 3])
            ax.tick_params(axis="y", length=0)
            ax.spines["left"].set_visible(False)
        ext.set_yticklabels(row_labels[::-1], fontsize=6.5)
        ext.set_xlim(-0.70, -0.40)
        ext.set_xticks([-0.6, -0.4])
        c.set_bounds(ext, x=(-0.7, -0.4))
        main.set_xlim(-0.32, 0.20)
        main.set_xticks([-0.2, -0.1, 0.0, 0.1, 0.2])
        c.set_bounds(main, x=(-0.3, 0.2))
        ext.spines["right"].set_visible(False)
        main.spines["left"].set_visible(False)
        main.set_yticklabels([])
        c.add_break_marks(ext, "right")
        c.add_break_marks(main, "left")
        main.axvline(0, color=c.GRAY, lw=0.6, ls=(0, (3, 2)), zorder=0)

        x_to_panel_mm = lambda value: scale_x * (
            (float(value) + 0.70) / 0.30 * 15.0
            if float(value) <= -0.40
            else 17.0 + (float(value) + 0.32) / 0.52 * 59.0
        )
        for row_index, (label, column) in enumerate(zip(row_labels, columns, strict=True)):
            values = target[column].dropna()
            y_center = 3.0 - row_index + 0.12
            marker_size = 2.4
            offsets = c.swarm_offsets(
                values.to_numpy(),
                x_to_mm=x_to_panel_mm,
                y_center=y_center,
                max_half_height_data=0.34,
                data_to_y_mm=(40 / 4) * scale_y,
                marker_diameter_pt=marker_size,
            )
            ys = y_center + np.asarray(offsets)
            ext_values = values.to_numpy() <= -0.40
            for ax, mask in ((ext, ext_values), (main, ~ext_values)):
                if not np.any(mask):
                    continue
                ax.plot(
                    values.to_numpy()[mask],
                    ys[mask],
                    linestyle="none",
                    marker="o",
                    markersize=marker_size,
                    markerfacecolor="#6E6E6E",
                    markeredgecolor="white",
                    markeredgewidth=0.4,
                    zorder=2,
                )
                c.register_axis_points(ax, values.to_numpy()[mask], ys[mask])
                c.register_swarm_points(
                    ax,
                    f"3c:{endpoint}:{label}",
                    values.to_numpy()[mask],
                    ys[mask],
                    marker_size,
                )
            for row_id, value in zip(target.loc[values.index, "target_id"], values, strict=True):
                c.record_value(3, "3c", label, row_id, column, value, f"3c_{endpoint}.csv")

            median = float(summary[f"median_{column}"])
            low = float(summary[f"ci_low_{column}"])
            high = float(summary[f"ci_high_{column}"])
            summary_y = 3.0 - row_index - 0.36
            summary_ax = ext if median <= -0.40 else main
            summary_ax.hlines(summary_y, low, high, color="black", lw=0.9, zorder=3)
            summary_ax.plot(median, summary_y, marker="o", color="black", ms=3.2, ls="none", zorder=4)
            c.register_axis_points(summary_ax, [low, median, high], [summary_y] * 3)
            c.record_value(3, "3c", f"{label} summary", "summary", f"median_{column}", median, f"3c_{endpoint}.csv")
            c.record_value(3, "3c", f"{label} summary", "summary", f"ci_low_{column}", low, f"3c_{endpoint}.csv")
            c.record_value(3, "3c", f"{label} summary", "summary", f"ci_high_{column}", high, f"3c_{endpoint}.csv")

        fig.text(14 / W, (y0 + 40.8) / H, endpoint, ha="left", va="bottom", fontsize=6.5, color=c.INK)
        fig.text(90 / W, (y0 + 40.8) / H, "← G better", ha="right", va="bottom", fontsize=6, color=c.MUTED)

    fig.text(
        52 / W,
        11.0 / H,
        "Δ log loss (G − comparator)",
        ha="center",
        va="top",
        fontsize=7,
        color=c.INK,
    )

    c.write_plotted_values(3)
    return fig
