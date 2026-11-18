"""Figure 4 (manuscript v0.7): what colours imply for other documents.

The hERG patent pair (former panel b) moved to Additional file 1; the status series use
black, grey and hatching, hues that the map does not use (review items S2, S4, S20)."""

from __future__ import annotations

import numpy as np
from matplotlib import patheffects
from matplotlib.colors import to_hex
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from chemgridmap.plotting import DEFAULT_COLORS

import common as c


W, H = 170.0, 130.0
SUPPORTED = "#1a1a1a"
UNCERTAIN = "#9a9a9a"
TOUCHED_EDGE = "#4d4d4d"
BAR_ALPHA = 0.72
CLASS_LABELS = {"inactive": "Inactive (≤ 6)", "medium": "Medium (6–7)", "active": "Active (≥ 7)"}
# Mantel–Haenszel supported-minus-uncertain differences standardised by class (time split 2015;
# rev7/analysis/s20_mh.json): estimate, dataset-bootstrap 95% interval, datasets, strata.
MH_2015 = {"IC50": (0.187, 0.088, 0.276, 19, 34), "Ki": (0.172, 0.093, 0.245, 19, 39)}


def _target_label(target: str, symbol: str, *, dagger: bool) -> str:
    marker = "†" if dagger else ""
    return f"{symbol}{marker} ({target})"


def _symbols() -> dict[str, str]:
    frame = c.load_csv("target_symbols.csv")
    return dict(zip(frame["target_id"], frame["symbol"], strict=True))


def _calibration_panel(
    fig: Figure,
    frame: object,
    rect: tuple[float, float, float, float],
    title: str,
    *,
    show_y: bool,
    show_y_label: bool,
    show_legend: bool,
) -> None:
    ax = c.add_axes_mm(fig, rect, (W, H))
    ax.set_label(f"4a:{title}")
    ax._qa_fig4a = True  # type: ignore[attr-defined]
    ax._qa_fig4a_endpoint = title  # type: ignore[attr-defined]
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.03, 1.03)
    ticks = np.arange(0.0, 1.01, 0.2)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.tick_params(axis="x", labelbottom=True)
    ax.tick_params(axis="y", labelleft=show_y)
    c.set_bounds(ax, x=(0, 1), y=(0, 1))
    ax.plot([0, 1], [0, 1], color=c.GRAY, lw=0.6, ls=(0, (3, 2)), zorder=0)
    if show_y_label:
        ax.set_ylabel("Observed agreement")

    visible = frame.loc[~frame["omit_in_figure"].astype(bool)].copy()
    plotted_records = []
    for display_class, marker in (
        ("inactive", "o"),
        ("medium", "s"),
        ("active", "^"),
    ):
        local = visible.loc[
            visible["display_class"].eq(display_class)
        ].sort_values("mean_p")
        if local.empty:
            continue
        color = DEFAULT_COLORS[display_class]
        ax.plot(
            local["mean_p"],
            local["agreement_rate"],
            color=color,
            lw=0.8,
            linestyle="-",
            marker=marker,
            markersize=3.2,
            markerfacecolor=color,
            markeredgecolor=color,
            markeredgewidth=0.5,
            label=CLASS_LABELS[display_class],
            zorder=3,
        )
        ax.vlines(
            local["mean_p"],
            local["wilson_ci_low"],
            local["wilson_ci_high"],
            color=color,
            lw=0.7,
            zorder=2,
        )
        c.register_axis_points(
            ax,
            np.r_[local["mean_p"], local["mean_p"], local["mean_p"]],
            np.r_[local["wilson_ci_low"], local["agreement_rate"], local["wilson_ci_high"]],
        )
        for index, row in enumerate(local.itertuples(index=False)):
            plotted_records.append(
                {
                    "endpoint": title,
                    "display_class": str(row.display_class),
                    "probability_bin": int(row.probability_bin),
                    "mean_p": float(row.mean_p),
                    "agreement_rate": float(row.agreement_rate),
                    "wilson_ci_low": float(row.wilson_ci_low),
                    "wilson_ci_high": float(row.wilson_ci_high),
                    "color": to_hex(color, keep_alpha=False).lower(),
                    "marker": marker,
                }
            )
            c.record_value(
                4,
                "4a",
                display_class.capitalize(),
                f"{title} bin {int(row.probability_bin)}",
                "mean_p",
                row.mean_p,
                "4a_calibration_by_class.csv",
            )
            c.record_value(
                4,
                "4a",
                display_class.capitalize(),
                f"{title} bin {int(row.probability_bin)}",
                "agreement_rate",
                row.agreement_rate,
                "4a_calibration_by_class.csv",
            )
            c.record_value(
                4,
                "4a",
                f"{display_class.capitalize()} Wilson 95% CI",
                f"{title} bin {int(row.probability_bin)}",
                "wilson_ci_low",
                row.wilson_ci_low,
                "4a_calibration_by_class.csv",
            )
            c.record_value(
                4,
                "4a",
                f"{display_class.capitalize()} Wilson 95% CI",
                f"{title} bin {int(row.probability_bin)}",
                "wilson_ci_high",
                row.wilson_ci_high,
                "4a_calibration_by_class.csv",
            )
    ax._qa_fig4a_records = plotted_records  # type: ignore[attr-defined]
    if show_legend:
        ax.legend(
            loc="upper left",
            bbox_to_anchor=(0.015, 0.99),
            frameon=False,
            fontsize=6,
            handlelength=1.3,
            handletextpad=0.35,
            borderaxespad=0,
            labelspacing=0.25,
        )
    c.mark_panel_title(fig.text(
        rect[0] / W,
        (rect[1] + rect[3] + 0.8) / H,
        title,
        ha="left",
        va="bottom",
        fontsize=6.5,
        color=c.INK,
    ))


def _row_panel(
    fig: Figure,
    filename: str,
    rect: tuple[float, float, float, float],
) -> None:
    frame = c.load_csv(filename)
    symbols = _symbols()
    ax = c.add_axes_mm(fig, rect, (W, H))
    count = len(frame)
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(count + 0.5, 0.5)
    ax.set_xticks([0.0, 0.5, 1.0])
    ax.set_xticks([0.25, 0.75], minor=True)
    ax.set_yticks(frame["dataset_rank"].astype(int))
    ax.set_yticklabels(
        [
            _target_label(
                target,
                symbols[target],
                dagger=target in {"CHEMBL240", "CHEMBL205"},
            )
            for target in frame["target_id"]
        ],
        fontsize=6,
    )
    ax.tick_params(axis="y", length=0, pad=2)
    ax.spines["left"].set_visible(False)
    c.set_bounds(ax, x=(0, 1))

    for _, row in frame.iterrows():
        rank = float(row["dataset_rank"])
        ax.axhline(rank, color=c.LIGHT_GRAY, lw=0.4, zorder=0)
        for status, dy, color in (
            ("supported", -0.22, SUPPORTED),
            ("uncertain", 0.22, UNCERTAIN),
        ):
            y = rank + dy
            low = float(row[f"{status}_ci_low"])
            high = float(row[f"{status}_ci_high"])
            value = float(row[f"{status}_agreement"])
            ax.hlines(y, low, high, color=color, lw=0.7, zorder=2)
            ax.plot(
                value,
                y,
                marker="o",
                markersize=3.0,
                markerfacecolor=color if status == "supported" else "white",
                markeredgecolor=color,
                markeredgewidth=0.5,
                linestyle="none",
                zorder=3,
            )
            c.register_axis_points(ax, [low, value, high], [y] * 3)
            c.record_value(
                4,
                "4c",
                status.capitalize(),
                row["target_id"],
                f"{status}_agreement",
                value,
                filename,
            )
            c.record_value(
                4,
                "4c",
                f"{status.capitalize()} 95% Wilson CI",
                row["target_id"],
                f"{status}_ci_low",
                low,
                filename,
            )
            c.record_value(
                4,
                "4c",
                f"{status.capitalize()} 95% Wilson CI",
                row["target_id"],
                f"{status}_ci_high",
                high,
                filename,
            )
    return ax


def build() -> Figure:
    c.init()
    fig = c.mm_figure(W, H)
    c.panel_letter(fig, "a", 1, 129, (W, H))
    c.panel_letter(fig, "b", 1, 73.5, (W, H))
    c.panel_letter(fig, "c", 101, 68, (W, H))

    calibration = c.load_csv("4a_calibration_by_class.csv")
    _calibration_panel(
        fig,
        calibration[calibration["endpoint"].eq("IC50")],
        (30, 80, 42, 42),
        "IC50",
        show_y=True,
        show_y_label=True,
        show_legend=True,
    )
    _calibration_panel(
        fig,
        calibration[calibration["endpoint"].eq("Ki")],
        (84, 80, 42, 42),
        "Ki",
        show_y=False,
        show_y_label=False,
        show_legend=False,
    )
    c.mark_panel_title(fig.text(
        78 / W,
        73.8 / H,
        "Predicted probability",
        ha="center",
        va="top",
        fontsize=7,
        color=c.INK,
    ))

    # 4c: 2015 time split
    # Gene symbols plus ChEMBL IDs need a wider left label area.
    _row_panel(fig, "4c_IC50.csv", (27.0, 15, 20, 42.5))
    ax_ki = _row_panel(fig, "4c_Ki.csv", (72.5, 15, 22, 42.5))
    ax_ki.set_xlabel("")
    c.mark_panel_title(fig.text(2 / W, 61.8 / H, "IC50 (18 datasets)", ha="left", va="bottom", fontsize=6.5, color=c.INK))
    for endpoint, x_text in (("IC50", 2.0), ("Ki", 48.0)):
        est, low, high, n_sets, n_strata = MH_2015[endpoint]
        c.mark_panel_title(fig.text(x_text / W, 58.6 / H, f"MH difference {est:.2f} ({low:.2f}–{high:.2f})",
                                    ha="left", va="bottom", fontsize=5.5, color=c.INK))
    c.mark_panel_title(fig.text(48 / W, 61.8 / H, "Ki (15 datasets)", ha="left", va="bottom", fontsize=6.5, color=c.INK))
    c.mark_panel_title(fig.text(
        56.5 / W,
        10.7 / H,
        "Agreement with later documents",
        ha="center",
        va="top",
        fontsize=7,
        color=c.INK,
    ))

    # 4d: three rule-P states grouped by past-document count.
    status_style = (
        ("P_supported", -0.32, SUPPORTED),
        ("uncertain_low_support", 0.0, UNCERTAIN),
        ("touched_all", 0.32, "white"),
    )
    for endpoint, x0, show_y in (("IC50", 108.0, True), ("Ki", 139.0, False)):
        frame = c.load_csv(f"4d_{endpoint}.csv")
        axd = c.add_axes_mm(fig, (x0, 15, 29, 43), (W, H))
        axd.set_xlim(-0.55, 1.55)
        axd.set_ylim(0, 1.75)
        axd.set_xticks([0, 1])
        axd.set_xticklabels(["1", "≥2"])
        axd.set_yticks(np.arange(0, 1.01, 0.2))
        axd.tick_params(axis="y", labelleft=show_y)
        c.set_bounds(axd, y=(0, 1))
        if show_y:
            axd.set_ylabel("Agreement")
        for _, row in frame.iterrows():
            x = 0.0 if row["earlier_units_group"] == "m_eq_1" else 1.0
            for state, dx, color in status_style:
                if row["state"] != state:
                    continue
                n_predictions = int(row["n_predictions"])
                if n_predictions < 8:
                    continue
                x_bar = x + dx
                height = float(row["agreement_rate"])
                low = float(row["agreement_ci_low"])
                high = float(row["agreement_ci_high"])
                touched = state == "touched_all"
                axd.bar(
                    x_bar,
                    height,
                    width=0.25,
                    color=color,
                    edgecolor=TOUCHED_EDGE if touched else "none",
                    hatch="//////" if touched else None,
                    linewidth=0.5 if touched else 0,
                    alpha=1.0 if touched else BAR_ALPHA,
                    zorder=1,
                )
                axd.patches[-1]._qa_fig4d_state = state  # type: ignore[attr-defined]
                axd.patches[-1]._qa_fig4d_color = color  # type: ignore[attr-defined]
                interval_color = TOUCHED_EDGE if touched else color
                interval = axd.vlines(
                    x_bar,
                    low,
                    high,
                    color=interval_color,
                    lw=0.8,
                    zorder=4,
                )
                caps = axd.hlines(
                    [low, high],
                    x_bar - 0.035,
                    x_bar + 0.035,
                    color=interval_color,
                    lw=0.8,
                    zorder=4,
                )
                interval._qa_fig4d_interval = True  # type: ignore[attr-defined]
                caps._qa_fig4d_interval = True  # type: ignore[attr-defined]
                mean = float(row["mean_probability"])
                axd.plot(
                    [x_bar - 0.15, x_bar + 0.15],
                    [mean, mean],
                    color="black",
                    lw=1.1,
                    solid_capstyle="butt",
                    path_effects=[patheffects.Stroke(linewidth=2.4, foreground="white"), patheffects.Normal()],
                    zorder=3,
                )
                c.register_axis_points(axd, [x_bar, x_bar], [low, high])
                c.register_axis_points(axd, [x_bar - 0.15, x_bar + 0.15], [mean, mean])
                series = {
                    "P_supported": "Supported",
                    "uncertain_low_support": "Uncertain low support",
                    "touched_all": "Touched by an offset pair",
                }[state]
                row_id = f"{endpoint} {row['earlier_units_group']}"
                c.record_value(
                    4,
                    "4d",
                    series,
                    row_id,
                    "agreement_rate",
                    height,
                    f"4d_{endpoint}.csv",
                )
                c.record_value(
                    4,
                    "4d",
                    f"{series} 95% CI",
                    row_id,
                    "agreement_ci_low",
                    low,
                    f"4d_{endpoint}.csv",
                )
                c.record_value(
                    4,
                    "4d",
                    f"{series} 95% CI",
                    row_id,
                    "agreement_ci_high",
                    high,
                    f"4d_{endpoint}.csv",
                )
                c.record_value(
                    4,
                    "4d",
                    f"{series} mean predicted probability",
                    row_id,
                    "mean_probability",
                    mean,
                    f"4d_{endpoint}.csv",
                )
                c.record_value(
                    4,
                    "4d",
                    f"{series} n",
                    row_id,
                    "n_predictions",
                    n_predictions,
                    f"4d_{endpoint}.csv",
                )
                n_datasets = int(row["n_datasets"])
                c.record_value(
                    4,
                    "4d",
                    f"{series} datasets",
                    row_id,
                    "n_datasets",
                    n_datasets,
                    f"4d_{endpoint}.csv",
                )
                count_label = axd.text(
                    x_bar,
                    1.02,
                    f"n = {n_predictions} · {n_datasets} sets",
                    ha="left",
                    va="center",
                    fontsize=6,
                    color=c.INK,
                    zorder=6,
                    rotation=90,
                    rotation_mode="anchor",
                )
                count_label._qa_fig4d_count = (  # type: ignore[attr-defined]
                    endpoint,
                    row["earlier_units_group"],
                    state,
                    n_predictions,
                    n_datasets,
                )
        axd._qa_fig4d = True  # type: ignore[attr-defined]
        axd._qa_fig4d_endpoint = endpoint  # type: ignore[attr-defined]
        axd._qa_fig4d_rows = frame[  # type: ignore[attr-defined]
            ["earlier_units_group", "state", "n_predictions", "n_datasets"]
        ].to_dict("records")
        c.mark_panel_title(fig.text(x0 / W, 58.8 / H, endpoint, ha="left", va="bottom", fontsize=6.5, color=c.INK))

    c.mark_panel_title(fig.text(
        139.5 / W,
        10.7 / H,
        "Earlier documents",
        ha="center",
        va="top",
        fontsize=7,
        color=c.INK,
    ))

    fig.legend(
        handles=[
            Line2D([0], [0], marker="o", linestyle="none", markersize=3.0, color=SUPPORTED),
            Line2D([0], [0], marker="o", linestyle="none", markersize=3.0, markerfacecolor="white",
                   markeredgecolor=UNCERTAIN, color=UNCERTAIN),
        ],
        labels=["Supported", "Uncertain"],
        loc="center left",
        bbox_to_anchor=(6 / W, 68.5 / H),
        frameon=False,
        fontsize=6.5,
        ncol=2,
        handlelength=1.2,
        handletextpad=0.4,
        columnspacing=1.2,
        borderaxespad=0,
    )
    legend_patches = [
        c.legend_patch(SUPPORTED, alpha=BAR_ALPHA),
        c.legend_patch(UNCERTAIN, alpha=BAR_ALPHA),
        c.legend_patch("white", edge=TOUCHED_EDGE, linewidth=0.5, alpha=1.0),
    ]
    legend_patches[2].set_hatch("//////")
    legend = fig.legend(
        handles=[
            *legend_patches,
            Line2D([0], [0], color="black", lw=1.1),
        ],
        labels=[
            "Supported",
            "Uncertain: low support",
            "Touched by an offset pair",
            "Mean predicted probability",
        ],
        loc="upper left",
        bbox_to_anchor=(112.5 / W, 68.7 / H),
        frameon=False,
        fontsize=6,
        ncol=2,
        handlelength=1.0,
        handletextpad=0.35,
        columnspacing=1.0,
        labelspacing=0.45,
        borderaxespad=0,
    )
    for patch, state in zip(
        legend.legend_handles[:3],
        ("P_supported", "uncertain_low_support", "touched_all"),
        strict=True,
    ):
        patch._qa_fig4d_legend_state = state  # type: ignore[attr-defined]

    c.write_plotted_values(4)
    return fig
