"""Additional file 1 figure: the hERG patent pair formerly shown as Fig. 4b (v0.6)."""

from __future__ import annotations

from matplotlib.figure import Figure

import common as c

W, H = 80.0, 66.0


def build() -> Figure:
    c.init()
    fig = c.mm_figure(W, H)
    # 4b: hERG patent pair
    pair = c.load_csv("4b_hERG_patent_pair.csv")
    ax = c.add_axes_mm(fig, (14, 12, 48, 48), (W, H))
    ax.set_xlim(4, 10)
    ax.set_ylim(4, 10)
    ax.set_xticks([4, 6, 8, 10])
    ax.set_yticks([4, 6, 8, 10])
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("pChEMBL, US-10131663-B2")
    ax.set_ylabel("pChEMBL, US-10590130-B2")
    ax.plot([4, 10], [4, 10], color=c.GRAY, lw=0.6, ls=(0, (3, 2)), zorder=0)
    median_offset = float(pair["offset_b_minus_a"].median())
    assert round(median_offset, 2) == 3.47
    ax.plot([4, 6.53], [7.47, 10], color="black", lw=0.6, zorder=1)
    ax.plot(
        pair["us_10131663_b2"],
        pair["us_10590130_b2"],
        linestyle="none",
        marker="o",
        markersize=3.4,
        markerfacecolor="none",
        markeredgecolor="black",
        markeredgewidth=0.7,
        zorder=2,
    )
    c.register_axis_points(ax, pair["us_10131663_b2"], pair["us_10590130_b2"])
    for _, row in pair.iterrows():
        c.record_value(4, "4b", "US-10590130-B2 vs US-10131663-B2", row["compound_id"], "us_10131663_b2", row["us_10131663_b2"], "4b_hERG_patent_pair.csv")
        c.record_value(4, "4b", "US-10590130-B2 vs US-10131663-B2", row["compound_id"], "us_10590130_b2", row["us_10590130_b2"], "4b_hERG_patent_pair.csv")
        c.record_value(4, "4b", "Patent offset", row["compound_id"], "offset_b_minus_a", row["offset_b_minus_a"], "4b_hERG_patent_pair.csv")

    spearman = float(pair[["us_10131663_b2", "us_10590130_b2"]].corr(method="spearman").iloc[0, 1])
    c.record_value(4, "4b", "Spearman rho", "summary", "spearman_rho", spearman, "4b_hERG_patent_pair.csv")
    ax.text(9.9, 8.6, "y = x", ha="right", va="top", fontsize=6.5, color=c.GRAY)
    ax.text(6.75, 9.7, "+3.47", ha="left", va="center", fontsize=6.5, color="black")
    ax.text(
        9.9,
        4.25,
        f"n = {len(pair)}\nSpearman ρ = {spearman:.2f}",
        ha="right",
        va="bottom",
        fontsize=6.5,
        color=c.INK,
    )

    return fig
