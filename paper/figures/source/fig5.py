"""Rebuild Figure 5: neighbourhoods and cores."""

from __future__ import annotations

from collections import Counter

import numpy as np
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

import common as c


W, H = 170.0, 150.0
LAYOUT_COLORS = {
    "projection": c.GRAY,
    "minimum_displacement": c.MUTED,
    "neighbor_refined": c.GREEN,
}
LAYOUT_BAR_COLORS = {
    "projection": c.shade(c.GRAY, 0.35),
    "minimum_displacement": "#B3B3B3",
    "neighbor_refined": c.shade(c.GREEN, 0.35),
}
LAYOUT_LABELS = {
    "projection": "Projection",
    "minimum_displacement": "Minimum-displacement grid",
    "neighbor_refined": "Refined grid",
}


def _recall_panel(fig: Figure) -> None:
    frame = c.load_csv("5a_recall.csv")
    fixed = c.load_csv("5a_fixedk_means.csv").set_index(["target_label", "layout"])
    ax = c.add_axes_mm(fig, (11, 106.5, 46, 36), (W, H))
    ax.set_xlim(-0.5, 2.5)
    ax.set_ylim(0, 1)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["CA II", "Thrombin", "hERG"])
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.set_ylabel("Neighbor recall (k = 10)")
    c.set_bounds(ax, y=(0, 1))

    layouts = ["projection", "minimum_displacement", "neighbor_refined"]
    offsets = [-0.27, 0.0, 0.27]
    for layout, offset in zip(layouts, offsets, strict=True):
        columns = [f"{layout}__r{i:02d}" for i in range(1, 11)]
        # Each column is one layout across ten seeds.
        for target_x, (_, row) in zip(range(3), frame.iterrows(), strict=True):
            values = row[columns].to_numpy(dtype=float)
            mean = float(values.mean())
            ax.bar(
                target_x + offset,
                mean,
                width=0.25,
                color=LAYOUT_BAR_COLORS[layout],
                edgecolor="none",
                zorder=1,
            )
            fixed_mean = float(fixed.loc[(row["Row_title"], layout), "fixed_recall_mean"])
            ax.hlines(
                fixed_mean,
                target_x + offset - 0.075,
                target_x + offset + 0.075,
                color="black",
                lw=1.0,
                zorder=4,
            )
            c.register_axis_points(
                ax,
                [target_x + offset - 0.075, target_x + offset + 0.075],
                [fixed_mean, fixed_mean],
            )
            seed_x = target_x + offset + np.linspace(-0.07, 0.07, 10)
            ax.plot(
                seed_x,
                values,
                linestyle="none",
                marker="o",
                markersize=1.8,
                markerfacecolor=LAYOUT_COLORS[layout],
                markeredgecolor="white",
                markeredgewidth=0.4,
                zorder=3,
            )
            c.register_axis_points(ax, [target_x + offset, target_x + offset], [0, mean])
            c.register_axis_points(ax, seed_x, values)
            for seed_index, value in enumerate(values, start=1):
                c.record_value(5, "5a", layout, f"{row['Row_title']} r{seed_index:02d}", columns[seed_index - 1], value, "5a_recall.csv")
            c.record_value(5, "5a", f"{layout} mean", row["Row_title"], "mean_recall", mean, "5a_recall.csv")
            c.record_value(
                5,
                "5a",
                f"{layout} fixed k mean",
                row["Row_title"],
                "fixed_recall_mean",
                fixed_mean,
                "5a_fixedk_means.csv",
            )

    handles = [
        Rectangle((0, 0), 1, 1, facecolor=LAYOUT_BAR_COLORS[key], edgecolor="none")
        for key in layouts
    ]
    handles.append(
        Line2D(
            [0],
            [0],
            color="black",
            lw=1.0,
            label="Fixed k = 10",
        )
    )
    ax.legend(
        handles,
        [LAYOUT_LABELS[key] for key in layouts] + ["Fixed k = 10"],
        loc="upper right",
        bbox_to_anchor=(0.992, 0.985),
        frameon=False,
        fontsize=6.5,
        handlelength=1.0,
        handleheight=0.8,
        handletextpad=0.35,
        labelspacing=0.25,
        borderaxespad=0,
    )


def _contrast_panel(fig: Figure) -> None:
    frame = c.load_csv("5b_layout_contrasts.csv")
    ax = c.add_axes_mm(fig, (85, 106.5, 25, 38), (W, H))
    ax.set_xlim(-0.05, 0.05)
    ax.set_xticks([-0.04, 0, 0.04])
    ax.set_xticks(np.arange(-0.04, 0.041, 0.02), minor=True)
    ax.set_ylim(-0.5, 14.5)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel(r"$\Delta e_\mathrm{grid}$ (other − refined)", fontsize=7)
    ax.xaxis.label._qa_fig5b_axis_title = True  # type: ignore[attr-defined]
    # Reserve space for group names and break the reference line only at the CA II title.
    ax.vlines(
        [0, 0],
        [-0.5, 14.5],
        [13.5, 14.5],
        color=c.GRAY,
        lw=0.6,
        linestyles=(0, (3, 2)),
        zorder=0,
    )

    comparison_labels = [
        ("minimum_displacement_minus_neighbor_refined", "Min. displacement", True),
        ("las_random_minus_neighbor_refined", "LAS, random init.", False),
        ("las_projection_minus_neighbor_refined", "LAS, projection init.", False),
        ("flas_random_minus_neighbor_refined", "FLAS", False),
    ]
    y_tick_positions = []
    y_tick_labels = []
    for group_index, target in enumerate(["CA II", "Thrombin", "hERG"]):
        group_y = 14 - group_index * 5
        fig.text(
            62 / W,
            ax.transData.transform((0, group_y))[1] / (fig.get_figheight() * fig.dpi),
            target,
            ha="left",
            va="center",
            fontsize=6.5,
            fontweight="bold",
            color=c.INK,
        )
        local = frame[frame["target_label"] == target].set_index("comparison")
        for row_index, (comparison, label, filled) in enumerate(comparison_labels, start=1):
            y = group_y - row_index
            row = local.loc[comparison]
            value = float(row["other_minus_refined_e_grid"])
            low = float(row["ci_low"])
            high = float(row["ci_high"])
            ax.hlines(y, low, high, color="black", lw=0.7, zorder=1)
            ax.plot(
                value,
                y,
                marker="o",
                markersize=3.0,
                markerfacecolor="black" if filled else "none",
                markeredgecolor="black",
                markeredgewidth=0.7,
                linestyle="none",
                zorder=2,
            )
            c.register_axis_points(ax, [low, value, high], [y] * 3)
            c.record_value(5, "5b", label, f"{target} {comparison}", "other_minus_refined_e_grid", value, "5b_layout_contrasts.csv")
            c.record_value(5, "5b", f"{label} 95% CI", f"{target} {comparison}", "ci_low", low, "5b_layout_contrasts.csv")
            c.record_value(5, "5b", f"{label} 95% CI", f"{target} {comparison}", "ci_high", high, "5b_layout_contrasts.csv")
            y_tick_positions.append(y)
            y_tick_labels.append(label)
    ax.set_yticks(y_tick_positions)
    ax.set_yticklabels(y_tick_labels, fontsize=6)
    ax.tick_params(axis="y", pad=2)

    ca_y = 14
    fig.text(
        110 / W,
        ax.transData.transform((0, ca_y))[1] / (fig.get_figheight() * fig.dpi),
        "Refined better →",
        ha="right",
        va="center",
        fontsize=6,
        color=c.MUTED,
    )
    fig.legend(
        handles=[
            c.ci_handle(marker="o", marker_size=3.0, color="black", hollow=False, line=False),
            c.ci_handle(marker="o", marker_size=3.0, color="black", hollow=True, line=False),
        ],
        labels=["Predeclared", "Post hoc"],
        loc="center left",
        bbox_to_anchor=(66 / W, 147.0 / H),
        frameon=False,
        fontsize=6.5,
        ncol=2,
        handlelength=1.1,
        handletextpad=0.35,
        columnspacing=1.1,
        borderaxespad=0,
    )


def _neighbor_grids(fig: Figure) -> str:
    neighbors = c.load_csv(
        "Figure_3C_CHEMBL81688_neighbors.csv",
        root=c.CGM_DIR / "paper/output/reliability_v1/figure_data",
    )
    coordinates = c.load_csv(
        "coordinates_seed_42.csv",
        root=c.CGM_DIR / "paper/output/reliability_v1/layouts/chembl205",
    )
    assert len(coordinates) == 889
    coordinate_index = coordinates.set_index("molecule_identity_key")
    query = neighbors[neighbors["is_query"]].iloc[0]
    smiles = str(query["canonical_smiles"])

    styles = {
        "Q": ("#1A1A1A", "none", "white", ""),
        "S": (c.BLUE, "none", "white", ""),
        "G": (c.ORANGE, "none", "white", ""),
        "F": ("#7A5195", "none", "white", ""),
    }
    letter_cap_mm = c.cap_height_mm(6)
    glyph_x_nudge_mm = {"S": 0.09, "G": 0.095}
    cell_mm = 22.7 / 9
    audits = {}
    for layout, x0, title in (
        ("minimum_displacement", 116.0, "Before refinement"),
        ("neighbor_refined", 142.5, "After refinement"),
    ):
        local = neighbors[neighbors["layout"] == layout]
        predicate = "minimum_displacement" if layout == "minimum_displacement" else "neighbor_refined"
        counts = Counter(local["figure_marker"])
        assert counts["Q"] == 1
        if predicate == "minimum_displacement":
            assert (counts["S"], counts["G"], counts["F"]) == (3, 9, 8)
            col_prefix, row_prefix = "grid_col", "grid_row"
        else:
            assert (counts["S"], counts["G"], counts["F"]) == (5, 7, 6)
            col_prefix, row_prefix = "refined_grid_col", "refined_grid_row"
        query_coordinate = coordinate_index.loc[query["molecule_identity_key"]]
        assert (int(query_coordinate[col_prefix]), int(query_coordinate[row_prefix])) == (11, 28)
        for _, row in local.iterrows():
            actual = coordinate_index.loc[row["molecule_identity_key"]]
            assert (int(actual[col_prefix]), int(actual[row_prefix])) == (
                int(row["grid_col"]),
                int(row["grid_row"]),
            )

        ax = c.add_axes_mm(fig, (x0, 120, 22.7, 22.7), (W, H))
        ax.set_xlim(6.5, 15.5)
        ax.set_ylim(23.5, 32.5)
        ax.set_aspect("equal", adjustable="box")
        ax.set_axis_off()
        ax.add_patch(
            Rectangle(
                (6.5, 23.5),
                9,
                9,
                facecolor="none",
                edgecolor="#BFBFBF",
                linewidth=0.4,
                clip_on=False,
            )
        )

        # Draw other entries in the neighbourhood window first, then overlay Q/S/G/F.
        in_window = coordinates[
            coordinates[col_prefix].between(7, 15)
            & coordinates[row_prefix].between(24, 32)
        ]
        for _, row in in_window.iterrows():
            col = int(row[col_prefix])
            y_row = int(row[row_prefix])
            ax.add_patch(
                Rectangle(
                    (col - 0.4, y_row - 0.4),
                    0.8,
                    0.8,
                    facecolor="#EDEDED",
                    edgecolor="none",
                    zorder=1,
                )
            )
            c.record_value(
                5,
                "5c",
                f"{predicate} other entries",
                row["molecule_identity_key"],
                f"{col_prefix}/{row_prefix}",
                col * 1000 + y_row,
                "coordinates_seed_42.csv",
            )

        for _, row in local.iterrows():
            marker = row["figure_marker"]
            col = int(row["grid_col"])
            y_row = int(row["grid_row"])
            face, edge, text_color, _ = styles[marker]
            ax.add_patch(
                Rectangle(
                    (col - 0.4, y_row - 0.4),
                    0.8,
                    0.8,
                    facecolor=face,
                    edgecolor=edge,
                    linewidth=0,
                    zorder=3,
                )
            )
            label = ax.text(
                col + glyph_x_nudge_mm.get(marker, 0.0) / cell_mm,
                y_row - 0.5 * letter_cap_mm / cell_mm,
                marker,
                ha="center",
                va="baseline",
                fontsize=6,
                fontweight="bold",
                color=text_color,
                zorder=4,
            )
            label._qa_cell_label = True  # type: ignore[attr-defined]
            label._qa_cell_letter = marker  # type: ignore[attr-defined]
            label._qa_cell_center_mm = (  # type: ignore[attr-defined]
                x0 + (col - 6.5) * cell_mm,
                120.0 + (y_row - 23.5) * cell_mm,
            )
            c.record_value(
                5,
                "5c",
                f"{predicate} {marker}",
                row["molecule_identity_key"],
                f"{col_prefix}/{row_prefix}",
                col * 1000 + y_row,
                "Figure_3C_CHEMBL81688_neighbors.csv",
            )

        c.mark_panel_title(fig.text(
            (x0 + 11.35) / W,
            143.3 / H,
            title,
            ha="center",
            va="bottom",
            fontsize=6.5,
            color=c.INK,
        ))
        fig.text(
            (x0 + 11.35) / W,
            119.2 / H,
            f"S {counts['S']}   G {counts['G']}   F {counts['F']}",
            ha="center",
            va="top",
            fontsize=6.5,
            color=c.INK,
        )
        audits[predicate] = {"S": counts["S"], "G": counts["G"], "F": counts["F"]}

    # 5c legend
    legend_rows = [
        ("Q", "white", "Query"),
        ("S", "white", "Shared"),
        ("G", "black", "Grid only"),
        ("F", "black", "Fingerprint only"),
        ("", "black", "Other entry"),
    ]
    for index, (marker, text_color, label) in enumerate(legend_rows):
        y = 113.5 - 2.9 * index
        face = "#EDEDED" if not marker else styles[marker][0]
        patch = Rectangle(
            (145 / W, (y - 1.25) / H),
            2.5 / W,
            2.5 / H,
            transform=fig.transFigure,
            facecolor=face,
            edgecolor="none",
            linewidth=0,
            clip_on=False,
        )
        fig.add_artist(patch)
        if marker:
            cell_text = fig.text(
                (146.25 + glyph_x_nudge_mm.get(marker, 0.0)) / W,
                (y - 0.5 * letter_cap_mm) / H,
                marker,
                ha="center",
                va="baseline",
                fontsize=6,
                fontweight="bold",
                color=styles[marker][2],
            )
            cell_text._qa_cell_label = True  # type: ignore[attr-defined]
            cell_text._qa_cell_letter = marker  # type: ignore[attr-defined]
            cell_text._qa_cell_center_mm = (146.25, y)  # type: ignore[attr-defined]
        c.mark_legend_text(fig.text(
            148.5 / W,
            y / H,
            label,
            ha="left",
            va="center",
            fontsize=6.5,
            color=c.INK,
        ))

    caption = fig.text(
        127.5 / W,
        100.3 / H,
        "CHEMBL81688 (Q)",
        ha="center",
        va="top",
        fontsize=6,
        color=c.MUTED,
    )
    caption._qa_fig5_c_caption = True  # type: ignore[attr-defined]
    return smiles


def _correlation_panel(fig: Figure) -> None:
    frame = c.load_csv("5d_layout_correlation.csv")
    sampled = frame[frame["row_type"] == "sampled_dataset"]
    references = frame[frame["row_type"] == "detailed_reference"]

    panels = {
        "IC50": (8, 16, 34, 70),
        "Ki": (47, 35, 34, 51),
        "reference": (47, 16, 34, 10.5),
    }
    for endpoint in ("IC50", "Ki"):
        local = sampled[sampled["endpoint"] == endpoint].sort_values("dataset_rank")
        ax = c.add_axes_mm(fig, panels[endpoint], (W, H))
        count = len(local)
        ax.set_xlim(-0.32, 0.22)
        ax.set_ylim(0.5, count + 0.5)
        ax.set_xticks([-0.2, 0, 0.2])
        ax.set_xticks([-0.3, -0.1, 0.1], minor=True)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        c.set_bounds(ax, x=(-0.3, 0.2))
        ax.axvline(0, color=c.GRAY, lw=0.6, ls=(0, (3, 2)), zorder=0)
        ax.hlines(
            local["dataset_rank"],
            local["rho_ci_low"],
            local["rho_ci_high"],
            color="#A6A6A6",
            lw=0.5,
            zorder=1,
        )
        ax.plot(
            local["rho"],
            local["dataset_rank"],
            linestyle="none",
            marker="o",
            markersize=2.2,
            markerfacecolor="#333333",
            markeredgecolor="#333333",
            markeredgewidth=0.4,
            zorder=2,
        )
        c.register_axis_points(ax, local["rho_ci_low"], local["dataset_rank"])
        c.register_axis_points(ax, local["rho_ci_high"], local["dataset_rank"])
        c.register_axis_points(ax, local["rho"], local["dataset_rank"])
        for _, row in local.iterrows():
            row_id = row["target_id"]
            c.record_value(5, "5d", endpoint, row_id, "rho", row["rho"], "5d_layout_correlation.csv")
            c.record_value(5, "5d", f"{endpoint} CI", row_id, "rho_ci_low", row["rho_ci_low"], "5d_layout_correlation.csv")
            c.record_value(5, "5d", f"{endpoint} CI", row_id, "rho_ci_high", row["rho_ci_high"], "5d_layout_correlation.csv")
        median = float(local["rho"].median())
        ax.axvline(median, color="black", lw=0.8, zorder=3)
        ax.text(
            median - 0.01,
            ax.get_ylim()[1],
            f"median {median:.3f}".replace("-", "−"),
            ha="right",
            va="top",
            fontsize=6,
            color=c.INK,
        )
        c.record_value(5, "5d", f"{endpoint} median", "summary", "median", median, "5d_layout_correlation.csv")
        c.mark_panel_title(fig.text(
            panels[endpoint][0] / W,
            86.8 / H,
            f"{endpoint} ({count} datasets)",
            ha="left",
            va="bottom",
            fontsize=6.5,
            color=c.INK,
        ))

    ax_ref = c.add_axes_mm(fig, panels["reference"], (W, H))
    ax_ref.set_xlim(-0.32, 0.22)
    ax_ref.set_ylim(0.5, 3.5)
    ax_ref.set_xticks([-0.2, 0, 0.2])
    ax_ref.set_xticks([-0.3, -0.1, 0.1], minor=True)
    ax_ref.yaxis.tick_right()
    ax_ref.tick_params(axis="y", length=0, pad=2, labelsize=6)
    ax_ref.spines["left"].set_visible(False)
    c.set_bounds(ax_ref, x=(-0.3, 0.2))
    ax_ref.axvline(0, color=c.GRAY, lw=0.6, ls=(0, (3, 2)), zorder=0)
    reference_order = ["CA II", "Thrombin", "hERG"]
    ref_plot = references.set_index("target_label").loc[reference_order].reset_index()
    ref_y = np.array([3, 2, 1])
    ax_ref.hlines(
        ref_y,
        ref_plot["rho_ci_low"],
        ref_plot["rho_ci_high"],
        color=c.MUTED,
        lw=0.7,
        zorder=1,
    )
    ax_ref.plot(
        ref_plot["rho"],
        ref_y,
        linestyle="none",
        marker="o",
        markersize=2.6,
        markerfacecolor="black",
        markeredgecolor="black",
        markeredgewidth=0.4,
        zorder=2,
    )
    ax_ref.set_yticks(ref_y)
    ax_ref.set_yticklabels(reference_order)
    c.register_axis_points(ax_ref, ref_plot["rho_ci_low"], ref_y)
    c.register_axis_points(ax_ref, ref_plot["rho_ci_high"], ref_y)
    c.register_axis_points(ax_ref, ref_plot["rho"], ref_y)
    for y, (_, row) in zip(ref_y, ref_plot.iterrows(), strict=True):
        c.record_value(5, "5d", "Detailed targets", row["target_label"], "rho", row["rho"], "5d_layout_correlation.csv")
        c.record_value(5, "5d", "Detailed targets CI", row["target_label"], "rho_ci_low", row["rho_ci_low"], "5d_layout_correlation.csv")
        c.record_value(5, "5d", "Detailed targets CI", row["target_label"], "rho_ci_high", row["rho_ci_high"], "5d_layout_correlation.csv")
    c.mark_panel_title(fig.text(47 / W, 27.3 / H, "Detailed targets (IC50)", ha="left", va="bottom", fontsize=6.5, color=c.INK))
    fig.text(3.5 / W, 51 / H, "Datasets, ordered by ρ", ha="center", va="center", rotation=90, fontsize=7, color=c.INK)
    fig.text(
        44.5 / W,
        11.5 / H,
        r"Spearman $\rho$ (recall, $e_\mathrm{grid} - e_\mathrm{fp}$)",
        ha="center",
        va="top",
        fontsize=7,
        color=c.INK,
    )


def _coverage_panel(fig: Figure) -> None:
    frame = c.load_csv("5e_core_coverage_by_group_size.csv")
    rules = ["most_frequent", "coverage", "fmcs"]
    titles = ["Most frequent\nscaffold", "Coverage\ntarget", "FMCS\n(10-s budget)"]
    strata = ["2", "3-4", "5-9", ">=10"]
    category_columns = [
        "coverage_ge_0_8",
        "coverage_gt_0_5_lt_0_8",
        "coverage_le_0_5",
        "no_core",
    ]
    category_colors = ["#08519C", "#3182BD", "#7FB2DA", "#D9D9D9"]
    upper_x = [104.0, 126.5, 149.0]
    for index, (rule, title, x0) in enumerate(zip(rules, titles, upper_x, strict=True)):
        local = frame[frame["rule"] == rule].set_index("size_stratum").loc[strata].reset_index()
        top = c.add_axes_mm(fig, (x0, 37.0, 19.5, 42.5), (W, H))
        bottom = c.add_axes_mm(fig, (x0, 16, 19.5, 18), (W, H))
        top.set_xlim(-0.6, 3.6)
        bottom.set_xlim(-0.6, 3.6)
        top.set_ylim(0, 1)
        bottom.set_ylim(14, 27)
        top.set_yticks([0, 0.5, 1])
        bottom.set_yticks([15, 20, 25])
        top.tick_params(axis="y", labelleft=index == 0)
        bottom.tick_params(axis="y", labelleft=index == 0)
        top.tick_params(axis="x", labelbottom=False)
        bottom.set_xticks([0, 1, 2, 3])
        bottom.set_xticklabels(["2", "3–4", "5–9", "≥10"], fontsize=6)
        if index == 0:
            top.set_ylabel("Fraction of groups")
            bottom.set_ylabel("Heavy atoms")

        x_positions = np.arange(4)
        bottom_values = np.zeros(4)
        for column, color in zip(category_columns, category_colors, strict=True):
            groups = local["groups"].to_numpy(dtype=float)
            counts = local[column].to_numpy(dtype=float)
            fractions = counts / groups
            top.bar(
                x_positions,
                fractions,
                bottom=bottom_values,
                width=0.7,
                color=color,
                edgecolor="white",
                linewidth=0.5,
                zorder=1,
            )
            for x, count, groups_n, fraction in zip(x_positions, counts, groups, fractions, strict=True):
                c.record_value(
                    5,
                    "5e",
                    f"{rule} {column}",
                    f"group {x + 2}",
                    "fraction",
                    fraction,
                    "5e_core_coverage_by_group_size.csv",
                )
                c.record_value(
                    5,
                    "5e",
                    f"{rule} {column}",
                    f"group {x + 2}",
                    "count",
                    count,
                    "5e_core_coverage_by_group_size.csv",
                )
            bottom_values += fractions

        for x, groups_n in zip(x_positions, local["groups"].to_numpy(dtype=float), strict=True):
            c.record_value(
                5,
                "5e",
                f"{rule} groups",
                f"group {x + 2}",
                "groups",
                groups_n,
                "5e_core_coverage_by_group_size.csv",
            )

        heavy = local["mean_core_heavy_atoms"].to_numpy(dtype=float)
        bottom.plot(x_positions, heavy, color=c.GRAY, lw=0.6, marker="o", markersize=2.8, markerfacecolor="black", markeredgecolor="black", markeredgewidth=0.4, zorder=2)
        c.register_axis_points(bottom, x_positions, heavy)
        for x, value in zip(x_positions, heavy, strict=True):
            c.record_value(5, "5e", f"{rule} heavy atoms", f"group {x + 2}", "mean_core_heavy_atoms", value, "5e_core_coverage_by_group_size.csv")
        title_text = c.mark_panel_title(
            fig.text(
                (x0 + 9.75) / W,
                80.3 / H,
                title,
                ha="center",
                va="bottom",
                fontsize=6.5,
                color=c.INK,
            )
        )
        title_text._qa_fig5e_title = True  # type: ignore[attr-defined]

    coverage_label = c.mark_legend_text(
        fig.text(101.5 / W, 92.6 / H, "Coverage:", ha="left", va="center", fontsize=6.5, color=c.INK)
    )
    coverage_label._qa_fig5e_legend = True  # type: ignore[attr-defined]
    legend_rows = (
        (115.0, 92.6, "≥0.8"),
        (132.0, 92.6, "> 0.5 and < 0.8"),
        (115.0, 88.4, "≤0.5"),
        (132.0, 88.4, "No ring scaffold"),
    )
    for (x, y, label), color in zip(legend_rows, category_colors, strict=True):
        fig.add_artist(
            Rectangle(
                (x / W, (y - 1.25) / H),
                2.5 / W,
                2.5 / H,
                transform=fig.transFigure,
                facecolor=color,
                edgecolor="none",
                clip_on=False,
            )
        )
        legend_text = c.mark_legend_text(
            fig.text((x + 3.3) / W, y / H, label, ha="left", va="center", fontsize=6.5, color=c.INK)
        )
        legend_text._qa_fig5e_legend = True  # type: ignore[attr-defined]
    fig.text(
        136.25 / W,
        11.5 / H,
        "Group size (members)",
        ha="center",
        va="top",
        fontsize=7,
        color=c.INK,
    )


def build() -> tuple[Figure, str]:
    c.init()
    fig = c.mm_figure(W, H)
    c.panel_letter(fig, "a", 1, 149, (W, H))
    c.panel_letter(fig, "b", 60, 149, (W, H))
    c.panel_letter(fig, "c", 112, 149, (W, H))
    c.panel_letter(fig, "d", 1, 94, (W, H))
    c.panel_letter(fig, "e", 97.0, 94.0, (W, H))
    _recall_panel(fig)
    _contrast_panel(fig)
    smiles = _neighbor_grids(fig)
    _correlation_panel(fig)
    _coverage_panel(fig)
    c.write_plotted_values(5)
    return fig, smiles
