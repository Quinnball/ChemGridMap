"""SVG, PNG and PDF rendering for molecule-resolved grid maps."""

from __future__ import annotations

import html
import os
import re
import textwrap
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import to_rgb
from matplotlib.figure import Figure
from matplotlib.patches import Polygon, Rectangle
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D


DEFAULT_COLORS = {
    "inactive": "#C44E52",
    "medium": "#E6A23C",
    "active": "#2A9D8F",
    "unlabelled": "#B8C2CC",
}

MOLECULE_TILE_OPACITY = 0.12
OVERVIEW_TILE_OPACITY = 0.58

# Keep these colors in sync with web/app.js.
MAP_COLOR_MODES = ("activity", "class_support", "neighbor_recall", "structure_group")
FADED_COLORS = {"active": "#b8e0db", "medium": "#f6dfb8",
                "inactive": "#ebc0c1", "unlabelled": "#dfe4e9"}
RELIABILITY_BINS = ((0, "#e8eef7"), (.2, "#bcd0ea"), (.4, "#86acd8"),
                    (.6, "#4f82bf"), (.8, "#1d5689"))


def _tile_color(row, colors, label_col, color_by):
    label = str(row[label_col]).lower()
    if color_by == "class_support" and row.get("class_status") == "uncertain":
        return FADED_COLORS.get(label, colors[label])
    if color_by == "neighbor_recall":
        value = row.get("neighbor_recall")
        if pd.isna(value):
            return DEFAULT_COLORS["unlabelled"]
        return next(color for low, color in reversed(RELIABILITY_BINS) if value >= low - 1e-9)
    if color_by == "structure_group":
        value = row.get("structure_group_color")
        return value if pd.notna(value) and value else DEFAULT_COLORS["unlabelled"]
    return colors[label]


def _color_legend(data, color_by):
    if color_by == "neighbor_recall":
        return "Representation neighbors retained in the grid", [
            (color, f"{round(low * 100)}–{round((RELIABILITY_BINS[i + 1][0] if i < 4 else 1) * 100)}%")
            for i, (low, color) in enumerate(RELIABILITY_BINS)]
    if color_by == "class_support":
        entries = []
        for label in ("inactive", "medium", "active"):
            entries.extend([(DEFAULT_COLORS[label], label), (FADED_COLORS[label], f"{label}: uncertain")])
        threshold = data.attrs.get("activity_confidence_fit", {}).get("confidence_threshold")
        threshold_text = f"{threshold:g}" if threshold is not None else "configured threshold"
        if data.get("offset_pair_touched", pd.Series(dtype=bool)).fillna(False).astype(bool).any():
            entries.append(("#172536", "offset document pair (bottom-left marker)"))
        return f"Activity class support (pale: probability below {threshold_text} or offset document pair)", entries
    if color_by == "structure_group":
        if "structure_group_color" not in data:
            return "Structure groups", [(DEFAULT_COLORS["unlabelled"], "No grouping")]
        entries = []
        for color, rows in data.groupby("structure_group_color", sort=False):
            groups = rows.structure_group.astype(str).drop_duplicates().tolist()
            entries.append((color, ", ".join(groups[:3]) + (", …" if len(groups) > 3 else "")))
        return "Structure groups (colors may repeat; IDs identify groups)", entries
    return "", []

EXTRA_COLORS = [
    "#2B6CB0",
    "#D97706",
    "#2F855A",
    "#C53030",
    "#6B46C1",
    "#00838F",
    "#805AD5",
    "#718096",
]


def _label_colors(
    labels: Iterable[str],
    color_map: Optional[Dict[str, str]] = None,
) -> Dict[str, str]:
    colors = dict(DEFAULT_COLORS)
    if color_map:
        colors.update({str(key).lower(): value for key, value in color_map.items()})
    unknown = sorted(
        {str(label).lower() for label in labels if str(label).lower() not in colors}
    )
    for index, label in enumerate(unknown):
        colors[label] = EXTRA_COLORS[index % len(EXTRA_COLORS)]
    return colors


def _svg_molecule_content(smiles: str, size: int) -> Optional[str]:
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    drawer = rdMolDraw2D.MolDraw2DSVG(size, size)
    options = drawer.drawOptions()
    options.clearBackground = False
    options.padding = 0.04
    options.bondLineWidth = 4.0
    options.minFontSize = 10
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    svg = re.sub(r"<\?xml[^>]*\?>", "", svg).strip()
    start = svg.find(">")
    end = svg.rfind("</svg>")
    if start < 0 or end < 0:
        return None
    return svg[start + 1 : end]


@lru_cache(maxsize=4096)
def _molecule_rgba(smiles: str, size: int = 220) -> Optional[np.ndarray]:
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    options = rdMolDraw2D.MolDrawOptions()
    options.padding = 0.04
    options.bondLineWidth = 4.0
    options.minFontSize = 10
    image = Draw.MolToImage(
        mol,
        size=(size, size),
        options=options,
    ).convert("RGBA")
    array = np.asarray(image).copy()
    rgb = array[:, :, :3]
    white = (rgb[:, :, 0] > 245) & (rgb[:, :, 1] > 245) & (rgb[:, :, 2] > 245)
    array[:, :, 3] = np.where(white, 0, 255).astype(np.uint8)
    return array


def render_svg(
    data: pd.DataFrame,
    output_path: Path,
    smiles_col: str = "canonical_smiles",
    label_col: str = "activity_class",
    color_map: Optional[Dict[str, str]] = None,
    tile_size: int = 100,
    molecule_margin: int = 5,
    draw_molecules: bool = True,
    color_by: str = "activity",
) -> Path:
    """Render a cropped, molecule-resolved vector SVG map."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    labels = data[label_col].astype(str).str.lower()
    colors = _label_colors(labels, color_map)

    min_col = int(data["grid_col"].min())
    max_col = int(data["grid_col"].max())
    min_row = int(data["grid_row"].min())
    max_row = int(data["grid_row"].max())
    width = (max_col - min_col + 1) * tile_size
    height = (max_row - min_row + 1) * tile_size
    legend_title, legend_entries = _color_legend(data, color_by)
    map_height = height
    if legend_entries:
        width = max(width, 620)
        height += 50 + 25 * ((len(legend_entries) + 2) // 3)

    parts = [
        '<svg xmlns="http://www.w3.org/2000/svg" '
        'width="{0}" height="{1}" viewBox="0 0 {0} {1}">'.format(width, height),
        '<rect width="100%" height="100%" fill="white"/>',
    ]
    for _, row in data.iterrows():
        column = int(row["grid_col"]) - min_col
        grid_row = int(row["grid_row"]) - min_row
        x = column * tile_size
        y = (max_row - min_row - grid_row) * tile_size
        label = str(row[label_col]).lower()
        color = _tile_color(row, colors, label_col, color_by)
        identity = row.get("molecule_identity_key", row.get("name", row[smiles_col]))
        fields = ["molecule_id", "activity_pchembl", "n_records", "pchembl_min",
                  "pchembl_max", "is_conflicted", "has_class_boundary_crossing", "source_rows"]
        title = str(identity) + " | " + "; ".join(
            "{}={}".format(field, row[field]) for field in fields if field in row.index
        )
        parts.append('<g data-molecule-id="{}" data-grid-row="{}" data-grid-col="{}">'.format(
            html.escape(str(identity), quote=True), int(row["grid_row"]), int(row["grid_col"])))
        parts.append("<title>{}</title>".format(html.escape(str(title))))
        parts.append(
            '<rect x="{x}" y="{y}" width="{size}" height="{size}" '
            'fill="{color}" fill-opacity="{opacity}" stroke="#FFFFFF" '
            'stroke-opacity="0.95" stroke-width="1" shape-rendering="crispEdges"/>'.format(
                x=x,
                y=y,
                size=tile_size,
                color=color,
                opacity=(MOLECULE_TILE_OPACITY if draw_molecules else OVERVIEW_TILE_OPACITY),
            )
        )
        if draw_molecules:
            molecule = _svg_molecule_content(
                str(row[smiles_col]), tile_size - 2 * molecule_margin
            )
            if molecule:
                parts.append(
                    '<g transform="translate({},{})">{}</g>'.format(
                        x + molecule_margin, y + molecule_margin, molecule
                    )
                )
        if pd.notna(row.get("is_conflicted")) and row.get("is_conflicted") == 1:
            corner = tile_size * 0.12
            parts.append(
                f'<path data-conflict-marker="true" d="M{x+tile_size-corner},{y+1} '
                f'L{x+tile_size-1},{y+1} L{x+tile_size-1},{y+corner} Z" fill="#172536"/>'
            )
        if color_by == "class_support" and bool(row.get("offset_pair_touched", False)):
            corner = tile_size * 0.12
            parts.append(
                f'<path data-offset-pair-marker="true" d="M{x+1},{y+tile_size-corner} '
                f'L{x+1},{y+tile_size-1} L{x+corner},{y+tile_size-1} Z" fill="#172536"/>'
            )
        parts.append("</g>")
    if legend_entries:
        parts.append('<g data-color-legend="{}" font-family="Arial, sans-serif" font-size="14" fill="#172536">'.format(color_by))
        parts.append(f'<text x="12" y="{map_height + 24}">{html.escape(legend_title)}</text>')
        for index, (color, text) in enumerate(legend_entries):
            x, y = 12 + (index % 3) * 200, map_height + 40 + (index // 3) * 25
            parts.append(f'<rect x="{x}" y="{y}" width="14" height="14" fill="{color}"/>')
            parts.append(f'<text x="{x + 20}" y="{y + 12}">{html.escape(text)}</text>')
        parts.append("</g>")
    parts.append("</svg>")
    output_path.write_text("\n".join(parts), encoding="utf-8", newline="\n")
    return output_path


def render_raster_and_pdf(
    data: pd.DataFrame,
    png_path: Path,
    pdf_path: Path,
    smiles_col: str = "canonical_smiles",
    label_col: str = "activity_class",
    color_map: Optional[Dict[str, str]] = None,
    draw_molecules: bool = True,
    color_by: str = "activity",
) -> None:
    """Render high-resolution PNG and PDF versions with Matplotlib."""
    labels = data[label_col].astype(str).str.lower()
    colors = _label_colors(labels, color_map)
    min_col = int(data["grid_col"].min())
    max_col = int(data["grid_col"].max())
    min_row = int(data["grid_row"].min())
    max_row = int(data["grid_row"].max())
    width_cells = max_col - min_col + 1
    height_cells = max_row - min_row + 1

    if draw_molecules:
        figure_size = (
            max(4.0, min(18.0, width_cells * 0.24)),
            max(4.0, min(18.0, height_cells * 0.24)),
        )
    else:
        longest = float(max(width_cells, height_cells))
        figure_size = (
            max(3.5, 9.0 * width_cells / longest),
            max(3.5, 9.0 * height_cells / longest),
        )
    legend_title, legend_entries = _color_legend(data, color_by)
    legend_columns = max(1, min(3, int(figure_size[0] / 2.6)))
    legend_title = textwrap.fill(legend_title, width=max(30, int(figure_size[0] * 12)))
    legend_height = (.6 + .22 * ((len(legend_entries) + legend_columns - 1) // legend_columns)
                     + .16 * legend_title.count("\n")) if legend_entries else 0
    if legend_entries:
        # Extend below the map so its original pixel positions and scale stay fixed.
        figure_size = (figure_size[0], figure_size[1] + legend_height)
    figure = Figure(figsize=figure_size, facecolor="white")
    FigureCanvasAgg(figure)
    legend_fraction = legend_height / figure_size[1]
    axis = figure.add_axes([0, legend_fraction, 1, 1 - legend_fraction])

    if draw_molecules:
        for _, row in data.iterrows():
            x = int(row["grid_col"]) - min_col
            y = int(row["grid_row"]) - min_row
            label = str(row[label_col]).lower()
            axis.add_patch(
                Rectangle(
                    (x - 0.5, y - 0.5),
                    1,
                    1,
                    facecolor=_tile_color(row, colors, label_col, color_by),
                    edgecolor="white",
                    linewidth=0.4,
                    alpha=MOLECULE_TILE_OPACITY,
                    zorder=1,
                )
            )
            if color_by == "class_support" and bool(row.get("offset_pair_touched", False)):
                axis.add_patch(Polygon(
                    ((x - .5, y - .38), (x - .5, y - .14), (x - .14, y - .5), (x - .38, y - .5)),
                    closed=True, facecolor="#172536", edgecolor="none", zorder=5,
                ))
            molecule = _molecule_rgba(str(row[smiles_col]))
            if molecule is not None:
                axis.imshow(
                    molecule,
                    extent=(x - 0.47, x + 0.47, y - 0.47, y + 0.47),
                    interpolation="bilinear",
                    zorder=3,
                )
    else:
        image = np.ones((height_cells, width_cells, 3), dtype=np.float32)
        for _, row in data.iterrows():
            x = int(row["grid_col"]) - min_col
            y = int(row["grid_row"]) - min_row
            label = str(row[label_col]).lower()
            rgb = np.asarray(to_rgb(_tile_color(row, colors, label_col, color_by)), dtype=np.float32)
            image[y, x] = OVERVIEW_TILE_OPACITY * rgb + (1.0 - OVERVIEW_TILE_OPACITY)
        axis.imshow(
            image,
            extent=(-0.5, width_cells - 0.5, -0.5, height_cells - 0.5),
            origin="lower",
            interpolation="nearest",
        )
        touched_rows = data[
            data.get("offset_pair_touched", pd.Series(False, index=data.index))
            .fillna(False).astype(bool)
        ] if color_by == "class_support" else data.iloc[0:0]
        for _, row in touched_rows.iterrows():
            x = int(row["grid_col"]) - min_col
            y = int(row["grid_row"]) - min_row
            axis.add_patch(Polygon(
                ((x - .5, y - .5), (x - .5, y - .1), (x - .1, y - .5)),
                closed=True, facecolor="#172536", edgecolor="none", zorder=5,
            ))

    axis.set_xlim(-0.55, width_cells - 0.45)
    axis.set_ylim(-0.55, height_cells - 0.45)
    axis.set_aspect("equal")
    axis.axis("off")
    if legend_entries:
        handles = [Rectangle((0, 0), 1, 1, facecolor=color, label=text) for color, text in legend_entries]
        figure.legend(handles=handles, loc="lower center", ncol=legend_columns, frameon=False,
                      title=legend_title, fontsize=8, title_fontsize=9)
    Path(png_path).parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(str(png_path), dpi=300, facecolor="white")
    figure.savefig(str(pdf_path), dpi=300, facecolor="white")


def render_grid_map(
    data: pd.DataFrame,
    output_prefix: Path,
    smiles_col: str = "canonical_smiles",
    label_col: str = "activity_class",
    color_map: Optional[Dict[str, str]] = None,
    detail: str = "auto",
    formats: Optional[Sequence[str]] = None,
    color_by: str = "activity",
) -> Dict[str, str]:
    """Export the same grid map as SVG, PNG and PDF."""
    if color_by not in MAP_COLOR_MODES:
        raise ValueError("Unknown map color mode: {}.".format(color_by))
    required_column = {"class_support": "class_status", "neighbor_recall": "neighbor_recall"}.get(color_by)
    if required_column and required_column not in data:
        raise ValueError(f"{color_by} coloring requires the {required_column} column.")
    normalized_detail = str(detail).strip().lower()
    if normalized_detail not in {"auto", "full", "overview"}:
        raise ValueError("render detail must be 'auto', 'full' or 'overview'.")
    if normalized_detail == "auto":
        normalized_detail = "full" if len(data) <= 2000 else "overview"
    draw_molecules = normalized_detail == "full"

    requested = tuple(formats) if formats is not None else ("svg", "png", "pdf")
    requested = tuple(str(value).strip().lower() for value in requested)
    invalid = sorted(set(requested).difference({"svg", "png", "pdf"}))
    if invalid:
        raise ValueError("Unknown output format(s): {}.".format(", ".join(invalid)))

    output_prefix = Path(output_prefix)
    svg_path = output_prefix.with_suffix(".svg")
    png_path = output_prefix.with_suffix(".png")
    pdf_path = output_prefix.with_suffix(".pdf")
    output_files: Dict[str, str] = {}
    if "svg" in requested:
        render_svg(
            data,
            svg_path,
            smiles_col=smiles_col,
            label_col=label_col,
            color_map=color_map,
            draw_molecules=draw_molecules,
            color_by=color_by,
        )
        output_files["svg"] = str(svg_path)
    if "png" in requested or "pdf" in requested:
        render_raster_and_pdf(
            data,
            png_path,
            pdf_path,
            smiles_col=smiles_col,
            label_col=label_col,
            color_map=color_map,
            draw_molecules=draw_molecules,
            color_by=color_by,
        )
        if "png" in requested:
            output_files["png"] = str(png_path)
        elif png_path.exists():
            png_path.unlink()
        if "pdf" in requested:
            output_files["pdf"] = str(pdf_path)
        elif pdf_path.exists():
            pdf_path.unlink()
    return output_files
