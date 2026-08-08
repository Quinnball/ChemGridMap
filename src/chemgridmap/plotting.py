"""SVG, PNG and PDF rendering for molecule-resolved grid maps."""

from __future__ import annotations

import html
import os
import re
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, Optional

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D


DEFAULT_COLORS = {
    "inactive": "#A50F15",
    "medium": "#F46D43",
    "active": "#00441B",
    "unlabelled": "#B8C2CC",
}

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
    options.padding = 0.08
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
    image = Draw.MolToImage(mol, size=(size, size)).convert("RGBA")
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
    molecule_margin: int = 8,
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
        color = colors[label]
        title = row.get("name", row.get("Compound", row[smiles_col]))
        parts.append("<g>")
        parts.append("<title>{}</title>".format(html.escape(str(title))))
        parts.append(
            '<rect x="{x}" y="{y}" width="{size}" height="{size}" '
            'fill="{color}" fill-opacity="0.52" stroke="#FFFFFF" '
            'stroke-width="1"/>'.format(
                x=x, y=y, size=tile_size, color=color
            )
        )
        molecule = _svg_molecule_content(
            str(row[smiles_col]), tile_size - 2 * molecule_margin
        )
        if molecule:
            parts.append(
                '<g transform="translate({},{})">{}</g>'.format(
                    x + molecule_margin, y + molecule_margin, molecule
                )
            )
        parts.append("</g>")
    parts.append("</svg>")
    output_path.write_text("\n".join(parts), encoding="utf-8")
    return output_path


def render_raster_and_pdf(
    data: pd.DataFrame,
    png_path: Path,
    pdf_path: Path,
    smiles_col: str = "canonical_smiles",
    label_col: str = "activity_class",
    color_map: Optional[Dict[str, str]] = None,
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

    figure = Figure(
        figsize=(
            max(4.0, min(18.0, width_cells * 0.24)),
            max(4.0, min(18.0, height_cells * 0.24)),
        ),
        facecolor="white",
    )
    FigureCanvasAgg(figure)
    axis = figure.add_axes([0, 0, 1, 1])

    for _, row in data.iterrows():
        x = int(row["grid_col"]) - min_col
        y = int(row["grid_row"]) - min_row
        label = str(row[label_col]).lower()
        axis.add_patch(
            Rectangle(
                (x - 0.5, y - 0.5),
                1,
                1,
                facecolor=colors[label],
                edgecolor="white",
                linewidth=0.4,
                alpha=0.52,
            )
        )
        molecule = _molecule_rgba(str(row[smiles_col]))
        if molecule is not None:
            axis.imshow(
                molecule,
                extent=(x - 0.43, x + 0.43, y - 0.43, y + 0.43),
                interpolation="bilinear",
                zorder=2,
            )

    axis.set_xlim(-0.55, width_cells - 0.45)
    axis.set_ylim(-0.55, height_cells - 0.45)
    axis.set_aspect("equal")
    axis.axis("off")
    Path(png_path).parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(str(png_path), dpi=300, facecolor="white")
    figure.savefig(str(pdf_path), dpi=300, facecolor="white")


def render_grid_map(
    data: pd.DataFrame,
    output_prefix: Path,
    smiles_col: str = "canonical_smiles",
    label_col: str = "activity_class",
    color_map: Optional[Dict[str, str]] = None,
) -> Dict[str, str]:
    """Export the same grid map as SVG, PNG and PDF."""
    output_prefix = Path(output_prefix)
    svg_path = output_prefix.with_suffix(".svg")
    png_path = output_prefix.with_suffix(".png")
    pdf_path = output_prefix.with_suffix(".pdf")
    render_svg(
        data,
        svg_path,
        smiles_col=smiles_col,
        label_col=label_col,
        color_map=color_map,
    )
    render_raster_and_pdf(
        data,
        png_path,
        pdf_path,
        smiles_col=smiles_col,
        label_col=label_col,
        color_map=color_map,
    )
    return {
        "svg": str(svg_path),
        "png": str(png_path),
        "pdf": str(pdf_path),
    }

