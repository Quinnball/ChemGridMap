"""Rebuild the README preview from saved coordinates, without rerunning the map.

Run with .[paper] installed. The bundled plotting table contains the original
CHEMBL205 seed-42 positions; the detail is the archived mixed 3 x 3 window.
"""
from __future__ import annotations

import argparse
import hashlib
import html
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import pandas as pd
import pymupdf
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
DATA = DOCS / "chemgridmap_demo_data.csv"
NS = "http://www.w3.org/2000/svg"
INK = "#172536"
MUTED = "#536579"
BLUE = "#2766A5"
COLORS = {"inactive": "#C44E52", "medium": "#E6A23C", "active": "#2A9D8F"}
PALE = {"inactive": "#FCF1F1", "medium": "#FFF8ED", "active": "#F0F9F7"}
FIELDS = ["molecule_identity_key", "canonical_smiles", "activity_pchembl",
          "activity_class", "n_records", "is_conflicted", "grid_col", "grid_row"]
WINDOW_COL, WINDOW_ROW = 13, 12
WIDTH, HEIGHT = 1440, 900
DETAIL_LEFT, DETAIL_TOP = 806, 182
DETAIL_STEP_X, DETAIL_STEP_Y = 196, 208
DETAIL_TILE_WIDTH, DETAIL_TILE_HEIGHT = 186, 198


def text(x, y, value, size=20, color=INK, bold=False, anchor="start"):
    return (f'<text x="{x}" y="{y}" font-family="Arial,Helvetica,sans-serif" '
            f'font-size="{size}" font-weight="{700 if bold else 400}" '
            f'fill="{color}" text-anchor="{anchor}">{html.escape(str(value))}</text>')


def rect(x, y, width, height, fill, stroke="none", sw=1):
    return (f'<rect x="{x}" y="{y}" width="{width}" height="{height}" '
            f'fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')


def molecule(smiles, x, y, width=170, height=148):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid structure: {smiles}")
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    options = drawer.drawOptions()
    options.clearBackground = False
    options.padding = 0.05
    options.bondLineWidth = 1.8
    options.minFontSize = 10
    options.maxFontSize = 15
    # Dark atom labels remain legible on the very pale activity fills.
    options.updateAtomPalette({7: (0.10, 0.23, 0.70), 8: (0.76, 0.13, 0.12),
                              9: (0.10, 0.40, 0.24), 16: (0.55, 0.37, 0.02),
                              17: (0.10, 0.40, 0.24)})
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    root = ET.fromstring(drawer.GetDrawingText())
    return (f'<g transform="translate({x},{y})">'
            + "".join(ET.tostring(child, encoding="unicode") for child in root)
            + "</g>")


def compose(data):
    if data.molecule_identity_key.duplicated().any() or data.duplicated(["grid_col", "grid_row"]).any():
        raise ValueError("The preview requires one molecule per occupied cell.")
    detail = data[data.grid_col.between(WINDOW_COL, WINDOW_COL + 2)
                  & data.grid_row.between(WINDOW_ROW, WINDOW_ROW + 2)]
    if len(detail) != 9:
        raise ValueError("The archived detail window must contain exactly nine cells.")

    parts = [f'<svg xmlns="{NS}" width="{WIDTH}" height="{HEIGHT}" '
             f'viewBox="0 0 {WIDTH} {HEIGHT}" role="img" aria-labelledby="title desc">',
             '<title id="title">ChemGridMap overview and enlarged molecular structures</title>',
             '<desc id="desc">889 CA II molecules occupy their archived grid cells. '
             'The blue box identifies the nine cells enlarged on the right, in the same order. '
             'Activity classes and median pChEMBL values are shown for each structure.</desc>',
             rect(0, 0, WIDTH, HEIGHT, "white"),
             rect(36, 34, 6, 65, "#2A9D8F"),
             text(58, 60, "From a chemical map to individual molecules", 32, bold=True),
             text(58, 94, "ChEMBL 37 / CHEMBL205 (CA II)  |  889 molecules  |  Morgan + UMAP", 20, MUTED),
             text(36, 148, "A  Complete map", 24, bold=True),
             text(806, 148, "B  The same nine cells, enlarged", 24, bold=True)]

    # Use integer cell indices for both panels so that the selection cannot drift.
    left, top, side = 36, 184, 660
    min_col, max_col = int(data.grid_col.min()), int(data.grid_col.max())
    min_row, max_row = int(data.grid_row.min()), int(data.grid_row.max())
    step = side / max(max_col - min_col + 1, max_row - min_row + 1)
    for row in data.itertuples():
        x, y = left + (row.grid_col - min_col) * step, top + (max_row - row.grid_row) * step
        parts += [f'<g data-molecule-id="{html.escape(row.molecule_identity_key)}">',
                  rect(x + .5, y + .5, step - 1, step - 1, COLORS[row.activity_class]), '</g>']
    box_x = left + (WINDOW_COL - min_col) * step
    box_y = top + (max_row - WINDOW_ROW - 2) * step
    parts += [rect(box_x - 3, box_y - 3, 3 * step + 6, 3 * step + 6, "none", "white", 7),
              rect(box_x - 3, box_y - 3, 3 * step + 6, 3 * step + 6, "none", BLUE, 3),
              rect(box_x - 4, box_y - 30, 25, 25, BLUE),
              text(box_x + 8.5, box_y - 11, "B", 18, "white", True, "middle")]
    # Share the detail-panel geometry with the line and arrowhead.
    start_x, start_y = box_x + 3 * step + 3, box_y + 1.5 * step
    end_x = DETAIL_LEFT
    end_y = DETAIL_TOP + (2 * DETAIL_STEP_Y + DETAIL_TILE_HEIGHT) / 2
    bend_x = (left + side + end_x) / 2
    head_base_x = end_x - 12
    parts += [f'<path id="detail-connector" d="M{start_x},{start_y} H{bend_x} '
              f'V{end_y} H{head_base_x}" fill="none" stroke="{BLUE}" '
              'stroke-width="2" stroke-linejoin="round"/>',
              f'<path id="detail-arrowhead" d="M{head_base_x},{end_y - 6} '
              f'L{end_x},{end_y} L{head_base_x},{end_y + 6} Z" fill="{BLUE}"/>']

    for row in detail.itertuples():
        x = DETAIL_LEFT + (row.grid_col - WINDOW_COL) * DETAIL_STEP_X
        y = DETAIL_TOP + (WINDOW_ROW + 2 - row.grid_row) * DETAIL_STEP_Y
        label = row.activity_class
        parts += [f'<g data-detail-molecule-id="{html.escape(row.molecule_identity_key)}" '
                  f'data-grid-col="{row.grid_col}" data-grid-row="{row.grid_row}">',
                  rect(x, y, DETAIL_TILE_WIDTH, DETAIL_TILE_HEIGHT, PALE[label], "#DCE3E8"),
                  rect(x, y, DETAIL_TILE_WIDTH, 4, COLORS[label]),
                  molecule(row.canonical_smiles, x + 8, y + 7),
                  text(x + 93, y + 169, row.molecule_identity_key, 15, bold=True, anchor="middle"),
                  text(x + 93, y + 190, f"{label.capitalize()}  |  pChEMBL {row.activity_pchembl:.2f}",
                       14, MUTED, anchor="middle"), '</g>']
    parts += [text(806, 833, "Structures, cell positions and median values are retained.", 18, MUTED),
              text(806, 860, "Median pChEMBL provides the activity annotation.", 17, MUTED)]

    for x, label, caption in [(36, "inactive", "Inactive: pChEMBL \u2264 6"),
                              (276, "medium", "Medium: 6 < pChEMBL < 7"),
                              (555, "active", "Active: pChEMBL \u2265 7")]:
        # The three legend items fit beneath the overview and its adjoining gap.
        parts += [rect(x, 870, 13, 13, COLORS[label]), text(x + 20, 882, caption, 16)]
    parts.append('</svg>')
    return "".join(parts), detail


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DATA)
    args = parser.parse_args()
    data = pd.read_csv(args.source)[FIELDS]
    DOCS.mkdir(exist_ok=True)
    if args.source.resolve() != DATA.resolve():
        data.to_csv(DATA, index=False)
        metadata = {"source": str(args.source.relative_to(ROOT)) if args.source.is_absolute()
                    else str(args.source), "source_sha256": hashlib.sha256(args.source.read_bytes()).hexdigest(),
                    "target": "CHEMBL205", "chembl_release": 37, "molecules": len(data),
                    "projection": "Morgan fingerprints / UMAP, archived seed 42",
                    "detail_window": {"min_col": WINDOW_COL, "min_row": WINDOW_ROW, "size": 3},
                    "selection": "Previously archived mixed window; illustrative, not a new result.",
                    "data_license": "ChEMBL attribution/share-alike terms; see validation_data/README.md"}
        (DOCS / "chemgridmap_demo_source.json").write_text(json.dumps(metadata, indent=2) + "\n")
    svg, detail = compose(data)
    svg_path = DOCS / "chemgridmap_demo.svg"
    svg_path.write_text(svg)
    with pymupdf.open(stream=svg.encode(), filetype="svg") as source:
        with pymupdf.open(stream=source.convert_to_pdf(), filetype="pdf") as pdf:
            pdf[0].get_pixmap(matrix=pymupdf.Matrix(3, 3), alpha=False).save(DOCS / "chemgridmap_demo.png")
    root = ET.fromstring(svg)
    assert not root.findall(f".//{{{NS}}}image"), "All structures must be native vector paths"
    assert len(root.findall(f".//{{{NS}}}g[@data-detail-molecule-id]")) == 9
    print(f"Rendered {len(data)} cells and {len(detail)} matching structures: {svg_path}")


if __name__ == "__main__":
    main()
