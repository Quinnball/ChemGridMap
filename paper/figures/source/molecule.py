"""Generate and embed vector RDKit molecules."""

from __future__ import annotations

import json
import os
from pathlib import Path
from tempfile import NamedTemporaryFile

import numpy as np
import pymupdf
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D


FONT_PATH = Path("/System/Library/Fonts/Supplemental/Arial.ttf")


def molecule_svg(smiles: str) -> tuple[str, dict[str, float]]:
    """Return SVG and printing measurements without changing the input molecule."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit cannot parse molecule: {smiles}")
    before = Chem.MolToSmiles(mol, isomericSmiles=True)
    rdDepictor.Compute2DCoords(mol)

    width_px, height_px = 900, 467
    drawer = rdMolDraw2D.MolDraw2DSVG(width_px, height_px)
    options = drawer.drawOptions()
    options.clearBackground = False
    options.fontFile = str(FONT_PATH)
    options.padding = 0.05
    options.updateAtomPalette({16: (0.60, 0.50, 0.00)})
    options.bondLineWidth = 8.25
    options.minFontSize = 76
    options.maxFontSize = 76
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    after = Chem.MolToSmiles(mol, isomericSmiles=True)
    if after != before:
        raise AssertionError("RDKit drawing changed canonical SMILES")

    coordinates = [drawer.GetDrawCoords(i) for i in range(mol.GetNumAtoms())]
    bond_lengths = []
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bond_lengths.append(
            ((coordinates[a].x - coordinates[b].x) ** 2
             + (coordinates[a].y - coordinates[b].y) ** 2) ** 0.5
        )

    canvas_mm = 27.0
    mm_per_px = canvas_mm / width_px
    measurements = {
        "median_bond_length_mm": float(round(float(np.median(bond_lengths)) * mm_per_px, 3)),
        "bond_line_width_pt": float(round(options.bondLineWidth * mm_per_px * 72 / 25.4, 3)),
        "atom_label_height_mm": float(round(options.maxFontSize * mm_per_px, 3)),
        "atom_label_height_pt": float(
            round(options.maxFontSize * mm_per_px * 72 / 25.4, 3)
        ),
    }
    return drawer.GetDrawingText(), measurements


def _coerce_molecule(value: str | Chem.Mol) -> Chem.Mol:
    """Accept SMILES or an RDKit molecule and return an independent copy."""
    if isinstance(value, Chem.Mol):
        return Chem.Mol(value)
    molecule = Chem.MolFromSmiles(value)
    if molecule is None:
        raise ValueError(f"RDKit cannot parse molecule: {value}")
    return molecule


def fixed_bond_molecule_svg(
    smiles: str | Chem.Mol,
    *,
    printed_bond_length_mm: float,
) -> tuple[str, dict[str, float]]:
    """Draw on an adaptive canvas and scale to a target printed bond length."""
    mol = _coerce_molecule(smiles)
    before = Chem.MolToSmiles(mol, isomericSmiles=True)
    highlight_atoms = []
    highlight_bonds = []
    if mol.HasProp("_f3_highlight_atoms"):
        highlight_atoms = json.loads(mol.GetProp("_f3_highlight_atoms"))
    if mol.HasProp("_f3_highlight_bonds"):
        highlight_bonds = json.loads(mol.GetProp("_f3_highlight_bonds"))
    highlight_color = (
        mol.GetProp("_f3_highlight_color")
        if mol.HasProp("_f3_highlight_color")
        else ""
    )
    rdDepictor.Compute2DCoords(mol)

    drawer = rdMolDraw2D.MolDraw2DSVG(-1, -1)
    options = drawer.drawOptions()
    options.clearBackground = False
    options.fontFile = str(FONT_PATH)
    options.padding = 0.05
    options.fixedBondLength = 100
    dummy_atoms = sum(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms())
    options.dummiesAreAttachments = dummy_atoms > 0
    options.updateAtomPalette({16: (0.60, 0.50, 0.00)})
    options.bondLineWidth = 3.3
    options.minFontSize = 28
    options.maxFontSize = 28
    if highlight_atoms:
        color = tuple(
            int(highlight_color[index:index + 2], 16) / 255
            for index in (1, 3, 5)
        )
        drawer.DrawMolecule(
            mol,
            highlightAtoms=highlight_atoms,
            highlightBonds=highlight_bonds,
            highlightAtomColors={index: color for index in highlight_atoms},
            highlightBondColors={index: color for index in highlight_bonds},
        )
    else:
        drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    after = Chem.MolToSmiles(mol, isomericSmiles=True)
    if after != before:
        raise AssertionError("RDKit drawing changed canonical SMILES")

    coordinates = [drawer.GetDrawCoords(i) for i in range(mol.GetNumAtoms())]
    bond_lengths_px = [
        (
            (coordinates[bond.GetBeginAtomIdx()].x - coordinates[bond.GetEndAtomIdx()].x) ** 2
            + (coordinates[bond.GetBeginAtomIdx()].y - coordinates[bond.GetEndAtomIdx()].y) ** 2
        )
        ** 0.5
        for bond in mol.GetBonds()
    ]
    median_bond_px = float(np.median(bond_lengths_px))
    mm_per_px = printed_bond_length_mm / median_bond_px
    width_mm = float(drawer.Width()) * mm_per_px
    height_mm = float(drawer.Height()) * mm_per_px
    measurements = {
        "heavy_atoms": float(mol.GetNumHeavyAtoms()),
        "dummy_atoms": float(dummy_atoms),
        "max_nitrogen_hydrogens": float(
            max(
                (
                    atom.GetTotalNumHs()
                    for atom in mol.GetAtoms()
                    if atom.GetAtomicNum() == 7
                ),
                default=0,
            )
        ),
        "canvas_width_px": float(drawer.Width()),
        "canvas_height_px": float(drawer.Height()),
        "median_bond_px": median_bond_px,
        "mm_per_px": mm_per_px,
        "width_mm": width_mm,
        "height_mm": height_mm,
        "printed_bond_length_mm": float(
            round(median_bond_px * mm_per_px, 6)
        ),
        "bond_line_width_mm": float(options.bondLineWidth * mm_per_px),
        "bond_line_width_pt": float(
            options.bondLineWidth * mm_per_px * 72 / 25.4
        ),
        "atom_label_height_mm": float(options.maxFontSize * mm_per_px),
        "atom_label_height_pt": float(
            options.maxFontSize * mm_per_px * 72 / 25.4
        ),
        "highlight_atom_count": float(len(highlight_atoms)),
        "highlight_bond_count": float(len(highlight_bonds)),
    }
    svg = drawer.GetDrawingText()
    if highlight_atoms:
        # Match highlight circle and fill colours; omit a sub-0.4 pt hairline.
        svg = svg.replace(f"stroke:{highlight_color};", "stroke:none;")
    return svg, measurements


def insert_molecule(pdf_path: Path, smiles: str, rect_mm_pt: tuple[float, float, float, float]) -> dict[str, float]:
    """Convert an RDKit SVG to vector PDF and place it in a specified millimetre rectangle."""
    svg, measurements = molecule_svg(smiles)
    with NamedTemporaryFile(suffix=".svg", delete=False) as stream:
        svg_path = Path(stream.name)
        stream.write(svg.encode("utf-8"))
    try:
        source_svg = pymupdf.open(svg_path, filetype="svg")
        source_pdf = pymupdf.open("pdf", source_svg.convert_to_pdf())
        pdf = pymupdf.open(pdf_path)
        page = pdf[0]
        x, y_top, width, height = rect_mm_pt
        points_per_mm = 72 / 25.4
        rect = pymupdf.Rect(
            x * points_per_mm,
            y_top * points_per_mm,
            (x + width) * points_per_mm,
            (y_top + height) * points_per_mm,
        )
        page.show_pdf_page(rect, source_pdf, 0, keep_proportion=True)
        with NamedTemporaryFile(suffix=".pdf", delete=False) as out:
            temp_path = Path(out.name)
        pdf.save(
            temp_path,
            garbage=4,
            deflate=True,
            no_new_id=True,
            reproducible=True,
        )
        pdf.close()
        source_pdf.close()
        source_svg.close()
        os.replace(temp_path, pdf_path)
    finally:
        svg_path.unlink(missing_ok=True)
    return measurements


def insert_fixed_bond_molecule(
    pdf_path: Path,
    smiles: str | Chem.Mol,
    *,
    center_x_mm: float,
    center_y_mm: float,
    canvas_height_mm: float,
    printed_bond_length_mm: float = 2.4,
) -> dict[str, float]:
    """Centre a fixed-bond-length RDKit vector molecule on the page."""
    svg, measurements = fixed_bond_molecule_svg(
        smiles,
        printed_bond_length_mm=printed_bond_length_mm,
    )
    with NamedTemporaryFile(suffix=".svg", delete=False) as stream:
        svg_path = Path(stream.name)
        stream.write(svg.encode("utf-8"))
    try:
        source_svg = pymupdf.open(svg_path, filetype="svg")
        source_pdf = pymupdf.open("pdf", source_svg.convert_to_pdf())
        pdf = pymupdf.open(pdf_path)
        page = pdf[0]
        points_per_mm = 72 / 25.4
        left_mm = center_x_mm - measurements["width_mm"] / 2
        top_mm = canvas_height_mm - center_y_mm - measurements["height_mm"] / 2
        rect = pymupdf.Rect(
            left_mm * points_per_mm,
            top_mm * points_per_mm,
            (left_mm + measurements["width_mm"]) * points_per_mm,
            (top_mm + measurements["height_mm"]) * points_per_mm,
        )
        page.show_pdf_page(rect, source_pdf, 0, keep_proportion=True)
        with NamedTemporaryFile(suffix=".pdf", delete=False) as out:
            temp_path = Path(out.name)
        pdf.save(
            temp_path,
            garbage=4,
            deflate=True,
            no_new_id=True,
            reproducible=True,
        )
        pdf.close()
        source_pdf.close()
        source_svg.close()
        os.replace(temp_path, pdf_path)
    finally:
        svg_path.unlink(missing_ok=True)
    return {
        **measurements,
        "center_x_mm": float(center_x_mm),
        "center_y_mm": float(center_y_mm),
    }
