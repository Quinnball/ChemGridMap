"""Figure 2: activity, support, and neighbourhood reliability layers for the CA II map."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon, Rectangle
from rdkit import Chem
from rdkit.Chem import rdFMCS
from chemgridmap.plotting import DEFAULT_COLORS, FADED_COLORS, RELIABILITY_BINS

import common as c


FIG_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = FIG_DIR / "data_fig2"
DEFAULT_DIR = DATA_DIR / "default"
FMCS_DIR = DATA_DIR / "fmcs"
RAW_INPUT = c.CGM_DIR / "validation_data/chembl37_chembl205/chembl37_chembl205_ic50_raw.csv"
FROZEN_COORDS = (
    c.CGM_DIR / "paper/output/reliability_v1/layouts/chembl205/coordinates_seed_42.csv"
)

EXPECTED_INPUT_SHA256 = (
    "e358cdc68574bd8ca238fe4fc955254ab434ebc26d1ffdac41ba644b90d18241"
)
OUTLINE = "#172536"
MAP_RECTS_MM = ((2.0, 58.0), (59.8, 58.0), (117.6, 58.0))
MAP_WIDTH_MM = 50.4
MAP_HEIGHT_MM = 60.0
SOURCE_HEIGHT_MM = 125.0
WIDE_MARKER_LEG = 0.62
TRIANGLE_INSET = 0.06
TRIANGLE_EDGEWIDTH = 0.3


def _tag_swatch(patch: Rectangle, color: str, panel: str) -> None:
    patch._qa_fig2_expected_fill = color.lower()  # type: ignore[attr-defined]
    patch._qa_fig2_color_panel = panel  # type: ignore[attr-defined]


def _ordinary_query_graph(pattern):
    """Build an ordinary molecular graph without implicit hydrogens from query atom and bond types."""
    molecule = Chem.RWMol()
    for atom in pattern.GetAtoms():
        ordinary = Chem.Atom(atom.GetAtomicNum())
        ordinary.SetNoImplicit(True)
        ordinary.SetIsAromatic(
            any(bond.GetIsAromatic() for bond in atom.GetBonds())
        )
        molecule.AddAtom(ordinary)
    for bond in pattern.GetBonds():
        molecule.AddBond(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            bond.GetBondType(),
        )
    result = molecule.GetMol()
    Chem.SanitizeMol(result)
    return result


def _fmcs_core_with_attachments(
    pattern,
    members: pd.DataFrame,
) -> tuple[Chem.Mol, dict[str, object]]:
    """Copy core hydrogen states from the first matching member and draw attachment points outside rings."""
    template = None
    member_match = None
    member_id = None
    for row in members.itertuples(index=False):
        molecule = Chem.MolFromSmiles(row.canonical_smiles)
        if molecule is None:
            continue
        match = molecule.GetSubstructMatch(pattern)
        if match:
            template = molecule
            member_match = match
            member_id = row.molecule_identity_key
            break
    if template is None or member_match is None:
        raise AssertionError("No G001 member matches the FMCS core")

    drawing = Chem.RWMol()
    for pattern_atom, member_idx in zip(
        pattern.GetAtoms(),
        member_match,
        strict=True,
    ):
        member_atom = template.GetAtomWithIdx(member_idx)
        atom = Chem.Atom(pattern_atom.GetAtomicNum())
        atom.SetFormalCharge(member_atom.GetFormalCharge())
        atom.SetIsAromatic(
            any(
                bond.GetBondType() == Chem.BondType.AROMATIC
                for bond in pattern_atom.GetBonds()
            )
        )
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(member_atom.GetTotalNumHs())
        drawing.AddAtom(atom)

    for bond in pattern.GetBonds():
        drawing.AddBond(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            bond.GetBondType(),
        )
        drawing_bond = drawing.GetBondBetweenAtoms(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
        )
        drawing_bond.SetIsAromatic(
            bond.GetBondType() == Chem.BondType.AROMATIC
        )

    matched_member_atoms = set(member_match)
    dummy_count = 0
    for core_idx, member_idx in enumerate(member_match):
        member_atom = template.GetAtomWithIdx(member_idx)
        for member_bond in member_atom.GetBonds():
            if member_bond.GetOtherAtomIdx(member_idx) in matched_member_atoms:
                continue
            dummy_idx = drawing.AddAtom(Chem.Atom(0))
            drawing.GetAtomWithIdx(dummy_idx).SetNoImplicit(True)
            drawing.AddBond(
                core_idx,
                dummy_idx,
                member_bond.GetBondType(),
            )
            dummy_count += 1

    molecule = drawing.GetMol()
    Chem.SanitizeMol(molecule)
    amide_n_idx = pattern.GetNumAtoms() - 1
    amide_n_hydrogens = molecule.GetAtomWithIdx(amide_n_idx).GetTotalNumHs()
    heavy_atoms = molecule.GetNumHeavyAtoms()
    if heavy_atoms != 22:
        raise AssertionError(f"Drawn FMCS core has heavy atoms {heavy_atoms} instead of 22")
    if amide_n_hydrogens != 0:
        raise AssertionError(f"FMCS amide N hydrogens {amide_n_hydrogens} instead of 0")
    encoded = Chem.MolToSmiles(molecule, isomericSmiles=True)
    if Chem.MolFromSmiles(encoded) is None:
        raise AssertionError("FMCS core with attachment points cannot be serialized as SMILES")
    return molecule, {
        "template_member": member_id,
        "heavy_atoms": heavy_atoms,
        "dummy_atoms": dummy_count,
        "amide_n_index": amide_n_idx,
        "amide_n_hydrogens": amide_n_hydrogens,
    }


def _add_triangle(ax, col: int, row: int) -> None:
    leg = WIDE_MARKER_LEG
    right = col + 0.5 - TRIANGLE_INSET
    top = row + 0.5 - TRIANGLE_INSET
    marker = Polygon(
        (
            (right, top),
            (right - leg, top),
            (right, top - leg),
        ),
        closed=True,
        facecolor=OUTLINE,
        edgecolor="white",
        linewidth=TRIANGLE_EDGEWIDTH,
        zorder=4,
    )
    marker.set_gid("wide-record-marker")
    marker._qa_fig2_triangle_leg = leg  # type: ignore[attr-defined]
    marker._qa_fig2_triangle_vertices = marker.get_xy().copy()  # type: ignore[attr-defined]
    marker._qa_fig2_triangle_edgewidth = TRIANGLE_EDGEWIDTH  # type: ignore[attr-defined]
    marker._qa_fig2_triangle_cell = (col, row)  # type: ignore[attr-defined]
    ax.add_patch(marker)


def _add_group_outline(ax, members: pd.DataFrame) -> None:
    occupied = {
        (int(row.grid_col), int(row.grid_row))
        for row in members.itertuples(index=False)
    }
    segments = []
    for col, row in sorted(occupied):
        if (col, row - 1) not in occupied:
            segments.append(((col - 0.5, row - 0.5), (col + 0.5, row - 0.5)))
        if (col, row + 1) not in occupied:
            segments.append(((col - 0.5, row + 0.5), (col + 0.5, row + 0.5)))
        if (col - 1, row) not in occupied:
            segments.append(((col - 0.5, row - 0.5), (col - 0.5, row + 0.5)))
        if (col + 1, row) not in occupied:
            segments.append(((col + 0.5, row - 0.5), (col + 0.5, row + 0.5)))
    ax.add_collection(
        LineCollection(
            segments,
            colors=OUTLINE,
            linewidths=0.8,
            zorder=5,
            capstyle="butt",
            joinstyle="miter",
        )
    )


def _legend_axis(fig):
    ax = c.add_axes_mm(fig, (0.0, 41.5, 170.0, 15.5), (170.0, SOURCE_HEIGHT_MM))
    ax.set_xlim(0.0, 170.0)
    ax.set_ylim(41.5, 57.0)
    ax.set_axis_off()
    return ax


def _legend_text(ax, x: float, y: float, text: str, **kwargs) -> None:
    c.mark_legend_text(
        ax.text(x, y, text, fontsize=6.5, color=c.INK, va="center", **kwargs)
    )


def _build_legends(fig, data: pd.DataFrame):
    ax = _legend_axis(fig)

    activity_rows = (
        ("Inactive (≤6)", "inactive", "a"),
        ("Medium (6–7)", "medium", "b"),
        ("Active (≥7)", "active", "c"),
    )
    for index, (label, key, marker) in enumerate(activity_rows):
        y = 53.3 - index * 3.1
        patch = Rectangle(
            (2.0, y - 1.25),
            2.5,
            2.5,
            facecolor=DEFAULT_COLORS[key],
            edgecolor="white",
            linewidth=0.4,
        )
        _tag_swatch(patch, DEFAULT_COLORS[key], "a")
        ax.add_patch(patch)
        _legend_text(ax, 5.3, y, label)
    y = 44.0
    legend_triangle = Polygon(
        ((3.25, y + 0.6), (2.0, y + 0.6), (3.25, y - 0.65)),
        closed=True,
        facecolor=OUTLINE,
        edgecolor="white",
        linewidth=TRIANGLE_EDGEWIDTH,
    )
    legend_triangle._qa_fig2_legend_triangle = True  # type: ignore[attr-defined]
    legend_triangle._qa_fig2_triangle_vertices = (  # type: ignore[attr-defined]
        legend_triangle.get_xy().copy()
    )
    legend_triangle._qa_fig2_triangle_edgewidth = TRIANGLE_EDGEWIDTH  # type: ignore[attr-defined]
    ax.add_patch(legend_triangle)
    _legend_text(ax, 5.3, y, "Wide record range")

    _legend_text(ax, 66.3, 55.5, "Supported", fontweight="bold", ha="center")
    _legend_text(ax, 84.5, 55.5, "Uncertain", fontweight="bold", ha="center")
    for index, label in enumerate(("inactive", "medium", "active")):
        y = 52.0 - index * 2.68
        _legend_text(ax, 63.6, y, label.capitalize(), ha="right")
        for x, source in ((65.1, DEFAULT_COLORS), (83.3, FADED_COLORS)):
            if label == "medium" and source is DEFAULT_COLORS:
                _legend_text(ax, x + 1.15, y, "—", ha="center")
            else:
                ax.add_patch(
                    Rectangle(
                        (x, y - 1.0),
                        2.3,
                        2.0,
                        facecolor=source[label],
                        edgecolor="white",
                        linewidth=0.4,
                    )
                )
                _tag_swatch(ax.patches[-1], source[label], "b")
    supported = int((data["class_status"] == "supported").sum())
    _legend_text(
        ax,
        63.6,
        43.8,
        f"Supported: {supported} of {len(data)}",
    )

    reliability_labels = ("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1")
    for index, (label, (_, color)) in enumerate(
        zip(reliability_labels, RELIABILITY_BINS, strict=True)
    ):
        x = 118.6 + (index % 3) * 16.0
        y = 53.0 - (index // 3) * 4.2
        patch = Rectangle(
            (x, y - 1.0),
            2.3,
            2.0,
            facecolor=color,
            edgecolor="white",
            linewidth=0.4,
        )
        patch._qa_fig2c_swatch = True  # type: ignore[attr-defined]
        _tag_swatch(patch, color, "c")
        ax.add_patch(patch)
        _legend_text(ax, x + 3.0, y, label)
    mean_recall = float(data["neighbor_recall"].mean())
    _legend_text(ax, 118.6, 45.0, f"Mean {mean_recall:.3f}")
    fixed_k_recall = float(data["neighbor_recall_fixed_k"].mean())
    _legend_text(ax, 118.6, 42.0, f"Fixed k {fixed_k_recall:.3f}")


def _validate_and_load() -> dict[str, object]:
    input_sha = hashlib.sha256(RAW_INPUT.read_bytes()).hexdigest()
    if input_sha != EXPECTED_INPUT_SHA256:
        raise AssertionError(f"Input SHA-256 mismatch: {input_sha}")

    data = pd.read_csv(
        DEFAULT_DIR / "caii_seed42_grid_coordinates.csv",
        float_precision="round_trip",
    )
    fmcs_data = pd.read_csv(
        FMCS_DIR / "caii_seed42_grid_coordinates.csv",
        float_precision="round_trip",
    )
    frozen = pd.read_csv(FROZEN_COORDS, float_precision="round_trip")
    default_groups = json.loads(
        (DEFAULT_DIR / "caii_seed42_structure_groups.json").read_text()
    )
    fmcs_groups = json.loads(
        (FMCS_DIR / "caii_seed42_structure_groups.json").read_text()
    )
    default_manifest = json.loads(
        (DEFAULT_DIR / "caii_seed42_run_manifest.json").read_text()
    )
    fmcs_manifest = json.loads(
        (FMCS_DIR / "caii_seed42_run_manifest.json").read_text()
    )
    members = pd.read_csv(
        DEFAULT_DIR / "caii_seed42_structure_members.csv",
        float_precision="round_trip",
    )

    checks: list[dict[str, object]] = []

    def check(name: str, expected, actual) -> None:
        checks.append(
            {
                "name": name,
                "expected": expected,
                "actual": actual,
                "passed": expected == actual,
            }
        )

    check("输入 SHA-256", EXPECTED_INPUT_SHA256, input_sha)
    check("默认运行 software_version", "0.5.0", default_manifest["software_version"])
    check("FMCS 运行 software_version", "0.5.0", fmcs_manifest["software_version"])
    check("默认运行条目数", 889, len(data))
    check("FMCS 运行条目数", 889, len(fmcs_data))
    frozen_indexed = frozen.set_index("molecule_identity_key")
    coordinate_matches = all(
        (
            int(frozen_indexed.loc[row.molecule_identity_key, "refined_grid_col"]),
            int(frozen_indexed.loc[row.molecule_identity_key, "refined_grid_row"]),
        )
        == (int(row.grid_col), int(row.grid_row))
        for row in data.itertuples(index=False)
    )
    check("默认坐标逐条匹配冻结 refined 坐标", True, coordinate_matches)
    check("supported 条目数", 506, int((data["class_status"] == "supported").sum()))
    check("平均 neighbor recall（三位）", 0.585, round(float(data["neighbor_recall"].mean()), 3))
    check("冲突条目数", 17, int((data["is_conflicted"] == 1).sum()))
    check("宽记录范围三角直角边（格子）", 0.62, WIDE_MARKER_LEG)

    default_by_id = {group["structure_group"]: group for group in default_groups["groups"]}
    fmcs_by_id = {group["structure_group"]: group for group in fmcs_groups["groups"]}
    g001 = default_by_id["G001"]
    g002 = default_by_id["G002"]
    g001_fmcs = fmcs_by_id["G001"]
    check("G001 成员数", 42, g001["size"])
    check("G001 核心重原子数", 12, g001["core_heavy_atoms"])
    check("G001 核心覆盖数", 42, g001["core_coverage_count"])
    check("G001 最常见骨架重原子数", 31, Chem.MolFromSmiles(g001["most_frequent_scaffold"]).GetNumHeavyAtoms())
    check("G001 最常见骨架覆盖数", 5, g001["most_frequent_scaffold_coverage_count"])
    check("G002 成员数", 28, g002["size"])
    check("G002 核心重原子数", 9, g002["core_heavy_atoms"])
    check("G002 核心覆盖数", 27, g002["core_coverage_count"])
    check("G002 exact scaffold 数", 17, g002["core_exact_scaffold_count"])
    check("G001 FMCS 核心重原子数", 22, g001_fmcs["core_heavy_atoms"])
    check("G001 FMCS 核心覆盖数", 39, g001_fmcs["core_coverage_count"])

    default_assignment = dict(zip(data["molecule_identity_key"], data["structure_group"], strict=True))
    fmcs_assignment = dict(
        zip(fmcs_data["molecule_identity_key"], fmcs_data["structure_group"], strict=True)
    )
    check("两次运行分组相同", True, default_assignment == fmcs_assignment)

    core_mol = Chem.MolFromSmiles(g001_fmcs["core_smiles"])
    pattern = Chem.MolFromSmarts(g001_fmcs["core_smarts"])
    if pattern is None:
        raise AssertionError("G001 FMCS core_smarts cannot be parsed")
    fmcs_core, fmcs_core_info = _fmcs_core_with_attachments(
        pattern,
        fmcs_data[fmcs_data["structure_group"] == "G001"],
    )
    check("G001 FMCS core_smiles 可解析", True, core_mol is not None)
    check("G001 FMCS core_smiles 重原子数", 22, core_mol.GetNumHeavyAtoms() if core_mol else None)
    check("G001 FMCS core_smarts 可解析", True, pattern is not None)
    member_pattern_matches = int(
        fmcs_data.loc[fmcs_data["structure_group"] == "G001", "canonical_smiles"]
        .map(lambda smiles: Chem.MolFromSmiles(smiles).HasSubstructMatch(pattern))
        .sum()
    )
    check("G001 FMCS SMARTS 匹配成员数", 39, member_pattern_matches)
    query_graph = _ordinary_query_graph(pattern)
    graph_mcs = rdFMCS.FindMCS(
        [core_mol, query_graph],
        ringMatchesRingOnly=True,
        completeRingsOnly=True,
    )
    graph_signature_match = (
        core_mol.GetNumHeavyAtoms() == pattern.GetNumAtoms()
        and graph_mcs.numAtoms == core_mol.GetNumHeavyAtoms()
        and member_pattern_matches == g001_fmcs["core_coverage_count"]
    )
    check("G001 FMCS core_smiles 与 SMARTS 普通图 MCS 原子数", 22, graph_mcs.numAtoms)
    check("G001 FMCS 绘制结构与 SMARTS 普通图一致", True, graph_signature_match)
    check("G001 FMCS 绘制结构重原子数", 22, fmcs_core.GetNumHeavyAtoms())
    check("G001 FMCS 绘制结构虚原子数", 2, fmcs_core_info["dummy_atoms"])
    check("G001 FMCS 绘制结构酰胺 N 氢数", 0, fmcs_core_info["amide_n_hydrogens"])

    g002_members = members[members["structure_group"] == "G002"]
    g002_member_ids = set(g002_members["molecule_identity_key"])
    core_member_ids = set(g002["core_member_ids"])
    noncore_ids = sorted(g002_member_ids - core_member_ids)
    representative_id = str(g002["representative_id"])
    representative_row = g002_members[
        g002_members["molecule_identity_key"] == representative_id
    ].iloc[0]
    noncore_row = g002_members[
        g002_members["molecule_identity_key"] == "CHEMBL85636"
    ].iloc[0]
    g002_core = Chem.MolFromSmiles(g002["core_smiles"])
    representative = Chem.MolFromSmiles(representative_row["canonical_smiles"])
    noncore = Chem.MolFromSmiles(noncore_row["canonical_smiles"])
    if g002_core is None or representative is None or noncore is None:
        raise AssertionError("G002 core or representative SMILES cannot be parsed")
    representative_matches = representative.GetSubstructMatches(g002_core)
    if len(representative_matches) != 1:
        raise AssertionError(
            f"CHEMBL81688 core match count is {len(representative_matches)} instead of 1"
        )
    representative_match = representative_matches[0]
    highlighted_bonds = []
    for query_bond in g002_core.GetBonds():
        first = representative_match[query_bond.GetBeginAtomIdx()]
        second = representative_match[query_bond.GetEndAtomIdx()]
        actual_bond = representative.GetBondBetweenAtoms(first, second)
        if actual_bond is None:
            raise AssertionError("CHEMBL81688 core match is missing a corresponding bond")
        highlighted_bonds.append(actual_bond.GetIdx())
    representative.SetProp(
        "_f3_highlight_atoms",
        json.dumps(list(representative_match)),
    )
    representative.SetProp(
        "_f3_highlight_bonds",
        json.dumps(highlighted_bonds),
    )
    representative.SetProp("_f3_highlight_color", "#CFE3F3")

    check("G002 representative_id", "CHEMBL81688", representative_id)
    check("G002 非核心成员", ["CHEMBL85636"], noncore_ids)
    check("CHEMBL81688 含 G002 核心", True, representative.HasSubstructMatch(g002_core))
    check("CHEMBL85636 不含 G002 核心", False, noncore.HasSubstructMatch(g002_core))
    check(
        "G002 代表分子 SMILES 取自 structure_members.csv",
        representative_row["canonical_smiles"],
        g002["representative_smiles"],
    )
    check("CHEMBL81688 高亮核心原子数", 9, len(representative_match))
    check(
        "CHEMBL81688 高亮核心键数",
        g002_core.GetNumBonds(),
        len(highlighted_bonds),
    )

    failed = [item for item in checks if not item["passed"]]
    if failed:
        raise AssertionError(f"Figure 2 data assertions failed: {failed}")

    return {
        "data": data,
        "fmcs_data": fmcs_data,
        "g001": g001,
        "g002": g002,
        "g001_fmcs": g001_fmcs,
        "fmcs_core": fmcs_core,
        "fmcs_core_info": fmcs_core_info,
        "g002_representative": representative,
        "g002_noncore_smiles": noncore_row["canonical_smiles"],
        "assertions": checks,
    }


def _write_plotted_values(data: pd.DataFrame) -> None:
    path = c.QA_DIR / "Figure_2_plotted_values.csv"
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
        source_file = "caii_seed42_grid_coordinates.csv"
        for row in data.itertuples(index=False):
            common_row = {
                "figure": 2,
                "row_identifier": row.molecule_identity_key,
                "source_file": source_file,
            }
            writer.writerow(
                {
                    **common_row,
                    "panel": "a",
                    "series": row.activity_class,
                    "variable": "activity_pchembl",
                    "value": float(row.activity_pchembl),
                }
            )
            writer.writerow(
                {
                    **common_row,
                    "panel": "b",
                    "series": row.class_status,
                    "variable": "class_confidence",
                    "value": float(row.class_confidence),
                }
            )
            writer.writerow(
                {
                    **common_row,
                    "panel": "c",
                    "series": "neighbor_recall",
                    "variable": "neighbor_recall",
                    "value": float(row.neighbor_recall),
                }
            )
        writer.writerow(
            {
                "figure": 2,
                "panel": "c",
                "series": "summary",
                "row_identifier": "mean",
                "variable": "neighbor_recall_mean",
                "source_file": source_file,
                "value": float(data["neighbor_recall"].mean()),
            }
        )
        writer.writerow(
            {
                "figure": 2,
                "panel": "c",
                "series": "summary",
                "row_identifier": "mean",
                "variable": "neighbor_recall_fixed_k_mean",
                "source_file": source_file,
                "value": float(data["neighbor_recall_fixed_k"].mean()),
            }
        )


def build():
    loaded = _validate_and_load()
    data: pd.DataFrame = loaded["data"]
    c.init()
    _write_plotted_values(data)
    fig = c.mm_figure(170.0, SOURCE_HEIGHT_MM)

    map_axes = []
    for x0, y0 in MAP_RECTS_MM:
        ax = c.add_axes_mm(
            fig,
            (x0, y0, MAP_WIDTH_MM, MAP_HEIGHT_MM),
            (170.0, SOURCE_HEIGHT_MM),
        )
        ax.set_xlim(1.5, 43.5)
        ax.set_ylim(-1.5, 48.5)
        ax.set_aspect("equal", adjustable="box")
        ax.set_axis_off()
        map_axes.append(ax)

    for panel_index, ax in enumerate(map_axes):
        panel_name = "abc"[panel_index]
        for row in data.itertuples(index=False):
            if panel_index == 0:
                if row.activity_class == "inactive":
                    color = DEFAULT_COLORS["inactive"]
                elif row.activity_class == "medium":
                    color = DEFAULT_COLORS["medium"]
                else:
                    color = DEFAULT_COLORS["active"]
            elif panel_index == 1:
                if row.class_status == "uncertain":
                    color = FADED_COLORS[row.activity_class]
                else:
                    color = DEFAULT_COLORS[row.activity_class]
            else:
                color = next(
                    bin_color
                    for low, bin_color in reversed(RELIABILITY_BINS)
                    if float(row.neighbor_recall) >= low - 1e-9
                )
            ax.add_patch(
                Rectangle(
                    (row.grid_col - 0.5, row.grid_row - 0.5),
                    1.0,
                    1.0,
                    facecolor=color,
                    edgecolor="white",
                    linewidth=0.4,
                    zorder=1,
                )
            )
            tile = ax.patches[-1]
            tile._qa_fig2_grid_panel = panel_name  # type: ignore[attr-defined]
            tile._qa_fig2_grid_row = int(row.grid_row)  # type: ignore[attr-defined]
            tile._qa_fig2_grid_col = int(row.grid_col)  # type: ignore[attr-defined]
            tile._qa_fig2_expected_fill = color.lower()  # type: ignore[attr-defined]

        if panel_index < 2:
            for row in data[data["is_conflicted"] == 1].itertuples(index=False):
                _add_triangle(ax, int(row.grid_col), int(row.grid_row))
                ax.patches[-1]._qa_fig2_triangle_panel = panel_name  # type: ignore[attr-defined]

    for group_id in ("G001", "G002"):
        members = data[data["structure_group"] == group_id]
        for ax in map_axes:
            _add_group_outline(ax, members)

    map_axes[0].text(
        2.3,
        36.6,
        "G002",
        fontsize=6.5,
        color=c.INK,
        ha="left",
        va="bottom",
        zorder=7,
    )
    map_axes[0].plot(
        (6.2, 6.8),
        (35.6, 33.3),
        color=OUTLINE,
        linewidth=0.4,
        solid_capstyle="butt",
        zorder=6,
    )
    map_axes[0].text(
        24.3,
        2.2,
        "G001",
        fontsize=6.5,
        color=c.INK,
        ha="left",
        va="center",
        zorder=7,
    )
    map_axes[0].plot(
        (22.55, 24.0),
        (2.2, 2.2),
        color=OUTLINE,
        linewidth=0.4,
        solid_capstyle="butt",
        zorder=6,
    )
    for ax in map_axes[1:]:
        ax.text(
            7.0,
            36.7,
            "G002",
            fontsize=6.5,
            color=c.INK,
            ha="left",
            va="bottom",
            zorder=7,
        )
        ax.plot(
            (7.5, 7.5),
            (36.1, 33.0),
            color=OUTLINE,
            linewidth=0.4,
            solid_capstyle="butt",
            zorder=6,
        )
        ax.text(
            24.3,
            2.2,
            "G001",
            fontsize=6.5,
            color=c.INK,
            ha="left",
            va="center",
            zorder=7,
        )
        ax.plot(
            (22.55, 24.0),
            (2.2, 2.2),
            color=OUTLINE,
            linewidth=0.4,
            solid_capstyle="butt",
            zorder=6,
        )

    titles = (
        (2.0, "Activity class"),
        (59.8, "Class support"),
        (117.6, "Neighbor recall (k = 10)"),
    )
    for x, title in titles:
        c.mark_panel_title(fig.text(
            x / 170.0,
            118.1 / SOURCE_HEIGHT_MM,
            title,
            fontsize=7,
            color=c.INK,
            ha="left",
            va="bottom",
        ))

    c.panel_letter(fig, "a", 1.0, 124.0, (170.0, SOURCE_HEIGHT_MM))
    c.panel_letter(fig, "b", 58.8, 124.0, (170.0, SOURCE_HEIGHT_MM))
    c.panel_letter(fig, "c", 116.6, 124.0, (170.0, SOURCE_HEIGHT_MM))
    c.panel_letter(fig, "d", 1.0, 41.0, (170.0, SOURCE_HEIGHT_MM))

    _build_legends(fig, data)
    fig.add_artist(
        Line2D(
            (94.5 / 170.0, 94.5 / 170.0),
            (4.0 / SOURCE_HEIGHT_MM, 41.0 / SOURCE_HEIGHT_MM),
            transform=fig.transFigure,
            color="#BFBFBF",
            linewidth=0.4,
            solid_capstyle="butt",
        )
    )

    g001 = loaded["g001"]
    g002 = loaded["g002"]
    g001_fmcs = loaded["g001_fmcs"]
    fmcs_core = loaded["fmcs_core"]
    c.mark_panel_title(fig.text(
        4.5 / 170.0,
        37.5 / SOURCE_HEIGHT_MM,
        f"G001 ({g001['size']} members)",
        fontsize=7,
        ha="left",
        va="baseline",
        color=c.INK,
    ))
    c.mark_panel_title(fig.text(
        98.0 / 170.0,
        37.5 / SOURCE_HEIGHT_MM,
        f"G002 ({g002['size']} members)",
        fontsize=7,
        ha="left",
        va="baseline",
        color=c.INK,
    ))

    structure_specs = (
        {
            "x": 17.0,
            "smiles": g001["most_frequent_scaffold"],
            "lines": (
                "Most frequent scaffold",
                f"{Chem.MolFromSmiles(g001['most_frequent_scaffold']).GetNumHeavyAtoms()} heavy atoms",
                f"in {g001['most_frequent_scaffold_coverage_count']} of {g001['size']} members",
            ),
            "slot": (2.0, 32.0),
        },
        {
            "x": 47.0,
            "smiles": g001["core_smiles"],
            "lines": (
                "Coverage-target core",
                f"{g001['core_heavy_atoms']} heavy atoms",
                f"in {g001['core_coverage_count']} of {g001['size']} members",
            ),
            "slot": (32.0, 62.0),
        },
        {
            "x": 77.0,
            "smiles": Chem.MolToSmiles(fmcs_core, isomericSmiles=True),
            "lines": (
                "FMCS core (10-s budget)",
                f"{g001_fmcs['core_heavy_atoms']} heavy atoms",
                f"in {g001_fmcs['core_coverage_count']} of {g001['size']} members",
                "wavy bonds: attachment points",
            ),
            "slot": (62.0, 92.0),
        },
        {
            "x": 108.5,
            "smiles": g002["core_smiles"],
            "lines": (
                "Coverage-target core",
                f"{g002['core_heavy_atoms']} heavy atoms",
                f"in {g002['core_coverage_count']} of {g002['size']} members",
            ),
            "slot": (96.0, 121.0),
        },
        {
            "x": 133.5,
            "smiles": loaded["g002_representative"],
            "lines": (
                "CHEMBL81688",
                "group center",
                "contains the core",
            ),
            "slot": (121.0, 146.0),
        },
        {
            "x": 158.5,
            "smiles": loaded["g002_noncore_smiles"],
            "lines": (
                "CHEMBL85636",
                "benzothiophene",
                "lacks the core",
            ),
            "slot": (146.0, 171.0),
        },
    )
    text_baselines_mm = (12.0, 9.0, 6.0, 3.0)
    for block_index, spec in enumerate(structure_specs):
        x = float(spec["x"])
        slot_x0, slot_x1 = spec["slot"]
        for line_index, text in enumerate(spec["lines"]):
            baseline = text_baselines_mm[line_index]
            label = fig.text(
                x / 170.0,
                baseline / SOURCE_HEIGHT_MM,
                text,
                fontsize=6.5,
                ha="center",
                va="baseline",
                color=c.INK,
            )
            label._qa_description_block = block_index  # type: ignore[attr-defined]
            label._qa_description_line = line_index  # type: ignore[attr-defined]
            label._qa_description_slot_mm = (slot_x0, slot_x1)  # type: ignore[attr-defined]
            if block_index == 2 and line_index == 3:
                label._qa_fig2_fmcs_wavy = True  # type: ignore[attr-defined]
    placements = [
        {**spec, "center_y": 25.0}
        for spec in structure_specs
    ]
    return fig, placements, loaded["assertions"]


def validate_supplementary(measurements: list[dict[str, float]]) -> dict[str, object]:
    expected = [31.0, 12.0, 22.0, 9.0, 17.0, 17.0]
    expected_dummy_atoms = [0, 0, 2, 0, 0, 0]
    heavy_atoms = [int(item["heavy_atoms"]) for item in measurements]
    dummy_atoms = [int(item["dummy_atoms"]) for item in measurements]
    bonds = [float(item["printed_bond_length_mm"]) for item in measurements]
    deviations = [abs(value - 2.4) / 2.4 for value in bonds]
    passed = (
        heavy_atoms == expected
        and dummy_atoms == expected_dummy_atoms
        and int(measurements[2]["max_nitrogen_hydrogens"]) == 0
        and len(measurements) == 6
        and max(deviations, default=1.0) <= 0.03
    )
    return {
        "name": "图 2 结构重原子数与统一键长",
        "passed": passed,
        "heavy_atoms": heavy_atoms,
        "dummy_atoms": dummy_atoms,
        "max_nitrogen_hydrogens": [
            int(item["max_nitrogen_hydrogens"]) for item in measurements
        ],
        "printed_bond_length_mm": bonds,
        "maximum_fractional_deviation": max(deviations, default=None),
        "measurements": measurements,
    }
