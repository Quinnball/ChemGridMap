"""Structure-defined groups and core-aligned comparisons, independent of layout."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import html
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, Geometry, rdBase
from rdkit.Chem import rdDepictor, rdFingerprintGenerator
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Cluster import Butina
from scipy.spatial.distance import cdist

from .cores import select_core
from .plotting import DEFAULT_COLORS, render_svg

GROUP_COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#8B6914", "#56B4E9",
                "#5B536E", "#B26249", "#658235", "#3B6C70", "#A95479", "#777777"]


def molecule_keys(data):
    for field in ("molecule_identity_key", "molecule_id", "name"):
        if field in data and data[field].notna().all() and data[field].astype(str).is_unique:
            return data[field].astype(str).tolist()
    if data.canonical_smiles.astype(str).is_unique:
        return data.canonical_smiles.astype(str).tolist()
    raise ValueError("Structure comparison needs unique molecular identifiers or unique structures.")


def fingerprint_distances(smiles, radius=2, n_bits=2048):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    mols = [Chem.MolFromSmiles(str(s)) for s in smiles]
    if any(m is None for m in mols):
        raise ValueError("Structure comparison encountered an invalid SMILES.")
    fps = [generator.GetFingerprint(m) for m in mols]
    distances = np.zeros((len(mols), len(mols)), dtype=np.float64)
    for i in range(1, len(mols)):
        d = 1 - np.asarray(DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i]))
        distances[i, :i] = d
        distances[:i, i] = d
    return mols, distances


def structural_groups(data, *, similarity=0.6, cluster_col=None, cluster_source=None,
                      radius=2, n_bits=2048, core_rule="most_frequent", core_min_coverage=0.8):
    """Butina groups use a center threshold, not an all-pairs similarity bound.

    ``core_rule="coverage"`` selects the largest member scaffold or scaffold-network
    ring fragment that matches at least ``core_min_coverage`` of the members.
    """
    if not 0 < similarity <= 1:
        raise ValueError("Group similarity must be in (0, 1].")
    if core_rule not in {"most_frequent", "coverage"}:
        raise ValueError("core_rule must be 'most_frequent' or 'coverage'.")
    if not 0 < core_min_coverage <= 1:
        raise ValueError("core_min_coverage must be in (0, 1].")
    if len(data) > 10000:
        raise ValueError("Exact grouping is limited to 10,000 molecules; compare an explicit subset.")
    frame = data.copy().reset_index(drop=True)
    keys = molecule_keys(frame)
    order = sorted(range(len(frame)), key=lambda i: (str(frame.canonical_smiles.iloc[i]), keys[i]))
    ordered = frame.iloc[order]
    mols, distances = fingerprint_distances(ordered.canonical_smiles, radius, n_bits)
    if cluster_col:
        if cluster_col not in frame or frame[cluster_col].isna().any() or not cluster_source:
            raise ValueError("External groups need a complete cluster column and an explicit source.")
        labels = ordered[cluster_col].astype(str).to_numpy()
        clusters = []
        for label in sorted(set(labels)):
            indices = np.flatnonzero(labels == label).tolist()
            center = indices[int(np.argmin(distances[np.ix_(indices, indices)].sum(axis=1)))]
            clusters.append((label, [center] + [i for i in indices if i != center]))
    else:
        found = Butina.ClusterData(distances, len(mols), 1 - similarity,
                                   isDistData=True, reordering=True)
        found = sorted(found, key=lambda c: (-len(c), str(ordered.canonical_smiles.iloc[c[0]]), keys[order[c[0]]]))
        clusters = [(f"G{i + 1:03d}", list(c)) for i, c in enumerate(found)]
    frame["comparison_id"] = keys
    frame["structure_group"] = ""
    frame["structure_group_color"] = ""
    frame["scaffold_smiles"] = ""
    summaries = []
    for rank, (label, indices) in enumerate(clusters):
        selected = [mols[i] for i in indices]
        scaffolds = [MurckoScaffold.MurckoScaffoldSmiles(mol=m, includeChirality=True) for m in selected]
        counts = Counter(s for s in scaffolds if s)
        ranked = sorted(counts, key=lambda s: (-counts[s], -Chem.MolFromSmiles(s).GetNumHeavyAtoms(), s))
        most_frequent = ranked[0] if ranked else ""
        frequent_pattern = Chem.MolFromSmiles(most_frequent) if most_frequent else None
        frequent_matched = [bool(frequent_pattern and m.HasSubstructMatch(frequent_pattern, useChirality=True))
                            for m in selected]
        if core_rule == "coverage" and most_frequent:
            choice = select_core(selected, min_coverage=core_min_coverage)
            core, matched = choice["core_smiles"], choice["matched"]
        else:
            choice, core, matched = None, most_frequent, frequent_matched
        covered = sum(matched)
        original_positions = [order[i] for i in indices]
        color = GROUP_COLORS[rank % len(GROUP_COLORS)]
        frame.loc[original_positions, "structure_group"] = label
        frame.loc[original_positions, "structure_group_color"] = color
        frame.loc[original_positions, "scaffold_smiles"] = scaffolds
        summaries.append({
            "structure_group": label, "size": len(indices), "color": color,
            "representative_id": keys[order[indices[0]]],
            "representative_smiles": str(ordered.canonical_smiles.iloc[indices[0]]),
            "representative_role": "group medoid" if cluster_col else "Butina center",
            "core_smiles": core, "core_exact_scaffold_count": counts.get(core, 0),
            "core_coverage_count": covered, "core_coverage_fraction": covered / len(indices),
            "core_member_ids": [keys[order[i]] for i, yes in zip(indices, matched) if yes],
            "core_rule": core_rule, "core_heavy_atoms": Chem.MolFromSmiles(core).GetNumHeavyAtoms() if core else 0,
            "core_coverage_target": core_min_coverage if core_rule == "coverage" else None,
            "core_target_met": None if choice is None else bool(choice["target_met"]),
            "most_frequent_scaffold": most_frequent,
            "most_frequent_scaffold_coverage_count": sum(frequent_matched),
            "distinct_scaffolds": len(counts), "acyclic_members": scaffolds.count(""),
            "core_status": "no ring scaffold" if not core else "selected scaffold has majority coverage" if covered > len(indices) / 2 else "selected scaffold has no majority coverage",
            "scaffold_counts": [{"smiles": s, "count": counts[s]} for s in ranked],
            "minimum_center_similarity": float(1 - distances[indices[0], indices].max()),
        })
    protocol = {"method": "external" if cluster_col else "RDKit Butina", "source": cluster_source if cluster_col else "normalized molecular structures",
                "cluster_column": cluster_col, "similarity_threshold": None if cluster_col else similarity,
                "distance": "1 - binary Morgan Tanimoto", "radius": radius, "n_bits": n_bits,
                "use_chirality_in_fingerprint": False, "reordering": True,
                "input_order": "canonical isomeric SMILES, then molecular identifier",
                "rdkit_version": rdBase.rdkitVersion, "scaffold": "typed Bemis-Murcko, chirality retained",
                "core_rule": ("largest member scaffold or scaffold-network ring fragment matching at least "
                              f"{core_min_coverage:.0%} of members; otherwise the highest-coverage candidate"
                              if core_rule == "coverage" else
                              "most frequent nonempty exact scaffold; ties: heavy atoms descending, then SMILES"),
                "core_min_coverage": core_min_coverage if core_rule == "coverage" else None,
                "coverage_rule": "explicit substructure matches with chirality; not inferred from cluster membership",
                "interpretation": "Groups are representation- and threshold-dependent, not chemical ground truth. Scaffold summaries do not explain every fingerprint bit. Colors may repeat; use group IDs.",
                "structure_sha256": hashlib.sha256(ordered.canonical_smiles.to_csv(index=False).encode()).hexdigest()}
    return frame, summaries, protocol


def group_layout_metrics(data, coordinates, *, k=10, space="grid"):
    """Report each group's same-group neighbor fraction and random-label baseline."""
    n = len(data)
    if n < 2:
        return pd.DataFrame()
    labels = data.structure_group.astype(str).to_numpy()
    distance = cdist(np.asarray(coordinates, float), np.asarray(coordinates, float))
    np.fill_diagonal(distance, np.inf)
    effective_k = min(k, n - 1)
    radius = np.partition(distance, effective_k - 1, axis=1)[:, effective_k - 1]
    neighbor = distance <= radius[:, None] + 1e-12
    fraction = (neighbor & (labels[:, None] == labels[None, :])).sum(axis=1) / neighbor.sum(axis=1)
    rows = []
    for group in sorted(set(labels)):
        members = labels == group
        size = int(members.sum())
        observed = float(fraction[members].mean())
        expected = (size - 1) / (n - 1)
        rows.append({"structure_group": group, "size": size, "space": space, "k": effective_k,
                     "same_group_fraction": observed, "random_label_expectation": expected,
                     "enrichment": observed / expected if expected else np.nan,
                     "excess_fraction": observed - expected,
                     "members_with_no_same_group_neighbor": int(np.sum(fraction[members] == 0))})
    return pd.DataFrame(rows)


def _svg_inner(svg):
    return svg[svg.index(">", svg.index("<svg")) + 1:svg.rfind("</svg>")]


def aligned_comparison(data, core_smiles, *, columns=3, cell_width=300, cell_height=255):
    """Align matching cores; highlight their complement, not an inferred causal change."""
    if len(data) == 0:
        raise ValueError("Select at least one molecule for comparison.")
    keys = molecule_keys(data)
    core = Chem.MolFromSmiles(core_smiles) if core_smiles else None
    if core is not None:
        rdDepictor.Compute2DCoords(core)
    mols, audits = [], []
    for key, row in zip(keys, data.to_dict("records")):
        mol = Chem.MolFromSmiles(str(row["canonical_smiles"]))
        original = Chem.MolToSmiles(mol, isomericSmiles=True)
        matches = mol.GetSubstructMatches(core, useChirality=True, uniquify=False, maxMatches=1024) if core else ()
        match = min(matches) if matches else ()
        if match:
            rdDepictor.GenerateDepictionMatching2DStructure(mol, core, list(enumerate(match)))
        else:
            rdDepictor.Compute2DCoords(mol)
        complement = sorted(set(range(mol.GetNumAtoms())) - set(match)) if match else []
        bonds = [b.GetIdx() for b in mol.GetBonds() if b.GetBeginAtomIdx() in complement or b.GetEndAtomIdx() in complement]
        audit = {"comparison_id": key, "core_smiles": core_smiles, "aligned": bool(match),
                 "core_matched": bool(match), "alignment_status": "core_aligned" if match else "core_absent",
                 "core_atom_indices": list(match), "highlight_atom_indices": complement,
                 "highlight_bond_indices": bonds, "match_count_capped_1024": len(matches),
                 "atom_match_rule": "lexicographically first chirality-aware match; symmetric matches may be equivalent",
                 "canonical_isomeric_smiles_before": original,
                 "canonical_isomeric_smiles_after": Chem.MolToSmiles(mol, isomericSmiles=True)}
        roundtrip = Chem.MolFromMolBlock(Chem.MolToMolBlock(mol))
        audit["molblock_roundtrip_preserves_structure"] = roundtrip is not None and Chem.MolToSmiles(roundtrip, isomericSmiles=True) == original
        if not audit["molblock_roundtrip_preserves_structure"]:
            # A template can contradict a molecule's E/Z geometry. Never export that alignment.
            mol = Chem.MolFromSmiles(original)
            rdDepictor.Compute2DCoords(mol)
            roundtrip = Chem.MolFromMolBlock(Chem.MolToMolBlock(mol))
            audit["aligned"] = False
            audit["alignment_status"] = "stereo_safe_unaligned"
            audit["molblock_roundtrip_preserves_structure"] = roundtrip is not None and Chem.MolToSmiles(roundtrip, isomericSmiles=True) == original
        audit["structure_preserved"] = audit["canonical_isomeric_smiles_after"] == original and audit["molblock_roundtrip_preserves_structure"]
        if not audit["structure_preserved"]:
            raise ValueError(f"Depiction structure check failed for {key}.")
        mols.append(mol)
        audits.append(audit)
    # One world-coordinate box gives every card the same scale and core position.
    positions = np.concatenate([m.GetConformer().GetPositions()[:, :2] for m in mols])
    lo, hi = positions.min(axis=0) - .7, positions.max(axis=0) + .7
    width, height = columns * cell_width, int(np.ceil(len(mols) / columns)) * cell_height
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">', '<rect width="100%" height="100%" fill="white"/>']
    for i, (mol, audit, (_, row)) in enumerate(zip(mols, audits, data.iterrows())):
        x, y = (i % columns) * cell_width, (i // columns) * cell_height
        label = str(row.get("activity_class", "unlabelled")).lower()
        color = DEFAULT_COLORS.get(label, DEFAULT_COLORS["unlabelled"])
        drawer = rdMolDraw2D.MolDraw2DSVG(cell_width - 10, cell_height - 62)
        options = drawer.drawOptions()
        options.clearBackground = False
        options.bondLineWidth = 1.8
        options.minFontSize = 10
        options.maxFontSize = 18
        options.highlightBondWidthMultiplier = 7
        drawer.SetScale(cell_width - 10, cell_height - 62, Geometry.Point2D(*lo), Geometry.Point2D(*hi))
        hc = (.98, .73, .40)
        drawer.DrawMolecule(mol, highlightAtoms=audit["highlight_atom_indices"],
                            highlightBonds=audit["highlight_bond_indices"],
                            highlightAtomColors={a: hc for a in audit["highlight_atom_indices"]},
                            highlightBondColors={b: hc for b in audit["highlight_bond_indices"]})
        audit["core_draw_coordinates"] = [[drawer.GetDrawCoords(a).x, drawer.GetDrawCoords(a).y] for a in audit["core_atom_indices"]]
        drawer.FinishDrawing()
        value = row.get("activity_pchembl", np.nan)
        value_text = f"pChEMBL {float(value):.2f}" if pd.notna(value) else "activity not supplied"
        text = f'{value_text} | {label}'
        if not audit["aligned"]:
            text += " | " + ("unaligned: stereo safety" if audit['core_matched'] else "core absent")
        parts.extend([f'<g data-comparison-id="{html.escape(audit["comparison_id"], quote=True)}" transform="translate({x},{y})">',
                      f'<rect x="3" y="3" width="{cell_width-6}" height="{cell_height-6}" rx="5" fill="white" stroke="#CAD4DB"/>',
                      f'<rect x="4" y="4" width="{cell_width-8}" height="4" fill="{color}"/>',
                      f'<g transform="translate(5,11)">{_svg_inner(drawer.GetDrawingText())}</g>',
                      f'<text x="12" y="{cell_height-31}" font-family="sans-serif" font-size="13" fill="#182B3B">{html.escape(audit["comparison_id"])}</text>',
                      f'<text x="12" y="{cell_height-13}" font-family="sans-serif" font-size="11" fill="#465965">{html.escape(text)}</text>', '</g>'])
    parts.append("</svg>")
    return "\n".join(parts), audits


def group_comparison(data, summary, limit=12):
    group = data[data.structure_group.eq(summary["structure_group"])].copy()
    group = group.sort_values(["scaffold_smiles", "canonical_smiles", "comparison_id"])
    # Include the representative first; no activity-based selection.
    representative = group[group.comparison_id.eq(summary["representative_id"])]
    group = pd.concat([representative, group[~group.comparison_id.eq(summary["representative_id"])]])
    selected = group.head(limit)
    svg, audit = aligned_comparison(selected, summary["core_smiles"])
    return selected, svg, audit


def save_group_outputs(data, summaries, protocol, output_dir, name="map"):
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    files = {}
    for suffix, contents in [("structure_groups", {"protocol": protocol, "groups": summaries})]:
        path = directory / f"{name}_{suffix}.json"
        path.write_text(json.dumps(contents, indent=2, allow_nan=False), encoding="utf-8")
        files[suffix] = str(path)
    path = directory / f"{name}_structure_members.csv"
    data.to_csv(path, index=False)
    files["structure_members"] = str(path)
    if {"grid_row", "grid_col"} <= set(data):
        palette = {s["structure_group"].lower(): s["color"] for s in summaries}
        path = directory / f"{name}_structure_groups.svg"
        render_svg(data, path, label_col="structure_group", color_map=palette,
                   draw_molecules=len(data) <= 2000)
        files["structure_svg"] = str(path)
    return files


def main(argv=None):
    parser = argparse.ArgumentParser(description="Compare structural groups in an existing coordinate table without changing the layout.")
    parser.add_argument("input", type=Path)
    parser.add_argument("-o", "--output-dir", type=Path, required=True)
    parser.add_argument("--similarity", type=float, default=.6)
    parser.add_argument("--cluster-col")
    parser.add_argument("--cluster-source")
    parser.add_argument("--group", help="Group ID to export; otherwise the largest group.")
    args = parser.parse_args(argv)
    data, groups, protocol = structural_groups(pd.read_csv(args.input), similarity=args.similarity,
                                                cluster_col=args.cluster_col, cluster_source=args.cluster_source)
    protocol["input_file_sha256"] = hashlib.sha256(args.input.read_bytes()).hexdigest()
    files = save_group_outputs(data, groups, protocol, args.output_dir)
    chosen = next((g for g in groups if g["structure_group"] == args.group), None) if args.group else groups[0]
    if chosen is None:
        parser.error("Unknown structure group.")
    selected, svg, audit = group_comparison(data, chosen)
    (args.output_dir / "comparison.svg").write_text(svg, encoding="utf-8")
    (args.output_dir / "atom_matches.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    selected.to_csv(args.output_dir / "comparison_members.csv", index=False)
    print(json.dumps({"groups": len(groups), "selected_group": chosen["structure_group"], "outputs": files}, indent=2))
    return 0
