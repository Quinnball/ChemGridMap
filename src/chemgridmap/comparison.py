"""Structure-defined groups and core-aligned comparisons, independent of layout."""

from __future__ import annotations

import argparse
import atexit
from collections import Counter
from concurrent.futures import as_completed
import hashlib
import html
import json
import os
from pathlib import Path
import re
import threading
import time

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, Geometry, rdBase
from rdkit.Chem import rdDepictor, rdFingerprintGenerator
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Cluster import Butina
from scipy.spatial.distance import cdist

from .cores import mcs_core, select_core
from .plotting import DEFAULT_COLORS, render_svg

GROUP_COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#8B6914", "#56B4E9",
                "#5B536E", "#B26249", "#658235", "#3B6C70", "#A95479", "#777777"]

# Small maps finish faster without starting workers. These are execution choices,
# not chemical parameters, and therefore stay out of the reproducible manifest.
# CA II subsets: 32 groups were slower with cold workers; 64 groups were faster
# both cold and warm. Keep smaller workloads serial even when a pool is available.
CORE_PARALLEL_MIN_GROUPS = 64
_core_executor = None
_core_executor_jobs = None
_core_executor_lock = threading.Lock()


def browser_core_jobs():
    """Browser-only worker default; library/CLI callers continue to default to one."""
    default = max(1, min(8, (os.cpu_count() or 1) - 1))
    try:
        return max(1, int(os.environ.get("CHEMGRIDMAP_JOBS", default)))
    except ValueError:
        return default


def shutdown_core_workers():
    """Finish submitted work and close this application's reusable loky pool."""
    global _core_executor, _core_executor_jobs
    with _core_executor_lock:
        executor, _core_executor, _core_executor_jobs = _core_executor, None, None
    if executor is not None:
        # Graceful shutdown also works in macOS sandboxes that disallow pgrep,
        # which loky's kill_workers process-tree traversal otherwise requires.
        executor.shutdown(wait=True)


atexit.register(shutdown_core_workers)


def _get_core_executor(jobs):
    global _core_executor, _core_executor_jobs
    from joblib.externals.loky import ProcessPoolExecutor
    from joblib.externals.loky.backend.context import get_context

    with _core_executor_lock:
        if _core_executor is not None and _core_executor_jobs != jobs:
            _core_executor.shutdown(wait=True)
            _core_executor = None
        if _core_executor is None:
            # loky starts a fresh interpreter, never fork a multithreaded server.
            # Limit nested numerical libraries inside the chemistry workers.
            _core_executor = ProcessPoolExecutor(
                max_workers=jobs, context=get_context("loky"),
                env={name: "1" for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
                                            "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")})
            _core_executor_jobs = jobs
        return _core_executor


def _summarize_core_payload(smiles, member_ids, core_rule, core_min_coverage, fmcs_timeout):
    """Workers receive only strings, IDs and numeric rule parameters."""
    selected = [Chem.MolFromSmiles(value) for value in smiles]
    return summarize_group_core(
        selected, member_ids, core_rule=core_rule, core_min_coverage=core_min_coverage,
        fmcs_timeout=fmcs_timeout, return_metrics=True)


def group_core_members(data, summaries):
    """Collect fixed groups once, preserving the interactive center-first order."""
    ordered = data.sort_values(["canonical_smiles", "comparison_id"])
    grouped = dict(tuple(ordered.groupby("structure_group", sort=False)))
    members = []
    for summary in summaries:
        frame = grouped[summary["structure_group"]]
        center = frame.comparison_id.eq(summary["representative_id"])
        members.append(pd.concat([frame[center], frame[~center]]))
    return members


def summarize_group_cores(groups, *, core_rule="coverage", core_min_coverage=0.8,
                          fmcs_timeout=1, jobs=1, progress=None):
    """Compute ordered ``(summary, metrics)`` pairs from ``(SMILES, IDs)`` groups.

    Each group's representative must be first, followed by the unchanged member
    order used by :func:`summarize_group_core`. Only plain data enters workers.
    Large groups are submitted first, then results return in original group order.
    FMCS is always serial because process load can change its timeout solution.
    ``progress(done, total)`` runs in the caller thread, never in a subprocess.
    """
    jobs = int(jobs)
    if jobs < 1:
        raise ValueError("jobs must be a positive integer.")
    payloads = [(list(smiles), list(ids)) for smiles, ids in groups]
    total = len(payloads)
    if progress is not None:
        progress(0, total)
    parallel = jobs > 1 and core_rule != "fmcs" and total >= CORE_PARALLEL_MIN_GROUPS
    results = [None] * total
    if parallel:
        executor = _get_core_executor(min(jobs, total))
        futures = {
            executor.submit(_summarize_core_payload, *payloads[index], core_rule,
                            core_min_coverage, fmcs_timeout): index
            for index in sorted(range(total), key=lambda index: -len(payloads[index][0]))}
        try:
            for done, future in enumerate(as_completed(futures), 1):
                results[futures[future]] = future.result()
                if progress is not None:
                    progress(done, total)
        except BaseException:
            for future in futures:
                future.cancel()
            raise
    else:
        for index, (smiles, ids) in enumerate(payloads):
            results[index] = _summarize_core_payload(
                smiles, ids, core_rule, core_min_coverage, fmcs_timeout)
            if progress is not None:
                progress(index + 1, total)
    return results


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


def summarize_group_core(selected, member_ids, *, core_rule="coverage", core_min_coverage=0.8,
                         fmcs_timeout=1, return_metrics=False):
    """Calculate one fixed group's core; build and interactive switching share this path.

    Members must be ordered with the representative first and the remaining
    members by canonical SMILES and identifier, as in the grouping protocol.
    Optional telemetry is kept outside legacy summaries to preserve run files.
    """
    if core_rule not in {"most_frequent", "coverage", "coverage_brics", "coverage_recap", "fmcs"}:
        raise ValueError("Unknown core rule.")
    started = time.perf_counter()
    scaffolds = [MurckoScaffold.MurckoScaffoldSmiles(mol=m, includeChirality=True) for m in selected]
    counts = Counter(s for s in scaffolds if s)
    ranked = sorted(counts, key=lambda s: (-counts[s], -Chem.MolFromSmiles(s).GetNumHeavyAtoms(), s))
    most_frequent = ranked[0] if ranked else ""
    frequent_pattern = Chem.MolFromSmiles(most_frequent) if most_frequent else None
    frequent_matched = [bool(frequent_pattern and m.HasSubstructMatch(frequent_pattern, useChirality=True))
                        for m in selected]
    if core_rule in {"coverage_brics", "coverage_recap"} or (core_rule == "coverage" and most_frequent):
        choice = select_core(selected, min_coverage=core_min_coverage,
                             include_brics=core_rule == "coverage_brics",
                             include_recap=core_rule == "coverage_recap")
        core, matched = choice["core_smiles"], choice["matched"]
    elif core_rule == "fmcs":
        choice = mcs_core(selected, min_coverage=core_min_coverage, timeout=int(fmcs_timeout))
        core, matched = choice["core_smiles"], [
            bool(pattern is not None and mol.HasSubstructMatch(pattern, useChirality=True))
            for mol in selected
        ] if (pattern := Chem.MolFromSmarts(choice["core_smarts"])) is not None else [False] * len(selected)
    else:
        choice, core, matched = None, most_frequent, frequent_matched
    covered = sum(matched)
    summary = {
        "core_smiles": core, "core_smiles_attached": choice.get("core_smiles_attached", "") if choice else "",
        "core_smiles_attached_source": (
            choice.get("core_smiles_attached_source", "") if choice else ""
        ),
        "core_exact_scaffold_count": counts.get(core, 0),
        "core_coverage_count": covered, "core_coverage_fraction": covered / len(selected),
        "core_member_ids": [key for key, yes in zip(member_ids, matched) if yes],
        "core_rule": core_rule,
        "core_heavy_atoms": (
            int(choice["core_heavy_atoms"]) if core_rule == "fmcs"
            else Chem.MolFromSmiles(core).GetNumHeavyAtoms() if core else 0
        ),
        "core_smarts": choice.get("core_smarts") if choice else None,
        "core_coverage_target": core_min_coverage if core_rule in {"coverage", "coverage_brics", "coverage_recap", "fmcs"} else None,
        "core_target_met": (
            bool(choice["target_met"]) if core_rule in {"coverage", "coverage_brics", "coverage_recap"} and choice
            else bool(choice["coverage_fraction"] >= core_min_coverage) if core_rule == "fmcs" and choice
            else None
        ),
        "fmcs_timed_out": bool(choice.get("timed_out", False)) if core_rule == "fmcs" and choice else False,
        "most_frequent_scaffold": most_frequent,
        "most_frequent_scaffold_coverage_count": sum(frequent_matched),
        "distinct_scaffolds": len(counts), "acyclic_members": scaffolds.count(""),
        "core_status": "no ring scaffold" if not core else "selected scaffold has majority coverage" if covered > len(selected) / 2 else "selected scaffold has no majority coverage",
        "scaffold_counts": [{"smiles": s, "count": counts[s]} for s in ranked],
    }
    metrics = {"candidates": choice.get("candidates", 0) if choice else len(counts),
               "seconds": time.perf_counter() - started}
    if core_rule == "coverage_brics":
        summary.update(candidates=metrics["candidates"],
                       brics_candidates=choice.get("brics_candidates", 0) if choice else 0,
                       brics_single_pass_members=choice.get("brics_single_pass_members", 0) if choice else 0)
    if core_rule == "coverage_recap":
        summary.update(candidates=metrics["candidates"],
                       recap_candidates=choice.get("recap_candidates", 0) if choice else 0,
                       recap_skipped_members=choice.get("recap_skipped_members", 0) if choice else 0)
    return (summary, metrics) if return_metrics else summary


def structural_groups(data, *, similarity=0.6, cluster_col=None, cluster_source=None,
                      radius=2, n_bits=2048, core_rule="most_frequent", core_min_coverage=0.8,
                      fmcs_timeout=1, jobs=1, progress=None):
    """Butina groups use a center threshold, not an all-pairs similarity bound.

    ``core_rule="coverage"`` selects the largest member scaffold or scaffold-network
    ring fragment that matches at least ``core_min_coverage`` of the members.

    Parameters
    ----------
    data : pandas.DataFrame
        Molecule table with canonical SMILES and structure metadata.
    similarity : float, default 0.6
        Morgan Tanimoto similarity threshold to a Butina center.
    cluster_col : str, optional
        Existing group-label column for external grouping.
    cluster_source : str, optional
        Provenance of external group labels.
    radius : int, default 2
        Morgan fingerprint radius.
    n_bits : int, default 2048
        Morgan fingerprint length.
    core_rule : {"coverage", "coverage_brics", "coverage_recap", "most_frequent", "fmcs"}, default "most_frequent"
        Group-core selection rule.
    core_min_coverage : float, default 0.8
        Coverage target for coverage and FMCS cores.
    fmcs_timeout : int, default 1
        FMCS search budget in seconds.
    jobs : int, default 1
        Core worker processes. Small maps and FMCS always use the serial path.
    progress : callable, optional
        Called with completed and total group counts during core calculation.

    Returns
    -------
    grouped : pandas.DataFrame
        Copy of ``data`` with group, scaffold and core-member columns.
    summaries : list of dict
        Per-group sizes, core choices, coverage and core evidence.
    protocol : dict
        Grouping and core-selection parameters and provenance.
    """
    if not 0 < similarity <= 1:
        raise ValueError("Group similarity must be in (0, 1].")
    if core_rule not in {"most_frequent", "coverage", "coverage_brics", "coverage_recap", "fmcs"}:
        raise ValueError("core_rule must be 'most_frequent', 'coverage', 'coverage_brics', 'coverage_recap' or 'fmcs'.")
    if int(fmcs_timeout) < 1:
        raise ValueError("fmcs_timeout must be a positive integer.")
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
    core_results = summarize_group_cores(
        [([str(ordered.canonical_smiles.iloc[i]) for i in indices],
          [keys[order[i]] for i in indices]) for _, indices in clusters],
        core_rule=core_rule, core_min_coverage=core_min_coverage,
        fmcs_timeout=fmcs_timeout, jobs=jobs, progress=progress)
    summaries = []
    for rank, (label, indices) in enumerate(clusters):
        selected = [mols[i] for i in indices]
        core_summary = core_results[rank][0]
        scaffolds = [MurckoScaffold.MurckoScaffoldSmiles(mol=m, includeChirality=True) for m in selected]
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
            **core_summary,
            "minimum_center_similarity": float(1 - distances[indices[0], indices].max()),
        })
    protocol = {"method": "external" if cluster_col else "RDKit Butina", "source": cluster_source if cluster_col else "normalized molecular structures",
                "cluster_column": cluster_col, "similarity_threshold": None if cluster_col else similarity,
                "distance": "1 - binary Morgan Tanimoto", "radius": radius, "n_bits": n_bits,
                "use_chirality_in_fingerprint": False, "reordering": True,
                "input_order": "canonical isomeric SMILES, then molecular identifier",
                "rdkit_version": rdBase.rdkitVersion, "scaffold": "typed Bemis-Murcko, chirality retained",
                "core_rule": ("FMCS matching at least "
                              f"{core_min_coverage:.0%} of members; timeout {int(fmcs_timeout)} s"
                              if core_rule == "fmcs" else
                              "largest member scaffold, scaffold-network ring fragment or BRICS fragment "
                              "(ring-containing, dummy atoms removed) matching at least "
                              f"{core_min_coverage:.0%} of members; otherwise the highest-coverage candidate"
                              if core_rule == "coverage_brics" else
                              "largest member scaffold, scaffold-network ring fragment or RECAP fragment "
                              "(ring-containing, dummy atoms removed) matching at least "
                              f"{core_min_coverage:.0%} of members; otherwise the highest-coverage candidate"
                              if core_rule == "coverage_recap" else
                              "largest member scaffold or scaffold-network ring fragment matching at least "
                              f"{core_min_coverage:.0%} of members; otherwise the highest-coverage candidate"
                              if core_rule == "coverage" else
                              "most frequent nonempty exact scaffold; ties: heavy atoms descending, then SMILES"),
                "core_min_coverage": core_min_coverage if core_rule in {"coverage", "coverage_brics", "coverage_recap", "fmcs"} else None,
                "fmcs_timeout": int(fmcs_timeout) if core_rule == "fmcs" else None,
                "coverage_rule": "explicit substructure matches with chirality; not inferred from cluster membership",
                "interpretation": "Groups are representation- and threshold-dependent, not chemical ground truth. Scaffold summaries do not explain every fingerprint bit. Colors may repeat; use group IDs.",
                "structure_sha256": hashlib.sha256(
                    ordered.canonical_smiles.to_csv(index=False, lineterminator="\n").encode("utf-8")
                ).hexdigest()}
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
    core = (Chem.MolFromSmarts(core_smiles) or Chem.MolFromSmiles(core_smiles)) if core_smiles else None
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
    svg, audit = aligned_comparison(selected, summary.get("core_smarts") or summary["core_smiles"])
    return selected, svg, audit


def save_group_outputs(data, summaries, protocol, output_dir, name="map"):
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    files = {}
    for suffix, contents in [("structure_groups", {"protocol": protocol, "groups": summaries})]:
        path = directory / f"{name}_{suffix}.json"
        path.write_text(json.dumps(contents, indent=2, allow_nan=False),
                        encoding="utf-8", newline="\n")
        files[suffix] = str(path)
    path = directory / f"{name}_structure_members.csv"
    identity_columns = [
        "molecule_identity_key", "molecule_id", "parent_molecule_id", "name",
        "comparison_id", "canonical_smiles", "smiles", "representation_smiles",
    ]
    grid_columns = ["grid_x", "grid_y", "grid_row", "grid_col"]
    grouping_columns = [
        "structure_group", "structure_group_color", "scaffold_smiles",
        "query_id", "measurement_status", "query_group_center_similarity", "contains_group_core",
    ]
    member_columns = [
        column
        for column in identity_columns + grid_columns + grouping_columns
        if column in data
    ]
    members = data[member_columns].copy()
    if "molecule_identity_key" not in members and "comparison_id" in members:
        members.insert(0, "molecule_identity_key", members["comparison_id"])
    if "comparison_id" in members:
        summary_frame = pd.DataFrame([
            {
                "structure_group": summary["structure_group"],
                "core_rule": summary.get("core_rule"),
                "core_smiles": summary.get("core_smiles"),
                "core_smarts": summary.get("core_smarts"),
                "core_coverage_target": summary.get("core_coverage_target"),
                "core_member_ids": summary.get("core_member_ids", []),
            }
            for summary in summaries
        ])
        members = members.merge(
            summary_frame, left_on="structure_group", right_on="structure_group", how="left"
        )
        members["core_member"] = [
            comparison_id in set(member_ids if isinstance(member_ids, list) else [])
            for comparison_id, member_ids in zip(
                members["comparison_id"], members["core_member_ids"]
            )
        ]
        members = members.drop(columns="core_member_ids")
    members.to_csv(path, index=False, encoding="utf-8", lineterminator="\n")
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
    data, groups, protocol = structural_groups(pd.read_csv(args.input, encoding="utf-8-sig"),
                                                similarity=args.similarity,
                                                cluster_col=args.cluster_col, cluster_source=args.cluster_source)
    protocol["input_file_sha256"] = hashlib.sha256(args.input.read_bytes()).hexdigest()
    files = save_group_outputs(data, groups, protocol, args.output_dir)
    chosen = next((g for g in groups if g["structure_group"] == args.group), None) if args.group else groups[0]
    if chosen is None:
        parser.error("Unknown structure group.")
    selected, svg, audit = group_comparison(data, chosen)
    (args.output_dir / "comparison.svg").write_text(svg, encoding="utf-8", newline="\n")
    (args.output_dir / "atom_matches.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8", newline="\n")
    selected.to_csv(args.output_dir / "comparison_members.csv", index=False,
                    encoding="utf-8", lineterminator="\n")
    print(json.dumps({"groups": len(groups), "selected_group": chosen["structure_group"], "outputs": files}, indent=2))
    return 0
