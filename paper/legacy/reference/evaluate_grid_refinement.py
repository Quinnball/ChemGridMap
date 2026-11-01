"""Evaluate fingerprint-neighbour grid refinement on all frozen manuscript layouts.

Refinement parameters were fixed on CA II seed 42; the other ten layouts are
held out. Besides the optimized quantity (k=10 Morgan-neighbour recall), the
script reports quantities the refinement does not optimize: other k, atom-pair
fingerprint neighbours, same-scaffold neighbours, trustworthiness, global
distance correlation, displacement and cross-seed stability.
"""
import itertools
import json
import time

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr
from sklearn.manifold import trustworthiness

from common import (RESULTS, TARGETS, _tie_inclusive_neighbor_indices, fingerprint_neighbors,
                    fingerprints, grid_xy, load_layout, mean_jaccard, per_molecule_recall,
                    space_neighbors)
from chemgridmap.refinement import REFINEMENT_DEFAULTS, refine_grid_layout

OUT = RESULTS / "grid_refinement"
OUT.mkdir(parents=True, exist_ok=True)


def lattice_side(layout):
    return int(np.rint(1 / np.diff(np.unique(layout.grid_x))[0])) + 1


def atom_pair_neighbors(target, k=10):
    gen = rdFingerprintGenerator.GetAtomPairGenerator(fpSize=2048)
    mols = [Chem.MolFromSmiles(s) for s in load_layout(target, 42).canonical_smiles]
    matrix = np.array([gen.GetFingerprintAsNumPy(m) for m in mols]).astype(bool)
    return _tie_inclusive_neighbor_indices(matrix, k, metric="jaccard")


def scaffolds(target):
    out = []
    for smiles in load_layout(target, 42).canonical_smiles:
        try:
            out.append(MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(smiles)))
        except Exception:  # pragma: no cover - defensive for unusual inputs
            out.append("")
    return np.array(out, dtype=object)


def same_scaffold_fraction(neighbors, scaffold):
    values = [np.mean(scaffold[row] == scaffold[i]) for i, row in enumerate(neighbors) if scaffold[i]]
    return float(np.mean(values))


def mean_abs_activity_difference(neighbors, values):
    return float(np.mean([np.mean(np.abs(values[row] - values[i])) for i, row in enumerate(neighbors)]))


rows, per_molecule, stability, layouts_out, diagnostics_all = [], [], [], [], {}
for target, (label, seeds) in TARGETS.items():
    fp_bool = fingerprints(target)
    fp = {k: fingerprint_neighbors(target, k) for k in (5, 10, 20)}
    ap = atom_pair_neighbors(target)
    scaffold = scaffolds(target)
    activity = load_layout(target, 42).activity_pchembl.to_numpy(float)
    fp_activity = mean_abs_activity_difference(fp[10], activity)
    refined_by_seed, base_by_seed = {}, {}
    for seed in seeds:
        layout = load_layout(target, seed)
        side = lattice_side(layout)
        cells = layout[["grid_col", "grid_row"]].to_numpy(int)
        projection = layout[["projection_x", "projection_y"]].to_numpy(float) * (side - 1)
        diag = {}
        t0 = time.perf_counter()
        refined = refine_grid_layout(cells, projection, fp[10], side, diagnostics=diag)
        seconds = time.perf_counter() - t0
        diagnostics_all[f"{target}_seed{seed}"] = diag
        layouts = {"minimum_displacement": cells / (side - 1), "neighbor_refined": refined / (side - 1)}
        base_by_seed[seed], refined_by_seed[seed] = layouts["minimum_displacement"], layouts["neighbor_refined"]
        proj01 = layout[["projection_x", "projection_y"]].to_numpy(float)
        proj_d = pdist(proj01)
        for name, xy in layouts.items():
            grid_nb = {k: space_neighbors(xy, k) for k in (5, 10, 20)}
            rec10 = per_molecule_recall(fp[10], grid_nb[10])
            rows.append({
                "target": label, "seed": seed, "layout": name, "n": len(xy), "side": side,
                "recall_k10": rec10.mean(),
                "precision_k10": float(np.mean([len(set(f) & set(g)) / len(g) for f, g in zip(fp[10], grid_nb[10])])),
                "recall_k5": per_molecule_recall(fp[5], grid_nb[5]).mean(),
                "recall_k20": per_molecule_recall(fp[20], grid_nb[20]).mean(),
                "atom_pair_recall_k10": per_molecule_recall(ap, grid_nb[10]).mean(),
                "same_scaffold_fraction": same_scaffold_fraction(grid_nb[10], scaffold),
                "trustworthiness_k10": float(trustworthiness(fp_bool, xy, n_neighbors=10, metric="jaccard")),
                "spearman_projection_vs_grid": float(spearmanr(proj_d, pdist(xy)).statistic),
                "mean_displacement_cells": float(np.linalg.norm(proj01 - xy, axis=1).mean() * (side - 1)),
                "grid_neighbor_abs_dpchembl": mean_abs_activity_difference(grid_nb[10], activity),
                "fingerprint_neighbor_abs_dpchembl": fp_activity,
                "refinement_seconds": seconds if name == "neighbor_refined" else 0.0,
                "moved_molecules": diag["moved_molecules"] if name == "neighbor_refined" else 0,
            })
            if seed == 42:
                for i, value in enumerate(rec10):
                    per_molecule.append({"target": label, "layout": name, "molecule_index": i,
                                         "molecule_identity_key": layout.molecule_identity_key[i],
                                         "fingerprint_neighbor_recall": value})
                frame = layout[["molecule_identity_key", "molecule_id", "activity_pchembl", "activity_class"]].copy()
                frame["target"], frame["layout"] = label, name
                frame["grid_col"] = np.rint(xy[:, 0] * (side - 1)).astype(int)
                frame["grid_row"] = np.rint(xy[:, 1] * (side - 1)).astype(int)
                frame["fingerprint_neighbor_recall"] = rec10
                layouts_out.append(frame)
        print(f"{label} seed {seed}: recall {rows[-2]['recall_k10']:.3f} -> {rows[-1]['recall_k10']:.3f} "
              f"({seconds:.1f}s)", flush=True)
    for a, b in itertools.combinations(seeds, 2):
        stability.append({"target": label, "seed_a": a, "seed_b": b,
                          "minimum_displacement": mean_jaccard(space_neighbors(base_by_seed[a]), space_neighbors(base_by_seed[b])),
                          "neighbor_refined": mean_jaccard(space_neighbors(refined_by_seed[a]), space_neighbors(refined_by_seed[b]))})

metrics = pd.DataFrame(rows)
metrics.to_csv(OUT / "layout_metrics.csv", index=False)
pd.DataFrame(stability).to_csv(OUT / "seed_stability.csv", index=False)
pd.DataFrame(per_molecule).to_csv(OUT / "per_molecule_recall_seed42.csv", index=False)
pd.concat(layouts_out).to_csv(OUT / "layouts_seed42.csv", index=False)
(OUT / "refinement_diagnostics.json").write_text(json.dumps(
    {"parameters": {k: list(v) if isinstance(v, tuple) else v for k, v in REFINEMENT_DEFAULTS.items()},
     "tuning_layout": "CA II seed 42", "runs": diagnostics_all}, indent=1, default=float))

summary = metrics.groupby(["target", "layout"]).mean(numeric_only=True).drop(columns=["seed", "n", "side"])
print(summary.round(3).T.to_string())
print(pd.DataFrame(stability).groupby("target")[["minimum_displacement", "neighbor_refined"]].mean().round(3))
