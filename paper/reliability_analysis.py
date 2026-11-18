"""Reproducible grid analyses for the reliability revision (T2a)."""
from __future__ import annotations

from datetime import datetime, timezone
from functools import lru_cache
import hashlib
import itertools
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr
from sklearn.manifold import trustworthiness

import chemgridmap
from chemgridmap.core import assign_to_grid, normalize_coordinates
from chemgridmap.projection import project_2d
from chemgridmap.provenance import environment_versions
from chemgridmap.refinement import REFINEMENT_DEFAULTS, refine_grid_layout
from chemgridmap.reliability import representation_neighbor_sets, grid_neighbor_sets, neighbor_recall
from chemgridmap.representations import morgan_fingerprints

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "paper/output/revision_v6"
OUTPUT = ROOT / "paper/output/reliability_v1"
SEEDS = (7, 21, 42, 84, 126, 168, 210, 252, 294, 336)
TARGETS = {"chembl205": "CA II", "chembl204": "Thrombin", "chembl240": "hERG"}
LAYOUTS = ("projection", "minimum_displacement", "neighbor_refined")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def dump_json(path, payload):
    Path(path).write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=str,
                                    allow_nan=False) + "\n", encoding="utf-8")


def write_manifest(directory, inputs, parameters):
    """Hash explicit inputs, current source files, and every output except self."""
    directory = Path(directory)
    sources = [p for base in (ROOT / "src", ROOT / "paper") for p in base.rglob("*")
               if p.is_file() and "output" not in p.parts
               and p.suffix in {".py", ".js", ".cjs", ".html"}]
    diff = subprocess.run(["git", "diff", "release-v0.4.0"], cwd=ROOT,
                          check=True, capture_output=True).stdout
    versions = environment_versions()
    versions.update({"chemgridmap.__version__": chemgridmap.__version__,
                     "chemgridmap.__file__": chemgridmap.__file__,
                     "thread_environment": {k: os.environ.get(k) for k in
                         ("NUMBA_NUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
                          "VECLIB_MAXIMUM_THREADS")}})
    def mapping(paths):
        return {str(p.resolve()): sha256(p) for p in sorted(set(map(Path, paths)))}
    manifest = {"created_at_utc": datetime.now(timezone.utc).isoformat(),
                "inputs_sha256": mapping(inputs), "parameters": parameters,
                "environment": versions,
                "code": {"git_diff_release_v0.4.0_sha256": hashlib.sha256(diff).hexdigest(),
                         "source_sha256": mapping(sources)},
                "outputs_sha256": mapping(p for p in directory.rglob("*") if p.is_file()
                    and p != directory / "manifest.json" and "__pycache__" not in p.parts)}
    dump_json(directory / "manifest.json", manifest)
    return manifest


@lru_cache(maxsize=3)
def load_target(target):
    data = pd.read_csv(FROZEN / target / "map_molecules.csv")
    if data.molecule_identity_key.tolist() != sorted(data.molecule_identity_key):
        raise ValueError(f"Frozen input ordering changed: {target}")
    return data, morgan_fingerprints(data.canonical_smiles.tolist())


def jaccard_distances(fp):
    bits = (np.asarray(fp) != 0).astype(np.float32)
    shared = np.rint(bits @ bits.T).astype(float)
    counts = bits.sum(axis=1).astype(float)
    union = counts[:, None] + counts[None, :] - shared
    return 1 - np.divide(shared, union, out=np.ones_like(shared), where=union > 0)


def mean_jaccard(a, b):
    return float(np.mean([len(set(x) & set(y)) / len(set(x) | set(y))
                          for x, y in zip(a, b)]))


def compare_frozen(current, frozen):
    if current.molecule_identity_key.tolist() != frozen.molecule_identity_key.tolist():
        raise ValueError("Frozen and regenerated identity order differ")
    cols = ["grid_col", "grid_row"]
    delta = np.abs(current[["projection_x", "projection_y"]].to_numpy() -
                   frozen[["projection_x", "projection_y"]].to_numpy())
    raw_delta = np.abs(current[["projection_x_raw", "projection_y_raw"]].to_numpy() -
                       frozen[["projection_x_raw", "projection_y_raw"]].to_numpy())
    mismatches = int(np.any(current[cols].to_numpy() != frozen[cols].to_numpy(), axis=1).sum())
    return {"n": len(current), "grid_mismatched_entries": mismatches,
            "projection_max_absolute_difference": float(delta.max()),
            "projection_raw_max_absolute_difference": float(raw_delta.max()),
            "passed": bool(mismatches == 0 and max(delta.max(), raw_delta.max()) <= 1e-5)}


def run_layouts(output_root=OUTPUT, seeds=SEEDS, targets=TARGETS):
    out = Path(output_root) / "layouts"
    out.mkdir(parents=True, exist_ok=True)
    timings, checks, inputs = [], [], []
    for target in targets:
        data, fp = load_target(target)
        inputs.append(FROZEN / target / "map_molecules.csv")
        neighbors = representation_neighbor_sets(fp, k=10)
        folder = out / target
        folder.mkdir(exist_ok=True)
        for seed in seeds:
            t0 = time.perf_counter()
            xy = project_2d(fp, method="umap", random_state=seed, umap_metric="jaccard")
            projection_seconds = time.perf_counter() - t0
            t0 = time.perf_counter()
            grid, side, _ = assign_to_grid(xy, method="dense")
            assignment_seconds = time.perf_counter() - t0
            frame = data.copy()
            frame[["projection_x_raw", "projection_y_raw"]] = xy
            frame[["projection_x", "projection_y"]] = normalize_coordinates(xy)
            frame[["grid_x", "grid_y"]] = grid
            cells = np.rint(grid * (side - 1)).astype(int)
            frame[["grid_col", "grid_row"]] = cells
            frame["grid_side"] = side
            frozen_path = FROZEN / target / f"coordinates_seed_{seed}.csv"
            if frozen_path.exists():
                inputs.append(frozen_path)
                check = {"target": target, "seed": seed,
                         **compare_frozen(frame, pd.read_csv(frozen_path))}
                checks.append(check)
                dump_json(out / "reproduction_check.json", checks)
                if not check["passed"]:
                    raise RuntimeError(f"Frozen layout reproduction failed: {check}")
            diagnostics = {}
            t0 = time.perf_counter()
            refined = refine_grid_layout(cells, normalize_coordinates(xy) * (side - 1),
                                         neighbors, side, random_state=seed, diagnostics=diagnostics)
            refinement_seconds = time.perf_counter() - t0
            frame[["refined_grid_col", "refined_grid_row"]] = refined
            frame.to_csv(folder / f"coordinates_seed_{seed}.csv", index=False)
            dump_json(folder / f"refinement_seed_{seed}.json", diagnostics)
            timings.append({"target": target, "seed": seed,
                            "projection_seconds": projection_seconds,
                            "assignment_seconds": assignment_seconds,
                            "refinement_seconds": refinement_seconds})
            pd.DataFrame(timings).to_csv(out / "timings.csv", index=False)
            print(f"layouts {target} seed={seed}: {projection_seconds:.1f}s + "
                  f"{assignment_seconds:.1f}s + {refinement_seconds:.1f}s", flush=True)
    headline = OUTPUT / "t1_acceptance_opus/headline/grid_refinement.csv"
    if headline.exists():
        inputs.append(headline)
    write_manifest(out, inputs, {"seeds": list(seeds), "targets": list(targets),
                   "occupancy": .4, "projection": "UMAP/Jaccard", "assignment": "dense",
                   "refinement": REFINEMENT_DEFAULTS, "refinement_random_state": "map seed",
                   "t1_bridge": str(headline), "t1_bridge_recomputed": False})


def layout_points(frame, name):
    if name == "projection":
        return frame[["projection_x", "projection_y"]].to_numpy(float)
    prefix = "refined_" if name == "neighbor_refined" else ""
    return frame[[prefix + "grid_col", prefix + "grid_row"]].to_numpy(float)


def evaluate_neighborhoods(fp, points, references=None, high_distances=None):
    references = references or {k: representation_neighbor_sets(fp, k=k) for k in (5, 10, 20)}
    neighbors = {k: grid_neighbor_sets(points, k=k) for k in (5, 10, 20)}
    recalls = {k: neighbor_recall(references[k], neighbors[k]) for k in (5, 10, 20)}
    high_distances = jaccard_distances(fp) if high_distances is None else high_distances
    result = {f"recall_k{k}": float(recalls[k].mean()) for k in recalls}
    result["precision_k10"] = float(neighbor_recall(neighbors[10], references[10]).mean())
    result["trustworthiness_k10"] = float(trustworthiness(high_distances, points,
                                                         n_neighbors=10, metric="precomputed"))
    return result, recalls[10], neighbors[10]


def run_neighborhoods(output_root=OUTPUT, seeds=SEEDS, targets=TARGETS):
    out = Path(output_root) / "neighborhoods"
    out.mkdir(parents=True, exist_ok=True)
    times = pd.read_csv(Path(output_root) / "layouts/timings.csv").set_index(["target", "seed"])
    rows, individual, inputs = [], [], [Path(output_root) / "layouts/timings.csv"]
    generator = rdFingerprintGenerator.GetAtomPairGenerator(fpSize=2048)
    for target in targets:
        data, fp = load_target(target)
        inputs.append(FROZEN / target / "map_molecules.csv")
        reference = {k: representation_neighbor_sets(fp, k=k) for k in (5, 10, 20)}
        high = jaccard_distances(fp)
        mols = [Chem.MolFromSmiles(s) for s in data.canonical_smiles]
        ap = np.array([generator.GetFingerprintAsNumPy(m) for m in mols])
        ap_neighbors = representation_neighbor_sets(ap, k=10)
        scaffolds = np.array([MurckoScaffold.MurckoScaffoldSmiles(mol=m) for m in mols])
        activity = data.activity_pchembl.to_numpy(float)
        def delta(neighbors):
            return float(np.mean([np.abs(activity[g] - activity[i]).mean()
                                  for i, g in enumerate(neighbors)]))
        for seed in seeds:
            path = Path(output_root) / "layouts" / target / f"coordinates_seed_{seed}.csv"
            inputs.append(path)
            frame = pd.read_csv(path)
            side = int(frame.grid_side.iloc[0])
            projection = layout_points(frame, "projection")
            for name in LAYOUTS:
                points = layout_points(frame, name)
                metrics, recall, neighbors = evaluate_neighborhoods(fp, points, reference, high)
                points_cells = points * (side - 1) if name == "projection" else points
                row = {"target": target, "target_label": TARGETS[target], "seed": seed,
                       "layout": name, "n": len(data), "grid_side": side, **metrics,
                       "atom_pair_recall_k10": float(neighbor_recall(ap_neighbors, neighbors).mean()),
                       "same_scaffold_fraction": float(np.mean([
                           (scaffolds[g] == scaffolds[i]).mean() for i, g in enumerate(neighbors)
                           if scaffolds[i]])),
                       "spearman_projection_vs_layout": float(spearmanr(pdist(projection), pdist(points)).statistic),
                       "mean_displacement_cells": float(np.linalg.norm(points_cells - projection * (side - 1), axis=1).mean()),
                       "neighbor_abs_dpchembl": delta(neighbors),
                       "fingerprint_neighbor_abs_dpchembl": delta(reference[10])}
                timing = times.loc[(target, seed)]
                row["projection_seconds"] = float(timing.projection_seconds)
                row["assignment_seconds"] = float(timing.assignment_seconds) if name != "projection" else 0.
                row["refinement_seconds"] = float(timing.refinement_seconds) if name == "neighbor_refined" else 0.
                row["total_layout_seconds"] = sum(row[x] for x in ("projection_seconds", "assignment_seconds", "refinement_seconds"))
                rows.append(row)
                if seed == 42:
                    individual.extend({"target": target, "seed": seed, "layout": name,
                                       "molecule_identity_key": key, "recall_k10": float(r)}
                                      for key, r in zip(data.molecule_identity_key, recall))
        print(f"neighborhoods {target}: complete", flush=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(out / "layout_metrics.csv", index=False)
    pd.DataFrame(individual).to_csv(out / "per_molecule_recall_seed42.csv", index=False)
    frame.groupby(["target", "layout"])[list(metrics)].agg(["mean", "median", "min", "max"]).to_csv(out / "summary.csv")
    write_manifest(out, inputs, {"k": [5, 10, 20], "ties": "relative 1e-9, absolute 1e-12",
        "same_scaffold_denominator": "focal entries with nonempty Murcko scaffold",
        "runtime_scope": "projection + assignment + refinement; excludes evaluation",
        "seeds": list(seeds), "targets": list(targets)})


def run_stability(output_root=OUTPUT, seeds=SEEDS, targets=TARGETS):
    out = Path(output_root) / "stability"
    out.mkdir(parents=True, exist_ok=True)
    rows, inputs = [], []
    for target in targets:
        neighbors = {}
        for seed in seeds:
            path = Path(output_root) / "layouts" / target / f"coordinates_seed_{seed}.csv"
            inputs.append(path)
            frame = pd.read_csv(path)
            for name in LAYOUTS:
                neighbors[seed, name] = grid_neighbor_sets(layout_points(frame, name), k=10)
        for a, b in itertools.combinations(seeds, 2):
            for name in LAYOUTS:
                rows.append({"target": target, "target_label": TARGETS[target], "layout": name,
                             "seed_a": a, "seed_b": b,
                             "neighbor_jaccard_k10": mean_jaccard(neighbors[a, name], neighbors[b, name])})
    frame = pd.DataFrame(rows)
    frame.to_csv(out / "seed_pair_jaccard.csv", index=False)
    frame.groupby(["target", "layout"]).neighbor_jaccard_k10.agg(["count", "mean", "median", "min", "max"]).to_csv(out / "summary.csv")
    write_manifest(out, inputs, {"seeds": list(seeds), "k": 10, "ties": "included"})
