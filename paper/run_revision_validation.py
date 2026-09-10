"""Reproduce the revision experiments from archived ChEMBL IC50 exports.

Run with the project environment: make paper
No outcome-dependent target, seed, or molecule selection is performed.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import itertools
import json
import os
from pathlib import Path
import time

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".cache/matplotlib"))
os.environ.setdefault("NUMBA_NUM_THREADS", "1")

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from sklearn.manifold import trustworthiness

from chemgridmap.chembl import curate_chembl_activity_data, save_chembl_curation
from chemgridmap.core import (assign_to_grid, compute_metrics, make_candidate_grid,
                             normalize_coordinates, _tie_inclusive_neighbor_indices)
from chemgridmap.plotting import render_svg
from chemgridmap.projection import project_2d
from chemgridmap.provenance import environment_versions
from chemgridmap.representations import morgan_fingerprints, rdkit_descriptors

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = {
    "targets": ["CHEMBL205", "CHEMBL204", "CHEMBL240"],
    "main_seeds": [7, 21, 42, 84, 126], "external_seeds": [7, 21, 42],
    "external_map_cap": 1200, "sampling_seed": 20260908,
    "occupancies": [0.90, 0.70, 0.50, 0.40, 0.30, 0.20],
    "default_occupancy": 0.40, "k": 10,
    "baseline_repeats": 5, "sparse_candidates": [64, 128, 192],
    "morgan_radius": 2, "morgan_bits": 2048,
    "umap_metric": "jaccard", "umap_n_neighbors": 15, "umap_min_dist": 0.1,
    "interpretation": "Descriptive software validation; seeds are not biological replicates.",
}


def dump(path, payload):
    Path(path).write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")


def evaluate(xy, grid, data):
    return compute_metrics(xy, xy, grid, values=data.activity_pchembl.to_numpy(),
                           labels=data.activity_class.to_numpy(), k=PROTOCOL["k"]
                           ).iloc[0].to_dict()


def coordinates(data, xy, grid, side):
    out = data.copy()
    out[["projection_x_raw", "projection_y_raw"]] = xy
    out[["projection_x", "projection_y"]] = normalize_coordinates(xy)
    out[["grid_x", "grid_y"]] = grid
    out["grid_col"] = np.rint(grid[:, 0] * (side - 1)).astype(int)
    out["grid_row"] = np.rint(grid[:, 1] * (side - 1)).astype(int)
    return out


def neighborhood_jaccard(a, b):
    aa = _tie_inclusive_neighbor_indices(a, PROTOCOL["k"])
    bb = _tie_inclusive_neighbor_indices(b, PROTOCOL["k"])
    return float(np.mean([len(set(x) & set(y)) / len(set(x) | set(y))
                          for x, y in zip(aa, bb)]))


def greedy(xy, order):
    cells, side, _ = make_candidate_grid(len(xy), target_occupancy=0.40)
    cost = cdist(normalize_coordinates(xy), cells, metric="sqeuclidean")
    available = np.ones(len(cells), dtype=bool)
    selected = np.empty(len(xy), dtype=int)
    for i in order:
        candidate = np.flatnonzero(available)
        selected[i] = candidate[np.argmin(cost[i, candidate])]
        available[selected[i]] = False
    return cells[selected], side


def assay_family(description):
    text = str(description).lower()
    if "nitrophen" in text or "4-npa" in text or "4npa" in text:
        return "Nitrophenyl substrate mentioned"
    if "ph-stat" in text or "ph stat" in text:
        return "pH-stat mentioned"
    if "stopped" in text or "co2" in text or "carbon dioxide" in text:
        return "CO2 / stopped-flow mentioned"
    return "Other / insufficient description"


def audit_case(curated, out):
    mol = curated.molecules
    selected = mol[mol.parent_molecule_id.fillna("").eq("CHEMBL20")].iloc[0]
    rows = curated.retained_records[
        curated.retained_records.molecule_identity_key.eq(selected.molecule_identity_key)
    ].copy()
    rows["description_group"] = rows.assay_description.map(assay_family)
    rows.to_csv(out / "acetazolamide_retained_records.csv", index=False)
    pd.DataFrame([selected]).to_csv(out / "acetazolamide_molecule.csv", index=False)
    summary = rows.groupby("description_group").pchembl_value.agg(["size", "median", "min", "max"])
    summary.to_csv(out / "acetazolamide_description_groups.csv")
    # This is a lexical audit, not a claim that descriptions exhaust assay conditions.
    dump(out / "acetazolamide_case.json", {
        "molecule_identity_key": selected.molecule_identity_key,
        "all_rows_resolve_to_molecule": bool(rows.molecule_identity_key.nunique() == 1),
        "row_count_matches_molecule": bool(len(rows) == selected.n_records),
        "median_matches": bool(np.isclose(rows.pchembl_value.median(), selected.activity_pchembl)),
        "unique_activity_ids": int(rows.activity_id.nunique()),
        "unique_assays": int(rows.assay_chembl_id.nunique()),
        "unique_documents": int(rows.document_chembl_id.nunique()),
    })


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, default=ROOT / "validation_data")
    parser.add_argument("--output", type=Path, default=ROOT / "paper/output/revision_v4")
    parser.add_argument("--targets", nargs="+", default=PROTOCOL["targets"])
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    dump(out / "protocol.json", PROTOCOL)
    all_baselines, all_seeds, all_sparse, all_stability = [], [], [], []
    inputs, counts, main_data = [], [], None
    for target in args.targets:
        matches = list(args.input_root.glob(f"*{target.lower()}*/*_ic50_raw.csv"))
        if len(matches) != 1:
            raise ValueError(f"Expected one raw CSV for {target}, found {matches}")
        path = matches[0]
        print(f"Curating {target}: {path.name}", flush=True)
        raw = pd.read_csv(path, low_memory=False)
        started = time.perf_counter()
        curated = curate_chembl_activity_data(raw, target_id=target)
        case = out / target.lower()
        save_chembl_curation(curated, case, target.lower())
        data = curated.molecules
        if target != "CHEMBL205" and len(data) > PROTOCOL["external_map_cap"]:
            data = data.sample(PROTOCOL["external_map_cap"], random_state=PROTOCOL["sampling_seed"])
        data = data.sort_values("molecule_identity_key").reset_index(drop=True)
        data.to_csv(case / "map_molecules.csv", index=False)
        counts.append({"target": target, "raw_records": len(raw),
                       "retained_records": len(curated.retained_records),
                       "curated_molecules": len(curated.molecules), "mapped_molecules": len(data),
                       "repeated_molecules": int((curated.molecules.n_records > 1).sum()),
                       "conflicted_molecules": int(curated.molecules.is_conflicted.sum()),
                       "structure_changed_molecules": int(curated.molecules.structure_changed.sum()),
                       "class_crossing_molecules": int(curated.molecules.has_class_boundary_crossing.sum()),
                       "curation_seconds": time.perf_counter() - started,
                       "identity_policy": curated.report["parameters"]["molecule_identity"]})
        inputs.append({"target": target, "file": str(path.relative_to(args.input_root)),
                       "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                       "chembl_status": json.loads((path.parent / "chembl_status.json").read_text())})
        fp = morgan_fingerprints(data.canonical_smiles.tolist())
        seeds = PROTOCOL["main_seeds"] if target == "CHEMBL205" else PROTOCOL["external_seeds"]
        layouts = {}
        for seed in seeds:
            print(f"  {target} UMAP seed {seed}; n={len(data)}", flush=True)
            t0 = time.perf_counter()
            xy = project_2d(fp, method="umap", random_state=seed, umap_metric="jaccard")
            projection_seconds = time.perf_counter() - t0
            diag = {}
            t0 = time.perf_counter()
            grid, side, _ = assign_to_grid(xy, method="dense", diagnostics=diag)
            assignment_seconds = time.perf_counter() - t0
            metrics = evaluate(xy, grid, data)
            base = {"target": target, "seed": seed, "n": len(data)}
            all_seeds.append({**base, **metrics, "projection_seconds": projection_seconds,
                              "assignment_seconds": assignment_seconds})
            all_baselines.append({**base, "method": "Minimum-cost", "repeat": 0, **metrics})
            coords = coordinates(data, xy, grid, side)
            coords.to_csv(case / f"coordinates_seed_{seed}.csv", index=False)
            layouts[seed] = (xy, grid)
            rng = np.random.default_rng(seed)
            for repeat in range(PROTOCOL["baseline_repeats"]):
                baseline, _ = greedy(xy, rng.permutation(len(xy)))
                all_baselines.append({**base, "method": "Greedy", "repeat": repeat,
                                      **evaluate(xy, baseline, data)})
                cells, _, _ = make_candidate_grid(len(xy), target_occupancy=0.4)
                baseline = cells[rng.choice(len(cells), len(xy), replace=False)]
                all_baselines.append({**base, "method": "Random", "repeat": repeat,
                                      **evaluate(xy, baseline, data)})
            for candidates in PROTOCOL["sparse_candidates"]:
                diagnostics = {}
                t0 = time.perf_counter()
                sparse, _, _ = assign_to_grid(xy, method="sparse", sparse_neighbors=candidates,
                                             diagnostics=diagnostics, sparse_adaptive=False)
                seconds = time.perf_counter() - t0
                mm = evaluate(xy, sparse, data)
                all_sparse.append({**base, "nearest_candidates": candidates, **mm,
                                   "assignment_seconds": seconds, **diagnostics,
                                   "dense_cost_ratio": mm["mean_squared_grid_displacement"] /
                                   metrics["mean_squared_grid_displacement"],
                                   "same_cell_fraction": float(np.all(sparse == grid, axis=1).mean())})
            if seed == 42:
                render_svg(coords, case / "molecular_grid.svg")
                render_svg(coords, case / "grid_overview.svg", draw_molecules=False)
                fp_bool = fp.astype(bool)
                dump(case / "upstream_fidelity.json", {
                    "fingerprint_to_projection_trustworthiness": float(trustworthiness(
                        fp_bool, xy, n_neighbors=10, metric="jaccard")),
                    "fingerprint_to_grid_trustworthiness": float(trustworthiness(
                        fp_bool, grid, n_neighbors=10, metric="jaccard")),
                    "n": len(data), "k": 10,
                    "note": "Secondary rank metric; fingerprint and grid distances contain ties."})
                if target == "CHEMBL205":
                    main_data = (data, fp, xy, grid, curated, raw)
                    audit_case(curated, case)
        for a, b in itertools.combinations(seeds, 2):
            all_stability.append({"target": target, "seed_a": a, "seed_b": b,
                                  "projection_neighbor_jaccard": neighborhood_jaccard(layouts[a][0], layouts[b][0]),
                                  "grid_neighbor_jaccard": neighborhood_jaccard(layouts[a][1], layouts[b][1])})
        # Checkpoint tables allow inspection while later targets are still running.
        for name, rows in [("datasets", counts), ("seeds", all_seeds), ("baselines", all_baselines),
                           ("sparse_calibration", all_sparse), ("seed_stability", all_stability)]:
            pd.DataFrame(rows).to_csv(out / f"{name}.csv", index=False)
    if main_data is not None:
        data, fp, xy, grid, curated, raw = main_data
        sensitivity = []
        for occupancy in PROTOCOL["occupancies"]:
            assigned, side, _ = assign_to_grid(xy, target_occupancy=occupancy, method="dense")
            sensitivity.append({"target_occupancy": occupancy, "actual_occupancy": len(data) / side ** 2,
                                "grid_side": side, **evaluate(xy, assigned, data)})
        pd.DataFrame(sensitivity).to_csv(out / "occupancy.csv", index=False)
        runtime = []
        for n in [100, 250, 500, len(data)]:
            for repeat in range(5):
                idx = np.sort(np.random.default_rng(repeat).choice(len(data), n, replace=False))
                t0 = time.perf_counter()
                assign_to_grid(xy[idx], method="dense")
                runtime.append({"n": n, "repeat": repeat, "seconds": time.perf_counter() - t0})
        pd.DataFrame(runtime).to_csv(out / "runtime.csv", index=False)
        original = curate_chembl_activity_data(raw, target_id="CHEMBL205", structure_policy="as-recorded")
        old = original.molecules.set_index("molecule_identity_key").loc[data.molecule_identity_key].reset_index()
        ablation = []
        old_fp = morgan_fingerprints(old.canonical_smiles.tolist())
        for seed in PROTOCOL["main_seeds"]:
            old_xy = project_2d(old_fp, method="umap", random_state=seed, umap_metric="jaccard")
            old_grid, _, _ = assign_to_grid(old_xy, method="dense")
            ablation.append({"seed": seed, "structure_policy": "as-recorded",
                             **evaluate(old_xy, old_grid, old)})
        pd.DataFrame(ablation).to_csv(out / "structure_ablation.csv", index=False)
        dump(out / "structure_audit.json", {
            "common_molecules": len(data),
            "changed_structures": int((old.canonical_smiles.values != data.canonical_smiles.values).sum()),
            "changed_activity_labels": int((old.activity_class.values != data.activity_class.values).sum()),
            "changed_fingerprints": int(np.any(old_fp != fp, axis=1).sum()),
            "old_multifragment_structures": int(old.canonical_smiles.str.contains(".", regex=False).sum()),
            "new_multifragment_structures": int(data.canonical_smiles.str.contains(".", regex=False).sum()),
        })
        alternatives = []
        for representation, method in [("Morgan", "pca"), ("Morgan", "tsne"), ("Descriptors", "umap")]:
            matrix = fp if representation == "Morgan" else rdkit_descriptors(data.canonical_smiles.tolist())
            projected = project_2d(matrix, method=method, random_state=42, umap_metric="euclidean")
            assigned, _, _ = assign_to_grid(projected, method="dense")
            alternatives.append({"representation": representation, "projection": method,
                                 **evaluate(projected, assigned, data)})
        pd.DataFrame(alternatives).to_csv(out / "alternative_inputs.csv", index=False)
    dump(out / "run_manifest.json", {"completed_at": datetime.now(timezone.utc).isoformat(),
                                    "inputs": inputs, "protocol": PROTOCOL,
                                    "environment": environment_versions(),
                                    "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    print(f"Validation complete: {out}", flush=True)


if __name__ == "__main__":
    main()
