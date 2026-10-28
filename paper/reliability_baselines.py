"""Six grid baselines with isolated timing and verifiable checkpoints."""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import lru_cache
import hashlib
import json
import math
import multiprocessing
from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr

from chemgridmap.core import assign_to_grid, normalize_coordinates
from chemgridmap.projection import project_2d
from chemgridmap.refinement import refine_grid_layout
from chemgridmap.reliability import representation_neighbor_sets
from grid_baselines import dpq_scores, flas, las, random_grid_layout
from reliability_analysis import (
    FROZEN, OUTPUT, ROOT, SEEDS, TARGETS, evaluate_neighborhoods,
    jaccard_distances, load_target, sha256, write_manifest, dump_json,
)

METHODS = ("projection_lap", "projection_lap_refined", "las_random",
           "las_projection", "flas_random", "random")
OCCUPANCIES = (0.4, 1.0)
ADDED_METRICS = ("dpq16_minus_tanimoto", "dpq2_minus_tanimoto", "dpq16_minus_euclidean",
                 "spearman_vs_projection", "mean_displacement_vs_projection_cells")
METRIC_COLUMNS = ("recall_k5", "recall_k10", "recall_k20", "trustworthiness_k10",
                  "dpq16_tanimoto", "dpq2_tanimoto", "dpq16_euclidean", "dpq2_euclidean")
PARAMETERS = {
    "methods": METHODS, "occupancy_targets": OCCUPANCIES,
    "feature_space": "Morgan radius=2, bits=2048, L2 normalized for LAS/FLAS",
    "filter": "separable uniform_filter1d, reflect; vacancy mask normalization",
    "empty_filter_window": "global feature mean", "lap_cost": "floating squared Euclidean",
    "radius_decay": 0.9, "random_radius_start": "max(H,W)/2-1",
    "projection_radius_start": "max(H,W)/10", "flas_nc": 49,
    "dpq": {"p": [16, 2], "primary": "Tanimoto/Jaccard",
            "secondary": "L2-normalized Morgan Euclidean",
            "tie_policy": "sorted: optimistic HD ordering; mean: pessimistic HD mean within exact 2D ties"},
    "runtime_seed": 42, "runtime_scope": "layout construction only, fingerprint preparation excluded",
    "runtime_refined": "includes representation neighbors, projection, LAP, refinement",
    "runtime_las_projection": "includes projection, LAP, LAS",
    "runtime_processes": "seed42 methods sequential, each in a fresh single-process worker",
    "other_seed_runtime": "not eligible for reported timings; cached projection initialization",
    "sources": ["doi:10.1111/cgf.14718", "https://github.com/Visual-Computing/LAS_FLAS"],
}


@lru_cache(maxsize=3)
def _prepared(target):
    data, fp = load_target(target)
    vectors = np.asarray(fp, dtype=float)
    norms = np.linalg.norm(vectors, axis=1, keepdims=True)
    vectors = np.divide(vectors, norms, out=np.zeros_like(vectors), where=norms != 0)
    high = jaccard_distances(fp)
    sparse = csr_matrix(vectors)
    sqnorm = np.einsum("ij,ij->i", vectors, vectors)
    euclidean = np.sqrt(np.maximum(sqnorm[:, None] + sqnorm[None, :] -
                                  2 * (sparse @ sparse.T).toarray(), 0))
    np.fill_diagonal(euclidean, 0)
    references = {k: representation_neighbor_sets(fp, k=k) for k in (5, 10, 20)}
    return data, fp, vectors, high, euclidean, references


def _measure_method(payload):
    target, seed, occupancy, method, output_root, timed = payload
    data, fp, vectors, high, euclidean, references = _prepared(target)
    layout_path = Path(output_root) / "layouts" / target / f"coordinates_seed_{seed}.csv"
    layout = pd.read_csv(layout_path)
    if data.molecule_identity_key.tolist() != layout.molecule_identity_key.tolist():
        raise ValueError("Baseline input identities do not match the layout order.")
    side = math.ceil(math.sqrt(len(data) / occupancy))
    projection_seconds = assignment_seconds = neighbors_seconds = algorithm_seconds = 0.0
    start = time.perf_counter()
    if method in ("projection_lap", "projection_lap_refined", "las_projection"):
        if timed:
            t0 = time.perf_counter()
            xy = project_2d(fp, method="umap", random_state=seed, umap_metric="jaccard")
            projection_seconds = time.perf_counter() - t0
        else:
            # Restore UMAP float32 precision so CSV loading does not change LAP ties.
            xy = layout[["projection_x_raw", "projection_y_raw"]].to_numpy(dtype=np.float32)
        t0 = time.perf_counter()
        grid, side, _ = assign_to_grid(xy, target_occupancy=occupancy, method="dense")
        assignment_seconds = time.perf_counter() - t0
        cells = np.rint(grid * (side - 1)).astype(int)
        if method == "projection_lap":
            algorithm_seconds = assignment_seconds
        elif method == "projection_lap_refined":
            # Time method-specific neighbor preparation separately from evaluation.
            t0 = time.perf_counter()
            neighbors = representation_neighbor_sets(fp, k=10)
            neighbors_seconds = time.perf_counter() - t0
            t0 = time.perf_counter()
            cells = refine_grid_layout(cells, normalize_coordinates(xy) * (side - 1),
                                       neighbors, side, random_state=seed)
            algorithm_seconds = time.perf_counter() - t0
        else:
            t0 = time.perf_counter()
            cells = las(vectors, (side, side), seed=seed, initial_cells=cells)
            algorithm_seconds = time.perf_counter() - t0
    else:
        t0 = time.perf_counter()
        if method == "las_random":
            cells = las(vectors, (side, side), seed=seed)
        elif method == "flas_random":
            cells = flas(vectors, (side, side), seed=seed)
        elif method == "random":
            cells = random_grid_layout(len(data), (side, side), seed=seed)
        else:
            raise ValueError(f"Unknown baseline: {method}")
        algorithm_seconds = time.perf_counter() - t0
    total_seconds = time.perf_counter() - start
    # Validate and evaluate after timing to exclude metric computation.
    if timed and method in ("projection_lap", "projection_lap_refined", "las_projection"):
        if not np.allclose(xy, layout[["projection_x_raw", "projection_y_raw"]], atol=1e-5, rtol=0):
            raise ValueError(f"Timed projection did not reproduce: {target} seed={seed}")
    if occupancy == .4 and method in ("projection_lap", "projection_lap_refined"):
        prefix = "refined_" if method == "projection_lap_refined" else ""
        np.testing.assert_array_equal(cells, layout[[prefix + "grid_col", prefix + "grid_row"]])
    metrics, _, _ = evaluate_neighborhoods(fp, cells, references, high)
    primary, secondary = dpq_scores(high, cells), dpq_scores(euclidean, cells)
    row = {"target": target, "target_label": TARGETS[target], "seed": seed,
           "method": method, "target_occupancy": occupancy,
           "actual_occupancy": len(data) / side ** 2, "n": len(data), "grid_side": side,
           **metrics, "dpq16_tanimoto": primary[16], "dpq2_tanimoto": primary[2],
           "dpq16_euclidean": secondary[16], "dpq2_euclidean": secondary[2],
           **_additional_metrics(high, euclidean, cells, layout, side),
           "runtime_reportable": bool(timed),
           "projection_recomputed": bool(timed and method in
               ("projection_lap", "projection_lap_refined", "las_projection")),
           "projection_seconds": projection_seconds, "assignment_seconds": assignment_seconds,
           "representation_neighbors_seconds": neighbors_seconds,
           "algorithm_seconds": algorithm_seconds, "total_seconds": total_seconds,
           "runtime_scope": "fresh_full_layout" if timed else "cached_projection_layout"}
    print(f"baselines {target} occupancy={occupancy} seed={seed} {method}: "
          f"recall10={row['recall_k10']:.4f}, DPQ16={primary[16]:.4f}, {total_seconds:.1f}s", flush=True)
    return row, cells


def _checkpoint_paths(output_root, target, occupancy, seed):
    parent = Path(output_root) / "baselines" / "checkpoints"
    stem = f"{target}_occupancy{occupancy:.1f}_seed{seed}"
    return parent / f"{stem}.csv", parent / f"{stem}.npz", parent / f"{stem}.json"


def _cache_key(output_root, target, occupancy, seed):
    inputs = [FROZEN / target / "map_molecules.csv",
              Path(output_root) / "layouts" / target / f"coordinates_seed_{seed}.csv"]
    sources = [Path(__file__), Path(__file__).with_name("grid_baselines.py"),
               Path(__file__).with_name("reliability_analysis.py")]
    sources.extend(sorted((ROOT / "src").rglob("*.py")))
    payload = {"target": target, "occupancy": occupancy, "seed": seed,
               "parameters": PARAMETERS, "files": {str(p): sha256(p) for p in inputs + sources}}
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def _cached_group(output_root, target, occupancy, seed, cache_key):
    csv_path, npz_path, metadata_path = _checkpoint_paths(output_root, target, occupancy, seed)
    if not all(p.exists() for p in (csv_path, npz_path, metadata_path)):
        return None
    metadata = json.loads(metadata_path.read_text())
    if metadata.get("cache_key") != cache_key or metadata.get("csv_sha256") != sha256(csv_path) \
            or metadata.get("coordinates_sha256") != sha256(npz_path):
        return None
    frame = pd.read_csv(csv_path)
    # Frozen checkpoint labels are translated in memory after their hashes pass.
    frame = frame.replace({"target_label": {"D3": "Thrombin"}})
    if len(frame) != len(METHODS) or set(frame.method) != set(METHODS):
        return None
    return frame


def _save_group(output_root, target, occupancy, seed, cache_key, results):
    csv_path, npz_path, metadata_path = _checkpoint_paths(output_root, target, occupancy, seed)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame([row for row, _ in results])
    frame.to_csv(csv_path, index=False)
    data, _ = load_target(target)
    np.savez_compressed(npz_path, molecule_identity_key=data.molecule_identity_key.to_numpy(str),
                        **{row["method"]: cells for row, cells in results})
    metadata = {"cache_key": cache_key, "csv_sha256": sha256(csv_path),
                "coordinates_sha256": sha256(npz_path)}
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    return frame


def _run_group(payload):
    output_root, target, occupancy, seed, cache_key = payload
    results = [_measure_method((target, seed, occupancy, method, output_root, False))
               for method in METHODS]
    return _save_group(output_root, target, occupancy, seed, cache_key, results)


def _additional_metrics(high, euclidean, cells, layout, side):
    """Compute pessimistic DPQ and correspondence to the saved projection."""
    primary = dpq_scores(high, cells, ties="mean")
    secondary = dpq_scores(euclidean, cells, powers=(16,), ties="mean")
    projection = layout[["projection_x", "projection_y"]].to_numpy(float)
    return {"dpq16_minus_tanimoto": primary[16], "dpq2_minus_tanimoto": primary[2],
            "dpq16_minus_euclidean": secondary[16],
            "spearman_vs_projection": float(spearmanr(pdist(cells), pdist(projection)).statistic),
            "mean_displacement_vs_projection_cells": float(
                np.linalg.norm(cells - projection * (side - 1), axis=1).mean())}


def supplement_baseline_metrics(output_root=OUTPUT, *, seeds=SEEDS, targets=TARGETS):
    """Append metrics without constructing layouts, timing, or writing checkpoints."""
    root = Path(output_root)
    out = root / "baselines"
    metric_path, manifest_path = out / "baseline_metrics.csv", out / "manifest.json"
    previous_manifest = json.loads(manifest_path.read_text())
    previous_hash = sha256(manifest_path)
    previous_metrics_hash = sha256(metric_path)
    inputs = [Path(p) for p in previous_manifest["inputs_sha256"]]
    for path in inputs:
        if sha256(path) != previous_manifest["inputs_sha256"][str(path)]:
            raise ValueError(f"Baseline input hash changed: {path}")
    frame = pd.read_csv(metric_path, float_precision="round_trip")
    old_columns = [column for column in frame if column not in ADDED_METRICS]
    original = frame[old_columns].copy(deep=True)
    keys = ["target", "target_occupancy", "seed", "method"]
    expected = {(target, occupancy, seed, method) for target in targets
                for occupancy in OCCUPANCIES for seed in seeds for method in METHODS}
    if len(frame) != len(expected) or set(frame[keys].itertuples(index=False, name=None)) != expected:
        raise ValueError("Baseline rows do not match the requested complete design.")
    audits, additions = [], {}
    for (target, occupancy, seed), group in frame.groupby(keys[:-1], sort=False):
        csv_path, npz_path, metadata_path = _checkpoint_paths(root, target, occupancy, seed)
        metadata = json.loads(metadata_path.read_text())
        if metadata["csv_sha256"] != sha256(csv_path) or metadata["coordinates_sha256"] != sha256(npz_path):
            raise ValueError(f"Checkpoint hash mismatch: {metadata_path}")
        saved = pd.read_csv(csv_path, float_precision="round_trip")
        # Compare current labels while preserving every frozen checkpoint byte.
        saved = saved.replace({"target_label": {"D3": "Thrombin"}})
        pd.testing.assert_frame_equal(
            group[old_columns].sort_values("method").reset_index(drop=True),
            saved[old_columns].sort_values("method").reset_index(drop=True),
            check_exact=False, atol=1e-14, rtol=0)
        layout_path = root / "layouts" / target / f"coordinates_seed_{seed}.csv"
        layout = pd.read_csv(layout_path)
        data, _, _, high, euclidean, _ = _prepared(target)
        identities = data.molecule_identity_key.to_numpy(str)
        np.testing.assert_array_equal(identities, layout.molecule_identity_key.to_numpy(str))
        inputs.extend([csv_path, npz_path, metadata_path, layout_path, FROZEN / target / "map_molecules.csv"])
        with np.load(npz_path, allow_pickle=False) as checkpoint:
            np.testing.assert_array_equal(identities, checkpoint["molecule_identity_key"])
            for index, row in group.iterrows():
                cells = checkpoint[row.method]
                side = int(row.grid_side)
                if cells.shape != (len(data), 2) or not np.isfinite(cells).all() or \
                        not np.equal(cells, np.rint(cells)).all() or \
                        np.any(cells < 0) or np.any(cells >= side) or len(np.unique(cells, axis=0)) != len(data):
                    raise ValueError(f"Invalid checkpoint coordinates: {npz_path}: {row.method}")
                recomputed = dpq_scores(high, cells, powers=(16,))[16]
                difference = abs(recomputed - row.dpq16_tanimoto)
                if not np.isclose(recomputed, row.dpq16_tanimoto, atol=1e-12, rtol=0):
                    raise ValueError(f"Checkpoint DPQ mismatch: {target}/{occupancy}/{seed}/{row.method}")
                audits.append({**{key: row[key] for key in keys},
                               "stored_dpq16_tanimoto": row.dpq16_tanimoto,
                               "recomputed_dpq16_tanimoto": recomputed, "absolute_difference": difference})
                additions[index] = _additional_metrics(high, euclidean, cells, layout, side)
        print(f"baseline metrics-only: {target}/{occupancy}/{seed} verified and evaluated", flush=True)
    for column in ADDED_METRICS:
        frame[column] = [additions[index][column] for index in frame.index]
    if not np.isfinite(frame[list(ADDED_METRICS)].to_numpy()).all():
        raise ValueError("Supplemental baseline metrics must be finite.")
    pd.testing.assert_frame_equal(frame[old_columns], original, check_exact=True)
    # Preserve existing summary values rather than re-rounding their inputs.
    summary = pd.read_csv(out / "summary.csv", header=[0, 1], index_col=[0, 1, 2],
                          float_precision="round_trip")
    extra_summary = frame.groupby(keys[:2] + ["method"])[list(ADDED_METRICS)].agg(["mean", "min", "max"])
    pd.testing.assert_index_equal(summary.index, extra_summary.index)
    for column in extra_summary:
        summary[column] = extra_summary[column]
    # Write only after every checkpoint and all old row values have been verified.
    frame.to_csv(metric_path, index=False)
    summary.to_csv(out / "summary.csv")
    pd.DataFrame(audits).to_csv(out / "checkpoint_dpq_verification.csv", index=False)
    audit = {"passed": True, "rows": len(frame), "checkpoint_groups": len(audits) // len(METHODS),
             "checkpoint_dpq16_max_absolute_difference": max(row["absolute_difference"] for row in audits),
             "original_columns_unchanged": old_columns, "added_columns": list(ADDED_METRICS),
             "layouts_recomputed": False, "timings_recomputed": False,
             "previous_baseline_metrics_sha256": previous_metrics_hash,
             "previous_manifest_sha256": previous_hash,
             "layout_generation_provenance": previous_manifest.get("layout_generation_provenance", {
                 "code": previous_manifest["code"], "parameters": previous_manifest["parameters"]})}
    dump_json(out / "metrics_only_verification.json", audit)
    manifest = write_manifest(out, inputs, {**PARAMETERS, "seeds": list(seeds), "targets": list(targets),
        "phase": "metrics-only", "completed_rows": len(frame), "expected_rows": len(expected),
        "projection_correspondence": "Spearman of all unique pairs; displacement from normalized projection * (grid_side - 1)",
        "checkpoint_verification": "metadata hashes, identities, old columns, sorted Tanimoto DPQ16 (absolute tolerance 1e-12)",
        "layouts_recomputed": False, "timings_recomputed": False})
    manifest["layout_generation_provenance"] = audit["layout_generation_provenance"]
    dump_json(manifest_path, manifest)
    return frame


def run_baselines(output_root=OUTPUT, *, seeds=SEEDS, targets=TARGETS, workers=1,
                  phase="all", resume=False):
    """Run isolated timings first, or append metrics from read-only checkpoints."""
    if phase == "metrics-only":
        return supplement_baseline_metrics(output_root, seeds=seeds, targets=targets)
    if phase not in ("all", "timing", "remaining"):
        raise ValueError("phase must be all, timing, remaining, or metrics-only.")
    output_root = Path(output_root)
    out = output_root / "baselines"
    out.mkdir(parents=True, exist_ok=True)
    jobs, frames = [], []
    for seed in sorted(seeds, key=lambda value: (value != 42, value)):
        if (phase == "timing" and seed != 42) or (phase == "remaining" and seed == 42):
            continue
        for target in targets:
            for occupancy in OCCUPANCIES:
                key = _cache_key(output_root, target, occupancy, seed)
                cached = _cached_group(output_root, target, occupancy, seed, key) if resume else None
                if cached is not None:
                    frames.append(cached)
                    print(f"baselines validated checkpoint: {target}/{occupancy}/{seed}", flush=True)
                elif seed == 42:
                    # Run methods sequentially, each in a fresh process.
                    results = []
                    with ProcessPoolExecutor(max_workers=1, max_tasks_per_child=1,
                                             mp_context=multiprocessing.get_context("spawn")) as executor:
                        for method in METHODS:
                            results.append(executor.submit(_measure_method,
                                (target, seed, occupancy, method, output_root, True)).result())
                    frames.append(_save_group(output_root, target, occupancy, seed, key, results))
                else:
                    jobs.append((output_root, target, occupancy, seed, key))
    if workers > 1 and jobs:
        with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("spawn")) as executor:
            futures = [executor.submit(_run_group, job) for job in jobs]
            for future in as_completed(futures):
                frames.append(future.result())
    else:
        frames.extend(_run_group(job) for job in jobs)
    # Merge only valid checkpoints after the two construction phases.
    frames = []
    inputs = []
    for target in targets:
        inputs.append(FROZEN / target / "map_molecules.csv")
        for seed in seeds:
            layout_path = output_root / "layouts" / target / f"coordinates_seed_{seed}.csv"
            if not layout_path.exists():
                continue
            inputs.append(layout_path)
            for occupancy in OCCUPANCIES:
                key = _cache_key(output_root, target, occupancy, seed)
                cached = _cached_group(output_root, target, occupancy, seed, key)
                if cached is not None:
                    frames.append(cached)
    if frames:
        frame = pd.concat(frames, ignore_index=True).sort_values(["target", "target_occupancy", "seed", "method"])
        frame.to_csv(out / "baseline_metrics.csv", index=False)
        metric_columns = list(METRIC_COLUMNS + ADDED_METRICS)
        frame.groupby(["target", "target_occupancy", "method"])[metric_columns].agg(
            ["mean", "min", "max"]).to_csv(out / "summary.csv")
        frame[frame.runtime_reportable].to_csv(out / "runtime_seed42.csv", index=False)
    write_manifest(out, inputs, {**PARAMETERS, "seeds": list(seeds), "targets": list(targets),
                                "phase": phase, "workers_other_seeds": workers,
                                "completed_rows": sum(len(f) for f in frames),
                                "expected_rows": len(seeds) * len(targets) * len(OCCUPANCIES) * len(METHODS)})
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
