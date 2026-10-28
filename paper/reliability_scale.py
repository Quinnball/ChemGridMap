"""Measure the full default CLI and refinement of the frozen hERG grid."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import re
import shlex
import subprocess
import sys
import tempfile
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUTPUT_ROOT = ROOT / "paper/output/reliability_v1"
RAW_INPUT = ROOT / "validation_data/chembl37_chembl240_large_scale/chembl37_chembl240_ic50_raw.csv"
FROZEN_GRID = ROOT / "paper/output/revision_v6/chembl240/full_grid_coordinates_cap2048.csv"
PRIOR_RESULT = ROOT.parent / (
    "project1/JCIM_review_innovation_20261024/results/grid_refinement/full_herg_scale_check.json")
SCALE_INPUTS = (RAW_INPUT, FROZEN_GRID, PRIOR_RESULT)
NAME = "full_herg_default"
CLI_RESOURCE_RUNNER = """import json, pathlib, resource, runpy, sys
resource_path = pathlib.Path(sys.argv.pop(1))
exit_code = 1
try:
    runpy.run_module('chemgridmap', run_name='__main__', alter_sys=True)
    exit_code = 0
except SystemExit as exc:
    exit_code = exc.code if isinstance(exc.code, int) else (0 if exc.code is None else 1)
    raise
finally:
    usage = resource.getrusage(resource.RUSAGE_SELF)
    resource_path.write_text(json.dumps({'cli_exit_code': exit_code,
        'peak_rss_bytes': usage.ru_maxrss, 'user_seconds': usage.ru_utime,
        'system_seconds': usage.ru_stime, 'rss_source': 'resource.getrusage(RUSAGE_SELF) on macOS'}))
"""


def dump(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def parse_time_l(text: str) -> dict:
    """Read macOS /usr/bin/time -l; maximum resident set size is in bytes."""
    times = re.search(r"([\d.]+)\s+real\s+([\d.]+)\s+user\s+([\d.]+)\s+sys", text)
    resident = re.search(r"(\d+)\s+maximum resident set size", text)
    if times is None:
        raise ValueError("Missing /usr/bin/time -l wall-time output.")
    return {"wall_seconds": float(times[1]), "user_seconds": float(times[2]),
            "system_seconds": float(times[3]),
            "peak_rss_bytes": int(resident[1]) if resident else None}


def cli_command(output: Path) -> list[str]:
    # runpy executes the real CLI; the wrapper only records process resources and
    # its actual exit status when macOS time -l cannot read kern.clockrate.
    return ["/usr/bin/time", "-l", str(ROOT / ".venv/bin/python"), "-B", "-c", CLI_RESOURCE_RUNNER,
            str(output / "cli_process_resource.json"),
            str(RAW_INPUT), "--input-format", "chembl", "--target-id", "CHEMBL240",
            "--projection", "umap", "--random-state", "42", "--assignment-method", "sparse",
            "--sparse-max-neighbors", "2048", "--output-dir", str(output / "cli"), "--name", NAME]


def _run_cli(output: Path) -> tuple[dict, pd.DataFrame, dict]:
    command = cli_command(output)
    environment = dict(os.environ)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    environment.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "chemgridmap_t2_mpl"))
    invocation = {"argv": command, "command": shlex.join(command), "cwd": str(ROOT),
                  "started_utc": datetime.now(timezone.utc).isoformat(),
                  "platform": platform.platform(), "processor": platform.processor(),
                  "timing_scope": "fresh CLI process, including curation, class support and all default exports",
                  "environment": {key: environment.get(key) for key in
                                  ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                                   "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS", "MPLCONFIGDIR",
                                   "PYTHONDONTWRITEBYTECODE")}}
    dump(output / "cli_invocation.json", invocation)
    print("scale: full default hERG CLI started", flush=True)
    (output / "cli_process_resource.json").unlink(missing_ok=True)
    start = time.perf_counter()
    with (output / "cli_stdout.log").open("w") as stdout, (output / "cli_stderr_time.log").open("w") as stderr:
        result = subprocess.run(command, cwd=ROOT, env=environment, stdout=stdout, stderr=stderr)
    invocation.update(returncode=result.returncode, wrapper_wall_seconds=time.perf_counter() - start,
                      finished_utc=datetime.now(timezone.utc).isoformat())
    dump(output / "cli_invocation.json", invocation)
    process_usage_path = output / "cli_process_resource.json"
    process_usage = json.loads(process_usage_path.read_text()) if process_usage_path.exists() else {}
    stderr_text = (output / "cli_stderr_time.log").read_text()
    time_blocked = "time: sysctl kern.clockrate: Operation not permitted" in stderr_text
    invocation.update(cli_exit_code=process_usage.get("cli_exit_code"),
                      time_l_memory_status="sandbox-blocked" if time_blocked else "measured")
    dump(output / "cli_invocation.json", invocation)
    if process_usage.get("cli_exit_code") != 0 or (result.returncode and not time_blocked):
        raise RuntimeError(f"Full hERG CLI exited {result.returncode}; see {output / 'cli_stderr_time.log'}")
    return _load_cli_outputs(output)


def _load_cli_outputs(output: Path) -> tuple[dict, pd.DataFrame, dict]:
    invocation = json.loads((output / "cli_invocation.json").read_text())
    process_usage = json.loads((output / "cli_process_resource.json").read_text())
    stderr_text = (output / "cli_stderr_time.log").read_text()
    if process_usage["cli_exit_code"] != 0:
        raise AssertionError("Cannot summarize an unsuccessful CLI process.")
    resource = parse_time_l(stderr_text)
    resource["time_l_real_seconds"] = resource["wall_seconds"]
    resource["wall_seconds"] = invocation["wrapper_wall_seconds"]
    resource["wall_time_source"] = "parent time.perf_counter; cross-checked against time -l real"
    resource["time_l_returncode"] = invocation["returncode"]
    resource["cli_exit_code"] = process_usage["cli_exit_code"]
    resource["time_l_memory_status"] = invocation["time_l_memory_status"]
    resource["peak_rss_source"] = "/usr/bin/time -l"
    if resource["peak_rss_bytes"] is None:
        resource["peak_rss_bytes"] = process_usage["peak_rss_bytes"]
        resource["peak_rss_source"] = process_usage["rss_source"]
    resource["peak_rss_mib"] = resource["peak_rss_bytes"] / 1024 ** 2
    manifest = json.loads((output / "cli" / f"{NAME}_run_manifest.json").read_text())
    data = pd.read_csv(output / "cli" / manifest["outputs"]["coordinates"], low_memory=False)
    parameters = manifest["parameters"]
    for key, expected in {"random_state": 42, "assignment_method": "sparse", "sparse_max_neighbors": 2048,
                          "grid_refinement": "neighbors", "core_rule": "coverage",
                          "core_min_coverage": 0.8, "activity_confidence": True,
                          "projection": "umap", "umap_metric": "jaccard"}.items():
        if parameters[key] != expected:
            raise AssertionError(f"Unexpected CLI parameter {key}: {parameters[key]!r}")
    if len(data) != 9581 or manifest["activity_confidence"] is None:
        raise AssertionError("Full hERG CLI must produce 9,581 entries with class support.")
    if manifest["structure_grouping"]["method"] != "RDKit Butina":
        raise AssertionError("Full default CLI must retain Butina grouping.")
    print(f"scale: full CLI completed in {resource['wall_seconds']:.2f} s", flush=True)
    return resource, data, manifest


def _refine_frozen(output: Path) -> tuple[dict, pd.DataFrame, list[np.ndarray]]:
    from chemgridmap.refinement import REFINEMENT_DEFAULTS, refine_grid_layout
    from chemgridmap.reliability import grid_neighbor_sets, neighbor_recall, representation_neighbor_sets
    from chemgridmap.representations import morgan_fingerprints

    data = pd.read_csv(FROZEN_GRID, low_memory=False)
    spacing = np.diff(np.unique(data.grid_x.to_numpy(float)))
    side = int(np.rint(1 / spacing.min())) + 1
    cells = data[["grid_col", "grid_row"]].to_numpy(int)
    if len(data) != 9581 or not np.array_equal(
            cells, np.rint(data[["grid_x", "grid_y"]].to_numpy(float) * (side - 1)).astype(int)):
        raise AssertionError("Frozen grid count or integer coordinate conversion changed.")
    projection = data[["projection_x", "projection_y"]].to_numpy(float) * (side - 1)
    start = time.perf_counter()
    fingerprints = morgan_fingerprints(data.canonical_smiles.tolist())
    reference = representation_neighbor_sets(fingerprints, k=10, metric="jaccard")
    fingerprint_seconds = time.perf_counter() - start
    before = neighbor_recall(reference, grid_neighbor_sets(cells, k=10))
    diagnostics = {}
    print("scale: frozen hERG refinement started", flush=True)
    start = time.perf_counter()
    refined = refine_grid_layout(cells, projection, reference, side, random_state=42, diagnostics=diagnostics)
    refinement_seconds = time.perf_counter() - start
    after = neighbor_recall(reference, grid_neighbor_sets(refined, k=10))
    prior = json.loads(PRIOR_RESULT.read_text())
    result = {"n": len(data), "side": side, "random_state": 42, "k": 10,
              "refinement_defaults": REFINEMENT_DEFAULTS,
              "fingerprint_neighbor_seconds": fingerprint_seconds,
              "refinement_seconds": refinement_seconds,
              "recall_minimum_displacement": float(before.mean()),
              "recall_neighbor_refined": float(after.mean()),
              "moved_molecules": diagnostics["moved_molecules"], "sweeps": len(diagnostics["sweeps"]),
              "mean_displacement_cells_before": float(np.linalg.norm(projection - cells, axis=1).mean()),
              "mean_displacement_cells_after": float(np.linalg.norm(projection - refined, axis=1).mean()),
              "prior_working_copy": prior,
              "recall_before_difference_from_prior": float(before.mean() - prior["recall_minimum_displacement"]),
              "recall_after_difference_from_prior": float(after.mean() - prior["recall_neighbor_refined"]),
              "diagnostics": diagnostics}
    coordinates = data[["molecule_identity_key", "molecule_id", "canonical_smiles"]].copy()
    coordinates[["minimum_grid_col", "minimum_grid_row"]] = cells
    coordinates[["refined_grid_col", "refined_grid_row"]] = refined
    coordinates["minimum_neighbor_recall_k10"] = before
    coordinates["refined_neighbor_recall_k10"] = after
    coordinates.to_csv(output / "frozen_refinement_coordinates.csv", index=False)
    dump(output / "frozen_refinement.json", result)
    print(f"scale: frozen refinement completed in {refinement_seconds:.2f} s", flush=True)
    return result, data, reference


def _structure_comparison(cli_data, frozen_data) -> dict:
    from rdkit import Chem
    from chemgridmap.representations import morgan_fingerprints

    differences = []
    for i in np.flatnonzero(cli_data.canonical_smiles.ne(frozen_data.canonical_smiles)):
        strings = [frozen_data.canonical_smiles.iloc[i], cli_data.canonical_smiles.iloc[i]]
        molecules = [Chem.MolFromSmiles(smiles) for smiles in strings]
        fp = morgan_fingerprints(strings).astype(bool)
        canonicalization = []
        for smiles in strings:
            history = [smiles]
            for _ in range(3):
                history.append(Chem.MolToSmiles(Chem.MolFromSmiles(history[-1])))
            canonicalization.append(history)
        differences.append({"molecule_identity_key": cli_data.molecule_identity_key.iloc[i],
                            "molecule_id": cli_data.molecule_id.iloc[i],
                            "frozen_canonical_smiles": strings[0], "cli_canonical_smiles": strings[1],
                            "mutual_substructure_match_with_chirality": bool(
                                molecules[0].HasSubstructMatch(molecules[1], useChirality=True) and
                                molecules[1].HasSubstructMatch(molecules[0], useChirality=True)),
                            "morgan_fingerprint_equal": bool(np.array_equal(fp[0], fp[1])),
                            "morgan_tanimoto": float((fp[0] & fp[1]).sum() / (fp[0] | fp[1]).sum()),
                            "repeated_canonicalization_from_frozen": canonicalization[0],
                            "repeated_canonicalization_from_cli": canonicalization[1]})
    return {"aligned_identity_keys_equal": bool(cli_data.molecule_identity_key.equals(frozen_data.molecule_identity_key)),
            "aligned_molecule_ids_equal": bool(cli_data.molecule_id.equals(frozen_data.molecule_id)),
            "n_entries": len(cli_data), "different_smiles_entries": len(differences),
            "mismatched_count": len(differences),
            "mismatched_molecule_ids": [row["molecule_id"] for row in differences],
            "different_fingerprint_entries": sum(not row["morgan_fingerprint_equal"] for row in differences),
            "differences": differences,
            "interpretation": "CLI and frozen-layout results are not a same-input paired baseline. Each is evaluated with its own exported structures and fingerprints; the observed non-idempotent normalization is a software/input-normalization defect candidate, not an equivalent-SMILES substitution."}


def _summarize_scale(output, resource, cli_data, cli_manifest, frozen, frozen_data) -> dict:
    from chemgridmap.reliability import grid_neighbor_sets, neighbor_recall, representation_neighbor_sets
    from chemgridmap.representations import morgan_fingerprints

    ordered = cli_data.set_index("molecule_identity_key").loc[frozen_data.molecule_identity_key].reset_index()
    comparison = _structure_comparison(ordered, frozen_data)
    comparison["original_identity_order_equal"] = bool(cli_data.molecule_identity_key.equals(frozen_data.molecule_identity_key))
    dump(output / "structure_comparison.json", comparison)
    # The real CLI prepares structures again. Its fingerprint reference must be
    # rebuilt from its own exported strings, even when entry identities match.
    cli_reference = representation_neighbor_sets(morgan_fingerprints(ordered.canonical_smiles.tolist()), k=10)
    cli_recall = neighbor_recall(cli_reference, grid_neighbor_sets(ordered[["grid_col", "grid_row"]].to_numpy(), k=10))
    recall_error = float(np.abs(cli_recall - ordered.neighbor_recall.to_numpy()).max())
    if recall_error > 1e-12:
        raise AssertionError(f"CLI recall independent recomputation differs by {recall_error}.")
    audit = ordered[["molecule_identity_key", "molecule_id"]].copy()
    audit["exported_neighbor_recall_k10"] = ordered.neighbor_recall.to_numpy()
    audit["recomputed_neighbor_recall_k10"] = cli_recall
    audit["absolute_difference"] = np.abs(cli_recall - ordered.neighbor_recall.to_numpy())
    audit.to_csv(output / "cli_recall_audit.csv", index=False)
    timings = cli_manifest["timings_seconds"]
    map_uninstrumented = timings["total_map_workflow"] - sum(
        value for key, value in timings.items() if key != "total_map_workflow")
    external_overhead = resource["wall_seconds"] - timings["total_map_workflow"]
    summary = {"target": "chembl240", "n_entries": len(cli_data), "seed": 42, "k": 10,
               "inputs": [str(path) for path in SCALE_INPUTS],
               "cli": {**resource, "recall_k10": float(cli_recall.mean()),
                       "timings_seconds": timings,
                       "map_uninstrumented_seconds": map_uninstrumented,
                       "outside_map_workflow_seconds": external_overhead,
                       "recall_recomputation_max_absolute_difference": recall_error,
                       "class_support_nonmissing_entries": int(cli_data.class_confidence.notna().sum()),
                       "structure_grouping": cli_manifest["structure_grouping"],
                       "parameters": cli_manifest["parameters"]},
               "frozen_refinement": {key: value for key, value in frozen.items() if key != "diagnostics"},
               "structure_comparison": comparison,
               "cli_vs_frozen_projection_max_absolute_difference": float(np.abs(
                   ordered[["projection_x", "projection_y"]].to_numpy() -
                   frozen_data[["projection_x", "projection_y"]].to_numpy()).max()),
               "timing_notes": [
                   "CLI wall time uses parent perf_counter, cross-checked with /usr/bin/time -l real.",
                   "Peak RSS uses /usr/bin/time -l when available; otherwise macOS resource.getrusage(RUSAGE_SELF), in bytes.",
                   "If time -l cannot read kern.clockrate in the sandbox, its exit 1 and stderr are retained separately from the CLI exit code.",
                   "CLI defaults include SVG/PNG/PDF, coverage cores, Butina groups and class support.",
                   "The manifest does not time grouping separately; map_uninstrumented includes grouping and other gaps.",
                   "outside_map_workflow includes startup, input/curation/class support, curation export and ZIP work.",
                   "Frozen refinement_seconds excludes fingerprint creation, neighbor search and recall evaluation."]}
    dump(output / "full_herg_summary.json", summary)
    rows = [{"target_id": "CHEMBL240", "n_entries": len(cli_data), "random_state": 42,
             "analysis": "cli_default_full_workflow", "input_layout": "raw_activity_csv",
             "recall_before_k10": None, "recall_after_k10": float(cli_recall.mean()),
             "wall_seconds": resource["wall_seconds"], "timing_scope": "full_cli",
             "peak_rss_bytes": resource["peak_rss_bytes"], "peak_rss_mib": resource["peak_rss_mib"]},
            {"target_id": "CHEMBL240", "n_entries": frozen["n"], "random_state": 42,
             "analysis": "frozen_layout_refinement", "input_layout": FROZEN_GRID.name,
             "recall_before_k10": frozen["recall_minimum_displacement"],
             "recall_after_k10": frozen["recall_neighbor_refined"],
             "wall_seconds": frozen["refinement_seconds"], "timing_scope": "refine_grid_layout_only",
             "peak_rss_bytes": None, "peak_rss_mib": None}]
    pd.DataFrame(rows).to_csv(output / "full_herg_scale.csv", index=False)
    stage_rows = [{"target_id": "CHEMBL240", "analysis": "cli_default_full_workflow",
                   "stage": key, "seconds": value, "source": "run_manifest"}
                  for key, value in timings.items()]
    stage_rows += [{"target_id": "CHEMBL240", "analysis": "cli_default_full_workflow",
                    "stage": key, "seconds": value, "source": source}
                   for key, value, source in [
                       ("full_cli_wall", resource["wall_seconds"], "parent time.perf_counter"),
                       ("time_l_real", resource["time_l_real_seconds"], "/usr/bin/time -l"),
                       ("map_uninstrumented", map_uninstrumented, "difference_not_isolated_stage"),
                       ("outside_map_workflow", external_overhead, "difference_not_isolated_stage")]]
    pd.DataFrame(stage_rows).to_csv(output / "full_herg_stage_timings.csv", index=False)
    return summary


def run_scale(output_root: Path | None = None) -> dict:
    """Run sequentially; the caller must avoid concurrent benchmark workloads."""
    output = (OUTPUT_ROOT if output_root is None else Path(output_root)) / "scale"
    output.mkdir(parents=True, exist_ok=True)
    resource, cli_data, cli_manifest = _run_cli(output)
    frozen, frozen_data, _ = _refine_frozen(output)
    return _summarize_scale(output, resource, cli_data, cli_manifest, frozen, frozen_data)


def summarize_scale(output_root: Path | None = None) -> dict:
    """Rebuild postprocessing from measured outputs without repeating benchmarks."""
    output = (OUTPUT_ROOT if output_root is None else Path(output_root)) / "scale"
    resource, cli_data, cli_manifest = _load_cli_outputs(output)
    frozen = json.loads((output / "frozen_refinement.json").read_text())
    frozen_data = pd.read_csv(FROZEN_GRID, low_memory=False)
    return _summarize_scale(output, resource, cli_data, cli_manifest, frozen, frozen_data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=OUTPUT_ROOT)
    parser.add_argument("--summarize-existing", action="store_true")
    args = parser.parse_args()
    (summarize_scale if args.summarize_existing else run_scale)(args.output_root)
