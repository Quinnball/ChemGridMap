"""Summarize neighbour colour contrasts and verify the Figure 1 CLI exports."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import time

os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "chemgridmap-t2b-mpl"))

import numpy as np
import pandas as pd

from chemgridmap.evidence import CLASSES
from chemgridmap.provenance import environment_versions
from chemgridmap.reliability import grid_neighbor_sets
from reliability_analysis import OUTPUT, ROOT, SEEDS, TARGETS, dump_json, layout_points, sha256

RAW_INPUT = ROOT / "validation_data/chembl37_chembl205/chembl37_chembl205_ic50_raw.csv"
CLI_NAME = "caii_seed42"
LAYOUTS = ("minimum_displacement", "neighbor_refined")
DEFINITIONS = ("eight_connected", "tie_aware_k10_union")
PROBABILITY_COLUMNS = [f"p_consensus_{name}" for name in CLASSES]


def eight_connected_pairs(cells):
    """Return unordered occupied-cell pairs at Chebyshev distance one."""
    cells = np.asarray(cells)
    if cells.ndim != 2 or cells.shape[1] != 2 or not np.equal(cells, np.rint(cells)).all():
        raise ValueError("Cells must be an integer N by 2 array.")
    occupied = {tuple(cell): i for i, cell in enumerate(cells)}
    if len(occupied) != len(cells):
        raise ValueError("Cells must have unique occupancy.")
    pairs = []
    for i, (x, y) in enumerate(cells):
        for dx, dy in ((1, -1), (1, 0), (1, 1), (0, 1)):
            j = occupied.get((x + dx, y + dy))
            if j is not None:
                pairs.append((min(i, j), max(i, j)))
    return np.asarray(sorted(pairs), dtype=int).reshape(-1, 2)


def tie_aware_neighbor_pairs(cells, k=10):
    """Take the symmetric union of tie-inclusive directed grid neighbours."""
    neighbors = grid_neighbor_sets(cells, k=k)
    pairs = {(min(i, int(j)), max(i, int(j)))
             for i, row in enumerate(neighbors) for j in row if i != j}
    return np.asarray(sorted(pairs), dtype=int).reshape(-1, 2)


def contrast_statistics(pairs, classes, probabilities, threshold=.8):
    """Sum independent consensus class disagreement probabilities over pairs."""
    pairs = np.asarray(pairs, dtype=int).reshape(-1, 2)
    probabilities = np.asarray(probabilities, dtype=float)
    if (not np.isfinite(probabilities).all() or (probabilities < 0).any()
            or not np.allclose(probabilities.sum(axis=1), 1, atol=1e-12, rtol=0)):
        raise ValueError("Consensus class probabilities must be finite and sum to one.")
    i, j = pairs.T
    different = np.asarray(classes)[i] != np.asarray(classes)[j]
    probability_different = 1 - (probabilities[i] * probabilities[j]).sum(axis=1)
    return {"neighbor_pairs": len(pairs),
            "displayed_contrasts": int(different.sum()),
            "expected_true_contrasts_among_displayed": float(probability_different[different].sum()),
            "supported_contrasts_p_ge_0.8": int((probability_different[different] >= threshold).sum()),
            "same_colour_pairs": int((~different).sum()),
            "expected_hidden_contrasts_same_colour": float(probability_different[~different].sum())}


def run_contrasts(output_root=OUTPUT, seeds=SEEDS, targets=TARGETS):
    """Read final annotations and frozen T2a layouts without changing either."""
    output_root = Path(output_root)
    adoption = json.loads((output_root / "activity/adoption.json").read_text())
    rows, inputs = [], [output_root / "activity/adoption.json"]
    for target in targets:
        annotation_path = output_root / "activity" / f"{target}_annotations.csv"
        annotations = pd.read_csv(annotation_path).set_index("molecule_identity_key", verify_integrity=True)
        inputs.append(annotation_path)
        for seed in seeds:
            layout_path = output_root / "layouts" / target / f"coordinates_seed_{seed}.csv"
            layout = pd.read_csv(layout_path)
            if set(layout.molecule_identity_key) != set(annotations.index):
                raise ValueError(f"Layout and annotation identities differ: {target} seed {seed}.")
            aligned = annotations.loc[layout.molecule_identity_key]
            inputs.append(layout_path)
            for name in LAYOUTS:
                cells = layout_points(layout, name)
                for definition, pairs in (
                        (DEFINITIONS[0], eight_connected_pairs(cells)),
                        (DEFINITIONS[1], tie_aware_neighbor_pairs(cells))):
                    rows.append({"target": target, "target_label": TARGETS[target],
                                 "seed": seed, "layout": name, "definition": definition,
                                 "n_entries": len(layout),
                                 **contrast_statistics(pairs, aligned.activity_class,
                                                       aligned[PROBABILITY_COLUMNS])})
    table = pd.DataFrame(rows)
    path = output_root / "activity/neighbor_contrasts.csv"
    table.to_csv(path, index=False)
    return {"inputs": inputs, "output": path, "rows": len(table),
            "parameters": {"adoption": adoption, "seeds": list(seeds),
                           "layouts": list(LAYOUTS), "definitions": list(DEFINITIONS),
                           "k": 10, "threshold": .8,
                           "probability_assumption": "independent consensus classes between entries"}}


def figure1_cli_command(directory):
    """Use the paper's UMAP projection and the CLI's default refinement."""
    return [str(ROOT / ".venv/bin/python"), "-B", "-m", "chemgridmap", str(RAW_INPUT),
            "--input-format", "chembl", "--target-id", "CHEMBL205", "--projection", "umap",
            "--random-state", "42", "--map-color-modes", "activity", "class_support", "neighbor_recall",
            "--output-dir", str(directory), "--name", CLI_NAME]


def _cli_signature(command, environment):
    source = {str(path.relative_to(ROOT)): sha256(path)
              for path in sorted((ROOT / "src/chemgridmap").rglob("*"))
              if path.is_file() and path.suffix in {".py", ".js", ".html", ".css"}}
    payload = {"argv": command, "raw_input_sha256": sha256(RAW_INPUT), "source_sha256": source,
               "versions": environment_versions(),
               "thread_environment": {key: environment[key] for key in (
                   "NUMBA_NUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")}}
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest(), payload


def valid_cli_cache(directory, signature):
    """Reuse exports only when the invocation and all saved file hashes match."""
    directory = Path(directory)
    path = directory / "cli_cache.json"
    if not path.exists():
        return False
    cache = json.loads(path.read_text())
    if cache.get("signature") != signature or cache.get("returncode") != 0 or not cache.get("outputs_sha256"):
        return False
    return all((directory / name).is_file() and sha256(directory / name) == digest
               for name, digest in cache["outputs_sha256"].items())


def prepare_figure1_cli(output_root=OUTPUT):
    """Run the real CLI once; signature-validated exports can be reused."""
    directory = Path(output_root) / "figure1_cli_run"
    directory.mkdir(parents=True, exist_ok=True)
    command = figure1_cli_command(directory)
    environment = dict(os.environ)
    for key in ("NUMBA_NUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        environment[key] = "1"
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    environment.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "chemgridmap-t2b-mpl"))
    signature, payload = _cli_signature(command, environment)
    if valid_cli_cache(directory, signature):
        print("Figure 1: verified CLI cache reused", flush=True)
        return directory
    invocation = {**payload, "signature": signature, "command": shlex.join(command), "cwd": str(ROOT),
                  "started_utc": datetime.now(timezone.utc).isoformat(),
                  "projection_note": "UMAP is explicit to reproduce T2a; the bare CLI projection default is PCA.",
                  "grid_refinement": "neighbors (CLI default)"}
    dump_json(directory / "cli_invocation.json", invocation)
    print("Figure 1: fresh CA II CLI run started", flush=True)
    start = time.perf_counter()
    with (directory / "cli_stdout.log").open("w") as stdout, (directory / "cli_stderr.log").open("w") as stderr:
        result = subprocess.run(command, cwd=ROOT, env=environment, stdout=stdout, stderr=stderr)
    invocation.update(returncode=result.returncode, wall_seconds=time.perf_counter() - start,
                      finished_utc=datetime.now(timezone.utc).isoformat())
    dump_json(directory / "cli_invocation.json", invocation)
    if result.returncode:
        raise RuntimeError(f"Figure 1 CLI exited {result.returncode}; see {directory / 'cli_stderr.log'}.")
    outputs = {path.name: sha256(path) for path in sorted(directory.glob(CLI_NAME + "*")) if path.is_file()}
    dump_json(directory / "cli_cache.json", {"signature": signature, "returncode": result.returncode,
                                           "outputs_sha256": outputs})
    return directory


def compare_figure1_entries(cli, layout, annotations, tolerance=1e-12):
    """Compare by identity rather than relying on export row order."""
    key = "molecule_identity_key"
    cli = cli.set_index(key, verify_integrity=True)
    layout = layout.set_index(key, verify_integrity=True)
    annotations = annotations.set_index(key, verify_integrity=True)
    if set(cli.index) != set(layout.index) or set(cli.index) != set(annotations.index):
        raise ValueError("CLI, T2a layout and annotation identity sets must agree.")
    layout = layout.loc[cli.index]
    annotations = annotations.loc[cli.index]
    comparison = pd.DataFrame(index=cli.index)
    for axis in ("col", "row"):
        comparison[f"cli_grid_{axis}"] = cli[f"grid_{axis}"]
        comparison[f"t2a_refined_grid_{axis}"] = layout[f"refined_grid_{axis}"]
    comparison["coordinates_match"] = np.all(
        cli[["grid_col", "grid_row"]].to_numpy() ==
        layout[["refined_grid_col", "refined_grid_row"]].to_numpy(), axis=1)
    for column in ["class_confidence", *PROBABILITY_COLUMNS]:
        comparison[f"cli_{column}"] = cli[column]
        comparison[f"annotation_{column}"] = annotations[column]
        comparison[f"{column}_absolute_difference"] = np.abs(cli[column] - annotations[column])
    numeric_columns = [f"{column}_absolute_difference" for column in ["class_confidence", *PROBABILITY_COLUMNS]]
    comparison["numeric_support_match"] = (comparison[numeric_columns] <= tolerance).all(axis=1)
    for column in ("activity_class", "class_status"):
        comparison[f"cli_{column}"] = cli[column]
        comparison[f"annotation_{column}"] = annotations[column]
        comparison[f"{column}_match"] = cli[column] == annotations[column]
    comparison["class_support_match"] = comparison[
        ["numeric_support_match", "activity_class_match", "class_status_match"]].all(axis=1)
    checks = {"n_entries": len(comparison), "expected_n_entries": 889,
              "grid_mismatched_entries": int((~comparison.coordinates_match).sum()),
              "class_support_mismatched_entries": int((~comparison.class_support_match).sum()),
              "class_confidence_max_absolute_difference": float(comparison.class_confidence_absolute_difference.max()),
              "consensus_probability_max_absolute_difference": float(comparison[
                  [f"{column}_absolute_difference" for column in PROBABILITY_COLUMNS]].to_numpy().max()),
              "support_absolute_tolerance": tolerance,
              "passed": bool(comparison.coordinates_match.all() and comparison.class_support_match.all())}
    return comparison.reset_index(), checks


def run_figure1_cli(output_root=OUTPUT):
    """Verify all 889 CA II CLI entries against the final analysis outputs."""
    output_root = Path(output_root)
    directory = prepare_figure1_cli(output_root)
    cli_path = directory / f"{CLI_NAME}_grid_coordinates.csv"
    layout_path = output_root / "layouts/chembl205/coordinates_seed_42.csv"
    annotation_path = output_root / "activity/chembl205_annotations.csv"
    comparison, checks = compare_figure1_entries(
        pd.read_csv(cli_path), pd.read_csv(layout_path), pd.read_csv(annotation_path))
    checks["passed"] = checks["passed"] and checks["n_entries"] == checks["expected_n_entries"]
    checks["inputs_sha256"] = {str(path): sha256(path) for path in (cli_path, layout_path, annotation_path)}
    checks["projection"] = "UMAP/Jaccard (explicit; bare CLI default is PCA)"
    comparison.to_csv(directory / "entry_comparison.csv", index=False)
    dump_json(directory / "consistency_check.json", checks)
    if not checks["passed"]:
        raise RuntimeError(f"Figure 1 CLI consistency failed: {checks}.")
    return checks


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps", nargs="+", choices=("contrasts", "prepare-cli", "figure1-cli"),
                        default=["contrasts", "figure1-cli"])
    args = parser.parse_args()
    for step in args.steps:
        result = {"contrasts": run_contrasts, "prepare-cli": prepare_figure1_cli,
                  "figure1-cli": run_figure1_cli}[step]()
        print(json.dumps(result, default=str), flush=True)
