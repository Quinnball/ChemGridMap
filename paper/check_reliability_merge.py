"""Reproduce T1 acceptance numbers using only the development package.

Run from the repository with its existing Python environment. Frozen inputs
are read only; all generated evidence goes into the requested output folder.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parent
FROZEN = ROOT / "paper/output/revision_v6"
REVIEW = PROJECT / "project1/JCIM_review_innovation_20261024"
CHECKS_ROOT = ROOT / "paper/output/reliability_v1/t1_checks"
FIGURE = PROJECT / "project1/JCIM_molecular_figure_revision_20261014_v2/figures/source_data/Figure_3_core_coverage.csv"
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / "src"))
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "chemgridmap_t1_mpl"))
os.environ.setdefault("NUMBA_NUM_THREADS", "1")

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

import chemgridmap
from chemgridmap.activity import activity_class
from chemgridmap.comparison import structural_groups
from chemgridmap.evidence import (CLASSES, class_probabilities, document_values,
                                 fit_between_unit_sd, predictive_parameters)
from chemgridmap.provenance import environment_versions
from chemgridmap.refinement import REFINEMENT_DEFAULTS, refine_grid_layout
from chemgridmap.reliability import grid_neighbor_sets, neighbor_recall, representation_neighbor_sets
from chemgridmap.representations import morgan_fingerprints

TARGETS = {"chembl205": ("CA II", (7, 21, 42, 84, 126)),
           "chembl204": ("Thrombin", (7, 21, 42)), "chembl240": ("hERG", (7, 21, 42))}
INPUTS = {}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source(path):
    path = Path(path).resolve()
    INPUTS[str(path)] = sha256(path)
    return path


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def read_csv(path):
    return pd.read_csv(source(path), low_memory=False)


def check(name, passed, **details):
    result = {"check": name, "passed": bool(passed), **details}
    print(json.dumps(result, ensure_ascii=False), flush=True)
    return result


def layout(target, seed):
    return read_csv(FROZEN / target / f"coordinates_seed_{seed}.csv").sort_values(
        "molecule_identity_key").reset_index(drop=True)


def package_code():
    diff = subprocess.run(["git", "diff", "release-v0.4.0", "--"], cwd=ROOT,
                          check=True, capture_output=True).stdout
    paths = sorted((ROOT / "src/chemgridmap").rglob("*.py")) + [Path(__file__).resolve()]
    paths += sorted((ROOT / "src/chemgridmap/web").glob("*")) + [ROOT / "pyproject.toml"]
    return {"tracked_diff_sha256": hashlib.sha256(diff).hexdigest(),
            "diff_command": "git diff release-v0.4.0 --",
            "source_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in paths if p.is_file()}}


def manifests(output, context, sections):
    # Only the numeric sections executed now belong to this writer. The UI
    # runner and prior smoke evidence retain their own parameters/provenance.
    folders = [output / section for section in sections if section in {"legacy", "default", "headline"}
               and (output / section).is_dir()]
    if folders and output != CHECKS_ROOT:
        folders.insert(0, output)
    for folder in reversed(folders):
        files = sorted(p for p in folder.rglob("*") if p.is_file() and p.name != "manifest.json")
        dump(folder / "manifest.json", {**context, "inputs_sha256": INPUTS,
             "outputs_sha256": {str(p.relative_to(folder)): sha256(p) for p in files}})


def run_cli(output, mode, extra):
    folder = output / mode
    folder.mkdir(parents=True, exist_ok=True)
    raw = source(ROOT / "validation_data/chembl37_chembl205/chembl37_chembl205_ic50_raw.csv")
    command = [sys.executable, "-m", "chemgridmap.cli", str(raw), "--input-format", "chembl",
               "--target-id", "CHEMBL205", "--projection", "umap", "--random-state", "42",
               "--name", "chembl205", "--output-dir", str(folder), "--output-formats", "svg", *extra]
    env = dict(os.environ, PYTHONPATH=str(ROOT / "src"), PYTHONDONTWRITEBYTECODE="1")
    started = time.perf_counter()
    result = subprocess.run(command, cwd=ROOT, env=env, capture_output=True, text=True)
    (folder / "cli.stdout.txt").write_text(result.stdout, encoding="utf-8")
    (folder / "cli.stderr.txt").write_text(result.stderr, encoding="utf-8")
    dump(folder / "command.json", {"argv": command, "returncode": result.returncode,
                                  "seconds": time.perf_counter() - started})
    if result.returncode:
        raise RuntimeError(f"CLI {mode} failed; see {folder / 'cli.stderr.txt'}")
    return folder, pd.read_csv(folder / "chembl205_grid_coordinates.csv", low_memory=False)


def legacy(output):
    old_manifest = json.loads(source(FROZEN / "run_manifest.json").read_text())
    archived = next(row for row in old_manifest["inputs"] if row["target"] == "CHEMBL205")
    raw = source(ROOT / "validation_data" / archived["file"])
    bundled = source(ROOT / "src/chemgridmap/data/chembl205_ic50.csv")
    source(ROOT / "paper/run_revision_validation.py")
    checks = [check("AC1.3_input_identity", sha256(raw) == archived["sha256"] == sha256(bundled),
                    input=str(raw), archived_sha256=archived["sha256"],
                    current_sha256=sha256(raw), bundled_sha256=sha256(bundled),
                    protocol="curated entries sorted by molecule_identity_key; UMAP/Jaccard seed 42; dense assignment")]
    folder, current = run_cli(output, "legacy", ["--grid-refinement", "none", "--core-rule",
                                                "most_frequent", "--no-activity-confidence"])
    frozen = layout("chembl205", 42)
    fields = ["molecule_identity_key", "canonical_smiles", "activity_pchembl", "grid_col", "grid_row"]
    aligned = frozen[fields].merge(current[fields], on="molecule_identity_key", how="outer",
                                   suffixes=("_frozen", "_current"), validate="one_to_one", indicator=True)
    same_input = (aligned["_merge"].eq("both") &
                  aligned.canonical_smiles_frozen.eq(aligned.canonical_smiles_current) &
                  np.isclose(aligned.activity_pchembl_frozen, aligned.activity_pchembl_current))
    same_cells = aligned.grid_col_frozen.eq(aligned.grid_col_current) & aligned.grid_row_frozen.eq(aligned.grid_row_current)
    aligned.to_csv(folder / "coordinate_comparison.csv", index=False)
    checks.append(check("AC1.3_coordinates", len(aligned) == 889 and same_input.all() and same_cells.all(),
                        entries=len(aligned), input_mismatches=int((~same_input).sum()),
                        cell_mismatches=int((~same_cells).sum())))
    _, groups, _ = structural_groups(current, core_rule="most_frequent")
    actual = pd.DataFrame([{"group": g["structure_group"], "size_current": g["size"],
                            "coverage_count_current": g["core_coverage_count"],
                            "exact_count_current": g["core_exact_scaffold_count"]}
                           for g in groups if g["size"] > 1])
    figure = read_csv(FIGURE).query("target == 'chembl205'")
    comparison = figure.merge(actual, on="group", how="outer", validate="one_to_one", indicator=True)
    matches = comparison["_merge"].eq("both")
    for field in ("size", "coverage_count", "exact_count"):
        matches &= comparison[field].eq(comparison[field + "_current"])
    comparison.to_csv(folder / "core_coverage_comparison.csv", index=False)
    checks.append(check("AC1.3_core_coverage", len(comparison) == 125 and matches.all(),
                        groups=len(comparison), mismatched_groups=int((~matches).sum())))
    return checks


def default_run(output):
    folder, data = run_cli(output, "default", [])
    report = json.loads((folder / "chembl205_curation_report.json").read_text())
    manifest = json.loads((folder / "chembl205_run_manifest.json").read_text())
    _, groups, _ = structural_groups(data, core_rule="coverage", core_min_coverage=0.8)
    group = next(g for g in groups if g["structure_group"] == "G002")
    case = data.loc[data.molecule_identity_key.eq("CHEMBL20")].iloc[0]
    benzothiazole = Chem.MolToSmiles(Chem.MolFromSmiles("c1ccc2scnc2c1"))
    checks = [check("AC1.4_counts", (report["input_rows"], report["retained_activity_records"],
                                    len(data)) == (2281, 1100, 889),
                    raw=report["input_rows"], retained=report["retained_activity_records"], entries=len(data)),
              check("AC1.4_G002", group["size"] == 28 and group["core_coverage_count"] == 27 and
                    Chem.MolToSmiles(Chem.MolFromSmiles(group["core_smiles"])) == benzothiazole,
                    core=group["core_smiles"], covered=group["core_coverage_count"], size=group["size"]),
              check("AC1.4_CHEMBL20", abs(float(case.p_consensus_active) - 0.61) <= 0.01,
                    p_active=float(case.p_consensus_active), evidence_units=int(case.n_evidence_units),
                    document_balanced_pchembl=float(case.document_balanced_pchembl), status=str(case.class_status))]
    expected = {"grid_refinement": "neighbors", "core_rule": "coverage", "core_min_coverage": 0.8,
                "activity_confidence": True, "merge_identical_documents": True, "map_color_modes": []}
    parameters = manifest["parameters"]
    differences = {key: {"expected": value, "actual": parameters.get(key)}
                   for key, value in expected.items() if parameters.get(key) != value}
    fit = report.get("activity_confidence", {})
    checks.append(check("AC1.4_manifest", not differences and "evidence_unit_policy" in fit,
                        expected_new_parameters=expected, parameter_differences=differences,
                        evidence_unit_policy=fit.get("evidence_unit_policy")))
    dump(folder / "case_summary.json", {"G002": group,
         "CHEMBL20": {k: case[k].item() if isinstance(case[k], np.generic) else case[k]
                      for k in ("molecule_identity_key", "n_records", "n_evidence_units", "p_consensus_active",
                                "document_balanced_pchembl", "class_status")}})
    return checks


def lodo_gaussian(target, records, folds=10, seed=20261024):
    # Match the review protocol: fit on other entries, then hold out each unit.
    doc = document_values(records, merge_identical=True)
    counts = doc.groupby("entry").size()
    multi = counts[counts >= 2].index.to_numpy()
    fold_of = dict(zip(np.random.default_rng(seed).permutation(multi), np.arange(len(multi)) % folds))
    rows = []
    for fold in range(folds):
        test = {entry for entry, assigned in fold_of.items() if assigned == fold}
        fit = fit_between_unit_sd(doc[~doc.entry.isin(test)])
        for entry in sorted(test):
            entry_docs = doc[doc.entry == entry]
            for _, held in entry_docs.iterrows():
                rest = entry_docs[entry_docs.document != held.document].document_value.to_numpy(float)
                own = max(np.var(rest, ddof=1), 1e-4) if len(rest) > 1 else 0.0
                loc, scale, df = predictive_parameters(np.array([rest.mean()]), np.array([own]),
                    np.array([float(len(rest))]), fit["tau"], fit["prior_df"])
                p = class_probabilities(loc, scale, df, 6.0, 7.0)[0]
                rows.append({"target": target, "fold": fold, "entry": entry,
                             "held_out_unit": held.documents, "remaining_units": len(rest),
                             "observed": CLASSES.index(activity_class(held.document_value, 6, 7)),
                             **{f"p{k}": float(p[k]) for k in range(3)}})
    return pd.DataFrame(rows)


def probability_scores(frame):
    p = frame[["p0", "p1", "p2"]].to_numpy(float)
    observed = frame.observed.to_numpy(int)
    confidence = p.max(axis=1)
    correct = (p.argmax(axis=1) == observed).astype(float)
    bins = np.minimum((confidence * 10).astype(int), 9)
    ece = sum(abs(correct[bins == b].mean() - confidence[bins == b].mean()) * (bins == b).mean()
              for b in range(10) if (bins == b).any())
    return {"predictions": len(frame), "entries": int(frame.entry.nunique()), "ece": float(ece),
            "accuracy": float(correct.mean()), "brier": float(((p - np.eye(3)[observed]) ** 2).sum(1).mean())}


def headline(output):
    folder = output / "headline"
    folder.mkdir(parents=True, exist_ok=True)
    reference = json.loads(source(REVIEW / "results/headline_numbers.json").read_text())
    # The read-only historical reference retains the legacy target key.
    reference["grid_refinement"] = {"Thrombin" if label == "D3" else label: section
                                    for label, section in reference["grid_refinement"].items()}
    for name in ("common.py", "evaluate_grid_refinement.py", "evaluate_core_selection.py",
                 "evaluate_activity_confidence.py"):
        source(REVIEW / "analysis" / name)
    rows, core_rows, validation, diagnostics, checks = [], [], [], {}, []
    for target, (label, seeds) in TARGETS.items():
        data = layout(target, 42)
        fingerprint = morgan_fingerprints(data.canonical_smiles.tolist())
        neighbors = representation_neighbor_sets(fingerprint, k=10, metric="jaccard")
        for seed in seeds:
            frozen = layout(target, seed)
            assert data.molecule_identity_key.equals(frozen.molecule_identity_key)
            side = int(np.rint(1 / np.diff(np.unique(frozen.grid_x))[0])) + 1
            cells = frozen[["grid_col", "grid_row"]].to_numpy(int)
            projection = frozen[["projection_x", "projection_y"]].to_numpy(float) * (side - 1)
            diag = {}
            started = time.perf_counter()
            refined = refine_grid_layout(cells, projection, neighbors, side, random_state=0, diagnostics=diag)
            row = {"target": label, "seed": seed, "refinement_seed": 0, "entries": len(data),
                   "recall_before": float(neighbor_recall(neighbors, grid_neighbor_sets(cells)).mean()),
                   "recall_after": float(neighbor_recall(neighbors, grid_neighbor_sets(refined)).mean()),
                   "seconds": time.perf_counter() - started}
            rows.append(row)
            diagnostics[f"{target}_seed_{seed}"] = diag
            saved = frozen[["molecule_identity_key", "grid_col", "grid_row"]].copy()
            saved[["refined_col", "refined_row"]] = refined
            saved.to_csv(folder / f"{target}_seed_{seed}_refined.csv", index=False)
            print(json.dumps(row), flush=True)
        old_groups = {g["structure_group"]: g for g in structural_groups(data, core_rule="most_frequent")[1]}
        groups = structural_groups(data, core_rule="coverage", core_min_coverage=0.8)[1]
        for group in groups:
            if group["size"] > 1:
                old = old_groups[group["structure_group"]]
                core_rows.append({"target": label, "group": group["structure_group"], "size": group["size"],
                                  "old_coverage_count": old["core_coverage_count"],
                                  "coverage_count": group["core_coverage_count"], "core": group["core_smiles"]})
        records = read_csv(FROZEN / target / f"{target}_retained_records.csv")
        validation.append(lodo_gaussian(label, records))
    metrics = pd.DataFrame(rows)
    metrics.to_csv(folder / "grid_refinement.csv", index=False)
    dump(folder / "refinement_diagnostics.json", diagnostics)
    for label, group in metrics.groupby("target", sort=False):
        actual = float(group.recall_after.mean())
        expected = reference["grid_refinement"][label]["recall_k10_neighbor_refined"]
        checks.append(check(f"headline_recall_{label}", abs(actual - expected) < 1e-12,
                            expected=expected, actual=actual, absolute_difference=abs(actual - expected),
                            frozen_layouts=len(group), refinement_seed=0))
    cores = pd.DataFrame(core_rows)
    cores.to_csv(folder / "core_coverage.csv", index=False)
    counts = {label: {"groups": len(group), "majority": int((group.coverage_count > group["size"] / 2).sum())}
              for label, group in cores.groupby("target", sort=False)}
    majority = int((cores.coverage_count > cores["size"] / 2).sum())
    checks.append(check("headline_core_coverage", len(cores) == 480 and majority == 470,
                        groups=len(cores), majority=majority, by_target=counts))
    predictions = pd.concat(validation, ignore_index=True)
    predictions.to_csv(folder / "lodo_gaussian_predictions.csv", index=False)
    scores = [{"target": name, **probability_scores(frame)}
              for name, frame in [("all", predictions), *predictions.groupby("target")]]
    pd.DataFrame(scores).to_csv(folder / "lodo_gaussian_scores.csv", index=False)
    expected = next(r for r in reference["activity_lodo_scores"] if r["target"] == "all" and r["method"] == "G")
    score = scores[0]
    checks.append(check("headline_gaussian_ECE", round(score["ece"], 4) == expected["ece"] and
                        score["predictions"] == 820 and score["entries"] == 289,
                        expected_rounded_ece=expected["ece"], **score))
    return checks


def application_check(output, mode):
    """Run AC1.6 and keep an environment block separate from a failed check."""
    folder = CHECKS_ROOT / mode
    folder.mkdir(parents=True, exist_ok=True)
    source(ROOT / "src/chemgridmap/app.py")
    source(ROOT / "src/chemgridmap/data/chembl205_ic50.csv")
    env = dict(os.environ, PYTHONPATH=str(ROOT / "src"), PYTHONDONTWRITEBYTECODE="1")
    if mode == "ui":
        if "CHEMGRIDMAP_NODE_PATH" in env:
            env["NODE_PATH"] = env["CHEMGRIDMAP_NODE_PATH"]
        command = ["node", str(source(ROOT / "paper/check_reliability_ui.cjs")), str(folder)]
    else:
        command = [sys.executable, "-m", "chemgridmap.app", "--smoke-test", "--output-dir", str(folder)]
    started = time.perf_counter()
    result = subprocess.run(command, cwd=ROOT, env=env, capture_output=True, text=True)
    if mode == "ui" and result.returncode and re.search(
            r"Cannot find module ['\"]playwright['\"]", result.stderr):
        hint = "Set CHEMGRIDMAP_NODE_PATH to the node_modules directory containing Playwright."
        result.stderr += "\n" + hint + "\n"
        print(hint, file=sys.stderr, flush=True)
    (output / f"{mode}_runner.stdout.txt").write_text(result.stdout, encoding="utf-8")
    (output / f"{mode}_runner.stderr.txt").write_text(result.stderr, encoding="utf-8")
    status = "passed" if result.returncode == 0 else "failed"
    if result.returncode and re.search(r"Operation not permitted|PermissionError|EPERM|sandbox.*denied",
                                       result.stdout + result.stderr, re.IGNORECASE):
        status = "sandbox-blocked"
    if mode == "ui" and result.returncode == 2:
        status = "sandbox-blocked"
    row = {"check": f"AC1.6_{mode}", "status": status,
           "passed": None if status == "sandbox-blocked" else status == "passed",
           "returncode": result.returncode, "argv": command, "output": str(folder),
           "seconds": time.perf_counter() - started}
    dump(output / f"{mode}_runner_command.json", row)
    print(json.dumps(row, ensure_ascii=False), flush=True)
    return [row]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=CHECKS_ROOT)
    parser.add_argument("--sections", nargs="+", choices=["legacy", "default", "headline", "smoke", "ui"],
                        default=["legacy", "default", "headline", "smoke", "ui"])
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    RDLogger.DisableLog("rdApp.*")
    assert Path(chemgridmap.__file__).resolve().is_relative_to(ROOT / "src"), "Wrong package imported"
    context = {"started_at": datetime.now(timezone.utc).isoformat(), "command": sys.argv,
               "parameters": {"sections": args.sections, "umap_seed": 42, "refinement_seed": 0,
                              "refinement": REFINEMENT_DEFAULTS, "core_min_coverage": 0.8,
                              "lodo_folds": 10, "lodo_seed": 20261024, "merge_identical": True},
               "environment": environment_versions(), "imported_package": chemgridmap.__file__,
               "package_version": chemgridmap.__version__, "code": package_code()}
    checks = []
    executed = []
    try:
        for section in args.sections:
            executed.append(section)
            if section in {"smoke", "ui"}:
                checks.extend(application_check(output, section))
            else:
                checks.extend({"legacy": legacy, "default": default_run, "headline": headline}[section](output))
    except Exception as exc:
        checks.append(check("execution", False, error=f"{type(exc).__name__}: {exc}"))
        raise
    finally:
        has_application_check = any(section in {"smoke", "ui"} for section in executed)
        dump(output / ("runner_checks.json" if has_application_check else "checks.json"), checks)
        context["finished_at"] = datetime.now(timezone.utc).isoformat()
        if has_application_check:
            context["inputs_sha256"] = INPUTS
            dump(output / "runner_context.json", context)
        manifests(output, context, executed)
    if any(item["passed"] is False for item in checks):
        return 1
    return 2 if any(item.get("status") == "sandbox-blocked" for item in checks) else 0


if __name__ == "__main__":
    raise SystemExit(main())
