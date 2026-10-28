"""Check T2b output hashes, required populations, and the pre-run file boundary."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from reliability_analysis import ROOT, OUTPUT, sha256, dump_json
from reliability_t2b import RESULT_MARKER, README_MARKER, REFERENCE, flatten


def latest_provenance(manifest):
    """Use the newest source record while retaining earlier run provenance."""
    return next(manifest[name] for name in
                ("software_checks_050", "portability_fix", "label_correction", "t2b_rework", "t2b")
                if name in manifest)


def audit_outputs(root=OUTPUT):
    root = Path(root)
    manifest = json.loads((root / "manifest.json").read_text())
    actual_paths = {str(p.relative_to(root)) for p in root.rglob("*") if p.is_file()
                    and p != root / "manifest.json" and p.name != ".DS_Store" and "__pycache__" not in p.parts}
    recorded_paths = set(manifest["outputs_sha256"])
    errors = [f"Root manifest missing outputs: {sorted(actual_paths - recorded_paths)}"] if actual_paths != recorded_paths else []
    verified = 0
    for relative, digest in manifest["outputs_sha256"].items():
        if not (root / relative).is_file() or sha256(root / relative) != digest:
            errors.append(f"Output hash mismatch: {relative}")
        verified += 1
    # Earlier provenance describes the accepted run, not the current source tree.
    current = latest_provenance(manifest)
    for section in ("inputs_sha256",):
        for name, digest in current[section].items():
            if sha256(name) != digest:
                errors.append(f"Current input hash mismatch: {name}")
    for name, digest in current["code"]["source_sha256"].items():
        if sha256(name) != digest:
            errors.append(f"Current source hash mismatch: {name}")
    for name, digest in current["reports_sha256"].items():
        if sha256(name) != digest:
            errors.append(f"Report hash mismatch: {name}")
    counts = {}
    for relative, expected in {"activity/mapped_entry_confidence.csv": 3289,
        "activity/mapped_entry_confidence_unmerged.csv": 3289, "activity/neighbor_contrasts.csv": 120,
        "cores/group_core_comparison.csv": 1440, "cores/core_target_sensitivity.csv": 2400,
        "cores/fmcs_recomputation.csv": 20, "cases/core_rule_examples.csv": 5,
        "cases/CHEMBL20_evidence_units.csv": 27, "figure_data/Figure_1ABC_map_seed42_refined.csv": 889,
        "figure_data/Figure_1D_G002_core.csv": 28, "figure_data/Figure_2D_supported_classes.csv": 9,
        "figure_data/Figure_2E_neighbor_contrasts.csv": 18, "figure_data/Figure_3D_core_coverage.csv": 960}.items():
        counts[relative] = len(pd.read_csv(root / relative))
        if counts[relative] != expected:
            errors.append(f"Wrong row count: {relative}: {counts[relative]} != {expected}")
    comparison = pd.read_csv(root / "activity/headline_comparison.csv")
    old = json.loads(REFERENCE.read_text())
    expected_headlines = {"/".join(path) for path, _ in flatten(old)}
    if set(comparison.headline_path) != expected_headlines or not comparison.headline_path.is_unique:
        errors.append("Historical headline comparison is incomplete or duplicated.")
    historical_scores = comparison[comparison.headline_path.str.startswith("activity_lodo_scores/")]
    differences = pd.to_numeric(historical_scores.difference, errors="coerce").dropna()
    if not differences.abs().le(5.00001e-5).all():
        errors.append("Historical LODO scores differ beyond reference rounding precision.")
    scores = pd.read_csv(root / "activity/lodo_scores.csv").set_index(["target", "method"])
    g = scores.loc[("all", "G")]
    if not (g.n_predictions == 820 and g.n_entries == 289 and np.isclose(g.ece, .0767537203890331, atol=1e-14)):
        errors.append("Gaussian baseline does not reproduce.")
    for relative in ("cores/core_checks.json", "cases/CHEMBL20_activity_summary.json",
                     "figure_data/t2b_checks.json", "figure1_cli_run/consistency_check.json",
                     "activity/protocol_checks.json"):
        if not json.loads((root / relative).read_text())["passed"]:
            errors.append(f"Check failed: {relative}")
    rework = audit_core_timing(root) if "t2b_rework" in manifest else None
    if rework:
        errors.extend(rework["errors"])
        counts.update(rework["required_row_counts"])
    return {"passed": not errors, "errors": errors, "verified_root_output_hashes": verified,
            "required_row_counts": counts, "historical_headline_leaves": len(expected_headlines),
            "core_timing_rework": rework,
            "t1_t2a_provenance_sections_retained": "t2a_rework" in manifest and "baseline_commit" in manifest}


def audit_core_timing(root):
    """Reconcile the new timing tables independently against group rows and the saved cache."""
    from reliability_cores import CACHE, EXPECTED_GROUPS, TARGETS

    root = Path(root)
    errors, counts = [], {}

    def require(condition, message):
        if not condition:
            errors.append(message)

    files = ("core_rule_timing.csv", "fmcs_full_recomputation.csv", "core_rule_timing_summary.csv",
             "core_timing_metadata.json", "core_timing_checks.json")
    for name in files:
        require((root / "cores" / name).is_file(), f"Missing core timing output: {name}")
    if errors:
        return {"passed": False, "errors": errors, "required_row_counts": counts}
    timing = pd.read_csv(root / "cores/core_rule_timing.csv")
    audit = pd.read_csv(root / "cores/fmcs_full_recomputation.csv", keep_default_na=False)
    summary = pd.read_csv(root / "cores/core_rule_timing_summary.csv")
    checks = json.loads((root / "cores/core_timing_checks.json").read_text())
    metadata = json.loads((root / "cores/core_timing_metadata.json").read_text())
    key = ["target", "group"]
    frozen = pd.read_csv(root / "cores/group_core_comparison.csv", keep_default_na=False)
    frozen = frozen.loc[frozen.rule.eq("fmcs")].set_index(key)
    keys = set(frozen.index)
    for name, table in (("core_rule_timing.csv", timing), ("fmcs_full_recomputation.csv", audit)):
        counts["cores/" + name] = len(table)
        require(len(table) == 480, f"Wrong row count: cores/{name}: {len(table)} != 480")
        require(not table.duplicated(key).any(), f"Duplicate target/group key: {name}")
        require(set(map(tuple, table[key].to_numpy())) == keys, f"Wrong group population: {name}")
        require(table.groupby("target").size().to_dict() == EXPECTED_GROUPS, f"Wrong target counts: {name}")
        require(table.target.map(TARGETS).eq(table.target_label).all(), f"Wrong target labels: {name}")
    if errors:
        return {"passed": False, "errors": errors, "required_row_counts": counts}
    timing, audit = timing.set_index(key).sort_index(), audit.set_index(key).sort_index()
    require(timing["size"].ge(2).all() and timing["size"].equals(audit["size"]), "Invalid nonsingleton group sizes.")
    require(timing["size"].equals(frozen["size"].sort_index()), "Group sizes differ from frozen comparison.")
    seconds = timing.filter(regex="_seconds$").to_numpy(dtype=float)
    require(np.isfinite(seconds).all() and (seconds >= 0).all(), "Core timing contains invalid seconds.")
    for name in ("same_cached_size", "same_cached_most_frequent_core", "same_cached_coverage_core"):
        require(timing[name].eq(True).all(), f"Core cache check failed: {name}")
    for left, right in (("fmcs_seconds", "recomputed_seconds"), ("history_mcs_seconds", "cached_seconds"),
                        ("fmcs_timed_out", "recomputed_timed_out"), ("history_mcs_timed_out", "cached_timed_out"),
                        ("fmcs_coverage_count", "recomputed_coverage_count")):
        require(np.allclose(timing[left], audit[right], atol=1e-12, rtol=1e-12), f"Timing/audit mismatch: {left}")
    cache = pd.read_csv(CACHE, keep_default_na=False)
    # The frozen cache retains its historical target label; normalize only in memory.
    cache["target"] = cache["target"].replace({"D3": "Thrombin"})
    cache = cache.set_index(["target", "group"])
    selected_cache = cache.loc[list(zip(timing.target_label, timing.index.get_level_values("group")))].copy()
    selected_cache.index = timing.index
    for left, right in (("history_aware_seconds", "aware_seconds"), ("history_mcs_seconds", "mcs_seconds"),
                        ("history_mcs_timed_out", "mcs_timed_out")):
        require(np.allclose(timing[left], selected_cache[right], atol=1e-12, rtol=1e-12), f"Historical cache mismatch: {left}")
    require(audit.cached_smarts.eq(selected_cache.mcs_core).all(), "Cached SMARTS differ from the historical cache.")
    require(audit.cached_coverage_count.eq(frozen.coverage_count.sort_index()).all(), "Cached coverage counts differ from frozen groups.")
    same = {}
    for label, field in (("smarts", "smarts"), ("heavy_atoms", "heavy_atoms"), ("ring_count", "ring_count"),
                         ("member_matches", "matched_member_ids"), ("coverage_count", "coverage_count"),
                         ("timeout_status", "timed_out")):
        same["same_" + label] = audit["cached_" + field].eq(audit["recomputed_" + field])
    same["same_coverage"] = pd.Series(np.isclose(audit.cached_coverage_fraction,
        audit.recomputed_coverage_fraction, atol=1e-14, rtol=0), index=audit.index)
    same = pd.DataFrame(same)
    for name in same:
        require(same[name].eq(audit[name]).all(), f"FMCS comparison flag mismatch: {name}")
    passed = same.all(axis=1)
    require(passed.eq(audit.passed).all(), "FMCS comparison passed flag mismatch.")
    timeout = audit.cached_timed_out | audit.recomputed_timed_out
    fmcs_counts = {"groups": len(audit), "exact_matches": int(passed.sum()), "different_groups": int((~passed).sum()),
                   "timeout_either_run_groups": int(timeout.sum()),
                   "differences_with_timeout": int((~passed & timeout).sum()),
                   "differences_without_timeout": int((~passed & ~timeout).sum()),
                   "cached_timeout_groups": int(audit.cached_timed_out.sum()),
                   "recomputed_timeout_groups": int(audit.recomputed_timed_out.sum()),
                   "difference_counts": {name: int((~same[name]).sum()) for name in same},
                   "non_timeout_reproduction_passed": bool(passed.loc[~timeout].all())}
    require(checks["fmcs_recomputation"] == fmcs_counts, "FMCS difference counts do not reconcile with audit rows.")
    require(fmcs_counts["non_timeout_reproduction_passed"], "FMCS differs for completed, non-timeout searches.")
    require(checks["passed"] and checks["groups_total"] == 480 and checks["nonsingleton_groups"] == EXPECTED_GROUPS,
            "Core timing checks failed or report incorrect populations.")
    require(metadata.get("workers") == 1 and metadata.get("platform") and metadata.get("cpu_model"),
            "Core timing machine metadata or single-process setting is missing.")
    rules = {"most_frequent": ("most_frequent_seconds", None), "coverage": ("coverage_seconds", None),
             "fmcs": ("fmcs_seconds", "fmcs_timed_out"), "historical_coverage": ("history_aware_seconds", None),
             "historical_fmcs": ("history_mcs_seconds", "history_mcs_timed_out")}
    expected = []
    for target, label in [*TARGETS.items(), ("overall", "Overall")]:
        group = timing if target == "overall" else timing.loc[target]
        for rule, (column, timeout_column) in rules.items():
            values = group[column]
            expected.append({"target": target, "target_label": label, "rule": rule, "groups": len(group),
                             "total_seconds": values.sum(), "median_seconds": values.median(),
                             "p90_seconds": values.quantile(.9), "max_seconds": values.max(),
                             "timeout_groups": group[timeout_column].sum() if timeout_column else 0})
    summary_key = ["target", "target_label", "rule"]
    expected = pd.DataFrame(expected).set_index(summary_key).sort_index()
    require(not summary.duplicated(summary_key).any(), "Duplicate timing summary keys.")
    summary = summary.set_index(summary_key).sort_index()
    require(summary.index.equals(expected.index), "Timing summary populations are incomplete.")
    if summary.index.equals(expected.index):
        require(np.allclose(summary[expected.columns], expected, atol=1e-10, rtol=1e-10),
                "Timing summary statistics do not reconcile with group rows.")
    si = pd.read_csv(root / "figure_data/SI_core_size_tradeoff.csv").set_index(summary_key)
    expected_si = expected.loc[expected.index.get_level_values("rule").isin(["most_frequent", "coverage", "fmcs"])]
    require(si.index.is_unique and set(si.index) == set(expected_si.index), "SI timing populations are incomplete or duplicated.")
    timing_fields = ["total_seconds", "median_seconds", "p90_seconds", "max_seconds", "timeout_groups"]
    for index in expected_si.index:
        if index not in si.index or not si.index.is_unique:
            continue
        for prefix, rule in (("timing_", index[2]), ("historical_", "historical_" + index[2])):
            fields = [prefix + name for name in timing_fields]
            if not set(fields) <= set(si.columns):
                require(False, f"SI timing columns are missing: {prefix}")
                break
            actual = si.loc[index, fields].to_numpy(dtype=float)
            if rule == "historical_most_frequent":
                require(np.isnan(actual).all(), "SI invents a historical most-frequent timing.")
            else:
                values = expected.loc[(index[0], index[1], rule), timing_fields].to_numpy(dtype=float)
                require(np.allclose(actual, values, atol=1e-10, rtol=1e-10), f"SI timing statistics mismatch: {index}, {prefix}")
    return {"passed": not errors, "errors": errors, "required_row_counts": counts,
            "fmcs_recomputation": fmcs_counts, "summary_rows": len(summary), "si_rows": len(si)}


def audit_scope(baseline, root=ROOT.parent, *, rework=False, reference_dir=None, external_changes=()):
    if rework:
        return audit_rework_scope(baseline, root, reference_dir=reference_dir, external_changes=external_changes)
    root = Path(root)
    before = json.loads(Path(baseline).read_text())
    changes, missing = [], []
    for relative, digest in before.items():
        path = root / relative
        if not path.exists():
            missing.append(relative)
        elif sha256(path) != digest:
            changes.append(relative)
    allowed = {"chemgridmap/paper/run_reliability_analysis.py", "chemgridmap/paper/test_reliability_analysis.py",
               "chemgridmap/paper/output/reliability_v1/RESULTS.md",
               "chemgridmap/paper/output/reliability_v1/figure_data/README.md",
               "chemgridmap/paper/output/reliability_v1/manifest.json"}
    # The calling executor may append its own log while this task is running.
    external_logs = [p for p in changes if p.startswith("project1/JCIM_reliability_revision_20261024/codex_tasks/T2/")
                     and p.endswith(".log")]
    unexpected = [p for p in changes if p not in allowed and p not in external_logs]
    current_files = set()
    for directory in (root / "chemgridmap", root / "project1"):
        for path in directory.rglob("*"):
            if path.is_file() and not any(part in {".git", ".venv", "node_modules", "__pycache__", ".pytest_cache", ".DS_Store"}
                                          for part in path.parts):
                current_files.add(str(path.relative_to(root)))
    new_files = sorted(current_files - set(before))
    new_code = {"chemgridmap/paper/" + name for name in (
        "reliability_activity.py", "reliability_cores.py", "reliability_contrasts.py", "reliability_t2b.py",
        "test_reliability_activity.py", "test_reliability_cores.py", "test_reliability_contrasts.py", "check_t2b_outputs.py")}
    output = "chemgridmap/paper/output/reliability_v1/"
    allowed_new_prefixes = [output + name + "/" for name in ("activity", "cores", "cases", "figure_data", "figure1_cli_run")]
    allowed_new_reports = {"project1/JCIM_reliability_revision_20261024/codex_tasks/T2/" + name
                           for name in ("T2b_REPORT.md", "T2_REPORT.md")}
    # These task drafts appeared concurrently and were not written by this task.
    external_task_drafts = {"project1/JCIM_reliability_revision_20261024/codex_tasks/T3/" + name
                           for name in ("T3a_task.md", "T3b_task.md")}
    external_new = [p for p in new_files if p in external_task_drafts]
    unexpected_new = [p for p in new_files if p not in new_code and p not in allowed_new_reports
                      and p not in external_task_drafts
                      and not any(p.startswith(prefix) for prefix in allowed_new_prefixes)]
    prefixes = {}
    for relative, marker in [("chemgridmap/paper/output/reliability_v1/RESULTS.md", RESULT_MARKER),
                             ("chemgridmap/paper/output/reliability_v1/figure_data/README.md", README_MARKER)]:
        prefix = (root / relative).read_text().split(marker, 1)[0]
        prefixes[relative] = hashlib.sha256(prefix.encode()).hexdigest() == before[relative]
    return {"passed": not unexpected and not unexpected_new and not missing and all(prefixes.values()),
            "baseline": str(baseline), "baseline_files": len(before), "changed_existing_files": changes,
            "unexpected_changes": unexpected, "missing_files": missing,
            "new_files": new_files, "unexpected_new_files": unexpected_new,
            "externally_created_task_drafts": [{"path": p, "sha256": sha256(root / p),
                                                "mtime_ns": (root / p).stat().st_mtime_ns}
                                               for p in external_new],
            "externally_appended_executor_logs": external_logs, "original_text_prefix_preserved": prefixes,
            "src_unchanged": not any(p.startswith("chemgridmap/src/") for p in changes),
            "tests_unchanged": not any(p.startswith("chemgridmap/tests/") for p in changes),
            "guard_unchanged": not any("codex_tasks/_guard/" in p for p in changes)}


def audit_rework_scope(baseline, root=ROOT.parent, *, reference_dir=None, external_changes=()):
    """Use complete pre-rework hashes and saved originals for the narrower R1–R3 scope."""
    root, baseline = Path(root), Path(baseline)
    reference_dir = Path(reference_dir) if reference_dir else baseline.parent
    before = json.loads(baseline.read_text())
    changes, missing = [], []
    for relative, digest in before.items():
        path = root / relative
        if not path.is_file():
            missing.append(relative)
        elif sha256(path) != digest:
            changes.append(relative)
    code = {"chemgridmap/paper/" + name for name in (
        "reliability_activity.py", "reliability_cores.py", "reliability_t2b.py",
        "run_reliability_analysis.py", "check_t2b_outputs.py",
        "test_reliability_activity.py", "test_reliability_cores.py", "test_reliability_t2b.py",
        "test_reliability_analysis.py", "test_check_t2b_outputs.py")}
    output = "chemgridmap/paper/output/reliability_v1/"
    reports = "project1/JCIM_reliability_revision_20261024/codex_tasks/T2/"
    allowed = code | {output + name for name in (
        "RESULTS.md", "figure_data/README.md", "figure_data/SI_core_size_tradeoff.csv",
        "manifest.json", "cores/manifest.json", "figure_data/t2b_manifest.json")}
    allowed.add(reports + "T2_REPORT.md")
    new_timing = {output + "cores/" + name for name in (
        "core_rule_timing.csv", "fmcs_full_recomputation.csv", "core_rule_timing_summary.csv",
        "core_timing_metadata.json", "core_timing_checks.json")}
    allowed_new = allowed | new_timing | {reports + "T2b_REWORK_REPORT.md"}
    ignored = {".git", ".venv", "node_modules", "__pycache__", ".pytest_cache", ".DS_Store"}
    current_files = {str(path.relative_to(root))
                     for directory in (root / "chemgridmap", root / "project1")
                     for path in directory.rglob("*")
                     if path.is_file() and not any(part in ignored for part in path.parts)}
    new_files = sorted(current_files - set(before))
    # Explicitly attributed coordinator edits never exempt software or frozen outputs.
    coordinator_prefix = "project1/JCIM_reliability_revision_20261024/codex_tasks/"
    invalid_external = [path for path in external_changes if not path.startswith(coordinator_prefix)
                        or "_guard" in Path(path).parts or ".." in Path(path).parts]
    external = [path for path in changes if path in external_changes and path not in invalid_external]
    unexpected = [path for path in changes if path not in allowed and path not in external]
    unexpected_new = [path for path in new_files if path not in allowed_new]
    prefixes, reference_hashes = {}, {}
    for relative, original_name, marker in (
            (output + "RESULTS.md", "RESULTS.md", RESULT_MARKER),
            (output + "figure_data/README.md", "README.md", README_MARKER)):
        original = reference_dir / original_name
        current_path = root / relative
        reference_hashes[original_name] = original.is_file() and sha256(original) == before.get(relative)
        if original.is_file() and current_path.is_file():
            old_text, new_text = original.read_text(), current_path.read_text()
            prefixes[relative] = (marker in old_text and marker in new_text and
                                  old_text.split(marker, 1)[0] == new_text.split(marker, 1)[0])
        else:
            prefixes[relative] = False
    tradeoff = output + "figure_data/SI_core_size_tradeoff.csv"
    original = reference_dir / "SI_core_size_tradeoff.csv"
    reference_hashes[original.name] = original.is_file() and sha256(original) == before.get(tradeoff)
    original_columns_preserved = False
    if original.is_file() and (root / tradeoff).is_file():
        old = pd.read_csv(original, dtype=str, keep_default_na=False)
        new = pd.read_csv(root / tradeoff, dtype=str, keep_default_na=False)
        # Exact cell strings also protect rounding, row order, and categorical labels.
        original_columns_preserved = (list(new.columns[:len(old.columns)]) == list(old.columns) and
                                      new.loc[:, old.columns].equals(old))
    changed_paths = changes + missing + new_files
    return {"passed": (not unexpected and not unexpected_new and not missing and not invalid_external and all(prefixes.values())
                       and all(reference_hashes.values()) and original_columns_preserved),
            "scope": "T2b rework R1–R3", "baseline": str(baseline), "baseline_files": len(before),
            "reference_dir": str(reference_dir), "reference_hashes_verified": reference_hashes,
            "changed_existing_files": changes, "unexpected_changes": unexpected, "missing_files": missing,
            "new_files": new_files, "unexpected_new_files": unexpected_new,
            "outside_scope_changes_observed": bool(missing or unexpected_new or any(path not in allowed for path in changes)),
            "invalid_external_change_paths": invalid_external,
            "externally_changed_files": [{"path": path, "before_sha256": before[path],
                                          "current_sha256": sha256(root / path),
                                          "mtime_ns": (root / path).stat().st_mtime_ns} for path in external],
            "original_text_prefix_preserved": prefixes,
            "si_original_columns_preserved": original_columns_preserved,
            "src_unchanged": not any(p.startswith("chemgridmap/src/") for p in changed_paths),
            "tests_unchanged": not any(p.startswith("chemgridmap/tests/") for p in changed_paths),
            "guard_unchanged": not any("codex_tasks/_guard/" in p for p in changed_paths)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--rework", action="store_true", help="Audit the narrower R1–R3 file boundary.")
    parser.add_argument("--reference-dir", type=Path, help="Saved original texts and SI table (default: baseline directory).")
    parser.add_argument("--external-change", action="append", default=[], metavar="RELATIVE_PATH",
                        help="Explicitly attribute a concurrent coordinator-file edit (repeatable; rework scope only).")
    parser.add_argument("--scope-output", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if (args.reference_dir or args.external_change or args.rework) and not (args.baseline and args.rework):
        parser.error("--rework, --reference-dir and --external-change require --baseline and --rework")
    if args.baseline:
        result = audit_scope(args.baseline, rework=args.rework, reference_dir=args.reference_dir,
                             external_changes=args.external_change)
        if args.scope_output:
            dump_json(args.scope_output, result)
    else:
        result = audit_outputs()
        if args.output:
            dump_json(args.output, result)
    print(json.dumps(result, indent=2))
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
