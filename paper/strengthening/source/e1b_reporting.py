"""根据 E1-B 已生成的表格写事实说明和可核验清单；不重跑计时。"""
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import platform
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
SOURCE = Path(__file__).resolve().parent
OUT = SOURCE.parent / "output/E1B"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def markdown_table(frame):
    """不增加 tabulate 依赖；CSV 中保留完整精度。"""
    def cell(value):
        if pd.isna(value):
            return "NA"
        if isinstance(value, float):
            return f"{value:.6f}"
        return str(value).replace("|", "\\|")
    rows = ["| " + " | ".join(frame.columns) + " |",
            "| " + " | ".join(["---"] * len(frame.columns)) + " |"]
    rows.extend("| " + " | ".join(cell(v) for v in row) + " |"
                for row in frame.itertuples(index=False, name=None))
    return "\n".join(rows)


def main():
    neighbor = json.loads((OUT / "neighbor_checks.json").read_text())
    cores = json.loads((OUT / "core_fmcs_provenance.json").read_text())
    scope = json.loads((OUT / "scope_check.json").read_text())
    rho = pd.read_csv(OUT / "neighbor_sar_summary.csv")
    errors = pd.read_csv(OUT / "layout_activity_error.csv")
    strata = pd.read_csv(OUT / "core_by_group_size.csv")
    budgets = pd.read_csv(OUT / "fmcs_budget_sweep.csv")
    assert len(rho) == 27 and rho.bootstrap_valid_count.eq(1000).all()
    assert neighbor["frozen_regression"]["passed"] and cores["summary"]["passed"]
    paired = errors.loc[errors.row_type.eq("paired_seed_difference")]
    completed = cores["summary"]["frozen_10s_completed"]
    regression = neighbor["frozen_regression"]
    metadata = cores["metadata"]
    lines = ["# E1-B results", "",
        "This analysis uses the frozen three-target data and coordinates. No layouts, activities, "
        "core-selection rules, or frozen results were changed. Full precision is retained in the CSV files.", "",
        "## 1. Neighborhood recall and displayed-activity error", "",
        "For each molecule, F is its tie-inclusive k=10 Morgan/Tanimoto neighborhood and G is its "
        "tie-inclusive k=10 Euclidean coordinate neighborhood; self is excluded. The relative and absolute "
        "tie tolerances are 1e-9 and 1e-12. The displayed activity y is the frozen pChEMBL median. "
        "Recall is |F intersect G|/|F|; e_grid and e_fp are the respective mean absolute activity "
        "differences. Layout-level means give equal weight to each molecule, not each edge.", "",
        "The prespecified correlation is Spearman rho(recall, e_grid - e_fp). Three main layouts were "
        "evaluated at seeds 7, 21, 42, 84, 126, 168, 210, 252, 294, and 336. Each of six baseline layouts "
        "was additionally evaluated once at occupancy 0.4, seed 42. The occupancy value on projection "
        "rows denotes the associated layout experiment; projection itself is not a lattice.", "",
        "The table reports seed-42 rho and percentile 95% CIs from 1,000 paired molecule bootstrap "
        "replicates, seed 20261031, with neighborhoods fixed. Ranks are recomputed in every replicate. "
        "The RNG is reset for each target/layout. All 27 CIs have 1,000 finite replicates. Main-layout "
        "medians and ranges use 10 seeds; baseline medians/ranges summarize only their single seed.", "",
        markdown_table(rho[["target_label", "layout", "rho_seed42", "rho_ci_low", "rho_ci_high",
                            "n_seeds", "rho_seed_median", "rho_seed_min", "rho_seed_max"]]), "",
        "Each target/seed/layout has five recall bins: [0,0.2), [0.2,0.4), [0.4,0.6), [0.6,0.8), "
        "and [0.8,1]. Bin n, mean e_grid and mean e_fp are in neighbor_sar_bins.csv (540 rows). "
        "Empty bins have n=0 and missing means. Seed-42 molecule-level values and neighborhood sizes "
        "are in neighbor_sar_by_molecule_seed42.csv.", "",
        "The projection_lap and projection_lap_refined baseline coordinates are exactly equal to "
        "minimum_displacement and neighbor_refined, respectively, at seed 42 for all three targets. "
        "Those baseline rows repeat the corresponding main-layout observations.", "",
        "## 2. Activity error across layouts", "",
        "Refinement did not use activity values. Each difference below is neighbor_refined minus "
        "minimum_displacement in mean e_grid. The point estimate averages 10 paired seed differences; "
        "the 95% percentile CI resamples those 10 differences 1,000 times with seed 20261031. "
        "Per-layout and per-seed means are in layout_activity_error.csv.", "",
        markdown_table(paired[["target_label", "n_seeds", "paired_difference_mean", "paired_ci_low", "paired_ci_high"]]), "",
        f"All {regression['rows']} frozen main-layout rows were reproduced. Maximum absolute differences: "
        f"e_grid {regression['neighbor_abs_dpchembl_max_absolute_difference']:.3g}; "
        f"e_fp {regression['fingerprint_neighbor_abs_dpchembl_max_absolute_difference']:.3g}; "
        f"mean recall {regression['recall_k10_max_absolute_difference']:.3g}. "
        f"The {regression['per_molecule_recall_rows']} frozen seed-42 per-molecule recalls also agree "
        f"(maximum absolute difference {regression['per_molecule_recall_max_absolute_difference']:.3g}).", "",
        "The molecule bootstrap does not model dependence induced by overlapping neighborhoods. "
        "Recall and e_grid - e_fp share F and G. The seed-based CI describes layout-seed variability "
        "for the fixed targets and activities. No correlation significance test or permutation p-value "
        "was calculated, and no significance claim is made.", "",
        "## 3. Core rules by group size", "",
        "The 480 groups comprise 196 groups of size 2, 142 of size 3-4, 94 of size 5-9, and 48 "
        "of size at least 10. All 480 full member-identity sets and 1,440 frozen rule rows were verified. "
        "Strict majority means coverage_count > size/2; target attainment means coverage >=0.8. "
        "No core is a categorical state and is excluded from mean coverage and mean heavy-atom count. "
        "Each of most_frequent and coverage has two No core groups. Per-target strata are also in the CSV.", "",
        markdown_table(strata.loc[strata.target.eq("overall"), ["size_stratum", "rule", "groups", "no_core_groups",
            "strict_majority_groups", "target_met_groups", "mean_coverage", "mean_core_heavy_atoms", "benzene_cores"]]), "",
        "Among the 48 groups of size at least 10, strict-majority coverage changes from 29 groups "
        "under most_frequent to 48 under coverage. Coverage of at least 0.8 changes from 13 to 47 groups.", "",
        "## 4. FMCS time budgets", "",
        "All 480 groups were recomputed at each timeout of 1, 2, and 5 integer seconds, sequentially "
        "in one process. The coverage target is 0.8. The unchanged mcs_core function uses "
        "ringMatchesRingOnly, completeRingsOnly, CompareElements and CompareOrder. A returned core is "
        "scored even if the search times out. The reported times are the function's internal wall clock: "
        "FindMCS, SMARTS parsing and member matching; preparation and exports are excluded.", "",
        markdown_table(budgets[["target_label", "timeout_seconds", "groups", "target_met_groups",
            "mean_coverage", "timeout_groups", "total_seconds", "median_seconds"]]), "",
        f"Current machine: {metadata['cpu_model']}; {metadata['platform']}; "
        f"{metadata['logical_cpu_count']} logical CPUs; Python {metadata['python_version']}; "
        f"RDKit {metadata['rdkit_version']}. Timing started at {metadata['started_at_utc']} "
        f"and finished at {metadata['finished_at_utc']}. Load averages at start/end were "
        f"{metadata['load_average_start']} / {metadata['load_average_finish']}. "
        "Other E1-B analysis jobs were paused; unrelated background activity was not controlled.", "",
        "The 10-second rows use the frozen full-recomputation coverage and timeout flags, paired with "
        "core_rule_timing.csv. These were checked against fmcs_full_recomputation.csv. Group-size "
        "stratification uses group_core_comparison.csv, whose FMCS queries originated in an earlier "
        "cache. The frozen coverage counts agree between those sources. The 10-second measurements "
        "were made in a separate run and are not contemporaneous with the new budgets.", "",
        f"After removing the {completed['excluded_timeout_groups']} timed-out groups from both rules "
        f"in the frozen 10-second timing table, {completed['groups']} paired groups remain. "
        f"FMCS total time is {completed['fmcs_total_seconds']:.9f} s and coverage-rule total time "
        f"for the same groups is {completed['coverage_total_seconds_same_groups']:.9f} s. "
        f"The descriptive ratio FMCS / coverage is {completed['fmcs_divided_by_coverage']:.9f}.", "",
        "## 5. Provenance", "",
        "manifest.json records explicit input and source SHA-256 hashes, output SHA-256 hashes, "
        "software versions, fixed parameters and timing metadata. neighbor_checks.json and "
        "core_fmcs_provenance.json retain the numerical checks. core_group_members.csv and "
        "fmcs_budget_by_group.csv retain the complete group membership and budget-level raw evidence.", "",
        "The whole-workspace scope check did not pass: concurrent changes were observed in Git "
        "metadata, a new chemgridmap_public worktree, and an external E1B run log. The CGM working-tree "
        "status remained identical at 45 lines, E1-A hashes were unchanged, and no changes outside "
        "the E1-B allowlist were observed in the frozen CGM non-Git files. The observations are "
        "retained in scope_check.json; attribution of the concurrent activity awaits coordinator review.", ""]
    (OUT / "RESULTS_E1B.md").write_text("\n".join(lines))
    inputs = dict(neighbor["inputs_sha256"])
    for name in cores["inputs"] + ["paper/strengthening/SELECTION_RULES.md"]:
        inputs[name] = sha256(ROOT / name)
    assert all(sha256(ROOT / name) == digest for name, digest in inputs.items())
    manifest = {"task": "E1-B", "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs_sha256": inputs,
        "source_sha256": {str(path.relative_to(ROOT)): sha256(path)
                          for path in sorted(SOURCE.glob("e1b_*.py")) + [SOURCE / "test_e1b.py"]},
        "parameters": {"neighbors": neighbor["parameters"], "cores": cores["parameters"]},
        "versions": {name: version(name) for name in ("numpy", "pandas", "scipy", "scikit-learn", "rdkit")},
        "python_version": platform.python_version(), "python_executable": sys.executable,
        "workspace_validation": {"whole_workspace_scope_passed": scope["passed"],
            "unexpected_change_count": len(scope["unexpected_changes"]),
            "cgm_status_identical": scope["cgm_status_identical"],
            "cgm_status_lines": scope["cgm_status_lines"],
            "e1a_sha256_unchanged": not scope["e1a_changed"],
            "scope_note": "Concurrent Git/worktree and external run-log changes retained for coordinator review."},
        "timing_metadata": metadata,
        "outputs_sha256": {path.name: sha256(path) for path in sorted(OUT.iterdir())
                           if path.is_file() and path.name != "manifest.json"}}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    print(f"REPORTING PASS: English results; {len(inputs)} input hashes, "
          f"{len(manifest['source_sha256'])} source hashes, {len(manifest['outputs_sha256'])} output hashes")


if __name__ == "__main__":
    main()
