"""T2b case, figure tables, result reconciliation, and collection provenance."""
from __future__ import annotations

from datetime import datetime, timezone
import csv
import json
import math
from statistics import NormalDist
from pathlib import Path
import subprocess
import sys
import hashlib

import numpy as np
import pandas as pd

import chemgridmap
from chemgridmap.evidence import document_values
from chemgridmap.provenance import environment_versions
from reliability_analysis import ROOT, FROZEN, OUTPUT, TARGETS, SEEDS, dump_json, sha256
from reliability_figure_data import markdown_table

REFERENCE = ROOT.parent / "project1/JCIM_review_innovation_20261024/results/headline_numbers.json"
RESULT_MARKER = "<!-- T2b generated results -->"
README_MARKER = "<!-- T2b generated column definitions -->"


def load_json(path):
    return json.loads(Path(path).read_text())


def append_section(path, marker, text):
    """Preserve the original prefix and replace only this generated section."""
    path = Path(path)
    original = path.read_text(encoding="utf-8") if path.exists() else ""
    prefix = original.split(marker, 1)[0]
    path.write_text(prefix + marker + "\n\n" + text.rstrip() + "\n", encoding="utf-8")


def run_activity_case(output_root=OUTPUT):
    root = Path(output_root)
    records_path = FROZEN / "chembl205/chembl205_retained_records.csv"
    records = pd.read_csv(records_path)
    selected = records.loc[records.molecule_identity_key.eq("CHEMBL20")]
    units = document_values(selected).sort_values(["document_value", "document"]).reset_index(drop=True)
    units.insert(0, "unit_rank", np.arange(1, len(units) + 1))
    units["n_documents"] = units.documents.str.split(";").map(len)
    units["is_merged"] = units.n_documents.gt(1)
    units["is_innocenti_2008"] = units.documents.str.split(";").map(lambda x: "CHEMBL1141531" in x)
    annotations_path = root / "activity/chembl205_annotations.csv"
    annotations = pd.read_csv(annotations_path)
    row = annotations.loc[annotations.molecule_identity_key.eq("CHEMBL20")].iloc[0]
    positions_path = root / "cases/G002_CHEMBL20_positions_seed42.csv"
    position = pd.read_csv(positions_path).set_index("molecule_identity_key").loc["CHEMBL20"]
    summary = {key: None if pd.isna(value) else value.item() if isinstance(value, np.generic) else value
               for key, value in row.items()}
    summary.update({"n_records": len(selected), "n_units": len(units), "seed": 42,
                    "refined_grid_col": int(position.refined_grid_col),
                    "refined_grid_row": int(position.refined_grid_row),
                    "lower_threshold": 6., "upper_threshold": 7.,
                    "default_model": load_json(root / "activity/adoption.json")["default_model"]})
    checks = {"records_54": len(selected) == 54, "units_27": len(units) == 27,
              "record_accounting": int(units.n_records.sum()) == len(selected),
              "median_7_185": bool(np.isclose(selected.pchembl_value.median(), 7.185)),
              "innocenti_records_25": int(units.loc[units.is_innocenti_2008, "n_records"].sum()) == 25}
    if summary["default_model"] == "G":
        checks.update({"balanced_value": bool(np.isclose(row.document_balanced_pchembl, 7.044074074074074)),
                       "active_probability": bool(np.isclose(row.p_consensus_active, .6132102836102352)),
                       "uncertain": row.class_status == "uncertain"})
    summary["checks"] = checks
    summary["passed"] = all(checks.values())
    if not summary["passed"]:
        raise ValueError(f"CHEMBL20 case checks failed: {checks}")
    units["lower_threshold"] = 6.
    units["upper_threshold"] = 7.
    units["document_balanced_pchembl"] = row.document_balanced_pchembl
    units.to_csv(root / "cases/CHEMBL20_evidence_units.csv", index=False)
    dump_json(root / "cases/CHEMBL20_activity_summary.json", summary)
    inputs = [records_path, annotations_path, positions_path, root / "activity/adoption.json"]
    write_new_manifest(root / "cases/t2b_manifest.json", inputs,
                       [root / "cases/CHEMBL20_evidence_units.csv", root / "cases/CHEMBL20_activity_summary.json",
                        root / "cases/core_rule_examples.csv"], {"scope": "T2b cases", "seed": 42})
    return summary


def write_new_manifest(path, inputs, outputs, parameters):
    """Write a new manifest without touching historical subdirectory manifests."""
    sources = sorted(p for base in (ROOT / "src", ROOT / "paper") for p in base.rglob("*.py")
                     if "output" not in p.parts and "__pycache__" not in p.parts)
    diff = subprocess.check_output(["git", "diff", "release-v0.4.0"], cwd=ROOT)
    result = {"collected_at_utc": datetime.now(timezone.utc).isoformat(),
              "inputs_sha256": {str(Path(p).resolve()): sha256(p) for p in inputs},
              "parameters": parameters,
              "environment": {**environment_versions(), "python_executable": sys.executable,
                              "chemgridmap.__version__": chemgridmap.__version__,
                              "chemgridmap.__file__": chemgridmap.__file__},
              "code": {"git_diff_release_v0.4.0_sha256": hashlib.sha256(diff).hexdigest(),
                       "source_sha256": {str(p.resolve()): sha256(p) for p in sources}},
              "outputs_sha256": {str(Path(p).resolve()): sha256(p) for p in outputs}}
    dump_json(path, result)
    return result


def flatten(value, path=()):
    if isinstance(value, dict):
        for key, child in value.items():
            yield from flatten(child, (*path, key))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            yield from flatten(child, (*path, str(index)))
    else:
        yield path, value


def reconcile_headlines(root):
    """Compare every historical headline leaf, keeping changed scopes explicit."""
    old = load_json(REFERENCE)
    # Preserve historical paths/values; translate labels only when querying current data.
    def current_label(label):
        return "Thrombin" if label == "D3" else label
    current, reasons, sources = {}, {}, {}

    def put(path, value, reason, source):
        current[path], reasons[path], sources[path] = value, reason, source

    metrics = pd.read_csv(root / "neighborhoods/layout_metrics.csv")
    stability = pd.read_csv(root / "stability/seed_pair_jaccard.csv")
    rename = {"spearman_projection_vs_grid": "spearman_projection_vs_layout",
              "grid_neighbor_abs_dpchembl": "neighbor_abs_dpchembl"}
    for label, section in old["grid_refinement"].items():
        frame = metrics[metrics.target_label.eq(current_label(label))]
        means = frame.groupby("layout").mean(numeric_only=True)
        for key in section:
            value = None
            if key == "seeds":
                value = frame.seed.nunique()
            elif key == "projection_recall":
                value = means.loc["projection", "recall_k10"]
            elif key == "fraction_of_grid_loss_recovered":
                p, b, a = means.loc[["projection", "minimum_displacement", "neighbor_refined"], "recall_k10"]
                value = (a - b) / (p - b)
            else:
                for layout in ("minimum_displacement", "neighbor_refined"):
                    suffix = "_" + layout
                    if key.endswith(suffix):
                        column = key[:-len(suffix)]
                        if column == "seed_jaccard":
                            value = stability.loc[stability.target_label.eq(current_label(label)) & stability.layout.eq(layout),
                                                  "neighbor_jaccard_k10"].mean()
                        else:
                            value = means.loc[layout, rename.get(column, column)]
            put(("grid_refinement", label, key), value,
                "T2a uses 10 map seeds and refinement seed equal to map seed; the reference used 5/3/3 maps and refinement seed 0. Runtime also varies by execution.",
                "neighborhoods/layout_metrics.csv; stability/seed_pair_jaccard.csv")
    scale = load_json(root / "scale/full_herg_summary.json")["frozen_refinement"]
    for key in old["grid_refinement_full_herg"]:
        put(("grid_refinement_full_herg", key), scale[key],
            "Same frozen full hERG layout; T2a rerun. Wall times vary by execution.", "scale/full_herg_summary.json")
    scores = pd.read_csv(root / "activity/lodo_scores.csv")
    for index, row in enumerate(old["activity_lodo_scores"]):
        target = {label: target for target, label in TARGETS.items()}.get(current_label(row["target"]), "all")
        found = scores[scores.method.eq(row["method"]) & scores.target.eq(target)].iloc[0]
        for key, value in row.items():
            put(("activity_lodo_scores", str(index), key), row[key] if key == "target" else found[key],
                "Same G/B0/MF/MP protocol; reference score table was rounded to four decimals.", "activity/lodo_scores.csv")
    annotations = pd.read_csv(root / "activity/mapped_entry_confidence.csv")
    predictions = pd.read_csv(root / "activity/lodo_predictions.csv")
    for label in old["activity_fits"]:
        frame = annotations[annotations.target_label.eq(current_label(label))]
        pred = predictions[predictions.target_label.eq(current_label(label)) & predictions.method.eq("G")]
        prob = pred[["p0", "p1", "p2"]].to_numpy()
        high = prob.max(axis=1) >= .8
        high_values = {"n": int(high.sum()), "mean_predicted": prob[high].max(axis=1).mean(),
                       "observed_agreement": np.mean(prob[high].argmax(axis=1) == pred.loc[high, "observed"])}
        for key, value in high_values.items():
            put(("activity_high_confidence_bin", label, key), value, "Same G LODO maximum-probability interval.",
                "activity/lodo_predictions.csv")
        records = pd.read_csv(FROZEN / next(t for t, name in TARGETS.items() if name == current_label(label)) /
                              f"{next(t for t, name in TARGETS.items() if name == current_label(label))}_retained_records.csv")
        from chemgridmap.evidence import fit_between_unit_sd
        fit = fit_between_unit_sd(document_values(records))
        for key in old["activity_fits"][label]:
            put(("activity_fits", label, key), fit[key], "Same full retained-record Gaussian fit.", "revision_v6 retained records")
    ace = annotations[annotations.molecule_identity_key.eq("CHEMBL20") & annotations.target.eq("chembl205")].iloc[0]
    for key, value in old["activity_acetazolamide"].items():
        put(("activity_acetazolamide", key), "CA II" if key == "target" else ace[key],
            "Final default annotation on the same retained records; probabilities describe the consensus.",
            "cases/CHEMBL20_activity_summary.json")
    for index, row in enumerate(old["activity_mapped_summary"]):
        frame = annotations[annotations.target_label.eq(current_label(row["target"]))]
        if row["scope"] == "one_unit":
            frame = frame[frame.n_evidence_units.eq(1)]
        elif row["scope"] == "multi_unit":
            frame = frame[frame.n_evidence_units.ge(2)]
        values = {"target": row["target"], "scope": row["scope"], "entries": len(frame),
                  "supported_fraction": frame.class_confidence.ge(.8).mean(),
                  "median_confidence": frame.class_confidence.median(),
                  "class_changed_by_document_balance": frame.class_changed_by_document_balance.sum()}
        values.update({"supported_" + c: frame.loc[frame.activity_class.eq(c), "class_confidence"].ge(.8).mean()
                       for c in ("inactive", "medium", "active")})
        for key in row:
            put(("activity_mapped_summary", str(index), key), values[key],
                "Same mapped population and support cutoff; reference table was rounded to four decimals.",
                "activity/mapped_entry_confidence.csv")
    contrasts = pd.read_csv(root / "activity/neighbor_contrasts.csv")
    # Secondary pairs retain the old k=10 definition; the main figure uses eight-adjacency.
    for index, row in enumerate(old["activity_neighbor_contrasts"]):
        subset = contrasts[(contrasts.target_label.eq(current_label(row["target"]))) & contrasts.layout.eq(row["layout"]) & contrasts.seed.eq(42)]
        subset = subset[subset.definition.eq("tie_aware_k10_union")]
        if len(subset):
            found = subset.iloc[0]
            for key in row:
                value = row[key] if key in {"target", "layout"} else found[key]
                put(("activity_neighbor_contrasts", str(index), key), value,
                    "Secondary k=10 symmetric union. Minimum-displacement seed 42 is unchanged; refined seed is now 42 rather than 0. Reference expected counts were rounded to one decimal.",
                    "activity/neighbor_contrasts.csv")
    cores = pd.read_csv(root / "cores/group_core_comparison.csv")
    for label, section in old["core_selection"].items():
        frame = cores[cores.target_label.eq(current_label(label))]
        for key in section:
            path = ("core_selection", label, key)
            if key.endswith("seconds"):
                put(path, None, "Historical whole-run timing was not repeated: FMCS uses the authorized cache with a 20-group recomputation audit; timing is not comparable.",
                    "cores/fmcs_recomputation.csv")
                continue
            if key == "groups":
                value = frame.group.nunique()
            elif key == "aware_same_as_current":
                pivot = frame[frame.rule.isin(["most_frequent", "coverage"])].pivot(index="group", columns="rule", values="core_smiles")
                value = (pivot.most_frequent.notna() & pivot.coverage.notna() &
                         pivot.most_frequent.eq(pivot.coverage)).sum()
            else:
                prefix, stat = key.split("_", 1)
                rule = {"current": "most_frequent", "aware": "coverage", "mcs": "fmcs"}[prefix]
                selected = frame[frame.rule.eq(rule)]
                value = {"majority": lambda: selected.strict_majority.sum(),
                         "ge_0.8": lambda: selected.target_met.sum(),
                         "median_heavy_atoms": lambda: selected.loc[selected.has_core, "core_heavy_atoms"].median(),
                         "benzene": lambda: selected.is_benzene.sum(),
                         "single_ring": lambda: selected.is_single_ring.sum()}[stat]()
            put(path, value, "Same 480 Butina groups. FMCS is a SMARTS query; its member-fragment SMILES is a display aid.",
                "cores/group_core_comparison.csv")
    rows = []
    for path, previous in flatten(old):
        if path not in current:
            raise ValueError(f"Unreconciled historical headline: {path}")
        value = current[path]
        numeric = isinstance(previous, (int, float)) and isinstance(value, (int, float, np.number))
        delta = float(value - previous) if numeric and np.isfinite(previous) and np.isfinite(value) else None
        same = previous == value if not numeric else bool(np.isclose(previous, value, atol=1e-12, rtol=1e-12))
        rows.append({"headline_path": "/".join(path), "reference_value": previous, "current_value": value,
                     "difference": delta, "status": "not_remeasured" if value is None else "same" if same else "changed",
                     "reason": reasons[path], "current_source": sources[path]})
    table = pd.DataFrame(rows)
    table.to_csv(root / "activity/headline_comparison.csv", index=False)
    return table


def run_t2b_figure_data(output_root=OUTPUT):
    root = Path(output_root)
    out = root / "figure_data"
    from reliability_contrasts import run_figure1_cli
    cli_check = run_figure1_cli(root)
    annotations = pd.read_csv(root / "activity/mapped_entry_confidence.csv")
    base = pd.read_csv(out / "Figure_1ABC_layout_seed42_refined.csv")
    columns = ["target", "molecule_identity_key", "activity_pchembl", "activity_class",
               "class_confidence", "class_status", "n_evidence_units", "discordant_document_pairs"]
    positions = base.merge(annotations[columns], on=["target", "molecule_identity_key"], validate="one_to_one")
    cores = pd.read_csv(root / "cores/group_core_comparison.csv")
    members = pd.read_csv(root / "cores/G001_G002_members.csv")
    g002 = members[members.group.eq("G002")].copy()
    units = pd.read_csv(root / "cases/CHEMBL20_evidence_units.csv")
    calibration = pd.read_csv(root / "activity/calibration.csv")
    sensitivity = pd.read_csv(root / "activity/threshold_sensitivity.csv")
    scores = pd.read_csv(root / "activity/lodo_scores.csv")
    rules = pd.read_csv(root / "activity/adoption_rules.csv")
    adoption = load_json(root / "activity/adoption.json")
    default = adoption["default_model"]
    main_calibration = calibration[calibration.method.isin([default, "B0"]) & ~calibration.target.eq("all")]
    support = sensitivity[sensitivity.evidence_policy.eq("merged") & sensitivity.scope.eq("all") &
                          sensitivity.threshold.eq(.8) & ~sensitivity.display_class.eq("all") & ~sensitivity.target.eq("all")].copy()
    contrasts = pd.read_csv(root / "activity/neighbor_contrasts.csv")
    refined = contrasts[contrasts.layout.eq("neighbor_refined") & contrasts.definition.eq("eight_connected")]
    stats = ["neighbor_pairs", "displayed_contrasts", "expected_true_contrasts_among_displayed",
             "supported_contrasts_p_ge_0.8", "same_colour_pairs", "expected_hidden_contrasts_same_colour"]
    long = refined.melt(id_vars=["target", "target_label", "seed"], value_vars=stats,
                        var_name="statistic", value_name="value")
    summary = long.groupby(["target", "target_label", "statistic"]).value.agg(
        seed_mean="mean", seed_min="min", seed_max="max", n_seeds="size").reset_index()
    main_contrasts = long[long.seed.eq(42)].merge(summary, on=["target", "target_label", "statistic"], validate="one_to_one")
    main_contrasts["definition"] = "eight_connected"
    main_contrasts["layout"] = "neighbor_refined"
    core_plot = cores[cores.rule.isin(["most_frequent", "coverage"])][
        ["target", "target_label", "group", "size", "rule", "has_core", "coverage_fraction", "strict_majority"]].copy()
    core_plot["prism_column"] = core_plot.target_label + " | " + core_plot.rule
    wide_rules = rules.pivot(index="method", columns="rule", values="passed").add_prefix("adoption_").reset_index()
    all_models = scores.merge(wide_rules, on="method", how="left", validate="many_to_one")
    all_models["selected_default"] = all_models.method.eq(default)
    all_models["eligible_for_adoption"] = all_models.method.isin(["O", "X"])
    tables = {
        "Figure_1ABC_map_seed42_refined.csv": positions,
        "Figure_1D_G002_core.csv": g002,
        "Figure_2B_CHEMBL20_evidence_units.csv": units,
        "Figure_2C_calibration.csv": main_calibration,
        "Figure_2D_supported_classes.csv": support,
        "Figure_2E_neighbor_contrasts.csv": main_contrasts,
        "Figure_3D_core_coverage.csv": core_plot,
        "SI_calibration_by_target.csv": calibration,
        "SI_threshold_sensitivity.csv": sensitivity,
        "SI_document_offsets.csv": pd.read_csv(root / "activity/document_offsets.csv"),
        "SI_core_target_sensitivity.csv": pd.read_csv(root / "cores/core_target_sensitivity_summary.csv"),
        "SI_core_size_tradeoff.csv": pd.read_csv(root / "cores/core_size_tradeoff.csv"),
        "SI_activity_models.csv": all_models}
    if len(positions) != 889 or len(g002) != 28 or len(units) != 27 or len(core_plot) != 960:
        raise ValueError("Figure row counts do not match their analysis populations.")
    if len(support) != 9 or len(main_contrasts) != 18:
        raise ValueError("Class-support or neighbor-contrast figure table is incomplete.")
    for name, table in tables.items():
        table.to_csv(out / name, index=False)
    write_figure_readme(root, tables)
    checks = {"passed": bool(cli_check["passed"]), "default_model": default,
              "row_counts": {name: len(table) for name, table in tables.items()},
              "figure1_cli_consistency": cli_check, "core_prism_columns": list(core_plot.prism_column.unique())}
    dump_json(out / "t2b_checks.json", checks)
    write_new_manifest(out / "t2b_manifest.json",
        [root / "activity/adoption.json", root / "activity/lodo_scores.csv", root / "activity/adoption_rules.csv",
         root / "activity/mapped_entry_confidence.csv", root / "activity/calibration.csv",
         root / "activity/threshold_sensitivity.csv", root / "activity/document_offsets.csv",
         root / "activity/neighbor_contrasts.csv", root / "cores/group_core_comparison.csv",
         root / "cores/G001_G002_members.csv", root / "cores/core_target_sensitivity_summary.csv",
         root / "cores/core_size_tradeoff.csv", root / "cases/CHEMBL20_evidence_units.csv",
         out / "Figure_1ABC_layout_seed42_refined.csv", root / "figure1_cli_run/consistency_check.json"],
        [out / name for name in tables] + [out / "t2b_checks.json", out / "README.md"],
        {"scope": "T2b", "default_model": default, "seed": 42, "primary_neighbor_definition": "eight_connected"})
    return checks


COLUMN_MEANINGS = {
    "groups": "Number of nonsingleton Butina groups",
    "no_core": "Groups with no nonempty candidate core, count",
    "strict_majority_groups": "Groups with core coverage strictly greater than one half, count",
    "target_met_groups": "Groups whose core reaches the specified coverage target, count",
    "mean_coverage": "Mean coverage fraction among groups with a core; No core is excluded",
    "median_coverage": "Median coverage fraction among groups with a core",
    "median_core_heavy_atoms": "Median core heavy-atom count among groups with a core",
    "median_core_ring_count": "Median RDKit core ring count among groups with a core",
    "benzene_cores": "Groups whose core is a single benzene ring, count",
    "single_ring_cores": "Groups whose core has exactly one ring, count",
    "most_frequent_core_smiles": "Selected most-frequent member scaffold SMILES",
    "coverage_core_smiles": "Selected coverage-target member scaffold or ring-fragment SMILES",
    "target": "Lowercase ChEMBL target ID; all or overall denotes the pooled analysis population",
    "target_label": "CA II, Thrombin, hERG, All, or Overall", "seed": "Map and refinement random seed",
    "molecule_identity_key": "Stable molecular entry identity", "molecule_id": "Display molecule identifier",
    "structure_group": "Butina group identifier", "group": "Butina group identifier",
    "scaffold_smiles": "Typed Bemis-Murcko scaffold SMILES",
    "refined_grid_col": "Optimized zero-based grid column, cells", "refined_grid_row": "Optimized zero-based grid row, cells",
    "is_G002": "Membership in CA II G002", "is_CHEMBL20": "Whether the entry is acetazolamide",
    "recall_k10": "Per-entry Morgan k=10 neighbor recall, including distance ties, dimensionless",
    "activity_pchembl": "Median retained-record pChEMBL, negative base-10 logarithm of molar IC50",
    "activity_class": "Displayed median class: inactive <=6, medium >6 and <7, active >=7",
    "class_confidence": "Consensus probability of the displayed class, dimensionless",
    "class_status": "supported if class_confidence >=0.8, otherwise uncertain",
    "n_evidence_units": "Number of evidence units after merging equal rounded document medians",
    "discordant_document_pairs": "Flagged document pairs sharing this entry, count",
    "canonical_smiles": "Canonical isomeric molecular SMILES",
    "matches_most_frequent_core": "Chirality-aware match to the most frequent scaffold core",
    "matches_coverage_core": "Chirality-aware match to the coverage-target core",
    "matches_fmcs_core": "Match to the authoritative FMCS SMARTS query",
    "matches_exact_scaffold": "Entry scaffold equals the selected most-frequent scaffold",
    "unit_rank": "One-based evidence-unit order, sorted by value then document ID",
    "entry": "Stable molecular entry identity", "document": "Representative document ID (first lexicographically)",
    "documents": "All constituent document IDs separated by semicolons",
    "document_value": "Evidence-unit pChEMBL: document median; equal rounded medians merged",
    "n_records": "Retained source records represented by the unit, count",
    "n_documents": "Constituent documents, count", "is_merged": "True when multiple documents form this unit",
    "is_innocenti_2008": "Unit contains CHEMBL1141531, including constituent documents",
    "lower_threshold": "Inactive boundary 6, pChEMBL", "upper_threshold": "Active boundary 7, pChEMBL",
    "document_balanced_pchembl": "Mean evidence-unit pChEMBL, each unit equally weighted",
    "method": "B0, G, MF, MP, O, X, or N; definitions in RESULTS section 2",
    "evaluation_set": "merged for common 820 predictions; unmerged for N's separate sensitivity population",
    "bin": "Zero-based five-bin index; 5 is the B0 point at probability 1",
    "bin_low": "Inclusive lower probability boundary", "bin_high": "Upper probability boundary; last bin includes 1",
    "n": "Number of predictions or entries in this row's explicitly defined population",
    "x": "Mean probability assigned to the B0 displayed class in this bin; B0 x=1",
    "y": "Observed fraction agreeing with the B0 displayed class in this bin",
    "plot_eligible": "Plot only if n>=8; all bins remain in the source table",
    "evidence_policy": "merged: equal rounded medians collapse; unmerged: one unit per document",
    "scope": "all entries, one_unit, or multi_unit under this evidence policy",
    "display_class": "Displayed class inactive/medium/active; all pools the classes",
    "threshold": "Consensus class-support cutoff, dimensionless",
    "supported_n": "Entries whose displayed class probability meets the cutoff, count",
    "supported_fraction": "supported_n divided by n; missing if n=0",
    "class_changed_by_document_balance": "Entries whose balanced-mean class differs from displayed median class, count",
    "discordant_document_entries": "Mapped entries with discordant_document_pairs >0, count",
    "statistic": "Neighbor-pair statistic; definitions in the T2b introduction above", "value": "Statistic at optimized seed 42",
    "seed_mean": "Mean statistic across all 10 seeds", "seed_min": "Minimum across all 10 seeds",
    "seed_max": "Maximum across all 10 seeds", "n_seeds": "Number of contributing seeds, 10",
    "definition": "eight_connected: occupied cells at Chebyshev distance 1; unordered unique pairs",
    "layout": "neighbor_refined: default neighbor-optimized layout",
    "size": "Group member count", "rule": "most_frequent, coverage, or fmcs",
    "has_core": "False denotes No core; do not interpret as zero coverage",
    "coverage_fraction": "Members matching the selected core divided by group size; missing for No core",
    "strict_majority": "Core matches strictly more than half the group",
    "prism_column": "Target and rule label for six Prism columns",
    "document_a": "Lexicographically first document in a pair", "document_b": "Second document in the pair",
    "shared_entries": "Entries represented by both documents, count",
    "median_offset": "Median document B minus A difference, pChEMBL",
    "mean_offset": "Mean B minus A difference, pChEMBL", "std_offset": "Sample SD of B minus A, pChEMBL",
    "sign_agreement": "Fraction of shared differences with the median's sign; zero does not agree",
    "systematic_offset": "True if shared_entries>=3, abs(median_offset)>=1, sign_agreement>=0.8",
    "coverage_target": "Required fraction of group members matched by a candidate core",
    "n_predictions": "Held-out evidence-unit predictions, count", "n_entries": "Distinct evaluated target-entry pairs, count",
    "accuracy": "Fraction of argmax class predictions matching the held-out class",
    "brier": "Mean sum over three classes of squared probability error, range 0 to 2",
    "log_loss": "Mean negative natural logarithm of held-out class probability, clipped below at 1e-6",
    "ece": "Maximum-probability calibration error using ten equal-width bins",
    "mean_confidence": "Mean maximum class probability",
    "high_confidence_n": "Predictions with maximum probability >=0.8, count",
    "high_confidence_mean_probability": "Mean maximum probability among high-confidence predictions",
    "high_confidence_accuracy": "Observed argmax accuracy among high-confidence predictions",
    "selected_default": "True for the model selected by all preregistered adoption rules",
    "eligible_for_adoption": "Only O and X are candidates; N uses a different evaluation set",
}


def write_figure_readme(root, tables):
    descriptions = dict(COLUMN_MEANINGS)
    descriptions.update({"adoption_" + rule: text for rule, text in {
        "overall_brier": "Candidate pooled Brier <= G, without tolerance",
        "overall_ece": "Candidate pooled ECE <= G, without tolerance",
        "herg_high_confidence_accuracy": "Candidate hERG high-confidence accuracy >= G +0.05",
        "caii_ece": "Candidate CA II ECE <= G +0.02",
        "thrombin_ece": "Candidate Thrombin ECE <= G +0.02"}.items()})
    text = ["## T2b finalized figure sources", "",
        "These files supplement the frozen T2a tables above. All probabilities, proportions and coverage values are dimensionless. "
        "No core is a category (has_core=False, coverage_fraction missing), not a measured zero. "
        "Empty groups/bins have missing fractions. Missing adoption flags mean the model is not an O/X candidate.", "",
        "Figure 2B uses pChEMBL boundaries 6 and 7 and the document-balanced value recorded as a constant column. "
        "Figure 2C and SI calibration use probability of the displayed B0 class; hide plot_eligible=False (n<8). "
        "ECE and high-confidence metrics instead use maximum probability. The B0 point is (1, B0 accuracy).", "",
        "Figure 2E counts unique unordered occupied-cell pairs at Chebyshev distance 1. neighbor_pairs is the total number of these pairs; "
        "same_colour_pairs counts pairs with equal displayed classes; displayed_contrasts counts different displayed classes; "
        "expected_true_contrasts_among_displayed sums 1 minus the dot product of the two consensus probability vectors over those pairs; "
        "supported_contrasts_p_ge_0.8 requires a displayed difference and probability of a true difference >=0.8. "
        "expected_hidden_contrasts_same_colour uses the same expectation for same-colour pairs. These expectations assume independence "
        "between entry consensus values and are model-based expected counts, not validated biological differences.", "",
        "Figure 3D contains six target/rule columns for Prism. FMCS and threshold tradeoffs are in SI. "
        "FMCS coverage is evaluated using its SMARTS query; a member-fragment SMILES is only a drawing aid. "
        "Figure 1 maps are in ../figure1_cli_run/: the command explicitly uses UMAP because the bare CLI projection default is PCA. "
        "Default refinement and support settings were retained and all 889 coordinates/support values were compared.", ""]
    source_map = {"Figure_1ABC_map_seed42_refined.csv": "T2a layout table joined to ../activity/mapped_entry_confidence.csv",
                  "Figure_1D_G002_core.csv": "../cores/G001_G002_members.csv",
                  "Figure_2B_CHEMBL20_evidence_units.csv": "../cases/CHEMBL20_evidence_units.csv",
                  "Figure_2C_calibration.csv": "../activity/calibration.csv, final default and B0 only",
                  "Figure_2D_supported_classes.csv": "../activity/threshold_sensitivity.csv, merged/all scope/0.8 cutoff",
                  "Figure_2E_neighbor_contrasts.csv": "../activity/neighbor_contrasts.csv, optimized/eight_connected",
                  "Figure_3D_core_coverage.csv": "../cores/group_core_comparison.csv, most_frequent and coverage",
                  "SI_calibration_by_target.csv": "../activity/calibration.csv",
                  "SI_threshold_sensitivity.csv": "../activity/threshold_sensitivity.csv",
                  "SI_document_offsets.csv": "../activity/document_offsets.csv; four counts in ../activity/document_pair_counts.csv",
                  "SI_core_target_sensitivity.csv": "../cores/core_target_sensitivity_summary.csv",
                  "SI_core_size_tradeoff.csv": "../cores/core_size_tradeoff.csv",
                  "SI_activity_models.csv": "../activity/lodo_scores.csv joined to ../activity/adoption_rules.csv"}
    for name, frame in tables.items():
        missing = set(frame) - set(descriptions)
        if missing:
            raise ValueError(f"Undocumented T2b figure columns in {name}: {sorted(missing)}")
        text += [f"### {name} ({len(frame)} rows)", "", "Source: " + source_map[name] + ".", "",
                 "| Column | Definition / unit |", "| --- | --- |"]
        text += [f"| {column} | {descriptions[column]} |" for column in frame]
        text.append("")
    append_section(root / "figure_data/README.md", README_MARKER, "\n".join(text))


def run_results(output_root=OUTPUT):
    root = Path(output_root)
    comparison = reconcile_headlines(root)
    scores = pd.read_csv(root / "activity/lodo_scores.csv")
    rules = pd.read_csv(root / "activity/adoption_rules.csv")
    adoption = load_json(root / "activity/adoption.json")
    support = pd.read_csv(root / "activity/threshold_sensitivity.csv")
    counts = pd.read_csv(root / "activity/document_pair_counts.csv")
    pairs = pd.read_csv(root / "activity/document_offsets.csv")
    contrasts = pd.read_csv(root / "activity/neighbor_contrasts.csv")
    tradeoff = pd.read_csv(root / "cores/core_size_tradeoff.csv")
    core_sensitivity = pd.read_csv(root / "cores/core_target_sensitivity_summary.csv")
    core_checks = load_json(root / "cores/core_checks.json")
    ace = load_json(root / "cases/CHEMBL20_activity_summary.json")
    figure_check = load_json(root / "figure_data/t2b_checks.json")
    primary = contrasts.query("seed == 42 and layout == 'neighbor_refined' and definition == 'eight_connected'")
    pooled = scores[scores.target.eq("all")]
    main_support = support.query("threshold == 0.8 and scope == 'all' and evidence_policy == 'merged'")
    summary_support = main_support[main_support.display_class.eq("all")]
    class_support = main_support[~main_support.display_class.eq("all") & ~main_support.target.eq("all")]
    lines = ["## 2. Evidence models and class support", "",
        "Sources: `activity/lodo_predictions.csv`, `lodo_scores.csv`, `model_fits.json`, `ridge_inner_cv.csv`, "
        "`ridge_document_offsets.csv`, `adoption_rules.csv`, `adoption.json` (all under `activity/`).", "",
        "The complete retained-record populations, not only the mapped subsets, fit each target's noise. "
        "Entries with at least two merged evidence units were randomized into ten folds (seed 20261024). "
        "Each fold fits only training-entry units, then leaves each test-entry unit out in turn. B0 is the class "
        "of the remaining-unit median with probability one; G is Gaussian; MF fits moderated-t variance likelihood; "
        "MP fits moderated-t predictive likelihood. Log loss clips probabilities at 1e-6; multiclass Brier sums "
        "three squared class errors; ECE uses ten equally wide maximum-probability bins.", "",
        "O fits unmerged repeated-entry document units with unpenalized entry means and ridge-penalized document offsets. "
        "Five-fold inner entry cross-fitting selects lambda from {0.1,0.3,1,3,10,30,100} using held-out-unit squared error. "
        "The corrected unmerged training units fit tau_O by the same pooled-SD method as G. A merged test unit uses its "
        "representative document (the lexicographically first constituent) for the prescribed offset, and unseen documents "
        "have offset zero. All constituent IDs remain in the exported prediction table. X flags document-pair endpoints "
        "within each training fold; a merged unit is removed if any constituent is flagged. It falls back to G if no "
        "remaining unit survives, without dropping the held-out unit. N retains equal-valued documents separately and "
        "has a different evaluation set; it cannot enter the adoption decision.", "",
        "### 2.1 Cross-fitted performance and adoption", "",
        markdown_table(pooled[["method", "n_predictions", "n_entries", "accuracy", "brier", "log_loss", "ece"]]), "",
        "G reproduces 820 predictions, 289 entries, ECE 0.0767537203890331 and Brier 0.5218617813580732. "
        "Full target-specific results are in `figure_data/SI_activity_models.csv`.", "",
        markdown_table(rules[["method", "rule", "G_value", "candidate_value", "required_boundary", "comparison", "passed"]]), "",
        f"**Decision: {adoption['decision']}** No tolerance was added to an adoption inequality. "
        "O improves Brier substantially but fails overall, CA II and Thrombin ECE limits. X improves hERG high-confidence "
        "agreement by at least 0.05 but fails overall ECE. These are negative adoption results, not missing models. "
        "O also increases pooled log loss from G's 0.938343 to 1.051252. "
        "The software default remains Gaussian, with no src/ or tests/ edits in T2b.", "",
        markdown_table(scores[scores.method.isin(["G", "O", "X"]) & ~scores.target.eq("all")][
            ["target_label", "method", "high_confidence_n", "high_confidence_mean_probability", "high_confidence_accuracy"]]), "",
        "Calibration figures use the model probability of B0's displayed class, not the maximum class probability. "
        "Five equal-width bins report x=mean probability, y=observed displayed-class agreement, and n; n<8 bins "
        "are retained but marked ineligible for plotting. B0 contributes its accuracy at x=1. These predictions "
        "evaluate a new evidence unit; mapped class support below instead concerns the consensus.", "",
        "### 2.2 Mapped annotations and sensitivity", "",
        "Sources: `activity/mapped_entry_confidence.csv`, `mapped_entry_confidence_unmerged.csv`, "
        "`threshold_sensitivity.csv`, `mapped_confidence_summary.csv`.", "",
        markdown_table(summary_support[["target_label", "n", "supported_n", "supported_fraction",
                                       "class_changed_by_document_balance", "discordant_document_entries"]]), "",
        markdown_table(class_support[["target_label", "display_class", "n", "supported_n", "supported_fraction"]]), "",
        "There are 3,289 mapped entries under each evidence policy. The complete support table separately reports "
        "all/one-unit/multi-unit populations and every displayed class for cutoffs 0.7, 0.8 and 0.9. "
        "The unmerged sensitivity below uses the same mapped entries but refits tau under the unmerged policy.", "",
        markdown_table(support.query("scope == 'all' and display_class == 'all' and target != 'all'")[
            ["target_label", "evidence_policy", "threshold", "n", "supported_n", "supported_fraction"]]), "",
        "No mapped medium-class entry meets the merged 0.8 support cutoff. This follows from this model and "
        "the one-unit-wide interval relative to the estimated noise; it is not a claim that intermediate activity is absent.", "",
        "### 2.3 Adjacent colour contrasts", "",
        "Source: `activity/neighbor_contrasts.csv` (120 rows: three targets, two layouts, ten seeds and two pair definitions). "
        "The primary definition is occupied cells at Chebyshev distance exactly one (eight-adjacency). "
        "The secondary definition is the symmetric union of tie-inclusive k=10 Euclidean grid-neighbor sets. "
        "Both use unique unordered pairs. Expected true differences among displayed contrasts sum "
        "1-sum_c p_i(c)p_j(c); supported contrasts additionally require that probability >=0.8. "
        "Same-colour expectations sum the same probability over displayed matches. These calculations assume "
        "independent entry consensus distributions; shared-document covariance is not modeled.", "",
        markdown_table(primary[["target_label", "neighbor_pairs", "displayed_contrasts",
            "expected_true_contrasts_among_displayed", "supported_contrasts_p_ge_0.8",
            "expected_hidden_contrasts_same_colour"]]), "",
        "These are optimized seed-42 primary-definition results. `figure_data/Figure_2E_neighbor_contrasts.csv` "
        "also supplies the ten-seed mean, minimum and maximum for every statistic; the 120-row table contains both definitions.", "",
        "### 2.4 Document-offset diagnostics", "",
        "Sources: `activity/document_offsets.csv`, `document_pair_counts.csv`, and each target's flagged-pair table. "
        "A pair is flagged with at least three shared entries, absolute median B-minus-A >=1 pChEMBL and same-sign fraction >=0.8. "
        "A/B are lexicographically ordered; zero differences do not count as same-sign; SD is the sample SD.", "",
        markdown_table(counts), "", markdown_table(pairs), "",
        "The hERG pair CHEMBL5725903 / CHEMBL5727624 has 30 shared entries and median offset 3.47 pChEMBL. "
        "The diagnostic flags shared mapped entries without deleting source records or relabeling displayed medians.", "",
        "## 3. Group-core rules", "",
        "Sources: `cores/group_core_comparison.csv`, `core_size_tradeoff.csv`, `core_target_sensitivity.csv`, "
        "`core_target_sensitivity_summary.csv`, `working_copy_reproduction.csv`, `fmcs_recomputation.csv` and `core_checks.json`.", "",
        "Butina groups are computed from mapped structures and do not depend on layout: 125 CA II, 152 Thrombin and 203 hERG "
        "nonsingleton groups (480 total). The long comparison contains 1,440 group/rule rows. Most-frequent scaffold, "
        "coverage target 0.8, and FMCS (threshold 0.8, completeRingsOnly) each report explicit coverage, size and ring counts.", "",
        markdown_table(tradeoff), "",
        "No core remains categorical: two groups under each scaffold rule have no core, with missing coverage. "
        "Means and medians exclude No core; strict-majority and target-met counts retain all 480 groups in their denominator. "
        "Coverage improves from 268/480 to 470/480 strictly-majority groups; 467/480 reach the coverage target. "
        "The cost is more generic cores: benzene cores increase from 10 to 30 and single-ring cores from 18 to 48.", "",
        f"All 480 cached FMCS results were checked against the current members; a fixed-seed sample of "
        f"{core_checks['fmcs_recomputed_groups']} groups was recomputed and matched. The full cache contains "
        f"{core_checks['fmcs_cached_timeout_groups']} timed-out results under the 10-second limit; those are budget-limited "
        "FMCS results, not certified maximum cores. The sample contains no timeout. SMARTS is authoritative for matching; "
        "core_smiles is a valid fragment from a matching member for drawing and cannot replace the query semantics.", "",
        "Coverage-target sensitivity (0.6, 0.7, 0.8, 0.9, 1.0):", "", markdown_table(core_sensitivity), "",
        "## 4. Cases and Figure 1 consistency", "",
        "Sources: `cases/CHEMBL20_evidence_units.csv`, `CHEMBL20_activity_summary.json`, `core_rule_examples.csv`; "
        "`cores/G001_G002_members.csv`; `figure1_cli_run/consistency_check.json` and `entry_comparison.csv`.", "",
        f"CHEMBL20 has {ace['n_records']} retained records and {ace['n_units']} units; median pChEMBL "
        f"{ace['activity_pchembl']:.6f}, document-balanced value {ace['document_balanced_pchembl']:.9f}, consensus scale "
        f"{ace['consensus_scale']:.9f}. Consensus probabilities (inactive/medium/active) are "
        f"{ace['p_consensus_inactive']:.12g} / {ace['p_consensus_medium']:.9f} / {ace['p_consensus_active']:.9f}. "
        f"New-document probabilities are {ace['p_new_document_inactive']:.9f} / {ace['p_new_document_medium']:.9f} / "
        f"{ace['p_new_document_active']:.9f}; status is {ace['class_status']}. Its seed-42 refined cell is "
        f"({ace['refined_grid_col']}, {ace['refined_grid_row']}). Innocenti 2008 (CHEMBL1141531) contributes 25/54 records "
        "to a single evidence unit. All unit IDs, constituent documents, record counts, merge flags and values are exported.", "",
        "G002 remains benzothiazole (`c1ccc2scnc2c1`) under both scaffold rules: core coverage 27/28, exact scaffold "
        "17/28, exception CHEMBL85636. G001 changes from a 31-heavy-atom core covering 5/42 to a 12-heavy-atom "
        "core covering 42/42. The five ranked alternative cases prioritize >=9-heavy-atom cores and then coverage gain, "
        "with all uncovered member IDs in `cases/core_rule_examples.csv`.", "",
        "The raw CA II CLI was run with `--projection umap --random-state 42 --map-color-modes activity class_support neighbor_recall`. "
        "UMAP is explicit because the bare CLI projection default is PCA. Default neighbor optimization remains enabled. "
        f"All {figure_check['figure1_cli_consistency']['n_entries']} entries match the T2a optimized seed-42 cells and "
        "the present class-support annotation exactly (maximum absolute support difference 0). "
        "T3 should use this run's exported activity, class_support and neighbor_recall maps.", "",
        "## Abstract headline numbers", "",
        "- 3,289 mapped entries; 820 merged-unit predictions from 289 repeated entries.",
        "- G remains the default: pooled Brier 0.521862 and ECE 0.076754; O and X fail preregistered ECE conditions.",
        "- Supported displayed classes: 506/889 CA II (56.9%), 528/1200 Thrombin (44.0%), 532/1200 hERG (44.3%); medium 0%.",
        "- Optimized seed-42 eight-adjacency displayed/supported contrasts: CA II 972/405, Thrombin 2013/752, hERG 1301/349.",
        "- Strict-majority group coverage: 268/480 to 470/480; 467/480 reach target 0.8; benzene cores 10 to 30.",
        "- CHEMBL20: 54 records, 27 units, balanced 7.044074, P(active consensus)=0.613210, uncertain.",
        "- T2a remains the source for ten-seed recall/stability, LAS/FLAS/DPQ and the 9,581-entry full hERG run; see section 1.", "",
        "## Item-by-item working-copy reconciliation", "",
        f"`activity/headline_comparison.csv` covers all {len(comparison)} leaves of "
        "`project1/JCIM_review_innovation_20261024/results/headline_numbers.json`, with old/current values, differences, "
        "status, reasons and source paths. Status counts:", "",
        markdown_table(comparison.groupby("status").size().rename("items").reset_index()), "",
        "Activity values and core counts reproduce the prior population; tiny differences arise from the reference's "
        "rounded tables. Grid means now use ten seeds and refinement seed equal to map seed. The secondary contrast "
        "comparison retains k=10; refined layouts change from refinement seed 0 to 42. Primary eight-adjacency figures "
        "have different denominators and must not be compared directly to old k=10 counts. Historical FMCS/core runtimes "
        "are explicitly not remeasured because cache reuse was authorized. The known full-hERG CHEMBL2035183 "
        "standardization limitation and sandbox-restricted native memory measurement remain as stated in section 1.", ""]
    append_section(root / "RESULTS.md", RESULT_MARKER, "\n".join(lines))
    if (root / "cores/core_rule_timing.csv").exists():
        update_rework_outputs(root)
        return refresh_rework_manifests(root)
    return refresh_manifests(root)


def refresh_manifests(root=OUTPUT):
    """Append current provenance while retaining earlier collection sections."""
    root = Path(root)
    for folder in ("activity", "cores"):
        path = root / folder / "manifest.json"
        if path.exists():
            old = load_json(path)
            inputs = list(old.get("inputs_sha256", {}))
            parameters = old.get("parameters", {})
        else:
            inputs = [FROZEN / target / "map_molecules.csv" for target in TARGETS]
            inputs += [ROOT.parent / "project1/JCIM_review_innovation_20261024/results/core_selection/group_core_comparison.csv"]
            parameters = {"scope": "T2b cores", "coverage_targets": [.6, .7, .8, .9, 1.], "fmcs_sample_seed": 20261024,
                          "fmcs_sample_size": 20, "fmcs_threshold": .8, "fmcs_complete_rings_only": True,
                          "fmcs_timeout_seconds": 10}
        write_new_manifest(path, inputs,
            [p for p in (root / folder).rglob("*") if p.is_file() and p != path], parameters)
    previous = load_json(root / "manifest.json")
    first_hash = previous.get("t2b", {}).get("previous_collection_sha256", sha256(root / "manifest.json"))
    inputs = {str(REFERENCE): sha256(REFERENCE)}
    task_root = ROOT.parent / "project1/JCIM_reliability_revision_20261024/codex_tasks"
    for path in (ROOT.parent / "docs/task.md", ROOT.parent / "AGENTS.md",
                 task_root / "00_方案与验收标准.md", task_root / "T2/T2_brief.md",
                 task_root / "T2/T2a_ACCEPTANCE.md", task_root / "T2/T2a_REWORK_REPORT.md"):
        inputs[str(path)] = sha256(path)
    for manifest_path in root.rglob("*manifest.json"):
        if manifest_path == root / "manifest.json":
            continue
        manifest = load_json(manifest_path)
        for path, digest in manifest.get("inputs_sha256", {}).items():
            p = Path(path)
            if not p.is_absolute():
                p = ROOT.parent / path if path.startswith(("chemgridmap/", "project1/")) else ROOT / path
            if p.is_file() and p.resolve() != (root / "manifest.json").resolve():
                inputs[str(p.resolve())] = sha256(p)
    sources = sorted(p for directory in (ROOT / "src", ROOT / "paper", ROOT / "tests")
                     for p in directory.rglob("*") if p.is_file() and "output" not in p.parts
                     and "__pycache__" not in p.parts and p.suffix in {".py", ".js", ".cjs", ".html"})
    diff = subprocess.check_output(["git", "diff", "release-v0.4.0"], cwd=ROOT)
    reports = ROOT.parent / "project1/JCIM_reliability_revision_20261024/codex_tasks/T2"
    previous["t2b"] = {"collected_at_utc": datetime.now(timezone.utc).isoformat(),
        "previous_collection_sha256": first_hash, "scope": "T2b final collection; T1 and T2a sections retained",
        "inputs_sha256": inputs,
        "parameters": {
            "activity": load_json(root / "activity/model_fits.json")["parameters"],
            "cores": load_json(root / "cores/manifest.json")["parameters"],
            "contrasts": {"seeds": list(SEEDS), "layouts": ["minimum_displacement", "neighbor_refined"],
                          "definitions": ["eight_connected", "tie_aware_k10_union"], "k": 10,
                          "support_threshold": .8, "probabilities": "consensus", "unordered_pairs": True},
            "figure_data": load_json(root / "figure_data/t2b_manifest.json")["parameters"],
            "t2a_historical_parameters": {name: load_json(root / name / "manifest.json").get("parameters", {})
                                          for name in ("layouts", "neighborhoods", "stability", "baselines", "scale")}},
        "environment": {**environment_versions(), "python_executable": sys.executable,
                        "chemgridmap.__version__": chemgridmap.__version__, "chemgridmap.__file__": chemgridmap.__file__},
        "code": {"git_diff_release_v0.4.0_sha256": hashlib.sha256(diff).hexdigest(),
                 "source_sha256": {str(p.resolve()): sha256(p) for p in sources}},
        "reports_sha256": {str(p): sha256(p) for p in (reports / "T2b_REPORT.md", reports / "T2_REPORT.md") if p.exists()},
        "historical_manifests": "Frozen T1/T2a run manifests retain original execution provenance. Current output hashes are below."}
    previous["outputs_sha256"] = {str(p.relative_to(root)): sha256(p) for p in sorted(root.rglob("*"))
                                  if p.is_file() and p != root / "manifest.json" and "__pycache__" not in p.parts
                                  and p.name != ".DS_Store"}
    dump_json(root / "manifest.json", previous)
    return previous["t2b"]


def replace_rework_block(text, name, content, before):
    """Replace one bounded report block, preserving all surrounding text."""
    start, end = f"<!-- T2b rework {name} -->", f"<!-- /T2b rework {name} -->"
    block = start + "\n\n" + content.rstrip() + "\n\n" + end + "\n\n"
    if start in text:
        left, rest = text.split(start, 1)
        _, right = rest.split(end, 1)
        return left + block + right.lstrip("\n")
    return text.replace(before, block + before, 1)


def append_timing_columns(root):
    """Append measured columns without changing the original CSV field strings."""
    summary = pd.read_csv(root / "cores/core_rule_timing_summary.csv", float_precision="round_trip")
    indexed = summary.set_index(["target", "rule"])
    stats = ["total_seconds", "median_seconds", "p90_seconds", "max_seconds", "timeout_groups"]
    path = root / "figure_data/SI_core_size_tradeoff.csv"
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        fields, rows = list(reader.fieldnames), list(reader)
    added = [prefix + stat for prefix in ("timing_", "historical_") for stat in stats]
    fields = [field for field in fields if field not in added] + added
    for row in rows:
        current = indexed.loc[(row["target"], row["rule"])]
        historical_rule = {"coverage": "historical_coverage", "fmcs": "historical_fmcs"}.get(row["rule"])
        historical = indexed.loc[(row["target"], historical_rule)] if historical_rule else None
        for stat in stats:
            row["timing_" + stat] = int(current[stat]) if stat == "timeout_groups" else float(current[stat])
            row["historical_" + stat] = ("" if historical is None else
                int(historical[stat]) if stat == "timeout_groups" else float(historical[stat]))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return summary


def core_timing_sections(root, summary):
    audit = pd.read_csv(root / "cores/fmcs_full_recomputation.csv")
    metadata = load_json(root / "cores/core_timing_metadata.json")
    timeout = audit.cached_timed_out | audit.recomputed_timed_out
    checks = ["same_smarts", "same_coverage_count", "same_coverage", "same_timeout_status",
              "same_member_matches", "same_heavy_atoms", "same_ring_count"]
    differences = pd.DataFrame([
        {"comparison": name, "different_groups": int((~audit[name]).sum()),
         "timeout_in_either_run": int((~audit[name] & timeout).sum()),
         "no_timeout_in_either_run": int((~audit[name] & ~timeout).sum())} for name in checks])
    platform_name = metadata["platform"]
    cpu = metadata["cpu_model"]
    table = summary[["target_label", "rule", "groups", "total_seconds", "median_seconds",
                     "p90_seconds", "max_seconds", "timeout_groups"]]
    lines = ["### 3.1 Sequential timing and full FMCS recomputation", "",
        "Sources: `cores/core_rule_timing.csv` (480 rows), `fmcs_full_recomputation.csv` (480 rows), "
        "`core_rule_timing_summary.csv`, `core_timing_metadata.json` and `core_timing_checks.json`. "
        "The timing columns are also appended to `figure_data/SI_core_size_tradeoff.csv`.", "",
        f"Machine: platform `{platform_name}`; CPU `{cpu}`. All 480 groups were processed sequentially in one process, "
        "with one timed call per rule per group. File reading, SMILES parsing, structural grouping and the separate "
        "cache-comparison audit are excluded. Most-frequent selection includes Murcko scaffold construction, "
        "frequency selection with the software tie-break and coverage matching. Coverage-target selection uses "
        "`chemgridmap.cores.select_core(min_coverage=0.8)`, including candidate generation and matching. FMCS uses "
        "`mcs_core(min_coverage=0.8, timeout=10)`, with ringMatchesRingOnly and completeRingsOnly enabled; "
        "its timer includes FindMCS, query parsing and coverage matching. The 10 s budget applies to the search; "
        "reported call times may exceed it. These are single-run wall times, not repeated benchmark estimates.", "",
        "All times below are seconds; p90 uses linear interpolation. `historical_coverage` and `historical_fmcs` "
        "are the read-only working-copy cache's `aware_seconds` and `mcs_seconds`, with the same rule-call boundaries. "
        "They were measured in a separate historical run whose machine/load context was not recorded in that CSV; "
        "they are contextual comparisons, not paired speedup experiments. Most-frequent historical timing is unavailable.", "",
        markdown_table(table), "",
        f"The historical cache has {int(audit.cached_timed_out.sum())} FMCS timeouts; the current full recomputation has "
        f"{int(audit.recomputed_timed_out.sum())}. {int(audit.passed.sum())}/480 groups match every audited field. "
        "Differences are counted separately below; a timeout in either run can return a different best-so-far query. "
        "Strict SMARTS equality, coverage, member identities and timeout status are distinct checks.", "",
        markdown_table(differences), "",
        "The earlier 20-group audit and all original core-size/coverage tables remain frozen evidence from T2b. "
        "The full-recomputation file records the current queries and coverage separately; added timing columns do not "
        "silently substitute the new FMCS solutions into the accepted structural tables. Timeout results do not certify "
        "a maximum common core. SMARTS remains authoritative for matching."]
    overall = summary[summary.target.eq("overall")].set_index("rule")
    headline = ("- Sequential three-rule timing over 480 groups (most-frequent / coverage target 0.8 / FMCS): "
        f"{overall.loc['most_frequent', 'total_seconds']:.6f} / {overall.loc['coverage', 'total_seconds']:.6f} / "
        f"{overall.loc['fmcs', 'total_seconds']:.6f} s; FMCS was recomputed for all 480 groups with "
        f"{int(audit.recomputed_timed_out.sum())} current timeouts (historical cache: {int(audit.cached_timed_out.sum())}).")
    return "\n".join(lines), headline


def update_rework_outputs(output_root=OUTPUT):
    """Refresh only the R1-R3 report and the authorized appended timing columns."""
    root = Path(output_root)
    summary = append_timing_columns(root)
    activity = rework_activity_sections(root)
    timing_text, timing_headline = core_timing_sections(root, summary)
    path = root / "RESULTS.md"
    text = path.read_text()
    text = replace_rework_block(text, "performance", activity["performance"],
                               "| method | rule | G_value |")
    text = replace_rework_block(text, "tau", activity["tau"], "### 2.2 Mapped annotations and sensitivity")
    text = text.replace("No mapped medium-class entry meets the merged 0.8 support cutoff. This follows from this model and "
        "the one-unit-wide interval relative to the estimated noise; it is not a claim that intermediate activity is absent.\n\n", "")
    text = replace_rework_block(text, "mapped_medium", activity["mapped_medium"], "### 2.3 Adjacent colour contrasts")
    text = replace_rework_block(text, "core timing", timing_text, "## 4. Cases and Figure 1 consistency")
    text = replace_rework_block(text, "headlines", "\n".join([*activity["headlines"], timing_headline]),
                               "## Item-by-item working-copy reconciliation")
    text = text.replace("Historical FMCS/core runtimes are explicitly not remeasured because cache reuse was authorized.",
        "The original historical reconciliation retains its cache-reuse status; section 3.1 now reports the separately measured full R2 timing and FMCS recomputation.")
    path.write_text(text)
    definitions = ["## T2b rework: appended core timing columns", "",
        "`SI_core_size_tradeoff.csv` retains every original row and field. The original structural columns describe "
        "the accepted cached cores; current recomputed FMCS queries/coverage are in `../cores/fmcs_full_recomputation.csv`.", "",
        "Timing summaries join by target and rule from `../cores/core_rule_timing_summary.csv`; `overall` pools all "
        "480 group calls. Platform, CPU and exact timing boundaries are in `../cores/core_timing_metadata.json`. "
        "One process executes all groups sequentially. Parsing/grouping, file I/O and the separate audit are excluded. "
        "FMCS uses threshold 0.8, completeRingsOnly, ringMatchesRingOnly and a 10 s search budget. "
        "Reported wall time includes post-search matching and can exceed 10 s.", "",
        "| Added column | Definition / unit |", "| --- | --- |"]
    meanings = {"total_seconds": "Sum of group call wall times / s", "median_seconds": "Median group wall time / s",
                "p90_seconds": "90th percentile of group wall time, linear interpolation / s",
                "max_seconds": "Maximum group wall time / s", "timeout_groups": "Number of timed-out groups / count"}
    for prefix, source in (("timing_", "Current sequential run"), ("historical_", "Read-only working-copy cache")):
        definitions += [f"| {prefix}{name} | {source}: {description} |" for name, description in meanings.items()]
    definitions += ["", "Historical coverage uses `aware_seconds`; historical FMCS uses `mcs_seconds`/`mcs_timed_out`. "
        "Historical most-frequent values are missing (not recorded), not zero. Only FMCS has a timeout budget; "
        "current scaffold-rule timeout counts are zero. Historical and current runs share rule-call boundaries but "
        "differ in run context, so timing ratios are descriptive only."]
    append_section(root / "figure_data/README.md", "<!-- T2b rework timing columns -->", "\n".join(definitions))


def refresh_rework_manifests(root=OUTPUT):
    """Record current rework provenance without relabeling frozen historical runs."""
    root = Path(root)
    inputs = [root / "activity/model_fits.json", root / "activity/lodo_scores.csv",
              root / "activity/mapped_entry_confidence.csv", root / "cores/core_size_tradeoff.csv",
              ROOT.parent / "project1/JCIM_review_innovation_20261024/results/core_selection/group_core_comparison.csv",
              ROOT.parent / "project1/JCIM_review_innovation_20261024/software/chemgridmap/src/chemgridmap/cores.py",
              ROOT.parent / "project1/JCIM_review_innovation_20261024/analysis/evaluate_core_selection.py",
              ROOT / "validation_data/chembl37_chembl204/download_manifest.json",
              ROOT / "validation_data/chembl37_chembl204/target_metadata.json",
              *[FROZEN / target / "map_molecules.csv" for target in TARGETS]]
    metadata = load_json(root / "cores/core_timing_metadata.json")
    parameters = {"scope": "T2b R1-R3 rework; frozen numeric activity/T2a outputs retained",
                  "core_timing": metadata, "support_threshold": .8, "medium_interval": [6., 7.]}
    for relative in ("cores/manifest.json", "figure_data/t2b_manifest.json"):
        path = root / relative
        previous = load_json(path)
        historical = previous.get("t2b_before_rework", previous)
        outputs = [p for p in path.parent.iterdir() if p.is_file() and p != path
                   and p.name != ".DS_Store"]
        current = write_new_manifest(path, inputs, outputs, parameters)
        current["t2b_before_rework"] = historical
        dump_json(path, current)
    manifest = load_json(root / "manifest.json")
    prior_digest = manifest.get("t2b_rework", {}).get("previous_collection_sha256", sha256(root / "manifest.json"))
    reports = ROOT.parent / "project1/JCIM_reliability_revision_20261024/codex_tasks/T2"
    task_inputs = [ROOT.parent / "docs/task.md", ROOT.parent / "AGENTS.md", reports / "T2b_ACCEPTANCE.md"]
    sources = sorted(p for directory in (ROOT / "src", ROOT / "paper", ROOT / "tests")
                     for p in directory.rglob("*") if p.is_file() and "output" not in p.parts
                     and "__pycache__" not in p.parts and p.suffix in {".py", ".js", ".cjs", ".html"})
    diff = subprocess.check_output(["git", "diff", "release-v0.4.0"], cwd=ROOT)
    manifest["t2b_rework"] = {"collected_at_utc": datetime.now(timezone.utc).isoformat(),
        "previous_collection_sha256": prior_digest, "parameters": parameters,
        "inputs_sha256": {str(p.resolve()): sha256(p) for p in [*inputs, *task_inputs]},
        "environment": {**environment_versions(), "python_executable": sys.executable},
        "code": {"git_diff_release_v0.4.0_sha256": hashlib.sha256(diff).hexdigest(),
                 "source_sha256": {str(p.resolve()): sha256(p) for p in sources}},
        "reports_sha256": {str(p): sha256(p) for p in (reports / "T2b_REWORK_REPORT.md", reports / "T2_REPORT.md") if p.exists()},
        "historical_provenance": "Existing T1/T2a/T2b root sections and frozen subdirectory manifests are retained; this section records the current rework."}
    manifest["outputs_sha256"] = {str(p.relative_to(root)): sha256(p) for p in sorted(root.rglob("*"))
        if p.is_file() and p != root / "manifest.json" and "__pycache__" not in p.parts and p.name != ".DS_Store"}
    dump_json(root / "manifest.json", manifest)
    return manifest["t2b_rework"]


def rework_activity_sections(root):
    root = Path(root)
    activity = root / "activity"
    fits = json.loads((activity / "model_fits.json").read_text())["full_population"]
    with (activity / "lodo_scores.csv").open() as handle:
        scores = {(r["target"], r["method"]): r for r in csv.DictReader(handle)}
    with (activity / "mapped_entry_confidence.csv").open() as handle:
        mapped = list(csv.DictReader(handle))
    order = [("chembl205", "CA II"), ("chembl204", "Thrombin"), ("chembl240", "hERG")]
    normal = NormalDist()
    limit = 0.5 / normal.inv_cdf(0.9)
    def f(value):
        return f"{float(value):.6f}"
    def change(target, metric):
        return f"{f(scores[target, 'G'][metric])} to {f(scores[target, 'O'][metric])}"
    n = int(scores["all", "G"]["n_predictions"])
    herg_n = int(scores["chembl240", "G"]["n_predictions"])
    share = f"{herg_n/n:.1%}"
    rows = ["| target | method | accuracy | brier | ece | n_predictions |",
            "| --- | --- | --- | --- | --- | --- |"]
    for target, label in order:
        for method in ("B0", "G", "O", "X"):
            row = scores[target, method]
            rows.append(f"| {label} | {method} | {f(row['accuracy'])} | {f(row['brier'])} | {f(row['ece'])} | {row['n_predictions']} |")
    performance = "\n\n".join([
        "Target-specific B0/G/O/X results from `activity/lodo_scores.csv`:",
        "Target-label correction (2026-10-24): CHEMBL204 is human thrombin (prothrombin; CHEMBL204). "
        "Early project files mislabeled it as the dopamine D3 receptor. Labels have been corrected to Thrombin; "
        "all numeric results remain unchanged.",
        "\n".join(rows),
        f"hERG contributes {herg_n}/{n} predictions ({share}), compared with {scores['chembl205', 'G']['n_predictions']} for CA II and {scores['chembl204', 'G']['n_predictions']} for Thrombin, so pooled performance is dominated by hERG. Relative to G, O improves hERG accuracy ({change('chembl240', 'accuracy')}), Brier ({change('chembl240', 'brier')}), and ECE ({change('chembl240', 'ece')}). For CA II, accuracy decreases ({change('chembl205', 'accuracy')}), while Brier ({change('chembl205', 'brier')}) and ECE ({change('chembl205', 'ece')}) worsen. Thrombin is mixed: accuracy increases ({change('chembl204', 'accuracy')}), but Brier ({change('chembl204', 'brier')}) and ECE ({change('chembl204', 'ece')}) worsen. Thus the data do not support a claim that all three Thrombin metrics deteriorate.",
        "A possible explanation for the poorer CA II/Thrombin probability scores is that O treats fitted document offsets as fixed: its predictive variance uses the corrected pooled SD but includes neither offset-estimation uncertainty nor additional variance for unseen-document offsets. This may lead to overconfident probabilities; it is an implementation-based explanation, not an experimentally established cause. The offset-adjusted call is in `paper/reliability_activity.py:168-169`, its Gaussian probability calculation in `paper/reliability_activity.py:119-124`, and the variance formula in `src/chemgridmap/evidence.py:184-192`.",
    ])
    rows = ["| target | policy | tau | n_entries | n_units |", "| --- | --- | --- | --- | --- |"]
    maxima, needed_units = [], []
    for target, label in order:
        for policy in ("default", "unmerged"):
            fit = fits[target][policy]
            rows.append(f"| {label} | {policy} | {f(fit['tau'])} | {fit['n_entries']} | {fit['n_units']} |")
        tau = fits[target]["default"]["tau"]
        maxima.append(2 * normal.cdf(0.5 / tau) - 1)
        needed_units.append(math.ceil((tau / limit)**2))
    tau = "\n\n".join([
        "Full-population between-unit SD estimates from `activity/model_fits.json`, `full_population.<target>.default/unmerged`:",
        "\n".join(rows),
        "Here tau (τ) is in pChEMBL log units; `default` merges equal-valued document units and `unmerged` retains them separately. The fits use all retained-record entries eligible to estimate within-entry between-unit variance, rather than only mapped entries. `n_entries` counts entries with at least two units and `n_units` counts their units; these are fit denominators, not all retained or mapped entries.",
        f"Under default G, the Gaussian consensus distribution used for class support has center μ equal to the document-balanced mean and SD s = τ/√m for m evidence units. A one-unit entry therefore has consensus SD τ. This must be distinguished from prediction of a new evidence unit, whose SD is τ√(1 + 1/m), hence √2τ for m = 1. For the medium interval 6 < x < 7, the maximum consensus probability over μ occurs at μ = 6.5 and is `2Φ(0.5/τ) − 1`: CA II {f(maxima[0])}, Thrombin {f(maxima[1])}, and hERG {f(maxima[2])}. Consequently, no one-unit entry can obtain medium support at the 0.8 threshold under these fitted models; this is a mathematical consequence, not a computation error.",
        f"For a multi-unit entry, the exact support condition is `Φ((7 − μ)/s) − Φ((6 − μ)/s) ≥ 0.8`. Such support is attainable for some μ if and only if `2Φ(0.5/s) − 1 ≥ 0.8`, equivalently `s ≤ 0.5/Φ⁻¹(0.9) = {f(limit)}` pChEMBL log units. This scale bound is necessary for any center; at μ = 6.5 it is also sufficient. At the bound, only that center attains 0.8. A smaller scale still requires μ to satisfy the exact probability condition. With s = τ/√m, the default fits require at least {needed_units[0]}, {needed_units[1]}, and {needed_units[2]} evidence units for CA II, Thrombin, and hERG, respectively, even at the optimal center. The single-unit impossibility therefore does not extend to all possible multi-unit entries.",
    ])
    rows = ["| target | displayed_medium_n | one_unit_n | multi_unit_n | minimum_consensus_scale | maximum_P_medium | supported_n |",
            "| --- | --- | --- | --- | --- | --- | --- |"]
    multi_counts = []
    for target, label in order:
        group = [r for r in mapped if r["target"] == target and r["activity_class"] == "medium"]
        one = sum(int(r["n_evidence_units"]) == 1 for r in group)
        multi = len(group) - one
        multi_counts.append(str(multi))
        minimum = min(float(r["consensus_scale"]) for r in group)
        maximum = max(float(r["p_consensus_medium"]) for r in group)
        supported = sum(float(r["p_consensus_medium"]) >= 0.8 for r in group)
        rows.append(f"| {label} | {len(group)} | {one} | {multi} | {f(minimum)} | {f(maximum)} | {supported} |")
    mapped_medium = "\n\n".join([
        "No mapped medium-class entry meets the merged 0.8 support cutoff. The actual medium-class entries in `activity/mapped_entry_confidence.csv` satisfy the following checks:",
        "\n".join(rows),
        f"All displayed-medium entries, including the {'/'.join(multi_counts)} multi-unit entries for CA II/Thrombin/hERG, have consensus scales above {f(limit)}. Their observed maximum medium probabilities are below 0.8, verifying zero support from the exported probabilities and from the necessary scale bound. This result concerns these fitted models and these mapped entries; it does not imply that intermediate activity is absent or that sufficiently precise multi-unit consensus estimates could never support the medium class.",
    ])
    default_values = ', '.join(f"{label} {f(fits[target]['default']['tau'])}" for target, label in order)
    unmerged_values = ', '.join(f(fits[target]['unmerged']['tau']) for target, _ in order)
    headlines = [
        f"- O improves accuracy, Brier and ECE on hERG, which supplies {herg_n}/{n} predictions ({share}); all three worsen on CA II. On Thrombin, accuracy improves but Brier and ECE worsen.",
        f"- Full-population default τ: {default_values} (unmerged: {unmerged_values}). One-unit maximum P(medium): {', '.join(f(v) for v in maxima)}; medium support at 0.8 requires consensus SD ≤{f(limit)} and a suitable center.",
    ]
    return {"performance": performance, "tau": tau, "mapped_medium": mapped_medium, "headlines": headlines}
