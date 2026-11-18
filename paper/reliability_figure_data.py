"""Export Prism tables and a numeric grid-results section without drawing figures."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from reliability_analysis import OUTPUT, SEEDS, TARGETS, write_manifest, dump_json


BASELINE_METRIC_COLUMNS = [
    "recall_k5", "recall_k10", "recall_k20", "precision_k10", "trustworthiness_k10",
    "dpq16_tanimoto", "dpq2_tanimoto", "dpq16_euclidean", "dpq2_euclidean",
    "dpq16_minus_tanimoto", "dpq2_minus_tanimoto", "dpq16_minus_euclidean",
    "spearman_vs_projection", "mean_displacement_vs_projection_cells",
]


def markdown_table(frame):
    def fmt(value):
        if isinstance(value, (float, np.floating)):
            return "—" if pd.isna(value) else f"{value:.6f}"
        return str(value)
    return "| " + " | ".join(map(str, frame.columns)) + " |\n| " + " | ".join(
        "---" for _ in frame.columns) + " |\n" + "\n".join(
            "| " + " | ".join(fmt(x) for x in row) + " |" for row in frame.itertuples(index=False, name=None))


def write_results(root, metrics, stability, baselines, scale, cases):
    summary = metrics.groupby(["target_label", "layout"], sort=False).agg(
        recall_mean=("recall_k10", "mean"), recall_min=("recall_k10", "min"),
        recall_max=("recall_k10", "max"), trustworthiness_mean=("trustworthiness_k10", "mean"),
        atom_pair_recall_mean=("atom_pair_recall_k10", "mean"),
        displacement_cells_mean=("mean_displacement_cells", "mean"),
        spearman_mean=("spearman_projection_vs_layout", "mean")).reset_index()
    stable = stability.groupby(["target_label", "layout"], sort=False).neighbor_jaccard_k10.agg(
        ["count", "mean", "median", "min", "max"]).reset_index()
    baseline = baselines.groupby(["target_label", "target_occupancy", "method"], sort=False).agg(
        recall_mean=("recall_k10", "mean"), recall_min=("recall_k10", "min"),
        recall_max=("recall_k10", "max"), dpq16_mean=("dpq16_tanimoto", "mean"),
        dpq16_min=("dpq16_tanimoto", "min"), dpq16_max=("dpq16_tanimoto", "max"),
        dpq16_minus_mean=("dpq16_minus_tanimoto", "mean"),
        dpq2_minus_mean=("dpq2_minus_tanimoto", "mean"),
        dpq16_minus_euclidean_mean=("dpq16_minus_euclidean", "mean"),
        trustworthiness_mean=("trustworthiness_k10", "mean"),
        spearman_vs_projection_mean=("spearman_vs_projection", "mean"),
        mean_displacement_vs_projection_cells=("mean_displacement_vs_projection_cells", "mean")).reset_index()
    timing = baselines.loc[baselines.runtime_reportable, ["target_label", "target_occupancy", "method",
        "projection_seconds", "assignment_seconds", "representation_neighbors_seconds", "algorithm_seconds", "total_seconds"]]
    bridge = pd.read_csv(root / "t1_acceptance_opus/headline/grid_refinement.csv")
    # Keep the frozen T1 bridge unchanged and translate its label only in memory.
    bridge["target"] = bridge.target.replace({"D3": "Thrombin"})
    old = bridge.groupby("target").recall_after.mean()
    lines = ["# RESULTS — 1. Grids and scale (T2a)", "",
        "This section covers T2a only; evidence models and overall core-rule comparisons are reserved for T2b. "
        "The three targets contain 889 mapped entries for CA II, 1,200 for thrombin (prothrombin; CHEMBL204), and 1,200 for hERG. Seeds: " +
        ", ".join(map(str, SEEDS)) + ". Ranges across seeds describe layout randomness, not biological replication.", "",
        "## 1.1 Headline values and neighborhood metrics", "",
        "Sources: `neighborhoods/layout_metrics.csv`, `figure_data/Figure_3A_neighbor_recall.csv`. "
        "Recall includes all neighbors tied at the kth distance and is averaged over entries.", ""]
    for target, label in TARGETS.items():
        subset = metrics[metrics.target.eq(target)]
        means = subset.groupby("layout").recall_k10.mean()
        before, after, projection = (means[x] for x in ("minimum_displacement", "neighbor_refined", "projection"))
        recovered = (after - before) / (projection - before)
        lines.append(f"- {label}: mean recall for projection / minimum displacement / neighbor refinement was "
                     f"{projection:.6f} / {before:.6f} / {after:.6f}; refinement increased recall by "
                     f"{after-before:.6f}, recovering {recovered:.1%} of the projection-to-grid recall decrease.")
    lines += ["", "Tradeoffs (means over the same 10 seeds):", ""]
    for target, label in TARGETS.items():
        means = metrics[metrics.target.eq(target)].groupby("layout").mean(numeric_only=True)
        before, after = means.loc["minimum_displacement"], means.loc["neighbor_refined"]
        lines.append(f"- {label}: mean displacement changed from {before.mean_displacement_cells:.6f} to "
                     f"{after.mean_displacement_cells:.6f} cells; Spearman correlation with the projection changed "
                     f"from {before.spearman_projection_vs_layout:.6f} to {after.spearman_projection_vs_layout:.6f}.")
    lines += ["", markdown_table(summary), "", "## 1.2 Stability across seeds", "",
        "Sources: `stability/seed_pair_jaccard.csv`, `figure_data/Figure_3B_seed_pair_jaccard.csv`. "
        "Each target and layout has 45 seed pairs. For each pair, neighbor-set Jaccard is computed per entry and averaged.", "",
        markdown_table(stable), "", "## 1.3 LAS / FLAS baselines", "",
        "Sources: `baselines/baseline_metrics.csv`, `baselines/summary.csv`, `figure_data/SI_las_baseline.csv`. "
        "Each group has 10 seeds, giving 360 rows. Target occupancy 1.0 still uses side length ceil(sqrt(N)), "
        "leaving a few vacancies; actual occupancy is recorded in the CSV files.", "",
        "DPQ16 and DPQ2 use Tanimoto distances unless marked euclidean. DPQ-minus uses the mean tie rule. "
        "Spearman correlation compares all pairwise layout distances with the same-seed UMAP projection. "
        "Displacement is the mean Euclidean distance from each assigned cell to its normalized projection "
        "coordinate multiplied by (grid_side - 1), in cells. All values in the following table are seed means "
        "except the explicitly labeled minima and maxima.", "",
        markdown_table(baseline), "", "Runtime values use only the original sequential, isolated single-process "
        "measurements for seed 42, in seconds; these measurements were retained without rerunning layouts or timing. "
        "total_seconds includes projection, LAP, neighbor preparation, and the algorithm required by each method, "
        "excluding shared fingerprint preparation and metric evaluation. algorithm_seconds is the final algorithm "
        "step (equal to assignment_seconds for plain LAP and not counted twice).", "",
        markdown_table(timing), "", "The following comparisons use means over 10 seeds for each target and occupancy:", ""]
    for (label, occupancy), subset in baseline.groupby(["target_label", "target_occupancy"]):
        best_recall = subset.loc[subset.recall_mean.idxmax()]
        best_dpq = subset.loc[subset.dpq16_mean.idxmax()]
        best_dpq_minus = subset.loc[subset.dpq16_minus_mean.idxmax()]
        best_trustworthiness = subset.loc[subset.trustworthiness_mean.idxmax()]
        best_spearman = subset.loc[subset.spearman_vs_projection_mean.idxmax()]
        least_displacement = subset.loc[subset.mean_displacement_vs_projection_cells.idxmin()]
        group_timing = timing[timing.target_label.eq(label) & timing.target_occupancy.eq(occupancy)]
        runtimes = "; ".join(f"{row.method} {row.total_seconds:.6f}" for row in group_timing.itertuples())
        lines.append(f"- {label}, target occupancy {occupancy:g}: highest recall, {best_recall.method} "
            f"({best_recall.recall_mean:.6f}); highest Tanimoto DPQ16, {best_dpq.method} ({best_dpq.dpq16_mean:.6f}); "
            f"highest Tanimoto DPQ-minus16, {best_dpq_minus.method} ({best_dpq_minus.dpq16_minus_mean:.6f}); "
            f"highest trustworthiness, {best_trustworthiness.method} ({best_trustworthiness.trustworthiness_mean:.6f}); "
            f"highest Spearman correlation with the projection, {best_spearman.method} "
            f"({best_spearman.spearman_vs_projection_mean:.6f}); lowest displacement from the projection, "
            f"{least_displacement.method} ({least_displacement.mean_displacement_vs_projection_cells:.6f} cells). "
            f"Seed 42 total runtime (seconds): {runtimes}.")
    lines += ["", "Method source: Barthel, Hezel, Jung, Schall, Computer Graphics Forum 42(1), 261–276 (2023), "
        "doi:10.1111/cgf.14718, equations 9–13 / algorithms 3–4; author reference implementation: "
        "https://github.com/Visual-Computing/LAS_FLAS. This independent implementation uses floating-point "
        "rectangular LAP (the reference implementation quantizes costs to integers). Vacancies are handled by "
        "U(mask·X)/U(mask), with the global mean for empty windows and reflect boundaries without wrapping.", "",
        "FLAS fixes the map at each radius, precomputes all N×M costs, and extracts submatrices for candidate windows. "
        "This is mathematically equivalent to computing costs per window but uses more memory; its timings do not "
        "directly represent the author reference implementation. Sparse Morgan inner products use CSR and have "
        "been independently tested against direct squared Euclidean distances.", "",
        "Morgan bit vectors are L2-normalized, giving squared Euclidean distance 2−2cos. With differing bit counts, "
        "the global distance ranking need not match Tanimoto. The primary DPQ uses Tanimoto; normalized Euclidean "
        "DPQ is also reported, both at p=16 and p=2. DPQ-minus is reported for Tanimoto at p=16 and p=2 and for "
        "normalized Euclidean distance at p=16.", "",
        "**DPQ tie rules:** DPQ_p (`ties=\"sorted\"`) orders neighbors tied in two-dimensional distance by "
        "high-dimensional distance, giving the optimistic estimate used by the reference implementation. "
        "DPQ-minus_p (`ties=\"mean\"`) replaces high-dimensional distances within each two-dimensional tie group "
        "by their mean, giving the pessimistic estimate. In the 20×20 linear color-gradient fixture, the original "
        "layout gives DPQ16=1 and DPQ-minus16=1; the seed 42 random "
        "layout gives DPQ16=0.4249241333, DPQ-minus16=0.0114411272, "
        "and DPQ2=0.0655529182. The <0.1 random-layout criterion applies to DPQ-minus16; "
        "optimistic DPQ16 is expected to be below 0.6 and above DPQ-minus16. Favoring the smallest "
        "high-dimensional distance within a two-dimensional tie raises early-rank gains, especially for p=16.", "",
        "## 1.4 Full hERG scale", "",
        "Sources: `scale/full_herg_summary.json`, `scale/full_herg_stage_timings.csv`, `figure_data/SI_full_herg_scale.csv`.", "",
        markdown_table(scale), ""]
    large = json.loads((root / "scale/full_herg_summary.json").read_text())
    cli = large["cli"]
    stage_frame = pd.read_csv(root / "scale/full_herg_stage_timings.csv")
    lines += [f"The full CLI processed {large['n_entries']:,} entries with recall {cli['recall_k10']:.6f}, "
              f"total elapsed time {cli['wall_seconds']:.3f} seconds, and peak RSS {cli['peak_rss_mib']:.3f} MiB. "
              f"RSS source: {cli['peak_rss_source']}; time -l memory measurement status: {cli['time_l_memory_status']}.", "",
        markdown_table(stage_frame[["stage", "seconds", "source"]]), "",
        "The run manifest does not separately record grouping/core calculations. map_uninstrumented is the "
        "remainder of total map time after recorded stages. outside_map_workflow also includes import, curation, "
        "class support, and ZIP export; neither remainder is an isolated stage measurement.", "",
        "**Input standardization difference:** The CLI and frozen layout have the same 9,581 identity keys, "
        "but the SMILES for CHEMBL2035183 changes during repeated RDKit standardization. The two parsed "
        "representations have Morgan Tanimoto similarity 0.7323943662 (one entry has a different fingerprint). "
        "The CLI and frozen layout therefore are not paired before/after layouts with identical inputs; each "
        "recall was computed independently from its own structures. `scale/structure_comparison.json` records "
        "the strings, bidirectional substructure comparisons, and three standardization rounds. This software "
        "input-standardization issue remains for follow-up; frozen src/ files were not modified.", "",
        "## 1.5 Cases and positions", "",
        "Sources: `cases/case_summary.json`, `cases/G002_CHEMBL20_positions_seed42.csv`, "
        "`figure_data/Figure_3C_CHEMBL81688_neighbors.csv`, `figure_data/Figure_1ABC_layout_seed42_refined.csv`.", ""]
    for example in cases["CHEMBL81688"]:
        lines.append(f"- CHEMBL81688 {example['layout']}: {example['fingerprint_neighbor_count']} fingerprint "
            f"neighbors, {example['grid_neighbor_count']} grid neighbors, {example['shared_neighbor_count']} shared "
            f"neighbors, recall {example['neighbor_recall']:.6f}; Q=query, S=shared, G=grid only, F=fingerprint only.")
    g = cases["G002_check"]
    lines += [f"- Under the coverage rule, G002 has {g['group_member_count']} members, identical to v2; core "
        f"`{g['core_smiles']}` covers {g['core_coverage_count']}/{g['group_member_count']}, with exact scaffold "
        f"matches in {g['core_exact_scaffold_count']}/{g['group_member_count']}.", "",
        "## 1.6 Link to the working copy and interpretation limits", "",
        "The T1 frozen-layout results for random_state=0 are quoted directly from "
        "`t1_acceptance_opus/headline/grid_refinement.csv` without recomputation: " +
        "; ".join(f"{label} {old[label]:.6f}" for label in TARGETS.values()) + ". "
        "The present analysis uses 10 map seeds per target, with refinement seeds equal to map seeds, "
        "so its means need not exactly match the earlier values.", "",
        "Neighbor refinement optimizes Morgan k=10 recall; increased recall is not independent validation of "
        "structural prediction. Atom-pair fingerprints, scaffolds, trustworthiness, displacement, and Spearman "
        "correlation provide additional views; full values are in the neighborhood tables. Lower global distance "
        "correlation or greater displacement is a method tradeoff. Grid proximity must not be equated directly "
        "with fingerprint similarity.", ""]
    (root / "RESULTS.md").write_text("\n".join(lines), encoding="utf-8")


def run_figure_data(output_root=OUTPUT):
    root = Path(output_root)
    out = root / "figure_data"
    out.mkdir(parents=True, exist_ok=True)
    paths = {"metrics": root / "neighborhoods/layout_metrics.csv",
             "individual": root / "neighborhoods/per_molecule_recall_seed42.csv",
             "stability": root / "stability/seed_pair_jaccard.csv",
             "baselines": root / "baselines/baseline_metrics.csv",
             "scale": root / "scale/full_herg_scale.csv",
             "case_neighbors": root / "cases/CHEMBL81688_neighbors_seed42.csv",
             "annotations": root / "cases/chembl205_annotations_seed42.csv"}
    data = {key: pd.read_csv(path, **({"float_precision": "round_trip"} if key == "baselines" else {}))
            for key, path in paths.items()}
    assert len(data["metrics"]) == 90
    assert len(data["stability"]) == 405
    assert len(data["baselines"]) == 360
    assert len(data["baselines"].loc[data["baselines"].runtime_reportable]) == 36
    for key, columns, expected in (
        ("metrics", ["target", "seed", "layout"], 90),
        ("stability", ["target", "layout", "seed_a", "seed_b"], 405),
        ("baselines", ["target", "seed", "target_occupancy", "method"], 360)):
        assert len(data[key][columns].drop_duplicates()) == expected
    assert len(data["individual"]) == 3289 * 3
    assert data["baselines"].loc[data["baselines"].runtime_reportable, "seed"].eq(42).all()
    assert np.isfinite(data["baselines"][BASELINE_METRIC_COLUMNS].to_numpy(dtype=float)).all()
    recall = data["metrics"][["target", "target_label", "layout", "seed", "n", "recall_k10"]].rename(
        columns={"n": "n_entries", "recall_k10": "mean_morgan_neighbor_recall_k10"})
    positions = data["annotations"].merge(data["individual"].query(
        "target == 'chembl205' and layout == 'neighbor_refined'")[["molecule_identity_key", "recall_k10"]],
        on="molecule_identity_key", validate="one_to_one")
    assert len(positions) == 889 and positions.is_G002.sum() == 28 and positions.is_CHEMBL20.sum() == 1
    tables = {"Figure_3A_neighbor_recall.csv": recall,
              "Figure_3B_seed_pair_jaccard.csv": data["stability"],
              "Figure_3C_CHEMBL81688_neighbors.csv": data["case_neighbors"],
              "Figure_1ABC_layout_seed42_refined.csv": positions,
              "SI_las_baseline.csv": data["baselines"], "SI_full_herg_scale.csv": data["scale"]}
    for name, frame in tables.items():
        frame.to_csv(out / name, index=False)
    descriptions = {
        "target": "Target ID (lowercase ChEMBL)", "target_id": "Target ID (uppercase ChEMBL)",
        "target_label": "Target short name",
        "layout": "Layout: projection / minimum_displacement / neighbor_refined",
        "seed": "Map and refinement random seed", "seed_a": "First map seed in the pair",
        "seed_b": "Second map seed in the pair", "random_state": "Random seed",
        "n": "Number of mapped entries", "n_entries": "Number of mapped entries",
        "mean_morgan_neighbor_recall_k10": "Mean per-entry Morgan recall at k=10, including distance ties",
        "neighbor_jaccard_k10": "Mean per-entry Jaccard of tied k=10 neighbor sets across two layouts",
        "molecule_identity_key": "Stable molecular identity key", "molecule_id": "Molecule display ID",
        "canonical_smiles": "Canonical isomeric SMILES", "query_id": "Query molecule ID",
        "is_query": "Whether the entry is query Q", "grid_col": "Grid column (zero-based, cells)",
        "grid_row": "Grid row (zero-based, cells)", "refined_grid_col": "Refined grid column (zero-based, cells)",
        "refined_grid_row": "Refined grid row (zero-based, cells)",
        "representation_neighbor": "Membership in the query's Morgan neighbor set",
        "grid_neighbor": "Membership in the query's grid neighbor set",
        "shared_neighbor": "Membership in both neighbor sets",
        "neighbor_relation": "query/shared/grid_only/representation_only",
        "figure_marker": "Q=query / S=shared / G=grid only / F=fingerprint only",
        "representation_distance_tanimoto": "Tanimoto distance from the query (1 - similarity)",
        "grid_distance_cells": "Euclidean distance from the query (cells)",
        "structure_group": "Butina group ID under the coverage rule",
        "scaffold_smiles": "Entry's Bemis–Murcko scaffold SMILES", "is_G002": "Membership in CA II group G002",
        "is_CHEMBL20": "Whether the entry is acetazolamide CHEMBL20",
        "recall_k10": "Morgan recall at k=10 with ties; per entry in Figure 1, mean in SI",
        "method": "Baseline: projection_lap/projection_lap_refined/las_random/las_projection/flas_random/random",
        "target_occupancy": "Target occupancy", "actual_occupancy": "Actual occupancy N/side²",
        "grid_side": "Square grid side length (cells)",
        "recall_k5": "Mean Morgan recall at k=5, including distance ties",
        "recall_k20": "Mean Morgan recall at k=20, including distance ties",
        "precision_k10": "Mean per-entry shared-neighbor count / grid-neighbor count",
        "trustworthiness_k10": "Trustworthiness from high-dimensional Jaccard distances to the layout, k=10",
        "dpq16_tanimoto": "Optimistic DPQ, sorted ties, Tanimoto high-dimensional distance, p=16",
        "dpq2_tanimoto": "Optimistic DPQ, sorted ties, Tanimoto high-dimensional distance, p=2",
        "dpq16_euclidean": "Optimistic DPQ, sorted ties, Euclidean distance on L2-normalized Morgan vectors, p=16",
        "dpq2_euclidean": "Optimistic DPQ, sorted ties, Euclidean distance on L2-normalized Morgan vectors, p=2",
        "dpq16_minus_tanimoto": "Pessimistic DPQ-minus, mean ties, Tanimoto high-dimensional distance, p=16",
        "dpq2_minus_tanimoto": "Pessimistic DPQ-minus, mean ties, Tanimoto high-dimensional distance, p=2",
        "dpq16_minus_euclidean": "Pessimistic DPQ-minus, mean ties, Euclidean distance on L2-normalized Morgan vectors, p=16",
        "spearman_vs_projection": "Spearman correlation of pairwise layout distances with the same-seed UMAP projection",
        "mean_displacement_vs_projection_cells": "Mean Euclidean distance from assigned cells to normalized same-seed projection coordinates multiplied by (grid_side - 1), in cells",
        "runtime_reportable": "True only for isolated sequential seed 42 measurements reportable as runtime",
        "projection_recomputed": "Whether the projection was recomputed for the original layout run",
        "projection_seconds": "UMAP time (seconds; zero when not recomputed)",
        "assignment_seconds": "LAP time (seconds)",
        "representation_neighbors_seconds": "Fingerprint-neighbor preparation time for refinement only (seconds)",
        "algorithm_seconds": "Final algorithm step (seconds); equals assignment_seconds for plain LAP",
        "total_seconds": "Total layout construction time (seconds), excluding shared fingerprint preparation and metric evaluation",
        "runtime_scope": "fresh_full_layout: independent timing; cached_projection_layout: nonreportable batch using saved projections when required",
        "analysis": "Scale-measurement workflow", "input_layout": "Raw CSV or frozen layout input",
        "recall_before_k10": "Recall before refinement at k=10 (blank when not separately exported by the full CLI)",
        "recall_after_k10": "Recall after refinement at k=10",
        "wall_seconds": "Elapsed time (seconds), with scope specified by timing_scope",
        "timing_scope": "full_cli: complete workflow; refine_grid_layout_only: refinement only",
        "peak_rss_bytes": "Peak resident memory (bytes); blank when not separately measured for refinement only",
        "peak_rss_mib": "Peak resident memory (MiB; 1 MiB = 1024² bytes)"}
    readme = ["# T2a Prism figure source data", "",
        "All files are tidy CSV in UTF-8; no figures are drawn. Boolean columns use True/False. "
        "Missing values mean unmeasured or not applicable, not zero.",
        "Neighbor sets retain distance ties and exclude self. pChEMBL-related data remain in the original "
        "neighborhoods tables. RESULTS.md reports the values from these files.", ""]
    for name, frame in tables.items():
        missing = set(frame.columns) - set(descriptions)
        if missing:
            raise ValueError(f"Undocumented figure columns: {name}: {missing}")
        readme += [f"## {name} ({len(frame)} rows)", "", "| Column | Meaning |", "| --- | --- |"]
        readme += [f"| {column} | {descriptions[column]} |" for column in frame]
        readme.append("")
    readme += ["Figure 1 currently contains only CA II layouts, group IDs, recall, and case markers. T2b will add "
        "activity classes and class support to the final Figure_1ABC_map_seed42_refined.csv.",
        "For SI runtime reporting, select runtime_reportable=True only. Other seeds contribute to quality means "
        "and ranges. All original timing values are retained when baseline metrics are supplemented. Where the "
        "sandbox blocks sysctl for scale time -l measurements, RSS is obtained using resource.getrusage; the "
        "original error and source are recorded in the scale JSON files.",
        "DPQ uses sorted high-dimensional distances within each two-dimensional tie group (optimistic); "
        "DPQ-minus replaces those distances by the tie-group mean (pessimistic). Projection correspondence "
        "uses layouts/<target>/coordinates_seed_<seed>.csv, scaled by (grid_side - 1).",
        "The full CLI and frozen-layout comparison each use their own structures for recall. Non-idempotent "
        "standardization of CHEMBL2035183 means these are not paired layouts with identical inputs. Stage "
        "timings are in ../scale/full_herg_stage_timings.csv.", ""]
    (out / "README.md").write_text("\n".join(readme), encoding="utf-8")
    cases = json.loads((root / "cases/case_summary.json").read_text())
    write_results(root, data["metrics"], data["stability"], data["baselines"], data["scale"], cases)
    dump_json(out / "checks.json", {"passed": True, "row_counts": {name: len(frame) for name, frame in tables.items()},
        "all_columns_documented": True, "G002_count": 28, "CHEMBL20_count": 1,
        "baseline_reportable_timings": 36, "baseline_metrics_finite": True,
        "baseline_metric_columns": BASELINE_METRIC_COLUMNS, "results_file": str(root / "RESULTS.md")})
    write_manifest(out, list(paths.values()) + [root / "cases/case_summary.json", root / "scale/full_herg_summary.json",
        root / "scale/full_herg_stage_timings.csv", root / "scale/structure_comparison.json",
        root / "t1_acceptance_opus/headline/grid_refinement.csv", root / "RESULTS.md"],
        {"scope": "T2a", "expected_baseline_rows": 360, "expected_stability_rows": 405})
