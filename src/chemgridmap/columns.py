"""Generate output-column and public-API reference documentation."""
from __future__ import annotations

import argparse
import inspect
from pathlib import Path
from typing import Dict

from . import __all__
import chemgridmap


def _entry(files, description, unit="text"):
    return {"files": tuple(files), "description": description, "unit_or_values": unit}


MAP_FILES = ("grid_coordinates.csv", "molecule_level.csv")
RECORD_FILES = ("retained_records.csv", "excluded_records.csv")


# Exact names and class-expanded patterns form the authoritative output schema.
COLUMNS: Dict[str, dict] = {
    "molecule_identity_key": _entry(MAP_FILES + RECORD_FILES + ("structure_members.csv", "identity_conflicts.csv"), "Stable entry key used for map rows and evidence joins."),
    "molecule_identity_source": _entry(MAP_FILES + RECORD_FILES, "How the entry identifier was obtained.", "parent ID or canonical SMILES fallback"),
    "canonical_smiles": _entry(MAP_FILES + RECORD_FILES + ("structure_members.csv",), "Canonical isomeric SMILES.", "SMILES"),
    "smiles": _entry(MAP_FILES + ("structure_members.csv",), "SMILES column consumed by the map workflow.", "SMILES"),
    "representation_smiles": _entry(MAP_FILES + RECORD_FILES + ("structure_members.csv",), "Structure used for fingerprint representation after the selected policy.", "SMILES"),
    "record_canonical_smiles": _entry(MAP_FILES + RECORD_FILES, "Unique canonical record structures joined with semicolons.", "SMILES list"),
    "structure_policy": _entry(MAP_FILES + RECORD_FILES, "Selected parent/as-recorded normalization policy."),
    "structure_changed": _entry(MAP_FILES + RECORD_FILES, "Whether parent normalization changed the record structure.", "0 or 1"),
    "parent_extraction_excluded": _entry(MAP_FILES + RECORD_FILES, "Whether parent extraction removed atoms.", "0 or 1"),
    "n_record_structures": _entry(MAP_FILES, "Number of distinct record structures for the entry.", "count"),
    "activity_pchembl": _entry(MAP_FILES, "Displayed median pChEMBL value.", "−log10 mol/L"),
    "activity_class": _entry(MAP_FILES, "Activity class from the displayed pChEMBL and thresholds.", "inactive, medium or active"),
    "n_records": _entry(MAP_FILES, "Retained source records for the entry.", "count"),
    "pchembl_min": _entry(MAP_FILES, "Minimum retained pChEMBL value.", "−log10 mol/L"),
    "pchembl_max": _entry(MAP_FILES, "Maximum retained pChEMBL value.", "−log10 mol/L"),
    "pchembl_median": _entry(MAP_FILES, "Median retained pChEMBL value.", "−log10 mol/L"),
    "pchembl_std": _entry(MAP_FILES, "Sample SD of retained pChEMBL values.", "−log10 mol/L"),
    "pchembl_range": _entry(MAP_FILES, "Maximum minus minimum pChEMBL.", "log units"),
    "has_inactive_measurement": _entry(MAP_FILES, "Whether any retained record is inactive.", "0 or 1"),
    "has_active_measurement": _entry(MAP_FILES, "Whether any retained record is active.", "0 or 1"),
    "crosses_clean_boundary": _entry(MAP_FILES, "Whether records cross an activity threshold.", "0 or 1"),
    "is_large_variation": _entry(MAP_FILES, "Whether record variation exceeds the configured range.", "0 or 1"),
    "is_conflicted": _entry(MAP_FILES, "Wide-record-range flag using all configured annotation rules.", "0 or 1"),
    "has_class_boundary_crossing": _entry(MAP_FILES, "Whether values occur on both sides of a class boundary.", "0 or 1"),
    "n_measurement_classes": _entry(MAP_FILES, "Number of observed activity classes.", "count"),
    "repeat_status": _entry(MAP_FILES, "Summary category for repeated activity records."),
    "wide_record_range": _entry(MAP_FILES, "Display alias for the wide-record-range flag.", "0 or 1"),
    "audit_lower_threshold": _entry(MAP_FILES, "Lower pChEMBL threshold saved with the row.", "−log10 mol/L"),
    "audit_upper_threshold": _entry(MAP_FILES, "Upper pChEMBL threshold saved with the row.", "−log10 mol/L"),
    "audit_conflict_range_threshold": _entry(MAP_FILES, "Wide-record-range threshold saved with the row.", "log units"),
    "source_rows": _entry(MAP_FILES, "0-based input rows contributing to the entry.", "semicolon-separated integers"),
    "molecule_id": _entry(MAP_FILES + ("structure_members.csv",), "Source molecule identifier(s)."),
    "parent_molecule_id": _entry(MAP_FILES + ("structure_members.csv",), "Source parent molecule identifier(s)."),
    "target_id": _entry(MAP_FILES, "Target ChEMBL identifier(s)."),
    "target_name": _entry(MAP_FILES, "Target display name(s)."),
    "assay_id": _entry(MAP_FILES, "Assay ChEMBL identifier(s)."),
    "document_id": _entry(MAP_FILES, "Source document identifier(s)."),
    "activity_id": _entry(MAP_FILES + RECORD_FILES, "Activity record identifier(s)."),
    "assay_type": _entry(MAP_FILES + RECORD_FILES, "Assay type(s) among retained records."),
    "n_assays": _entry(MAP_FILES, "Distinct retained assays.", "count"),
    "n_documents": _entry(MAP_FILES, "Distinct source documents.", "count"),
    "n_evidence_units": _entry(MAP_FILES, "Document-balanced units after optional identical-value merging.", "count"),
    "top_unit_documents": _entry(MAP_FILES, "Documents represented by the unit with most records."),
    "top_unit_record_fraction": _entry(MAP_FILES, "Share of records in the largest evidence unit.", "fraction 0–1"),
    "document_balanced_pchembl": _entry(MAP_FILES, "Mean of document-level evidence units.", "−log10 mol/L"),
    "document_balanced_class": _entry(MAP_FILES, "Activity class of the document-balanced estimate."),
    "class_changed_by_document_balance": _entry(MAP_FILES, "Whether balancing changes the displayed class.", "0 or 1"),
    "consensus_scale": _entry(MAP_FILES, "Posterior scale of the entry consensus.", "log units"),
    "p_consensus_{class}": _entry(MAP_FILES, "Posterior probability that the consensus is in each activity class.", "fraction 0–1"),
    "p_new_document_{class}": _entry(MAP_FILES, "Predictive probability for a new document in each class.", "fraction 0–1"),
    "class_confidence": _entry(MAP_FILES, "Posterior consensus probability for the displayed class.", "fraction 0–1"),
    "class_reproducibility": _entry(MAP_FILES, "Predictive probability for the displayed class in a new document.", "fraction 0–1"),
    "discordant_document_pairs": _entry(MAP_FILES, "Flagged document pairs containing both sides of the entry.", "count"),
    "discordant_documents": _entry(MAP_FILES, "Documents in flagged pairs that share this entry."),
    "offset_pair_touched": _entry(MAP_FILES, "Whether this entry occurs in any document of a flagged offset pair.", "boolean"),
    "status_reason": _entry(MAP_FILES, "Support/display decision: supported, low_support or offset_pair."),
    "class_status": _entry(MAP_FILES, "Displayed support status.", "supported or uncertain"),
    "source_row": _entry(("grid_coordinates.csv",) + RECORD_FILES, "0-based source row carried through curation.", "integer"),
    "projection_x_raw": _entry(MAP_FILES, "Unscaled first projection coordinate.", "projection units"),
    "projection_y_raw": _entry(MAP_FILES, "Unscaled second projection coordinate.", "projection units"),
    "projection_x": _entry(MAP_FILES, "Normalized first projection coordinate.", "fraction 0–1"),
    "projection_y": _entry(MAP_FILES, "Normalized second projection coordinate.", "fraction 0–1"),
    "grid_x": _entry(MAP_FILES + ("structure_members.csv",), "Continuous horizontal grid coordinate.", "grid units"),
    "grid_y": _entry(MAP_FILES + ("structure_members.csv",), "Continuous vertical grid coordinate.", "grid units"),
    "grid_col": _entry(MAP_FILES + ("structure_members.csv",), "Zero-based horizontal grid column.", "integer"),
    "grid_row": _entry(MAP_FILES + ("structure_members.csv",), "Zero-based vertical grid row.", "integer"),
    "neighbor_recall": _entry(MAP_FILES, "Tie-inclusive representation-neighbor recall on the grid.", "fraction 0–1"),
    "neighbor_shared_count": _entry(MAP_FILES, "Neighbors present in both representation and grid sets.", "count"),
    "representation_neighbor_count": _entry(MAP_FILES, "Representation neighbors within the kth-distance shell.", "count"),
    "grid_neighbor_count": _entry(MAP_FILES, "Grid neighbors within the kth-distance shell.", "count"),
    "neighbor_recall_fixed_k": _entry(MAP_FILES, "Expected recall with random tie-breaking at fixed k.", "fraction 0–1"),
    "neighbor_kth_similarity": _entry(MAP_FILES, "Tanimoto similarity of the kth Morgan fingerprint neighbor.", "similarity 0–1"),
    "neighbor_kth_distance": _entry(MAP_FILES, "Euclidean distance to the kth non-fingerprint representation neighbor.", "representation distance"),
    "comparison_id": _entry(MAP_FILES + ("structure_members.csv",), "Identifier used in structure comparisons."),
    "structure_group": _entry(MAP_FILES + ("structure_members.csv",), "Structure-group identifier."),
    "structure_group_color": _entry(MAP_FILES + ("structure_members.csv",), "Group layer color.", "hex RGB"),
    "scaffold_smiles": _entry(MAP_FILES + ("structure_members.csv",), "Bemis–Murcko scaffold SMILES.", "SMILES"),
    "name": _entry(("structure_members.csv",), "Supplied molecule name when present."),
    "core_member": _entry(("structure_members.csv",), "Whether the member is covered by the group core.", "boolean"),
    "core_rule": _entry(("structure_members.csv",), "Core-selection rule.", "coverage, most_frequent or fmcs"),
    "core_smiles": _entry(("structure_members.csv",), "Selected group core SMILES.", "SMILES"),
    "core_smarts": _entry(("structure_members.csv",), "FMCS core SMARTS when available.", "SMARTS"),
    "core_coverage_target": _entry(("structure_members.csv",), "Requested core-coverage target.", "fraction 0–1"),
    "map_name": _entry(("metrics.csv",), "Map/output prefix."),
    "representation_type": _entry(("metrics.csv",), "Representation used for projection.", "morgan, descriptors or user_coordinates"),
    "projection_method": _entry(("metrics.csv",), "Two-dimensional projection method.", "pca, tsne, umap or coordinates"),
    "n_samples": _entry(("metrics.csv",), "Mapped molecules.", "count"),
    "k": _entry(("metrics.csv",), "Neighbor count used for reliability.", "integer"),
    "trustworthiness_k": _entry(("metrics.csv",), "Neighborhood size used by trustworthiness.", "integer"),
    "trustworthiness_sample_size": _entry(("metrics.csv",), "Maximum deterministic trustworthiness sample.", "count"),
    "trustworthiness_scope": _entry(("metrics.csv",), "Whether trustworthiness was exhaustive or sampled."),
    "annotation_neighborhood_policy": _entry(("metrics.csv",), "Policy for kth-distance ties."),
    "projection_to_grid_trustworthiness": _entry(("metrics.csv",), "Projection-to-grid trustworthiness.", "fraction 0–1"),
    "projection_grid_knn_overlap": _entry(("metrics.csv",), "Mean kNN overlap between projection and grid.", "fraction 0–1"),
    "projection_grid_tie_aware_recall": _entry(("metrics.csv",), "Tie-aware projection-to-grid recall.", "fraction 0–1"),
    "projection_grid_tie_aware_precision": _entry(("metrics.csv",), "Tie-aware projection-to-grid precision.", "fraction 0–1"),
    "projection_tie_neighborhood_size": _entry(("metrics.csv",), "Average tie-inclusive projection neighborhood size.", "count"),
    "grid_tie_neighborhood_size": _entry(("metrics.csv",), "Average tie-inclusive grid neighborhood size.", "count"),
    "mean_grid_displacement": _entry(("metrics.csv",), "Mean projection-to-grid displacement.", "grid units"),
    "mean_squared_grid_displacement": _entry(("metrics.csv",), "Mean squared projection-to-grid displacement.", "squared grid units"),
    "p95_grid_displacement": _entry(("metrics.csv",), "95th percentile projection-to-grid displacement.", "grid units"),
    "representation_knn_value_abs_diff": _entry(("metrics.csv",), "Mean absolute kNN activity difference in representation space.", "pChEMBL"),
    "projection_knn_value_abs_diff": _entry(("metrics.csv",), "Mean absolute kNN activity difference in projection space.", "pChEMBL"),
    "grid_knn_value_abs_diff": _entry(("metrics.csv",), "Mean absolute kNN activity difference on the grid.", "pChEMBL"),
    "projection_label_purity": _entry(("metrics.csv",), "Projection kNN label purity.", "fraction 0–1"),
    "grid_label_purity": _entry(("metrics.csv",), "Grid kNN label purity.", "fraction 0–1"),
    "projection_distance_metric": _entry(("metrics.csv",), "Metric used by the 2D projection."),
    "projection_fallback": _entry(("metrics.csv",), "Reason PCA replaced UMAP, or blank."),
    "representation_distance_metric": _entry(("metrics.csv",), "Metric used for representation neighbors.", "jaccard or euclidean"),
    "grid_assignment_method": _entry(("metrics.csv",), "Assignment solver selected for this map."),
    "grid_candidate_cells": _entry(("metrics.csv",), "Candidate lattice cells.", "count"),
    "grid_side": _entry(("metrics.csv",), "Grid side length.", "cells"),
    "grid_sizing_strategy": _entry(("metrics.csv",), "How grid size was selected."),
    "grid_target_occupancy": _entry(("metrics.csv",), "Requested cell occupancy.", "fraction 0–1"),
    "grid_actual_occupancy": _entry(("metrics.csv",), "Occupied fraction of candidate cells.", "fraction 0–1"),
    "grid_padding_cells": _entry(("metrics.csv",), "Explicit padding used for grid sizing, if any.", "cells"),
    "coordinate_scaling": _entry(("metrics.csv",), "Projection normalization mode."),
    "grid_refinement": _entry(("metrics.csv",), "Grid refinement mode."),
    "representation_grid_tie_aware_recall": _entry(("metrics.csv",), "Mean tie-aware neighbor recall.", "fraction 0–1"),
    "representation_grid_fixed_k_recall": _entry(("metrics.csv",), "Mean fixed-k neighbor recall.", "fraction 0–1"),
    "assignment_solver": _entry(("metrics.csv",), "Assignment algorithm implementation."),
    "assignment_cost_dtype": _entry(("metrics.csv",), "Numeric dtype of assignment costs."),
    "assignment_candidate_edges": _entry(("metrics.csv",), "Candidate edges supplied to the sparse solver.", "count"),
    "assignment_fallback_edges_added": _entry(("metrics.csv",), "Fallback assignment edges added.", "count"),
    "assignment_fallback_edges_used": _entry(("metrics.csv",), "Fallback assignment edges selected.", "count"),
    "document_a": _entry(("document_offsets.csv",), "First document in the ordered pair."),
    "document_b": _entry(("document_offsets.csv",), "Second document in the ordered pair."),
    "shared_entries": _entry(("document_offsets.csv",), "Entries shared by both documents.", "count"),
    "median_offset": _entry(("document_offsets.csv",), "Median B-minus-A log-unit offset.", "log units"),
    "mean_offset": _entry(("document_offsets.csv",), "Mean B-minus-A log-unit offset.", "log units"),
    "std_offset": _entry(("document_offsets.csv",), "Sample SD of paired offsets.", "log units"),
    "sign_agreement": _entry(("document_offsets.csv",), "Fraction of differences sharing the median sign.", "fraction 0–1"),
    "systematic_offset": _entry(("document_offsets.csv",), "Whether the pair meets offset criteria.", "boolean"),
    "touched_entry_count": _entry(("document_offsets.csv",), "Size of the union of entries in the two documents.", "count"),
    "total_evaluated_pairs": _entry(("document_offsets.csv",), "Pairs meeting the minimum-shared threshold.", "count"),
    "total_observed_pairs": _entry(("document_offsets.csv",), "Pairs sharing at least one entry.", "count"),
    "total_possible_pairs": _entry(("document_offsets.csv",), "All document pairs.", "count"),
    "exclusion_reason": _entry(RECORD_FILES, "Semicolon-separated curation exclusion reasons."),
    "original_smiles": _entry(RECORD_FILES, "Original SMILES field.", "SMILES"),
    "activity_pchembl_raw": _entry(RECORD_FILES, "Source pChEMBL value used for aggregation.", "−log10 mol/L"),
    "action_type": _entry(RECORD_FILES, "ChEMBL activity action type."),
    "activity_comment": _entry(RECORD_FILES, "ChEMBL activity comment."),
    "activity_properties": _entry(RECORD_FILES, "ChEMBL activity properties field."),
    "assay_chembl_id": _entry(RECORD_FILES, "ChEMBL assay identifier."),
    "assay_description": _entry(RECORD_FILES, "Assay description."),
    "assay_variant_accession": _entry(RECORD_FILES, "Assay variant accession."),
    "assay_variant_mutation": _entry(RECORD_FILES, "Assay variant mutation."),
    "bao_endpoint": _entry(RECORD_FILES, "BioAssay Ontology endpoint identifier."),
    "bao_format": _entry(RECORD_FILES, "BioAssay Ontology format identifier."),
    "bao_label": _entry(RECORD_FILES, "BioAssay Ontology format label."),
    "data_validity_comment": _entry(RECORD_FILES, "ChEMBL data-validity comment."),
    "data_validity_description": _entry(RECORD_FILES, "ChEMBL data-validity description."),
    "document_chembl_id": _entry(RECORD_FILES, "Source document identifier."),
    "document_journal": _entry(RECORD_FILES, "Source journal."),
    "document_year": _entry(RECORD_FILES, "Source publication year.", "year"),
    "modality": _entry(RECORD_FILES, "ChEMBL molecule modality."),
    "molecule_chembl_id": _entry(RECORD_FILES, "Source molecule identifier."),
    "molecule_pref_name": _entry(RECORD_FILES, "Preferred molecule name."),
    "parent_molecule_chembl_id": _entry(RECORD_FILES, "Source parent identifier."),
    "pchembl_value": _entry(RECORD_FILES, "ChEMBL pChEMBL value.", "−log10 mol/L"),
    "potential_duplicate": _entry(RECORD_FILES, "ChEMBL duplicate-citation flag.", "boolean"),
    "qudt_units": _entry(RECORD_FILES, "QUDT unit identifier."),
    "record_id": _entry(RECORD_FILES, "Source record identifier."),
    "relation": _entry(RECORD_FILES, "Unstandardized relation."),
    "src_id": _entry(RECORD_FILES, "ChEMBL source-data identifier."),
    "standard_flag": _entry(RECORD_FILES, "ChEMBL standard-value flag.", "boolean"),
    "standard_relation": _entry(RECORD_FILES, "Standardized relation."),
    "standard_text_value": _entry(RECORD_FILES, "Standardized text value."),
    "standard_type": _entry(RECORD_FILES, "Standardized activity type."),
    "standard_units": _entry(RECORD_FILES, "Standardized activity units."),
    "standard_upper_value": _entry(RECORD_FILES, "Standardized upper value."),
    "standard_value": _entry(RECORD_FILES, "Standardized activity value."),
    "target_chembl_id": _entry(RECORD_FILES, "Target ChEMBL identifier."),
    "target_organism": _entry(RECORD_FILES, "Target organism."),
    "target_pref_name": _entry(RECORD_FILES, "Preferred target name."),
    "target_tax_id": _entry(RECORD_FILES, "Target taxonomy identifier."),
    "text_value": _entry(RECORD_FILES, "Source text value."),
    "toid": _entry(RECORD_FILES, "ChEMBL target or document object identifier."),
    "type": _entry(RECORD_FILES, "Activity value type."),
    "units": _entry(RECORD_FILES, "Unstandardized activity units."),
    "uo_units": _entry(RECORD_FILES, "Unit Ontology identifier."),
    "upper_value": _entry(RECORD_FILES, "Unstandardized upper value."),
    "value": _entry(RECORD_FILES, "Unstandardized activity value."),
    "ligand_efficiency_bei": _entry(RECORD_FILES, "Binding efficiency index."),
    "ligand_efficiency_le": _entry(RECORD_FILES, "Ligand efficiency."),
    "ligand_efficiency_lle": _entry(RECORD_FILES, "Lipophilic ligand efficiency."),
    "ligand_efficiency_sei": _entry(RECORD_FILES, "Surface efficiency index."),
    "ligand_efficiency": _entry(RECORD_FILES, "Combined ligand efficiency value."),
    "action_type_action_type": _entry(RECORD_FILES, "Raw action-type field from the export."),
    "action_type_description": _entry(RECORD_FILES, "Action-type description."),
    "action_type_parent_type": _entry(RECORD_FILES, "Parent action type."),
    "structures": _entry(("identity_conflicts.csv",), "Normalized structures for the conflicting parent.", "semicolon-separated SMILES"),
    "record_count": _entry(("identity_conflicts.csv",), "Records associated with one conflicting identity.", "count"),
}


def describe(columns, files):
    """Return one dictionary entry per emitted column, expanding class patterns."""
    rows = []
    for column in columns:
        candidates = [column]
        for activity_class in ("inactive", "medium", "active"):
            if column in {f"p_consensus_{activity_class}", f"p_new_document_{activity_class}"}:
                candidates.append(column.replace(activity_class, "{class}"))
        key = next((name for name in candidates if name in COLUMNS), None)
        if key is None:
            raise KeyError(f"No column description for {files}:{column}")
        item = COLUMNS[key]
        if files not in item["files"]:
            raise KeyError(f"{key} is not registered for {files}")
        rows.append({
            "file": files,
            "column": column,
            **{name: item[name] for name in ("description", "unit_or_values")},
        })
    return rows


def output_columns_markdown(frames=None) -> str:
    """Render the schema; when frames are supplied, validate them first."""
    rows = []
    for files, entry in COLUMNS.items():
        if files in {"structures", "record_count"}:
            listed_files = ("identity_conflicts.csv",)
        else:
            listed_files = entry["files"]
        for filename in listed_files:
            rows.append((filename, files, entry["description"], entry["unit_or_values"]))
    if frames:
        for filename, frame in frames.items():
            describe(frame.columns, filename)
    lines = ["# Output columns", ""]
    for filename in sorted(set(row[0] for row in rows)):
        lines.extend([f"## `{filename}`", "", "| Column | Meaning | Unit or values |", "|---|---|---|"])
        for file_name, column, description, unit in rows:
            if file_name == filename:
                lines.append(f"| `{column}` | {description} | {unit} |")
        lines.append("")
    return "\n".join(lines)


def api_markdown() -> str:
    """Render signatures and the first docstring paragraph for every public name."""
    lines = ["# Public API", ""]
    for name in __all__:
        value = getattr(chemgridmap, name)
        try:
            signature = str(inspect.signature(value))
        except (TypeError, ValueError):
            signature = "()"
        summary = (inspect.getdoc(value) or "").split("\n\n", 1)[0]
        lines.extend([f"## `{name}`", "", "```python", f"{name}{signature}", "```", "", summary, ""])
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", type=Path, help="Write OUTPUT_COLUMNS.md.")
    parser.add_argument("--write-api", type=Path, help="Write API.md.")
    args = parser.parse_args(argv)
    if args.write:
        args.write.parent.mkdir(parents=True, exist_ok=True)
        args.write.write_text(output_columns_markdown(), encoding="utf-8", newline="\n")
    if args.write_api:
        args.write_api.parent.mkdir(parents=True, exist_ok=True)
        args.write_api.write_text(api_markdown(), encoding="utf-8", newline="\n")
    if not args.write and not args.write_api:
        parser.error("Choose --write and/or --write-api.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
