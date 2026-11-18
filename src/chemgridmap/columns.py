"""Generate output-column and public-API reference documentation."""
from __future__ import annotations

import argparse
import inspect
from pathlib import Path
from typing import Dict

from . import __all__
import chemgridmap


def _entry(files, description, unit="text", by_file=None):
    return {
        "files": tuple(files),
        "description": description,
        "unit_or_values": unit,
        "by_file": by_file or {},
    }


def _description(entry, filename):
    """Use a file-specific meaning when a column has multiple definitions."""
    return entry["by_file"].get(filename, entry["description"])


MAP_FILES = ("grid_coordinates.csv", "molecule_level.csv")
RECORD_FILES = ("retained_records.csv", "excluded_records.csv")


# Exact names and class-expanded patterns form the authoritative output schema.
COLUMNS: Dict[str, dict] = {
    "molecule_identity_key": _entry(MAP_FILES + RECORD_FILES + ("structure_members.csv", "identity_conflicts.csv", "structure_variants.csv"), "Stable entry key used for map rows and evidence joins."),
    "molecule_identity_source": _entry(MAP_FILES + RECORD_FILES, "How the entry identifier was obtained.", "parent ID or canonical SMILES fallback"),
    "canonical_smiles": _entry(MAP_FILES + RECORD_FILES + ("structure_members.csv",), "Canonical isomeric SMILES.", "SMILES"),
    "smiles": _entry(
        MAP_FILES + ("structure_members.csv", "structure_variants.csv"),
        "SMILES column consumed by the map workflow.",
        "SMILES",
        {"structure_variants.csv": "Parent-structure SMILES writing observed for this identity."},
    ),
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
    "class_confidence": _entry(MAP_FILES, "Class support: posterior probability that the consensus (document-balanced mean) lies in the displayed class. The column keeps its 0.5 name for compatibility; it is not the chance that a new document agrees (see class_reproducibility).", "fraction 0–1"),
    "class_reproducibility": _entry(MAP_FILES, "Class reproducibility: predictive probability that one new document reports the displayed class; the quantity checked against later documents in the evaluation.", "fraction 0–1"),
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
    "scaffold_smiles": _entry(
        MAP_FILES + ("structure_members.csv",),
        (
            "Bemis–Murcko scaffold SMILES. Charge-separated groups can leave "
            "representation artifacts after removal of the negative oxygen: "
            "sulfoxides may appear as `[SH+]` and N-oxides as `[nH+]`. These "
            "scaffold strings are not neutral structures. Sulfones retain "
            "`S(=O)(=O)`, so sulfoxides and sulfones yield different scaffolds."
        ),
        "SMILES",
    ),
    "name": _entry(("structure_members.csv",), "Supplied molecule name when present."),
    "core_member": _entry(("structure_members.csv",), "Whether the member is covered by the group core.", "boolean"),
    "core_rule": _entry(("structure_members.csv",), "Core-selection rule.", "coverage, most_frequent or fmcs"),
    "core_smiles": _entry(
        ("structure_members.csv",),
        (
            "Selected group core SMILES. Charge-separated groups can leave "
            "representation artifacts after removal of the negative oxygen: "
            "sulfoxides may appear as `[SH+]` and N-oxides as `[nH+]`. These "
            "core strings are not neutral structures. Sulfones retain "
            "`S(=O)(=O)`, so sulfoxides and sulfones yield different cores."
        ),
        "SMILES",
    ),
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
    "record_count": _entry(
        ("identity_conflicts.csv", "structure_variants.csv"),
        "Records associated with one conflicting identity.",
        "count",
        {"structure_variants.csv": "Total retained records for this identity."},
    ),
    "normalized_smiles": _entry(
        ("structure_variants.csv",),
        "Standardized structure form for this SMILES writing.",
        "SMILES",
    ),
    "structure_record_count": _entry(
        ("structure_variants.csv",),
        "Retained records using this SMILES writing.",
        "count",
    ),
    "standard_inchikey": _entry(
        ("structure_variants.csv",),
        "Standard InChIKey shared by the writing variants for this identity.",
        "InChIKey",
    ),
    "selected": _entry(
        ("structure_variants.csv",),
        "Whether this writing normalizes to the map's selected structure (the normalized form with the most records).",
        "boolean",
    ),
}


# Optional own-data fields are emitted only when their feature is requested.
COLUMNS.update({
    "n_measured": _entry(("metrics.csv",), "Measured reference molecules used for activity/class diagnostics in a query run.", "count"),
    "n_query": _entry(("metrics.csv",), "Unmeasured molecules included in structural/layout diagnostics.", "count"),
    "raw_activity_value": _entry(MAP_FILES + RECORD_FILES, "Original activity value before unit conversion; the input column is also retained.", "input unit"),
    "raw_activity_unit": _entry(MAP_FILES + RECORD_FILES, "Original per-row unit, or the selected common unit, before conversion.", "M, mM, uM, µM, μM, nM, pM"),
    "sigma_doc": _entry(MAP_FILES, "Assumed between-unit SD, only emitted for assumed_sigma.", "log units"),
    "sigma_doc_degrees_of_freedom": _entry(MAP_FILES, "Empty for assumed_sigma because the SD was not estimated.", "empty"),
    "measurement_status": _entry(MAP_FILES + ("structure_members.csv",), "Whether the molecule belongs to measured reference data or the unmeasured query set.", "measured or not_measured"),
    "query_id": _entry(MAP_FILES + ("query_check.csv", "query_neighbors.csv", "structure_members.csv"), "User query identifier or generated Q00001-style identifier; empty for measured rows."),
    "input_smiles": _entry(("query_check.csv",), "Original query SMILES before structure preparation.", "SMILES"),
    "status": _entry(("query_check.csv",), "Query inclusion decision, one row per input row.", "mapped, invalid_smiles, duplicate_query or in_reference"),
    "first_query_id": _entry(("query_check.csv",), "First query identifier with the same canonical structure; empty otherwise."),
    "reference_id": _entry(("query_check.csv",), "Measured molecule with the same canonical structure; empty otherwise."),
    "rank": _entry(("query_neighbors.csv",), "Similarity rank among measured neighbors, including every tie at the fifth neighbor.", "integer starting at 1"),
    "measured_id": _entry(("query_neighbors.csv",), "Measured reference molecule identifier."),
    "tanimoto": _entry(("query_neighbors.csv",), "Morgan fingerprint Tanimoto using the map radius and bit count, including for descriptor maps.", "similarity 0–1"),
    "nearest_measured_id": _entry(MAP_FILES, "Most similar measured molecule in Morgan fingerprint space; no activity is inferred."),
    "nearest_measured_similarity": _entry(MAP_FILES, "Morgan Tanimoto to the nearest measured molecule.", "similarity 0–1"),
    "query_group_center_similarity": _entry(MAP_FILES + ("structure_members.csv",), "Morgan Tanimoto to the assigned reference group center; empty when unassigned.", "similarity 0–1"),
    "contains_group_core": _entry(MAP_FILES + ("structure_members.csv",), "Whether the unmeasured molecule contains its assigned reference group's core; empty when no core/group applies.", "boolean or empty"),
})
for _column in ("source_row", "canonical_smiles"):
    COLUMNS[_column]["files"] += ("query_check.csv",)
for _column in ("activity_class", "activity_pchembl", "class_confidence", "class_status", "n_evidence_units"):
    COLUMNS[_column]["files"] += ("query_neighbors.csv",)
COLUMNS["name"]["files"] += MAP_FILES
COLUMNS["activity_class"]["unit_or_values"] = "inactive, medium, active; not_measured for unmeasured map rows"

# U8 新文件保留既有检查/近邻列契约；参照文件 schema 不变。
for _entry_value in COLUMNS.values():
    _files = _entry_value["files"]
    for _legacy, _current in (("query_check.csv", "new_molecules_check.csv"),
                             ("query_neighbors.csv", "new_molecules_neighbors.csv")):
        if _legacy in _files:
            _entry_value["files"] += (_current,)
    if any(_file in _files for _file in MAP_FILES + ("query_check.csv", "structure_members.csv")):
        _entry_value["files"] += ("new_molecules.csv",)
COLUMNS["tanimoto"]["files"] += ("new_molecules.csv",)
COLUMNS["exclusion_reason"]["files"] += ("new_molecules_excluded_records.csv",)
COLUMNS.update({
    "label": _entry(("new_molecules.csv",), "User label, used only for display and filtering."),
    "is_new": _entry(("new_molecules.csv",), "Identifies the current new-molecule overlay.", "boolean"),
    "nearest_reference_id": _entry(("new_molecules.csv",), "Placement target: highest Morgan Tanimoto, then smallest reference ID."),
    "placement_distance": _entry(("new_molecules.csv",), "Euclidean distance from the target reference cell.", "grid cells"),
    "placement_recall": _entry(("new_molecules.csv",), "Share of the k = 10 most similar reference molecules (ties retained) that are among the 10 reference cells nearest to the placed cell.", "fraction 0–1"),
    "placement_recall_shared": _entry(("new_molecules.csv",), "Number of those similar reference molecules that are also grid neighbours of the placed cell.", "count"),
    "placement_recall_fingerprint_count": _entry(("new_molecules.csv",), "Number of similar reference molecules counted, including ties at the tenth similarity.", "count"),
    "neighbor_evidence": _entry(("new_molecules.csv",), "Evidence from similar tie-inclusive top-five reference neighbors; never a prediction.", "none, insufficient, active, medium, inactive, mixed"),
    "similar_reference_count": _entry(("new_molecules.csv",), "Top-five neighbors including ties meeting group_similarity.", "count"),
    "supported_reference_count": _entry(("new_molecules.csv",), "Supported similar neighbors; zero when reference support was not assessed.", "count"),
    "supported_reference_classes": _entry(("new_molecules.csv",), "Classes of supported similar neighbors; empty when support was not assessed.", "semicolon-separated classes"),
    "evidence_reference_count": _entry(("new_molecules.csv",), "Neighbors contributing to evidence; uses observed classes when support was not assessed.", "count"),
    "evidence_reference_classes": _entry(("new_molecules.csv",), "Classes contributing to neighbor evidence, including observed-class fallback.", "semicolon-separated classes"),
    "reference_support_assessed": _entry(("new_molecules.csv",), "Whether reference class support was assessed.", "boolean"),
    "neighbor_comparison": _entry(("new_molecules.csv",), "Own measured class compared with consistent neighbor evidence only.", "agrees, disagrees, or empty"),
    "user_median": _entry(("new_molecules.csv",), "Own retained-record median, when available.", "pActivity"),
    "reference_median": _entry(("new_molecules.csv",), "Median of the structurally identical reference molecule.", "pActivity"),
    "median_difference": _entry(("new_molecules.csv",), "Own median minus structurally identical reference median.", "log units"),
})


# Optional rule comparisons have their own schema; existing export text stays fixed.
CORE_RULE_FILES = ("core_rules/<stem>_cores_<rule>.csv",
                   "core_rules/group_<group-hash>_core_comparison.csv")
for _column in ("structure_group", "core_rule", "core_smiles", "core_coverage_target"):
    COLUMNS[_column]["files"] += CORE_RULE_FILES
COLUMNS["core_rule"]["by_file_units"] = {
    _file: "coverage, coverage_brics, coverage_recap, most_frequent or fmcs; the comparison CSV adds a generic_murcko_supplementary row" for _file in CORE_RULE_FILES}
COLUMNS.update({
    "size": _entry(CORE_RULE_FILES, "Number of original reference members in the group.", "count"),
    "core_smiles_attached": _entry(CORE_RULE_FILES, "Core depiction with attachment points generated from a matching member.", "SMILES"),
    "core_heavy_atoms": _entry(CORE_RULE_FILES, "Number of non-hydrogen atoms in the selected core, excluding attachment dummies.", "count"),
    "core_coverage_count": _entry(CORE_RULE_FILES, "Reference members matching the core with chirality enabled.", "count"),
    "core_coverage_fraction": _entry(CORE_RULE_FILES, "Matching reference members divided by group size.", "fraction 0–1"),
    "core_target_met": _entry(CORE_RULE_FILES, "Whether the coverage target is met; empty for most_frequent.", "boolean or empty"),
    "candidates": _entry(CORE_RULE_FILES, "Candidates evaluated; FMCS does not enumerate this candidate pool.", "count"),
    "seconds": _entry(CORE_RULE_FILES, "Elapsed core calculation time for this group.", "seconds"),
    "brics_candidates": _entry(CORE_RULE_FILES, "BRICS candidates absent from the original scaffold/ring-fragment pool; only for coverage_brics.", "count"),
    "brics_single_pass_members": _entry(CORE_RULE_FILES, "Members with more than 10 BRICS bonds, decomposed with singlePass=True; only for coverage_brics.", "count"),
    "recap_candidates": _entry(CORE_RULE_FILES, "RECAP candidates absent from the original scaffold/ring-fragment pool; only for coverage_recap.", "count"),
    "recap_skipped_members": _entry(CORE_RULE_FILES, "Members with more than 80 heavy atoms that were not RECAP-decomposed (they still count in coverage); only for coverage_recap.", "count"),
    "fmcs_timed_out": _entry((CORE_RULE_FILES[1],), "FMCS reached its fixed 10-second comparison budget; the returned result may be incomplete.", "boolean"),
})


def _unit(entry, filename):
    return entry.get("by_file_units", {}).get(filename, entry["unit_or_values"])


# Optional manual exports extend only their own schemas, preserving all legacy dictionaries.
from .app import MANUAL_COLUMNS, FEEDBACK_COLUMNS

_MANUAL_DESCRIPTIONS = {
    "value": ("Exact manually entered value; empty for a qualitative class.", "entered unit"),
    "unit": ("Unit of the entered exact value.", "pActivity, M, mM, µM, nM, pM"),
    "measured_value": ("Exact entered value; uploaded multiple-record aggregates use pActivity.", "measured_unit"),
    "measured_unit": ("Unit of the measured value.", "pActivity, M, mM, µM, nM, pM"),
    "measured_class": ("Measured class, from an exact value or qualitative user choice.", "active, medium, inactive"),
    "class_source": ("Origin of the measured class.", "value, user_class"),
    "result_source": ("Origin of the measured result.", "upload, manual, recorded"),
    "batch": ("User-supplied batch or publication used as an evidence unit.", "text"),
    "note": ("User-supplied note.", "text"),
    "entered_at": ("Time the molecule was manually added.", "ISO 8601 with timezone"),
    "recorded_at": ("Last manual result save time; empty for uploaded results.", "ISO 8601 with timezone"),
    "nearest_reference_similarity": ("Morgan Tanimoto similarity to the placement reference.", "0 to 1"),
    "contains_group_core_rule": ("Rule used for the saved new-molecule core match; API inspection uses the current map rule.", "coverage, coverage_brics, most_frequent, fmcs"),
}
for _filename, _names in (("manual_entries.csv", MANUAL_COLUMNS),
                           ("neighbor_evidence_feedback.csv", FEEDBACK_COLUMNS)):
    for _name in _names:
        if _name not in COLUMNS:
            _meaning, _units = _MANUAL_DESCRIPTIONS[_name]
            COLUMNS[_name] = _entry((), _meaning, _units)
        COLUMNS[_name]["files"] += (_filename,)
        if _name == "neighbor_comparison":
            COLUMNS[_name]["by_file_units"] = {**COLUMNS[_name].get("by_file_units", {}),
                                               _filename: "agrees, disagrees, no_comparison"}
for _name in _MANUAL_DESCRIPTIONS:
    if _name not in {"value", "unit", "entered_at"}:
        if _name not in COLUMNS:
            _meaning, _units = _MANUAL_DESCRIPTIONS[_name]
            COLUMNS[_name] = _entry((), _meaning, _units)
        COLUMNS[_name]["files"] += ("new_molecules.csv",)


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
            "description": _description(item, files),
            "unit_or_values": _unit(item, files),
        })
    return rows


def output_columns_markdown(frames=None) -> str:
    """Render the schema; when frames are supplied, validate them first."""
    rows = []
    for files, entry in COLUMNS.items():
        if files == "structures":
            listed_files = ("identity_conflicts.csv",)
        else:
            listed_files = entry["files"]
        for filename in listed_files:
            rows.append((filename, files, _description(entry, filename), _unit(entry, filename)))
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
    lines.extend([
        "## Optional core-rule outputs", "",
        "When building with `coverage_brics`, the existing structure-member tables also use that `core_rule` value.",
        "Their legacy column-dictionary text remains unchanged for byte-identical default exports.",
        "Map-page switches write only `core_rules/<stem>_cores_<rule>.csv` and the matching `.json`.",
        "The JSON records the rule, build coverage target, RDKit version, each group's `structure_sha256`,",
        "the original run manifest SHA-256, and total elapsed seconds. Structure hashes use the grouping",
        "protocol's sorted canonical-structure serialization. Existing run files, manifest and ZIP are unchanged.",
        "The comparison CSV uses the first 12 hexadecimal digits of the group identifier SHA-256 for `<group-hash>`.",
        "BRICS and RECAP columns are present for their rule exports and in the five-rule comparison table.",
        "The comparison table ends with a `generic_murcko_supplementary` row: the most frequent generic Murcko",
        "scaffold (all atoms carbon, all bonds single) and the members sharing it exactly. It is supplementary",
        "information, not a selectable core rule, and its target columns are empty.",
        "Comparison ZIPs for non-build rules also live under `core_rules/` and include the rule in their names.", "",
        "## Optional own-data outputs and metadata", "",
        "Reference filenames retain the run's `<name>_` prefix and never change on insertion.",
        "New results live in `new_molecules/<input-stem>/`, with numeric suffixes for existing names.",
        "`new_molecules_check.csv` and `new_molecules_neighbors.csv` retain the legacy query column contracts.",
        "`new_molecules.csv` contains all checked molecules; placement fields are empty for excluded rows.",
        "Measured input also writes `new_molecules_excluded_records.csv`: every original input column,",
        "followed by `exclusion_reason`; the header is retained when no records are excluded.",
        "Group overlay summaries add `n_query_members` and `n_new_members`; saved reference groups stay unchanged.",
        "Generic reference records with unit conversion or assumed sigma also export",
        "`retained_records.csv`, `molecule_level.csv` and `curation_report.json`.", "",
        "Unit-conversion reports contain `excluded_reasons` counts and the selected unit/column.",
        "Fits and manifests identify assumed SDs with `method=assumed_sigma`, `assumed_sigma`,",
        "and `evidence_unit_source` (the supplied column or `one unit per compound`).",
        "With no evidence-unit column, `document_offset_detection` is `not applicable`.",
        "The existing automatic `literature_fallback` remains a separate method.", "",
        "The separate new-molecule manifest records the input and reference-manifest hashes, parameters and timing.",
        "Reference layout, metrics, diagnostics and original ZIP contain reference molecules only.",
        "New-molecule neighbors summarize observed reference evidence, never predictions.",
        "Legacy query-named schema sections below describe the retained low-level query functions.", "",
        '## Manual-entry and measured-result feedback outputs',
        '',
        '`manual_entries.csv` has one row per manually added molecule. Its SHA-256 and row count',
        'are recorded in `new_molecules_manifest.json`. `neighbor_evidence_feedback.csv` includes',
        'all currently measured new molecules, including uploaded results. Saving, modifying or',
        'clearing a manual result rewrites this table and the current collection ZIP. An explicit',
        'feedback download can also create it for an unchanged uploaded collection. File-only',
        'insertion retains its original outputs until a manual action or feedback download.',
        '',
        'Class-only results use `class_source=user_class` and leave numeric activity empty.',
        '`no_comparison` covers absent, insufficient or mixed neighbor evidence. These tables',
        'record observations and their agreement with reference evidence; no activity is predicted.',
        'Timestamps use ISO 8601 with a timezone; unknown uploaded recording times remain blank.',
        'Browser refresh retains the current session. After restarting, **Open a previous',
        'result** restores the saved active collection, manual entries and measured results.',
        '',
        'Exact CSV column order:',
        '',
        "- `manual_entries.csv`: " + ", ".join(f"`{name}`" for name in MANUAL_COLUMNS) + ".",
        "- `neighbor_evidence_feedback.csv`: " + ", ".join(f"`{name}`" for name in FEEDBACK_COLUMNS) + ".", "",
        "## Own-data CLI options", "",
        "| Option | Meaning |", "|---|---|",
        "| `--value-unit {p,M,mM,uM,nM,pM}` | Default `p` preserves existing pActivity values; concentrations become −log10 mol/L. |",
        "| `--unit-column COL` | Per-row concentration units, also accepting µM and μM; case sensitive, whitespace trimmed. Mutually exclusive with explicit `--value-unit`. |",
        "| `--assumed-sigma S` | Finite positive SD in log units for record-table support; incompatible with molecule-only input or disabled support. |",
        "| `--query-file PATH` | Insert unmeasured CSV/TSV molecules into a completed Morgan reference grid, up to 1,000 after checking. Measured new molecules use the browser. |",
        "| `--query-smiles-column COL` | Override automatic SMILES-column detection. |",
        "| `--query-id-column COL` | Query IDs, unique within the input; default generated IDs preserve input order. |", "",
    ])
    lines.extend([
        '',
        '## Session cache (browser runs)',
        '',
        'A successful browser build additionally writes `session/reference.pkl.gz`,',
        '`session/state.pkl.gz` and `session/session_manifest.json`. The reference blob',
        'contains the complete built result and source-lookup data; the state blob retains',
        'current core-rule caches, comparisons, selected rule, output paths and active',
        'new-molecule collection. The active collection stores `session/overlay.pkl.gz`',
        'and its own `session_manifest.json` under `new_molecules/<set>/`. Input frames are',
        'reconstructed from the byte-preserved run `input.csv`.',
        '',
        'The JSON sidecar records SHA-256 hashes for snapshot files, `source_tree_digest`,',
        'Python/pandas/numpy/RDKit versions, run-manifest SHA-256 and creation time. Both',
        'snapshots and existing scientific outputs are verified before opening; changed',
        'CSV files, missing files or incompatible versions require regeneration. Stored',
        'paths are relative, so a complete run folder can move to another output root.',
        'Opening is read-only and restores cached rules and the active collection without',
        'recomputing; core/feedback operations subsequently update their session state.',
        '',
        'These compressed caches occupy additional disk space and are excluded from the',
        'existing run manifest and ZIP. They are not new scientific result columns.',
        'Deleting `session/` preserves exports but prevents restoring that session; delete',
        'the entire run folder when the result is no longer needed. The import-page list',
        "only scans this launcher's output root and shows its latest 20 compatible results.",
        '',
        'Browser core computation uses parallel loky workers from 64 groups, with a default of',
        '`min(8, os.cpu_count() - 1)` (at least 1), overridden by `CHEMGRIDMAP_JOBS`.',
        'Library `jobs=1` and CLI `--jobs 1` remain serial by default and do not read that',
        'environment variable. FMCS is always serial',
        'because load changes its timeout solution. No parallel setting changes output',
        'column order, core selection or pre-existing exported files. Browser-only marker',
        'sizes/tooltips do not change the exported `map_with_new_molecules.svg`.',
        '',
    ])
    return "\n".join(lines)


def api_markdown() -> str:
    """Render signatures and complete docstrings for every public name."""
    lines = ["# Public API", ""]
    for name in __all__:
        value = getattr(chemgridmap, name)
        try:
            signature = str(inspect.signature(value))
        except (TypeError, ValueError):
            signature = "()"
        documentation = inspect.getdoc(value) or ""
        lines.extend([
            f"## `{name}`", "", "```python", f"{name}{signature}", "```", "",
            documentation, "",
        ])
    from .cores import brics_candidate_cores, generic_scaffold_summary, recap_candidate_cores, select_core
    from .comparison import summarize_group_core
    from .inspection import build_source_index
    from .app import AppSession
    for value in (brics_candidate_cores, recap_candidate_cores, generic_scaffold_summary, select_core, summarize_group_core, build_source_index,
                  AppSession.source_members, AppSession.set_core_rule, AppSession.group_cores):
        name = value.__module__ + "." + value.__qualname__
        lines.extend([f"## `{name}`", "", "```python", f"{name}{inspect.signature(value)}", "```", "",
                      inspect.getdoc(value) or "", ""])
    lines.extend([
        "## Local browser endpoints", "",
        "All endpoints retain the local session-token requirement; errors return HTTP 400.",
        "`GET /api/source-members?source=ID&source=ID` returns source metadata, distinct molecule counts",
        "and the union of mapped members (`grid_row`, `grid_col`, molecule ID and identity key).",
        "Add `scope=new` to limit the index to the current inserted molecules and their batches.",
        "`GET /api/inspect` adds `map_molecule_count` to retained sources and a `source_union_molecule_count`;",
        "inserted molecules expose `new_molecule_sources`. Possible source-lookup matches are excluded.",
        "`POST /api/core-rule` accepts `{\"rule\": \"coverage_brics\"}` (or `coverage_recap`, `coverage`, `most_frequent`).",
        "It returns the current rule, build rule, group summaries, downloadable files and cache status.",
        "FMCS can be restored only when it was the build rule; global switching never starts new FMCS jobs.",
        "`GET /api/group-cores?group=ID` returns five core summaries with SVGs, timing and the comparison CSV,",
        "plus `generic` (supplementary generic Murcko scaffold, member count and SVG); `GET /api/group` also returns `generic`;",
        "FMCS receives a fixed 10-second budget. Existing group membership and map coordinates are unchanged.", "",
    ])
    lines.extend([
        '## Manual new molecules and measured results',
        '',
        'These local endpoints use the same session token as the existing new-molecule',
        'API. A completed Morgan reference map is required. Validation failures return',
        'HTTP 400 with an `error`; rejected previews return an explicit status.',
        '',
        '| Method and path | Payload or response |',
        '|---|---|',
        '| `POST /api/new-molecules/preview` | `{"smiles":"CCF"}`; returns `status`, `canonical_smiles`, `svg`, and `reference_id` or `query_id` when matched. Status is `ok`, `invalid`, `in_reference` or `duplicate_new`. |',
        '| `POST /api/new-molecules/add` | `smiles`, optional `query_id`, `label`, `note`, and `measurement`; appends to the current collection. |',
        '| `POST /api/new-molecules/result` | `query_id` and `measurement`; records or updates a manual measured result. Uploaded measured results are read-only. |',
        '| `POST /api/new-molecules/result/clear` | `{"query_id":"M001"}`; removes a manual result and restores the unmeasured state. |',
        '| `GET /api/new-molecules/feedback` | Downloads `neighbor_evidence_feedback.csv`; also creates the export on demand for an otherwise unchanged uploaded collection. |',
        '',
        'An exact value uses this shared measurement object:',
        '',
        '```json',
        '{"mode":"value","value":100,"unit":"nM","batch":"R1","note":"Exact measurement"}',
        '```',
        '',
        'A qualitative or censored result uses:',
        '',
        '```json',
        '{"mode":"class","measured_class":"inactive","batch":"R2","note":"> 10 µM"}',
        '```',
        '',
        'Supported units are `pActivity`, `M`, `mM`, `µM`, `nM` and `pM`; values must be',
        'finite and strictly positive. `units.convert_activity_values` performs conversion,',
        "and the reference map's thresholds determine the class. Class-only input accepts",
        '`active`, `medium` or `inactive`, leaves pActivity empty and sets',
        '`class_source=user_class`. Batch and note are optional free text.',
        '',
        'Successful add/result/clear calls return `status=ok`, `query_id`, `row`, `point`',
        'and `job`. Add can instead return the same rejected status as preview. IDs default',
        'to the first unused `M001`, `M002`, etc.; supplied IDs must be unique in the current',
        'new-molecule set. Normalization and group assignment use the file insertion path.',
        'Only the new cell is appended in the browser; existing reference and new cells',
        "remain fixed. Results change only that molecule's measurement and comparison.",
        '',
        '`GET /api/new-molecules` adds `directory` and',
        '`replacement={total,manual,recorded,directory}` for the replacement confirmation.',
        'Its rows expose `result_source` (`upload`, `manual`, `recorded` or empty),',
        '`measured_value`, `measured_unit`, `measured_class`, `class_source`, `batch`,',
        '`note` and `recorded_at`. The report adds `manual_count` and',
        '`feedback_counts={recorded,agrees,disagrees,no_comparison}`. The inspector adds',
        '`measured_result`, including editability and the current neighbor comparison.',
        'The session remains available after browser refresh. After restarting, opening a',
        'saved result restores the active collection, manual additions and recorded results.',
        '',
        "Hand-entered changes update the collection's CSVs, SVG, manifest and ZIP.",
        '`manual_entries.csv` records manual additions and its SHA-256 and row count are',
        'stored in the manifest. `neighbor_evidence_feedback.csv` covers all currently',
        'measured new molecules, including uploaded results. Original uploaded bytes and',
        'reference exports are never rewritten. File-only insertion retains its previous',
        'output contract until a manual action or explicit feedback download is requested.',
        '',
    ])
    lines.extend([
        '',
        '## Saved sessions and core progress',
        '',
        'Browser computation uses separate loky processes for `most_frequent`, `coverage`,',
        '`coverage_brics` and `coverage_recap` when there are at least 64 groups; smaller runs stay serial.',
        'Defaults use',
        '`min(8, os.cpu_count() - 1)`, at least 1; `CHEMGRIDMAP_JOBS` overrides the browser',
        'worker count only. `build_grid_map(..., jobs=1)` and CLI `--jobs 1` remain serial',
        'by default and do not read this environment variable;',
        'explicit `jobs>1` requires a guarded Python main entry point on spawn platforms.',
        'FMCS remains serial even when `jobs>1` because parallel load changes timeout',
        'solutions. Returned group ordering and scientific values are unchanged.',
        '',
        '- `GET /api/core-progress` returns `{running, rule, done, total}`. Poll every second',
        '  while the existing `POST /api/core-rule` request is pending; the map stays usable.',
        '- `GET /api/previous-runs` returns `{runs, unavailable}`. `runs` contains up to 20',
        '  newest compatible runs under the current launcher output root. Each row has',
        '  `run_name`, `created_at`, `filename`, `reference_map`, `molecules`, `core_rule`',
        '  and integer `new_molecules`. `unavailable` counts `missing_snapshot`,',
        '  `incompatible`, and `modified` runs without listing individual old runs.',
        '- `POST /api/open-run` accepts only `{"run_name":"directory-name"}` from that list.',
        '  Absolute paths and traversal are rejected, as is opening during a build. The',
        '  response has the same session envelope as `GET /api/session`, plus',
        '  `restored={new_set,total,manual,recorded}` for the restoration status message.',
        '  The Python method `AppSession.open_run(name)` returns the job snapshot plus',
        '  `restored`. Opening restores original input/settings, cached core rules,',
        '  last-selected rule, active inserted collection, manual additions and measured',
        "  results without recalculating or rewriting the run's existing files. Switching",
        '  to an already computed rule returns `cached: true`. A new build creates a new run.',
        '',
        'Snapshot failures return HTTP 400 with `error`, `code`, `files`, and `details`.',
        'The client translates `session_missing`, `session_version`, `session_modified`,',
        '`session_unsafe`, `session_invalid`, `run_not_listed`, and `map_running`; affected',
        'relative filenames remain visible. Invalid files or incompatible versions require',
        'regeneration. No partial state is applied.',
        '',
        '`session/` stores compressed pickle state and a `session_manifest.json` with',
        'SHA-256 hashes, `source_tree_digest`, Python/pandas/numpy/RDKit versions and the',
        'run-manifest hash/creation time. Restore verifies both snapshot and scientific',
        'output hashes before restricted deserialization. Paths are relative to the run,',
        'so copying the whole folder to another output root preserves reopening support.',
        'Active new-molecule collections have their own `session/`, updated after insertion',
        'or manual changes. These private caches consume additional disk space and never',
        'enter existing export ZIPs or run manifests. Deleting `session/` disables reopening',
        'without removing exports; delete the entire run folder when no longer needed.',
        '',
    ])
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
