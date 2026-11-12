# Output columns

## `core_rules/<stem>_cores_<rule>.csv`

| Column | Meaning | Unit or values |
|---|---|---|
| `structure_group` | Structure-group identifier. | text |
| `core_rule` | Core-selection rule. | coverage, coverage_brics, most_frequent or fmcs |
| `core_smiles` | Selected group core SMILES. Charge-separated groups can leave representation artifacts after removal of the negative oxygen: sulfoxides may appear as `[SH+]` and N-oxides as `[nH+]`. These core strings are not neutral structures. Sulfones retain `S(=O)(=O)`, so sulfoxides and sulfones yield different cores. | SMILES |
| `core_coverage_target` | Requested core-coverage target. | fraction 0–1 |
| `size` | Number of original reference members in the group. | count |
| `core_smiles_attached` | Core depiction with attachment points generated from a matching member. | SMILES |
| `core_heavy_atoms` | Number of non-hydrogen atoms in the selected core, excluding attachment dummies. | count |
| `core_coverage_count` | Reference members matching the core with chirality enabled. | count |
| `core_coverage_fraction` | Matching reference members divided by group size. | fraction 0–1 |
| `core_target_met` | Whether the coverage target is met; empty for most_frequent. | boolean or empty |
| `candidates` | Candidates evaluated; FMCS does not enumerate this candidate pool. | count |
| `seconds` | Elapsed core calculation time for this group. | seconds |
| `brics_candidates` | BRICS candidates absent from the original scaffold/ring-fragment pool; only for coverage_brics. | count |
| `brics_single_pass_members` | Members with more than 10 BRICS bonds, decomposed with singlePass=True; only for coverage_brics. | count |

## `core_rules/group_<group-hash>_core_comparison.csv`

| Column | Meaning | Unit or values |
|---|---|---|
| `structure_group` | Structure-group identifier. | text |
| `core_rule` | Core-selection rule. | coverage, coverage_brics, most_frequent or fmcs |
| `core_smiles` | Selected group core SMILES. Charge-separated groups can leave representation artifacts after removal of the negative oxygen: sulfoxides may appear as `[SH+]` and N-oxides as `[nH+]`. These core strings are not neutral structures. Sulfones retain `S(=O)(=O)`, so sulfoxides and sulfones yield different cores. | SMILES |
| `core_coverage_target` | Requested core-coverage target. | fraction 0–1 |
| `size` | Number of original reference members in the group. | count |
| `core_smiles_attached` | Core depiction with attachment points generated from a matching member. | SMILES |
| `core_heavy_atoms` | Number of non-hydrogen atoms in the selected core, excluding attachment dummies. | count |
| `core_coverage_count` | Reference members matching the core with chirality enabled. | count |
| `core_coverage_fraction` | Matching reference members divided by group size. | fraction 0–1 |
| `core_target_met` | Whether the coverage target is met; empty for most_frequent. | boolean or empty |
| `candidates` | Candidates evaluated; FMCS does not enumerate this candidate pool. | count |
| `seconds` | Elapsed core calculation time for this group. | seconds |
| `brics_candidates` | BRICS candidates absent from the original scaffold/ring-fragment pool; only for coverage_brics. | count |
| `brics_single_pass_members` | Members with more than 10 BRICS bonds, decomposed with singlePass=True; only for coverage_brics. | count |
| `fmcs_timed_out` | FMCS reached its fixed 10-second comparison budget; the returned result may be incomplete. | boolean |

## `document_offsets.csv`

| Column | Meaning | Unit or values |
|---|---|---|
| `document_a` | First document in the ordered pair. | text |
| `document_b` | Second document in the ordered pair. | text |
| `shared_entries` | Entries shared by both documents. | count |
| `median_offset` | Median B-minus-A log-unit offset. | log units |
| `mean_offset` | Mean B-minus-A log-unit offset. | log units |
| `std_offset` | Sample SD of paired offsets. | log units |
| `sign_agreement` | Fraction of differences sharing the median sign. | fraction 0–1 |
| `systematic_offset` | Whether the pair meets offset criteria. | boolean |
| `touched_entry_count` | Size of the union of entries in the two documents. | count |
| `total_evaluated_pairs` | Pairs meeting the minimum-shared threshold. | count |
| `total_observed_pairs` | Pairs sharing at least one entry. | count |
| `total_possible_pairs` | All document pairs. | count |

## `excluded_records.csv`

| Column | Meaning | Unit or values |
|---|---|---|
| `molecule_identity_key` | Stable entry key used for map rows and evidence joins. | text |
| `molecule_identity_source` | How the entry identifier was obtained. | parent ID or canonical SMILES fallback |
| `canonical_smiles` | Canonical isomeric SMILES. | SMILES |
| `representation_smiles` | Structure used for fingerprint representation after the selected policy. | SMILES |
| `record_canonical_smiles` | Unique canonical record structures joined with semicolons. | SMILES list |
| `structure_policy` | Selected parent/as-recorded normalization policy. | text |
| `structure_changed` | Whether parent normalization changed the record structure. | 0 or 1 |
| `parent_extraction_excluded` | Whether parent extraction removed atoms. | 0 or 1 |
| `activity_id` | Activity record identifier(s). | text |
| `assay_type` | Assay type(s) among retained records. | text |
| `source_row` | 0-based source row carried through curation. | integer |
| `exclusion_reason` | Semicolon-separated curation exclusion reasons. | text |
| `original_smiles` | Original SMILES field. | SMILES |
| `activity_pchembl_raw` | Source pChEMBL value used for aggregation. | −log10 mol/L |
| `action_type` | ChEMBL activity action type. | text |
| `activity_comment` | ChEMBL activity comment. | text |
| `activity_properties` | ChEMBL activity properties field. | text |
| `assay_chembl_id` | ChEMBL assay identifier. | text |
| `assay_description` | Assay description. | text |
| `assay_variant_accession` | Assay variant accession. | text |
| `assay_variant_mutation` | Assay variant mutation. | text |
| `bao_endpoint` | BioAssay Ontology endpoint identifier. | text |
| `bao_format` | BioAssay Ontology format identifier. | text |
| `bao_label` | BioAssay Ontology format label. | text |
| `data_validity_comment` | ChEMBL data-validity comment. | text |
| `data_validity_description` | ChEMBL data-validity description. | text |
| `document_chembl_id` | Source document identifier. | text |
| `document_journal` | Source journal. | text |
| `document_year` | Source publication year. | year |
| `modality` | ChEMBL molecule modality. | text |
| `molecule_chembl_id` | Source molecule identifier. | text |
| `molecule_pref_name` | Preferred molecule name. | text |
| `parent_molecule_chembl_id` | Source parent identifier. | text |
| `pchembl_value` | ChEMBL pChEMBL value. | −log10 mol/L |
| `potential_duplicate` | ChEMBL duplicate-citation flag. | boolean |
| `qudt_units` | QUDT unit identifier. | text |
| `record_id` | Source record identifier. | text |
| `relation` | Unstandardized relation. | text |
| `src_id` | ChEMBL source-data identifier. | text |
| `standard_flag` | ChEMBL standard-value flag. | boolean |
| `standard_relation` | Standardized relation. | text |
| `standard_text_value` | Standardized text value. | text |
| `standard_type` | Standardized activity type. | text |
| `standard_units` | Standardized activity units. | text |
| `standard_upper_value` | Standardized upper value. | text |
| `standard_value` | Standardized activity value. | text |
| `target_chembl_id` | Target ChEMBL identifier. | text |
| `target_organism` | Target organism. | text |
| `target_pref_name` | Preferred target name. | text |
| `target_tax_id` | Target taxonomy identifier. | text |
| `text_value` | Source text value. | text |
| `toid` | ChEMBL target or document object identifier. | text |
| `type` | Activity value type. | text |
| `units` | Unstandardized activity units. | text |
| `uo_units` | Unit Ontology identifier. | text |
| `upper_value` | Unstandardized upper value. | text |
| `value` | Unstandardized activity value. | text |
| `ligand_efficiency_bei` | Binding efficiency index. | text |
| `ligand_efficiency_le` | Ligand efficiency. | text |
| `ligand_efficiency_lle` | Lipophilic ligand efficiency. | text |
| `ligand_efficiency_sei` | Surface efficiency index. | text |
| `ligand_efficiency` | Combined ligand efficiency value. | text |
| `action_type_action_type` | Raw action-type field from the export. | text |
| `action_type_description` | Action-type description. | text |
| `action_type_parent_type` | Parent action type. | text |
| `raw_activity_value` | Original activity value before unit conversion; the input column is also retained. | input unit |
| `raw_activity_unit` | Original per-row unit, or the selected common unit, before conversion. | M, mM, uM, µM, μM, nM, pM |

## `grid_coordinates.csv`

| Column | Meaning | Unit or values |
|---|---|---|
| `molecule_identity_key` | Stable entry key used for map rows and evidence joins. | text |
| `molecule_identity_source` | How the entry identifier was obtained. | parent ID or canonical SMILES fallback |
| `canonical_smiles` | Canonical isomeric SMILES. | SMILES |
| `smiles` | SMILES column consumed by the map workflow. | SMILES |
| `representation_smiles` | Structure used for fingerprint representation after the selected policy. | SMILES |
| `record_canonical_smiles` | Unique canonical record structures joined with semicolons. | SMILES list |
| `structure_policy` | Selected parent/as-recorded normalization policy. | text |
| `structure_changed` | Whether parent normalization changed the record structure. | 0 or 1 |
| `parent_extraction_excluded` | Whether parent extraction removed atoms. | 0 or 1 |
| `n_record_structures` | Number of distinct record structures for the entry. | count |
| `activity_pchembl` | Displayed median pChEMBL value. | −log10 mol/L |
| `activity_class` | Activity class from the displayed pChEMBL and thresholds. | inactive, medium, active; not_measured for unmeasured map rows |
| `n_records` | Retained source records for the entry. | count |
| `pchembl_min` | Minimum retained pChEMBL value. | −log10 mol/L |
| `pchembl_max` | Maximum retained pChEMBL value. | −log10 mol/L |
| `pchembl_median` | Median retained pChEMBL value. | −log10 mol/L |
| `pchembl_std` | Sample SD of retained pChEMBL values. | −log10 mol/L |
| `pchembl_range` | Maximum minus minimum pChEMBL. | log units |
| `has_inactive_measurement` | Whether any retained record is inactive. | 0 or 1 |
| `has_active_measurement` | Whether any retained record is active. | 0 or 1 |
| `crosses_clean_boundary` | Whether records cross an activity threshold. | 0 or 1 |
| `is_large_variation` | Whether record variation exceeds the configured range. | 0 or 1 |
| `is_conflicted` | Wide-record-range flag using all configured annotation rules. | 0 or 1 |
| `has_class_boundary_crossing` | Whether values occur on both sides of a class boundary. | 0 or 1 |
| `n_measurement_classes` | Number of observed activity classes. | count |
| `repeat_status` | Summary category for repeated activity records. | text |
| `wide_record_range` | Display alias for the wide-record-range flag. | 0 or 1 |
| `audit_lower_threshold` | Lower pChEMBL threshold saved with the row. | −log10 mol/L |
| `audit_upper_threshold` | Upper pChEMBL threshold saved with the row. | −log10 mol/L |
| `audit_conflict_range_threshold` | Wide-record-range threshold saved with the row. | log units |
| `source_rows` | 0-based input rows contributing to the entry. | semicolon-separated integers |
| `molecule_id` | Source molecule identifier(s). | text |
| `parent_molecule_id` | Source parent molecule identifier(s). | text |
| `target_id` | Target ChEMBL identifier(s). | text |
| `target_name` | Target display name(s). | text |
| `assay_id` | Assay ChEMBL identifier(s). | text |
| `document_id` | Source document identifier(s). | text |
| `activity_id` | Activity record identifier(s). | text |
| `assay_type` | Assay type(s) among retained records. | text |
| `n_assays` | Distinct retained assays. | count |
| `n_documents` | Distinct source documents. | count |
| `n_evidence_units` | Document-balanced units after optional identical-value merging. | count |
| `top_unit_documents` | Documents represented by the unit with most records. | text |
| `top_unit_record_fraction` | Share of records in the largest evidence unit. | fraction 0–1 |
| `document_balanced_pchembl` | Mean of document-level evidence units. | −log10 mol/L |
| `document_balanced_class` | Activity class of the document-balanced estimate. | text |
| `class_changed_by_document_balance` | Whether balancing changes the displayed class. | 0 or 1 |
| `consensus_scale` | Posterior scale of the entry consensus. | log units |
| `p_consensus_{class}` | Posterior probability that the consensus is in each activity class. | fraction 0–1 |
| `p_new_document_{class}` | Predictive probability for a new document in each class. | fraction 0–1 |
| `class_confidence` | Posterior consensus probability for the displayed class. | fraction 0–1 |
| `class_reproducibility` | Predictive probability for the displayed class in a new document. | fraction 0–1 |
| `discordant_document_pairs` | Flagged document pairs containing both sides of the entry. | count |
| `discordant_documents` | Documents in flagged pairs that share this entry. | text |
| `offset_pair_touched` | Whether this entry occurs in any document of a flagged offset pair. | boolean |
| `status_reason` | Support/display decision: supported, low_support or offset_pair. | text |
| `class_status` | Displayed support status. | supported or uncertain |
| `source_row` | 0-based source row carried through curation. | integer |
| `projection_x_raw` | Unscaled first projection coordinate. | projection units |
| `projection_y_raw` | Unscaled second projection coordinate. | projection units |
| `projection_x` | Normalized first projection coordinate. | fraction 0–1 |
| `projection_y` | Normalized second projection coordinate. | fraction 0–1 |
| `grid_x` | Continuous horizontal grid coordinate. | grid units |
| `grid_y` | Continuous vertical grid coordinate. | grid units |
| `grid_col` | Zero-based horizontal grid column. | integer |
| `grid_row` | Zero-based vertical grid row. | integer |
| `neighbor_recall` | Tie-inclusive representation-neighbor recall on the grid. | fraction 0–1 |
| `neighbor_shared_count` | Neighbors present in both representation and grid sets. | count |
| `representation_neighbor_count` | Representation neighbors within the kth-distance shell. | count |
| `grid_neighbor_count` | Grid neighbors within the kth-distance shell. | count |
| `neighbor_recall_fixed_k` | Expected recall with random tie-breaking at fixed k. | fraction 0–1 |
| `neighbor_kth_similarity` | Tanimoto similarity of the kth Morgan fingerprint neighbor. | similarity 0–1 |
| `neighbor_kth_distance` | Euclidean distance to the kth non-fingerprint representation neighbor. | representation distance |
| `comparison_id` | Identifier used in structure comparisons. | text |
| `structure_group` | Structure-group identifier. | text |
| `structure_group_color` | Group layer color. | hex RGB |
| `scaffold_smiles` | Bemis–Murcko scaffold SMILES. Charge-separated groups can leave representation artifacts after removal of the negative oxygen: sulfoxides may appear as `[SH+]` and N-oxides as `[nH+]`. These scaffold strings are not neutral structures. Sulfones retain `S(=O)(=O)`, so sulfoxides and sulfones yield different scaffolds. | SMILES |
| `name` | Supplied molecule name when present. | text |
| `raw_activity_value` | Original activity value before unit conversion; the input column is also retained. | input unit |
| `raw_activity_unit` | Original per-row unit, or the selected common unit, before conversion. | M, mM, uM, µM, μM, nM, pM |
| `sigma_doc` | Assumed between-unit SD, only emitted for assumed_sigma. | log units |
| `sigma_doc_degrees_of_freedom` | Empty for assumed_sigma because the SD was not estimated. | empty |
| `measurement_status` | Whether the molecule belongs to measured reference data or the unmeasured query set. | measured or not_measured |
| `query_id` | User query identifier or generated Q00001-style identifier; empty for measured rows. | text |
| `nearest_measured_id` | Most similar measured molecule in Morgan fingerprint space; no activity is inferred. | text |
| `nearest_measured_similarity` | Morgan Tanimoto to the nearest measured molecule. | similarity 0–1 |
| `query_group_center_similarity` | Morgan Tanimoto to the assigned reference group center; empty when unassigned. | similarity 0–1 |
| `contains_group_core` | Whether the unmeasured molecule contains its assigned reference group's core; empty when no core/group applies. | boolean or empty |

## `identity_conflicts.csv`

| Column | Meaning | Unit or values |
|---|---|---|
| `molecule_identity_key` | Stable entry key used for map rows and evidence joins. | text |
| `structures` | Normalized structures for the conflicting parent. | semicolon-separated SMILES |
| `record_count` | Records associated with one conflicting identity. | count |

## `manual_entries.csv`

| Column | Meaning | Unit or values |
|---|---|---|
| `canonical_smiles` | Canonical isomeric SMILES. | SMILES |
| `value` | Unstandardized activity value. | text |
| `query_id` | User query identifier or generated Q00001-style identifier; empty for measured rows. | text |
| `input_smiles` | Original query SMILES before structure preparation. | SMILES |
| `label` | User label, used only for display and filtering. | text |
| `unit` | Unit of the entered exact value. | pActivity, M, mM, µM, nM, pM |
| `measured_class` | Measured class, from an exact value or qualitative user choice. | active, medium, inactive |
| `class_source` | Origin of the measured class. | value, user_class |
| `batch` | User-supplied batch or publication used as an evidence unit. | text |
| `note` | User-supplied note. | text |
| `entered_at` | Time the molecule was manually added. | ISO 8601 with timezone |

## `metrics.csv`

| Column | Meaning | Unit or values |
|---|---|---|
| `map_name` | Map/output prefix. | text |
| `representation_type` | Representation used for projection. | morgan, descriptors or user_coordinates |
| `projection_method` | Two-dimensional projection method. | pca, tsne, umap or coordinates |
| `n_samples` | Mapped molecules. | count |
| `k` | Neighbor count used for reliability. | integer |
| `trustworthiness_k` | Neighborhood size used by trustworthiness. | integer |
| `trustworthiness_sample_size` | Maximum deterministic trustworthiness sample. | count |
| `trustworthiness_scope` | Whether trustworthiness was exhaustive or sampled. | text |
| `annotation_neighborhood_policy` | Policy for kth-distance ties. | text |
| `projection_to_grid_trustworthiness` | Projection-to-grid trustworthiness. | fraction 0–1 |
| `projection_grid_knn_overlap` | Mean kNN overlap between projection and grid. | fraction 0–1 |
| `projection_grid_tie_aware_recall` | Tie-aware projection-to-grid recall. | fraction 0–1 |
| `projection_grid_tie_aware_precision` | Tie-aware projection-to-grid precision. | fraction 0–1 |
| `projection_tie_neighborhood_size` | Average tie-inclusive projection neighborhood size. | count |
| `grid_tie_neighborhood_size` | Average tie-inclusive grid neighborhood size. | count |
| `mean_grid_displacement` | Mean projection-to-grid displacement. | grid units |
| `mean_squared_grid_displacement` | Mean squared projection-to-grid displacement. | squared grid units |
| `p95_grid_displacement` | 95th percentile projection-to-grid displacement. | grid units |
| `representation_knn_value_abs_diff` | Mean absolute kNN activity difference in representation space. | pChEMBL |
| `projection_knn_value_abs_diff` | Mean absolute kNN activity difference in projection space. | pChEMBL |
| `grid_knn_value_abs_diff` | Mean absolute kNN activity difference on the grid. | pChEMBL |
| `projection_label_purity` | Projection kNN label purity. | fraction 0–1 |
| `grid_label_purity` | Grid kNN label purity. | fraction 0–1 |
| `projection_distance_metric` | Metric used by the 2D projection. | text |
| `projection_fallback` | Reason PCA replaced UMAP, or blank. | text |
| `representation_distance_metric` | Metric used for representation neighbors. | jaccard or euclidean |
| `grid_assignment_method` | Assignment solver selected for this map. | text |
| `grid_candidate_cells` | Candidate lattice cells. | count |
| `grid_side` | Grid side length. | cells |
| `grid_sizing_strategy` | How grid size was selected. | text |
| `grid_target_occupancy` | Requested cell occupancy. | fraction 0–1 |
| `grid_actual_occupancy` | Occupied fraction of candidate cells. | fraction 0–1 |
| `grid_padding_cells` | Explicit padding used for grid sizing, if any. | cells |
| `coordinate_scaling` | Projection normalization mode. | text |
| `grid_refinement` | Grid refinement mode. | text |
| `representation_grid_tie_aware_recall` | Mean tie-aware neighbor recall. | fraction 0–1 |
| `representation_grid_fixed_k_recall` | Mean fixed-k neighbor recall. | fraction 0–1 |
| `assignment_solver` | Assignment algorithm implementation. | text |
| `assignment_cost_dtype` | Numeric dtype of assignment costs. | text |
| `assignment_candidate_edges` | Candidate edges supplied to the sparse solver. | count |
| `assignment_fallback_edges_added` | Fallback assignment edges added. | count |
| `assignment_fallback_edges_used` | Fallback assignment edges selected. | count |
| `n_measured` | Measured reference molecules used for activity/class diagnostics in a query run. | count |
| `n_query` | Unmeasured molecules included in structural/layout diagnostics. | count |

## `molecule_level.csv`

| Column | Meaning | Unit or values |
|---|---|---|
| `molecule_identity_key` | Stable entry key used for map rows and evidence joins. | text |
| `molecule_identity_source` | How the entry identifier was obtained. | parent ID or canonical SMILES fallback |
| `canonical_smiles` | Canonical isomeric SMILES. | SMILES |
| `smiles` | SMILES column consumed by the map workflow. | SMILES |
| `representation_smiles` | Structure used for fingerprint representation after the selected policy. | SMILES |
| `record_canonical_smiles` | Unique canonical record structures joined with semicolons. | SMILES list |
| `structure_policy` | Selected parent/as-recorded normalization policy. | text |
| `structure_changed` | Whether parent normalization changed the record structure. | 0 or 1 |
| `parent_extraction_excluded` | Whether parent extraction removed atoms. | 0 or 1 |
| `n_record_structures` | Number of distinct record structures for the entry. | count |
| `activity_pchembl` | Displayed median pChEMBL value. | −log10 mol/L |
| `activity_class` | Activity class from the displayed pChEMBL and thresholds. | inactive, medium, active; not_measured for unmeasured map rows |
| `n_records` | Retained source records for the entry. | count |
| `pchembl_min` | Minimum retained pChEMBL value. | −log10 mol/L |
| `pchembl_max` | Maximum retained pChEMBL value. | −log10 mol/L |
| `pchembl_median` | Median retained pChEMBL value. | −log10 mol/L |
| `pchembl_std` | Sample SD of retained pChEMBL values. | −log10 mol/L |
| `pchembl_range` | Maximum minus minimum pChEMBL. | log units |
| `has_inactive_measurement` | Whether any retained record is inactive. | 0 or 1 |
| `has_active_measurement` | Whether any retained record is active. | 0 or 1 |
| `crosses_clean_boundary` | Whether records cross an activity threshold. | 0 or 1 |
| `is_large_variation` | Whether record variation exceeds the configured range. | 0 or 1 |
| `is_conflicted` | Wide-record-range flag using all configured annotation rules. | 0 or 1 |
| `has_class_boundary_crossing` | Whether values occur on both sides of a class boundary. | 0 or 1 |
| `n_measurement_classes` | Number of observed activity classes. | count |
| `repeat_status` | Summary category for repeated activity records. | text |
| `wide_record_range` | Display alias for the wide-record-range flag. | 0 or 1 |
| `audit_lower_threshold` | Lower pChEMBL threshold saved with the row. | −log10 mol/L |
| `audit_upper_threshold` | Upper pChEMBL threshold saved with the row. | −log10 mol/L |
| `audit_conflict_range_threshold` | Wide-record-range threshold saved with the row. | log units |
| `source_rows` | 0-based input rows contributing to the entry. | semicolon-separated integers |
| `molecule_id` | Source molecule identifier(s). | text |
| `parent_molecule_id` | Source parent molecule identifier(s). | text |
| `target_id` | Target ChEMBL identifier(s). | text |
| `target_name` | Target display name(s). | text |
| `assay_id` | Assay ChEMBL identifier(s). | text |
| `document_id` | Source document identifier(s). | text |
| `activity_id` | Activity record identifier(s). | text |
| `assay_type` | Assay type(s) among retained records. | text |
| `n_assays` | Distinct retained assays. | count |
| `n_documents` | Distinct source documents. | count |
| `n_evidence_units` | Document-balanced units after optional identical-value merging. | count |
| `top_unit_documents` | Documents represented by the unit with most records. | text |
| `top_unit_record_fraction` | Share of records in the largest evidence unit. | fraction 0–1 |
| `document_balanced_pchembl` | Mean of document-level evidence units. | −log10 mol/L |
| `document_balanced_class` | Activity class of the document-balanced estimate. | text |
| `class_changed_by_document_balance` | Whether balancing changes the displayed class. | 0 or 1 |
| `consensus_scale` | Posterior scale of the entry consensus. | log units |
| `p_consensus_{class}` | Posterior probability that the consensus is in each activity class. | fraction 0–1 |
| `p_new_document_{class}` | Predictive probability for a new document in each class. | fraction 0–1 |
| `class_confidence` | Posterior consensus probability for the displayed class. | fraction 0–1 |
| `class_reproducibility` | Predictive probability for the displayed class in a new document. | fraction 0–1 |
| `discordant_document_pairs` | Flagged document pairs containing both sides of the entry. | count |
| `discordant_documents` | Documents in flagged pairs that share this entry. | text |
| `offset_pair_touched` | Whether this entry occurs in any document of a flagged offset pair. | boolean |
| `status_reason` | Support/display decision: supported, low_support or offset_pair. | text |
| `class_status` | Displayed support status. | supported or uncertain |
| `projection_x_raw` | Unscaled first projection coordinate. | projection units |
| `projection_y_raw` | Unscaled second projection coordinate. | projection units |
| `projection_x` | Normalized first projection coordinate. | fraction 0–1 |
| `projection_y` | Normalized second projection coordinate. | fraction 0–1 |
| `grid_x` | Continuous horizontal grid coordinate. | grid units |
| `grid_y` | Continuous vertical grid coordinate. | grid units |
| `grid_col` | Zero-based horizontal grid column. | integer |
| `grid_row` | Zero-based vertical grid row. | integer |
| `neighbor_recall` | Tie-inclusive representation-neighbor recall on the grid. | fraction 0–1 |
| `neighbor_shared_count` | Neighbors present in both representation and grid sets. | count |
| `representation_neighbor_count` | Representation neighbors within the kth-distance shell. | count |
| `grid_neighbor_count` | Grid neighbors within the kth-distance shell. | count |
| `neighbor_recall_fixed_k` | Expected recall with random tie-breaking at fixed k. | fraction 0–1 |
| `neighbor_kth_similarity` | Tanimoto similarity of the kth Morgan fingerprint neighbor. | similarity 0–1 |
| `neighbor_kth_distance` | Euclidean distance to the kth non-fingerprint representation neighbor. | representation distance |
| `comparison_id` | Identifier used in structure comparisons. | text |
| `structure_group` | Structure-group identifier. | text |
| `structure_group_color` | Group layer color. | hex RGB |
| `scaffold_smiles` | Bemis–Murcko scaffold SMILES. Charge-separated groups can leave representation artifacts after removal of the negative oxygen: sulfoxides may appear as `[SH+]` and N-oxides as `[nH+]`. These scaffold strings are not neutral structures. Sulfones retain `S(=O)(=O)`, so sulfoxides and sulfones yield different scaffolds. | SMILES |
| `name` | Supplied molecule name when present. | text |
| `raw_activity_value` | Original activity value before unit conversion; the input column is also retained. | input unit |
| `raw_activity_unit` | Original per-row unit, or the selected common unit, before conversion. | M, mM, uM, µM, μM, nM, pM |
| `sigma_doc` | Assumed between-unit SD, only emitted for assumed_sigma. | log units |
| `sigma_doc_degrees_of_freedom` | Empty for assumed_sigma because the SD was not estimated. | empty |
| `measurement_status` | Whether the molecule belongs to measured reference data or the unmeasured query set. | measured or not_measured |
| `query_id` | User query identifier or generated Q00001-style identifier; empty for measured rows. | text |
| `nearest_measured_id` | Most similar measured molecule in Morgan fingerprint space; no activity is inferred. | text |
| `nearest_measured_similarity` | Morgan Tanimoto to the nearest measured molecule. | similarity 0–1 |
| `query_group_center_similarity` | Morgan Tanimoto to the assigned reference group center; empty when unassigned. | similarity 0–1 |
| `contains_group_core` | Whether the unmeasured molecule contains its assigned reference group's core; empty when no core/group applies. | boolean or empty |

## `neighbor_evidence_feedback.csv`

| Column | Meaning | Unit or values |
|---|---|---|
| `canonical_smiles` | Canonical isomeric SMILES. | SMILES |
| `activity_pchembl` | Displayed median pChEMBL value. | −log10 mol/L |
| `query_id` | User query identifier or generated Q00001-style identifier; empty for measured rows. | text |
| `label` | User label, used only for display and filtering. | text |
| `nearest_reference_id` | Placement target: highest Morgan Tanimoto, then smallest reference ID. | text |
| `neighbor_evidence` | Evidence from similar tie-inclusive top-five reference neighbors; never a prediction. | none, insufficient, active, medium, inactive, mixed |
| `similar_reference_count` | Top-five neighbors including ties meeting group_similarity. | count |
| `evidence_reference_classes` | Classes contributing to neighbor evidence, including observed-class fallback. | semicolon-separated classes |
| `neighbor_comparison` | Own measured class compared with consistent neighbor evidence only. | agrees, disagrees, no_comparison |
| `measured_class` | Measured class, from an exact value or qualitative user choice. | active, medium, inactive |
| `class_source` | Origin of the measured class. | value, user_class |
| `batch` | User-supplied batch or publication used as an evidence unit. | text |
| `note` | User-supplied note. | text |
| `result_source` | Origin of the measured result. | upload, manual, recorded |
| `measured_value` | Exact entered value; uploaded multiple-record aggregates use pActivity. | measured_unit |
| `measured_unit` | Unit of the measured value. | pActivity, M, mM, µM, nM, pM |
| `nearest_reference_similarity` | Morgan Tanimoto similarity to the placement reference. | 0 to 1 |
| `recorded_at` | Last manual result save time; empty for uploaded results. | ISO 8601 with timezone |

## `new_molecules.csv`

| Column | Meaning | Unit or values |
|---|---|---|
| `molecule_identity_key` | Stable entry key used for map rows and evidence joins. | text |
| `molecule_identity_source` | How the entry identifier was obtained. | parent ID or canonical SMILES fallback |
| `canonical_smiles` | Canonical isomeric SMILES. | SMILES |
| `smiles` | SMILES column consumed by the map workflow. | SMILES |
| `representation_smiles` | Structure used for fingerprint representation after the selected policy. | SMILES |
| `record_canonical_smiles` | Unique canonical record structures joined with semicolons. | SMILES list |
| `structure_policy` | Selected parent/as-recorded normalization policy. | text |
| `structure_changed` | Whether parent normalization changed the record structure. | 0 or 1 |
| `parent_extraction_excluded` | Whether parent extraction removed atoms. | 0 or 1 |
| `n_record_structures` | Number of distinct record structures for the entry. | count |
| `activity_pchembl` | Displayed median pChEMBL value. | −log10 mol/L |
| `activity_class` | Activity class from the displayed pChEMBL and thresholds. | inactive, medium, active; not_measured for unmeasured map rows |
| `n_records` | Retained source records for the entry. | count |
| `pchembl_min` | Minimum retained pChEMBL value. | −log10 mol/L |
| `pchembl_max` | Maximum retained pChEMBL value. | −log10 mol/L |
| `pchembl_median` | Median retained pChEMBL value. | −log10 mol/L |
| `pchembl_std` | Sample SD of retained pChEMBL values. | −log10 mol/L |
| `pchembl_range` | Maximum minus minimum pChEMBL. | log units |
| `has_inactive_measurement` | Whether any retained record is inactive. | 0 or 1 |
| `has_active_measurement` | Whether any retained record is active. | 0 or 1 |
| `crosses_clean_boundary` | Whether records cross an activity threshold. | 0 or 1 |
| `is_large_variation` | Whether record variation exceeds the configured range. | 0 or 1 |
| `is_conflicted` | Wide-record-range flag using all configured annotation rules. | 0 or 1 |
| `has_class_boundary_crossing` | Whether values occur on both sides of a class boundary. | 0 or 1 |
| `n_measurement_classes` | Number of observed activity classes. | count |
| `repeat_status` | Summary category for repeated activity records. | text |
| `wide_record_range` | Display alias for the wide-record-range flag. | 0 or 1 |
| `audit_lower_threshold` | Lower pChEMBL threshold saved with the row. | −log10 mol/L |
| `audit_upper_threshold` | Upper pChEMBL threshold saved with the row. | −log10 mol/L |
| `audit_conflict_range_threshold` | Wide-record-range threshold saved with the row. | log units |
| `source_rows` | 0-based input rows contributing to the entry. | semicolon-separated integers |
| `molecule_id` | Source molecule identifier(s). | text |
| `parent_molecule_id` | Source parent molecule identifier(s). | text |
| `target_id` | Target ChEMBL identifier(s). | text |
| `target_name` | Target display name(s). | text |
| `assay_id` | Assay ChEMBL identifier(s). | text |
| `document_id` | Source document identifier(s). | text |
| `activity_id` | Activity record identifier(s). | text |
| `assay_type` | Assay type(s) among retained records. | text |
| `n_assays` | Distinct retained assays. | count |
| `n_documents` | Distinct source documents. | count |
| `n_evidence_units` | Document-balanced units after optional identical-value merging. | count |
| `top_unit_documents` | Documents represented by the unit with most records. | text |
| `top_unit_record_fraction` | Share of records in the largest evidence unit. | fraction 0–1 |
| `document_balanced_pchembl` | Mean of document-level evidence units. | −log10 mol/L |
| `document_balanced_class` | Activity class of the document-balanced estimate. | text |
| `class_changed_by_document_balance` | Whether balancing changes the displayed class. | 0 or 1 |
| `consensus_scale` | Posterior scale of the entry consensus. | log units |
| `p_consensus_{class}` | Posterior probability that the consensus is in each activity class. | fraction 0–1 |
| `p_new_document_{class}` | Predictive probability for a new document in each class. | fraction 0–1 |
| `class_confidence` | Posterior consensus probability for the displayed class. | fraction 0–1 |
| `class_reproducibility` | Predictive probability for the displayed class in a new document. | fraction 0–1 |
| `discordant_document_pairs` | Flagged document pairs containing both sides of the entry. | count |
| `discordant_documents` | Documents in flagged pairs that share this entry. | text |
| `offset_pair_touched` | Whether this entry occurs in any document of a flagged offset pair. | boolean |
| `status_reason` | Support/display decision: supported, low_support or offset_pair. | text |
| `class_status` | Displayed support status. | supported or uncertain |
| `source_row` | 0-based source row carried through curation. | integer |
| `projection_x_raw` | Unscaled first projection coordinate. | projection units |
| `projection_y_raw` | Unscaled second projection coordinate. | projection units |
| `projection_x` | Normalized first projection coordinate. | fraction 0–1 |
| `projection_y` | Normalized second projection coordinate. | fraction 0–1 |
| `grid_x` | Continuous horizontal grid coordinate. | grid units |
| `grid_y` | Continuous vertical grid coordinate. | grid units |
| `grid_col` | Zero-based horizontal grid column. | integer |
| `grid_row` | Zero-based vertical grid row. | integer |
| `neighbor_recall` | Tie-inclusive representation-neighbor recall on the grid. | fraction 0–1 |
| `neighbor_shared_count` | Neighbors present in both representation and grid sets. | count |
| `representation_neighbor_count` | Representation neighbors within the kth-distance shell. | count |
| `grid_neighbor_count` | Grid neighbors within the kth-distance shell. | count |
| `neighbor_recall_fixed_k` | Expected recall with random tie-breaking at fixed k. | fraction 0–1 |
| `neighbor_kth_similarity` | Tanimoto similarity of the kth Morgan fingerprint neighbor. | similarity 0–1 |
| `neighbor_kth_distance` | Euclidean distance to the kth non-fingerprint representation neighbor. | representation distance |
| `comparison_id` | Identifier used in structure comparisons. | text |
| `structure_group` | Structure-group identifier. | text |
| `structure_group_color` | Group layer color. | hex RGB |
| `scaffold_smiles` | Bemis–Murcko scaffold SMILES. Charge-separated groups can leave representation artifacts after removal of the negative oxygen: sulfoxides may appear as `[SH+]` and N-oxides as `[nH+]`. These scaffold strings are not neutral structures. Sulfones retain `S(=O)(=O)`, so sulfoxides and sulfones yield different scaffolds. | SMILES |
| `name` | Supplied molecule name when present. | text |
| `core_member` | Whether the member is covered by the group core. | boolean |
| `core_rule` | Core-selection rule. | coverage, most_frequent or fmcs |
| `core_smiles` | Selected group core SMILES. Charge-separated groups can leave representation artifacts after removal of the negative oxygen: sulfoxides may appear as `[SH+]` and N-oxides as `[nH+]`. These core strings are not neutral structures. Sulfones retain `S(=O)(=O)`, so sulfoxides and sulfones yield different cores. | SMILES |
| `core_smarts` | FMCS core SMARTS when available. | SMARTS |
| `core_coverage_target` | Requested core-coverage target. | fraction 0–1 |
| `raw_activity_value` | Original activity value before unit conversion; the input column is also retained. | input unit |
| `raw_activity_unit` | Original per-row unit, or the selected common unit, before conversion. | M, mM, uM, µM, μM, nM, pM |
| `sigma_doc` | Assumed between-unit SD, only emitted for assumed_sigma. | log units |
| `sigma_doc_degrees_of_freedom` | Empty for assumed_sigma because the SD was not estimated. | empty |
| `measurement_status` | Whether the molecule belongs to measured reference data or the unmeasured query set. | measured or not_measured |
| `query_id` | User query identifier or generated Q00001-style identifier; empty for measured rows. | text |
| `input_smiles` | Original query SMILES before structure preparation. | SMILES |
| `status` | Query inclusion decision, one row per input row. | mapped, invalid_smiles, duplicate_query or in_reference |
| `first_query_id` | First query identifier with the same canonical structure; empty otherwise. | text |
| `reference_id` | Measured molecule with the same canonical structure; empty otherwise. | text |
| `tanimoto` | Morgan fingerprint Tanimoto using the map radius and bit count, including for descriptor maps. | similarity 0–1 |
| `nearest_measured_id` | Most similar measured molecule in Morgan fingerprint space; no activity is inferred. | text |
| `nearest_measured_similarity` | Morgan Tanimoto to the nearest measured molecule. | similarity 0–1 |
| `query_group_center_similarity` | Morgan Tanimoto to the assigned reference group center; empty when unassigned. | similarity 0–1 |
| `contains_group_core` | Whether the unmeasured molecule contains its assigned reference group's core; empty when no core/group applies. | boolean or empty |
| `label` | User label, used only for display and filtering. | text |
| `is_new` | Identifies the current new-molecule overlay. | boolean |
| `nearest_reference_id` | Placement target: highest Morgan Tanimoto, then smallest reference ID. | text |
| `placement_distance` | Euclidean distance from the target reference cell. | grid cells |
| `neighbor_evidence` | Evidence from similar tie-inclusive top-five reference neighbors; never a prediction. | none, insufficient, active, medium, inactive, mixed |
| `similar_reference_count` | Top-five neighbors including ties meeting group_similarity. | count |
| `supported_reference_count` | Supported similar neighbors; zero when reference support was not assessed. | count |
| `supported_reference_classes` | Classes of supported similar neighbors; empty when support was not assessed. | semicolon-separated classes |
| `evidence_reference_count` | Neighbors contributing to evidence; uses observed classes when support was not assessed. | count |
| `evidence_reference_classes` | Classes contributing to neighbor evidence, including observed-class fallback. | semicolon-separated classes |
| `reference_support_assessed` | Whether reference class support was assessed. | boolean |
| `neighbor_comparison` | Own measured class compared with consistent neighbor evidence only. | agrees, disagrees, or empty |
| `user_median` | Own retained-record median, when available. | pActivity |
| `reference_median` | Median of the structurally identical reference molecule. | pActivity |
| `median_difference` | Own median minus structurally identical reference median. | log units |
| `measured_class` | Measured class, from an exact value or qualitative user choice. | active, medium, inactive |
| `class_source` | Origin of the measured class. | value, user_class |
| `batch` | User-supplied batch or publication used as an evidence unit. | text |
| `note` | User-supplied note. | text |
| `result_source` | Origin of the measured result. | upload, manual, recorded |
| `measured_value` | Exact entered value; uploaded multiple-record aggregates use pActivity. | measured_unit |
| `measured_unit` | Unit of the measured value. | pActivity, M, mM, µM, nM, pM |
| `nearest_reference_similarity` | Morgan Tanimoto similarity to the placement reference. | 0 to 1 |
| `recorded_at` | Last manual result save time; empty for uploaded results. | ISO 8601 with timezone |
| `contains_group_core_rule` | Rule used for the saved new-molecule core match; API inspection uses the current map rule. | coverage, coverage_brics, most_frequent, fmcs |

## `new_molecules_check.csv`

| Column | Meaning | Unit or values |
|---|---|---|
| `canonical_smiles` | Canonical isomeric SMILES. | SMILES |
| `source_row` | 0-based source row carried through curation. | integer |
| `query_id` | User query identifier or generated Q00001-style identifier; empty for measured rows. | text |
| `input_smiles` | Original query SMILES before structure preparation. | SMILES |
| `status` | Query inclusion decision, one row per input row. | mapped, invalid_smiles, duplicate_query or in_reference |
| `first_query_id` | First query identifier with the same canonical structure; empty otherwise. | text |
| `reference_id` | Measured molecule with the same canonical structure; empty otherwise. | text |

## `new_molecules_excluded_records.csv`

| Column | Meaning | Unit or values |
|---|---|---|
| `exclusion_reason` | Semicolon-separated curation exclusion reasons. | text |

## `new_molecules_neighbors.csv`

| Column | Meaning | Unit or values |
|---|---|---|
| `activity_pchembl` | Displayed median pChEMBL value. | −log10 mol/L |
| `activity_class` | Activity class from the displayed pChEMBL and thresholds. | inactive, medium, active; not_measured for unmeasured map rows |
| `n_evidence_units` | Document-balanced units after optional identical-value merging. | count |
| `class_confidence` | Posterior consensus probability for the displayed class. | fraction 0–1 |
| `class_status` | Displayed support status. | supported or uncertain |
| `query_id` | User query identifier or generated Q00001-style identifier; empty for measured rows. | text |
| `rank` | Similarity rank among measured neighbors, including every tie at the fifth neighbor. | integer starting at 1 |
| `measured_id` | Measured reference molecule identifier. | text |
| `tanimoto` | Morgan fingerprint Tanimoto using the map radius and bit count, including for descriptor maps. | similarity 0–1 |

## `query_check.csv`

| Column | Meaning | Unit or values |
|---|---|---|
| `canonical_smiles` | Canonical isomeric SMILES. | SMILES |
| `source_row` | 0-based source row carried through curation. | integer |
| `query_id` | User query identifier or generated Q00001-style identifier; empty for measured rows. | text |
| `input_smiles` | Original query SMILES before structure preparation. | SMILES |
| `status` | Query inclusion decision, one row per input row. | mapped, invalid_smiles, duplicate_query or in_reference |
| `first_query_id` | First query identifier with the same canonical structure; empty otherwise. | text |
| `reference_id` | Measured molecule with the same canonical structure; empty otherwise. | text |

## `query_neighbors.csv`

| Column | Meaning | Unit or values |
|---|---|---|
| `activity_pchembl` | Displayed median pChEMBL value. | −log10 mol/L |
| `activity_class` | Activity class from the displayed pChEMBL and thresholds. | inactive, medium, active; not_measured for unmeasured map rows |
| `n_evidence_units` | Document-balanced units after optional identical-value merging. | count |
| `class_confidence` | Posterior consensus probability for the displayed class. | fraction 0–1 |
| `class_status` | Displayed support status. | supported or uncertain |
| `query_id` | User query identifier or generated Q00001-style identifier; empty for measured rows. | text |
| `rank` | Similarity rank among measured neighbors, including every tie at the fifth neighbor. | integer starting at 1 |
| `measured_id` | Measured reference molecule identifier. | text |
| `tanimoto` | Morgan fingerprint Tanimoto using the map radius and bit count, including for descriptor maps. | similarity 0–1 |

## `retained_records.csv`

| Column | Meaning | Unit or values |
|---|---|---|
| `molecule_identity_key` | Stable entry key used for map rows and evidence joins. | text |
| `molecule_identity_source` | How the entry identifier was obtained. | parent ID or canonical SMILES fallback |
| `canonical_smiles` | Canonical isomeric SMILES. | SMILES |
| `representation_smiles` | Structure used for fingerprint representation after the selected policy. | SMILES |
| `record_canonical_smiles` | Unique canonical record structures joined with semicolons. | SMILES list |
| `structure_policy` | Selected parent/as-recorded normalization policy. | text |
| `structure_changed` | Whether parent normalization changed the record structure. | 0 or 1 |
| `parent_extraction_excluded` | Whether parent extraction removed atoms. | 0 or 1 |
| `activity_id` | Activity record identifier(s). | text |
| `assay_type` | Assay type(s) among retained records. | text |
| `source_row` | 0-based source row carried through curation. | integer |
| `exclusion_reason` | Semicolon-separated curation exclusion reasons. | text |
| `original_smiles` | Original SMILES field. | SMILES |
| `activity_pchembl_raw` | Source pChEMBL value used for aggregation. | −log10 mol/L |
| `action_type` | ChEMBL activity action type. | text |
| `activity_comment` | ChEMBL activity comment. | text |
| `activity_properties` | ChEMBL activity properties field. | text |
| `assay_chembl_id` | ChEMBL assay identifier. | text |
| `assay_description` | Assay description. | text |
| `assay_variant_accession` | Assay variant accession. | text |
| `assay_variant_mutation` | Assay variant mutation. | text |
| `bao_endpoint` | BioAssay Ontology endpoint identifier. | text |
| `bao_format` | BioAssay Ontology format identifier. | text |
| `bao_label` | BioAssay Ontology format label. | text |
| `data_validity_comment` | ChEMBL data-validity comment. | text |
| `data_validity_description` | ChEMBL data-validity description. | text |
| `document_chembl_id` | Source document identifier. | text |
| `document_journal` | Source journal. | text |
| `document_year` | Source publication year. | year |
| `modality` | ChEMBL molecule modality. | text |
| `molecule_chembl_id` | Source molecule identifier. | text |
| `molecule_pref_name` | Preferred molecule name. | text |
| `parent_molecule_chembl_id` | Source parent identifier. | text |
| `pchembl_value` | ChEMBL pChEMBL value. | −log10 mol/L |
| `potential_duplicate` | ChEMBL duplicate-citation flag. | boolean |
| `qudt_units` | QUDT unit identifier. | text |
| `record_id` | Source record identifier. | text |
| `relation` | Unstandardized relation. | text |
| `src_id` | ChEMBL source-data identifier. | text |
| `standard_flag` | ChEMBL standard-value flag. | boolean |
| `standard_relation` | Standardized relation. | text |
| `standard_text_value` | Standardized text value. | text |
| `standard_type` | Standardized activity type. | text |
| `standard_units` | Standardized activity units. | text |
| `standard_upper_value` | Standardized upper value. | text |
| `standard_value` | Standardized activity value. | text |
| `target_chembl_id` | Target ChEMBL identifier. | text |
| `target_organism` | Target organism. | text |
| `target_pref_name` | Preferred target name. | text |
| `target_tax_id` | Target taxonomy identifier. | text |
| `text_value` | Source text value. | text |
| `toid` | ChEMBL target or document object identifier. | text |
| `type` | Activity value type. | text |
| `units` | Unstandardized activity units. | text |
| `uo_units` | Unit Ontology identifier. | text |
| `upper_value` | Unstandardized upper value. | text |
| `value` | Unstandardized activity value. | text |
| `ligand_efficiency_bei` | Binding efficiency index. | text |
| `ligand_efficiency_le` | Ligand efficiency. | text |
| `ligand_efficiency_lle` | Lipophilic ligand efficiency. | text |
| `ligand_efficiency_sei` | Surface efficiency index. | text |
| `ligand_efficiency` | Combined ligand efficiency value. | text |
| `action_type_action_type` | Raw action-type field from the export. | text |
| `action_type_description` | Action-type description. | text |
| `action_type_parent_type` | Parent action type. | text |
| `raw_activity_value` | Original activity value before unit conversion; the input column is also retained. | input unit |
| `raw_activity_unit` | Original per-row unit, or the selected common unit, before conversion. | M, mM, uM, µM, μM, nM, pM |

## `structure_members.csv`

| Column | Meaning | Unit or values |
|---|---|---|
| `molecule_identity_key` | Stable entry key used for map rows and evidence joins. | text |
| `canonical_smiles` | Canonical isomeric SMILES. | SMILES |
| `smiles` | SMILES column consumed by the map workflow. | SMILES |
| `representation_smiles` | Structure used for fingerprint representation after the selected policy. | SMILES |
| `molecule_id` | Source molecule identifier(s). | text |
| `parent_molecule_id` | Source parent molecule identifier(s). | text |
| `grid_x` | Continuous horizontal grid coordinate. | grid units |
| `grid_y` | Continuous vertical grid coordinate. | grid units |
| `grid_col` | Zero-based horizontal grid column. | integer |
| `grid_row` | Zero-based vertical grid row. | integer |
| `comparison_id` | Identifier used in structure comparisons. | text |
| `structure_group` | Structure-group identifier. | text |
| `structure_group_color` | Group layer color. | hex RGB |
| `scaffold_smiles` | Bemis–Murcko scaffold SMILES. Charge-separated groups can leave representation artifacts after removal of the negative oxygen: sulfoxides may appear as `[SH+]` and N-oxides as `[nH+]`. These scaffold strings are not neutral structures. Sulfones retain `S(=O)(=O)`, so sulfoxides and sulfones yield different scaffolds. | SMILES |
| `name` | Supplied molecule name when present. | text |
| `core_member` | Whether the member is covered by the group core. | boolean |
| `core_rule` | Core-selection rule. | coverage, most_frequent or fmcs |
| `core_smiles` | Selected group core SMILES. Charge-separated groups can leave representation artifacts after removal of the negative oxygen: sulfoxides may appear as `[SH+]` and N-oxides as `[nH+]`. These core strings are not neutral structures. Sulfones retain `S(=O)(=O)`, so sulfoxides and sulfones yield different cores. | SMILES |
| `core_smarts` | FMCS core SMARTS when available. | SMARTS |
| `core_coverage_target` | Requested core-coverage target. | fraction 0–1 |
| `measurement_status` | Whether the molecule belongs to measured reference data or the unmeasured query set. | measured or not_measured |
| `query_id` | User query identifier or generated Q00001-style identifier; empty for measured rows. | text |
| `query_group_center_similarity` | Morgan Tanimoto to the assigned reference group center; empty when unassigned. | similarity 0–1 |
| `contains_group_core` | Whether the unmeasured molecule contains its assigned reference group's core; empty when no core/group applies. | boolean or empty |

## `structure_variants.csv`

| Column | Meaning | Unit or values |
|---|---|---|
| `molecule_identity_key` | Stable entry key used for map rows and evidence joins. | text |
| `smiles` | Parent-structure SMILES writing observed for this identity. | SMILES |
| `record_count` | Total retained records for this identity. | count |
| `normalized_smiles` | Standardized structure form for this SMILES writing. | SMILES |
| `structure_record_count` | Retained records using this SMILES writing. | count |
| `standard_inchikey` | Standard InChIKey shared by the writing variants for this identity. | InChIKey |
| `selected` | Whether this writing normalizes to the map's selected structure (the normalized form with the most records). | boolean |

## Optional core-rule outputs

When building with `coverage_brics`, the existing structure-member tables also use that `core_rule` value.
Their legacy column-dictionary text remains unchanged for byte-identical default exports.
Map-page switches write only `core_rules/<stem>_cores_<rule>.csv` and the matching `.json`.
The JSON records the rule, build coverage target, RDKit version, each group's `structure_sha256`,
the original run manifest SHA-256, and total elapsed seconds. Structure hashes use the grouping
protocol's sorted canonical-structure serialization. Existing run files, manifest and ZIP are unchanged.
The comparison CSV uses the first 12 hexadecimal digits of the group identifier SHA-256 for `<group-hash>`.
BRICS columns are present for BRICS rule exports and in the four-rule comparison table.
Comparison ZIPs for non-build rules also live under `core_rules/` and include the rule in their names.

## Optional own-data outputs and metadata

Reference filenames retain the run's `<name>_` prefix and never change on insertion.
New results live in `new_molecules/<input-stem>/`, with numeric suffixes for existing names.
`new_molecules_check.csv` and `new_molecules_neighbors.csv` retain the legacy query column contracts.
`new_molecules.csv` contains all checked molecules; placement fields are empty for excluded rows.
Measured input also writes `new_molecules_excluded_records.csv`: every original input column,
followed by `exclusion_reason`; the header is retained when no records are excluded.
Group overlay summaries add `n_query_members` and `n_new_members`; saved reference groups stay unchanged.
Generic reference records with unit conversion or assumed sigma also export
`retained_records.csv`, `molecule_level.csv` and `curation_report.json`.

Unit-conversion reports contain `excluded_reasons` counts and the selected unit/column.
Fits and manifests identify assumed SDs with `method=assumed_sigma`, `assumed_sigma`,
and `evidence_unit_source` (the supplied column or `one unit per compound`).
With no evidence-unit column, `document_offset_detection` is `not applicable`.
The existing automatic `literature_fallback` remains a separate method.

The separate new-molecule manifest records the input and reference-manifest hashes, parameters and timing.
Reference layout, metrics, diagnostics and original ZIP contain reference molecules only.
New-molecule neighbors summarize observed reference evidence, never predictions.
Legacy query-named schema sections below describe the retained low-level query functions.

## Manual-entry and measured-result feedback outputs

`manual_entries.csv` has one row per manually added molecule. Its SHA-256 and row count
are recorded in `new_molecules_manifest.json`. `neighbor_evidence_feedback.csv` includes
all currently measured new molecules, including uploaded results. Saving, modifying or
clearing a manual result rewrites this table and the current collection ZIP. An explicit
feedback download can also create it for an unchanged uploaded collection. File-only
insertion retains its original outputs until a manual action or feedback download.

Class-only results use `class_source=user_class` and leave numeric activity empty.
`no_comparison` covers absent, insufficient or mixed neighbor evidence. These tables
record observations and their agreement with reference evidence; no activity is predicted.
Timestamps use ISO 8601 with a timezone; unknown uploaded recording times remain blank.
Browser refresh retains the current session. After restarting, **Open a previous
result** restores the saved active collection, manual entries and measured results.

Exact CSV column order:

- `manual_entries.csv`: `query_id`, `input_smiles`, `canonical_smiles`, `value`, `unit`, `measured_class`, `class_source`, `batch`, `label`, `note`, `entered_at`.
- `neighbor_evidence_feedback.csv`: `query_id`, `canonical_smiles`, `label`, `result_source`, `measured_value`, `measured_unit`, `activity_pchembl`, `measured_class`, `class_source`, `neighbor_evidence`, `neighbor_comparison`, `nearest_reference_id`, `nearest_reference_similarity`, `similar_reference_count`, `evidence_reference_classes`, `batch`, `note`, `recorded_at`.

## Own-data CLI options

| Option | Meaning |
|---|---|
| `--value-unit {p,M,mM,uM,nM,pM}` | Default `p` preserves existing pActivity values; concentrations become −log10 mol/L. |
| `--unit-column COL` | Per-row concentration units, also accepting µM and μM; case sensitive, whitespace trimmed. Mutually exclusive with explicit `--value-unit`. |
| `--assumed-sigma S` | Finite positive SD in log units for record-table support; incompatible with molecule-only input or disabled support. |
| `--query-file PATH` | Insert unmeasured CSV/TSV molecules into a completed Morgan reference grid, up to 1,000 after checking. Measured new molecules use the browser. |
| `--query-smiles-column COL` | Override automatic SMILES-column detection. |
| `--query-id-column COL` | Query IDs, unique within the input; default generated IDs preserve input order. |


## Session cache (browser runs)

A successful browser build additionally writes `session/reference.pkl.gz`,
`session/state.pkl.gz` and `session/session_manifest.json`. The reference blob
contains the complete built result and source-lookup data; the state blob retains
current core-rule caches, comparisons, selected rule, output paths and active
new-molecule collection. The active collection stores `session/overlay.pkl.gz`
and its own `session_manifest.json` under `new_molecules/<set>/`. Input frames are
reconstructed from the byte-preserved run `input.csv`.

The JSON sidecar records SHA-256 hashes for snapshot files, `source_tree_digest`,
Python/pandas/numpy/RDKit versions, run-manifest SHA-256 and creation time. Both
snapshots and existing scientific outputs are verified before opening; changed
CSV files, missing files or incompatible versions require regeneration. Stored
paths are relative, so a complete run folder can move to another output root.
Opening is read-only and restores cached rules and the active collection without
recomputing; core/feedback operations subsequently update their session state.

These compressed caches occupy additional disk space and are excluded from the
existing run manifest and ZIP. They are not new scientific result columns.
Deleting `session/` preserves exports but prevents restoring that session; delete
the entire run folder when the result is no longer needed. The import-page list
only scans this launcher's output root and shows its latest 20 compatible results.

Browser core computation uses parallel loky workers from 64 groups, with a default of
`min(8, os.cpu_count() - 1)` (at least 1), overridden by `CHEMGRIDMAP_JOBS`.
Library `jobs=1` and CLI `--jobs 1` remain serial by default and do not read that
environment variable. FMCS is always serial
because load changes its timeout solution. No parallel setting changes output
column order, core selection or pre-existing exported files. Browser-only marker
sizes/tooltips do not change the exported `map_with_new_molecules.svg`.

## Offset display default

`offset_display` defaults to `none`: `offset_pair_touched` and
`document_offsets.csv` are still written, and touched entries retain the
support-layer marker, but `class_status` and `status_reason` follow support
alone. Use `offset_display="uncertain"` or `--offset-display uncertain` to
show touched entries as uncertain.

## Audit export subsets

`selected_molecule.csv` contains one row from `grid_coordinates.csv` with exactly the same columns.

`source_records.csv` contains filtered rows from `retained_records.csv` with exactly the same columns.

`comparison_members.csv` contains up to 12 rows from `grid_coordinates.csv` for the selected
structure group, with exactly the same columns.
