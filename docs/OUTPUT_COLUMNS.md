# Output columns

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
| `activity_class` | Activity class from the displayed pChEMBL and thresholds. | inactive, medium or active |
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
| `scaffold_smiles` | Bemis–Murcko scaffold SMILES. | SMILES |

## `identity_conflicts.csv`

| Column | Meaning | Unit or values |
|---|---|---|
| `molecule_identity_key` | Stable entry key used for map rows and evidence joins. | text |
| `structures` | Normalized structures for the conflicting parent. | semicolon-separated SMILES |
| `record_count` | Records associated with one conflicting identity. | count |

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
| `activity_class` | Activity class from the displayed pChEMBL and thresholds. | inactive, medium or active |
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
| `scaffold_smiles` | Bemis–Murcko scaffold SMILES. | SMILES |

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
| `scaffold_smiles` | Bemis–Murcko scaffold SMILES. | SMILES |
| `name` | Supplied molecule name when present. | text |
| `core_member` | Whether the member is covered by the group core. | boolean |
| `core_rule` | Core-selection rule. | coverage, most_frequent or fmcs |
| `core_smiles` | Selected group core SMILES. | SMILES |
| `core_smarts` | FMCS core SMARTS when available. | SMARTS |
| `core_coverage_target` | Requested core-coverage target. | fraction 0–1 |
