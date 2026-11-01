# Reproduction map

Run from the repository root after restoring the validation data. Upstream E1A,
E1B, E2 and E4 stages are restored archive inputs in this source tree; the make
targets below rebuild the post hoc chain and figures from those outputs.

## Manuscript

| Location | Output | Script | Make target |
| --- | --- | --- | --- |
| Abstract | `paper/strengthening/output/E1A/`, `E1B/`, `E2/` | Restored E1A, E1B and E2 analyses | Archived input; no make target |
| Introduction, scale examples | `paper/output/reliability_v1/` | Restored frozen layout outputs | Archived input; no make target |
| Implementation, workflow overview | Hand-composed Figure 1 source | No runnable source in this tree | Not reproduced |
| Implementation, evidence model | `paper/strengthening/output/E1A/` | `e1a_evidence_audit.py` | No target; run the script directly |
| Implementation, layouts | `paper/strengthening/output/E1B/`, `E6/` | `e1b_neighbors.py`, `e6_analysis.py` | `paper-e6` |
| Implementation, cores | `paper/strengthening/output/E1B/` | `e1b_cores.py` | No target; run the script directly |
| Implementation, viewer | `paper/figures/data_fig6/` | `fig6_capture.cjs` | `paper-figure6-capture` |
| Testing and reproducibility | `paper/output/current/`, `paper/strengthening/output/E10/` | Package tests, `e10_analysis.py` | `make test`, `paper-e10` |
| Evaluation design, data and models | `paper/strengthening/data/chembl37_api/`, `output/E2/` | Download and E2 analyses | Archived input; no make target |
| Results, annotated CA II map | `paper/figures/data_fig2/`, `paper/figures/data_final/` | ChemGridMap CLI, `fig2.py` | `paper-fig2-data`, `paper-figures` |
| Results, median-class disagreement | `paper/strengthening/output/E1A/metrics.csv`, `E2/metrics_by_target.csv` | E1A and E2 analyses | Archived input; no make target |
| Results, evidence-model evaluation | `paper/strengthening/output/E2/`, `E7/` | E2, `e7_analysis.py` | `paper-e7` |
| Results, dependence and calibration | `paper/strengthening/output/E7/`, `E9/`, `E12b/` | `e7_analysis.py`, `e9_analysis.py`, `e12b_analysis.py` | `paper-e7 paper-e9 paper-e12b` |
| Results, offset-document display rule | `paper/strengthening/output/E5c/`, `E5d/`, `E9/`, `E12a/` | E5c, E5d, E9 and E12a scripts | `paper-e5c paper-e5d paper-e9 paper-e12a` |
| Results, new-target validation | `paper/strengthening/output/E11/`, `E12a/` | `e11_analysis.py`, `e12a_rule_calibration.py` | `paper-e11 paper-e12a` |
| Results, counting and between-document SD | `paper/strengthening/output/E1A/`, `E2/` | Restored E1A and E2 analyses | Archived input; no make target |
| Results, single-document cells | `paper/strengthening/output/E1A/`, `E2/` | Restored E1A and E2 analyses | Archived input; no make target |
| Results, later documents and time split | `paper/strengthening/output/E2/`, `E7/`, `E9/`, `E12b/` | E2, E7, E9 and E12b scripts | `paper-e7 paper-e9 paper-e12b` |
| Results, consensus and flagged pairs | `paper/strengthening/output/E1A/`, `E2/`, `E4/` | Restored analyses | Archived input; no make target |
| Results, neighborhoods | `paper/strengthening/output/E6/`, `E10/`, `paper/figures/data_final/` | E6, E10, `data_final.py`, `fig5.py` | `paper-e6 paper-e10 paper-figures-data paper-figures` |
| Results, cores | `paper/strengthening/output/E1B/`, `T1/`, `paper/figures/data_final/` | E1-B, `t1_scale.py`, `data_final.py`, `fig5.py` | `paper-t1 paper-figures-data paper-figures` |
| Results, software comparison | External comparison audit | No migrated runnable source | Not reproduced |
| Results, performance | `paper/strengthening/output/T1/` | `t1_scale.py` | `paper-t1` |
| Conclusions | Outputs cited by the Results sections above | Same as cited sections | Same as cited sections |
| Fig. 1 | Hand-composed workflow source | No generating script | Not reproduced |
| Fig. 2a, activity layer | `paper/figures/data_fig2/default/` | CLI, `fig2.py` | `paper-fig2-data paper-figures` |
| Fig. 2b, class-support layer | `paper/figures/data_fig2/default/` | CLI, `fig2.py` | `paper-fig2-data paper-figures` |
| Fig. 2c, neighbor-recall layer | `paper/figures/data_fig2/default/` | CLI, `fig2.py` | `paper-fig2-data paper-figures` |
| Fig. 2d, G001/G002 cores | `paper/figures/data_fig2/{default,fmcs}/` | CLI, `fig2.py` | `paper-fig2-data paper-figures` |
| Fig. 3a, median-class disagreement | `paper/figures/data_final/3a_B0_disagreement.csv` | `data_final.py`, `fig3.py` | `paper-figures-data paper-figures` |
| Fig. 3b, AUC and comparators | `paper/figures/data_final/3b_*.csv` | `data_final.py`, `fig3.py` | `paper-figures-data paper-figures` |
| Fig. 3c, paired log-loss differences | `paper/figures/data_final/3c_*.csv` | `data_final.py`, `fig3.py` | `paper-figures-data paper-figures` |
| Fig. 4a, calibration curves | `paper/figures/data_final/4a_calibration_curves.csv` | `data_final.py`, `fig4.py` | `paper-figures-data paper-figures` |
| Fig. 4b, hERG patent pair | `paper/figures/data_final/4b_hERG_patent_pair.csv` | `data_final.py`, `fig4.py` | `paper-figures-data paper-figures` |
| Fig. 4c, later-document agreement | `paper/figures/data_final/4c_*.csv` | `data_final.py`, `fig4.py` | `paper-figures-data paper-figures` |
| Fig. 4d, display-rule groups | `paper/figures/data_final/4d_*.csv` | `data_final.py`, `fig4.py` | `paper-figures-data paper-figures` |
| Fig. 5a, detailed-target recall | `paper/figures/data_final/5a_recall.csv` | `data_final.py`, `fig5.py` | `paper-figures-data paper-figures` |
| Fig. 5b, layout contrasts | `paper/figures/data_final/5b_layout_contrasts.csv` | `data_final.py`, `fig5.py` | `paper-figures-data paper-figures` |
| Fig. 5c, minimum/refined-window query | `paper/output/reliability_v1/`, `fig5.py` | Frozen layouts, `fig5.py` | `paper-figures` |
| Fig. 5d, 77-dataset correlation | `paper/figures/data_final/5d_layout_correlation.csv` | `e2c_analysis.py`, `e6_analysis.py`, `data_final.py`, `fig5.py` | `paper-e2c paper-e6 paper-figures-data paper-figures` |
| Fig. 5e, core coverage by size | `paper/figures/data_final/5e_core_coverage_by_group_size.csv` | `data_final.py`, `fig5.py` | `paper-figures-data paper-figures` |
| Fig. 6a, grid view | `paper/figures/data_fig6/` | `fig6_capture.cjs`, `fig6.py` | `paper-figure6-capture paper-figures` |
| Fig. 6b, inspector | `paper/figures/data_fig6/` | `fig6_capture.cjs`, `fig6.py` | `paper-figure6-capture paper-figures` |
| Fig. 6c, neighbor view | `paper/figures/data_fig6/` | `fig6_capture.cjs`, `fig6.py` | `paper-figure6-capture paper-figures` |
| Fig. 6d, source records | `paper/figures/data_fig6/` | `fig6_capture.cjs`, `fig6.py` | `paper-figure6-capture paper-figures` |
| Table 1 | External comparison audit | No migrated runnable source | Not reproduced |
| Table 2 | `paper/strengthening/output/E2/hypothesis_tests.csv` | Restored E2 analysis | Archived input; no make target |
| Table 3 | `paper/strengthening/output/E2/hypothesis_tests.csv`, `E11/validation_tests.csv` | Restored E2, `e11_analysis.py` | `paper-e11` |

## Supplementary information

| Location | Output | Script | Make target |
| --- | --- | --- | --- |
| S1, data and curation | `paper/output/revision_v6/`, `output/E2/target_inclusion.csv`, `output/E12b/censored_records.csv` | `make paper`, restored E2, `e12b_analysis.py` | `make paper paper-e12b` |
| Table S1.1 | `paper/output/revision_v6/*/*_curation_report.json` | Frozen package curation | `make paper` |
| Table S1.2 | `paper/output/revision_v6/*/*_curation_report.json` | Frozen package curation | `make paper` |
| Table S1.3 | `output/E12b/censored_records.csv`, `output/E2/target_inclusion.csv` | `e12b_analysis.py`, restored E2 | `paper-e12b` |
| S2, random target sample | `data/chembl37_api/selection_log_*.csv`, `output/E2/target_inclusion.csv` | Target download and restored E2 | Archived input; no make target |
| Table S2.1 | `data/chembl37_api/selection_log_IC50.csv`, `selection_log_Ki.csv` | `download_chembl37_sample_parallel.py` | No target; run the downloader |
| Table S2.2 | `data/chembl37_api/selected_targets.csv`, `output/E2/target_inclusion.csv` | Target download and restored E2 | Archived input; no make target |
| Table S2.3 | `output/E2/target_inclusion.csv` | Restored E2 analysis | Archived input; no make target |
| Table S2.4 | `data/chembl37_api/download_manifest.json`, `output/E12b/selection_summary.csv` | Downloader and `e12b_analysis.py` | `paper-e12b` |
| S3, evidence model and design | `output/E1A/metrics.csv`, `output/E5d/fair_comparison_summary.csv` | Restored E1A, `e5d_analysis.py` | `paper-e5d` |
| Table S3.1 | `output/E1A/metrics.csv`, `output/E5d/fair_comparison_summary.csv` | Restored E1A, `e5d_analysis.py` | `paper-e5d` |
| Table S3.2 | Historical model archive from the earlier manuscript | The archived model was not migrated into this repository | Not reproduced |
| S4, three detailed targets | `output/E1A/` | Restored E1A analysis | Archived input; no make target |
| Table S4.1 | `output/E1A/metrics.csv` | Restored E1A analysis | Archived input; no make target |
| Table S4.2 | `output/E1A/bootstrap_ci.csv` | Restored E1A analysis | Archived input; no make target |
| Table S4.3 | `output/E1A/paired_differences.csv` | Restored E1A analysis | Archived input; no make target |
| Table S4.4 | `output/E1A/tau_influence.csv` | Restored E1A analysis | Archived input; no make target |
| Table S4.5 | `output/E1A/pseudo_consensus.csv` | Restored E1A analysis | Archived input; no make target |
| Table S4.6 | `output/E1A/count_flip_summary.csv` | Restored E1A analysis | Archived input; no make target |
| Table S4.7 | `output/E1A/single_unit_thresholds.csv` | Restored E1A analysis | Archived input; no make target |
| Table S4.8 | `output/E1A/metrics.csv` | Restored E1A analysis | Archived input; no make target |
| Table S4.9 | `output/E1A/count_flip_summary.csv` | Restored E1A analysis | Archived input; no make target |
| S5, 52-dataset results | `output/E2/`, `E7/`, `E12b/` | Restored E2, E7 and E12b scripts | `paper-e7 paper-e12b` |
| Table S5.1 | `output/E2/metrics_by_target.csv`, `forest_plot_data.csv` | Restored E2 analysis | Archived input; no make target |
| Table S5.2 | `output/E2/h1_h2_summary.csv` | Restored E2 analysis | Archived input; no make target |
| Table S5.3 | `output/E2/hypothesis_tests.csv` | Restored E2 analysis | Archived input; no make target |
| Table S5.4 | `output/E7/hypotheses_flag_clean.csv` | `e7_analysis.py` | `paper-e7` |
| Table S5.5 | `output/E7/permutation_by_target.csv` | `e7_analysis.py` | `paper-e7` |
| Table S5.6 | `output/E2/sigma_doc_by_target.csv` | Restored E2 analysis | Archived input; no make target |
| Table S5.7 | `output/E2/agreement_by_display_status.csv` | Restored E2 analysis | Archived input; no make target |
| Table S5.8 | `output/E12b/scores_h4_h7.csv` | `e12b_analysis.py` | `paper-e12b` |
| S6, dependence and intervals | `output/E7/`, `E12b/` | `e7_analysis.py`, `e12b_analysis.py` | `paper-e7 paper-e12b` |
| Table S6.1 | `output/E7/two_stage_bootstrap.csv` | `e7_analysis.py` | `paper-e7` |
| Table S6.2 | `output/E12b/re_h5_h6_alternatives.csv` | `e12b_analysis.py` | `paper-e12b` |
| Table S6.3 | `output/E7/stratified_auc_pooled.csv` | `e7_analysis.py` | `paper-e7` |
| Table S6.4 | `output/E7/shared_document_cluster_tests.csv` | `e7_analysis.py` | `paper-e7` |
| Table S6.5 | `output/E7/exclude_nonrandom_ki_hypotheses.csv` | `e7_analysis.py` | `paper-e7` |
| Table S6.6 | `output/E7/entry_equal_weight_h3.csv` | `e7_analysis.py` | `paper-e7` |
| Table S6.7 | `output/E7/random_effects.csv` | `e7_analysis.py` | `paper-e7` |
| S7, comparators | `output/E5d/fair_comparison_summary.csv`, `E4/e4a_comparison_summary.csv` | `e5d_analysis.py`, restored E4 | `paper-e5d` |
| Table S7.1 | `output/E5d/fair_comparison_summary.csv` | `e5d_analysis.py` | `paper-e5d` |
| Table S7.2 | `output/E4/e4a_comparison_summary.csv` | Restored E4 analysis | Archived input; no make target |
| S8, calibration and display rule | `output/E12a/`, `E9/`, `E5d/`, `E12b/`, `figures/data_final/` | E5d, E9, E12a, E12b, `data_final.py` | `paper-e5d paper-e9 paper-e12a paper-e12b paper-figures-data` |
| Table S8.1 | `output/E12a/calibration_all_designs.csv` | `e12a_rule_calibration.py` | `paper-e12a` |
| Table S8.2 | `figures/data_final/si_4a_without_CHEMBL5113.csv` | `data_final.py` | `paper-figures-data` |
| Table S8.3 | `output/E12a/rule_tests.csv` | `e12a_rule_calibration.py` | `paper-e12a` |
| Table S8.4 | `output/E12a/rule_tests.csv` | `e12a_rule_calibration.py` | `paper-e12a` |
| Table S8.5 | `output/E5d/saturation_loss.csv` | `e5d_analysis.py` | `paper-e5d` |
| Table S8.6 | `output/E9/detailed_maps_rule_p.csv` | `e9_analysis.py` | `paper-e9` |
| Table S8.7 | `output/E9/time_split_by_earlier_units.csv` | `e9_analysis.py` | `paper-e9` |
| Table S8.8 | `output/E9/old_dl_diagnostic.csv`, `output/E12a/calibration_all_designs.csv` | `e9_analysis.py`, `e12a_rule_calibration.py` | `paper-e9 paper-e12a` |
| Table S8.9 | `output/E12a/rule_states_by_design.csv` | `e12a_rule_calibration.py` | `paper-e12a` |
| Table S8.10 | `output/E12a/lodo_named.csv` | `e12a_rule_calibration.py` | `paper-e12a` |
| Table S8.11 | `output/E12a/lodo.csv` | `e12a_rule_calibration.py` | `paper-e12a` |
| Table S8.12 | `output/E12a/time_split_rule.csv` | `e12a_rule_calibration.py` | `paper-e12a` |
| Table S8.13 | `output/E12a/full_map_touched_units.csv` | `e12a_rule_calibration.py` | `paper-e12a` |
| Table S8.14 | `output/E12b/ece_reference.csv` | `e12b_analysis.py` | `paper-e12b` |
| Table S8.15 | `output/E5d/map_impact_summary.csv` | `e5d_analysis.py` | `paper-e5d` |
| S9, post hoc model alternatives | `output/E5/`, `E5b/`, `E5c/`, `E9/` | E5, E5b, E5c, `e9_analysis.py` | `paper-e5 paper-e5b paper-e5c paper-e9` |
| Table S9.1 | `output/E5c/metrics_pooled.csv` | `e5c_analysis.py` | `paper-e5c` |
| Table S9.2 | `output/E5{,b,c}/metrics_pooled.csv` | E5, E5b and E5c scripts | `paper-e5 paper-e5b paper-e5c` |
| Table S9.3 | `output/E5b/decision_rule_numbers.csv`, `E5c/decision_rule_numbers.csv` | E5b and E5c scripts | `paper-e5b paper-e5c` |
| Table S9.4 | `output/E9/rule2_recomputed.csv` | `e9_analysis.py` | `paper-e9` |
| Table S9.5 | `output/E9/decision_by_summary.csv` | `e9_analysis.py` | `paper-e9` |
| S10, document and near-duplicate folds | `output/E5/` | `e5_analysis.py` | `paper-e5` |
| Table S10.1 | `output/E5/metrics_by_dataset.csv` | `e5_analysis.py` | `paper-e5` |
| Table S10.2 | `output/E5/document_folds_paired.csv` | `e5_analysis.py` | `paper-e5` |
| Table S10.3 | `output/E5/near_duplicate_oof_paired.csv` | `e5_analysis.py` | `paper-e5` |
| S11, time split | `output/E2/`, `E7/`, `E9/`, `E12b/`, `E4/` | E2, E7, E9, E12b, restored E4 | `paper-e7 paper-e9 paper-e12b` |
| Table S11.1 | `output/E7/time_split_by_dataset.csv` | `e7_analysis.py` | `paper-e7` |
| Table S11.2 | `output/E2/time_split.csv` | Restored E2 analysis | Archived input; no make target |
| Table S11.3 | `output/E12b/posthoc_tests_two_sided.csv` | `e12b_analysis.py` | `paper-e12b` |
| Table S11.4 | `output/E9/time_split_orexin_tests.csv` | `e9_analysis.py` | `paper-e9` |
| Table S11.5 | `output/E4/time_split_posthoc.csv` | Restored E4 analysis | Archived input; no make target |
| Table S11.6 | `output/E12b/time_split_populations.csv` | `e12b_analysis.py` | `paper-e12b` |
| S12, flagged document pairs | `output/E2/`, `E1A/`, `E4/` | Restored E2, E1A and E4 analyses | Archived input; no make target |
| Table S12.1 | `output/E2/flagged_pairs_census.csv` | Restored E2 analysis | Archived input; no make target |
| Table S12.2 | `output/E1A/flagged_pairs_profile.csv` | Restored E1A analysis | Archived input; no make target |
| Table S12.3 | `output/E4/flag_null.csv` | Restored E4 analysis | Archived input; no make target |
| Table S12.4 | `output/E4/flag_clean_trainfold.csv` | Restored E4 analysis | Archived input; no make target |
| S13, grids and layout baselines | `output/E2c/`, `E6/`, `E7/`, `E10/` | E2c, E6, E7, E10 scripts | `paper-e2c paper-e6 paper-e7 paper-e10` |
| Table S13.1 | `output/E6/recall_fixedk_summary.csv` | `e6_analysis.py` | `paper-e6` |
| Table S13.2 | `output/E6/recall_fixedk_77.csv` | `e6_analysis.py` | `paper-e6` |
| Table S13.3 | `output/E10/l1_l2_change_summary.csv` | `e10_analysis.py` | `paper-e10` |
| Table S13.4 | `output/E10/reproducibility.csv`, `connected_controls.csv` | `e10_analysis.py` | `paper-e10` |
| Table S13.5 | `output/E7/locality_lambda.csv` | `e7_analysis.py` | `paper-e7` |
| Table S13.6 | `output/E6/methods_refinement.md` | `e6_analysis.py` | `paper-e6` |
| Table S13.7 | `output/E6/methods_las_flas.md` | `e6_analysis.py` | `paper-e6` |
| S14, cores | `output/E1B/`, `T1/`, `E12c/` | Restored E1B, `t1_scale.py`, archived E12c | `paper-t1`; E12c archived |
| Table S14.1 | `output/E1B/core_by_group_size.csv` | Restored E1B analysis | Archived input; no make target |
| Table S14.2 | `output/E1B/fmcs_budget_sweep.csv` | Restored E1B analysis | Archived input; no make target |
| Table S14.3 | `output/T1/core_timing_summary.csv` | `t1_scale.py` | `paper-t1` |
| Table S14.4 | `output/E12c/core_stereo_settings.csv` | Archived E12c analysis | Archived input; no make target |
| S15, G002 case study | `paper/output/reliability_v1/cases/` | `reliability_cases.py` via `run_reliability_analysis.py` | `paper/run_reliability_analysis.py --steps cases` |
| Table S15.1 | `paper/output/reliability_v1/cases/case_summary.json` | `reliability_cases.py` via `run_reliability_analysis.py` | `paper/run_reliability_analysis.py --steps cases` |
| Table S15.2 | `paper/output/reliability_v1/cases/G002_CHEMBL20_positions_seed42.csv` | `reliability_cases.py` via `run_reliability_analysis.py` | `paper/run_reliability_analysis.py --steps cases` |
| S16, software and performance | CI record, `output/T1/` | CI, `t1_scale.py` | `paper-t1` |
| Table S16.1 | CI run record | External CI | Not reproduced |
| Table S16.2 | `output/T1/scaling_summary.csv` | `t1_scale.py` | `paper-t1` |
| S17, interface views | `paper/figures/data_fig6/` | `fig6_capture.cjs` | `paper-figure6-capture` |
| Table S17.1 | `paper/figures/data_fig6/` | `fig6_capture.cjs` | `paper-figure6-capture` |
| S18, reproduction | All stage `manifest.json` files | Stage scripts and `manifest_utils.py` | `paper-posthoc` |
| Table S18.1 | `paper/strengthening/output/*/manifest.json`, `T1/*metadata.json` | Stage scripts | `paper-posthoc` |
| S19, sources of Table 1 | External comparison audit | No migrated runnable source | Not reproduced |
| Table S19.1 | External comparison audit | No migrated runnable source | Not reproduced |
| Table S19.2 | Archived task comparison | No migrated runnable source | Not reproduced |
| S20, planned versus reported | `output/E7/planned_vs_reported.csv`, `known_deviations.csv`, E2/E11 outputs | E7, restored E2, E11 | `paper-e7 paper-e11` |
| Table S20.1 | `output/E7/planned_vs_reported.csv` | `e7_analysis.py` | `paper-e7` |
| Table S20.2 | `output/E7/known_deviations.csv` | `e7_analysis.py` | `paper-e7` |
| S21, prospective validation | `output/E11/`, `output/E12a/e11_sensitivity.csv` | E11 and E12a scripts | `paper-e11 paper-e12a` |
| Table S21.1 | `output/E11/target_inclusion.csv` | `e11_analysis.py` | `paper-e11` |
| Table S21.2 | `output/E11/validation_tests.csv` | `e11_analysis.py` | `paper-e11` |
| Table S21.3 | `output/E11/descriptives.csv` | `e11_analysis.py` | `paper-e11` |
| Table S21.4 | `output/E12a/e11_sensitivity.csv` | `e12a_rule_calibration.py` | `paper-e12a` |
| Table S21.5 | `output/E11/map_impact_summary.csv` | `e11_analysis.py` | `paper-e11` |

## Archived runs

- All layout results are the 0.5 deterministic path.
- `output/as_executed_pre05/E2c/` and `output/as_executed_pre05/E6/` contain
  the originally executed spectral-path outputs. E10 uses them only as the
  comparison arm.
- The archived B1c rows were regenerated from source. E5's archived executable
  had rewritten that cache by pooling all datasets; the source-defined
  per-dataset result is authoritative.

## Dependency order

E2c → E6 → E7 → E10 → E12b. The E5 family runs
E5 → E5b → E5c; E9 follows the E5 family; E11 follows E5c.
