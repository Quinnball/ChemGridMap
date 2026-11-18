# Post hoc analyses for manuscript v0.7 (release 0.6.0)

These scripts answer reviewer items from the November 2026 review. All of them are post hoc: they
were written after the predeclared analyses (E2, E2c, E11, E15 and E28), whose plans had
hashes recorded, and do not change those results.

| Script | Review item | Input | Output |
| --- | --- | --- | --- |
| `s20_mh.py` | S20 | `strengthening/output/E5b/predictions_e5b.csv` | `s20_mh.json`, `mh_strata.csv` |
| `s20_per_dataset_tests.py` | S20 | `s20_mh.json` | `s20_tests.json` |
| `e28_reference_tables.py`, `run_reference_table.sh` | R11, R12, S18, S19 | The E28 raw ChEMBL exports | One reference table per dataset under `E28_REF` |
| `e28_posthoc.py` | R11, R12, S18, S21, S22 | E28 run (`E28_RUN`) and `E28_REF` | Tables under `E28_OUT` |
| `e28_supported_count.py` | R2 Q3 | `compounds_annotated_<sample>.csv.gz` from `e28_posthoc.py` | `q3_supported_count_<sample>.csv` |
| `e28_s19_sensitivity.py` | S19 | `E28_RUN`, `E28_REF` | `s19_summary_<sample>.csv` |
| `e28_ki_components.py` | R12 | `E28_REF` | Ki datasets joined by shared later documents |
| `s35_binding_only.py` | S35 | The 52 evaluation downloads and `E2/target_inclusion.csv` | `s35_by_dataset.csv`, `s35_summary.json` |
| `r15_bounds.py` | R15 | Full IC50 exports of CA II and thrombin in `validation_data/` | `r15_summary.json`, prediction tables |
| `../cvd_separations.py` | S43 | The rendered map colors | `cvd_separations.{csv,md}` |

`e28_s19_sensitivity.py` recomputes neighbor evidence with similarity thresholds 0.5, 0.6 and 0.7,
k = 3, 5 and 10, and count-based Morgan fingerprints. The default setting (2,048-bit Morgan, radius 2,
threshold 0.6, k = 5, ties kept) must reproduce the E28 labels of the hashed-plan run; the script checks this.

Figures for v0.7 are built by `paper/figures/source/build_v07.py` (Figures 4 and 6 and the Additional
file 1 hERG-pair figure), `paper/figures/source/export_fig1.cjs` (Figure 1, from
`Figure_1_workflow.svg`) and `paper/figures/source/fig7_capture.mjs` followed by `fig7_compose.py`
(Figure 7). See `paper/REPRODUCE.md`.
