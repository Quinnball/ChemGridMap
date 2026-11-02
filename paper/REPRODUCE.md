# Reproduction map

Run every command from the repository root after restoring the validation data.
The commands use the locked environment described in `paper/README.md`.

## Validation data

The validation package is archived at DOI `{{validation data DOI}}`. It contains
the following directories at paths relative to the repository root:

| Restored path | Contents and producer |
| --- | --- |
| `paper/strengthening/data/chembl37_api/` | Raw ChEMBL 37 activity downloads and selection manifests used by E1-A, E1-B, E2, E2b, E4, E12b, E15 and E16. |
| `paper/strengthening/data/chembl37_api_validation/` | Raw validation-target downloads used by E11. |
| `paper/strengthening/output/` | Archived outputs from E1-A through E16, including `as_executed_pre05/E2c/` and `as_executed_pre05/E6/`. |
| `paper/strengthening/qa/` | Caches and QA records from the same stages, including the E2 curation cache. |
| `paper/output/revision_v6/` | Frozen three-target curation, saved layouts and earlier manuscript outputs. |
| `paper/output/reliability_v1/` | Frozen reliability layouts, activity, core, case and scale inputs consumed by later analyses. |
| `paper/figures/` | Archived Figure 2–6 inputs, browser captures and figure source tables. |

Restore without overwriting local files:

```bash
curl -L --fail --output /tmp/chemgridmap-validation-data.zip \
  "https://doi.org/{{validation data DOI}}"
unzip -n /tmp/chemgridmap-validation-data.zip -d .
make paper
```

`make paper` restores the smaller 0.3.0 layout export from
`validation_data/chemgridmap-0.3.0-validation.zip` when it is absent. The
validation ZIP supplies the ChEMBL 37 API inputs and the upstream outputs that
are intentionally not committed.

## Stage commands

Each stage has one Make target:

| Stage or result | Command |
| --- | --- |
| Current three-target layouts, validation and package tests | `make paper` |
| E1-A activity evidence | `make paper-e1a` |
| E1-B neighbors and cores | `make paper-e1b` |
| E2 curation, predictions and primary metrics | `make paper-e2` |
| E2b tests and E2b R1 time split | `make paper-e2b` |
| E2c sampled-target layouts | `make paper-e2c` |
| E4a-E4e chance, fold, occupancy and release checks | `make paper-e4` |
| E5 through E5d | `make paper-e5 paper-e5b paper-e5c paper-e5d` |
| E6 fixed-k diagnostics | `make paper-e6` |
| E7 dependency and reporting checks | `make paper-e7` |
| E9 calibration synthesis | `make paper-e9` |
| E10 disconnected-layout comparison | `make paper-e10` |
| E11 new-target validation | `make paper-e11` |
| E12a rule calibration | `make paper-e12a` |
| E12b secondary scores and descriptions | `make paper-e12b` |
| E12c core stereochemistry sensitivity | `make paper-e12c` |
| E15 leak-free reassessment | `make paper-e15` |
| E16 identity sensitivity | `make paper-e16` |
| E22 calibration and agreement by displayed class | `make paper-e22` |
| Detailed-target `reliability_v1` outputs | See the dedicated commands below |
| T1 timing and scaling | `make paper-t1` |
| Fig. 2 map inputs and all manuscript figures | `make paper-fig2-data paper-figures` |
| Fig. 2–6 from restored archived inputs | `make paper-figures-archived` |
| Fig. 6 browser captures | `make paper-figure6-capture` |
| Current structural comparison figures | `make paper-structure` |
| Whole reproducible chain, including figures | `make paper-all` |

`make paper-e2b` runs the accepted default E2b analysis and then
`e2b_analysis.py --r1-outputs`. The R1 step is required for the reported
25-dataset, 1,230-prediction time split and for
`hypothesis_tests_by_stratum.csv`. It does not change the default CLI behavior.

## Manuscript and supplementary map

| Location | Archived output | Runnable command |
| --- | --- | --- |
| Abstract: held-out disagreement, support and agreement | `strengthening/output/E1A/`, `E2/` | `make paper-e1a paper-e2` |
| Abstract: later-document and time-split statements | `strengthening/output/E2/`, `E7/`, `E9/`, `E12b/` | `make paper-e2b paper-e7 paper-e9 paper-e12b` |
| Introduction and scale examples | `paper/output/reliability_v1/` | See the `reliability_v1` commands below |
| Implementation: workflow overview | Hand-composed Figure 1 source | Not reproducible; see below |
| Implementation: evidence model | `strengthening/output/E1A/` | `make paper-e1a` |
| Implementation: layouts | `strengthening/output/E1B/`, `E6/` | `make paper-e1b paper-e6` |
| Implementation: cores | `strengthening/output/E1B/`, `T1/`, `E12c/` | `make paper-e1b paper-t1 paper-e12c` |
| Implementation: viewer | `figures/data_fig6/` | `make paper-figure6-capture paper-figures` |
| Evaluation design, data and models | `strengthening/data/chembl37_api/`, `output/E2/` | `make paper-e2` |
| Results: annotated CA II map | `figures/data_fig2/`, `figures/data_final/` | `make paper-fig2-data paper-figures` |
| Results: evidence-model evaluation | `output/E2/`, `E7/` | `make paper-e2 paper-e7` |
| Results: calibration, offsets and display rule | `output/E5c/`, `E5d/`, `E9/`, `E12a/` | `make paper-e5c paper-e5d paper-e9 paper-e12a` |
| Results: new-target validation | `output/E11/`, `E12a/` | `make paper-e11 paper-e12a` |
| Results: neighborhoods and disconnected layouts | `output/E2c/`, `E6/`, `E10/` | `make paper-e2c paper-e6 paper-e10` |
| Results: cores and performance | `output/E1B/`, `T1/`, `figures/data_final/` | `make paper-e1b paper-t1 paper-figures-data` |
| Results: offset association with overconfidence | `output/E15/q5_offset_association.csv` | `make paper-e15` |
| Results: display rule in held-out literature | `output/E15/q1_rule_by_design.csv`, `q1_tests.csv`, `e11_document_overlap.csv` | `make paper-e15` |
| Results: descriptive audit of offset pairs | `output/E15/q3_*` | `make paper-e15` |
| Results: corrected standard-InChIKey identity check | `output/E16/` | `make paper-e16` |
| Results: software comparison | External comparison audit | Not reproducible; see below |
| Fig. 1 | Hand-composed workflow source | Not reproducible; see below |
| Fig. 2a-d | `figures/data_fig2/{default,fmcs}/` | `make paper-fig2-data paper-figures` |
| Fig. 3a-c | `figures/data_final/3*.csv` | `make paper-figures-data paper-figures` |
| Fig. 4a-d | `figures/data_final/4*.csv` | `make paper-figures-data paper-figures` |
| Fig. 5a-e | `figures/data_final/5*.csv`, `paper/output/reliability_v1/` | `make paper-figures-data paper-figures` |
| Fig. 6a-d | `figures/data_fig6/` | `make paper-figure6-capture paper-figures` |
| Table 1 | External comparison audit | Not reproducible; see below |
| Table 2 | `strengthening/output/E2/hypothesis_tests.csv` | `make paper-e2` |
| Table 3 | `strengthening/output/E2/`, `E11/`, `E12a/` | `make paper-e2b paper-e11 paper-e12a` |
| Table 3 note | `strengthening/output/E2/`, `E11/`, `E12a/`, `E15/` | `make paper-e2b paper-e11 paper-e12a paper-e15` |
| S1-S4: curation, random sample and three detailed targets | `output/E1A/`, `E2/`, `revision_v6/` | `make paper paper-e1a paper-e2` |
| S2.1: target-selection download logs | `data/chembl37_api/selection_log_*.csv` | `.venv/bin/python paper/strengthening/source/download_chembl37_sample.py` |
| S2.5-S2.6: reviewed identity variants and true conflicts | `output/E16/reviewed_variants.csv`, `all_parent_conflicts.csv` | `make paper-e16` |
| S5-S6: 52 datasets, dependence and intervals | `output/E2/`, `E7/`, `E12b/` | `make paper-e2b paper-e7 paper-e12b` |
| S7-S9: comparators, calibration and model alternatives | `output/E4/`, `E5/`, `E5b/`, `E5c/`, `E5d/`, `E9/`, `E12a/` | `make paper-e4 paper-e5 paper-e5b paper-e5c paper-e5d paper-e9 paper-e12a` |
| S8.16-S8.22: leak-free calibration, display-rule checks and offset audit | `output/E15/` | `make paper-e15` |
| S10: document and near-duplicate folds | `output/E5/` | `make paper-e5` |
| S11: time split | `output/E2/`, `E7/`, `E9/`, `E12b/`, `E4/` | `make paper-e2b paper-e4 paper-e7 paper-e9 paper-e12b` |
| S12: flagged document pairs | `output/E2/`, `E1A/`, `E4/` | `make paper-e1a paper-e2 paper-e4` |
| S13: grids, layouts and fixed-k baselines | `output/E2c/`, `E6/`, `E7/`, `E10/` | `make paper-e2c paper-e6 paper-e7 paper-e10` |
| S14: cores and stereochemistry | `output/E1B/`, `T1/`, `E12c/` | `make paper-e1b paper-t1 paper-e12c` |
| S15: G002 case study | `paper/output/reliability_v1/cases/` | `.venv/bin/python paper/run_reliability_analysis.py --steps cases` |
| S16.2-S16.3: scaling summaries | `output/T1/` | `make paper-t1` |
| S17: interface views | `figures/data_fig6/` | `make paper-figure6-capture paper-figures` |
| S18: manifests | All `strengthening/output/*/manifest.json` | `make paper-all` |
| S19: sources of Table 1 | External comparison audit | Not reproducible; see below |
| S20: planned versus reported | `output/E7/`, `E2/`, `E11/` | `make paper-e2b paper-e7 paper-e11` |
| S21: prospective validation | `output/E11/`, `E12a/` | `make paper-e11 paper-e12a` |
| S21.6-S21.8: leak-free prospective checks and identity sensitivity | `output/E15/`, `E16/` | `make paper-e15 paper-e16` |
| S22: offset-document pair audit | `output/E15/q3_*` | `make paper-e15` |
| E22 class reliability | `output/E22/`; script `e22_class_reliability.py` | `make paper-e22`; main text Table 4 and the calibration section; Additional file 1, S8, S11 and S21 |

For S3.2 and S16.1, see the explicit exceptions below.

The `paper/output/reliability_v1/` tree is generated independently of
`make paper`. The commands from S18 are:

```bash
.venv/bin/python -B paper/run_reliability_analysis.py \
  --steps layouts neighborhoods stability --workers 1
.venv/bin/python -B paper/run_reliability_analysis.py \
  --steps baselines --baseline-phase all --workers 1
.venv/bin/python -B paper/run_reliability_analysis.py \
  --steps activity cores cores-timing cases figure-data results --workers 1
.venv/bin/python -B paper/run_reliability_analysis.py \
  --steps scale --workers 1
```

## Archived runs

A clean rerun used a `git archive` of public commit `969704f`. Its main chain ran
from 06:39 to 11:45 (+0800) on 2026-11-02 and stopped at E9, whose check of the
three detailed maps relied on the earlier default display. After E9 was made to
request the display rule explicitly (`offset_display="uncertain"`), E9 and every
later stage were run in the same checkout from 14:55 to 15:27.

The checkout restored only inputs that no analysis stage produces: the ChEMBL
37 activity downloads and validation downloads, the frozen `revision_v6` and
`reliability_v1` data, and `as_executed_pre05/`. Stage outputs, caches, manifests
and figures were all regenerated. The comparison scripts and the file-by-file
results are in `paper/strengthening/qa/E24/` of the validation data.

The following counts come from `paper/strengthening/qa/E24/summary.csv`. Manifests,
logs and timing metadata are excluded from the judged counts. “Last bits” means
rows, columns and text agree, but at least one numeric field is not byte-equal.

| Stage | Byte-identical | Last bits only (maximum absolute difference) | Expected different | Missing/new |
| --- | ---: | ---: | ---: | ---: |
| E1-A | 11 | — | 0 | 0/0 |
| E1-B | 6 | 1 (`3.8`) | 3 | 3/0 |
| E2 | 20 | — | 1 | 0/0 |
| E2b | — | — | — | — |
| E2c | 4 | — | 0 | 0/0 |
| E4 | 10 | — | 0 | 1/0 |
| E5 | 28 | — | 0 | 0/0 |
| E5b | 20 | 1 (`5.55e-16`) | 0 | 1/0 |
| E5c | 5 | 19 (`6.84e-10`) | 1 | 1/0 |
| E5d | 12 | — | 0 | 1/0 |
| E6 | 6 | — | 0 | 0/0 |
| E7 | 23 | — | 0 | 0/0 |
| E9 | 13 | 1 (`2.84e-14`) | 1 | 0/0 |
| E10 | 7 | — | 0 | 0/0 |
| E11 | 1 | 5 (`3.97e-23`) | 1 | 0/0 |
| E12a | 0 | 9 (`2.22e-16`) | 0 | 0/0 |
| E12b | 9 | — | 0 | 0/0 |
| E12c | 1 | — | 0 | 0/0 |
| E15 | 14 | 1 (`0`) | 0 | 0/0 |
| E16 | 10 | — | 0 | 0/0 |
| E22 | 6 | 2 (`2.22e-16`) | 0 | 0/0 |
| T1 | 7 | 1 (`4.01`) | 201 | 8/1 |
| Figure 2 data | 32 | — | 14 | 2/0 |
| Figure 6 data | 4 | — | 2 | 6/0 |
| Figure 3–5 data | 16 | 2 (`0`) | 1 | 0/0 |

E2b writes its outputs into the E2 and E7 directories, so it has no row of its own. The large T1 difference count is explained by the old archive's
column schema, logs, timing fields and default-display metadata; it is not 201
independent scientific disagreements.

The observed differences fall into the following expected classes:

- **Timing.** T1's wall-clock fields and log files differ by run. Across the six
  fixed sizes, the recorded group timings differ by at most `4.01` in the
  compared summary, and individual `core_timing_by_group.csv` times differ by at
  most `0.411`. E1-B has one 2 s FMCS group change from timed out to complete
  (64 to 63); T1 has two such 2 s groups change from timed out to complete. The
  saved group cores themselves are unchanged.
- **Output format.** The archived T1 runs were written by an earlier software
  release. Their scale outputs lack later fields such as `wide_record_range`,
  `offset_pair_touched`, `status_reason`, `neighbor_recall_fixed_k` and
  `touched_entry_count`, so those CSV schemas cannot be byte-compared with the
  current outputs.
- **Text and paths.** E2 has 86 and E11 has 65 non-numeric cells different in
  `target_inclusion.csv`; these are checkout paths and current identity-failure
  wording, with no numeric difference. E11's five numeric CSV comparisons have
  a maximum absolute difference of `3.97e-23`.
- **Default display.** Figure 2's `document_offsets.csv` now lists all four
  evaluated pairs instead of no rows. The current default also removes “or an
  offset document pair” from the class-support layers and from the Figure 6
  legend. Figure 6 now reads
  `Pale: uncertain (class support below 0.8).`
- **Jitter.** Figure 3 has one scatter arrangement displaced by the floating
  last-bit differences in the figure input tables; the plotted input values
  otherwise agree. This is the only positional change in the pixel comparison of Figs. 2–5.
- **Serialization and metadata.** Manifests, run logs, timing records, archives,
  PDFs and PNG/SVG metadata also differ where they contain times, paths, hashes
  or software provenance. The scientific CSV comparison above is the judged
  result.

The sealed E2 and E11 plans still fix the software version and identity rule.
Their archived runs used the pre-0.5 `normalized-structure` rule, and E15 must
use that rule because R5 verifies the E11 execution. Current public commands make
the choice explicit; the default remains `inchikey`, used by E16.

## Archived provenance

The commit fields were not recorded in the archived manifests. The versions
below come from the per-stage manifests; `not recorded` means the manifest did
not contain that package version.

| Stage or figure | Manifest | Commit | Recorded package version |
| --- | --- | --- | --- |
| E1-A | `paper/strengthening/output/E1A/manifest.json` | not recorded | not recorded (Python 3.12.14) |
| E1-B | `paper/strengthening/output/E1B/manifest.json` | not recorded | not recorded (Python dependencies recorded) |
| E2 | `paper/strengthening/output/E2/manifest.json` | not recorded | not recorded; sealed pre-0.5 run |
| E2c | `paper/strengthening/output/E2c/manifest.json` | not recorded | not recorded |
| E4 | no `manifest.json` | not recorded | not recorded |
| E5 | `paper/strengthening/output/E5/manifest.json` | not recorded | chemgridmap 0.5.0 |
| E5b-E5d | `paper/strengthening/output/E5b/manifest.json`, `E5c/`, `E5d/` | not recorded | not recorded (Python 3.12.14) |
| E6 | `paper/strengthening/output/E6/manifest.json` | not recorded | chemgridmap 0.5.0 |
| E7 | `paper/strengthening/output/E7/manifest.json` | not recorded | chemgridmap 0.5.0 |
| E9 | `paper/strengthening/output/E9/manifest.json` | not recorded | chemgridmap 0.5.0 |
| E10 | `paper/strengthening/output/E10/manifest.json` | not recorded | chemgridmap 0.5.0 |
| E11 | `paper/strengthening/output/E11/manifest.json` | not recorded | not recorded; sealed pre-0.5 run |
| E12a | `paper/strengthening/output/E12a/manifest.json` | not recorded | not recorded (Python 3.12.14) |
| E12b-E12c | `paper/strengthening/output/E12b/manifest.json`, `E12c/` | not recorded | not recorded; source hash 047c269f... |
| E15 | `paper/strengthening/output/E15/manifest.json` | not recorded | chemgridmap 0.5.0 |
| E16 | `paper/strengthening/output/E16/manifest.json` | not recorded | chemgridmap 0.5.0 |
| T1 | `paper/strengthening/output/T1/input_manifest.json`, `timing_metadata.json` | not recorded | not recorded |
| Fig. 2 inputs | no figure manifest | clean rerun at `969704f` | not recorded; current default CLI |
| Fig. 3-5 data | `paper/figures/data_final/data_manifest.json` | not recorded | producer `paper/figures/source/data_final.py` |
| Fig. 6 captures | `paper/figures/data_fig6/capture_report.json` | clean rerun at `969704f` | not recorded; current default settings, 3× device scale |

## Not reproducible

- Table 1 and S19: the software comparison was manually curated from external
  sources; no migrated runnable comparison dataset or script exists.
- Fig. 1: the workflow diagram is a hand-composed manuscript source.
- S3.2: the historical model belongs to an earlier manuscript and its model
  archive was not migrated into this repository.
- S16.1: the CI run record is external execution evidence, not a scientific
  output produced by the local analysis scripts.

The following archived files are not regenerated by the current stage scripts:

- `RESULTS_E1B.md`, `RESULTS_E5b.md`, `RESULTS_E5c.md`, `RESULTS_E5d.md` and
  `RESULTS_T1.md`: text reports written when the stages were first run.
- E1-B `scope_check.json` and `manifest.json`: no current script writes
  `scope_check.json`, while `e1b_reporting.py` requires it.
- E4 `flag_clean_trainfold_R0_bug.csv`: a retained record from an earlier
  erroneous run.
- T1 `ui/` screenshots and `logs/run_10001_*`: manual UI evidence and historical
  failure/retry logs.
- Figure 6 `attempt1_long_inspector/`: a retained failed capture attempt.
- Figure 2 `*_cli.log`: historical CLI logs.

## Dependency order

E1-A → E1-B and E2 → E2b → E5 → E5b → E5c → E5d → E9.
E2c → E6 → E10. E2, E2c and E4 → E7 → E12b.
E5c and E11 → E12a → E15 → E16 and E22. E1-B, T1 and E12b → E12c.
