# Paper analyses

ChemGridMap 0.5.0 includes the software, fixed validation data and the analysis
sources used for the reliability revision. Scripts write results under
`paper/`; they do not modify the installed package or the original data.

## Environment

From the repository root:

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements-paper.lock.txt
.venv/bin/python -m pip install -e ".[dev,paper]"
```

Python 3.12 and the locked environment are recommended. Public ChEMBL data
remain subject to CC BY-SA 3.0.

## Restore fixed inputs

The 0.3.0 validation ZIP in `validation_data/` supplies the three frozen raw
exports and saved layouts. Restore it only when the local targets are absent:

```bash
unzip -n validation_data/chemgridmap-0.3.0-validation.zip -d .
```

The `-n` option preserves existing files. The archive includes `SHA256.json`
for checking the restored contents. For a source checkout without that ZIP,
obtain the same public ChEMBL 37 activities from the ChEMBL API and select them
with `paper/strengthening/SELECTION_RULES.md`.

Regenerate the three-target layouts, record audit, figures and baseline tests:

```bash
make paper
```

This is equivalent to `make paper-current`. Outputs go to
`paper/output/current`; the frozen v6 inputs and earlier figures remain intact.

The remaining manuscript checks and figure builds are:

```bash
make paper-structure
make paper-legacy
make paper-application
make paper-record-tasks
.venv/bin/python paper/run_reliability_analysis.py \
  --steps layouts neighborhoods stability baselines scale cases figure-data
```

The revision validators, figure scripts and output auditors are listed in the
`Makefile`. Each writes a separate output directory and does not replace the
frozen baseline.

## Multi-target ChEMBL 37 data

The E1/E2 scripts use downloaded ChEMBL 37 activity files and do not commit
those large CSVs. The selection criteria, seeds, endpoint counts and inclusion
rule are recorded in
[`strengthening/SELECTION_RULES.md`](strengthening/SELECTION_RULES.md).
The original ChEMBL 37 API downloads will also be archived with the validation
dataset on Zenodo because the live ChEMBL API will switch to a newer release.

Download the selected targets and endpoints from the ChEMBL API:

```bash
.venv/bin/python paper/strengthening/source/download_chembl37_sample.py
```

For a parallel download:

```bash
.venv/bin/python paper/strengthening/source/download_chembl37_sample_parallel.py
```

The downloader writes to `paper/strengthening/data/chembl37_api/` and records
the ChEMBL release, source URLs, counts and download outcomes. A validation
package may be used instead only if its files follow the same column schema and
selection rules; verify its manifest before analysis.

## E1-A: activity evidence

Generate the three-target held-out predictions, document-unit diagnostics,
support strata, calibration metrics and bootstrap intervals:

```bash
make paper-e1a
```

`--regression-only` checks the API-derived records against the frozen v6
curation without rebuilding outputs.

Main outputs are under `paper/strengthening/output/E1A/`. Detail tables and
fold fits are under `paper/strengthening/qa/E1A/`.

## E1-B: neighbors and cores

Run the neighbor-preservation analysis after E1-A and the frozen reliability
layouts exist:

```bash
make paper-e1b
```

The target also runs the scaffold, coverage-aware and FMCS core comparison,
including the 1, 2, 5 and frozen 10 second budgets.

Outputs are under `paper/strengthening/output/E1B/`. The FMCS timeout is part of
each group record; timed-out searches are not treated as deterministic.

## E2: multi-target evidence evaluation

Build the target-level curation, ten-fold entry predictions, comparators and
bootstrap summaries:

```bash
make paper-e2
```

Build hypothesis tests, time splits, negative controls, pseudo-consensus
calibration, sigma-document summaries and forest-plot data:

```bash
make paper-e2b
```

The target runs the default E2b analysis and then
`e2b_analysis.py --r1-outputs`. The R1 step produces the reported 25-dataset,
1,230-prediction time split and `hypothesis_tests_by_stratum.csv`; default CLI
behavior is unchanged when the script is called without the flag.

E2b follows the pre-declared H7 rule: Gaussian (`G`) remains the default when
the difference from `exclude-flagged` (`X`) is not significant at both
endpoints; X remains an explicit software option.

Outputs are under `paper/strengthening/output/E2/`.

## E2c: sampled-target layout evaluation

Run the pre-declared 80-target-endpoint evaluation of projection,
minimum-displacement and neighbor-refined layouts:

```bash
make paper-e2c
```

Outputs are under `paper/strengthening/output/E2c/`. E6 and E7 reuse the
projection metrics and tests from this run.

## E4: chance, folds, occupancy and release checks

Run E4a-E4e in dependency order:

```bash
make paper-e4
```

Outputs are under `paper/strengthening/output/E4/` and
`paper/strengthening/qa/E4/`.

## E5: evidence-model post hoc analyses

Builds the entry-, document-, near-duplicate- and time-split model
comparisons, paired tests, random-effects summaries, pseudo-consensus
diagnostics and calibration summaries:

```bash
make paper-e5
```

Depends on E2. Outputs go to `paper/strengthening/output/E5/`; intermediate
fold and bootstrap tables go to `paper/strengthening/qa/E5/`. The archived
single-host run took about 8 minutes.

## E5b: cross-dataset calibration

Fits the cross-dataset calibration models and evaluates grid-posterior,
non-parametric and decision-rule variants:

```bash
make paper-e5b
```

Depends on E5, E2 and the two CA II map runs from `paper-fig2-data` (the make
target builds them first). Outputs are under
`paper/strengthening/output/E5b/` and `paper/strengthening/qa/E5b/`. The
isolated single-host run took about 68 minutes.

On macOS, the non-parametric-prior EM uses the Accelerate-compiled float32
kernel built by `e5b_models.py:ensure_native_em_kernel`. Other platforms use
the float64 NumPy implementation. The archived outputs were generated on
macOS, so candidate models that depend on a non-parametric prior (the G-NP
series in E5b and E5c) can differ slightly on another platform. The Gaussian
model G and its predictions are unaffected.

## E5c: offset-aware support

Repeats the support analysis with document-offset-aware candidate models:

```bash
make paper-e5c
```

Depends on E5b and the CA II figure input. Outputs are under
`paper/strengthening/output/E5c/` and `paper/strengthening/qa/E5c/`. The
isolated single-host run took about 8 minutes.

## E5d: display-rule validation

Evaluates the offset-aware display rule, fair comparators and map-level effects:

```bash
make paper-e5d
```

Depends on E5c and E2. Outputs are under
`paper/strengthening/output/E5d/` and `paper/strengthening/qa/E5d/`. The
isolated single-host run took under 1 minute.

## E6: fixed-k layout diagnostics

Builds fixed-k recall and excess-error summaries from the saved reliability
layouts and E2c results:

```bash
make paper-e6
```

Depends on E2c and E4. Outputs are under
`paper/strengthening/output/E6/` and `paper/strengthening/qa/E6/`. The migrated
prepare, regression-diagnostic and analysis sequence took about 40 minutes on
the isolated host.

## E7: reporting and dependency sensitivity

Builds two-stage bootstrap summaries, stratified AUC, shared-document clusters,
negative controls and the planned-versus-reported audit:

```bash
make paper-e7
```

Depends on E2, E2c and E4. Outputs are under
`paper/strengthening/output/E7/` and `paper/strengthening/qa/E7/`. The isolated
single-host run took under 1 minute.

## E9: calibration synthesis

Recomputes the calibration decisions, time splits and effect diagnostics used
by the final calibration discussion:

```bash
make paper-e9
```

Depends on E5, E5b, E5c and E5d. Outputs are under
`paper/strengthening/output/E9/` and `paper/strengthening/qa/E9/`. The isolated
single-host run took under 1 minute.

## E10: disconnected-layout comparison

Compares the 0.5 deterministic reconnect path with the originally executed
spectral path. The script reads UMAP's fitted `reducer.graph_`, checks that
exactly the disconnected datasets changed, reruns their projections twice, and
compares the archived and current L1/L2 results without recomputing the tests:

```bash
make paper-e10
```

Depends on E2c and E6. The pre-0.5 E2c and E6 files are restored under
`paper/strengthening/output/as_executed_pre05/`. Outputs are under
`paper/strengthening/output/E10/` and `paper/strengthening/qa/E10/`.

## E11: new-target validation

Checks the downloaded validation targets, reproduces the E5c regression and
runs the prospective display-rule evaluation:

```bash
make paper-e11
```

Depends on E5c and the validation API data under
`paper/strengthening/data/chembl37_api_validation/`. Outputs are under
`paper/strengthening/output/E11/` and `paper/strengthening/qa/E11/`. The isolated
single-host run took about 6 minutes.

## E12a: rule calibration

Runs leave-one-dataset-out calibration, component tests and the Figure 4d source
table:

```bash
make paper-e12a
```

Depends on E5c, E5d, E7, E9 and E11. Outputs are under
`paper/strengthening/output/E12a/` and `paper/strengthening/qa/E12a/`. The
isolated single-host run took about 1 minute.

## E12b: secondary scores and descriptions

Builds the alternative score, random-effects, time-population, core and data
description tables:

```bash
make paper-e12b
```

Depends on E2, E4, E6, E7, E9 and E1-B. Outputs are under
`paper/strengthening/output/E12b/` and `paper/strengthening/qa/E12b/`. The
isolated single-host run took about 2 minutes.

## E12c: core stereochemistry sensitivity

Run the ported stereochemistry sensitivity analysis with the repository
package:

```bash
make paper-e12c
```

Depends on E1-B, T1 and E12b. The output is
`paper/strengthening/output/E12c/core_stereo_settings.csv`.

## E15: leak-free display-rule reassessment and offset audit

Recomputes display rule P without held-out-document pairs, uses strict
document folds, and audits the S1-S6 features of systematic offset pairs:

```bash
make paper-e15
```

Depends on E2, E5, E5c, E7, E9, E11 and E12a. Outputs are under
`paper/strengthening/output/E15/` and `paper/strengthening/qa/E15/`. The
entry-fold and near-duplicate requests are entry-grouped designs; strict
document folds remove every merged unit containing a held-fold document.
`E15_REUSE_FRAMES=1` may reuse the QA cache during local debugging, but a clean
reproduction should be run without it.

## T1: timing and scaling

Runs the core-rule timing and default-CLI scaling analysis:

```bash
make paper-t1
```

Depends on E1-B and the ChEMBL 37 inputs. Outputs are under
`paper/strengthening/output/T1/`. Timing is machine-dependent and is not a
deterministic scientific output.

## Figure data and figures

Build the figure 3–5 source tables after all dependent analyses:

```bash
make paper-figures-data
```

Generate the two CA II map variants used by Figure 2:

```bash
make paper-fig2-data
```

Figure 6 needs the browser screenshots in `paper/figures/data_fig6/`; regenerate
them with `make paper-figure6-capture` when a local Playwright browser is
available. Then build the vector and raster figures:

```bash
make paper-figures
```

Figure 1 is a manually composed workflow and is not generated by this source
set.

## Run all analyses and figures

The dependency-ordered full target is:

```bash
make paper-all
```

`make paper-posthoc` remains available for the post hoc chain without the
upstream download/curation and figure targets.

## Tests and figure checks

```bash
.venv/bin/python -m pytest -q
.venv/bin/python -m pytest -q \
  paper/strengthening/source/test_e1a.py \
  paper/strengthening/source/test_e1b.py \
  paper/strengthening/source/test_e2.py \
  paper/strengthening/source/test_e2c.py
```

The package tests under `tests/` and the analysis-script tests under `paper/`
are reported as separate groups. In the E17 locked environment they produced
`208 passed, 6 skipped` and `207 passed, 32 skipped`, respectively. Data-dependent
paper tests skip when the validation package is absent.

## Testing with restored validation data

Restore the validation package (DOI `{{validation data DOI}}`) into
`paper/strengthening/data/` and `paper/output/`, then run the analysis and
figure source tests:

```bash
unzip -n validation_data/chemgridmap-0.3.0-validation.zip -d .
.venv/bin/python -m pytest paper/strengthening/source paper/figures/source
```

Tests skip archived-output regressions when their required data or outputs are
not present. The validation package also supplies the ChEMBL 37 API files and
upstream outputs needed to run every stage independently.

The E5, E5b, E5c, E5d, E6, E7, E9, E10 and E11 scripts write a shared manifest
through `paper/strengthening/source/manifest_utils.py`. It records stage and
analysis role, timing, software versions, parameters, repository-relative input
and source paths, and hashes of every stage output. It contains no task hash or
absolute path. E2c, E12a and E12b retain their existing manifest formats.

The layout dependency chain is E2c → E6 → E7 → E10 → E12b.

The main figure builder uses frozen outputs and native RDKit drawings. If the
optional caption source is absent from a public checkout, the relayout script
prints a warning, omits that Markdown file, and still builds the figures.

## Interpretation boundaries

The grid map is a representation-dependent visualization. Grid proximity is
not inferred as chemical or biological similarity. Group membership is defined
before projection. Document merging and offsets are diagnostics, not assay
harmonization. Timing and FMCS results may vary with machine load; seed
variation is reported rather than treated as biological replication.
