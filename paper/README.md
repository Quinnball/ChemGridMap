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
.venv/bin/python paper/strengthening/source/e1a_evidence_audit.py
```

`--regression-only` checks the API-derived records against the frozen v6
curation without rebuilding outputs.

Main outputs are under `paper/strengthening/output/E1A/`. Detail tables and
fold fits are under `paper/strengthening/qa/E1A/`.

## E1-B: neighbors and cores

Run the neighbor-preservation analysis after E1-A and the frozen reliability
layouts exist:

```bash
.venv/bin/python paper/strengthening/source/e1b_neighbors.py
```

Run the scaffold, coverage-aware and FMCS core comparison, including the 1, 2,
5 and frozen 10 second budgets:

```bash
.venv/bin/python paper/strengthening/source/e1b_cores.py
```

Outputs are under `paper/strengthening/output/E1B/`. The FMCS timeout is part of
each group record; timed-out searches are not treated as deterministic.

## E2: multi-target evidence evaluation

Build the target-level curation, ten-fold entry predictions, comparators and
bootstrap summaries:

```bash
.venv/bin/python paper/strengthening/source/e2_analysis.py
```

Build hypothesis tests, time splits, negative controls, pseudo-consensus
calibration, sigma-document summaries and forest-plot data:

```bash
.venv/bin/python paper/strengthening/source/e2b_analysis.py
```

E2b follows the pre-declared H7 rule: Gaussian (`G`) remains the default when
the difference from `exclude-flagged` (`X`) is not significant at both
endpoints; X remains an explicit software option.

Outputs are under `paper/strengthening/output/E2/`.

## E2c: sampled-target layout evaluation

Run the pre-declared 80-target-endpoint evaluation of projection,
minimum-displacement and neighbor-refined layouts:

```bash
.venv/bin/python paper/strengthening/source/e2c_analysis.py
```

Outputs are under `paper/strengthening/output/E2c/`.

## Tests and figure checks

```bash
.venv/bin/python -m pytest -q
.venv/bin/python -m pytest -q \
  paper/strengthening/source/test_e1a.py \
  paper/strengthening/source/test_e1b.py \
  paper/strengthening/source/test_e2.py \
  paper/strengthening/source/test_e2c.py
```

The main figure builder uses frozen outputs and native RDKit drawings. If the
optional caption source is absent from a public checkout, the relayout script
prints a warning, omits that Markdown file, and still builds the figures.

## Interpretation boundaries

The grid map is a representation-dependent visualization. Grid proximity is
not inferred as chemical or biological similarity. Group membership is defined
before projection. Document merging and offsets are diagnostics, not assay
harmonization. Timing and FMCS results may vary with machine load; seed
variation is reported rather than treated as biological replication.
