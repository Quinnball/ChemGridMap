# Manuscript analyses

## Current candidate: one entry

From a development checkout installed with `.[dev,umap,paper]`, run:

```bash
make paper
```

This is an alias for `make paper-current`. It restores the public v0.3.0
validation ZIP after SHA-256 checks, refusing to replace differing local files;
recurates all three raw exports; checks both ID and cell queries for all 3,289
mapped molecules and ten altered-input cases; runs the untouched web CSV through
the CLI with UMAP; reconciles 25 linked CA II records with the transcription of
Innocenti et al.'s Tables 2-5; and runs pytest. New evidence lives only in
`paper/output/current`. The source check matches additive and concentration
before comparing endpoints, units and IC50 values. It records pChEMBL arithmetic
differences without overwriting database annotations. The publisher PDF is not
redistributed; source locators and its digest are in `paper/data/`.

The saved v6 projections and negative geometry findings are not optimized or
silently replaced. `claim_checks.json` records the current source reconciliation
alongside the computational checks. Author approval remains a separate step.
The manuscript evaluates record traceability and layout behavior; it makes no
participant-level usability claim and requires no participant data to reproduce.

## Optional task prototype (outside the manuscript)

The earlier exploratory task code is retained but is not run by `make paper`.
Generate its preview explicitly with `.venv/bin/python paper/prepare_tasks.py`.

Open `paper/output/current/task_materials/task_preview.html` to preview the
controlled tasks. Before collecting data, a human must review the answer key,
consent/institutional requirements, participant count and stopping rule, and
provide equal practice. The two modes have identical evidence access. The test
is of display mode within this interface, not superiority to complete external
applications. Responses are downloaded locally, not sent to a server.

After real observations have been collected:

```bash
.venv/bin/python paper/analyze_tasks.py /path/to/responses_*.json \
  --answer-key paper/output/current/task_materials/answer_key.json \
  --output /path/to/private-study-analysis
```

The scorer refuses empty observations and reports participant-level descriptive
results without inferential or general superiority claims. Keep participant
responses out of the public repository unless separately reviewed and authorized.

## Archived evidence and legacy commands

The `0.3.0` application uses archived ChEMBL 37 IC50 records for CHEMBL205,
CHEMBL204 and CHEMBL240. The local reproduction source archive includes all
three raw CSV files. Analysis outputs are written to `paper/output/revision_v4`.
The original rc2 layout manifests remain unchanged; v5 adds record-retrieval
validation, fixed-molecule seed checks and full fingerprint-to-grid metrics.
The v6 follow-up adds a full metadata-preserving pandas task baseline and
an assay-condition audit. Earlier numerical layout results are not replaced.

Download the validation archive from the v0.3.0 release and extract it at the
repository root. It supplies the public source CSVs and saved layouts:

```bash
python paper/run_record_tasks.py --source paper/output/revision_v6
python paper/build_revision_figures.py --source paper/output/revision_v6 --output paper/output/revision_v6/figures
```

Both retrieval routes receive the same coordinate and retained-record tables.
The pandas baseline explicitly joins by molecular identity and verifies the
record links and median. It is not a deliberately information-poor control.
All 3,289 mapped entries are queried by ID and cell; three altered-input
controls test error handling. This evaluates correctness, not user speed.
The acetazolamide follow-up uses archived assay descriptions and the metadata
and abstract for DOI 10.1016/j.bmcl.2008.02.008. It does not claim a new
experimental finding or independently transcribed full-text tables.

Use a Python 3.12 environment and the supplied dependency snapshot:

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements-paper.lock.txt
.venv/bin/python -m pip install -e '.[dev,umap,paper]'
make test
make paper
```

`run_revision_validation.py` performs record auditing, parent normalization,
multi-seed projection, assignment controls, occupancy sensitivity and the
as-recorded ablation. `validate_adaptive_assignment.py` calibrates adaptive
sparse matching against the same dense layouts. `build_revision_figures.py`
uses the saved results and native RDKit structures, without generated data.
`audit_revision_outputs.py` checks saved counts, medians, grid cells and SVG
identifiers. A real public-CLI check, when present, is compared with the
manuscript seed-42 layout as well.

The optional full hERG scale check is separate:

```bash
.venv/bin/python paper/validate_adaptive_assignment.py --large
.venv/bin/python paper/extend_large_candidate_check.py
.venv/bin/python paper/build_revision_figures.py
```

The 1,024-candidate cap does not converge for the saved 9,581-molecule case;
the separate 2,048-cap check is retained rather than silently replacing it.
Timing excludes representation, UMAP, curation and rendering unless a field
explicitly states otherwise. These are software runs, not biological replicates.

After reproducing the baseline and generating its figures, run `make paper-application`
to create `paper/output/revision_v5` without replacing the baseline. The new
`run_application_validation.py` tracks molecular identities across saved seeds
and verifies the same acetazolamide source records through ID and cell selection.
The real website download is supplied separately in the review archive. After
placing it at `paper/output/revision_v5/web_export_map/chembl205_web_raw.csv`, run
the documented ChEMBL-to-map command with name `chembl205_web` in that directory.
`scripts/check_clean_install.py --python PATH_TO_NEW_ENV_PYTHON` then validates
an independently installed wheel, runs pytest, maps that CSV and retrieves its
records. It records missing web-export identifiers rather than filling them in.

## Historical analyses

`run_comparisons.py`, `run_large_scale_validation.py`,
`benchmark_sparse_neighbors.py` and `profile_large_scale_stages.py` are retained
for previous versions. The 770-molecule external-vector demo and the earlier
9,628-entry hERG results are not the current manuscript experiments. No MPNN
model is trained or evaluated in this software note.
