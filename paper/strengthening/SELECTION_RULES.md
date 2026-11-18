# Target selection rules for the multi-target evaluation (pre-declared)

- **Declared:** 2026-10-28, by Opus, before any activity data for these targets were downloaded.
- **Purpose:** provide a random, auditable sample of ChEMBL 37 targets for evaluating the evidence model beyond CA II, thrombin and hERG.
- **Changes:** any change after download must be recorded below with its reason and date. The original rules stay in place.

## Source

- ChEMBL web services, `https://www.ebi.ac.uk/chembl/api/data`.
- The release must report `chembl_db_version == "ChEMBL_37"`; otherwise stop.
- The full SQLite dump was not used: the measured transfer rate (24–57 KB/s) would take more than a day for 5.8 GB.

## Candidate targets

All targets with `target_type = SINGLE PROTEIN` and `organism = Homo sapiens`, excluding CHEMBL204, CHEMBL205 and CHEMBL240, which are already analysed.

## Eligibility for download

The count is the API `total_count` of activities with `pchembl_value` not null and the given `standard_type`.

| Endpoint | Eligible if count ≥ |
|---|---|
| IC50 | 1,000 |
| Ki | 1,000 |

## Sampling

- Candidates are sorted by ChEMBL ID, then shuffled with `numpy.random.default_rng(seed).permutation`. The seed is 20261028 for IC50 and 20261029 for Ki.
- Walk the shuffled list in order, query each candidate's count, and select it if eligible. Stop at **50 IC50 targets** and **30 Ki targets**, or when the list is exhausted.
- Every queried candidate is logged with its count, eligible or not.
- In addition, Ki data are downloaded for CHEMBL205, CHEMBL204 and CHEMBL240, so the original three targets also get a second endpoint.

## Download

- All activities of each selected target and endpoint that have a pChEMBL value. Only the fields used by curation, evidence units and time splits are kept:
  - identifiers: activity, molecule, parent molecule, document, assay and target IDs;
  - values: standard type, relation, value and units; pChEMBL;
  - quality flags: data-validity comment and potential-duplicate flag;
  - context: document year, assay type and source ID;
  - structure: canonical SMILES.
- Each file must contain exactly `total_count` records; otherwise retry, and record the outcome.

## Analysis inclusion (applied later, also pre-declared)

- A target–endpoint enters the evaluation if it has **≥ 30 entries with ≥ 2 evidence units** after the software's default curation.
- Targets that fail are reported with the reason. None are silently dropped.
- Curation, evidence-unit definition, thresholds (6 and 7), models and metrics are identical to the three-target analysis, unless the task sheet for that stage records a pre-declared change.

## Amendments

(none)

### Amendment 1 (2026-10-28, before any analysis)

Also download **IC50** for the three original targets (CHEMBL205, CHEMBL204, CHEMBL240) in the same field format. Purposes:

1. Regression check: the API-format records must reproduce the frozen retained-record populations after curation.
2. Provide `document_year` and `assay_type` for hERG, which the frozen hERG file lacks, for the time-split analysis.

These files are not part of the random sample.

### Amendment 2 (2026-10-28 21:10, execution only; rules unchanged)

The sequential count queries ran at about 6 s per target. 446 targets had been queried, 24 of them eligible for IC50, so at that rate the selection would take more than 6 hours. The run was stopped, and the procedure was re-executed with `source/download_chembl37_sample_parallel.py`, which:

- queries counts concurrently in batches of 60, using 6 threads;
- reuses the counts already logged by the sequential run;
- selects targets by walking the same sealed shuffled order and taking the first eligible ones.

The selected set is therefore identical to what the sequential procedure would produce. Downloads also run concurrently, with 4 threads. No eligibility threshold, seed, target count or field set was changed. No activity data of sampled targets had been downloaded or inspected before this amendment.
