# ChEMBL 37 validation data

This directory contains provenance bundles used to validate the direct ChEMBL
CSV workflow and the scalable grid-assignment path reported in the accompanying
manuscript.

- `chembl37_chembl205/chembl37_chembl205_ic50_raw.csv` is the flattened
  CHEMBL205 IC50 activity export retrieved from the ChEMBL 37 REST API.
- `chembl37_chembl205/chembl_status.json` records the ChEMBL service release
  metadata observed during retrieval.
- `chembl37_chembl204/chembl37_chembl204_ic50_raw.csv` contains the D3 IC50
  export, with release and retrieval metadata in the same directory.
- `chembl37_chembl205/historical_vs_chembl37_comparison.json` summarizes the
  comparison with the archived 966-record project table.
- `chembl37_chembl240_large_scale/` contains compact release, target, curation,
  fidelity, sparse-candidate sensitivity, stage-profile, and benchmark
  metadata for the historical 9,628-entry validation, alongside the archived
  raw CSV. The revised parent-normalized analysis yields 9,581 entries and
  writes new results to `paper/output/revision_v4/chembl240/`; it does not
  overwrite the historical summaries.
- The complete ChEMBL 37 KCNH2/hERG IC50 input, retained-record,
  molecule-level, coordinate, and generated-map bundle is distributed as
  `chembl37_chembl240_large_scale.zip` with the GitHub `v0.2.0` release.

The reproduction archive includes the flattened raw CSV files and available
release/target metadata. Some local directories additionally retain API pages;
these pages are not required to rerun the archived CSV analysis.
Generated maps can be regenerated with the commands documented in the project
README.

## Data license and attribution

ChEMBL data are provided by EMBL-EBI under the Creative Commons
Attribution-ShareAlike 3.0 Unported license (CC BY-SA 3.0):
https://creativecommons.org/licenses/by-sa/3.0/. The ChEMBL interface
documentation states the same data license:
https://chembl.gitbook.io/chembl-interface-documentation/about. The derived CSV,
audit, coordinate, and metric files in this directory are redistributed under
those terms. Cite the current ChEMBL paper and identify ChEMBL 37 when reusing
these files. The ChemGridMap source code remains separately licensed under the
BSD 3-Clause License.
