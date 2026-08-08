# ChEMBL 37 validation data

This directory contains the compact provenance bundle used to validate the
direct ChEMBL CSV workflow reported in the accompanying manuscript.

- `chembl37_chembl205/chembl37_chembl205_ic50_raw.csv` is the flattened
  CHEMBL205 IC50 activity export retrieved from the ChEMBL 37 REST API.
- `chembl37_chembl205/chembl_status.json` records the ChEMBL service release
  metadata observed during retrieval.
- `chembl37_chembl205/historical_vs_chembl37_comparison.json` summarizes the
  comparison with the archived 966-record project table.

The paginated API responses and generated map files are intentionally excluded
from version control because they can be regenerated from the CSV with the
command documented in the project README.
