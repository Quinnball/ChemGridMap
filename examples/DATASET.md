# ChEMBL205 molecule-level example

`chembl205_raw_activity_legacy.csv` and `chembl205_embedding_umap.csv` are the
activity-record and molecule-level examples used with ChemGridMap.

## Contents

- 966 legacy activity records
- 963 records with valid SMILES
- 770 unique canonical molecules after median aggregation
- molecule-level pChEMBL values
- inactive, medium and active labels
- 128 external embedding columns (`emb_000` to `emb_127`)
- precomputed UMAP coordinates (`projection_x`, `projection_y`)
- precomputed grid coordinates retained from the manuscript run

The activity classes are:

- inactive: pChEMBL <= 6
- medium: 6 < pChEMBL < 7
- active: pChEMBL >= 7

## Preparation

The local source workbook contained 966 IC50 activity records for the
CHEMBL205 target. Three rows had no usable SMILES. The ChEMBL importer parses
structures with RDKit, groups records by canonical isomeric SMILES, and uses
the median pChEMBL as the molecule-level value, yielding 770 unique molecules.
The incorrect protein-sequence column present in the historical workbook was
not used and is not included in the distributed CSV.

The 128-dimensional columns are an externally generated learned molecular
representation. ChemGridMap does not train that model; it reads the vectors as
an optional input representation.

The workbook metadata records a creation date of 16 September 2025. The exact
ChEMBL release number and full activity-audit columns were not retained in the
historical workbook. The bundled raw file is therefore labelled `legacy`, and
its curation report warns which checks cannot be repeated. A publication
archive should use a newly downloaded full ChEMBL CSV with all audit columns.

## Source and license

The underlying structures and bioactivity records were obtained from ChEMBL,
which makes its data available under the Creative Commons
Attribution-ShareAlike 3.0 Unported License:

https://www.ebi.ac.uk/chembl/

This processed example is distributed under the same CC BY-SA 3.0 terms.
Please retain this attribution when redistributing the table.

Recommended ChEMBL citation:

Zdrazil B, Felix E, Hunter F, et al. The ChEMBL Database in 2023: a drug
discovery platform spanning multiple bioactivity data types and time periods.
Nucleic Acids Research. 2024;52(D1):D1180-D1192.
https://doi.org/10.1093/nar/gkad1004
