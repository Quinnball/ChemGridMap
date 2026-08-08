# ChemGridMap

ChemGridMap turns a ChEMBL activity export into an inspectable
one-molecule-per-cell chemical map. A single command can validate the downloaded
records, aggregate repeated measurements to molecule-level median pChEMBL
values, calculate a molecular representation and two-dimensional projection,
assign every molecule to a unique grid cell, and export the final map.

The same program also accepts already curated molecule tables, external learned
embedding columns, and existing two-dimensional coordinates. Representation,
projection, curation reports, and final grid coordinates remain separate in the
outputs so that the transformation can be audited.

![ChemGridMap example using the 770-molecule project dataset](docs/chemgridmap_demo.png)

## Installation

Create a clean environment and install the package from this directory:

```bash
conda create -n chemgridmap -c conda-forge python=3.11 rdkit
conda activate chemgridmap
pip install -e .
```

Install UMAP support when it is needed:

```bash
pip install -e ".[umap]"
```

## ChEMBL CSV to chemical map

On the ChEMBL target page, filter the activities to the endpoint of interest,
include the pChEMBL and activity-quality columns, and download the results as a
CSV. For example, a CHEMBL205 IC50 export can be processed directly:

```bash
chemgridmap chembl205_activities.csv \
  --input-format chembl \
  --target-id CHEMBL205 \
  --activity-type IC50 \
  --representation morgan \
  --projection pca \
  --output-dir chembl205_map \
  --name chembl205
```

PCA works with the base installation. For UMAP:

```bash
pip install -e ".[umap]"
chemgridmap chembl205_activities.csv \
  --input-format chembl \
  --target-id CHEMBL205 \
  --activity-type IC50 \
  --representation morgan \
  --projection umap \
  --output-dir chembl205_map \
  --name chembl205
```

The ChEMBL importer recognizes common web-export and API column names,
including `Molecule ChEMBL ID`, `Smiles`, `pChEMBL Value`,
`Target ChEMBL ID`, `Standard Type`, `Standard Relation`, `Standard Units`,
`Data Validity Comment`, and `Potential Duplicate`.

The default curation policy:

- retains the requested target and activity type;
- requires finite pChEMBL values;
- retains exact (`=`) activities in nM when those columns are present;
- retains blank or `Manually validated` data-validity comments;
- removes rows flagged by ChEMBL as potential duplicate citations;
- removes missing or invalid SMILES;
- uses RDKit canonical isomeric SMILES as molecule identity without salt stripping;
- aggregates repeated measurements by median pChEMBL;
- records the number, minimum, maximum, standard deviation, and range of measurements;
- flags repeated molecules when the pChEMBL range is greater than 1.0 or measurements
  span the inactive (`<=6`) and active (`>=7`) cores.

Use `--keep-potential-duplicates` only when retaining ChEMBL duplicate-citation
flags is an intentional analysis decision. If an audit column is absent, the
program continues but records a warning in the curation report.

## Output files

A ChEMBL run writes:

- `*_molecule_level.csv`: curated one-row-per-molecule table;
- `*_retained_records.csv`: activity records used for aggregation;
- `*_excluded_records.csv`: excluded rows and explicit reasons;
- `*_curation_report.json`: detected columns, parameters, counts, and warnings;
- `*_grid_coordinates.csv`: normalized projection and assigned grid cells;
- `*_metrics.csv`: projection-to-grid fidelity and local annotation metrics;
- `*.svg`, `*.png`, and `*.pdf`: molecule-resolved chemical maps.

## Bundled example

The bundled example is the 770-molecule ChEMBL-derived dataset used in the
accompanying manuscript. It contains molecule-level pChEMBL values, activity
classes, a 128-dimensional external embedding, and the corresponding
precomputed UMAP coordinates.

```bash
chemgridmap examples/chembl205_embedding_umap.csv \
  --smiles-col smiles \
  --value-col activity_pchembl \
  --label-col activity_class \
  --x-col projection_x \
  --y-col projection_y \
  --grid-padding 20 \
  --output-dir gridmap_output \
  --name chembl205_embedding_umap
```

The original 966-row project table is also supplied as
`examples/chembl205_raw_activity_legacy.csv`. It can be processed end to end:

```bash
chemgridmap examples/chembl205_raw_activity_legacy.csv \
  --input-format chembl \
  --target-id CHEMBL205 \
  --activity-type IC50 \
  --representation morgan \
  --projection pca \
  --output-dir gridmap_output \
  --name chembl205_direct
```

This legacy file predates retention of the complete ChEMBL export metadata, so
the audit report correctly warns that target, relation, units, validity, and
potential-duplicate fields cannot be rechecked from that file.

## ChEMBL 37 validation

The direct-download workflow was independently tested against the official
ChEMBL 37 REST API export for CHEMBL205 IC50 activities. The current release
returned 2,281 raw records; the default audit retained 1,100 records and
generated 892 unique molecules without curation warnings. All 966 legacy rows
and all 770 legacy structures were recovered from the current release.

The archived validation input and outputs are under
`validation_data/chembl37_chembl205/`. The UMAP example can be reproduced with:

```bash
chemgridmap validation_data/chembl37_chembl205/chembl37_chembl205_ic50_raw.csv \
  --input-format chembl \
  --target-id CHEMBL205 \
  --activity-type IC50 \
  --representation morgan \
  --projection umap \
  --grid-padding 20 \
  --output-dir validation_data/chembl37_chembl205/chemgridmap_output_umap \
  --name chembl37_chembl205_ic50_umap
```

## Curated tables and existing coordinates

Existing coordinates can be passed directly:

```bash
chemgridmap molecules.csv \
  --smiles-col smiles \
  --x-col umap_x \
  --y-col umap_y \
  --label-col activity_class \
  --output-dir gridmap_output
```

## Python

```python
import pandas as pd
from chemgridmap import build_grid_map

data = pd.read_csv("examples/chembl205_embedding_umap.csv")
result = build_grid_map(
    data,
    output_dir="gridmap_output",
    name="chembl205_embedding_umap",
    value_col="activity_pchembl",
    label_col="activity_class",
    x_col="projection_x",
    y_col="projection_y",
    grid_padding=20,
)

print(result.metrics)
```

To recompute the two-dimensional layout from the supplied 128-dimensional
embedding:

```bash
pip install -e ".[umap]"
chemgridmap examples/chembl205_embedding_umap.csv \
  --value-col activity_pchembl \
  --label-col activity_class \
  --representation embedding \
  --projection umap \
  --grid-padding 20 \
  --output-dir gridmap_output \
  --name chembl205_recomputed_umap
```

## Input notes

The default input is a CSV file with one molecule per row and a `smiles`
column. Repeated canonical structures raise an error by default because the
generic input route should not silently decide how measurements are aggregated.
Use `--duplicate-policy first` only when keeping the first row is scientifically
appropriate. Use `--input-format chembl` for the explicit, audited ChEMBL
median-aggregation workflow.

Activity colors are recognized automatically for labels named `inactive`,
`medium`, and `active`. Other labels receive a categorical palette. If a
continuous value is supplied without a label column, the default activity
classes are `<= 6`, `6-7`, and `>= 7`; both thresholds can be changed from the
command line.

Grid assignment removes overlap but is not lossless. Inspect
projection-to-grid trustworthiness, exact k-nearest-neighbor overlap, and the
continuous projection coordinates together with the final map.
`--grid-padding` controls the number of unused candidate cells around the
square-root grid size. Small examples usually need a value of 1-3; the
manuscript experiment used 20 for 770 molecules.

## Paper example

The script in `paper/run_comparisons.py` rebuilds the projection and
representation comparisons from the included molecule-level table:

```bash
python paper/run_comparisons.py \
  examples/chembl205_embedding_umap.csv \
  --output paper/output
```

The raw legacy activity table, molecule-level example, and provenance notes are
provided in `examples/`.

## License

ChemGridMap source code is released under the BSD 3-Clause License. The
ChEMBL-derived example table is provided separately under the ChEMBL
Attribution-ShareAlike terms described in `examples/DATASET.md`.
