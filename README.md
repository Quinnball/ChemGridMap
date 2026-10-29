# ChemGridMap

[![Tests](https://github.com/Quinnball/ChemGridMap/actions/workflows/tests.yml/badge.svg)](https://github.com/Quinnball/ChemGridMap/actions/workflows/tests.yml)
[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD--3--Clause-blue.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-3776AB.svg)](https://www.python.org/)

ChemGridMap turns molecular structures and pActivity measurements into an
auditable one-molecule-per-cell chemical map. It preserves source records,
separates structure groups from map position, and exposes activity support,
neighbor recall and group-core coverage as distinct reliability layers.
Use it from Python, the command line or a local browser interface.

![Three reliability layers and core-aligned structural comparison](docs/structure_comparison.png)

## Install

Python 3.10 or newer is required; Python 3.12 is recommended. UMAP is installed
with the core package.

After the 0.5.0 PyPI release:

```bash
python -m pip install chemgridmap
```

For the current GitHub source:

```bash
python -m pip install "git+https://github.com/Quinnball/ChemGridMap.git"
```

For development:

```bash
git clone https://github.com/Quinnball/ChemGridMap.git
cd ChemGridMap
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
```

On Windows PowerShell, activate the environment with
`.venv\Scripts\Activate.ps1`.

## Quick start

Download the activities for ChEMBL target CA II (`CHEMBL205`), include
`pChEMBL Value`, `Standard Relation`, `Standard Units`, document and
quality-flag fields, then run:

```bash
chemgridmap chembl205_ic50.csv \
  --input-format chembl \
  --target-id CHEMBL205 \
  --output-dir caii_map
```

The 0.5.0 defaults are UMAP, neighbor-refined grid placement,
coverage-target group cores, and Gaussian activity-class support. To launch the
same workflow in the local browser:

```bash
chemgridmap-app
```

The browser keeps files on the computer, lets you choose CSV columns, displays
all three reliability layers, and exports the source-record audit with the map.

The package also contains a CA II smoke test:

```bash
python -m chemgridmap.app --smoke-test --output-dir caii_smoke
```

This uses UMAP by default and verifies that CHEMBL20 resolves to its retained
source records and median pChEMBL value.

## Reading the three reliability layers

**Activity evidence.** The main activity layer keeps the observed
molecule-level median and class. A faded/`uncertain` cell means that the
document-balanced evidence does not reach the default 0.8 support threshold;
it is not a missing measurement and is not silently reclassified. Inspect the
cell to compare document counts, evidence units, class probabilities and
document-offset diagnostics.

**Neighbor recall.** The neighbor-recall layer reports how many
fingerprint-space neighbors also appear among the selected cell's tie-aware
grid neighbors. A dark cell has better local neighborhood agreement; a pale
cell is a warning that the 2D layout does not preserve that local fingerprint
neighborhood. Grid adjacency by itself is not evidence of chemical or
biological similarity.

**Core coverage.** For each representation-defined structure group, the default
core is the largest member scaffold or scaffold-network ring fragment that
matches at least 80% of members. The comparison reports the actual coverage
count and fraction. `coverage`, `most_frequent` and optional `fmcs` cores are
selectable; an FMCS timeout is recorded because a partial search can depend on
machine speed and load.

To reproduce the earlier 0.4 layout and core policy explicitly:

```bash
chemgridmap chembl205_ic50.csv --input-format chembl \
  --projection pca --grid-refinement none --core-rule most_frequent \
  --no-activity-confidence --output-dir legacy_0_4
```

## Generic activity records

Non-ChEMBL CSVs can be used directly when they contain a compound ID, SMILES,
pActivity value (`-log10 M`) and source-document ID:

```csv
compound_id,smiles,p_activity,document_id,relation
CMP-001,CCO,6.3,DOC-A,=
CMP-001,OCC,5.9,DOC-B,=
CMP-002,c1ccccc1,7.2,DOC-A,=
```

```bash
chemgridmap records.csv \
  --compound-column compound_id \
  --smiles-column smiles \
  --value-column p_activity \
  --document-column document_id \
  --relation-column relation \
  --output-dir generic_map
```

Rows are grouped by compound ID and the displayed pActivity is the median.
With a document column, class support is calculated by default with the same
evidence-unit and variance code used for ChEMBL input. `--evidence-model
gaussian` is the default; `exclude-flagged` is the explicit alternative that
removes evidence units containing documents in systematic offset pairs and
falls back to Gaussian support where no unit remains. Use
`--no-merge-identical-documents` to keep re-reported identical values as
separate evidence units.

A curated one-row-per-molecule table remains supported through `--smiles-col`,
`--value-col` and `--label-col`.

## Outputs

A run writes:

- `<name>_grid_coordinates.csv`: molecule annotations, evidence support,
  projection coordinates and grid cells;
- `<name>_metrics.csv`: projection/grid fidelity and neighbor metrics;
- `<name>.svg`, `.png` and `.pdf`: the main activity-colored map;
- optional `<name>_class_support`, `_neighbor_recall` and
  `_structure_group` layers;
- structure-group summaries and member tables;
- an inspector-compatible source-record table for ChEMBL input;
- `<name>_run_manifest.json`: parameters, software versions, diagnostics and
  SHA-256 digests;
- `<name>_results.zip`: the files from one run.

For additional map layers:

```bash
chemgridmap chembl205_ic50.csv --input-format chembl \
  --map-color-modes activity class_support neighbor_recall structure_group \
  --output-dir caii_layers
```

## Paper reproduction

The commands for the main application analyses, E1-A/E1-B reliability analyses,
E2 multi-target evidence evaluation and E2c layout evaluation are in
[`paper/README.md`](paper/README.md). The analysis sources are included under
[`paper/strengthening/`](paper/strengthening/), with target selection rules in
[`paper/strengthening/SELECTION_RULES.md`](paper/strengthening/SELECTION_RULES.md).

## Citation

After the accompanying paper is published, use the preferred citation recorded
in [`CITATION.cff`](CITATION.cff). Until then, cite the ChemGridMap software
version used for the analysis:

> Qiu, Zhenming. ChemGridMap 0.5.0.

## License

ChemGridMap source code is released under the BSD 3-Clause License; see
[`LICENSE`](LICENSE). ChEMBL-derived data are distributed under the ChEMBL
data license, CC BY-SA 3.0. Check the source snapshot and attribution before
redistributing derived tables or maps.
