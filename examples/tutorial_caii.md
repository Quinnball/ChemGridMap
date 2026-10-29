# CA II tutorial: CSV to audited map

This tutorial uses a ChEMBL 37 CA II (`CHEMBL205`) IC50 activity download. It
assumes ChemGridMap 0.5.0 and a CSV containing structures, pChEMBL values,
target/endpoint fields, documents and quality flags.

## 1. Build the map

```bash
chemgridmap chembl205_ic50.csv \
  --input-format chembl \
  --target-id CHEMBL205 \
  --output-dir caii_map \
  --name caii \
  --map-color-modes activity class_support neighbor_recall structure_group
```

The command writes the source audit, coordinate and metric tables, the main map
and the reliability layers. The 0.5.0 defaults are UMAP, neighbor refinement,
coverage cores and Gaussian activity support.

## 2. Read the layers

Open `caii_map/caii.svg` first. Cell color is the observed activity class;
cell opacity indicates class-support status. A faded cell is `uncertain`, not
unmeasured.

Then open:

- `caii_map/caii_neighbor_recall.svg` to check how many fingerprint neighbors
  are also grid neighbors;
- `caii_map/caii_structure_group.svg` to inspect representation-defined groups;
- `caii_map/caii_structure_groups.json` to read core coverage and, when selected,
  FMCS timeout status.

Grid proximity is a layout property. It is not proof of structural or
biological similarity.

## 3. Open one molecule's source records

For the packaged example, CHEMBL20 has 54 retained records and median pChEMBL
7.185:

```bash
chemgridmap inspect caii_map/caii_grid_coordinates.csv \
  --records caii_map/caii_retained_records.csv \
  --molecule-id CHEMBL20 \
  --output-dir caii_record_audit
```

Read `caii_record_audit/audit_summary.json` for the verified counts, median,
annotation rules, document warnings and run binding. The source rows are in
`source_records.csv`.

## 4. Use the browser

```bash
chemgridmap-app
```

Drop the same CSV into the page, keep the 0.5.0 defaults or adjust them under
Advanced settings, generate the map, click a cell, and open the source-record
audit. The browser uses the same curation, evidence and grid code as the CLI.
