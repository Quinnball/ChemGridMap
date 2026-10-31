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
coverage cores, Gaussian activity support, a 0.8 support level and the
offset-aware display rule.

## 2. Read the layers

Open `caii_map/caii.svg` first. Cell color is the observed activity class;
the main map is the activity layer. Class-support status is shown in the
separate `caii_class_support` layer. There, a faded cell is `uncertain`, not
unmeasured. Pale cells can be below the support threshold or touched by an
offset document pair; a bottom-left marker identifies the latter. Full-detail
cells use an opaque colored frame, a white inset and the molecular structure.
The same files at overview scale show solid cell colors.

Then open:

- `caii_map/caii_neighbor_recall.svg` to check how many fingerprint neighbors
  are also grid neighbors; the coordinates CSV also contains
  `neighbor_recall_fixed_k`, which breaks distance ties at random with total
  weight k, and `neighbor_kth_similarity` for the kth fingerprint neighbor;
- `caii_map/caii_structure_group.svg` to inspect representation-defined groups;
- `caii_map/caii_structure_groups.json` to read core coverage and, when selected,
  FMCS timeout status.

Grid proximity is a layout property. It is not proof of structural or
biological similarity.

Support level, offset display, minimum shared entries, minimum log-unit
difference and sign agreement can be changed from the CLI or browser. The fit
and run manifest record all five values. Set offset display to `none` to keep
the support-only status while retaining the offset marker and diagnostics.

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
annotation rules, offset document pairs and their touched-entry counts, the
configured support/offset rules, σ_doc value and degrees of freedom, the kth
neighbor value and run binding. The source rows are in `source_records.csv`.

## 4. Use the browser

```bash
chemgridmap-app
```

Drop the same CSV into the page, keep the 0.5.0 defaults or adjust them under
Advanced settings, generate the map, click a cell, and open the source-record
audit. The browser uses the same curation, evidence and grid code as the CLI.

See [`docs/OUTPUT_COLUMNS.md`](../docs/OUTPUT_COLUMNS.md) for every output
column and [`docs/API.md`](../docs/API.md) for the public Python API.
