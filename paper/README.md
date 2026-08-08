# Manuscript example

`run_comparisons.py` rebuilds the six configurations reported in the
accompanying manuscript:

- Morgan fingerprint with PCA, t-SNE and UMAP
- Morgan fingerprint, RDKit descriptors and an external embedding with UMAP

The input must be the prepared molecule-level table, with one canonical
molecule per row. The expected columns are:

- `smiles`
- `activity_pchembl`
- optional `activity_class`
- `emb_000`, `emb_001`, ... for the external embedding comparison

The prepared table is included in `examples`. Run:

```bash
python paper/run_comparisons.py \
  examples/chembl205_embedding_umap.csv \
  --output paper/output
```

The original 966-record assay download is not included. This repository begins
from the 770-molecule table used by the visualization method. See
`examples/DATASET.md` for its preparation and licensing notes.
