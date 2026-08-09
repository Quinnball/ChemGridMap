# Manuscript analyses

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

`run_large_scale_validation.py` retrieves the official ChEMBL REST API records
for CHEMBL240, performs the full curation audit, builds the approximately
10,000-molecule map, and records timing, memory, environment, and fidelity
metadata. The archived input can be reused without network access:

```bash
python paper/run_large_scale_validation.py \
  --skip-download \
  --timing-repeats 5 \
  --hardware-label "Apple M2 Max, 12 cores, 32 GB"
```

`benchmark_sparse_neighbors.py` records sensitivity to the number of sparse
candidate cells. `profile_large_scale_stages.py` profiles the main processing
stages. The 966-record historical example, the compact ChEMBL 37 validation
metadata, and the smaller CHEMBL205 bundle are included under `examples/` and
`validation_data/`. The complete CHEMBL240 bundle is attached to the `v0.2.0`
GitHub release; provenance and licensing notes are documented in the
corresponding README files.
