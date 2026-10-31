# Public API

## `ChemblCurationResult`

```python
ChemblCurationResult(molecules: 'pd.DataFrame', retained_records: 'pd.DataFrame', excluded_records: 'pd.DataFrame', report: 'Dict[str, object]', identity_conflicts: 'pd.DataFrame' = <factory>) -> None
```

Tables and audit information produced from a ChEMBL activity export.

## `GridMapResult`

```python
GridMapResult(data: 'pd.DataFrame', representation: 'np.ndarray', projection: 'np.ndarray', grid: 'np.ndarray', metrics: 'pd.DataFrame', preprocessing: 'Dict[str, int]', feature_names: 'List[str]' = <factory>, output_files: 'Dict[str, str]' = <factory>, structure_groups: 'List[Dict]' = <factory>, grouping_protocol: 'Dict' = <factory>) -> None
```

Data and files produced by one grid-map run.

## `annotate_activity_confidence`

```python
annotate_activity_confidence(molecules: 'pd.DataFrame', records: 'pd.DataFrame', *, lower: 'float' = 6.0, upper: 'float' = 7.0, confidence_threshold: 'float' = 0.8, model: 'str' = 'gaussian', fit: 'Optional[Dict[str, object]]' = None, entry: 'str' = 'molecule_identity_key', value: 'str' = 'pchembl_value', document: 'str' = 'document_chembl_id', displayed: 'str' = 'activity_pchembl', merge_identical: 'bool' = True, offset_display: 'str' = 'uncertain', offset_min_shared: 'int' = 3, offset_threshold: 'float' = 1.0, offset_sign_agreement: 'float' = 0.8)
```

Add document-balanced estimates and class confidence to entry-level annotations.

## `assign_to_grid`

```python
assign_to_grid(points_2d: 'np.ndarray', padding: 'Optional[int]' = None, target_occupancy: 'float' = 0.4, coordinate_scaling: 'str' = 'isotropic', method: 'str' = 'auto', dense_max_pairs: 'int' = 20000000, sparse_neighbors: 'int' = 128, diagnostics: 'Optional[Dict]' = None, sparse_adaptive: 'bool' = True, sparse_max_neighbors: 'int' = 1024) -> 'Tuple[np.ndarray, int, int]'
```

Assign each projected point to a unique dense- or sparse-matched cell.

## `build_grid_map`

```python
build_grid_map(data: 'pd.DataFrame', output_dir: 'Optional[Path]' = None, name: 'str' = 'grid_map', smiles_col: 'str' = 'smiles', value_col: 'Optional[str]' = None, label_col: 'Optional[str]' = None, representation: 'str' = 'morgan', projection: 'str' = 'umap', umap_metric: 'str' = 'auto', x_col: 'Optional[str]' = None, y_col: 'Optional[str]' = None, embedding_prefix: 'str' = 'emb_', radius: 'int' = 2, n_bits: 'int' = 2048, lower_threshold: 'float' = 6.0, upper_threshold: 'float' = 7.0, duplicate_policy: 'str' = 'error', random_state: 'int' = 42, k: 'int' = 10, trustworthiness_sample_size: 'int' = 3000, grid_padding: 'Optional[int]' = None, grid_occupancy: 'float' = 0.4, coordinate_scaling: 'str' = 'isotropic', assignment_method: 'str' = 'auto', dense_max_pairs: 'int' = 20000000, sparse_neighbors: 'int' = 128, render_detail: 'str' = 'auto', output_formats: 'Optional[List[str]]' = None, color_map: 'Optional[Dict[str, str]]' = None, provenance: 'Optional[Dict]' = None, sparse_adaptive: 'bool' = True, sparse_max_neighbors: 'int' = 1024, grouping: 'str' = 'butina', group_similarity: 'float' = 0.6, cluster_col: 'Optional[str]' = None, cluster_source: 'Optional[str]' = None, grid_refinement: 'str' = 'neighbors', core_rule: 'str' = 'coverage', core_min_coverage: 'float' = 0.8, fmcs_timeout: 'int' = 1, compound_col: 'Optional[str]' = None, document_col: 'Optional[str]' = None, relation_col: 'Optional[str]' = None, evidence_model: 'str' = 'gaussian', merge_identical_documents: 'bool' = True, confidence_threshold: 'float' = 0.8, offset_display: 'str' = 'uncertain', offset_min_shared: 'int' = 3, offset_threshold: 'float' = 1.0, offset_sign_agreement: 'float' = 0.8, map_color_modes: 'Optional[List[str]]' = None, map_color_formats: 'Optional[Sequence[str]]' = None) -> 'GridMapResult'
```

Build and optionally save a one-molecule-per-cell grid chemical map.

## `canonicalize_smiles`

```python
canonicalize_smiles(smiles: 'object') -> 'Optional[str]'
```

Return an RDKit canonical SMILES, or None for an invalid structure.

## `class_probabilities`

```python
class_probabilities(center, scale, df, lower: 'float', upper: 'float') -> 'np.ndarray'
```

Probabilities of the inactive (<=lower), medium and active (>=upper) classes.

## `curate_chembl_activity_data`

```python
curate_chembl_activity_data(data: 'pd.DataFrame', target_id: 'Optional[str]' = None, activity_type: 'str' = 'IC50', lower_threshold: 'float' = 6.0, upper_threshold: 'float' = 7.0, conflict_range_threshold: 'float' = 1.0, exclude_potential_duplicates: 'bool' = True, molecule_identity: 'str' = 'auto', structure_policy: 'str' = 'auto', assay_types: 'Optional[List[str]]' = None, identity_conflicts: 'str' = 'error') -> 'ChemblCurationResult'
```

Convert a ChEMBL activity export into an audited molecule-level table.

## `document_pair_offsets`

```python
document_pair_offsets(records: 'pd.DataFrame', *, entry: 'str' = 'molecule_identity_key', value: 'str' = 'pchembl_value', document: 'str' = 'document_chembl_id', min_shared: 'int' = 3, offset_threshold: 'float' = 1.0, sign_agreement: 'float' = 0.8) -> 'pd.DataFrame'
```

Diagnose document offsets from the median value of each shared entry.

## `document_values`

```python
document_values(records: 'pd.DataFrame', *, entry: 'str' = 'molecule_identity_key', value: 'str' = 'pchembl_value', document: 'str' = 'document_chembl_id', merge_identical: 'bool' = True, decimals: 'int' = 2) -> 'pd.DataFrame'
```

Collapse source records to evidence units: one median value per entry and document.

## `fit_between_unit_sd`

```python
fit_between_unit_sd(doc: 'pd.DataFrame', *, min_entries: 'int' = 5, fallback_sd: 'Optional[float]' = 0.64) -> 'Dict[str, float]'
```

Pooled SD between evidence units of the same entry (Gaussian model).

## `grid_neighbor_sets`

```python
grid_neighbor_sets(cells: 'np.ndarray', k: 'int' = 10) -> 'List[np.ndarray]'
```

Tie-inclusive k-nearest-neighbour sets on the lattice (Euclidean cell distance).

## `mcs_core`

```python
mcs_core(mols: 'Sequence[Chem.Mol]', *, min_coverage: 'float' = 0.8, timeout: 'int' = 10) -> 'Dict'
```

Maximum common substructure present in at least ``min_coverage`` of members (comparison method).

## `neighbor_recall`

```python
neighbor_recall(reference: 'List[np.ndarray]', displayed: 'List[np.ndarray]') -> 'np.ndarray'
```

Fraction of each molecule's reference neighbours that are also displayed neighbours.

## `neighbor_recall_fixed_k`

```python
neighbor_recall_fixed_k(reference_weights: 'Sequence[Mapping[int, float]]', displayed_weights: 'Sequence[Mapping[int, float]]', *, k: 'int' = 10) -> 'np.ndarray'
```

Expected recall after independently breaking each row's distance ties.

## `normalize_coordinates`

```python
normalize_coordinates(points: 'np.ndarray', mode: 'str' = 'isotropic') -> 'np.ndarray'
```

Translate and scale 2D coordinates into a unit-square frame.

## `inspect_entry`

```python
inspect_entry(coordinates, records, *, molecule_id=None, cell=None, parameters=None, activity_confidence=None, k=None)
```

Select one entry by exact ID or zero-based (row, column), then verify its records.

## `representation_neighbor_sets`

```python
representation_neighbor_sets(matrix: 'np.ndarray', k: 'int' = 10, metric: 'str' = 'jaccard', batch_size: 'int' = 1024) -> 'List[np.ndarray]'
```

Tie-inclusive k-nearest-neighbour sets, excluding self.

## `select_core`

```python
select_core(mols: 'Sequence[Chem.Mol]', *, min_coverage: 'float' = 0.8) -> 'Dict'
```

Largest candidate core matching at least ``min_coverage`` of the members.

## `structural_groups`

```python
structural_groups(data, *, similarity=0.6, cluster_col=None, cluster_source=None, radius=2, n_bits=2048, core_rule='most_frequent', core_min_coverage=0.8, fmcs_timeout=1)
```

Butina groups use a center threshold, not an all-pairs similarity bound.
