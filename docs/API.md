# Public API

## `ChemblCurationResult`

```python
ChemblCurationResult(molecules: 'pd.DataFrame', retained_records: 'pd.DataFrame', excluded_records: 'pd.DataFrame', report: 'Dict[str, object]', identity_conflicts: 'pd.DataFrame' = <factory>, structure_variants: 'pd.DataFrame' = <factory>) -> None
```

Tables and audit information produced from a ChEMBL activity export.

## `GridMapResult`

```python
GridMapResult(data: 'pd.DataFrame', representation: 'np.ndarray', projection: 'np.ndarray', grid: 'np.ndarray', metrics: 'pd.DataFrame', preprocessing: 'Dict[str, int]', feature_names: 'List[str]' = <factory>, output_files: 'Dict[str, str]' = <factory>, structure_groups: 'List[Dict]' = <factory>, grouping_protocol: 'Dict' = <factory>) -> None
```

Data and files produced by one grid-map run.

## `annotate_activity_confidence`

```python
annotate_activity_confidence(molecules: 'pd.DataFrame', records: 'pd.DataFrame', *, lower: 'float' = 6.0, upper: 'float' = 7.0, confidence_threshold: 'float' = 0.8, model: 'str' = 'gaussian', fit: 'Optional[Dict[str, object]]' = None, entry: 'str' = 'molecule_identity_key', value: 'str' = 'pchembl_value', document: 'str' = 'document_chembl_id', displayed: 'str' = 'activity_pchembl', merge_identical: 'bool' = True, offset_display: 'str' = 'none', offset_min_shared: 'int' = 3, offset_threshold: 'float' = 1.0, offset_sign_agreement: 'float' = 0.8)
```

Add document-balanced estimates and class confidence to entry-level annotations.

Displayed classes, medians and source records are unchanged. Added columns
describe how strongly the evidence units support the displayed class.

Parameters
----------
molecules : pandas.DataFrame
    One row per molecular entry, containing ``entry`` and ``displayed``.
records : pandas.DataFrame
    Source records with entry, activity-value and document columns.
lower : float, default 6.0
    Upper pActivity bound of the inactive class.
upper : float, default 7.0
    Lower pActivity bound of the active class.
confidence_threshold : float, default 0.8
    Minimum class support used for a supported display.
model : {"gaussian", "exclude-flagged", "moderated_t"}, default "gaussian"
    Evidence model. ``exclude-flagged`` removes units containing documents
    in systematic offset pairs; ``moderated_t`` is available to Python callers.
fit : dict, optional
    Precomputed variance fit to reuse.
entry : str, default "molecule_identity_key"
    Entry identifier column in ``molecules`` and ``records``.
value : str, default "pchembl_value"
    Activity-value column in ``records``.
document : str, default "document_chembl_id"
    Source-document column in ``records``.
displayed : str, default "activity_pchembl"
    Displayed activity column whose class receives support values.
merge_identical : bool, default True
    Merge documents reporting the same entry value into one evidence unit.
offset_display : {"none", "uncertain"}, default "none"
    Keep support-only display or show touched entries as uncertain.
offset_min_shared : int, default 3
    Minimum shared entries for an evaluated document pair.
offset_threshold : float, default 1.0
    Minimum absolute median offset in log units.
offset_sign_agreement : float, default 0.8
    Minimum fraction of paired differences sharing the median sign.

Returns
-------
annotated : pandas.DataFrame
    Copy of ``molecules`` with evidence, support and offset columns.
fit : dict
    Sigma-document estimate, thresholds, diagnostics, evaluated document
    pairs and warnings.

## `assign_to_grid`

```python
assign_to_grid(points_2d: 'np.ndarray', padding: 'Optional[int]' = None, target_occupancy: 'float' = 0.4, coordinate_scaling: 'str' = 'isotropic', method: 'str' = 'auto', dense_max_pairs: 'int' = 20000000, sparse_neighbors: 'int' = 128, diagnostics: 'Optional[Dict]' = None, sparse_adaptive: 'bool' = True, sparse_max_neighbors: 'int' = 1024) -> 'Tuple[np.ndarray, int, int]'
```

Assign each projected point to a unique dense- or sparse-matched cell.

## `build_grid_map`

```python
build_grid_map(data: 'pd.DataFrame', output_dir: 'Optional[Path]' = None, name: 'str' = 'grid_map', smiles_col: 'str' = 'smiles', value_col: 'Optional[str]' = None, label_col: 'Optional[str]' = None, representation: 'str' = 'morgan', projection: 'str' = 'umap', umap_metric: 'str' = 'auto', x_col: 'Optional[str]' = None, y_col: 'Optional[str]' = None, embedding_prefix: 'str' = 'emb_', radius: 'int' = 2, n_bits: 'int' = 2048, lower_threshold: 'float' = 6.0, upper_threshold: 'float' = 7.0, duplicate_policy: 'str' = 'error', random_state: 'int' = 42, k: 'int' = 10, trustworthiness_sample_size: 'int' = 3000, grid_padding: 'Optional[int]' = None, grid_occupancy: 'float' = 0.4, coordinate_scaling: 'str' = 'isotropic', assignment_method: 'str' = 'auto', dense_max_pairs: 'int' = 20000000, sparse_neighbors: 'int' = 128, render_detail: 'str' = 'auto', output_formats: 'Optional[List[str]]' = None, color_map: 'Optional[Dict[str, str]]' = None, provenance: 'Optional[Dict]' = None, sparse_adaptive: 'bool' = True, sparse_max_neighbors: 'int' = 1024, grouping: 'str' = 'butina', group_similarity: 'float' = 0.6, cluster_col: 'Optional[str]' = None, cluster_source: 'Optional[str]' = None, grid_refinement: 'str' = 'neighbors', core_rule: 'str' = 'coverage', core_min_coverage: 'float' = 0.8, fmcs_timeout: 'int' = 1, compound_col: 'Optional[str]' = None, document_col: 'Optional[str]' = None, relation_col: 'Optional[str]' = None, evidence_model: 'str' = 'gaussian', merge_identical_documents: 'bool' = True, confidence_threshold: 'float' = 0.8, offset_display: 'str' = 'none', offset_min_shared: 'int' = 3, offset_threshold: 'float' = 1.0, offset_sign_agreement: 'float' = 0.8, map_color_modes: 'Optional[List[str]]' = None, map_color_formats: 'Optional[Sequence[str]]' = None) -> 'GridMapResult'
```

Build and optionally save a one-molecule-per-cell grid chemical map.

``grid_refinement="neighbors"`` refines the minimum-displacement assignment
so that representation-space neighbours sit closer on the lattice. Every
run reports per-molecule neighbour reliability (``neighbor_recall``).
Extra color layers use ``map_color_formats`` (SVG/PNG); ``None`` keeps SVG
plus PNG when enabled by ``output_formats``. Main-map formats are unchanged.

Parameters
----------
data : pandas.DataFrame
    Molecule table, or generic activity-record table when ``compound_col``
    is supplied.
output_dir : pathlib.Path, optional
    Directory for map files and the run manifest. No files are written when
    omitted.
name : str, default "grid_map"
    Output filename prefix and map name.
smiles_col : str, default "smiles"
    SMILES column consumed by the map workflow.
value_col : str, optional
    Numeric activity column.
label_col : str, optional
    Activity-class label column.
representation : {"morgan", "descriptors", "embedding"}, default "morgan"
    Molecular representation used for projection and neighbour diagnostics.
projection : {"pca", "tsne", "umap"}, default "umap"
    Two-dimensional projection method.
umap_metric : {"auto", "euclidean", "jaccard", "cosine"}, default "auto"
    UMAP distance; auto uses Jaccard for Morgan fingerprints.
x_col, y_col : str, optional
    Existing coordinate columns supplied together instead of a projection.
embedding_prefix : str, default "emb_"
    Prefix of numeric embedding columns.
radius : int, default 2
    Morgan fingerprint radius.
n_bits : int, default 2048
    Morgan fingerprint length.
lower_threshold : float, default 6.0
    Upper pActivity bound of the inactive class.
upper_threshold : float, default 7.0
    Lower pActivity bound of the active class.
duplicate_policy : {"error", "first", "keep"}, default "error"
    Handling of repeated canonical structures.
random_state : int, default 42
    Seed for projection and deterministic tie-breaking.
k : int, default 10
    Neighbour count used for reliability diagnostics.
trustworthiness_sample_size : int, default 3000
    Maximum deterministic sample used for trustworthiness.
grid_padding : int, optional
    Explicit padding added to each grid side.
grid_occupancy : float, default 0.40
    Target fraction of occupied candidate cells when padding is omitted.
coordinate_scaling : {"isotropic", "independent"}, default "isotropic"
    Projection normalization mode.
assignment_method : {"auto", "dense", "sparse"}, default "auto"
    Grid-assignment solver.
dense_max_pairs : int, default 20000000
    Pair-count limit used when choosing dense assignment.
sparse_neighbors : int, default 128
    Initial candidate cells per molecule for sparse assignment.
render_detail : {"auto", "full", "overview"}, default "auto"
    Structure rendering detail.
output_formats : sequence of {"svg", "png", "pdf"}, optional
    Main-map formats.
color_map : dict, optional
    Overrides for activity- or reliability-layer colors.
provenance : dict, optional
    Additional source and parameter metadata for the run manifest.
sparse_adaptive : bool, default True
    Expand sparse candidates when initial candidates are insufficient.
sparse_max_neighbors : int, default 1024
    Candidate cap for adaptive sparse assignment.
grouping : {"butina", "external", "none"}, default "butina"
    Structure-grouping method.
group_similarity : float, default 0.6
    Morgan Tanimoto threshold to a Butina center.
cluster_col, cluster_source : str, optional
    External group labels and their provenance.
grid_refinement : {"neighbors", "none"}, default "neighbors"
    Neighbour-guided refinement mode.
core_rule : {"coverage", "most_frequent", "fmcs"}, default "coverage"
    Group-core selection rule.
core_min_coverage : float, default 0.8
    Coverage target for coverage and FMCS cores.
fmcs_timeout : int, default 1
    FMCS search budget in seconds.
compound_col : str, optional
    Compound identifier for generic activity records.
document_col : str, optional
    Source-document column for generic activity support.
relation_col : str, optional
    Generic relation column; only exact ``=`` records are retained.
evidence_model : {"gaussian", "exclude-flagged"}, default "gaussian"
    Activity-support model for generic records.
merge_identical_documents : bool, default True
    Merge identical entry values across documents as one evidence unit.
confidence_threshold : float, default 0.8
    Minimum class support for a supported display.
offset_display : {"none", "uncertain"}, default "none"
    Keep support-only display or show touched entries as uncertain.
offset_min_shared : int, default 3
    Minimum shared entries for an evaluated document pair.
offset_threshold : float, default 1.0
    Minimum absolute median offset in log units.
offset_sign_agreement : float, default 0.8
    Minimum sign agreement for a systematic offset.
map_color_modes : sequence of str, optional
    Additional reliability layers to export.
map_color_formats : sequence of {"svg", "png"}, optional
    Formats for additional reliability layers.

Returns
-------
GridMapResult
    Prepared data, representation, projection, grid, metrics, groups and
    output-file paths.

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
curate_chembl_activity_data(data: 'pd.DataFrame', target_id: 'Optional[str]' = None, activity_type: 'str' = 'IC50', lower_threshold: 'float' = 6.0, upper_threshold: 'float' = 7.0, conflict_range_threshold: 'float' = 1.0, exclude_potential_duplicates: 'bool' = True, molecule_identity: 'str' = 'auto', structure_policy: 'str' = 'auto', assay_types: 'Optional[List[str]]' = None, identity_conflicts: 'str' = 'error', identity_rule: 'str' = 'inchikey') -> 'ChemblCurationResult'
```

Convert a ChEMBL activity export into an audited molecule-level table.

Non-null pChEMBL values are used as supplied by ChEMBL. When the related
metadata columns are present, the importer also enforces exact relations,
nM standard units, accepted validity flags, the requested target, and the
requested activity type.

Parameters
----------
data : pandas.DataFrame
    Raw ChEMBL activity export.
target_id : str, optional
    Expected Target ChEMBL ID. When omitted, a single target in the input
    is used; multiple targets require an explicit value.
activity_type : str, default "IC50"
    Standard activity type retained.
lower_threshold : float, default 6.0
    Upper pChEMBL bound of the inactive class.
upper_threshold : float, default 7.0
    Lower pChEMBL bound of the active class.
conflict_range_threshold : float, default 1.0
    Activity range above which repeated records are marked as wide.
exclude_potential_duplicates : bool, default True
    Remove records flagged as potential duplicate citations.
molecule_identity : {"auto", "parent-id", "canonical-smiles"}, default "auto"
    Entry-identity policy.
structure_policy : {"auto", "parent", "as-recorded"}, default "auto"
    Parent-structure normalization policy.
assay_types : sequence of str, optional
    Assay-type values retained.
identity_conflicts : {"error", "exclude"}, default "error"
    Stop on conflicting identities or exclude and report them.
identity_rule : {"inchikey", "normalized-structure"}, default "inchikey"
    Identity-conflict rule. ``normalized-structure`` reproduces analyses
    predating version 0.5.

Returns
-------
ChemblCurationResult
    Molecule-level annotations, retained and excluded records, conflicts,
    structure variants and an audit report.

## `document_pair_offsets`

```python
document_pair_offsets(records: 'pd.DataFrame', *, entry: 'str' = 'molecule_identity_key', value: 'str' = 'pchembl_value', document: 'str' = 'document_chembl_id', min_shared: 'int' = 3, offset_threshold: 'float' = 1.0, sign_agreement: 'float' = 0.8) -> 'pd.DataFrame'
```

Diagnose document offsets from the median value of each shared entry.

Documents A/B are sorted by ID; offsets are B minus A. The standard
deviation is the sample SD (zero for a single shared entry). Agreement is
the fraction with the median offset's sign; zero differences do not agree.
Only pairs with at least ``min_shared`` entries are returned. Counts of
all possible, observed (at least one shared entry), and evaluated pairs
are retained in ``attrs``.

## `document_values`

```python
document_values(records: 'pd.DataFrame', *, entry: 'str' = 'molecule_identity_key', value: 'str' = 'pchembl_value', document: 'str' = 'document_chembl_id', merge_identical: 'bool' = True, decimals: 'int' = 2) -> 'pd.DataFrame'
```

Collapse source records to evidence units: one median value per entry and document.

With ``merge_identical``, documents reporting the same entry value (at the
reported precision) form one unit.

## `fit_between_unit_sd`

```python
fit_between_unit_sd(doc: 'pd.DataFrame', *, min_entries: 'int' = 5, fallback_sd: 'Optional[float]' = 0.64) -> 'Dict[str, object]'
```

Pooled SD between evidence units of the same entry (Gaussian model).

## `prepare_generic_activity_records`

```python
prepare_generic_activity_records(data: 'pd.DataFrame', *, compound_column: 'str', value_column: 'str', document_column: 'Optional[str]', smiles_column: 'str', relation_column: 'Optional[str]' = None, lower: 'float' = 6.0, upper: 'float' = 7.0, conflict_range_threshold: 'float' = 1.0, confidence: 'bool' = True, confidence_threshold: 'float' = 0.8, merge_identical: 'bool' = True, evidence_model: 'str' = 'gaussian', offset_display: 'str' = 'none', offset_min_shared: 'int' = 3, offset_threshold: 'float' = 1.0, offset_sign_agreement: 'float' = 0.8) -> 'Tuple[pd.DataFrame, pd.DataFrame, Optional[Dict[str, object]], Dict[str, object]]'
```

Curate a generic activity-record CSV into the molecule-level map schema.

Each compound ID is one entry. Values must already be pActivity values
(-log10 M). If a document column is supplied, class support is calculated
by default with the same evidence functions used by the ChEMBL path.

Parameters
----------
data : pandas.DataFrame
    Generic activity-record table.
compound_column : str
    Compound identifier column.
value_column : str
    pActivity value column.
document_column : str, optional
    Source-document column; required when ``confidence`` is true.
smiles_column : str
    Molecular structure column.
relation_column : str, optional
    Relation column; only exact ``=`` records are retained.
lower : float, default 6.0
    Upper pActivity bound of the inactive class.
upper : float, default 7.0
    Lower pActivity bound of the active class.
conflict_range_threshold : float, default 1.0
    Activity range above which repeated records are marked as wide.
confidence : bool, default True
    Estimate evidence-unit class support.
confidence_threshold : float, default 0.8
    Minimum class support for a supported display.
merge_identical : bool, default True
    Merge documents reporting the same entry value into one evidence unit.
evidence_model : {"gaussian", "exclude-flagged"}, default "gaussian"
    Activity-support model.
offset_display : {"none", "uncertain"}, default "none"
    Keep support-only display or show touched entries as uncertain.
offset_min_shared : int, default 3
    Minimum shared entries for an evaluated document pair.
offset_threshold : float, default 1.0
    Minimum absolute median offset in log units.
offset_sign_agreement : float, default 0.8
    Minimum sign agreement for a systematic offset.

Returns
-------
molecules : pandas.DataFrame
    One curated row per compound.
records : pandas.DataFrame
    Retained source records with normalized activity fields.
fit : dict or None
    Activity-confidence diagnostics, or ``None`` when disabled.
report : dict
    Input counts, excluded records and warnings.

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

Each input row assigns weights to the complete tie-inclusive kth-distance
shell. Strictly closer neighbours receive weight 1, tied boundary
neighbours share the remaining weight, and each row sums to ``k`` (or
``n - 1`` for datasets smaller than ``k + 1``).

## `normalize_coordinates`

```python
normalize_coordinates(points: 'np.ndarray', mode: 'str' = 'isotropic') -> 'np.ndarray'
```

Translate and scale 2D coordinates into a unit-square frame.

Isotropic scaling uses one common factor for both axes, so distances and
angles are preserved up to translation and uniform scale. Independent
scaling is retained for reproducing maps generated before version 0.3.

## `inspect_entry`

```python
inspect_entry(coordinates, records, *, molecule_id=None, cell=None, parameters=None, activity_confidence=None, k=None)
```

Select one entry by exact ID or zero-based (row, column), then verify its records.

Parameters
----------
coordinates : pandas.DataFrame
    Grid-coordinate table from one completed run.
records : pandas.DataFrame
    Retained source-record table from the same run.
molecule_id : str, optional
    Exact molecule identity or source identifier to select.
cell : tuple of int, optional
    Zero-based ``(row, column)`` cell to select.
parameters : dict, optional
    Curation annotation thresholds for verification.
activity_confidence : dict, optional
    Activity-confidence fit from the run manifest or curation report.
k : int, optional
    Neighbour count used to label the kth-neighbour diagnostics.

Returns
-------
selected : pandas.DataFrame
    One-row coordinate subset for the selected entry.
evidence : pandas.DataFrame
    Source records linked to the selected map row.
summary : dict
    Verified annotations, confidence diagnostics, offset pairs and warnings.

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

Candidates meeting the target are ranked by heavy atoms, then coverage;
if none meets it, the highest-coverage candidate is used. Remaining ties
are broken by canonical SMILES.

## `structural_groups`

```python
structural_groups(data, *, similarity=0.6, cluster_col=None, cluster_source=None, radius=2, n_bits=2048, core_rule='most_frequent', core_min_coverage=0.8, fmcs_timeout=1)
```

Butina groups use a center threshold, not an all-pairs similarity bound.

``core_rule="coverage"`` selects the largest member scaffold or scaffold-network
ring fragment that matches at least ``core_min_coverage`` of the members.

Parameters
----------
data : pandas.DataFrame
    Molecule table with canonical SMILES and structure metadata.
similarity : float, default 0.6
    Morgan Tanimoto similarity threshold to a Butina center.
cluster_col : str, optional
    Existing group-label column for external grouping.
cluster_source : str, optional
    Provenance of external group labels.
radius : int, default 2
    Morgan fingerprint radius.
n_bits : int, default 2048
    Morgan fingerprint length.
core_rule : {"coverage", "most_frequent", "fmcs"}, default "most_frequent"
    Group-core selection rule.
core_min_coverage : float, default 0.8
    Coverage target for coverage and FMCS cores.
fmcs_timeout : int, default 1
    FMCS search budget in seconds.

Returns
-------
grouped : pandas.DataFrame
    Copy of ``data`` with group, scaffold and core-member columns.
summaries : list of dict
    Per-group sizes, core choices, coverage and core evidence.
protocol : dict
    Grouping and core-selection parameters and provenance.
