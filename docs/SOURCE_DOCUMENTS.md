# Offline source-document helpers

```python
from chemgridmap.inspection import (
    load_documents, normalize_doi, reference_link, build_sources,
    document_coverage, enrich_source_records,
    build_source_index,
)

load_documents() -> dict
normalize_doi(value) -> str
reference_link(value) -> dict | None
build_sources(records, document_column, reference_column=None) -> list[dict]
document_coverage(sources) -> dict
enrich_source_records(records, document_column=None) -> pandas.DataFrame
build_source_index(molecules, records, document_column, reference_column=None) -> dict
```

`load_documents` lazily reads the bundled ChEMBL 37 document table and caches
the mapping from document ID to CSV fields. Treat the returned table as
read-only. The source metadata, download time, SHA-256 and CC BY-SA 3.0 license
are in `data/chembl37_documents.json`; no network request is made.

`normalize_doi` strips whitespace and DOI resolver/`doi:` prefixes while
preserving case. `reference_link` recognizes DOI values, `PMID:` or one to
nine digit PubMed IDs, and HTTP(S) URLs. It returns plain-text `kind`, `label`
and `url` strings, or `None` for other text and unsupported protocols.

`build_sources` groups retained records by the evidence-unit column, or by
nonempty reference values when `document_column=None`. Each source contains
`id`, `doc_type`, `journal`, `year`, `n_records`, `links` and `references`.
Links include DOI, PubMed, Google Patents and the ChEMBL document page where
available, plus recognized user references. References retain their supplied
`value`, including plain text without a link. Results sort by ascending year
then ID, with missing years last. Neither source records nor calculations are
modified. Renderers must escape displayed text and use `target="_blank"` and
`rel="noopener noreferrer"` when displaying the links.

`document_coverage` returns `total`, `linked`, `unlinked`, `doi_or_pubmed`,
`patent` and `chembl_only` counts. Patent sources take precedence over DOI or
PubMed sources so the document categories do not overlap. `linked` also
counts safe user reference URLs.

`enrich_source_records` returns a copy with `document_doi`,
`document_pubmed_id`, `document_patent_id` and `document_url` appended in that
order. Unknown document IDs receive empty values. If no document column is
specified, `document_chembl_id` is used when present. The helper is intended
for inspector CSV and evidence-ZIP downloads, separate from map output files.
Existing same-named input fields stay in their original positions; the appended
fields then have repeated CSV headings. The on-screen record table keeps its original schema.

## Highlight molecules from retained sources

`build_source_index` groups the **retained records** with exactly the same
column choice and `_reference_text(value).strip()` normalization as
`build_sources`: ChEMBL document ID, a selected generic document/batch column,
or the reference-column fallback. A molecule belongs to a source if at least
one retained record has that identifier. Counts include only distinct mapped
reference molecules; excluded records never contribute. The source metadata
(`doc_type`, journal, year and links) comes from `build_sources`.

The browser builds this index lazily once per reference map. Its
`GET /api/source-members?source=ID` response contains `sources`,
`map_molecule_count` and `members` with grid coordinates and identifiers.
Repeated `source` parameters form a union. Unknown sources or a missing map
return HTTP 400. Inspector sources gain `map_molecule_count`, and the inspector
also supplies `source_union_molecule_count` for the shared-source button.

Source highlighting dims other molecules and uses a dark dashed outline in
the grid and projection. It replaces structure-group highlighting, preserves
the selected-cell outline and provides an explicit clear button. Rebuilding,
inserting or removing new molecules clears it.

For inserted molecules, `scope=new` uses only the current inserted population.
Measured batches use retained measurement records; unmeasured batches use the
selected uploaded batch/document column. Duplicate structures count once and
reference matches do not become new molecules. The new-molecule inspector
exposes `new_molecule_sources`; replacing/removing the upload invalidates that
separate index. No original export is rewritten by either kind of highlight.

The possible-source lookup described below identifies structural matches,
not retained-record provenance, and therefore never supplies highlight buttons
or contributes to this index.

## Optional source lookup file

For generic activity records without linked sources, the browser accepts a
separate **Source lookup: ChEMBL export (optional)** file. Choose an activity
export for the same target. It must contain recognizable SMILES and ChEMBL
document ID columns. Matching uses rows with valid structures and document IDs
of the form `CHEMBL` followed by digits. No extra dataset is bundled or fetched.

Both input structures pass through `canonicalize_smiles` and
`parent_structure(...)[0]`; each distinct SMILES is standardized once per map
build. Salts match their parent structures, but stereoisomers remain distinct.
The result lists every matching document, using the existing bundled document
metadata and safe DOI, PubMed, patent and ChEMBL links. **Same value** marks a
document only when a lookup record's `pchembl_value` is within 0.005 of one of
the molecule's retained activity records after conversion to pChEMBL. Without a
pChEMBL column, no document receives this mark. Marked documents sort first,
then by ascending year and document ID, with missing years last.

**Possible sources** means that the documents measured the same structure;
matching does not establish where the user's measurements originated. Target
and activity type do not filter matches. Data audit warns when the lookup file
contains more than one target, because possible sources may concern other
targets. The inspector shows this section only for measured molecules whose
own sources have no links, including a message when no structure matches.

Matching runs once during map generation and usually takes a few seconds for a
single-target export. It runs locally without uploading data. Lookup results,
counts and the file name remain in the browser session; neither the lookup file
nor its name is added to map exports, manifests or `app_options.json`. Loading
or removing a lookup file discards the current map. Changing measured input
clears the lookup file. ChEMBL input ignores the optional file and hides its
control, because its own document IDs already identify its sources.
