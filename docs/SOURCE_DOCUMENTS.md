# Offline source-document helpers

```python
from chemgridmap.inspection import (
    load_documents, normalize_doi, reference_link, build_sources,
    document_coverage, enrich_source_records,
)

load_documents() -> dict
normalize_doi(value) -> str
reference_link(value) -> dict | None
build_sources(records, document_column, reference_column=None) -> list[dict]
document_coverage(sources) -> dict
enrich_source_records(records, document_column=None) -> pandas.DataFrame
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

