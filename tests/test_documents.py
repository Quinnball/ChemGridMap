"""Offline document lookup, safe links and source grouping."""
import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path

import pandas as pd
import pytest

from chemgridmap import curate_chembl_activity_data
from chemgridmap.inspection import (
    build_sources, document_coverage, enrich_source_records, load_documents,
    normalize_doi, reference_link,
)


DATA = Path(__file__).resolve().parents[1] / "src" / "chemgridmap" / "data"
DOCUMENT_SHA256 = "436bbdf1ce80433ebe382adcba1f1456912a50ac43747ec74f2546286d58eaef"


@pytest.mark.parametrize("value", [
    "10.1021/Jm00131a011", " doi:10.1021/Jm00131a011 ",
    "DOI: 10.1021/Jm00131a011", "https://doi.org/10.1021/Jm00131a011",
    "http://dx.doi.org/10.1021/Jm00131a011",
])
def test_doi_prefixes_keep_case(value):
    assert normalize_doi(value) == "10.1021/Jm00131a011"
    assert reference_link(value) == {
        "kind": "doi", "label": "DOI 10.1021/Jm00131a011",
        "url": "https://doi.org/10.1021/Jm00131a011",
    }


def test_doi_url_escapes_suffix_and_preserves_slashes():
    assert reference_link('10.1234/A B/C?D#E%"<x>')["url"] == (
        "https://doi.org/10.1234/A%20B/C%3FD%23E%25%22%3Cx%3E")


@pytest.mark.parametrize("value", ["PMID:2810337", "pmid: 2810337", "2810337", 2810337, 2810337.0])
def test_pubmed_references(value):
    assert reference_link(value) == {
        "kind": "pubmed", "label": "PubMed 2810337",
        "url": "https://pubmed.ncbi.nlm.nih.gov/2810337/",
    }


@pytest.mark.parametrize("value", [
    "javascript:alert(1)", "data:text/html,<script>alert(1)</script>",
    "file:///tmp/paper.pdf", "ftp://example.org/paper", "//example.org/paper",
    "1234567890", "plain reference", "", None,
])
def test_unsafe_or_unknown_references_do_not_become_links(value):
    assert reference_link(value) is None


@pytest.mark.parametrize("value", ["https://example.org/paper?a=1&b=2", "http://example.org/paper"])
def test_http_references(value):
    assert reference_link(value) == {"kind": "url", "label": value, "url": value}


def test_explicit_pmid_prefix_accepts_a_numeric_identifier():
    assert reference_link("PMID:1234567890")["url"] == "https://pubmed.ncbi.nlm.nih.gov/1234567890/"


def test_bundled_documents_have_exact_provenance_and_full_unique_inventory():
    path = DATA / "chembl37_documents.csv.gz"
    assert path.stat().st_size == 997982
    assert hashlib.sha256(path.read_bytes()).hexdigest() == DOCUMENT_SHA256
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 101100
    assert len({row["document_chembl_id"] for row in rows}) == len(rows)
    assert Counter(row["doc_type"] for row in rows) == {
        "PUBLICATION": 93488, "PATENT": 7145, "DATASET": 465, "BOOK": 2,
    }
    assert sum(bool(row["doi"]) for row in rows) == 93787
    assert sum(bool(row["pubmed_id"]) for row in rows) == 88134
    assert sum(bool(row["patent_id"]) for row in rows) == 7145
    assert sum(not any(row[field] for field in ("doi", "pubmed_id", "patent_id"))
               for row in rows) == 115
    metadata = json.loads((DATA / "chembl37_documents.json").read_text(encoding="utf-8"))
    assert metadata["sha256"] == DOCUMENT_SHA256
    assert metadata["rows"] == len(rows)
    assert metadata["chembl_db_version"] == "ChEMBL_37"
    assert "CC BY-SA 3.0" in metadata["license"]
    document = load_documents()["CHEMBL1124628"]
    assert document["doi"] == "10.1021/jm00131a011"
    assert document["pubmed_id"] == "2810337"


@pytest.fixture(scope="module")
def caii_records():
    curated = curate_chembl_activity_data(
        pd.read_csv(DATA / "chembl205_ic50.csv"), target_id="CHEMBL205", activity_type="IC50")
    return curated.retained_records


def test_caii_chembl20_has_all_27_sources(caii_records):
    selected = caii_records[caii_records.molecule_identity_key.eq("CHEMBL20")]
    sources = build_sources(selected, "document_chembl_id")
    assert len(sources) == 27
    assert sum(source["n_records"] for source in sources) == len(selected)
    publications = [source for source in sources if source["doc_type"] == "PUBLICATION"]
    assert publications
    assert all(any(link["kind"] == "doi" for link in source["links"])
               for source in publications)
    drugmatrix = next(source for source in sources if source["id"] == "CHEMBL1909046")
    assert drugmatrix["links"] == [{
        "kind": "chembl", "label": "ChEMBL CHEMBL1909046",
        "url": "https://www.ebi.ac.uk/chembl/explore/document/CHEMBL1909046",
    }]
    assert sources[-1] == drugmatrix
    ordered = [(int(source["year"]), source["id"]) for source in sources if source["year"]]
    assert ordered == sorted(ordered)


def test_caii_patents_link_to_google_patents(caii_records):
    sources = build_sources(caii_records, "document_chembl_id")
    patents = [source for source in sources if source["doc_type"] == "PATENT"]
    assert patents
    assert all(any(link["url"].startswith("https://patents.google.com/patent/")
                   for link in source["links"]) for source in patents)


def test_patent_number_loses_hyphens_only():
    document = next(row for row in load_documents().values()
                    if row["patent_id"] == "US-10174040-B2")
    sources = build_sources(pd.DataFrame({"doc": [document["document_chembl_id"]]}), "doc")
    link = next(link for link in sources[0]["links"] if link["kind"] == "patent")
    assert link["url"] == "https://patents.google.com/patent/US10174040B2"
    assert link["label"] == "Patent US-10174040-B2"


def test_combined_document_and_lab_evidence_preserves_unlinked_units():
    records = pd.DataFrame({"run_id": ["R3", "CHEMBL1124628", "R1", "R2", "R1"]})
    sources = build_sources(records, "run_id")
    assert [source["id"] for source in sources] == ["CHEMBL1124628", "R1", "R2", "R3"]
    assert [source["n_records"] for source in sources] == [1, 2, 1, 1]
    assert sources[0]["links"]
    assert all(not source["links"] for source in sources[1:])
    assert document_coverage(sources) == {
        "total": 4, "doi_or_pubmed": 1, "patent": 0, "chembl_only": 0,
        "linked": 1, "unlinked": 3,
    }


def test_reference_values_group_by_evidence_unit_and_preserve_plain_text():
    invalid = ' javascript:alert("paper") '
    records = pd.DataFrame({
        "run_id": ["R1", "R1", "R1", "R1", "R2", "R2", "R3"],
        "reference": ["10.1234/Paper", "10.1234/Paper", "PMID:2810337", invalid,
                      "https://example.org/paper", "data:text/plain,paper", "no reference"],
    })
    original = records.copy(deep=True)
    sources = build_sources(records, "run_id", "reference")
    assert len(sources) == 3
    assert [source["n_records"] for source in sources] == [4, 2, 1]
    assert [reference["value"] for reference in sources[0]["references"]] == [
        "10.1234/Paper", "PMID:2810337", invalid]
    assert [link["kind"] for link in sources[0]["links"]] == ["doi", "pubmed"]
    assert sources[0]["references"][-1]["url"] == ""
    assert [link["kind"] for link in sources[1]["links"]] == ["url"]
    assert document_coverage(sources)["linked"] == 2
    pd.testing.assert_frame_equal(records, original)


def test_without_evidence_column_sources_group_by_reference():
    records = pd.DataFrame({"reference": ["10.1234/Paper", "10.1234/Paper", "2810337",
                                          "javascript:alert(1)", "", None]})
    sources = build_sources(records, None, "reference")
    assert [source["id"] for source in sources] == ["10.1234/Paper", "2810337", "javascript:alert(1)"]
    assert [source["n_records"] for source in sources] == [2, 1, 1]
    assert document_coverage(sources)["linked"] == 2
    assert build_sources(records, None) == []


def test_document_links_and_equivalent_reference_urls_are_not_duplicated():
    records = pd.DataFrame({"doc": ["CHEMBL1124628"], "reference": ["doi:10.1021/jm00131a011"]})
    source = build_sources(records, "doc", "reference")[0]
    assert len(source["links"]) == 3
    assert len(source["references"]) == 1


def test_patent_coverage_categories_do_not_overlap():
    document = next(row for row in load_documents().values()
                    if row["patent_id"] and row["doi"])
    sources = build_sources(pd.DataFrame({"doc": [document["document_chembl_id"], "CHEMBL1909046"]}), "doc")
    coverage = document_coverage(sources)
    assert coverage == {"total": 2, "doi_or_pubmed": 0, "patent": 1,
                        "chembl_only": 1, "linked": 2, "unlinked": 0}


def test_source_record_download_preserves_original_columns_and_appends_document_fields():
    records = pd.DataFrame({"source_row": [9, 12], "run_id": ["CHEMBL1124628", "R1"],
                            "pActivity": [6.2, 7.1], "notes": ["first", "second"]})
    result = enrich_source_records(records, "run_id")
    assert list(result.columns) == list(records.columns) + [
        "document_doi", "document_pubmed_id", "document_patent_id", "document_url"]
    pd.testing.assert_frame_equal(result[records.columns], records)
    assert result.loc[0, "document_doi"] == "10.1021/jm00131a011"
    assert result.loc[0, "document_pubmed_id"] == "2810337"
    assert result.loc[0, "document_url"] == "https://www.ebi.ac.uk/chembl/explore/document/CHEMBL1124628"
    assert result.iloc[1, -4:].tolist() == ["", "", "", ""]
    assert len(records.columns) == 4
