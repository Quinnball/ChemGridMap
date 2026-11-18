"""来源链接只扩充检查器，不改变实测记录的计算。"""

import io
import json
import zipfile
from email.message import Message
from types import SimpleNamespace

import pandas as pd
import pytest

from chemgridmap.app import AppSession, DEMO, LocalHandler, validate_options
from chemgridmap.inspection import build_sources


LINK_COLUMNS = ["document_doi", "document_pubmed_id", "document_patent_id", "document_url"]
OPTIONS = {
    "input_format": "molecule", "compound_column": "compound_id", "smiles_col": "smiles",
    "value_col": "activity", "document_column": "run_id", "projection": "pca",
    "grouping": "none", "output_formats": ["svg"],
}


def run_session(root, frame, options):
    session = AppSession(root)
    session.load(frame.to_csv(index=False).encode(), "records.csv")
    session.start(options, background=False)
    assert session.job["state"] == "done", session.job
    return session


def reference_frame():
    references = ["doi:10.1021/jm00131a011", "PMID:2810337", "https://example.org/article?q=1&x=2",
                  "javascript:alert(1)", "data:text/html,<h1>bad</h1>", "10.1021/jm00131a011"]
    return pd.DataFrame([
        {"compound_id": f"C{i}", "smiles": smiles, "run_id": f"R{j // 2 + 1}",
         "activity": 5.5 + .2 * i + .03 * j, "reference": reference}
        for i, smiles in enumerate(["CCO", "CCN", "CCC", "CCF", "CCCl"], 1)
        for j, reference in enumerate(references)
    ])


@pytest.fixture(scope="module")
def caii_session(tmp_path_factory):
    session = AppSession(tmp_path_factory.mktemp("u5_caii"))
    session.load(DEMO.read_bytes(), DEMO.name)
    session.start({"projection": "pca", "grouping": "none", "output_formats": ["svg"]},
                  background=False)
    assert session.job["state"] == "done", session.job
    return session


def test_caii_sources_complete_sorted_and_downloads_appended(caii_session):
    session = caii_session
    detail, records = session.inspect(molecule_id="CHEMBL20")
    sources = detail["sources"]
    assert len(sources) == 27
    assert sum(source["n_records"] for source in sources) == detail["total_records"]
    assert sources == sorted(sources, key=lambda source: (
        int(source["year"]) if source["year"] else float("inf"), source["id"]))
    for source in sources:
        links = {link["kind"]: link["url"] for link in source["links"]}
        assert links["chembl"].endswith("/" + source["id"])
        if source["doc_type"] == "PUBLICATION":
            assert links["doi"].startswith("https://doi.org/")
        elif source["doc_type"] == "PATENT":
            assert links["patent"].startswith("https://patents.google.com/patent/")
        if source["id"] == "CHEMBL1909046":
            assert list(links) == ["chembl"]
    original = session.curation.retained_records
    assert records.columns.tolist() == original.columns.tolist() + LINK_COLUMNS
    pd.testing.assert_frame_equal(records[original.columns], original.loc[records.index])
    with zipfile.ZipFile(session.files[detail["audit_file"]]) as archive:
        assert archive.read("source_records.csv") == records.to_csv(index=False, lineterminator="\n").encode()
        assert pd.read_csv(io.BytesIO(archive.read("source_records.csv"))).columns[-4:].tolist() == LINK_COLUMNS
    handler = object.__new__(LocalHandler)
    handler.server = SimpleNamespace(session=session, origin="http://127.0.0.1:8765")
    handler.path = "/api/evidence.csv?id=CHEMBL20"
    handler.headers = Message()
    handler.headers["Host"] = "127.0.0.1:8765"
    handler.headers["X-ChemGridMap-Token"] = session.token
    response = {}
    handler.send = lambda data, content_type="application/json", **kwargs: response.update(
        data=data, content_type=content_type, **kwargs)
    handler.do_GET()
    assert response["content_type"] == "text/csv"
    assert response["data"] == records.to_csv(index=False, lineterminator="\n").encode("utf-8-sig")
    coverage = session.job["report"]["document_coverage"]
    assert coverage["kind"] == "chembl"
    assert coverage["total"] == original[session.curation.report["detected_columns"]["document_id"]].nunique()
    assert coverage["linked"] == coverage["total"]
    saved_report = json.loads(session.files["map_curation_report.json"].read_text(encoding="utf-8"))
    assert "document_coverage" not in saved_report


def test_combined_document_ids_and_laboratory_units(tmp_path):
    frame = reference_frame()
    frame["run_id"] = ["CHEMBL1124628", "R1", "R2", "R3", "CHEMBL1909046", "R1"] * 5
    session = run_session(tmp_path, frame, OPTIONS)
    detail, records = session.inspect(molecule_id="C1")
    sources = {source["id"]: source for source in detail["sources"]}
    assert set(sources) == {"CHEMBL1124628", "CHEMBL1909046", "R1", "R2", "R3"}
    assert sources["CHEMBL1124628"]["links"]
    assert sources["CHEMBL1909046"]["links"]
    assert all(not sources[key]["links"] for key in ("R1", "R2", "R3"))
    assert records.loc[records.run_id.eq("CHEMBL1124628"), "document_doi"].item() == "10.1021/jm00131a011"


def test_reference_option_keeps_calculations_and_groups_values(tmp_path):
    frame = reference_frame()
    plain = run_session(tmp_path / "plain", frame, OPTIONS)
    selected = run_session(tmp_path / "selected", frame, {**OPTIONS, "reference_column": "reference"})
    pd.testing.assert_frame_equal(plain.result.data, selected.result.data)
    assert plain.activity_confidence == selected.activity_confidence
    assert "reference_column" not in json.loads(plain.files["app_options.json"].read_text(encoding="utf-8"))
    assert json.loads(selected.files["app_options.json"].read_text(encoding="utf-8"))["reference_column"] == "reference"
    detail, _ = selected.inspect(molecule_id="C1")
    assert [source["id"] for source in detail["sources"]] == ["R1", "R2", "R3"]
    assert [len(source["references"]) for source in detail["sources"]] == [2, 2, 2]
    assert selected.job["report"]["document_coverage"]["linked"] == 3
    links = [link for source in detail["sources"] for link in source["links"]]
    assert {link["kind"] for link in links} == {"doi", "pubmed", "url"}
    assert all(link["url"].startswith(("http://", "https://")) for link in links)
    assert len(build_sources(frame[frame.compound_id.eq("C1")], None, "reference")) == 6


def test_reference_without_evidence_unit_does_not_create_units(tmp_path):
    options = {**OPTIONS, "document_column": None, "assumed_sigma": .64,
               "reference_column": "reference"}
    session = run_session(tmp_path, reference_frame(), options)
    detail, _ = session.inspect(molecule_id="C1")
    assert len(detail["sources"]) == 6
    assert detail["molecule"]["n_evidence_units"] == 1
    assert session.job["report"]["document_coverage"]["linked"] == 4


def test_numeric_reference_column_preserves_text_and_missing_values(tmp_path):
    frame = reference_frame()
    frame["reference"] = ["0002810337", "2810337", None, "123", None, None] * 5
    session = run_session(tmp_path, frame, {**OPTIONS, "reference_column": "reference"})
    detail, _ = session.inspect(molecule_id="C1")
    values = [reference["value"] for source in detail["sources"] for reference in source["references"]]
    assert "0002810337" in values
    assert "2810337" in values
    assert "2810337.0" not in values


def test_plain_reference_text_and_existing_download_fields_are_preserved(tmp_path):
    frame = reference_frame()
    frame["reference"] = ["NA", "NULL", "N/A", "", "doi:10.1021/jm00131a011", "nan"] * 5
    frame["document_doi"] = "original input value"
    session = run_session(tmp_path, frame, {**OPTIONS, "reference_column": "reference"})
    detail, records = session.inspect(molecule_id="C1")
    references = [item for source in detail["sources"] for item in source["references"]]
    assert {item["value"] for item in references} == {"NA", "NULL", "N/A", "nan", "doi:10.1021/jm00131a011"}
    assert all(not item["url"] for item in references if item["value"] != "doi:10.1021/jm00131a011")
    assert records.columns[-4:].tolist() == LINK_COLUMNS
    assert records.loc[:, "document_doi"].iloc[:, 0].eq("original input value").all()
    assert records.loc[:, "document_doi"].iloc[:, 1].eq("").all()
    assert detail["records"][0]["document_doi"] == "original input value"
    assert records.to_csv(index=False).splitlines()[0].split(",").count("document_doi") == 2


@pytest.mark.parametrize("value", [None, ""])
def test_unselected_reference_is_omitted_from_options(value):
    assert validate_options({**OPTIONS, "reference_column": value}) == validate_options(OPTIONS)


def test_reference_option_requires_generic_records_and_existing_column(tmp_path):
    with pytest.raises(ValueError, match="generic activity records"):
        validate_options({"reference_column": "ref"})
    session = AppSession(tmp_path)
    session.load(reference_frame().to_csv(index=False).encode(), "records.csv")
    with pytest.raises(ValueError, match="reference column was not found"):
        session.start({**OPTIONS, "reference_column": "missing"}, background=False)
