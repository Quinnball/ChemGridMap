"""来源索引只使用保留记录，且按已放置的不同分子计数。"""

import hashlib

import pandas as pd
import pytest

from chemgridmap.app import AppSession, DEMO
from chemgridmap.inspection import build_source_index, build_sources
from test_app_u8_session import api


OPTIONS = {"input_format": "molecule", "compound_column": "id", "smiles_col": "smiles",
           "value_col": "value", "document_column": "batch", "relation_column": "relation",
           "projection": "pca", "grouping": "none", "grid_refinement": "none",
           "output_formats": ["svg"]}
RAW = (b"id,smiles,value,batch,relation,reference\n"
       b"A,CCO,6, R1 ,=,00123\nA,CCO,6.1,R2,=,00123\n"
       b"B,CCN,7,R1,=,00123\nB,CCN,8,R1,=,00123\n"
       b"C,CCC,5,R3,=,00456\nC,CCC,7,R1,<,00123\n"
       b"D,c1ccccc1,8,R2,=,00456\n")


def build_session(root, options=None, raw=RAW):
    session = AppSession(root)
    session.load(raw, "records.csv")
    session.start(options or OPTIONS, background=False)
    assert session.job["state"] == "done", session.job
    return session


def test_source_index_normalizes_like_build_sources_and_counts_distinct_map_entries():
    molecules = pd.DataFrame({"molecule_identity_key": ["A", "B"], "molecule_id": ["A", "B"],
                              "grid_row": [0, 1], "grid_col": [2, 3]})
    records = pd.DataFrame({"molecule_identity_key": ["A", "A", "B", "X", "A"],
                            "source": [123., " 123 ", "123", "other", None]})
    for document, reference in (("source", None), (None, "source")):
        index = build_source_index(molecules, records, document, reference)
        assert list(index) == ["123"]
        assert index["123"]["map_molecule_count"] == 2
        expected = next(source for source in build_sources(records, document, reference)
                        if source["id"] == "123")
        assert {key: index["123"][key] for key in expected} == expected
        assert index["123"]["members"] == molecules.to_dict(orient="records")


def test_generic_sources_retained_union_cached_and_inspection_counts(tmp_path, monkeypatch):
    session = build_session(tmp_path)
    hashes = {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in session.files.values()}
    result = api(session, "/api/source-members?source=R1&source=R2&source=R1")
    assert result["status"] == 200
    payload = result["data"]
    assert payload["map_molecule_count"] == 3
    assert {row["molecule_id"] for row in payload["members"]} == {"A", "B", "D"}
    assert [source["map_molecule_count"] for source in payload["sources"]] == [2, 2]
    assert hashes == {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in session.files.values()}
    monkeypatch.setattr("chemgridmap.app.build_source_index", lambda *args: pytest.fail("来源索引重复构建"))
    detail, _ = session.inspect(molecule_id="A")
    assert detail["source_union_molecule_count"] == 3
    assert [source["map_molecule_count"] for source in detail["sources"]] == [2, 2]
    assert session.source_members(["R3"])["map_molecule_count"] == 1


def test_reference_fallback_uses_original_text(tmp_path):
    session = build_session(tmp_path, {**OPTIONS, "document_column": None, "reference_column": "reference"})
    detail, _ = session.inspect(molecule_id="A")
    assert detail["sources"][0]["id"] == "00123"
    assert detail["sources"][0]["map_molecule_count"] == 2
    assert session.source_members(["00456"])["map_molecule_count"] == 2


def test_source_endpoint_errors_and_authentication(tmp_path):
    empty = AppSession(tmp_path / "empty")
    response = api(empty, "/api/source-members?source=R1")
    assert response["status"] == 400 and response["data"]["error"] == "Build a map first."
    session = build_session(tmp_path / "map")
    for path in ("/api/source-members?source=missing", "/api/source-members",
                 "/api/source-members?source=R1&scope=invalid"):
        assert api(session, path)["status"] == 400
    assert api(session, "/api/source-members?source=R1", token=False)["status"] == 403


@pytest.mark.parametrize("measured", [False, True])
def test_new_batch_counts_only_new_molecules_and_invalidates_on_replacement(tmp_path, measured):
    session = build_session(tmp_path)
    raw = (b"id,smiles,value,batch,relation\nN1,CCF,6,R1,=\nN2,CCCl,7,R1,=\n"
           b"CTRL,CCO,6,R1,=\nBAD,bad_smiles,8,R1,=\nN3,CCBr,8,R2,<\n")
    session.load_query(raw, "new.csv")
    options = {"smiles_col": "smiles", "id_col": "id", "document_col": "batch"}
    if measured:
        options.update(activity_col="value", relation_col="relation")
    session.insert_new(options)
    detail, _ = session.inspect(molecule_id="N1")
    assert [(source["id"], source["map_molecule_count"]) for source in detail["new_molecule_sources"]] == [("R1", 2)]
    payload = api(session, "/api/source-members?source=R1&scope=new")["data"]
    assert payload["scope"] == "new"
    assert {row["molecule_identity_key"] for row in payload["members"]} == {"query:N1", "query:N2"}
    assert session.source_members(["R1"])["map_molecule_count"] == 2
    assert api(session, "/api/source-members?source=R2&scope=new")["status"] == (400 if measured else 200)
    session.load_query(b"id,smiles,batch\nN4,CCF,R1\n", "replacement.csv")
    session.insert_new({"smiles_col": "smiles", "id_col": "id", "document_col": "batch"})
    assert session.source_members(["R1"], "new")["map_molecule_count"] == 1
    session.clear_query()
    assert api(session, "/api/source-members?source=R1&scope=new")["status"] == 400
    session.start(OPTIONS, background=False)
    assert session.source_members(["R1"])["map_molecule_count"] == 2


def test_caii_chembl18_source_counts_match_independent_retained_record_join(tmp_path):
    session = build_session(tmp_path, {"projection": "pca", "grouping": "none",
                                     "output_formats": ["svg"]}, DEMO.read_bytes())
    records = session.curation.retained_records
    document = session.curation.report["detected_columns"]["document_id"]
    detail, _ = session.inspect(molecule_id="CHEMBL18")
    mapped = set(session.result.data.molecule_identity_key.astype(str))
    union = set()
    for source in detail["sources"]:
        expected = set(records.loc[records[document].astype(str).str.strip().eq(source["id"]),
                                   "molecule_identity_key"].astype(str)) & mapped
        actual = session.source_members([source["id"]])
        assert source["map_molecule_count"] == actual["map_molecule_count"] == len(expected)
        assert {member["molecule_identity_key"] for member in actual["members"]} == expected
        union |= expected
    assert detail["source_union_molecule_count"] == len(union)
    assert session.source_members([source["id"] for source in detail["sources"]])["map_molecule_count"] == len(union)
