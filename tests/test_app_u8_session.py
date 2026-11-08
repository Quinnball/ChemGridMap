"""固定参照图上的新分子会话、HTTP 数据契约与独立文件。"""

from email.message import Message
import hashlib
import io
import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from chemgridmap.app import AppSession, LocalHandler


REFERENCE = b"molecule_id,smiles,value\nM1,CCO,5\nM2,CCN,7.2\nM3,CCC,6.5\nM4,c1ccccc1,8\n"
OPTIONS = {"input_format": "molecule", "value_col": "value", "projection": "pca",
           "grouping": "none", "grid_refinement": "none", "render_detail": "overview",
           "output_formats": ["svg"]}
QUERIES = b"id,smiles,label\nN1,CCF,round1\nN2,CCCl,round2\nREF,OCC,control\nDUP,FCC,round1\nBAD,bad_smiles,invalid\n"


def api(session, path, content=b"", *, method="GET", token=True):
    handler = object.__new__(LocalHandler)
    handler.server = SimpleNamespace(session=session, origin="http://127.0.0.1:8765")
    handler.path = path
    handler.headers = Message()
    handler.headers["Host"] = "127.0.0.1:8765"
    handler.headers["Content-Length"] = str(len(content))
    handler.headers["X-Filename"] = "new.tsv"
    if token:
        handler.headers["X-ChemGridMap-Token"] = session.token
    handler.rfile = io.BytesIO(content)
    response = {}

    def send(data, content_type="application/json", status=200, **kwargs):
        response.update(status=status, data=data if isinstance(data, bytes) else json.loads(json.dumps(data, allow_nan=False)),
                        content_type=content_type, **kwargs)

    handler.send = send
    getattr(handler, "do_" + method)()
    return response


@pytest.fixture
def session(tmp_path):
    app = AppSession(tmp_path)
    app.load(REFERENCE, "reference.csv")
    app.start(OPTIONS, background=False)
    assert app.job["state"] == "done", app.job
    return app


def insert(session, raw=QUERIES, name="batch.csv", **options):
    session.load_query(raw, name)
    return session.insert_new({"smiles_col": "smiles", "id_col": "id", **options})


def test_insertion_keeps_every_reference_file_and_does_not_rebuild(session, monkeypatch):
    files = {name: path.read_bytes() for name, path in session.files.items()}
    original = session.result.data.copy(deep=True)
    original_representation = session.result.representation.copy()
    monkeypatch.setattr("chemgridmap.app.build_grid_map", lambda *args, **kwargs: pytest.fail("重新建图"))
    insert(session, label_col="label")
    pd.testing.assert_frame_equal(session.result.data, original)
    assert (session.result.representation == original_representation).all()
    assert all(session.files[name].read_bytes() == content for name, content in files.items())
    assert len(session.points()) == len(original) + 2
    assert session.new_result.report["neighbor_evidence_counts"]
    assert not session.new_result.report["reference_support_assessed"]
    for payload in (session.snapshot(), session.new_molecules_info(), session.points()):
        json.dumps(payload, allow_nan=False)


def test_replace_remove_refresh_rebuild_and_reference_replacement(session):
    insert(session, label_col="label")
    old_path = session.files["new_molecules_input.csv"]
    old_hash = hashlib.sha256(old_path.read_bytes()).hexdigest()
    restored = api(session, "/api/session")["data"]["new_molecules"]
    assert restored["input"]["filename"] == "batch.csv"
    assert restored["options"]["label_column"] == "label"
    assert len(restored["rows"]) == 5
    insert(session, b"id,smiles\nN3,CCBr\n", name="batch.csv")
    assert session.files["new_molecules_input.csv"].parent.name == "batch-2"
    assert hashlib.sha256(old_path.read_bytes()).hexdigest() == old_hash
    assert session.new_result.data.query_id.tolist() == ["N3"]
    old_run = session.job["output_dir"]
    session.start(OPTIONS, background=False)
    assert session.job["state"] == "done", session.job
    assert session.job["output_dir"] != old_run
    assert session.new_result.data.query_id.tolist() == ["N3"]
    assert session.files["new_molecules_input.csv"].read_bytes() == b"id,smiles\nN3,CCBr\n"
    current_path = session.files["new_molecules_input.csv"]
    session.clear_query()
    assert current_path.is_file() and old_path.is_file()
    assert session.new_result is None and session.new_molecules_info() is None
    assert session.result is not None and len(session.points()) == 4
    assert "query_molecules" not in session.job["report"]
    assert not any(name.startswith("new_molecules") for name in session.files)
    insert(session)
    session.load(REFERENCE, "another-reference.csv")
    assert session.result is None and session.new_result is None
    assert session.query_frame is None and session.pending_query is None


def test_upload_is_staged_until_successful_replacement(session):
    insert(session)
    saved = session.new_molecules_info()
    session.load_query(b"id,structure\nnew,CCBr\n", "bad-columns.csv")
    with pytest.raises(ValueError, match="column"):
        session.insert_new({"smiles_col": "missing"})
    assert session.new_molecules_info() == saved
    assert session.query_filename == "batch.csv"
    assert session.new_result.data.query_id.tolist() == ["N1", "N2"]


def test_http_new_molecule_workflow_and_strict_json(session):
    response = api(session, "/api/new-molecules/input", QUERIES, method="POST")
    assert response["status"] == 200 and response["data"]["rows"] == 5
    payload = {"smiles_col": "smiles", "id_col": "id", "label_col": "label"}
    assert api(session, "/api/new-molecules", json.dumps(payload).encode(), method="POST")["status"] == 200
    rows = api(session, "/api/new-molecules")["data"]["rows"]
    assert {row["status"] for row in rows} == {"mapped", "in_reference", "duplicate_query", "invalid_smiles"}
    assert api(session, "/api/points")["status"] == 200
    assert api(session, "/api/inspect?id=N1")["status"] == 200
    assert api(session, "/api/inspect?id=REF")["status"] == 200
    assert api(session, "/api/new-molecules/clear", method="POST")["status"] == 200
    assert api(session, "/api/new-molecules")["data"]["rows"] == []


@pytest.mark.parametrize("path", ["/api/new-molecules/input", "/api/new-molecules", "/api/new-molecules/clear"])
def test_new_endpoints_keep_session_authentication(session, path):
    assert api(session, path, b"{}", method="POST", token=False)["status"] == 403


def test_new_input_requires_existing_morgan_map(session, tmp_path):
    empty = AppSession(tmp_path / "empty")
    empty.load(REFERENCE, "ref.csv")
    assert api(empty, "/api/new-molecules/input", QUERIES, method="POST")["status"] == 400
    session.job["options"]["representation"] = "descriptors"
    response = api(session, "/api/new-molecules/input", QUERIES, method="POST")
    assert response["status"] == 400 and "Morgan" in response["data"]["error"]


def test_own_measurements_records_units_reference_difference_and_inspection(session):
    raw = ("id,smiles,value,unit,relation,batch,label\n"
           "N1,CCF,100,nM,=,R1,lab\nN1,FCC,0.1,uM,=,R2,lab\n"
           "N1,CCF,0.1,µM,=,R3,lab\nN1,CCF,1000,nM,<,R1,lab\n"
           "CTRL,OCC,1000,nM,=,R1,control\nCTRL,CCO,1,uM,=,R2,control\n").encode()
    insert(session, raw, activity_col="value", unit_col="unit", relation_col="relation",
           document_col="batch", label_col="label")
    detail, records = session.inspect(molecule_id="N1")
    assert detail["molecule"]["measurement_status"] == "measured"
    assert detail["molecule"]["activity_pchembl"] == pytest.approx(7)
    assert detail["molecule"]["n_records"] == 3
    assert len(records) == 3
    assert set(records.evidence_unit) == {"R1", "R2", "R3"}
    assert set(records.raw_activity_unit) == {"nM", "uM", "µM"}
    assert session.new_result.report["excluded_activity_records"] == 1
    filename = "new_molecules_excluded_records.csv"
    assert filename in api(session, "/api/new-molecules")["data"]["files"]
    assert filename in api(session, "/api/session")["data"]["job"]["files"]
    downloaded = api(session, "/api/file?name=" + filename)
    assert downloaded["status"] == 200 and downloaded["filename"] == filename
    assert downloaded["data"] == session.files[filename].read_bytes()
    exported = pd.read_csv(io.BytesIO(downloaded["data"]), dtype=str)
    assert exported.columns.tolist() == ["id", "smiles", "value", "unit", "relation", "batch", "label", "exclusion_reason"]
    assert exported["relation"].tolist() == ["<"]
    assert exported["exclusion_reason"].tolist() == ["non_exact_relation"]
    matched, _ = session.inspect(molecule_id="CTRL")
    assert matched["molecule"]["molecule_id"] == "M1"
    own = matched["new_molecule_measurement"]
    assert own == {"own_median": 6., "reference_median": 5., "median_difference": 1., "own_n_records": 2}
    assert session.job["report"]["query_reference_differences"] == {
        "median": 1., "min": 1., "max": 1., "count": 1}
    for path in ("/api/session", "/api/new-molecules", "/api/points", "/api/inspect?id=N1", "/api/inspect?id=CTRL"):
        assert api(session, path)["status"] == 200


def test_new_id_collision_keeps_new_and_reference_selections_distinct(session):
    insert(session, b"id,smiles\nM1,CCF\nMATCH,OCC\nDUP,CCO\n")
    new, _ = session.inspect(molecule_id="M1")
    assert new["molecule"]["molecule_identity_key"] == "query:M1"
    assert new["molecule"]["canonical_smiles"] == "CCF"
    reference = session.result.data[session.result.data.molecule_id.eq("M1")].iloc[0]
    by_cell, _ = session.inspect(row=reference.grid_row, col=reference.grid_col)
    assert by_cell["molecule"]["canonical_smiles"] == "CCO"
    for identity in ("MATCH", "DUP"):
        matched, _ = session.inspect(molecule_id=identity)
        assert matched["molecule"]["canonical_smiles"] == "CCO"
        assert matched["query_match"]["reference_id"] == "M1"


def test_existing_group_and_reference_audit_exports_survive_insertion(session):
    session.start({**OPTIONS, "grouping": "butina"}, background=False)
    group = session.result.structure_groups[0]["structure_group"]
    comparison = session.compare_group(group)
    inspected, _ = session.inspect(molecule_id="M1")
    paths = {name: session.files[name] for name in (comparison["archive"], inspected["audit_file"])}
    originals = {name: path.read_bytes() for name, path in paths.items()}
    insert(session)
    assert session.compare_group(group) == comparison
    assert all(session.files[name] == path for name, path in paths.items())
    assert all(paths[name].read_bytes() == content for name, content in originals.items())
    session.clear_query()
    assert session.compare_group(group) == comparison
    assert all(session.files[name] == path for name, path in paths.items())
    assert all(paths[name].read_bytes() == content for name, content in originals.items())


def test_generic_reference_span_threshold_is_reused_by_new_measurements(session):
    session.start({**OPTIONS, "compound_column": "molecule_id", "conflict_range_threshold": 2},
                  background=False)
    assert session.job["state"] == "done", session.job
    assert session.result.insertion_parameters["conflict_range_threshold"] == 2
    insert(session, b"id,smiles,value\nN1,CCF,5\nN1,CCF,6.5\nN2,CCCl,5\nN2,CCCl,7.1\n",
           activity_col="value")
    own = session.new_result.data.set_index("query_id")
    assert own.loc["N1", "pchembl_range"] == pytest.approx(1.5)
    assert own.loc["N1", "wide_record_range"] == 0
    assert own.loc["N2", "pchembl_range"] == pytest.approx(2.1)
    assert own.loc["N2", "wide_record_range"] == 1
