"""Browser inputs use the existing units, evidence and unmeasured-molecule API."""

from email.message import Message
import io
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from chemgridmap.app import AppSession, DEMO, LocalHandler, MAX_UPLOAD, json_rows, validate_options
from chemgridmap.chembl import curate_chembl_activity_data


QUERY_DEMO = Path(__file__).with_name("data") / "caii_queries.csv"
SMALL = b'molecule_id,smiles,value\nM1,CCO,5\nM2,CCN,7\nM3,CCC,8\n'


def api(session, path, content=b"", *, method="POST", headers=None, token=True):
    """Exercise the real handler without opening a socket in restricted CI."""
    handler = object.__new__(LocalHandler)
    handler.server = SimpleNamespace(session=session, origin="http://127.0.0.1:8765")
    handler.path = path
    merged = {"Host": "127.0.0.1:8765", "Content-Length": str(len(content))}
    if token:
        merged["X-ChemGridMap-Token"] = session.token
    merged.update(headers or {})
    handler.headers = Message()
    for key, value in merged.items():
        handler.headers[key] = value
    handler.rfile = io.BytesIO(content)
    response = {}

    def send(data, content_type="application/json", status=200, filename=None):
        # Match production serialization, including its ban on non-finite numbers.
        response.update(status=status, payload=json.loads(json.dumps(data, allow_nan=False)))

    handler.send = send
    getattr(handler, f"do_{method}")()
    return response["status"], response["payload"]


@pytest.mark.parametrize("path", ["/api/query-input", "/api/query-input/clear"])
@pytest.mark.parametrize("headers,token", [({}, False), ({"Host": "example.com"}, True),
                                         ({"Origin": "https://example.com"}, True)])
def test_query_input_requires_session_host_and_origin(tmp_path, path, headers, token):
    session = AppSession(tmp_path)
    session.load(SMALL, "measured.csv")
    assert api(session, path, b"smiles\nCCF\n", headers=headers, token=token)[0] == 403
    assert session.query_frame is None


@pytest.mark.parametrize("delimiter", [",", ";", "\t", "|"])
def test_query_upload_delimiters_clear_and_main_input_replacement(tmp_path, delimiter):
    session = AppSession(tmp_path)
    session.load(SMALL, "measured.csv")
    content = (f"query_id{delimiter}SMILES\n0001{delimiter}CCF\n0002{delimiter}CCCl\n").encode("utf-8-sig")
    status, info = api(session, "/api/query-input", content, headers={"X-Filename": "queries.tsv"})
    assert status == 200
    assert info["filename"] == "queries.tsv" and info["rows"] == 2
    assert info["columns"] == ["query_id", "SMILES"] and info["detected"]["smiles"] == "SMILES"
    assert session.query_frame.query_id.tolist() == ["0001", "0002"]
    assert api(session, "/api/session", method="GET")[1]["query_input"] == info
    assert api(session, "/api/query-input/clear")[1] == {"query_input": None}
    assert session.query_frame is None and not session.query_raw and not session.query_filename
    assert "query_input" not in api(session, "/api/session", method="GET")[1]
    api(session, "/api/query-input", content)
    assert api(session, "/api/input", SMALL)[0] == 200
    assert session.query_frame is None and not session.query_raw and not session.query_filename


@pytest.mark.parametrize("length", [-1, MAX_UPLOAD + 1])
def test_query_upload_size_limit_before_reading(tmp_path, length):
    session = AppSession(tmp_path)
    session.load(SMALL, "measured.csv")
    status, detail = api(session, "/api/query-input", headers={"Content-Length": str(length)})
    assert status == 400 and "50 MB" in detail["error"]
    assert session.query_frame is None


@pytest.mark.parametrize("content", [b"", b"smiles\n", b"smiles\n\xff\n"])
def test_query_upload_rejects_empty_or_non_utf8_data(tmp_path, content):
    session = AppSession(tmp_path)
    session.load(SMALL, "measured.csv")
    assert api(session, "/api/query-input", content)[0] == 400
    assert session.query_frame is None


def test_running_map_cannot_replace_or_clear_queries(tmp_path):
    session = AppSession(tmp_path)
    session.load(SMALL, "measured.csv")
    session.load_query(b"smiles\nCCF\n", "queries.csv")
    session.update(state="running")
    for path in ("/api/query-input", "/api/query-input/clear"):
        assert api(session, path, b"smiles\nCCCl\n")[0] == 400
    assert session.query_frame.smiles.tolist() == ["CCF"]


@pytest.mark.parametrize("options", [
    {"value_unit": "NM"}, {"value_unit": "nM", "unit_column": "units"},
    {"value_unit": "p", "unit_column": "units"}, {"unit_column": ""},
    {"assumed_sigma": 0}, {"assumed_sigma": -1}, {"assumed_sigma": float("nan")},
    {"assumed_sigma": float("inf")}, {"assumed_sigma": .64, "activity_confidence": False},
    {"input_format": "molecule", "assumed_sigma": .64},
    {"query_smiles_col": []}, {"query_id_col": ""},
])
def test_new_options_reject_invalid_combinations(options):
    with pytest.raises(ValueError):
        validate_options(options)


def test_new_options_are_conditional_and_columns_are_checked(tmp_path):
    new = {"value_unit", "unit_column", "assumed_sigma", "query_smiles_col", "query_id_col"}
    assert not (new & validate_options({}).keys())
    assert validate_options({"input_format": "molecule", "unit_column": "units"})["unit_column"] == "units"
    session = AppSession(tmp_path)
    session.load(SMALL, "measured.csv")
    with pytest.raises(ValueError, match="query_id_col"):
        session.start({"query_id_col": "query_id"}, background=False)
    session.load_query(b"smiles\nCCF\n", "queries.csv")
    with pytest.raises(ValueError, match="query_smiles_col"):
        session.start({"query_smiles_col": "missing"}, background=False)
    with pytest.raises(ValueError, match="unit column"):
        session.start({"input_format": "molecule", "unit_column": "missing"}, background=False)


@pytest.fixture(scope="module")
def caii_sessions(tmp_path_factory):
    root = tmp_path_factory.mktemp("app_caii_queries")
    sessions = []
    for name, query in (("reference", False), ("with_queries", True)):
        session = AppSession(root / name)
        session.load(DEMO.read_bytes(), DEMO.name)
        options = {"target_id": "CHEMBL205", "projection": "pca", "render_detail": "overview",
                   "output_formats": ["svg"]}
        if query:
            session.load_query(QUERY_DEMO.read_bytes(), QUERY_DEMO.name)
            options.update(query_smiles_col="smiles", query_id_col="query_id")
        session.start(options, background=False)
        assert session.job["state"] == "done", session.job
        sessions.append(session)
    return sessions


def test_caii_query_points_inspector_and_export_agree(caii_sessions):
    baseline, session = caii_sessions
    status, points = api(session, "/api/points", method="GET")
    assert status == 200 and len(points) == 934
    assert sum(point["measurement_status"] == "not_measured" for point in points) == 45
    assert all("measurement_status" not in point for point in api(baseline, "/api/points", method="GET")[1])
    report = session.job["report"]["query_molecules"]
    assert {key: report[key] for key in ("mapped", "invalid_smiles", "duplicate_query", "in_reference")} == {
        "mapped": 45, "invalid_smiles": 3, "duplicate_query": 2, "in_reference": 5}
    assert report["structure_policy"] == session.curation.report["parameters"]["structure_policy"]
    assert {"new_molecules_check.csv", "new_molecules_neighbors.csv"} <= session.files.keys()
    before = set(Path(session.job["output_dir"]).iterdir())
    status, detail = api(session, "/api/inspect?id=ANALOG_001", method="GET")
    assert status == 200 and detail["molecule"]["measurement_status"] == "not_measured"
    assert not ({"audit_file", "neighbors"} & detail.keys())
    assert detail["records"] == []
    assert detail["summary"]["reference_support_assessed"]
    assert set(Path(session.job["output_dir"]).iterdir()) == before
    assert detail["molecule"].get("activity_pchembl") is None
    assert detail["molecule"].get("class_confidence") is None
    saved = pd.read_csv(session.files["new_molecules_neighbors.csv"])
    expected = saved[saved.query_id.eq("ANALOG_001")]
    assert detail["query_neighbors"] == json_rows(expected)
    row = detail["molecule"]
    by_cell = api(session, f'/api/inspect?row={row["grid_row"]}&col={row["grid_col"]}', method="GET")[1]
    assert by_cell == detail


def test_query_addition_preserves_measured_audit_values(caii_sessions):
    baseline, session = caii_sessions
    before, before_records = baseline.inspect(molecule_id="CHEMBL20")
    after, after_records = session.inspect(molecule_id="CHEMBL20")
    pd.testing.assert_frame_equal(before_records, after_records)
    assert session.activity_confidence == baseline.activity_confidence
    # Insertion leaves the complete reference geometry and audit digests unchanged.
    assert before["summary"] == after["summary"]
    for column, value in before["molecule"].items():
        assert after["molecule"][column] == value, column


def test_own_nm_records_with_assumed_sigma_and_generic_id_lookup(tmp_path):
    curated = curate_chembl_activity_data(pd.read_csv(DEMO), target_id="CHEMBL205")
    records = curated.retained_records
    own = pd.DataFrame({"molecule_identity_key": records.molecule_identity_key,
                        "smiles": records.representation_smiles,
                        "IC50_nM": np.power(10., 9. - records.activity_pchembl_raw.to_numpy(float)),
                        "relation": "="})
    session = AppSession(tmp_path)
    session.load(own.to_csv(index=False).encode(), "caii_own_data_nm.csv")
    session.start({"input_format": "molecule", "compound_column": "molecule_identity_key",
                   "smiles_col": "smiles", "value_col": "IC50_nM", "value_unit": "nM",
                   "relation_column": "relation", "assumed_sigma": .64, "projection": "pca",
                   "render_detail": "overview", "output_formats": ["svg"]}, background=False)
    assert session.job["state"] == "done", session.job
    assert len(session.result.data) == 889
    assert {"supported", "uncertain"} <= set(session.result.data.class_status)
    status, detail = api(session, "/api/inspect?id=CHEMBL20", method="GET")
    assert status == 200
    summary = detail["summary"]
    assert summary["method"] == "assumed_sigma" and summary["sigma_doc"] == summary["assumed_sigma"] == .64
    assert summary["evidence_unit_source"] == "one unit per compound"
    assert "sigma_doc_degrees_of_freedom" not in summary
    assert session.job["report"]["activity_confidence"]["method"] == "assumed_sigma"
    expected = curated.molecules.set_index("molecule_identity_key").activity_pchembl
    actual = session.result.data.set_index("molecule_identity_key").activity_pchembl
    np.testing.assert_allclose(actual.loc[expected.index], expected, rtol=0, atol=1e-12)


@pytest.mark.parametrize("use_column", [False, True])
def test_molecule_concentrations_and_generic_ids(tmp_path, use_column):
    session = AppSession(tmp_path)
    session.load(b'molecule_id,smiles,value,units\nA1,CCO,1000,nM\nA2,CCN,10,nM\nA3,CCC,1,nM\n', "molecules.csv")
    units = {"unit_column": "units"} if use_column else {"value_unit": "nM"}
    session.start({"input_format": "molecule", "value_col": "value", "projection": "pca",
                   "output_formats": ["svg"], **units}, background=False)
    assert session.job["state"] == "done", session.job
    status, detail = api(session, "/api/inspect?id=A1", method="GET")
    assert status == 200 and detail["molecule"]["value"] == 6


@pytest.mark.parametrize("document", [None, "batch"])
def test_chembl_assumed_sigma_uses_actual_evidence_column(tmp_path, document):
    session = AppSession(tmp_path)
    session.load(b'Smiles,pChEMBL Value,Molecule ChEMBL ID,batch\nCCO,5,CHEMBL1,A\nCCN,7,CHEMBL2,B\nCCC,8,CHEMBL3,C\n', "records.csv")
    session.start({"projection": "pca", "document_column": document, "assumed_sigma": .64,
                   "output_formats": ["svg"]}, background=False)
    assert session.job["state"] == "done", session.job
    detail, _ = session.inspect(molecule_id="CHEMBL1")
    assert detail["summary"]["evidence_unit_source"] == (document or "one unit per compound")
    assert detail["summary"]["sigma_doc"] == .64 and detail["summary"]["method"] == "assumed_sigma"
    assert session.job["report"]["activity_confidence"]["evidence_unit_source"] == (document or "one unit per compound")


def test_generic_data_audit_keeps_missing_evidence_warning(tmp_path):
    session = AppSession(tmp_path)
    session.load(SMALL, "generic.csv")
    session.start({"input_format": "molecule", "compound_column": "molecule_id",
                   "value_col": "value", "projection": "pca", "output_formats": ["svg"]},
                  background=False)
    assert session.job["state"] == "done", session.job
    assert any("no document identifiers" in warning and "assumed-sigma" in warning
               for warning in session.job["report"]["warnings"])
