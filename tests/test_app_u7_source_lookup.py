"""来源查找只进入本机会话与检查器，不改变计算或导出。"""

from email.message import Message
import io
from types import SimpleNamespace
import zipfile

import pandas as pd
import pytest

from chemgridmap import app as app_module
from chemgridmap.app import AppSession, LocalHandler, MAX_UPLOAD


MEASURED = (b"compound,smiles,activity,run_id,reference\n"
            b"C1,CCO,6,CHEMBL1124628,doi:10.1021/jm00131a011\n"
            b"C1,OCC,7,CHEMBL1124628,doi:10.1021/jm00131a011\n"
            b"C2,CCN,5,R2,notebook\nC3,CCC,8,R3,notebook\n")
LOOKUP = (b"canonical_smiles,document_chembl_id,pchembl_value,target_chembl_id\n"
          b"OCC,CHEMBL1124628,6,CHEMBL205\n"
          b"CCN,CHEMBL1909046,5,CHEMBL206\n")
OPTIONS = {"input_format": "molecule", "compound_column": "compound",
           "smiles_col": "smiles", "value_col": "activity", "assumed_sigma": .64,
           "projection": "pca", "grouping": "none", "grid_refinement": "none",
           "render_detail": "overview", "output_formats": ["svg"]}


def request(session, path, body=b"", *, method="POST", token=True, length=None):
    """调用真实处理器，不开监听端口。"""
    handler = object.__new__(LocalHandler)
    handler.server = SimpleNamespace(session=session, origin="http://127.0.0.1:8765")
    handler.path, handler.rfile = path, io.BytesIO(body)
    handler.headers = Message()
    handler.headers["Host"] = "127.0.0.1:8765"
    handler.headers["Content-Length"] = str(len(body) if length is None else length)
    handler.headers["X-Filename"] = "sources.csv"
    if token:
        handler.headers["X-ChemGridMap-Token"] = session.token
    response = {}
    handler.send = lambda data, content_type="application/json", **kwargs: response.update(
        data=data, **kwargs)
    getattr(handler, "do_" + method)()
    return response


def run_session(root, *, lookup=True, **options):
    session = AppSession(root)
    session.load(MEASURED, "measured.csv")
    if lookup:
        session.load_source(LOOKUP, "sources.csv")
    session.start({**OPTIONS, **options}, background=False)
    assert session.job["state"] == "done", session.job
    return session


@pytest.fixture(scope="module")
def source_session(tmp_path_factory):
    return run_session(tmp_path_factory.mktemp("source_lookup"))


def test_source_upload_clear_session_and_replace_measured_input(tmp_path):
    session = AppSession(tmp_path)
    with pytest.raises(ValueError, match="measured data first"):
        session.load_source(LOOKUP, "sources.csv")
    session.load(MEASURED, "measured.csv")
    response = request(session, "/api/source-input", LOOKUP)
    assert response["data"] == {"filename": "sources.csv", "rows": 2,
                                "columns": pd.read_csv(io.BytesIO(LOOKUP)).columns.tolist()}
    payload = request(session, "/api/session", method="GET")["data"]
    assert payload["source_input"] == response["data"]
    session.load(MEASURED, "replacement.csv")
    assert session.source_input_info() is None and session.possible_sources is None
    session.load_source(LOOKUP, "folder\\sources.csv")
    assert session.source_filename == "sources.csv"
    assert request(session, "/api/source-input/clear")["data"] == {"source_input": None}
    assert "source_input" not in request(session, "/api/session", method="GET")["data"]


@pytest.mark.parametrize("action", ["load", "clear"])
def test_source_changes_discard_map_and_keep_saved_files(tmp_path, action):
    session = run_session(tmp_path)
    saved = session.files["map.svg"]
    if action == "load":
        session.load_source(LOOKUP, "replacement.csv")
    else:
        session.clear_source()
    assert session.job["state"] == "ready"
    assert session.result is session.generic_records is session.possible_sources is None
    assert not session.files and saved.exists()


@pytest.mark.parametrize("path", ["/api/source-input", "/api/source-input/clear"])
def test_source_endpoints_require_token_and_enforce_upload_limit(tmp_path, path):
    session = AppSession(tmp_path)
    session.load(MEASURED, "measured.csv")
    assert request(session, path, LOOKUP, token=False)["status"] == 403
    rejected = request(session, path, length=MAX_UPLOAD + 1)
    assert rejected == {"status": 400, "data": {"error": "Maximum upload size is 50 MB."}}
    assert session.source_frame is None


def test_invalid_upload_and_running_job_preserve_source_input(tmp_path):
    session = AppSession(tmp_path)
    session.load(MEASURED, "measured.csv")
    session.load_source(LOOKUP, "sources.csv")
    bad = request(session, "/api/source-input", b"smiles\nCCO\n")
    assert bad == {"status": 400, "data": {"error":
        "The source lookup file needs SMILES and ChEMBL document ID columns, as in a ChEMBL activity export."}}
    assert session.source_filename == "sources.csv"
    session.update(state="running")
    with pytest.raises(ValueError, match="running"):
        session.load_source(LOOKUP, "replacement.csv")
    with pytest.raises(ValueError, match="running"):
        session.clear_source()
    assert session.source_filename == "sources.csv"


def test_source_audit_and_inspector_include_possible_and_missing_sources(source_session):
    session = source_session
    assert session.job["report"]["source_lookup"] == {
        "file": "sources.csv", "total": 3, "found": 2, "same": 2, "target_count": 2}
    detail, _ = session.inspect(molecule_id="C1")
    assert detail["sources"] == []
    assert detail["source_lookup_file"] == "sources.csv"
    source, = detail["possible_sources"]
    assert source["id"] == "CHEMBL1124628" and source["same_values"] == ["6.00"]
    empty, _ = session.inspect(molecule_id="C3")
    assert empty["possible_sources"] == [] and empty["source_lookup_file"] == "sources.csv"


@pytest.mark.parametrize("selected", [{"document_column": "run_id"},
                                     {"reference_column": "reference"}])
def test_own_linked_sources_hide_lookup_section(tmp_path, selected):
    session = run_session(tmp_path, **selected)
    linked, _ = session.inspect(molecule_id="C1")
    assert any(source["links"] for source in linked["sources"])
    assert "possible_sources" not in linked and "source_lookup_file" not in linked
    unlinked, _ = session.inspect(molecule_id="C2")
    assert unlinked["possible_sources"]


def test_unmeasured_molecule_has_no_possible_sources(tmp_path):
    session = AppSession(tmp_path)
    session.load(MEASURED, "measured.csv")
    session.load_source(LOOKUP, "sources.csv")
    session.load_query(b"query_id,smiles\nQ1,CCF\n", "queries.csv")
    session.start({**OPTIONS, "query_id_col": "query_id"}, background=False)
    assert session.job["state"] == "done", session.job
    detail, _ = session.inspect(molecule_id="Q1")
    assert "possible_sources" not in detail and "source_lookup_file" not in detail
    assert session.job["report"]["source_lookup"]["total"] == 3


@pytest.mark.parametrize("mode", ["chembl", "molecule_table", "no_file"])
def test_matching_is_not_called_without_generic_source_input(tmp_path, monkeypatch, mode):
    session = AppSession(tmp_path)
    session.load(MEASURED, "measured.csv")
    if mode != "no_file":
        session.load_source(LOOKUP, "sources.csv")
    monkeypatch.setattr(app_module, "match_source_lookup", lambda *args: pytest.fail("unexpected lookup"))
    options = {**OPTIONS}
    if mode == "chembl":
        session.load(b"Smiles,pChEMBL Value,Molecule ChEMBL ID\nCCO,6,CHEMBL20\nCCN,5,CHEMBL200\nCCC,8,CHEMBL300\n", "chembl.csv")
        session.load_source(LOOKUP, "sources.csv")
        options.update(input_format="chembl", compound_column=None)
    elif mode == "molecule_table":
        session.load(b"smiles,activity\nCCO,6\nCCN,5\nCCC,8\n", "molecules.csv")
        session.load_source(LOOKUP, "sources.csv")
        options.update(compound_column=None, assumed_sigma=None)
    session.start(options, background=False)
    assert session.job["state"] == "done", session.job
    assert session.possible_sources is None and "source_lookup" not in session.job["report"]
    detail, _ = session.inspect(row=int(session.result.data.iloc[0].grid_row),
                                col=int(session.result.data.iloc[0].grid_col))
    assert "possible_sources" not in detail and "source_lookup_file" not in detail


def test_lookup_runs_once_during_build_and_never_leaks_into_exports(tmp_path, monkeypatch):
    calls, progress = [], []
    original = app_module.match_source_lookup
    def tracked(*args):
        calls.append(1)
        return original(*args)
    monkeypatch.setattr(app_module, "match_source_lookup", tracked)
    session = AppSession(tmp_path)
    session.load(MEASURED, "measured.csv")
    session.load_source(LOOKUP, "distinct-source-lookup.csv")
    update = session.update
    def track_update(**values):
        progress.append(values.get("message"))
        update(**values)
    monkeypatch.setattr(session, "update", track_update)
    session.start(OPTIONS, background=False)
    assert session.job["state"] == "done", session.job
    session.inspect(molecule_id="C1")
    session.inspect(molecule_id="C2")
    assert calls == [1]
    assert "Matching structures to the source lookup file..." in progress
    for path in session.files.values():
        contents = [path.read_bytes()]
        if path.suffix == ".zip":
            with zipfile.ZipFile(path) as archive:
                contents = [archive.read(name) for name in archive.namelist()]
        for content in contents:
            assert b"distinct-source-lookup.csv" not in content
            assert b"possible_sources" not in content and b'"source_lookup"' not in content


def test_u6_reference_table_builds_indexes_once(tmp_path, monkeypatch):
    from chemgridmap.inspection import build_sources
    frame = pd.DataFrame({"molecule_identity_key": [f"C{i}" for i in range(100)],
                          "molecule_id": [f"C{i}" for i in range(100)],
                          "canonical_smiles": ["CCO"] * 100})
    records = pd.DataFrame({"molecule_identity_key": [f"C{i}" for i in range(100)],
                            "doc": ["CHEMBL1124628"] * 100})
    checks = pd.DataFrame({"status": ["in_reference"] * 100,
                          "query_id": [f"Q{i}" for i in range(100)],
                          "reference_id": [f"C{i}" for i in range(100)]})
    expected_source = build_sources(records.iloc[:1], "doc")[0]
    count = []
    groupby = pd.DataFrame.groupby
    def track_groups(self, by=None, *args, **kwargs):
        if isinstance(by, str) and by == "molecule_identity_key":
            count.append(1)
        return groupby(self, by, *args, **kwargs)
    monkeypatch.setattr(pd.DataFrame, "groupby", track_groups)
    monkeypatch.setattr(app_module, "_select_molecule", lambda *args: pytest.fail("repeated scan"))
    result = SimpleNamespace(data=frame, query_check=checks)
    rows = AppSession(tmp_path)._query_reference_matches(result, None, records, {"document_column": "doc"})
    assert count == [1]
    assert rows == [{"query_id": f"Q{i}", "reference_id": f"C{i}", "molecule_id": f"C{i}",
                     "parent_molecule_id": "", "source_count": 1,
                     "linked_source_count": int(bool(expected_source["links"]))} for i in range(100)]
