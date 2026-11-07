"""未测分子按已有结构检查结果跳转，并复用已测分子的来源文献。"""

from email.message import Message
import json
import os
from pathlib import Path
from types import SimpleNamespace
from urllib.parse import quote

import pandas as pd
import pytest

from chemgridmap.app import AppSession, DEMO
from chemgridmap.app import LocalHandler
from chemgridmap.queries import QUERY_CHECK_COLUMNS


OPTIONS = {"projection": "pca", "grouping": "none", "grid_refinement": "none",
           "render_detail": "overview", "output_formats": ["svg"]}
INPUTS = Path(os.environ.get("CHEMGRIDMAP_U6_INPUT_DIR", "/private/tmp/chemgridmap_u6/inputs"))


def inspect_api(session, query_id):
    """通过真实 GET 处理器检查响应，不在沙箱中打开监听端口。"""
    handler = object.__new__(LocalHandler)
    handler.server = SimpleNamespace(session=session, origin="http://127.0.0.1:8765")
    handler.path = "/api/inspect?id=" + quote(query_id)
    handler.headers = Message()
    handler.headers["Host"] = "127.0.0.1:8765"
    handler.headers["X-ChemGridMap-Token"] = session.token
    response = {}
    handler.send = lambda data, content_type="application/json", **kwargs: response.update(
        data=json.loads(json.dumps(data, allow_nan=False)), **kwargs)
    handler.do_GET()
    return response


def run_query_session(root, measured, query, *, query_id_col="query_id", **options):
    session = AppSession(root)
    session.load(measured, "measured.csv")
    session.load_query(query, "queries.csv")
    session.start({**OPTIONS, "query_id_col": query_id_col, **options}, background=False)
    assert session.job["state"] == "done", session.job
    return session


@pytest.fixture(scope="module")
def demo_query_session(tmp_path_factory):
    # 仓库中的既有小样例与任务指定的 55 行演示文件逐字节相同。
    source = INPUTS / "demo_unmeasured_molecules_caii.csv"
    if not source.is_file():
        source = Path(__file__).with_name("data") / "caii_queries.csv"
    return run_query_session(tmp_path_factory.mktemp("u6_demo"), DEMO.read_bytes(),
                             source.read_bytes())


@pytest.fixture(scope="module")
def simulated_query_session(tmp_path_factory):
    source = INPUTS / "sim_generated_molecules.csv"
    if not source.is_file():
        pytest.skip("Set CHEMGRIDMAP_U6_INPUT_DIR to copied optional simulated candidates.")
    return run_query_session(tmp_path_factory.mktemp("u6_simulated"), DEMO.read_bytes(),
                             source.read_bytes(), query_id_col="gen_id")


def assert_reference_matches(session, expected_count):
    checks = session.result.query_check
    matched = checks[checks.status.eq("in_reference")]
    assert len(matched) == expected_count
    rows = session.job["report"]["query_reference_matches"]
    assert [row["query_id"] for row in rows] == matched.query_id.tolist()
    by_id = {row["query_id"]: row for row in rows}
    for check in matched.itertuples():
        response = inspect_api(session, check.query_id)
        assert response.get("status", 200) == 200
        detail = response["data"]
        assert detail["query_match"] == {
            "query_id": check.query_id, "status": "in_reference",
            "reference_id": check.reference_id, "first_query_id": "",
        }
        molecule = detail["molecule"]
        assert molecule["molecule_identity_key"] == check.reference_id
        assert molecule["molecule_id"].startswith("CHEMBL")
        assert molecule["measurement_status"] == "measured"
        assert detail["sources"]
        direct, records = session.inspect(molecule_id=check.reference_id)
        assert "query_match" not in direct
        assert detail["sources"] == direct["sources"]
        assert detail["records"] == direct["records"]
        assert len(records) == detail["total_records"]
        assert by_id[check.query_id] == {
            "query_id": check.query_id, "reference_id": check.reference_id,
            "molecule_id": molecule["molecule_id"],
            "parent_molecule_id": molecule.get("parent_molecule_id", ""),
            "source_count": len(direct["sources"]),
            "linked_source_count": sum(bool(source["links"]) for source in direct["sources"]),
        }
    # 浏览器扩充字段不能流入既有报告和 query_check 导出列。
    saved = json.loads(session.files["map_curation_report.json"].read_text(encoding="utf-8"))
    assert "query_reference_matches" not in saved
    assert checks.columns.tolist() == QUERY_CHECK_COLUMNS
    exported = pd.read_csv(session.files["map_query_check.csv"], keep_default_na=False)
    assert exported.columns.tolist() == QUERY_CHECK_COLUMNS
    pd.testing.assert_frame_equal(checks.fillna(""), exported, check_dtype=False)


def test_demo_all_five_matches_open_reference_evidence_and_audit_counts(demo_query_session):
    assert len(demo_query_session.query_frame) == 55
    assert_reference_matches(demo_query_session, 5)


def test_simulated_all_eight_matches_open_reference_evidence_and_audit_counts(simulated_query_session):
    assert len(simulated_query_session.query_frame) == 313
    assert_reference_matches(simulated_query_session, 8)


@pytest.mark.parametrize("suffix", ["001", "002"])
def test_demo_duplicate_opens_first_mapped_query(demo_query_session, suffix):
    session = demo_query_session
    duplicate = "DUPLICATE_" + suffix
    first_id = "ANALOG_" + suffix
    detail, records = session.inspect(molecule_id=duplicate)
    first, _ = session.inspect(molecule_id=first_id)
    assert detail["query_match"] == {
        "query_id": duplicate, "status": "duplicate_query", "reference_id": "",
        "first_query_id": first_id,
    }
    assert detail["molecule"] == first["molecule"]
    assert detail["query_neighbors"] == first["query_neighbors"]
    assert records.empty


@pytest.mark.parametrize("suffix", ["001", "002", "003"])
def test_demo_invalid_smiles_has_specific_api_error(demo_query_session, suffix):
    query_id = "INVALID_" + suffix
    response = inspect_api(demo_query_session, query_id)
    assert response["status"] == 400
    assert response["data"]["error"] == f"Unmeasured molecule {query_id} has invalid SMILES."


def test_measured_ids_and_grid_cells_do_not_include_query_match(demo_query_session):
    session = demo_query_session
    for molecule_id in ("CHEMBL20", "ANALOG_001"):
        detail, _ = session.inspect(molecule_id=molecule_id)
        assert "query_match" not in detail
        molecule = detail["molecule"]
        by_cell, _ = session.inspect(row=molecule["grid_row"], col=molecule["grid_col"])
        assert "query_match" not in by_cell


def test_duplicate_of_reference_match_preserves_both_query_checks(tmp_path):
    measured = b"molecule_id,smiles,value\nM1,CCO,5\nM2,CCN,7\nM3,CCC,8\n"
    session = run_query_session(
        tmp_path, measured, b"query_id,smiles\nFIRST,OCC\nREPEAT,CCO\nNEW,CCF\n",
        input_format="molecule", value_col="value")
    original_check = session.result.query_check.copy(deep=True)
    first, _ = session.inspect(molecule_id="FIRST")
    duplicate, _ = session.inspect(molecule_id="REPEAT")
    assert first["molecule"]["molecule_id"] == duplicate["molecule"]["molecule_id"] == "M1"
    assert first["query_match"] == {
        "query_id": "FIRST", "status": "in_reference", "reference_id": "M1", "first_query_id": "",
    }
    assert duplicate["query_match"] == {
        "query_id": "REPEAT", "status": "duplicate_query", "reference_id": "M1",
        "first_query_id": "FIRST",
    }
    pd.testing.assert_frame_equal(original_check, session.result.query_check)
    assert session.job["report"]["query_reference_matches"] == [{
        "query_id": "FIRST", "reference_id": "M1", "molecule_id": "M1",
        "parent_molecule_id": "", "source_count": 0, "linked_source_count": 0,
    }]


@pytest.mark.parametrize("document_column,with_references", [
    (None, False), (None, True), ("run_id", True),
])
def test_generic_sources_count_real_links_and_support_missing_documents(
        tmp_path, document_column, with_references):
    frame = pd.DataFrame([
        {"compound": compound, "smiles": smiles, "activity": 5 + index + .1 * repeat,
         "run_id": document, "reference": reference}
        for index, (compound, smiles) in enumerate([("C1", "CCO"), ("C2", "CCN"), ("C3", "CCC")])
        for repeat, (document, reference) in enumerate([
            ("D1", "doi:10.1021/jm00131a011"), ("D1", "doi:10.1021/jm00131a011"),
            ("D2", "local notebook"), ("D3", "https://example.org/study"),
        ])
    ])
    options = {"input_format": "molecule", "compound_column": "compound", "value_col": "activity",
               "document_column": document_column}
    if with_references:
        options["reference_column"] = "reference"
    session = run_query_session(tmp_path, frame.to_csv(index=False).encode(),
                                b"query_id,smiles\nQ1,OCC\n", **options)
    detail, _ = session.inspect(molecule_id="Q1")
    direct, _ = session.inspect(molecule_id="C1")
    row, = session.job["report"]["query_reference_matches"]
    assert row["source_count"] == (3 if document_column else 0)
    assert row["linked_source_count"] == (2 if document_column else 0)
    if document_column:
        assert row["source_count"] == len(detail["sources"])
        assert row["linked_source_count"] == sum(bool(source["links"]) for source in detail["sources"])
    assert detail["sources"] == direct["sources"]


def test_stereoisomer_is_not_resolved_to_measured_structure(tmp_path):
    measured = b"molecule_id,smiles,value\nM1,C[C@H](O)F,5\nM2,CCN,7\nM3,CCC,8\n"
    session = run_query_session(
        tmp_path, measured, b"query_id,smiles\nOPPOSITE,C[C@@H](O)F\n",
        input_format="molecule", value_col="value")
    detail, _ = session.inspect(molecule_id="OPPOSITE")
    assert detail["molecule"]["measurement_status"] == "not_measured"
    assert detail["molecule"]["molecule_id"] == "OPPOSITE"
    assert "query_match" not in detail
    assert session.job["report"]["query_reference_matches"] == []


def test_reference_match_without_measured_identifier_uses_canonical_structure(tmp_path):
    session = run_query_session(
        tmp_path, b"smiles,value\nCCO,5\nCCN,7\nCCC,8\n",
        b"query_id,smiles\nMATCH,OCC\n", input_format="molecule", value_col="value")
    detail, _ = session.inspect(molecule_id="MATCH")
    assert detail["molecule"]["canonical_smiles"] == "CCO"
    assert detail["query_match"] == {
        "query_id": "MATCH", "status": "in_reference", "reference_id": "CCO", "first_query_id": "",
    }
    row, = session.job["report"]["query_reference_matches"]
    assert row["reference_id"] == row["molecule_id"] == "CCO"
    assert row["source_count"] == row["linked_source_count"] == 0
