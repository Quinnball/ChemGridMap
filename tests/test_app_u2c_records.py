"""Generic-record browser audits preserve the source rows without changing exports."""

from email.message import Message
import io
import json
import math
from types import SimpleNamespace

import pandas as pd
import pytest

from chemgridmap.app import AppSession, LocalHandler


GENERIC_NOTE = (
    "Generic activity records: the listed rows are the retained records for this compound; "
    "values were converted to pActivity."
)
OPTIONS = {
    "input_format": "molecule", "compound_column": "compound_id", "smiles_col": "smiles",
    "value_col": "ic50", "unit_column": "unit", "relation_column": "relation",
    "document_column": "run_id", "projection": "pca", "grouping": "none",
    "output_formats": ["svg"],
}


def assay_frame():
    rows = []
    for index, smiles in enumerate(("CCO", "CCN", "CCC", "CCF", "CCCl")):
        for run, offset in (("R1", 0), ("R2", .03), ("R3", 1.2)):
            unit = "nM" if run != "R2" else "µM"
            scale = 9 if unit == "nM" else 6
            rows.append({"note": "retained", "smiles": smiles, "relation": "=",
                         "ic50": 10. ** (scale - (5.7 + .3 * index + offset)),
                         "unit": unit, "compound_id": f"C{index + 1}", "run_id": run})
    rows.append({**rows[0], "smiles": "CCBr", "compound_id": "C6"})
    rows.extend([{**rows[0], "relation": ">", "note": "censored"},
                 {**rows[0], "ic50": None, "note": "failed"}])
    return pd.DataFrame(rows)


def build_session(directory, frame, options=None):
    session = AppSession(directory)
    session.load(frame.to_csv(index=False).encode(), "own_records.csv")
    session.start(options or OPTIONS, background=False)
    assert session.job["state"] == "done", session.job
    return session


def evidence_csv(session, compound):
    # Invoke the production HTTP handler without binding a socket.
    handler = object.__new__(LocalHandler)
    handler.server = SimpleNamespace(session=session, origin="http://127.0.0.1:8765")
    handler.path = f"/api/evidence.csv?id={compound}"
    handler.headers = Message()
    handler.headers["Host"] = "127.0.0.1:8765"
    handler.headers["X-ChemGridMap-Token"] = session.token
    response = {}
    handler.send = lambda data, content_type="application/json", status=200, filename=None: response.update(
        data=data, content_type=content_type, status=status, filename=filename)
    handler.do_GET()
    assert response["status"] == 200
    assert response["content_type"] == "text/csv"
    assert response["filename"] == "source_records.csv"
    return pd.read_csv(io.BytesIO(response["data"]))


@pytest.fixture(scope="module")
def generic_session(tmp_path_factory):
    return build_session(tmp_path_factory.mktemp("u2c_records"), assay_frame())


def test_generic_records_counts_median_columns_offsets_and_csv(generic_session):
    session = generic_session
    detail, records = session.inspect(molecule_id="C1")
    expected_columns = ["source_row", "compound_id", "run_id", "ic50", "unit", "relation",
                        "activity_pchembl_raw", "note", "smiles"]
    original_columns = expected_columns.copy()
    expected_columns += ["document_doi", "document_pubmed_id", "document_patent_id", "document_url"]
    assert records.columns.tolist() == expected_columns
    assert detail["record_columns"] == original_columns
    assert list(detail["records"][0]) == original_columns
    assert detail["total_records"] == detail["molecule"]["n_records"] == len(records) == 3
    assert records.source_row.tolist() == [0, 1, 2]
    assert set(records.unit) == {"nM", "µM"}
    assert records.note.tolist() == ["retained"] * 3
    assert records.activity_pchembl_raw.median() == pytest.approx(
        detail["molecule"]["activity_pchembl"], rel=0, abs=1e-9)
    summary = detail["summary"]
    assert summary["record_count_verified"] and summary["median_verified"]
    assert summary["input_format"] == "generic_activity"
    assert summary["interpretation"] == GENERIC_NOTE
    assert summary["warnings"] == [
        "Offset document pair: R1 vs R3, median difference 1.2 log units, n shared 5; "
        "the pair touches 6 entries. This entry has records from R1, R3; "
        "the display is unchanged because the offset rule is off.",
        "Offset document pair: R2 vs R3, median difference 1.17 log units, n shared 5; "
        "the pair touches 5 entries. This entry has records from R2, R3; "
        "the display is unchanged because the offset rule is off.",
    ]
    assert not {"n_assays", "n_documents"} & summary.keys()
    downloaded = evidence_csv(session, "C1")
    pd.testing.assert_frame_equal(downloaded[original_columns], records[original_columns].reset_index(drop=True))
    assert downloaded.columns.tolist() == expected_columns
    assert downloaded[expected_columns[-4:]].isna().all().all()
    # A compound represented only in R1 must not inherit the unrelated R2/R3 pair.
    single, _ = session.inspect(molecule_id="C6")
    assert len(single["summary"]["warnings"]) == 1
    assert "R1 vs R3" in single["summary"]["warnings"][0]
    assert "records from R1;" in single["summary"]["warnings"][0]


def test_generic_job_shows_all_evaluated_pairs_without_mutating_fit(generic_session):
    session = generic_session
    confidence = session.job["report"]["activity_confidence"]
    assert [(pair["document_a"], pair["document_b"], pair["systematic_offset"])
            for pair in confidence["document_offsets_top10"]] == [
        ("R1", "R3", True), ("R2", "R3", True), ("R1", "R2", False)]
    assert confidence["document_pair_count"] == 3
    assert confidence["discordant_document_pair_count"] == 2
    assert len(session.activity_confidence["document_offsets_top10"]) == 2
    manifest = json.loads(session.files["map_run_manifest.json"].read_text(encoding="utf-8"))
    assert manifest["activity_confidence"] == session.activity_confidence
    assert manifest["source"]["generic_activity"]["activity_confidence"] == session.activity_confidence


@pytest.mark.parametrize("missing", [False, True])
def test_offset_hints_use_normalized_evidence_ids_but_display_original_values(tmp_path, missing):
    frame = assay_frame()
    frame.loc[frame.run_id.eq("R1"), "run_id"] = None if missing else " R1 "
    session = build_session(tmp_path, frame)
    detail, records = session.inspect(molecule_id="C6")
    assert len(detail["summary"]["warnings"]) == 1
    assert ("unknown_document" if missing else "R1") in detail["summary"]["warnings"][0]
    assert pd.isna(records.run_id.iloc[0]) if missing else records.run_id.iloc[0] == " R1 "


@pytest.mark.parametrize("change,message", [("count", "Record counts"), ("median", "Recalculated activity")])
def test_generic_audit_rejects_mismatched_map(generic_session, change, message):
    session = generic_session
    original = session.result.data.copy()
    index = session.result.data.index[session.result.data.molecule_id.eq("C1")][0]
    column = "n_records" if change == "count" else "activity_pchembl"
    session.result.data.loc[index, column] += 1
    try:
        with pytest.raises(ValueError, match=message):
            session.inspect(molecule_id="C1")
    finally:
        session.result.data = original


def test_generic_records_are_not_truncated_at_100(tmp_path):
    frame = assay_frame()
    frame = pd.concat([frame, pd.concat([frame.iloc[[0]]] * 101)], ignore_index=True)
    session = build_session(tmp_path, frame)
    detail, records = session.inspect(molecule_id="C1")
    assert detail["total_records"] == len(detail["records"]) == len(records) == 104
    assert len(evidence_csv(session, "C1")) == 104


@pytest.mark.parametrize("document", [None, "run_id"])
def test_p_only_generic_job_keeps_reports_and_inspects_legacy_records(tmp_path, document):
    frame = assay_frame()
    frame["activity"] = frame["ic50"].map(lambda value: -math.log10(value)
                                         if pd.notna(value) else value)
    frame["activity"] += frame.unit.map({"nM": 9, "µM": 6})
    frame = frame.drop(columns=["ic50", "unit"])
    options = {key: value for key, value in OPTIONS.items() if key != "unit_column"}
    options.update(value_col="activity", document_column=document)
    session = build_session(tmp_path, frame, options)
    report = session.job["report"]
    assert {"generic_activity", "activity_confidence"} <= report.keys()
    assert (report["activity_confidence"] is None) == (document is None)
    assert "excluded_reasons" not in report["generic_activity"]
    assert "activity_pchembl_raw" not in session.generic_records
    detail, records = session.inspect(molecule_id="C1")
    assert detail["summary"]["median_verified"]
    assert "activity_pchembl_raw" in records
    assert "activity_pchembl_raw" not in session.generic_records


def test_chembl_job_pair_preview_does_not_mutate_curation_or_inspection(tmp_path):
    frame = assay_frame()
    frame = frame[frame.note.eq("retained")].copy()
    frame["pChEMBL Value"] = frame.ic50.map(lambda value: -math.log10(value))
    frame["pChEMBL Value"] += frame.unit.map({"nM": 9, "µM": 6})
    frame = frame.rename(columns={"compound_id": "Molecule ChEMBL ID", "smiles": "Smiles",
                                  "run_id": "Document ChEMBL ID"})
    session = build_session(tmp_path, frame, {"input_format": "chembl", "projection": "pca",
                                            "grouping": "none", "output_formats": ["svg"]})
    assert len(session.job["report"]["activity_confidence"]["document_offsets_top10"]) == 3
    assert len(session.curation.report["activity_confidence"]["document_offsets_top10"]) == 2
    saved = json.loads(session.files["map_curation_report.json"].read_text(encoding="utf-8"))
    assert saved == session.curation.report
    detail, records = session.inspect(molecule_id="C1")
    assert "record_columns" not in detail
    assert "input_format" not in detail["summary"]
    assert len(detail["summary"]["document_offsets"]) == 2
    assert detail["summary"]["warnings"][0].startswith("Offset document pair: R1 vs R3,")
    assert "molecule_identity_key" in records


def test_molecule_table_keeps_empty_record_audit(tmp_path):
    frame = pd.DataFrame({"molecule_id": ["M1", "M2", "M3"], "smiles": ["CCO", "CCN", "CCC"],
                          "value": [5., 7., 8.]})
    session = build_session(tmp_path, frame, {"input_format": "molecule", "value_col": "value",
                                            "projection": "pca", "grouping": "none",
                                            "output_formats": ["svg"]})
    detail, records = session.inspect(molecule_id="M1")
    assert records.empty and detail["records"] == [] and detail["total_records"] == 0
    assert detail["summary"]["interpretation"] == "Molecule-table input: no assay-record audit is available."
    assert not detail["summary"]["median_verified"]
    assert "record_columns" not in detail


@pytest.mark.parametrize("action", ["load", "load_query", "clear_query", "rebuild"])
def test_generic_retained_records_reset_with_map_state(tmp_path, monkeypatch, action):
    session = AppSession(tmp_path)
    session.load(assay_frame().to_csv(index=False).encode(), "own.csv")
    session.generic_records = pd.DataFrame({"source_row": [0]})
    if action == "load":
        session.load(b"smiles\nCCO\n", "replacement.csv")
    elif action == "load_query":
        session.load_query(b"smiles\nCCF\n", "queries.csv")
    elif action == "clear_query":
        session.clear_query()
    else:
        monkeypatch.setattr(session, "run", lambda options: None)
        session.start(OPTIONS, background=False)
    assert session.generic_records is None
