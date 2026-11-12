"""U9b manual entry and measured-result contracts on a fixed reference map."""
import hashlib
import json
from pathlib import Path
import zipfile

import pandas as pd
import pytest

from chemgridmap.app import AppSession, MANUAL_COLUMNS, FEEDBACK_COLUMNS, normalize_measurement
from chemgridmap.queries import insert_new_molecules
from test_app_u8_session import api


@pytest.fixture
def session(tmp_path):
    app = AppSession(tmp_path)
    app.load(b"molecule_id,smiles,value\nR1,CCO,8\nR2,CCN,8\nR3,CCC,8\nR4,c1ccccc1,8\n", "reference.csv")
    app.start({"input_format": "molecule", "value_col": "value", "projection": "pca",
               "grouping": "butina", "group_similarity": .1, "grid_refinement": "none",
               "render_detail": "overview", "output_formats": ["svg"]}, background=False)
    assert app.job["state"] == "done", app.job
    return app


def add(session, smiles="CCF", **fields):
    return session.add_new({"smiles": smiles, **fields})


def upload(session, raw=b"id,smiles\nN1,CCCl\nN2,CCBr\n", **options):
    session.load_query(raw, "original.csv")
    session.insert_new({"id_col": "id", **options})


def test_preview_all_states_and_structure_policy(session):
    assert session.preview_new({"smiles": "bad"})["status"] == "invalid"
    match = session.preview_new({"smiles": "OCC"})
    assert match["status"] == "in_reference" and match["reference_id"] == "R1"
    preview = session.preview_new({"smiles": "FCC"})
    assert preview["status"] == "ok" and preview["canonical_smiles"] == "CCF" and "<svg" in preview["svg"]
    add(session)
    assert session.preview_new({"smiles": "FCC"})["query_id"] == "M001"
    assert add(session, "FCC")["status"] == "duplicate_new"
    assert add(session, "CCO")["status"] == "in_reference"
    assert add(session, "bad")["status"] == "invalid"
    assert len(session.new_result.data) == 1
    session.result.data["structure_policy"] = "parent"
    assert session.preview_new({"smiles": "CCF.[Na+]"})["status"] == "duplicate_new"


def test_append_preserves_reference_existing_rows_and_single_file_evidence(session):
    upload(session)
    reference = session.result.data.copy(deep=True)
    old = session.new_result.data.copy(deep=True)
    raw = session.files["new_molecules_input.csv"].read_bytes()
    reference_files = {key: path.read_bytes() for key, path in session.reference_files.items()}
    response = add(session, label="next", note="note")
    assert response["query_id"] == "M001"
    pd.testing.assert_frame_equal(session.result.data, reference)
    pd.testing.assert_frame_equal(session.new_result.data.iloc[:len(old)].reset_index(drop=True), old)
    assert session.files["new_molecules_input.csv"].read_bytes() == raw
    assert {key: path.read_bytes() for key, path in session.reference_files.items()} == reference_files
    single = insert_new_molecules(session.result, pd.DataFrame({"query_id": ["M001"], "smiles": ["CCF"], "label": ["next"]}),
                                  id_col="query_id", label_column="label")
    actual = session.new_result.data.iloc[[-1]].reset_index(drop=True)
    skip = {"source_row", "grid_row", "grid_col", "grid_x", "grid_y", "placement_distance"}
    columns = [column for column in single.data if column not in skip]
    pd.testing.assert_frame_equal(actual[columns], single.data[columns])
    pd.testing.assert_frame_equal(session.new_result.neighbors[session.new_result.neighbors.query_id.eq("M001")].reset_index(drop=True), single.neighbors)


def test_five_one_by_one_same_evidence_and_groups_as_batch(session):
    smiles = ["CCF", "CCCl", "CCBr", "CCCO", "CCCN"]
    frame = pd.DataFrame({"query_id": [f"M{i:03d}" for i in range(1, 6)], "smiles": smiles})
    batch = insert_new_molecules(session.result, frame, id_col="query_id")
    for value in smiles:
        add(session, value)
    columns = ["query_id", "neighbor_evidence", "neighbor_comparison", "structure_group", "contains_group_core",
               "nearest_reference_id", "tanimoto", "similar_reference_count", "evidence_reference_classes"]
    pd.testing.assert_frame_equal(session.new_result.data[columns], batch.data[columns])


def test_id_collision_auto_sequence_limit_and_http(session):
    assert add(session, query_id="M002")["query_id"] == "M002"
    assert add(session, "CCCl")["query_id"] == "M001"
    with pytest.raises(ValueError, match="unique"):
        add(session, "CCBr", query_id="M002")
    assert add(session, "CCBr")["query_id"] == "M003"
    session.new_result.data = pd.concat([session.new_result.data] * 334, ignore_index=True).iloc[:1000]
    response = api(session, "/api/new-molecules/add", json.dumps({"smiles": "CCCO"}).encode(), method="POST")
    assert response["status"] == 400 and "1,000" in response["data"]["error"]


@pytest.mark.parametrize("value,unit,expected,category", [(100, "nM", 7, "active"), (1000, "nM", 6, "inactive"),
    (.2, "µM", 6.698970004336019, "medium"), (1e-7, "M", 7, "active"), (.0001, "mM", 7, "active"),
    (100000, "pM", 7, "active"), (6.5, "pActivity", 6.5, "medium")])
def test_shared_numeric_conversion_and_thresholds(value, unit, expected, category):
    measured = normalize_measurement({"mode": "value", "value": value, "unit": unit}, lower=6, upper=7)
    assert measured["activity_pchembl"] == pytest.approx(expected)
    assert measured["measured_class"] == category


@pytest.mark.parametrize("value,unit", [("bad", "nM"), (">10", "µM"), (0, "nM"), (-1, "pActivity"),
                                         ("nan", "nM"), ("inf", "nM"), (5, "bad")])
def test_invalid_exact_measurements(value, unit):
    with pytest.raises(ValueError):
        normalize_measurement({"mode": "value", "value": value, "unit": unit}, lower=6, upper=7)


def test_result_class_edit_remove_sources_counts_and_reference_immutable(session):
    upload(session)
    original = session.new_result.data.copy(deep=True)
    neighbors = session.new_result.neighbors.copy(deep=True)
    reference = session.result.data.copy(deep=True)
    result = session.record_result({"query_id": "N1", "measurement": {"mode": "class", "measured_class": "inactive", "batch": "B1", "note": "retest"}})
    assert result["row"]["neighbor_comparison"] == "disagrees"
    assert result["row"]["class_source"] == "user_class" and result["row"]["activity_pchembl"] is None
    assert session.new_molecules_info()["report"]["feedback_counts"] == {"recorded": 1, "agrees": 0, "disagrees": 1, "no_comparison": 0}
    detail, _ = session.inspect(molecule_id="N1")
    assert detail["measured_result"]["editable"] and detail["measured_result"]["note"] == "retest"
    assert session.source_members(["B1"], "new")["map_molecule_count"] == 1
    result = session.record_result({"query_id": "N1", "measurement": {"value": 100, "unit": "nM"}})
    assert result["row"]["neighbor_comparison"] == "agrees"
    session.clear_result({"query_id": "N1"})
    pd.testing.assert_frame_equal(session.new_result.data[original.columns], original)
    assert not session._source_index("new")
    pd.testing.assert_frame_equal(session.new_result.neighbors, neighbors)
    pd.testing.assert_frame_equal(session.result.data, reference)


@pytest.mark.parametrize("evidence", ["none", "mixed", "insufficient"])
def test_noncomparable_evidence_is_not_invented(session, evidence):
    add(session)
    for table in (session.new_result.data, session.new_result.all_molecules):
        table.loc[:, "neighbor_evidence"] = evidence
    result = session.record_result({"query_id": "M001", "measurement": {"mode": "class", "measured_class": "active"}})
    assert result["row"]["neighbor_comparison"] == "no_comparison"
    assert result["row"]["neighbor_evidence"] == evidence


def test_uploaded_measured_readonly_feedback_includes_upload_and_raw_preserved(session):
    raw = b"id,smiles,value,unit,batch\r\nUP,CCCl,100,nM,Lab A\r\n"
    upload(session, raw, value_col="value", unit_column="unit", document_column="batch")
    assert "neighbor_evidence_feedback.csv" not in session.files
    for operation in (session.record_result, session.clear_result):
        with pytest.raises(ValueError, match="read-only"):
            operation({"query_id": "UP", "measurement": {"measured_class": "inactive"}})
    assert not session.inspect(molecule_id="UP")[0]["measured_result"]["editable"]
    add(session, measurement={"mode": "value", "value": 200, "unit": "nM", "batch": "B2", "note": "manual"})
    assert session.files["new_molecules_input.csv"].read_bytes() == raw
    feedback = pd.read_csv(session.files["neighbor_evidence_feedback.csv"], keep_default_na=False)
    assert list(feedback) == FEEDBACK_COLUMNS and len(feedback) == 2
    assert feedback.result_source.tolist() == ["upload", "manual"]
    assert feedback.measured_unit.tolist() == ["nM", "nM"]
    manual = pd.read_csv(session.files["manual_entries.csv"], keep_default_na=False)
    assert list(manual) == MANUAL_COLUMNS and manual.query_id.tolist() == ["M001"]
    assert "+00:00" in manual.entered_at.iloc[0]
    manifest = json.loads(session.files["new_molecules_manifest.json"].read_text(encoding="utf-8"))
    assert manifest["manual_entries"]["rows"] == 1
    assert manifest["manual_entries"]["sha256"] == hashlib.sha256(session.files["manual_entries.csv"].read_bytes()).hexdigest()
    with zipfile.ZipFile(session.files["new_molecules_results.zip"]) as archive:
        for name in ("manual_entries.csv", "neighbor_evidence_feedback.csv", "new_molecules_input.csv"):
            assert archive.read(name) == session.files[name].read_bytes()
    session.clear_result({"query_id": "M001"})
    assert len(pd.read_csv(session.files["neighbor_evidence_feedback.csv"])) == 1


def test_http_session_refresh_manual_directory_svg_marker_and_schema(session):
    result = api(session, "/api/new-molecules/add", json.dumps({"smiles": "CCF", "measurement": {
        "mode": "class", "measured_class": "inactive"}}).encode(), method="POST")
    assert result["status"] == 200 and result["data"]["status"] == "ok"
    assert Path(session.new_molecules_info()["directory"]).name == "manual"
    refreshed = api(session, "/api/session")["data"]["new_molecules"]
    assert refreshed["rows"][0]["measured_class"] == "inactive"
    assert 'data-neighbor-mismatch="true"' in session.files["map_with_new_molecules.svg"].read_text(encoding="utf-8")
    assert b"neighbor_evidence" in api(session, "/api/new-molecules/feedback")["data"]
    from chemgridmap.columns import describe
    assert len(describe(MANUAL_COLUMNS, "manual_entries.csv")) == len(MANUAL_COLUMNS)
    assert len(describe(FEEDBACK_COLUMNS, "neighbor_evidence_feedback.csv")) == len(FEEDBACK_COLUMNS)
    session.set_core_rule("most_frequent")
    detail, _ = session.inspect(molecule_id="M001")
    assert detail["molecule"]["contains_group_core_rule"] == "most_frequent"
    session.clear_query()
    assert not session.manual_entries and not session.measured_results


def test_numeric_manual_nonposition_fields_equal_single_file_and_clear_all_values(session):
    measurement = {"mode": "value", "value": 100, "unit": "nM", "batch": "B1"}
    add(session, measurement=measurement, label="trial")
    source = pd.DataFrame({"query_id": ["M001"], "smiles": ["CCF"], "label": ["trial"], "value": [100.], "batch": ["B1"]})
    single = insert_new_molecules(session.result, source, id_col="query_id", label_column="label",
                                  value_col="value", value_unit="nM", document_column="batch")
    ignore = {"source_row", "grid_row", "grid_col", "grid_x", "grid_y", "placement_distance", "user_median"}
    columns = [field for field in single.data if field not in ignore]
    pd.testing.assert_frame_equal(session.new_result.data[columns], single.data[columns])
    session.record_result({"query_id": "M001", "measurement": {"value": 1000, "unit": "nM", "batch": "B2"}})
    detail, records = session.inspect(molecule_id="M001")
    assert detail["measured_result"]["measured_class"] == "inactive"
    assert records.raw_activity_value.tolist() == [1000.]
    assert session.source_members(["B2"], "new")["map_molecule_count"] == 1
    session.record_result({"query_id": "M001", "measurement": {"mode": "class", "measured_class": "medium"}})
    assert session.new_result.records.empty
    session.clear_result({"query_id": "M001"})
    row = session.new_result.data.iloc[0]
    assert row.measurement_status == "not_measured" and pd.isna(row.activity_pchembl)
    assert pd.isna(row.pchembl_median) and pd.isna(row.n_records)


def test_removing_result_keeps_molecule_note_and_adding_preserves_uploaded_comparison(session):
    add(session, note="original molecule note")
    session.record_result({"query_id": "M001", "measurement": {"measured_class": "inactive", "note": "result note"}})
    assert session.measured_result_info("M001")["note"] == "result note"
    assert session.measured_result_info("M001")["result_source"] == "recorded"
    session.clear_result({"query_id": "M001"})
    assert session.new_molecules_info()["rows"][0]["note"] == "original molecule note"
    assert session.manual_entries[0]["note"] == "original molecule note"
    upload(session, b"id,smiles,value\nUP,CCCl,7\n", value_col="value")
    for table in (session.new_result.data, session.new_result.all_molecules):
        table.loc[:, "neighbor_evidence"] = "mixed"
        table.loc[:, "neighbor_comparison"] = ""
    original = session.new_result.data.copy(deep=True)
    add(session)
    pd.testing.assert_frame_equal(session.new_result.data.iloc[:1][original.columns].reset_index(drop=True), original, check_dtype=False)
    exported = pd.read_csv(session.files["new_molecules.csv"], keep_default_na=False)
    assert exported.iloc[0].neighbor_comparison == ""
    feedback = pd.read_csv(session.files["neighbor_evidence_feedback.csv"], keep_default_na=False)
    assert feedback.iloc[0].neighbor_comparison == "no_comparison"


def test_reference_preview_coordinates_disambiguate_colliding_query_ids(session):
    add(session, query_id="R1")
    preview = session.preview_new({"smiles": "CCO"})
    assert preview["reference_id"] == "R1"
    detail, _ = session.inspect(row=preview["grid_row"], col=preview["grid_col"])
    assert detail["molecule"]["molecule_id"] == "R1"
    assert detail["molecule"]["canonical_smiles"] == "CCO"


def test_numeric_append_keeps_every_original_uploaded_field_and_neighbor_check(session):
    upload(session)
    originals = {name: getattr(session.new_result, name).copy(deep=True)
                 for name in ("data", "all_molecules", "check", "neighbors")}
    add(session, measurement={"mode": "value", "value": 50, "unit": "nM", "batch": "manual assay"})
    for name, original in originals.items():
        current = getattr(session.new_result, name)
        pd.testing.assert_frame_equal(current.iloc[:len(original)][original.columns].reset_index(drop=True),
                                      original, check_dtype=False)
    for name in ("data", "all_molecules"):
        old_rows = getattr(session.new_result, name).iloc[:len(originals[name])]
        assert old_rows.measurement_status.eq("not_measured").all()
        assert old_rows.activity_class.eq("not_measured").all()
        for field in ("activity_pchembl", "n_records", "pchembl_median", "class_status"):
            if field in old_rows and field not in originals[name]:
                assert old_rows[field].isna().all()
    assert not session.measurement_rows().iloc[:2].result_source.any()


@pytest.mark.parametrize("unit,stored", [("p", "pActivity"), ("pActivity", "pActivity"),
    ("uM", "µM"), ("µM", "µM"), ("μM", "µM"), ("nM", "nM")])
def test_manual_unit_aliases_are_exported_with_display_names(unit, stored):
    result = normalize_measurement({"value": 1, "unit": unit}, lower=6, upper=7)
    assert result["measured_unit"] == stored
    assert result["activity_pchembl"] == pytest.approx(1 if stored == "pActivity" else 6 if stored == "µM" else 9)


@pytest.mark.parametrize("document_column", ["run_id", "batch"])
def test_manual_numeric_source_has_one_record_under_uploaded_source_column(session, document_column):
    raw = f"id,smiles,value,{document_column}\nUP,CCCl,7,uploaded assay\n".encode()
    upload(session, raw, value_col="value", document_column=document_column)
    add(session, measurement={"value": 50, "unit": "nM", "batch": "manual assay"})
    index = session._source_index("new")
    assert set(index) == {"uploaded assay", "manual assay"}
    assert index["uploaded assay"]["n_records"] == 1
    assert index["manual assay"]["n_records"] == 1
    assert index["manual assay"]["map_molecule_count"] == 1
    assert session.source_members(["manual assay"], "new")["map_molecule_count"] == 1
    detail, _ = session.inspect(molecule_id="M001")
    assert [(source["id"], source["n_records"]) for source in detail["new_molecule_sources"]] == [("manual assay", 1)]
    session.clear_result({"query_id": "M001"})
    assert set(session._source_index("new")) == {"uploaded assay"}
    upload(session, f"id,smiles,{document_column}\nN1,CCCl,original batch\n".encode(),
           document_column=document_column)
    session.record_result({"query_id": "N1", "measurement": {"value": 100, "unit": "nM", "batch": "retest"}})
    assert set(session._source_index("new")) == {"retest"}
    assert session._source_index("new")["retest"]["n_records"] == 1
    session.clear_result({"query_id": "N1"})
    assert set(session._source_index("new")) == {"original batch"}
    assert session._source_index("new")["original batch"]["n_records"] == 1
