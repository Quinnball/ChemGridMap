import json
from pathlib import Path

import pandas as pd
import pytest

from chemgridmap import curate_chembl_activity_data, inspect_entry
from chemgridmap.cli import main


@pytest.fixture
def case():
    raw = pd.DataFrame({"Smiles": ["CCO"] * 3 + ["CCN"],
                        "Molecule ChEMBL ID": ["CHEMBL20"] * 3 + ["CHEMBL200"],
                        "pChEMBL Value": [5.2, 7.2, 8.0, 5.0]})
    curated = curate_chembl_activity_data(raw)
    grid = curated.molecules.copy()
    grid["grid_row"] = [2, 3]
    grid["grid_col"] = [4, 5]
    return grid, curated.retained_records


def test_id_and_cell_recover_the_same_measurements(case):
    grid, records = case
    row, evidence, summary = inspect_entry(grid, records, molecule_id="chembl20")
    _, by_cell, _ = inspect_entry(grid, records, cell=(int(row.iloc[0].grid_row), int(row.iloc[0].grid_col)))
    pd.testing.assert_frame_equal(evidence, by_cell)
    assert summary["n_records"] == 3
    assert summary["median_pchembl"] == 7.2
    assert summary["is_conflicted"]
    assert summary["n_assays"] is None
    assert summary["offset_pair_touched"] is None
    assert summary["status_reason"] is None


def test_inspector_reports_kth_fingerprint_neighbor(case):
    grid, records = case
    grid = grid.copy()
    grid["neighbor_kth_similarity"] = [0.42, 0.51]
    _, _, summary = inspect_entry(grid, records, molecule_id="CHEMBL20", k=10)
    assert summary["neighbor_kth_label"] == "10th fingerprint neighbor: Tanimoto 0.510"


@pytest.mark.parametrize(("offset_display", "expected_rule"), [
    ("uncertain", "This entry has records from DOC_A, so it is shown as uncertain."),
    ("none", "This entry has records from DOC_A; the display is unchanged because the offset rule is off."),
])
def test_inspector_reports_each_touching_offset_pair_with_present_documents(
        offset_display, expected_rule):
    records = pd.DataFrame({
        "source_row": [0],
        "molecule_identity_key": ["CHEMBL4"],
        "activity_pchembl_raw": [6.5],
        "document_chembl_id": ["DOC_A"],
    })
    grid = pd.DataFrame({
        "molecule_identity_key": ["CHEMBL4"], "source_rows": ["0"],
        "n_records": [1], "activity_pchembl": [6.5],
        "grid_row": [0], "grid_col": [0],
        "offset_pair_touched": [True], "status_reason": ["offset_pair"],
    })
    activity_confidence = {
        "confidence_threshold": 0.8,
        "offset_display": offset_display,
        "document_offsets": [{
            "document_a": "DOC_A", "document_b": "DOC_B",
            "median_offset": 3.0, "shared_entries": 3, "touched_entry_count": 7,
        }],
    }

    _, _, summary = inspect_entry(
        grid, records, molecule_id="CHEMBL4",
        activity_confidence=activity_confidence,
    )
    assert len(summary["document_offsets"]) == 1
    assert summary["status_reason"] == "offset_pair"
    assert summary["confidence_threshold"] == 0.8
    assert (
        "Offset document pair: DOC_A vs DOC_B, median difference 3 log units, "
        "n shared 3; the pair touches 7 entries. "
        + expected_rule
    ) in summary["warnings"]

    _, _, unavailable = inspect_entry(
        grid, records, molecule_id="CHEMBL4"
    )
    assert unavailable["document_offsets"] == []
    assert "Offset document pair details are unavailable for this map." in unavailable["warnings"]


def test_id_is_exact_not_substring(case):
    with pytest.raises(ValueError, match="0 entries"):
        inspect_entry(*case, molecule_id="CHEMBL2")


def test_smiles_identity_keys_remain_case_sensitive():
    raw = pd.DataFrame({"Smiles": ["c1ccccc1", "C1CCCCC1"], "pChEMBL Value": [5.0, 8.0]})
    curated = curate_chembl_activity_data(raw)
    grid = curated.molecules.copy()
    grid["grid_row"], grid["grid_col"] = [0, 1], [0, 0]
    assert grid.molecule_identity_key.str.upper().nunique() == 1
    for identity in grid.molecule_identity_key:
        chosen, records, summary = inspect_entry(grid, curated.retained_records, molecule_id=identity)
        assert len(chosen) == 1 and len(records) == 1
        assert summary["molecule_identity_key"] == identity


def test_empty_cell_is_not_replaced_by_nearest_molecule(case):
    with pytest.raises(ValueError, match="0 entries"):
        inspect_entry(*case, cell=(0, 0))


def test_missing_source_row_is_rejected(case):
    grid, records = case
    with pytest.raises(ValueError, match="counts"):
        inspect_entry(grid, records.drop(index=0), molecule_id="CHEMBL20")


def test_changed_measurements_are_rejected(case):
    grid, records = case
    records.loc[1, "activity_pchembl_raw"] = 6.0
    with pytest.raises(ValueError, match="activity"):
        inspect_entry(grid, records, molecule_id="CHEMBL20")


def test_different_record_identity_is_rejected(case):
    grid, records = case
    records.loc[0, "molecule_identity_key"] = "different"
    with pytest.raises(ValueError, match="identity"):
        inspect_entry(grid, records, molecule_id="CHEMBL20")


def test_inspect_cli_exports_auditable_native_svg(case, tmp_path):
    grid, records = case
    grid.to_csv(tmp_path / "grid.csv", index=False)
    records.to_csv(tmp_path / "records.csv", index=False)
    assert main(["inspect", str(tmp_path / "grid.csv"), "--records", str(tmp_path / "records.csv"),
                 "--molecule-id", "CHEMBL20", "--output-dir", str(tmp_path / "audit")]) == 0
    summary = json.loads(
        (tmp_path / "audit/audit_summary.json").read_text(encoding="utf-8")
    )
    assert summary["median_verified"] and summary["record_count_verified"]
    svg = (tmp_path / "audit/selected_molecule.svg").read_text(encoding="utf-8")
    assert "data-conflict-marker" in svg and "<image" not in svg
    assert len(summary["inputs"]["coordinates"]["sha256"]) == 64
