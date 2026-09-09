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
    summary = json.loads((tmp_path / "audit/audit_summary.json").read_text())
    assert summary["median_verified"] and summary["record_count_verified"]
    svg = (tmp_path / "audit/selected_molecule.svg").read_text()
    assert "data-conflict-marker" in svg and "<image" not in svg
    assert len(summary["inputs"]["coordinates"]["sha256"]) == 64
