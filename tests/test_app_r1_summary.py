"""Browser record counts describe the uploaded records, before molecule aggregation."""

import json

import pandas as pd
import pytest

from chemgridmap.app import AppSession


@pytest.mark.parametrize("value_unit", ["p", "nM"])
@pytest.mark.parametrize("with_queries", [False, True])
def test_generic_summary_counts_input_records_retained_records_and_wide_molecules(
        tmp_path, value_unit, with_queries):
    # C1 has three wide-range records, C2 two narrow records, and C3 one record.
    # The seventh source record is censored and must not enter the retained count.
    frame = pd.DataFrame({
        "compound": ["C1", "C1", "C1", "C2", "C2", "C3", "C4"],
        "smiles": ["CCO", "CCO", "CCO", "CCN", "CCN", "CCC", "CCCC"],
        "activity": [5., 7., 8., 6.2, 6.3, 7.5, 9.],
        "relation": ["=", "=", "=", "=", "=", "=", ">"],
    })
    if value_unit == "nM":
        frame["activity"] = 10. ** (9. - frame["activity"])
    session = AppSession(tmp_path)
    session.load(frame.to_csv(index=False).encode(), "own_records.csv")
    options = {
        "input_format": "molecule", "compound_column": "compound",
        "smiles_col": "smiles", "value_col": "activity", "relation_column": "relation",
        "projection": "pca", "grouping": "none", "output_formats": ["svg"],
    }
    if value_unit != "p":
        options.update(value_unit=value_unit, assumed_sigma=.64)
    if with_queries:
        session.load_query(b"query_id,smiles\nQ1,CCF\nQ2,CCCl\n", "queries.csv")
        options.update(query_smiles_col="smiles", query_id_col="query_id")
    session.start(options, background=False)
    assert session.job["state"] == "done", session.job
    report = session.job["report"]
    assert report["input_rows"] == 7
    assert report["retained_activity_records"] == 6
    assert report["repeated_conflicted_molecules"] == 1
    assert session.job["molecules"] == (5 if with_queries else 3)
    if with_queries:
        assert report["query_molecules"]["mapped"] == 2
    else:
        assert "query_molecules" not in report
    # Updating the in-memory job summary must not overwrite map preprocessing or exports.
    assert session.result.preprocessing["input_rows"] == 3
    manifest = json.loads(session.files["map_run_manifest.json"].read_text(encoding="utf-8"))
    exported_report = manifest["source"]["generic_activity"]
    assert exported_report["input_rows"] == 7
    assert exported_report["retained_activity_records"] == 6
    assert "repeated_conflicted_molecules" not in exported_report


def test_generic_summary_reports_zero_wide_ranges(tmp_path):
    session = AppSession(tmp_path)
    session.load(b"compound,smiles,activity\nC1,CCO,5\nC2,CCN,7\nC3,CCC,8\n", "own_records.csv")
    session.start({
        "input_format": "molecule", "compound_column": "compound", "value_col": "activity",
        "projection": "pca", "grouping": "none", "output_formats": ["svg"],
    }, background=False)
    assert session.job["state"] == "done", session.job
    report = session.job["report"]
    assert report["input_rows"] == report["retained_activity_records"] == 3
    assert report["repeated_conflicted_molecules"] == 0


def test_chembl_summary_preserves_the_existing_curation_report(tmp_path):
    session = AppSession(tmp_path)
    session.load(
        b"Molecule ChEMBL ID,Smiles,pChEMBL Value\n"
        b"CHEMBL1,CCO,5\nCHEMBL1,CCO,8\nCHEMBL2,CCN,7\nCHEMBL3,CCC,6\n",
        "chembl.csv",
    )
    session.start({"projection": "pca", "grouping": "none", "output_formats": ["svg"]},
                  background=False)
    assert session.job["state"] == "done", session.job
    assert session.job["report"] is session.curation.report
    assert session.job["report"]["input_rows"] == 4
    assert session.job["report"]["retained_activity_records"] == 4
    assert session.job["report"]["repeated_conflicted_molecules"] == 1
