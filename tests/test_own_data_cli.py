"""Public CLI paths for own measurements and unmeasured structures."""
import json
import zipfile

import numpy as np
import pandas as pd
import pytest

from chemgridmap.cli import build_parser, main
from chemgridmap.columns import describe
from chemgridmap.inspection import inspect_entry


@pytest.mark.parametrize("arguments", [
    ["--value-unit", "p", "--unit-column", "unit"],
    ["--assumed-sigma", "0"],
    ["--assumed-sigma", "nan"],
    ["--assumed-sigma", "inf"],
])
def test_new_cli_options_reject_conflicts_and_invalid_sigma(arguments):
    with pytest.raises(SystemExit):
        build_parser().parse_args(arguments)


def test_molecule_only_cli_rejects_assumed_sigma(tmp_path, capsys):
    source = tmp_path / "molecules.csv"
    pd.DataFrame({"smiles": ["CCO", "CCN"], "value": [6, 7]}).to_csv(source, index=False)
    with pytest.raises(SystemExit):
        main([str(source), "--value-col", "value", "--assumed-sigma", "0.64"])
    assert "record-table" in capsys.readouterr().err


def test_units_assumed_sigma_queries_cli_audit(tmp_path):
    source = tmp_path / "records.csv"
    pd.DataFrame({
        "compound": ["A", "A", "B", "C", "D", "E", "F", "G"],
        "smiles": ["CCO", "OCC", "CCN", "CCC", "COC", "CCF", "CCCl", "CCBr"],
        "IC50": [1000, 10, 100, 1, 0, "bad", 5, 10],
        "unit": ["nM", "nM", "nM", " µM ", "nM", "nM", "NM", "nM"],
        "relation": ["=", "=", "=", "=", "=", "=", "=", ">"],
    }).to_csv(source, index=False)
    query = tmp_path / "queries.tsv"
    pd.DataFrame({"query_id": ["001", "002", "003", "004"],
                  "SMILES": ["CCCO", "OCCC", "invalid", "CCN"]}).to_csv(query, sep="\t", index=False)
    output = tmp_path / "out"
    assert main([
        str(source), "--compound-column", "compound", "--value-column", "IC50",
        "--unit-column", "unit", "--relation-column", "relation", "--assumed-sigma", "0.64",
        "--query-file", str(query), "--query-id-column", "query_id",
        "--projection", "pca", "--grouping", "none", "--render-detail", "overview",
        "--output-formats", "svg", "--map-color-modes", "activity", "class_support",
        "--output-dir", str(output), "--name", "own",
    ]) == 0
    coordinates = pd.read_csv(output / "own_grid_coordinates.csv")
    records = pd.read_csv(output / "own_retained_records.csv")
    manifest = json.loads((output / "own_run_manifest.json").read_text(encoding="utf-8"))
    report = json.loads((output / "own_curation_report.json").read_text(encoding="utf-8"))
    fit = manifest["activity_confidence"]
    assert fit["method"] == "assumed_sigma"
    assert fit["assumed_sigma"] == 0.64
    assert fit["evidence_unit_source"] == "one unit per compound"
    assert fit["document_offset_detection"] == "not applicable"
    assert report["excluded_reasons"]["nonpositive_activity_value"] == 1
    assert report["excluded_reasons"]["invalid_activity_value"] == 1
    assert report["excluded_reasons"]["unrecognized_activity_unit"] == 1
    assert report["excluded_reasons"]["non_exact_relation"] == 1
    np.testing.assert_allclose(records.activity_pchembl_raw, [6, 8, 7, 6], atol=1e-12)
    assert {"raw_activity_value", "raw_activity_unit"} <= set(records.columns)
    _, _, inspection = inspect_entry(coordinates, records, molecule_id="A", activity_confidence=fit)
    assert inspection["median_verified"] is True
    unmeasured = coordinates[coordinates.measurement_status.eq("not_measured")]
    assert len(unmeasured) == 1
    assert unmeasured.class_confidence.isna().all()
    assert unmeasured.activity_pchembl.isna().all()
    check = pd.read_csv(output / "own_query_check.csv", dtype={"query_id": str})
    assert check.query_id.tolist() == ["001", "002", "003", "004"]
    assert set(check.status) == {"mapped", "duplicate_query", "invalid_smiles", "in_reference"}
    neighbors = pd.read_csv(output / "own_query_neighbors.csv")
    describe(check.columns, "query_check.csv")
    describe(neighbors.columns, "query_neighbors.csv")
    with zipfile.ZipFile(output / "own_results.zip") as bundle:
        assert {"own_query_check.csv", "own_query_neighbors.csv", "own_retained_records.csv"} <= set(bundle.namelist())


def test_chembl_assumed_sigma_without_document_column(tmp_path):
    source = tmp_path / "chembl.csv"
    pd.DataFrame({
        "molecule_chembl_id": ["CHEMBL1", "CHEMBL2", "CHEMBL3"],
        "canonical_smiles": ["CCO", "CCN", "CCC"],
        "pchembl_value": [6.5, 7.2, 5.5],
        "standard_type": ["IC50"] * 3, "standard_relation": ["="] * 3,
        "standard_units": ["nM"] * 3, "target_chembl_id": ["CHEMBL205"] * 3,
    }).to_csv(source, index=False)
    output = tmp_path / "out"
    assert main([str(source), "--input-format", "chembl", "--assumed-sigma", "0.64",
                 "--projection", "pca", "--grouping", "none", "--output-formats", "svg",
                 "--render-detail", "overview", "--output-dir", str(output)]) == 0
    manifest = json.loads((output / "grid_map_run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["activity_confidence"]["method"] == "assumed_sigma"
    assert manifest["activity_confidence"]["document_offset_detection"] == "not applicable"


def test_single_molecule_units_preserve_exclusions_in_manifest(tmp_path):
    source = tmp_path / "molecules.csv"
    pd.DataFrame({"smiles": ["CCO", "CCN", "CCC", "CCF"],
                  "IC50": [10, 100, 0, "bad"]}).to_csv(source, index=False)
    assert main([str(source), "--value-column", "IC50", "--value-unit", "nM",
                 "--projection", "pca", "--grouping", "none", "--render-detail", "overview",
                 "--output-formats", "svg", "--output-dir", str(tmp_path / "out")]) == 0
    manifest = json.loads((tmp_path / "out/grid_map_run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["parameters"]["value_unit"] == "nM"
    assert manifest["preprocessing"]["excluded_reasons"]["nonpositive_activity_value"] == 1
    assert manifest["preprocessing"]["excluded_reasons"]["invalid_activity_value"] == 1
    coordinates = pd.read_csv(tmp_path / "out/grid_map_grid_coordinates.csv")
    np.testing.assert_allclose(coordinates.IC50, [8, 7], atol=1e-12)
    assert coordinates.raw_activity_value.tolist() == [10, 100]
