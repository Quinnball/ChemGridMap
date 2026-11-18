"""U8d 的参照说明、提醒边界及导出契约。"""

import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from chemgridmap import build_grid_map, insert_new_molecules
from chemgridmap.app import AppSession
from chemgridmap.cli import main
from chemgridmap.queries import _insertion_warnings


@pytest.mark.parametrize("none,placed,median,warn", [
    (9, 10, .34, True),       # 恰好 90% 仍提醒。
    (89, 100, .34, False),
    (10, 10, .35, False),     # 恰好 0.35 不提醒。
    (10, 10, .349999, True),
    (9, 10, .36, False),
    (0, 0, None, False),
])
def test_structure_warning_boundaries(none, placed, median, warn):
    queries = pd.DataFrame({"neighbor_evidence": ["none"] * none + ["active"] * (placed - none),
                            "tanimoto": [median] * placed})
    report = _insertion_warnings(queries, pd.DataFrame())
    assert report["none_share"] == (none / placed if placed else 0)
    assert report["median_nearest_tanimoto"] == median
    assert report["warnings"] == ([{"code": "structural_mismatch", "placed": placed,
                                    "none": none, "median": median}] if warn else [])
    json.dumps(report, allow_nan=False)


@pytest.mark.parametrize("column", ["Target ChEMBL ID", "target_chembl_id", "TARGET CHEMBL ID"])
def test_target_warning_uses_detected_column_and_all_uploaded_rows(column):
    # 即使没有已放置的新分子，上传的不同靶点仍应提醒；空值和同靶点不提醒。
    queries = pd.DataFrame(columns=["neighbor_evidence", "tanimoto"])
    source = pd.DataFrame({column: [" CHEMBL205 ", "chembl240", "CHEMBL204", "CHEMBL240", "", None]})
    report = _insertion_warnings(queries, source, "chembl205")
    assert report["warnings"] == [{"code": "target_mismatch", "targets": ["CHEMBL204", "CHEMBL240"],
                                   "reference": "CHEMBL205"}]
    assert _insertion_warnings(queries, source.iloc[:1], "CHEMBL205")["warnings"] == []
    assert _insertion_warnings(queries, source, None)["warnings"] == []
    assert _insertion_warnings(queries, source.rename(columns={column: "other"}), "CHEMBL205")["warnings"] == []


@pytest.mark.parametrize("name_column,names,expected", [
    ("Target Name", ["Thrombin", "Carbonic anhydrase 2", None], "CHEMBL205 · Carbonic anhydrase 2 · IC50"),
    ("Target pref name", ["Thrombin", "Carbonic anhydrase 2", "Carbonic anhydrase 2"],
     "CHEMBL205 · Carbonic anhydrase 2 · IC50"),
    ("Target Name", ["Thrombin", None, " "], "CHEMBL205 · IC50"),
    ("unrecognized name", ["Thrombin", "Carbonic anhydrase 2", None], "CHEMBL205 · IC50"),
])
def test_reference_description_uses_only_selected_target_and_name_fallback(tmp_path, name_column, names, expected):
    session = AppSession(tmp_path)
    source = pd.DataFrame({"Target ChEMBL ID": ["CHEMBL204", "chembl205 ", "CHEMBL205"], name_column: names})
    session.load(source.to_csv(index=False).encode(), "reference.csv")
    session.result = SimpleNamespace(data=pd.DataFrame({"molecule_id": ["a", "b"]}))
    session.curation = SimpleNamespace(report={"parameters": {"target_id": "CHEMBL205"}})
    session.job = {"state": "done", "options": {"input_format": "chembl", "target_id": None, "activity_type": "IC50"}}
    assert session.snapshot()["reference_map"] == {"description": expected, "count": 2, "target_id": "CHEMBL205"}
    assert "reference_map" not in session.job


def test_generic_description_and_missing_target_column(tmp_path):
    session = AppSession(tmp_path)
    session.load(b"smiles,pActivity\nCCO,5\nCCN,7\n", "lab.csv")
    assert "reference_map" not in session.snapshot()
    session.result = SimpleNamespace(data=session.frame)
    session.job = {"state": "done", "options": {"input_format": "molecule", "value_col": "pActivity"}}
    assert session.snapshot()["reference_map"] == {"description": "lab.csv · pActivity", "count": 2, "target_id": None}
    session.job["options"] = {"input_format": "chembl", "target_id": "CHEMBL205", "activity_type": "IC50"}
    session.frame["Target Name"] = "An unverified name"
    assert session.snapshot()["reference_map"]["description"] == "CHEMBL205 · IC50"


@pytest.fixture
def reference_frame():
    return pd.DataFrame({"molecule_id": ["A", "B", "C", "D"],
                         "smiles": ["CCO", "CCN", "CCC", "c1ccccc1"], "value": [5, 7.2, 6.5, 8]})


def test_insert_report_manifest_and_reference_immutability(tmp_path, reference_frame):
    reference_frame["target_id"] = "CHEMBL205"
    reference = build_grid_map(reference_frame, value_col="value", projection="pca", grouping="none",
                               grid_refinement="none", output_dir=tmp_path / "ref", output_formats=["svg"],
                               render_detail="overview")
    original = reference.data.copy(deep=True)
    original_files = {path: Path(path).read_bytes() for path in reference.output_files.values()}
    inserted = insert_new_molecules(reference, pd.DataFrame({"smiles": ["[Xe]"], "Target ChEMBL ID": ["CHEMBL204"]}),
                                    output_dir=tmp_path / "insert")
    assert inserted.report["warnings"] == [
        {"code": "structural_mismatch", "placed": 1, "none": 1, "median": 0.0},
        {"code": "target_mismatch", "targets": ["CHEMBL204"], "reference": "CHEMBL205"}]
    manifest = json.loads(Path(inserted.output_files["new_molecules_manifest"]).read_text(encoding="utf-8"))
    for field in ("warnings", "none_share", "median_nearest_tanimoto"):
        assert manifest[field] == manifest["report"][field] == inserted.report[field]
    assert inserted.report["none_share"] == 1.0
    assert inserted.report["median_nearest_tanimoto"] == 0.0
    pd.testing.assert_frame_equal(reference.data, original)
    assert all(Path(path).read_bytes() == content for path, content in original_files.items())


def test_library_query_path_uses_selected_target_provenance(tmp_path, reference_frame):
    result = build_grid_map(reference_frame, value_col="value", projection="pca", grouping="none",
                            grid_refinement="none", output_dir=tmp_path, output_formats=["svg"],
                            render_detail="overview", provenance={"curation_parameters": {"target_id": "CHEMBL205"}},
                            query_data=pd.DataFrame({"smiles": ["CCF"], "target_chembl_id": ["CHEMBL204"]}))
    assert {"code": "target_mismatch", "targets": ["CHEMBL204"], "reference": "CHEMBL205"} in result.query_report["warnings"]


def test_session_and_cli_pass_explicit_target_without_target_column(tmp_path, reference_frame):
    source = reference_frame.rename(columns={"value": "pchembl_value", "molecule_id": "molecule_chembl_id"})
    raw = source.to_csv(index=False).encode()
    query = b"smiles,target_chembl_id\nCCF,CHEMBL204\n"
    session = AppSession(tmp_path / "app")
    session.load(raw, "reference.csv")
    session.start({"target_id": "CHEMBL205", "projection": "pca", "grouping": "none", "grid_refinement": "none",
                   "render_detail": "overview", "output_formats": ["svg"]}, background=False)
    assert session.job["state"] == "done", session.job
    session.load_query(query, "new.csv")
    snapshot = session.insert_new()
    expected = {"code": "target_mismatch", "targets": ["CHEMBL204"], "reference": "CHEMBL205"}
    assert expected in snapshot["report"]["new_molecules"]["warnings"]
    assert snapshot["reference_map"]["description"] == "CHEMBL205 · IC50"
    (tmp_path / "reference.csv").write_bytes(raw)
    (tmp_path / "new.csv").write_bytes(query)
    assert main([str(tmp_path / "reference.csv"), "--input-format", "chembl", "--target-id", "CHEMBL205",
                 "--projection", "pca", "--grouping", "none", "--grid-refinement", "none", "--render-detail", "overview",
                 "--output-formats", "svg", "--output-dir", str(tmp_path / "cli"),
                 "--query-file", str(tmp_path / "new.csv")]) == 0
    manifest = json.loads(next((tmp_path / "cli").rglob("new_molecules_manifest.json")).read_text(encoding="utf-8"))
    assert expected in manifest["warnings"]
