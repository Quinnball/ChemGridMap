"""用外部输入副本核对 U8；不将模拟数据新增到软件包。"""

import hashlib
import json
import os
from pathlib import Path
from unittest.mock import patch
import zipfile

import numpy as np
import pandas as pd
import pytest

from chemgridmap.app import AppSession
from chemgridmap.columns import describe
from chemgridmap.queries import QUERY_CHECK_COLUMNS, QUERY_NEIGHBOR_COLUMNS


INPUTS = Path(os.environ.get("CHEMGRIDMAP_U8_INPUT_DIR", "/private/tmp/chemgridmap_u8/inputs"))


def file_hashes(directory):
    return {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in directory.iterdir() if path.is_file()}


def check_insertion_files(inserted, raw, reference_directory):
    expected = {"new_molecules_input.csv", "new_molecules_check.csv", "new_molecules.csv",
                "new_molecules_neighbors.csv", "map_with_new_molecules.svg",
                "new_molecules_manifest.json", "new_molecules_results.zip"}
    if inserted.report["measured"]:
        expected.add("new_molecules_excluded_records.csv")
    paths = [Path(value) for value in inserted.output_files.values()]
    directory = paths[0].parent
    assert {path.name for path in directory.iterdir()} == expected
    assert (directory / "new_molecules_input.csv").read_bytes() == raw
    assert inserted.check.columns.tolist() == QUERY_CHECK_COLUMNS
    assert inserted.neighbors.columns.tolist() == QUERY_NEIGHBOR_COLUMNS
    describe(inserted.check.columns, "new_molecules_check.csv")
    describe(inserted.neighbors.columns, "new_molecules_neighbors.csv")
    describe(inserted.all_molecules.columns, "new_molecules.csv")
    manifest = json.loads((directory / "new_molecules_manifest.json").read_text(encoding="utf-8"))
    json.dumps(manifest, allow_nan=False)
    assert manifest["input_sha256"] == hashlib.sha256(raw).hexdigest()
    assert Path(manifest["reference_run_directory"]) == reference_directory
    reference_manifest = reference_directory / "map_run_manifest.json"
    assert manifest["reference_manifest_sha256"] == hashlib.sha256(reference_manifest.read_bytes()).hexdigest()
    with zipfile.ZipFile(directory / "new_molecules_results.zip") as archive:
        assert set(archive.namelist()) == expected - {"new_molecules_results.zip"}
        for name in archive.namelist():
            assert archive.read(name) == (directory / name).read_bytes()
    return directory


@pytest.fixture(scope="module")
def real_session(tmp_path_factory):
    names = ["chembl205_ic50.csv", "demo_unmeasured_molecules_caii.csv",
             "sim_generated_molecules.csv", "sim_lab_assay_results.csv"]
    if not all((INPUTS / name).is_file() for name in names):
        pytest.skip("Set CHEMGRIDMAP_U8_INPUT_DIR to copied U8 inputs.")
    session = AppSession(tmp_path_factory.mktemp("u8_real_inputs"))
    session.load((INPUTS / names[0]).read_bytes(), names[0])
    session.start({"projection": "pca", "grid_refinement": "none",
                   "render_detail": "overview", "output_formats": ["svg"]}, background=False)
    assert session.job["state"] == "done", session.job
    assert len(session.result.data) == 889
    return session


@pytest.mark.parametrize("filename,id_column", [
    ("demo_unmeasured_molecules_caii.csv", "query_id"),
    ("sim_generated_molecules.csv", "gen_id"),
])
def test_real_candidates_preserve_all_reference_files_and_place_nearest(real_session, filename, id_column):
    session = real_session
    directory = Path(session.job["output_dir"])
    original_files = file_hashes(directory)
    reference = session.result.data.copy(deep=True)
    raw = (INPUTS / filename).read_bytes()
    session.load_query(raw, filename)
    # 禁止插入流程偷偷重算任何参照地图。
    with patch("chemgridmap.app.build_grid_map", side_effect=AssertionError("reference rebuilt")), \
         patch("chemgridmap.pipeline.project_2d", side_effect=AssertionError("reference reprojected")), \
         patch("chemgridmap.pipeline.compute_representation", side_effect=AssertionError("reference recomputed")):
        session.insert_new({"smiles_col": "smiles", "id_col": id_column})
    inserted = session.new_result
    pd.testing.assert_frame_equal(reference, session.result.data)
    assert file_hashes(directory) == original_files
    check_insertion_files(inserted, raw, directory)
    json.dumps(session.snapshot(), allow_nan=False)
    json.dumps(session.new_molecules_info(), allow_nan=False)
    json.dumps(session.points(), allow_nan=False)
    if filename.startswith("demo_"):
        assert {key: inserted.report[key] for key in ("mapped", "in_reference", "duplicate_query", "invalid_smiles")} == {
            "mapped": 45, "in_reference": 5, "duplicate_query": 2, "invalid_smiles": 3}
    else:
        assert inserted.report["input_rows"] == 313
        assert inserted.report["in_reference"] == 8
    assert sum(inserted.report["neighbor_evidence_counts"].values()) == len(inserted.data)
    # 独立穷举每一步所有空格，核对距离和行列并列规则。
    settings = session.result.insertion_parameters
    available = {(row, col) for row in range(settings["grid_rows"])
                 for col in range(settings["grid_cols"])} - set(zip(reference.grid_row, reference.grid_col))
    ref = reference.set_index("molecule_identity_key")
    for row in inserted.data.sort_values(["tanimoto", "query_id"], ascending=[False, True]).itertuples():
        target = ref.loc[row.nearest_reference_id]
        nearest = min(available, key=lambda cell: ((cell[0]-target.grid_row)**2 +
                                                   (cell[1]-target.grid_col)**2, cell[0], cell[1]))
        assert (row.grid_row, row.grid_col) == nearest
        assert row.placement_distance == pytest.approx(np.hypot(nearest[0]-target.grid_row,
                                                              nearest[1]-target.grid_col))
        available.remove(nearest)


def test_simulated_measurements_units_controls_and_offsets(real_session):
    session = real_session
    raw = (INPUTS / "sim_lab_assay_results.csv").read_bytes()
    directory = Path(session.job["output_dir"])
    before = file_hashes(directory)
    session.load_query(raw, "sim_lab_assay_results.csv")
    session.insert_new({"smiles_col": "smiles", "id_col": "compound_id", "activity_col": "ic50",
                        "unit_col": "unit", "relation_col": "relation", "document_col": "run_id",
                        "label_col": "design_id"})
    inserted = session.new_result
    assert file_hashes(directory) == before
    check_insertion_files(inserted, raw, directory)
    records = inserted.records
    assert len(records) == 104
    expected = -np.log10(pd.to_numeric(records.ic50) * records.unit.map({"nM": 1e-9, "uM": 1e-6, "µM": 1e-6}))
    np.testing.assert_allclose(records.activity_pchembl_raw, expected, rtol=0, atol=1e-12)
    assert inserted.report["excluded_activity_records"] == 30
    assert inserted.report["excluded_reasons"]["non_exact_relation"] == 28
    assert inserted.report["excluded_reasons"]["invalid_activity_value"] == 2
    controls = [("CTRL-AAZ", "CHEMBL20", 0.0383953867693414),
                ("CTRL-MZA", "CHEMBL19", 0.16401760455646652),
                ("CTRL-EZA", "CHEMBL18", -0.08240688999363854)]
    rows = inserted.all_molecules.set_index("query_id")
    for own, reference_id, difference in controls:
        assert rows.loc[own, "status"] == "in_reference"
        assert rows.loc[own, "reference_id"] == reference_id
        assert rows.loc[own, "median_difference"] == pytest.approx(difference, abs=1e-12)
        detail, _ = session.inspect(molecule_id=own)
        json.dumps(detail, allow_nan=False)
    pairs = {(row["document_a"], row["document_b"]): row for row in inserted.report["document_pair_offsets"]}
    expected_offsets = {("R1", "R2"): (-0.025295339700091013, False),
                        ("R1", "R3"): (1.1754694791501894, True),
                        ("R2", "R3"): (1.1029064260568795, True)}
    assert pairs.keys() == expected_offsets.keys()
    for pair, (offset, flagged) in expected_offsets.items():
        assert pairs[pair]["median_offset"] == pytest.approx(offset, abs=1e-12)
        assert pairs[pair]["systematic_offset"] is flagged
    json.dumps(session.new_molecules_info(), allow_nan=False)
