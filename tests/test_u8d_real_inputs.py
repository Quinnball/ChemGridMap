"""用临时目录中的真实输入副本验证参照图错配提醒。"""

import hashlib
import json
import os
from pathlib import Path

import pytest

from chemgridmap.app import AppSession


INPUTS = Path(os.environ.get("CHEMGRIDMAP_U8D_INPUT_DIR", "/private/tmp/chemgridmap_u8d/inputs"))


@pytest.fixture(scope="module")
def reference_sessions(tmp_path_factory):
    sessions = {}

    def get_session(filename, target_id):
        source = INPUTS / filename
        if not source.is_file():
            pytest.skip("Set CHEMGRIDMAP_U8D_INPUT_DIR to copied U8d inputs.")
        if filename not in sessions:
            session = AppSession(tmp_path_factory.mktemp("u8d_" + target_id))
            session.load(source.read_bytes(), filename)
            session.start({"target_id": target_id, "projection": "pca",
                           "grid_refinement": "none", "render_detail": "overview",
                           "output_formats": ["svg"]}, background=False)
            assert session.job["state"] == "done", session.job
            sessions[filename] = session
        return sessions[filename]

    return get_session


@pytest.mark.parametrize("reference_file,target_id,query_file,id_column,expected_warning", [
    ("chembl205_ic50.csv", "CHEMBL205", "sim_generated_molecules.csv", "gen_id", False),
    ("chembl205_ic50.csv", "CHEMBL205", "demo_unmeasured_molecules_caii.csv", "query_id", False),
    ("Thrombin_CHEMBL204_IC50.csv", "CHEMBL204", "demo_unmeasured_molecules_caii.csv", "query_id", True),
])
def test_real_candidates_warn_only_for_mismatched_reference(
        reference_sessions, reference_file, target_id, query_file, id_column, expected_warning):
    source = INPUTS / query_file
    if not source.is_file():
        pytest.skip("Set CHEMGRIDMAP_U8D_INPUT_DIR to copied U8d inputs.")
    session = reference_sessions(reference_file, target_id)
    reference_directory = Path(session.job["output_dir"])
    before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
              for p in reference_directory.iterdir() if p.is_file()}
    session.load_query(source.read_bytes(), query_file)
    session.insert_new({"smiles_col": "smiles", "id_col": id_column})
    inserted = session.new_result
    report = inserted.report
    placed = len(inserted.data)
    none = int(inserted.data.neighbor_evidence.eq("none").sum())
    median = float(inserted.data.tanimoto.median())
    assert placed > 0
    assert report["none_share"] == pytest.approx(none / placed)
    assert report["median_nearest_tanimoto"] == pytest.approx(median)
    if expected_warning:
        assert none / placed >= .9 and median < .35
        assert len(report["warnings"]) == 1
        warning = report["warnings"][0]
        assert {key: warning[key] for key in ("code", "placed", "none")} == {
            "code": "structural_mismatch", "placed": placed, "none": none}
        assert warning["median"] == pytest.approx(median)
    else:
        assert report["warnings"] == []
    output_directory = Path(next(iter(inserted.output_files.values()))).parent
    manifest = json.loads((output_directory / "new_molecules_manifest.json").read_text(encoding="utf-8"))
    for key in ("warnings", "none_share", "median_nearest_tanimoto"):
        assert manifest[key] == report[key]
        assert manifest["report"][key] == report[key]
    after = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
             for p in reference_directory.iterdir() if p.is_file()}
    assert before == after
    json.dumps(session.new_molecules_info(), allow_nan=False)
