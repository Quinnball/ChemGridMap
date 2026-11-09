"""U8d 双语 Node VM 契约；不启动浏览器。"""
import json
from pathlib import Path
import shutil
import subprocess

import pytest

from chemgridmap.app import ASSETS, AppSession


def run_contract(scenario, language, payload=None):
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node 不可用；真实浏览器验收由沙箱外另行完成。")
    command = [node, str(Path(__file__).with_name("u8d_web_contract.cjs")), str(ASSETS), scenario, language]
    if payload is not None:
        command.append(str(payload))
    result = subprocess.run(command, text=True, capture_output=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize("language", ["en", "zh"])
@pytest.mark.parametrize("scenario", ["workflow", "zeroPlaced", "counts", "legends", "reference", "warnings"])
def test_u8d_web_contract(scenario, language):
    run_contract(scenario, language)


def test_u8d_web_consumes_snapshot_and_structured_warnings(tmp_path):
    app = AppSession(tmp_path / "app")
    app.load(("molecule_chembl_id,canonical_smiles,pchembl_value,target_chembl_id,target_pref_name,standard_type,standard_relation\n"
              "R1,CCO,5,CHEMBL205,Carbonic anhydrase 2,IC50,=\n"
              "R2,CCN,7.2,CHEMBL205,Carbonic anhydrase 2,IC50,=\n"
              "R3,CCC,6.5,CHEMBL205,Carbonic anhydrase 2,IC50,=\n"
              "R4,c1ccccc1,8,CHEMBL205,Carbonic anhydrase 2,IC50,=\n").encode(), "reference.csv")
    app.start({"input_format": "chembl", "target_id": "CHEMBL205", "activity_type": "IC50",
               "projection": "pca", "grouping": "none", "grid_refinement": "none",
               "render_detail": "overview", "output_formats": ["svg"]}, background=False)
    assert app.job["state"] == "done", app.job
    app.load_query(b"id,smiles,target_chembl_id\nS1,C[Si](C)(C)C,CHEMBL204\n", "silicon.csv")
    app.insert_new({"smiles_col": "smiles", "id_col": "id"})
    payload = tmp_path / "snapshot.json"
    payload.write_text(json.dumps({"job": app.snapshot(), "new": app.new_molecules_info(),
                                   "points": app.points()}, allow_nan=False), encoding="utf-8")
    for language in ["en", "zh"]:
        run_contract("payload", language, payload)
