"""U8 的 Node VM 交互契约；不启动真实浏览器。"""
from pathlib import Path
import shutil
import subprocess

import pytest
from chemgridmap.app import ASSETS


@pytest.mark.parametrize("scenario", ["dialogs", "workflow", "list", "csv", "layers", "inspector", "audit", "projection", "outlines", "u8b_formatting", "u8b_markers", "u8b_summary_exports", "u8c_match_notes", "u8c_neighbor_counts"])
def test_u8_web_contract(scenario):
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node 不可用，真实浏览器验收仍须单独执行。")
    result = subprocess.run([node, str(Path(__file__).with_name("u8_web_contract.cjs")), str(ASSETS), scenario],
                            text=True, capture_output=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr


def test_u8_web_consumes_actual_appsession_payload(tmp_path):
    """实际 AppSession JSON 驱动 Node DOM，捕捉前后端字段偏差。"""
    import json
    from chemgridmap.app import AppSession

    node = shutil.which("node")
    if node is None:
        pytest.skip("Node 不可用。")
    app = AppSession(tmp_path / "app")
    app.load(b"molecule_id,smiles,value\nR1,CCO,5\nR2,CCN,7.2\nR3,CCC,6.5\nR4,c1ccccc1,8\n", "reference.csv")
    app.start({"input_format": "molecule", "value_col": "value", "projection": "pca",
               "grouping": "none", "grid_refinement": "none", "render_detail": "overview",
               "output_formats": ["svg"]}, background=False)
    assert app.job["state"] == "done", app.job
    app.load_query(b"id,smiles,value,units,relation,batch,label\nN1,CCF,10,nM,=,R1,round1\nN1,CCF,0.02,uM,=,R2,round1\nN2,CCCl,300,nM,=,R1,round2\nCTRL,OCC,100,nM,=,R1,control\nBAD,bad_smiles,1,nM,=,R1,bad\n", "actual.csv")
    app.insert_new({"smiles_col": "smiles", "id_col": "id", "activity_col": "value",
                    "unit_col": "units", "relation_col": "relation", "document_col": "batch", "label_col": "label"})
    detail, _ = app.inspect(molecule_id="N1")
    payload = tmp_path / "actual-payload.json"
    payload.write_text(json.dumps({"job": app.job, "new": app.new_molecules_info(), "points": app.points(), "detail": detail}, allow_nan=False), encoding="utf-8")
    result = subprocess.run([node, str(Path(__file__).with_name("u8_web_contract.cjs")), str(ASSETS), "payload", str(payload)],
                            text=True, capture_output=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
