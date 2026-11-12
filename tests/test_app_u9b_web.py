"""U9b 双语浏览器逻辑契约；真实浏览器排版由外部验收。"""
from pathlib import Path
import shutil
import subprocess

import pytest
from chemgridmap.app import ASSETS


@pytest.mark.parametrize("scenario", ["errors", "controls", "preview_cache", "preview", "validation", "add", "measured", "list_audit", "replacement", "markers", "mutation_highlight", "cropped_patch", "patch"])
def test_u9b_web_contract(scenario):
    node = shutil.which("node")
    if node is None:
        pytest.skip("未安装 Node；不能执行 VM 契约。")
    result = subprocess.run([node, str(Path(__file__).with_name("u9b_web_contract.cjs")), str(ASSETS), scenario],
                            text=True, capture_output=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr


def test_u9b_web_uses_actual_session_payload(tmp_path):
    """以真实新增、已测上传和反馈结果，校验前后端字段与本地化。"""
    import json
    from chemgridmap.app import AppSession

    node = shutil.which("node")
    if node is None:
        pytest.skip("未安装 Node。")
    app = AppSession(tmp_path / "app")
    app.load(b"molecule_id,smiles,value\nR1,CCO,5\nR2,CCN,7.2\nR3,CCC,6.5\nR4,c1ccccc1,8\n", "reference.csv")
    app.start({"input_format": "molecule", "value_col": "value", "projection": "pca",
               "grouping": "none", "grid_refinement": "none", "render_detail": "overview",
               "output_formats": ["svg"]}, background=False)
    assert app.job["state"] == "done", app.job
    app.load_query(b"id,smiles,value\nN1,CCCl,8\nN2,CCBr,\n", "uploaded.csv")
    app.insert_new({"smiles_col": "smiles", "id_col": "id", "activity_col": "value", "unit": "p"})
    app.add_new({"smiles": "CCF", "measurement": {"mode": "value", "value": 15, "unit": "nM", "batch": "B1", "note": "manual"}})
    app.record_result({"query_id": "N2", "measurement": {"mode": "class", "measured_class": "inactive", "batch": "B2", "note": "feedback"}})
    payload = {"job": app.job, "new": app.new_molecules_info(), "points": app.points(),
               "details": [app.inspect(molecule_id=identity)[0] for identity in ("N1", "N2", "M001")]}
    saved = tmp_path / "payload.json"
    saved.write_text(json.dumps(payload, allow_nan=False), encoding="utf-8")
    result = subprocess.run([node, str(Path(__file__).with_name("u9b_web_contract.cjs")), str(ASSETS), "payload", str(saved)],
                            text=True, capture_output=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
