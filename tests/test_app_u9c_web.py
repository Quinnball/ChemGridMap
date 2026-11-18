"""U9c 双语 Node VM 契约，不替代真实浏览器排版验收。"""
from pathlib import Path
import shutil
import subprocess

import pytest
from chemgridmap.app import ASSETS


@pytest.mark.parametrize("scenario", ["progress", "list", "restore", "markers", "controls"])
def test_u9c_web_contract(scenario):
    node = shutil.which("node")
    if node is None:
        pytest.skip("未安装 Node，不能执行 VM 契约。")
    result = subprocess.run(
        [node, str(Path(__file__).with_name("u9c_web_contract.cjs")), str(ASSETS), scenario],
        text=True, capture_output=True, timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr
