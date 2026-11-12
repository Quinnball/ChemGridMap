"""U9c-R1 前端生产函数契约；原生浏览器排版/按键仍需浏览器验收。"""
import os
from pathlib import Path
import shutil
import subprocess

import pytest

from chemgridmap.app import ASSETS


@pytest.mark.parametrize("scenario", ["language", "redraw", "time", "list", "layout"])
def test_u9c_r1_web_contract(scenario):
    node = shutil.which("node")
    if node is None:
        pytest.skip("未安装 Node，不能执行 VM 契约。")
    result = subprocess.run(
        [node, str(Path(__file__).with_name("u9c_r1_web_contract.cjs")), str(ASSETS), scenario],
        env={**os.environ, "TZ": "Asia/Taipei"},
        text=True, capture_output=True, timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr
