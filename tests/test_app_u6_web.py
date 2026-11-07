"""U6 结构匹配界面契约；使用 Node VM，不启动真实浏览器。"""

from pathlib import Path
import shutil
import subprocess

import pytest

from chemgridmap.app import ASSETS


@pytest.mark.parametrize("scenario", ["matched", "duplicates", "audit", "limit", "no_sources", "language", "invalid"])
def test_u6_web_contract(scenario):
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node is unavailable; browser acceptance remains a separate check.")
    script = Path(__file__).with_name("u6_web_contract.cjs")
    result = subprocess.run([node, str(script), str(ASSETS), scenario],
                            text=True, capture_output=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
