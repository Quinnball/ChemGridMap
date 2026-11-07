"""U7 来源查找与 U6 收尾契约；使用 Node VM，不启动真实浏览器。"""

from pathlib import Path
import shutil
import subprocess

import pytest

from chemgridmap.app import ASSETS


@pytest.mark.parametrize("scenario", ["sources", "audit", "controls", "errors", "language", "u6_layout"])
def test_u7_web_contract(scenario):
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node is unavailable; browser acceptance remains a separate check.")
    script = Path(__file__).with_name("u7_web_contract.cjs")
    result = subprocess.run([node, str(script), str(ASSETS), scenario],
                            text=True, capture_output=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
