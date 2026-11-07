"""U5 来源文献与中文排版契约；使用 Node VM，不启动真实浏览器。"""

from pathlib import Path
import shutil
import subprocess

import pytest

from chemgridmap.app import ASSETS


@pytest.mark.parametrize("scenario", ["sources", "references", "coverage", "punctuation"])
def test_u5_web_contract(scenario):
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node is unavailable; browser acceptance remains a separate check.")
    script = Path(__file__).with_name("u5_web_contract.cjs")
    result = subprocess.run([node, str(script), str(ASSETS), scenario],
                            text=True, capture_output=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
