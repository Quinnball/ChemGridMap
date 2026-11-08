"""U8b 双语标签及模态提示框契约；真实浏览器遮挡验收在沙箱外完成。"""

from pathlib import Path
import shutil
import subprocess

import pytest

from chemgridmap.app import ASSETS


@pytest.mark.parametrize("scenario", ["labels", "dialogHelp"])
def test_u8b_dialog_contract(scenario):
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node is unavailable; browser acceptance remains a separate check.")
    script = Path(__file__).with_name("u8b_dialog_contract.cjs")
    result = subprocess.run([node, str(script), str(ASSETS), scenario],
                            text=True, capture_output=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
