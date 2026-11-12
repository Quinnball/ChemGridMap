"""U9a 中英文 Node VM 契约；真实浏览器验收由沙箱外执行。"""
from pathlib import Path
import shutil
import subprocess

import pytest

from chemgridmap.app import ASSETS


@pytest.mark.parametrize("language", ["en", "zh"])
@pytest.mark.parametrize("scenario", ["sources", "source_highlight", "layers", "new_batches",
                                      "rules", "busy", "comparison", "stale_comparison", "language", "errors", "reset"])
def test_u9a_web_contract(scenario, language):
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node 不可用；保留沙箱外浏览器验收。")
    result = subprocess.run(
        [node, str(Path(__file__).with_name("u9a_web_contract.cjs")), str(ASSETS), scenario, language],
        text=True, capture_output=True, timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr
