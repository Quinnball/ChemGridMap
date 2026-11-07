"""U4 界面契约；仅使用 Node 内置模块，不启动真实浏览器。"""

from pathlib import Path
from email.message import Message
import shutil
import subprocess
from types import SimpleNamespace

import pytest

from chemgridmap.app import ASSETS, LocalHandler


@pytest.mark.parametrize("scenario", ["dictionary", "static", "help", "language", "dynamic", "server"])
def test_u4_web_contract(scenario):
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node is unavailable; browser acceptance remains a separate check.")
    script = Path(__file__).with_name("u4_web_contract.cjs")
    result = subprocess.run([node, str(script), str(ASSETS), scenario],
                            text=True, capture_output=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr


def test_i18n_static_route_returns_original_dictionary_bytes_without_socket():
    handler = object.__new__(LocalHandler)
    handler.server = SimpleNamespace(session=object(), origin="http://127.0.0.1:8765")
    handler.path = "/i18n.json"
    handler.headers = Message()
    handler.headers["Host"] = "127.0.0.1:8765"
    response = {}
    handler.send = lambda data, content_type="application/json", status=200, filename=None: response.update(
        data=data, content_type=content_type, status=status)
    handler.do_GET()
    assert response["status"] == 200
    assert response["content_type"] == "application/json"
    assert response["data"] == (ASSETS / "i18n.json").read_bytes()
