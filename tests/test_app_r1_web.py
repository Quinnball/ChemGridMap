"""R1 browser contracts, using Python and optional Node built-ins only."""

from pathlib import Path
import re
import shutil
import subprocess

import pytest

from chemgridmap.app import ASSETS


def test_unmeasured_outline_stays_screen_sized_below_all_highlights():
    stylesheet = (ASSETS / "app.css").read_text(encoding="utf-8")
    rules = list(re.finditer(r"([^{}]+)\{([^{}]+)\}", stylesheet))
    base = next(rule for rule in rules if ".unmeasured-outline>rect" in rule[1])
    assert "#map-stage circle.unmeasured-outline" in base[1]
    assert re.search(r"stroke:\s*#8A99A6\s*;", base[2])
    assert re.search(r"stroke-width:\s*calc\(1\.5\s*\*\s*var\(--zoom-inverse,\s*1\)\)\s*;", base[2])
    assert "!important" not in base[2]
    for selector in (".group-highlight>rect", ".neighbor-shared>rect", ".neighbor-grid_only>rect",
                     ".neighbor-representation_only>rect", "g.selected-cell>rect"):
        highlight = next(rule for rule in rules if selector in rule[1])
        assert highlight.start() > base.start()
        assert re.search(r"stroke:[^;]+!important", highlight[2])
        assert re.search(r"stroke-width:[^;]+!important", highlight[2])


@pytest.mark.parametrize("scenario", ["neighbors", "summary", "groups", "plural", "outlines",
                                      "audit", "input_reload", "support_legend", "generic_inspector"])
def test_r1_browser_behavior_with_node_builtins(scenario):
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node is unavailable; the Python browser contracts remain active.")
    script = Path(__file__).with_name("r1_web_contract.cjs")
    result = subprocess.run([node, str(script), str(ASSETS / "app.js"), scenario],
                            text=True, capture_output=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
