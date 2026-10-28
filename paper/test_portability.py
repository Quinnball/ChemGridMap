"""Exercise portable runner environments without regenerating frozen outputs."""
import copy
import importlib
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from check_t2b_outputs import latest_provenance


@pytest.mark.parametrize("newest", ["software_checks_050", "portability_fix", "label_correction",
                                   "t2b_rework", "t2b"])
def test_latest_provenance_prefers_newest_available_section(newest):
    names = ["software_checks_050", "portability_fix", "label_correction", "t2b_rework", "t2b"]
    manifest = {name: {"code": {"source_sha256": {"source.py": name}}}
                for name in reversed(names[names.index(newest):])}
    before = copy.deepcopy(manifest)
    assert latest_provenance(manifest) is manifest[newest]
    assert manifest == before


@pytest.fixture
def application_runner(tmp_path, monkeypatch):
    runner = importlib.import_module("check_reliability_merge")
    monkeypatch.setattr(runner, "CHECKS_ROOT", tmp_path / "checks")
    monkeypatch.setattr(runner, "INPUTS", {})
    output = tmp_path / "runner"
    output.mkdir()
    return runner, output


@pytest.mark.parametrize("custom,existing,expected", [
    ("custom_modules", "existing_modules", "custom_modules"),
    (None, "existing_modules", "existing_modules"),
    (None, None, None),
    ("", "existing_modules", ""),
])
def test_ui_node_path_override_and_inheritance(application_runner, monkeypatch,
                                               custom, existing, expected):
    runner, output = application_runner
    for name, value in (("CHEMGRIDMAP_NODE_PATH", custom), ("NODE_PATH", existing)):
        if value is None:
            monkeypatch.delenv(name, raising=False)
        else:
            monkeypatch.setenv(name, value)
    captured = {}

    def fake_run(command, **kwargs):
        captured.update(kwargs)
        captured["command"] = command
        return subprocess.CompletedProcess(command, 0, "ui completed\n", "")

    monkeypatch.setattr(runner.subprocess, "run", fake_run)
    result = runner.application_check(output, "ui")
    assert captured["env"].get("NODE_PATH") == expected
    assert ("NODE_PATH" in captured["env"]) == (expected is not None)
    assert os.environ.get("NODE_PATH") == existing
    assert captured["command"][0] == "node"
    assert Path(captured["command"][-1]) == runner.CHECKS_ROOT / "ui"
    assert result[0]["status"] == "passed"
    assert (output / "ui_runner.stdout.txt").read_text() == "ui completed\n"


@pytest.mark.parametrize("stderr,needs_hint", [
    ("Error: Cannot find module 'playwright'\n", True),
    ('Error: Cannot find module "playwright"\n', True),
    ("Error: UI assertion failed\n", False),
])
def test_ui_missing_playwright_reports_configuration_hint(application_runner, monkeypatch,
                                                         capsys, stderr, needs_hint):
    runner, output = application_runner
    monkeypatch.setattr(runner.subprocess, "run", lambda command, **kwargs:
                        subprocess.CompletedProcess(command, 1, "", stderr))
    result = runner.application_check(output, "ui")
    stored = (output / "ui_runner.stderr.txt").read_text()
    terminal = capsys.readouterr().err
    assert stored.startswith(stderr)
    assert ("CHEMGRIDMAP_NODE_PATH" in stored) is needs_hint
    assert ("CHEMGRIDMAP_NODE_PATH" in terminal) is needs_hint
    assert result[0]["status"] == "failed"
    assert result[0]["passed"] is False


@pytest.mark.parametrize("module,directory", [
    ("run_reliability_analysis", "chemgridmap-t2a-mpl"),
    ("reliability_contrasts", "chemgridmap-t2b-mpl"),
])
@pytest.mark.parametrize("explicit", [False, True])
def test_module_mplconfigdir_uses_system_temp_and_preserves_user_value(tmp_path, module,
                                                                    directory, explicit):
    root = Path(__file__).resolve().parents[1]
    temporary = tmp_path / "system_temp"
    temporary.mkdir()
    environment = dict(os.environ, TMPDIR=str(temporary), PYTHONDONTWRITEBYTECODE="1",
                       PYTHONPATH=os.pathsep.join((str(root / "src"), str(root / "paper"))))
    environment.pop("MPLCONFIGDIR", None)
    expected = temporary / directory
    if explicit:
        expected = tmp_path / "user_mpl_config"
        environment["MPLCONFIGDIR"] = str(expected)
    code = ("import importlib, json, os; "
            f"importlib.import_module({module!r}); "
            "print(json.dumps(os.environ['MPLCONFIGDIR']))")
    result = subprocess.run([sys.executable, "-B", "-c", code], cwd=tmp_path,
                            env=environment, capture_output=True, text=True, check=True)
    assert json.loads(result.stdout.strip().splitlines()[-1]) == str(expected)


@pytest.mark.parametrize("explicit", [False, True])
def test_figure_cli_mplconfigdir_default_and_inheritance(tmp_path, monkeypatch, explicit):
    runner = importlib.import_module("reliability_contrasts")
    temporary = tmp_path / "system_temp"
    temporary.mkdir()
    monkeypatch.setattr(runner.tempfile, "gettempdir", lambda: str(temporary))
    monkeypatch.setattr(runner, "_cli_signature", lambda command, environment: ("test", {}))
    monkeypatch.setattr(runner, "valid_cli_cache", lambda directory, signature: False)
    monkeypatch.delenv("MPLCONFIGDIR", raising=False)
    expected = temporary / "chemgridmap-t2b-mpl"
    if explicit:
        expected = tmp_path / "user_mpl_config"
        monkeypatch.setenv("MPLCONFIGDIR", str(expected))
    captured = {}

    def fake_run(command, **kwargs):
        captured.update(kwargs)
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(runner.subprocess, "run", fake_run)
    output = runner.prepare_figure1_cli(tmp_path)
    assert output == tmp_path / "figure1_cli_run"
    assert captured["env"]["MPLCONFIGDIR"] == str(expected)
    assert os.environ.get("MPLCONFIGDIR") == (str(expected) if explicit else None)
