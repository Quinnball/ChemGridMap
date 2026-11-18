"""Group parallelism preserves chemistry, order, serial defaults and pool lifetime."""

import inspect
import json
import os
from pathlib import Path
import subprocess
import sys
import threading
from concurrent.futures import Future

import pandas as pd
import pytest
from rdkit import Chem

from chemgridmap import build_grid_map
from chemgridmap import comparison
from chemgridmap.cli import build_parser


def payloads():
    smiles = ["CCNC(=O)c1ccc(O)cc1", "CCCNC(=O)c1ccc(O)cc1",
              "CCCCNC(=O)c1ccc(O)cc1", "Oc1ccccc1"]
    return [(smiles * count, [f"G{group}-M{index}" for index in range(4 * count)])
            for group, count in enumerate((4, 7, 5, 6))]


def stable(result):
    return [(summary, {key: value for key, value in metrics.items() if key != "seconds"})
            for summary, metrics in result]


@pytest.mark.parametrize("rule", ["most_frequent", "coverage", "coverage_brics"])
def test_parallel_matches_serial_and_reports_parent_progress(rule, monkeypatch):
    # Exercise real workers on a small fixture; the measured production cutoff
    # is covered separately so this test need not repeat dozens of cheap groups.
    monkeypatch.setattr(comparison, "CORE_PARALLEL_MIN_GROUPS", 2)
    groups = payloads()
    expected = comparison.summarize_group_cores(groups, core_rule=rule)
    progress = []
    actual = comparison.summarize_group_cores(
        groups, core_rule=rule, jobs=2,
        progress=lambda done, total: progress.append((done, total, os.getpid(), threading.get_ident())))
    assert stable(actual) == stable(expected)
    assert progress == [(done, len(groups), os.getpid(), threading.get_ident())
                        for done in range(len(groups) + 1)]


def test_parallel_submits_plain_data_largest_first_and_restores_order(monkeypatch):
    monkeypatch.setattr(comparison, "CORE_PARALLEL_MIN_GROUPS", 2)
    groups, submitted = payloads(), []
    class Executor:
        def submit(self, function, *arguments):
            smiles, ids, rule, target, timeout = arguments
            assert all(isinstance(value, str) for value in smiles + ids)
            assert isinstance(rule, str) and isinstance(target, float) and isinstance(timeout, int)
            submitted.append(len(smiles))
            future = Future()
            future.set_result(function(*arguments))
            return future
    monkeypatch.setattr(comparison, "_get_core_executor", lambda jobs: Executor())
    expected = comparison.summarize_group_cores(groups)
    actual = comparison.summarize_group_cores(groups, jobs=8)
    assert submitted == sorted((len(smiles) for smiles, _ in groups), reverse=True)
    assert stable(actual) == stable(expected)


@pytest.mark.parametrize("rule,groups", [("fmcs", payloads()), ("coverage", payloads()[:1]),
                                       ("coverage", payloads() * 2),
                                       ("coverage", [(["c1ccccc1"], [str(i)]) for i in range(4)])])
def test_fmcs_and_small_workloads_never_start_workers(monkeypatch, rule, groups):
    if rule == "fmcs":
        monkeypatch.setattr(comparison, "CORE_PARALLEL_MIN_GROUPS", 2)
    monkeypatch.setattr(comparison, "_get_core_executor", lambda *_: pytest.fail("must stay serial"))
    expected = comparison.summarize_group_cores(groups, core_rule=rule)
    assert stable(comparison.summarize_group_cores(groups, core_rule=rule, jobs=8)) == stable(expected)


def test_library_cli_defaults_ignore_browser_environment(monkeypatch):
    monkeypatch.setenv("CHEMGRIDMAP_JOBS", "6")
    assert comparison.browser_core_jobs() == 6
    assert inspect.signature(comparison.summarize_group_cores).parameters["jobs"].default == 1
    assert inspect.signature(comparison.structural_groups).parameters["jobs"].default == 1
    assert inspect.signature(build_grid_map).parameters["jobs"].default == 1
    assert build_parser().parse_args([]).jobs == 1
    assert build_parser().parse_args(["--jobs", "4"]).jobs == 4
    with pytest.raises(SystemExit):
        build_parser().parse_args(["--jobs", "0"])
    with pytest.raises(ValueError, match="positive integer"):
        comparison.summarize_group_cores([], jobs=0)


@pytest.mark.parametrize("cpus,expected", [(None, 1), (1, 1), (2, 1), (6, 5), (12, 8), (32, 8)])
def test_browser_default_worker_count(monkeypatch, cpus, expected):
    monkeypatch.delenv("CHEMGRIDMAP_JOBS", raising=False)
    monkeypatch.setattr(comparison.os, "cpu_count", lambda: cpus)
    assert comparison.browser_core_jobs() == expected
    monkeypatch.setenv("CHEMGRIDMAP_JOBS", "invalid")
    assert comparison.browser_core_jobs() == expected


@pytest.mark.parametrize("rule", ["most_frequent", "coverage", "coverage_brics"])
def test_build_group_summary_and_member_order_unchanged(rule, monkeypatch):
    monkeypatch.setattr(comparison, "CORE_PARALLEL_MIN_GROUPS", 2)
    rows = [{"canonical_smiles": Chem.MolToSmiles(Chem.MolFromSmiles(smiles)),
             "molecule_id": key, "series": f"S{index}"}
            for index, (structures, keys) in enumerate(payloads())
            for smiles, key in zip(structures, keys)]
    data = pd.DataFrame(rows).sample(frac=1, random_state=8)
    options = dict(cluster_col="series", cluster_source="test input", core_rule=rule)
    expected_frame, expected_groups, expected_protocol = comparison.structural_groups(data, **options)
    actual_frame, actual_groups, actual_protocol = comparison.structural_groups(data, jobs=2, **options)
    pd.testing.assert_frame_equal(actual_frame, expected_frame)
    assert actual_groups == expected_groups and actual_protocol == expected_protocol
    for summary, members in zip(actual_groups, comparison.group_core_members(actual_frame, actual_groups)):
        expected = actual_frame[actual_frame.structure_group.eq(summary["structure_group"])].copy()
        expected = expected.sort_values(["canonical_smiles", "comparison_id"])
        center = expected.comparison_id.eq(summary["representative_id"])
        expected = pd.concat([expected[center], expected[~center]])
        pd.testing.assert_frame_equal(members, expected)


def test_pool_reuse_and_shutdown_leave_no_worker_processes(monkeypatch):
    monkeypatch.setattr(comparison, "CORE_PARALLEL_MIN_GROUPS", 2)
    comparison.summarize_group_cores(payloads(), jobs=2)
    first = comparison._get_core_executor(2)
    assert first is comparison._get_core_executor(2)
    processes = list(first._processes.values())
    assert len(processes) == 2 and all(process.is_alive() for process in processes)
    comparison.shutdown_core_workers()
    assert all(not process.is_alive() for process in processes)
    assert comparison._core_executor is None
    comparison.shutdown_core_workers()


def _process_exists(pid):
    """Return whether ``pid`` is a running process, without signalling it."""
    if sys.platform == "win32":
        # os.kill(pid, 0) would call TerminateProcess on Windows.
        import ctypes
        kernel32 = ctypes.windll.kernel32
        handle = kernel32.OpenProcess(0x1000, False, pid)  # PROCESS_QUERY_LIMITED_INFORMATION
        if not handle:
            return False
        try:
            code = ctypes.c_ulong()
            kernel32.GetExitCodeProcess(handle, ctypes.byref(code))
            return code.value == 259  # STILL_ACTIVE
        finally:
            kernel32.CloseHandle(handle)
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    return True


def _run_child(arguments, tmp_path):
    """Run a child interpreter with output in files, so that orphaned workers
    holding inherited pipes cannot block the test past its timeout."""
    out, err = tmp_path / "child.out", tmp_path / "child.err"
    with open(out, "w", encoding="utf-8") as stdout, open(err, "w", encoding="utf-8") as stderr:
        returncode = subprocess.run([sys.executable, *arguments], stdout=stdout, stderr=stderr,
                                    timeout=60).returncode
    return returncode, out.read_text(encoding="utf-8"), err.read_text(encoding="utf-8")


def test_normal_process_exit_reaps_workers(tmp_path):
    script = """
import json, os
from chemgridmap.comparison import _get_core_executor
executor = _get_core_executor(2)
assert executor.submit(os.getpid).result() != os.getpid()
print(json.dumps(list(executor._processes)), flush=True)
"""
    returncode, stdout, stderr = _run_child(["-c", script], tmp_path)
    assert returncode == 0, stderr
    worker_pids = json.loads(stdout)
    assert len(worker_pids) == 2
    assert not any(_process_exists(pid) for pid in worker_pids)


# Ctrl+C reaches the app as KeyboardInterrupt on every platform. SIGTERM is a
# POSIX stop request; on Windows os.kill(..., SIGTERM) is TerminateProcess,
# which no handler can intercept, so that case exists only on POSIX.
APP_STOPS = {"ctrl_c": "raise KeyboardInterrupt",
             "sigterm": "os.kill(os.getpid(), signal.SIGTERM)"}
if sys.platform == "win32":
    del APP_STOPS["sigterm"]


@pytest.mark.parametrize("stop", sorted(APP_STOPS))
def test_app_stop_reaps_workers_without_binding_a_port(tmp_path, stop):
    script = f"""
import json, os, signal, sys
from chemgridmap import app, comparison
worker_pids = []
class Server:
    origin = 'http://127.0.0.1'
    def __init__(self, *arguments):
        executor = comparison._get_core_executor(2)
        executor.submit(os.getpid).result()
        worker_pids.extend(executor._processes)
    def serve_forever(self):
        {APP_STOPS[stop]}
    def server_close(self):
        pass
app.LocalServer = Server
assert app.main(['--no-browser', '--output-dir', sys.argv[1]]) == 0
print(json.dumps(worker_pids), flush=True)
"""
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    returncode, stdout, stderr = _run_child(["-c", script, str(output_dir)], tmp_path)
    assert returncode == 0, stderr
    worker_pids = json.loads(stdout.splitlines()[-1])
    assert len(worker_pids) == 2
    assert not any(_process_exists(pid) for pid in worker_pids)


def test_empty_groups_report_complete_without_starting_workers():
    progress = []
    assert comparison.summarize_group_cores([], jobs=8, progress=lambda *value: progress.append(value)) == []
    assert progress == [(0, 0)]


@pytest.fixture(scope="module")
def caii_core_groups():
    from chemgridmap.chembl import curate_chembl_activity_data
    source = Path(comparison.__file__).parent / "data" / "chembl205_ic50.csv"
    curated = curate_chembl_activity_data(pd.read_csv(source), target_id="CHEMBL205")
    frame, summaries, _ = comparison.structural_groups(curated.molecules, core_rule="most_frequent")
    return [(members.canonical_smiles.tolist(), members.comparison_id.tolist())
            for members in comparison.group_core_members(frame, summaries)]


@pytest.mark.parametrize("rule", ["coverage", "coverage_brics"])
def test_packaged_caii_parallel_matches_serial(caii_core_groups, rule):
    expected = comparison.summarize_group_cores(caii_core_groups, core_rule=rule)
    actual = comparison.summarize_group_cores(caii_core_groups, core_rule=rule, jobs=4)
    assert stable(actual) == stable(expected)
