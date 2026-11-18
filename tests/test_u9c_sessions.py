"""Saved sessions restore complete state without computation or export mutation."""
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import pickle
import shutil
import subprocess
import sys
from urllib.parse import urlencode

import pandas as pd
import pytest

from chemgridmap.app import AppSession, SessionSnapshotError, SessionUnpickler
from test_app_u8_session import api
from test_app_u7_source_lookup import MEASURED, LOOKUP, OPTIONS as RECORD_OPTIONS

OPTIONS = {"input_format": "molecule", "value_col": "value", "projection": "pca",
           "grouping": "butina", "group_similarity": .1, "grid_refinement": "none",
           "render_detail": "overview", "output_formats": ["svg"]}
REFERENCE = b"molecule_id,smiles,value\nR1,CCO,8\nR2,CCN,8\nR3,CCC,8\nR4,c1ccccc1,8\n"


def hashes(root):
    return {path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in root.rglob("*") if path.is_file()}


def build(root, records=False):
    session = AppSession(root)
    session.load(MEASURED if records else REFERENCE, "reference.csv")
    if records:
        session.load_source(LOOKUP, "sources.csv")
    session.start({**(RECORD_OPTIONS if records else OPTIONS), "grouping": "butina"}, background=False)
    assert session.job["state"] == "done", session.job
    return session


def exercise(session):
    session.set_core_rule("coverage_brics")
    session.group_cores("G001")
    session.load_query(b"id,smiles\nN1,CCCl\nN2,CCBr\n", "batch.csv")
    session.insert_new({"id_col": "id"})
    session.add_new({"smiles": "CCF", "query_id": "MANUAL"})
    session.record_result({"query_id": "N1", "measurement": {"mode": "class", "measured_class": "inactive"}})
    session.record_result({"query_id": "MANUAL", "measurement": {"value": 100, "unit": "nM", "batch": "test"}})


def collect(session, reference_id="R1"):
    paths = ["/api/session", "/api/points", "/api/groups", "/api/group?group=G001",
             "/api/group-cores?group=G001", "/api/inspect?id=" + reference_id, "/api/inspect?id=N1",
             "/api/source-members?source=test&scope=new",
             "/api/new-molecules"]
    source_ids = list(session._source_index())
    if source_ids:
        paths.append("/api/source-members?" + urlencode({"source": source_ids[0]}))
    output = {}
    for path in paths:
        response = api(session, path)
        assert response["status"] == 200, response
        output[path] = response["data"]
    output["feedback"] = api(session, "/api/new-molecules/feedback")["data"].decode()
    return json.loads(json.dumps(output, allow_nan=False))


@pytest.mark.parametrize("records", [False, True])
def test_complete_roundtrip_in_fresh_process_and_read_only_open(tmp_path, records):
    session = build(tmp_path, records)
    exercise(session)
    reference_id = "C1" if records else "R1"
    expected = collect(session, reference_id)
    root = Path(session.job["output_dir"])
    before = hashes(root)
    target = tmp_path / "restored.json"
    script = """
import json, sys
from pathlib import Path
from chemgridmap.app import AppSession
from test_u9c_sessions import collect, hashes
session=AppSession(sys.argv[1]); root=Path(sys.argv[1])/sys.argv[2]
before=hashes(root); session.open_run(sys.argv[2]); assert hashes(root)==before
Path(sys.argv[3]).write_text(json.dumps(collect(session,sys.argv[4]),allow_nan=False))
assert session.set_core_rule('coverage_brics')['cached']
assert session.group_cores('G001')['cores']
"""
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([str(Path("src").resolve()), str(Path("tests").resolve())]))
    result = subprocess.run([sys.executable, "-c", script, str(tmp_path), root.name, str(target), reference_id],
                            env=env, capture_output=True, text=True, timeout=90)
    assert result.returncode == 0, result.stdout + result.stderr
    assert json.loads(target.read_text(encoding="utf-8")) == expected
    # Only subsequent cached-rule selection updates session state. Existing exports are unchanged.
    after = hashes(root)
    assert {key: after[key] for key in before if "session/" not in key and not key.startswith("audit_")} == {
        key: value for key, value in before.items() if "session/" not in key and not key.startswith("audit_")}


def test_moved_directory_and_future_mutations_stay_inside_new_root(tmp_path):
    session = build(tmp_path / "original")
    exercise(session)
    source = Path(session.job["output_dir"])
    destination = tmp_path / "moved" / source.name
    destination.parent.mkdir()
    shutil.move(source, destination)
    before = hashes(destination)
    restored = AppSession(destination.parent)
    restored.open_run(destination.name)
    assert hashes(destination) == before
    assert restored.job["output_dir"] == str(destination)
    assert all(path.is_relative_to(destination) for path in restored.files.values())
    assert all(Path(path).is_relative_to(destination) for path in restored.result.output_files.values())
    assert all(Path(path).is_relative_to(destination) for path in restored.new_result.output_files.values())
    restored.clear_result({"query_id": "N1"})
    assert restored.measured_result_info("N1")["measurement_status"] == "not_measured"
    assert not source.exists()
    restored.clear_query()
    again = AppSession(destination.parent)
    again.open_run(destination.name)
    assert again.new_result is None and again.manual_entries == [] and again.measured_results == {}
    assert again.previous_runs()["runs"][0]["new_molecules"] == 0


@pytest.mark.parametrize("kind", ["output", "snapshot", "version", "source", "deleted"])
def test_tampering_is_rejected_with_translatable_reason_and_no_partial_restore(tmp_path, kind):
    session = build(tmp_path)
    root = Path(session.job["output_dir"])
    path = root / "session" / "session_manifest.json"
    if kind in {"version", "source"}:
        manifest = json.loads(path.read_text(encoding="utf-8"))
        if kind == "version":
            manifest["versions"]["python"] = "0.0.0"
        else:
            manifest["source_tree_digest"]["source_sha256"] = "changed"
        path.write_text(json.dumps(manifest))
    elif kind == "snapshot":
        (root / "session/state.pkl.gz").write_bytes(b"changed")
    elif kind == "deleted":
        session.files["map_grid_coordinates.csv"].unlink()
    else:
        with session.files["map_grid_coordinates.csv"].open("ab") as stream:
            stream.write(b"\n")
    other = AppSession(tmp_path)
    other.load(REFERENCE, "live.csv")
    before = other.input_info(), other.snapshot()
    response = api(other, "/api/open-run", json.dumps({"run_name": root.name}).encode(), method="POST")
    assert response["status"] == 400
    assert response["data"]["code"] == ("session_version" if kind in {"source", "version"} else "session_modified")
    assert "regenerate" in response["data"]["error"]
    if kind not in {"source", "version"}:
        assert response["data"]["files"]
    assert (other.input_info(), other.snapshot()) == before


def test_restricted_unpickler_rejects_executable_globals_even_with_matching_hash(tmp_path):
    marker = tmp_path / "must-not-exist"
    class Payload:
        def __reduce__(self):
            return os.system, ("touch " + str(marker),)
    raw = pickle.dumps(Payload())
    with pytest.raises(SessionSnapshotError, match="session_unsafe"):
        SessionUnpickler(io.BytesIO(raw)).load()
    for dangerous in (pd.read_pickle, __import__("numpy").load, eval):
        with pytest.raises(SessionSnapshotError, match="session_unsafe"):
            SessionUnpickler(io.BytesIO(pickle.dumps(dangerous))).load()
    session = build(tmp_path / "runs")
    root = Path(session.job["output_dir"])
    path = root / "session/reference.pkl.gz"
    path.write_bytes(gzip.compress(raw))
    manifest_path = root / "session/session_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["files"]["session/reference.pkl.gz"] = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(SessionSnapshotError, match="session_unsafe"):
        AppSession(session.output_root).open_run(root.name)
    assert not marker.exists()


@pytest.mark.parametrize("run_name", ["../outside", "/tmp", "..", ".", "run/../other", "run\\other", "missing"])
def test_path_traversal_and_unknown_run_are_rejected(tmp_path, run_name):
    with pytest.raises(SessionSnapshotError, match="run_not_listed"):
        AppSession(tmp_path).open_run(run_name)


def test_building_cannot_open_and_new_endpoints_require_authentication(tmp_path):
    session = build(tmp_path)
    name = Path(session.job["output_dir"]).name
    for method, path in [("GET", "/api/previous-runs"), ("GET", "/api/core-progress"), ("POST", "/api/open-run")]:
        response = api(session, path, b"{}", method=method, token=False)
        assert response["status"] == 403
    session.update(state="running")
    response = api(session, "/api/open-run", json.dumps({"run_name": name}).encode(), method="POST")
    assert response["status"] == 400 and response["data"]["code"] == "map_running"


def test_list_is_recent_twenty_with_unavailable_counts_and_open_envelope(tmp_path):
    session = build(tmp_path)
    root = Path(session.job["output_dir"])
    for i in range(22):
        shutil.copytree(root, tmp_path / f"copied-{i:02}")
    (tmp_path / "old-run").mkdir()
    broken = tmp_path / "copied-00/session/session_manifest.json"
    manifest = json.loads(broken.read_text(encoding="utf-8")); manifest["versions"]["python"] = "old"
    broken.write_text(json.dumps(manifest))
    response = api(session, "/api/previous-runs")
    assert response["status"] == 200
    listing = response["data"]
    assert len(listing["runs"]) == 20
    assert listing["unavailable"] == {"missing_snapshot": 1, "incompatible": 1, "modified": 0}
    assert listing["runs"] == sorted(listing["runs"], key=lambda row: (row["created_at"], row["run_name"]), reverse=True)
    row = listing["runs"][0]
    assert row["filename"] == "reference.csv" and row["molecules"] == 4 and row["core_rule"] == "coverage"
    opened = api(session, "/api/open-run", json.dumps({"run_name": row["run_name"]}).encode(), method="POST")
    assert opened["status"] == 200
    assert opened["data"]["job"]["state"] == "done" and opened["data"]["input"]["filename"] == "reference.csv"
    assert opened["data"]["restored"]["total"] == 0


def test_snapshot_files_are_not_in_export_manifest_or_zip(tmp_path):
    import zipfile
    session = build(tmp_path)
    manifest = json.loads(Path(session.result.output_files["run_manifest"]).read_text(encoding="utf-8"))
    assert not any("session" in name for name in manifest["outputs"].values())
    with zipfile.ZipFile(session.files["chemgridmap_results.zip"]) as archive:
        assert not any(name.startswith("session/") for name in archive.namelist())


@pytest.mark.parametrize("kind", ["missing_metadata", "missing_active", "bad_job", "bad_result"])
def test_malformed_verified_state_does_not_partially_restore(tmp_path, kind):
    session = build(tmp_path)
    root = Path(session.job["output_dir"])
    manifest_path = root / "session/session_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if kind == "missing_metadata":
        del manifest["summary"]
    else:
        name = "reference.pkl.gz" if kind == "bad_result" else "state.pkl.gz"
        blob = root / "session" / name
        state = SessionUnpickler(io.BytesIO(gzip.decompress(blob.read_bytes()))).load()
        if kind == "missing_active":
            del state["active_set"]
        elif kind == "bad_job":
            state["job"] = None
        else:
            state["result"] = None
        blob.write_bytes(gzip.compress(pickle.dumps(state, protocol=5)))
        manifest["files"]["session/" + name] = hashlib.sha256(blob.read_bytes()).hexdigest()
    manifest_path.write_text(json.dumps(manifest))
    other = AppSession(tmp_path)
    other.load(REFERENCE, "existing.csv")
    before = other.input_info(), other.snapshot(), other.result
    with pytest.raises(SessionSnapshotError, match="session_invalid"):
        other.open_run(root.name)
    assert (other.input_info(), other.snapshot(), other.result) == before
    if kind == "missing_metadata":
        assert other.previous_runs()["unavailable"]["modified"] == 1


def test_open_uses_saved_objects_without_rebuilding_or_recurating(tmp_path, monkeypatch):
    from chemgridmap import app as app_module
    session = build(tmp_path)
    exercise(session)
    name = Path(session.job["output_dir"]).name
    def unexpected(*args, **kwargs):
        pytest.fail("Opening a saved result must not recalculate it")
    for method in ("build_grid_map", "curate_chembl_activity_data", "summarize_group_cores", "insert_new_molecules"):
        monkeypatch.setattr(app_module, method, unexpected)
    restored = AppSession(tmp_path)
    restored.open_run(name)
    assert restored.points() == session.points()
    assert restored.set_core_rule("coverage_brics")["cached"]
