"""The optional comparator must not become a runtime or test dependency."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "paper"))
from run_viewer_tasks import choose_cases, git_blob_sha


def test_case_selection_is_deterministic_and_not_activity_ranked():
    frame = pd.DataFrame([
        ["B", 5, 1, 9.0], ["A", 5, 1, 2.0], ["C", 3, 0, 7.0],
        ["D", 1, 0, 9.0], ["E", 1, 0, 1.0],
    ], columns=["molecule_identity_key", "n_records", "is_conflicted", "activity_pchembl"])
    expected = [("flagged_repeat", "A"), ("unflagged_repeat", "C"), ("singleton", "D")]
    for seed in (1, 42):
        assert [(group, row.molecule_identity_key) for group, row in choose_cases(frame.sample(frac=1, random_state=seed))] == expected


def test_git_blob_hash_includes_header(tmp_path):
    path = tmp_path / "empty"
    path.write_bytes(b"")
    assert git_blob_sha(path) == "e69de29bb2d1d6434b8b29ae775ad8c2e48c5391"
