"""CLI 有无插入必须逐字节保留同一组参照输出。"""

import json
from pathlib import Path
import zipfile

import pandas as pd

from chemgridmap.cli import main


def manifest_without_allowed_variation(raw):
    value = json.loads(raw)
    for key in ("software_source_sha256", "git_commit", "git_dirty", "timings_seconds"):
        value.pop(key, None)
    return value


def compare_reference_directories(before, after):
    names = {path.name for path in Path(before).iterdir() if path.is_file()}
    assert names == {path.name for path in Path(after).iterdir() if path.is_file()}
    for name in names:
        a, b = (Path(before) / name).read_bytes(), (Path(after) / name).read_bytes()
        if name.endswith("_run_manifest.json"):
            assert manifest_without_allowed_variation(a) == manifest_without_allowed_variation(b)
        elif name.endswith(".zip"):
            with zipfile.ZipFile(Path(before) / name) as left, zipfile.ZipFile(Path(after) / name) as right:
                assert left.namelist() == right.namelist()
                for member in left.namelist():
                    x, y = left.read(member), right.read(member)
                    if member.endswith("_run_manifest.json"):
                        assert manifest_without_allowed_variation(x) == manifest_without_allowed_variation(y)
                    else:
                        assert x == y, member
        else:
            assert a == b, name


def test_cli_query_file_keeps_all_reference_outputs_byte_identical(tmp_path, monkeypatch):
    monkeypatch.setenv("SOURCE_DATE_EPOCH", "1793923200")
    source = tmp_path / "references.csv"
    pd.DataFrame({"molecule_id": ["R1", "R2", "R3", "R4"],
                  "smiles": ["CCO", "CCN", "CCC", "CCCC"],
                  "activity": [5.5, 6.5, 7.5, 8.0]}).to_csv(source, index=False)
    query = tmp_path / "new.tsv"
    query.write_bytes(b"id\tsmiles\r\nN1\tCCF\r\nN2\tCCCl\r\n")
    before, after = tmp_path / "before", tmp_path / "after"
    shared = [str(source), "--value-column", "activity", "--projection", "pca",
              "--grouping", "none", "--grid-refinement", "none", "--output-formats", "svg",
              "--render-detail", "overview"]
    assert main([*shared, "--output-dir", str(before)]) == 0
    assert main([*shared, "--output-dir", str(after), "--query-file", str(query),
                 "--query-id-column", "id"]) == 0
    compare_reference_directories(before, after)
    new_directory = after / "new_molecules" / "new"
    assert (new_directory / "new_molecules_input.csv").read_bytes() == query.read_bytes()
    placed = pd.read_csv(new_directory / "new_molecules.csv")
    assert placed.status.eq("mapped").all() and len(placed) == 2
