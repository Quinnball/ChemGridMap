"""核对 R1 写入范围；禁止进入后台下载子目录。"""
from pathlib import Path
import hashlib
import json
import os
import subprocess

ROOT = Path(__file__).resolve().parents[3]
REV = "paper/strengthening"
QA = ROOT / REV / "qa/E1A"
EXCLUDED = REV + "/data/chembl37_api"


def main():
    before = json.loads((QA / "workspace_R1_before.json").read_text())
    result_path = QA / "workspace_R1_changes.json"
    if not result_path.exists():
        result_path.write_text("{}\n")
    status = subprocess.check_output(["git", "-C", str(ROOT), "status",
        "--porcelain=v1", "--untracked-files=all"], text=True)
    (QA / "cgm_status_R1_after.txt").write_text(status)
    after = {}
    for base, dirs, files in os.walk(ROOT):
        dirs[:] = [d for d in dirs if (Path(base) / d).relative_to(ROOT).as_posix() != EXCLUDED]
        for name in files:
            path = Path(base) / name
            try:
                stat = path.stat()
            except FileNotFoundError:
                continue
            after[str(path.relative_to(ROOT))] = [stat.st_size, stat.st_mtime_ns]
    changes = {"created": sorted(set(after) - set(before)),
               "modified": sorted(p for p in set(before) & set(after) if before[p] != after[p]),
               "deleted": sorted(set(before) - set(after))}
    def allowed(path):
        return (path.startswith(REV + "/output/E1A/") or
                path.startswith(REV + "/qa/E1A/") or
                (path.startswith(REV + "/source/e1a_") and path.endswith(".py")) or
                path in [REV + "/source/test_e1a.py",
                         REV + "/output/E1A/E1A_REPORT.md"])
    unallowed = [p for names in changes.values() for p in names if not allowed(p)]
    # 不自动豁免范围外变更；保留实际证据供报告说明。
    metadata = ROOT / REV / "data/flagged_pairs_metadata.json"
    changes.update({"unexpected_changes": unallowed,
        "all_observed_workspace_changes_within_allowlist": not unallowed,
        "cgm_status_lines": len(status.splitlines()),
        "cgm_status_identical": status == (QA / "cgm_status_R1_before.txt").read_text(),
        "metadata_sha256_unchanged": hashlib.sha256(metadata.read_bytes()).hexdigest() ==
            (QA / "metadata_R1_before.sha256").read_text().strip(),
        "data_changes_excluding_download": [p for p in unallowed if p.startswith(REV + "/data/")],
        "excluded_without_traversal": EXCLUDED,
        "snapshot_method": "size and mtime_ns; independent guard and input SHA-256 checks also required"})
    result_path.write_text(json.dumps(changes, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(changes, indent=2, ensure_ascii=False))
    assert not unallowed, "存在范围外变化，需如实调查报告"
    assert changes["cgm_status_identical"]
    assert changes["metadata_sha256_unchanged"]
    print("SCOPE PASS: allowed paths only; workspace status identical; non-download data unchanged")


if __name__ == "__main__":
    main()
