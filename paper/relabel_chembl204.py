"""T2c: relabel existing outputs without rerunning analyses; retain a replayable audit."""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import zipfile

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "paper/output/reliability_v1"
DIRECTORIES = {"activity", "baselines", "cases", "cores", "figure1_cli_run",
               "figure_data", "neighborhoods", "stability", "layouts", "scale"}
RULES = {"d3_ece": "thrombin_ece", "adoption_d3_ece": "adoption_thrombin_ece"}
LABEL_COLUMNS = {"target_label", "prism_column", "rule"}
LATER_ROOT_SECTIONS = {"portability_fix", "software_checks_050"}
SOURCE_RECORDS = ("software_checks_050", "portability_fix", "label_correction")
NOTICE = ("Label correction (2026-10-24): CHEMBL204 is human thrombin "
          "(prothrombin; CHEMBL204). Early project files incorrectly labeled it as the "
          "dopamine D3 receptor. The target label has been corrected to Thrombin; "
          "all numerical results are unchanged.")


def digest(data):
    return hashlib.sha256(data).hexdigest()


def file_hash(path):
    with Path(path).open("rb") as stream:
        digest = hashlib.sha256()
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
        return digest.hexdigest()


def latest_source_record(manifest):
    """Return the newest section that records the current source tree."""
    for name in SOURCE_RECORDS:
        section = manifest.get(name)
        if section and "source_sha256" in section.get("code", {}):
            return section
    raise AssertionError("No source hash record found")


def json_bytes(value):
    return (json.dumps(value, ensure_ascii=False, indent=2) + "\n").encode()


def csv_rows(data):
    return list(csv.reader(io.StringIO(data.decode("utf-8"), newline="")))


def csv_bytes(rows, newline="\n"):
    stream = io.StringIO(newline="")
    csv.writer(stream, lineterminator=newline).writerows(rows)
    return stream.getvalue().encode()


def relabel_cell(column, value):
    if column == "target_label" and value == "D3":
        return "Thrombin"
    if column == "prism_column" and (value == "D3" or value.startswith("D3 | ")):
        return "Thrombin" + value[2:]
    if column == "rule":
        return RULES.get(value, value)
    return value


def assert_csv_preserved(before, after):
    old, new = csv_rows(before), csv_rows(after)
    assert len(old) == len(new), "CSV row count changed"
    if not old:
        return
    assert new[0] == [RULES.get(name, name) for name in old[0]], "Unexpected CSV header edit"
    for old_row, new_row in zip(old[1:], new[1:]):
        assert len(old_row) == len(new_row) == len(old[0]), "CSV width changed"
        for name, a, b in zip(old[0], old_row, new_row):
            if name in LABEL_COLUMNS:
                assert b == relabel_cell(name, a), f"Unexpected label edit: {name}"
            else:
                assert a == b, f"Non-label cell changed: {name}"


def relabel_csv(data):
    rows = csv_rows(data)
    if not rows:
        return data, {"changed_columns": [], "changed_cells": 0, "changed_headers": 0}
    header = rows[0]
    new = [[RULES.get(name, name) for name in header]]
    changed, columns = 0, set()
    for row in rows[1:]:
        updated = [relabel_cell(name, value) for name, value in zip(header, row)]
        for name, a, b in zip(header, row, updated):
            if a != b:
                changed += 1
                columns.add(name)
        new.append(updated)
    headers = sum(a != b for a, b in zip(header, new[0]))
    columns.update(a for a, b in zip(header, new[0]) if a != b)
    # Unchanged files remain byte-identical. CSV parsing never converts numbers.
    result = csv_bytes(new, "\r\n" if b"\r\n" in data else "\n") if changed or headers else data
    assert_csv_preserved(data, result)
    return result, {"changed_columns": sorted(columns), "changed_cells": changed,
                    "changed_headers": headers, "rows": len(rows) - 1,
                    "non_label_cells_identical": True}


def relabel_json_value(value, path=(), edits=None):
    """Only known label fields and the specified rule identifiers may change."""
    edits = [] if edits is None else edits
    if isinstance(value, dict):
        result = {}
        for key, child in value.items():
            new_key = RULES.get(key, key)
            assert new_key not in result, "Renamed JSON key collides"
            if new_key != key:
                edits.append({"path": list(path + (key,)), "kind": "key", "before": key, "after": new_key})
            result[new_key] = relabel_json_value(child, path + (key,), edits)
        return result
    if isinstance(value, list):
        return [relabel_json_value(child, path + (i,), edits) for i, child in enumerate(value)]
    new = value
    if path and path[-1] == "target_label" and value == "D3":
        new = "Thrombin"
    elif path and path[-1] == "rule" and isinstance(value, str):
        new = RULES.get(value, value)
    if new != value:
        edits.append({"path": list(path), "kind": "value", "before": value, "after": new})
    return new


def assert_json_preserved(before, after, edits):
    # Undo only the recorded label edits, then compare the entire JSON tree.
    restored = copy.deepcopy(after)
    for edit in reversed(edits):
        if edit["kind"] == "key":
            assert edit["path"][-1] == edit["before"] and RULES.get(edit["before"]) == edit["after"]
        else:
            assert ((edit["path"][-1] == "target_label" and edit["before"] == "D3" and edit["after"] == "Thrombin")
                    or (edit["path"][-1] == "rule" and RULES.get(edit["before"]) == edit["after"]))
        node = restored
        for key in edit["path"][:-1]:
            node = node[RULES.get(key, key) if isinstance(key, str) else key]
        key = edit["path"][-1]
        if edit["kind"] == "key":
            node[key] = node.pop(edit["after"])
        else:
            node[key] = edit["before"]
    assert restored == before, "Non-label JSON field changed"


def relabel_markdown(data, relative):
    lines = data.decode().splitlines(keepends=True)
    result, edits = [], []
    if relative == "figure_data/README.md":
        result.append(NOTICE + "\n\n")
        edits.append({"path": relative, "line_before": 0, "before": "", "after": NOTICE,
                      "context_before": "", "context_after": lines[0].rstrip()})
    for i, line in enumerate(lines):
        if line.startswith("Target-label issue requiring coordinated correction:"):
            new = NOTICE + "\n"
        else:
            new = re.sub(r"\bD3\b", "Thrombin", line)
            for old, replacement in RULES.items():
                new = new.replace(old, replacement)
        result.append(new)
        if new != line:
            edits.append({"path": relative, "line_before": i + 1, "before": line.rstrip("\r\n"),
                          "after": new.rstrip("\r\n"),
                          "context_before": lines[i - 1].rstrip() if i else "",
                          "context_after": lines[i + 1].rstrip() if i + 1 < len(lines) else ""})
    return "".join(result).encode(), edits


def output_files(root):
    return sorted(p for p in root.rglob("*") if p.is_file() and p.name != ".DS_Store"
                  and "__pycache__" not in p.parts)


def allowed(path, root):
    parts = path.relative_to(root).parts
    return parts[0] in DIRECTORIES and "checkpoints" not in parts


def residual_rows(root):
    result = []
    for path in output_files(root):
        relative = path.relative_to(root).as_posix()
        if relative.startswith("label_correction/") or path.suffix not in {".csv", ".json", ".md", ".log", ".txt", ".py", ".xml", ".cjs", ".js"}:
            continue
        for i, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
            if not re.search(r"\bD3\b", line):
                continue
            if relative.split("/")[0] in {"t1_checks", "t1_acceptance_opus", "t1_rework"} or "checkpoints" in path.parts:
                reason = "冻结的历史验收记录或基线 checkpoint"
            elif relative == "activity/headline_comparison.csv":
                reason = "冻结参考的路径及历史字符串对照；不是当前 target_label"
            elif "BL21 D3 strain" in line:
                reason = "原始 assay_description 中的菌株名，禁止更改自由文本"
            elif relative == "figure_data/t2b_checks.json":
                reason = "T2b 历史检查中的组合列名；不符合 JSON 精确值 D3 的更改范围"
            elif path.suffix == ".log":
                reason = "已执行分析的历史日志，未重跑分析"
            elif "dopamine D3 receptor" in line:
                reason = "更正说明中引用历史误标名称"
            else:
                raise AssertionError(f"Unclassified D3: {relative}:{i}")
            result.append([relative, i, reason, line])
    return result


def verify(root, *, check_sources=True):
    root = Path(root)
    audit = json.loads((root / "label_correction/relabel_audit.json").read_text())
    assert file_hash(root / "label_correction/originals.zip") == audit["originals_sha256"]
    with zipfile.ZipFile(root / "label_correction/originals.zip") as originals:
        for entry in audit["files"]:
            path = root / entry["path"]
            assert file_hash(path) == entry["after_sha256"], f"Post-correction hash mismatch: {path}"
            if entry["before_sha256"] == entry["after_sha256"]:
                continue
            before = originals.read(entry["path"])
            assert digest(before) == entry["before_sha256"]
            after = path.read_bytes()
            if entry["kind"] == "csv":
                assert_csv_preserved(before, after)
            elif entry["kind"] == "json":
                assert_json_preserved(json.loads(before), json.loads(after), entry["edits"])
            elif entry["kind"] == "manifest":
                old, new = json.loads(before), json.loads(after)
                assert {k: v for k, v in new.items() if k not in {"outputs_sha256", "label_correction"}} == {
                    k: v for k, v in old.items() if k != "outputs_sha256"}
            elif entry["kind"] == "markdown":
                assert relabel_markdown(before, entry["path"])[0] == after
        current = json.loads((root / "manifest.json").read_text())
        old = json.loads(originals.read("manifest.json"))
        root_fields = {"outputs_sha256", "label_correction"} | LATER_ROOT_SECTIONS
        assert {k: v for k, v in current.items() if k not in root_fields} == {
            k: v for k, v in old.items() if k != "outputs_sha256"}
    paths = {p.relative_to(root).as_posix() for p in output_files(root) if p != root / "manifest.json"}
    assert paths == set(current["outputs_sha256"]), "Root output inventory mismatch"
    for relative, expected in current["outputs_sha256"].items():
        assert file_hash(root / relative) == expected, f"Output hash mismatch: {relative}"
    for entry in audit["files"]:
        if entry["kind"] != "manifest":
            continue
        manifest = json.loads((root / entry["path"]).read_text())
        for name, expected in manifest.get("outputs_sha256", {}).items():
            assert file_hash(name) == expected, f"Submanifest hash mismatch: {name}"
    if check_sources:
        for name, expected in latest_source_record(current)["code"]["source_sha256"].items():
            assert file_hash(name) == expected, f"Source hash mismatch: {name}"
    return {"passed": True, **audit["summary"], "root_outputs_verified": len(paths),
            "root_manifest_sha256": file_hash(root / "manifest.json")}


def run(root=OUTPUT):
    root = Path(root).resolve()
    if (root / "label_correction/relabel_audit.json").exists():
        return verify(root)  # Rerunning is a read-only verification; retain the original audit.
    correction = root / "label_correction"
    correction.mkdir(exist_ok=True)
    if (correction / "originals.zip").exists():
        raise RuntimeError("Incomplete prior run: restore originals.zip before retrying; original audit must not be overwritten.")
    files = [p for p in output_files(root) if allowed(p, root) or p == root / "RESULTS.md"]
    manifests = [p for p in files if p.name.endswith("manifest.json")]
    records, markdown_edits, backups = [], [], {"manifest.json": (root / "manifest.json").read_bytes()}
    updates = {}
    for path in files:
        if path in manifests:
            continue
        relative = path.relative_to(root).as_posix()
        before = path.read_bytes()
        after, detail, kind = before, {}, "unchanged"
        if path.suffix == ".csv":
            after, detail = relabel_csv(before)
            kind = "csv"
        elif path.suffix == ".json":
            old, edits = json.loads(before), []
            new = relabel_json_value(old, edits=edits)
            assert_json_preserved(old, new, edits)
            after = json_bytes(new) if edits else before
            detail = {"edits": edits, "changed_cells": len(edits), "non_label_fields_identical": True}
            kind = "json"
        elif relative in {"RESULTS.md", "figure_data/README.md"}:
            after, edits = relabel_markdown(before, relative)
            markdown_edits.extend(edits)
            detail = {"changed_lines": len(edits)}
            kind = "markdown"
        if after != before:
            backups[relative] = before
            updates[path] = after
        records.append({"path": relative, "kind": kind, "before_sha256": digest(before),
                        "after_sha256": digest(after), **detail})
    # Save recoverable originals before the first mutation, including all manifests.
    for path in manifests:
        backups[path.relative_to(root).as_posix()] = path.read_bytes()
        for name in json.loads(path.read_text()).get("outputs_sha256", {}):
            assert Path(name).is_file(), f"Missing manifest output: {name}"
    with zipfile.ZipFile(correction / "originals.zip", "x", zipfile.ZIP_DEFLATED) as archive:
        for name, data in backups.items():
            archive.writestr(name, data)
    for path, data in updates.items():
        path.write_bytes(data)
    fields = ["path", "line_before", "before", "after", "context_before", "context_after"]
    (correction / "markdown_replacements.csv").write_bytes(csv_bytes([fields] + [[row[k] for k in fields] for row in markdown_edits]))
    (correction / "residual_d3.csv").write_bytes(csv_bytes([["path", "line", "reason", "context"], *residual_rows(root)]))
    source_paths = sorted(p for directory in (ROOT / "src", ROOT / "paper", ROOT / "tests")
                          for p in directory.rglob("*") if p.is_file() and "output" not in p.parts
                          and "__pycache__" not in p.parts and p.suffix in {".py", ".js", ".cjs", ".html"})
    source_paths.append(ROOT / "validation_data/README.md")
    provenance = {"date": "2026-10-24", "scope": "CHEMBL204 labels only; no analyses rerun",
                  "target_labels": {"chembl205": "CA II", "chembl204": "Thrombin", "chembl240": "hERG"},
                  "audit": "label_correction/relabel_audit.json", "numerical_results_unchanged": True,
                  "code": {"source_sha256": {str(p): file_hash(p) for p in source_paths}}}
    pending = set(manifests)
    while pending:
        ready = [p for p in sorted(pending) if not {
            Path(name) for name in json.loads(p.read_text()).get("outputs_sha256", {})} & pending]
        assert ready, "Manifest output dependency cycle"
        for path in ready:
            before = path.read_bytes()
            manifest = json.loads(before)
            if "outputs_sha256" in manifest:
                manifest["outputs_sha256"] = {name: file_hash(name) for name in manifest["outputs_sha256"]}
            manifest["label_correction"] = {"date": provenance["date"], "scope": provenance["scope"],
                                           "previous_manifest_sha256": digest(before),
                                           "historical_fields_retained": True}
            after = json_bytes(manifest)
            relative = path.relative_to(root).as_posix()
            backups[relative] = before
            path.write_bytes(after)
            records.append({"path": relative, "kind": "manifest", "before_sha256": digest(before),
                            "after_sha256": digest(after), "changed_keys": ["outputs_sha256", "label_correction"],
                            "historical_fields_identical": True})
            pending.remove(path)
    summary = {"csv_files_checked": sum(r["kind"] == "csv" for r in records),
               "json_files_checked": sum(r["kind"] == "json" for r in records),
               "csv_cells_changed": sum(r.get("changed_cells", 0) for r in records if r["kind"] == "csv"),
               "csv_headers_changed": sum(r.get("changed_headers", 0) for r in records),
               "json_fields_changed": sum(r.get("changed_cells", 0) for r in records if r["kind"] == "json"),
               "markdown_lines_changed": len(markdown_edits), "manifests_refreshed": len(manifests) + 1}
    audit = {"passed": True, "summary": summary, "files": records,
             "root_manifest_before_sha256": digest(backups["manifest.json"]),
             "root_manifest_note": "The root is the outer hash envelope; its final SHA is printed and recorded in T2c_REPORT.md to avoid a self-hash cycle.",
             "originals_sha256": file_hash(correction / "originals.zip")}
    (correction / "relabel_audit.json").write_bytes(json_bytes(audit))
    root_manifest = json.loads(backups["manifest.json"])
    previous = root_manifest.get("t2b_rework", root_manifest["t2b"])
    provenance["inputs_sha256"] = {name: file_hash(name) for name in previous["inputs_sha256"]}
    provenance["reports_sha256"] = previous["reports_sha256"].copy()
    provenance["previous_collection_sha256"] = digest(backups["manifest.json"])
    root_manifest["label_correction"] = provenance
    root_manifest["outputs_sha256"] = {p.relative_to(root).as_posix(): file_hash(p) for p in output_files(root)
                                       if p != root / "manifest.json"}
    (root / "manifest.json").write_bytes(json_bytes(root_manifest))
    return verify(root)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true", help="Only replay and verify the saved correction audit.")
    parser.add_argument("--refresh-sources", action="store_true",
                        help="After output verification, record reviewed source-only fixes without repeating relabeling.")
    args = parser.parse_args()
    if args.refresh_sources:
        verify(args.root, check_sources=False)
        path = args.root / "manifest.json"
        manifest = json.loads(path.read_text())
        section = manifest["label_correction"]
        previous = section["code"]["source_sha256"]
        updates = {name: {"before": expected, "after": file_hash(name)} for name, expected in previous.items()
                   if file_hash(name) != expected}
        if updates:
            section.setdefault("source_hash_updates", []).append(updates)
            section["code"]["source_sha256"] = {name: file_hash(name) for name in previous}
            path.write_bytes(json_bytes(manifest))
    print(json.dumps(verify(args.root) if args.check else run(args.root), indent=2))


if __name__ == "__main__":
    main()
