#!/usr/bin/env python3
"""Select ChEMBL 37 targets by the pre-declared rules and download their activities.

Run by the coordinator (network required). Writes only under
strengthening/data/chembl37_api/. The rules are in SELECTION_RULES.md in that
directory. Re-running resumes: completed activity files are verified, not
downloaded again.
"""
from __future__ import annotations

import csv
import hashlib
import json
import time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import numpy as np

API = "https://www.ebi.ac.uk/chembl/api/data"
OUT = Path(__file__).resolve().parents[1] / "data/chembl37_api"
EXCLUDED = {"CHEMBL204", "CHEMBL205", "CHEMBL240"}
PLAN = {"IC50": {"seed": 20261028, "n": 50}, "Ki": {"seed": 20261029, "n": 30}}
MIN_RECORDS = 1000
FIELDS = ("activity_id", "molecule_chembl_id", "parent_molecule_chembl_id", "canonical_smiles",
          "standard_type", "standard_relation", "standard_value", "standard_units", "pchembl_value",
          "data_validity_comment", "potential_duplicate", "document_chembl_id", "document_year",
          "assay_chembl_id", "assay_type", "target_chembl_id", "src_id")


def fetch(endpoint, params, attempts=6):
    url = f"{API}/{endpoint}.json?{urlencode(params)}"
    for attempt in range(attempts):
        try:
            with urlopen(Request(url, headers={"User-Agent": "ChemGridMap/0.5 (research)"}), timeout=180) as r:
                return json.load(r)
        except Exception as error:  # Retries are logged; the last error is raised.
            print(f"retry {attempt + 1}/{attempts} {endpoint}: {error}", flush=True)
            time.sleep(min(60, 2 ** attempt * 3))
    raise RuntimeError(f"failed: {url}")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, payload):
    Path(path).write_text(json.dumps(payload, indent=1, ensure_ascii=False) + "\n")


def candidates():
    path = OUT / "human_single_protein_targets.csv"
    if path.exists():
        return [row["target_chembl_id"] for row in csv.DictReader(path.open())]
    rows, offset = [], 0
    while True:
        page = fetch("target", {"target_type": "SINGLE PROTEIN", "organism": "Homo sapiens",
                                "only": "target_chembl_id,pref_name", "limit": 1000, "offset": offset})
        rows += page["targets"]
        offset += 1000
        if not page["page_meta"]["next"]:
            break
    rows.sort(key=lambda r: r["target_chembl_id"])
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["target_chembl_id", "pref_name"])
        writer.writeheader()
        writer.writerows(rows)
    return [r["target_chembl_id"] for r in rows]


def count(target, standard_type):
    page = fetch("activity", {"target_chembl_id": target, "standard_type": standard_type,
                              "pchembl_value__isnull": "false", "limit": 1, "only": "activity_id"})
    return int(page["page_meta"]["total_count"])


def select(all_targets, standard_type):
    log_path = OUT / f"selection_log_{standard_type}.csv"
    done = {}
    if log_path.exists():
        done = {r["target_chembl_id"]: r for r in csv.DictReader(log_path.open())}
    plan = PLAN[standard_type]
    pool = [t for t in all_targets if t not in EXCLUDED]
    order = [pool[i] for i in np.random.default_rng(plan["seed"]).permutation(len(pool))]
    selected, rows = [], []
    for rank, target in enumerate(order, 1):
        if len(selected) >= plan["n"]:
            break
        n = int(done[target]["count"]) if target in done else count(target, standard_type)
        eligible = n >= MIN_RECORDS
        rows.append({"rank": rank, "target_chembl_id": target, "count": n, "eligible": eligible})
        if eligible:
            selected.append(target)
        if target not in done:
            with log_path.open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=["rank", "target_chembl_id", "count", "eligible"])
                writer.writeheader()
                writer.writerows(rows)
            time.sleep(0.2)
    print(f"{standard_type}: queried {len(rows)}, selected {len(selected)}", flush=True)
    return selected, {r["target_chembl_id"]: r["count"] for r in rows}


def download(target, standard_type, expected):
    folder = OUT / "activities"
    folder.mkdir(exist_ok=True)
    path = folder / f"{target}_{standard_type}.csv"
    if path.exists():
        with path.open() as f:
            n = sum(1 for _ in f) - 1
        if n == expected:
            return {"file": path.name, "records": n, "sha256": sha(path), "status": "verified_existing"}
    rows, offset = [], 0
    while offset < expected:
        page = fetch("activity", {"target_chembl_id": target, "standard_type": standard_type,
                                  "pchembl_value__isnull": "false", "only": ",".join(FIELDS),
                                  "limit": 1000, "offset": offset, "order_by": "activity_id"})
        rows += page["activities"]
        offset += 1000
        if not page["page_meta"]["next"]:
            break
    tmp = path.with_suffix(".part")
    with tmp.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(FIELDS), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    status = "complete" if len(rows) == expected else f"count_mismatch_expected_{expected}"
    tmp.replace(path)
    print(f"{target} {standard_type}: {len(rows)}/{expected} {status}", flush=True)
    return {"file": path.name, "records": len(rows), "expected": expected, "sha256": sha(path), "status": status}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    status = fetch("status", {})
    write_json(OUT / "api_status.json", {**status, "retrieved_at": time.strftime("%Y-%m-%dT%H:%M:%S%z")})
    if status.get("chembl_db_version") != "ChEMBL_37":
        raise SystemExit(f"Release is {status.get('chembl_db_version')}, not ChEMBL_37; stopping.")
    all_targets = candidates()
    plan_rows, counts = [], {}
    for standard_type in PLAN:
        chosen, observed = select(all_targets, standard_type)
        counts[standard_type] = observed
        plan_rows += [(t, standard_type, observed[t]) for t in chosen]
    for target in sorted(EXCLUDED):
        plan_rows.append((target, "Ki", count(target, "Ki")))
    with (OUT / "selected_targets.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["target_chembl_id", "standard_type", "pchembl_records", "role"])
        for target, standard_type, n in plan_rows:
            writer.writerow([target, standard_type, n, "original_target_second_endpoint" if target in EXCLUDED else "sampled"])
    manifest = {"rules": "SELECTION_RULES.md", "rules_sha256": sha(OUT / "SELECTION_RULES.md"),
                "candidates": len(all_targets), "files": []}
    for target, standard_type, n in plan_rows:
        manifest["files"].append({"target": target, "standard_type": standard_type,
                                  **download(target, standard_type, n)})
        write_json(OUT / "download_manifest.json", manifest)
    manifest["finished_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    write_json(OUT / "download_manifest.json", manifest)
    bad = [f for f in manifest["files"] if not str(f["status"]).startswith(("complete", "verified"))]
    print(f"done: {len(manifest['files'])} files, {sum(f['records'] for f in manifest['files'])} records, problems {len(bad)}")


if __name__ == "__main__":
    main()
