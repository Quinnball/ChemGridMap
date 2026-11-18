#!/usr/bin/env python3
"""E11: continue the sealed sampling order to select new validation targets.

The plan is paper/strengthening/plans/E11_validation_plan_sealed.md. The
order, seeds, eligibility rule, fields and file checks reuse
download_chembl37_sample.py unchanged. Querying resumes after the last rank of
the first selection. Output is written under chembl37_api_validation/.
"""
from __future__ import annotations

import csv
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

import download_chembl37_sample as base

ORIGINAL = base.OUT
OUT = ORIGINAL.parent / "chembl37_api_validation"
START = {"IC50": 1015, "Ki": 2105}  # first rank not queried by the first selection
N_NEW = {"IC50": 30, "Ki": 30}
THREADS = 6
BATCH = 60


def sealed_order(all_targets, standard_type):
    pool = [t for t in all_targets if t not in base.EXCLUDED]
    seed = base.PLAN[standard_type]["seed"]
    return [pool[i] for i in np.random.default_rng(seed).permutation(len(pool))]


def check_against_first_log(order, standard_type):
    """The reconstructed order must reproduce every rank of the first selection log."""
    log = list(csv.DictReader((ORIGINAL / f"selection_log_{standard_type}.csv").open()))
    for row in log:
        if order[int(row["rank"]) - 1] != row["target_chembl_id"]:
            raise SystemExit(f"{standard_type}: order differs from the first log at rank {row['rank']}")
    if int(log[-1]["rank"]) != START[standard_type] - 1:
        raise SystemExit(f"{standard_type}: the first log ends at rank {log[-1]['rank']}")


def select_new(order, standard_type):
    log_path = OUT / f"selection_log_{standard_type}.csv"
    rows, selected, position = [], [], START[standard_type] - 1
    with ThreadPoolExecutor(THREADS) as executor:
        while len(selected) < N_NEW[standard_type] and position < len(order):
            batch = order[position:position + BATCH]
            counts = list(executor.map(lambda t: base.count(t, standard_type), batch))
            for offset, (target, n) in enumerate(zip(batch, counts)):
                if len(selected) >= N_NEW[standard_type]:
                    break
                eligible = n >= base.MIN_RECORDS
                rows.append({"rank": position + offset + 1, "target_chembl_id": target,
                             "count": n, "eligible": eligible})
                if eligible:
                    selected.append((target, n))
            position += len(batch)
            with log_path.open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=["rank", "target_chembl_id", "count", "eligible"])
                writer.writeheader()
                writer.writerows(rows)
            print(f"{standard_type}: walked {len(rows)}, selected {len(selected)}", flush=True)
    return selected


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    status = base.fetch("status", {})
    base.write_json(OUT / "api_status.json", {**status, "retrieved_at": time.strftime("%Y-%m-%dT%H:%M:%S%z")})
    if status.get("chembl_db_version") != "ChEMBL_37":
        raise SystemExit(f"Release is {status.get('chembl_db_version')}, not ChEMBL_37; stopping.")
    all_targets = base.candidates()  # the target list saved by the first run
    plan_rows = []
    for standard_type in base.PLAN:
        order = sealed_order(all_targets, standard_type)
        check_against_first_log(order, standard_type)
        plan_rows += [(t, standard_type, n) for t, n in select_new(order, standard_type)]
    with (OUT / "selected_targets.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["target_chembl_id", "standard_type", "pchembl_records", "role"])
        writer.writerows([(t, s, n, "validation") for t, s, n in plan_rows])
    base.OUT = OUT  # download() writes to base.OUT / "activities"
    manifest = {"plan": "paper/strengthening/plans/E11_validation_plan_sealed.md",
                "rules": "SELECTION_RULES.md", "rules_sha256": base.sha(ORIGINAL / "SELECTION_RULES.md"),
                "candidates": len(all_targets), "start_rank": START, "files": [],
                "started_at": time.strftime("%Y-%m-%dT%H:%M:%S%z")}
    with ThreadPoolExecutor(4) as executor:
        manifest["files"] = list(executor.map(
            lambda row: {"target": row[0], "standard_type": row[1], **base.download(row[0], row[1], row[2])},
            plan_rows))
    manifest["finished_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    base.write_json(OUT / "download_manifest.json", manifest)
    bad = [f for f in manifest["files"] if not str(f["status"]).startswith(("complete", "verified"))]
    print(f"done: {len(manifest['files'])} files, {sum(f['records'] for f in manifest['files'])} records, "
          f"problems {len(bad)}", flush=True)


if __name__ == "__main__":
    main()
