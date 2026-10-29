#!/usr/bin/env python3
"""Parallel execution of the pre-declared selection and download.

Same rules and order as download_chembl37_sample.py (Amendment 2 in
SELECTION_RULES.md). Count queries run concurrently in batches, but the
selection takes the first eligible targets in the sealed shuffled order, so
the chosen set is identical to the sequential procedure. Counts already
logged by the sequential run are reused.
"""
from __future__ import annotations

import csv
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

import download_chembl37_sample as base

THREADS = 6
BATCH = 60


def select_parallel(all_targets, standard_type):
    log_path = base.OUT / f"selection_log_{standard_type}.csv"
    known = {}
    if log_path.exists():
        known = {r["target_chembl_id"]: int(r["count"]) for r in csv.DictReader(log_path.open())}
    plan = base.PLAN[standard_type]
    pool = [t for t in all_targets if t not in base.EXCLUDED]
    order = [pool[i] for i in np.random.default_rng(plan["seed"]).permutation(len(pool))]
    rows, selected, position = [], [], 0
    with ThreadPoolExecutor(THREADS) as executor:
        while len(selected) < plan["n"] and position < len(order):
            batch = order[position:position + BATCH]
            missing = [t for t in batch if t not in known]
            for target, n in zip(missing, executor.map(lambda t: base.count(t, standard_type), missing)):
                known[target] = n
            for offset, target in enumerate(batch):
                if len(selected) >= plan["n"]:
                    break
                n = known[target]
                eligible = n >= base.MIN_RECORDS
                rows.append({"rank": position + offset + 1, "target_chembl_id": target, "count": n, "eligible": eligible})
                if eligible:
                    selected.append(target)
            position += len(batch)
            with log_path.open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=["rank", "target_chembl_id", "count", "eligible"])
                writer.writeheader()
                writer.writerows(rows)
            print(f"{standard_type}: walked {len(rows)}, selected {len(selected)}", flush=True)
    return selected, {r["target_chembl_id"]: r["count"] for r in rows}


def main():
    status = base.fetch("status", {})
    base.write_json(base.OUT / "api_status_parallel_run.json", {**status, "retrieved_at": time.strftime("%Y-%m-%dT%H:%M:%S%z")})
    if status.get("chembl_db_version") != "ChEMBL_37":
        raise SystemExit("Release is not ChEMBL_37; stopping.")
    all_targets = base.candidates()
    plan_rows = []
    for standard_type in base.PLAN:
        chosen, observed = select_parallel(all_targets, standard_type)
        plan_rows += [(t, standard_type, observed[t]) for t in chosen]
    for target in sorted(base.EXCLUDED):
        plan_rows.append((target, "Ki", base.count(target, "Ki")))
    with (base.OUT / "selected_targets.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["target_chembl_id", "standard_type", "pchembl_records", "role"])
        for target, standard_type, n in plan_rows:
            writer.writerow([target, standard_type, n, "original_target_second_endpoint" if target in base.EXCLUDED else "sampled"])
    manifest = {"rules": "SELECTION_RULES.md", "rules_sha256": base.sha(base.OUT / "SELECTION_RULES.md"),
                "candidates": len(all_targets), "execution": "parallel (Amendment 2)", "files": []}
    with ThreadPoolExecutor(4) as executor:
        results = list(executor.map(lambda row: {"target": row[0], "standard_type": row[1], **base.download(row[0], row[1], row[2])}, plan_rows))
    manifest["files"] = results
    manifest["finished_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    base.write_json(base.OUT / "download_manifest.json", manifest)
    bad = [f for f in results if not str(f["status"]).startswith(("complete", "verified"))]
    print(f"done: {len(results)} files, {sum(f['records'] for f in results)} records, problems {len(bad)}", flush=True)


if __name__ == "__main__":
    main()
