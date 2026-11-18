"""E28: run the predeclared aggregation for every set and collect descriptive quantities.

Reads the per-dataset outputs written by ``u10_neighbor_validation.py dataset`` and writes,
for each set, the outputs of ``aggregate`` plus ``descriptive.json`` (category fractions,
pooled agreement per consistent class, insertion times and placement distances). Usage::

    python paper/u10_e28_report.py --runs /path/e28 --runs-2010 /path/e28_2010 --output /path/e28_out
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import u10_neighbor_validation as v

SETS = {
    "primary": ("primary",),
    "all52": ("primary", "secondary"),
    "replication": ("replication",),
}


def _link_set(runs, parts, target):
    target.mkdir(parents=True, exist_ok=True)
    for part in parts:
        for folder in sorted((runs / part).glob("*/")):
            link = target / folder.name
            if not link.exists():
                link.symlink_to(folder.resolve())
    return target


def insertion_stats(root):
    seconds, sizes, distances, statuses = [], [], [], {}
    for manifest in sorted(root.glob("*/later_*/new_molecules/*/new_molecules_manifest.json")):
        info = json.loads(manifest.read_text(encoding="utf-8"))
        table = pd.read_csv(manifest.parent / "new_molecules.csv", usecols=lambda c: c in {"placement_distance", "query_id"})
        seconds.append(float(info.get("timings_seconds", {}).get("through_svg_and_tables", np.nan)))
        sizes.append(int(len(table)))
        if "placement_distance" in table:
            distances.extend(pd.to_numeric(table.placement_distance, errors="coerce").dropna().tolist())
    for path in sorted(root.glob("*/dataset.json")):
        info = json.loads(path.read_text(encoding="utf-8"))
        for key, value in (info.get("later_status_counts") or {}).items():
            statuses[key] = statuses.get(key, 0) + int(value)
    frame = pd.DataFrame({"seconds": seconds, "placed": sizes})
    full = frame[frame.placed >= 900]
    return {
        "batches": int(len(frame)),
        "batch_seconds_median": float(frame.seconds.median()) if len(frame) else None,
        "batches_900_plus": int(len(full)),
        "batch_900_plus_seconds_median": float(full.seconds.median()) if len(full) else None,
        "batch_900_plus_seconds_range": [float(full.seconds.min()), float(full.seconds.max())] if len(full) else None,
        "seconds_per_1000_compounds_median": float((frame.seconds / frame.placed * 1000).median()) if len(frame) else None,
        "placement_distance_median": float(np.median(distances)) if distances else None,
        "placement_distance_q25_q75": [float(np.quantile(distances, 0.25)), float(np.quantile(distances, 0.75))] if distances else None,
        "placement_distance_n": int(len(distances)),
        "later_status_counts": statuses,
    }


def descriptive(out):
    metrics = pd.read_csv(out / "metrics_by_dataset.csv")
    pooled = pd.read_csv(out / "pooled_by_category.csv")
    inclusion = pd.read_csv(out / "dataset_inclusion.csv")
    result = {"datasets": int(len(inclusion)), "included": int(inclusion.included.fillna(False).astype(bool).sum())}
    if "reason" not in inclusion:
        inclusion["reason"] = ""
    result["excluded"] = inclusion.loc[~inclusion.included.fillna(False).astype(bool), ["dataset", "reason"]].to_dict(orient="records")
    result["undated_raw_records"] = int(inclusion.undated_raw_records.sum())
    result["identity_conflicts_excluded"] = inclusion.loc[
        (inclusion.get("identity_conflicts_full") == "exclude") | (inclusion.get("identity_conflicts_reference") == "exclude"), "dataset"].tolist()
    by_endpoint = {}
    categories = ("none", "insufficient", "mixed", *v.CLASSES)
    for endpoint, part in metrics.groupby("endpoint"):
        total = part.n.sum()
        entry = {"datasets": int(len(part)), "later_compounds": int(total)}
        for category in categories:
            n = int(part[f"n_{category}"].sum())
            entry[f"{category}_n"] = n
            entry[f"{category}_pooled_fraction"] = n / total if total else None
            entry[f"{category}_median_fraction"] = float((part[f"n_{category}"] / part.n).median())
        entry["consistent_n"] = int(part.n_consistent.sum())
        entry["consistent_pooled_fraction"] = entry["consistent_n"] / total if total else None
        entry["consistent_median_fraction"] = float((part.n_consistent / part.n).median())
        entry["datasets_without_consistent"] = int((part.n_consistent == 0).sum())
        cons = part[part.n_consistent >= v.MIN_CONSISTENT]
        entry["n1_datasets"] = int(len(cons))
        entry["n1_agreement_median"] = float(cons.consistent_agreement.median()) if len(cons) else None
        entry["n1_chance_median"] = float(cons.consistent_chance.median()) if len(cons) else None
        pooled_part = pooled[(pooled.endpoint == endpoint) & pooled.neighbor_evidence.isin(v.CLASSES)]
        agree = int(pooled_part.agree.sum())
        n = int(pooled_part.n.sum())
        entry["consistent_pooled_agreement"] = agree / n if n else None
        entry["consistent_pooled_agreement_ci"] = list(v.wilson(agree, n)) if n else None
        entry["per_class"] = pooled_part[["neighbor_evidence", "n", "agree", "agreement", "ci_low", "ci_high", "datasets"]].to_dict(orient="records")
        measured = pooled[pooled.endpoint == endpoint][[f"measured_{c}" for c in v.CLASSES]].sum()
        entry["measured_distribution"] = {c: int(measured[f"measured_{c}"]) for c in v.CLASSES}
        for name in ("consistent", "unsettled", "none"):
            column = f"nearest_{name}_agreement"
            count = f"nearest_{name}_n"
            pooled_rate = (part[column] * part[count]).sum() / part[count].sum() if part[count].sum() else None
            entry[f"nearest_{name}_pooled"] = float(pooled_rate) if pooled_rate is not None else None
            entry[f"nearest_{name}_n"] = int(part[count].sum())
        by_endpoint[endpoint] = entry
    result["by_endpoint"] = by_endpoint
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", type=Path, required=True)
    parser.add_argument("--runs-2010", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    report = {}
    sets = dict(SETS)
    roots = {name: (args.runs, parts) for name, parts in sets.items()}
    if args.runs_2010 and (args.runs_2010 / "primary").exists():
        roots["primary_cutoff2010"] = (args.runs_2010, ("primary",))
    for name, (runs, parts) in roots.items():
        if not all((runs / part).exists() for part in parts):
            continue
        root = _link_set(runs, parts, args.output / "_links" / name)
        out = args.output / name
        summary = v.aggregate(root, out)
        entry = {"tests": summary["tests"], **descriptive(out), "insertion": insertion_stats(root)}
        (out / "descriptive.json").write_text(json.dumps(entry, indent=2, default=float), encoding="utf-8")
        report[name] = entry
    (args.output / "report.json").write_text(json.dumps(report, indent=2, default=float), encoding="utf-8")
    print(json.dumps({k: {"tests": r["tests"], "included": r["included"]} for k, r in report.items()}, indent=2, default=float))


if __name__ == "__main__":
    main()
