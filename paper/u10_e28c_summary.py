"""E28 part C: pool the per-dataset core-rule comparisons of the sampled datasets.

Reads ``core_rules_by_group.csv``, ``generic_scaffold_by_group.csv`` and ``dataset.json``
written by ``u10_core_rules.py --raw`` for each dataset and writes the descriptive
quantities listed in section 6 of the E28 plan, pooled over all groups and by endpoint.
No hypothesis tests. Usage::

    python paper/u10_e28c_summary.py --root /path/e28c --output /path/e28_out/cores
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from u10_core_rules import RULES

SIZE_ORDER = ("2", "3-4", "5-9", ">=10")


def summarize(table, generic):
    keyed = table.set_index(["target", "endpoint", "group", "rule"])
    # Groups without a ring scaffold have no core under any rule; they stay in the counts (as empty cores).
    wide = keyed.core_smiles.unstack("rule").fillna("")
    heavy = keyed.core_heavy_atoms.unstack("rule").fillna(0)
    out = {"groups": int(len(wide)), "datasets": int(table[["target", "endpoint"]].drop_duplicates().shape[0]),
           "by_rule": {}, "by_size_class": {}, "changes": {}, "generic": {}}
    for rule in RULES:
        part = table[table.rule.eq(rule)]
        out["by_rule"][rule] = {
            "target_met": int(part.target_met.sum()), "strict_majority": int(part.strict_majority.sum()),
            "mean_heavy_atoms": round(float(part.core_heavy_atoms.mean()), 2),
            "benzene_cores": int(part.is_benzene.sum()), "total_seconds": round(float(part.seconds.sum()), 1),
            "median_ms": round(float(part.seconds.median() * 1000), 2)}
        out["by_size_class"][rule] = {
            label: {"groups": int(len(sub)), "mean_heavy_atoms": round(float(sub.core_heavy_atoms.mean()), 2),
                    "target_met": int(sub.target_met.sum()), "strict_majority": int(sub.strict_majority.sum())}
            for label, sub in part.groupby("size_class") }
    for rule, base in (("coverage_brics", "coverage"), ("coverage_recap", "coverage"), ("coverage_recap", "coverage_brics")):
        out["changes"][f"{rule}_vs_{base}"] = {
            "changed": int((wide[rule] != wide[base]).sum()), "larger": int((heavy[rule] > heavy[base]).sum()),
            "smaller": int((heavy[rule] < heavy[base]).sum())}
    out["groups_without_ring_core"] = int((wide["coverage"] == "").sum())
    out["changes"]["brics_equals_recap"] = int(((wide["coverage_brics"] == wide["coverage_recap"]) & (wide["coverage"] != "")).sum())
    out["generic"] = {
        "generic_majority": int(generic.generic_majority.sum()),
        "generic_groups_more_than_typed": int(generic.generic_groups_more_than_typed.sum()),
        "typed_exact_majority": int((generic.most_frequent_scaffold_exact_count > generic["size"] / 2).sum())}
    skipped = table[table.rule.eq("coverage_recap")].recap_skipped_members.fillna("").astype(str)
    out["recap_groups_with_skipped_members"] = int(pd.to_numeric(skipped, errors="coerce").fillna(0).gt(0).sum())
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    args.output.mkdir(parents=True, exist_ok=True)
    tables, generics, infos = [], [], []
    for folder in sorted(p for p in args.root.iterdir() if p.is_dir()):
        if not (folder / "dataset.json").exists():
            infos.append({"dataset": folder.name, "status": "failed"})
            continue
        info = json.loads((folder / "dataset.json").read_text(encoding="utf-8"))
        infos.append({"dataset": folder.name, "status": "ok", **info})
        if info.get("groups", 0) == 0:
            continue
        table = pd.read_csv(folder / "core_rules_by_group.csv")
        generic = pd.read_csv(folder / "generic_scaffold_by_group.csv")
        table["endpoint"] = info["endpoint"]
        generic["endpoint"] = info["endpoint"]
        tables.append(table)
        generics.append(generic)
    table, generic = pd.concat(tables, ignore_index=True), pd.concat(generics, ignore_index=True)
    table.to_csv(args.output / "core_rules_by_group.csv", index=False, lineterminator="\n")
    generic.to_csv(args.output / "generic_scaffold_by_group.csv", index=False, lineterminator="\n")
    pd.DataFrame(infos).to_csv(args.output / "datasets.csv", index=False, lineterminator="\n")
    summary = {"all": summarize(table, generic),
               "by_endpoint": {endpoint: summarize(table[table.endpoint.eq(endpoint)], generic[generic.endpoint.eq(endpoint)])
                               for endpoint in sorted(table.endpoint.unique())},
               "datasets_subsampled": int(sum(1 for i in infos if i.get("entries", 0) > i.get("mapped_entries", 0))),
               "datasets_failed": [i["dataset"] for i in infos if i["status"] == "failed"]}
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary["all"], indent=2))


if __name__ == "__main__":
    main()
