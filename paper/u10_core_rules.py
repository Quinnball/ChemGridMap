"""U10 (post hoc): compare BRICS and RECAP candidate pools and the generic Murcko scaffold.

Uses the same 480 nonsingleton groups of the three detailed targets as the frozen
core comparison (``reliability_cores.run_cores``): Butina groups at similarity 0.6
on ``paper/output/revision_v6/<target>/map_molecules.csv`` sorted by identity key.

Writes, per group and rule, the selected core, its heavy atoms, coverage and
whether the 0.8 target is met, and per group the generic Murcko scaffold and the
number of members sharing it. Usage::

    python paper/u10_core_rules.py --input paper/output/revision_v6 --output paper/output/u10_cores
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, rdBase

from chemgridmap.comparison import structural_groups, summarize_group_core
from chemgridmap.cores import generic_scaffold_summary

TARGETS = {"chembl205": "CA II", "chembl204": "Thrombin", "chembl240": "hERG"}
EXPECTED_GROUPS = {"chembl205": 125, "chembl204": 152, "chembl240": 203}
RULES = ("most_frequent", "coverage", "coverage_brics", "coverage_recap")
SIZE_CLASSES = ((2, 2, "2"), (3, 4, "3-4"), (5, 9, "5-9"), (10, 10 ** 9, ">=10"))


def size_class(size):
    return next(label for low, high, label in SIZE_CLASSES if low <= size <= high)


def group_rows(frame, groups, target, label):
    """Rows for every rule in every nonsingleton group, plus generic-scaffold rows."""
    rows, generic_rows = [], []
    for group in groups:
        group_id = group["structure_group"]
        members = frame.loc[frame.structure_group.eq(group_id)]
        ids = members.molecule_identity_key.astype(str).tolist()
        mols = [Chem.MolFromSmiles(s) for s in members.canonical_smiles]
        common = {"target": target, "target_label": label, "group": group_id,
                  "size": len(mols), "size_class": size_class(len(mols))}
        for rule in RULES:
            started = time.perf_counter()
            summary = summarize_group_core(mols, ids, core_rule=rule, core_min_coverage=0.8)
            seconds = time.perf_counter() - started
            core = summary["core_smiles"]
            rows.append({**common, "rule": rule, "core_smiles": core,
                         "core_heavy_atoms": summary["core_heavy_atoms"],
                         "coverage_count": summary["core_coverage_count"],
                         "coverage_fraction": summary["core_coverage_fraction"],
                         "strict_majority": summary["core_coverage_count"] > len(mols) / 2,
                         "target_met": summary["core_coverage_count"] / len(mols) >= 0.8 - 1e-12,
                         "is_benzene": core == "c1ccccc1",
                         "extra_candidates": summary.get("brics_candidates", summary.get("recap_candidates", "")),
                         "recap_skipped_members": summary.get("recap_skipped_members", ""),
                         "seconds": seconds})
        generic = generic_scaffold_summary(mols)
        exact_count = int(group["core_exact_scaffold_count"])
        generic_rows.append({**common, **{k: v for k, v in generic.items() if k != "size"},
                             "most_frequent_scaffold_exact_count": exact_count,
                             "generic_groups_more_than_typed": generic["generic_scaffold_count"] > exact_count,
                             "generic_majority": generic["generic_scaffold_count"] > len(mols) / 2})
    return rows, generic_rows


SUBSET = 1200
SUBSET_SEED = 42


def run_dataset(raw_path, target, endpoint, output_root):
    """E28 part C: the same comparison on one sampled dataset (activity-independent subset)."""
    from chemgridmap import curate_chembl_activity_data
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(raw_path, low_memory=False, encoding="utf-8-sig")
    try:
        curated = curate_chembl_activity_data(raw, target_id=target, activity_type=endpoint)
    except ValueError:
        curated = curate_chembl_activity_data(raw, target_id=target, activity_type=endpoint,
                                              identity_conflicts="exclude")
    data = curated.molecules.sort_values("molecule_identity_key").reset_index(drop=True)
    entries = len(data)
    if entries > SUBSET:
        chosen = np.sort(np.random.default_rng(SUBSET_SEED).choice(entries, SUBSET, replace=False))
        data = data.iloc[chosen].reset_index(drop=True)
    frame, groups, _ = structural_groups(data, similarity=0.6, core_rule="most_frequent")
    groups = [g for g in groups if g["size"] >= 2]
    rows, generic_rows = group_rows(frame, groups, target, endpoint)
    pd.DataFrame(rows).to_csv(output_root / "core_rules_by_group.csv", index=False, lineterminator="\n")
    pd.DataFrame(generic_rows).to_csv(output_root / "generic_scaffold_by_group.csv", index=False, lineterminator="\n")
    info = {"target": target, "endpoint": endpoint, "entries": int(entries), "mapped_entries": int(len(data)),
            "groups": int(len(groups))}
    (output_root / "dataset.json").write_text(json.dumps(info, indent=2), encoding="utf-8")
    return info


def run(input_root, output_root):
    input_root, output_root = Path(input_root), Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    rows, generic_rows = [], []
    for target, label in TARGETS.items():
        data = pd.read_csv(input_root / target / "map_molecules.csv", low_memory=False)
        data = data.sort_values("molecule_identity_key").reset_index(drop=True)
        frame, groups, _ = structural_groups(data, similarity=0.6, core_rule="most_frequent")
        groups = [g for g in groups if g["size"] >= 2]
        if len(groups) != EXPECTED_GROUPS[target]:
            raise ValueError(f"Unexpected nonsingleton group count for {target}: {len(groups)}.")
        for group in groups:
            group_id = group["structure_group"]
            members = frame.loc[frame.structure_group.eq(group_id)]
            ids = members.molecule_identity_key.astype(str).tolist()
            mols = [Chem.MolFromSmiles(s) for s in members.canonical_smiles]
            common = {"target": target, "target_label": label, "group": group_id,
                      "size": len(mols), "size_class": size_class(len(mols))}
            for rule in RULES:
                started = time.perf_counter()
                summary = summarize_group_core(mols, ids, core_rule=rule, core_min_coverage=0.8)
                seconds = time.perf_counter() - started
                core = summary["core_smiles"]
                rows.append({**common, "rule": rule, "core_smiles": core,
                             "core_heavy_atoms": summary["core_heavy_atoms"],
                             "coverage_count": summary["core_coverage_count"],
                             "coverage_fraction": summary["core_coverage_fraction"],
                             "strict_majority": summary["core_coverage_count"] > len(mols) / 2,
                             "target_met": summary["core_coverage_count"] / len(mols) >= 0.8 - 1e-12,
                             "is_benzene": core == "c1ccccc1",
                             "extra_candidates": summary.get("brics_candidates", summary.get("recap_candidates", "")),
                             "recap_skipped_members": summary.get("recap_skipped_members", ""),
                             "seconds": seconds})
            generic = generic_scaffold_summary(mols)
            most_frequent = next(r for r in reversed(rows) if r["rule"] == "most_frequent")
            exact_count = int(group["core_exact_scaffold_count"])
            generic_rows.append({**common, **{k: v for k, v in generic.items() if k != "size"},
                                 "most_frequent_scaffold_exact_count": exact_count,
                                 "generic_groups_more_than_typed": generic["generic_scaffold_count"] > exact_count,
                                 "generic_majority": generic["generic_scaffold_count"] > len(mols) / 2,
                                 "most_frequent_core": most_frequent["core_smiles"]})
    table, generic_table = pd.DataFrame(rows), pd.DataFrame(generic_rows)
    table.to_csv(output_root / "core_rules_by_group.csv", index=False, lineterminator="\n")
    generic_table.to_csv(output_root / "generic_scaffold_by_group.csv", index=False, lineterminator="\n")

    wide = table.pivot_table(index=["target", "group"], columns="rule", values="core_smiles", aggfunc="first")
    heavy = table.pivot_table(index=["target", "group"], columns="rule", values="core_heavy_atoms", aggfunc="first")
    summary = {"groups": int(len(wide)), "by_rule": {}, "changes": {}, "generic": {}, "by_size_class": {}}
    for rule in RULES:
        part = table[table.rule.eq(rule)]
        summary["by_rule"][rule] = {
            "target_met": int(part.target_met.sum()), "strict_majority": int(part.strict_majority.sum()),
            "mean_heavy_atoms": round(float(part.core_heavy_atoms.mean()), 2),
            "benzene_cores": int(part.is_benzene.sum()), "total_seconds": round(float(part.seconds.sum()), 2),
            "median_ms": round(float(part.seconds.median() * 1000), 2)}
        summary["by_size_class"][rule] = {
            label: {"groups": int(len(sub)), "mean_heavy_atoms": round(float(sub.core_heavy_atoms.mean()), 2),
                    "target_met": int(sub.target_met.sum())}
            for label, sub in part.groupby("size_class", sort=False)}
    for rule, base in (("coverage_brics", "coverage"), ("coverage_recap", "coverage"),
                       ("coverage_recap", "coverage_brics")):
        changed = wide[rule] != wide[base]
        larger = heavy[rule] > heavy[base]
        summary["changes"][f"{rule}_vs_{base}"] = {
            "changed": int(changed.sum()), "larger": int(larger.sum()),
            "smaller": int((heavy[rule] < heavy[base]).sum())}
    summary["changes"]["brics_equals_recap"] = int((wide["coverage_brics"] == wide["coverage_recap"]).sum())
    summary["generic"] = {
        "generic_majority": int(generic_table.generic_majority.sum()),
        "generic_target_met": int((generic_table.generic_scaffold_fraction >= 0.8 - 1e-12).sum()),
        "generic_groups_more_than_typed": int(generic_table.generic_groups_more_than_typed.sum()),
        "typed_exact_majority": int((generic_table.most_frequent_scaffold_exact_count
                                     > generic_table["size"] / 2).sum())}
    summary["rdkit_version"] = rdBase.rdkitVersion
    (output_root / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, help="Frozen revision_v6 folder (three detailed targets).")
    parser.add_argument("--raw", type=Path, help="One sampled ChEMBL export (E28 part C).")
    parser.add_argument("--target")
    parser.add_argument("--endpoint", default="IC50")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.raw:
        print(json.dumps(run_dataset(args.raw, args.target, args.endpoint, args.output), indent=2))
    else:
        print(json.dumps(run(args.input, args.output), indent=2))


if __name__ == "__main__":
    main()
