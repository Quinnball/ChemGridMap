"""Compare core-selection rules on all 480 nonsingleton groups of the manuscript.

Rules: (1) current most-frequent Bemis-Murcko scaffold, (2) coverage-aware
largest scaffold-network candidate with coverage >= 0.8, (3) RDKit FMCS with
threshold 0.8 (complete rings only; comparison method).
"""
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

from common import RESULTS, TARGETS, load_layout
from chemgridmap.comparison import structural_groups
from chemgridmap.cores import mcs_core, select_core

RDLogger.DisableLog("rdApp.*")
OUT = RESULTS / "core_selection"
OUT.mkdir(parents=True, exist_ok=True)
# Historical snapshot: original absolute path replaced; see v2_source_data/ in this folder.
FIG = "<research workspace>/figures/source_data"

previous = OUT / "group_core_comparison.csv"
cached_mcs = (pd.read_csv(previous).set_index(["target", "group"]) if previous.exists() else None)
MCS_COLUMNS = ["mcs_core", "mcs_coverage", "mcs_heavy_atoms", "mcs_timed_out", "mcs_seconds"]
sensitivity = []
rows = []
for target, (label, _) in TARGETS.items():
    data = load_layout(target, 42)
    frame, summaries, _ = structural_groups(data, similarity=0.6)
    keyed = frame.set_index("comparison_id")
    for summary in summaries:
        if summary["size"] < 2:
            continue
        members = keyed[keyed.structure_group == summary["structure_group"]]
        mols = [Chem.MolFromSmiles(s) for s in members.canonical_smiles]
        median_atoms = float(np.median([m.GetNumHeavyAtoms() for m in mols]))
        current_atoms = Chem.MolFromSmiles(summary["core_smiles"]).GetNumHeavyAtoms() if summary["core_smiles"] else 0
        aware = select_core(mols, min_coverage=0.8)
        for threshold in (0.6, 0.7, 0.8, 0.9, 1.0):
            choice = select_core(mols, min_coverage=threshold)
            sensitivity.append({"target": label, "group": summary["structure_group"], "target_coverage": threshold,
                                "coverage": choice["coverage_fraction"], "heavy_atoms": choice["core_heavy_atoms"],
                                "target_met": choice["target_met"], "member_heavy_atoms": median_atoms})
        key = (label, summary["structure_group"])
        if cached_mcs is not None and key in cached_mcs.index:
            cached = cached_mcs.loc[key]
            mcs = {"core_smarts": cached.mcs_core, "coverage_fraction": cached.mcs_coverage,
                   "core_heavy_atoms": int(cached.mcs_heavy_atoms), "timed_out": bool(cached.mcs_timed_out),
                   "seconds": float(cached.mcs_seconds)}
        else:
            mcs = mcs_core(mols, min_coverage=0.8)
        rows.append({
            "target": label, "group": summary["structure_group"], "size": summary["size"],
            "median_member_heavy_atoms": median_atoms,
            "current_core": summary["core_smiles"], "current_coverage": summary["core_coverage_fraction"] if summary["core_smiles"] else np.nan,
            "current_heavy_atoms": current_atoms,
            "aware_core": aware["core_smiles"], "aware_coverage": aware["coverage_fraction"] if aware["core_smiles"] else np.nan,
            "aware_heavy_atoms": aware["core_heavy_atoms"], "aware_target_met": aware["target_met"],
            "aware_candidates": aware["candidates"], "aware_seconds": aware["seconds"],
            "mcs_core": mcs["core_smarts"], "mcs_coverage": mcs["coverage_fraction"], "mcs_heavy_atoms": mcs["core_heavy_atoms"],
            "mcs_timed_out": mcs["timed_out"], "mcs_seconds": mcs["seconds"],
        })
    print(f"{label}: {sum(r['target'] == label for r in rows)} nonsingleton groups", flush=True)

table = pd.DataFrame(rows)
sens = pd.DataFrame(sensitivity)
sens.to_csv(OUT / "coverage_target_sensitivity.csv", index=False)
print(sens.groupby(["target", "target_coverage"]).agg(groups=("group", "size"), met=("target_met", "sum"),
      median_heavy=("heavy_atoms", "median"), median_coverage=("coverage", "median")).to_string())
table.to_csv(OUT / "group_core_comparison.csv", index=False)

# Reproduction check against the manuscript's Figure 3A source data.
fig = pd.read_csv(f"{FIG}/Figure_3_core_coverage.csv")
print("Figure 3A columns:", list(fig.columns)[:8])

summary = []
for label, g in table.groupby("target", sort=False):
    for rule in ("current", "aware", "mcs"):
        cov = g[f"{rule}_coverage"]
        atoms = g[f"{rule}_heavy_atoms"]
        summary.append({"target": label, "rule": rule, "groups": len(g),
                        "no_core": int(cov.isna().sum() if rule != "mcs" else (atoms == 0).sum()),
                        "coverage_gt_0.5": int((cov > 0.5).sum()), "coverage_ge_0.8": int((cov >= 0.8).sum()),
                        "median_coverage": float(cov.median()), "median_core_heavy_atoms": float(atoms.median()),
                        "median_core_fraction_of_member": float((atoms / g.median_member_heavy_atoms).median()),
                        "seconds_total": float(g[f"{rule}_seconds"].sum()) if rule != "current" else np.nan})
summary = pd.DataFrame(summary)
summary.to_csv(OUT / "core_rule_summary.csv", index=False)
print(summary.round(3).to_string(index=False))
print(table[table.group == "G002"][["target", "current_core", "current_coverage", "aware_core", "aware_coverage", "mcs_core", "mcs_coverage"]].to_string(index=False))
