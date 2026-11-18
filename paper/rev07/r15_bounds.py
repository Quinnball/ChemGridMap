"""R15 (post hoc, exploratory): bound records and the agreement of inactive displays.

The evaluation downloads requested only records with a pChEMBL value, so they hold no bounds such as
"> 10,000 nM". The full IC50 exports of CA II (CHEMBL205) and thrombin (CHEMBL204) in
validation_data/ keep every relation. For these two development targets this script reports

1. the share of bound records among the records that place a compound in the inactive range
   (exact pChEMBL <= 6, or a bound that implies pActivity <= 6), and the same for the active range;
2. how many entries of each displayed class have bound records from documents that report no exact
   value for that entry;
3. held-out agreement with exact documents (each evidence unit of an entry with at least two units held
   out in turn; displayed class = class of the median of the remaining records, as model B0 of E1-A),
   and the same agreement recomputed with bound documents added as held-out documents: a bound that
   implies the displayed class counts as agreement, a bound that excludes it as disagreement, and a
   bound compatible with more than one class is left out.

Bounds pass every curation rule except the relation (and the pChEMBL value, which ChEMBL gives only for
exact relations): nM units, no data-validity comment other than "Manually validated", no
potential-duplicate flag and a valid structure. They are linked to entries by parent identifier.
Usage: PYTHONPATH=src python paper/rev07/r15_bounds.py <output folder>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from chemgridmap import curate_chembl_activity_data, document_values

ROOT = Path(__file__).resolve().parents[2]
TARGETS = {
    "CHEMBL205": ROOT / "validation_data/chembl37_chembl205/chembl37_chembl205_ic50_raw.csv",
    "CHEMBL204": ROOT / "validation_data/chembl37_chembl204/chembl37_chembl204_ic50_raw.csv",
}
LOWER, UPPER = 6.0, 7.0
BOUND_RELATIONS = {">", ">=", ">>", "<", "<=", "<<"}
DISQUALIFYING = ("data_validity_flag", "chembl_potential_duplicate", "invalid_or_missing_smiles",
                 "non_nm_standard_units", "wrong_target", "wrong_activity_type")


def klass(p):
    return np.where(p <= LOWER, "inactive", np.where(p >= UPPER, "active", "medium"))


def bound_outcome(displayed, relation, p_bound):
    """'agree', 'disagree' or None (the bound is compatible with more than one class)."""
    upper_bound = relation.startswith(">")  # IC50 > V  ->  pActivity < p_bound
    if upper_bound:
        if displayed == "inactive":
            return "agree" if p_bound <= LOWER else None
        if displayed == "medium":
            return "disagree" if p_bound <= LOWER else None
        return "disagree" if p_bound < UPPER or (p_bound == UPPER and relation != ">=") else None
    if displayed == "active":
        return "agree" if p_bound >= UPPER else None
    if displayed == "medium":
        return "disagree" if p_bound >= UPPER else None
    return "disagree" if p_bound > LOWER or (p_bound == LOWER and relation != "<=") else None


def bounds_table(excluded):
    rel = excluded.standard_relation.fillna("").astype(str).str.strip().str.strip("'\"")
    reasons = excluded.exclusion_reason.fillna("").astype(str)
    keep = rel.isin(BOUND_RELATIONS) & reasons.str.contains("non_exact_standard_relation")
    for reason in DISQUALIFYING:
        keep &= ~reasons.str.contains(reason)
    b = excluded.loc[keep].copy()
    b["relation"] = rel[keep]
    value = pd.to_numeric(b.standard_value, errors="coerce")
    b = b[np.isfinite(value) & (value > 0) & b.standard_units.astype(str).str.strip().eq("nM")].copy()
    b["p_bound"] = -np.log10(pd.to_numeric(b.standard_value) * 1e-9)
    b["entry"] = b.parent_molecule_chembl_id.fillna(b.molecule_chembl_id).astype(str)
    b["document"] = b.document_chembl_id.astype(str)
    return b


def exact_loo(records):
    """Hold out each evidence unit of entries with at least two units (model B0)."""
    units = document_values(records, entry="molecule_identity_key")
    rows = []
    for entry, group in units.groupby("entry"):
        if len(group) < 2:
            continue
        recs = records[records.molecule_identity_key.eq(entry)]
        for _, unit in group.iterrows():
            docs = set(unit.documents.split(";"))
            remaining = recs[~recs.document_chembl_id.astype(str).isin(docs)]
            displayed = str(klass(np.array([remaining.pchembl_value.astype(float).median()]))[0])
            observed = str(klass(np.array([unit.document_value]))[0])
            rows.append({"entry": entry, "document": unit.document, "displayed": displayed,
                         "outcome": "agree" if observed == displayed else "disagree", "kind": "exact"})
    return pd.DataFrame(rows)


def analyse(target, path):
    raw = pd.read_csv(path, low_memory=False)
    cur = curate_chembl_activity_data(raw, target_id=target, activity_type="IC50")
    records = cur.retained_records.copy()
    records["pchembl_value"] = pd.to_numeric(records.pchembl_value)
    bounds = bounds_table(cur.excluded_records)
    entry_of_parent = records.groupby("parent_molecule_chembl_id").molecule_identity_key.first()
    bounds["entry"] = bounds.entry.map(entry_of_parent).fillna("unmapped:" + bounds.entry)
    mapped = ~bounds.entry.str.startswith("unmapped:")

    # 1. share of bounds among records in the inactive and the active range
    upper = bounds.relation.str.startswith(">")
    inactive_bounds = int((upper & (bounds.p_bound <= LOWER)).sum())
    active_bounds = int((~upper & (bounds.p_bound >= UPPER)).sum())
    exact_inactive = int((records.pchembl_value <= LOWER).sum())
    exact_active = int((records.pchembl_value >= UPPER).sum())
    only_bounds = bounds.loc[~mapped].entry.nunique()

    # 2. entries with bounds from documents without an exact value for them
    med = records.groupby("molecule_identity_key").pchembl_value.median()
    displayed = pd.Series(klass(med.values), index=med.index)
    exact_docs = records.groupby("molecule_identity_key").document_chembl_id.agg(lambda s: set(s.astype(str)))
    b = bounds[mapped].copy()
    b = b[[doc not in exact_docs.get(e, set()) for e, doc in zip(b.entry, b.document)]]
    b["displayed"] = b.entry.map(displayed)
    b["outcome"] = [bound_outcome(d, r, p) for d, r, p in zip(b.displayed, b.relation, b.p_bound)]

    def per_document(g):
        outcomes = set(g.outcome.dropna())
        return outcomes.pop() if len(outcomes) == 1 else None
    bdocs = b.groupby(["entry", "document"]).apply(per_document, include_groups=False).rename("outcome").reset_index()
    bdocs["displayed"] = bdocs.entry.map(displayed)
    bdocs["kind"] = "bound"

    # 3. held-out agreement, exact units only and with bound documents added
    loo = exact_loo(records)
    summary = {"target": target, "raw_records": int(len(raw)), "retained_records": int(len(records)),
               "entries": int(len(med)), "bound_records_usable": int(len(bounds)),
               "bound_records_linked": int(mapped.sum()), "compounds_with_bounds_only": int(only_bounds),
               "inactive_range": {"exact": exact_inactive, "bounds": inactive_bounds,
                                  "bound_share": inactive_bounds / max(1, exact_inactive + inactive_bounds)},
               "active_range": {"exact": exact_active, "bounds": active_bounds,
                                "bound_share": active_bounds / max(1, exact_active + active_bounds)},
               "exact_b0_disagreement": float((loo.outcome == "disagree").mean()) if len(loo) else None,
               "exact_predictions": int(len(loo)), "by_class": {}}
    for cls in ("inactive", "medium", "active"):
        ent = displayed.index[displayed.eq(cls)]
        with_bounds = bdocs[bdocs.displayed.eq(cls)]
        ex = loo[loo.displayed.eq(cls)]
        bd = with_bounds[with_bounds.outcome.notna()]
        n_ex, a_ex = len(ex), int((ex.outcome == "agree").sum())
        n_bd, a_bd = len(bd), int((bd.outcome == "agree").sum())
        summary["by_class"][cls] = {
            "entries": int(len(ent)),
            "entries_with_bound_documents": int(with_bounds.entry.nunique()),
            "bound_documents": int(len(with_bounds)),
            "bound_documents_informative": int(n_bd), "bound_documents_agree": a_bd,
            "bound_documents_disagree": int(n_bd - a_bd),
            "exact_predictions": int(n_ex), "exact_agree": a_ex,
            "exact_agreement": a_ex / n_ex if n_ex else None,
            "bound_agreement": a_bd / n_bd if n_bd else None,
            "combined_agreement": (a_ex + a_bd) / (n_ex + n_bd) if n_ex + n_bd else None,
        }
    return summary, pd.concat([loo, bdocs], ignore_index=True)


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    summaries = []
    for target, path in TARGETS.items():
        summary, table = analyse(target, path)
        table.to_csv(out / f"r15_predictions_{target}.csv", index=False)
        summaries.append(summary)
        print(json.dumps(summary, indent=1))
    (out / "r15_summary.json").write_text(json.dumps(summaries, indent=1) + "\n")


if __name__ == "__main__":
    main(sys.argv[1])
