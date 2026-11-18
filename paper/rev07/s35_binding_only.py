"""S35 (post hoc): assay-type composition, and H1 and sigma_doc with binding assays only.

For each of the 52 evaluation datasets this script curates the downloaded records twice with unchanged
rules, once with all assay types (as in the evaluation) and once with binding assays only (ChEMBL assay
type B), and reports for each version the number of entries with two or more evidence units, sigma_doc
and the H1 disagreement rate: each evidence unit of an entry with at least two units is held out in
turn and compared with the class of the median of the remaining records (model B0). The predeclared
analyses used the earlier identity rule, so both versions use it here (--identity-rule
normalized-structure). Datasets with fewer than 30 multi-unit entries after the binding filter are
reported but not summarized.
Usage: PYTHONPATH=src python paper/rev07/s35_binding_only.py <folder with <target>_<endpoint>.csv> <inclusion.csv> <output folder>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from chemgridmap import curate_chembl_activity_data, document_values, fit_between_unit_sd

LOWER, UPPER = 6.0, 7.0


def klass(p):
    return np.where(p <= LOWER, 0, np.where(p >= UPPER, 2, 1))


def b0_disagreement(records):
    units = document_values(records, entry="molecule_identity_key")
    multi = units.groupby("entry").filter(lambda g: len(g) >= 2)
    if multi.empty:
        return np.nan, 0, 0
    recs = records[records.molecule_identity_key.isin(multi.entry.unique())]
    by_entry = {e: g for e, g in recs.groupby("molecule_identity_key")}
    disagree = total = 0
    for entry, group in multi.groupby("entry"):
        r = by_entry[entry]
        docs = r.document_chembl_id.astype(str)
        values = r.pchembl_value.astype(float)
        for _, unit in group.iterrows():
            remaining = values[~docs.isin(unit.documents.split(";"))]
            total += 1
            disagree += int(klass(remaining.median()) != klass(unit.document_value))
    return disagree / total, total, multi.entry.nunique()


def one(raw, target, endpoint, assay_types=None):
    cur = curate_chembl_activity_data(raw, target_id=target, activity_type=endpoint, assay_types=assay_types,
                                      identity_rule="normalized-structure", identity_conflicts="exclude")
    records = cur.retained_records.copy()
    records["pchembl_value"] = pd.to_numeric(records.pchembl_value)
    h1, predictions, multi = b0_disagreement(records)
    sigma = fit_between_unit_sd(document_values(records, entry="molecule_identity_key"), fallback_sd=None) \
        if multi >= 5 else {"tau": np.nan}
    composition = records.assay_type.fillna("missing").value_counts(normalize=True).to_dict()
    return {"retained_records": int(len(records)), "multi_unit_entries": int(multi), "predictions": int(predictions),
            "h1": float(h1), "sigma_doc": float(sigma["tau"]), "composition": composition}


def main(folder, inclusion, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    inc = pd.read_csv(inclusion)
    inc = inc[inc.included.astype(str).str.lower().eq("true")]
    rows = []
    for _, row in inc.iterrows():
        path = Path(folder) / f"{row.target}_{row.endpoint}.csv"
        raw = pd.read_csv(path, low_memory=False)
        both = one(raw, row.target, row.endpoint)
        binding = one(raw, row.target, row.endpoint, assay_types=["B"])
        rows.append({"target": row.target, "endpoint": row.endpoint,
                     "share_B": both["composition"].get("B", 0.0), "share_F": both["composition"].get("F", 0.0),
                     "share_A": both["composition"].get("A", 0.0),
                     "all_multi_unit_entries": both["multi_unit_entries"], "all_predictions": both["predictions"],
                     "all_h1": both["h1"], "all_sigma_doc": both["sigma_doc"],
                     "binding_multi_unit_entries": binding["multi_unit_entries"],
                     "binding_predictions": binding["predictions"], "binding_h1": binding["h1"],
                     "binding_sigma_doc": binding["sigma_doc"]})
        print(rows[-1], flush=True)
    table = pd.DataFrame(rows)
    table.to_csv(out / "s35_by_dataset.csv", index=False)
    summary = {}
    for endpoint, t in table.groupby("endpoint"):
        q = t[t.binding_multi_unit_entries >= 30]
        dh = q.binding_h1 - q.all_h1
        ds = q.binding_sigma_doc - q.all_sigma_doc
        summary[endpoint] = {
            "datasets": int(len(t)), "median_share_binding": float(t.share_B.median()),
            "datasets_all_binding": int((t.share_B >= 0.999).sum()),
            "datasets_qualifying_binding_only": int(len(q)),
            "median_h1_all": float(t.all_h1.median()), "median_sigma_all": float(t.all_sigma_doc.median()),
            "median_h1_all_qualifying": float(q.all_h1.median()), "median_h1_binding": float(q.binding_h1.median()),
            "median_sigma_all_qualifying": float(q.all_sigma_doc.median()),
            "median_sigma_binding": float(q.binding_sigma_doc.median()),
            "median_h1_change": float(dh.median()), "median_sigma_change": float(ds.median()),
            "h1_change_p_two_sided": float(stats.wilcoxon(dh[dh != 0]).pvalue) if (dh != 0).sum() >= 5 else None,
            "sigma_change_p_two_sided": float(stats.wilcoxon(ds[ds != 0]).pvalue) if (ds != 0).sum() >= 5 else None,
        }
    (out / "s35_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:4])
