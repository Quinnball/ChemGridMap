"""S20 (post hoc): supported-versus-uncertain agreement difference, stratified by dataset and displayed class.

Mantel-Haenszel risk difference with the Sato variance, a dataset bootstrap (2,000 resamples of
datasets), the crude difference, per-class differences and per-dataset class-standardised differences.
Input: strengthening/output/E5b/predictions_e5b.csv (model G rows, same-value merge basis).
Usage: python paper/rev07/s20_mh.py <predictions_e5b.csv> <output folder>
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def mh(t):
    a, n1, c, n0 = t.k_supported.values, t.n_supported.values, t.k_uncertain.values, t.n_uncertain.values
    N = n1 + n0
    W = n1 * n0 / N
    rd = np.sum(W * (a / n1 - c / n0)) / W.sum()
    P = (n1**2 * c - n0**2 * a + n1 * n0 * (n0 - n1) / 2) / N**2
    Q = (a * (n0 - c) + c * (n1 - a)) / (2 * N)
    var = (rd * P.sum() + Q.sum()) / W.sum() ** 2
    return rd, np.sqrt(max(var, 0))


def main(predictions, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(predictions, low_memory=False)
    d = d[d.method.eq("G") & d.merge_basis.eq("same_value")].copy()
    d["design"] = d.split_design + d.cutoff_year.fillna(0).astype(int).astype(str).str.replace("^0$", "", regex=True)
    g = d.groupby(["design", "target", "endpoint", "displayed", "display_status"]).agg(n=("agree", "size"), k=("agree", "sum")).reset_index()
    w = g.pivot_table(index=["design", "target", "endpoint", "displayed"], columns="display_status", values=["n", "k"], fill_value=0)
    w.columns = [f"{a}_{b}" for a, b in w.columns]
    w = w.reset_index()
    w = w[(w.n_supported > 0) & (w.n_uncertain > 0)]
    w.to_csv(out / "mh_strata.csv", index=False)
    rng = np.random.default_rng(20261118)
    rows = []
    for (des, ep), t in w.groupby(["design", "endpoint"]):
        rd, se = mh(t)
        crude = t.k_supported.sum() / t.n_supported.sum() - t.k_uncertain.sum() / t.n_uncertain.sum()
        targets = t.target.unique()
        bs = []
        for _ in range(2000):
            s = rng.choice(targets, len(targets))
            bs.append(mh(pd.concat([t[t.target.eq(x)] for x in s]))[0])
        per = {}
        for cl, tc in t.groupby("displayed"):
            r, _s = mh(tc)
            per[int(cl)] = (round(r, 3), int(tc.target.nunique()), int(tc.n_supported.sum() + tc.n_uncertain.sum()))
        pdv = [round(mh(tt)[0], 4) for _, tt in t.groupby("target") if tt.n_supported.sum() >= 5 and tt.n_uncertain.sum() >= 5]
        rows.append(dict(design=des, endpoint=ep, mh_rd=round(rd, 4), sato_lo=round(rd - 1.96 * se, 4), sato_hi=round(rd + 1.96 * se, 4),
                         boot_lo=round(float(np.percentile(bs, 2.5)), 4), boot_hi=round(float(np.percentile(bs, 97.5)), 4),
                         crude=round(crude, 4), strata=int(len(t)), datasets=int(t.target.nunique()), per_class=per, per_dataset=pdv))
    (out / "s20_mh.json").write_text(json.dumps(rows, indent=1) + "\n")
    for r in rows:
        print(r)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
