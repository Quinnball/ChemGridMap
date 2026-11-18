"""Post hoc E28 analyses for the v0.7 revision (R11, R12, S18, S19, S21, S22).

Inputs: the E28 run outputs (hashed-plan run) (compounds.csv, new_molecules_neighbors.csv) and the
reference tables rebuilt by reconstruct.py. All analyses are post hoc.
"""
import json, os, sys, glob
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
RDLogger.DisableLog("rdApp.*")
CL = ("inactive", "medium", "active")
HERE = Path(__file__).resolve().parent
# E28_RUN: E28 run folder (hashed-plan run); E28_REF: tables from e28_reference_tables.py; E28_OUT: result folder.
RUN = Path(os.environ.get("E28_RUN", HERE / "e28_run")); REF = Path(os.environ.get("E28_REF", HERE / "e28_ref")); OUT = Path(os.environ.get("E28_OUT", HERE / "output" / "e28"))
OUT.mkdir(parents=True, exist_ok=True)
BANDS = [(0.6, 0.7, "0.6-0.7"), (0.7, 0.8, "0.7-0.8"), (0.8, 0.9, "0.8-0.9"), (0.9, 1.0, "0.9-<1"), (1.0, 1.01, "1")]

def band(t):
    if t >= 1 - 1e-12: return "1"
    for lo, hi, name in BANDS[:4]:
        if lo - 1e-12 <= t < hi - 1e-12: return name
    return "<0.6"

_scaf = {}
def murcko(smi):
    if smi in _scaf: return _scaf[smi]
    m = Chem.MolFromSmiles(smi); s = ""
    if m is not None:
        try: s = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m))
        except Exception: s = ""
    _scaf[smi] = s; return s

_flat = {}
def flat(smi):
    if smi in _flat: return _flat[smi]
    m = Chem.MolFromSmiles(smi); s = ""
    if m is not None:
        Chem.RemoveStereochemistry(m); s = Chem.MolToSmiles(m)
    _flat[smi] = s; return s

def load(sample):
    rows, nbrs, refs = [], [], {}
    for d in sorted((RUN / sample).glob("*/compounds.csv")):
        name = d.parent.name
        rp = REF / sample / name / "reference.csv"
        if not rp.exists():
            print("missing reference", sample, name); continue
        c = pd.read_csv(d); c["dataset"] = name
        n = pd.concat([pd.read_csv(f) for f in sorted(d.parent.glob("later_*/new_molecules/*/new_molecules_neighbors.csv"))], ignore_index=True)
        n["dataset"] = name
        r = pd.read_csv(rp); r["dataset"] = name
        ld = pd.read_csv(REF / sample / name / "later_documents.csv")
        c = c.merge(ld, on="query_id", how="left")
        rows.append(c); nbrs.append(n); refs[name] = r
    return pd.concat(rows, ignore_index=True), pd.concat(nbrs, ignore_index=True), refs

def evidence(nb, sim=0.6, k=5, support=True):
    """Neighbor evidence from a tie-inclusive neighbor table (rank = competition rank)."""
    close = nb[(nb["rank"] <= k) & (nb.tanimoto >= sim - 1e-12)] if "rank" in nb else nb[nb.tanimoto >= sim - 1e-12]
    out = {}
    for q, g in close.groupby(["dataset", "query_id"]):
        sup = g[g.class_status.eq("supported")] if support else g
        cls = sorted(set(sup.activity_class.dropna()) & set(CL))
        out[q] = "insufficient" if not cls else cls[0] if len(cls) == 1 else "mixed"
    return out

def n1_dataset(t):
    """Observed minus chance agreement among consistent compounds of one dataset."""
    c = t[t.ev.isin(CL)]
    if len(c) == 0: return np.nan, 0, np.nan, np.nan
    obs = float(c.measured_class.eq(c.ev).mean())
    pm = c.measured_class.value_counts(normalize=True); pc = c.ev.value_counts(normalize=True)
    exp = float(sum(pm.get(x, 0) * pc.get(x, 0) for x in CL))
    return obs - exp, len(c), obs, exp

def n1_test(df, min_n=10, alternative="two-sided"):
    vals = []
    for ds, t in df.groupby("dataset"):
        d, n, o, e = n1_dataset(t)
        if n >= min_n: vals.append(d)
    v = np.array(vals); nz = v[v != 0]
    p = float(stats.wilcoxon(nz, alternative=alternative).pvalue) if len(nz) >= 5 else np.nan
    return dict(datasets=len(v), positive=int((v > 0).sum()), zeros=int((v == 0).sum()),
                median=float(np.median(v)) if len(v) else np.nan, p=p)

def pooled(df):
    c = df[df.ev.isin(CL)]
    return dict(n=int(len(c)), agree=int(c.measured_class.eq(c.ev).sum()),
                agreement=float(c.measured_class.eq(c.ev).mean()) if len(c) else np.nan)

def two_stage(df, value, reps=2000, seed=20261118):
    rng = np.random.default_rng(seed); groups = [g for _, g in df.groupby("dataset")]
    est = []
    for _ in range(reps):
        pick = rng.integers(0, len(groups), len(groups)); num = den = 0
        for i in pick:
            g = groups[i]; idx = rng.integers(0, len(g), len(g)); s = g.iloc[idx]
            num += value(s).sum(); den += len(s)
        est.append(num / den if den else np.nan)
    return float(np.nanpercentile(est, 2.5)), float(np.nanpercentile(est, 97.5))

def main(sample="primary"):
    comp, nb, refs = load(sample)
    comp = comp[comp.measured_class.isin(CL)].copy()
    key = list(zip(comp.dataset, comp.query_id))
    # validation: evidence recomputed from the saved neighbor table equals the saved label
    base = evidence(nb)
    comp["ev"] = [base.get(k, "none") for k in key]
    comp["ev_check"] = comp.ev.eq(comp.neighbor_evidence)
    report = {"sample": sample, "compounds": int(len(comp)), "evidence_reproduced": float(comp.ev_check.mean())}
    # reference lookups
    reftab = pd.concat(refs.values(), ignore_index=True).rename(columns={"molecule_identity_key": "measured_id"})
    reftab["measured_id"] = reftab.measured_id.astype(str)
    nb["measured_id"] = nb.measured_id.astype(str)
    nb = nb.merge(reftab[["dataset", "measured_id", "canonical_smiles", "documents_all", "class_status"]].rename(
        columns={"canonical_smiles": "ref_smiles", "class_status": "ref_status_rebuilt"}), on=["dataset", "measured_id"], how="left")
    report["neighbor_status_rebuilt_match"] = float((nb.class_status == nb.ref_status_rebuilt).mean())
    report["neighbor_smiles_found"] = float(nb.ref_smiles.notna().mean())
    comp["band"] = comp.tanimoto.map(band)
    comp["endpoint"] = comp.endpoint.astype(str)
    # --- R11: nearest-similarity bands and stereo-only identity
    near = nb.sort_values(["dataset", "query_id", "rank", "measured_id"]).groupby(["dataset", "query_id"]).head(1)
    near = near.set_index(["dataset", "query_id"])
    q_smiles = dict(zip(key, comp.canonical_smiles))
    comp["nearest_same_constitution"] = [
        (flat(q_smiles[k]) == flat(near.loc[k, "ref_smiles"])) if (k in near.index and isinstance(near.loc[k, "ref_smiles"], str)) else False
        for k in key]
    r11 = []
    for ep, e in comp.groupby("endpoint"):
        for b in ["<0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-<1", "1"]:
            t = e[e.band.eq(b)]
            p = pooled(t); test = n1_test(t, min_n=10)
            r11.append(dict(endpoint=ep, band=b, compounds=int(len(t)),
                            share=float(len(t) / len(e)), consistent=p["n"], consistent_share=float(p["n"] / len(t)) if len(t) else np.nan,
                            agreement=p["agreement"], same_constitution=int(t.nearest_same_constitution.sum()),
                            n1_datasets=test["datasets"], n1_median=test["median"], n1_positive=test["positive"], n1_p=test["p"]))
    r11 = pd.DataFrame(r11); r11.to_csv(OUT / f"r11_bands_{sample}.csv", index=False)
    excl = []
    for ep, e in comp.groupby("endpoint"):
        for label, mask in (("all", np.ones(len(e), bool)), ("without_T1", ~e.band.eq("1").values),
                            ("without_T>=0.9", ~e.band.isin(["0.9-<1", "1"]).values),
                            ("without_same_constitution", ~e.nearest_same_constitution.values)):
            t = e[mask]; p = pooled(t); test = n1_test(t, min_n=10)
            chance = []
            for ds, g in t.groupby("dataset"):
                d, n, o, x = n1_dataset(g)
                if n: chance.append((n, o, x))
            pooled_chance = sum(n * x for n, o, x in chance) / sum(n for n, o, x in chance)
            med_agree = float(np.median([o for n, o, x in chance if n >= 10]))
            med_chance = float(np.median([x for n, o, x in chance if n >= 10]))
            excl.append(dict(endpoint=ep, subset=label, compounds=int(len(t)), consistent=p["n"], agreement=p["agreement"],
                             pooled_chance=pooled_chance, median_agreement=med_agree, median_chance=med_chance, **{f"n1_{k}": v for k, v in test.items()}))
    excl = pd.DataFrame(excl); excl.to_csv(OUT / f"r11_exclusions_{sample}.csv", index=False)
    # --- R12: series continuation
    sup = nb[(nb.tanimoto >= 0.6 - 1e-12) & nb.class_status.eq("supported")].copy()
    sup["ref_scaffold"] = sup.ref_smiles.map(lambda s: murcko(s) if isinstance(s, str) else "")
    qscaf = {k: murcko(s) for k, s in q_smiles.items()}
    qdocs = {k: set(str(d).split(";")) if isinstance(d, str) else set() for k, d in zip(key, comp.documents)}
    scaf_flag, doc_flag = {}, {}
    for (ds, q), g in sup.groupby(["dataset", "query_id"]):
        k = (ds, q)
        scaf_flag[k] = bool(qscaf.get(k)) and any(qscaf.get(k) == s for s in g.ref_scaffold)
        docs = set().union(*[set(str(x).split(";")) for x in g.documents_all.dropna()]) if g.documents_all.notna().any() else set()
        doc_flag[k] = bool(qdocs.get(k, set()) & docs)
    comp["series_scaffold"] = [scaf_flag.get(k, False) for k in key]
    comp["series_document"] = [doc_flag.get(k, False) for k in key]
    comp["series_any"] = comp.series_scaffold | comp.series_document
    r12 = []
    for ep, e in comp.groupby("endpoint"):
        c = e[e.ev.isin(CL)]
        for label, col in (("scaffold", "series_scaffold"), ("document", "series_document"), ("either", "series_any")):
            for val in (True, False):
                t = c[c[col].eq(val)]; test = n1_test(e[e[col].eq(val) | ~e.ev.isin(CL)][lambda x: x.ev.isin(CL)], min_n=10)
                r12.append(dict(endpoint=ep, indicator=label, continuation=val, consistent=int(len(t)),
                                share=float(len(t) / len(c)), agreement=float(t.measured_class.eq(t.ev).mean()) if len(t) else np.nan,
                                active_share=float(t.ev.eq("active").mean()) if len(t) else np.nan,
                                **{f"n1_{k}": v for k, v in test.items()}))
    r12 = pd.DataFrame(r12); r12.to_csv(OUT / f"r12_series_{sample}.csv", index=False)
    # --- S18: N2 within similarity bands and nearest class; support filter at matched similarity
    comp["near_class"] = comp.nearest_reference_class
    near_status = near["class_status"].to_dict()
    comp["near_status"] = [near_status.get(k, "") for k in key]
    comp["near_agree"] = comp.measured_class.eq(comp.near_class)
    comp["group"] = np.where(comp.ev.isin(CL), "consistent", np.where(comp.ev.isin(["insufficient", "mixed"]), "unsettled", "none"))
    def mh(df, a_mask, b_mask, strata):
        num = den = 0.0; used = 0
        for _, g in df.groupby(strata):
            a = g[a_mask(g)]; b = g[b_mask(g)]
            if len(a) and len(b):
                w = len(a) * len(b) / (len(a) + len(b)); num += w * (a.near_agree.mean() - b.near_agree.mean()); den += w; used += 1
        return (num / den if den else np.nan), used
    s18 = []
    for ep, e in comp.groupby("endpoint"):
        e6 = e[e.tanimoto >= 0.6 - 1e-12]
        cons = lambda g: g.group.eq("consistent"); uns = lambda g: g.group.eq("unsettled")
        crude = e6[cons(e6)].near_agree.mean() - e6[uns(e6)].near_agree.mean()
        mb, ub = mh(e6, cons, uns, ["dataset", "band"])
        mc, uc = mh(e6, cons, uns, ["dataset", "band", "near_class"])
        sp = lambda g: g.near_status.eq("supported"); un = lambda g: ~g.near_status.eq("supported")
        crude_s = e6[sp(e6)].near_agree.mean() - e6[un(e6)].near_agree.mean()
        ms, us = mh(e6, sp, un, ["dataset", "band"])
        msc, usc = mh(e6, sp, un, ["dataset", "band", "near_class"])
        # per-dataset band-standardised N2
        per = []
        for ds, g in e6.groupby("dataset"):
            v, u = mh(g, cons, uns, ["band"])
            if not np.isnan(v) and g[cons(g)].shape[0] >= 5 and g[uns(g)].shape[0] >= 5: per.append(v)
        per = np.array(per); nz = per[per != 0]
        s18.append(dict(endpoint=ep, n2_crude=crude, n2_mh_band=mb, n2_mh_band_strata=ub, n2_mh_band_class=mc, n2_mh_band_class_strata=uc,
                        n2_band_per_dataset_n=len(per), n2_band_per_dataset_median=float(np.median(per)) if len(per) else np.nan,
                        n2_band_per_dataset_positive=int((per > 0).sum()), n2_band_p=float(stats.wilcoxon(nz).pvalue) if len(nz) >= 5 else np.nan,
                        support_crude=crude_s, support_mh_band=ms, support_mh_band_strata=us, support_mh_band_class=msc,
                        nearest_supported_agreement=float(e6[sp(e6)].near_agree.mean()), nearest_uncertain_agreement=float(e6[un(e6)].near_agree.mean()),
                        nearest_supported_n=int(sp(e6).sum()), nearest_uncertain_n=int(un(e6).sum())))
    s18 = pd.DataFrame(s18); s18.to_csv(OUT / f"s18_{sample}.csv", index=False)
    # --- S22: zeros, lift over baselines on the consistent compounds
    s22 = []
    for ep, e in comp.groupby("endpoint"):
        c = e[e.ev.isin(CL)]
        maj = {}
        for ds, r in refs.items():
            maj[ds] = r.activity_class.value_counts().idxmax()
        later_maj = e.groupby("dataset").measured_class.agg(lambda s: s.value_counts().idxmax()).to_dict()
        s22.append(dict(endpoint=ep, consistent=int(len(c)), evidence_agreement=float(c.measured_class.eq(c.ev).mean()),
                        nearest_class_agreement=float(c.measured_class.eq(c.near_class).mean()),
                        reference_majority_agreement=float((c.measured_class == c.dataset.map(maj)).mean()),
                        later_majority_agreement=float((c.measured_class == c.dataset.map(later_maj)).mean()),
                        active_rule_agreement=float(c.measured_class.eq("active").mean()),
                        n1=n1_test(e, min_n=10, alternative="greater")))
    pd.DataFrame(s22).to_csv(OUT / f"s22_{sample}.csv", index=False)
    # --- S21: two-stage intervals of pooled agreement of consistent evidence, by evidence class
    s21 = []
    for ep, e in comp.groupby("endpoint"):
        c = e[e.ev.isin(CL)].copy(); c["hit"] = c.measured_class.eq(c.ev)
        for cls in ("all",) + CL:
            t = c if cls == "all" else c[c.ev.eq(cls)]
            if len(t) == 0: continue
            lo, hi = two_stage(t, lambda s: s.hit)
            per = t.groupby("dataset").hit.agg(["mean", "size"]); per = per[per["size"] >= 10]
            s21.append(dict(endpoint=ep, evidence=cls, n=int(len(t)), datasets=int(t.dataset.nunique()), agreement=float(t.hit.mean()),
                            two_stage_lo=lo, two_stage_hi=hi, median_per_dataset=float(per["mean"].median()) if len(per) else np.nan,
                            datasets_min10=int(len(per))))
        none_share = e.ev.eq("none").mean()
        s21.append(dict(endpoint=ep, evidence="none_share", n=int(len(e)), datasets=int(e.dataset.nunique()), agreement=float(none_share),
                        two_stage_lo=two_stage(e, lambda s: s.ev.eq("none"))[0], two_stage_hi=two_stage(e, lambda s: s.ev.eq("none"))[1],
                        median_per_dataset=float(e.groupby("dataset").ev.agg(lambda s: s.eq("none").mean()).median()), datasets_min10=np.nan))
    pd.DataFrame(s21).to_csv(OUT / f"s21_{sample}.csv", index=False)
    comp.drop(columns=["ev_check"]).to_csv(OUT / f"compounds_annotated_{sample}.csv.gz", index=False)
    nb.to_csv(OUT / f"neighbors_{sample}.csv.gz", index=False)
    report.update(r11_rows=len(r11))
    (OUT / f"report_{sample}.json").write_text(json.dumps(report, indent=1, default=float))
    print(json.dumps(report, default=float))

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "primary")
