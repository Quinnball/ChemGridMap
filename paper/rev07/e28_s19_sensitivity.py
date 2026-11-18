"""S19 sensitivity (post hoc): neighbor evidence with other similarity thresholds, numbers of
neighbors and count-based Morgan fingerprints. Recomputes the tie-inclusive neighbor table
from the rebuilt reference entries; the default (bits, 0.6, k=5) must reproduce the labels of the E28 run."""
import os, sys, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
RDLogger.DisableLog("rdApp.*")
CL = ("inactive", "medium", "active")
HERE = Path(__file__).resolve().parent
# E28_RUN: E28 run folder (hashed-plan run); E28_REF: tables from e28_reference_tables.py; E28_OUT: result folder.
RUN = Path(os.environ.get("E28_RUN", HERE / "e28_run")); REF = Path(os.environ.get("E28_REF", HERE / "e28_ref")); OUT = Path(os.environ.get("E28_OUT", HERE / "output" / "e28"))
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def category(order, sims, status, classes, k, thr):
    if not len(order): return "none"
    cutoff = sims[order[min(k, len(order)) - 1]]
    top = [i for i in order if sims[i] >= cutoff - 1e-12]
    close = [i for i in top if sims[i] >= thr - 1e-12]
    if not close: return "none"
    sup = {classes[i] for i in close if status[i] == "supported"} & set(CL)
    return "insufficient" if not sup else next(iter(sup)) if len(sup) == 1 else "mixed"

def run(sample):
    rows = []
    for d in sorted((RUN / sample).glob("*/compounds.csv")):
        name = d.parent.name; rp = REF / sample / name / "reference.csv"
        if not rp.exists(): print("missing", name); continue
        ref = pd.read_csv(rp); comp = pd.read_csv(d)
        keys = ref.molecule_identity_key.astype(str).tolist()
        mols = [Chem.MolFromSmiles(s) for s in ref.canonical_smiles]
        bits = [gen.GetFingerprint(m) for m in mols]; counts = [gen.GetCountFingerprint(m) for m in mols]
        status = ref.class_status.tolist(); classes = ref.activity_class.tolist()
        rank_key = np.argsort(np.array(keys))  # tie order by identifier
        pos = np.empty(len(keys), int); pos[rank_key] = np.arange(len(keys))
        for q in comp.itertuples():
            m = Chem.MolFromSmiles(q.canonical_smiles)
            out = {"dataset": name, "endpoint": q.endpoint, "query_id": q.query_id, "measured_class": q.measured_class, "saved": q.neighbor_evidence}
            for fp_name, refs, fp in (("bits", bits, gen.GetFingerprint(m)), ("counts", counts, gen.GetCountFingerprint(m))):
                sims = np.asarray(DataStructs.BulkTanimotoSimilarity(fp, refs))
                order = np.lexsort((pos, -sims))
                for k in (3, 5, 10):
                    for thr in (0.5, 0.6, 0.7):
                        if fp_name == "counts" and (k != 5): continue
                        out[f"{fp_name}_k{k}_t{thr}"] = category(order, sims, status, classes, k, thr)
            rows.append(out)
    df = pd.DataFrame(rows); df.to_csv(OUT / f"s19_labels_{sample}.csv.gz", index=False)
    return df

def summarize(df, sample):
    df = df[df.measured_class.isin(CL)]
    res = []
    for col in [c for c in df.columns if c.startswith(("bits_", "counts_"))]:
        for ep, e in df.groupby("endpoint"):
            c = e[e[col].isin(CL)]
            vals = []
            for ds, g in e.groupby("dataset"):
                gc = g[g[col].isin(CL)]
                if len(gc) >= 10:
                    pm = gc.measured_class.value_counts(normalize=True); pc = gc[col].value_counts(normalize=True)
                    vals.append(gc.measured_class.eq(gc[col]).mean() - sum(pm.get(x, 0) * pc.get(x, 0) for x in CL))
            v = np.array(vals); nz = v[v != 0]
            res.append(dict(setting=col, endpoint=ep, none=float(e[col].eq("none").mean()),
                            unsettled=float(e[col].isin(["insufficient", "mixed"]).mean()), consistent=float(e[col].isin(CL).mean()),
                            consistent_n=int(len(c)), agreement=float(c.measured_class.eq(c[col]).mean()),
                            active_agreement=float(c[c[col].eq("active")].measured_class.eq("active").mean()),
                            inactive_agreement=float(c[c[col].eq("inactive")].measured_class.eq("inactive").mean()) if c[col].eq("inactive").any() else np.nan,
                            n1_datasets=len(v), n1_positive=int((v > 0).sum()), n1_median=float(np.median(v)) if len(v) else np.nan,
                            n1_p_two_sided=float(stats.wilcoxon(nz).pvalue) if len(nz) >= 5 else np.nan))
    r = pd.DataFrame(res); r.to_csv(OUT / f"s19_summary_{sample}.csv", index=False)
    check = float((df["bits_k5_t0.6"] == df.saved).mean())
    print("default reproduces E28 labels:", check)
    return r

if __name__ == "__main__":
    s = sys.argv[1] if len(sys.argv) > 1 else "primary"
    summarize(run(s), s)
