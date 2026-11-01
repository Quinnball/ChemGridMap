"""Validate and apply evidence-aware activity confidence.

1. Fit the document-variance model on each full curated target population.
2. Leave-one-document-out prediction, cross-fitted over 10 entry folds:
   for every entry with >= 2 documents, predict the display class of each
   held-out document from the remaining documents. Compared methods:
   B0 current display (class of the remaining-record median, probability 1),
   G Gaussian with the pooled between-unit SD (default), MF/MP moderated t
   models fitted by variance likelihood or leave-one-unit-out predictive likelihood.
3. Annotate the 3289 mapped entries and analyse colour contrasts between
   grid neighbours on the seed-42 layouts (minimum-displacement and refined).
"""
import json

import numpy as np
import pandas as pd

from common import LAYOUTS, RESULTS, TARGETS, load_layout, space_neighbors
from chemgridmap.activity import activity_class
from chemgridmap.evidence import (CLASSES, annotate_activity_confidence, class_probabilities,
                                  document_values, fit_between_unit_sd, fit_document_variance,
                                  fit_predictive_variance, predictive_parameters)

OUT = RESULTS / "activity_confidence"
OUT.mkdir(parents=True, exist_ok=True)
LOWER, UPPER = 6.0, 7.0
cls = lambda v: CLASSES.index(activity_class(v, LOWER, UPPER))  # noqa: E731


def full_records(target):
    return pd.read_csv(LAYOUTS / target / f"{target}_retained_records.csv", low_memory=False)


def lodo(target, records, folds=10, seed=20261024):
    doc = document_values(records)
    counts = doc.groupby("entry").size()
    multi = counts[counts >= 2].index.to_numpy()
    rng = np.random.default_rng(seed)
    fold_of = dict(zip(rng.permutation(multi), np.arange(len(multi)) % folds))
    rows = []
    for fold in range(folds):
        test = {e for e, f in fold_of.items() if f == fold}
        train = doc[~doc.entry.isin(test)]
        fits = {"G": fit_between_unit_sd(train), "MF": fit_document_variance(train), "MP": fit_predictive_variance(train)}
        for entry in sorted(test):
            e_docs = doc[doc.entry == entry]
            e_recs = records[records.molecule_identity_key == entry]
            for _, held in e_docs.iterrows():
                rest = e_docs[e_docs.document != held.document].document_value.to_numpy(float)
                rest_records = e_recs[~e_recs.document_chembl_id.astype(str).isin(held.documents.split(";"))].pchembl_value
                m = len(rest)
                own = max(np.var(rest, ddof=1), 1e-4) if m > 1 else 0.0
                row = {"target": target, "fold": fold, "entry": entry, "held_out_unit": held.documents,
                       "remaining_units": m, "observed": cls(held.document_value)}
                probs = {"B0": np.eye(3)[cls(float(np.median(rest_records)))]}
                for name, fit in fits.items():
                    loc, scale, df = predictive_parameters(np.array([rest.mean()]), np.array([own]), np.array([float(m)]),
                                                           fit["tau"], fit["prior_df"])
                    probs[name] = class_probabilities(loc, scale, df, LOWER, UPPER)[0]
                for name, p in probs.items():
                    row.update({f"{name}_p{k}": p[k] for k in range(3)})
                rows.append(row)
    return pd.DataFrame(rows)


def scores(frame, method):
    p = frame[[f"{method}_p{k}" for k in range(3)]].to_numpy(float)
    y = np.eye(3)[frame.observed.to_numpy(int)]
    conf, pred = p.max(1), p.argmax(1)
    correct = (pred == frame.observed.to_numpy()).astype(float)
    bins = np.minimum((conf * 10).astype(int), 9)
    ece = sum(abs(correct[bins == b].mean() - conf[bins == b].mean()) * (bins == b).mean()
              for b in range(10) if (bins == b).any())
    return {"method": method, "n_predictions": len(frame), "n_entries": frame.entry.nunique(),
            "accuracy": correct.mean(), "brier": float(((p - y) ** 2).sum(1).mean()),
            "log_loss": float(-np.log(np.clip((p * y).sum(1), 1e-6, 1)).mean()),
            "mean_confidence": float(conf.mean()), "ece": float(ece)}


validation, fits, annotated = [], {}, []
for target, (label, _) in TARGETS.items():
    records = full_records(target)
    fits[label] = fit_between_unit_sd(document_values(records))
    frame = lodo(target, records)
    frame["target"] = label
    validation.append(frame)
    mapped = load_layout(target, 42)
    out, _ = annotate_activity_confidence(mapped, records, lower=LOWER, upper=UPPER, fit=fits[label])
    out["target"] = label
    annotated.append(out)
    print(label, {k: round(v, 3) if isinstance(v, float) else v for k, v in fits[label].items()}, flush=True)

validation = pd.concat(validation, ignore_index=True)
validation.to_csv(OUT / "lodo_predictions.csv", index=False)
score_rows = [{"target": t, **scores(f, m)} for t, f in [("all", validation), *validation.groupby("target")]
              for m in ("B0", "G", "MF", "MP")]
score_table = pd.DataFrame(score_rows)
score_table.to_csv(OUT / "lodo_scores.csv", index=False)
print(score_table.round(3).to_string(index=False))

# Reliability diagram for the class shown from the remaining documents (M vs B1).
rel = []
for method in ("G", "MF", "MP"):
    p = validation[[f"{method}_p{k}" for k in range(3)]].to_numpy(float)
    shown = validation[[f"B0_p{k}" for k in range(3)]].to_numpy(float).argmax(1)
    p_shown = p[np.arange(len(p)), shown]
    hit = (shown == validation.observed.to_numpy()).astype(float)
    edges = np.linspace(0, 1, 6)
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (p_shown >= lo) & (p_shown < hi if hi < 1 else p_shown <= hi)
        if sel.any():
            rel.append({"method": method, "bin_low": lo, "bin_high": hi, "n": int(sel.sum()),
                        "mean_predicted": float(p_shown[sel].mean()), "observed_frequency": float(hit[sel].mean())})
pd.DataFrame(rel).to_csv(OUT / "reliability_shown_class.csv", index=False)
print(pd.DataFrame(rel).round(3).to_string(index=False))

annotated = pd.concat(annotated, ignore_index=True)
keep = ["target", "molecule_identity_key", "molecule_id", "activity_pchembl", "activity_class", "n_records",
        "n_documents", "n_evidence_units", "top_unit_documents", "top_unit_record_fraction", "document_balanced_pchembl",
        "document_balanced_class", "class_changed_by_document_balance", "consensus_scale",
        *[f"p_consensus_{c}" for c in CLASSES], "class_confidence", "class_reproducibility", "class_status"]
annotated[keep].to_csv(OUT / "mapped_entry_confidence.csv", index=False)

summary = []
for label, g in annotated.groupby("target", sort=False):
    for scope, sel in [("all", np.ones(len(g), bool)), ("one_unit", g.n_evidence_units == 1),
                       ("multi_unit", g.n_evidence_units >= 2)]:
        h = g[sel]
        summary.append({"target": label, "scope": scope, "entries": len(h),
                        "supported_fraction": float((h.class_status == "supported").mean()),
                        "median_confidence": float(h.class_confidence.median()),
                        "class_changed_by_document_balance": int(h.class_changed_by_document_balance.sum()),
                        **{f"supported_{c}": float((h[h.activity_class == c].class_status == "supported").mean())
                           if (h.activity_class == c).any() else np.nan for c in CLASSES}})
summary = pd.DataFrame(summary)
summary.to_csv(OUT / "mapped_confidence_summary.csv", index=False)
print(summary.round(3).to_string(index=False))

# Colour contrasts between grid neighbours on the seed-42 layouts.
layouts = pd.read_csv(RESULTS / "grid_refinement/layouts_seed42.csv")
contrast_rows = []
for target, (label, _) in TARGETS.items():
    g = annotated[annotated.target == label].reset_index(drop=True)
    probs = g[[f"p_consensus_{c}" for c in CLASSES]].to_numpy(float)
    shown = g.activity_class.map(CLASSES.index).to_numpy()
    for layout_name, lay in layouts[layouts.target == label].groupby("layout"):
        lay = lay.set_index("molecule_identity_key").loc[g.molecule_identity_key]
        nb = space_neighbors(lay[["grid_col", "grid_row"]].to_numpy(float))
        pairs = {(min(i, j), max(i, j)) for i, row in enumerate(nb) for j in row}
        i, j = np.array(sorted(pairs)).T
        differ = 1.0 - (probs[i] * probs[j]).sum(1)
        shown_diff = shown[i] != shown[j]
        contrast_rows.append({"target": label, "layout": layout_name, "neighbor_pairs": len(i),
                              "displayed_contrasts": int(shown_diff.sum()),
                              "expected_true_contrasts_among_displayed": float(differ[shown_diff].sum()),
                              "supported_contrasts_p_ge_0.8": int((differ[shown_diff] >= 0.8).sum()),
                              "expected_hidden_contrasts_same_colour": float(differ[~shown_diff].sum())})
contrasts = pd.DataFrame(contrast_rows)
contrasts.to_csv(OUT / "neighbor_contrasts_seed42.csv", index=False)
print(contrasts.round(1).to_string(index=False))

case = annotated[annotated.molecule_id == "CHEMBL20"].iloc[0]
(OUT / "fits_and_case.json").write_text(json.dumps({
    "fits": fits, "acetazolamide": {k: (float(case[k]) if isinstance(case[k], (float, np.floating, np.integer, int)) else str(case[k]))
                                     for k in keep if k in case}}, indent=1))
print(json.dumps({k: str(case[k]) for k in keep}, indent=0))
