"""Cross-fit evidence models and annotate the frozen mapped populations."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import spsolve

from chemgridmap.activity import activity_class
from chemgridmap.evidence import (
    CLASSES, annotate_activity_confidence, class_probabilities, document_pair_offsets,
    document_values, fit_between_unit_sd, fit_document_variance,
    fit_predictive_variance, predictive_parameters,
)
from reliability_analysis import FROZEN, OUTPUT, TARGETS, dump_json, write_manifest

LOWER, UPPER = 6.0, 7.0
SEED = 20261024
LAMBDAS = (0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0)
METHODS = ("B0", "G", "MF", "MP", "O", "X", "N")


def _json_safe(value):
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, np.generic):
        value = value.item()
    return None if isinstance(value, float) and not np.isfinite(value) else value


def entry_folds(doc, folds=10, seed=SEED):
    """Assign repeated entries exactly as the reference Gaussian protocol."""
    counts = doc.groupby("entry").size()
    multi = counts[counts >= 2].index.to_numpy()
    return dict(zip(np.random.default_rng(seed).permutation(multi), np.arange(len(multi)) % folds))


def repeated_units(doc):
    return doc[doc.groupby("entry").entry.transform("size") >= 2].copy()


def fit_document_ridge(doc, penalty):
    """Fit unpenalized entry means and penalized document offsets by centering."""
    doc = repeated_units(doc)
    if doc.empty:
        raise ValueError("Document ridge needs entries with at least two units.")
    entry_codes, entries = pd.factorize(doc.entry, sort=True)
    document_codes, documents = pd.factorize(doc.document, sort=True)
    n = len(doc)
    entry_design = sparse.csr_matrix((np.ones(n), (np.arange(n), entry_codes)), shape=(n, len(entries)))
    document_design = sparse.csr_matrix((np.ones(n), (np.arange(n), document_codes)), shape=(n, len(documents)))
    counts = np.asarray(entry_design.sum(axis=0)).ravel()
    cross = entry_design.T @ document_design
    normal = (document_design.T @ document_design - cross.T @ sparse.diags(1.0 / counts) @ cross
              + penalty * sparse.eye(len(documents), format="csc"))
    values = doc.document_value.to_numpy(float)
    entry_means = np.asarray(entry_design.T @ values).ravel() / counts
    rhs = document_design.T @ (values - entry_means[entry_codes])
    delta = spsolve(normal.tocsc(), rhs)
    means = np.asarray(entry_design.T @ (values - delta[document_codes])).ravel() / counts
    return {"penalty": float(penalty), "offsets": dict(zip(documents, delta)),
            "entry_means": dict(zip(entries, means)), "n_entries": len(entries),
            "n_documents": len(documents), "n_units": n}


def corrected_values(doc, offsets):
    """Use each unit's representative document; unseen offsets are zero."""
    return doc.document_value.to_numpy(float) - doc.document.map(offsets).fillna(0).to_numpy(float)


def select_ridge_penalty(doc, folds=5, seed=SEED):
    """Select lambda by entry-cross-fitted leave-one-unit-out squared error."""
    doc = repeated_units(doc)
    assignments = entry_folds(doc, folds, seed)
    rows = []
    for fold in range(folds):
        test_entries = {e for e, f in assignments.items() if f == fold}
        train = doc[~doc.entry.isin(test_entries)]
        test = doc[doc.entry.isin(test_entries)]
        for penalty in LAMBDAS:
            fitted = fit_document_ridge(train, penalty)
            adjusted = pd.Series(corrected_values(test, fitted["offsets"]), index=test.index)
            sums = adjusted.groupby(test.entry).transform("sum")
            counts = adjusted.groupby(test.entry).transform("size")
            # Rest-unit means are computed only after the held unit is removed.
            rest = (sums - adjusted) / (counts - 1)
            predicted = rest.to_numpy() + test.document.map(fitted["offsets"]).fillna(0).to_numpy()
            squared = (test.document_value.to_numpy(float) - predicted) ** 2
            rows.append({"inner_fold": fold, "lambda": penalty, "n_predictions": len(squared),
                         "squared_error_sum": float(squared.sum()), "mse": float(squared.mean())})
    cv = pd.DataFrame(rows)
    totals = cv.groupby("lambda")[["squared_error_sum", "n_predictions"]].sum()
    totals["mse"] = totals.squared_error_sum / totals.n_predictions
    selected = float(totals.mse.idxmin())
    fitted = fit_document_ridge(doc, selected)
    adjusted = doc.copy()
    adjusted["document_value"] = corrected_values(doc, fitted["offsets"])
    fitted["variance_fit"] = fit_between_unit_sd(adjusted)
    fitted["inner_entry_folds"] = assignments
    return fitted, cv


def flagged_documents(records):
    table = document_pair_offsets(records)
    flagged = table[table.systematic_offset.astype(bool)]
    return set(flagged.document_a) | set(flagged.document_b)


def unflagged_units(doc, flagged):
    """Exclude a merged unit if any constituent document is flagged."""
    keep = doc.documents.map(lambda text: not (set(text.split(";")) & flagged))
    return doc.loc[keep]


def gaussian_prediction(values, fit, held_offset=0.0):
    values = np.asarray(values, dtype=float)
    own = max(values.var(ddof=1), 1e-4) if len(values) > 1 else 0.0
    center, scale, df = predictive_parameters([values.mean() + held_offset], [own],
                                              [len(values)], fit["tau"], fit["prior_df"])
    return class_probabilities(center, scale, df, LOWER, UPPER)[0]


def excluded_prediction(rest_doc, flagged, fit, fallback):
    remaining = unflagged_units(rest_doc, flagged)
    if remaining.empty:
        return np.asarray(fallback), 0, True
    return gaussian_prediction(remaining.document_value.to_numpy(float), fit), len(remaining), False


def lodo_models(target, records):
    """Evaluate all models; N retains its separate unmerged evaluation set."""
    doc = document_values(records)
    raw_doc = document_values(records, merge_identical=False)
    assignments = entry_folds(doc)
    rows, cv_rows, offset_rows, provenance = [], [], [], []
    by_entry = {entry: frame for entry, frame in records.groupby("molecule_identity_key")}
    for fold in range(10):
        test = {entry for entry, assigned in assignments.items() if assigned == fold}
        train = doc[~doc.entry.isin(test)]
        raw_train = raw_doc[~raw_doc.entry.isin(test)]
        fits = {"G": fit_between_unit_sd(train), "MF": fit_document_variance(train),
                "MP": fit_predictive_variance(train)}
        ridge, cv = select_ridge_penalty(raw_train)
        flagged = flagged_documents(records[~records.molecule_identity_key.isin(test)])
        fits["X"] = fit_between_unit_sd(unflagged_units(train, flagged))
        cv_rows.append(cv.assign(target=target, outer_fold=fold))
        offset_rows.extend({"target": target, "outer_fold": fold, "document": d, "delta": v,
                            "lambda": ridge["penalty"]} for d, v in ridge["offsets"].items())
        provenance.append({"fold": fold, "test_entries": sorted(test), "fits": fits,
                           "O": {k: v for k, v in ridge.items() if k != "entry_means"},
                           "X_flagged_documents": sorted(flagged)})
        for entry in sorted(test):
            units = doc[doc.entry == entry]
            source_records = by_entry[entry]
            for held in units.itertuples(index=False):
                rest_doc = units[units.document != held.document]
                rest = rest_doc.document_value.to_numpy(float)
                rest_records = source_records.loc[~source_records.document_chembl_id.astype(str).isin(
                    held.documents.split(";")), "pchembl_value"]
                displayed = CLASSES.index(activity_class(float(np.median(rest_records)), LOWER, UPPER))
                observed = CLASSES.index(activity_class(held.document_value, LOWER, UPPER))
                probabilities = {"B0": np.eye(3)[displayed]}
                probabilities.update({m: gaussian_prediction(rest, fits[m]) for m in ("G", "MF", "MP")})
                probabilities["O"] = gaussian_prediction(corrected_values(rest_doc, ridge["offsets"]),
                    ridge["variance_fit"], ridge["offsets"].get(held.document, 0.0))
                probabilities["X"], x_used, x_fallback = excluded_prediction(
                    rest_doc, flagged, fits["X"], probabilities["G"])
                for method, p in probabilities.items():
                    rows.append({"target": target, "target_label": TARGETS[target], "method": method,
                        "evaluation_set": "merged", "fold": fold, "entry": entry,
                        "held_out_unit": held.documents, "held_out_document": held.document,
                        "held_out_value": held.document_value, "remaining_units": len(rest),
                        "used_units": x_used if method == "X" else len(rest),
                        "fallback_to_G": method == "X" and x_fallback,
                        "observed": observed, "displayed": displayed,
                        **{f"p{k}": float(p[k]) for k in range(3)}})
        print(f"activity {target} fold={fold}: lambda={ridge['penalty']:g}", flush=True)
    n_assignments = entry_folds(raw_doc)
    n_fits = []
    for fold in range(10):
        test = {entry for entry, assigned in n_assignments.items() if assigned == fold}
        fit = fit_between_unit_sd(raw_doc[~raw_doc.entry.isin(test)])
        n_fits.append({"fold": fold, "test_entries": sorted(test), "fit": fit})
        for entry in sorted(test):
            units = raw_doc[raw_doc.entry == entry]
            source_records = by_entry[entry]
            for held in units.itertuples(index=False):
                rest = units.loc[units.document != held.document, "document_value"].to_numpy(float)
                rest_records = source_records.loc[source_records.document_chembl_id.astype(str) != held.document,
                                                   "pchembl_value"]
                p = gaussian_prediction(rest, fit)
                rows.append({"target": target, "target_label": TARGETS[target], "method": "N",
                    "evaluation_set": "unmerged", "fold": fold, "entry": entry,
                    "held_out_unit": held.documents, "held_out_document": held.document,
                    "held_out_value": held.document_value, "remaining_units": len(rest), "used_units": len(rest),
                    "fallback_to_G": False,
                    "observed": CLASSES.index(activity_class(held.document_value, LOWER, UPPER)),
                    "displayed": CLASSES.index(activity_class(float(np.median(rest_records)), LOWER, UPPER)),
                    **{f"p{k}": float(p[k]) for k in range(3)}})
    return pd.DataFrame(rows), {"outer_folds": provenance, "N_outer_folds": n_fits}, pd.concat(cv_rows), offset_rows


def probability_scores(frame):
    p = frame[["p0", "p1", "p2"]].to_numpy(float)
    observed = frame.observed.to_numpy(int)
    confidence = p.max(axis=1)
    correct = (p.argmax(axis=1) == observed).astype(float)
    bins = np.minimum((confidence * 10).astype(int), 9)
    ece = sum(abs(correct[bins == b].mean() - confidence[bins == b].mean()) * (bins == b).mean()
              for b in range(10) if (bins == b).any())
    high = confidence >= .8
    return {"n_predictions": len(frame), "n_entries": int(frame.entry.nunique()),
            "accuracy": float(correct.mean()), "brier": float(((p - np.eye(3)[observed]) ** 2).sum(1).mean()),
            "log_loss": float(-np.log(np.clip(p[np.arange(len(p)), observed], 1e-6, 1)).mean()),
            "ece": float(ece), "mean_confidence": float(confidence.mean()),
            "high_confidence_n": int(high.sum()),
            "high_confidence_mean_probability": float(confidence[high].mean()) if high.any() else np.nan,
            "high_confidence_accuracy": float(correct[high].mean()) if high.any() else np.nan}


def calibration_bins(frame):
    """Bin probability of the displayed class, including empty bins."""
    p = frame[["p0", "p1", "p2"]].to_numpy(float)
    displayed = frame.displayed.to_numpy(int)
    confidence = p[np.arange(len(p)), displayed]
    correct = (displayed == frame.observed.to_numpy(int)).astype(float)
    bins = np.minimum((confidence * 5).astype(int), 4)
    rows = []
    for b in range(5):
        selected = bins == b
        n = int(selected.sum())
        rows.append({"bin": b, "bin_low": b / 5, "bin_high": (b + 1) / 5, "n": n,
                     "x": float(confidence[selected].mean()) if n else np.nan,
                     "y": float(correct[selected].mean()) if n else np.nan, "plot_eligible": n >= 8})
    return rows


def evaluate_predictions(predictions):
    metrics, calibration = [], []
    for target, frame in [("all", predictions), *predictions.groupby("target", sort=False)]:
        for method, selected in frame.groupby("method", sort=False):
            base = {"target": target, "target_label": TARGETS.get(target, "All"), "method": method,
                    "evaluation_set": "unmerged" if method == "N" else "merged"}
            metrics.append({**base, **probability_scores(selected)})
            if method == "B0":
                calibration.append({**base, "bin": 5, "bin_low": 1., "bin_high": 1.,
                    "x": 1., "y": float((selected.observed == selected.displayed).mean()),
                    "n": len(selected), "plot_eligible": len(selected) >= 8})
            else:
                calibration.extend({**base, **row} for row in calibration_bins(selected))
    return pd.DataFrame(metrics), pd.DataFrame(calibration)


def adoption_decision(metrics):
    """Apply every preregistered inequality without a tolerance."""
    indexed = metrics.set_index(["target", "method"])
    rows, candidates = [], []
    for method in ("O", "X"):
        rules = [("overall_brier", "all", "brier", 0., "le"),
                 ("overall_ece", "all", "ece", 0., "le"),
                 ("herg_high_confidence_accuracy", "chembl240", "high_confidence_accuracy", .05, "ge"),
                 ("caii_ece", "chembl205", "ece", .02, "le"),
                 ("thrombin_ece", "chembl204", "ece", .02, "le")]
        passes = []
        for name, target, metric, margin, direction in rules:
            baseline = float(indexed.loc[(target, "G"), metric])
            candidate = float(indexed.loc[(target, method), metric])
            boundary = baseline + margin
            passed = bool(candidate <= boundary if direction == "le" else candidate >= boundary)
            passes.append(passed)
            rows.append({"method": method, "rule": name, "target": target, "metric": metric,
                "G_value": baseline, "candidate_value": candidate, "difference": candidate - baseline,
                "required_boundary": boundary, "comparison": "<=" if direction == "le" else ">=", "passed": passed})
        if all(passes):
            candidates.append(method)
    chosen = min(candidates, key=lambda m: indexed.loc[("all", m), "brier"]) if candidates else "G"
    return {"default_model": chosen, "eligible_models": candidates, "N_excluded_from_adoption": True,
            "rules": rows, "decision": "Retain G; neither candidate meets all rules." if chosen == "G"
            else f"Adopt {chosen}; all rules pass and overall Brier is lowest among eligible candidates."}, pd.DataFrame(rows)


def _apply_adjusted_annotation(base, units, fit):
    """Update consensus quantities for an adjusted set of evidence units."""
    moments = units.groupby("entry").document_value.agg(["mean", "size"]).reindex(base.molecule_identity_key)
    scale = fit["tau"] / np.sqrt(moments["size"].to_numpy(float))
    predictive_scale = fit["tau"] * np.sqrt(1 + 1 / moments["size"].to_numpy(float))
    mean = moments["mean"].to_numpy(float)
    consensus = class_probabilities(mean, scale, np.inf, LOWER, UPPER)
    predictive = class_probabilities(mean, predictive_scale, np.inf, LOWER, UPPER)
    out = base.copy()
    out["document_balanced_pchembl"] = mean
    out["document_balanced_class"] = [activity_class(v, LOWER, UPPER) for v in mean]
    out["class_changed_by_document_balance"] = (out.document_balanced_class != out.activity_class).astype(int)
    out["consensus_scale"] = scale
    out["model_units_used"] = moments["size"].to_numpy(int)
    for k, name in enumerate(CLASSES):
        out[f"p_consensus_{name}"] = consensus[:, k]
        out[f"p_new_document_{name}"] = predictive[:, k]
    index = np.array([CLASSES.index(c) for c in out.activity_class])
    out["class_confidence"] = consensus[np.arange(len(out)), index]
    out["class_reproducibility"] = predictive[np.arange(len(out)), index]
    out["class_status"] = np.where(out.class_confidence >= .8, "supported", "uncertain")
    return out


def annotate_target(mapped, records, model, merge_identical=True):
    """Fit on complete retained records; annotate only frozen mapped entries."""
    doc = document_values(records, merge_identical=merge_identical)
    gaussian_fit = fit_between_unit_sd(doc)
    out, fit = annotate_activity_confidence(mapped, records, fit=gaussian_fit,
                                           merge_identical=merge_identical)
    out["model_units_used"] = out.n_evidence_units
    out["fallback_to_G"] = False
    if model == "O":
        ridge, cv = select_ridge_penalty(document_values(records, merge_identical=False))
        adjusted = doc.copy()
        adjusted["document_value"] = corrected_values(doc, ridge["offsets"])
        out = _apply_adjusted_annotation(out, adjusted, ridge["variance_fit"])
        fit.update({"model": "document_offset", "ridge": ridge, "inner_cv": cv.to_dict("records")})
    elif model == "X":
        flagged = flagged_documents(records)
        remaining = unflagged_units(doc, flagged)
        x_fit = fit_between_unit_sd(remaining)
        usable = out.molecule_identity_key.isin(remaining.entry)
        updated = _apply_adjusted_annotation(out.loc[usable], remaining, x_fit)
        out.loc[usable, updated.columns] = updated
        out.loc[~usable, "fallback_to_G"] = True
        fit.update({"model": "exclude_offset_documents", "excluded_fit": x_fit,
                    "flagged_documents": sorted(flagged)})
    elif model != "G":
        raise ValueError(f"Unknown annotation model: {model}")
    out["model"] = model
    out["evidence_policy"] = "merged" if merge_identical else "unmerged"
    return out, fit


def support_summary(frame):
    rows = []
    for target, data in [("all", frame), *frame.groupby("target", sort=False)]:
        for policy, group in data.groupby("evidence_policy", sort=False):
            for scope, selected in [("all", group), ("one_unit", group[group.n_evidence_units == 1]),
                                    ("multi_unit", group[group.n_evidence_units >= 2])]:
                for display in ("all", *CLASSES):
                    subset = selected if display == "all" else selected[selected.activity_class == display]
                    for threshold in (.7, .8, .9):
                        supported = int((subset.class_confidence >= threshold).sum())
                        rows.append({"target": target, "target_label": TARGETS.get(target, "All"),
                            "evidence_policy": policy, "scope": scope, "display_class": display,
                            "threshold": threshold, "n": len(subset), "supported_n": supported,
                            "supported_fraction": supported / len(subset) if len(subset) else np.nan,
                            "class_changed_by_document_balance": int(subset.class_changed_by_document_balance.sum()),
                            "discordant_document_entries": int((subset.discordant_document_pairs > 0).sum())})
    return pd.DataFrame(rows)


def run_activity(output_root=OUTPUT):
    out = Path(output_root) / "activity"
    out.mkdir(parents=True, exist_ok=True)
    predictions, cv_tables, offset_rows, provenance, inputs = [], [], [], {}, []
    records_by_target = {}
    for target in TARGETS:
        path = FROZEN / target / f"{target}_retained_records.csv"
        inputs.append(path)
        records = pd.read_csv(path, low_memory=False)
        records_by_target[target] = records
        frame, fits, cv, offsets = lodo_models(target, records)
        predictions.append(frame)
        provenance[target] = fits
        cv_tables.append(cv)
        offset_rows.extend(offsets)
    predictions = pd.concat(predictions, ignore_index=True)
    metrics, calibration = evaluate_predictions(predictions)
    decision, rules = adoption_decision(metrics)
    baseline = metrics.query("target == 'all' and method == 'G'").iloc[0]
    if not (baseline.n_predictions == 820 and baseline.n_entries == 289
            and abs(baseline.ece - .0767537203890331) < 1e-12
            and abs(baseline.brier - .5218617813580732) < 1e-12):
        raise RuntimeError("Gaussian cross-validation does not reproduce the frozen baseline.")
    predictions.to_csv(out / "lodo_predictions.csv", index=False)
    metrics.to_csv(out / "lodo_scores.csv", index=False)
    calibration.to_csv(out / "calibration.csv", index=False)
    rules.to_csv(out / "adoption_rules.csv", index=False)
    pd.concat(cv_tables, ignore_index=True).to_csv(out / "ridge_inner_cv.csv", index=False)
    pd.DataFrame(offset_rows).to_csv(out / "ridge_document_offsets.csv", index=False)
    dump_json(out / "adoption.json", _json_safe(decision))
    print(f"activity default_model={decision['default_model']}", flush=True)
    annotated, n_annotated, full_fits, pair_tables, count_rows = [], [], {}, [], []
    for target, records in records_by_target.items():
        path = FROZEN / target / "map_molecules.csv"
        inputs.append(path)
        mapped = pd.read_csv(path, low_memory=False)
        frame, fit = annotate_target(mapped, records, decision["default_model"])
        separate, n_fit = annotate_target(mapped, records, decision["default_model"], merge_identical=False)
        frame = frame.assign(target=target, target_label=TARGETS[target])
        separate = separate.assign(target=target, target_label=TARGETS[target])
        frame.to_csv(out / f"{target}_annotations.csv", index=False)
        annotated.append(frame)
        n_annotated.append(separate)
        full_fits[target] = {"default": fit, "unmerged": n_fit}
        offsets = document_pair_offsets(records)
        flagged = offsets[offsets.systematic_offset.astype(bool)].copy()
        counts = {"target": target, "target_label": TARGETS[target],
            "all_possible_pairs": offsets.attrs["document_pair_possible_count"],
            "shared_at_least_one_pairs": offsets.attrs["observed_document_pair_count"],
            "shared_threshold_pairs": offsets.attrs["document_pair_count"], "flagged_pairs": len(flagged),
            "mapped_discordant_entries": int((frame.discordant_document_pairs > 0).sum())}
        count_rows.append(counts)
        flagged["target"], flagged["target_label"] = target, TARGETS[target]
        flagged.to_csv(out / f"{target}_document_offsets.csv", index=False)
        pair_tables.append(flagged)
    annotated = pd.concat(annotated, ignore_index=True)
    n_annotated = pd.concat(n_annotated, ignore_index=True)
    if len(annotated) != 3289 or len(n_annotated) != 3289:
        raise RuntimeError("Mapped annotation count must be 3289 for each evidence policy.")
    annotated.to_csv(out / "mapped_entry_confidence.csv", index=False)
    n_annotated.to_csv(out / "mapped_entry_confidence_unmerged.csv", index=False)
    summaries = support_summary(pd.concat([annotated, n_annotated], ignore_index=True))
    summaries.to_csv(out / "threshold_sensitivity.csv", index=False)
    summaries.query("threshold == 0.8").to_csv(out / "mapped_confidence_summary.csv", index=False)
    pairs = pd.concat(pair_tables, ignore_index=True)
    pairs.to_csv(out / "document_offsets.csv", index=False)
    pd.DataFrame(count_rows).to_csv(out / "document_pair_counts.csv", index=False)
    pair = pairs.query("target == 'chembl240' and document_a == 'CHEMBL5725903' and document_b == 'CHEMBL5727624'")
    if len(pair) != 1 or pair.iloc[0].shared_entries != 30 or abs(pair.iloc[0].median_offset - 3.47) > 1e-12:
        raise RuntimeError("The preregistered hERG document pair was not reproduced.")
    parameters = {"lower": LOWER, "upper": UPPER, "outer_folds": 10, "inner_folds": 5, "seed": SEED,
        "lambda_grid": LAMBDAS, "lambda_objective": "mean squared held-out unit error",
        "inner_split": "entry-level cross-fitting followed by leave-one-unmerged-unit-out prediction",
        "O_training": "unmerged document units of repeated entries; unpenalized entry means",
        "O_tau": "G pooled SD on training unmerged document values minus fitted document offsets",
        "O_merged_unit_offset": "representative document: lexicographically first constituent document",
        "O_unseen_document_offset": 0.,
        "X_flagged_document": "either endpoint of a default-threshold flagged document pair",
        "X_merged_unit_exclusion": "exclude if any constituent document is flagged",
        "X_empty_prediction": "fallback to G; retain every held-out unit in evaluation",
        "N": "separate unmerged evaluation set; not eligible for adoption",
        "calibration": "probability of B0 displayed class; five equal-width bins; n>=8 to plot",
        "ECE": "maximum class probability; ten equal-width bins",
        "default_model": decision["default_model"]}
    dump_json(out / "model_fits.json", _json_safe({"parameters": parameters,
        "cross_validation": provenance, "full_population": full_fits}))
    write_manifest(out, inputs, parameters)
    return decision


if __name__ == "__main__":
    run_activity()
