"""E1-A rank statistics and entry-cluster bootstrap without refitting frozen models."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import rankdata


def auc_mann_whitney(correct, probability):
    correct = np.asarray(correct, dtype=bool)
    n1, n0 = correct.sum(), (~correct).sum()
    if not n1 or not n0:
        return np.nan
    ranks = rankdata(probability, method="average")
    return float((ranks[correct].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def cluster_weights(frame, rng, repetitions):
    """Sample each target's original entry count and apply the same multiplicity to rows of an entry."""
    weights = np.zeros((repetitions, len(frame)), dtype=float)
    for target in sorted(frame.target.unique()):
        positions = np.flatnonzero(frame.target.to_numpy() == target)
        codes, entries = pd.factorize(frame.iloc[positions].entry, sort=True)
        draws = rng.integers(0, len(entries), size=(repetitions, len(entries)))
        counts = np.zeros((repetitions, len(entries)), dtype=int)
        np.add.at(counts, (np.arange(repetitions)[:, None], draws), 1)
        weights[:, positions] = counts[:, codes]
    return weights


def weighted_scores(frame, weights):
    """Multiplicity-weighted frozen metrics; average-rank AUC is equivalent under repetition."""
    p = frame[["p0", "p1", "p2"]].to_numpy(float)
    observed = frame.observed.to_numpy(int)
    displayed = frame.displayed.to_numpy(int)
    confidence = p.max(axis=1)
    correct = (p.argmax(axis=1) == observed).astype(float)
    # The R1 negative control permutes binary agreement without changing classes or probabilities.
    agrees = (frame.agree.to_numpy(float) if "agree" in frame
              else (displayed == observed).astype(float))
    shown = p[np.arange(len(p)), displayed]
    total = weights.sum(axis=1)

    def mean(values):
        return (weights * values).sum(axis=1) / total

    high = confidence >= .8
    high_n = weights[:, high].sum(axis=1)
    result = {
        "accuracy": mean(correct),
        "brier": mean(((p - np.eye(3)[observed]) ** 2).sum(axis=1)),
        "log_loss": mean(-np.log(np.clip(p[np.arange(len(p)), observed], 1e-6, 1))),
        "mean_confidence": mean(confidence), "agreement_rate": mean(agrees),
        "mean_p_displayed": mean(shown), "calibration_bias": mean(shown - agrees),
        "high_confidence_n": high_n,
    }
    for name, values in [("high_confidence_mean_probability", confidence),
                         ("high_confidence_accuracy", correct)]:
        result[name] = np.divide((weights[:, high] * values[high]).sum(axis=1), high_n,
                                out=np.full(len(weights), np.nan), where=high_n > 0)
    bins = np.minimum((confidence * 10).astype(int), 9)
    result["ece"] = sum(np.abs((weights[:, bins == b] * (correct - confidence)[bins == b]).sum(axis=1))
                        for b in range(10)) / total
    # Each tied positive-negative pair contributes 1/2, equivalent to the Mann-Whitney rank sum.
    order = np.argsort(shown, kind="stable")
    starts = np.r_[0, np.flatnonzero(np.diff(shown[order]) != 0) + 1]
    positives = np.add.reduceat(weights[:, order] * agrees[order], starts, axis=1)
    negatives = np.add.reduceat(weights[:, order] * (1 - agrees[order]), starts, axis=1)
    numerator = (positives * (np.cumsum(negatives, axis=1) - .5 * negatives)).sum(axis=1)
    denominator = positives.sum(axis=1) * negatives.sum(axis=1)
    result["auc"] = np.divide(numerator, denominator, out=np.full(len(weights), np.nan), where=denominator > 0)
    return result


def percentile_ci(values):
    finite = np.asarray(values)[np.isfinite(values)]
    if not len(finite):
        return np.nan, np.nan, 0
    low, high = np.percentile(finite, [2.5, 97.5])
    return float(low), float(high), len(finite)


def capped_tau(doc, cap=4):
    moments = doc.groupby("entry").document_value.agg(["size", "var"])
    moments = moments[moments["size"] >= 2]
    weights = np.minimum(moments["size"] - 1, cap)
    return float(np.sqrt((weights * moments["var"]).sum() / weights.sum()))
