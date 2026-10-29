"""E1-A 的秩统计与按条目整群抽样；不拟合或修改冻结模型。"""
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
    """每个靶点分别抽取其原有条目数；相同条目的各行使用相同重数。"""
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
    """冻结指标的重数加权形式；AUC 在重复观测的平均秩下等价。"""
    p = frame[["p0", "p1", "p2"]].to_numpy(float)
    observed = frame.observed.to_numpy(int)
    displayed = frame.displayed.to_numpy(int)
    confidence = p.max(axis=1)
    correct = (p.argmax(axis=1) == observed).astype(float)
    # R1 负控只置换二分类一致指标，不改动三分类 observed 或概率。
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
    # 同分块内每个正负配对贡献 1/2；块间贡献与 Mann–Whitney 秩和完全相同。
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
