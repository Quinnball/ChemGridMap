"""E4a 事后分析：比较对象、显示类别校准和类别支持度。"""
from __future__ import annotations

import gzip
import hashlib
import json
import pickle
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[4]
REV = ROOT / "project1/JCIM_reliability_revision_20261024"
S = REV / "strengthening"
sys.path[:0] = [
    str(S / "source"),
    str(ROOT / "chemgridmap/paper"),
    str(ROOT / "chemgridmap/src"),
]

import numpy as np
import pandas as pd
from scipy import stats

from chemgridmap.evidence import class_probabilities, document_values
from e1a_statistics import (
    auc_mann_whitney,
    cluster_weights,
    percentile_ci,
)

E2 = S / "output/E2"
OUT = S / "output/E4"
QA = S / "qa/E4"
POPULATION = "all"
ENDPOINTS = ("IC50", "Ki")
BOOT_N = 2000
BOOTSTRAP_SEED = 20261103
CURVE_BINS = 5
CURVE_MIN_PREDICTIONS = 8
ORIGINAL_KI_TARGETS = ("CHEMBL205", "CHEMBL204")
ALIGN = ["fold", "entry", "held_out_unit"]


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, float_precision="round_trip", low_memory=False)


def load_records() -> dict[tuple[str, str], pd.DataFrame]:
    path = S / "qa/E2/curated_records.pkl.gz"
    return pickle.loads(gzip.decompress(path.read_bytes()))


def displayed_ece(
    score,
    agree,
    weights: np.ndarray | None = None,
) -> float:
    """显示类别的 10 箱 ECE；权重无需对每行归一。"""
    score = np.asarray(score, dtype=float)
    agree = np.asarray(agree, dtype=float)
    if weights is None:
        weights = np.ones((1, len(score)), dtype=float)
    weights = np.asarray(weights, dtype=float)
    bins = np.minimum((score * 10).astype(int), 9)
    total = weights.sum(axis=1)
    numerator = np.zeros(len(weights), dtype=float)
    for index in range(10):
        chosen = bins == index
        numerator += np.abs(
            (weights[:, chosen] * (agree[chosen] - score[chosen])).sum(axis=1)
        )
    return np.divide(
        numerator,
        total,
        out=np.full(len(weights), np.nan),
        where=total > 0,
    )


def bootstrap_binary_auc(agree, score, frame, weights) -> np.ndarray:
    """在固定目标条目重数下计算二分类 AUC。"""
    agree = np.asarray(agree, dtype=bool)
    score = np.asarray(score, dtype=float)
    order = np.argsort(score, kind="stable")
    starts = np.r_[0, np.flatnonzero(np.diff(score[order]) != 0) + 1]
    positives = np.add.reduceat(
        weights[:, order] * agree[order], starts, axis=1
    )
    negatives = np.add.reduceat(
        weights[:, order] * (1 - agree[order]), starts, axis=1
    )
    numerator = (
        positives * (np.cumsum(negatives, axis=1) - 0.5 * negatives)
    ).sum(axis=1)
    denominator = positives.sum(axis=1) * negatives.sum(axis=1)
    return np.divide(
        numerator,
        denominator,
        out=np.full(len(weights), np.nan),
        where=denominator > 0,
    )


def auc_cluster_ci(agree, score, frame) -> tuple[float, float, float, int]:
    point = auc_mann_whitney(agree, score)
    if len(frame) == 0:
        return point, np.nan, np.nan, 0
    local = frame[["target", "entry"]].copy()
    weights = cluster_weights(local, np.random.default_rng(BOOTSTRAP_SEED), BOOT_N)
    draws = bootstrap_binary_auc(agree, score, local, weights)
    low, high, valid = percentile_ci(draws)
    return point, low, high, valid


def bootstrap_median_ci(values) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan, np.nan
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    draws = np.empty(BOOT_N, dtype=float)
    for start in range(0, BOOT_N, 100):
        end = min(start + 100, BOOT_N)
        indices = rng.integers(0, len(values), size=(end - start, len(values)))
        draws[start:end] = np.median(values[indices], axis=1)
    return tuple(float(value) for value in np.percentile(draws, [2.5, 97.5]))


def one_sided_wilcoxon(values) -> tuple[float, float, int]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0 or np.allclose(values, 0.0):
        return np.nan, 1.0, len(values)
    result = stats.wilcoxon(values, alternative="greater", zero_method="wilcox")
    return float(result.statistic), float(result.pvalue), len(values)


def distance_score(center: float, displayed: int) -> float:
    """离显示类别边界的距离；越大表示类别内位置越深。"""
    if displayed == 0:
        return float(6.0 - center)
    if displayed == 2:
        return float(center - 7.0)
    if displayed == 1:
        return float(min(center - 6.0, 7.0 - center))
    raise ValueError(f"未知显示类别：{displayed}")


def distance_scores(center, displayed) -> np.ndarray:
    return np.array(
        [
            distance_score(float(value), int(shown))
            for value, shown in zip(center, displayed)
        ],
        dtype=float,
    )


def shown_probability(center, scale, shown) -> np.ndarray:
    probabilities = class_probabilities(
        np.atleast_1d(center),
        np.atleast_1d(scale),
        np.inf,
        6.0,
        7.0,
    )
    shown = np.broadcast_to(np.asarray(shown, dtype=int), len(probabilities))
    return probabilities[np.arange(len(probabilities)), shown]


def g_with_support(
    predictions: pd.DataFrame,
    records_by_key: dict[tuple[str, str], pd.DataFrame],
) -> pd.DataFrame:
    """从同一组剩余单元重建 G 中心和显示类别支持度。"""
    g = predictions.loc[
        (predictions.population == POPULATION) & (predictions.method == "G")
    ].copy()
    g["remaining_center"] = np.nan
    for (target, endpoint), target_rows in g.groupby(
        ["target", "endpoint"], sort=True
    ):
        doc = document_values(
            records_by_key[(target, endpoint)], merge_identical=True
        )
        units = dict(tuple(doc.groupby("entry", sort=True)))
        for index, row in target_rows.iterrows():
            available = units[row.entry]
            held = set(str(row.held_out_unit).split(";"))
            remaining = available.loc[
                ~available.documents.map(
                    lambda value: bool(set(str(value).split(";")) & held)
                )
            ]
            if len(remaining) != int(row.m):
                raise ValueError(
                    f"剩余单元数不一致：{target}/{endpoint}/{row.entry}"
                )
            g.at[index, "remaining_center"] = float(
                remaining.document_value.mean()
            )
    g["display_support"] = shown_probability(
        g.remaining_center.to_numpy(float),
        g.tau.to_numpy(float) / np.sqrt(g.m.to_numpy(float)),
        g.displayed.to_numpy(int),
    )
    return g


def comparator_table(
    predictions: pd.DataFrame,
    g: pd.DataFrame,
    comparators: pd.DataFrame,
    metrics: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for (target, endpoint), target_g in g.groupby(
        ["target", "endpoint"], sort=True
    ):
        target_g = target_g.sort_values(ALIGN).reset_index(drop=True)
        b0 = predictions.loc[
            (predictions.target == target)
            & (predictions.endpoint == endpoint)
            & (predictions.population == POPULATION)
            & (predictions.method == "B0")
        ].sort_values(ALIGN).reset_index(drop=True)
        fixed = predictions.loc[
            (predictions.target == target)
            & (predictions.endpoint == endpoint)
            & (predictions.population == POPULATION)
            & (predictions.method == "G-fixed")
        ].sort_values(ALIGN).reset_index(drop=True)
        if not target_g[ALIGN].equals(fixed[ALIGN]):
            raise ValueError(f"G-fixed 留出键不一致：{target}/{endpoint}")

        b_sd = b0.loc[b0.remaining_records >= 2].reset_index(drop=True)
        archived = comparators.loc[
            (comparators.target == target)
            & (comparators.endpoint == endpoint)
            & (comparators.population == POPULATION)
            & (comparators.model == "B-sd")
        ]
        if len(archived) != 1:
            raise ValueError(f"缺少比较对象：{target}/{endpoint}/B-sd")
        archived = archived.iloc[0]
        scores = {
            "G": (target_g.agree, target_g.p_displayed, target_g),
            "D": (
                target_g.agree,
                distance_scores(
                    target_g.remaining_center,
                    target_g.displayed,
                ),
                target_g,
            ),
            "G-fixed": (fixed.agree, fixed.p_displayed, fixed),
            "B-sd": (b_sd.agree, b_sd.b_sd, b_sd),
        }
        for model, (agree, score, frame) in scores.items():
            point, low, high, valid = auc_cluster_ci(agree, score, frame)
            if model == "B-sd":
                both_missing = np.isnan(point) and np.isnan(archived.auc)
                same_value = np.isclose(
                    point, archived.auc, rtol=0, atol=1e-12
                )
                if not (both_missing or same_value):
                    raise ValueError(f"B-sd AUC 未复现：{target}/{endpoint}")
            rows.append({
                "target": target,
                "endpoint": endpoint,
                "population": POPULATION,
                "comparator": model,
                "n_predictions": len(frame),
                "n_entries": int(frame.entry.nunique()),
                "auc": point,
                "auc_ci_low": low,
                "auc_ci_high": high,
                "ci_fully_above_half": bool(low > 0.5),
                "bootstrap_valid": valid,
                "bootstrap_repetitions": BOOT_N,
                "bootstrap_seed": BOOTSTRAP_SEED,
                "analysis_role": "post hoc",
            })

    result = pd.DataFrame(rows)
    expected = metrics.loc[
        (metrics.population == POPULATION)
        & (metrics.stratum == "all")
        & metrics.model.isin(("G", "G-fixed")),
        ["target", "endpoint", "model", "auc"],
    ].copy()
    actual = result.loc[
        result.comparator.isin(("G", "G-fixed")),
        ["target", "endpoint", "comparator", "auc"],
    ].rename(columns={"comparator": "model"})
    check = expected.merge(
        actual,
        on=["target", "endpoint", "model"],
        suffixes=("_archive", "_recomputed"),
        validate="one_to_one",
    )
    if len(check) != len(expected):
        raise ValueError("G/G-fixed 与 metrics_by_target.csv 的键不一致")
    error = float(
        (check.auc_archive - check.auc_recomputed).abs().max()
    )
    if error > 1e-12:
        raise ValueError(f"G/G-fixed AUC 未复现，最大误差 {error}")
    result.attrs["maximum_metric_auc_error"] = error
    return result


def comparison_summary(
    comparators: pd.DataFrame,
    metrics: pd.DataFrame,
    forest: pd.DataFrame,
) -> tuple[pd.DataFrame, dict]:
    rows = []
    for endpoint in ENDPOINTS:
        endpoint_frame = comparators.loc[comparators.endpoint == endpoint]
        for comparator in ("G", "D", "G-fixed", "B-sd"):
            selected = endpoint_frame.loc[endpoint_frame.comparator == comparator]
            rows.append({
                "row_type": "auc_ci_above_half",
                "endpoint": endpoint,
                "comparison": comparator,
                "n_datasets": len(selected),
                "n_auc_valid": int(selected.auc.notna().sum()),
                "n_fully_above_half": int(selected.ci_fully_above_half.sum()),
                "bootstrap_repetitions": BOOT_N,
                "bootstrap_seed": BOOTSTRAP_SEED,
                "analysis_role": "post hoc",
            })

        wide = endpoint_frame.pivot(
            index="target", columns="comparator", values="auc"
        )
        for comparator in ("D", "G-fixed"):
            differences = (wide["G"] - wide[comparator]).dropna()
            statistic, p_value, n_targets = one_sided_wilcoxon(differences)
            ci_low, ci_high = bootstrap_median_ci(differences)
            rows.append({
                "row_type": "paired_auc_difference",
                "endpoint": endpoint,
                "comparison": f"G-{comparator}",
                "n_datasets": n_targets,
                "median_difference": float(np.median(differences)),
                "ci_low": ci_low,
                "ci_high": ci_high,
                "wilcoxon_statistic": statistic,
                "raw_p": p_value,
                "alternative": "greater",
                "multiple_correction": "none",
                "bootstrap_repetitions": BOOT_N,
                "bootstrap_seed": BOOTSTRAP_SEED,
                "analysis_role": "post hoc",
            })

    g = metrics.loc[
        (metrics.population == POPULATION)
        & (metrics.stratum == "all")
        & (metrics.model == "G")
    ]
    ki = g.loc[g.endpoint == "Ki"].copy()
    excluded = ki.loc[~ki.target.isin(ORIGINAL_KI_TARGETS)]
    statistic, p_value, n_targets = one_sided_wilcoxon(excluded.auc - 0.5)
    rows.append({
        "row_type": "ki_excluding_original_targets_h3",
        "endpoint": "Ki",
        "comparison": "G_AUC_vs_0.5",
        "n_datasets": n_targets,
        "median_auc": float(excluded.auc.median()),
        "wilcoxon_statistic": statistic,
        "raw_p": p_value,
        "alternative": "greater",
        "excluded_targets": ";".join(ORIGINAL_KI_TARGETS),
        "analysis_role": "post hoc",
    })
    result = pd.DataFrame(rows)

    forest_counts = {}
    checks = {}
    for endpoint in ENDPOINTS:
        archived = forest.loc[
            (forest.population == POPULATION) & (forest.endpoint == endpoint)
        ]
        archived_count = int((archived.g_auc_ci_low > 0.5).sum())
        computed = int(
            (
                comparators.loc[
                    (comparators.endpoint == endpoint)
                    & (comparators.comparator == "G"),
                    "ci_fully_above_half",
                ]
            ).sum()
        )
        if computed != archived_count:
            raise ValueError(
                f"G AUC CI 计数未与 forest 一致：{endpoint}"
            )
        forest_counts[endpoint] = {
            "n_datasets": len(archived),
            "g_ci_fully_above_half": computed,
        }
        checks[f"g_ci_count_matches_forest_{endpoint}"] = True
    return result, {
        "forest_g_ci_counts": forest_counts,
        "checks": checks,
    }


def calibration_outputs(
    predictions: pd.DataFrame,
    g: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    b0 = predictions.loc[
        (predictions.population == POPULATION) & (predictions.method == "B0")
    ].copy()
    detail_rows = []
    pooled_rows = []
    support_rows = []
    for endpoint in ENDPOINTS:
        for model, frame in (
            ("G", g.loc[g.endpoint == endpoint]),
            ("B0", b0.loc[b0.endpoint == endpoint]),
        ):
            for target, target_frame in frame.groupby("target", sort=True):
                detail_rows.append({
                    "row_type": "target_displayed_ece",
                    "endpoint": endpoint,
                    "target": target,
                    "model": model,
                    "metric": "displayed_ece",
                    "n_predictions": len(target_frame),
                    "n_entries": int(target_frame.entry.nunique()),
                    "estimate": float(
                        displayed_ece(
                            target_frame.p_displayed,
                            target_frame.agree,
                        )[0]
                    ),
                    "analysis_role": "post hoc",
                })
            endpoint_frame = frame.reset_index(drop=True)
            weights = cluster_weights(
                endpoint_frame[["target", "entry"]],
                np.random.default_rng(BOOTSTRAP_SEED),
                BOOT_N,
            )
            draws = displayed_ece(
                endpoint_frame.p_displayed,
                endpoint_frame.agree,
                weights,
            )
            low, high, valid = percentile_ci(draws)
            pooled_rows.append({
                "row_type": "pooled_displayed_ece",
                "endpoint": endpoint,
                "target": "",
                "model": model,
                "metric": "displayed_ece",
                "n_predictions": len(endpoint_frame),
                "n_entries": int(
                    endpoint_frame[["target", "entry"]].drop_duplicates().shape[0]
                ),
                "n_targets": int(endpoint_frame.target.nunique()),
                "estimate": float(
                    displayed_ece(
                        endpoint_frame.p_displayed,
                        endpoint_frame.agree,
                    )[0]
                ),
                "ci_low": low,
                "ci_high": high,
                "bootstrap_valid": valid,
                "bootstrap_repetitions": BOOT_N,
                "bootstrap_seed": BOOTSTRAP_SEED,
                "analysis_role": "post hoc",
            })

        support = g.loc[g.endpoint == endpoint].reset_index(drop=True)
        for target, target_frame in support.groupby("target", sort=True):
            point, low, high, valid = auc_cluster_ci(
                target_frame.agree,
                target_frame.display_support,
                target_frame,
            )
            support_rows.append({
                "row_type": "target_support_auc",
                "endpoint": endpoint,
                "target": target,
                "model": "display_support",
                "metric": "support_auc",
                "n_predictions": len(target_frame),
                "n_entries": int(target_frame.entry.nunique()),
                "estimate": point,
                "ci_low": low,
                "ci_high": high,
                "bootstrap_valid": valid,
                "bootstrap_repetitions": BOOT_N,
                "bootstrap_seed": BOOTSTRAP_SEED,
                "analysis_role": "post hoc",
            })
        point, low, high, valid = auc_cluster_ci(
            support.agree,
            support.display_support,
            support,
        )
        support_rows.append({
            "row_type": "pooled_support_auc",
            "endpoint": endpoint,
            "target": "",
            "model": "display_support",
            "metric": "support_auc",
            "n_predictions": len(support),
            "n_entries": int(
                support[["target", "entry"]].drop_duplicates().shape[0]
            ),
            "n_targets": int(support.target.nunique()),
            "estimate": point,
            "ci_low": low,
            "ci_high": high,
            "bootstrap_valid": valid,
            "bootstrap_repetitions": BOOT_N,
            "bootstrap_seed": BOOTSTRAP_SEED,
            "analysis_role": "post hoc",
        })

        support_bins = np.minimum(
            (support.display_support * 10).astype(int), 9
        )
        for index in range(10):
            selected = support.loc[support_bins == index]
            support_rows.append({
                "row_type": "support_bin",
                "endpoint": endpoint,
                "target": "",
                "model": "display_support",
                "metric": "mean_support_vs_agreement",
                "n_predictions": len(selected),
                "n_entries": int(
                    selected[["target", "entry"]].drop_duplicates().shape[0]
                ) if len(selected) else 0,
                "bin": index + 1,
                "bin_low": index / 10,
                "bin_high": (index + 1) / 10,
                "mean_score": (
                    float(selected.display_support.mean())
                    if len(selected)
                    else np.nan
                ),
                "agreement_rate": (
                    float(selected.agree.mean())
                    if len(selected)
                    else np.nan
                ),
                "analysis_role": "post hoc",
            })

    calibration = pd.DataFrame(detail_rows + pooled_rows + support_rows)
    calibration = calibration.sort_values(
        ["endpoint", "row_type", "model", "target", "bin"],
        na_position="last",
    ).reset_index(drop=True)

    curve_rows = []
    for endpoint in ENDPOINTS:
        for model, frame in (
            ("G", g.loc[g.endpoint == endpoint]),
            ("B0", b0.loc[b0.endpoint == endpoint]),
        ):
            bin_index = np.minimum(
                (frame.p_displayed * CURVE_BINS).astype(int),
                CURVE_BINS - 1,
            )
            for index in range(CURVE_BINS):
                selected = frame.loc[bin_index == index]
                eligible = len(selected) >= CURVE_MIN_PREDICTIONS
                curve_rows.append({
                    "endpoint": endpoint,
                    "model": model,
                    "bin": index + 1,
                    "bin_low": index / CURVE_BINS,
                    "bin_high": (index + 1) / CURVE_BINS,
                    "n_predictions": len(selected),
                    "mean_probability": (
                        float(selected.p_displayed.mean())
                        if eligible
                        else np.nan
                    ),
                    "observed_agreement": (
                        float(selected.agree.mean())
                        if eligible
                        else np.nan
                    ),
                    "eligible": eligible,
                    "minimum_predictions": CURVE_MIN_PREDICTIONS,
                    "analysis_role": "post hoc",
                })
    curve = pd.DataFrame(curve_rows)

    qa_details = g.loc[
        :,
        [
            "target",
            "endpoint",
            "fold",
            "entry",
            "held_out_unit",
            "m",
            "displayed",
            "observed",
            "agree",
            "tau",
            "remaining_center",
            "p_displayed",
            "display_support",
        ],
    ].copy()
    qa_details["analysis_role"] = "post hoc"
    return calibration, curve, qa_details


def write_outputs() -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    predictions = read_csv(E2 / "predictions.csv")
    metrics = read_csv(E2 / "metrics_by_target.csv")
    comparators_archive = read_csv(E2 / "comparators.csv")
    forest = read_csv(E2 / "forest_plot_data.csv")
    inclusion = read_csv(E2 / "target_inclusion.csv")
    included = inclusion.loc[inclusion.included.astype(bool)]
    if len(included) != 52:
        raise ValueError(f"纳入数据集不是 52 个，而是 {len(included)} 个")

    records = load_records()
    g = g_with_support(predictions, records)
    comparator = comparator_table(
        predictions,
        g,
        comparators_archive,
        metrics,
    )
    summary, summary_checks = comparison_summary(
        comparator,
        metrics,
        forest,
    )
    calibration, curve, qa_details = calibration_outputs(predictions, g)

    comparator.to_csv(OUT / "comparators_posthoc.csv", index=False)
    calibration.to_csv(OUT / "calibration_displayed.csv", index=False)
    curve.to_csv(OUT / "calibration_curve_displayed.csv", index=False)
    summary.to_csv(OUT / "e4a_comparison_summary.csv", index=False)
    qa_details.to_csv(QA / "e4a_g_support_details.csv", index=False)

    expected_forest = {
        "IC50": {"n_datasets": 29, "g_ci_fully_above_half": 24},
        "Ki": {"n_datasets": 23, "g_ci_fully_above_half": 19},
    }
    observed_forest = summary_checks["forest_g_ci_counts"]
    if observed_forest != expected_forest:
        raise ValueError(
            f"forest CI 计数不符合任务书：{observed_forest}"
        )
    ki = summary.loc[summary.row_type == "ki_excluding_original_targets_h3"].iloc[0]
    if len(ki) and (
        abs(float(ki.median_auc) - 0.747) > 5e-4
        or float(ki.raw_p) > 5e-7
    ):
        raise ValueError(
            f"Ki 去除原靶点结果不符：median={ki.median_auc}, p={ki.raw_p}"
        )

    checks = {
        "included_datasets": len(included),
        "maximum_metric_auc_error": comparator.attrs["maximum_metric_auc_error"],
        "forest_g_ci_counts": observed_forest,
        "ki_excluding_original_targets": {
            "n_datasets": int(ki.n_datasets),
            "median_auc": float(ki.median_auc),
            "wilcoxon_p": float(ki.raw_p),
        },
        "output_rows": {
            "comparators_posthoc.csv": len(comparator),
            "calibration_displayed.csv": len(calibration),
            "calibration_curve_displayed.csv": len(curve),
            "e4a_comparison_summary.csv": len(summary),
            "e4a_g_support_details.csv": len(qa_details),
        },
        "input_sha256": {
            name: hashlib.sha256((E2 / name).read_bytes()).hexdigest()
            for name in (
                "predictions.csv",
                "metrics_by_target.csv",
                "comparators.csv",
                "forest_plot_data.csv",
                "target_inclusion.csv",
            )
        },
        "seeds": {
            "bootstrap": BOOTSTRAP_SEED,
            "bootstrap_repetitions": BOOT_N,
        },
        "checks": summary_checks["checks"],
        "passed": True,
    }
    (QA / "e4a_checks.json").write_text(
        json.dumps(checks, indent=2, ensure_ascii=False) + "\n"
    )
    return checks


if __name__ == "__main__":
    result = write_outputs()
    print(json.dumps(result, indent=2, ensure_ascii=False))
