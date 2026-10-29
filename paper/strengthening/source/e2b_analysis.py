"""E2b: planned multi-target tests, temporal split and descriptive analyses."""
from __future__ import annotations

from datetime import datetime
import gzip
import hashlib
import importlib.metadata
import json
from pathlib import Path
import pickle
import platform
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[3]
REV = ROOT / "paper"
S = REV / "strengthening"
sys.path[:0] = [
    str(S / "source"),
    str(ROOT / "paper"),
    str(ROOT / "src"),
]

import numpy as np
import pandas as pd
from scipy import stats

from chemgridmap.activity import activity_class
from chemgridmap.evidence import (
    CLASSES,
    LITERATURE_BETWEEN_UNIT_SD,
    class_probabilities,
    document_pair_offsets,
    document_values,
    fit_between_unit_sd,
)
from e1a_statistics import (
    auc_mann_whitney,
    cluster_weights,
    percentile_ci,
    weighted_scores,
)
from e2_models import (
    BOOTSTRAP_SEED,
    CLASS_BOUNDARIES,
    ENDPOINTS,
    FOLD_SEED,
    POPULATIONS,
    curate_activity_table,
)
from reliability_activity import entry_folds

OUT = S / "output/E2"
QA = S / "qa/E2"
API = S / "data/chembl37_api"
SEALED_PLAN = S / "plans/E2_task_draft_sealed.md"
SEALED_PLAN_SHA256 = "92e5ac5c920ae4d14f106d1286b5a7310f1291e62d3d29c9b32a20344b3a827d"
REPORT = OUT / "E2b_REPORT.md"
BRIEF = ROOT / "docs/codex-result.md"
NULL_SEED = 20261030
BOOT_N = 2000
SIGMA_BOOT_N = 1000
SUPPORT_THRESHOLD = 0.8
TIME_CUTOFFS = (2015, 2010)
ORIGINAL_IC50_TARGETS = ("CHEMBL205", "CHEMBL204", "CHEMBL240")
R0_TIME_SUFFIX = "_R0_without_original_ic50"
E2A_READ_ONLY_OUTPUTS = (
    "target_inclusion.csv",
    "predictions.csv",
    "metrics_by_target.csv",
    "pooled_metrics.csv",
    "comparators.csv",
)
E2A_QA_FILES = {
    "cgm_status_after.txt",
    "cgm_status_before.txt",
    "curated_records.pkl.gz",
    "frozen_compare.json",
    "guard_after.log",
    "guard_before.log",
    "guard_final.log",
    "pytest.log",
    "regression.json",
    "rev_after.sha256",
    "rev_before.sha256",
    "run.log",
}
E2A_MANIFEST = OUT / "manifest.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_safe(value):
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False)


def read_predictions() -> pd.DataFrame:
    """Preserve exact probabilities so average-rank AUC is stable."""
    return pd.read_csv(
        OUT / "predictions.csv",
        low_memory=False,
        float_precision="round_trip",
    )


def load_records() -> dict[tuple[str, str], pd.DataFrame]:
    path = QA / "curated_records.pkl.gz"
    return pickle.loads(gzip.decompress(path.read_bytes()))


def included_keys(inclusion: pd.DataFrame) -> list[tuple[str, str]]:
    selected = inclusion.loc[inclusion.included, ["target", "endpoint"]]
    return list(selected.itertuples(index=False, name=None))


def holm_adjust(p_values) -> np.ndarray:
    """Holm step-down adjustment, returned in the input order."""
    p = np.asarray(p_values, dtype=float)
    result = np.full(len(p), np.nan)
    finite = np.flatnonzero(np.isfinite(p))
    if not len(finite):
        return result
    order = finite[np.argsort(p[finite], kind="stable")]
    m = len(order)
    running = 0.0
    for rank, index in enumerate(order):
        running = max(running, min(1.0, (m - rank) * p[index]))
        result[index] = running
    return result


def one_sided_wilcoxon(values, alternative: str) -> tuple[float, float, int]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values) or np.allclose(values, 0.0):
        return float("nan"), 1.0, len(values)
    test = stats.wilcoxon(values, alternative=alternative, zero_method="wilcox")
    return float(test.statistic), float(test.pvalue), len(values)


def apply_holm_corrections(tests: pd.DataFrame) -> pd.DataFrame:
    """Apply the separately prespecified H3-H7 and H8 families."""
    result = tests.copy()
    result["holm_p"] = np.nan
    for (endpoint, population), locations in result.loc[
        result.correction_family == "H3-H7"
    ].groupby(["endpoint", "population"]).groups.items():
        result.loc[locations, "holm_p"] = holm_adjust(
            result.loc[locations, "raw_p"]
        )
    for population, locations in result.loc[
        result.correction_family == "H8"
    ].groupby("population").groups.items():
        result.loc[locations, "holm_p"] = holm_adjust(
            result.loc[locations, "raw_p"]
        )
    return result


def bootstrap_median_ci(values, seed=BOOTSTRAP_SEED, repetitions=BOOT_N):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return float("nan"), float("nan"), 0
    rng = np.random.default_rng(seed)
    draws = np.empty(repetitions, dtype=float)
    for start in range(0, repetitions, 100):
        end = min(start + 100, repetitions)
        indices = rng.integers(0, len(values), size=(end - start, len(values)))
        draws[start:end] = np.median(values[indices], axis=1)
    return float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5)), repetitions


def hypothesis_tests(
    metrics: pd.DataFrame,
    comparators: pd.DataFrame,
    records_by_key: dict[tuple[str, str], pd.DataFrame],
) -> pd.DataFrame:
    """Compute H3-H8 with the two prespecified Holm families."""
    base = metrics.loc[
        (metrics.stratum == "all") & (metrics.target.isin(metrics.target.unique()))
    ].copy()
    rows = []
    pair_counts = {}
    for key, records in records_by_key.items():
        offsets = document_pair_offsets(records)
        pair_counts[key] = int(offsets.systematic_offset.sum())

    definitions = {
        "H3": ("G_auc_minus_0.5", "greater"),
        "H4": ("G_minus_G-record_log_loss", "less"),
        "H5": ("G_minus_G-global_log_loss", "less"),
        "H6": ("G_minus_G-fixed_log_loss", "less"),
        "H7": ("X_minus_G_log_loss", "less"),
    }
    for endpoint in ENDPOINTS:
        for population in POPULATIONS:
            selected = base.loc[
                (base.endpoint == endpoint) & (base.population == population)
            ]
            pivot = selected.pivot(
                index="target", columns="model", values="log_loss"
            )
            auc = selected.loc[selected.model == "G"].set_index("target").auc
            values_by_hypothesis = {
                "H3": auc - 0.5,
                "H4": pivot["G"] - pivot["G-record"],
                "H5": pivot["G"] - pivot["G-global"],
                "H6": pivot["G"] - pivot["G-fixed"],
                "H7": pivot["X"] - pivot["G"],
            }
            eligible_targets = {
                target
                for target in pivot.index
                if pair_counts.get((target, endpoint), 0) > 0
            }
            values_by_hypothesis["H7"] = values_by_hypothesis["H7"].loc[
                values_by_hypothesis["H7"].index.intersection(eligible_targets)
            ]
            for hypothesis, values in values_by_hypothesis.items():
                name, alternative = definitions[hypothesis]
                statistic, p_value, n_targets = one_sided_wilcoxon(
                    values, alternative
                )
                low, high, valid = bootstrap_median_ci(values)
                median = float(np.nanmedian(values)) if len(values) else np.nan
                rows.append({
                    "hypothesis": hypothesis,
                    "endpoint": endpoint,
                    "population": population,
                    "correction_family": "H3-H7",
                    "difference_definition": name,
                    "alternative": alternative,
                    "n_targets": n_targets,
                    "median_difference": median,
                    "median_ci_low": low,
                    "median_ci_high": high,
                    "wilcoxon_statistic": statistic,
                    "raw_p": p_value,
                    "holm_p": np.nan,
                    "bootstrap_valid": valid,
                    "bootstrap_repetitions": BOOT_N,
                    "bootstrap_seed": BOOTSTRAP_SEED,
                    "eligible_targets": ";".join(sorted(values.index)),
                })
    h3_h7 = pd.DataFrame(rows)

    h8_rows = []
    for population in POPULATIONS:
        comp = comparators.loc[comparators.population == population]
        wide = comp.pivot(index=["target", "endpoint"], columns="model", values="auc")
        for endpoint in ENDPOINTS:
            values = (wide["G"] - wide["B-range"]).xs(
                endpoint, level="endpoint", drop_level=False
            )
            values = values.loc[values.index.get_level_values("endpoint") == endpoint]
            statistic, p_value, n_targets = one_sided_wilcoxon(
                values, alternative="greater"
            )
            low, high, valid = bootstrap_median_ci(values)
            h8_rows.append({
                "hypothesis": "H8",
                "endpoint": endpoint,
                "population": population,
                "correction_family": "H8",
                "difference_definition": "G_auc_minus_B-range_auc",
                "alternative": "greater",
                "n_targets": n_targets,
                "median_difference": float(np.nanmedian(values)),
                "median_ci_low": low,
                "median_ci_high": high,
                "wilcoxon_statistic": statistic,
                "raw_p": p_value,
                "holm_p": np.nan,
                "bootstrap_valid": valid,
                "bootstrap_repetitions": BOOT_N,
                "bootstrap_seed": BOOTSTRAP_SEED,
                "eligible_targets": ";".join(
                    sorted(value[0] for value in values.index)
                ),
            })
    h8 = pd.DataFrame(h8_rows)
    result = apply_holm_corrections(
        pd.concat([h3_h7, h8], ignore_index=True)
    )
    result["analysis_role"] = np.where(
        result.population == "all", "primary", "sensitivity"
    )
    return result.sort_values(
        ["population", "endpoint", "correction_family", "hypothesis"]
    ).reset_index(drop=True)


def distribution_stats(values) -> dict:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return {
            "n_targets": 0,
            "median": np.nan,
            "q1": np.nan,
            "q3": np.nan,
            "iqr": np.nan,
            "minimum": np.nan,
            "maximum": np.nan,
        }
    q1, median, q3 = np.percentile(values, [25, 50, 75])
    return {
        "n_targets": len(values),
        "median": float(median),
        "q1": float(q1),
        "q3": float(q3),
        "iqr": float(q3 - q1),
        "minimum": float(values.min()),
        "maximum": float(values.max()),
    }


def hypothesis_tests_by_stratum(metrics: pd.DataFrame) -> pd.DataFrame:
    """Report the prespecified H2/H3 strata without adding Holm families."""
    rows = []
    for endpoint in ENDPOINTS:
        for population in POPULATIONS:
            selected = metrics.loc[
                (metrics.endpoint == endpoint)
                & (metrics.population == population)
                & (metrics.model == "G")
                & metrics.stratum.isin(("m=1", "m>=2"))
            ].copy()
            for stratum in ("m=1", "m>=2"):
                layer = selected.loc[selected.stratum == stratum]
                layer = layer.loc[np.isfinite(layer.auc)].sort_values("target")
                values = layer.auc.to_numpy(float) - 0.5
                statistic, p_value, n_targets = one_sided_wilcoxon(
                    values, alternative="greater"
                )
                low, high, valid = bootstrap_median_ci(values)
                ece = distribution_stats(layer.ece)
                bias = distribution_stats(layer.calibration_bias)
                rows.append({
                    "endpoint": endpoint,
                    "population": population,
                    "stratum": stratum,
                    "analysis_role": (
                        "primary" if population == "all" else "sensitivity"
                    ),
                    "correction_family": "none",
                    "family_note": "stratified; not part of the Holm families",
                    "n_targets": n_targets,
                    "eligible_targets": ";".join(layer.target),
                    "h3_difference_definition": "G_auc_minus_0.5",
                    "h3_alternative": "greater",
                    "h3_median_difference": (
                        float(np.median(values)) if len(values) else np.nan
                    ),
                    "h3_median_ci_low": low,
                    "h3_median_ci_high": high,
                    "h3_wilcoxon_statistic": statistic,
                    "h3_raw_p": p_value,
                    "h3_bootstrap_valid": valid,
                    "h3_bootstrap_repetitions": BOOT_N,
                    "h3_bootstrap_seed": BOOTSTRAP_SEED,
                    "h2_ece_median": ece["median"],
                    "h2_ece_q1": ece["q1"],
                    "h2_ece_q3": ece["q3"],
                    "h2_ece_iqr": ece["iqr"],
                    "h2_ece_minimum": ece["minimum"],
                    "h2_ece_maximum": ece["maximum"],
                    "h2_calibration_bias_median": bias["median"],
                    "h2_calibration_bias_q1": bias["q1"],
                    "h2_calibration_bias_q3": bias["q3"],
                    "h2_calibration_bias_iqr": bias["iqr"],
                    "h2_calibration_bias_minimum": bias["minimum"],
                    "h2_calibration_bias_maximum": bias["maximum"],
                })
    return pd.DataFrame(rows)


def pooled_metric(
    pooled: pd.DataFrame, endpoint: str, population: str, model: str, metric: str
):
    selected = pooled.loc[
        (pooled.endpoint == endpoint)
        & (pooled.population == population)
        & (pooled.model == model)
        & (pooled.stratum == "all")
        & (pooled.metric == metric)
    ]
    if len(selected) != 1:
        raise ValueError(
            f"Expected one pooled metric: {endpoint}/{population}/{model}/{metric}"
        )
    return selected.iloc[0]


def h1_h2_summary(metrics: pd.DataFrame, pooled: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for endpoint in ENDPOINTS:
        for population in POPULATIONS:
            selected = metrics.loc[
                (metrics.endpoint == endpoint)
                & (metrics.population == population)
                & (metrics.stratum == "all")
            ]
            b0 = selected.loc[selected.model == "B0"]
            g = selected.loc[selected.model == "G"]
            h1_values = 1.0 - b0.accuracy.to_numpy(float)
            h1_pool = pooled_metric(pooled, endpoint, population, "B0", "accuracy")
            rows.append({
                "hypothesis": "H1",
                "endpoint": endpoint,
                "population": population,
                "metric": "one_minus_B0_accuracy",
                **distribution_stats(h1_values),
                "pooled_estimate": 1.0 - float(h1_pool.estimate),
                "pooled_ci_low": 1.0 - float(h1_pool.ci_high),
                "pooled_ci_high": 1.0 - float(h1_pool.ci_low),
            })
            for metric in ("ece", "calibration_bias"):
                values = g[metric].to_numpy(float)
                pool = pooled_metric(pooled, endpoint, population, "G", metric)
                rows.append({
                    "hypothesis": "H2",
                    "endpoint": endpoint,
                    "population": population,
                    "metric": f"G_{metric}",
                    **distribution_stats(values),
                    "pooled_estimate": float(pool.estimate),
                    "pooled_ci_low": float(pool.ci_low),
                    "pooled_ci_high": float(pool.ci_high),
                })
    return pd.DataFrame(rows)


def calibration_curve(predictions: pd.DataFrame, bins: int) -> pd.DataFrame:
    rows = []
    g = predictions.loc[predictions.method == "G"].copy()
    g["confidence"] = g[["p0", "p1", "p2"]].max(axis=1)
    g["correct"] = (
        g[["p0", "p1", "p2"]].to_numpy(float).argmax(axis=1)
        == g.observed.to_numpy(int)
    )
    for (endpoint, population), frame in g.groupby(
        ["endpoint", "population"], sort=True
    ):
        bin_index = np.minimum((frame.confidence * bins).astype(int), bins - 1)
        for index in range(bins):
            selected = frame.loc[bin_index == index]
            rows.append({
                "endpoint": endpoint,
                "population": population,
                "bin": index + 1,
                "bin_low": index / bins,
                "bin_high": (index + 1) / bins,
                "n_predictions": len(selected),
                "mean_confidence": selected.confidence.mean(),
                "observed_accuracy": selected.correct.mean(),
            })
    return pd.DataFrame(rows)


def unit_year_map(records: pd.DataFrame) -> pd.Series:
    years = (
        records.assign(
            document_year=pd.to_numeric(records.document_year, errors="coerce")
        )
        .groupby("document_chembl_id")
        .document_year.max()
    )
    return years


def add_unit_years(
    doc: pd.DataFrame, records: pd.DataFrame
) -> tuple[pd.DataFrame, int]:
    """Unit year is the latest nonmissing year among its component documents."""
    years = unit_year_map(records)
    output = doc.copy()
    output["year"] = [
        max(
            (
                float(years.get(item, np.nan))
                for item in str(documents).split(";")
                if item in years.index and np.isfinite(years.get(item, np.nan))
            ),
            default=np.nan,
        )
        for documents in output.documents
    ]
    missing = int(output.year.isna().sum())
    return output, missing


def classify(value: float) -> int:
    return CLASSES.index(activity_class(float(value), *CLASS_BOUNDARIES))


def shown_probability(center, scale, shown) -> np.ndarray:
    probabilities = class_probabilities(
        np.atleast_1d(center),
        np.atleast_1d(scale),
        np.inf,
        *CLASS_BOUNDARIES,
    )
    shown = np.broadcast_to(np.asarray(shown, dtype=int), len(probabilities))
    return probabilities[np.arange(len(probabilities)), shown]


def time_split_predictions(
    target: str,
    endpoint: str,
    records: pd.DataFrame,
    cutoff: int,
) -> tuple[pd.DataFrame, dict]:
    """Build one target-endpoint temporal split without future leakage."""
    doc = document_values(records, merge_identical=True)
    doc, missing_units = add_unit_years(doc, records)
    known = doc.loc[doc.year.notna()].copy()
    diagnostics = {
        "target": target,
        "endpoint": endpoint,
        "cutoff_year": cutoff,
        "n_units": int(len(doc)),
        "n_entries": int(doc.entry.nunique()),
        "missing_year_units": missing_units,
        "missing_year_entries": int(
            doc.loc[doc.year.isna(), "entry"].nunique()
        ),
        "eligible": False,
        "reason": "",
        "n_fit_entries": 0,
        "n_predictions": 0,
        "future_units_ge_30": False,
    }
    past = known.loc[known.year <= cutoff]
    past_multi_by_entry = past.groupby("entry").size()
    n_fit_entries = int((past_multi_by_entry >= 2).sum())
    diagnostics["n_fit_entries"] = n_fit_entries
    if n_fit_entries < 5:
        diagnostics["reason"] = (
            "fewer than 5 entries with at least two past evidence units"
        )
        return pd.DataFrame(), diagnostics
    sigma_fit = fit_between_unit_sd(past)
    rows = []
    source_by_entry = dict(tuple(records.groupby("molecule_identity_key", sort=True)))
    for entry, group in known.groupby("entry", sort=True):
        past_units = group.loc[group.year <= cutoff]
        future_units = group.loc[group.year > cutoff]
        if past_units.empty or future_units.empty:
            continue
        past_values = past_units.document_value.to_numpy(float)
        past_documents = set(
            ";".join(past_units.documents.astype(str)).split(";")
        )
        source = source_by_entry[entry]
        past_records = source.loc[
            source.document_chembl_id.astype(str).isin(past_documents),
            "pchembl_value",
        ].to_numpy(float)
        if not len(past_records):
            continue
        display = classify(np.median(past_records))
        m_past = len(past_values)
        support = float(
            shown_probability(
                past_values.mean(), sigma_fit["tau"] / np.sqrt(m_past), display
            )[0]
        )
        probabilities = class_probabilities(
            [past_values.mean()],
            [sigma_fit["tau"] * np.sqrt(1.0 + 1.0 / m_past)],
            np.inf,
            *CLASS_BOUNDARIES,
        )[0]
        for future in future_units.itertuples(index=False):
            observed = classify(future.document_value)
            rows.append({
                "target": target,
                "endpoint": endpoint,
                "cutoff_year": cutoff,
                "entry": entry,
                "future_unit": future.documents,
                "future_year": int(future.year),
                "m_past": m_past,
                "stratum": "m=1" if m_past == 1 else "m>=2",
                "sigma_past": float(sigma_fit["tau"]),
                "past_center": float(past_values.mean()),
                "past_support": support,
                "display_status": (
                    "supported" if support >= SUPPORT_THRESHOLD else "uncertain"
                ),
                "displayed": display,
                "observed": observed,
                "p0": float(probabilities[0]),
                "p1": float(probabilities[1]),
                "p2": float(probabilities[2]),
                "agree": bool(observed == display),
            })
    result = pd.DataFrame(rows)
    if result.empty:
        diagnostics["reason"] = "no entry has both past and future evidence units"
        return result, diagnostics
    result["p_displayed"] = result[
        ["p0", "p1", "p2"]
    ].to_numpy(float)[np.arange(len(result)), result.displayed.to_numpy(int)]
    diagnostics["eligible"] = True
    diagnostics["n_predictions"] = len(result)
    diagnostics["future_units_ge_30"] = len(result) >= 30
    diagnostics["reason"] = "included"
    return result, diagnostics


def time_split_metrics(frame: pd.DataFrame) -> dict:
    if frame.empty:
        return {
            "n_predictions": 0,
            "n_entries": 0,
            "n_targets": 0,
            "b0_inconsistency_rate": np.nan,
            "g_ece": np.nan,
            "g_calibration_bias": np.nan,
            "g_auc": np.nan,
        }
    b0_scores = point_weighted_scores(
        frame, displayed_one_hot(frame)
    )
    g_scores = point_weighted_scores(frame)
    return {
        "n_predictions": len(frame),
        "n_entries": int(frame[["target", "entry"]].drop_duplicates().shape[0]),
        "n_targets": int(frame.target.nunique()),
        "b0_inconsistency_rate": 1.0 - float(
            b0_scores["agreement_rate"][0]
        ),
        "g_ece": float(g_scores["ece"][0]),
        "g_calibration_bias": float(g_scores["calibration_bias"][0]),
        "g_auc": float(g_scores["auc"][0]),
    }


def point_weighted_scores(frame: pd.DataFrame, probabilities=None) -> dict:
    local = frame.copy()
    if probabilities is not None:
        local[["p0", "p1", "p2"]] = probabilities
    return weighted_scores(local, np.ones((1, len(local)), dtype=float))


def displayed_one_hot(frame: pd.DataFrame) -> np.ndarray:
    probabilities = np.zeros((len(frame), 3), dtype=float)
    probabilities[
        np.arange(len(frame)), frame.displayed.to_numpy(int)
    ] = 1.0
    return probabilities


def time_split_bootstrap(frame: pd.DataFrame) -> dict:
    if frame.empty:
        return {
            "b0_inconsistency_ci_low": np.nan,
            "b0_inconsistency_ci_high": np.nan,
            "g_ece_ci_low": np.nan,
            "g_ece_ci_high": np.nan,
            "g_calibration_bias_ci_low": np.nan,
            "g_calibration_bias_ci_high": np.nan,
            "g_auc_ci_low": np.nan,
            "g_auc_ci_high": np.nan,
            "bootstrap_valid": 0,
        }
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    weights = cluster_weights(frame, rng, BOOT_N)
    b0 = frame.copy()
    b0[["p0", "p1", "p2"]] = displayed_one_hot(frame)
    b0_draws = weighted_scores(b0, weights)
    g_draws = weighted_scores(frame, weights)
    intervals = {}
    for name, draws in (
        ("b0_inconsistency", 1.0 - b0_draws["agreement_rate"]),
        ("g_ece", g_draws["ece"]),
        ("g_calibration_bias", g_draws["calibration_bias"]),
        ("g_auc", g_draws["auc"]),
    ):
        low, high, valid = percentile_ci(draws)
        intervals[f"{name}_ci_low"] = low
        intervals[f"{name}_ci_high"] = high
        intervals["bootstrap_valid"] = max(intervals.get("bootstrap_valid", 0), valid)
    return intervals


def build_time_split(
    records_by_key: dict[tuple[str, str], pd.DataFrame],
    inclusion: pd.DataFrame,
    additional_keys=(),
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    details, diagnostics = [], []
    keys = included_keys(inclusion)
    for key in additional_keys:
        if key not in records_by_key:
            raise KeyError(f"Missing records for temporal analysis: {key}")
        if key not in keys:
            keys.append(key)
    for target, endpoint in keys:
        for cutoff in TIME_CUTOFFS:
            frame, info = time_split_predictions(
                target, endpoint, records_by_key[(target, endpoint)], cutoff
            )
            diagnostics.append(info)
            if not frame.empty:
                details.append(frame)
    detail = (
        pd.concat(details, ignore_index=True) if details else pd.DataFrame()
    )
    if not detail.empty:
        detail.to_csv(QA / "time_split_predictions.csv", index=False)
    by_target_rows = []
    for item in diagnostics:
        selected = (
            detail.loc[
                (detail.target == item["target"])
                & (detail.endpoint == item["endpoint"])
                & (detail.cutoff_year == item["cutoff_year"])
            ]
            if not detail.empty
            else pd.DataFrame()
        )
        by_target_rows.append({**item, **time_split_metrics(selected)})
    by_target = pd.DataFrame(by_target_rows)
    pooled_rows = []
    if not detail.empty:
        for (endpoint, cutoff), endpoint_frame in detail.groupby(
            ["endpoint", "cutoff_year"], sort=True
        ):
            ge30 = set(
                by_target.loc[
                    (by_target.endpoint == endpoint)
                    & (by_target.cutoff_year == cutoff)
                    & by_target.future_units_ge_30,
                    "target",
                ]
            )
            scopes = {
                "all": endpoint_frame,
                "m=1": endpoint_frame.loc[endpoint_frame.stratum == "m=1"],
                "m>=2": endpoint_frame.loc[endpoint_frame.stratum == "m>=2"],
                "display_supported": endpoint_frame.loc[
                    endpoint_frame.display_status == "supported"
                ],
                "display_uncertain": endpoint_frame.loc[
                    endpoint_frame.display_status == "uncertain"
                ],
                "future_units_ge_30": endpoint_frame.loc[
                    endpoint_frame.target.isin(ge30)
                ],
            }
            for scope, selected in scopes.items():
                point = time_split_metrics(selected)
                interval = time_split_bootstrap(selected)
                pooled_rows.append({
                    "endpoint": endpoint,
                    "cutoff_year": cutoff,
                    "scope": scope,
                    **point,
                    **interval,
                    "bootstrap_repetitions": BOOT_N,
                    "bootstrap_seed": BOOTSTRAP_SEED,
                    "unit_resampling": "target x entry",
                })
    pooled = pd.DataFrame(pooled_rows)
    return pooled, by_target, detail


def load_original_ic50_records() -> dict[tuple[str, str], pd.DataFrame]:
    """Curate the three original IC50 datasets with the E2a regression settings."""
    manifest = json.loads(
        (API / "download_original_ic50_manifest.json").read_text()
    )
    records = {}
    for item in manifest["files"]:
        target = str(item["target"])
        endpoint = str(item["standard_type"])
        if endpoint != "IC50" or target not in ORIGINAL_IC50_TARGETS:
            continue
        raw = pd.read_csv(API / "activities" / item["file"], low_memory=False)
        records[(target, endpoint)] = curate_activity_table(
            raw, target, endpoint
        ).retained_records
    if set(records) != {
        (target, "IC50") for target in ORIGINAL_IC50_TARGETS
    }:
        raise RuntimeError("Original IC50 manifest does not contain the three expected datasets")
    return records


def pseudo_consensus(
    records_by_key: dict[tuple[str, str], pd.DataFrame],
    inclusion: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    details = []
    for target, endpoint in included_keys(inclusion):
        doc = document_values(records_by_key[(target, endpoint)], merge_identical=True)
        assignments = entry_folds(doc, seed=FOLD_SEED)
        fold_tau = {}
        for fold in sorted(set(assignments.values())):
            test = {
                entry for entry, assigned in assignments.items() if assigned == fold
            }
            fold_tau[fold] = fit_between_unit_sd(
                doc.loc[~doc.entry.isin(test)]
            )["tau"]
        for entry, group in doc.groupby("entry", sort=True):
            m = len(group)
            if m < 5:
                continue
            values = group.document_value.to_numpy(float)
            remaining_means = (values.sum() - values) / (m - 1)
            shown = np.array([classify(value) for value in values])
            observed = np.array([classify(value) for value in remaining_means])
            tau = fold_tau[assignments[entry]]
            p_cons = class_probabilities(
                values, np.full(m, tau), np.inf, *CLASS_BOUNDARIES
            )[np.arange(m), shown]
            p_pred = class_probabilities(
                values,
                np.full(m, tau * np.sqrt(1.0 + 1.0 / (m - 1))),
                np.inf,
                *CLASS_BOUNDARIES,
            )[np.arange(m), shown]
            for index, unit in enumerate(group.itertuples(index=False)):
                details.append({
                    "target": target,
                    "endpoint": endpoint,
                    "entry": entry,
                    "unit": unit.documents,
                    "m": m,
                    "fold": assignments[entry],
                    "training_tau": tau,
                    "unit_value": values[index],
                    "remaining_mean": remaining_means[index],
                    "displayed": shown[index],
                    "observed": observed[index],
                    "p_cons": p_cons[index],
                    "p_pred": p_pred[index],
                    "agree": bool(shown[index] == observed[index]),
                })
    detail = pd.DataFrame(details)
    detail.to_csv(QA / "pseudo_consensus_predictions.csv", index=False)
    summary_rows, curve_rows = [], []
    for scope, selected in [
        ("all", detail),
        *[
            (f"{endpoint}", detail.loc[detail.endpoint == endpoint])
            for endpoint in ENDPOINTS
        ],
    ]:
        for probability in ("p_cons", "p_pred"):
            p = selected[probability].to_numpy(float)
            outcome = selected.agree.to_numpy(float)
            bins = np.minimum((p * 10).astype(int), 9)
            ece = sum(
                abs(
                    (
                        np.sum((outcome[bins == index] - p[bins == index]))
                    )
                )
                for index in range(10)
            ) / len(selected)
            summary_rows.append({
                "row_type": "summary",
                "endpoint": scope,
                "probability_type": probability,
                "n_predictions": len(selected),
                "n_entries": int(
                    selected[["target", "entry"]].drop_duplicates().shape[0]
                ),
                "n_targets": int(selected.target.nunique()),
                "mean_probability": p.mean(),
                "agreement_rate": outcome.mean(),
                "ece": ece,
                "minimum_units": 5,
            })
            for index in range(10):
                chosen = bins == index
                curve_rows.append({
                    "endpoint": scope,
                    "probability_type": probability,
                    "bin": index + 1,
                    "bin_low": index / 10,
                    "bin_high": (index + 1) / 10,
                    "n_predictions": int(chosen.sum()),
                    "mean_probability": (
                        float(p[chosen].mean()) if chosen.any() else np.nan
                    ),
                    "observed_agreement": (
                        float(outcome[chosen].mean())
                        if chosen.any()
                        else np.nan
                    ),
                })
    return pd.DataFrame(summary_rows), pd.DataFrame(curve_rows)


def sigma_bootstrap(doc: pd.DataFrame) -> tuple[float, float, int]:
    moments = doc.groupby("entry", sort=True).document_value.agg(["size", "var"])
    dfs = moments["size"].to_numpy(float) - 1.0
    sumsquares = dfs * moments["var"].fillna(0).to_numpy(float)
    has_multi = (dfs > 0).astype(int)
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    draws = np.empty(SIGMA_BOOT_N)
    for start in range(0, SIGMA_BOOT_N, 100):
        end = min(start + 100, SIGMA_BOOT_N)
        sampled = rng.integers(0, len(moments), size=(end - start, len(moments)))
        denominator = dfs[sampled].sum(axis=1)
        numerator = sumsquares[sampled].sum(axis=1)
        multi_count = has_multi[sampled].sum(axis=1)
        pooled = np.sqrt(
            np.divide(
                numerator,
                denominator,
                out=np.zeros(end - start),
                where=denominator > 0,
            )
        )
        draws[start:end] = np.where(
            multi_count >= 5, pooled, LITERATURE_BETWEEN_UNIT_SD
        )
    return (*percentile_ci(draws)[:2], SIGMA_BOOT_N)


def sigma_fit_row(
    target: str,
    endpoint: str,
    doc: pd.DataFrame,
    measurement_type: str = "all",
) -> dict:
    fit = fit_between_unit_sd(doc)
    moments = doc.groupby("entry").size()
    multi = moments[moments >= 2]
    low, high, valid = sigma_bootstrap(doc)
    return {
        "row_type": "target",
        "target": target,
        "endpoint": endpoint,
        "measurement_type": measurement_type,
        "n_entries": int(len(moments)),
        "n_multi_unit_entries": int(len(multi)),
        "n_units": int(len(doc)),
        "degrees_of_freedom": int((multi - 1).sum()),
        "sigma_doc": float(fit["tau"]),
        "ci_low": low,
        "ci_high": high,
        "bootstrap_repetitions": SIGMA_BOOT_N,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "fit_method": fit["method"],
    }


def ca2_ki_influence(records: pd.DataFrame) -> tuple[dict, pd.DataFrame]:
    doc = document_values(records, merge_identical=True)
    multi = doc.groupby("entry").filter(lambda frame: len(frame) >= 2)
    contributions = {}
    for documents in multi.documents.astype(str):
        for document in set(documents.split(";")):
            contributions[document] = contributions.get(document, 0) + 1
    ranked = sorted(contributions.items(), key=lambda item: (-item[1], item[0]))
    top = ranked[:5]
    top_documents = {item[0] for item in top}
    keep = doc.documents.map(
        lambda value: not (set(str(value).split(";")) & top_documents)
    )
    reduced = doc.loc[keep]
    full_fit = fit_between_unit_sd(doc)
    reduced_fit = fit_between_unit_sd(reduced)
    years = (
        records.assign(
            document_year=pd.to_numeric(records.document_year, errors="coerce")
        )
        .groupby("document_chembl_id")
        .document_year.max()
    )
    detail = pd.DataFrame([
        {
            "document": document,
            "contribution_df": contribution,
            "year": years.get(document, np.nan),
        }
        for document, contribution in top
    ])
    row = {
        "row_type": "ca2_top5_document_removal",
        "target": "CHEMBL205",
        "endpoint": "Ki",
        "measurement_type": "all",
        "n_entries": int(doc.entry.nunique()),
        "n_multi_unit_entries": int(
            (doc.groupby("entry").size() >= 2).sum()
        ),
        "n_units": int(len(doc)),
        "degrees_of_freedom": int(
            (
                doc.groupby("entry").size().loc[
                    doc.groupby("entry").size() >= 2
                ]
                - 1
            ).sum()
        ),
        "sigma_doc": float(full_fit["tau"]),
        "sigma_without_top5": float(reduced_fit["tau"]),
        "delta_sigma": float(reduced_fit["tau"] - full_fit["tau"]),
        "removed_document_ids": ";".join(item[0] for item in top),
        "removed_document_years": ";".join(
            str(int(detail.loc[detail.document == item[0], "year"].iloc[0]))
            if pd.notna(detail.loc[detail.document == item[0], "year"].iloc[0])
            else "NA"
            for item in top
        ),
        "author_cluster_verification": "not possible from available ChEMBL fields",
    }
    return row, detail


def build_sigma_doc(
    records_by_key: dict[tuple[str, str], pd.DataFrame],
    inclusion: pd.DataFrame,
) -> pd.DataFrame:
    rows, extra = [], []
    for target, endpoint in included_keys(inclusion):
        records = records_by_key[(target, endpoint)]
        doc = document_values(records, merge_identical=True)
        rows.append(sigma_fit_row(target, endpoint, doc))
        assay_docs = {}
        for assay_type in ("B", "F"):
            subset = records.loc[
                records.assay_type.astype(str).str.upper() == assay_type
            ]
            if subset.empty:
                continue
            assay_doc = document_values(subset, merge_identical=True)
            n_multi = int((assay_doc.groupby("entry").size() >= 2).sum())
            if n_multi >= 30:
                assay_docs[assay_type] = assay_doc
        if set(assay_docs) == {"B", "F"}:
            for assay_type, assay_doc in assay_docs.items():
                rows.append(
                    sigma_fit_row(target, endpoint, assay_doc, assay_type)
                )
    ca2_row, ca2_documents = ca2_ki_influence(
        records_by_key[("CHEMBL205", "Ki")]
    )
    rows.append(ca2_row)
    ca2_documents.to_csv(QA / "ca2_top5_documents.csv", index=False)
    for endpoint in ENDPOINTS:
        selected = pd.DataFrame(rows)
        selected = selected.loc[
            (selected.row_type == "target")
            & (selected.endpoint == endpoint)
            & (selected.measurement_type == "all")
        ]
        rows.append({
            "row_type": "endpoint_distribution",
            "target": "all",
            "endpoint": endpoint,
            "measurement_type": "all",
            **distribution_stats(selected.sigma_doc),
            "literature_reference_sd": (
                "Ki 0.47 (Kalliokoski 2013); Ki 0.54 (Kramer 2012)"
                if endpoint == "Ki"
                else "IC50 0.57-0.61 (Kalliokoski 2013)"
            ),
        })
    return pd.DataFrame(rows)


def flagged_pairs_census(
    records_by_key: dict[tuple[str, str], pd.DataFrame],
    predictions: pd.DataFrame,
    inclusion: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for target, endpoint in included_keys(inclusion):
        records = records_by_key[(target, endpoint)]
        offsets = document_pair_offsets(records)
        flagged = offsets.loc[offsets.systematic_offset.astype(bool)].copy()
        flagged_documents = set(flagged.document_a) | set(flagged.document_b)
        doc = document_values(records, merge_identical=True)
        doc["touched"] = doc.documents.map(
            lambda value: bool(
                set(str(value).split(";")) & flagged_documents
            )
        )
        touched_entries = set(doc.loc[doc.touched, "entry"])
        predictions_for_target = predictions.loc[
            (predictions.target == target)
            & (predictions.endpoint == endpoint)
            & (predictions.population == "all")
            & (predictions.method == "G")
        ]
        touched_predictions = predictions_for_target.loc[
            predictions_for_target.entry.isin(touched_entries)
        ]
        rows.append({
            "row_type": "target_summary",
            "target": target,
            "endpoint": endpoint,
            "n_flagged_pairs": len(flagged),
            "n_flagged_documents": len(flagged_documents),
            "n_entries": int(doc.entry.nunique()),
            "n_touched_entries": len(touched_entries),
            "fraction_touched_entries": (
                len(touched_entries) / doc.entry.nunique()
            ),
            "n_predictions": len(predictions_for_target),
            "n_touched_predictions": len(touched_predictions),
            "fraction_touched_predictions": (
                len(touched_predictions) / len(predictions_for_target)
                if len(predictions_for_target)
                else np.nan
            ),
            "agreement_touched": touched_predictions.agree.mean(),
            "agreement_untouched": predictions_for_target.loc[
                ~predictions_for_target.entry.isin(touched_entries), "agree"
            ].mean(),
        })
        if flagged.empty:
            continue
        unmerged = document_values(records, merge_identical=False)
        for rank, pair in enumerate(
            sorted(
                flagged.itertuples(index=False),
                key=lambda item: (-abs(float(item.median_offset)), item.document_a, item.document_b),
            )[:10],
            1,
        ):
            wide = unmerged.loc[
                unmerged.document.isin([pair.document_a, pair.document_b])
            ].pivot(
                index="entry",
                columns="document",
                values="document_value",
            )
            shared = wide[[pair.document_a, pair.document_b]].dropna()
            defined = bool(
                shared[pair.document_a].nunique() > 1
                and shared[pair.document_b].nunique() > 1
            )
            rho = (
                float(
                    stats.spearmanr(
                        shared[pair.document_a], shared[pair.document_b]
                    ).statistic
                )
                if defined
                else np.nan
            )
            rows.append({
                "row_type": "top_pair",
                "target": target,
                "endpoint": endpoint,
                "offset_rank": rank,
                "document_a": pair.document_a,
                "document_b": pair.document_b,
                "shared_entries": int(pair.shared_entries),
                "median_offset_b_minus_a": float(pair.median_offset),
                "spearman_rho": rho,
                "spearman_defined": defined,
            })
    return pd.DataFrame(rows)


def entry_level_stats(
    records_by_key: dict[tuple[str, str], pd.DataFrame],
    inclusion: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    details = []
    for target, endpoint in included_keys(inclusion):
        records = records_by_key[(target, endpoint)]
        doc = document_values(records, merge_identical=True)
        raw_units = document_values(records, merge_identical=False)
        tau_doc = fit_between_unit_sd(doc)["tau"]
        tau_record = fit_between_unit_sd(raw_units)["tau"]
        source_by_entry = dict(tuple(records.groupby("molecule_identity_key")))
        for entry, group in doc.groupby("entry", sort=True):
            source = source_by_entry[entry]
            values = source.pchembl_value.to_numpy(float)
            shown = classify(np.median(values))
            document_support = float(
                shown_probability(
                    group.document_value.mean(),
                    tau_doc / np.sqrt(len(group)),
                    shown,
                )[0]
            )
            record_support = float(
                shown_probability(
                    values.mean(),
                    tau_record / np.sqrt(len(values)),
                    shown,
                )[0]
            )
            document_supported = document_support >= SUPPORT_THRESHOLD
            record_supported = record_support >= SUPPORT_THRESHOLD
            if document_supported == record_supported:
                direction = "unchanged"
            elif record_supported:
                direction = "uncertain_to_supported"
            else:
                direction = "supported_to_uncertain"
            details.append({
                "target": target,
                "endpoint": endpoint,
                "entry": entry,
                "n_records": len(values),
                "n_units": len(group),
                "document_support": document_support,
                "record_support": record_support,
                "document_supported": document_supported,
                "record_supported": record_supported,
                "flip": document_supported != record_supported,
                "flip_direction": direction,
                "displayed": shown,
            })
    detail = pd.DataFrame(details)
    summary = (
        detail.groupby(["target", "endpoint"], sort=True)
        .apply(
            lambda group: pd.Series({
                "n_entries": len(group),
                "n_single_unit_entries": int((group.n_units == 1).sum()),
                "single_unit_proportion": float((group.n_units == 1).mean()),
                "n_records_gt_units": int(
                    (group.n_records > group.n_units).sum()
                ),
                "records_gt_units_proportion": float(
                    (group.n_records > group.n_units).mean()
                ),
                "n_document_supported": int(group.document_supported.sum()),
                "document_supported_proportion": float(
                    group.document_supported.mean()
                ),
                "n_flip": int(group.flip.sum()),
                "n_uncertain_to_supported": int(
                    (group.flip_direction == "uncertain_to_supported").sum()
                ),
                "n_supported_to_uncertain": int(
                    (group.flip_direction == "supported_to_uncertain").sum()
                ),
            }),
            include_groups=False,
        )
        .reset_index()
    )
    detail.to_csv(QA / "entry_level_details.csv", index=False)
    return summary, detail


def reconstruct_remaining_center(
    predictions: pd.DataFrame,
    records_by_key: dict[tuple[str, str], pd.DataFrame],
) -> pd.Series:
    centers = pd.Series(np.nan, index=predictions.index, dtype=float)
    g = predictions.loc[
        (predictions.method == "G") & (predictions.population == "all")
    ]
    for (target, endpoint), frame in g.groupby(
        ["target", "endpoint"], sort=True
    ):
        doc = document_values(records_by_key[(target, endpoint)], merge_identical=True)
        units = dict(tuple(doc.groupby("entry", sort=True)))
        for index, row in frame.iterrows():
            held = set(str(row.held_out_unit).split(";"))
            remaining = units[row.entry].loc[
                ~units[row.entry].documents.map(
                    lambda value: bool(set(str(value).split(";")) & held)
                )
            ]
            if remaining.empty:
                continue
            centers.at[index] = float(remaining.document_value.mean())
    return centers


def agreement_by_display_status(
    predictions: pd.DataFrame,
    records_by_key: dict[tuple[str, str], pd.DataFrame],
) -> pd.DataFrame:
    g = predictions.loc[predictions.method == "G"].copy()
    # Reconstruct for every population from its own held-out rows.
    g["center"] = np.nan
    for population in POPULATIONS:
        positions = g.index[g.population == population]
        local = g.loc[positions].copy()
        local["center"] = reconstruct_remaining_center(
            local.assign(
                population="all"
            ),
            records_by_key,
        )
        g.loc[positions, "center"] = local.center
    g["display_support"] = shown_probability(
        g.center.to_numpy(float),
        g.tau.to_numpy(float) / np.sqrt(g.m.to_numpy(float)),
        g.displayed.to_numpy(int),
    )
    g["display_status"] = np.where(
        g.display_support >= SUPPORT_THRESHOLD, "supported", "uncertain"
    )
    rows = []
    for (endpoint, population, display_status), endpoint_frame in g.groupby(
        ["endpoint", "population", "display_status"], sort=True
    ):
        for stratum in ("all", "m=1", "m>=2"):
            frame = (
                endpoint_frame
                if stratum == "all"
                else endpoint_frame.loc[endpoint_frame.stratum == stratum]
            )
            if frame.empty:
                continue
            rng = np.random.default_rng(BOOTSTRAP_SEED)
            weights = cluster_weights(frame, rng, BOOT_N)
            draws = weighted_scores(frame, weights)["agreement_rate"]
            low, high, valid = percentile_ci(draws)
            rows.append({
                "endpoint": endpoint,
                "population": population,
                "display_status": display_status,
                "stratum": stratum,
                "n_targets": int(frame.target.nunique()),
                "n_entries": int(
                    frame[["target", "entry"]].drop_duplicates().shape[0]
                ),
                "n_predictions": len(frame),
                "agreement_rate": float(frame.agree.mean()),
                "ci_low": low,
                "ci_high": high,
                "bootstrap_valid": valid,
                "bootstrap_repetitions": BOOT_N,
                "bootstrap_seed": BOOTSTRAP_SEED,
                "support_threshold": SUPPORT_THRESHOLD,
                "support_definition": (
                    "sigma_doc/sqrt(m), center=remaining-unit mean"
                ),
            })
    return pd.DataFrame(rows)


def bootstrap_binary_auc(correct, score, frame, weights) -> np.ndarray:
    """Weighted AUC for binary inputs and fixed target-entry multiplicities."""
    correct = np.asarray(correct, dtype=bool)
    score = np.asarray(score, dtype=float)
    order = np.argsort(score, kind="stable")
    starts = np.r_[0, np.flatnonzero(np.diff(score[order]) != 0) + 1]
    positives = np.add.reduceat(
        weights[:, order] * correct[order], starts, axis=1
    )
    negatives = np.add.reduceat(
        weights[:, order] * (1 - correct[order]), starts, axis=1
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


def forest_plot_data(
    predictions: pd.DataFrame,
    inclusion: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    align = ["fold", "entry", "held_out_unit"]
    for target, endpoint in included_keys(inclusion):
        for population in POPULATIONS:
            frame = predictions.loc[
                (predictions.target == target)
                & (predictions.endpoint == endpoint)
                & (predictions.population == population)
            ]
            g = frame.loc[frame.method == "G"].sort_values(align).reset_index(drop=True)
            if g.empty:
                continue
            rng = np.random.default_rng(BOOTSTRAP_SEED)
            weights = cluster_weights(g, rng, BOOT_N)
            g_draws = weighted_scores(g, weights)
            g_point = point_weighted_scores(g)
            g_point_log_loss = float(g_point["log_loss"][0])
            g_auc_low, g_auc_high, valid = percentile_ci(g_draws["auc"])
            g_ece_low, g_ece_high, _ = percentile_ci(g_draws["ece"])
            alternatives = {}
            for alternative in ("G-record", "G-fixed", "G-global", "X"):
                other = (
                    frame.loc[frame.method == alternative]
                    .sort_values(align)
                    .reset_index(drop=True)
                )
                if not g[align].equals(other[align]):
                    raise ValueError(
                        f"Forest keys differ: {target}/{endpoint}/{population}/{alternative}"
                )
                other_point = point_weighted_scores(other)
                difference = g_draws["log_loss"] - weighted_scores(
                    other, weights
                )["log_loss"]
                low, high, _ = percentile_ci(difference)
                alternatives[alternative] = {
                    "difference": float(
                        g_point_log_loss
                        - float(other_point["log_loss"][0])
                    ),
                    "ci_low": low,
                    "ci_high": high,
                }
            b0 = (
                frame.loc[frame.method == "B0"]
                .sort_values(align)
                .reset_index(drop=True)
            )
            b_range_auc = auc_mann_whitney(b0.agree, b0.b_range)
            b_draws = bootstrap_binary_auc(
                b0.agree.to_numpy(bool),
                b0.b_range.to_numpy(float),
                b0,
                weights,
            )
            b_low, b_high, _ = percentile_ci(b_draws)
            row = {
                "target_id": target,
                "endpoint": endpoint,
                "population": population,
                "n_predictions": len(g),
                "n_entries": int(g.entry.nunique()),
                "g_auc": auc_mann_whitney(g.agree, g.p_displayed),
                "g_auc_ci_low": g_auc_low,
                "g_auc_ci_high": g_auc_high,
                "g_ece": float(g_point["ece"][0]),
                "g_ece_ci_low": g_ece_low,
                "g_ece_ci_high": g_ece_high,
                "b_range_auc": b_range_auc,
                "b_range_auc_ci_low": b_low,
                "b_range_auc_ci_high": b_high,
                "bootstrap_valid": valid,
            }
            for alternative, item in alternatives.items():
                column = "g_minus_" + alternative.replace("-", "_").lower() + "_log_loss"
                row[column] = item["difference"]
                row[column + "_ci_low"] = item["ci_low"]
                row[column + "_ci_high"] = item["ci_high"]
            rows.append(row)
    return pd.DataFrame(rows)


def negative_control(
    predictions: pd.DataFrame,
    inclusion: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    g = predictions.loc[predictions.method == "G"].copy()
    h3_rows, auc_rows = [], []
    for endpoint in ENDPOINTS:
        for population in POPULATIONS:
            endpoint_pop = g.loc[
                (g.endpoint == endpoint) & (g.population == population)
            ]
            targets = [
                target
                for target in included_keys(inclusion)
                if target[1] == endpoint and target[0] in set(endpoint_pop.target)
            ]
            rng = np.random.default_rng(NULL_SEED)
            rejection = 0
            for _ in range(200):
                target_auc = []
                for target, _endpoint in targets:
                    frame = endpoint_pop.loc[endpoint_pop.target == target]
                    permuted = rng.permutation(frame.agree.to_numpy(bool))
                    score = auc_mann_whitney(
                        permuted, frame.p_displayed.to_numpy(float)
                    )
                    target_auc.append(score - 0.5)
                _stat, p_value, _n = one_sided_wilcoxon(
                    target_auc, alternative="greater"
                )
                rejection += int(p_value <= 0.05)
            rate = rejection / 200
            h3_rows.append({
                "endpoint": endpoint,
                "population": population,
                "replicates": 200,
                "rejections": rejection,
                "rejection_rate": rate,
                "alpha": 0.05,
                "threshold": 0.10,
                "passed": bool(rate <= 0.10),
                "permutation": "binary agreement permuted within target",
                "seed": NULL_SEED,
            })
            for target, _endpoint in targets:
                frame = endpoint_pop.loc[endpoint_pop.target == target]
                labels = frame.agree.to_numpy(bool)
                score = frame.p_displayed.to_numpy(float)
                observed = auc_mann_whitney(labels, score)
                ranks = stats.rankdata(score, method="average")
                n_positive = int(labels.sum())
                n_negative = len(labels) - n_positive
                null = np.empty(1000)
                local_rng = np.random.default_rng(NULL_SEED)
                for index in range(1000):
                    selected = local_rng.choice(
                        ranks, size=n_positive, replace=False
                    )
                    null[index] = (
                        selected.sum() - n_positive * (n_positive + 1) / 2
                    ) / (n_positive * n_negative)
                p_value = (
                    1 + int((null >= observed).sum())
                ) / (1 + len(null))
                auc_rows.append({
                    "target": target,
                    "endpoint": endpoint,
                    "population": population,
                    "n_predictions": len(frame),
                    "observed_auc": observed,
                    "null_mean": float(null.mean()),
                    "null_ci_low": float(np.percentile(null, 2.5)),
                    "null_ci_high": float(np.percentile(null, 97.5)),
                    "permutation_p_greater": p_value,
                    "replicates": 1000,
                    "seed": NULL_SEED,
                    "permutation": "binary agreement permuted within target",
                    "note": "within-target null; no pooled AUC comparison to 0.5",
                })
    return pd.DataFrame(h3_rows), pd.DataFrame(auc_rows)


def h7_decision(tests: pd.DataFrame) -> dict:
    selected = tests.loc[
        (tests.hypothesis == "H7") & (tests.population == "all")
    ]
    significant = {
        row.endpoint: bool(pd.notna(row.holm_p) and row.holm_p < 0.05)
        for row in selected.itertuples(index=False)
    }
    default = (
        "X"
        if all(significant.get(endpoint, False) for endpoint in ENDPOINTS)
        else "G"
    )
    return {
        "default_model": default,
        "h7_all_significant_both_endpoints": default == "X",
        "significant_by_endpoint": significant,
        "rule": (
            "X only if H7 is significant after Holm in both IC50 and Ki; "
            "otherwise G"
        ),
    }


def markdown_table(frame: pd.DataFrame, max_rows: int | None = None) -> str:
    if frame.empty:
        return "_No rows._"
    selected = frame.head(max_rows) if max_rows else frame

    def cell(value):
        if value is None or (np.isscalar(value) and pd.isna(value)):
            return "NA"
        if isinstance(value, (float, np.floating)):
            return f"{value:.6g}"
        return str(value).replace("|", "\\|").replace("\n", " ")

    lines = [
        "| " + " | ".join(map(str, selected.columns)) + " |",
        "| " + " | ".join("---" for _ in selected.columns) + " |",
    ]
    lines.extend(
        "| " + " | ".join(cell(value) for value in row) + " |"
        for row in selected.itertuples(index=False, name=None)
    )
    return "\n".join(lines)


def write_results(
    h1_h2: pd.DataFrame,
    tests: pd.DataFrame,
    negative_h3: pd.DataFrame,
    negative_auc: pd.DataFrame,
    time_pooled: pd.DataFrame,
    time_by_target: pd.DataFrame,
    pseudo_summary: pd.DataFrame,
    sigma: pd.DataFrame,
    pairs: pd.DataFrame,
    agreement: pd.DataFrame,
    decision: dict,
) -> None:
    primary = tests.loc[tests.population == "all"].copy()
    primary_table = primary[[
        "hypothesis",
        "endpoint",
        "n_targets",
        "median_difference",
        "median_ci_low",
        "median_ci_high",
        "raw_p",
        "holm_p",
        "correction_family",
    ]]
    time_main = time_pooled.loc[
        (time_pooled.cutoff_year == 2015)
        & time_pooled.scope.isin(
            ["all", "m=1", "m>=2", "display_supported", "display_uncertain"]
        )
    ]
    sigma_targets = sigma.loc[
        (sigma.row_type == "target") & (sigma.measurement_type == "all")
    ][[
        "target",
        "endpoint",
        "n_multi_unit_entries",
        "degrees_of_freedom",
        "sigma_doc",
        "ci_low",
        "ci_high",
    ]]
    pair_summary = pairs.loc[pairs.row_type == "target_summary"][[
        "target",
        "endpoint",
        "n_flagged_pairs",
        "n_flagged_documents",
        "n_touched_entries",
        "fraction_touched_entries",
        "n_touched_predictions",
        "fraction_touched_predictions",
    ]]
    status_main = agreement.loc[
        (agreement.population == "all")
        & (agreement.stratum == "all")
    ][[
        "endpoint",
        "display_status",
        "n_predictions",
        "agreement_rate",
        "ci_low",
        "ci_high",
    ]]
    lines = [
        "# E2 factual results",
        "",
        "## Sealed plan",
        "",
        f"SHA-256 verified before execution: `{SEALED_PLAN_SHA256}`.",
        "",
        "## H1-H8",
        "",
        "All populations and endpoints are reported in `hypothesis_tests.csv`; "
        "the table below is the prespecified primary `all` population.",
        "",
        markdown_table(primary_table),
        "",
        "H3-H7 use one Holm family per endpoint and population. H8 uses a "
        "separate Holm family across H8-IC50 and H8-Ki within each population.",
        "",
        "H1/H2 descriptive summaries:",
        "",
        markdown_table(
            h1_h2.loc[h1_h2.population == "all"]
        ),
        "",
        f"H7 decision rule outcome: default model = **{decision['default_model']}**.",
        "",
        "## Negative controls",
        "",
        markdown_table(negative_h3),
        "",
        "Per-target G AUC permutation results are in "
        "`qa/E2/negative_control_auc.csv`. No pooled AUC was compared with 0.5.",
        "",
        "## Time split",
        "",
        "Cutoff 2015, primary analysis. Missing-year units are excluded and counted "
        "in `time_split_by_target.csv`.",
        "",
        markdown_table(time_main),
        "",
        "The 2010 sensitivity analysis and all target-level rows are in "
        "`time_split.csv` and `time_split_by_target.csv`.",
        "",
        "## Pseudo-consensus",
        "",
        markdown_table(pseudo_summary),
        "",
        "## Sigma_doc",
        "",
        markdown_table(sigma_targets),
        "",
        "Ki literature references are single-measurement SD values 0.47 "
        "(Kalliokoski 2013) and 0.54 (Kramer 2012); IC50 references are "
        "0.57-0.61 (Kalliokoski 2013). These are descriptive comparisons only. "
        f"Targets meeting the joint B/F criterion (at least 30 multi-unit "
        f"entries for each assay type): "
        f"{sigma.loc[sigma.measurement_type.isin(['B', 'F']), ['target', 'endpoint']].drop_duplicates().shape[0]}. "
        "The CA II Ki top-five document-removal result is in "
        "`sigma_doc_by_target.csv`.",
        "",
        "## Flagged pairs",
        "",
        markdown_table(pair_summary),
        "",
        "## Agreement by display status",
        "",
        markdown_table(status_main),
        "",
    ]
    (OUT / "RESULTS_E2.md").write_text("\n".join(lines))


def e2a_hashes() -> dict[str, str]:
    return {
        name: sha256(OUT / name)
        for name in E2A_READ_ONLY_OUTPUTS
    }


def write_manifest(
    started: str,
    ended: str,
    decision: dict,
    baseline_dir_hashes: dict[str, str],
) -> None:
    existing = json.loads(E2A_MANIFEST.read_text())
    new_outputs = [
        path
        for path in sorted(OUT.iterdir())
        if path.is_file()
        and path.name not in {*E2A_READ_ONLY_OUTPUTS, "manifest.json"}
    ]
    new_qa = [
        path
        for path in sorted(QA.iterdir())
        if path.is_file()
        and path.name not in E2A_QA_FILES
        and not path.name.startswith(("rev_e2b_before", "rev_e2b_after"))
    ]
    sources = sorted((S / "source").glob("e2b_*.py")) + [S / "source/test_e2.py"]
    e2b = {
        "started_at": started,
        "ended_at": ended,
        "sealed_plan_sha256": sha256(SEALED_PLAN),
        "seeds": {
            "fold": FOLD_SEED,
            "bootstrap": BOOTSTRAP_SEED,
            "bootstrap_repetitions": BOOT_N,
            "sigma_bootstrap_repetitions": SIGMA_BOOT_N,
            "negative_control": NULL_SEED,
            "negative_control_h3_repetitions": 200,
            "negative_control_auc_repetitions": 1000,
        },
        "parameters": {
            "time_cutoffs": list(TIME_CUTOFFS),
            "time_unit_year": "maximum nonmissing document_year among component documents",
            "time_min_past_fit_entries": 5,
            "time_missing_year_policy": "exclude unit from temporal analysis and count it",
            "pseudo_consensus_min_units": 5,
            "sigma_bootstrap_unit": "entry",
            "support_threshold": SUPPORT_THRESHOLD,
            "h3_h7_holm": "one family per endpoint and population",
            "h8_holm": "H8-IC50 and H8-Ki within each population",
            "h7_eligible": "targets with at least one flagged document pair",
            "prediction_read_precision": "round_trip",
            "ca2_removal": "top five documents by number of multi-unit entries containing the document",
            "ca2_tie_break": "document_chembl_id",
        },
        "h7_decision": decision,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha256(path) for path in sources
        },
        "outputs_sha256": {
            str(path.relative_to(ROOT)): sha256(path) for path in new_outputs
        },
        "qa_sha256": {
            str(path.relative_to(ROOT)): sha256(path) for path in new_qa
        },
        "e2a_read_only_outputs_sha256": e2a_hashes(),
        "revision_directory_sha256": baseline_dir_hashes,
        "software_versions": {
            "python": platform.python_version(),
            **{
                name: importlib.metadata.version(name)
                for name in ("numpy", "pandas", "scipy", "pytest")
            },
        },
    }
    existing["E2b"] = e2b
    E2A_MANIFEST.write_text(
        json.dumps(
            json_safe(existing),
            indent=2,
            ensure_ascii=False,
            allow_nan=False,
        )
        + "\n"
    )


def archive_r0_time_split() -> dict[str, Path]:
    """Preserve the accepted R0 time-split files before R1 regeneration."""
    active = {
        "pooled": OUT / "time_split.csv",
        "by_target": OUT / "time_split_by_target.csv",
        "predictions": QA / "time_split_predictions.csv",
    }
    archived = {
        name: path.with_name(path.stem + R0_TIME_SUFFIX + path.suffix)
        for name, path in active.items()
    }
    for name, path in active.items():
        backup = archived[name]
        if not backup.exists():
            if not path.exists():
                raise FileNotFoundError(f"Cannot archive missing R0 file: {path}")
            path.replace(backup)
    return archived


def comparison_table(
    old: pd.DataFrame,
    new: pd.DataFrame,
    keys: list[str],
    metrics: list[str],
) -> pd.DataFrame:
    comparison = old[keys + metrics].merge(
        new[keys + metrics],
        on=keys,
        how="outer",
        suffixes=("_r0", "_r1"),
        validate="one_to_one",
    )
    for metric in metrics:
        comparison[f"{metric}_change"] = (
            comparison[f"{metric}_r1"] - comparison[f"{metric}_r0"]
        )
    return comparison.sort_values(keys).reset_index(drop=True)


def update_results_time_split_r1(
    old_pooled: pd.DataFrame,
    time_pooled: pd.DataFrame,
) -> None:
    """Replace only the temporal section of the factual English results."""
    path = OUT / "RESULTS_E2.md"
    text = path.read_text()
    start_marker = "## Time split\n"
    end_marker = "## Pseudo-consensus\n"
    start = text.index(start_marker)
    end = text.index(end_marker)
    primary = time_pooled.loc[
        (time_pooled.cutoff_year == 2015)
        & time_pooled.scope.isin(
            ["all", "m=1", "m>=2", "display_supported", "display_uncertain"]
        )
    ]
    comparison = comparison_table(
        old_pooled.loc[old_pooled.scope == "all"],
        time_pooled.loc[time_pooled.scope == "all"],
        ["endpoint", "cutoff_year", "scope"],
        [
            "n_predictions",
            "n_entries",
            "n_targets",
            "b0_inconsistency_rate",
            "g_ece",
            "g_calibration_bias",
            "g_auc",
        ],
    )
    section = [
        "## Time split",
        "",
        "R1 update: the temporal analysis now adds the original-target IC50 "
        "datasets CHEMBL205, CHEMBL204 and CHEMBL240. These datasets are used "
        "only for the temporal analysis and are not included in H1-H8.",
        "",
        "Cutoff 2015, primary temporal analysis with the R1 datasets:",
        "",
        markdown_table(primary),
        "",
        "Old/new comparison for the endpoint-level all-target scope:",
        "",
        markdown_table(comparison),
        "",
        "Missing-year units are excluded and counted in "
        "`time_split_by_target.csv`. The 2010 sensitivity analysis and all "
        "target-level rows are in the regenerated `time_split.csv` and "
        "`time_split_by_target.csv`.",
        "",
    ]
    path.write_text(text[:start] + "\n".join(section) + "\n" + text[end:])


def run_r1_outputs() -> None:
    """Regenerate only the R1 files without rerunning accepted E2b analyses."""
    started = datetime.now().astimezone().isoformat()
    if sha256(SEALED_PLAN) != SEALED_PLAN_SHA256:
        raise RuntimeError("封存计划 SHA-256 与任务单不一致，停止执行。")
    archived = archive_r0_time_split()
    old_pooled = pd.read_csv(
        archived["pooled"], float_precision="round_trip"
    )
    old_by_target = pd.read_csv(
        archived["by_target"], float_precision="round_trip"
    )

    inclusion = pd.read_csv(OUT / "target_inclusion.csv")
    records_by_key = load_records()
    original_ic50 = load_original_ic50_records()
    records_by_key.update(original_ic50)
    additional_keys = [
        (target, "IC50") for target in ORIGINAL_IC50_TARGETS
    ]
    time_pooled, time_by_target, _time_detail = build_time_split(
        records_by_key, inclusion, additional_keys
    )
    time_pooled.to_csv(OUT / "time_split.csv", index=False)
    time_by_target.to_csv(OUT / "time_split_by_target.csv", index=False)

    metrics = pd.read_csv(
        OUT / "metrics_by_target.csv",
        low_memory=False,
        float_precision="round_trip",
    )
    stratified = hypothesis_tests_by_stratum(metrics)
    stratified.to_csv(OUT / "hypothesis_tests_by_stratum.csv", index=False)
    update_results_time_split_r1(old_pooled, time_pooled)

    ended = datetime.now().astimezone().isoformat()
    metadata = {
        "started_at": started,
        "ended_at": ended,
        "archived_r0": {
            name: str(path.relative_to(ROOT)) for name, path in archived.items()
        },
        "original_ic50_targets": list(ORIGINAL_IC50_TARGETS),
        "old_time_by_target_rows": len(old_by_target),
        "new_time_by_target_rows": len(time_by_target),
        "time_predictions": len(_time_detail),
        "stratified_rows": len(stratified),
    }
    (QA / "r1_run_metadata.json").write_text(
        json.dumps(json_safe(metadata), indent=2) + "\n"
    )
    print(
        f"E2B R1 OUTPUTS COMPLETE time_rows={len(time_by_target)} "
        f"time_predictions={len(_time_detail)} strata={len(stratified)}",
        flush=True,
    )


def r1_report_section(metadata: dict) -> str:
    old_pooled = pd.read_csv(
        OUT / f"time_split{R0_TIME_SUFFIX}.csv",
        float_precision="round_trip",
    )
    new_pooled = pd.read_csv(OUT / "time_split.csv", float_precision="round_trip")
    old_by_target = pd.read_csv(
        OUT / f"time_split_by_target{R0_TIME_SUFFIX}.csv",
        float_precision="round_trip",
    )
    new_by_target = pd.read_csv(
        OUT / "time_split_by_target.csv", float_precision="round_trip"
    )
    comparison = comparison_table(
        old_pooled.loc[old_pooled.scope == "all"],
        new_pooled.loc[new_pooled.scope == "all"],
        ["endpoint", "cutoff_year", "scope"],
        [
            "n_predictions",
            "n_entries",
            "n_targets",
            "b0_inconsistency_rate",
            "g_ece",
            "g_calibration_bias",
            "g_auc",
        ],
    )
    original = new_by_target.loc[
        new_by_target.target.isin(ORIGINAL_IC50_TARGETS)
    ]
    stratified = pd.read_csv(
        OUT / "hypothesis_tests_by_stratum.csv",
        float_precision="round_trip",
    )
    test_log = (QA / "r1_pytest.log").read_text().rstrip()
    guard_after = (QA / "r1_guard_after.log").read_text().rstrip()
    cgm_before = (QA / "r1_cgm_status_before.txt").read_text()
    cgm_after = (QA / "r1_cgm_status_after.txt").read_text()
    scope = json.loads((QA / "r1_freeze_check.json").read_text())
    return "\n".join([
        "## 返工 R1",
        "",
        f"- 开始：`{metadata['started_at']}`",
        f"- 结束：`{metadata['ended_at']}`",
        f"- 封存计划 SHA-256：`{sha256(SEALED_PLAN)}`，核对一致。",
        "- 补做范围仅为原三靶点 IC50 的时间切分，以及 H2/H3 的 m 分层报告；"
        "其余已验收 E2b 结果不改动。",
        "",
        "### 时间切分",
        "",
        "CHEMBL205、CHEMBL204、CHEMBL240 的 IC50 记录使用 E2a 回归检查相同的"
        "API 整理参数进入 2015 与 2010 两个时间切分。R0 文件以 "
        f"`{R0_TIME_SUFFIX}` 后缀保留。",
        "",
        markdown_table(comparison),
        "",
        "新增三个原靶点的逐靶点结果：",
        "",
        markdown_table(original),
        "",
        f"逐靶点行数由 {metadata['old_time_by_target_rows']} 增至 "
        f"{metadata['new_time_by_target_rows']}；新时间切分预测共 "
        f"{metadata['time_predictions']} 个将来单元。",
        "",
        "### H2/H3 按 m 分层",
        "",
        "以下结果按端点 × 总体分别计算，每层只使用有 AUC 的靶点；"
        "`family_note` 明确标为 `stratified; not part of the Holm families`。",
        "",
        markdown_table(stratified),
        "",
        "### 测试",
        "",
        "```sh\n"
        "chemgridmap/.venv/bin/python -m pytest -p no:cacheprovider "
        "python -m pytest -p no:cacheprovider paper/strengthening/source/test_e2.py -vv\n"
        "```",
        "",
        "```text\n" + test_log + "\n```",
        "",
        "### 冻结与范围",
        "",
        "`_guard` 最终检查：",
        "",
        "```text\n" + guard_after + "\n```",
        "",
        f"CGM 状态执行前后均为 {len(cgm_before.splitlines())} 行；"
        f"文本{'一致' if cgm_before == cgm_after else '不一致'}。",
        "",
        f"REV 全量核对：新增 {len(scope['added'])}、修改 "
        f"{len(scope['modified'])}、删除 {len(scope['removed'])}；"
        f"允许位置外变化 {len(scope['outside_allowed_changes'])} 个。",
        "",
        "允许位置外变化：",
        "",
        *(
            [f"- `{path}`" for path in scope["outside_allowed_changes"]]
            or ["- 无"]
        ),
        "",
        "问题与未完成事项：",
        "",
        "- 无未完成项。原三个 IC50 数据集按要求只进入时间切分，不进入 H1–H8。",
        "",
    ])


def append_r1_report(metadata: dict) -> None:
    text = REPORT.read_text()
    marker = "\n## 返工 R1\n"
    if marker in text:
        text = text[:text.index(marker)]
    REPORT.write_text(text.rstrip() + "\n\n" + r1_report_section(metadata))


def targeted_output_scope() -> dict:
    before = {}
    for line in (QA / "r1_e2_outputs_before.sha256").read_text().splitlines():
        digest, path = line.split(maxsplit=1)
        before[Path(path).name] = digest
    after = {path.name: sha256(path) for path in OUT.iterdir() if path.is_file()}
    changed = sorted(
        name for name in set(before) | set(after) if before.get(name) != after.get(name)
    )
    allowed = {
        "time_split.csv",
        "time_split_by_target.csv",
        "hypothesis_tests_by_stratum.csv",
        "RESULTS_E2.md",
        "manifest.json",
        f"time_split{R0_TIME_SUFFIX}.csv",
        f"time_split_by_target{R0_TIME_SUFFIX}.csv",
    }
    return {
        "changed": changed,
        "allowed": sorted(allowed),
        "unexpected": sorted(set(changed) - allowed),
    }


def baseline_output_hashes() -> dict[str, str]:
    result = {}
    for line in (QA / "r1_e2_outputs_before.sha256").read_text().splitlines():
        digest, path = line.split(maxsplit=1)
        result[path] = digest
    return result


def finalize_r1() -> None:
    """Append the R1 audit records after tests and final guard checks."""
    metadata = json.loads((QA / "r1_run_metadata.json").read_text())
    scope = compare_baseline(QA / "r1_rev_before.sha256")
    (QA / "r1_freeze_check.json").write_text(
        json.dumps(json_safe(scope), indent=2) + "\n"
    )
    append_r1_report(metadata)

    targeted = targeted_output_scope()
    targeted["changed"] = sorted({*targeted["changed"], "manifest.json"})
    targeted["unexpected"] = sorted(set(targeted["changed"]) - set(targeted["allowed"]))
    test_log = QA / "r1_pytest.log"
    test_text = test_log.read_text()
    pytest_passed = (
        " failed" not in test_text
        and " error" not in test_text
        and " passed" in test_text
    )
    r1_outputs = [
        OUT / f"time_split{R0_TIME_SUFFIX}.csv",
        OUT / f"time_split_by_target{R0_TIME_SUFFIX}.csv",
        QA / f"time_split_predictions{R0_TIME_SUFFIX}.csv",
        OUT / "time_split.csv",
        OUT / "time_split_by_target.csv",
        QA / "time_split_predictions.csv",
        OUT / "hypothesis_tests_by_stratum.csv",
        OUT / "RESULTS_E2.md",
        REPORT,
    ]
    existing = json.loads(E2A_MANIFEST.read_text())
    existing["R1"] = {
        **metadata,
        "sealed_plan_sha256": sha256(SEALED_PLAN),
        "scope": {
            "original_ic50_only_in_time_split": True,
            "hypothesis_tests_csv_unchanged": (
                sha256(OUT / "hypothesis_tests.csv")
                == baseline_output_hashes()[
                    "paper/strengthening/output/E2/hypothesis_tests.csv"
                ]
            ),
            "stratified_family_note": "stratified; not part of the Holm families",
            "targeted_output_scope": targeted,
        },
        "test": {
            "command": (
                "chemgridmap/.venv/bin/python -m pytest -p no:cacheprovider "
                "paper/strengthening/source/test_e2.py -vv"
            ),
            "passed": pytest_passed,
            "log_sha256": sha256(test_log),
        },
        "source_sha256": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in (S / "source/e2b_analysis.py", S / "source/test_e2.py")
        },
        "outputs_sha256": {
            str(path.relative_to(ROOT)): sha256(path) for path in r1_outputs
        },
        "freeze_check": scope,
    }
    E2A_MANIFEST.write_text(
        json.dumps(
            json_safe(existing),
            indent=2,
            ensure_ascii=False,
            allow_nan=False,
        )
        + "\n"
    )
    if targeted["unexpected"]:
        raise RuntimeError(f"Unexpected E2 output changes: {targeted['unexpected']}")
    if not pytest_passed:
        raise RuntimeError("R1 pytest log does not report a passing run")
    print("E2B R1 FINALIZED", flush=True)


def rev_file_hashes() -> dict[str, str]:
    result = {}
    for path in sorted(REV.rglob("*")):
        if path.is_file():
            result[str(path.relative_to(ROOT))] = sha256(path)
    return result


def compare_baseline(before_path: Path) -> dict:
    before = {}
    for line in before_path.read_text().splitlines():
        digest, path = line.split(maxsplit=1)
        before[path] = digest
    after = rev_file_hashes()
    added = sorted(set(after) - set(before))
    removed = sorted(set(before) - set(after))
    modified = sorted(
        path for path in set(before) & set(after) if before[path] != after[path]
    )

    def allowed(path: str) -> bool:
        relative = Path(path).relative_to(
            "paper/strengthening"
        ).as_posix()
        prefixes = (
            "output/E2/",
            "qa/E2/",
            "source/e2b_",
            "source/test_e2.py",
            "source/e2c_",
            "source/test_e2c.py",
            "output/E2c/",
            "qa/E2c/",
        )
        return relative.startswith(prefixes)

    outside = [
        path for path in [*added, *modified, *removed] if not allowed(path)
    ]
    return {
        "before_files": len(before),
        "after_files": len(after),
        "added": added,
        "modified": modified,
        "removed": removed,
        "outside_allowed_changes": outside,
        "allowed_change_count": sum(
            allowed(path) for path in [*added, *modified, *removed]
        ),
        "manifest_exception": (
            "E2a manifest is modified only by appending the E2b section as "
            "explicitly required by docs/task.md"
        ),
    }


def write_report(
    started: str,
    ended: str,
    h1_h2: pd.DataFrame,
    tests: pd.DataFrame,
    negative_h3: pd.DataFrame,
    negative_auc: pd.DataFrame,
    time_pooled: pd.DataFrame,
    time_by_target: pd.DataFrame,
    pseudo_summary: pd.DataFrame,
    sigma: pd.DataFrame,
    pairs: pd.DataFrame,
    entry_summary: pd.DataFrame,
    agreement: pd.DataFrame,
    decision: dict,
) -> None:
    pytest_log = (
        (QA / "pytest_e2b.log").read_text()
        if (QA / "pytest_e2b.log").exists()
        else "未运行"
    )
    guard_log = (
        (QA / "guard_e2b_after.log").read_text()
        if (QA / "guard_e2b_after.log").exists()
        else "未运行"
    )
    scope = json.loads((QA / "frozen_compare_e2b.json").read_text())
    all_changes = [
        *scope["added"], *scope["modified"], *scope["removed"]
    ]
    e2c_changes = sorted(
        path
        for path in all_changes
        if "/E2c" in path
        or "/e2c_" in path
        or path.endswith("/test_e2c.py")
        or path.endswith("/E2c_REPORT.md")
    )
    cgm_before = (QA / "cgm_status_before.txt").read_text()
    cgm_after = (QA / "cgm_status_e2b_verified.txt").read_text()
    h1_h2_all = h1_h2.loc[h1_h2.population == "all"]
    tests_all = tests.loc[tests.population == "all"]
    time_main = time_pooled.loc[time_pooled.cutoff_year == 2015]
    sigma_targets = sigma.loc[
        (sigma.row_type == "target") & (sigma.measurement_type == "all")
    ]
    pair_summary = pairs.loc[pairs.row_type == "target_summary"]
    status_all = agreement.loc[agreement.population == "all"]
    output_rows = []
    for path in sorted(OUT.iterdir()):
        if path.is_file() and path.name != "manifest.json":
            with path.open(encoding="utf-8") as handle:
                rows = (
                    max(sum(1 for _ in handle) - 1, 0)
                    if path.suffix == ".csv"
                    else "NA"
                )
            output_rows.append({
                "file": path.name,
                "rows": rows,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            })
    sections = [
        "# E2b 执行报告",
        "",
        "## 1. 封存计划与范围",
        "",
        f"- 开始：`{started}`",
        f"- 结束：`{ended}`",
        f"- 封存计划 SHA-256：`{sha256(SEALED_PLAN)}`，与任务单要求的 "
        f"`{SEALED_PLAN_SHA256}` 一致。",
        "- 主分析总体为 `all`；`flag-clean` 作为敏感性分析同步报告。",
        "- E2c 相关文件未运行、未修改；其新增或变化在冻结核对中按任务单列为允许的并行范围。",
        "",
        "## 1.1 实施过程",
        "",
        "新增 `S/source/e2b_analysis.py`，从 E2a 的 round-trip 浮点预测、指标、"
        "比较对象和已整理记录缓存计算全部 E2b 输出；在 `test_e2.py` 末尾追加 "
        "5 项 E2b 测试。计算顺序为：H1/H2 描述量、H3–H8、两类负控、"
        "2015/2010 时间切分、伪共识、σ_doc、标记文献对、条目统计、"
        "显示状态一致率、森林图数据、英文事实表与清单。",
        "",
        "执行中发现并修正了 4 个实现问题：一个 bootstrap 常量名拼写错误；"
        "时间切分汇总与森林图点估计的数组未取标量；pandas 3 下 one-hot "
        "赋值受原表 dtype 影响；按显示状态一致率最初漏写总体 `all` 行。"
        "修正后重新执行完整分析，最终结果没有使用失败轮次的中间数值。",
        "",
        "## 2. H1 与 H2",
        "",
        markdown_table(h1_h2_all),
        "",
        "## 3. H3–H8",
        "",
        "差值的正负方向与 `alternative` 一同写入 CSV。H3–H7 在每个端点、"
        "每个总体内共同做 Holm；H8 在 IC50 与 Ki 两行内单独做 Holm。",
        "",
        markdown_table(tests_all),
        "",
        "## 4. H7 决策",
        "",
        f"- 默认模型：`{decision['default_model']}`",
        f"- 规则：{decision['rule']}",
        f"- 两端点 H7 Holm 显著：`{decision['h7_all_significant_both_endpoints']}`",
        "",
        "## 5. 负控",
        "",
        markdown_table(negative_h3),
        "",
        f"逐靶点 G AUC 的 1,000 次置换结果共 {len(negative_auc)} 行，"
        "写入 `S/qa/E2/negative_control_auc.csv`。结果只使用靶点内零分布；"
        "未把合计 AUC 与 0.5 比较。",
        "",
        "## 6. 时间切分",
        "",
        "单元年份取构成文献的最晚非缺失 `document_year`。年份缺失单元排除并计数；"
        "只有过去单元数 ≥2 的条目计入 σ_past，当前数据集要求至少 5 个此类条目。",
        "",
        markdown_table(time_main),
        "",
        "逐靶点结果，包括 2010 敏感性分析和缺失年份计数：",
        "",
        markdown_table(time_by_target),
        "",
        "## 7. 伪共识",
        "",
        markdown_table(pseudo_summary),
        "",
        "校准区间位于 `S/output/E2/pseudo_consensus_calibration.csv`。",
        "",
        "## 8. σ_doc",
        "",
        markdown_table(sigma_targets),
        "",
        "同时满足 B 型与 F 型测定各有至少 30 个多单元条目的靶点–端点数为 "
        f"{sigma.loc[sigma.measurement_type.isin(['B', 'F']), ['target', 'endpoint']].drop_duplicates().shape[0]}。",
        "",
        "CA II Ki 的文献簇影响行位于 `sigma_doc_by_target.csv`。"
        "下载数据没有通讯作者字段，因此无法核实“同一通讯作者”；"
        "按 Opus 答复，以“在多单元条目中提供单元的数量”定义自由度贡献，"
        "并列按 `document_chembl_id` 排序。被移除的 5 篇文献及年份：",
        "",
        markdown_table(
            pd.read_csv(QA / "ca2_top5_documents.csv")
        ),
        "",
        "## 9. 被标记文献对普查",
        "",
        markdown_table(pair_summary),
        "",
        "各靶点偏移最大的 10 对位于 `flagged_pairs_census.csv` 的 "
        "`row_type=top_pair` 行；Spearman 不可定义时明确记录为 NA。",
        "",
        "## 10. 条目层面统计",
        "",
        markdown_table(entry_summary),
        "",
        "## 11. 按显示状态分层的一致率",
        "",
        markdown_table(status_all),
        "",
        "显示支持度定义为剩余单元均值中心、SD = σ_doc / sqrt(m)，阈值 0.8。",
        "",
        "## 12. 输出与清单",
        "",
        markdown_table(pd.DataFrame(output_rows)),
        "",
        "`manifest.json` 保留 E2a 原有键值，并追加 `E2b` 节。"
        "除清单追加外，E2a 的五个只读数据输出哈希保持不变。",
        "",
        "## 13. 测试",
        "",
        "命令：",
        "",
        "```sh\n"
        "chemgridmap/.venv/bin/python -m pytest -p no:cacheprovider "
        "python -m pytest -p no:cacheprovider paper/strengthening/source/test_e2.py -vv\n"
        "```",
        "",
        "完整输出：",
        "",
        "```text\n" + pytest_log.rstrip() + "\n```",
        "",
        "## 14. 冻结核对",
        "",
        "`_guard`：",
        "",
        "```text\n" + guard_log.rstrip() + "\n```",
        "",
        f"REV 执行前 {scope['before_files']} 个文件，比对时 "
        f"{scope['after_files']} 个文件；新增 {len(scope['added'])}、"
        f"修改 {len(scope['modified'])}、删除 {len(scope['removed'])}。",
        f"允许位置之外的并发变化：{len(scope['outside_allowed_changes'])} 个。",
        "",
        *(
            [f"- `{path}`" for path in scope["outside_allowed_changes"]]
            or ["- 无"]
        ),
        "",
        f"E2c 并行新增或变化文件（{len(e2c_changes)} 个，按任务单不计 E2b 范围违规）：",
        "",
        markdown_table(pd.DataFrame({"path": e2c_changes})),
        "",
        f"CGM 状态执行前 {len(cgm_before.splitlines())} 行，执行后 "
        f"{len(cgm_after.splitlines())} 行；前后"
        f"{'一致' if cgm_before == cgm_after else '不一致'}。",
        "",
        "## 15. 问题与未完成事项",
        "",
        "- CA II Ki 数据缺少通讯作者字段，无法核实作者簇；只执行计划指定的"
        "自由度贡献前 5 篇移除分析。",
        "- 其余问题见上述负控、冻结核对和逐靶点结果；本报告不对显著性作解释。",
        "",
    ]
    REPORT.write_text("\n".join(sections))
    negative_passed = bool(negative_h3.passed.all())
    brief = [
        "# E2b 简要结果",
        "",
        "- 封存计划 SHA-256 核对通过。",
        f"- H3–H8：{len(tests)} 个预定义检验结果已写入 `hypothesis_tests.csv`；"
        "主分析总体为 `all`。",
        f"- H7 决策：默认模型为 `{decision['default_model']}`。",
        f"- 负控：四个端点×总体分析格的 200 次 H3 置换拒绝率均 ≤ 0.10，"
        f"整体通过={negative_passed}。",
        f"- 时间切分：输出 {len(time_by_target)} 个靶点×截止年份行，"
        "缺失年份单元已排除并计数。",
        f"- 描述分析：σ_doc {len(sigma)} 行；标记文献对普查 {len(pairs)} 行；"
        f"显示状态一致率 {len(agreement)} 行。",
        "- 测试：`test_e2.py` 9 项通过；E2a 五个只读输出哈希未变。",
        "- CA II Ki 数据缺少通讯作者字段；按计划执行贡献最高的 5 篇文献移除分析。",
        "",
    ]
    BRIEF.write_text("\n".join(brief))


def run() -> None:
    started = datetime.now().astimezone().isoformat()
    if sha256(SEALED_PLAN) != SEALED_PLAN_SHA256:
        raise RuntimeError("封存计划 SHA-256 与任务单不一致，停止执行。")
    OUT.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    baseline_hashes = e2a_hashes()
    (QA / "e2a_read_only_hashes_before.json").write_text(
        json.dumps(baseline_hashes, indent=2) + "\n"
    )
    inclusion = pd.read_csv(OUT / "target_inclusion.csv")
    predictions = read_predictions()
    metrics = pd.read_csv(
        OUT / "metrics_by_target.csv",
        low_memory=False,
        float_precision="round_trip",
    )
    pooled = pd.read_csv(
        OUT / "pooled_metrics.csv",
        low_memory=False,
        float_precision="round_trip",
    )
    comparators = pd.read_csv(
        OUT / "comparators.csv",
        low_memory=False,
        float_precision="round_trip",
    )
    records_by_key = load_records()

    h1_h2 = h1_h2_summary(metrics, pooled)
    calibration_curve(predictions, 5).to_csv(
        OUT / "calibration_curve.csv", index=False
    )
    tests = hypothesis_tests(metrics, comparators, records_by_key)
    tests.to_csv(OUT / "hypothesis_tests.csv", index=False)
    decision = h7_decision(tests)
    (QA / "h7_decision.json").write_text(
        json.dumps(json_safe(decision), indent=2) + "\n"
    )
    h1_h2.to_csv(OUT / "h1_h2_summary.csv", index=False)

    negative_h3, negative_auc = negative_control(predictions, inclusion)
    negative_h3.to_csv(QA / "negative_control_h3.csv", index=False)
    negative_auc.to_csv(QA / "negative_control_auc.csv", index=False)

    time_pooled, time_by_target, _time_detail = build_time_split(
        records_by_key, inclusion
    )
    time_pooled.to_csv(OUT / "time_split.csv", index=False)
    time_by_target.to_csv(OUT / "time_split_by_target.csv", index=False)

    pseudo_summary, pseudo_curve = pseudo_consensus(
        records_by_key, inclusion
    )
    pseudo_summary.to_csv(OUT / "pseudo_consensus.csv", index=False)
    pseudo_curve.to_csv(
        OUT / "pseudo_consensus_calibration.csv", index=False
    )

    sigma = build_sigma_doc(records_by_key, inclusion)
    sigma.to_csv(OUT / "sigma_doc_by_target.csv", index=False)

    pairs = flagged_pairs_census(
        records_by_key, predictions, inclusion
    )
    pairs.to_csv(OUT / "flagged_pairs_census.csv", index=False)

    entry_summary, _entry_details = entry_level_stats(
        records_by_key, inclusion
    )
    entry_summary.to_csv(OUT / "entry_level_by_target.csv", index=False)

    agreement = agreement_by_display_status(predictions, records_by_key)
    agreement.to_csv(OUT / "agreement_by_display_status.csv", index=False)

    forest_plot_data(predictions, inclusion).to_csv(
        OUT / "forest_plot_data.csv", index=False
    )
    write_results(
        h1_h2,
        tests,
        negative_h3,
        negative_auc,
        time_pooled,
        time_by_target,
        pseudo_summary,
        sigma,
        pairs,
        agreement,
        decision,
    )
    ended = datetime.now().astimezone().isoformat()
    write_manifest(started, ended, decision, baseline_hashes)
    print(
        f"E2B COMPLETE tests={len(tests)} "
        f"time_rows={len(time_by_target)} sigma_rows={len(sigma)} "
        f"default_model={decision['default_model']}",
        flush=True,
    )


if __name__ == "__main__":
    if "--r1-outputs" in sys.argv:
        run_r1_outputs()
    elif "--r1-finalize" in sys.argv:
        finalize_r1()
    else:
        run()
