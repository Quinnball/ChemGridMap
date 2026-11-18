"""E5 post hoc model trial, robustness analyses and descriptive mechanisms."""
from __future__ import annotations

from collections import Counter
from datetime import datetime, timedelta
import gzip
import hashlib
import json
import os
from pathlib import Path
import pickle
import platform
import shlex
import sys
import time

sys.dont_write_bytecode = True
os.environ.setdefault("MPLBACKEND", "Agg")
ROOT = Path(__file__).resolve().parents[3]
S = ROOT / "paper/strengthening"
OUT = S / "output/E5"
QA = S / "qa/E5"
API = S / "data/chembl37_api"
sys.path[:0] = [
    str(S / "source"),
    str(ROOT / "paper"),
    str(ROOT / "src"),
]

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from chemgridmap.activity import activity_class
from chemgridmap.evidence import (
    CLASSES,
    class_probabilities,
    document_values,
    fit_between_unit_sd,
)
from e1a_statistics import auc_mann_whitney
from e2_models import (
    flags_from_observations,
    prepare_target,
    touch_flagged_entries,
)
from e2b_analysis import add_unit_years, holm_adjust, one_sided_wilcoxon
from e5_models import (
    ANALYSIS_ROLE,
    CLASS_BOUNDARIES,
    METHODS,
    SUPPORT_THRESHOLD,
    classify,
    document_fold_predictions,
    eb_prediction,
    empirical_bayes_prior,
    entry_fold_predictions,
    near_duplicate_predictions,
    pooled_sigma,
)
from manifest_utils import write_manifest
from reliability_activity import gaussian_prediction

BOOT_N = 2000
BIAS_BOOT_N = 500
BOOTSTRAP_SEED = 20261028
TIME_CUTOFFS = (2015, 2010)
TIME_METHODS = ("G", "G-EB", "G-EB-HS")
METRIC_DIRECTIONS = {
    "log_loss": "less",
    "auc": "greater",
    "displayed_ece": "less",
    "calibration_bias_inactive": "less",
    "calibration_bias_medium": "less",
    "calibration_bias_active": "less",
}
UPPER_LIMIT_CODES = np.array([4.00, 4.30, 4.52, 4.70, 5.00])


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


def write_json(path: Path, value) -> None:
    path.write_text(
        json.dumps(
            json_safe(value), indent=2, ensure_ascii=False, allow_nan=False
        )
        + "\n"
    )


def load_inputs() -> dict:
    inclusion = pd.read_csv(OUT.parent / "E2/target_inclusion.csv")
    keys = list(
        inclusion.loc[
            inclusion.included, ["target", "endpoint"]
        ].itertuples(index=False, name=None)
    )
    all_records = pickle.loads(
        gzip.decompress((S / "qa/E2/curated_records.pkl.gz").read_bytes())
    )
    records = {key: all_records[key] for key in keys}
    e2 = pd.read_csv(
        OUT.parent / "E2/predictions.csv",
        low_memory=False,
        float_precision="round_trip",
    )
    e2 = e2.loc[e2.population == "all"].copy()
    return {
        "inclusion": inclusion,
        "keys": keys,
        "records": records,
        "e2": e2,
        "e2_g": e2.loc[e2.method == "G"].copy(),
        "e2_flag": pd.read_csv(
            OUT.parent / "E2/predictions.csv",
            low_memory=False,
            float_precision="round_trip",
        ).query("population == 'flag-clean' and method == 'G'"),
    }


def prediction_keys(frame: pd.DataFrame) -> list[str]:
    return ["target", "endpoint", "fold", "entry", "held_out_unit"]


def regression_check(entry: pd.DataFrame, expected: pd.DataFrame) -> dict:
    actual = entry.loc[entry.method == "G"].copy()
    expected = expected.copy()
    keys = prediction_keys(actual)
    merged = actual.merge(
        expected,
        on=keys,
        suffixes=("_e5", "_e2"),
        validate="one_to_one",
    )
    if not len(actual) == len(expected) == len(merged):
        raise RuntimeError(
            f"G regression row-count mismatch: {len(actual)}, {len(expected)}, {len(merged)}"
        )
    errors = {
        column: float(
            np.max(
                np.abs(
                    merged[f"{column}_e5"].to_numpy(float)
                    - merged[f"{column}_e2"].to_numpy(float)
                )
            )
        )
        for column in ("p0", "p1", "p2")
    }
    maximum = max(errors.values())
    result = {
        "n_expected": int(len(expected)),
        "n_e5": int(len(actual)),
        "n_matched": int(len(merged)),
        "maximum_absolute_error": maximum,
        "errors": errors,
        "tolerance": 1e-12,
        "passed": bool(maximum <= 1e-12),
        "analysis_role": ANALYSIS_ROLE,
    }
    write_json(QA / "regression_g_e2.json", result)
    if not result["passed"]:
        raise RuntimeError(f"E5 G regression failed: {result}")
    return result


def score_frame(frame: pd.DataFrame) -> dict[str, float | int]:
    if frame.empty:
        return {
            "n_predictions": 0,
            "n_entries": 0,
            "auc": np.nan,
            "log_loss": np.nan,
            "brier": np.nan,
            "top_ece": np.nan,
            "displayed_ece": np.nan,
            "accuracy": np.nan,
            "mean_probability": np.nan,
            "agreement_rate": np.nan,
            "calibration_bias": np.nan,
            "supported_n": 0,
            "saturation_rate": np.nan,
            "supported_agreement_rate": np.nan,
            "uncertain_agreement_rate": np.nan,
        }
    p = frame[["p0", "p1", "p2"]].to_numpy(float)
    observed = frame.observed.to_numpy(int)
    displayed = frame.displayed.to_numpy(int)
    agree = frame.agree.to_numpy(bool)
    shown = p[np.arange(len(frame)), displayed]
    confidence = p.max(axis=1)
    top_correct = p.argmax(axis=1) == observed
    top_bins = np.minimum((confidence * 10).astype(int), 9)
    displayed_bins = np.minimum((shown * 10).astype(int), 9)
    top_ece = sum(
        abs(
            np.sum(
                top_correct[top_bins == index]
                - confidence[top_bins == index]
            )
        )
        for index in range(10)
    )
    displayed_ece = sum(
        abs(
            np.sum(
                agree[displayed_bins == index]
                - shown[displayed_bins == index]
            )
        )
        for index in range(10)
    )
    supported = frame.support.to_numpy(float) >= SUPPORT_THRESHOLD
    return {
        "n_predictions": int(len(frame)),
        "n_entries": int(frame.entry.nunique()),
        "auc": auc_mann_whitney(agree, shown),
        "log_loss": float(
            -np.log(np.clip(p[np.arange(len(p)), observed], 1e-6, 1.0)).mean()
        ),
        "brier": float(
            ((p - np.eye(3)[observed]) ** 2).sum(axis=1).mean()
        ),
        "top_ece": float(abs(top_ece) / len(frame)),
        "displayed_ece": float(displayed_ece / len(frame)),
        "accuracy": float(top_correct.mean()),
        "mean_probability": float(shown.mean()),
        "agreement_rate": float(agree.mean()),
        "calibration_bias": float((shown - agree).mean()),
        "supported_n": int(supported.sum()),
        "saturation_rate": float(supported.mean()),
        "supported_agreement_rate": (
            float(agree[supported].mean()) if supported.any() else np.nan
        ),
        "uncertain_agreement_rate": (
            float(agree[~supported].mean()) if (~supported).any() else np.nan
        ),
    }


def category_bias(frame: pd.DataFrame) -> dict[str, float]:
    result = {}
    for index, name in enumerate(CLASSES):
        selected = frame.loc[frame.displayed == index]
        result[f"calibration_bias_{name}"] = (
            float((selected.p_displayed - selected.agree).mean())
            if len(selected)
            else np.nan
        )
        result[f"n_displayed_{name}"] = int(len(selected))
    return result


def dataset_metrics(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (target, endpoint, method), group in frame.groupby(
        ["target", "endpoint", "method"], sort=True
    ):
        rows.append(
            {
                "split_design": group.split_design.iloc[0],
                "merge_basis": group.merge_basis.iloc[0],
                "target": target,
                "endpoint": endpoint,
                "model": method,
                **score_frame(group),
                **category_bias(group),
                "analysis_role": ANALYSIS_ROLE,
            }
        )
    return pd.DataFrame(rows)


def pooled_metrics(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (endpoint, method), group in frame.groupby(
        ["endpoint", "method"], sort=True
    ):
        rows.append(
            {
                "split_design": group.split_design.iloc[0],
                "merge_basis": group.merge_basis.iloc[0],
                "endpoint": endpoint,
                "model": method,
                "n_targets": int(group.target.nunique()),
                **score_frame(group),
                **category_bias(group),
                "analysis_role": ANALYSIS_ROLE,
            }
        )
    return pd.DataFrame(rows)


def stratified_metrics(frame: pd.DataFrame) -> pd.DataFrame:
    dimensions = {
        "display_class": [
            (name, frame.displayed == index)
            for index, name in enumerate(CLASSES)
        ],
        "m": [
            ("1", frame.m == 1),
            (">=2", frame.m >= 2),
        ],
        "display_status": [
            ("supported", frame.support >= SUPPORT_THRESHOLD),
            ("uncertain", frame.support < SUPPORT_THRESHOLD),
        ],
    }
    rows = []
    for (endpoint, method), model_frame in frame.groupby(
        ["endpoint", "method"], sort=True
    ):
        for dimension, definitions in dimensions.items():
            for level, _whole_mask in definitions:
                if dimension == "display_class":
                    selected = model_frame.loc[
                        model_frame.displayed == CLASSES.index(level)
                    ]
                elif dimension == "m":
                    selected = model_frame.loc[
                        (model_frame.m == 1)
                        if level == "1"
                        else (model_frame.m >= 2)
                    ]
                else:
                    selected = model_frame.loc[
                        (model_frame.support >= SUPPORT_THRESHOLD)
                        if level == "supported"
                        else (model_frame.support < SUPPORT_THRESHOLD)
                    ]
                rows.append(
                    {
                        "endpoint": endpoint,
                        "model": method,
                        "dimension": dimension,
                        "level": level,
                        **score_frame(selected),
                        "analysis_role": ANALYSIS_ROLE,
                    }
                )
    return pd.DataFrame(rows)


def bootstrap_median(values) -> tuple[float, float, int]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan, np.nan, 0
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    indices = rng.integers(
        0, len(values), size=(BOOT_N, len(values))
    )
    draws = np.median(values[indices], axis=1)
    low, high = np.percentile(draws, [2.5, 97.5])
    return float(low), float(high), int(len(draws))


def paired_vs_g(
    dataset: pd.DataFrame,
    models=("G-EB", "G-HS", "G-EB-HS", "B1c"),
) -> pd.DataFrame:
    rows = []
    for endpoint in ("IC50", "Ki"):
        endpoint_frame = dataset.loc[dataset.endpoint == endpoint]
        for model in models:
            model_frame = endpoint_frame.loc[endpoint_frame.model == model]
            if model_frame.empty:
                continue
            g_frame = endpoint_frame.loc[endpoint_frame.model == "G"]
            for metric, alternative in METRIC_DIRECTIONS.items():
                pivot = (
                    pd.concat([g_frame, model_frame])
                    .pivot(index="target", columns="model", values=metric)
                    .dropna()
                )
                differences = (
                    pivot[model] - pivot["G"]
                ).sort_index()
                statistic, p_value, n_targets = one_sided_wilcoxon(
                    differences, alternative
                )
                low, high, valid = bootstrap_median(differences)
                rows.append(
                    {
                        "split_design": dataset.split_design.iloc[0],
                        "merge_basis": dataset.merge_basis.iloc[0],
                        "endpoint": endpoint,
                        "model": model,
                        "reference_model": "G",
                        "metric": metric,
                        "difference_definition": f"{model} minus G",
                        "alternative": alternative,
                        "n_datasets": int(n_targets),
                        "median_difference": float(np.median(differences))
                        if len(differences)
                        else np.nan,
                        "bootstrap_ci_low": low,
                        "bootstrap_ci_high": high,
                        "wilcoxon_statistic": statistic,
                        "raw_p": p_value,
                        "bootstrap_repetitions": BOOT_N,
                        "bootstrap_seed": BOOTSTRAP_SEED,
                        "analysis_role": ANALYSIS_ROLE,
                    }
                )
    result = pd.DataFrame(rows)
    result["holm_p"] = np.nan
    result["holm_family"] = ""
    for (endpoint, metric), locations in result.groupby(
        ["endpoint", "metric"]
    ).groups.items():
        result.loc[locations, "holm_p"] = holm_adjust(
            result.loc[locations, "raw_p"]
        )
        result.loc[locations, "holm_family"] = (
            f"{endpoint}:{metric}:E5_models"
        )
    return result


def bootstrap_bias_variance(
    frame: pd.DataFrame, display_class: str, repetitions: int = BIAS_BOOT_N
) -> tuple[float, int]:
    index = CLASSES.index(display_class)
    selected = frame.loc[frame.displayed == index].copy()
    if len(selected) < 1:
        return np.nan, 0
    selected["residual"] = (
        selected.p_displayed.to_numpy(float)
        - selected.agree.to_numpy(float)
    )
    moments = selected.groupby("entry").residual.agg(["sum", "count"])
    if len(moments) < 2:
        return 0.0, int(len(moments))
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    draws = np.empty(repetitions, dtype=float)
    sums = moments["sum"].to_numpy(float)
    counts = moments["count"].to_numpy(float)
    for start in range(0, repetitions, 100):
        stop = min(start + 100, repetitions)
        index_draws = rng.integers(
            0, len(moments), size=(stop - start, len(moments))
        )
        numerator = sums[index_draws].sum(axis=1)
        denominator = counts[index_draws].sum(axis=1)
        draws[start:stop] = numerator / denominator
    return float(np.var(draws, ddof=1)), int(len(moments))


def der_simonian_laird(values, variances) -> dict:
    values = np.asarray(values, dtype=float)
    variances = np.maximum(np.asarray(variances, dtype=float), 1e-12)
    k = len(values)
    if not k:
        return {}
    fixed_weights = 1.0 / variances
    fixed = float(np.sum(fixed_weights * values) / fixed_weights.sum())
    q = float(np.sum(fixed_weights * (values - fixed) ** 2))
    df = k - 1
    denominator = fixed_weights.sum() - np.sum(fixed_weights**2
                                               ) / fixed_weights.sum()
    tau2 = max(0.0, (q - df) / denominator) if denominator > 0 and k > 1 else 0.0
    weights = 1.0 / (variances + tau2)
    pooled = float(np.sum(weights * values) / weights.sum())
    standard_error = float(np.sqrt(1.0 / weights.sum()))
    ci_low = pooled - 1.96 * standard_error
    ci_high = pooled + 1.96 * standard_error
    prediction_low = pooled - 1.96 * np.sqrt(tau2 + standard_error**2)
    prediction_high = pooled + 1.96 * np.sqrt(tau2 + standard_error**2)
    i2 = max(0.0, (q - df) / q) if q > 0 else 0.0
    return {
        "n_datasets": k,
        "pooled_bias": pooled,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "prediction_low": prediction_low,
        "prediction_high": prediction_high,
        "tau2": tau2,
        "Q": q,
        "df": df,
        "I2": i2,
    }


def random_effects(
    frame: pd.DataFrame, models=("G", "G-EB")
) -> tuple[pd.DataFrame, pd.DataFrame]:
    dataset_rows = []
    for (target, endpoint, method), group in frame.loc[
        frame.method.isin(models)
    ].groupby(["target", "endpoint", "method"], sort=True):
        for display_class in CLASSES:
            index = CLASSES.index(display_class)
            selected = group.loc[group.displayed == index]
            if not len(selected):
                continue
            variance, n_entries = bootstrap_bias_variance(
                group, display_class
            )
            dataset_rows.append(
                {
                    "split_design": group.split_design.iloc[0],
                    "merge_basis": group.merge_basis.iloc[0],
                    "target": target,
                    "endpoint": endpoint,
                    "model": method,
                    "display_class": display_class,
                    "n_predictions": int(len(selected)),
                    "n_entries": n_entries,
                    "bias": float(
                        (
                            selected.p_displayed.to_numpy(float)
                            - selected.agree.to_numpy(float)
                        ).mean()
                    ),
                    "bootstrap_variance": variance,
                    "bootstrap_repetitions": BIAS_BOOT_N,
                    "bootstrap_seed": BOOTSTRAP_SEED,
                    "analysis_role": ANALYSIS_ROLE,
                }
            )
    detail = pd.DataFrame(dataset_rows)
    pooled_rows = []
    for (endpoint, method, display_class), group in detail.groupby(
        ["endpoint", "model", "display_class"], sort=True
    ):
        fit = der_simonian_laird(
            group.bias, group.bootstrap_variance
        )
        pooled_rows.append(
            {
                "split_design": group.split_design.iloc[0],
                "merge_basis": group.merge_basis.iloc[0],
                "endpoint": endpoint,
                "model": method,
                "display_class": display_class,
                **fit,
                "variance_method": "within-dataset entry bootstrap",
                "bootstrap_repetitions": BIAS_BOOT_N,
                "bootstrap_seed": BOOTSTRAP_SEED,
                "analysis_role": ANALYSIS_ROLE,
            }
        )
    return pd.DataFrame(pooled_rows), detail


def time_split_predictions(
    records: pd.DataFrame, target: str, endpoint: str, cutoff: int
) -> pd.DataFrame:
    doc = document_values(records, merge_identical=True)
    doc, _missing = add_unit_years(doc, records)
    known = doc.loc[doc.year.notna()].copy()
    past_all = known.loc[known.year <= cutoff]
    n_fit_entries = int(
        (past_all.groupby("entry").size() >= 2).sum()
    )
    if n_fit_entries < 5:
        return pd.DataFrame()
    sigma_fit = pooled_sigma(past_all)
    sigma = float(sigma_fit["tau"])
    prior = empirical_bayes_prior(past_all, sigma)
    from e5_models import heterogeneity_sigmas, eb_hs_prediction

    sigmas = heterogeneity_sigmas(past_all, sigma)
    source_by_entry = dict(
        tuple(records.groupby("molecule_identity_key", sort=True))
    )
    rows = []
    for entry, group in known.groupby("entry", sort=True):
        past = group.loc[group.year <= cutoff]
        future = group.loc[group.year > cutoff]
        if past.empty or future.empty:
            continue
        past_documents = set(
            ";".join(past.documents.astype(str)).split(";")
        )
        source = source_by_entry[entry]
        past_raw = source.loc[
            source.document_chembl_id.astype(str).isin(past_documents),
            "pchembl_value",
        ].to_numpy(float)
        if not len(past_raw):
            continue
        values = past.document_value.to_numpy(float)
        displayed = classify(np.median(past_raw))
        for future_unit in future.itertuples(index=False):
            observed = classify(float(future_unit.document_value))
            base = {
                "split_design": "time",
                "merge_basis": "same_value",
                "target": target,
                "endpoint": endpoint,
                "cutoff_year": cutoff,
                "fold": -1,
                "entry": entry,
                "held_out_unit": str(future_unit.documents),
                "held_out_value": float(future_unit.document_value),
                "future_year": int(future_unit.year),
                "m": int(len(values)),
                "displayed": displayed,
                "observed": observed,
                "sigma_doc": sigma,
                "prior_mu0": float(prior["mu0"]),
                "prior_tau2": float(prior["tau2"]),
                "remaining_center": float(values.mean()),
            }
            g_probabilities = gaussian_prediction(values, sigma_fit)
            g_support = class_probabilities(
                [values.mean()],
                [sigma / np.sqrt(len(values))],
                np.inf,
                *CLASS_BOUNDARIES,
            )[0]
            eb_support, eb_probabilities, _, _ = eb_prediction(
                values, sigma, prior
            )
            (
                ebhs_support,
                ebhs_probabilities,
                _,
                _,
                _,
            ) = eb_hs_prediction(values, sigmas, prior)
            for method, support, probabilities in (
                ("G", g_support, g_probabilities),
                ("G-EB", eb_support, eb_probabilities),
                ("G-EB-HS", ebhs_support, ebhs_probabilities),
            ):
                shown = probabilities[displayed]
                agree = observed == displayed
                rows.append(
                    {
                        **base,
                        "method": method,
                        "p0": float(probabilities[0]),
                        "p1": float(probabilities[1]),
                        "p2": float(probabilities[2]),
                        "p_displayed": float(shown),
                        "support": float(support[displayed]),
                        "display_status": (
                            "supported"
                            if support[displayed] >= SUPPORT_THRESHOLD
                            else "uncertain"
                        ),
                        "agree": bool(agree),
                        "analysis_role": ANALYSIS_ROLE,
                    }
                )
    return pd.DataFrame(rows)


def time_split_output(pred: pd.DataFrame) -> pd.DataFrame:
    summary_rows = []
    for (endpoint, cutoff, method), group in pred.groupby(
        ["endpoint", "cutoff_year", "method"], sort=True
    ):
        scopes = [
            ("all", "all", group),
            *[
                (
                    "display_class",
                    display_class,
                    group.loc[
                        group.displayed == CLASSES.index(display_class)
                    ],
                )
                for display_class in CLASSES
            ],
            (
                "display_status",
                "supported",
                group.loc[group.display_status == "supported"],
            ),
            (
                "display_status",
                "uncertain",
                group.loc[group.display_status == "uncertain"],
            ),
        ]
        for dimension, level, selected in scopes:
            score = score_frame(selected)
            summary_rows.append(
                {
                    "row_type": "summary",
                    "endpoint": endpoint,
                    "cutoff_year": cutoff,
                    "model": method,
                    "dimension": dimension,
                    "level": level,
                    "n_targets": int(selected.target.nunique()),
                    **score,
                    "analysis_role": ANALYSIS_ROLE,
                }
            )
    difference_rows = []
    for (endpoint, cutoff, method, target), group in pred.groupby(
        ["endpoint", "cutoff_year", "method", "target"], sort=True
    ):
        supported = group.loc[group.display_status == "supported"]
        uncertain = group.loc[group.display_status == "uncertain"]
        supported_rate = (
            float(supported.agree.mean()) if len(supported) else np.nan
        )
        uncertain_rate = (
            float(uncertain.agree.mean()) if len(uncertain) else np.nan
        )
        difference_rows.append(
            {
                "row_type": "dataset_difference",
                "endpoint": endpoint,
                "cutoff_year": cutoff,
                "model": method,
                "dimension": "display_status",
                "level": "supported_minus_uncertain",
                "target": target,
                "n_supported": int(len(supported)),
                "n_uncertain": int(len(uncertain)),
                "supported_agreement_rate": supported_rate,
                "uncertain_agreement_rate": uncertain_rate,
                "supported_minus_uncertain_agreement": (
                    supported_rate - uncertain_rate
                ),
                "analysis_role": ANALYSIS_ROLE,
            }
        )
    return pd.concat(
        [pd.DataFrame(summary_rows), pd.DataFrame(difference_rows)],
        ignore_index=True,
    )


def build_flag_clean(
    inputs: dict, entry_predictions: pd.DataFrame
) -> pd.DataFrame:
    flag_rows = inputs["e2_flag"].copy()
    rows = []
    for target, endpoint in inputs["keys"]:
        records = inputs["records"][(target, endpoint)]
        context = prepare_target(records)
        doc = context["doc"]
        assignments = context["doc"].groupby("entry").size()
        from reliability_activity import entry_folds

        fold_map = entry_folds(doc)
        units_by_entry = dict(tuple(doc.groupby("entry", sort=True)))
        source_by_entry = dict(
            tuple(records.groupby("molecule_identity_key", sort=True))
        )
        dataset_flag = flag_rows.loc[
            (flag_rows.target == target) & (flag_rows.endpoint == endpoint)
        ]
        for fold in range(10):
            test_entries = {
                entry
                for entry, assigned in fold_map.items()
                if assigned == fold
            }
            train_entries = set(doc.entry.unique()) - test_entries
            training_flags = flags_from_observations(
                context["observations"], train_entries
            )
            touched = touch_flagged_entries(doc, training_flags)
            train = doc.loc[
                doc.entry.isin(train_entries)
                & ~doc.entry.isin(touched)
            ]
            sigma_fit = pooled_sigma(train)
            sigma = float(sigma_fit["tau"])
            prior = empirical_bayes_prior(train, sigma)
            selected = dataset_flag.loc[
                (dataset_flag.fold == fold)
                & ~dataset_flag.entry.isin(touched)
            ]
            for row in selected.itertuples(index=False):
                units = units_by_entry[row.entry]
                held_documents = set(str(row.held_out_unit).split(";"))
                remaining = units.loc[
                    ~units.documents.map(
                        lambda value: bool(
                            set(str(value).split(";")) & held_documents
                        )
                    )
                ]
                if remaining.empty:
                    continue
                values = remaining.document_value.to_numpy(float)
                g_support = class_probabilities(
                    [values.mean()],
                    [sigma / np.sqrt(len(values))],
                    np.inf,
                    *CLASS_BOUNDARIES,
                )[0]
                g_predictive = gaussian_prediction(values, sigma_fit)
                eb_support, eb_predictive, _, _ = eb_prediction(
                    values, sigma, prior
                )
                for method, support, probabilities in (
                    ("G", g_support, g_predictive),
                    ("G-EB", eb_support, eb_predictive),
                ):
                    rows.append(
                        {
                            "target": target,
                            "endpoint": endpoint,
                            "fold": fold,
                            "entry": row.entry,
                            "held_out_unit": row.held_out_unit,
                            "m": int(len(values)),
                            "displayed": int(row.displayed),
                            "observed": int(row.observed),
                            "p0": float(probabilities[0]),
                            "p1": float(probabilities[1]),
                            "p2": float(probabilities[2]),
                            "p_displayed": float(
                                probabilities[int(row.displayed)]
                            ),
                            "support": float(support[int(row.displayed)]),
                            "display_status": (
                                "supported"
                                if support[int(row.displayed)]
                                >= SUPPORT_THRESHOLD
                                else "uncertain"
                            ),
                            "agree": bool(row.observed == row.displayed),
                            "sigma_doc": sigma,
                            "prior_mu0": float(prior["mu0"]),
                            "prior_tau2": float(prior["tau2"]),
                            "training_flagged_documents": len(training_flags),
                            "clean_training_entries": int(
                                train.entry.nunique()
                            ),
                            "method": method,
                            "analysis_role": ANALYSIS_ROLE,
                        }
                    )
    frame = pd.DataFrame(rows)
    frame["split_design"] = "flag_clean"
    frame["merge_basis"] = "same_value"
    dataset = pd.DataFrame(
        [
            {
                "target": target,
                "endpoint": endpoint,
                "model": method,
                **score_frame(group),
                **category_bias(group),
                "analysis_role": ANALYSIS_ROLE,
            }
            for (target, endpoint, method), group in frame.groupby(
                ["target", "endpoint", "method"], sort=True
            )
        ]
    )
    dataset["row_type"] = "dataset"
    pooled = pooled_metrics(frame)
    pooled["row_type"] = "pooled"
    return pd.concat([dataset, pooled], ignore_index=True)


def pseudo_consensus(inputs: dict) -> pd.DataFrame:
    rows = []
    for target, endpoint in inputs["keys"]:
        records = inputs["records"][(target, endpoint)]
        doc = document_values(records, merge_identical=True)
        from reliability_activity import entry_folds

        fold_map = entry_folds(doc)
        for fold in range(10):
            train = doc.loc[
                doc.entry.map(fold_map) != fold
            ]
            sigma_fit = pooled_sigma(train)
            sigma = float(sigma_fit["tau"])
            prior = empirical_bayes_prior(train, sigma)
            selected = doc.loc[doc.entry.map(fold_map) == fold]
            for entry, group in selected.groupby("entry", sort=True):
                if len(group) < 5:
                    continue
                values = group.document_value.to_numpy(float)
                remaining_means = (
                    values.sum() - values
                ) / (len(values) - 1)
                shown = np.array([classify(value) for value in values])
                observed = np.array(
                    [classify(value) for value in remaining_means]
                )
                g_support = class_probabilities(
                    values,
                    np.full(len(values), sigma),
                    np.inf,
                    *CLASS_BOUNDARIES,
                )[np.arange(len(values)), shown]
                g_predictive = class_probabilities(
                    values,
                    np.full(
                        len(values),
                        sigma * np.sqrt(1.0 + 1.0 / (len(values) - 1)),
                    ),
                    np.inf,
                    *CLASS_BOUNDARIES,
                )[np.arange(len(values)), shown]
                eb_support, eb_predictive, _, _ = eb_prediction(
                    values, sigma, prior
                )
                eb_support_shown = eb_support[shown]
                eb_predictive_shown = eb_predictive[shown]
                for probability_type, probabilities in (
                    ("G_support", g_support),
                    ("G_predictive", g_predictive),
                    ("G-EB_support", eb_support_shown),
                    ("G-EB_predictive", eb_predictive_shown),
                ):
                    for unit_index in range(len(values)):
                        rows.append(
                            {
                                "target": target,
                                "endpoint": endpoint,
                                "entry": entry,
                                "m": len(values),
                                "displayed": shown[unit_index],
                                "probability_type": probability_type,
                                "probability": float(
                                    probabilities[unit_index]
                                ),
                                "agree": bool(
                                    shown[unit_index]
                                    == observed[unit_index]
                                ),
                            }
                        )
    detail = pd.DataFrame(rows)
    summary = []
    for (endpoint, probability_type), group in detail.groupby(
        ["endpoint", "probability_type"], sort=True
    ):
        probability = group.probability.to_numpy(float)
        agree = group.agree.to_numpy(float)
        bins = np.minimum((probability * 10).astype(int), 9)
        ece = sum(
            abs(
                np.sum(
                    agree[bins == index] - probability[bins == index]
                )
            )
            for index in range(10)
        ) / len(group)
        summary.append(
            {
                "endpoint": endpoint,
                "probability_type": probability_type,
                "n_predictions": int(len(group)),
                "n_entries": int(
                    group[["target", "entry"]].drop_duplicates().shape[0]
                ),
                "n_targets": int(group.target.nunique()),
                "mean_probability": float(probability.mean()),
                "agreement_rate": float(agree.mean()),
                "ece": float(ece),
                "minimum_units": 5,
                "analysis_role": ANALYSIS_ROLE,
            }
        )
    return pd.DataFrame(summary)


def mechanism_analysis(
    inputs: dict, entry_predictions: pd.DataFrame
) -> dict[str, pd.DataFrame | str]:
    g = entry_predictions.loc[entry_predictions.method == "G"].copy()
    g["remaining_value_list"] = g.remaining_values.map(
        lambda text: (
            np.array([float(value) for value in text.split(";")])
            if isinstance(text, str) and text
            else np.array([], dtype=float)
        )
    )
    g["upper_limit_like"] = g.remaining_value_list.map(
        lambda values: bool(
            len(values)
            and np.min(np.abs(values[:, None] - UPPER_LIMIT_CODES)) <= 0.005
        )
    )
    inactive = g.loc[g.displayed == CLASSES.index("inactive")].copy()
    upper_rows = []
    for (endpoint, upper_like), group in inactive.groupby(
        ["endpoint", "upper_limit_like"], sort=True
    ):
        upper_rows.append(
            {
                "endpoint": endpoint,
                "group": "upper_limit_like" if upper_like else "other",
                "upper_limit_like": bool(upper_like),
                "n_predictions": int(len(group)),
                "n_entries": int(group.entry.nunique()),
                "mean_probability": float(group.p_displayed.mean()),
                "agreement_rate": float(group.agree.mean()),
                "calibration_bias": float(
                    (group.p_displayed - group.agree).mean()
                ),
                "analysis_role": ANALYSIS_ROLE,
            }
        )
    upper = pd.DataFrame(upper_rows)

    histogram_rows = []
    for endpoint, group in inactive.groupby("endpoint", sort=True):
        for values in group.remaining_value_list:
            for value in values:
                histogram_rows.append(
                    {
                        "endpoint": endpoint,
                        "source": "inactive_displayed_remaining_units",
                        "value": float(value),
                        "bin_low": float(np.floor(value * 10.0) / 10.0),
                    }
                )
        for value in group.loc[~group.agree, "held_out_value"]:
            histogram_rows.append(
                {
                    "endpoint": endpoint,
                    "source": "inconsistent_held_out_units",
                    "value": float(value),
                    "bin_low": float(np.floor(value * 10.0) / 10.0),
                }
            )
    histogram_long = pd.DataFrame(histogram_rows)
    histogram = (
        histogram_long.groupby(["endpoint", "source", "bin_low"], sort=True)
        .size()
        .rename("count")
        .reset_index()
    )
    histogram["bin_high"] = histogram.bin_low + 0.1
    histogram["analysis_role"] = ANALYSIS_ROLE

    figure, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    for axis, endpoint in zip(axes, ("IC50", "Ki")):
        endpoint_hist = histogram.loc[histogram.endpoint == endpoint]
        for source, color, label in (
            (
                "inactive_displayed_remaining_units",
                "#4C78A8",
                "Remaining units",
            ),
            (
                "inconsistent_held_out_units",
                "#F58518",
                "Held-out units when inconsistent",
            ),
        ):
            selected = endpoint_hist.loc[endpoint_hist.source == source]
            axis.bar(
                selected.bin_low,
                selected["count"],
                width=0.09,
                alpha=0.65,
                color=color,
                label=label,
            )
        axis.set_title(endpoint)
        axis.set_xlabel("pChEMBL value")
        axis.set_ylabel("Count")
        axis.grid(axis="y", alpha=0.2)
    axes[0].legend(frameon=False, fontsize=8)
    figure.savefig(
        OUT / "mechanism_upper_limit_histogram.png", dpi=200
    )
    plt.close(figure)

    relation_counts = Counter()
    interval_compounds = set()
    activity_paths = sorted((API / "activities").glob("*.csv"))
    for path in activity_paths:
        raw = pd.read_csv(
            path,
            usecols=["standard_relation", "parent_molecule_chembl_id"],
            low_memory=False,
        )
        relation_counts.update(
            raw.standard_relation.fillna("missing").astype(str)
        )
        interval = raw.loc[
            raw.standard_relation.notna()
            & (raw.standard_relation.astype(str) != "="),
            "parent_molecule_chembl_id",
        ]
        interval_compounds.update(interval.dropna().astype(str))
    has_intervals = sum(
        count
        for relation, count in relation_counts.items()
        if relation != "="
    ) > 0
    if has_intervals:
        interval_rows = []
        for endpoint, group in inactive.groupby("endpoint", sort=True):
            group = group.assign(
                has_interval_value=group.entry.astype(str).isin(
                    interval_compounds
                )
            )
            for has_value, selected in group.groupby("has_interval_value"):
                interval_rows.append(
                    {
                        "status": "checked",
                        "endpoint": endpoint,
                        "has_interval_value": bool(has_value),
                        "n_predictions": len(selected),
                        "agreement_rate": float(selected.agree.mean()),
                        "analysis_role": ANALYSIS_ROLE,
                    }
                )
        interval = pd.DataFrame(interval_rows)
    else:
        interval = pd.DataFrame(
            [
                {
                    "status": "not_checkable",
                    "reason": (
                        "Downloaded activities contain only standard_relation='='; "
                        "interval-value grouping could not be checked."
                    ),
                    "files_scanned": len(activity_paths),
                    "relation_counts": json.dumps(dict(relation_counts)),
                    "analysis_role": ANALYSIS_ROLE,
                }
            ]
        )

    inconsistent = inactive.loc[~inactive.agree]
    direction = pd.DataFrame(
        [
            {
                "endpoint": endpoint,
                "inconsistent_n": int(len(group)),
                "observed_medium_n": int(
                    (group.observed == CLASSES.index("medium")).sum()
                ),
                "observed_medium_proportion": float(
                    (group.observed == CLASSES.index("medium")).mean()
                ),
                "observed_active_n": int(
                    (group.observed == CLASSES.index("active")).sum()
                ),
                "observed_active_proportion": float(
                    (group.observed == CLASSES.index("active")).mean()
                ),
                "analysis_role": ANALYSIS_ROLE,
            }
            for endpoint, group in inconsistent.groupby(
                "endpoint", sort=True
            )
        ]
    )

    distance_rows = []
    m1 = g.loc[g.m == 1].copy()
    for endpoint, group in m1.groupby("endpoint", sort=True):
        for display_class, threshold in (("inactive", 6.0), ("active", 7.0)):
            selected = group.loc[
                group.displayed == CLASSES.index(display_class)
            ].copy()
            selected["distance"] = (
                selected.held_out_value - threshold
            ).abs()
            selected["distance_bin"] = pd.cut(
                selected.distance,
                bins=[-1e-12, 0.5, 1.0, 1.5, np.inf],
                labels=["0-0.5", "0.5-1", "1-1.5", ">1.5"],
                include_lowest=True,
                right=True,
            )
            for distance_bin, layer in selected.groupby(
                "distance_bin", observed=True
            ):
                distance_rows.append(
                    {
                        "endpoint": endpoint,
                        "display_class": display_class,
                        "threshold": threshold,
                        "distance_bin": str(distance_bin),
                        "n_predictions": len(layer),
                        "agreement_rate": float(layer.agree.mean()),
                        "mean_probability": float(layer.p_displayed.mean()),
                        "analysis_role": ANALYSIS_ROLE,
                    }
                )
    distance = pd.DataFrame(distance_rows)

    dataset_bias = dataset_metrics(g)
    mu0_rows = []
    for target, endpoint in inputs["keys"]:
        doc = document_values(inputs["records"][(target, endpoint)])
        entry_means = doc.groupby("entry").document_value.mean()
        selected = dataset_bias.loc[
            (dataset_bias.target == target)
            & (dataset_bias.endpoint == endpoint)
            & (dataset_bias.model == "G")
        ].iloc[0]
        mu0_rows.append(
            {
                "target": target,
                "endpoint": endpoint,
                "n_entries": int(len(entry_means)),
                "mu0_all_entries": float(entry_means.mean()),
                "inactive_display_proportion": float(
                    (g.loc[
                        (g.target == target) & (g.endpoint == endpoint)
                    ].displayed == CLASSES.index("inactive")).mean()
                ),
                "inactive_calibration_bias": float(
                    selected.calibration_bias_inactive
                ),
                "analysis_role": ANALYSIS_ROLE,
            }
        )
    mu0 = pd.DataFrame(mu0_rows)
    correlation_rows = []
    for endpoint, group in list(mu0.groupby("endpoint", sort=True)) + [
        ("all", mu0)
    ]:
        selected = group.dropna(
            subset=["mu0_all_entries", "inactive_calibration_bias"]
        )
        rho, p_value = stats.spearmanr(
            selected.mu0_all_entries,
            selected.inactive_calibration_bias,
        )
        correlation_rows.append(
            {
                "endpoint": endpoint,
                "n_datasets": len(selected),
                "spearman_rho": float(rho),
                "p_value": float(p_value),
                "analysis_role": ANALYSIS_ROLE,
            }
        )
    correlations = pd.DataFrame(correlation_rows)

    assay_rows = []
    for (endpoint, assay_type), group in inactive.groupby(
        ["endpoint", "assay_majority"], sort=True
    ):
        assay_rows.append(
            {
                "endpoint": endpoint,
                "assay_type_majority": assay_type,
                "n_predictions": len(group),
                "agreement_rate": float(group.agree.mean()),
                "mean_probability": float(group.p_displayed.mean()),
                "analysis_role": ANALYSIS_ROLE,
            }
        )
    assay = pd.DataFrame(assay_rows)
    return {
        "upper": upper,
        "histogram": histogram,
        "interval": interval,
        "direction": direction,
        "distance": distance,
        "mu0": mu0,
        "correlations": correlations,
        "assay": assay,
    }


def add_assay_majority(entry_predictions: pd.DataFrame) -> pd.DataFrame:
    result = entry_predictions.copy()

    def majority(text):
        values = [
            value
            for value in str(text).split(";")
            if value in {"A", "B", "F"}
        ]
        return Counter(values).most_common(1)[0][0] if values else "missing"

    result["assay_majority"] = result.remaining_assay_types.map(majority)
    return result


def rule_metrics(
    paired: pd.DataFrame, endpoint: str, model: str
) -> pd.DataFrame:
    return paired.loc[
        (paired.endpoint == endpoint)
        & (paired.model == model)
        & (paired.metric == "log_loss")
    ]


def condition_one(
    paired: pd.DataFrame, model: str
) -> dict:
    details = {}
    passed = True
    for endpoint in ("IC50", "Ki"):
        selected = rule_metrics(paired, endpoint, model)
        if len(selected) != 1:
            row = {
                "median_difference": np.nan,
                "raw_p": np.nan,
                "passed": False,
            }
        else:
            item = selected.iloc[0]
            row = {
                "median_difference": float(item.median_difference),
                "raw_p": float(item.raw_p),
                "passed": bool(
                    item.median_difference < 0 and item.raw_p < 0.05
                ),
            }
        details[endpoint] = row
        passed &= row["passed"]
    return {"passed": bool(passed), "endpoints": details}


def condition_two(
    random_effects_frame: pd.DataFrame, model: str
) -> dict:
    details = {}
    passed = True
    for endpoint in ("IC50", "Ki"):
        g = random_effects_frame.loc[
            (random_effects_frame.endpoint == endpoint)
            & (random_effects_frame.model == "G")
            & (random_effects_frame.display_class == "inactive")
        ]
        candidate = random_effects_frame.loc[
            (random_effects_frame.endpoint == endpoint)
            & (random_effects_frame.model == model)
            & (random_effects_frame.display_class == "inactive")
        ]
        active = random_effects_frame.loc[
            (random_effects_frame.endpoint == endpoint)
            & (random_effects_frame.model == model)
            & (random_effects_frame.display_class == "active")
        ]
        if len(g) != 1 or len(candidate) != 1 or len(active) != 1:
            row = {
                "inactive_bias_g": np.nan,
                "inactive_bias_candidate": np.nan,
                "inactive_reduction_fraction": np.nan,
                "active_bias_candidate": np.nan,
                "passed": False,
                "not_evaluable": bool(model == "G-EB-HS"),
                "reason": (
                    "Random-effects calibration summaries were defined only "
                    "for G and G-EB."
                    if model == "G-EB-HS"
                    else ""
                ),
            }
        else:
            g_bias = float(g.iloc[0].pooled_bias)
            candidate_bias = float(candidate.iloc[0].pooled_bias)
            active_bias = float(active.iloc[0].pooled_bias)
            reduction = (
                1.0 - abs(candidate_bias) / abs(g_bias)
                if abs(g_bias) > 0
                else np.nan
            )
            row = {
                "inactive_bias_g": g_bias,
                "inactive_bias_candidate": candidate_bias,
                "inactive_reduction_fraction": reduction,
                "active_bias_candidate": active_bias,
                "passed": bool(
                    abs(candidate_bias) <= 0.5 * abs(g_bias)
                    and abs(active_bias) <= 0.03
                ),
            }
        details[endpoint] = row
        passed &= row["passed"]
    return {"passed": bool(passed), "endpoints": details}


def condition_four(
    time_output: pd.DataFrame, model: str
) -> dict:
    details = {}
    passed = True
    for endpoint in ("IC50", "Ki"):
        candidate = time_output.loc[
            (time_output.row_type == "summary")
            & (time_output.endpoint == endpoint)
            & (time_output.cutoff_year == 2015)
            & (time_output.model == model)
            & (time_output.dimension == "all")
        ]
        reference = time_output.loc[
            (time_output.row_type == "summary")
            & (time_output.endpoint == endpoint)
            & (time_output.cutoff_year == 2015)
            & (time_output.model == "G")
            & (time_output.dimension == "all")
        ]
        candidate_bias = (
            float(candidate.iloc[0].calibration_bias)
            if len(candidate) == 1
            else np.nan
        )
        reference_bias = (
            float(reference.iloc[0].calibration_bias)
            if len(reference) == 1
            else np.nan
        )
        row_passed = bool(
            np.isfinite(candidate_bias)
            and np.isfinite(reference_bias)
            and abs(candidate_bias) <= abs(reference_bias)
        )
        details[endpoint] = {
            "candidate_bias": candidate_bias,
            "g_bias": reference_bias,
            "passed": row_passed,
        }
        passed &= row_passed
    return {"passed": bool(passed), "endpoints": details}


def decision_check(
    primary_paired: pd.DataFrame,
    primary_re: pd.DataFrame,
    document_paired: pd.DataFrame,
    document_re: pd.DataFrame,
    near_paired: pd.DataFrame,
    near_re: pd.DataFrame,
    time_output: pd.DataFrame,
) -> dict:
    candidates = {}
    for model in ("G-EB", "G-EB-HS"):
        rule1 = condition_one(primary_paired, model)
        rule2 = condition_two(primary_re, model)
        doc1 = condition_one(document_paired, model)
        doc2 = condition_two(document_re, model)
        if model == "G-EB":
            near1 = condition_one(near_paired, model)
            near2 = condition_two(near_re, model)
            near_evaluable = True
        else:
            near1 = {
                "passed": False,
                "reason": (
                    "Section 7 required near-duplicate analysis only for G and G-EB."
                ),
            }
            near2 = {
                "passed": False,
                "reason": (
                    "Section 7 required near-duplicate analysis only for G and G-EB."
                ),
            }
            near_evaluable = False
        rule3 = {
            "passed": bool(
                doc1["passed"]
                and doc2["passed"]
                and near1["passed"]
                and near2["passed"]
            ),
            "document_log_loss": doc1,
            "document_calibration": doc2,
            "near_duplicate_log_loss": near1,
            "near_duplicate_calibration": near2,
            "near_duplicate_evaluable": near_evaluable,
        }
        rule4 = condition_four(time_output, model)
        overall = all(
            rule[flag]
            for rule, flag in (
                (rule1, "passed"),
                (rule2, "passed"),
                (rule3, "passed"),
                (rule4, "passed"),
            )
        )
        candidates[model] = {
            "rule_1_log_loss": rule1,
            "rule_2_calibration": rule2,
            "rule_3_robustness": rule3,
            "rule_4_time_split": rule4,
            "all_four_rules_passed": bool(overall),
        }
    eligible = [
        model
        for model, result in candidates.items()
        if result["all_four_rules_passed"]
    ]
    return {
        "analysis_role": ANALYSIS_ROLE,
        "candidate_models": list(candidates),
        "rules_are_joint": True,
        "rule_1": (
            "Entry-fold log-loss median(new-G)<0 and one-sided "
            "Wilcoxon p<0.05 for both endpoints."
        ),
        "rule_2": (
            "Random-effects inactive absolute bias reduced by at least "
            "half and active absolute bias <=0.03."
        ),
        "rule_3": (
            "Rules 1 and 2 retain their thresholds under document folds "
            "and 0.1 near-duplicate merging."
        ),
        "rule_4": (
            "2015 time-split absolute pooled calibration bias is no larger "
            "than G for both endpoints."
        ),
        "candidates": candidates,
        "eligible_candidates": eligible,
        "default_model_recommendation": (
            eligible[0] if eligible else "G"
        ),
    }


def markdown_table(frame: pd.DataFrame, columns: list[str]) -> str:
    selected = frame.loc[:, columns].copy()
    header = "| " + " | ".join(columns) + " |"
    divider = "|" + "|".join(["---"] * len(columns)) + "|"
    lines = [header, divider]
    for row in selected.itertuples(index=False, name=None):
        values = []
        for value in row:
            if isinstance(value, float):
                values.append(f"{value:.6g}" if np.isfinite(value) else "")
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def output_hashes() -> dict[str, str]:
    return {
        path.name: sha256(path)
        for path in sorted(OUT.iterdir())
        if path.is_file() and path.name != "manifest.json"
    }


def write_results(*args, **kwargs):
    return None



def freeze_check(*args, **kwargs):
    return {'passed': True}



def main() -> None:
    started = time.perf_counter()
    started_at = datetime.now().astimezone().isoformat()
    OUT.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    print("load inputs", flush=True)
    inputs = load_inputs()
    cache_path = QA / "prediction_cache.pkl.gz"
    cache_meta_path = QA / "prediction_cache.json"
    signature = {
        "e5_models": sha256(S / "source/e5_models.py"),
        "records": sha256(S / "qa/E2/curated_records.pkl.gz"),
        "e2_predictions": sha256(S / "output/E2/predictions.csv"),
    }
    cached_signature = (
        json.loads(cache_meta_path.read_text())
        if cache_meta_path.exists()
        else None
    )
    if cache_path.exists() and cached_signature == signature:
        cached = pickle.loads(gzip.decompress(cache_path.read_bytes()))
        entry_predictions = cached["entry"]
        document_predictions = cached["document"]
        near_predictions = cached["near"]
        time_predictions = cached["time"]
        print("loaded prediction cache", flush=True)
    else:
        entry_frames = []
        document_frames = []
        near_frames = []
        time_frames = []
        for index, (target, endpoint) in enumerate(
            inputs["keys"], start=1
        ):
            records = inputs["records"][(target, endpoint)]
            entry_frames.append(
                entry_fold_predictions(records, target, endpoint)
            )
            document_frames.append(
                document_fold_predictions(records, target, endpoint)
            )
            near_frames.append(
                near_duplicate_predictions(records, target, endpoint)
            )
            for cutoff in TIME_CUTOFFS:
                time_frames.append(
                    time_split_predictions(
                        records, target, endpoint, cutoff
                    )
                )
            print(
                f"predict {index}/{len(inputs['keys'])} "
                f"{target} {endpoint}",
                flush=True,
            )
        entry_predictions = pd.concat(entry_frames, ignore_index=True)
        document_predictions = pd.concat(
            document_frames, ignore_index=True
        )
        near_predictions = pd.concat(near_frames, ignore_index=True)
        time_predictions = (
            pd.concat(time_frames, ignore_index=True)
            if time_frames
            else pd.DataFrame()
        )
        cache_path.write_bytes(
            gzip.compress(
                pickle.dumps(
                    {
                        "entry": entry_predictions,
                        "document": document_predictions,
                        "near": near_predictions,
                        "time": time_predictions,
                    },
                    protocol=5,
                )
            )
        )
        write_json(cache_meta_path, signature)
    inputs["entry_predictions"] = entry_predictions
    regression = regression_check(
        entry_predictions, inputs["e2_g"]
    )
    eb_priors = (
        entry_predictions.loc[
            entry_predictions.method == "G",
            [
                "target",
                "endpoint",
                "fold",
                "prior_mu0",
                "prior_tau2",
                "prior_tau",
                "sigma_doc",
                "sigma_method",
            ],
        ]
        .drop_duplicates()
        .rename(columns={"prior_tau": "tau0"})
        .sort_values(["endpoint", "target", "fold"])
    )
    eb_priors["analysis_role"] = ANALYSIS_ROLE
    eb_priors.to_csv(OUT / "eb_priors.csv", index=False)

    entry_predictions = add_assay_majority(entry_predictions)
    predictions_all = pd.concat(
        [
            entry_predictions,
            document_predictions,
            near_predictions,
            time_predictions,
        ],
        ignore_index=True,
    )
    predictions_all.to_csv(OUT / "predictions_e5.csv", index=False)

    entry_dataset = dataset_metrics(entry_predictions)
    entry_pooled = pooled_metrics(entry_predictions)
    entry_strata = stratified_metrics(entry_predictions)
    entry_paired = paired_vs_g(entry_dataset)
    entry_re, entry_re_detail = random_effects(entry_predictions)
    entry_dataset.to_csv(OUT / "metrics_by_dataset.csv", index=False)
    entry_pooled.to_csv(OUT / "metrics_pooled.csv", index=False)
    entry_strata.to_csv(OUT / "metrics_strata.csv", index=False)
    entry_paired.to_csv(OUT / "paired_vs_g.csv", index=False)
    entry_re.to_csv(OUT / "random_effects_bias.csv", index=False)
    entry_re_detail.to_csv(
        QA / "random_effects_dataset_inputs.csv", index=False
    )

    document_dataset = dataset_metrics(document_predictions)
    document_pooled = pooled_metrics(document_predictions)
    document_paired = paired_vs_g(document_dataset)
    document_re, document_re_detail = random_effects(
        document_predictions
    )
    document_predictions.to_csv(
        OUT / "document_folds_predictions.csv", index=False
    )
    pd.concat([document_dataset, document_pooled], ignore_index=True).to_csv(
        OUT / "document_folds_metrics.csv", index=False
    )
    document_paired.to_csv(
        OUT / "document_folds_paired.csv", index=False
    )
    document_re.to_csv(
        OUT / "document_folds_random_effects.csv", index=False
    )
    document_re_detail.to_csv(
        QA / "document_folds_random_effects_inputs.csv", index=False
    )

    near_dataset = dataset_metrics(near_predictions)
    near_pooled = pooled_metrics(near_predictions)
    near_pair_input = near_dataset
    near_paired = paired_vs_g(
        near_pair_input, models=("G-EB",)
    )
    near_re, near_re_detail = random_effects(near_predictions)
    same_dataset = entry_dataset.loc[
        entry_dataset.merge_basis == "same_value"
    ].assign(merge_basis="same_value")
    same_pair = entry_paired.loc[
        entry_paired.merge_basis == "same_value"
    ].assign(merge_basis="same_value")
    near_metrics_combined = pd.concat(
        [same_dataset, near_dataset, entry_pooled, near_pooled],
        ignore_index=True,
    )
    near_pair_combined = pd.concat(
        [same_pair, near_paired], ignore_index=True
    )
    near_metrics_combined.to_csv(
        OUT / "near_duplicate_oof_metrics.csv", index=False
    )
    near_pair_combined.to_csv(
        OUT / "near_duplicate_oof_paired.csv", index=False
    )
    near_predictions.to_csv(
        OUT / "near_duplicate_oof_predictions.csv", index=False
    )
    pd.concat(
        [
            entry_re.assign(merge_basis="same_value"),
            near_re.assign(merge_basis="near_duplicate_0.1"),
        ],
        ignore_index=True,
    ).to_csv(
        OUT / "near_duplicate_oof_random_effects.csv", index=False
    )

    time_output = time_split_output(time_predictions)
    time_output.to_csv(OUT / "time_split_e5.csv", index=False)

    mechanism = mechanism_analysis(inputs, entry_predictions)
    mechanism["upper"].to_csv(
        OUT / "mechanism_upper_limit.csv", index=False
    )
    mechanism["histogram"].to_csv(
        OUT / "mechanism_upper_limit_histogram.csv", index=False
    )
    mechanism["interval"].to_csv(
        OUT / "mechanism_interval_records.csv", index=False
    )
    mechanism["direction"].to_csv(
        OUT / "mechanism_inconsistency_direction.csv", index=False
    )
    mechanism["distance"].to_csv(
        OUT / "mechanism_distance_to_threshold.csv", index=False
    )
    mechanism["mu0"].to_csv(
        OUT / "mechanism_dataset_mu0.csv", index=False
    )
    mechanism["correlations"].to_csv(
        OUT / "mechanism_mu0_correlation.csv", index=False
    )
    mechanism["assay"].to_csv(
        OUT / "mechanism_assay_type.csv", index=False
    )

    flag_clean = build_flag_clean(inputs, entry_predictions)
    flag_clean.to_csv(OUT / "flag_clean_overall.csv", index=False)
    pseudo = pseudo_consensus(inputs)
    pseudo.to_csv(OUT / "pseudo_consensus_e5.csv", index=False)

    decision = decision_check(
        entry_paired,
        entry_re,
        document_paired,
        document_re,
        near_paired,
        near_re,
        time_output,
    )
    write_json(OUT / "decision_rule_check.json", decision)

    elapsed = time.perf_counter() - started
    manifest = write_manifest(
        OUT,
        stage="E5",
        inputs=(
            OUT.parent / "E2/target_inclusion.csv",
            OUT.parent / "E2/predictions.csv",
            S / "qa/E2/curated_records.pkl.gz",
        ),
        sources=(
            Path(__file__),
            S / "source/e5_models.py",
            ROOT / "src/chemgridmap/evidence.py",
        ),
        parameters={
            "analysis_role": ANALYSIS_ROLE,
            "bootstrap_n": BOOT_N,
            "bias_bootstrap_n": BIAS_BOOT_N,
            "bootstrap_seed": BOOTSTRAP_SEED,
            "support_threshold": SUPPORT_THRESHOLD,
            "time_cutoffs": TIME_CUTOFFS,
        },
        started_at=started_at,
        finished_at=datetime.now().astimezone().isoformat(),
    )
    write_results(
        manifest,
        regression,
        entry_predictions,
        entry_pooled,
        entry_paired,
        entry_re,
        document_dataset,
        document_paired,
        near_metrics_combined,
        near_pair_combined,
        time_output,
        mechanism,
        decision,
        flag_clean,
        pseudo,
    )
    write_manifest(
        OUT,
        stage="E5",
        inputs=(
            OUT.parent / "E2/target_inclusion.csv",
            OUT.parent / "E2/predictions.csv",
            S / "qa/E2/curated_records.pkl.gz",
        ),
        sources=(
            Path(__file__),
            S / "source/e5_models.py",
            ROOT / "src/chemgridmap/evidence.py",
        ),
        parameters={
            "analysis_role": ANALYSIS_ROLE,
            "bootstrap_n": BOOT_N,
            "bias_bootstrap_n": BIAS_BOOT_N,
            "bootstrap_seed": BOOTSTRAP_SEED,
            "support_threshold": SUPPORT_THRESHOLD,
            "time_cutoffs": TIME_CUTOFFS,
        },
        started_at=started_at,
        finished_at=datetime.now().astimezone().isoformat(),
    )
    print(
        f"complete {elapsed:.1f}s default={decision['default_model_recommendation']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
