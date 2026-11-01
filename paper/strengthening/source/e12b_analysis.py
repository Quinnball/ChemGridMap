"""E12b: post hoc scores, tests, core settings, and data descriptions."""
from __future__ import annotations

from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import pickle
import platform
import sys
import time

sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parents[3]
S = ROOT / "paper/strengthening"
OUT = S / "output/E12b"
QA = S / "qa/E12b"
CGM = ROOT
FROZEN_PACKAGE = CGM
E2 = S / "output/E2"
E2_QA = S / "qa/E2"
E4 = S / "output/E4"
E6 = S / "output/E6"
E7 = S / "output/E7"
E9 = S / "output/E9"
E1B = S / "output/E1B"
FIGURE_DATA = ROOT / "paper/figures/data_final"
API = S / "data/chembl37_api"

sys.path[:0] = [str(S / "source"), str(CGM / "src"), str(CGM / "paper")]
import chemgridmap

if not Path(chemgridmap.__file__).resolve().is_relative_to(FROZEN_PACKAGE):
    raise RuntimeError(f"chemgridmap imported outside frozen copy: {chemgridmap.__file__}")

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy import stats

from chemgridmap.evidence import document_values
from e1a_statistics import percentile_ci, weighted_scores
from e2b_analysis import add_unit_years, load_original_ic50_records
from e4b_analysis import build_temporal_frame
from e7_analysis import (
    holm_adjust,
    metric_draws,
    one_sided_wilcoxon,
    point_metric_values,
    two_stage_weights,
)

ENDPOINTS = ("IC50", "Ki")
BOOT_N = 2000
BOOT_SEED = 20261103
SIMULATION_SEED = 20261113
NEAR_DUPLICATE_GAP = 0.1
SUPPORT_THRESHOLD = 0.8
DEVELOPMENT_TARGETS = {"CHEMBL204", "CHEMBL205", "CHEMBL240"}
ER_TARGET = "CHEMBL206"
ATR_TARGET = "CHEMBL5024"
SIZE_ORDER = ("2", "3-4", "5-9", ">=10")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(frame: pd.DataFrame, name: str, *, qa: bool = False) -> pd.DataFrame:
    output = frame.copy()
    output["analysis_role"] = "post hoc"
    path = (QA if qa else OUT) / name
    output.to_csv(path, index=False, float_format="%.17g", na_rep="NA")
    return output


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


def read_predictions() -> pd.DataFrame:
    return pd.read_csv(E2 / "predictions.csv", low_memory=False, float_precision="round_trip")


def load_records() -> dict[tuple[str, str], pd.DataFrame]:
    records = pickle.loads(gzip.decompress((E2_QA / "curated_records.pkl.gz").read_bytes()))
    records.update(load_original_ic50_records())
    return records


def two_sided_wilcoxon(values) -> tuple[float, float, int]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values) or np.allclose(values, 0):
        return float("nan"), 1.0, len(values)
    result = stats.wilcoxon(values, alternative="two-sided", zero_method="wilcox")
    return float(result.statistic), float(result.pvalue), len(values)


def population_holm(frame: pd.DataFrame, family: str, family_column="correction_family") -> pd.DataFrame:
    result = frame.copy()
    result["holm_p_two_sided"] = np.nan
    locations = result.index[result[family_column] == family]
    result.loc[locations, "holm_p_two_sided"] = holm_adjust(
        result.loc[locations, "raw_p_two_sided"]
    )
    return result


def score_probabilities(frame: pd.DataFrame, score: str) -> np.ndarray:
    probabilities = frame[["p0", "p1", "p2"]].to_numpy(float)
    observed = frame.observed.to_numpy(int)
    one_hot = np.eye(3)[observed]
    if score == "brier":
        return ((probabilities - one_hot) ** 2).sum(axis=1)
    clip = 1e-6 if score == "log_loss_1e-6" else 1e-15
    probability = probabilities[np.arange(len(frame)), observed]
    return -np.log(np.clip(probability, clip, 1.0))


def dataset_scores(predictions: pd.DataFrame, score: str) -> pd.DataFrame:
    selected = predictions.loc[predictions.population == "all"].copy()
    selected["score_value"] = score_probabilities(selected, score)
    return (
        selected.groupby(["target", "endpoint", "method"], sort=True)
        .score_value.mean()
        .unstack("method")
    )


def build_scores_h4_h7(predictions: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    source = pd.read_csv(E2 / "hypothesis_tests.csv")
    flags = pd.read_csv(E2 / "flagged_pairs_census.csv")
    flagged_keys = set(
        flags.loc[
            (flags.row_type == "target_summary") & (flags.n_flagged_pairs > 0),
            ["target", "endpoint"],
        ].itertuples(index=False, name=None)
    )
    definitions = {
        "H4": ("G_minus_G-record_log_loss", "G", "G-record", -1),
        "H5": ("G_minus_G-global_log_loss", "G", "G-global", -1),
        "H6": ("G_minus_G-fixed_log_loss", "G", "G-fixed", -1),
        "H7": ("X_minus_G_log_loss", "X", "G", -1),
    }
    rows = []
    regression = {}
    for endpoint in ENDPOINTS:
        endpoint_source = source.loc[
            (source.endpoint == endpoint) & (source.population == "all")
        ].set_index("hypothesis")
        h3_p = float(endpoint_source.loc["H3", "raw_p"])
        for score in ("log_loss_1e-6", "log_loss_1e-15", "brier"):
            wide = dataset_scores(predictions, score)
            for hypothesis, (definition, left, right, favorable_sign) in definitions.items():
                values = (wide[left] - wide[right]).copy()
                values = values.loc[values.index.get_level_values("endpoint") == endpoint]
                if hypothesis == "H7":
                    values = values.loc[
                        values.index.isin(
                            [key for key in values.index if key in flagged_keys]
                        )
                    ]
                values = values.dropna()
                statistic, one_p, count = one_sided_wilcoxon(values, "less")
                two_stat, two_p, _ = two_sided_wilcoxon(values)
                rows.append({
                    "score": score,
                    "endpoint": endpoint,
                    "hypothesis": hypothesis,
                    "difference_definition": definition,
                    "n_datasets": count,
                    "median_difference": float(values.median()) if count else np.nan,
                    "mean_difference": float(values.mean()) if count else np.nan,
                    "g_advantaged_datasets": int((favorable_sign * values > 0).sum()),
                    "original_direction": "less",
                    "one_sided_statistic": statistic,
                    "one_sided_p": one_p,
                    "holm_p_h3_h7": np.nan,
                    "two_sided_statistic": two_stat,
                    "raw_p_two_sided": two_p,
                    "holm_p_h4_h7": np.nan,
                    "correction_family_predeclared": "H3-H7",
                    "correction_family_two_sided": "H4-H7",
                })
            locations = [
                index for index, row in enumerate(rows)
                if row["score"] == score and row["endpoint"] == endpoint
            ]
            pvalues = np.r_[
                h3_p,
                [rows[index]["one_sided_p"] for index in locations],
            ]
            adjusted = holm_adjust(pvalues)
            for index, value in zip(locations, adjusted[1:]):
                rows[index]["holm_p_h3_h7"] = value
            two_adjusted = holm_adjust([rows[index]["raw_p_two_sided"] for index in locations])
            for index, value in zip(locations, two_adjusted):
                rows[index]["holm_p_h4_h7"] = value

    result = pd.DataFrame(rows)
    original = result.loc[
        result.score == "log_loss_1e-6"
    ].set_index(["endpoint", "hypothesis"])
    for endpoint in ENDPOINTS:
        for hypothesis in definitions:
            expected = source.loc[
                (source.endpoint == endpoint)
                & (source.population == "all")
                & (source.hypothesis == hypothesis)
            ].iloc[0]
            actual = original.loc[(endpoint, hypothesis)]
            regression[f"{endpoint}:{hypothesis}"] = {
                "median_absolute_error": abs(actual.median_difference - expected.median_difference),
                "one_sided_p_absolute_error": abs(actual.one_sided_p - expected.raw_p),
                "holm_absolute_error": abs(actual.holm_p_h3_h7 - expected.holm_p),
            }
    maximum = {
        key: max(item[key] for item in regression.values())
        for key in ("median_absolute_error", "one_sided_p_absolute_error", "holm_absolute_error")
    }
    if maximum["median_absolute_error"] > 1e-12 or maximum["one_sided_p_absolute_error"] > 1e-12 or maximum["holm_absolute_error"] > 1e-12:
        raise RuntimeError(f"H4-H7 regression failed: {maximum}")
    return result, {"rows": regression, "maximum_errors": maximum, "passed": True}


def dl_details(values, variances) -> dict:
    values = np.asarray(values, float)
    variances = np.maximum(np.asarray(variances, float), 1e-12)
    fixed_weights = 1.0 / variances
    fixed = float(np.sum(fixed_weights * values) / fixed_weights.sum())
    q = float(np.sum(fixed_weights * (values - fixed) ** 2))
    c = float(fixed_weights.sum() - np.sum(fixed_weights**2) / fixed_weights.sum())
    tau2 = max(0.0, (q - (len(values) - 1)) / c) if c > 0 else 0.0
    weights = 1.0 / (variances + tau2)
    estimate = float(np.sum(weights * values) / weights.sum())
    se = float(np.sqrt(1.0 / weights.sum()))
    z = float(stats.norm.ppf(0.975))
    return {
        "n_datasets": len(values),
        "estimate": estimate,
        "ci_low": estimate - z * se,
        "ci_high": estimate + z * se,
        "tau2": tau2,
        "q": q,
        "i_squared": max(0.0, (q - (len(values) - 1)) / q * 100.0) if q > 0 else 0.0,
        "max_weight": float(weights.max() / weights.sum()),
        "max_weight_target": "",
        "weights": weights,
    }


def build_random_effect_alternatives() -> tuple[pd.DataFrame, dict]:
    source = pd.read_csv(E7 / "random_effects_dataset_inputs.csv")
    reference = pd.read_csv(E7 / "random_effects.csv")
    rows = []
    for endpoint in ENDPOINTS:
        sensitive = ATR_TARGET if endpoint == "Ki" else ER_TARGET
        for hypothesis in ("H4", "H5", "H6"):
            selected = source.loc[
                (source.endpoint == endpoint)
                & (source.hypothesis == hypothesis)
            ].sort_values("target").reset_index(drop=True)
            values = selected.point_estimate.to_numpy(float)
            variances = selected.bootstrap_variance.to_numpy(float)
            targets = selected.target.to_numpy(str)

            dl = dl_details(values, variances)
            dl["max_weight_target"] = targets[int(np.argmax(dl["weights"]))]
            raw_weights = 1.0 / np.maximum(variances, 1e-12)
            fixed = float(np.sum(raw_weights * values) / raw_weights.sum())
            fixed_se = float(np.sqrt(1.0 / raw_weights.sum()))
            hksj_se = float(
                np.sqrt(
                    np.sum(dl["weights"] * (values - dl["estimate"]) ** 2)
                    / ((len(values) - 1) * dl["weights"].sum())
                )
            )
            hksj_half = float(stats.t.ppf(0.975, len(values) - 1) * hksj_se)
            mean_t_half = float(
                stats.t.ppf(0.975, len(values) - 1)
                * values.std(ddof=1)
                / np.sqrt(len(values))
            )
            keep = targets != sensitive
            excluded_dl = dl_details(values[keep], variances[keep])
            excluded_dl["max_weight_target"] = targets[keep][
                int(np.argmax(excluded_dl["weights"]))
            ]
            methods = [
                {
                    "method": "DerSimonian-Laird",
                    "excluded_target": "",
                    **{key: dl[key] for key in ("estimate", "ci_low", "ci_high", "tau2", "q", "i_squared", "max_weight", "max_weight_target")},
                },
                {
                    "method": "fixed-effect",
                    "excluded_target": "",
                    "estimate": fixed,
                    "ci_low": fixed - 1.96 * fixed_se,
                    "ci_high": fixed + 1.96 * fixed_se,
                    "tau2": 0.0,
                    "q": dl["q"],
                    "i_squared": dl["i_squared"],
                    "max_weight": float(raw_weights.max() / raw_weights.sum()),
                    "max_weight_target": targets[int(np.argmax(raw_weights))],
                },
                {
                    "method": "DerSimonian-Laird-HKSJ",
                    "excluded_target": "",
                    "estimate": dl["estimate"],
                    "ci_low": dl["estimate"] - hksj_half,
                    "ci_high": dl["estimate"] + hksj_half,
                    "tau2": dl["tau2"],
                    "q": dl["q"],
                    "i_squared": dl["i_squared"],
                    "max_weight": dl["max_weight"],
                    "max_weight_target": dl["max_weight_target"],
                },
                {
                    "method": "simple-mean-t",
                    "excluded_target": "",
                    "estimate": float(values.mean()),
                    "ci_low": float(values.mean() - mean_t_half),
                    "ci_high": float(values.mean() + mean_t_half),
                    "tau2": float(values.var(ddof=1)),
                    "q": dl["q"],
                    "i_squared": dl["i_squared"],
                    "max_weight": np.nan,
                    "max_weight_target": "",
                },
                {
                    "method": "DerSimonian-Laird-without-sensitive-dataset",
                    "excluded_target": sensitive,
                    **{key: excluded_dl[key] for key in ("estimate", "ci_low", "ci_high", "tau2", "q", "i_squared", "max_weight", "max_weight_target")},
                },
            ]
            for item in methods:
                rows.append({
                    "endpoint": endpoint,
                    "hypothesis": hypothesis,
                    "n_datasets": len(values) - (1 if item["excluded_target"] else 0),
                    **item,
                })
    result = pd.DataFrame(rows)

    checks = {}
    for hypothesis in ("H5", "H6"):
        expected = reference.loc[
            (reference.endpoint == "Ki") & (reference.hypothesis == hypothesis)
        ].iloc[0]
        actual = result.loc[
            (result.endpoint == "Ki")
            & (result.hypothesis == hypothesis)
            & (result.method == "DerSimonian-Laird")
        ].iloc[0]
        checks[hypothesis] = {
            "estimate": float(actual.estimate),
            "ci_low": float(actual.ci_low),
            "ci_high": float(actual.ci_high),
            "maximum_absolute_error": float(
                max(
                    abs(actual.estimate - expected.random_estimate),
                    abs(actual.ci_low - expected.ci_low),
                    abs(actual.ci_high - expected.ci_high),
                )
            ),
        }
        if checks[hypothesis]["maximum_absolute_error"] > 1e-12:
            raise RuntimeError(f"Random-effects regression failed for {hypothesis}")
    return result, {"checks": checks, "passed": True}


def test_summary(
    name: str,
    section: str,
    source_file: str,
    values,
    *,
    family: str,
    units: str = "dataset",
    note: str = "",
    alternative: str = "greater",
) -> dict:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    statistic, one_p, count = one_sided_wilcoxon(values, alternative)
    two_stat, two_p, _ = two_sided_wilcoxon(values)
    return {
        "test_name": name,
        "manuscript_section": section,
        "source_file": source_file,
        "unit": units,
        "n_units": count,
        "median_difference": float(np.median(values)) if count else np.nan,
        "positive_units": int((values > 0).sum()),
        "zero_units": int((values == 0).sum()),
        "negative_units": int((values < 0).sum()),
        "one_sided_statistic": statistic,
        "one_sided_p": one_p,
        "one_sided_alternative": alternative,
        "two_sided_statistic": two_stat,
        "raw_p_two_sided": two_p,
        "correction_family": family,
        "note": note,
    }


def component_values(members: pd.DataFrame, values: pd.Series) -> pd.DataFrame:
    frame = members.copy()
    frame = frame.merge(
        values.rename("dataset_value").reset_index(),
        on=["target", "endpoint"],
        how="inner",
        validate="many_to_one",
    )
    return (
        frame.groupby(["endpoint", "component"], sort=True)
        .dataset_value.mean()
        .reset_index()
    )


def build_posthoc_tests(records: dict[tuple[str, str], pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    cluster_members = pd.read_csv(E7 / "shared_document_cluster_members.csv")
    cluster_reference = pd.read_csv(E7 / "shared_document_cluster_tests.csv")
    comparators = pd.read_csv(E2 / "comparators.csv")
    comparators = comparators.loc[comparators.population == "all"]
    wide = comparators.pivot(
        index=["target", "endpoint"], columns="model", values="auc"
    )
    wide["h3_difference"] = wide.G - 0.5
    wide["h8_difference"] = wide.G - wide["B-range"]
    cluster_values = {}
    for measure in ("h3_difference", "h8_difference"):
        cluster_values[measure] = component_values(cluster_members, wide[measure])

    for hypothesis, measure in (("H3", "h3_difference"), ("H8", "h8_difference")):
        for endpoint in ENDPOINTS:
            values = cluster_values[measure].loc[
                cluster_values[measure].endpoint == endpoint, "dataset_value"
            ]
            row = test_summary(
                f"{hypothesis} shared-document components",
                "Results: dependence checks",
                "S/output/E7/shared_document_cluster_members.csv; S/output/E2/comparators.csv",
                values,
                family="dependence",
                units="shared-document component",
            )
            expected = cluster_reference.loc[
                (cluster_reference.hypothesis == hypothesis)
                & (cluster_reference.endpoint == endpoint)
            ].iloc[0]
            if abs(row["one_sided_p"] - expected.raw_p) > 1e-12:
                raise RuntimeError(f"Shared-cluster regression failed: {hypothesis} {endpoint}")
            rows.append(row)

    time_by_dataset = pd.read_csv(E7 / "time_split_by_dataset.csv")
    time_2015 = time_by_dataset.loc[time_by_dataset.cutoff_year == 2015].copy()
    time_2015["difference"] = (
        time_2015.supported_agreement_rate - time_2015.uncertain_agreement_rate
    )
    eligible = time_2015.loc[
        (time_2015.supported_n_predictions >= 5)
        & (time_2015.uncertain_n_predictions >= 5)
        & time_2015.difference.notna()
    ].copy()
    for endpoint in ENDPOINTS:
        values = eligible.loc[eligible.endpoint == endpoint, "difference"]
        rows.append(test_summary(
            f"2015 supported minus uncertain ({endpoint}, development targets included)",
            "Results: agreement with later documents",
            "S/output/E7/time_split_by_dataset.csv",
            values,
            family="time split",
        ))
        without = eligible.loc[
            (eligible.endpoint == endpoint)
            & ~eligible.target.isin(DEVELOPMENT_TARGETS),
            "difference",
        ]
        rows.append(test_summary(
            f"2015 supported minus uncertain ({endpoint}, development targets excluded)",
            "Results: agreement with later documents",
            "S/output/E7/time_split_by_dataset.csv",
            without,
            family="time split",
        ))

    ki_members = cluster_members.loc[cluster_members.endpoint == "Ki"]
    ki_time = eligible.loc[eligible.endpoint == "Ki", ["target", "endpoint", "difference"]]
    ki_frame = ki_members.merge(ki_time, on=["target", "endpoint"], validate="one_to_one")
    ki_component = ki_frame.groupby(["endpoint", "component"], sort=True).difference.mean()
    rows.append(test_summary(
        "2015 supported minus uncertain (Ki, shared-document components)",
        "Results: dependence checks",
        "S/output/E7/time_split_by_dataset.csv; S/output/E7/shared_document_cluster_members.csv",
        ki_component,
        family="time split",
        units="shared-document component",
    ))

    latest = build_temporal_frame(records, pd.read_csv(E2 / "target_inclusion.csv"), "latest")
    near = latest.loc[latest.minimum_past_gap > NEAR_DUPLICATE_GAP].copy()
    write_csv(
        latest[["target", "endpoint", "cutoff_year", "entry", "future_unit", "future_value", "minimum_past_gap"]],
        "time_split_near_duplicate_gaps.csv",
        qa=True,
    )
    near_2015 = near.loc[near.cutoff_year == 2015].copy()
    near_2015["difference"] = (
        near_2015.display_status.eq("supported").groupby(
            [near_2015.target, near_2015.endpoint]
        ).transform("mean")
    )
    near_dataset_rows = []
    for (target, endpoint), group in near_2015.groupby(["target", "endpoint"], sort=True):
        supported = group.loc[group.display_status == "supported"]
        uncertain = group.loc[group.display_status == "uncertain"]
        if len(supported) < 5 or len(uncertain) < 5:
            continue
        near_dataset_rows.append({
            "target": target,
            "endpoint": endpoint,
            "n_predictions": len(group),
            "supported_n": len(supported),
            "uncertain_n": len(uncertain),
            "difference": supported.agree.mean() - uncertain.agree.mean(),
        })
    near_dataset = pd.DataFrame(near_dataset_rows)
    write_csv(near_dataset, "time_split_near_duplicate_by_dataset.csv", qa=True)
    for endpoint in ENDPOINTS:
        values = near_dataset.loc[near_dataset.endpoint == endpoint, "difference"]
        rows.append(test_summary(
            f"2015 supported minus uncertain after removing near-duplicate future units ({endpoint})",
            "Results: agreement with later documents",
            "qa/E12b/time_split_near_duplicate_by_dataset.csv; derived with e4b_analysis.build_temporal_frame",
            values,
            family="time split",
        ))

    fixed = pd.read_csv(E6 / "recall_fixedk_77.csv")
    fixed = fixed.loc[fixed.layout.isin(["neighbor_refined", "minimum_displacement"])].copy()
    l1 = fixed.loc[
        fixed.layout == "neighbor_refined", ["target", "endpoint", "rho_fixed"]
    ].rename(columns={"rho_fixed": "difference"})
    pivot = fixed.pivot(
        index=["target", "endpoint"], columns="layout", values="fixed_e_grid"
    ).reset_index()
    pivot["difference"] = pivot.neighbor_refined - pivot.minimum_displacement
    l2 = pivot[["target", "endpoint", "difference"]]
    reference_fixed = pd.read_csv(E6 / "l1_l2_fixedk.csv")
    for test_id, data, label in (
        ("L1", l1, "fixed-k recall correlation"),
        ("L2", l2, "fixed-k activity-neighbor error reduction"),
    ):
        for endpoint in ENDPOINTS:
            values = data.loc[data.endpoint == endpoint, "difference"]
            row = test_summary(
                f"{label} ({endpoint})",
                "Results: neighborhoods and layout",
                "S/output/E6/recall_fixedk_77.csv",
                values,
                family="layout",
                alternative="less",
            )
            expected = reference_fixed.loc[
                (reference_fixed.test_id == test_id)
                & (reference_fixed.endpoint == endpoint)
                & (reference_fixed.population == "all77")
                & (reference_fixed.status == "complete_all77")
            ].iloc[0]
            if abs(row["one_sided_p"] - expected.raw_p) > 1e-12 or row["n_units"] != expected.n_targets:
                raise RuntimeError(f"Fixed-k regression failed: {test_id} {endpoint}")
            rows.append(row)

    result = pd.DataFrame(rows)
    for family in ("dependence", "time split"):
        locations = result.index[result.correction_family == family]
        result.loc[locations, "holm_p_two_sided"] = holm_adjust(
            result.loc[locations, "raw_p_two_sided"]
        )
    for endpoint in ENDPOINTS:
        locations = result.index[
            (result.correction_family == "layout") & (result.test_name.str.endswith(f"({endpoint})"))
        ]
        result.loc[locations, "holm_p_two_sided"] = holm_adjust(
            result.loc[locations, "raw_p_two_sided"]
        )
    return result, near_dataset


def top_label_ece_from_arrays(confidence, correct, group_codes=None) -> tuple[float, float | None]:
    bins = np.minimum((confidence * 10).astype(int), 9)
    pooled = sum(
        abs(correct[bins == index].sum() - confidence[bins == index].sum())
        for index in range(10)
    ) / len(confidence)
    if group_codes is None:
        return float(pooled), None
    count = int(group_codes.max()) + 1
    encoded = group_codes * 10 + bins
    correct_sum = np.bincount(encoded, weights=correct, minlength=count * 10)
    confidence_sum = np.bincount(encoded, weights=confidence, minlength=count * 10)
    group_n = np.bincount(group_codes, minlength=count)
    ece = (
        np.abs(correct_sum - confidence_sum).reshape(count, 10).sum(axis=1)
        / group_n
    )
    return float(pooled), float(np.median(ece))


def simulation_median_ece(confidence, group_codes, draws: np.ndarray) -> np.ndarray:
    count = int(group_codes.max()) + 1
    bins = np.minimum((confidence * 10).astype(int), 9)
    encoded = group_codes * 10 + bins
    one_hot = np.zeros((len(confidence), count * 10), dtype=float)
    one_hot[np.arange(len(confidence)), encoded] = 1.0
    correct_sum = draws @ one_hot
    confidence_sum = one_hot.T @ confidence
    group_n = np.bincount(group_codes, minlength=count)
    ece = (
        np.abs(correct_sum - confidence_sum).reshape(len(draws), count, 10).sum(axis=2)
        / group_n
    )
    return np.median(ece, axis=1)


def build_ece_reference(predictions: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    selected = predictions.loc[
        (predictions.method == "G") & (predictions.population == "all")
    ].copy()
    rng = np.random.default_rng(SIMULATION_SEED)
    rows = []
    regression = {}
    for endpoint in ENDPOINTS:
        frame = selected.loc[selected.endpoint == endpoint].reset_index(drop=True)
        probabilities = frame[["p0", "p1", "p2"]].to_numpy(float)
        observed = frame.observed.to_numpy(int)
        confidence = probabilities.max(axis=1)
        correct = (probabilities.argmax(axis=1) == observed).astype(float)
        group_codes, groups = pd.factorize(frame.target, sort=True)
        pooled, median = top_label_ece_from_arrays(confidence, correct, group_codes)
        draws = rng.binomial(
            1,
            np.broadcast_to(confidence, (1000, len(confidence))),
            size=(1000, len(confidence)),
        ).astype(float)
        simulated_medians = simulation_median_ece(confidence, group_codes, draws)
        simulated_pooled = np.array([
            sum(
                abs(draw[bins].sum() - confidence[bins].sum())
                for bins in [np.minimum((confidence * 10).astype(int), 9) == index for index in range(10)]
            ) / len(confidence)
            for draw in draws
        ])
        rows.append({
            "endpoint": endpoint,
            "n_predictions": len(frame),
            "n_datasets": len(groups),
            "simulation_repetitions": 1000,
            "simulation_seed": SIMULATION_SEED,
            "outcome_probability": "max(p0,p1,p2), matching top-label ECE",
            "observed_pooled_ece": pooled,
            "observed_dataset_median_ece": median,
            "simulated_pooled_ece_mean": float(simulated_pooled.mean()),
            "simulated_dataset_median_ece_mean": float(simulated_medians.mean()),
            "simulated_dataset_median_ece_q025": float(np.percentile(simulated_medians, 2.5)),
            "simulated_dataset_median_ece_q975": float(np.percentile(simulated_medians, 97.5)),
        })
        regression[endpoint] = {
            "pooled": pooled,
            "dataset_median": median,
        }
    result = pd.DataFrame(rows)
    if abs(regression["IC50"]["pooled"] - 0.0435160441444185) > 1e-12 or abs(regression["Ki"]["pooled"] - 0.0250259849085929) > 1e-12:
        raise RuntimeError("Top-label ECE regression failed")
    return result, {"observed": regression, "passed": True}


def bootstrap_time_population(frame: pd.DataFrame) -> dict:
    point = point_metric_values(frame)
    weights = two_stage_weights(frame, BOOT_SEED, BOOT_N)
    draws = metric_draws(frame, weights)
    result = {}
    for name in ("supported_agreement_rate", "uncertain_agreement_rate"):
        low, high, valid = percentile_ci(draws[name])
        result[f"{name}_ci_low"] = low
        result[f"{name}_ci_high"] = high
        result[f"{name}_valid"] = valid
    return {**point, **result}


def build_time_populations() -> tuple[pd.DataFrame, dict]:
    population_frames = {
        "planned_including_development_IC50": pd.read_csv(
            E2_QA / "time_split_predictions.csv",
            float_precision="round_trip",
        ),
        "52_evaluation_datasets": pd.read_csv(
            E2_QA / "time_split_predictions_R0_without_original_ic50.csv",
            float_precision="round_trip",
        ),
    }
    rows = []
    checks = {}
    for population, data in population_frames.items():
        for cutoff in (2015, 2010):
            for endpoint in ENDPOINTS:
                frame = data.loc[
                    (data.cutoff_year == cutoff) & (data.endpoint == endpoint)
                ].copy()
                boot = bootstrap_time_population(frame)
                single = frame.loc[
                    (frame.display_status == "supported") & (frame.m_past == 1)
                ]
                rows.append({
                    "cutoff_year": cutoff,
                    "endpoint": endpoint,
                    "population": population,
                    "n_predictions": len(frame),
                    "n_datasets": int(frame.target.nunique()),
                    "supported_agreement_rate": boot["supported_agreement_rate"],
                    "supported_ci_low": boot["supported_agreement_rate_ci_low"],
                    "supported_ci_high": boot["supported_agreement_rate_ci_high"],
                    "uncertain_agreement_rate": boot["uncertain_agreement_rate"],
                    "uncertain_ci_low": boot["uncertain_agreement_rate_ci_low"],
                    "uncertain_ci_high": boot["uncertain_agreement_rate_ci_high"],
                    "g_auc": boot["g_auc"],
                    "single_previous_document_supported_n": len(single),
                    "single_previous_document_agreement_rate": (
                        float(single.agree.mean()) if len(single) else np.nan
                    ),
                    "single_previous_document_mean_probability": (
                        float(single.p_displayed.mean()) if len(single) else np.nan
                    ),
                    "bootstrap_method": "two_stage: dataset then entry",
                    "bootstrap_repetitions": BOOT_N,
                    "bootstrap_seed": BOOT_SEED,
                })
    result = pd.DataFrame(rows)
    planned_2015_ic50 = result.loc[
        (result.population == "planned_including_development_IC50")
        & (result.cutoff_year == 2015)
        & (result.endpoint == "IC50")
    ].iloc[0]
    e7_boot = pd.read_csv(E7 / "two_stage_bootstrap.csv")
    expected = e7_boot.loc[
        (e7_boot.analysis_section == "time_split_2015")
        & (e7_boot.endpoint == "IC50")
        & (e7_boot.population == "all")
        & (e7_boot.bootstrap_method == "two_stage")
    ].iloc[0]
    checks["planned_2015_IC50"] = {
        "n_predictions": int(planned_2015_ic50.n_predictions),
        "supported": float(planned_2015_ic50.supported_agreement_rate),
        "ci_low": float(planned_2015_ic50.supported_ci_low),
        "ci_high": float(planned_2015_ic50.supported_ci_high),
        "e7_maximum_absolute_error": float(
            max(
                abs(planned_2015_ic50.supported_agreement_rate - expected.supported_agreement_rate_estimate),
                abs(planned_2015_ic50.supported_ci_low - expected.supported_agreement_rate_ci_low),
                abs(planned_2015_ic50.supported_ci_high - expected.supported_agreement_rate_ci_high),
            )
        ),
    }
    if checks["planned_2015_IC50"]["e7_maximum_absolute_error"] > 1e-12:
        raise RuntimeError("Planned time-split bootstrap does not reproduce E7")
    expected_52 = result.loc[
        (result.population == "52_evaluation_datasets")
        & (result.cutoff_year == 2015)
        & (result.endpoint == "IC50")
    ].iloc[0]
    checks["single_previous_document_n"] = {
        "planned_2015_IC50": int(planned_2015_ic50.single_previous_document_supported_n),
        "52_evaluation_2015_IC50": int(expected_52.single_previous_document_supported_n),
    }
    if checks["single_previous_document_n"] != {
        "planned_2015_IC50": 225,
        "52_evaluation_2015_IC50": 206,
    }:
        raise RuntimeError("Unexpected single-previous-document counts")
    return result, {"checks": checks, "passed": True}


def coverage_bucket(value: float) -> str:
    if not np.isfinite(value):
        return "no_core"
    if value >= 0.8 - 1e-12:
        return "coverage_ge_0_8"
    if value > 0.5 + 1e-12:
        return "coverage_gt_0_5_lt_0_8"
    return "coverage_le_0_5"


def core_pattern(core: str, rule: str):
    if not core or pd.isna(core):
        return None
    pattern = Chem.MolFromSmarts(core) if rule == "fmcs" else Chem.MolFromSmiles(core)
    if pattern is None:
        raise ValueError(f"Could not parse {rule} core: {core}")
    return pattern


def build_cores() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    cache = pd.read_csv(
        ROOT / "paper/output/reliability_v1/cores/group_core_comparison.csv",
        low_memory=False,
    )
    members = pd.read_csv(E1B / "core_group_members.csv")
    members["size_stratum"] = pd.Categorical(
        members.size_stratum, categories=SIZE_ORDER, ordered=True
    )
    molecule_tables = {}
    for target in ("chembl205", "chembl204", "chembl240"):
        data = pd.read_csv(
            ROOT / f"paper/output/revision_v6/{target}/map_molecules.csv",
            low_memory=False,
        )
        molecule_tables[target] = dict(
            zip(
                data.molecule_identity_key.astype(str),
                data.canonical_smiles.astype(str),
            )
        )
    detail_rows = []
    keyed_members = {
        (row.target, row.group): str(row.member_ids).split("|")
        for row in members.itertuples(index=False)
    }
    for row in cache.itertuples(index=False):
        ids = keyed_members[(row.target, row.group)]
        mols = [Chem.MolFromSmiles(molecule_tables[row.target][key]) for key in ids]
        if any(mol is None for mol in mols):
            raise ValueError(f"Invalid detail-target SMILES in {row.target}/{row.group}")
        core_string = row.core_smarts if row.rule == "fmcs" else row.core_smiles
        pattern = core_pattern(core_string, row.rule)
        for setting, use_chirality in (
            ("original", row.rule != "fmcs"),
            ("unified_use_chirality", True),
        ):
            coverage_count = (
                sum(bool(mol.HasSubstructMatch(pattern, useChirality=use_chirality)) for mol in mols)
                if pattern is not None
                else 0
            )
            detail_rows.append({
                "setting": setting,
                "target": row.target,
                "group": row.group,
                "size": int(row.size),
                "size_stratum": members.loc[
                    (members.target == row.target) & (members.group == row.group),
                    "size_stratum",
                ].iloc[0],
                "rule": row.rule,
                "core_string": core_string,
                "use_chirality": bool(use_chirality),
                "has_core": bool(row.has_core),
                "coverage_count": int(coverage_count),
                "coverage_fraction": coverage_count / len(mols) if pattern is not None else np.nan,
                "core_heavy_atoms": (
                    int(row.core_heavy_atoms) if bool(row.has_core) else np.nan
                ),
                "is_benzene": bool(row.is_benzene),
            })
    detail = pd.DataFrame(detail_rows)

    grouped_rows = []
    for (setting, rule, stratum), frame in detail.groupby(
        ["setting", "rule", "size_stratum"], sort=True, observed=True
    ):
        buckets = frame.coverage_fraction.map(coverage_bucket)
        grouped_rows.append({
            "setting": setting,
            "rule": rule,
            "size_stratum": stratum,
            "coverage_ge_0_8": int((buckets == "coverage_ge_0_8").sum()),
            "coverage_gt_0_5_lt_0_8": int((buckets == "coverage_gt_0_5_lt_0_8").sum()),
            "coverage_le_0_5": int((buckets == "coverage_le_0_5").sum()),
            "no_ring_scaffold_groups": int((buckets == "no_core").sum()),
            "groups": len(frame),
            "mean_core_heavy_atoms": float(
                frame.loc[frame.has_core, "core_heavy_atoms"].mean()
            ),
            "benzene_groups": int(frame.is_benzene.sum()),
        })
    chirality = pd.DataFrame(grouped_rows).sort_values(
        ["setting", "rule", "size_stratum"]
    )

    figure = pd.read_csv(FIGURE_DATA / "5e_core_coverage_by_group_size.csv")
    original = chirality.loc[chirality.setting == "original"].reset_index(drop=True)
    maximum_error = 0.0
    for row in figure.itertuples(index=False):
        actual = original.loc[
            (original.rule == row.rule) & (original.size_stratum == row.size_stratum)
        ].iloc[0]
        for expected, observed in (
            (row.coverage_ge_0_8, actual.coverage_ge_0_8),
            (row.coverage_gt_0_5_lt_0_8, actual.coverage_gt_0_5_lt_0_8),
            (row.coverage_le_0_5, actual.coverage_le_0_5),
            (row.no_core, actual.no_ring_scaffold_groups),
            (row.groups, actual.groups),
        ):
            maximum_error = max(maximum_error, abs(float(expected) - float(observed)))
        maximum_error = max(
            maximum_error,
            abs(float(row.mean_core_heavy_atoms) - float(actual.mean_core_heavy_atoms)),
        )
    if maximum_error > 1e-12:
        raise RuntimeError(f"Core Figure 5e regression failed: {maximum_error}")

    unified = detail.loc[detail.setting == "unified_use_chirality"]
    cost_rows = []
    for rule, frame in unified.groupby("rule", sort=True):
        for stratum in SIZE_ORDER:
            selected = frame.loc[frame.size_stratum == stratum]
            cost_rows.append({
                "setting": "unified_use_chirality",
                "rule": rule,
                "size_stratum": stratum,
                "groups": len(selected),
                "mean_core_heavy_atoms": float(
                    selected.loc[selected.has_core, "core_heavy_atoms"].mean()
                ),
                "benzene_groups": int(selected.is_benzene.sum()),
                "no_ring_scaffold_groups": int((~selected.has_core).sum()),
            })
    costs = pd.DataFrame(cost_rows)
    checks = {
        "figure_5e_maximum_absolute_error": maximum_error,
        "original_benzene_counts": original.groupby("rule").benzene_groups.sum().to_dict(),
        "unified_changed_group_core_matches": int(
            (
                detail.loc[detail.setting == "unified_use_chirality"]
                .set_index(["target", "group", "rule"])
                .coverage_count
                != detail.loc[detail.setting == "original"]
                .set_index(["target", "group", "rule"])
                .coverage_count
            ).sum()
        ),
        "passed": True,
    }
    return chirality, costs, checks


def relation_series(raw: pd.DataFrame) -> pd.Series:
    return (
        raw.standard_relation.fillna("")
        .astype(str)
        .str.strip()
        .str.replace(r"^'([=<>~]+)'$", r"\1", regex=True)
    )


def build_censored_records(records: dict[tuple[str, str], pd.DataFrame]) -> tuple[pd.DataFrame, dict]:
    inclusion = pd.read_csv(E2 / "target_inclusion.csv")
    included = inclusion.loc[inclusion.included].copy()
    rows = []
    for item in included.itertuples(index=False):
        input_path = API / "activities" / Path(str(item.input_file)).name
        raw = pd.read_csv(input_path, low_memory=False)
        relation = relation_series(raw)
        non_exact = int((relation != "=").sum())
        retained = len(records[(item.target, item.endpoint)])
        other = len(raw) - non_exact - retained
        if other < 0:
            raise ValueError(f"Negative other-exclusion count: {item.target} {item.endpoint}")
        rows.append({
            "population": "52_evaluation_datasets",
            "target": item.target,
            "endpoint": item.endpoint,
            "input_file": str(input_path.relative_to(ROOT)),
            "records_total": len(raw),
            "records_retained": retained,
            "excluded_non_exact_relation": non_exact,
            "excluded_other_reasons": other,
        })

    original_manifest = json.loads(
        (API / "download_original_ic50_manifest.json").read_text()
    )
    manifest_by_name = {item["file"]: item for item in original_manifest["files"]}
    for target in ("CHEMBL205", "CHEMBL204", "CHEMBL240"):
        report_path = (
            ROOT
            / f"paper/output/revision_v6/{target.lower()}/{target.lower()}_curation_report.json"
        )
        report = json.loads(report_path.read_text())
        file_name = f"{target}_IC50.csv"
        total = int(report["input_rows"])
        retained = int(report["retained_activity_records"])
        non_exact = int(report["exclusion_counts"].get("non_exact_standard_relation", 0))
        rows.append({
            "population": "detailed_IC50",
            "target": target,
            "endpoint": "IC50",
            "input_file": f"paper/output/revision_v6/{target.lower()}/{target.lower()}_curation_report.json",
            "records_total": total,
            "records_retained": retained,
            "excluded_non_exact_relation": non_exact,
            "excluded_other_reasons": total - retained - non_exact,
            "e2_api_file": str((API / "activities" / file_name).relative_to(ROOT)),
            "e2_api_records": int(manifest_by_name[file_name]["records"]),
        })
    result = pd.DataFrame(rows)
    ca2 = result.loc[
        (result.population == "detailed_IC50") & (result.target == "CHEMBL205")
    ].iloc[0]
    checks = {
        "ca2_records_total": int(ca2.records_total),
        "ca2_records_retained": int(ca2.records_retained),
        "ca2_excluded_non_exact_relation": int(ca2.excluded_non_exact_relation),
        "rows": len(result),
        "passed": bool(
            int(ca2.records_total) == 2281
            and int(ca2.records_retained) == 1100
            and int(ca2.excluded_non_exact_relation) == 1056
        ),
    }
    if not checks["passed"]:
        raise RuntimeError("CA II censoring regression failed")
    return result, checks


def build_selection_summary() -> tuple[pd.DataFrame, dict]:
    inclusion = pd.read_csv(E2 / "target_inclusion.csv")
    sampled = inclusion.loc[inclusion.role == "sampled"].copy()
    rows = []
    for endpoint in ENDPOINTS:
        frame = sampled.loc[sampled.endpoint == endpoint]
        insufficient = frame.reason.eq(
            "fewer than 30 entries with two or more evidence units"
        )
        identity_conflict = frame.reason.str.startswith("curation failure:")
        other = ~(frame.included | insufficient | identity_conflict)
        rows.append({
            "endpoint": endpoint,
            "sampled_target_endpoints": len(frame),
            "included": int(frame.included.sum()),
            "excluded_fewer_than_30_multi_unit_entries": int(insufficient.sum()),
            "excluded_parent_identity_conflict": int(identity_conflict.sum()),
            "excluded_other": int(other.sum()),
        })
    totals = {
        "endpoint": "Total",
        "sampled_target_endpoints": int(sum(row["sampled_target_endpoints"] for row in rows)),
        "included": int(sum(row["included"] for row in rows)),
        "excluded_fewer_than_30_multi_unit_entries": int(sum(row["excluded_fewer_than_30_multi_unit_entries"] for row in rows)),
        "excluded_parent_identity_conflict": int(sum(row["excluded_parent_identity_conflict"] for row in rows)),
        "excluded_other": int(sum(row["excluded_other"] for row in rows)),
    }
    rows.append(totals)
    result = pd.DataFrame(rows)
    checks = {
        "sampled": totals["sampled_target_endpoints"],
        "included": totals["included"],
        "fewer_than_30": totals["excluded_fewer_than_30_multi_unit_entries"],
        "parent_identity_conflict": totals["excluded_parent_identity_conflict"],
        "passed": totals["sampled_target_endpoints"] == 80 and totals["included"] == 50,
    }
    if not checks["passed"]:
        raise RuntimeError("Selection regression failed")
    return result, checks


def markdown_table(frame: pd.DataFrame) -> str:
    def cell(value):
        if pd.isna(value):
            return "NA"
        if isinstance(value, (float, np.floating)):
            return f"{value:.6g}"
        return str(value).replace("|", "\\|")

    lines = [
        "| " + " | ".join(frame.columns) + " |",
        "| " + " | ".join("---" for _ in frame.columns) + " |",
    ]
    lines.extend(
        "| " + " | ".join(cell(value) for value in row) + " |"
        for row in frame.itertuples(index=False, name=None)
    )
    return "\n".join(lines)


def write_outputs() -> dict:
    OUT.mkdir(parents=True, exist_ok=False)
    QA.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    started_at = datetime.now(timezone.utc).isoformat()
    predictions = read_predictions()
    records = load_records()

    scores, scores_regression = build_scores_h4_h7(predictions)
    random_effects, random_effects_regression = build_random_effect_alternatives()
    posthoc, near_dataset = build_posthoc_tests(records)
    ece, ece_regression = build_ece_reference(predictions)
    time_population, time_population_regression = build_time_populations()
    chirality, costs, core_regression = build_cores()
    censored, censoring_regression = build_censored_records(records)
    selection, selection_regression = build_selection_summary()

    outputs = {
        "scores_h4_h7.csv": write_csv(scores, "scores_h4_h7.csv"),
        "re_h5_h6_alternatives.csv": write_csv(random_effects, "re_h5_h6_alternatives.csv"),
        "posthoc_tests_two_sided.csv": write_csv(posthoc, "posthoc_tests_two_sided.csv"),
        "ece_reference.csv": write_csv(ece, "ece_reference.csv"),
        "time_split_populations.csv": write_csv(time_population, "time_split_populations.csv"),
        "core_chirality.csv": write_csv(chirality, "core_chirality.csv"),
        "core_costs.csv": write_csv(costs, "core_costs.csv"),
        "censored_records.csv": write_csv(censored, "censored_records.csv"),
        "selection_summary.csv": write_csv(selection, "selection_summary.csv"),
    }
    checks = {
        "scores_h4_h7": scores_regression,
        "random_effects": random_effects_regression,
        "ece": ece_regression,
        "time_split_populations": time_population_regression,
        "cores": core_regression,
        "censored_records": censoring_regression,
        "selection": selection_regression,
    }
    elapsed = time.perf_counter() - started
    finished_at = datetime.now(timezone.utc).isoformat()

    manifest = {
        "started_at_utc": started_at,
        "finished_at_utc": finished_at,
        "elapsed_seconds": elapsed,
        "python": sys.version,
        "platform": platform.platform(),
        "chemgridmap_file": str(Path(chemgridmap.__file__).resolve()),
        "chemgridmap_frozen_sha256": sha256(Path(chemgridmap.__file__).resolve()),
        "checks": checks,
        "outputs_sha256": {
            path.name: sha256(path) for path in sorted(OUT.iterdir()) if path.is_file()
        },
        "qa_sha256": {
            path.name: sha256(path) for path in sorted(QA.iterdir()) if path.is_file()
        },
    }
    (OUT / "manifest.json").write_text(
        json.dumps(json_safe(manifest), indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    )
    (QA / "timing.json").write_text(
        json.dumps(
            json_safe({
                "started_at_utc": started_at,
                "finished_at_utc": finished_at,
                "elapsed_seconds": elapsed,
            }),
            indent=2,
            ensure_ascii=False,
            allow_nan=False,
        )
        + "\n"
    )
    return {"outputs": outputs, "checks": checks, "elapsed_seconds": elapsed}


if __name__ == "__main__":
    result = write_outputs()
    print(json.dumps(json_safe({
        "outputs": sorted(result["outputs"]),
        "checks": {name: value.get("passed") for name, value in result["checks"].items()},
        "elapsed_seconds": result["elapsed_seconds"],
    }), indent=2))
