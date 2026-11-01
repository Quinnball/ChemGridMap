"""E22: post-hoc reliability by displayed activity class."""
from __future__ import annotations

from datetime import datetime
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[3]
S = ROOT / "paper/strengthening"
OUT = S / "output/E22"
DATA = S / "data/chembl37_api_validation"
ACTIVITIES = DATA / "activities"
E2 = S / "output/E2"
E5C = S / "output/E5c"
E9 = S / "output/E9"
E11 = S / "output/E11"
E15 = S / "output/E15"

sys.path[:0] = [
    str(S / "source"),
    str(ROOT / "paper"),
    str(ROOT / "src"),
]

import numpy as np
import pandas as pd
from scipy import stats

from chemgridmap.evidence import CLASSES
from e15_leakfree import load_e2_records, wilcoxon_details
from e2_models import curate_activity_table
from e7_analysis import holm_adjust, two_stage_weights
from e9_analysis import random_effect_details
from manifest_utils import write_manifest

ANALYSIS_ROLE = "post hoc"
ENDPOINTS = ("IC50", "Ki")
DESIGNS = ("entry", "document", "near_duplicate", "time_2015", "E11")
BOOT_SEED = 20261105
BOOT_REPETITIONS = 2000
SUPPORT_THRESHOLD = 0.8
MIN_STATE_PREDICTIONS = 5
MIN_BIAS_PREDICTIONS = 10
TEST_CATEGORIES = ("inactive", "active", "inactive+active")
BIAS_CATEGORIES = tuple(CLASSES)
BIAS_OUTPUT_CATEGORIES = (*BIAS_CATEGORIES, "all")
BIAS_HOLM_FAMILY = "inactive, medium, active x IC50, Ki"


def _design_mask(frame: pd.DataFrame, design: str) -> np.ndarray:
    if design == "entry":
        return (
            frame.split_design.eq("entry")
            & frame.merge_basis.eq("same_value")
        ).to_numpy(bool)
    if design == "document":
        return frame.split_design.eq("document").to_numpy(bool)
    if design == "near_duplicate":
        return (
            frame.split_design.eq("entry")
            & frame.merge_basis.eq("near_duplicate_0.1")
        ).to_numpy(bool)
    if design == "time_2015":
        return (
            frame.split_design.eq("time")
            & frame.cutoff_year.eq(2015)
        ).to_numpy(bool)
    raise ValueError(f"Unknown E5c design: {design}")


def add_support_status(frame: pd.DataFrame, source: str) -> pd.DataFrame:
    """Use saved support for E5c and saved G status for E11."""
    result = frame.copy()
    if source == "E5c":
        supported = result.support.to_numpy(float) >= SUPPORT_THRESHOLD
    elif source == "E11":
        supported = result.g_status.eq("supported").to_numpy(bool)
    else:
        raise ValueError(f"Unknown support source: {source}")
    result["e22_status"] = np.where(supported, "supported", "uncertain")
    return result


def load_designs() -> dict[str, pd.DataFrame]:
    columns = [
        "split_design",
        "merge_basis",
        "target",
        "endpoint",
        "entry",
        "displayed",
        "p_displayed",
        "support",
        "agree",
        "method",
        "cutoff_year",
    ]
    predictions = pd.read_csv(
        E5C / "predictions_e5c.csv",
        usecols=columns,
        low_memory=False,
        float_precision="round_trip",
    )
    predictions = predictions.loc[predictions.method == "G"].copy()
    designs = {
        name: add_support_status(
            predictions.loc[_design_mask(predictions, name)].reset_index(drop=True),
            "E5c",
        )
        for name in DESIGNS[:-1]
    }
    e11 = pd.read_csv(
        E11 / "predictions_e11.csv",
        low_memory=False,
        float_precision="round_trip",
    )
    designs["E11"] = add_support_status(
        e11.loc[e11.method == "G"].reset_index(drop=True),
        "E11",
    )
    return designs


def class_mask(frame: pd.DataFrame, category: str) -> np.ndarray:
    displayed = frame.displayed.to_numpy(int)
    if category == "all":
        return np.ones(len(frame), dtype=bool)
    if category == "inactive+active":
        return displayed != CLASSES.index("medium")
    return displayed == CLASSES.index(category)


def bin_composition(predictions: pd.DataFrame) -> pd.DataFrame:
    """Describe the datasets contributing to each class and probability bin."""
    rows = []
    for endpoint in ENDPOINTS:
        endpoint_frame = predictions.loc[predictions.endpoint == endpoint].copy()
        endpoint_frame["probability_bin"] = np.clip(
            np.floor(endpoint_frame.p_displayed * 5 + 1e-12).astype(int),
            0,
            4,
        )
        for display_index, display_class in enumerate(CLASSES):
            selected = endpoint_frame.loc[endpoint_frame.displayed == display_index]
            for probability_bin in range(5):
                group = selected.loc[
                    selected.probability_bin == probability_bin
                ]
                n_predictions = len(group)
                dataset_counts = group.groupby("target", sort=True).size()
                if n_predictions:
                    largest_dataset, largest_count = min(
                        dataset_counts.items(),
                        key=lambda item: (-int(item[1]), str(item[0])),
                    )
                    largest_fraction = int(largest_count) / n_predictions
                else:
                    largest_dataset = ""
                    largest_fraction = np.nan
                rows.append(
                    {
                        "analysis_role": ANALYSIS_ROLE,
                        "design": "entry",
                        "endpoint": endpoint,
                        "display_class": display_class,
                        "probability_bin": probability_bin,
                        "n_predictions": n_predictions,
                        "n_datasets": int(len(dataset_counts)),
                        "largest_dataset": str(largest_dataset),
                        "largest_dataset_fraction": largest_fraction,
                        "agreement_rate": (
                            float(group.agree.mean()) if n_predictions else np.nan
                        ),
                    }
                )
    return pd.DataFrame(rows)


def validate_bin_composition(composition: pd.DataFrame) -> None:
    """Check the archival bin counts and the reviewer-relevant Ki medium row."""
    if len(composition) != 30:
        raise RuntimeError(
            f"Expected 30 figure 4a bin composition rows, found {len(composition)}"
        )
    for endpoint in ENDPOINTS:
        for display_class in CLASSES:
            selected = composition.loc[
                composition.endpoint.eq(endpoint)
                & composition.display_class.eq(display_class)
            ]
            if selected.probability_bin.tolist() != [0, 1, 2, 3, 4]:
                raise RuntimeError(
                    f"Incomplete probability bins for {endpoint} {display_class}"
                )
    ki_medium_top = row_for(
        composition,
        endpoint="Ki",
        display_class="medium",
        probability_bin=4,
    )
    if (
        int(ki_medium_top.n_predictions) != 960
        or int(ki_medium_top.n_datasets) != 1
        or ki_medium_top.largest_dataset != "CHEMBL5024"
        or not np.isclose(
            float(ki_medium_top.largest_dataset_fraction),
            1.0,
            atol=1e-12,
            rtol=0,
        )
        or not np.isclose(
            float(ki_medium_top.agreement_rate),
            0.9958333333333333,
            atol=1e-12,
            rtol=0,
        )
    ):
        raise RuntimeError("Ki medium highest bin composition mismatch")


def cell_metrics(
    frame: pd.DataFrame,
    weights: np.ndarray,
    mask: np.ndarray,
) -> dict[str, float | int]:
    """Match E9's two-stage cell bootstrap: resample the full endpoint frame."""
    mask = np.asarray(mask, dtype=bool)
    selected = frame.loc[mask]
    n = len(selected)
    if not n:
        return {
            "n_predictions": 0,
            "n_datasets": 0,
            "mean_probability": np.nan,
            "agreement_rate": np.nan,
            "agreement_ci_low": np.nan,
            "agreement_ci_high": np.nan,
            "bias": np.nan,
            "bias_ci_low": np.nan,
            "bias_ci_high": np.nan,
            "valid_repetitions": 0,
        }
    selected_weights = weights[:, mask]
    denominator = selected_weights.sum(axis=1)
    agreement = np.divide(
        (selected_weights * selected.agree.to_numpy(float)).sum(axis=1),
        denominator,
        out=np.full(len(weights), np.nan),
        where=denominator > 0,
    )
    residual = selected.p_displayed.to_numpy(float) - selected.agree.to_numpy(float)
    bias = np.divide(
        (selected_weights * residual).sum(axis=1),
        denominator,
        out=np.full(len(weights), np.nan),
        where=denominator > 0,
    )
    valid = np.isfinite(agreement) & np.isfinite(bias)
    agreement_ci = (
        np.percentile(agreement[valid], [2.5, 97.5])
        if valid.any()
        else (np.nan, np.nan)
    )
    bias_ci = (
        np.percentile(bias[valid], [2.5, 97.5])
        if valid.any()
        else (np.nan, np.nan)
    )
    mean_probability = float(selected.p_displayed.mean())
    agreement_rate = float(selected.agree.mean())
    return {
        "n_predictions": n,
        "n_datasets": int(selected.target.nunique()),
        "mean_probability": mean_probability,
        "agreement_rate": agreement_rate,
        "agreement_ci_low": float(agreement_ci[0]),
        "agreement_ci_high": float(agreement_ci[1]),
        "bias": mean_probability - agreement_rate,
        "bias_ci_low": float(bias_ci[0]),
        "bias_ci_high": float(bias_ci[1]),
        "valid_repetitions": int(valid.sum()),
    }


def reliability_by_design(designs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for design in DESIGNS:
        frame = designs[design]
        for endpoint in ENDPOINTS:
            endpoint_frame = frame.loc[frame.endpoint == endpoint].reset_index(drop=True)
            weights = two_stage_weights(
                endpoint_frame,
                seed=BOOT_SEED,
                repetitions=BOOT_REPETITIONS,
            )
            for display_class in (*CLASSES, "all"):
                display_mask = class_mask(endpoint_frame, display_class)
                for status in ("supported", "uncertain"):
                    mask = display_mask & endpoint_frame.e22_status.eq(status).to_numpy(bool)
                    rows.append(
                        {
                            "analysis_role": ANALYSIS_ROLE,
                            "design": design,
                            "endpoint": endpoint,
                            "display_class": display_class,
                            "status": status,
                            **cell_metrics(endpoint_frame, weights, mask),
                            "bootstrap_method": "two_stage",
                            "bootstrap_repetitions": BOOT_REPETITIONS,
                            "bootstrap_seed": BOOT_SEED,
                        }
                    )
    return pd.DataFrame(rows)


def dataset_state_differences(
    frame: pd.DataFrame,
    category: str,
) -> tuple[np.ndarray, int, int]:
    selected = frame.loc[class_mask(frame, category)]
    values = []
    n_supported = 0
    n_uncertain = 0
    for _, group in selected.groupby("target", sort=True):
        supported = group.loc[group.e22_status == "supported"]
        uncertain = group.loc[group.e22_status == "uncertain"]
        if (
            len(supported) < MIN_STATE_PREDICTIONS
            or len(uncertain) < MIN_STATE_PREDICTIONS
        ):
            continue
        values.append(
            float(supported.agree.mean()) - float(uncertain.agree.mean())
        )
        n_supported += len(supported)
        n_uncertain += len(uncertain)
    return np.asarray(values, dtype=float), n_supported, n_uncertain


def reliability_tests(designs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for design in DESIGNS:
        start = len(rows)
        frame = designs[design]
        for endpoint in ENDPOINTS:
            endpoint_frame = frame.loc[frame.endpoint == endpoint].reset_index(drop=True)
            for category in TEST_CATEGORIES:
                values, n_supported, n_uncertain = dataset_state_differences(
                    endpoint_frame,
                    category,
                )
                detail = wilcoxon_details(values)
                rows.append(
                    {
                        "analysis_role": ANALYSIS_ROLE,
                        "design": design,
                        "endpoint": endpoint,
                        "category": category,
                        "n_datasets": detail["n_units"],
                        "n_supported_predictions": n_supported,
                        "n_uncertain_predictions": n_uncertain,
                        "positive_datasets": detail["positive_units"],
                        "zero_datasets": detail["zero_units"],
                        "negative_datasets": detail["negative_units"],
                        "median_difference": detail["median_difference"],
                        "two_sided_p": detail["two_sided_p"],
                        "wilcoxon_method": detail["wilcoxon_method"],
                        "method_reason": detail["method_reason"],
                    }
                )
        family = rows[start:]
        for row, adjusted in zip(family, holm_adjust([r["two_sided_p"] for r in family])):
            row["holm_two_sided_p"] = float(adjusted)
    return pd.DataFrame(rows)


def bias_dataset_counts(
    designs: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    q4 = pd.read_csv(E15 / "q4_calibration_dataset_level.csv")
    q4 = q4.loc[q4.row_type == "summary"]
    rows = []
    for design in ("entry", "time_2015"):
        start = len(rows)
        for endpoint in ENDPOINTS:
            endpoint_frame = designs[design].loc[
                designs[design].endpoint == endpoint
            ].reset_index(drop=True)
            selected_full = None
            for display_class in BIAS_OUTPUT_CATEGORIES:
                selected = endpoint_frame.loc[class_mask(endpoint_frame, display_class)]
                if display_class == "inactive":
                    selected_full = selected
                biases = (
                    selected.groupby("target")
                    .apply(
                        lambda group: float(
                            (group.p_displayed - group.agree).mean()
                        ),
                        include_groups=False,
                    )
                    .loc[lambda values: selected.groupby("target").size() >= MIN_BIAS_PREDICTIONS]
                )
                values = biases.to_numpy(float)
                match = q4.loc[
                    q4.design.eq(design)
                    & q4.endpoint.eq(endpoint)
                    & q4.display_class.eq(display_class)
                ]
                if len(match) != 1:
                    raise RuntimeError(
                        f"E15 Q4 summary row is not unique: {design} {endpoint} {display_class}"
                    )
                p_value = float(match.iloc[0].dataset_bias_two_sided_p_min10)
                row = {
                    "analysis_role": ANALYSIS_ROLE,
                    "design": design,
                    "endpoint": endpoint,
                    "display_class": display_class,
                    "n_datasets": int(len(values)),
                    "positive_datasets": int((values > 0).sum()),
                    "zero_datasets": int((values == 0).sum()),
                    "negative_datasets": int((values < 0).sum()),
                    "bias_min": float(values.min()) if len(values) else np.nan,
                    "bias_max": float(values.max()) if len(values) else np.nan,
                    "bias_median": float(np.median(values)) if len(values) else np.nan,
                    "two_sided_p": p_value,
                    "wilcoxon_method": str(match.iloc[0].dataset_wilcoxon_method),
                    "p_value_source": "output/E15/q4_calibration_dataset_level.csv",
                    "holm_family": BIAS_HOLM_FAMILY,
                }
                if design == "entry" and display_class == "inactive":
                    re = random_effect_details(selected_full, 20)
                    row.update(
                        {
                            "random_effects_n_datasets": re["n_datasets"],
                            "random_effects_bias": re["random_estimate"],
                            "random_effects_ci_low": re["ci_low"],
                            "random_effects_ci_high": re["ci_high"],
                            "random_effects_i_squared": re["i_squared"],
                            "random_effects_pi_low": re["pi_low"],
                            "random_effects_pi_high": re["pi_high"],
                        }
                    )
                rows.append(row)
        family = [
            row for row in rows[start:]
            if row["display_class"] in BIAS_CATEGORIES
        ]
        for row, adjusted in zip(family, holm_adjust([r["two_sided_p"] for r in family])):
            row["holm_two_sided_p"] = float(adjusted)
    return pd.DataFrame(rows)


def dataset_bias_test(selected: pd.DataFrame) -> dict[str, float | int | str]:
    sizes = selected.groupby("target").size()
    eligible_targets = sizes.index[sizes >= MIN_BIAS_PREDICTIONS]
    eligible = selected.loc[selected.target.isin(eligible_targets)]
    biases = (
        eligible.groupby("target")
        .apply(
            lambda group: float((group.p_displayed - group.agree).mean()),
            include_groups=False,
        )
        .to_numpy(float)
    )
    detail = wilcoxon_details(biases)
    return {
        "bias_dataset_n": detail["n_units"],
        "bias_positive_datasets": detail["positive_units"],
        "bias_zero_datasets": detail["zero_units"],
        "bias_negative_datasets": detail["negative_units"],
        "bias_dataset_median": detail["median_difference"],
        "bias_dataset_two_sided_p": detail["two_sided_p"],
        "bias_dataset_wilcoxon_method": detail["wilcoxon_method"],
    }


def e11_calibration(designs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    e11 = designs["E11"]
    rows = []
    for scope, excluded in (("all", set()), ("without_CHEMBL4792", {"CHEMBL4792"})):
        scoped = e11.loc[~e11.target.isin(excluded)].reset_index(drop=True)
        for endpoint in ENDPOINTS:
            endpoint_frame = scoped.loc[scoped.endpoint == endpoint].reset_index(drop=True)
            weights = two_stage_weights(
                endpoint_frame,
                seed=BOOT_SEED,
                repetitions=BOOT_REPETITIONS,
            )
            for display_class in (*CLASSES, "all"):
                mask = class_mask(endpoint_frame, display_class)
                selected = endpoint_frame.loc[mask]
                rows.append(
                    {
                        "analysis_role": ANALYSIS_ROLE,
                        "scope": scope,
                        "endpoint": endpoint,
                        "display_class": display_class,
                        **cell_metrics(endpoint_frame, weights, mask),
                        **dataset_bias_test(selected),
                        "bootstrap_method": "two_stage",
                        "bootstrap_repetitions": BOOT_REPETITIONS,
                        "bootstrap_seed": BOOT_SEED,
                    }
                )
    return pd.DataFrame(rows)


def one_sided_wilcoxon_details(values) -> dict:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    nonzero = np.abs(values[values != 0])
    ties = len(np.unique(nonzero)) < len(nonzero)
    method = "exact" if len(values) <= 25 and not ties else "approx"
    if not len(values) or np.allclose(values, 0):
        return {
            "n_datasets": len(values),
            "positive_datasets": int((values > 0).sum()),
            "zero_datasets": int((values == 0).sum()),
            "negative_datasets": int((values < 0).sum()),
            "median_difference": float(np.median(values)) if len(values) else np.nan,
            "one_sided_p": 1.0,
            "wilcoxon_method": "not_applicable",
        }
    result = stats.wilcoxon(
        values,
        alternative="greater",
        zero_method="wilcox",
        method=method,
    )
    return {
        "n_datasets": len(values),
        "positive_datasets": int((values > 0).sum()),
        "zero_datasets": int((values == 0).sum()),
        "negative_datasets": int((values < 0).sum()),
        "median_difference": float(np.median(values)),
        "one_sided_p": float(result.pvalue),
        "wilcoxon_method": method,
    }


def e11_default_v2(designs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    e11 = designs["E11"]
    rows = []
    for endpoint in ENDPOINTS:
        endpoint_frame = e11.loc[e11.endpoint == endpoint]
        values = []
        for _, group in endpoint_frame.groupby("target", sort=True):
            supported = group.loc[group.e22_status == "supported"]
            uncertain = group.loc[group.e22_status == "uncertain"]
            if (
                len(supported) < MIN_STATE_PREDICTIONS
                or len(uncertain) < MIN_STATE_PREDICTIONS
            ):
                continue
            values.append(float(supported.agree.mean()) - float(uncertain.agree.mean()))
        rows.append(
            {
                "analysis_role": ANALYSIS_ROLE,
                "endpoint": endpoint,
                "alternative": "greater",
                **one_sided_wilcoxon_details(values),
            }
        )
    return pd.DataFrame(rows)


def e11_document_overlap() -> pd.DataFrame:
    inclusion = pd.read_csv(E11 / "target_inclusion.csv")
    included = inclusion.loc[inclusion.included, ["target", "endpoint", "retained_records"]]
    if len(included) != 23:
        raise RuntimeError(f"Expected 23 E11 datasets, found {len(included)}")
    evaluation_records = load_e2_records()
    evaluation_documents = {
        str(document)
        for records in evaluation_records.values()
        for document in records.document_chembl_id.dropna().unique()
    }
    documents_by_key = {}
    for row in included.itertuples(index=False):
        raw = pd.read_csv(
            ACTIVITIES / f"{row.target}_{row.endpoint}.csv",
            low_memory=False,
        )
        curated = curate_activity_table(raw, row.target, row.endpoint).retained_records
        found = len(curated)
        expected = int(row.retained_records)
        if found != expected:
            raise RuntimeError(
                f"Retained record mismatch for {row.target} {row.endpoint}: "
                f"{found} != {expected}"
            )
        documents_by_key[(row.target, row.endpoint)] = set(
            curated.document_chembl_id.dropna().astype(str)
        )
    rows = []
    for endpoint in ENDPOINTS:
        keys = [key for key in documents_by_key if key[1] == endpoint]
        documents = set().union(*(documents_by_key[key] for key in keys))
        rows.append(
            {
                "analysis_role": ANALYSIS_ROLE,
                "endpoint": endpoint,
                "n_datasets": len(keys),
                "e11_documents": len(documents),
                "evaluation_documents": len(evaluation_documents),
                "shared_documents": len(documents & evaluation_documents),
                "shared_fraction_of_e11": (
                    len(documents & evaluation_documents) / len(documents)
                    if documents
                    else np.nan
                ),
            }
        )
    all_documents = set().union(*documents_by_key.values())
    rows.append(
        {
            "analysis_role": ANALYSIS_ROLE,
            "endpoint": "combined",
            "n_datasets": len(documents_by_key),
            "e11_documents": len(all_documents),
            "evaluation_documents": len(evaluation_documents),
            "shared_documents": len(all_documents & evaluation_documents),
            "shared_fraction_of_e11": (
                len(all_documents & evaluation_documents) / len(all_documents)
                if all_documents
                else np.nan
            ),
        }
    )
    return pd.DataFrame(rows)


def time_split_auc_summary() -> pd.DataFrame:
    time_split = pd.read_csv(E2 / "time_split_by_target.csv")
    time_split = time_split.loc[
        time_split.cutoff_year.eq(2015)
        & time_split.eligible.astype(bool)
        & time_split.g_auc.notna()
    ].copy()
    inclusion = pd.read_csv(E2 / "target_inclusion.csv")
    evaluation_keys = set(
        zip(
            inclusion.loc[inclusion.included, "target"].astype(str),
            inclusion.loc[inclusion.included, "endpoint"].astype(str),
        )
    )
    time_split["in_evaluation_52"] = [
        (str(target), str(endpoint)) in evaluation_keys
        for target, endpoint in zip(time_split.target, time_split.endpoint)
    ]
    rows = []
    for scope, selected in (
        ("all", time_split),
        ("evaluation_52", time_split.loc[time_split.in_evaluation_52]),
    ):
        for endpoint in ENDPOINTS:
            values = selected.loc[selected.endpoint == endpoint, "g_auc"].to_numpy(float)
            rows.append(
                {
                    "analysis_role": ANALYSIS_ROLE,
                    "scope": scope,
                    "endpoint": endpoint,
                    "n_datasets": len(values),
                    "auc_median": float(np.median(values)) if len(values) else np.nan,
                    "auc_q1": float(np.percentile(values, 25)) if len(values) else np.nan,
                    "auc_q3": float(np.percentile(values, 75)) if len(values) else np.nan,
                    "cutoff_year": 2015,
                }
            )
    return pd.DataFrame(rows)


def row_for(
    frame: pd.DataFrame,
    **conditions,
) -> pd.Series:
    selected = frame
    for column, value in conditions.items():
        selected = selected.loc[selected[column] == value]
    if len(selected) != 1:
        raise RuntimeError(f"Expected one row for {conditions}, found {len(selected)}")
    return selected.iloc[0]


def validate_expected(
    reliability: pd.DataFrame,
    tests: pd.DataFrame,
    bias_counts: pd.DataFrame,
    e11_calibration_rows: pd.DataFrame,
    e11_v2: pd.DataFrame,
    overlap: pd.DataFrame,
) -> dict:
    e9 = pd.read_csv(E9 / "calibration_summary.csv")
    for status in ("supported", "uncertain"):
        actual = row_for(
            reliability,
            design="entry",
            endpoint="IC50",
            display_class="all",
            status=status,
        )
        expected = row_for(
            e9,
            design="entry",
            endpoint="IC50",
            population=f"G_{status}",
            display_class="all",
        )
        for column in (
            "n",
            "agreement_rate",
            "agreement_ci_low",
            "agreement_ci_high",
        ):
            actual_column = {
                "n": "n_predictions",
                "agreement_rate": "agreement_rate",
                "agreement_ci_low": "agreement_ci_low",
                "agreement_ci_high": "agreement_ci_high",
            }[column]
            if not np.isclose(
                float(actual[actual_column]),
                float(expected[column]),
                atol=1e-12,
                rtol=0,
            ):
                raise RuntimeError(f"E9 aggregate mismatch: {status} {column}")
    for endpoint in ENDPOINTS:
        actual = row_for(
            reliability,
            design="entry",
            endpoint=endpoint,
            display_class="all",
            status="supported",
        )
        expected = row_for(
            e9,
            design="entry",
            endpoint=endpoint,
            population="G_supported",
            display_class="all",
        )
        if not np.isclose(
            float(actual.agreement_ci_low),
            float(expected.agreement_ci_low),
            atol=1e-12,
            rtol=0,
        ) or not np.isclose(
            float(actual.agreement_ci_high),
            float(expected.agreement_ci_high),
            atol=1e-12,
            rtol=0,
        ):
            raise RuntimeError(f"E9 entry supported CI mismatch: {endpoint}")

    entry_ic50_supported_inactive = row_for(
        reliability,
        design="entry",
        endpoint="IC50",
        display_class="inactive",
        status="supported",
    )
    if int(entry_ic50_supported_inactive.n_predictions) != 615:
        raise RuntimeError("Entry IC50 supported-inactive count mismatch")
    time_ki_supported_inactive = row_for(
        reliability,
        design="time_2015",
        endpoint="Ki",
        display_class="inactive",
        status="supported",
    )
    if int(time_ki_supported_inactive.n_predictions) != 190:
        raise RuntimeError("Time Ki supported-inactive count mismatch")

    e11_ic50_inactive = row_for(
        e11_calibration_rows,
        scope="all",
        endpoint="IC50",
        display_class="inactive",
    )
    e11_ki_inactive = row_for(
        e11_calibration_rows,
        scope="all",
        endpoint="Ki",
        display_class="inactive",
    )
    e11_ki_reduced = row_for(
        e11_calibration_rows,
        scope="without_CHEMBL4792",
        endpoint="Ki",
        display_class="inactive",
    )
    checks = {
        "e11_ic50_inactive_n": int(e11_ic50_inactive.n_predictions),
        "e11_ic50_inactive_positive_datasets": int(
            e11_ic50_inactive.bias_positive_datasets
        ),
        "e11_ki_inactive_n": int(e11_ki_inactive.n_predictions),
        "e11_ki_without_4792_n": int(e11_ki_reduced.n_predictions),
        "v2_ic50_datasets": int(row_for(e11_v2, endpoint="IC50").n_datasets),
        "v2_ic50_positive": int(row_for(e11_v2, endpoint="IC50").positive_datasets),
        "v2_ki_datasets": int(row_for(e11_v2, endpoint="Ki").n_datasets),
        "v2_ki_positive": int(row_for(e11_v2, endpoint="Ki").positive_datasets),
        "overlap_ic50": row_for(overlap, endpoint="IC50").shared_documents,
        "overlap_ic50_n": row_for(overlap, endpoint="IC50").e11_documents,
        "overlap_ki": row_for(overlap, endpoint="Ki").shared_documents,
        "overlap_ki_n": row_for(overlap, endpoint="Ki").e11_documents,
    }
    expected = {
        "e11_ic50_inactive_n": 447,
        "e11_ic50_inactive_positive_datasets": 7,
        "e11_ki_inactive_n": 500,
        "e11_ki_without_4792_n": 202,
        "v2_ic50_datasets": 12,
        "v2_ic50_positive": 11,
        "v2_ki_datasets": 11,
        "v2_ki_positive": 11,
        "overlap_ic50": 621,
        "overlap_ic50_n": 2353,
        "overlap_ki": 1165,
        "overlap_ki_n": 1629,
    }
    for name, value in expected.items():
        if checks[name] != value:
            raise RuntimeError(f"E22 expected value mismatch {name}: {checks[name]} != {value}")
    if not (0.03 <= float(row_for(
        bias_counts,
        design="entry",
        endpoint="IC50",
        display_class="inactive",
    ).holm_two_sided_p) <= 0.05):
        raise RuntimeError("Entry IC50 inactive Holm result is outside the expected range")
    medium_ki = row_for(
        bias_counts,
        design="entry",
        endpoint="Ki",
        display_class="medium",
    )
    if not (0.0001 <= float(medium_ki.two_sided_p) <= 0.0002):
        raise RuntimeError("Entry Ki medium raw P is outside the expected range")
    for design in ("entry", "time_2015"):
        design_rows = bias_counts.loc[bias_counts.design == design]
        if len(design_rows) != 8:
            raise RuntimeError(
                f"Expected 8 bias count rows for {design}, found {len(design_rows)}"
            )
        family_rows = design_rows.loc[
            design_rows.display_class.isin(BIAS_CATEGORIES)
        ]
        all_rows = design_rows.loc[design_rows.display_class == "all"]
        if len(family_rows) != 6 or family_rows.holm_two_sided_p.isna().any():
            raise RuntimeError(f"Incorrect Holm family for {design}")
        if len(all_rows) != 2 or all_rows.holm_two_sided_p.notna().any():
            raise RuntimeError(f"all rows must remain outside the Holm family for {design}")
        if not design_rows.holm_family.eq(BIAS_HOLM_FAMILY).all():
            raise RuntimeError(f"Incorrect Holm family description for {design}")
    if not (len(tests) == 30):
        raise RuntimeError(f"Expected 30 class reliability tests, found {len(tests)}")
    return checks


def main() -> None:
    started_at = datetime.now().astimezone().isoformat()
    designs = load_designs()
    reliability = reliability_by_design(designs)
    tests = reliability_tests(designs)
    bias_counts = bias_dataset_counts(designs)
    e11_calibration_rows = e11_calibration(designs)
    e11_v2 = e11_default_v2(designs)
    overlap = e11_document_overlap()
    auc = time_split_auc_summary()
    composition = bin_composition(designs["entry"])
    validate_bin_composition(composition)
    checks = validate_expected(
        reliability,
        tests,
        bias_counts,
        e11_calibration_rows,
        e11_v2,
        overlap,
    )
    OUT.mkdir(parents=True, exist_ok=True)
    outputs = {
        "class_reliability_by_design.csv": reliability,
        "class_reliability_tests.csv": tests,
        "class_bias_dataset_counts.csv": bias_counts,
        "e11_class_calibration.csv": e11_calibration_rows,
        "e11_v2_default_display.csv": e11_v2,
        "e11_document_overlap_analyzed.csv": overlap,
        "time_split_auc_by_dataset_summary.csv": auc,
        "fig4a_bin_composition.csv": composition,
    }
    for name, frame in outputs.items():
        frame.to_csv(OUT / name, index=False)
    activity_inputs = [
        ACTIVITIES / f"{target}_{endpoint}.csv"
        for target, endpoint in pd.read_csv(E11 / "target_inclusion.csv")
        .loc[lambda frame: frame.included, ["target", "endpoint"]]
        .itertuples(index=False, name=None)
    ]
    write_manifest(
        OUT,
        stage="E22",
        inputs=(
            E5C / "predictions_e5c.csv",
            E11 / "predictions_e11.csv",
            E11 / "target_inclusion.csv",
            E2 / "time_split_by_target.csv",
            E2 / "target_inclusion.csv",
            E15 / "q4_calibration_dataset_level.csv",
            E9 / "calibration_summary.csv",
            DATA / "download_manifest.json",
            DATA / "selected_targets.csv",
            *activity_inputs,
        ),
        sources=(
            Path(__file__),
            S / "source/test_e22.py",
            S / "source/e7_analysis.py",
            S / "source/e15_leakfree.py",
            S / "source/e9_analysis.py",
            S / "source/e2_models.py",
            S / "source/manifest_utils.py",
            ROOT / "src/chemgridmap/evidence.py",
        ),
        parameters={
            "analysis_role": ANALYSIS_ROLE,
            "bootstrap_method": "two_stage",
            "bootstrap_repetitions": BOOT_REPETITIONS,
            "bootstrap_seed": BOOT_SEED,
            "support_threshold": SUPPORT_THRESHOLD,
            "minimum_state_predictions": MIN_STATE_PREDICTIONS,
            "minimum_bias_dataset_predictions": MIN_BIAS_PREDICTIONS,
            "reliability_test_categories": TEST_CATEGORIES,
            "bias_count_categories": BIAS_OUTPUT_CATEGORIES,
            "bias_holm_categories": BIAS_CATEGORIES,
            "bias_holm_family": BIAS_HOLM_FAMILY,
            "fig4a_probability_bins": "five equal-width bins per displayed class",
            "fig4a_bin_rule": "clip(floor(p_displayed * 5 + 1e-12), 0, 4)",
            "fig4a_largest_dataset": "dataset with the most predictions in the bin",
        },
        started_at=started_at,
        finished_at=datetime.now().astimezone().isoformat(),
    )
    print(
        f"E22 outputs={len(outputs)} reliability_rows={len(reliability)} "
        f"checks={checks}",
        flush=True,
    )


if __name__ == "__main__":
    main()
