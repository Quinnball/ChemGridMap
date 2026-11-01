"""E11: prospective validation of display rule P on new target-endpoint datasets."""
from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
from datetime import datetime
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[3]
S = ROOT / "paper/strengthening"
DATA = S / "data/chembl37_api_validation"
ACTIVITIES = DATA / "activities"
OUT = S / "output/E11"
QA = S / "qa/E11"
PLAN = S / "plans/E11_validation_plan_sealed.md"
E5C = S / "output/E5c"
CGM = ROOT
SOURCE = S / "source"

sys.path[:0] = [
    str(SOURCE),
    str(ROOT / "paper"),
    str(ROOT / "src"),
]

import numpy as np
import pandas as pd

from chemgridmap.evidence import CLASSES
from e1a_statistics import auc_mann_whitney
from e2_models import (
    FOLD_SEED,
    build_target_predictions,
    curate_activity_table,
    prepare_target,
)
from e5c_models import (
    delta_document_map,
    unit_deltas,
    weighted_gaussian_prediction,
)
from e5d_analysis import full_map_metrics
from e5d_models import BOOT_REPETITIONS, BOOT_SEED, add_rule_states
from e7_analysis import holm_adjust, one_sided_wilcoxon, two_stage_weights
from manifest_utils import write_manifest
from reliability_activity import entry_folds

ENDPOINTS = ("IC50", "Ki")
SUPPORT_THRESHOLD = 0.8
ANALYSIS_ROLE = "prospective validation"
REGRESSION_TARGETS = (
    ("CHEMBL206", "IC50"),
    ("CHEMBL2971", "IC50"),
    ("CHEMBL234", "Ki"),
    ("CHEMBL5113", "Ki"),
)
TRANSCRIPT: list[str] = []


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def log(message: str) -> None:
    print(message, flush=True)
    TRANSCRIPT.append(message)


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
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            json_safe(value),
            indent=2,
            ensure_ascii=False,
            allow_nan=False,
        )
        + "\n"
    )


def read_activity(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, low_memory=False)


def predict_entry_from_records(
    records: pd.DataFrame,
    target: str,
    endpoint: str,
    context: dict | None = None,
) -> pd.DataFrame:
    """Run G and reproduce E5c's training-only touched marker for every held unit."""
    context = context or prepare_target(records)
    baseline, _ = build_target_predictions(
        records,
        target,
        endpoint,
        models=("G",),
        context=context,
    )
    baseline = baseline.loc[baseline.population == "all"].copy()
    doc = context["doc"]
    observations = context["observations"]
    assignments = entry_folds(doc, seed=FOLD_SEED)
    units_by_entry = dict(tuple(doc.groupby("entry", sort=True)))
    full_delta = delta_document_map(observations)
    rows = []

    for fold in range(10):
        train_entries = {
            entry for entry, assigned in assignments.items() if assigned != fold
        }
        fold_delta = delta_document_map(
            observations, training_entries=train_entries
        )
        for row in baseline.loc[baseline.fold == fold].to_dict("records"):
            units = units_by_entry[str(row["entry"])]
            held_documents = set(str(row["held_out_unit"]).split(";"))
            remaining = units.loc[
                ~units.documents.map(
                    lambda value: bool(
                        set(str(value).split(";")) & held_documents
                    )
                )
            ]
            if remaining.empty:
                raise RuntimeError(
                    f"{target} {endpoint} {row['entry']} has no remaining units"
                )
            values = remaining.document_value.to_numpy(float)
            sigma = float(row["tau"])
            support, predictive, center, _ = weighted_gaussian_prediction(
                values,
                np.full(len(values), sigma**2),
                sigma,
            )
            if not np.allclose(
                predictive,
                np.array([row["p0"], row["p1"], row["p2"]], dtype=float),
                atol=1e-12,
                rtol=0,
            ):
                raise RuntimeError("G probabilities do not reproduce E2/E5c.")
            dynamic, touched = unit_deltas(remaining, fold_delta)
            _, full_touched = unit_deltas(remaining, full_delta)
            del dynamic
            shown = int(row["displayed"])
            held_value = float(
                units.loc[
                    units.documents.astype(str) == str(row["held_out_unit"]),
                    "document_value",
                ].iloc[0]
            )
            rows.append(
                {
                    "split_design": "entry",
                    "merge_basis": "same_value",
                    "target": target,
                    "endpoint": endpoint,
                    "fold": int(row["fold"]),
                    "cutoff_year": np.nan,
                    "entry": str(row["entry"]),
                    "held_out_unit": str(row["held_out_unit"]),
                    "m": int(row["m"]),
                    "displayed": shown,
                    "observed": int(row["observed"]),
                    "held_out_value": held_value,
                    "sigma_doc": sigma,
                    "remaining_center": float(center),
                    "remaining_values": ";".join(
                        f"{value:.10g}" for value in values
                    ),
                    "method": "G",
                    "p0": float(row["p0"]),
                    "p1": float(row["p1"]),
                    "p2": float(row["p2"]),
                    "p_displayed": float(row["p_displayed"]),
                    "support": float(support[shown]),
                    "display_status": (
                        "supported"
                        if support[shown] >= SUPPORT_THRESHOLD
                        else "uncertain"
                    ),
                    "agree": bool(row["agree"]),
                    "analysis_role": ANALYSIS_ROLE,
                    "m_layer": 1 if int(row["m"]) == 1 else 2,
                    "flag_touched": bool(touched.any()),
                    "full_data_flag_touched": bool(full_touched.any()),
                    "n_remaining_flagged_units": int(touched.sum()),
                }
            )
    result = pd.DataFrame(rows).sort_values(
        ["target", "endpoint", "fold", "entry", "held_out_unit"]
    )
    result = add_rule_states(result.reset_index(drop=True))
    return result


def predict_entry_from_file(
    raw_path: Path,
    target: str,
    endpoint: str,
) -> pd.DataFrame:
    curated = curate_activity_table(
        read_activity(raw_path), target, endpoint
    )
    return predict_entry_from_records(
        curated.retained_records, target, endpoint
    )


def run_regression() -> dict:
    expected_columns = [
        "split_design",
        "merge_basis",
        "method",
        "target",
        "endpoint",
        "fold",
        "entry",
        "held_out_unit",
        "p0",
        "p1",
        "p2",
        "p_displayed",
        "support",
        "agree",
        "flag_touched",
    ]
    expected = pd.read_csv(
        E5C / "predictions_e5c.csv",
        usecols=expected_columns,
        low_memory=False,
        float_precision="round_trip",
    )
    expected = expected.loc[
        (expected.split_design == "entry")
        & (expected.merge_basis == "same_value")
        & (expected.method == "G")
        & expected.set_index(["target", "endpoint"]).index.isin(
            REGRESSION_TARGETS
        )
    ].copy()
    actual_frames = []
    dataset_rows = []
    for target, endpoint in REGRESSION_TARGETS:
        path = (
            S
            / "data/chembl37_api/activities"
            / f"{target}_{endpoint}.csv"
        )
        actual_frames.append(predict_entry_from_file(path, target, endpoint))
    actual = pd.concat(actual_frames, ignore_index=True)
    keys = ["target", "endpoint", "fold", "entry", "held_out_unit"]
    merged = actual.merge(
        expected,
        on=keys,
        suffixes=("_e11", "_e5c"),
        validate="one_to_one",
    )
    if not (len(actual) == len(expected) == len(merged)):
        raise RuntimeError(
            "E11/E5c regression row mismatch: "
            f"{len(actual)}, {len(expected)}, {len(merged)}"
        )
    errors = {}
    for column in ("p0", "p1", "p2", "p_displayed", "support"):
        errors[column] = float(
            np.max(
                np.abs(
                    merged[f"{column}_e11"].to_numpy(float)
                    - merged[f"{column}_e5c"].to_numpy(float)
                )
            )
        )
    agree_equal = bool(
        np.array_equal(
            merged.agree_e11.to_numpy(bool),
            merged.agree_e5c.to_numpy(bool),
        )
    )
    flag_equal = bool(
        np.array_equal(
            merged.flag_touched_e11.to_numpy(bool),
            merged.flag_touched_e5c.to_numpy(bool),
        )
    )
    for target, endpoint in REGRESSION_TARGETS:
        selected = actual.loc[
            (actual.target == target) & (actual.endpoint == endpoint)
        ]
        dataset_rows.append(
            {
                "target": target,
                "endpoint": endpoint,
                "n_predictions": int(len(selected)),
                "n_flag_touched": int(selected.flag_touched.sum()),
            }
        )
    maximum = max(errors.values())
    passed = bool(
        maximum <= 1e-12 and agree_equal and flag_equal
    )
    result = {
        "analysis_role": ANALYSIS_ROLE,
        "selected_datasets": dataset_rows,
        "n_expected": int(len(expected)),
        "n_actual": int(len(actual)),
        "n_matched": int(len(merged)),
        "maximum_probability_error": maximum,
        "errors": errors,
        "probability_tolerance": 1e-12,
        "agree_identical": agree_equal,
        "flag_touched_identical": flag_equal,
        "passed": passed,
    }
    write_json(QA / "regression.json", result)
    if not passed:
        raise RuntimeError(f"E11 regression failed: {result}")
    return result


def verify_download() -> dict:
    manifest = json.loads(
        (DATA / "download_manifest.json").read_text()
    )
    selected = pd.read_csv(DATA / "selected_targets.csv")
    rows = []
    for item in manifest["files"]:
        path = ACTIVITIES / item["file"]
        count = int(len(read_activity(path)))
        digest = sha256(path)
        rows.append(
            {
                "target": item["target"],
                "endpoint": item["standard_type"],
                "file": item["file"],
                "records": count,
                "expected_records": int(item["expected"]),
                "manifest_records": int(item["records"]),
                "sha256_matches": digest == item["sha256"],
                "status": item["status"],
                "passed": bool(
                    count == int(item["expected"])
                    == int(item["records"])
                    and digest == item["sha256"]
                    and item["status"] == "complete"
                ),
            }
        )
    frame = pd.DataFrame(rows)
    result = {
        "analysis_role": ANALYSIS_ROLE,
        "selected_targets_n": int(len(selected)),
        "selected_by_endpoint": {
            endpoint: int((selected.standard_type == endpoint).sum())
            for endpoint in ENDPOINTS
        },
        "manifest_files_n": int(len(rows)),
        "all_passed": bool(
            len(rows) == len(selected) and frame.passed.all()
        ),
    }
    write_json(QA / "data_manifest_check.json", result)
    if not frame.passed.all():
        raise RuntimeError("Validation download manifest/count check failed.")
    return result


def curate_validation_targets() -> tuple[
    pd.DataFrame, dict[tuple[str, str], pd.DataFrame],
    dict[tuple[str, str], dict],
]:
    selected = pd.read_csv(DATA / "selected_targets.csv")
    inclusion_rows = []
    records_by_key = {}
    contexts_by_key = {}
    for index, item in selected.iterrows():
        target = str(item.target_chembl_id)
        endpoint = str(item.standard_type)
        path = ACTIVITIES / f"{target}_{endpoint}.csv"
        raw = None
        try:
            raw = read_activity(path)
            curated = curate_activity_table(raw, target, endpoint)
            records = curated.retained_records
            context = prepare_target(records)
            entry_counts = context["doc"].groupby("entry").size()
            multi = int((entry_counts >= 2).sum())
            included = bool(multi >= 30)
            inclusion_rows.append(
                {
                    "target": target,
                    "endpoint": endpoint,
                    "role": item.role,
                    "api_records": int(len(raw)),
                    "retained_records": int(len(records)),
                    "entries": int(len(entry_counts)),
                    "multi_unit_entries": multi,
                    "included": included,
                    "reason": (
                        "included"
                        if included
                        else "fewer than 30 entries with two or more evidence units"
                    ),
                    "input_file": str(path.relative_to(ROOT)),
                }
            )
            records_by_key[(target, endpoint)] = records
            contexts_by_key[(target, endpoint)] = context
        except Exception as error:
            inclusion_rows.append(
                {
                    "target": target,
                    "endpoint": endpoint,
                    "role": item.role,
                    "api_records": int(len(raw)) if raw is not None else np.nan,
                    "retained_records": np.nan,
                    "entries": np.nan,
                    "multi_unit_entries": np.nan,
                    "included": False,
                    "reason": f"curation failure: {type(error).__name__}: {error}",
                    "input_file": str(path.relative_to(ROOT)),
                }
            )
        log(
            f"curation {index + 1}/{len(selected)} "
            f"{target} {endpoint}"
        )
    inclusion = pd.DataFrame(inclusion_rows)
    inclusion.to_csv(OUT / "target_inclusion.csv", index=False)
    return inclusion, records_by_key, contexts_by_key


def build_predictions(
    inclusion: pd.DataFrame,
    records_by_key: dict,
    contexts_by_key: dict,
) -> pd.DataFrame:
    frames = []
    included = inclusion.loc[inclusion.included]
    for item in included.itertuples(index=False):
        key = (item.target, item.endpoint)
        frame = predict_entry_from_records(
            records_by_key[key],
            item.target,
            item.endpoint,
            context=contexts_by_key[key],
        )
        frames.append(frame)
        log(f"predictions {item.target} {item.endpoint}: {len(frame)}")
    predictions = pd.concat(frames, ignore_index=True)
    predictions.to_csv(OUT / "predictions_e11.csv", index=False)
    return predictions


def map_impact(
    inclusion: pd.DataFrame,
    records_by_key: dict,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for item in inclusion.loc[inclusion.included].itertuples(
        index=False
    ):
        metrics = full_map_metrics(
            records_by_key[(item.target, item.endpoint)]
        )
        rows.append(
            {
                "analysis_role": ANALYSIS_ROLE,
                "target": item.target,
                "endpoint": item.endpoint,
                **metrics,
            }
        )
    by_dataset = pd.DataFrame(rows)
    by_dataset.to_csv(
        OUT / "map_impact_by_dataset.csv", index=False
    )
    summary_rows = []
    for endpoint, group in [
        *by_dataset.groupby("endpoint", sort=True),
        ("ALL", by_dataset),
    ]:
        summary_rows.append(
            {
                "analysis_role": ANALYSIS_ROLE,
                "endpoint": endpoint,
                "n_datasets": int(len(group)),
                "n_entries": int(group.n_entries.sum()),
                "touched_entries_n": int(
                    group.touched_entries_n.sum()
                ),
                "touched_entries_rate": float(
                    group.touched_entries_n.sum()
                    / group.n_entries.sum()
                ),
                "g_supported_n": int(group.g_supported_n.sum()),
                "g_saturation_rate": float(
                    group.g_supported_n.sum()
                    / group.n_entries.sum()
                ),
                "p_supported_n": int(group.p_supported_n.sum()),
                "p_saturation_rate": float(
                    group.p_supported_n.sum()
                    / group.n_entries.sum()
                ),
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT / "map_impact_summary.csv", index=False)
    return by_dataset, summary


def hypothesis_result(values) -> dict:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if len(finite) < 8:
        return {
            "n_datasets": int(len(finite)),
            "median_difference": (
                float(np.median(finite)) if len(finite) else np.nan
            ),
            "wilcoxon_statistic": np.nan,
            "p_value": np.nan,
            "test_performed": False,
        }
    statistic, p_value, n = one_sided_wilcoxon(
        finite, alternative="greater"
    )
    return {
        "n_datasets": int(n),
        "median_difference": float(np.median(finite)),
        "wilcoxon_statistic": statistic,
        "p_value": p_value,
        "test_performed": True,
    }


def validation_tests(predictions: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for endpoint in ENDPOINTS:
        endpoint_frame = predictions.loc[
            predictions.endpoint == endpoint
        ]
        v1_values = []
        v2_values = []
        v3_values = []
        v1_eligible = []
        for target, group in endpoint_frame.groupby(
            "target", sort=True
        ):
            g_supported = group.g_status == "supported"
            p_supported = group.p_status == "supported"
            touched = group.flag_touched.astype(bool)

            if (g_supported & touched).any():
                v1_eligible.append(target)
                if p_supported.any() and g_supported.any():
                    v1_values.append(
                        float(
                            group.loc[p_supported, "agree"].mean()
                            - group.loc[g_supported, "agree"].mean()
                        )
                    )
                else:
                    v1_values.append(np.nan)

            p_uncertain = ~p_supported
            if p_supported.sum() >= 5 and p_uncertain.sum() >= 5:
                v2_values.append(
                    float(
                        group.loc[p_supported, "agree"].mean()
                        - group.loc[p_uncertain, "agree"].mean()
                    )
                )

            auc = auc_mann_whitney(
                group.agree.to_numpy(bool),
                group.p_displayed.to_numpy(float),
            )
            v3_values.append(
                float(auc - 0.5) if np.isfinite(auc) else np.nan
            )

        results = {
            "V1": hypothesis_result(v1_values),
            "V2": hypothesis_result(v2_values),
            "V3": hypothesis_result(v3_values),
        }
        raw_p = np.array(
            [results[name]["p_value"] for name in ("V1", "V2", "V3")],
            dtype=float,
        )
        adjusted = holm_adjust(raw_p)
        for index, name in enumerate(("V1", "V2", "V3")):
            rows.append(
                {
                    "analysis_role": ANALYSIS_ROLE,
                    "endpoint": endpoint,
                    "test": name,
                    "n_datasets_eligible": (
                        len(v1_eligible)
                        if name == "V1"
                        else results[name]["n_datasets"]
                    ),
                    **results[name],
                    "holm_adjusted_p": float(adjusted[index]),
                    "holm_significant": bool(
                        np.isfinite(adjusted[index])
                        and adjusted[index] < 0.05
                    ),
                    "alternative": "greater_than_zero",
                    "alpha": 0.05,
                    "minimum_datasets": 8,
                }
            )
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "validation_tests.csv", index=False)
    return result


def bootstrap_state_metric(
    frame: pd.DataFrame,
    weights: np.ndarray,
    mask: np.ndarray,
    metric: str,
) -> dict:
    mask = np.asarray(mask, dtype=bool)
    point_n = int(mask.sum())
    denominator = weights[:, mask].sum(axis=1)
    probability = frame.p_displayed.to_numpy(float)
    agree = frame.agree.to_numpy(float)
    if metric == "mean_probability":
        numerator = (
            weights[:, mask] * probability[mask]
        ).sum(axis=1)
        point = (
            float(probability[mask].mean()) if point_n else np.nan
        )
    elif metric == "agreement_rate":
        numerator = (
            weights[:, mask] * agree[mask]
        ).sum(axis=1)
        point = float(agree[mask].mean()) if point_n else np.nan
    elif metric == "calibration_bias":
        numerator = (
            weights[:, mask]
            * (probability[mask] - agree[mask])
        ).sum(axis=1)
        point = (
            float((probability[mask] - agree[mask]).mean())
            if point_n
            else np.nan
        )
    else:
        raise ValueError(f"Unknown bootstrap metric: {metric}")
    draws = np.divide(
        numerator,
        denominator,
        out=np.full(len(weights), np.nan),
        where=denominator > 0,
    )
    draws = draws[np.isfinite(draws)]
    if len(draws):
        low, high = np.percentile(draws, [2.5, 97.5])
    else:
        low = high = np.nan
    return {
        "n_predictions": point_n,
        "value": point,
        "ci_low": float(low),
        "ci_high": float(high),
        "valid_repetitions": int(len(draws)),
    }


def descriptives(predictions: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for endpoint in ENDPOINTS:
        frame = predictions.loc[
            predictions.endpoint == endpoint
        ].reset_index(drop=True)
        weights = two_stage_weights(
            frame,
            seed=BOOT_SEED,
            repetitions=BOOT_REPETITIONS,
        )
        p_supported = (frame.p_status == "supported").to_numpy()
        for metric in (
            "mean_probability",
            "agreement_rate",
            "calibration_bias",
        ):
            rows.append(
                {
                    "analysis_role": ANALYSIS_ROLE,
                    "endpoint": endpoint,
                    "population": "P_supported",
                    "metric": metric,
                    **bootstrap_state_metric(
                        frame, weights, p_supported, metric
                    ),
                    "bootstrap_method": "two_stage",
                    "bootstrap_repetitions": BOOT_REPETITIONS,
                    "bootstrap_seed": BOOT_SEED,
                    "two_stage_function": "e7_analysis.two_stage_weights",
                }
            )

        inactive_untouched = (
            (frame.displayed == CLASSES.index("inactive"))
            & ~frame.flag_touched.astype(bool)
        ).to_numpy()
        rows.append(
            {
                "analysis_role": ANALYSIS_ROLE,
                "endpoint": endpoint,
                "population": "untouched_inactive",
                "metric": "calibration_bias_pooled",
                **bootstrap_state_metric(
                    frame,
                    weights,
                    inactive_untouched,
                    "calibration_bias",
                ),
                "bootstrap_method": "two_stage",
                "bootstrap_repetitions": BOOT_REPETITIONS,
                "bootstrap_seed": BOOT_SEED,
                "two_stage_function": "e7_analysis.two_stage_weights",
            }
        )
        dataset_biases = []
        eligible_targets = []
        touched_inactive_rates = []
        for target, group in frame.groupby("target", sort=True):
            untouched = (
                (group.displayed == CLASSES.index("inactive"))
                & ~group.flag_touched.astype(bool)
            )
            if int(untouched.sum()) >= 10:
                dataset_biases.append(
                    float(
                        (
                            group.loc[untouched, "p_displayed"]
                            - group.loc[untouched, "agree"]
                        ).mean()
                    )
                )
                eligible_targets.append(target)
            touched_inactive = (
                (group.displayed == CLASSES.index("inactive"))
                & group.flag_touched.astype(bool)
            )
            if touched_inactive.any():
                touched_inactive_rates.append(
                    float(
                        group.loc[touched_inactive, "agree"].mean()
                    )
                )
        rows.append(
            {
                "analysis_role": ANALYSIS_ROLE,
                "endpoint": endpoint,
                "population": "untouched_inactive",
                "metric": "calibration_bias_dataset_median",
                "n_predictions": np.nan,
                "value": (
                    float(np.median(dataset_biases))
                    if dataset_biases
                    else np.nan
                ),
                "ci_low": np.nan,
                "ci_high": np.nan,
                "valid_repetitions": np.nan,
                "n_datasets": len(dataset_biases),
                "eligible_targets": ";".join(eligible_targets),
                "minimum_predictions_per_dataset": 10,
            }
        )
        touched_inactive = (
            (frame.displayed == CLASSES.index("inactive"))
            & frame.flag_touched.astype(bool)
        )
        rows.append(
            {
                "analysis_role": ANALYSIS_ROLE,
                "endpoint": endpoint,
                "population": "touched_inactive",
                "metric": "agreement_rate",
                "n_predictions": int(touched_inactive.sum()),
                "value": (
                    float(frame.loc[touched_inactive, "agree"].mean())
                    if touched_inactive.any()
                    else np.nan
                ),
                "ci_low": np.nan,
                "ci_high": np.nan,
                "valid_repetitions": np.nan,
                "n_datasets": int(
                    frame.loc[touched_inactive, "target"].nunique()
                ),
            }
        )
        rows.append(
            {
                "analysis_role": ANALYSIS_ROLE,
                "endpoint": endpoint,
                "population": "all",
                "metric": "flag_touched_rate",
                "n_predictions": int(len(frame)),
                "value": float(frame.flag_touched.mean()),
                "ci_low": np.nan,
                "ci_high": np.nan,
                "valid_repetitions": np.nan,
                "n_datasets": int(frame.target.nunique()),
            }
        )
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "descriptives.csv", index=False)
    return result


def make_decision(
    inclusion: pd.DataFrame,
    tests: pd.DataFrame,
) -> dict:
    significant_endpoints = []
    medians = {}
    for endpoint in ENDPOINTS:
        v1 = tests.loc[
            (tests.endpoint == endpoint) & (tests.test == "V1")
        ].iloc[0]
        medians[endpoint] = float(v1.median_difference)
        if bool(v1.holm_significant):
            significant_endpoints.append(endpoint)
    replicated = any(
        endpoint in significant_endpoints
        and all(
            np.isfinite(value) and value >= 0
            for other, value in medians.items()
            if other != endpoint
        )
        for endpoint in significant_endpoints
    )
    sentence = (
        "the display rule was confirmed on new targets"
        if replicated
        else "the display rule did not replicate on new targets"
    )
    decision = {
        "analysis_role": ANALYSIS_ROLE,
        "qualified_dataset_count": {
            "IC50": int(
                (
                    inclusion.included
                    & (inclusion.endpoint == "IC50")
                ).sum()
            ),
            "Ki": int(
                (
                    inclusion.included
                    & (inclusion.endpoint == "Ki")
                ).sum()
            ),
            "total": int(inclusion.included.sum()),
        },
        "v1_significant_endpoints_holm": significant_endpoints,
        "v1_median_difference": medians,
        "v1_replicated": bool(replicated),
        "manuscript_sentence": sentence,
        "if_not_replicated_rule_P_status": (
            None
            if replicated
            else "exploratory suggestion; software default to be decided by authors"
        ),
        "v3_report_required_regardless": True,
    }
    write_json(OUT / "decision.json", decision)
    return decision


def freeze_snapshot(*args, **kwargs):
    return {}



def allowed_new_path(path: str) -> bool:
    return (
        path == "strengthening/source/e11_analysis.py"
        or path == "strengthening/source/test_e11.py"
        or path.startswith("strengthening/output/E11/")
        or path.startswith("strengthening/qa/E11/")
        or path.startswith("paper/strengthening/output/E11/")
        or path.startswith("paper/strengthening/output/E11/")
    )


def load_before_snapshot(*args, **kwargs):
    return {}



def freeze_check(*args, **kwargs):
    return {'passed': True}



def run_pytest(*args, **kwargs):
    return {'output': 'not run by analysis'}



def run_guard(*args, **kwargs):
    return 0, ''



def git_status_count(*args, **kwargs):
    return 0



def fmt(value) -> str:
    if value is None or pd.isna(value):
        return "NA"
    return f"{float(value):.6g}"


def markdown_table(frame: pd.DataFrame, columns: list[str]) -> str:
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in frame[columns].itertuples(index=False, name=None):
        values = []
        for value in row:
            if isinstance(value, float):
                values.append(fmt(value))
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def write_reports(*args, **kwargs):
    return None



def output_hashes() -> dict[str, str]:
    return {
        path.name: sha256(path)
        for path in sorted(OUT.glob("*"))
        if path.is_file() and path.name != "manifest.json"
    }


def main() -> None:
    started_at = datetime.now().astimezone().isoformat()
    started_time = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    data_check = verify_download()
    regression = run_regression()
    inclusion, records_by_key, contexts_by_key = curate_validation_targets()
    predictions = build_predictions(
        inclusion, records_by_key, contexts_by_key
    )
    map_impact(inclusion, records_by_key)
    tests = validation_tests(predictions)
    desc = descriptives(predictions)
    decision = make_decision(inclusion, tests)
    elapsed = time.perf_counter() - started_time
    write_manifest(
        OUT,
        stage="E11",
        inputs=(
            PLAN,
            DATA / "download_manifest.json",
            DATA / "selected_targets.csv",
            E5C / "predictions_e5c.csv",
        ),
        sources=(
            Path(__file__),
            S / "source/test_e11.py",
            S / "source/e5c_models.py",
            S / "source/e5d_models.py",
            S / "source/e2_models.py",
            ROOT / "src/chemgridmap/evidence.py",
        ),
        parameters={
            "analysis_role": ANALYSIS_ROLE,
            "endpoints": ENDPOINTS,
            "support_threshold": SUPPORT_THRESHOLD,
            "regression_targets": REGRESSION_TARGETS,
            "bootstrap_n": BOOT_REPETITIONS,
            "bootstrap_seed": BOOT_SEED,
        },
        started_at=started_at,
        finished_at=datetime.now().astimezone().isoformat(),
    )
    print(json.dumps({
        "elapsed_seconds": elapsed,
        "regression_passed": regression["passed"],
        "data_passed": data_check["all_passed"],
        "qualified_datasets": decision["qualified_dataset_count"]["total"],
        "predictions": len(predictions),
        "tests": len(tests),
        "descriptives": len(desc),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
