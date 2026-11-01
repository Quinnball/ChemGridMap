"""E12a: post hoc rule-state, sensitivity, and calibration analyses."""
from __future__ import annotations

from datetime import datetime
import gzip
import hashlib
import json
import os
from pathlib import Path
import pickle
import platform
import sys
import time
from tempfile import gettempdir

sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parents[3]
S = ROOT / "paper/strengthening"
OUT = S / "output/E12a"
QA = S / "qa/E12a"
E2 = S / "output/E2"
E2_QA = S / "qa/E2"
E5C = S / "output/E5c"
E5D = S / "output/E5d"
E7 = S / "output/E7"
E9 = S / "output/E9"
E11 = S / "output/E11"
CGM = ROOT
sys.path[:0] = [str(S / "source"), str(CGM / "src"), str(CGM / "paper")]
import chemgridmap  # noqa: E402

CHEMGRIDMAP_FILE = Path(chemgridmap.__file__).resolve()
if CGM not in CHEMGRIDMAP_FILE.parents:
    raise RuntimeError(f"chemgridmap imported outside repository: {CHEMGRIDMAP_FILE}")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from chemgridmap.activity import activity_class  # noqa: E402
from chemgridmap.evidence import (  # noqa: E402
    CLASSES,
    document_values,
    fit_between_unit_sd,
)
from e2_models import (  # noqa: E402
    document_pair_observations,
    flags_from_observations,
    touch_flagged_entries,
)
from e5c_models import weighted_gaussian_prediction  # noqa: E402
from e5d_models import SUPPORT_THRESHOLD, add_rule_states  # noqa: E402
from e7_analysis import (  # noqa: E402
    conditional_agreement,
    connected_components,
    holm_adjust,
    one_sided_wilcoxon,
    two_stage_weights,
)
import e9_analysis as e9  # noqa: E402

ANALYSIS_ROLE = "post hoc"
ENDPOINTS = ("IC50", "Ki")
DESIGN_SPECS = (
    ("entry", "entry", "same_value", None),
    ("document", "document", "same_value", None),
    ("near_duplicate", "entry", "near_duplicate_0.1", None),
    ("time_2015", "time", "same_value", 2015),
    ("time_2010", "time", "same_value", 2010),
)
DESIGN_LOOKUP = {
    name: (split, merge, cutoff)
    for name, split, merge, cutoff in DESIGN_SPECS
}
DEVELOPMENT_TARGETS = {"CHEMBL204", "CHEMBL205", "CHEMBL240"}
TARGET_LABELS = {
    "CHEMBL204": "thrombin",
    "CHEMBL205": "CA II",
    "CHEMBL240": "hERG",
    "CHEMBL5113": "orexin receptor 1",
    "CHEMBL5024": "orexin receptor 2",
    "CHEMBL4792": "orexin receptor 2",
    "CHEMBL2971": "JAK2",
    "CHEMBL206": "estrogen receptor",
    "CHEMBL325": "HDAC1",
    "CHEMBL1163125": "BRD4",
}
STATE_MASKS = (
    "G_supported",
    "P_supported",
    "demoted",
    "touched_all",
    "touched_uncertain",
    "uncertain_low_support",
    "untouched",
    "all",
)
CALIBRATION_POPULATIONS = (
    "all",
    "untouched",
    "touched",
    "G_supported",
    "P_supported",
)
BOOT_REPETITIONS = 2000
E9_BOOT_SEED = 20261105
NEW_BOOT_SEED = 20261112


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_safe(value):
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def write_json(path: Path, value) -> None:
    path.write_text(
        json.dumps(json_safe(value), indent=2, ensure_ascii=False, allow_nan=False)
        + "\n"
    )


def with_n_alias(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    if "n_predictions" in result and "n" not in result:
        result["n"] = result["n_predictions"]
    return result


def write_csv(frame: pd.DataFrame, name: str) -> None:
    with_n_alias(frame).to_csv(OUT / name, index=False)


def dev_targets_in(frame: pd.DataFrame) -> str:
    return ";".join(
        sorted(set(frame.target.astype(str)) & DEVELOPMENT_TARGETS)
    )


def design_frame(predictions: pd.DataFrame, design: str) -> pd.DataFrame:
    split_design, merge_basis, cutoff = DESIGN_LOOKUP[design]
    selected = predictions.loc[
        (predictions.split_design == split_design)
        & (predictions.merge_basis == merge_basis)
        & (predictions.method == "G")
    ].copy()
    if cutoff is not None:
        selected = selected.loc[selected.cutoff_year == cutoff]
    return add_rule_states(selected.reset_index(drop=True))


def load_designs() -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    usecols = [
        "split_design",
        "merge_basis",
        "target",
        "endpoint",
        "fold",
        "cutoff_year",
        "entry",
        "held_out_unit",
        "m",
        "m_layer",
        "displayed",
        "method",
        "p_displayed",
        "support",
        "agree",
        "flag_touched",
    ]
    predictions = pd.read_csv(
        E5C / "predictions_e5c.csv",
        usecols=usecols,
        low_memory=False,
        float_precision="round_trip",
    )
    designs = {
        name: design_frame(predictions, name) for name in DESIGN_LOOKUP
    }
    return predictions, designs


def state_mask(frame: pd.DataFrame, state: str) -> np.ndarray:
    support = frame.support.to_numpy(float) >= SUPPORT_THRESHOLD
    touched = frame.flag_touched.to_numpy(bool)
    p_status = frame.p_status.to_numpy()
    masks = {
        "G_supported": support,
        "P_supported": p_status == "supported",
        "demoted": support & touched,
        "touched_all": touched,
        "touched_uncertain": p_status == "uncertain_touched",
        "uncertain_low_support": p_status == "uncertain_low_support",
        "untouched": ~touched,
        "all": np.ones(len(frame), dtype=bool),
    }
    return masks[state]


def bootstrap_cell(
    frame: pd.DataFrame,
    weights: np.ndarray,
    mask: np.ndarray,
) -> dict[str, float | int]:
    """Agreement/bias intervals matching E9, plus mean-probability interval."""
    mask = np.asarray(mask, dtype=bool)
    n = int(mask.sum())
    if not n:
        return {
            "mean_probability_ci_low": np.nan,
            "mean_probability_ci_high": np.nan,
            "agreement_ci_low": np.nan,
            "agreement_ci_high": np.nan,
            "bias_ci_low": np.nan,
            "bias_ci_high": np.nan,
            "valid_repetitions": 0,
        }
    base = e9.bootstrap_cell(frame, weights, mask)
    denominator = weights[:, mask].sum(axis=1)
    p_values = frame.p_displayed.to_numpy(float)[mask]
    probability_draws = np.divide(
        (weights[:, mask] * p_values).sum(axis=1),
        denominator,
        out=np.full(len(weights), np.nan),
        where=denominator > 0,
    )
    finite = np.isfinite(probability_draws)
    if finite.any():
        probability_bounds = np.percentile(
            probability_draws[finite], [2.5, 97.5]
        )
    else:
        probability_bounds = (np.nan, np.nan)
    return {
        "mean_probability_ci_low": float(probability_bounds[0]),
        "mean_probability_ci_high": float(probability_bounds[1]),
        **base,
    }


def point_summary(frame: pd.DataFrame, mask: np.ndarray) -> dict[str, float | int]:
    selected = frame.loc[np.asarray(mask, dtype=bool)]
    n = len(selected)
    mean_probability = float(selected.p_displayed.mean()) if n else np.nan
    agreement = float(selected.agree.mean()) if n else np.nan
    return {
        "n_predictions": n,
        "n_datasets": int(selected.target.nunique()) if n else 0,
        "mean_probability": mean_probability,
        "agreement_rate": agreement,
        "bias": mean_probability - agreement if n else np.nan,
    }


def rule_states(designs: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, dict]:
    rows = []
    metadata = {}
    for design, frame in designs.items():
        metadata[design] = {}
        for endpoint in ENDPOINTS:
            endpoint_frame = frame.loc[frame.endpoint == endpoint].reset_index(drop=True)
            weights = two_stage_weights(
                endpoint_frame,
                seed=E9_BOOT_SEED,
                repetitions=BOOT_REPETITIONS,
            )
            metadata[design][endpoint] = {
                "n_predictions": len(endpoint_frame),
                "n_datasets": int(endpoint_frame.target.nunique()),
                "weights": weights,
            }
            for state in STATE_MASKS:
                mask = state_mask(endpoint_frame, state)
                point = point_summary(endpoint_frame, mask)
                boot = bootstrap_cell(endpoint_frame, weights, mask)
                selected = endpoint_frame.loc[mask]
                rows.append({
                    "analysis_role": ANALYSIS_ROLE,
                    "design": design,
                    "endpoint": endpoint,
                    "state": state,
                    **point,
                    **boot,
                    "is_development_target": bool(
                        set(selected.target.astype(str)) & DEVELOPMENT_TARGETS
                    ),
                    "development_targets": dev_targets_in(selected),
                    "bootstrap_method": "two_stage",
                    "bootstrap_repetitions": BOOT_REPETITIONS,
                    "bootstrap_seed": E9_BOOT_SEED,
                    "two_stage_function": "e7_analysis.two_stage_weights",
                })
    return pd.DataFrame(rows), metadata


def two_sided_wilcoxon(values) -> tuple[float, float, int]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values) or np.allclose(values, 0):
        return np.nan, 1.0, len(values)
    result = stats.wilcoxon(values, alternative="two-sided", zero_method="wilcox")
    return float(result.statistic), float(result.pvalue), len(values)


def wilcoxon_record(values) -> dict:
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    statistic, two_p, n = two_sided_wilcoxon(finite)
    _, one_p, _ = one_sided_wilcoxon(finite, alternative="greater")
    return {
        "n_units": n,
        "median_difference": float(np.median(finite)) if n else np.nan,
        "positive_units": int((finite > 0).sum()),
        "zero_units": int((finite == 0).sum()),
        "negative_units": int((finite < 0).sum()),
        "wilcoxon_statistic": statistic,
        "two_sided_p": two_p,
        "one_sided_p": one_p,
    }


def dataset_differences(frame: pd.DataFrame, test: str) -> pd.DataFrame:
    rows = []
    for target, group in frame.groupby("target", sort=True):
        g_supported = group.g_status == "supported"
        p_supported = group.p_status == "supported"
        demoted = g_supported & group.flag_touched.astype(bool)
        touched_uncertain = group.p_status == "uncertain_touched"
        low_support = group.p_status == "uncertain_low_support"
        if test == "rise":
            eligible = bool(demoted.any()) and bool(g_supported.any()) and bool(
                p_supported.any()
            )
            supported_n = int(g_supported.sum())
            comparison_n = int(p_supported.sum())
        elif test == "retained_minus_demoted":
            eligible = int(p_supported.sum()) >= 5 and int(demoted.sum()) >= 5
            supported_n = int(p_supported.sum())
            comparison_n = int(demoted.sum())
        elif test == "touched_minus_low_support":
            eligible = (
                int(touched_uncertain.sum()) >= 5
                and int(low_support.sum()) >= 5
            )
            supported_n = int(touched_uncertain.sum())
            comparison_n = int(low_support.sum())
        else:
            raise ValueError(test)
        difference = np.nan
        if eligible:
            if test == "rise":
                difference = float(
                    group.loc[p_supported, "agree"].mean()
                    - group.loc[g_supported, "agree"].mean()
                )
            elif test == "retained_minus_demoted":
                difference = float(
                    group.loc[p_supported, "agree"].mean()
                    - group.loc[demoted, "agree"].mean()
                )
            else:
                difference = float(
                    group.loc[touched_uncertain, "agree"].mean()
                    - group.loc[low_support, "agree"].mean()
                )
        rows.append({
            "target": target,
            "eligible": eligible,
            "difference": difference,
            "left_n": supported_n,
            "right_n": comparison_n,
        })
    return pd.DataFrame(rows)


def component_values(
    dataset_values: pd.DataFrame,
    membership: pd.DataFrame,
) -> pd.DataFrame:
    eligible = dataset_values.loc[dataset_values.eligible]
    mapping = membership[["component", "target"]].drop_duplicates()
    joined = eligible.merge(mapping, on="target", how="inner", validate="many_to_one")
    rows = []
    for component, group in joined.groupby("component", sort=True):
        values = group.difference.to_numpy(float)
        values = values[np.isfinite(values)]
        if not len(values):
            continue
        rows.append({
            "component": component,
            "difference": float(values.mean()),
            "n_eligible_datasets": int(len(values)),
        })
    return pd.DataFrame(rows)


def rule_tests(
    designs: dict[str, pd.DataFrame],
    cluster_members: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for design, frame in designs.items():
        for endpoint in ENDPOINTS:
            endpoint_frame = frame.loc[frame.endpoint == endpoint]
            endpoint_members = cluster_members.loc[
                cluster_members.endpoint == endpoint
            ]
            for test in (
                "rise",
                "retained_minus_demoted",
                "touched_minus_low_support",
            ):
                differences = dataset_differences(endpoint_frame, test)
                for unit, values in (
                    ("dataset", differences),
                    ("shared_document_component", component_values(
                        differences, endpoint_members
                    )),
                ):
                    result = wilcoxon_record(values.difference)
                    rows.append({
                        "analysis_role": ANALYSIS_ROLE,
                        "design": design,
                        "endpoint": endpoint,
                        "test": test,
                        "unit": unit,
                        **result,
                        "is_development_target": bool(
                            set(endpoint_frame.target.astype(str))
                            & DEVELOPMENT_TARGETS
                        ),
                        "development_targets": dev_targets_in(endpoint_frame),
                        "alternative": "greater_than_zero",
                        "zero_method": "wilcox",
                    })
    result = pd.DataFrame(rows)
    result["holm_two_sided_p"] = np.nan
    for endpoint, indices in result.groupby("endpoint").groups.items():
        result.loc[indices, "holm_two_sided_p"] = holm_adjust(
            result.loc[indices, "two_sided_p"]
        )
    return result


def lodo_point_quantities(frame: pd.DataFrame) -> list[tuple[str, np.ndarray]]:
    inactive = frame.displayed.to_numpy() == CLASSES.index("inactive")
    touched = frame.flag_touched.to_numpy(bool)
    support = frame.support.to_numpy(float) >= SUPPORT_THRESHOLD
    g_supported = frame.g_status.to_numpy() == "supported"
    p_supported = frame.p_status.to_numpy() == "supported"
    demoted = g_supported & touched
    return [
        ("inactive_bias_all", inactive),
        ("inactive_bias_untouched", inactive & ~touched),
        ("inactive_bias_touched", inactive & touched),
        ("G_supported_agreement", g_supported),
        ("P_supported_agreement_and_probability", p_supported),
        ("demoted_agreement", demoted),
        ("touched_uncertain_agreement", frame.p_status.to_numpy() == "uncertain_touched"),
        ("uncertain_low_support_agreement", frame.p_status.to_numpy() == "uncertain_low_support"),
    ]


def lodo_rows(
    frame: pd.DataFrame,
    design: str,
    endpoint: str,
    scenario: str,
    omitted: set[str],
    named: bool,
) -> list[dict]:
    rows = []
    weights = None
    if named:
        weights = two_stage_weights(
            frame, seed=NEW_BOOT_SEED, repetitions=BOOT_REPETITIONS
        )
    for quantity, mask in lodo_point_quantities(frame):
        point = point_summary(frame, mask)
        boot = (
            bootstrap_cell(frame, weights, mask)
            if named
            else {
                "mean_probability_ci_low": np.nan,
                "mean_probability_ci_high": np.nan,
                "agreement_ci_low": np.nan,
                "agreement_ci_high": np.nan,
                "bias_ci_low": np.nan,
                "bias_ci_high": np.nan,
                "valid_repetitions": 0,
            }
        )
        rows.append({
            "analysis_role": ANALYSIS_ROLE,
            "design": design,
            "endpoint": endpoint,
            "scenario": scenario,
            "omitted_targets": ";".join(sorted(omitted)),
            "quantity": quantity,
            **point,
            **boot,
            "test_n_units": np.nan,
            "two_sided_p": np.nan,
            "is_development_target": bool(
                set(frame.target.astype(str)) & DEVELOPMENT_TARGETS
            ),
            "development_targets": dev_targets_in(frame),
            "bootstrap_method": "two_stage" if named else "not_applicable",
            "bootstrap_repetitions": BOOT_REPETITIONS if named else 0,
            "bootstrap_seed": NEW_BOOT_SEED if named else np.nan,
            "two_stage_function": (
                "e7_analysis.two_stage_weights" if named else ""
            ),
        })
    rise = dataset_differences(frame, "rise")
    result = wilcoxon_record(rise.difference)
    rows.append({
        "analysis_role": ANALYSIS_ROLE,
        "design": design,
        "endpoint": endpoint,
        "scenario": scenario,
        "omitted_targets": ";".join(sorted(omitted)),
        "quantity": "rise_test",
        "n_predictions": np.nan,
        "n_datasets": int(len(frame.target.unique())),
        "mean_probability": np.nan,
        "agreement_rate": np.nan,
        "bias": np.nan,
        "mean_probability_ci_low": np.nan,
        "mean_probability_ci_high": np.nan,
        "agreement_ci_low": np.nan,
        "agreement_ci_high": np.nan,
        "bias_ci_low": np.nan,
        "bias_ci_high": np.nan,
        "valid_repetitions": 0,
        "test_n_units": result["n_units"],
        "median_difference": result["median_difference"],
        "positive_units": result["positive_units"],
        "zero_units": result["zero_units"],
        "negative_units": result["negative_units"],
        "two_sided_p": result["two_sided_p"],
        "one_sided_p": result["one_sided_p"],
        "is_development_target": bool(
            set(frame.target.astype(str)) & DEVELOPMENT_TARGETS
        ),
        "development_targets": dev_targets_in(frame),
        "bootstrap_method": "not_applicable",
        "bootstrap_repetitions": 0,
        "bootstrap_seed": np.nan,
        "two_stage_function": "",
    })
    return rows


def build_lodo(
    designs: dict[str, pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    all_rows = []
    named_rows = []
    top_ic50 = {}
    for design in ("entry", "document"):
        for endpoint in ENDPOINTS:
            full_frame = designs[design].loc[
                designs[design].endpoint == endpoint
            ].reset_index(drop=True)
            targets = sorted(full_frame.target.unique())
            all_rows.extend(
                lodo_rows(
                    full_frame,
                    design,
                    endpoint,
                    "none",
                    set(),
                    named=False,
                )
            )
            for target in targets:
                selected = full_frame.loc[
                    full_frame.target != target
                ].reset_index(drop=True)
                all_rows.extend(
                    lodo_rows(
                        selected,
                        design,
                        endpoint,
                        f"omit_{target}",
                        {target},
                        named=False,
                    )
                )

            if endpoint == "Ki":
                scenarios = (
                    ("none", set()),
                    ("without_CHEMBL5113", {"CHEMBL5113"}),
                    ("without_CHEMBL5024", {"CHEMBL5024"}),
                    ("without_CHEMBL5113_and_CHEMBL5024", {
                        "CHEMBL5113",
                        "CHEMBL5024",
                    }),
                )
            else:
                demoted = (
                    (full_frame.g_status == "supported")
                    & full_frame.flag_touched.astype(bool)
                )
                counts = demoted.groupby(full_frame.target).sum().sort_values(
                    ascending=False, kind="stable"
                )
                top_targets = list(counts.head(4).index)
                top_ic50[design] = [
                    {
                        "target": target,
                        "target_label": TARGET_LABELS.get(target, ""),
                        "demoted_predictions": int(counts.loc[target]),
                    }
                    for target in top_targets
                ]
                scenarios = [("none", set())]
                scenarios.extend(
                    (f"omit_{target}", {target}) for target in top_targets
                )
                scenarios.append(("omit_top4_demoted", set(top_targets)))
            for scenario, omitted in scenarios:
                selected = full_frame.loc[
                    ~full_frame.target.isin(omitted)
                ].reset_index(drop=True)
                named_rows.extend(
                    lodo_rows(
                        selected,
                        design,
                        endpoint,
                        scenario,
                        omitted,
                        named=True,
                    )
                )
    return pd.DataFrame(all_rows), pd.DataFrame(named_rows), top_ic50


def weighted_dataset_biases(selected: pd.DataFrame) -> np.ndarray:
    return (
        selected.groupby("target", sort=True)
        .apply(
            lambda group: float(
                (group.p_displayed - group.agree).mean()
            ),
            include_groups=False,
        )
        .to_numpy(float)
    )


def calibration_all_designs(
    designs: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    rows = []
    for design, frame in designs.items():
        for endpoint in ENDPOINTS:
            endpoint_frame = frame.loc[frame.endpoint == endpoint].reset_index(drop=True)
            weights = two_stage_weights(
                endpoint_frame,
                seed=E9_BOOT_SEED,
                repetitions=BOOT_REPETITIONS,
            )
            masks = {
                "all": np.ones(len(endpoint_frame), dtype=bool),
                "untouched": ~endpoint_frame.flag_touched.to_numpy(bool),
                "touched": endpoint_frame.flag_touched.to_numpy(bool),
                "G_supported": (
                    endpoint_frame.g_status.to_numpy() == "supported"
                ),
                "P_supported": (
                    endpoint_frame.p_status.to_numpy() == "supported"
                ),
            }
            for population in CALIBRATION_POPULATIONS:
                for display_name in (*CLASSES, "all"):
                    display_mask = (
                        np.ones(len(endpoint_frame), dtype=bool)
                        if display_name == "all"
                        else endpoint_frame.displayed.to_numpy()
                        == CLASSES.index(display_name)
                    )
                    mask = masks[population] & display_mask
                    selected = endpoint_frame.loc[mask]
                    point = point_summary(endpoint_frame, mask)
                    boot = bootstrap_cell(endpoint_frame, weights, mask)
                    stats10 = e9.dataset_bias_stats(selected, 10)
                    stats1 = e9.dataset_bias_stats(selected, 1)
                    random_effects = e9.random_effect_details(selected, 20)
                    dataset_biases = weighted_dataset_biases(selected)
                    rows.append({
                        "analysis_role": ANALYSIS_ROLE,
                        "design": design,
                        "endpoint": endpoint,
                        "population": population,
                        "display_class": display_name,
                        **point,
                        **boot,
                        "dataset_median_min10_n": stats10["dataset_median_n"],
                        "dataset_bias_median_min10": stats10["bias_median"],
                        "dataset_bias_q1_min10": stats10["bias_q1"],
                        "dataset_bias_q3_min10": stats10["bias_q3"],
                        "positive_datasets_min10": stats10["positive_datasets"],
                        "zero_datasets_min10": stats10["zero_datasets"],
                        "negative_datasets_min10": stats10["negative_datasets"],
                        "dataset_bias_two_sided_p_min10": stats10[
                            "wilcoxon_two_sided_p"
                        ],
                        "dataset_median_all_n": stats1["dataset_median_n"],
                        "dataset_bias_median_all": stats1["bias_median"],
                        "dataset_bias_q1_all": stats1["bias_q1"],
                        "dataset_bias_q3_all": stats1["bias_q3"],
                        "positive_datasets_all": stats1["positive_datasets"],
                        "zero_datasets_all": stats1["zero_datasets"],
                        "negative_datasets_all": stats1["negative_datasets"],
                        "dataset_bias_two_sided_p_all": stats1[
                            "wilcoxon_two_sided_p"
                        ],
                        "dataset_bias_simple_mean": (
                            float(dataset_biases.mean())
                            if len(dataset_biases)
                            else np.nan
                        ),
                        "random_effects_n_datasets": random_effects["n_datasets"],
                        "random_effects_bias": random_effects["random_estimate"],
                        "random_effects_ci_low": random_effects["ci_low"],
                        "random_effects_ci_high": random_effects["ci_high"],
                        "random_effects_tau2": random_effects["tau2"],
                        "random_effects_i_squared": random_effects[
                            "i_squared"
                        ],
                        "random_effects_pi_low": random_effects["pi_low"],
                        "random_effects_pi_high": random_effects["pi_high"],
                        "is_development_target": bool(
                            set(selected.target.astype(str))
                            & DEVELOPMENT_TARGETS
                        ),
                        "development_targets": dev_targets_in(selected),
                        "minimum_dataset_predictions_median": 10,
                        "minimum_dataset_predictions_random_effects": 20,
                        "bootstrap_method": "two_stage",
                        "bootstrap_repetitions": BOOT_REPETITIONS,
                        "bootstrap_seed": E9_BOOT_SEED,
                        "two_stage_function": "e7_analysis.two_stage_weights",
                    })
    return pd.DataFrame(rows)


def time_split_rule(
    design_results: pd.DataFrame,
    designs: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    rows = design_results.loc[
        design_results.design.isin(("time_2015", "time_2010"))
    ].copy()
    extra = {
        "P_supported_inactive_n_predictions": [],
        "P_supported_inactive_n_datasets": [],
        "P_supported_inactive_mean_probability": [],
        "P_supported_inactive_agreement_rate": [],
        "P_supported_inactive_bias": [],
        "P_supported_inactive_agreement_ci_low": [],
        "P_supported_inactive_agreement_ci_high": [],
        "P_supported_inactive_bias_ci_low": [],
        "P_supported_inactive_bias_ci_high": [],
        "P_supported_inactive_without_CHEMBL5113_n_predictions": [],
        "P_supported_inactive_without_CHEMBL5113_mean_probability": [],
        "P_supported_inactive_without_CHEMBL5113_agreement_rate": [],
        "P_supported_inactive_without_CHEMBL5113_bias": [],
        "P_supported_inactive_without_CHEMBL5113_agreement_ci_low": [],
        "P_supported_inactive_without_CHEMBL5113_agreement_ci_high": [],
        "P_supported_inactive_without_CHEMBL5113_bias_ci_low": [],
        "P_supported_inactive_without_CHEMBL5113_bias_ci_high": [],
    }
    for row in rows.itertuples(index=False):
        frame = designs[row.design]
        endpoint_frame = frame.loc[frame.endpoint == row.endpoint].reset_index(drop=True)
        mask = (
            endpoint_frame.p_status.to_numpy() == "supported"
        ) & (
            endpoint_frame.displayed.to_numpy() == CLASSES.index("inactive")
        )
        weights = two_stage_weights(
            endpoint_frame,
            seed=E9_BOOT_SEED,
            repetitions=BOOT_REPETITIONS,
        )
        point = point_summary(endpoint_frame, mask)
        boot = bootstrap_cell(endpoint_frame, weights, mask)
        for key, value in (
            ("n_predictions", point["n_predictions"]),
            ("n_datasets", point["n_datasets"]),
            ("mean_probability", point["mean_probability"]),
            ("agreement_rate", point["agreement_rate"]),
            ("bias", point["bias"]),
            ("agreement_ci_low", boot["agreement_ci_low"]),
            ("agreement_ci_high", boot["agreement_ci_high"]),
            ("bias_ci_low", boot["bias_ci_low"]),
            ("bias_ci_high", boot["bias_ci_high"]),
        ):
            extra[f"P_supported_inactive_{key}"].append(value)

        reduced = endpoint_frame.loc[
            endpoint_frame.target != "CHEMBL5113"
        ].reset_index(drop=True)
        reduced_mask = (
            reduced.p_status.to_numpy() == "supported"
        ) & (
            reduced.displayed.to_numpy() == CLASSES.index("inactive")
        )
        reduced_weights = two_stage_weights(
            reduced,
            seed=NEW_BOOT_SEED,
            repetitions=BOOT_REPETITIONS,
        )
        reduced_point = point_summary(reduced, reduced_mask)
        reduced_boot = bootstrap_cell(
            reduced, reduced_weights, reduced_mask
        )
        for key, value in (
            ("n_predictions", reduced_point["n_predictions"]),
            ("mean_probability", reduced_point["mean_probability"]),
            ("agreement_rate", reduced_point["agreement_rate"]),
            ("bias", reduced_point["bias"]),
            ("agreement_ci_low", reduced_boot["agreement_ci_low"]),
            ("agreement_ci_high", reduced_boot["agreement_ci_high"]),
            ("bias_ci_low", reduced_boot["bias_ci_low"]),
            ("bias_ci_high", reduced_boot["bias_ci_high"]),
        ):
            extra[
                f"P_supported_inactive_without_CHEMBL5113_{key}"
            ].append(value)
    for column, values in extra.items():
        rows[column] = values
    return rows.reset_index(drop=True)


def fig4d_three_groups(designs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    time_frame = designs["time_2015"]
    rows = []
    for endpoint in ENDPOINTS:
        endpoint_frame = time_frame.loc[
            time_frame.endpoint == endpoint
        ].copy()
        endpoint_frame = endpoint_frame.loc[
            ~(
                endpoint_frame.target.isin(DEVELOPMENT_TARGETS)
                & endpoint_frame.flag_touched.notna()
                & (endpoint_frame.endpoint == "IC50")
            )
        ].reset_index(drop=True)
        weights = two_stage_weights(
            endpoint_frame,
            seed=NEW_BOOT_SEED,
            repetitions=BOOT_REPETITIONS,
        )
        for m_name, m_mask in (
            ("m_eq_1", endpoint_frame.m_layer.to_numpy() == 1),
            ("m_ge_2", endpoint_frame.m_layer.to_numpy() >= 2),
        ):
            for state, state_values in (
                ("P_supported", endpoint_frame.p_status.to_numpy() == "supported"),
                (
                    "uncertain_low_support",
                    endpoint_frame.p_status.to_numpy() == "uncertain_low_support",
                ),
                (
                    "touched_all",
                    endpoint_frame.flag_touched.to_numpy(bool),
                ),
            ):
                mask = m_mask & state_values
                selected = endpoint_frame.loc[mask]
                point = point_summary(endpoint_frame, mask)
                boot = bootstrap_cell(endpoint_frame, weights, mask)
                rows.append({
                    "analysis_role": ANALYSIS_ROLE,
                    "design": "time_2015",
                    "endpoint": endpoint,
                    "earlier_units_group": m_name,
                    "state": state,
                    **point,
                    **boot,
                    "is_development_target": bool(
                        set(selected.target.astype(str)) & DEVELOPMENT_TARGETS
                    ),
                    "development_targets": dev_targets_in(selected),
                    "population_note": (
                        "52-dataset time split; development-target IC50 "
                        "predictions excluded"
                    ),
                    "bootstrap_method": "two_stage",
                    "bootstrap_repetitions": BOOT_REPETITIONS,
                    "bootstrap_seed": NEW_BOOT_SEED,
                    "two_stage_function": "e7_analysis.two_stage_weights",
                })
    return pd.DataFrame(rows)


def shared_components_for_predictions(
    predictions: pd.DataFrame,
    endpoint: str,
) -> pd.DataFrame:
    docs = {}
    selected_predictions = predictions.loc[predictions.endpoint == endpoint]
    for (target, item_endpoint), group in selected_predictions.groupby(
        ["target", "endpoint"], sort=True
    ):
        documents = set()
        for value in group.held_out_unit.astype(str):
            documents.update(part for part in value.split(";") if part)
        docs[(target, item_endpoint)] = documents
    nodes = sorted(key for key in docs if key[1] == endpoint)
    edges = []
    for number, pair_a in enumerate(nodes):
        for pair_b in nodes[number + 1:]:
            docs_a, docs_b = docs[pair_a], docs[pair_b]
            denominator = min(len(docs_a), len(docs_b))
            overlap = len(docs_a & docs_b) / denominator if denominator else 0.0
            if overlap >= 0.2:
                edges.append((pair_a, pair_b))
    rows = []
    for number, component in enumerate(
        connected_components(nodes, edges), start=1
    ):
        component_id = f"{endpoint}_C{number:02d}"
        for target in sorted(target for target, _ in component):
            rows.append({
                "endpoint": endpoint,
                "component": component_id,
                "target": target,
                "component_size": len(component),
                "threshold": 0.2,
            })
    return pd.DataFrame(rows, columns=[
        "endpoint", "component", "target", "component_size", "threshold"
    ])


def e11_sensitivity() -> tuple[pd.DataFrame, pd.DataFrame]:
    usecols = [
        "target",
        "endpoint",
        "entry",
        "held_out_unit",
        "displayed",
        "p_displayed",
        "support",
        "agree",
        "flag_touched",
    ]
    predictions = pd.read_csv(
        E11 / "predictions_e11.csv",
        usecols=usecols,
        low_memory=False,
        float_precision="round_trip",
    )
    predictions = add_rule_states(predictions.reset_index(drop=True))
    rows = []
    component_frames = []
    for endpoint in ENDPOINTS:
        endpoint_frame = predictions.loc[
            predictions.endpoint == endpoint
        ].reset_index(drop=True)
        v1_values = dataset_differences(endpoint_frame, "rise")
        for item in v1_values.itertuples(index=False):
            rows.append({
                "analysis_role": ANALYSIS_ROLE,
                "section": "dataset_v1",
                "endpoint": endpoint,
                "scenario": "all",
                "unit": "dataset",
                "target": item.target,
                "n_predictions": np.nan,
                "n_datasets": 1,
                "value": item.difference,
                "two_sided_p": np.nan,
                "is_development_target": item.target in DEVELOPMENT_TARGETS,
                "development_targets": (
                    item.target if item.target in DEVELOPMENT_TARGETS else ""
                ),
            })
        test_specs = [("all", set())]
        if endpoint == "Ki":
            test_specs.append(("without_CHEMBL4792", {"CHEMBL4792"}))
        for scenario, omitted in test_specs:
            selected = endpoint_frame.loc[
                ~endpoint_frame.target.isin(omitted)
            ].reset_index(drop=True)
            differences = dataset_differences(selected, "rise")
            result = wilcoxon_record(differences.difference)
            rows.append({
                "analysis_role": ANALYSIS_ROLE,
                "section": "v1_test",
                "endpoint": endpoint,
                "scenario": scenario,
                "unit": "dataset",
                "n_datasets": result["n_units"],
                "value": result["median_difference"],
                "positive_units": result["positive_units"],
                "zero_units": result["zero_units"],
                "negative_units": result["negative_units"],
                "two_sided_p": result["two_sided_p"],
                "one_sided_p": result["one_sided_p"],
                "is_development_target": bool(
                    set(selected.target.astype(str)) & DEVELOPMENT_TARGETS
                ),
                "development_targets": dev_targets_in(selected),
            })
            for target in sorted(selected.target.unique()):
                reduced = selected.loc[
                    selected.target != target
                ].reset_index(drop=True)
                lodo_result = wilcoxon_record(
                    dataset_differences(reduced, "rise").difference
                )
                rows.append({
                    "analysis_role": ANALYSIS_ROLE,
                    "section": "v1_lodo",
                    "endpoint": endpoint,
                    "scenario": scenario,
                    "unit": "dataset",
                    "omitted_targets": target,
                    "n_datasets": lodo_result["n_units"],
                    "value": lodo_result["median_difference"],
                    "positive_units": lodo_result["positive_units"],
                    "zero_units": lodo_result["zero_units"],
                    "negative_units": lodo_result["negative_units"],
                    "two_sided_p": lodo_result["two_sided_p"],
                    "one_sided_p": lodo_result["one_sided_p"],
                    "is_development_target": bool(
                        set(reduced.target.astype(str)) & DEVELOPMENT_TARGETS
                    ),
                    "development_targets": dev_targets_in(reduced),
                })

            weights = two_stage_weights(
                selected,
                seed=NEW_BOOT_SEED,
                repetitions=BOOT_REPETITIONS,
            )
            state_masks = {
                "P_supported": selected.p_status.to_numpy() == "supported",
                "demoted": (
                    (selected.g_status.to_numpy() == "supported")
                    & selected.flag_touched.to_numpy(bool)
                ),
                "touched_uncertain": (
                    selected.p_status.to_numpy() == "uncertain_touched"
                ),
                "uncertain_low_support": (
                    selected.p_status.to_numpy() == "uncertain_low_support"
                ),
            }
            for state, mask in state_masks.items():
                point = point_summary(selected, mask)
                boot = bootstrap_cell(selected, weights, mask)
                rows.append({
                    "analysis_role": ANALYSIS_ROLE,
                    "section": "state_summary",
                    "endpoint": endpoint,
                    "scenario": scenario,
                    "unit": "prediction",
                    "state": state,
                    **point,
                    **boot,
                    "is_development_target": bool(
                        set(selected.target.astype(str)) & DEVELOPMENT_TARGETS
                    ),
                    "development_targets": dev_targets_in(selected),
                    "bootstrap_method": "two_stage",
                    "bootstrap_repetitions": BOOT_REPETITIONS,
                    "bootstrap_seed": NEW_BOOT_SEED,
                    "two_stage_function": "e7_analysis.two_stage_weights",
                })

        membership = shared_components_for_predictions(predictions, endpoint)
        component_frames.append(membership)
        shared_values = component_values(v1_values, membership)
        shared_result = wilcoxon_record(shared_values.difference)
        rows.append({
            "analysis_role": ANALYSIS_ROLE,
            "section": "shared_component_v1",
            "endpoint": endpoint,
            "scenario": "all",
            "unit": "shared_document_component",
            "n_datasets": 1,
            "value": shared_result["median_difference"],
            "test_n_units": shared_result["n_units"],
            "positive_units": shared_result["positive_units"],
            "zero_units": shared_result["zero_units"],
            "negative_units": shared_result["negative_units"],
            "two_sided_p": shared_result["two_sided_p"],
            "one_sided_p": shared_result["one_sided_p"],
            "is_development_target": bool(
                set(endpoint_frame.target.astype(str)) & DEVELOPMENT_TARGETS
            ),
            "development_targets": dev_targets_in(endpoint_frame),
        })
    return pd.DataFrame(rows), pd.concat(component_frames, ignore_index=True)


def full_map_touched_units() -> tuple[pd.DataFrame, dict]:
    records = pickle.loads(
        gzip.decompress((E2_QA / "curated_records.pkl.gz").read_bytes())
    )
    inclusion = pd.read_csv(E2 / "target_inclusion.csv")
    rows = []
    for item in inclusion.loc[inclusion.included].itertuples(index=False):
        target = str(item.target)
        endpoint = str(item.endpoint)
        target_records = records[(target, endpoint)]
        doc = document_values(target_records, merge_identical=True)
        observations = document_pair_observations(target_records)
        flags = flags_from_observations(observations)
        touched_entries = touch_flagged_entries(doc, flags)
        fit = fit_between_unit_sd(doc)
        sigma = float(fit["tau"])
        record_class = (
            target_records.assign(
                pchembl_value=pd.to_numeric(
                    target_records.pchembl_value, errors="coerce"
                )
            )
            .groupby("molecule_identity_key", sort=True)
            .pchembl_value.median()
            .map(lambda value: activity_class(float(value), 6.0, 7.0))
        )
        evidence_units = doc.groupby("entry", sort=True).size()
        g_supported_entries = set()
        for entry, group in doc.groupby("entry", sort=True):
            values = group.document_value.to_numpy(float)
            support = weighted_gaussian_prediction(
                values,
                np.full(len(values), sigma**2),
                sigma,
            )[0]
            display_class = str(record_class.loc[entry])
            if float(support[CLASSES.index(display_class)]) >= SUPPORT_THRESHOLD:
                g_supported_entries.add(entry)
        demoted_entries = g_supported_entries & touched_entries
        entry_count = len(evidence_units)

        def split_count(entries: set[str]) -> tuple[int, int]:
            n_single = sum(int(evidence_units.loc[entry] == 1) for entry in entries)
            return n_single, len(entries) - n_single

        touched_single, touched_multi = split_count(touched_entries)
        demoted_single, demoted_multi = split_count(demoted_entries)
        rows.append({
            "analysis_role": ANALYSIS_ROLE,
            "target": target,
            "endpoint": endpoint,
            "n_entries": entry_count,
            "single_unit_entries_n": int((evidence_units == 1).sum()),
            "single_unit_entries_rate": (
                float((evidence_units == 1).mean()) if entry_count else np.nan
            ),
            "touched_entries_n": len(touched_entries),
            "touched_single_unit_n": touched_single,
            "touched_multiple_unit_n": touched_multi,
            "demoted_entries_n": len(demoted_entries),
            "demoted_single_unit_n": demoted_single,
            "demoted_multiple_unit_n": demoted_multi,
            "g_supported_entries_n": len(g_supported_entries),
            "is_development_target": target in DEVELOPMENT_TARGETS,
            "development_targets": (
                target if target in DEVELOPMENT_TARGETS else ""
            ),
        })
    by_dataset = pd.DataFrame(rows)
    total = {
        "analysis_role": ANALYSIS_ROLE,
        "target": "ALL",
        "endpoint": "ALL",
        "n_entries": int(by_dataset.n_entries.sum()),
        "single_unit_entries_n": int(by_dataset.single_unit_entries_n.sum()),
        "single_unit_entries_rate": float(
            by_dataset.single_unit_entries_n.sum() / by_dataset.n_entries.sum()
        ),
        "touched_entries_n": int(by_dataset.touched_entries_n.sum()),
        "touched_single_unit_n": int(by_dataset.touched_single_unit_n.sum()),
        "touched_multiple_unit_n": int(by_dataset.touched_multiple_unit_n.sum()),
        "demoted_entries_n": int(by_dataset.demoted_entries_n.sum()),
        "demoted_single_unit_n": int(by_dataset.demoted_single_unit_n.sum()),
        "demoted_multiple_unit_n": int(by_dataset.demoted_multiple_unit_n.sum()),
        "g_supported_entries_n": int(by_dataset.g_supported_entries_n.sum()),
        "is_development_target": True,
        "development_targets": ";".join(sorted(DEVELOPMENT_TARGETS & set(by_dataset.target))),
    }
    result = pd.concat([by_dataset, pd.DataFrame([total])], ignore_index=True)

    existing = pd.read_csv(E5D / "map_impact_by_dataset.csv")
    merged = by_dataset.merge(
        existing[["target", "endpoint", "touched_entries_n", "g_supported_n"]],
        on=["target", "endpoint"],
        suffixes=("_new", "_e5d"),
        validate="one_to_one",
    )
    regression = {
        "max_abs_touched_entries_difference": int(
            np.max(np.abs(
                merged.touched_entries_n_new - merged.touched_entries_n_e5d
            ))
        ),
        "max_abs_g_supported_difference": int(
            np.max(np.abs(
                merged.g_supported_entries_n - merged.g_supported_n
            ))
        ),
        "n_datasets": int(len(merged)),
    }
    if regression["max_abs_touched_entries_difference"] != 0:
        raise RuntimeError(f"Touched-entry regression failed: {regression}")
    if regression["max_abs_g_supported_difference"] != 0:
        raise RuntimeError(f"Supported-entry regression failed: {regression}")
    return result, regression


def calibration_regression(calibration: pd.DataFrame) -> dict:
    reference = pd.read_csv(E9 / "calibration_summary.csv")
    keys = ["design", "endpoint", "population", "display_class"]
    mine = calibration.loc[
        calibration.population.isin(("all", "untouched", "touched", "G_supported", "P_supported"))
    ][[
        *keys,
        "n_predictions",
        "n_datasets",
        "mean_probability",
        "agreement_rate",
        "bias",
        "agreement_ci_low",
        "agreement_ci_high",
        "bias_ci_low",
        "bias_ci_high",
        "dataset_bias_median_min10",
        "random_effects_n_datasets",
        "random_effects_bias",
        "random_effects_ci_low",
        "random_effects_ci_high",
        "random_effects_tau2",
        "random_effects_i_squared",
        "random_effects_pi_low",
        "random_effects_pi_high",
    ]].rename(columns={"n_predictions": "n"})
    reference = reference.rename(columns={
        "mean_p": "mean_probability",
        "bias_median": "dataset_bias_median_min10",
    })
    reference = reference[keys + [
        column for column in mine.columns if column not in keys
    ]]
    merged = mine.merge(reference, on=keys, suffixes=("_new", "_e9"), validate="one_to_one")
    errors = {}
    for column in mine.columns:
        if column in keys:
            continue
        left = pd.to_numeric(merged[f"{column}_new"], errors="coerce")
        right = pd.to_numeric(merged[f"{column}_e9"], errors="coerce")
        both_nan = left.isna() & right.isna()
        difference = (left - right).abs()
        difference.loc[both_nan] = 0.0
        if difference.isna().any():
            raise RuntimeError(f"Non-comparable calibration values: {column}")
        errors[column] = float(difference.max())
    return {
        "n_overlapping_rows": int(len(merged)),
        "max_absolute_errors": errors,
        "max_absolute_error": float(max(errors.values())),
    }


def rule_test_regression(rule_output: pd.DataFrame) -> dict:
    expected = {
        ("document", "IC50"): (15, 15, 3.05e-5),
        ("document", "Ki"): (11, 9, 0.0049),
        ("entry", "IC50"): (14, 12, 0.0012),
        ("entry", "Ki"): (11, 7, 0.062),
    }
    rows = []
    for (design, endpoint), (n, positive, p_value) in expected.items():
        row = rule_output.loc[
            (rule_output.design == design)
            & (rule_output.endpoint == endpoint)
            & (rule_output.test == "rise")
            & (rule_output.unit == "dataset")
        ].iloc[0]
        if int(row.n_units) != n or int(row.positive_units) != positive:
            raise RuntimeError(
                f"Rise regression count failed: {design} {endpoint}: {row.to_dict()}"
            )
        if not np.isclose(row.one_sided_p, p_value, rtol=0.03, atol=2e-5):
            raise RuntimeError(
                f"Rise regression P failed: {design} {endpoint}: "
                f"{row.one_sided_p} vs {p_value}"
            )
        rows.append({
            "design": design,
            "endpoint": endpoint,
            "n_units": int(row.n_units),
            "positive_units": int(row.positive_units),
            "one_sided_p": float(row.one_sided_p),
        })
    return {"rows": rows}


def lodo_extrema(lodo: pd.DataFrame) -> dict:
    result = {}
    for (design, endpoint, quantity), group in lodo.groupby(
        ["design", "endpoint", "quantity"], sort=True
    ):
        if quantity == "rise_test":
            values = group.two_sided_p
            metric = "two_sided_p"
        else:
            values = group.bias if quantity.startswith("inactive_bias") else group.agreement_rate
            metric = "bias" if quantity.startswith("inactive_bias") else "agreement_rate"
        finite = values.replace([np.inf, -np.inf], np.nan).dropna()
        if finite.empty:
            continue
        min_index = finite.idxmin()
        max_index = finite.idxmax()
        key = f"{design}|{endpoint}|{quantity}"
        result[key] = {
            "metric": metric,
            "minimum": float(lodo.loc[min_index, metric]),
            "minimum_omitted_targets": lodo.loc[min_index, "omitted_targets"],
            "maximum": float(lodo.loc[max_index, metric]),
            "maximum_omitted_targets": lodo.loc[max_index, "omitted_targets"],
        }
    return result


def input_hashes() -> dict[str, str]:
    paths = [
        E5C / "predictions_e5c.csv",
        E5D / "map_impact_by_dataset.csv",
        E7 / "shared_document_cluster_members.csv",
        E9 / "calibration_summary.csv",
        E11 / "predictions_e11.csv",
        E2 / "target_inclusion.csv",
        E2_QA / "curated_records.pkl.gz",
    ]
    return {
        str(path.relative_to(ROOT)): sha256(path)
        for path in paths
    }


def main() -> None:
    started_at = datetime.now().astimezone()
    started_perf = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(Path(gettempdir()) / "e12a_mpl"))

    predictions, designs = load_designs()
    states, _ = rule_states(designs)
    write_csv(states, "rule_states_by_design.csv")

    cluster_members = pd.read_csv(E7 / "shared_document_cluster_members.csv")
    tests = rule_tests(designs, cluster_members)
    write_csv(tests, "rule_tests.csv")

    lodo, lodo_named, top_ic50 = build_lodo(designs)
    write_csv(lodo, "lodo.csv")
    write_csv(lodo_named, "lodo_named.csv")

    calibration = calibration_all_designs(designs)
    write_csv(calibration, "calibration_all_designs.csv")

    time_rule = time_split_rule(states, designs)
    write_csv(time_rule, "time_split_rule.csv")
    fig4d = fig4d_three_groups(designs)
    write_csv(fig4d, "fig4d_three_groups.csv")

    e11, e11_components = e11_sensitivity()
    write_csv(e11, "e11_sensitivity.csv")
    e11_components.to_csv(QA / "e11_shared_document_components.csv", index=False)

    full_map, map_regression = full_map_touched_units()
    write_csv(full_map, "full_map_touched_units.csv")

    calibration_check = calibration_regression(calibration)
    rule_check = rule_test_regression(tests)
    summary = {
        "analysis_role": ANALYSIS_ROLE,
        "started_at": started_at.isoformat(),
        "elapsed_seconds": time.perf_counter() - started_perf,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "frozen_chemgridmap": str(CHEMGRIDMAP_FILE),
        "loaded_prediction_rows": int(len(predictions)),
        "design_counts": {
            design: {
                endpoint: int(
                    (frame.endpoint == endpoint).sum()
                )
                for endpoint in ENDPOINTS
            }
            for design, frame in designs.items()
        },
        "rule_test_regression": rule_check,
        "calibration_regression": calibration_check,
        "full_map_regression": map_regression,
        "top_ic50_demoted": top_ic50,
        "lodo_extrema": lodo_extrema(lodo),
        "outputs": sorted(
            path.name for path in OUT.glob("*.csv")
        ),
    }
    write_json(QA / "summary.json", summary)
    write_json(OUT / "manifest.json", {
        **summary,
        "inputs": input_hashes(),
        "bootstrap_repetitions": BOOT_REPETITIONS,
        "e9_overlap_bootstrap_seed": E9_BOOT_SEED,
        "new_bootstrap_seed": NEW_BOOT_SEED,
        "analysis_role": ANALYSIS_ROLE,
        "network_used": False,
        "git_write_operations": False,
    })
    print(json.dumps(json_safe(summary), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
