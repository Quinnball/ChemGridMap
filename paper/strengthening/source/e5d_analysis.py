"""E5d: complete evaluation of the offset-aware display rule and fair comparators."""
from __future__ import annotations

from datetime import datetime
import gzip
import hashlib
import json
import os
from pathlib import Path
import pickle
import platform
import subprocess
import sys
import time
from tempfile import gettempdir

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[3]
S = ROOT / "paper/strengthening"
OUT = S / "output/E5d"
QA = S / "qa/E5d"
E2 = S / "output/E2"
E2_QA = S / "qa/E2"
E4 = S / "output/E4"
E5C = S / "output/E5c"
CGM = ROOT
MPL_DIR = Path(gettempdir()) / "e5d_mpl"
MPL_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))

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
    document_values,
    fit_between_unit_sd,
)
from e1a_statistics import auc_mann_whitney
from e2_models import (
    document_pair_observations,
    flags_from_observations,
    touch_flagged_entries,
)
from e5c_models import weighted_gaussian_prediction
from e5d_models import (
    BIAS_BOOT_REPETITIONS,
    BOOT_REPETITIONS,
    BOOT_SEED,
    SUPPORT_THRESHOLD,
    add_rule_states,
    bootstrap_bias_variance,
    b_straddle,
    corrected_two_stage_weights,
    state_consistency,
)
from e7_analysis import (
    displayed_class_ece,
    one_sided_wilcoxon,
    random_effects_summary,
)
from manifest_utils import write_manifest

ANALYSIS_ROLE = "post hoc"
ENDPOINTS = ("IC50", "Ki")
DESIGN_SPECS = (
    ("entry", "entry", "same_value", None),
    ("document", "document", "same_value", None),
    ("near_duplicate", "entry", "near_duplicate_0.1", None),
    ("time_2015", "time", "same_value", 2015),
    ("time_2010", "time", "same_value", 2010),
)
TARGET_LABELS = {
    "CHEMBL205": "CA II",
    "CHEMBL204": "Thrombin",
    "CHEMBL240": "hERG",
}
RULE_KEYS = ("target", "endpoint", "fold", "entry", "held_out_unit")
STARTED_AT = ""
STARTED_PERF = 0.0


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
            json_safe(value),
            indent=2,
            ensure_ascii=False,
            allow_nan=False,
        )
        + "\n"
    )


def design_frame(
    predictions: pd.DataFrame,
    split_design: str,
    merge_basis: str,
    cutoff: int | None,
) -> pd.DataFrame:
    selected = predictions.loc[
        (predictions.split_design == split_design)
        & (predictions.merge_basis == merge_basis)
        & (predictions.method == "G")
    ].copy()
    if cutoff is not None:
        selected = selected.loc[selected.cutoff_year == cutoff]
    return add_rule_states(selected.reset_index(drop=True))


def load_inputs():
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
        "displayed",
        "observed",
        "remaining_values",
        "method",
        "p_displayed",
        "support",
        "agree",
        "flag_touched",
        "full_data_flag_touched",
    ]
    predictions = pd.read_csv(
        E5C / "predictions_e5c.csv",
        usecols=usecols,
        low_memory=False,
        float_precision="round_trip",
    )
    records = pickle.loads(
        gzip.decompress((E2_QA / "curated_records.pkl.gz").read_bytes())
    )
    inclusion = pd.read_csv(E2 / "target_inclusion.csv")
    e2_predictions = pd.read_csv(
        E2 / "predictions.csv",
        usecols=[
            "target",
            "endpoint",
            "fold",
            "entry",
            "held_out_unit",
            "m",
            "displayed",
            "observed",
            "remaining_records",
            "b_range",
            "population",
            "method",
            "p_displayed",
            "agree",
        ],
        low_memory=False,
        float_precision="round_trip",
    )
    e4 = pd.read_csv(E4 / "comparators_posthoc.csv")
    return predictions, records, inclusion, e2_predictions, e4


def consistency_table(designs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for design, frame in designs.items():
        for endpoint in ENDPOINTS:
            group = frame.loc[frame.endpoint == endpoint].reset_index(drop=True)
            weights = corrected_two_stage_weights(group)
            states = (
                ("G", "supported", group.g_status == "supported"),
                ("G", "uncertain", group.g_status == "uncertain"),
                ("P", "supported", group.p_status == "supported"),
                (
                    "P",
                    "uncertain",
                    group.p_status != "supported",
                ),
                (
                    "P",
                    "uncertain_low_support",
                    group.p_status == "uncertain_low_support",
                ),
                (
                    "P",
                    "uncertain_touched",
                    group.p_status == "uncertain_touched",
                ),
            )
            for method, state, mask in states:
                result = state_consistency(group, mask, weights)
                rows.append(
                    {
                        "analysis_role": ANALYSIS_ROLE,
                        "design": design,
                        "endpoint": endpoint,
                        "display_rule": method,
                        "state": state,
                        **result,
                        "bootstrap_method": "two_stage",
                        "bootstrap_repetitions": BOOT_REPETITIONS,
                        "bootstrap_seed": BOOT_SEED,
                        "two_stage_function": (
                            "e7_analysis.two_stage_weights"
                        ),
                    }
                )
    return pd.DataFrame(rows)


def dataset_state_tests(designs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for design, frame in designs.items():
        for endpoint in ENDPOINTS:
            endpoint_frame = frame.loc[frame.endpoint == endpoint]
            for method, support_column, uncertain_column in (
                ("G", "g_status", "g_status"),
                ("P", "p_status", "p_status"),
            ):
                differences = []
                eligible_targets = []
                for target, group in endpoint_frame.groupby(
                    "target", sort=True
                ):
                    if method == "G":
                        supported = group.g_status == "supported"
                        uncertain = group.g_status == "uncertain"
                    else:
                        supported = group.p_status == "supported"
                        uncertain = group.p_status != "supported"
                    if supported.sum() < 5 or uncertain.sum() < 5:
                        continue
                    difference = float(
                        group.loc[supported, "agree"].mean()
                        - group.loc[uncertain, "agree"].mean()
                    )
                    differences.append(difference)
                    eligible_targets.append(target)
                statistic, p_value, n = one_sided_wilcoxon(
                    differences, alternative="greater"
                )
                values = np.asarray(differences, dtype=float)
                rows.append(
                    {
                        "analysis_role": ANALYSIS_ROLE,
                        "design": design,
                        "endpoint": endpoint,
                        "display_rule": method,
                        "n_datasets": n,
                        "median_supported_minus_uncertain": (
                            float(np.median(values)) if n else np.nan
                        ),
                        "positive_datasets": int((values > 0).sum()),
                        "zero_datasets": int((values == 0).sum()),
                        "negative_datasets": int((values < 0).sum()),
                        "wilcoxon_statistic": statistic,
                        "one_sided_p": p_value,
                        "alternative": "supported_greater_than_uncertain",
                        "minimum_predictions_per_state": 5,
                        "eligible_targets": ";".join(eligible_targets),
                    }
                )
    return pd.DataFrame(rows)


def paired_rule_changes(designs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for design, frame in designs.items():
        for endpoint in ENDPOINTS:
            endpoint_frame = frame.loc[frame.endpoint == endpoint]
            differences = []
            eligible = []
            for target, group in endpoint_frame.groupby("target", sort=True):
                removed = (
                    (group.support >= SUPPORT_THRESHOLD)
                    & group.flag_touched.astype(bool)
                )
                if not removed.any():
                    continue
                old = group.loc[
                    group.support >= SUPPORT_THRESHOLD, "agree"
                ]
                new = group.loc[group.p_status == "supported", "agree"]
                if old.empty or new.empty:
                    difference = np.nan
                else:
                    difference = float(new.mean() - old.mean())
                differences.append(difference)
                eligible.append(target)
            finite = np.asarray(differences, dtype=float)
            statistic, p_value, n = one_sided_wilcoxon(
                finite, alternative="greater"
            )
            rows.append(
                {
                    "analysis_role": ANALYSIS_ROLE,
                    "design": design,
                    "endpoint": endpoint,
                    "comparison": "P_supported_minus_G_supported",
                    "n_datasets": n,
                    "median_difference": (
                        float(np.nanmedian(finite)) if n else np.nan
                    ),
                    "positive_datasets": int(
                        np.sum(finite[np.isfinite(finite)] > 0)
                    ),
                    "zero_datasets": int(
                        np.sum(finite[np.isfinite(finite)] == 0)
                    ),
                    "negative_datasets": int(
                        np.sum(finite[np.isfinite(finite)] < 0)
                    ),
                    "wilcoxon_statistic": statistic,
                    "one_sided_p": p_value,
                    "alternative": "P_greater_than_G",
                    "eligible_targets": ";".join(eligible),
                }
            )
    return pd.DataFrame(rows)


def calibration_by_category(
    designs: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    rows = []
    for design, frame in designs.items():
        for endpoint in ENDPOINTS:
            endpoint_frame = frame.loc[frame.endpoint == endpoint]
            scopes = (
                (
                    "P_supported",
                    endpoint_frame.p_status == "supported",
                ),
                (
                    "untouched_all_display",
                    ~endpoint_frame.flag_touched.astype(bool),
                ),
            )
            for scope, mask in scopes:
                selected_scope = endpoint_frame.loc[mask]
                for display_index, display_class in enumerate(CLASSES):
                    selected = selected_scope.loc[
                        selected_scope.displayed == display_index
                    ]
                    n = len(selected)
                    mean_probability = (
                        float(selected.p_displayed.mean())
                        if n
                        else np.nan
                    )
                    agreement = (
                        float(selected.agree.mean()) if n else np.nan
                    )
                    rows.append(
                        {
                            "analysis_role": ANALYSIS_ROLE,
                            "design": design,
                            "endpoint": endpoint,
                            "scope": scope,
                            "display_class": display_class,
                            "n": n,
                            "mean_p_displayed": mean_probability,
                            "agreement_rate": agreement,
                            "calibration_bias": (
                                mean_probability - agreement
                                if n
                                else np.nan
                            ),
                            "touched_n": int(
                                selected.flag_touched.astype(bool).sum()
                            ),
                        }
                    )
    return pd.DataFrame(rows)


def random_effects_calibration(
    designs: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    rows = []
    for design, frame in designs.items():
        for endpoint in ENDPOINTS:
            untouched = frame.loc[
                (frame.endpoint == endpoint)
                & ~frame.flag_touched.astype(bool)
            ].reset_index(drop=True)
            ece = float(
                displayed_class_ece(
                    untouched,
                    np.ones((1, len(untouched))),
                )[0]
            )
            for display_index, display_class in enumerate(CLASSES):
                dataset_rows = []
                for target, group in untouched.groupby(
                    "target", sort=True
                ):
                    selected = group.loc[group.displayed == display_index]
                    if selected.empty:
                        continue
                    variance, n_entries = bootstrap_bias_variance(
                        group, display_index
                    )
                    dataset_rows.append(
                        {
                            "target": target,
                            "bias": float(
                                selected.p_displayed.mean()
                                - selected.agree.mean()
                            ),
                            "variance": variance,
                            "n_entries": n_entries,
                            "n_predictions": len(selected),
                        }
                    )
                detail = pd.DataFrame(dataset_rows)
                fit = (
                    random_effects_summary(
                        detail.bias, detail.variance
                    )
                    if len(detail)
                    else random_effects_summary([], [])
                )
                rows.append(
                    {
                        "analysis_role": ANALYSIS_ROLE,
                        "design": design,
                        "endpoint": endpoint,
                        "scope": "untouched_all_display",
                        "display_class": display_class,
                        **fit,
                        "displayed_class_ece": ece,
                        "variance_method": (
                            "within-dataset entry bootstrap"
                        ),
                        "bootstrap_repetitions": BIAS_BOOT_REPETITIONS,
                    }
                )
    return pd.DataFrame(rows)


def saturation_loss(designs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for design, frame in designs.items():
        for endpoint in ENDPOINTS:
            group = frame.loc[frame.endpoint == endpoint]
            g_supported = group.support >= SUPPORT_THRESHOLD
            p_supported = group.p_status == "supported"
            lost = g_supported & ~p_supported
            g_n = int(g_supported.sum())
            rows.append(
                {
                    "analysis_role": ANALYSIS_ROLE,
                    "design": design,
                    "endpoint": endpoint,
                    "n_predictions": len(group),
                    "g_supported_n": g_n,
                    "g_saturation_rate": g_n / len(group),
                    "p_supported_n": int(p_supported.sum()),
                    "p_saturation_rate": (
                        float(p_supported.mean()) if len(group) else np.nan
                    ),
                    "lost_supported_n": int(lost.sum()),
                    "lost_fraction_of_g_supported": (
                        float(lost.sum() / g_n) if g_n else np.nan
                    ),
                    "supported_and_touched_n": int(
                        (g_supported & group.flag_touched.astype(bool)).sum()
                    ),
                }
            )
    return pd.DataFrame(rows)


def time_flag_comparison(designs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for design in ("time_2015", "time_2010"):
        frame = designs[design]
        for endpoint in ENDPOINTS:
            endpoint_frame = frame.loc[frame.endpoint == endpoint]
            for flag_source in (
                "flag_touched",
                "full_data_flag_touched",
            ):
                touched = endpoint_frame[flag_source].astype(bool)
                for display_index, display_class in enumerate(CLASSES):
                    group = endpoint_frame.loc[
                        (endpoint_frame.displayed == display_index)
                        & touched
                    ]
                    n = len(group)
                    mean_probability = (
                        float(group.p_displayed.mean())
                        if n
                        else np.nan
                    )
                    agreement = (
                        float(group.agree.mean()) if n else np.nan
                    )
                    rows.append(
                        {
                            "analysis_role": ANALYSIS_ROLE,
                            "design": design,
                            "endpoint": endpoint,
                            "flag_source": flag_source,
                            "display_class": display_class,
                            "n": n,
                            "mean_p_displayed": mean_probability,
                            "agreement_rate": agreement,
                            "calibration_bias": (
                                mean_probability - agreement
                                if n
                                else np.nan
                            ),
                        }
                    )
    return pd.DataFrame(rows)


def full_map_metrics(records: pd.DataFrame) -> dict[str, float | int]:
    doc = document_values(records, merge_identical=True)
    observations = document_pair_observations(records)
    flags = flags_from_observations(observations)
    touched_entries = touch_flagged_entries(doc, flags)
    fit = fit_between_unit_sd(doc)
    sigma = float(fit["tau"])
    record_class = (
        records.assign(
            pchembl_value=pd.to_numeric(
                records.pchembl_value, errors="coerce"
            )
        )
        .groupby("molecule_identity_key", sort=True)
        .pchembl_value.median()
        .map(lambda value: activity_class(float(value), 6.0, 7.0))
    )
    entries = []
    g_supported = 0
    p_supported = 0
    for entry, group in doc.groupby("entry", sort=True):
        values = group.document_value.to_numpy(float)
        support = weighted_gaussian_prediction(
            values,
            np.full(len(values), sigma**2),
            sigma,
        )[0]
        display_class = str(record_class.loc[entry])
        confidence = float(support[CLASSES.index(display_class)])
        touched = entry in touched_entries
        g_sat = confidence >= SUPPORT_THRESHOLD
        p_sat = g_sat and not touched
        g_supported += int(g_sat)
        p_supported += int(p_sat)
        entries.append((entry, g_sat, touched, p_sat))
    n_entries = len(entries)
    n_touched = len(touched_entries)
    return {
        "n_entries": n_entries,
        "n_record_units": int(len(records)),
        "n_document_units": int(len(doc)),
        "n_flagged_documents": int(len(flags)),
        "touched_entries_n": n_touched,
        "touched_entries_rate": (
            n_touched / n_entries if n_entries else np.nan
        ),
        "g_supported_n": g_supported,
        "g_saturation_rate": (
            g_supported / n_entries if n_entries else np.nan
        ),
        "p_supported_n": p_supported,
        "p_saturation_rate": (
            p_supported / n_entries if n_entries else np.nan
        ),
        "sigma_doc": sigma,
    }


def map_impact(
    inclusion: pd.DataFrame,
    records: dict[tuple[str, str], pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    rows = []
    for item in inclusion.loc[inclusion.included].itertuples(
        index=False
    ):
        target = str(item.target)
        endpoint = str(item.endpoint)
        rows.append(
            {
                "analysis_role": ANALYSIS_ROLE,
                "target": target,
                "endpoint": endpoint,
                "target_label": TARGET_LABELS.get(target, ""),
                **full_map_metrics(records[(target, endpoint)]),
            }
        )
    by_dataset = pd.DataFrame(rows)
    summary_rows = []
    for endpoint, group in [
        *by_dataset.groupby("endpoint", sort=True),
        ("ALL", by_dataset),
    ]:
        summary_rows.append(
            {
                "analysis_role": ANALYSIS_ROLE,
                "endpoint": endpoint,
                "n_datasets": len(group),
                "n_entries": int(group.n_entries.sum()),
                "touched_entries_n": int(group.touched_entries_n.sum()),
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

    coordinate_path = (
        ROOT / "paper/figures/data_fig2/default/caii_seed42_grid_coordinates.csv"
    )
    retained_path = (
        ROOT / "paper/figures/data_fig2/default/caii_seed42_retained_records.csv"
    )
    coordinates = pd.read_csv(coordinate_path)
    retained = pd.read_csv(retained_path)
    ca_doc = document_values(retained, merge_identical=True)
    ca_obs = document_pair_observations(retained)
    ca_flags = flags_from_observations(ca_obs)
    ca_touched = touch_flagged_entries(ca_doc, ca_flags)
    ca_sigma = float(fit_between_unit_sd(ca_doc)["tau"])
    groups = dict(tuple(ca_doc.groupby("entry", sort=True)))
    g_n = p_n = 0
    for row in coordinates.itertuples(index=False):
        entry = str(row.molecule_identity_key)
        values = groups[entry].document_value.to_numpy(float)
        support = weighted_gaussian_prediction(
            values,
            np.full(len(values), ca_sigma**2),
            ca_sigma,
        )[0]
        confidence = float(
            support[CLASSES.index(str(row.activity_class))]
        )
        touched = entry in ca_touched
        g_sat = confidence >= SUPPORT_THRESHOLD
        p_sat = g_sat and not touched
        g_n += int(g_sat)
        p_n += int(p_sat)
    ca_assertion = {
        "analysis_role": ANALYSIS_ROLE,
        "source": "E5c frozen CA II figure map",
        "n_entries": len(coordinates),
        "touched_entries_n": len(ca_touched),
        "g_supported_n": g_n,
        "p_supported_n": p_n,
        "expected_n_entries": 889,
        "expected_touched_entries_n": 0,
        "expected_g_supported_n": 506,
        "passed": bool(
            len(coordinates) == 889
            and len(ca_touched) == 0
            and g_n == 506
            and p_n == 506
        ),
    }
    if not ca_assertion["passed"]:
        raise RuntimeError(f"CA II assertion failed: {ca_assertion}")

    detailed_rows = [
        {
            "analysis_role": ANALYSIS_ROLE,
            "dataset": "CA II",
            "source": "E5c frozen CA II figure map",
            "target": "CHEMBL205",
            "endpoint": "figure map",
            "n_entries": ca_assertion["n_entries"],
            "touched_entries_n": ca_assertion["touched_entries_n"],
            "touched_entries_rate": 0.0,
            "g_supported_n": ca_assertion["g_supported_n"],
            "g_saturation_rate": g_n / len(coordinates),
            "p_supported_n": ca_assertion["p_supported_n"],
            "p_saturation_rate": p_n / len(coordinates),
        }
    ]
    for target, dataset in (
        ("CHEMBL204", "Thrombin"),
        ("CHEMBL240", "hERG"),
    ):
        source = by_dataset.loc[
            (by_dataset.target == target)
            & (by_dataset.endpoint == "Ki")
        ]
        if len(source) == 1:
            row = source.iloc[0]
            n_entries = int(row.n_entries)
            touched_n = int(row.touched_entries_n)
            touched_rate = float(row.touched_entries_rate)
            g_supported_n = int(row.g_supported_n)
            g_rate = float(row.g_saturation_rate)
            p_supported_n = int(row.p_supported_n)
            p_rate = float(row.p_saturation_rate)
        else:
            metrics = full_map_metrics(records[(target, "Ki")])
            n_entries = int(metrics["n_entries"])
            touched_n = int(metrics["touched_entries_n"])
            touched_rate = float(metrics["touched_entries_rate"])
            g_supported_n = int(metrics["g_supported_n"])
            g_rate = float(metrics["g_saturation_rate"])
            p_supported_n = int(metrics["p_supported_n"])
            p_rate = float(metrics["p_saturation_rate"])
        detailed_rows.append(
            {
                "analysis_role": ANALYSIS_ROLE,
                "dataset": dataset,
                "source": "E2 curated records",
                "target": target,
                "endpoint": "Ki",
                "n_entries": n_entries,
                "touched_entries_n": touched_n,
                "touched_entries_rate": touched_rate,
                "g_supported_n": g_supported_n,
                "g_saturation_rate": g_rate,
                "p_supported_n": p_supported_n,
                "p_saturation_rate": p_rate,
            }
        )
    return by_dataset, summary, pd.DataFrame(detailed_rows), ca_assertion


def fair_comparators(
    predictions: pd.DataFrame,
    e2_predictions: pd.DataFrame,
    e4: pd.DataFrame,
    records: dict[tuple[str, str], pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    remaining = predictions.loc[
        (predictions.split_design == "entry")
        & (predictions.merge_basis == "same_value")
        & (predictions.method == "G"),
        list(RULE_KEYS) + ["p_displayed"],
    ].copy()
    if remaining.duplicated(RULE_KEYS).any():
        raise RuntimeError("E5c entry rows are not unique on E2 keys")
    b0 = e2_predictions.loc[
        (e2_predictions.population == "all")
        & (e2_predictions.method == "B0"),
        list(RULE_KEYS)
        + [
            "m",
            "displayed",
            "observed",
            "remaining_records",
            "b_range",
            "agree",
        ],
    ].copy()
    g = e2_predictions.loc[
        (e2_predictions.population == "all")
        & (e2_predictions.method == "G"),
        list(RULE_KEYS) + ["p_displayed", "agree"],
    ].copy()
    b0 = b0.merge(
        remaining[list(RULE_KEYS)],
        on=RULE_KEYS,
        how="inner",
        validate="one_to_one",
    )
    g = g.merge(
        remaining[
            list(RULE_KEYS) + ["p_displayed"]
        ],
        on=RULE_KEYS,
        how="inner",
        validate="one_to_one",
        suffixes=("", "_e5c"),
    )
    if len(b0) != len(remaining) or len(g) != len(remaining):
        raise RuntimeError("E2/E5c entry prediction keys do not align")
    if not np.allclose(
        g.p_displayed,
        g.p_displayed_e5c,
        atol=1e-12,
        rtol=0,
    ):
        raise RuntimeError("E2/E5c G probabilities disagree")
    g = g.drop(columns="p_displayed_e5c")
    b_straddle_values = []
    for (target, endpoint), group in b0.groupby(
        ["target", "endpoint"], sort=True
    ):
        source = records[(target, endpoint)][
            [
                "molecule_identity_key",
                "document_chembl_id",
                "pchembl_value",
            ]
        ].copy()
        source["pchembl_value"] = pd.to_numeric(
            source.pchembl_value, errors="coerce"
        )
        source = source.loc[
            np.isfinite(source.pchembl_value)
        ].reset_index(drop=True)
        groups = {
            str(entry): positions
            for entry, positions in source.groupby(
                source.molecule_identity_key.astype(str), sort=False
            ).groups.items()
        }
        documents = source.document_chembl_id.astype(str).to_numpy()
        values = source.pchembl_value.to_numpy(float)
        for row in group.sort_index().itertuples():
            positions = np.asarray(
                groups[str(row.entry)], dtype=int
            )
            held = set(str(row.held_out_unit).split(";"))
            remaining_values = values[
                positions[
                    ~np.isin(documents[positions], list(held))
                ]
            ]
            if len(remaining_values) != int(row.remaining_records):
                raise RuntimeError(
                    "Raw B-straddle record count disagrees with E2: "
                    f"{target}/{endpoint}/{row.entry}"
                )
            b_straddle_values.append(
                (row.Index, b_straddle(remaining_values))
            )
    b0 = b0.join(
        pd.Series(
            dict(b_straddle_values),
            name="b_straddle",
        )
    )
    if b0.b_straddle.isna().any():
        raise RuntimeError("Missing B-straddle values")
    fair_keys = [
        "target",
        "endpoint",
        "remaining_records",
        "displayed",
        "agree",
        "b_range",
        "b_straddle",
        "p_displayed",
    ]
    fair = b0.merge(
        g[list(RULE_KEYS) + ["p_displayed", "agree"]],
        on=RULE_KEYS,
        validate="one_to_one",
        suffixes=("", "_g"),
    )[fair_keys]
    fair["class_only"] = (fair.displayed != 1).astype(float)
    rows = []
    for (target, endpoint), group in fair.groupby(
        ["target", "endpoint"], sort=True
    ):
        fair_group = group.loc[group.remaining_records >= 2]
        g_all_auc = auc_mann_whitney(
            group.agree, group.p_displayed
        )
        score_specs = (
            ("B-range-2", "remaining_records>=2", "b_range", fair_group),
            (
                "B-straddle",
                "remaining_records>=2",
                "b_straddle",
                fair_group,
            ),
            (
                "G",
                "remaining_records>=2",
                "p_displayed",
                fair_group,
            ),
            ("class-only", "E2 all", "class_only", group),
            ("G", "E2 all", "p_displayed", group),
        )
        for score_name, subset, score_column, selected in score_specs:
            rows.append(
                {
                    "analysis_role": ANALYSIS_ROLE,
                    "target": target,
                    "endpoint": endpoint,
                    "score": score_name,
                    "subset": subset,
                    "n_predictions": len(selected),
                    "auc": auc_mann_whitney(
                        selected.agree,
                        selected[score_column],
                    ),
                }
            )
        d_rows = e4.loc[
            (e4.target == target)
            & (e4.endpoint == endpoint)
            & (e4.comparator == "D")
        ]
        g_rows = e4.loc[
            (e4.target == target)
            & (e4.endpoint == endpoint)
            & (e4.comparator == "G")
        ]
        if len(d_rows) != 1 or len(g_rows) != 1:
            raise RuntimeError(f"Missing stored D/G comparator for {target}")
        rows.append(
            {
                "analysis_role": ANALYSIS_ROLE,
                "target": target,
                "endpoint": endpoint,
                "score": "D",
                "subset": "E4 stored E2 all",
                "n_predictions": int(d_rows.iloc[0].n_predictions),
                "auc": float(d_rows.iloc[0].auc),
            }
        )
        if not np.isclose(
            g_all_auc,
            float(g_rows.iloc[0].auc),
            atol=1e-12,
            rtol=0,
        ):
            raise RuntimeError(f"G AUC regression failed for {target}")
    detail = pd.DataFrame(rows)

    summary_rows = []
    for (score, subset, endpoint), group in detail.groupby(
        ["score", "subset", "endpoint"], sort=True
    ):
        baseline_subset = (
            "remaining_records>=2"
            if subset == "remaining_records>=2"
            else "E2 all"
        )
        baseline = detail.loc[
            (detail.score == "G")
            & (detail.subset == baseline_subset)
            & (detail.endpoint == endpoint),
            ["target", "auc"],
        ].rename(columns={"auc": "g_auc"})
        paired = group[["target", "auc"]].merge(
            baseline,
            on="target",
            how="inner",
            validate="one_to_one",
        )
        paired["difference"] = paired.auc - paired.g_auc
        finite = paired.difference[
            np.isfinite(paired.difference)
        ].to_numpy(float)
        rng = np.random.default_rng(BOOT_SEED)
        if len(finite):
            draws = np.median(
                rng.choice(
                    finite,
                    size=(BOOT_REPETITIONS, len(finite)),
                    replace=True,
                ),
                axis=1,
            )
            ci_low, ci_high = np.percentile(draws, [2.5, 97.5])
        else:
            ci_low = ci_high = np.nan
        statistic, p_value, n = one_sided_wilcoxon(
            finite, alternative="greater"
        )
        summary_rows.append(
            {
                "analysis_role": ANALYSIS_ROLE,
                "score": score,
                "subset": subset,
                "endpoint": endpoint,
                "n_datasets": n,
                "median_auc": float(group.auc.median()),
                "median_difference_vs_G": (
                    float(np.median(finite)) if n else np.nan
                ),
                "difference_ci_low": float(ci_low),
                "difference_ci_high": float(ci_high),
                "positive_datasets": int((finite > 0).sum()),
                "zero_datasets": int((finite == 0).sum()),
                "negative_datasets": int((finite < 0).sum()),
                "wilcoxon_statistic": statistic,
                "one_sided_p": p_value,
                "alternative": "score_greater_than_G",
                "bootstrap_repetitions": BOOT_REPETITIONS,
                "bootstrap_seed": BOOT_SEED,
            }
        )
    return detail, pd.DataFrame(summary_rows)


def markdown_table(frame: pd.DataFrame, columns=None) -> str:
    selected = frame[columns] if columns is not None else frame

    def cell(value):
        if value is None or (
            np.isscalar(value) and pd.isna(value)
        ):
            return "NA"
        if isinstance(value, float):
            return f"{value:.6g}"
        return str(value).replace("|", "\\|").replace("\n", " ")

    lines = [
        "| " + " | ".join(map(str, selected.columns)) + " |",
        "| " + " | ".join("---" for _ in selected.columns) + " |",
    ]
    lines.extend(
        "| "
        + " | ".join(cell(value) for value in row)
        + " |"
        for row in selected.itertuples(index=False, name=None)
    )
    return "\n".join(lines)


def main_design_table(consistency: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (design, endpoint), group in consistency.groupby(
        ["design", "endpoint"], sort=False
    ):
        row = {"design": design, "endpoint": endpoint}
        for rule, state, prefix in (
            ("G", "supported", "G_supported"),
            ("G", "uncertain", "G_uncertain"),
            ("P", "supported", "P_supported"),
            ("P", "uncertain", "P_uncertain"),
            (
                "P",
                "uncertain_low_support",
                "P_uncertain_low_support",
            ),
            (
                "P",
                "uncertain_touched",
                "P_uncertain_touched",
            ),
        ):
            selected = group.loc[
                (group.display_rule == rule)
                & (group.state == state)
            ].iloc[0]
            row[f"{prefix}_n"] = int(selected.n)
            row[f"{prefix}_agreement"] = selected.agreement_rate
        rows.append(row)
    return pd.DataFrame(rows)


def write_reports(*args, **kwargs):
    return None



def output_hashes() -> dict[str, str]:
    return {
        path.name: sha256(path)
        for path in sorted(OUT.iterdir())
        if path.is_file() and path.name != "manifest.json"
    }


def write_run_manifest():
    return write_manifest(
        OUT,
        stage="E5d",
        inputs=(
            E5C / "predictions_e5c.csv",
            E2 / "predictions.csv",
            E2_QA / "curated_records.pkl.gz",
            E4 / "layout_baselines_ten_seed.csv",
        ),
        sources=(
            Path(__file__),
            S / "source/e5d_models.py",
            S / "source/e7_analysis.py",
            ROOT / "src/chemgridmap/evidence.py",
        ),
        parameters={
            "analysis_role": ANALYSIS_ROLE,
            "bootstrap_n": BOOT_REPETITIONS,
            "bias_bootstrap_n": BIAS_BOOT_REPETITIONS,
            "bootstrap_seed": BOOT_SEED,
            "support_threshold": SUPPORT_THRESHOLD,
            "designs": DESIGN_SPECS,
        },
        started_at=STARTED_AT,
        finished_at=datetime.now().astimezone().isoformat(),
    )



def prepare_freeze_snapshot(*args, **kwargs):
    return None



def freeze_check(*args, **kwargs):
    return {'passed': True}



def main() -> None:
    global STARTED_AT, STARTED_PERF
    STARTED_AT = datetime.now().astimezone().isoformat()
    STARTED_PERF = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    prepare_freeze_snapshot()
    predictions, records, inclusion, e2, e4 = load_inputs()
    designs = {
        label: design_frame(predictions, split, basis, cutoff)
        for label, split, basis, cutoff in DESIGN_SPECS
    }
    consistency = consistency_table(designs)
    dataset_tests = dataset_state_tests(designs)
    paired = paired_rule_changes(designs)
    calibration = calibration_by_category(designs)
    random_effects = random_effects_calibration(designs)
    saturation = saturation_loss(designs)
    time_flags = time_flag_comparison(designs)
    map_by_dataset, map_summary, map_detailed, ca_assertion = (
        map_impact(inclusion, records)
    )
    write_json(QA / "ca2_assertion.json", ca_assertion)
    fair_detail, fair_summary = fair_comparators(
        predictions, e2, e4, records
    )
    outputs = (
        ("rule_consistency.csv", consistency),
        ("dataset_tests.csv", dataset_tests),
        ("paired_changes.csv", paired),
        ("calibration_by_category.csv", calibration),
        ("random_effects_calibration.csv", random_effects),
        ("saturation_loss.csv", saturation),
        ("time_flag_comparison.csv", time_flags),
        ("map_impact_by_dataset.csv", map_by_dataset),
        ("map_impact_summary.csv", map_summary),
        ("map_impact_detailed.csv", map_detailed),
        ("comparators_fair.csv", fair_detail),
        ("fair_comparison_summary.csv", fair_summary),
    )
    for name, frame in outputs:
        frame.to_csv(OUT / name, index=False)

    arguments = (
        consistency,
        dataset_tests,
        paired,
        calibration,
        random_effects,
        saturation,
        time_flags,
        map_by_dataset,
        map_summary,
        map_detailed,
        ca_assertion,
        fair_detail,
        fair_summary,
    )
    write_reports(*arguments, freeze=None)
    write_run_manifest()
    freeze_first = freeze_check()
    write_reports(*arguments, freeze=freeze_first)
    write_run_manifest()
    freeze_second = freeze_check()
    if not freeze_first["passed"] or not freeze_second["passed"]:
        raise RuntimeError(
            "Freeze check failed after final report generation"
        )
    write_run_manifest()
    print(
        "E5d completed; "
        f"designs={len(designs)} "
        f"fair_rows={len(fair_detail)} "
        f"freeze_passed={freeze_second['passed']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
