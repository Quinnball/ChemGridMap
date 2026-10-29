"""E4b: post hoc temporal-split sensitivity analyses."""
from __future__ import annotations

import gzip
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

from chemgridmap.evidence import (
    CLASSES,
    class_probabilities,
    document_values,
    fit_between_unit_sd,
)
from e2_models import CLASS_BOUNDARIES
from e2b_analysis import (
    load_original_ic50_records,
    one_sided_wilcoxon,
    time_split_metrics,
)

E2_OUTPUT = S / "output/E2"
E2_QA = S / "qa/E2"
OUT = S / "output/E4"
OUTPUT = OUT / "time_split_posthoc.csv"
BASELINE = E2_QA / "time_split_predictions.csv"
RECORDS = E2_QA / "curated_records.pkl.gz"
INCLUSION = E2_OUTPUT / "target_inclusion.csv"

TIME_CUTOFFS = (2015, 2010)
ENDPOINTS = ("IC50", "Ki")
ORIGINAL_IC50_TARGETS = ("CHEMBL205", "CHEMBL204", "CHEMBL240")
SUPPORT_THRESHOLD = 0.8
NEAR_DUPLICATE_LOG10_GAP = 0.1
_RECORDS_CACHE = None
_TEMPORAL_FRAME_CACHE = {}
_POSTHOC_CACHE = None

COLUMNS = [
    "analysis_role",
    "analysis_section",
    "variant",
    "record_type",
    "cutoff_year",
    "endpoint",
    "target",
    "sample_filter",
    "display_status",
    "m_past",
    "n_predictions",
    "n_entries",
    "n_targets",
    "n_datasets",
    "supported_n",
    "supported_agreement_rate",
    "uncertain_n",
    "uncertain_agreement_rate",
    "dataset_difference",
    "median_dataset_difference",
    "positive_difference_datasets",
    "wilcoxon_statistic",
    "wilcoxon_p_one_sided",
    "mean_p_displayed",
    "calibration_bias",
    "auc",
    "dataset_auc_median",
    "dataset_auc_q1",
    "dataset_auc_q3",
    "dataset_auc_iqr",
]


def load_records() -> dict[tuple[str, str], pd.DataFrame]:
    global _RECORDS_CACHE
    if _RECORDS_CACHE is not None:
        return _RECORDS_CACHE
    records = pickle.loads(gzip.decompress(RECORDS.read_bytes()))
    records.update(load_original_ic50_records())
    _RECORDS_CACHE = records
    return records


def temporal_keys(
    inclusion: pd.DataFrame,
    records: dict[tuple[str, str], pd.DataFrame],
) -> list[tuple[str, str]]:
    keys = list(
        inclusion.loc[inclusion.included, ["target", "endpoint"]]
        .itertuples(index=False, name=None)
    )
    for target in ORIGINAL_IC50_TARGETS:
        key = (target, "IC50")
        if key not in keys:
            if key not in records:
                raise KeyError(f"Missing temporal records: {key}")
            keys.append(key)
    return keys


def unit_year_map(records: pd.DataFrame, year_mode: str) -> pd.Series:
    """Return a document-year map; unit years aggregate these document years."""
    if year_mode not in {"latest", "earliest"}:
        raise ValueError(f"Unknown year mode: {year_mode}")
    years = records.assign(
        document_year=pd.to_numeric(records.document_year, errors="coerce")
    ).groupby("document_chembl_id").document_year
    return years.max() if year_mode == "latest" else years.min()


def add_unit_years(
    documents: pd.DataFrame, records: pd.DataFrame, year_mode: str
) -> pd.DataFrame:
    years = unit_year_map(records, year_mode)
    output = documents.copy()
    output["year"] = [
        (
            min(
                (
                    float(years[item])
                    for item in str(unit_documents).split(";")
                    if item in years.index and np.isfinite(years[item])
                ),
                default=np.nan,
            )
            if year_mode == "earliest"
            else max(
                (
                    float(years[item])
                    for item in str(unit_documents).split(";")
                    if item in years.index and np.isfinite(years[item])
                ),
                default=np.nan,
            )
        )
        for unit_documents in output.documents
    ]
    return output


def classify(value: float) -> int:
    from chemgridmap.activity import activity_class

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


def temporal_predictions(
    target: str,
    endpoint: str,
    records: pd.DataFrame,
    cutoff: int,
    year_mode: str,
) -> pd.DataFrame:
    """Rebuild E2b temporal predictions with configurable unit-year aggregation."""
    known = document_values(records, merge_identical=True)
    known = add_unit_years(known, records, year_mode)
    known = known.loc[known.year.notna()].copy()
    past = known.loc[known.year <= cutoff]
    past_multi_by_entry = past.groupby("entry").size()
    if int((past_multi_by_entry >= 2).sum()) < 5:
        return pd.DataFrame()

    sigma_fit = fit_between_unit_sd(past)
    source_by_entry = dict(
        tuple(records.groupby("molecule_identity_key", sort=True))
    )
    rows = []
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
                past_values.mean(),
                sigma_fit["tau"] / np.sqrt(m_past),
                display,
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
                "p_displayed": float(probabilities[display]),
                "future_value": float(future.document_value),
                "minimum_past_gap": float(
                    np.min(np.abs(past_values - future.document_value))
                ),
            })
    return pd.DataFrame(rows)


def build_temporal_frame(
    records: dict[tuple[str, str], pd.DataFrame],
    inclusion: pd.DataFrame,
    year_mode: str,
) -> pd.DataFrame:
    if year_mode in _TEMPORAL_FRAME_CACHE:
        return _TEMPORAL_FRAME_CACHE[year_mode]
    frames = []
    for target, endpoint in temporal_keys(inclusion, records):
        for cutoff in TIME_CUTOFFS:
            frame = temporal_predictions(
                target,
                endpoint,
                records[(target, endpoint)],
                cutoff,
                year_mode,
            )
            if not frame.empty:
                frames.append(frame)
    output = (
        pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    )
    _TEMPORAL_FRAME_CACHE[year_mode] = output
    return output


def base_row(section: int, variant: str, record_type: str, **values) -> dict:
    return {
        "analysis_role": "post hoc",
        "analysis_section": section,
        "variant": variant,
        "record_type": record_type,
        **values,
    }


def dataset_comparisons(
    frame: pd.DataFrame, variant: str
) -> list[dict]:
    """Build eligible dataset contrasts and one paired test per scope/cutoff."""
    rows = []
    for cutoff in TIME_CUTOFFS:
        cutoff_frame = frame.loc[frame.cutoff_year == cutoff]
        for endpoint in (*ENDPOINTS, "ALL"):
            selected = (
                cutoff_frame
                if endpoint == "ALL"
                else cutoff_frame.loc[cutoff_frame.endpoint == endpoint]
            )
            eligible = []
            eligible_predictions = 0
            for (target, dataset_endpoint), dataset in selected.groupby(
                ["target", "endpoint"], sort=True
            ):
                supported = dataset.loc[
                    dataset.display_status == "supported"
                ]
                uncertain = dataset.loc[
                    dataset.display_status == "uncertain"
                ]
                if len(supported) < 5 or len(uncertain) < 5:
                    continue
                supported_rate = float(supported.agree.mean())
                uncertain_rate = float(uncertain.agree.mean())
                difference = supported_rate - uncertain_rate
                eligible.append(difference)
                eligible_predictions += len(dataset)
                rows.append(base_row(
                    1,
                    variant,
                    "dataset_comparison",
                    cutoff_year=cutoff,
                    endpoint=endpoint,
                    target=target,
                    n_predictions=len(dataset),
                    supported_n=len(supported),
                    supported_agreement_rate=supported_rate,
                    uncertain_n=len(uncertain),
                    uncertain_agreement_rate=uncertain_rate,
                    dataset_difference=difference,
                ))
            statistic, p_value, n_datasets = one_sided_wilcoxon(
                eligible, alternative="greater"
            )
            rows.append(base_row(
                1,
                variant,
                "paired_summary",
                cutoff_year=cutoff,
                endpoint=endpoint,
                n_predictions=eligible_predictions,
                n_datasets=n_datasets,
                median_dataset_difference=(
                    float(np.median(eligible)) if eligible else np.nan
                ),
                positive_difference_datasets=int(
                    np.sum(np.asarray(eligible) > 0)
                ),
                wilcoxon_statistic=statistic,
                wilcoxon_p_one_sided=p_value,
            ))
    return rows


def metric_row(
    section: int,
    variant: str,
    cutoff: int,
    endpoint: str,
    sample_filter: str,
    selected: pd.DataFrame,
    display_status: str | None = None,
    m_past: str | None = None,
) -> dict:
    metrics = time_split_metrics(selected)
    supported = selected.loc[selected.display_status == "supported"]
    uncertain = selected.loc[selected.display_status == "uncertain"]
    return base_row(
        section,
        variant,
        "m_past_stratum" if m_past else "pooled_summary",
        cutoff_year=cutoff,
        endpoint=endpoint,
        sample_filter=sample_filter,
        display_status=display_status,
        m_past=m_past,
        n_predictions=metrics["n_predictions"],
        n_entries=metrics["n_entries"],
        n_targets=metrics["n_targets"],
        supported_n=len(supported),
        supported_agreement_rate=(
            float(supported.agree.mean()) if len(supported) else np.nan
        ),
        uncertain_n=len(uncertain),
        uncertain_agreement_rate=(
            float(uncertain.agree.mean()) if len(uncertain) else np.nan
        ),
        mean_p_displayed=(
            float(selected.p_displayed.mean()) if len(selected) else np.nan
        ),
        calibration_bias=metrics["g_calibration_bias"],
        auc=metrics["g_auc"],
    )


def pooled_after_reference_exclusion(
    frame: pd.DataFrame, variant: str
) -> list[dict]:
    """Section 2: remove CHEMBL5113/Ki and summarize pooled results."""
    rows = []
    filtered = frame.loc[
        ~((frame.target == "CHEMBL5113") & (frame.endpoint == "Ki"))
    ]
    for cutoff in TIME_CUTOFFS:
        cutoff_frame = filtered.loc[filtered.cutoff_year == cutoff]
        for endpoint in (*ENDPOINTS, "ALL"):
            selected = (
                cutoff_frame
                if endpoint == "ALL"
                else cutoff_frame.loc[cutoff_frame.endpoint == endpoint]
            )
            for status in ("all", "supported", "uncertain"):
                status_frame = (
                    selected
                    if status == "all"
                    else selected.loc[selected.display_status == status]
                )
                rows.append(metric_row(
                    2,
                    variant,
                    cutoff,
                    endpoint,
                    "exclude_CHEMBL5113_Ki",
                    status_frame,
                    display_status=status,
                ))
    return rows


def auc_distribution(frame: pd.DataFrame, variant: str) -> list[dict]:
    """Section 3: dataset AUC distribution together with pooled AUC."""
    rows = []
    for cutoff in TIME_CUTOFFS:
        cutoff_frame = frame.loc[frame.cutoff_year == cutoff]
        for endpoint in (*ENDPOINTS, "ALL"):
            selected = (
                cutoff_frame
                if endpoint == "ALL"
                else cutoff_frame.loc[cutoff_frame.endpoint == endpoint]
            )
            for label, sample_filter, sample in (
                (
                    "all datasets",
                    "all_datasets",
                    selected,
                ),
                (
                    "excluding CHEMBL5113/Ki",
                    "exclude_CHEMBL5113_Ki",
                    selected.loc[
                        ~(
                            (selected.target == "CHEMBL5113")
                            & (selected.endpoint == "Ki")
                        )
                    ],
                ),
            ):
                del label
                dataset_auc = []
                for _, dataset in sample.groupby(
                    ["target", "endpoint"], sort=True
                ):
                    value = float(time_split_metrics(dataset)["g_auc"])
                    if np.isfinite(value):
                        dataset_auc.append(value)
                if dataset_auc:
                    q1, median, q3 = np.percentile(
                        np.asarray(dataset_auc), [25, 50, 75]
                    )
                else:
                    q1 = median = q3 = np.nan
                pooled = time_split_metrics(sample)
                rows.append(base_row(
                    3,
                    variant,
                    "auc_summary",
                    cutoff_year=cutoff,
                    endpoint=endpoint,
                    sample_filter=sample_filter,
                    n_predictions=pooled["n_predictions"],
                    n_datasets=len(dataset_auc),
                    calibration_bias=pooled["g_calibration_bias"],
                    auc=pooled["g_auc"],
                    dataset_auc_median=float(median),
                    dataset_auc_q1=float(q1),
                    dataset_auc_q3=float(q3),
                    dataset_auc_iqr=float(q3 - q1),
                ))
    return rows


def m_past_strata(frame: pd.DataFrame, variant: str) -> list[dict]:
    """Section 4: display reliability by one versus multiple past units."""
    rows = []
    for cutoff in TIME_CUTOFFS:
        cutoff_frame = frame.loc[frame.cutoff_year == cutoff]
        for endpoint in (*ENDPOINTS, "ALL"):
            selected = (
                cutoff_frame
                if endpoint == "ALL"
                else cutoff_frame.loc[cutoff_frame.endpoint == endpoint]
            )
            for m_past in ("m=1", "m>=2"):
                stratum = selected.loc[selected.stratum == m_past]
                for status in ("supported", "uncertain"):
                    status_frame = stratum.loc[
                        stratum.display_status == status
                    ]
                    rows.append(metric_row(
                        4,
                        variant,
                        cutoff,
                        endpoint,
                        "all_datasets",
                        status_frame,
                        display_status=status,
                        m_past=m_past,
                    ))
    return rows


def build_posthoc() -> pd.DataFrame:
    global _POSTHOC_CACHE
    if _POSTHOC_CACHE is not None:
        return _POSTHOC_CACHE.copy(deep=True)
    inclusion = pd.read_csv(INCLUSION, float_precision="round_trip")
    records = load_records()
    baseline = pd.read_csv(BASELINE, float_precision="round_trip")
    latest_rebuilt = build_temporal_frame(records, inclusion, "latest")

    rows = []
    rows.extend(dataset_comparisons(baseline, "latest_unit_year"))
    rows.extend(pooled_after_reference_exclusion(
        baseline, "latest_unit_year"
    ))
    rows.extend(auc_distribution(baseline, "latest_unit_year"))
    rows.extend(m_past_strata(baseline, "latest_unit_year"))

    earliest = build_temporal_frame(records, inclusion, "earliest")
    rows.extend(dataset_comparisons(earliest, "earliest_unit_year"))
    rows.extend(pooled_after_reference_exclusion(
        earliest, "earliest_unit_year"
    ))

    deduplicated = baseline.merge(
        latest_rebuilt[
            [
                "target",
                "endpoint",
                "cutoff_year",
                "entry",
                "future_unit",
                "minimum_past_gap",
            ]
        ],
        on=["target", "endpoint", "cutoff_year", "entry", "future_unit"],
        validate="one_to_one",
    )
    deduplicated = deduplicated.loc[
        deduplicated.minimum_past_gap > NEAR_DUPLICATE_LOG10_GAP
    ].drop(columns="minimum_past_gap")
    rows.extend(dataset_comparisons(
        deduplicated, "drop_near_duplicate_future_units"
    ))
    rows.extend(pooled_after_reference_exclusion(
        deduplicated, "drop_near_duplicate_future_units"
    ))

    output = pd.DataFrame(rows, columns=COLUMNS)
    output = output.sort_values(
        [
            "analysis_section",
            "variant",
            "cutoff_year",
            "endpoint",
            "record_type",
            "target",
            "sample_filter",
            "m_past",
            "display_status",
        ],
        na_position="first",
        kind="stable",
    ).reset_index(drop=True)
    _POSTHOC_CACHE = output
    return output.copy(deep=True)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    output = build_posthoc()
    output.to_csv(OUTPUT, index=False)
    print(f"Wrote {len(output)} rows to {OUTPUT}")


if __name__ == "__main__":
    main()
