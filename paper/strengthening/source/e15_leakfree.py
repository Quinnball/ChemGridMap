"""E15 leak-free display-rule evaluation (post hoc, sealed plan)."""
from __future__ import annotations

from datetime import datetime
import gzip
import hashlib
import json
import os
from pathlib import Path
import pickle
import sys
import time

sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parents[3]
S = ROOT / "paper/strengthening"
SOURCE = S / "source"
OUT = S / "output/E15"
QA = S / "qa/E15"
E2 = S / "output/E2"
E5 = S / "output/E5"
E5C = S / "output/E5c"
E9 = S / "output/E9"
E11 = S / "output/E11"
E12A = S / "output/E12a"
VALIDATION_DATA = S / "data/chembl37_api_validation"

sys.path[:0] = [str(SOURCE), str(ROOT / "paper"), str(ROOT / "src")]

import numpy as np
import pandas as pd
from scipy import stats

from chemgridmap.activity import activity_class
from chemgridmap.evidence import CLASSES, class_probabilities, document_values
from e1a_statistics import auc_mann_whitney
from e2_models import FOLD_SEED, curate_activity_table, document_pair_observations
from e5_analysis import score_frame
from e5_models import (
    CLASS_BOUNDARIES,
    SUPPORT_THRESHOLD,
    document_fold_map,
    pooled_sigma,
)
from e5c_models import (
    delta_document_map,
    pair_observations_from_units,
    unit_deltas,
)
from e7_analysis import holm_adjust, two_stage_weights
from e9_analysis import bootstrap_cell, random_effect_details
from manifest_utils import write_manifest
from reliability_activity import entry_folds, gaussian_prediction

ANALYSIS_ROLE = "post hoc"
ENDPOINTS = ("IC50", "Ki")
BOOT_REPETITIONS = 2000
BOOT_SEED = 20261105
OREXIN_TARGETS = {"CHEMBL5113", "CHEMBL4792"}
DESIGNS = (
    "entry",
    "strict_document",
    "near_duplicate",
    "time_2015",
    "E11",
)
DESIGN_LABELS = {
    "entry": "条目分折",
    "strict_document": "严格文献分折",
    "near_duplicate": "近重复合并",
    "time_2015": "时间切分 2015",
    "E11": "E11 新靶点",
}

STARTED_AT = datetime.now().astimezone().isoformat(timespec="seconds")
STARTED_PERF = time.perf_counter()


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
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(json_safe(value), indent=2, ensure_ascii=False, allow_nan=False)
        + "\n",
        encoding="utf-8",
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def split_documents(value) -> set[str]:
    return {part for part in str(value).split(";") if part}


def pair_key(document_a: str, document_b: str) -> tuple[str, str]:
    return tuple(sorted((str(document_a), str(document_b))))


def strict_fold_masks(
    units: pd.DataFrame,
    fold: int,
) -> tuple[np.ndarray, np.ndarray, set[str]]:
    """Return training/held masks and every document named in the held fold."""
    held_mask = units.document_fold.to_numpy(int) == int(fold)
    held_documents = {
        document
        for value in units.loc[held_mask, "documents"]
        for document in split_documents(value)
    }
    touches_held = units.documents.map(
        lambda value: bool(split_documents(value) & held_documents)
    ).to_numpy(bool)
    return ~touches_held, held_mask, held_documents


def load_e2_records() -> dict[tuple[str, str], pd.DataFrame]:
    path = S / "qa/E2/curated_records.pkl.gz"
    records = pickle.loads(gzip.decompress(path.read_bytes()))
    inclusion = pd.read_csv(E2 / "target_inclusion.csv")
    keys = {
        (str(row.target), str(row.endpoint))
        for row in inclusion.loc[inclusion.included].itertuples(index=False)
    }
    return {key: records[key] for key in keys}


def archived_rule_frame() -> pd.DataFrame:
    frame = pd.read_csv(
        E5C / "predictions_e5c.csv",
        low_memory=False,
        float_precision="round_trip",
    )
    return frame.loc[frame.method == "G"].copy()


def archived_e11_frame() -> pd.DataFrame:
    frame = pd.read_csv(
        E11 / "predictions_e11.csv",
        low_memory=False,
        float_precision="round_trip",
    )
    return frame.loc[frame.method == "G"].copy()


def training_entry_folds(units: pd.DataFrame) -> dict[str, int]:
    original = units.loc[:, ["entry", "document", "document_value", "documents"]]
    return entry_folds(original, seed=FOLD_SEED)


def filtered_observations(
    observations: pd.DataFrame,
    *,
    training_entries: set[str] | None = None,
    held_documents: set[str] | None = None,
    excluded_pairs: set[tuple[str, str]] | None = None,
) -> pd.DataFrame:
    selected = observations
    if training_entries is not None:
        selected = selected.loc[selected.entry.astype(str).isin(training_entries)]
    if held_documents:
        keep = np.fromiter(
            (
                not (
                    split_documents(document_a)
                    | split_documents(document_b)
                )
                & held_documents
                for document_a, document_b in zip(
                    selected.document_a.to_numpy(object),
                    selected.document_b.to_numpy(object),
                )
            ),
            dtype=bool,
            count=len(selected),
        )
        selected = selected.loc[keep]
    if excluded_pairs:
        keys = pd.Series(
            [
                pair_key(document_a, document_b)
                for document_a, document_b in zip(
                    selected.document_a.to_numpy(object),
                    selected.document_b.to_numpy(object),
                )
            ],
            index=selected.index,
            dtype=object,
        )
        selected = selected.loc[~keys.isin(excluded_pairs)]
    return selected


def leakfree_delta_for_row(
    units: pd.DataFrame,
    observations: pd.DataFrame,
    row: pd.Series,
    *,
    replace_withheld_documents: bool,
    ignore_held_document_pairs: bool,
    excluded_pairs: set[tuple[str, str]] | None = None,
) -> dict[str, float]:
    held_documents = split_documents(row.held_out_unit)
    if replace_withheld_documents:
        keep = units.documents.map(
            lambda value: not (split_documents(value) & held_documents)
        )
        training = units.loc[keep]
    else:
        training = units.loc[
            units.entry_fold != int(row.fold)
        ]
    selected = filtered_observations(
        observations,
        training_entries=set(training.entry.astype(str)),
        held_documents=held_documents if ignore_held_document_pairs else None,
        excluded_pairs=excluded_pairs,
    )
    return delta_document_map(selected)


def remaining_units_for_row(
    group: pd.DataFrame,
    row: pd.Series,
    *,
    removed_documents: set[str] | None = None,
) -> pd.DataFrame:
    documents = split_documents(row.held_out_unit)
    if removed_documents:
        documents |= removed_documents
    keep = group.documents.map(
        lambda value: not (split_documents(value) & documents)
    )
    return group.loc[keep]


def recompute_archived_flag_touched(
    records_by_key: dict[tuple[str, str], pd.DataFrame],
    predictions: pd.DataFrame,
) -> pd.DataFrame:
    """R1: independently rebuild the archived training-fold touched marker."""
    rows = []
    for (target, endpoint), dataset in predictions.groupby(
        ["target", "endpoint"], sort=True
    ):
        records = records_by_key[(target, endpoint)]
        original = document_values(records, merge_identical=True)
        if dataset.split_design.iloc[0] == "entry":
            units = original.copy()
            units["entry_fold"] = (
                units.entry.map(training_entry_folds(units)).fillna(-1).astype(int)
            )
            observations = document_pair_observations(records)
            groups = dict(tuple(units.groupby("entry", sort=True)))
            fold_delta = {}
            for fold in sorted(dataset.fold.unique()):
                training = units.loc[units.entry_fold != int(fold)]
                selected = filtered_observations(
                    observations,
                    training_entries=set(training.entry.astype(str)),
                )
                fold_delta[int(fold)] = delta_document_map(selected)
            for item in dataset.to_dict("records"):
                row = pd.Series(item)
                delta = fold_delta[int(row.fold)]
                remaining = remaining_units_for_row(
                    groups[str(row.entry)], row
                )
                _, touched = unit_deltas(remaining, delta)
                rows.append(
                    {
                        "target": target,
                        "endpoint": endpoint,
                        "fold": int(row.fold),
                        "entry": str(row.entry),
                        "held_out_unit": str(row.held_out_unit),
                        "observed_flag_touched": bool(row.flag_touched),
                        "recomputed_flag_touched": bool(touched.any()),
                    }
                )
        elif dataset.split_design.iloc[0] == "document":
            units = original.copy()
            units["document_fold"] = document_fold_map(units)
            groups = dict(tuple(units.groupby("entry", sort=True)))
            observations = document_pair_observations(records)
            for fold, fold_predictions in dataset.groupby("fold", sort=True):
                train = units.loc[units.document_fold != int(fold)]
                selected = filtered_observations(
                    observations,
                    training_entries=set(train.entry.astype(str)),
                )
                delta = delta_document_map(selected)
                for item in fold_predictions.to_dict("records"):
                    row = pd.Series(item)
                    remaining = remaining_units_for_row(
                        groups[str(row.entry)], row
                    )
                    remaining = remaining.loc[
                        remaining.document_fold != int(row.fold)
                    ]
                    _, touched = unit_deltas(remaining, delta)
                    rows.append(
                        {
                            "target": target,
                            "endpoint": endpoint,
                            "fold": int(row.fold),
                            "entry": str(row.entry),
                            "held_out_unit": str(row.held_out_unit),
                            "observed_flag_touched": bool(row.flag_touched),
                            "recomputed_flag_touched": bool(touched.any()),
                        }
                    )
        else:
            raise ValueError(
                f"R1 received unexpected archived design: {dataset.split_design.iloc[0]}"
            )
    result = pd.DataFrame(rows)
    result["matches"] = (
        result.observed_flag_touched == result.recomputed_flag_touched
    )
    return result


def delta_document_map_from_units_local(units: pd.DataFrame) -> dict[str, float]:
    observations = pair_observations_from_units(units)
    return delta_document_map(observations, expand_documents=True)


def apply_entry_like_leakfree_flags(
    records_by_key: dict[tuple[str, str], pd.DataFrame],
    predictions: pd.DataFrame,
    *,
    source_label: str,
    excluded_pairs_by_dataset: dict[
        tuple[str, str], set[tuple[str, str]]
    ] | None = None,
) -> pd.DataFrame:
    excluded_pairs_by_dataset = excluded_pairs_by_dataset or {}
    outputs = []
    for (target, endpoint), dataset in predictions.groupby(
        ["target", "endpoint"], sort=True
    ):
        records = records_by_key[(target, endpoint)]
        original = document_values(records, merge_identical=True)
        units = original.copy()
        units["entry_fold"] = (
            units.entry.map(training_entry_folds(units)).fillna(-1).astype(int)
        )
        observations = document_pair_observations(records)
        groups = dict(tuple(units.groupby("entry", sort=True)))
        excluded = excluded_pairs_by_dataset.get((target, endpoint), set())
        fold_observations = {}
        delta_cache = {}
        for fold in sorted(dataset.fold.unique()):
            training = units.loc[units.entry_fold != int(fold)]
            selected = filtered_observations(
                observations,
                training_entries=set(training.entry.astype(str)),
                excluded_pairs=excluded,
            )
            fold_observations[int(fold)] = selected
        rows = []
        for item in dataset.to_dict("records"):
            row = pd.Series(item)
            cache_key = (int(row.fold), str(row.held_out_unit))
            if cache_key not in delta_cache:
                delta_cache[cache_key] = delta_document_map(
                    filtered_observations(
                        fold_observations[int(row.fold)],
                        held_documents=split_documents(row.held_out_unit),
                    )
                )
            delta = delta_cache[cache_key]
            remaining = remaining_units_for_row(
                groups[str(row.entry)], row
            )
            _, touched = unit_deltas(remaining, delta)
            rows.append(
                {
                    **item,
                    "design": source_label,
                    "leakfree_flag_touched": bool(touched.any()),
                    "n_leakfree_remaining_flagged_units": int(touched.sum()),
                    "pair_exclusion_rule": (
                        "training entries; pairs containing the held-out "
                        "document are ignored"
                    ),
                }
            )
        outputs.append(pd.DataFrame(rows))
    return pd.concat(outputs, ignore_index=True)


def apply_near_duplicate_leakfree_flags(
    records_by_key: dict[tuple[str, str], pd.DataFrame],
    predictions: pd.DataFrame,
) -> pd.DataFrame:
    from e4c_analysis import near_duplicate_units

    outputs = []
    for (target, endpoint), dataset in predictions.groupby(
        ["target", "endpoint"], sort=True
    ):
        records = records_by_key[(target, endpoint)]
        original = document_values(records, merge_identical=True)
        assignments = entry_folds(original, seed=FOLD_SEED)
        units = near_duplicate_units(records)
        units["entry_fold"] = (
            units.entry.map(assignments).fillna(-1).astype(int)
        )
        observations = pair_observations_from_units(units)
        groups = dict(tuple(units.groupby("entry", sort=True)))
        fold_observations = {}
        delta_cache = {}
        for fold in sorted(dataset.fold.unique()):
            training = units.loc[units.entry_fold != int(fold)]
            selected = filtered_observations(
                observations,
                training_entries=set(training.entry.astype(str)),
            )
            fold_observations[int(fold)] = selected
        rows = []
        for item in dataset.to_dict("records"):
            row = pd.Series(item)
            cache_key = (int(row.fold), str(row.held_out_unit))
            if cache_key not in delta_cache:
                delta_cache[cache_key] = delta_document_map(
                    filtered_observations(
                        fold_observations[int(row.fold)],
                        held_documents=split_documents(row.held_out_unit),
                    ),
                    expand_documents=True,
                )
            delta = delta_cache[cache_key]
            group = groups[str(row.entry)]
            remaining = remaining_units_for_row(group, row)
            _, touched = unit_deltas(remaining, delta)
            rows.append(
                {
                    **item,
                    "design": "near_duplicate",
                    "leakfree_flag_touched": bool(touched.any()),
                    "n_leakfree_remaining_flagged_units": int(touched.sum()),
                    "pair_exclusion_rule": (
                        "training entries; pairs containing the held-out "
                        "document are ignored"
                    ),
                }
            )
        outputs.append(pd.DataFrame(rows))
    return pd.concat(outputs, ignore_index=True)


def document_fold_predictions_e15(
    records: pd.DataFrame,
    target: str,
    endpoint: str,
    *,
    strict_leak_free: bool,
) -> pd.DataFrame:
    """Return model G in document folds; off reproduces E5 exactly."""
    from e5_models import classify

    units = document_values(records, merge_identical=True)
    units["document_fold"] = document_fold_map(units)
    groups = dict(tuple(units.groupby("entry", sort=True)))
    source_by_entry = dict(
        tuple(records.groupby("molecule_identity_key", sort=True))
    )
    rows = []
    for fold in range(10):
        if strict_leak_free:
            keep_training, held_mask, held_fold_documents = strict_fold_masks(
                units, fold
            )
            train = units.loc[keep_training].copy()
        else:
            held_mask = units.document_fold.to_numpy(int) == fold
            held_fold_documents = {
                document
                for value in units.loc[held_mask, "documents"]
                for document in split_documents(value)
            }
            train = units.loc[~held_mask].copy()
        sigma_fit = pooled_sigma(train)
        sigma = float(sigma_fit["tau"])
        delta_map = delta_document_map_from_units_local(train)
        test = units.loc[held_mask]
        for held in test.itertuples(index=False):
            group = groups[str(held.entry)]
            if strict_leak_free:
                remaining = group.loc[
                    group.documents.map(
                        lambda value: not (
                            split_documents(value) & held_fold_documents
                        )
                    )
                ]
            else:
                held_documents = split_documents(held.documents)
                remaining = group.loc[
                    group.documents.map(
                        lambda value: not (
                            split_documents(value) & held_documents
                        )
                    )
                ].loc[lambda frame: frame.document_fold != fold]
            if remaining.empty:
                continue
            source = source_by_entry[str(held.entry)]
            documents = set(
                ";".join(remaining.documents.astype(str)).split(";")
            )
            documents.discard("")
            raw = source.loc[
                source.document_chembl_id.astype(str).isin(documents)
            ]
            if raw.empty:
                continue
            values = remaining.document_value.to_numpy(float)
            displayed = classify(float(np.median(raw.pchembl_value)))
            observed = classify(float(held.document_value))
            predictive = gaussian_prediction(values, sigma_fit)
            support = class_probabilities(
                [values.mean()],
                [sigma / np.sqrt(len(values))],
                np.inf,
                *CLASS_BOUNDARIES,
            )[0]
            _, touched = unit_deltas(remaining, delta_map)
            shown = int(displayed)
            rows.append(
                {
                    "split_design": "document",
                    "merge_basis": "same_value",
                    "target": target,
                    "endpoint": endpoint,
                    "fold": fold,
                    "entry": str(held.entry),
                    "held_out_unit": str(held.documents),
                    "m": int(len(values)),
                    "displayed": shown,
                    "observed": observed,
                    "p0": float(predictive[0]),
                    "p1": float(predictive[1]),
                    "p2": float(predictive[2]),
                    "p_displayed": float(predictive[shown]),
                    "support": float(support[shown]),
                    "agree": bool(observed == shown),
                    "flag_touched": bool(touched.any()),
                    "sigma_doc": sigma,
                    "design": (
                        "strict_document"
                        if strict_leak_free
                        else "document"
                    ),
                    "leakfree_flag_touched": bool(touched.any()),
                    "n_leakfree_remaining_flagged_units": int(touched.sum()),
                    "pair_exclusion_rule": (
                        "fitted statistics and pairs use training units only; "
                        "units containing any held-fold document are removed"
                    ),
                }
            )
    return pd.DataFrame(rows).sort_values(
        ["endpoint", "fold", "entry", "held_out_unit"]
    ).reset_index(drop=True)


def archived_document_g() -> pd.DataFrame:
    frame = pd.read_csv(
        E5 / "document_folds_predictions.csv",
        low_memory=False,
        float_precision="round_trip",
    )
    return frame.loc[frame.method == "G"].reset_index(drop=True)


def curate_e11(records_by_key: dict) -> tuple[dict, pd.DataFrame]:
    inclusion = pd.read_csv(E11 / "target_inclusion.csv")
    selected = pd.read_csv(VALIDATION_DATA / "selected_targets.csv")
    curated = {}
    r5_rows = []
    for item in selected.itertuples(index=False):
        target = str(item.target_chembl_id)
        endpoint = str(item.standard_type)
        raw_path = VALIDATION_DATA / f"activities/{target}_{endpoint}.csv"
        raw = pd.read_csv(raw_path, low_memory=False)
        expected = inclusion.loc[
            (inclusion.target == target)
            & (inclusion.endpoint == endpoint)
        ]
        expected_n = (
            int(expected.retained_records.iloc[0])
            if len(expected) and np.isfinite(expected.retained_records.iloc[0])
            else np.nan
        )
        try:
            # E15 R5 checks the E11 execution under its sealed pre-0.5 rule.
            result = curate_activity_table(
                raw,
                target,
                endpoint,
                identity_rule="normalized-structure",
            )
            retained_records = result.retained_records
            records_by_key[(target, endpoint)] = retained_records
            curated[(target, endpoint)] = retained_records
            recomputed_n = int(len(retained_records))
            failure = ""
        except Exception as error:
            retained_records = None
            recomputed_n = np.nan
            failure = f"{type(error).__name__}: {error}"
            if np.isfinite(expected_n):
                raise
        r5_rows.append(
            {
                "target": target,
                "endpoint": endpoint,
                "recomputed_retained_records": recomputed_n,
                "archived_retained_records": expected_n,
                "matches": bool(
                    (not np.isfinite(expected_n) and not np.isfinite(recomputed_n))
                    or (
                        np.isfinite(expected_n)
                        and recomputed_n == expected_n
                    )
                ),
                "curation_error": failure,
            }
        )
    r5 = pd.DataFrame(r5_rows)
    if not r5.matches.all():
        raise RuntimeError("R5 E11 curation counts do not match the archive")
    return curated, r5


def full_data_flagged_pairs(
    records: pd.DataFrame,
) -> set[tuple[str, str]]:
    observations = document_pair_observations(records)
    flagged = set()
    for (a, b), group in observations.groupby(
        ["document_a", "document_b"], sort=True
    ):
        if len(group) < 3:
            continue
        median = float(np.median(group.difference))
        sign = float(np.sign(median))
        agreement = float(np.mean(group.difference * sign > 0))
        if abs(median) >= 1.0 and agreement >= 0.8:
            flagged.add(pair_key(a, b))
    return flagged


def build_design_frame(
    archived: pd.DataFrame,
    records_by_key: dict,
    *,
    design: str,
    excluded_pairs_by_dataset=None,
) -> pd.DataFrame:
    if design in {"entry", "E11"}:
        base = archived.loc[
            (archived.split_design == "entry")
            & (archived.merge_basis == "same_value")
        ].reset_index(drop=True)
        result = apply_entry_like_leakfree_flags(
            records_by_key,
            base,
            source_label=design,
            excluded_pairs_by_dataset=excluded_pairs_by_dataset,
        )
    elif design == "near_duplicate":
        base = archived.loc[
            (archived.split_design == "entry")
            & (archived.merge_basis == "near_duplicate_0.1")
        ].reset_index(drop=True)
        result = apply_near_duplicate_leakfree_flags(
            records_by_key, base
        )
    elif design == "time_2015":
        result = archived.loc[
            (archived.split_design == "time")
            & (archived.merge_basis == "same_value")
            & (archived.cutoff_year == 2015)
        ].copy().reset_index(drop=True)
        result["design"] = design
        result["leakfree_flag_touched"] = result.flag_touched.astype(bool)
        result["n_leakfree_remaining_flagged_units"] = result[
            "n_remaining_flagged_units"
        ]
        result["pair_exclusion_rule"] = (
            "only documents at or before the 2015 cutoff are used"
        )
    else:
        raise ValueError(design)
    result["g_supported"] = result.support.to_numpy(float) >= SUPPORT_THRESHOLD
    result["p_supported"] = (
        result.g_supported & ~result.leakfree_flag_touched.astype(bool)
    )
    result["demoted"] = (
        result.g_supported & result.leakfree_flag_touched.astype(bool)
    )
    result["uncertain"] = ~result.g_supported
    return result


def build_strict_document_frame(records_by_key: dict) -> pd.DataFrame:
    frames = []
    for number, ((target, endpoint), records) in enumerate(
        sorted(records_by_key.items()), start=1
    ):
        frames.append(
            document_fold_predictions_e15(
                records,
                target,
                endpoint,
                strict_leak_free=True,
            )
        )
        if number % 10 == 0:
            print(f"strict document folds: {number}/{len(records_by_key)}", flush=True)
    result = pd.concat(frames, ignore_index=True)
    result["g_supported"] = result.support.to_numpy(float) >= SUPPORT_THRESHOLD
    result["p_supported"] = (
        result.g_supported & ~result.leakfree_flag_touched.astype(bool)
    )
    result["demoted"] = (
        result.g_supported & result.leakfree_flag_touched.astype(bool)
    )
    result["uncertain"] = ~result.g_supported
    return result


def wilcoxon_details(values) -> dict:
    """Plan rule: exact for n<=25 without ties, otherwise normal approximation."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    nonzero = np.abs(values[values != 0])
    ties = len(np.unique(nonzero)) < len(nonzero)
    n = len(values)
    if not n or np.allclose(values, 0):
        return {
            "n_units": n,
            "median_difference": float(np.median(values)) if n else np.nan,
            "positive_units": int((values > 0).sum()),
            "zero_units": int((values == 0).sum()),
            "negative_units": int((values < 0).sum()),
            "wilcoxon_statistic": np.nan,
            "two_sided_p": 1.0,
            "wilcoxon_method": "not_applicable_all_zero" if n else "not_applicable_empty",
            "method_reason": "no non-zero paired differences",
        }
    method = "exact" if n <= 25 and not ties else "approx"
    result = stats.wilcoxon(
        values,
        alternative="two-sided",
        zero_method="wilcox",
        method=method,
    )
    return {
        "n_units": n,
        "median_difference": float(np.median(values)),
        "positive_units": int((values > 0).sum()),
        "zero_units": int((values == 0).sum()),
        "negative_units": int((values < 0).sum()),
        "wilcoxon_statistic": float(result.statistic),
        "two_sided_p": float(result.pvalue),
        "wilcoxon_method": method,
        "method_reason": (
            "n<=25 and no tied non-zero absolute differences"
            if method == "exact"
            else "n>25 or tied non-zero absolute differences"
        ),
    }


def frame_weights(frame: pd.DataFrame) -> np.ndarray:
    return two_stage_weights(
        frame,
        seed=BOOT_SEED,
        repetitions=BOOT_REPETITIONS,
    )


def cell_metrics(frame: pd.DataFrame, mask: np.ndarray) -> dict:
    mask = np.asarray(mask, dtype=bool)
    selected = frame.loc[mask]
    weights = frame_weights(frame)
    boot = bootstrap_cell(frame, weights, mask)
    return {
        "n_predictions": int(mask.sum()),
        "n_datasets": int(selected.target.nunique()) if len(selected) else 0,
        "agreement_rate": float(selected.agree.mean()) if len(selected) else np.nan,
        "agreement_ci_low": boot["agreement_ci_low"],
        "agreement_ci_high": boot["agreement_ci_high"],
        "mean_probability": (
            float(selected.p_displayed.mean()) if len(selected) else np.nan
        ),
        "mean_support": (
            float(selected.support.mean()) if len(selected) else np.nan
        ),
        "bias": (
            float((selected.p_displayed - selected.agree).mean())
            if len(selected)
            else np.nan
        ),
        "bias_ci_low": boot["bias_ci_low"],
        "bias_ci_high": boot["bias_ci_high"],
        "valid_repetitions": boot["valid_repetitions"],
    }


def dataset_rise_values(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for target, group in frame.groupby("target", sort=True):
        demoted = group.demoted.to_numpy(bool)
        g_supported = group.g_supported.to_numpy(bool)
        p_supported = group.p_supported.to_numpy(bool)
        eligible = bool(demoted.any() and g_supported.any() and p_supported.any())
        rise = (
            float(
                group.loc[p_supported, "agree"].mean()
                - group.loc[g_supported, "agree"].mean()
            )
            if eligible
            else np.nan
        )
        rows.append(
            {
                "target": target,
                "eligible": eligible,
                "rise": rise,
                "g_supported_n": int(g_supported.sum()),
                "p_supported_n": int(p_supported.sum()),
                "demoted_n": int(demoted.sum()),
            }
        )
    return pd.DataFrame(rows)


def q1(frames: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    tests = []
    for design in DESIGNS:
        frame = frames[design]
        for endpoint in ENDPOINTS:
            selected = frame.loc[frame.endpoint == endpoint].reset_index(drop=True)
            base = cell_metrics(selected, selected.g_supported.to_numpy(bool))
            rule = cell_metrics(selected, selected.p_supported.to_numpy(bool))
            demoted = cell_metrics(selected, selected.demoted.to_numpy(bool))
            rows.append(
                {
                    "analysis_role": ANALYSIS_ROLE,
                    "design": design,
                    "design_label": DESIGN_LABELS[design],
                    "endpoint": endpoint,
                    "g_supported_n": base["n_predictions"],
                    "g_supported_n_datasets": base["n_datasets"],
                    "g_supported_agreement": base["agreement_rate"],
                    "g_supported_agreement_ci_low": base["agreement_ci_low"],
                    "g_supported_agreement_ci_high": base["agreement_ci_high"],
                    "p_supported_n": rule["n_predictions"],
                    "p_supported_n_datasets": rule["n_datasets"],
                    "p_supported_agreement": rule["agreement_rate"],
                    "p_supported_agreement_ci_low": rule["agreement_ci_low"],
                    "p_supported_agreement_ci_high": rule["agreement_ci_high"],
                    "supported_rise": (
                        rule["agreement_rate"] - base["agreement_rate"]
                    ),
                    "demoted_n": demoted["n_predictions"],
                    "demoted_n_datasets": demoted["n_datasets"],
                    "demoted_rate": (
                        demoted["n_predictions"] / base["n_predictions"]
                        if base["n_predictions"]
                        else np.nan
                    ),
                    "demoted_agreement": demoted["agreement_rate"],
                    "demoted_agreement_ci_low": demoted["agreement_ci_low"],
                    "demoted_agreement_ci_high": demoted["agreement_ci_high"],
                    "bootstrap_method": "two_stage",
                    "bootstrap_repetitions": BOOT_REPETITIONS,
                    "bootstrap_seed": BOOT_SEED,
                }
            )
            differences = dataset_rise_values(selected)
            test = wilcoxon_details(
                differences.loc[differences.eligible, "rise"]
            )
            tests.append(
                {
                    "analysis_role": ANALYSIS_ROLE,
                    "design": design,
                    "design_label": DESIGN_LABELS[design],
                    "endpoint": endpoint,
                    "test": "rise",
                    "alternative": "two_sided",
                    "zero_method": "wilcox",
                    **test,
                }
            )
    tests = pd.DataFrame(tests)
    tests["holm_two_sided_p"] = holm_adjust(tests.two_sided_p)
    results = pd.DataFrame(rows)
    results = results.merge(
        tests[
            [
                "design",
                "endpoint",
                "n_units",
                "median_difference",
                "two_sided_p",
                "holm_two_sided_p",
                "wilcoxon_method",
            ]
        ],
        on=["design", "endpoint"],
        how="left",
    )
    return results, tests


def q2(frames: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    summaries = []
    tests = []
    comparisons = []
    for design in DESIGNS:
        frame = frames[design]
        for endpoint in ENDPOINTS:
            selected = frame.loc[frame.endpoint == endpoint].reset_index(drop=True)
            masks = {
                "St": (
                    selected.g_supported.to_numpy(bool)
                    & selected.leakfree_flag_touched.to_numpy(bool)
                ),
                "Su": (
                    selected.g_supported.to_numpy(bool)
                    & ~selected.leakfree_flag_touched.to_numpy(bool)
                ),
                "U": ~selected.g_supported.to_numpy(bool),
            }
            for state, mask in masks.items():
                summaries.append(
                    {
                        "analysis_role": ANALYSIS_ROLE,
                        "design": design,
                        "design_label": DESIGN_LABELS[design],
                        "endpoint": endpoint,
                        "state": state,
                        **cell_metrics(selected, mask),
                        "bootstrap_method": "two_stage",
                        "bootstrap_repetitions": BOOT_REPETITIONS,
                        "bootstrap_seed": BOOT_SEED,
                    }
                )
            rows = []
            for target, group in selected.groupby("target", sort=True):
                st = (
                    group.g_supported.to_numpy(bool)
                    & group.leakfree_flag_touched.to_numpy(bool)
                )
                su = (
                    group.g_supported.to_numpy(bool)
                    & ~group.leakfree_flag_touched.to_numpy(bool)
                )
                eligible = bool(st.sum() >= 5 and su.any())
                difference = (
                    float(group.loc[st, "agree"].mean() - group.loc[su, "agree"].mean())
                    if eligible
                    else np.nan
                )
                rows.append(
                    {
                        "target": target,
                        "eligible": eligible,
                        "difference": difference,
                        "St_n": int(st.sum()),
                        "Su_n": int(su.sum()),
                    }
                )
            detail = pd.DataFrame(rows)
            test = wilcoxon_details(
                detail.loc[detail.eligible, "difference"]
            )
            tests.append(
                {
                    "analysis_role": ANALYSIS_ROLE,
                    "design": design,
                    "design_label": DESIGN_LABELS[design],
                    "endpoint": endpoint,
                    "test": "St_minus_Su",
                    "alternative": "two_sided",
                    "zero_method": "wilcox",
                    **test,
                }
            )
            weights = frame_weights(selected)
            st_mask = masks["St"]
            u_mask = masks["U"]
            denominator = weights[:, st_mask].sum(axis=1) * weights[:, u_mask].sum(axis=1)
            numerator = (
                (weights[:, st_mask] * selected.agree.to_numpy(float)[st_mask]).sum(axis=1)
                * weights[:, u_mask].sum(axis=1)
                - weights[:, st_mask].sum(axis=1)
                * (weights[:, u_mask] * selected.agree.to_numpy(float)[u_mask]).sum(axis=1)
            )
            draws = np.divide(
                numerator,
                denominator,
                out=np.full(len(weights), np.nan),
                where=denominator > 0,
            )
            finite = draws[np.isfinite(draws)]
            st_rate = (
                float(selected.loc[st_mask, "agree"].mean())
                if st_mask.any()
                else np.nan
            )
            u_rate = (
                float(selected.loc[u_mask, "agree"].mean())
                if u_mask.any()
                else np.nan
            )
            comparisons.append(
                {
                    "analysis_role": ANALYSIS_ROLE,
                    "design": design,
                    "endpoint": endpoint,
                    "st_agreement": st_rate,
                    "u_agreement": u_rate,
                    "st_minus_u": st_rate - u_rate,
                    "ci_low": float(np.percentile(finite, 2.5)) if len(finite) else np.nan,
                    "ci_high": float(np.percentile(finite, 97.5)) if len(finite) else np.nan,
                    "valid_repetitions": int(len(finite)),
                    "bootstrap_method": "two_stage",
                    "bootstrap_repetitions": BOOT_REPETITIONS,
                    "bootstrap_seed": BOOT_SEED,
                }
            )
    tests = pd.DataFrame(tests)
    tests["holm_two_sided_p"] = holm_adjust(tests.two_sided_p)
    summary = pd.DataFrame(summaries).merge(
        tests[
            [
                "design",
                "endpoint",
                "n_units",
                "median_difference",
                "two_sided_p",
                "holm_two_sided_p",
                "wilcoxon_method",
            ]
        ],
        on=["design", "endpoint"],
        how="left",
    )
    return summary, tests, pd.DataFrame(comparisons)


def q4(
    frames: dict[str, pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    detail_rows = []
    summary_rows = []
    for design in ("entry", "time_2015"):
        frame = frames[design]
        for endpoint in ENDPOINTS:
            endpoint_frame = frame.loc[
                frame.endpoint == endpoint
            ].reset_index(drop=True)
            weights = frame_weights(endpoint_frame)
            for class_name in (*CLASSES, "all"):
                mask = (
                    np.ones(len(endpoint_frame), dtype=bool)
                    if class_name == "all"
                    else endpoint_frame.displayed.to_numpy(int)
                    == CLASSES.index(class_name)
                )
                selected = endpoint_frame.loc[mask]
                for target, group in selected.groupby("target", sort=True):
                    detail_rows.append(
                        {
                            "row_type": "dataset",
                            "design": design,
                            "endpoint": endpoint,
                            "display_class": class_name,
                            "target": target,
                            "n_predictions": int(len(group)),
                            "agreement_rate": float(group.agree.mean()),
                            "mean_probability": float(group.p_displayed.mean()),
                            "bias": float(
                                (group.p_displayed - group.agree).mean()
                            ),
                            "is_development_target": target
                            in {"CHEMBL204", "CHEMBL205", "CHEMBL240"},
                        }
                    )
                boot = bootstrap_cell(endpoint_frame, weights, mask)
                dataset_biases = (
                    selected.groupby("target")
                    .apply(
                        lambda group: float(
                            (group.p_displayed - group.agree).mean()
                        ),
                        include_groups=False,
                    )
                    .to_numpy(float)
                )
                detail = pd.DataFrame(
                    {
                        "bias": dataset_biases,
                        "n": selected.groupby("target").size().to_numpy(),
                    }
                )
                eligible = detail.loc[detail.n >= 10, "bias"].to_numpy(float)
                median_test = wilcoxon_details(eligible)
                re = random_effect_details(selected, 20)
                reduced = endpoint_frame.loc[
                    endpoint_frame.target != "CHEMBL5113"
                ].reset_index(drop=True)
                reduced_mask = (
                    np.ones(len(reduced), dtype=bool)
                    if class_name == "all"
                    else reduced.displayed.to_numpy(int)
                    == CLASSES.index(class_name)
                )
                reduced_weights = frame_weights(reduced)
                reduced_boot = bootstrap_cell(
                    reduced, reduced_weights, reduced_mask
                )
                reduced_selected = reduced.loc[reduced_mask]
                summary_rows.append(
                    {
                        "row_type": "summary",
                        "design": design,
                        "endpoint": endpoint,
                        "display_class": class_name,
                        "target": "ALL",
                        "n_predictions": int(mask.sum()),
                        "n_datasets": int(selected.target.nunique()),
                        "agreement_rate": (
                            float(selected.agree.mean()) if len(selected) else np.nan
                        ),
                        "agreement_ci_low": boot["agreement_ci_low"],
                        "agreement_ci_high": boot["agreement_ci_high"],
                        "mean_probability": (
                            float(selected.p_displayed.mean())
                            if len(selected)
                            else np.nan
                        ),
                        "bias": (
                            float((selected.p_displayed - selected.agree).mean())
                            if len(selected)
                            else np.nan
                        ),
                        "bias_ci_low": boot["bias_ci_low"],
                        "bias_ci_high": boot["bias_ci_high"],
                        "dataset_median_min10_n": len(eligible),
                        "dataset_bias_median_min10": median_test[
                            "median_difference"
                        ],
                        "dataset_bias_two_sided_p_min10": median_test[
                            "two_sided_p"
                        ],
                        "dataset_wilcoxon_method": median_test[
                            "wilcoxon_method"
                        ],
                        "random_effects_n_datasets": re["n_datasets"],
                        "random_effects_bias": re["random_estimate"],
                        "random_effects_ci_low": re["ci_low"],
                        "random_effects_ci_high": re["ci_high"],
                        "random_effects_tau2": re["tau2"],
                        "without_CHEMBL5113_n": int(reduced_mask.sum()),
                        "without_CHEMBL5113_bias": (
                            float(
                                (
                                    reduced_selected.p_displayed
                                    - reduced_selected.agree
                                ).mean()
                            )
                            if len(reduced_selected)
                            else np.nan
                        ),
                        "without_CHEMBL5113_bias_ci_low": reduced_boot[
                            "bias_ci_low"
                        ],
                        "without_CHEMBL5113_bias_ci_high": reduced_boot[
                            "bias_ci_high"
                        ],
                        "bootstrap_method": "two_stage",
                        "bootstrap_repetitions": BOOT_REPETITIONS,
                        "bootstrap_seed": BOOT_SEED,
                    }
                )
    detail = pd.DataFrame(detail_rows)
    summary = pd.DataFrame(summary_rows)
    combined = pd.concat([summary, detail], ignore_index=True, sort=False)
    return combined, summary


def q5(frames: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for design in ("entry", "strict_document"):
        frame = frames[design]
        for endpoint in ENDPOINTS:
            selected = frame.loc[frame.endpoint == endpoint].reset_index(drop=True)
            inactive = selected.displayed.to_numpy(int) == CLASSES.index("inactive")
            touched = selected.leakfree_flag_touched.to_numpy(bool)
            for state, mask in (
                ("touched_inactive", inactive & touched),
                ("untouched_inactive", inactive & ~touched),
                ("all_inactive", inactive),
            ):
                rows.append(
                    {
                        "analysis_role": ANALYSIS_ROLE,
                        "design": design,
                        "endpoint": endpoint,
                        "state": state,
                        **cell_metrics(selected, mask),
                        "bootstrap_method": "two_stage",
                        "bootstrap_repetitions": BOOT_REPETITIONS,
                        "bootstrap_seed": BOOT_SEED,
                    }
                )
    return pd.DataFrame(rows)


def strict_document_effect(
    old: pd.DataFrame,
    strict: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for endpoint in (*ENDPOINTS, "ALL"):
        old_frame = old if endpoint == "ALL" else old.loc[old.endpoint == endpoint]
        new_frame = strict if endpoint == "ALL" else strict.loc[strict.endpoint == endpoint]
        old_score = score_frame(old_frame)
        new_score = score_frame(new_frame)
        for metric in ("auc", "log_loss", "supported_agreement_rate"):
            old_value = float(old_score[metric])
            new_value = float(new_score[metric])
            rows.append(
                {
                    "analysis_role": ANALYSIS_ROLE,
                    "endpoint": endpoint,
                    "metric": metric,
                    "original_document_fold": old_value,
                    "strict_document_fold": new_value,
                    "difference_strict_minus_original": new_value - old_value,
                    "original_n_predictions": int(len(old_frame)),
                    "strict_n_predictions": int(len(new_frame)),
                }
            )
    return pd.DataFrame(rows)


def original_document_training_flags(
    records_by_key: dict,
) -> pd.DataFrame:
    """R3 reviewer convention: archived G, flags from training documents only."""
    archived = archived_rule_frame()
    document = archived.loc[
        (archived.split_design == "document")
        & (archived.merge_basis == "same_value")
    ].copy()
    outputs = []
    for (target, endpoint), dataset in document.groupby(
        ["target", "endpoint"], sort=True
    ):
        units = document_values(
            records_by_key[(target, endpoint)], merge_identical=True
        )
        units["document_fold"] = document_fold_map(units)
        groups = dict(tuple(units.groupby("entry", sort=True)))
        all_observations = pair_observations_from_units(units)
        fold_documents = {
            fold: {
                document
                for value in units.loc[
                    units.document_fold == fold, "documents"
                ]
                for document in split_documents(value)
            }
            for fold in range(10)
        }
        for fold, fold_frame in dataset.groupby("fold", sort=True):
            selected_pairs = filtered_observations(
                all_observations,
                held_documents=fold_documents[int(fold)],
            )
            delta = delta_document_map(
                selected_pairs,
                expand_documents=True,
            )
            for item in fold_frame.to_dict("records"):
                row = pd.Series(item)
                remaining = remaining_units_for_row(
                    groups[str(row.entry)], row
                )
                remaining = remaining.loc[
                    remaining.document_fold != int(row.fold)
                ]
                _, touched = unit_deltas(remaining, delta)
                outputs.append(
                    {
                        **item,
                        "training_document_flag_touched": bool(touched.any()),
                    }
                )
    result = pd.DataFrame(outputs)
    result["g_supported"] = result.support >= SUPPORT_THRESHOLD
    result["p_supported"] = (
        result.g_supported & ~result.training_document_flag_touched
    )
    return result


def e11_v1(
    frame: pd.DataFrame,
    eligibility_frame: pd.DataFrame | None = None,
) -> dict:
    rows = []
    eligible_by_target = None
    if eligibility_frame is not None:
        eligible_by_target = set(
            eligibility_frame.loc[
                (eligibility_frame.support >= SUPPORT_THRESHOLD)
                & eligibility_frame.flag_touched.astype(bool),
                "target",
            ].astype(str)
        )
    for target, group in frame.groupby("target", sort=True):
        g_supported = group.g_supported.to_numpy(bool)
        p_supported = group.p_supported.to_numpy(bool)
        eligible = (
            str(target) in eligible_by_target
            if eligible_by_target is not None
            else (
                (g_supported & group.leakfree_flag_touched.to_numpy(bool)).any()
                and g_supported.any()
                and p_supported.any()
            )
        )
        if eligible:
            if not (g_supported.any() and p_supported.any()):
                rows.append(np.nan)
                continue
            rows.append(
                float(
                    group.loc[p_supported, "agree"].mean()
                    - group.loc[g_supported, "agree"].mean()
                )
            )
    values = np.asarray(rows, dtype=float)
    if not len(values) or np.allclose(values, 0):
        p_value = 1.0
        statistic = np.nan
        method = "not_applicable"
    else:
        method = "exact" if len(values) <= 25 and len(np.unique(np.abs(values))) == len(values) else "approx"
        test = stats.wilcoxon(
            values,
            alternative="greater",
            zero_method="wilcox",
            method=method,
        )
        p_value = float(test.pvalue)
        statistic = float(test.statistic)
    return {
        "n_datasets": len(values),
        "median_difference": float(np.median(values)) if len(values) else np.nan,
        "positive_datasets": int((values > 0).sum()),
        "zero_datasets": int((values == 0).sum()),
        "negative_datasets": int((values < 0).sum()),
        "wilcoxon_statistic": statistic,
        "one_sided_p": p_value,
        "wilcoxon_method": method,
    }


def regression_checks(
    frames: dict[str, pd.DataFrame],
    records_by_key: dict,
    e11_records: dict,
    r5: pd.DataFrame,
) -> tuple[dict, dict]:
    archived = archived_rule_frame()
    entry = archived.loc[
        (archived.split_design == "entry")
        & (archived.merge_basis == "same_value")
    ]
    document = archived.loc[
        (archived.split_design == "document")
        & (archived.merge_basis == "same_value")
    ]
    r1_entry = recompute_archived_flag_touched(records_by_key, entry)
    r1_document = recompute_archived_flag_touched(records_by_key, document)
    r1_e11 = recompute_archived_flag_touched(
        e11_records, archived_e11_frame()
    )
    r1 = {
        "entry_rows": int(len(r1_entry)),
        "entry_matches": int(r1_entry.matches.sum()),
        "entry_all_match": bool(r1_entry.matches.all()),
        "document_rows": int(len(r1_document)),
        "document_matches": int(r1_document.matches.sum()),
        "document_all_match": bool(r1_document.matches.all()),
        "e11_rows": int(len(r1_e11)),
        "e11_matches": int(r1_e11.matches.sum()),
        "e11_all_match": bool(r1_e11.matches.all()),
        "all_match": bool(
            r1_entry.matches.all()
            and r1_document.matches.all()
            and r1_e11.matches.all()
        ),
    }
    if not r1["all_match"]:
        raise RuntimeError(f"R1 failed: {r1}")

    from e5d_models import add_rule_states

    archived_states = pd.read_csv(E12A / "rule_states_by_design.csv")
    r2_rows = []
    for design_name, design_frame in (
        ("entry", entry),
        ("document", document),
    ):
        states = add_rule_states(design_frame.copy())
        for endpoint in ENDPOINTS:
            selected = states.loc[
                (states.endpoint == endpoint)
                & (states.p_status == "supported")
            ]
            expected = archived_states.loc[
                (archived_states.design == design_name)
                & (archived_states.endpoint == endpoint)
                & (archived_states.state == "P_supported")
            ].iloc[0]
            value = float(selected.agree.mean())
            r2_rows.append(
                {
                    "design": design_name,
                    "endpoint": endpoint,
                    "recomputed_supported_agreement": value,
                    "archived_supported_agreement": float(
                        expected.agreement_rate
                    ),
                    "absolute_difference": abs(
                        value - float(expected.agreement_rate)
                    ),
                    "passes_1e_12": bool(
                        abs(value - float(expected.agreement_rate)) <= 1e-12
                    ),
                }
            )
    r2 = pd.DataFrame(r2_rows)
    if not r2.passes_1e_12.all():
        raise RuntimeError("R2 supported-agreement regression failed")

    r3_rows = []
    document_training = original_document_training_flags(records_by_key)
    for endpoint in ENDPOINTS:
        entry_frame = frames["entry"]
        endpoint_entry = entry_frame.loc[entry_frame.endpoint == endpoint]
        entry_rule = float(
            endpoint_entry.loc[
                endpoint_entry.p_supported, "agree"
            ].mean()
        )
        endpoint_document = document_training.loc[
            document_training.endpoint == endpoint
        ]
        document_rule = float(
            endpoint_document.loc[
                endpoint_document.p_supported, "agree"
            ].mean()
        )
        old_e11 = archived_e11_frame()
        e11_test = e11_v1(
            frames["E11"].loc[frames["E11"].endpoint == endpoint],
            eligibility_frame=old_e11.loc[
                old_e11.endpoint == endpoint
            ],
        )
        r3_rows.append(
            {
                "endpoint": endpoint,
                "entry_supported_agreement": entry_rule,
                "training_document_supported_agreement": document_rule,
                **e11_test,
            }
        )
    r3 = pd.DataFrame(r3_rows)

    old_document = archived_document_g()
    replay_frames = []
    for (target, endpoint), records in sorted(records_by_key.items()):
        replay_frames.append(
            document_fold_predictions_e15(
                records,
                target,
                endpoint,
                strict_leak_free=False,
            )
        )
    replay = pd.concat(replay_frames, ignore_index=True)
    keys = ["target", "endpoint", "fold", "entry", "held_out_unit"]
    columns = [
        "displayed",
        "observed",
        "p0",
        "p1",
        "p2",
        "p_displayed",
        "support",
        "agree",
    ]
    old = old_document.sort_values(keys).reset_index(drop=True)
    new = replay.sort_values(keys).reset_index(drop=True)
    key_match = old[keys].equals(new[keys])
    numeric_match = bool(
        np.allclose(
            old[columns].to_numpy(float),
            new[columns].to_numpy(float),
            rtol=0,
            atol=1e-12,
        )
    )
    r4 = {
        "archived_rows": int(len(old)),
        "replayed_rows": int(len(new)),
        "key_match": bool(key_match),
        "numeric_match_atol_1e_12": numeric_match,
        "passed": bool(len(old) == len(new) and key_match and numeric_match),
    }
    if not r4["passed"]:
        raise RuntimeError(f"R4 replay failed: {r4}")
    r5_summary = {
        "rows": int(len(r5)),
        "matches": int(r5.matches.sum()),
        "passed": bool(r5.matches.all()),
    }
    checks = {"R1": r1, "R2": r2_rows, "R3": r3_rows, "R4": r4, "R5": r5_summary}
    return checks, {
        "r1_entry": r1_entry,
        "r1_document": r1_document,
        "r1_e11": r1_e11,
        "r2": r2,
        "r3": r3,
        "r5": r5,
    }


def decision(
    q1_tests: pd.DataFrame,
    q2_tests: pd.DataFrame,
    q2_comparisons: pd.DataFrame,
    q3_decision: dict,
) -> dict:
    q1_endpoints = {}
    for endpoint in ENDPOINTS:
        endpoint_tests = q1_tests.loc[q1_tests.endpoint == endpoint]
        strict = endpoint_tests.loc[
            endpoint_tests.design == "strict_document"
        ].iloc[0]
        time = endpoint_tests.loc[
            endpoint_tests.design == "time_2015"
        ].iloc[0]
        decline = endpoint_tests.loc[
            (endpoint_tests.holm_two_sided_p < 0.05)
            & (endpoint_tests.median_difference < 0)
        ]
        passed = bool(
            strict.median_difference > 0
            and strict.holm_two_sided_p < 0.05
            and time.median_difference >= 0
            and decline.empty
        )
        q1_endpoints[endpoint] = {
            "confirmed": passed,
            "strict_document_median_rise": float(
                strict.median_difference
            ),
            "strict_document_holm_p": float(strict.holm_two_sided_p),
            "time_2015_median_rise": float(time.median_difference),
            "significant_declines": int(len(decline)),
        }

    q2_endpoints = {}
    for endpoint in ENDPOINTS:
        endpoint_tests = q2_tests.loc[q2_tests.endpoint == endpoint]
        candidates = endpoint_tests.loc[
            endpoint_tests.design.isin(("entry", "strict_document"))
            & (endpoint_tests.median_difference < 0)
            & (endpoint_tests.holm_two_sided_p < 0.05)
        ].sort_values("holm_two_sided_p")
        time = endpoint_tests.loc[
            endpoint_tests.design == "time_2015"
        ].iloc[0]
        reasons = []
        for candidate in candidates.itertuples(index=False):
            comparison = q2_comparisons.loc[
                (q2_comparisons.design == candidate.design)
                & (q2_comparisons.endpoint == endpoint)
            ].iloc[0]
            supports = bool(
                comparison.st_minus_u > 0
                and comparison.ci_low > 0
            )
            reasons.append(
                {
                    "design": candidate.design,
                    "dataset_median_St_minus_Su": float(
                        candidate.median_difference
                    ),
                    "holm_p": float(candidate.holm_two_sided_p),
                    "time_2015_median_St_minus_Su": float(
                        time.median_difference
                    ),
                    "St_minus_U": float(comparison.st_minus_u),
                    "St_minus_U_ci_low": float(comparison.ci_low),
                    "St_minus_U_ci_high": float(comparison.ci_high),
                    "supported": bool(
                        time.median_difference <= 0 and supports
                    ),
                }
            )
        q2_endpoints[endpoint] = {
            "supported": any(item["supported"] for item in reasons),
            "candidates": reasons,
        }
    q1_confirmed = any(item["confirmed"] for item in q1_endpoints.values())
    q2_confirmed = any(item["supported"] for item in q2_endpoints.values())
    q3_label = q3_decision["decision"]
    if q1_confirmed:
        row = "Q1 confirmed: rule P is a validated display rule"
    elif q2_confirmed and q3_label in {"strong", "weak"}:
        row = (
            "Q1 not confirmed; Q2 supported; Q3 strong/weak: use three-state "
            "display and report it as a data-audit tool"
        )
    elif q2_confirmed and q3_label == "none":
        row = (
            "Q1 not confirmed; Q2 supported; Q3 none: three-state display is "
            "descriptive, not a primary innovation"
        )
    elif q3_label == "strong":
        row = (
            "Q1 not confirmed; Q2 unsupported; Q3 strong: treat offset pairs "
            "as a data-audit tool and remove default downgrading"
        )
    else:
        row = (
            "Q1 not confirmed; Q2 unsupported; Q3 weak/none: rule P is only "
            "an optional transparent display with no claimed benefit"
        )
    return {
        "analysis_role": ANALYSIS_ROLE,
        "Q1": {
            "decision": "confirmed" if q1_confirmed else "not_confirmed",
            "endpoints": q1_endpoints,
        },
        "Q2": {
            "decision": "supported" if q2_confirmed else "unsupported",
            "endpoints": q2_endpoints,
        },
        "Q3": q3_decision,
        "section_4_row": row,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    print("loading 52-dataset records", flush=True)
    records_by_key = load_e2_records()
    e11_cache = QA / "curated_e11.pkl.gz"
    if e11_cache.exists():
        print("loading cached E11 curation", flush=True)
        e11_records = pickle.loads(gzip.decompress(e11_cache.read_bytes()))
        r5 = pd.read_csv(QA / "r5.csv")
    else:
        print("re-curating E11 records (R5)", flush=True)
        e11_records, r5 = curate_e11({})
        e11_cache.write_bytes(
            gzip.compress(pickle.dumps(e11_records, protocol=5))
        )
        r5.to_csv(QA / "r5.csv", index=False)
    print("loading archived predictions", flush=True)

    archived = archived_rule_frame()
    design_cache = QA / "design_frames.pkl.gz"
    if os.environ.get("E15_REUSE_FRAMES") == "1" and design_cache.exists():
        print("loading cached design frames", flush=True)
        design_frames = pickle.loads(
            gzip.decompress(design_cache.read_bytes())
        )
    else:
        design_frames = {}
        print("building entry-fold leak-free flags", flush=True)
        design_frames["entry"] = build_design_frame(
            archived, records_by_key, design="entry"
        )
        print("building strict document-fold predictions", flush=True)
        design_frames["strict_document"] = build_strict_document_frame(
            records_by_key
        )
        print("building near-duplicate leak-free flags", flush=True)
        design_frames["near_duplicate"] = build_design_frame(
            archived, records_by_key, design="near_duplicate"
        )
        print("building 2015 time-split frame", flush=True)
        design_frames["time_2015"] = build_design_frame(
            archived, records_by_key, design="time_2015"
        )
        print("building E11 leak-free flags", flush=True)
        design_frames["E11"] = build_design_frame(
            archived_e11_frame(), e11_records, design="E11"
        )
        design_cache.write_bytes(
            gzip.compress(pickle.dumps(design_frames, protocol=5))
        )

    print("running R1-R5 regression checks", flush=True)
    checks, check_tables = regression_checks(
        design_frames, records_by_key, e11_records, r5
    )
    for name, table in check_tables.items():
        table.to_csv(QA / f"{name}.csv", index=False)
    (QA / "regression_checks.json").write_text(
        json.dumps(json_safe(checks), indent=2, ensure_ascii=False)
        + "\n",
        encoding="utf-8",
    )
    print("R1-R5 passed", flush=True)

    q1_table, q1_tests = q1(design_frames)
    q1_table.to_csv(OUT / "q1_rule_by_design.csv", index=False)
    q1_tests.to_csv(OUT / "q1_tests.csv", index=False)

    without_orexin = {
        design: frame.loc[~frame.target.isin(OREXIN_TARGETS)].reset_index(drop=True)
        for design, frame in design_frames.items()
    }
    orexin_q1, orexin_tests = q1(without_orexin)
    pd.concat(
        [
            orexin_q1.assign(row_type="summary"),
            orexin_tests.assign(row_type="test"),
        ],
        ignore_index=True,
        sort=False,
    ).to_csv(OUT / "q1_without_orexin.csv", index=False)

    evaluation_pairs = {}
    all_evaluation_pairs = set()
    for key, records in records_by_key.items():
        flagged = full_data_flagged_pairs(records)
        evaluation_pairs[key] = flagged
        all_evaluation_pairs |= flagged
    e11_new_pair_frame = build_design_frame(
        archived_e11_frame(),
        e11_records,
        design="E11",
        excluded_pairs_by_dataset={
            key: all_evaluation_pairs for key in e11_records
        },
    )
    e11_only_frames = {**design_frames, "E11": e11_new_pair_frame}
    e11_q1, e11_tests = q1(e11_only_frames)
    e11_q1 = e11_q1.loc[e11_q1.design == "E11"]
    e11_tests = e11_tests.loc[e11_tests.design == "E11"]
    pd.concat(
        [
            e11_q1.assign(row_type="summary"),
            e11_tests.assign(row_type="test"),
        ],
        ignore_index=True,
        sort=False,
    ).to_csv(OUT / "q1_e11_new_pairs.csv", index=False)

    evaluation_documents = {
        document
        for records in records_by_key.values()
        for document in records.document_chembl_id.dropna().astype(str)
    }
    overlap_rows = []
    for endpoint in ENDPOINTS:
        e11_documents = {
            document
            for key, records in e11_records.items()
            if key[1] == endpoint
            for document in records.document_chembl_id.dropna().astype(str)
        }
        overlap_rows.append(
            {
                "scope": "endpoint",
                "endpoint": endpoint,
                "e11_documents": len(e11_documents),
                "evaluation_documents": len(evaluation_documents),
                "shared_documents": len(e11_documents & evaluation_documents),
                "shared_fraction_of_e11": (
                    len(e11_documents & evaluation_documents)
                    / len(e11_documents)
                    if e11_documents
                    else np.nan
                ),
            }
        )
    all_e11_documents = {
        document
        for records in e11_records.values()
        for document in records.document_chembl_id.dropna().astype(str)
    }
    overlap_rows.append(
        {
            "scope": "all",
            "endpoint": "ALL",
            "e11_documents": len(all_e11_documents),
            "evaluation_documents": len(evaluation_documents),
            "shared_documents": len(all_e11_documents & evaluation_documents),
            "shared_fraction_of_e11": (
                len(all_e11_documents & evaluation_documents)
                / len(all_e11_documents)
            ),
        }
    )
    pd.DataFrame(overlap_rows).to_csv(
        OUT / "e11_document_overlap.csv", index=False
    )

    q2_table, q2_tests, q2_comparisons = q2(design_frames)
    q2_table.to_csv(OUT / "q2_three_state.csv", index=False)
    q2_tests.to_csv(OUT / "q2_tests.csv", index=False)
    q2_comparisons.to_csv(
        QA / "q2_St_minus_U_bootstrap.csv", index=False
    )

    q4_combined, q4_summary = q4(design_frames)
    q4_combined.to_csv(
        OUT / "q4_calibration_dataset_level.csv", index=False
    )
    q4_summary.to_csv(QA / "q4_summary.csv", index=False)

    q5(design_frames).to_csv(
        OUT / "q5_offset_association.csv", index=False
    )

    old_document = archived_document_g()
    strict_effect = strict_document_effect(
        old_document, design_frames["strict_document"]
    )
    strict_effect.to_csv(
        OUT / "strict_document_fold_effect.csv", index=False
    )

    from e15_offset_audit import run as run_offset_audit

    q3_decision = run_offset_audit(
        records_by_key=records_by_key,
        e11_records=e11_records,
        out_dir=OUT,
    )
    final_decision = decision(
        q1_tests, q2_tests, q2_comparisons, q3_decision
    )
    write_json(OUT / "decision.json", final_decision)

    finished_at = datetime.now().astimezone().isoformat(timespec="seconds")
    input_paths = [
        S / "qa/E2/curated_records.pkl.gz",
        E5C / "predictions_e5c.csv",
        E5 / "document_folds_predictions.csv",
        E11 / "predictions_e11.csv",
        E11 / "target_inclusion.csv",
        VALIDATION_DATA / "selected_targets.csv",
        E12A / "rule_states_by_design.csv",
        E9 / "calibration_summary.csv",
    ]
    source_paths = [
        SOURCE / "e15_leakfree.py",
        SOURCE / "e15_offset_audit.py",
        SOURCE / "e5_models.py",
        SOURCE / "e5c_models.py",
        SOURCE / "e2_models.py",
        SOURCE / "e7_analysis.py",
        SOURCE / "e9_analysis.py",
        SOURCE / "e11_analysis.py",
        SOURCE / "e12a_rule_calibration.py",
        SOURCE / "manifest_utils.py",
    ]
    write_manifest(
        OUT,
        stage="E15",
        inputs=input_paths,
        sources=source_paths,
        parameters={
            "analysis_role": ANALYSIS_ROLE,
            "bootstrap_repetitions": BOOT_REPETITIONS,
            "bootstrap_seed": BOOT_SEED,
            "wilcoxon_zero_method": "wilcox",
            "wilcoxon_method_rule": (
                "exact when n<=25 and no tied non-zero absolute values; "
                "otherwise approx"
            ),
            "designs": list(DESIGNS),
            "ordinary_holm_family_size": 10,
        },
        started_at=STARTED_AT,
        finished_at=finished_at,
    )
    elapsed = time.perf_counter() - STARTED_PERF
    print(f"E15 complete in {elapsed:.1f} seconds", flush=True)


if __name__ == "__main__":
    main()
