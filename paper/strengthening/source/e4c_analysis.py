"""E4c：文献对标记零分布、训练折 flag-clean 与近重复合并。"""
from __future__ import annotations

import gzip
from collections import defaultdict
from itertools import combinations
from pathlib import Path
import pickle
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[4]
CGM = ROOT / "chemgridmap"
REV = ROOT / "project1/JCIM_reliability_revision_20261024"
S = REV / "strengthening"
OUT = S / "output/E4"
QA = S / "qa/E4"
sys.path[:0] = [str(S / "source"), str(CGM / "paper"), str(CGM / "src")]

import numpy as np
import pandas as pd

from chemgridmap.activity import activity_class
from chemgridmap.evidence import CLASSES, document_values, fit_between_unit_sd
from e1a_statistics import auc_mann_whitney
from e2_models import (
    FOLD_SEED,
    document_pair_observations,
    flags_from_observations,
    touch_flagged_entries,
)
from reliability_activity import entry_folds, gaussian_prediction

ROLE = "post hoc"
NULL_SEED = 20261103
NULL_REPETITIONS = 200
NEAR_DUPLICATE_THRESHOLD = 0.1
SUPPORT_THRESHOLD = 0.8
LOWER, UPPER = 6.0, 7.0


def load_inputs() -> dict[str, object]:
    """只读 E1/E2 已验收输出与整理记录缓存。"""
    predictions = pd.read_csv(
        S / "output/E2/predictions.csv", float_precision="round_trip"
    )
    keys = (
        predictions[["target", "endpoint"]]
        .drop_duplicates()
        .sort_values(["endpoint", "target"])
        .itertuples(index=False, name=None)
    )
    records = pickle.loads(
        gzip.decompress((S / "qa/E2/curated_records.pkl.gz").read_bytes())
    )
    return {
        "records": {key: records[key] for key in keys},
        "predictions": predictions,
        "census": pd.read_csv(
            S / "output/E2/flagged_pairs_census.csv",
            float_precision="round_trip",
        ),
        "sigma": pd.read_csv(
            S / "output/E2/sigma_doc_by_target.csv",
            float_precision="round_trip",
        ),
        "e1a_predictions": pd.read_csv(
            S / "output/E1A/predictions.csv", float_precision="round_trip"
        ),
        "herg_records": pd.read_csv(
            CGM / "paper/output/revision_v6/chembl240/chembl240_retained_records.csv",
            low_memory=False,
            float_precision="round_trip",
        ),
    }


def shared_entry_stratum(n_shared: int) -> str:
    n_shared = int(n_shared)
    if n_shared in (3, 4):
        return str(n_shared)
    return ">=5"


def pair_structure(records: pd.DataFrame) -> dict[tuple[str, str], np.ndarray]:
    """生成可重复模拟的文献对索引；仅保留至少有 3 个共享条目的文献对。"""
    units = document_values(records, merge_identical=False).sort_values(
        ["entry", "document"]
    )
    pairs: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    for _, group in units.groupby("entry", sort=True):
        documents = group.document.astype(str).to_numpy()
        positions = group.index.to_numpy()
        for left, right in combinations(range(len(documents)), 2):
            pairs[(documents[left], documents[right])].append(
                (int(positions[left]), int(positions[right]))
            )
    return {
        pair: np.asarray(indices, dtype=int)
        for pair, indices in sorted(pairs.items())
        if len(indices) >= 3
    }


def simulate_null_flags(
    structure: dict[tuple[str, str], np.ndarray],
    n_units: int,
    sigma_doc: float,
    rng: np.random.Generator,
    repetitions: int = NULL_REPETITIONS,
) -> dict[str, np.ndarray]:
    """按原条目×文档结构模拟 Gaussian 单元值，并复用 CGM 默认标记规则。"""
    values = rng.normal(0.0, float(sigma_doc), size=(repetitions, n_units))
    counts = {
        "total": np.zeros(repetitions, dtype=int),
        "3": np.zeros(repetitions, dtype=int),
        "4": np.zeros(repetitions, dtype=int),
        ">=5": np.zeros(repetitions, dtype=int),
    }
    for indices in structure.values():
        differences = values[:, indices[:, 1]] - values[:, indices[:, 0]]
        medians = np.median(differences, axis=1)
        signs = np.sign(medians)
        agreement = np.mean(differences * signs[:, None] > 0, axis=1)
        flagged = (np.abs(medians) >= 1.0) & (agreement >= 0.8)
        counts["total"] += flagged
        counts[shared_entry_stratum(len(indices))] += flagged
    return counts


def flag_null_analysis(
    records: dict[tuple[str, str], pd.DataFrame],
    census: pd.DataFrame,
    sigma: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, object]]:
    summaries = census.loc[census.row_type == "target_summary"].set_index(
        ["target", "endpoint"]
    )
    pairs = census.loc[census.row_type == "top_pair"]
    sigma_rows = sigma.loc[sigma.row_type == "target"].set_index(
        ["target", "endpoint"]
    )
    rng = np.random.default_rng(NULL_SEED)
    rows = []
    stratum_expected: dict[tuple[str, str], list[float]] = defaultdict(list)
    for target, endpoint in records:
        frame = records[(target, endpoint)]
        structure = pair_structure(frame)
        sigma_doc = float(sigma_rows.loc[(target, endpoint), "sigma_doc"])
        counts = simulate_null_flags(
            structure, len(document_values(frame, merge_identical=False)), sigma_doc, rng
        )
        for stratum in ("3", "4", ">=5"):
            stratum_expected[(endpoint, stratum)].append(float(counts[stratum].mean()))
        actual = int(summaries.loc[(target, endpoint), "n_flagged_pairs"])
        rows.append(
            {
                "row_type": "dataset",
                "target": target,
                "endpoint": endpoint,
                "n_units": int(len(document_values(frame, merge_identical=False))),
                "n_candidate_pairs": len(structure),
                "sigma_doc": sigma_doc,
                "simulation_repetitions": NULL_REPETITIONS,
                "seed": NULL_SEED,
                "expected_false_flagged_pairs": float(counts["total"].mean()),
                "actual_flagged_pairs": actual,
                "analysis_role": ROLE,
            }
        )
    detail = pd.DataFrame(rows)

    summary_rows = []
    for endpoint in ("IC50", "Ki"):
        selected = detail.loc[detail.endpoint == endpoint]
        endpoint_pairs = pairs.loc[pairs.endpoint == endpoint]
        for stratum in ("3", "4", ">=5"):
            pair_values = endpoint_pairs.loc[
                endpoint_pairs.shared_entries.map(shared_entry_stratum) == stratum
            ]
            summary_rows.append(
                {
                    "row_type": "shared_stratum_total",
                    "target": "",
                    "endpoint": endpoint,
                    "shared_entries_stratum": stratum,
                    "simulation_repetitions": NULL_REPETITIONS,
                    "seed": NULL_SEED,
                    "expected_false_flagged_pairs": float(
                        sum(stratum_expected[(endpoint, stratum)])
                    ),
                    "actual_flagged_pairs": int(len(pair_values)),
                    "analysis_role": ROLE,
                }
            )
        summary_rows.append(
            {
                "row_type": "endpoint_total",
                "target": "",
                "endpoint": endpoint,
                "n_datasets": int(len(selected)),
                "simulation_repetitions": NULL_REPETITIONS,
                "seed": NULL_SEED,
                "expected_false_flagged_pairs": float(
                    selected.expected_false_flagged_pairs.sum()
                ),
                "actual_flagged_pairs": int(selected.actual_flagged_pairs.sum()),
                "analysis_role": ROLE,
            }
        )
    summary_rows.append(
        {
            "row_type": "grand_total",
            "target": "",
            "endpoint": "",
            "n_datasets": len(detail),
            "simulation_repetitions": NULL_REPETITIONS,
            "seed": NULL_SEED,
            "expected_false_flagged_pairs": float(
                detail.expected_false_flagged_pairs.sum()
            ),
            "actual_flagged_pairs": int(detail.actual_flagged_pairs.sum()),
            "analysis_role": ROLE,
        }
    )
    result = pd.concat([detail, pd.DataFrame(summary_rows)], ignore_index=True)
    endpoint_totals = result.loc[result.row_type == "endpoint_total"].set_index(
        "endpoint"
    )
    checks = {
        endpoint: {
            "expected": float(endpoint_totals.loc[endpoint, "expected_false_flagged_pairs"]),
            "actual": int(endpoint_totals.loc[endpoint, "actual_flagged_pairs"]),
        }
        for endpoint in ("IC50", "Ki")
    }
    if checks["IC50"]["actual"] != 39 or checks["Ki"]["actual"] != 36:
        raise ValueError(f"E2 标记对总数与原验收结果不一致：{checks}")
    return result, checks


def score_g(frame: pd.DataFrame) -> dict[str, float]:
    """与 E2 一致的 G AUC 与高置信度准确率定义。"""
    probabilities = frame[["p0", "p1", "p2"]].to_numpy(float)
    observed = frame.observed.to_numpy(int)
    confidence = probabilities.max(axis=1)
    correct = probabilities.argmax(axis=1) == observed
    high = confidence >= SUPPORT_THRESHOLD
    return {
        "auc": auc_mann_whitney(
            frame.observed.to_numpy(int) == frame.displayed.to_numpy(int),
            frame.p_displayed.to_numpy(float),
        ),
        "high_confidence_n": int(high.sum()),
        "high_confidence_accuracy": (
            float(correct[high].mean()) if high.any() else np.nan
        ),
    }


def training_fold_g(
    records: pd.DataFrame,
    b0: pd.DataFrame,
    target: str,
    endpoint: str,
    source: str,
) -> pd.DataFrame:
    """只用当前折训练条目识别标记文献，并重算保留条目的 G 预测。"""
    doc = document_values(records, merge_identical=True)
    observations = document_pair_observations(records)
    assignments = entry_folds(doc, seed=FOLD_SEED)
    units_by_entry = dict(tuple(doc.groupby("entry", sort=True)))
    rows = []
    for fold in range(10):
        test_entries = {
            entry for entry, assigned in assignments.items() if assigned == fold
        }
        train_entries = set(doc.entry.unique()) - test_entries
        training_flags = flags_from_observations(observations, train_entries)
        train_doc = doc.loc[doc.entry.isin(train_entries)]
        touched_entries = touch_flagged_entries(doc, training_flags)
        clean_train = train_doc.loc[
            ~train_doc.entry.isin(touched_entries)
        ]
        fit = fit_between_unit_sd(clean_train)
        fold_predictions = b0.loc[b0.fold == fold]
        for held in fold_predictions.itertuples(index=False):
            if int(assignments[held.entry]) != fold:
                raise ValueError("E4c 分折与 E1/E2 冻结分折不匹配。")
            if held.entry in touched_entries:
                continue
            held_documents = set(str(held.held_out_unit).split(";"))
            units = units_by_entry[held.entry]
            remaining = units.loc[
                ~units.documents.map(
                    lambda value: bool(
                        set(str(value).split(";")) & held_documents
                    )
                )
            ]
            if remaining.empty:
                raise ValueError("训练折 flag-clean 留出后没有剩余证据单元。")
            probabilities = gaussian_prediction(
                remaining.document_value.to_numpy(float), fit
            )
            rows.append(
                {
                    "source": source,
                    "target": target,
                    "endpoint": endpoint,
                    "fold": fold,
                    "entry": held.entry,
                    "held_out_unit": held.held_out_unit,
                    "observed": int(held.observed),
                    "displayed": int(held.displayed),
                    "p0": float(probabilities[0]),
                    "p1": float(probabilities[1]),
                    "p2": float(probabilities[2]),
                    "p_displayed": float(probabilities[int(held.displayed)]),
                    "training_flagged_documents": len(training_flags),
                    "clean_training_entries": int(clean_train.entry.nunique()),
                }
            )
    result = pd.DataFrame(rows)
    result["analysis_role"] = ROLE
    return result


def flag_clean_trainfold_analysis(inputs: dict[str, object]) -> pd.DataFrame:
    predictions = inputs["predictions"]
    e1a = inputs["e1a_predictions"]
    rows = []

    e1a_original = e1a.loc[
        (e1a.target == "chembl240")
        & (e1a.population == "flag-clean")
        & (e1a.method == "G")
    ]
    e1a_b0 = e1a.loc[
        (e1a.target == "chembl240")
        & (e1a.population == "all")
        & (e1a.method == "B0")
    ]
    e1a_new = training_fold_g(
        inputs["herg_records"], e1a_b0, "chembl240", "IC50", "E1A"
    )
    original_scores = score_g(e1a_original)
    new_scores = score_g(e1a_new)
    rows.append(
        {
            "row_type": "dataset",
            "source": "E1A",
            "target": "chembl240",
            "endpoint": "IC50",
            "n_predictions_original": len(e1a_original),
            "n_predictions_trainfold": len(e1a_new),
            "auc_original_flag_clean": original_scores["auc"],
            "auc_trainfold_flag_clean": new_scores["auc"],
            "auc_difference": new_scores["auc"] - original_scores["auc"],
            "high_confidence_n_original": original_scores["high_confidence_n"],
            "high_confidence_accuracy_original": original_scores[
                "high_confidence_accuracy"
            ],
            "high_confidence_n_trainfold": new_scores["high_confidence_n"],
            "high_confidence_accuracy_trainfold": new_scores[
                "high_confidence_accuracy"
            ],
            "high_confidence_accuracy_difference": new_scores[
                "high_confidence_accuracy"
            ]
            - original_scores["high_confidence_accuracy"],
            "analysis_role": ROLE,
        }
    )

    for (target, endpoint), records in inputs["records"].items():
        original = predictions.loc[
            (predictions.target == target)
            & (predictions.endpoint == endpoint)
            & (predictions.population == "flag-clean")
            & (predictions.method == "G")
        ]
        b0 = predictions.loc[
            (predictions.target == target)
            & (predictions.endpoint == endpoint)
            & (predictions.population == "all")
            & (predictions.method == "B0")
        ]
        new = training_fold_g(records, b0, target, endpoint, "E2")
        original_scores = score_g(original)
        new_scores = score_g(new)
        rows.append(
            {
                "row_type": "dataset",
                "source": "E2",
                "target": target,
                "endpoint": endpoint,
                "n_predictions_original": len(original),
                "n_predictions_trainfold": len(new),
                "auc_original_flag_clean": original_scores["auc"],
                "auc_trainfold_flag_clean": new_scores["auc"],
                "auc_difference": new_scores["auc"] - original_scores["auc"],
                "high_confidence_n_original": original_scores[
                    "high_confidence_n"
                ],
                "high_confidence_accuracy_original": original_scores[
                    "high_confidence_accuracy"
                ],
                "high_confidence_n_trainfold": new_scores["high_confidence_n"],
                "high_confidence_accuracy_trainfold": new_scores[
                    "high_confidence_accuracy"
                ],
                "high_confidence_accuracy_difference": new_scores[
                    "high_confidence_accuracy"
                ]
                - original_scores["high_confidence_accuracy"],
                "analysis_role": ROLE,
            }
        )
    return pd.DataFrame(rows)


def near_duplicate_units(records: pd.DataFrame) -> pd.DataFrame:
    """把同一目中相邻值差不超过 0.1 的文献单元做单链接聚类。"""
    units = document_values(records, merge_identical=False).sort_values(
        ["entry", "document_value", "document"]
    )
    rows = []
    for entry, group in units.groupby("entry", sort=True):
        values = group.document_value.to_numpy(float)
        starts = np.r_[
            0,
            np.flatnonzero(np.diff(values) > NEAR_DUPLICATE_THRESHOLD + 1e-12)
            + 1,
            len(values),
        ]
        for start, stop in zip(starts[:-1], starts[1:]):
            selected = group.iloc[start:stop]
            rows.append(
                {
                    "entry": entry,
                    "document": str(selected.document.iloc[0]),
                    "document_value": float(selected.document_value.mean()),
                    "n_records": int(selected.n_records.sum()),
                    "documents": ";".join(selected.document.astype(str)),
                }
            )
    return pd.DataFrame(rows)


def near_duplicate_g(
    records: pd.DataFrame,
    b0: pd.DataFrame,
    target: str,
    endpoint: str,
) -> tuple[pd.DataFrame, dict[str, float]]:
    units = near_duplicate_units(records)
    fit = fit_between_unit_sd(units)
    units_by_entry = dict(tuple(units.groupby("entry", sort=True)))
    records_by_entry = dict(tuple(records.groupby("molecule_identity_key", sort=True)))
    entry_fold_rows = b0[["entry", "fold"]].drop_duplicates()
    if entry_fold_rows.entry.duplicated().any():
        raise ValueError("一个条目被分配到了多个 E2 冻结折。")
    rows = []
    for held_entry in entry_fold_rows.itertuples(index=False):
        available = units_by_entry[held_entry.entry]
        if len(available) < 2:
            continue
        source = records_by_entry[held_entry.entry]
        for held_index, held_cluster in available.iterrows():
            held_documents = set(str(held_cluster.documents).split(";"))
            remaining_units = available.drop(index=held_index)
            remaining_records = source.loc[
                ~source.document_chembl_id.astype(str).isin(held_documents),
                "pchembl_value",
            ].to_numpy(float)
            if not len(remaining_records):
                continue
            observed = CLASSES.index(
                activity_class(float(held_cluster.document_value), LOWER, UPPER)
            )
            displayed = CLASSES.index(
                activity_class(float(np.median(remaining_records)), LOWER, UPPER)
            )
            probabilities = gaussian_prediction(
                remaining_units.document_value.to_numpy(float), fit
            )
            rows.append(
                {
                    "target": target,
                    "endpoint": endpoint,
                    "fold": int(held_entry.fold),
                    "entry": held_entry.entry,
                    "held_out_unit": held_cluster.documents,
                    "observed": observed,
                    "displayed": displayed,
                    "p0": float(probabilities[0]),
                    "p1": float(probabilities[1]),
                    "p2": float(probabilities[2]),
                    "p_displayed": float(probabilities[displayed]),
                }
            )
    result = pd.DataFrame(rows)
    result["agree"] = result.observed == result.displayed
    p = result[["p0", "p1", "p2"]].to_numpy(float)
    observed = result.observed.to_numpy(int)
    scores = {
        "auc": auc_mann_whitney(
            result.agree.to_numpy(bool), result.p_displayed.to_numpy(float)
        ),
        "log_loss": float(
            -np.log(np.clip(p[np.arange(len(p)), observed], 1e-6, 1.0)).mean()
        ),
        "sigma_doc": float(fit["tau"]),
        "n_units": int(len(units)),
        "n_predictions": int(len(result)),
    }
    return result, scores


def near_duplicate_analysis(inputs: dict[str, object]) -> pd.DataFrame:
    predictions = inputs["predictions"]
    sigma = inputs["sigma"].loc[inputs["sigma"].row_type == "target"].set_index(
        ["target", "endpoint"]
    )
    rows = []
    focus_targets = {"CHEMBL5024", "CHEMBL244", "CHEMBL2000"}
    for (target, endpoint), records in inputs["records"].items():
        original = predictions.loc[
            (predictions.target == target)
            & (predictions.endpoint == endpoint)
            & (predictions.population == "all")
            & (predictions.method == "G")
        ]
        b0 = predictions.loc[
            (predictions.target == target)
            & (predictions.endpoint == endpoint)
            & (predictions.population == "all")
            & (predictions.method == "B0")
        ]
        _, near_scores = near_duplicate_g(records, b0, target, endpoint)
        p = original[["p0", "p1", "p2"]].to_numpy(float)
        observed = original.observed.to_numpy(int)
        original_auc = auc_mann_whitney(
            original.observed.to_numpy(int)
            == original.displayed.to_numpy(int),
            original.p_displayed.to_numpy(float),
        )
        original_log_loss = float(
            -np.log(np.clip(p[np.arange(len(p)), observed], 1e-6, 1.0)).mean()
        )
        original_sigma = float(sigma.loc[(target, endpoint), "sigma_doc"])
        rows.append(
            {
                "target": target,
                "endpoint": endpoint,
                "report_priority": "focus" if target in focus_targets else "",
                "n_units_original": int(
                    len(document_values(records, merge_identical=False))
                ),
                "n_units_near_duplicate": near_scores["n_units"],
                "n_predictions_original": int(len(original)),
                "n_predictions_near_duplicate": near_scores["n_predictions"],
                "sigma_doc_original": original_sigma,
                "sigma_doc_near_duplicate": near_scores["sigma_doc"],
                "sigma_doc_difference": near_scores["sigma_doc"] - original_sigma,
                "auc_original": original_auc,
                "auc_near_duplicate": near_scores["auc"],
                "auc_difference": near_scores["auc"] - original_auc,
                "log_loss_original": original_log_loss,
                "log_loss_near_duplicate": near_scores["log_loss"],
                "log_loss_difference": near_scores["log_loss"]
                - original_log_loss,
                "sigma_fit_method": "full-data pooled between near-duplicate units",
                "analysis_role": ROLE,
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["report_priority", "endpoint", "target"], ascending=[False, True, True]
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    inputs = load_inputs()
    null, null_checks = flag_null_analysis(
        inputs["records"], inputs["census"], inputs["sigma"]
    )
    flag_clean = flag_clean_trainfold_analysis(inputs)
    near = near_duplicate_analysis(inputs)
    null.to_csv(OUT / "flag_null.csv", index=False)
    flag_clean.to_csv(OUT / "flag_clean_trainfold.csv", index=False)
    near.to_csv(OUT / "near_duplicate_merge.csv", index=False)
    focus = near.loc[near.report_priority == "focus"]
    print(
        "FLAG NULL",
        null_checks,
        "FLAG CLEAN datasets=",
        len(flag_clean),
        "NEAR DUPLICATES datasets=",
        len(near),
        flush=True,
    )
    print(focus.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
