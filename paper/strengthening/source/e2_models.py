"""E2a 多靶点整理、折外预测与现有做法比较对象。"""
from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[3]
S = ROOT / "paper/strengthening"
sys.path[:0] = [
    str(S / "source"),
    str(ROOT / "paper"),
    str(ROOT / "src"),
]

import numpy as np
import pandas as pd

from chemgridmap.activity import activity_class
from chemgridmap.chembl import curate_chembl_activity_data
from chemgridmap.evidence import (
    CLASSES,
    LITERATURE_BETWEEN_UNIT_SD,
    class_probabilities,
    document_pair_offsets,
    document_values,
    fit_between_unit_sd,
)
from e1a_models import b1_confidence, record_prediction
from e1a_statistics import auc_mann_whitney
from reliability_activity import (
    entry_folds,
    excluded_prediction,
    flagged_documents,
    gaussian_prediction,
    unflagged_units,
)

ENDPOINTS = ("IC50", "Ki")
MODELS = ("B0", "B1", "G", "G-record", "G-fixed", "G-global", "X")
POPULATIONS = ("all", "flag-clean")
STRATA = ("all", "m=1", "m>=2")
FOLD_SEED = 20261024
BOOTSTRAP_SEED = 20261028
CLASS_BOUNDARIES = (6.0, 7.0)
CURATION_DEFAULTS = {
    "activity_type": None,
    "lower_threshold": 6.0,
    "upper_threshold": 7.0,
    "conflict_range_threshold": 1.0,
    "exclude_potential_duplicates": True,
    "molecule_identity_requested": "auto",
    "structure_policy": "auto",
    "assay_types": None,
}


def curation_settings(target_id: str, endpoint: str) -> dict:
    """原三靶点沿用冻结报告；其余靶点使用封存计划指定的默认整理。"""
    frozen = ROOT / f"chemgridmap/paper/output/revision_v6/{target_id.lower()}"
    report_path = frozen / f"{target_id.lower()}_curation_report.json"
    if endpoint == "IC50" and report_path.exists():
        params = json.loads(report_path.read_text())["parameters"]
        return {
            "activity_type": params["activity_type"],
            "lower_threshold": params["lower_threshold"],
            "upper_threshold": params["upper_threshold"],
            "conflict_range_threshold": params["conflict_range_threshold"],
            "exclude_potential_duplicates": params["exclude_potential_duplicates"],
            "molecule_identity_requested": params["molecule_identity_requested"],
            "structure_policy": params["structure_policy"],
            "assay_types": params["assay_types"],
        }
    return {**CURATION_DEFAULTS, "activity_type": endpoint}


def curate_activity_table(raw: pd.DataFrame, target_id: str, endpoint: str):
    """API 列由 CGM 的规范化列名检测直接识别，不在此处改写字段时间格式。"""
    params = curation_settings(target_id, endpoint)
    return curate_chembl_activity_data(
        raw,
        target_id=target_id,
        activity_type=params["activity_type"],
        lower_threshold=params["lower_threshold"],
        upper_threshold=params["upper_threshold"],
        conflict_range_threshold=params["conflict_range_threshold"],
        exclude_potential_duplicates=params["exclude_potential_duplicates"],
        molecule_identity=params["molecule_identity_requested"],
        structure_policy=params["structure_policy"],
        assay_types=params["assay_types"],
    )


def clean_training(doc: pd.DataFrame, excluded_entries: set[str]) -> pd.DataFrame:
    """flag-clean 拟合时整条移除完整数据中涉及标记文献的条目。"""
    return doc.loc[~doc.entry.isin(excluded_entries)].copy()


def touch_flagged_entries(doc: pd.DataFrame, flagged: set[str]) -> set[str]:
    touched = doc.documents.map(
        lambda value: bool(set(str(value).split(";")) & flagged)
    )
    return set(doc.loc[touched, "entry"])


def document_pair_observations(records: pd.DataFrame) -> pd.DataFrame:
    """一次生成文献对共享条目观测，供每个训练折重新判定标记文献。"""
    doc = document_values(records, merge_identical=False)
    rows = []
    for entry, group in doc.groupby("entry", sort=True):
        values = group[["document", "document_value"]].itertuples(index=False, name=None)
        for (document_a, value_a), (document_b, value_b) in combinations(values, 2):
            rows.append((document_a, document_b, entry, float(value_b - value_a)))
    return pd.DataFrame(
        rows, columns=["document_a", "document_b", "entry", "difference"]
    )


def flags_from_observations(
    observations: pd.DataFrame, training_entries: set[str] | None = None
) -> set[str]:
    """按 CGM 默认阈值重新找出训练条目中的系统性偏移文献。"""
    selected = observations
    if training_entries is not None:
        selected = selected.loc[selected.entry.isin(training_entries)]
    flagged = set()
    for _, group in selected.groupby(["document_a", "document_b"], sort=True):
        if len(group) < 3:
            continue
        median = float(np.median(group.difference))
        sign = float(np.sign(median))
        agreement = float(np.mean(group.difference.to_numpy(float) * sign > 0))
        if abs(median) >= 1.0 and agreement >= 0.8:
            flagged.update((group.document_a.iloc[0], group.document_b.iloc[0]))
    return flagged


def assert_flag_helper_matches_cgm(records: pd.DataFrame) -> None:
    """启动时核对快速文献对路径与原 CGM 默认输出一致。"""
    expected = flagged_documents(records)
    actual = flags_from_observations(document_pair_observations(records))
    if actual != expected:
        raise RuntimeError("E2 文献对标记与 CGM 默认实现不一致。")


def pooled_global_fit(
    docs_by_target: dict[str, pd.DataFrame],
    excluded_entries: dict[str, set[str]],
    target: str,
    population: str,
) -> dict:
    """留一靶点合并其它纳入靶点的多单元条目，排除当前靶点全部数据。"""
    frames = []
    donors = []
    for donor, doc in docs_by_target.items():
        if donor == target:
            continue
        selected = (
            clean_training(doc, excluded_entries.get(donor, set()))
            if population == "flag-clean"
            else doc
        )
        if selected.empty:
            continue
        selected = selected.copy()
        selected["entry"] = donor + "::" + selected.entry.astype(str)
        frames.append(selected)
        donors.append(donor)
    if not frames:
        return {
            "tau": LITERATURE_BETWEEN_UNIT_SD,
            "prior_df": float("inf"),
            "model": "gaussian",
            "method": "no_other_included_target_fallback",
            "n_entries": 0,
            "n_units": 0,
            "donor_targets": [],
        }
    fit = fit_between_unit_sd(pd.concat(frames, ignore_index=True))
    fit["donor_targets"] = donors
    return fit


def build_global_fits(
    docs_by_target: dict[str, pd.DataFrame],
    excluded_entries: dict[str, set[str]],
) -> dict[str, dict[str, dict]]:
    fits = {}
    for target in sorted(docs_by_target):
        fits[target] = {
            population: pooled_global_fit(
                docs_by_target, excluded_entries, target, population
            )
            for population in POPULATIONS
        }
    return fits


def comparator_scores(records: pd.DataFrame, held_documents: str) -> tuple[float, float, int]:
    """只使用剔除留出文献后的原始记录计算 range 与 SD 分数。"""
    held = set(str(held_documents).split(";"))
    remaining = records.loc[
        ~records.document_chembl_id.astype(str).isin(held), "pchembl_value"
    ].to_numpy(float)
    if not len(remaining):
        raise ValueError("留出文献后没有剩余原始记录。")
    range_score = float(float(np.max(remaining) - np.min(remaining)) <= 1.0)
    sd_score = (
        float(-np.std(remaining, ddof=1))
        if len(remaining) >= 2
        else float("nan")
    )
    return range_score, sd_score, len(remaining)


def _fit_for_population(
    doc: pd.DataFrame,
    population: str,
    excluded_entries: set[str],
) -> pd.DataFrame:
    return clean_training(doc, excluded_entries) if population == "flag-clean" else doc


def _prediction_row(
    base: dict,
    method: str,
    probabilities,
    tau: float,
    used_units: int,
    fallback: bool,
) -> dict:
    probabilities = np.asarray(probabilities, dtype=float)
    return {
        **base,
        "method": method,
        "tau": float(tau),
        "used_units": int(used_units),
        "fallback_to_G": bool(fallback),
        "p0": float(probabilities[0]),
        "p1": float(probabilities[1]),
        "p2": float(probabilities[2]),
    }


def prepare_target(records: pd.DataFrame) -> dict:
    """缓存一个靶点不随折变化的文档单元、文献对与 flag-clean 条目。"""
    doc = document_values(records, merge_identical=True)
    raw_units = document_values(
        records, document="activity_id", merge_identical=False
    )
    observations = document_pair_observations(records)
    full_flags = flags_from_observations(observations)
    return {
        "doc": doc,
        "raw_units": raw_units,
        "observations": observations,
        "full_flags": full_flags,
        "clean_entries": touch_flagged_entries(doc, full_flags),
    }


def build_target_predictions(
    records: pd.DataFrame,
    target: str,
    endpoint: str,
    global_fits: dict[str, dict] | None = None,
    models: tuple[str, ...] = MODELS,
    context: dict | None = None,
) -> tuple[pd.DataFrame, set[str]]:
    """构造一个靶点的所有折外预测；总体长表稍后统一补齐 B1。"""
    context = context or prepare_target(records)
    doc = context["doc"]
    raw_units = context["raw_units"]
    observations = context["observations"]
    clean_entries = context["clean_entries"]
    assignments = entry_folds(doc, seed=FOLD_SEED)
    units_by_entry = dict(tuple(doc.groupby("entry", sort=True)))
    records_by_entry = dict(tuple(records.groupby("molecule_identity_key", sort=True)))
    rows = []

    for fold in range(10):
        test_entries = {entry for entry, assigned in assignments.items() if assigned == fold}
        train_entries = set(doc.entry.unique()) - test_entries
        g_train = doc.loc[~doc.entry.isin(test_entries)]
        record_train = raw_units.loc[~raw_units.entry.isin(test_entries)]
        fold_flags = flags_from_observations(observations, train_entries)
        x_train = unflagged_units(g_train, fold_flags)
        g_fit = fit_between_unit_sd(g_train)
        record_fit = fit_between_unit_sd(record_train)
        x_fit = fit_between_unit_sd(x_train)
        clean_g_train = clean_training(g_train, clean_entries)
        clean_record_train = clean_training(record_train, clean_entries)
        clean_x_train = clean_training(x_train, clean_entries)
        clean_g_fit = fit_between_unit_sd(clean_g_train)
        clean_record_fit = fit_between_unit_sd(clean_record_train)
        clean_x_fit = fit_between_unit_sd(clean_x_train)
        fixed_fit = {**g_fit, "tau": LITERATURE_BETWEEN_UNIT_SD}
        global_fit_all = (global_fits or {}).get(target, {}).get("all", fixed_fit)
        global_fit_clean = (
            (global_fits or {}).get(target, {}).get("flag-clean", fixed_fit)
        )

        for entry in sorted(test_entries):
            units = units_by_entry[entry]
            source_records = records_by_entry[entry]
            for held in units.itertuples(index=False):
                held_documents = str(held.documents)
                held_set = set(held_documents.split(";"))
                remaining_units = units.loc[
                    ~units.documents.map(
                        lambda value: bool(set(str(value).split(";")) & held_set)
                    )
                ]
                remaining_raw = source_records.loc[
                    ~source_records.document_chembl_id.astype(str).isin(held_set),
                    "pchembl_value",
                ].to_numpy(float)
                if remaining_units.empty or not len(remaining_raw):
                    raise ValueError("留出一个文献单元后没有剩余证据。")
                rest = remaining_units.document_value.to_numpy(float)
                observed = CLASSES.index(
                    activity_class(float(held.document_value), *CLASS_BOUNDARIES)
                )
                displayed = CLASSES.index(
                    activity_class(float(np.median(remaining_raw)), *CLASS_BOUNDARIES)
                )
                range_score, sd_score, remaining_n = comparator_scores(
                    source_records, held_documents
                )
                base = {
                    "target": target,
                    "endpoint": endpoint,
                    "fold": int(fold),
                    "entry": entry,
                    "held_out_unit": held_documents,
                    "m": int(len(rest)),
                    "displayed": displayed,
                    "observed": observed,
                    "remaining_records": remaining_n,
                    "b_range": range_score,
                    "b_sd": sd_score,
                }
                all_probabilities = {
                    "B0": np.eye(3)[displayed],
                    "G": gaussian_prediction(rest, g_fit),
                    "G-record": record_prediction(
                        source_records, held_documents, record_fit["tau"]
                    )[0],
                    "G-fixed": gaussian_prediction(rest, fixed_fit),
                    "G-global": gaussian_prediction(rest, global_fit_all),
                }
                all_probabilities["X"], x_used, x_fallback = excluded_prediction(
                    remaining_units, fold_flags, x_fit, all_probabilities["G"]
                )
                clean_available = entry not in clean_entries
                clean_probabilities = {
                    "B0": np.eye(3)[displayed],
                    "G": gaussian_prediction(rest, clean_g_fit),
                    "G-record": record_prediction(
                        source_records, held_documents, clean_record_fit["tau"]
                    )[0],
                    "G-fixed": gaussian_prediction(rest, fixed_fit),
                    "G-global": gaussian_prediction(rest, global_fit_clean),
                }
                clean_probabilities["X"], clean_x_used, clean_x_fallback = excluded_prediction(
                    remaining_units,
                    fold_flags,
                    clean_x_fit,
                    clean_probabilities["G"],
                )
                for method in models:
                    if method == "B1":
                        continue
                    if method not in all_probabilities:
                        raise ValueError(f"未知模型：{method}")
                    tau_all = {
                        "G": g_fit["tau"],
                        "G-record": record_fit["tau"],
                        "G-fixed": fixed_fit["tau"],
                        "G-global": global_fit_all["tau"],
                        "X": g_fit["tau"],
                    }.get(method, np.nan)
                    rows.append(
                        _prediction_row(
                            {**base, "population": "all", "flag_clean": clean_available},
                            method,
                            all_probabilities[method],
                            tau_all,
                            x_used if method == "X" else len(rest),
                            x_fallback if method == "X" else False,
                        )
                    )
                    if not clean_available:
                        continue
                    tau_clean = {
                        "G": clean_g_fit["tau"],
                        "G-record": clean_record_fit["tau"],
                        "G-fixed": fixed_fit["tau"],
                        "G-global": global_fit_clean["tau"],
                        "X": clean_g_fit["tau"],
                    }.get(method, np.nan)
                    rows.append(
                        _prediction_row(
                            {**base, "population": "flag-clean", "flag_clean": True},
                            method,
                            clean_probabilities[method],
                            tau_clean,
                            clean_x_used if method == "X" else len(rest),
                            clean_x_fallback if method == "X" else False,
                        )
                    )

    output = pd.DataFrame(rows)
    if output.empty:
        return output, clean_entries

    # B1 的同靶点折外一致率来自 B0 全部测试预测，不先按条目平均。
    if "B1" in models:
        b0 = output.loc[
            (output.population == "all") & (output.method == "B0")
        ].copy()
        b1_rows = []
        for fold in range(10):
            if not (b0.fold == fold).any():
                continue
            confidence = b1_confidence(b0, target, fold)
            selected = b0.loc[b0.fold == fold].copy()
            for row in selected.to_dict("records"):
                probabilities = np.full(3, (1.0 - confidence) / 2.0)
                probabilities[int(row["displayed"])] = confidence
                b1_rows.append(
                    _prediction_row(
                        {key: value for key, value in row.items()
                         if key not in {"method", "tau", "used_units", "fallback_to_G",
                                        "p0", "p1", "p2"}},
                        "B1",
                        probabilities,
                        np.nan,
                        row["m"],
                        False,
                    )
                )
                if row["flag_clean"]:
                    b1_rows.append(
                        _prediction_row(
                            {
                                **{
                                    key: value for key, value in row.items()
                                    if key not in {"method", "tau", "used_units",
                                                   "fallback_to_G", "p0", "p1", "p2"}
                                },
                                "population": "flag-clean",
                            },
                            "B1",
                            probabilities,
                            np.nan,
                            row["m"],
                            False,
                        )
                    )
        output = pd.concat([output, pd.DataFrame(b1_rows)], ignore_index=True)
    output["p_displayed"] = output[
        ["p0", "p1", "p2"]
    ].to_numpy(float)[np.arange(len(output)), output.displayed.to_numpy(int)]
    output["agree"] = output.observed == output.displayed
    output["stratum"] = np.where(output.m == 1, "m=1", "m>=2")
    probabilities = output[["p0", "p1", "p2"]].to_numpy(float)
    if not np.allclose(probabilities.sum(axis=1), 1.0, atol=1e-12, rtol=0):
        raise ValueError("E2 模型概率之和不为 1。")
    if "G-global" in models:
        expected = pd.MultiIndex.from_frame(
            output.loc[
                (output.population == "all") & (output.method == "B0"),
                ["target", "endpoint", "fold", "entry", "held_out_unit", "m"],
            ]
        ).sort_values()
        for method in models:
            actual = pd.MultiIndex.from_frame(
                output.loc[
                    (output.population == "all") & (output.method == method),
                    ["target", "endpoint", "fold", "entry", "held_out_unit", "m"],
                ]
            ).sort_values()
            if not actual.equals(expected):
                raise ValueError(f"{method} 的留出键与 B0 不一致。")
    return (
        output.sort_values(
            ["target", "endpoint", "population", "method", "fold", "entry",
             "held_out_unit"]
        ).reset_index(drop=True),
        clean_entries,
    )


def comparator_aucs(predictions: pd.DataFrame) -> pd.DataFrame:
    """按靶点、总体汇总 B-range、B-sd 与相同子集上的 G AUC。"""
    rows = []
    for (target, endpoint, population), frame in predictions.groupby(
        ["target", "endpoint", "population"], sort=True
    ):
        b0 = frame.loc[frame.method == "B0"]
        g = frame.loc[frame.method == "G"]
        keys = ["fold", "entry", "held_out_unit"]
        g_sd = g.drop(columns=["remaining_records"]).merge(
            b0[keys + ["remaining_records"]], on=keys, validate="one_to_one"
        )
        g_sd = g_sd.loc[g_sd.remaining_records >= 2]
        b_sd = b0.loc[b0.remaining_records >= 2]
        for model, frame_for_n, labels, scores in (
            ("B-range", b0, b0.agree.to_numpy(bool), b0.b_range.to_numpy(float)),
            ("B-sd", b_sd, b_sd.agree.to_numpy(bool), b_sd.b_sd.to_numpy(float)),
            ("G", g, g.agree.to_numpy(bool), g.p_displayed.to_numpy(float)),
            (
                "G (B-sd subset)",
                g_sd,
                g_sd.agree.to_numpy(bool),
                g_sd.p_displayed.to_numpy(float),
            ),
        ):
            rows.append({
                "target": target,
                "endpoint": endpoint,
                "population": population,
                "model": model,
                "n_predictions": int(len(frame_for_n)),
                "n_entries": int(frame_for_n.entry.nunique()),
                "auc": auc_mann_whitney(labels, scores),
            })
    return pd.DataFrame(rows)
