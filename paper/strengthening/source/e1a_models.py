"""E1-A 的六个模型与事后 flag-clean 总体；冻结预测只复制、不改写。"""
from __future__ import annotations

import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / "paper"), str(ROOT / "src")]

import numpy as np
import pandas as pd

from chemgridmap.activity import activity_class
from chemgridmap.evidence import (
    CLASSES, LITERATURE_BETWEEN_UNIT_SD, class_probabilities,
    document_pair_offsets, document_values, fit_between_unit_sd,
)
from reliability_activity import entry_folds, flagged_documents, unflagged_units

MODELS = ("B0", "B1", "G", "G-record", "G-fixed", "X")
ROW_KEY = ["target", "fold", "entry", "held_out_unit"]
NOTE = "post hoc sensitivity analysis motivated by the flagged hERG pair"
FOLD_SEED = 20261024


def b1_confidence(b0, target, fold):
    """同靶点折外条目的全部 B0 预测一致率，不按条目再平均。"""
    selected = b0[b0.target == target]
    if "method" in selected:
        selected = selected[selected.method == "B0"]
    if "evaluation_set" in selected:
        selected = selected[selected.evaluation_set == "merged"]
    held_entries = set(selected.loc[selected.fold == fold, "entry"])
    train = selected[~selected.entry.isin(held_entries)]
    if train.empty:
        raise ValueError(f"B1 没有折外预测：target={target}, fold={fold}")
    return float((train.observed == train.displayed).mean())


def record_prediction(source_records, held_documents, tau):
    """去除留出单元的全部文献记录，以剩余原始记录预测文献单元。"""
    held = set(str(held_documents).split(";"))
    remaining = source_records.loc[
        ~source_records.document_chembl_id.astype(str).isin(held), "pchembl_value"
    ].to_numpy(float)
    if not len(remaining):
        raise ValueError("G-record 留出文献后没有剩余记录")
    scale = float(tau) * np.sqrt(1.0 + 1.0 / len(remaining))
    probabilities = class_probabilities([remaining.mean()], [scale], [np.inf], 6.0, 7.0)[0]
    return probabilities, len(remaining)


def _fit_row(target, fold, method, fit, train, test):
    """保留实际拟合样本与固定折，供后续输出和审计复核。"""
    return {
        "target": target, "fold": fold, "method": method,
        "tau": fit["tau"], "prior_df": fit["prior_df"],
        "fit_method": fit["method"], "n_entries": fit["n_entries"],
        "n_units": fit["n_units"],
        "degrees_of_freedom": fit.get("degrees_of_freedom", np.nan),
        "training_entries": int(train.entry.nunique()),
        "test_entries": len(test),
        "test_entry_ids": json.dumps(sorted(test), ensure_ascii=False),
    }


def build_predictions(records, frozen):
    """返回无 population 重复的预测、标记条目集合、折拟合表及剔除统计。"""
    frozen = frozen[frozen.evaluation_set == "merged"].copy()
    base = frozen[frozen.method.isin(["B0", "G", "X"])].copy()
    if base.duplicated(["method", *ROW_KEY]).any():
        raise ValueError("冻结预测的模型与留出键重复")
    b0 = base[base.method == "B0"]
    expected_keys = pd.MultiIndex.from_frame(b0[ROW_KEY]).sort_values()
    for method in ("G", "X"):
        actual_keys = pd.MultiIndex.from_frame(base.loc[base.method == method, ROW_KEY]).sort_values()
        if not actual_keys.equals(expected_keys):
            raise ValueError(f"冻结 {method} 与 B0 的留出键不一致")

    additions, fit_rows, flag_rows = [], [], []
    flagged_entries, row_metadata, fold_tau = {}, {}, {}
    for target, source in records.items():
        doc = document_values(source)
        rec_doc = document_values(source, document="activity_id", merge_identical=False)
        assignments = entry_folds(doc, seed=FOLD_SEED)
        source_by_entry = dict(tuple(source.groupby("molecule_identity_key")))
        units_by_entry = dict(tuple(doc.groupby("entry")))

        # flag-clean 用完整数据识别文献对，并剔除触及任一构成文献的整个条目。
        pairs = document_pair_offsets(source)
        pairs = pairs[pairs.systematic_offset.astype(bool)]
        flagged = set(pairs.document_a) | set(pairs.document_b)
        touched = doc.documents.map(lambda text: bool(set(text.split(";")) & flagged))
        excluded = set(doc.loc[touched, "entry"])
        flagged_entries[target] = excluded
        target_b0 = b0[b0.target == target]
        removed = target_b0.entry.isin(excluded)
        flag_rows.append({
            "target": target, "n_flagged_pairs": len(pairs),
            "n_flagged_documents": len(flagged), "n_flagged_entries": len(excluded),
            "n_all_entries": int(doc.entry.nunique()),
            "n_prediction_entries": int(target_b0.entry.nunique()),
            "n_excluded_prediction_entries": int(target_b0.loc[removed, "entry"].nunique()),
            "n_predictions_all": len(target_b0), "n_predictions_excluded": int(removed.sum()),
            "n_predictions_flag_clean": int((~removed).sum()), "sensitivity_note": NOTE,
        })

        for fold in range(10):
            test = {entry for entry, assigned in assignments.items() if assigned == fold}
            train = doc[~doc.entry.isin(test)]
            record_train = rec_doc[~rec_doc.entry.isin(test)]
            g_fit = fit_between_unit_sd(train)
            record_fit = fit_between_unit_sd(record_train)
            x_flags = flagged_documents(source[~source.molecule_identity_key.isin(test)])
            x_train = unflagged_units(train, x_flags)
            x_fit = fit_between_unit_sd(x_train)
            fixed_fit = {**g_fit, "tau": LITERATURE_BETWEEN_UNIT_SD,
                         "method": "fixed_literature_between_unit_sd"}
            for method, fit, fitting_units in (
                ("G", g_fit, train), ("G-record", record_fit, record_train),
                ("G-fixed", fixed_fit, train), ("X", x_fit, x_train),
            ):
                fit_rows.append(_fit_row(target, fold, method, fit, fitting_units, test))
                fold_tau[(target, fold, method)] = fit["tau"]

            fold_b0 = target_b0[target_b0.fold == fold]
            if fold_b0.empty:
                continue
            confidence = b1_confidence(b0, target, fold)
            fit_rows.append({
                "target": target, "fold": fold, "method": "B1", "tau": np.nan,
                "prior_df": np.nan, "fit_method": "out_of_fold_B0_agreement",
                "b1_confidence": confidence,
                "training_entries": int(target_b0.loc[~target_b0.entry.isin(test), "entry"].nunique()),
                "test_entries": len(test),
                "test_entry_ids": json.dumps(sorted(test), ensure_ascii=False),
                "n_training_predictions": int((~target_b0.entry.isin(test)).sum()),
            })
            for held in fold_b0.to_dict("records"):
                entry = held["entry"]
                if assignments[entry] != fold:
                    raise ValueError("冻结预测与 seed 20261024 的文献条目折不一致")
                units = units_by_entry[entry]
                remaining_units = units[units.documents != held["held_out_unit"]]
                if len(remaining_units) != int(held["remaining_units"]):
                    raise ValueError("冻结留出单元与当前文献合并结果不一致")
                source_records = source_by_entry[entry]
                held_documents = set(held["held_out_unit"].split(";"))
                remaining_records = source_records.loc[
                    ~source_records.document_chembl_id.astype(str).isin(held_documents), "pchembl_value"
                ].to_numpy(float)
                displayed = CLASSES.index(activity_class(float(np.median(remaining_records)), 6.0, 7.0))
                if displayed != int(held["displayed"]):
                    raise ValueError("显示类别不等于剩余原始记录中位数的类别")
                row_metadata[(target, fold, entry, held["held_out_unit"])] = len(remaining_records)

                b1_p = np.full(3, (1.0 - confidence) / 2.0)
                b1_p[displayed] = confidence
                record_p, record_n = record_prediction(source_records, held["held_out_unit"], record_fit["tau"])
                rest_values = remaining_units.document_value.to_numpy(float)
                fixed_p = class_probabilities(
                    [rest_values.mean()],
                    [LITERATURE_BETWEEN_UNIT_SD * np.sqrt(1.0 + 1.0 / len(rest_values))],
                    [np.inf], 6.0, 7.0,
                )[0]
                for method, probabilities, tau, used in (
                    ("B1", b1_p, np.nan, len(rest_values)),
                    ("G-record", record_p, record_fit["tau"], record_n),
                    ("G-fixed", fixed_p, LITERATURE_BETWEEN_UNIT_SD, len(rest_values)),
                ):
                    additions.append({
                        **held, "method": method, "tau": tau, "used_units": used,
                        "fallback_to_G": False,
                        **{f"p{k}": float(probabilities[k]) for k in range(3)},
                    })

    # B0/G/X 的冻结字段原样沿用；tau 只添加溯源信息，X 回退时标注实际 G 的 tau。
    base["tau"] = [
        np.nan if row.method == "B0" else fold_tau[
            (row.target, row.fold, "G" if row.method == "X" and row.fallback_to_G else row.method)
        ]
        for row in base.itertuples(index=False)
    ]
    predictions = pd.concat([base, pd.DataFrame(additions)], ignore_index=True)
    predictions["remaining_records"] = [
        row_metadata[(row.target, row.fold, row.entry, row.held_out_unit)]
        for row in predictions.itertuples(index=False)
    ]
    probabilities = predictions[["p0", "p1", "p2"]].to_numpy(float)
    predictions["p_displayed"] = probabilities[
        np.arange(len(predictions)), predictions.displayed.to_numpy(int)
    ]
    predictions["agree"] = predictions.observed == predictions.displayed
    predictions["flag_clean"] = [
        row.entry not in flagged_entries[row.target] for row in predictions.itertuples(index=False)
    ]
    predictions["stratum"] = np.where(predictions.remaining_units == 1, "m=1", "m>=2")
    predictions["sensitivity_note"] = NOTE
    predictions = predictions.sort_values(["target", "method", "fold", "entry", "held_out_unit"]).reset_index(drop=True)
    if not np.allclose(probabilities.sum(axis=1), 1.0, atol=1e-12, rtol=0):
        raise ValueError("模型概率之和不为 1")
    for method in MODELS:
        keys = pd.MultiIndex.from_frame(predictions.loc[predictions.method == method, ROW_KEY]).sort_values()
        if not keys.equals(expected_keys):
            raise ValueError(f"{method} 与冻结 B0 的留出键不一致")
    return predictions, flagged_entries, pd.DataFrame(fit_rows), pd.DataFrame(flag_rows)
