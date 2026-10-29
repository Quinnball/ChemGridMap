"""E1-A 3.5–3.9：保持预定定义的描述分析，不修改冻结输入。"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from reliability_activity import entry_folds
from chemgridmap.activity import activity_class
from chemgridmap.evidence import (CLASSES, LITERATURE_BETWEEN_UNIT_SD,
    class_probabilities, document_pair_offsets, document_values, fit_between_unit_sd)
from e1a_statistics import capped_tau, percentile_ci


TARGETS = ("chembl205", "chembl204", "chembl240")
NOTE = "post hoc sensitivity analysis motivated by the flagged hERG pair"
BOOT_SEED, BOOT_N, FOLD_SEED = 20261028, 2000, 20261024
LOWER, UPPER, SUPPORT, DF_CAP = 6.0, 7.0, 0.8, 4


def _shown_probability(center, scale, shown):
    p = class_probabilities(np.atleast_1d(center), scale, np.inf, LOWER, UPPER)
    shown = np.broadcast_to(np.asarray(shown, dtype=int), len(p))
    return p[np.arange(len(p)), shown]


def _tau_bootstrap(doc, rng):
    """所有条目（含单单元）均被抽样；重复抽中的条目分别贡献自由度。"""
    moments = doc.groupby("entry", sort=True).document_value.agg(["size", "var"])
    degrees = moments["size"].to_numpy(int) - 1
    sumsquares = degrees * moments["var"].fillna(0).to_numpy(float)
    repeated = (degrees > 0).astype(int)
    result = np.empty(BOOT_N)
    # 分批仅限制内存；随机整数流与一次生成全部抽样一致。
    for start in range(0, BOOT_N, 100):
        end = min(start + 100, BOOT_N)
        draws = rng.integers(0, len(moments), size=(end - start, len(moments)))
        denominator = degrees[draws].sum(axis=1)
        numerator = sumsquares[draws].sum(axis=1)
        multi_count = repeated[draws].sum(axis=1)
        pooled = np.sqrt(np.divide(numerator, denominator,
            out=np.zeros(end - start), where=denominator > 0))
        # 与冻结 fit_between_unit_sd 默认的至少五个重复条目规则一致。
        result[start:end] = np.where(multi_count >= 5, pooled, LITERATURE_BETWEEN_UNIT_SD)
    return result


def tau_influence(records, documents, flagged_entries):
    rows, loo_all, draws_all = [], [], []
    rng = np.random.default_rng(BOOT_SEED)
    for target in TARGETS:
        doc = documents[target]
        fitted = fit_between_unit_sd(doc)
        clean = doc[~doc.entry.isin(flagged_entries[target])]
        clean_fit = fit_between_unit_sd(clean)
        moments = doc.groupby("entry", sort=True).document_value.agg(["size", "var"])
        multi = moments[moments["size"] >= 2]
        degrees = int((multi["size"] - 1).sum())
        draws = _tau_bootstrap(doc, rng)
        low, high, valid = percentile_ci(draws)
        rows.append({"row_type": "summary", "target": target,
            "n_entries": len(moments), "n_multi_unit_entries": len(multi),
            "n_units": len(doc), "degrees_of_freedom": degrees,
            "tau_full": fitted["tau"], "tau_flag_clean": clean_fit["tau"],
            "flag_clean_n_entries": clean.entry.nunique(),
            "excluded_n_entries": len(flagged_entries[target]),
            "tau_cap": capped_tau(doc, DF_CAP), "df_cap": DF_CAP,
            "tau_ci_low": low, "tau_ci_high": high, "n_bootstrap_valid": valid,
            "n_bootstrap": BOOT_N, "bootstrap_seed": BOOT_SEED,
            "bootstrap_population": "all retained entries, including single-unit entries",
            "tau_fit_method": fitted["method"], "flag_clean_fit_method": clean_fit["method"]})
        local = []
        for entry, item in multi.iterrows():
            without = fit_between_unit_sd(doc[doc.entry != entry])
            local.append({"row_type": "leave_one_entry", "target": target, "entry": entry,
                "n_units": int(item["size"]), "entry_df": int(item["size"] - 1),
                "degrees_of_freedom_fraction": (item["size"] - 1) / degrees,
                "tau_full": fitted["tau"], "tau_without_entry": without["tau"],
                "delta_tau": without["tau"] - fitted["tau"],
                "is_flagged_entry": entry in flagged_entries[target]})
        ranked = sorted(local, key=lambda item: (-abs(item["delta_tau"]), item["entry"]))
        for rank, item in enumerate(ranked, 1):
            item["absolute_delta_rank"] = rank
        loo_all.extend(ranked)
        rows.extend(ranked[:5])
        draws_all.extend({"target": target, "replicate": i + 1, "tau": value,
                          "bootstrap_seed": BOOT_SEED} for i, value in enumerate(draws))
    return pd.DataFrame(rows), pd.DataFrame(loo_all), pd.DataFrame(draws_all)


def pseudo_consensus(documents):
    units = []
    for target in TARGETS:
        doc = documents[target]
        assignments = entry_folds(doc, seed=FOLD_SEED)
        folds = {}
        for fold in sorted(set(assignments.values())):
            test = {entry for entry, assigned in assignments.items() if assigned == fold}
            folds[fold] = fit_between_unit_sd(doc[~doc.entry.isin(test)])["tau"]
        for entry, group in doc.groupby("entry", sort=True):
            m = len(group)
            if m < 5:
                continue
            values = group.document_value.to_numpy(float)
            remaining_means = (values.sum() - values) / (m - 1)
            shown = np.array([CLASSES.index(activity_class(value, LOWER, UPPER)) for value in values])
            observed = np.array([CLASSES.index(activity_class(value, LOWER, UPPER))
                                 for value in remaining_means])
            tau = folds[assignments[entry]]
            p_cons = _shown_probability(values, tau, shown)
            p_pred = _shown_probability(values, tau * np.sqrt(1 + 1 / (m - 1)), shown)
            for i, unit in enumerate(group.itertuples(index=False)):
                units.append({"row_type": "unit", "target": target, "entry": entry,
                    "unit": unit.documents, "entry_units": m, "fold": assignments[entry],
                    "training_tau": tau, "unit_value": values[i],
                    "remaining_mean": remaining_means[i], "displayed": int(shown[i]),
                    "displayed_class": CLASSES[shown[i]], "observed": int(observed[i]),
                    "observed_class": CLASSES[observed[i]], "p_cons": p_cons[i],
                    "p_pred": p_pred[i], "agree": bool(shown[i] == observed[i])})
    detail = pd.DataFrame(units)
    summaries = []
    for target in ("all", *TARGETS):
        selected = detail if target == "all" else detail[detail.target == target]
        summaries.append({"row_type": "summary", "target": target, "n_units": len(selected),
            "n_entries": len(selected[["target", "entry"]].drop_duplicates()),
            "mean_p_cons": selected.p_cons.mean(), "mean_p_pred": selected.p_pred.mean(),
            "n_agree": int(selected.agree.sum()), "agreement_rate": selected.agree.mean(),
            "analysis_scope": "descriptive only; entries with at least five merged units"})
    return pd.concat([pd.DataFrame(summaries), detail], ignore_index=True)


def count_flips(records, documents, mapped):
    rows = []
    for target in TARGETS:
        source = records[target]
        rec_units = document_values(source, document="activity_id", merge_identical=False)
        tau_rec = fit_between_unit_sd(rec_units)["tau"]
        tau = fit_between_unit_sd(documents[target])["tau"]
        moments = source.groupby("molecule_identity_key").pchembl_value.agg(["size", "mean", "median"])
        unit_n = documents[target].groupby("entry").size()
        doc_counts = source.groupby(["molecule_identity_key", "document_chembl_id"]).size()
        max_doc = doc_counts.groupby(level=0).max()
        selected = mapped[mapped.target == target]
        for item in selected.itertuples(index=False):
            entry = item.molecule_identity_key
            moment = moments.loc[entry]
            shown = CLASSES.index(item.activity_class)
            p_document = float(_shown_probability([item.document_balanced_pchembl],
                tau / np.sqrt(unit_n.loc[entry]), shown)[0])
            assert abs(p_document - item.class_confidence) <= 1e-12, "完整 G 支持度与冻结注释不一致"
            p_record = float(_shown_probability([moment["mean"]],
                tau_rec / np.sqrt(moment["size"]), shown)[0])
            supported_document = bool(item.class_confidence >= SUPPORT)
            supported_record = bool(p_record >= SUPPORT)
            rows.append({"target": target, "entry": entry,
                "n_records": int(moment["size"]), "n_units": int(unit_n.loc[entry]),
                "max_document_record_fraction": max_doc.loc[entry] / moment["size"],
                "document_support": item.class_confidence, "record_support": p_record,
                "recomputed_document_support": p_document,
                "document_supported": supported_document, "record_supported": supported_record,
                "flip": supported_document != supported_record,
                "flip_direction": ("uncertain_to_supported" if supported_record else "supported_to_uncertain")
                    if supported_document != supported_record else "unchanged",
                "displayed": shown, "displayed_class": item.activity_class,
                "median": item.activity_pchembl, "retained_record_median": moment["median"],
                "record_mean": moment["mean"], "document_mean": item.document_balanced_pchembl,
                "tau_document": tau, "tau_record": tau_rec})
    detail = pd.DataFrame(rows)
    summaries = []
    for target in ("all", *TARGETS):
        selected = detail if target == "all" else detail[detail.target == target]
        summaries.append({"target": target, "n_mapped_entries": len(selected),
            "n_flip": int(selected.flip.sum()),
            "n_uncertain_to_supported": int((selected.flip_direction == "uncertain_to_supported").sum()),
            "n_supported_to_uncertain": int((selected.flip_direction == "supported_to_uncertain").sum()),
            "n_document_supported": int(selected.document_supported.sum()),
            "n_record_supported": int(selected.record_supported.sum()),
            "support_threshold": SUPPORT})
    return detail[detail.flip].copy(), pd.DataFrame(summaries), detail


def flagged_pairs_profile(records, predictions, flagged_entries):
    path = Path(__file__).resolve().parents[1] / "data/flagged_pairs_metadata.json"
    metadata = json.loads(path.read_text())
    g = predictions[predictions.method == "G"].copy()
    assert not g.duplicated(["target", "fold", "entry", "held_out_unit"]).any()
    g["agree"] = g.observed == g.displayed
    rows, shared_union = [], {target: set() for target in TARGETS}
    expected = {(pair["target"], pair["document_a"], pair["document_b"]) for pair in metadata["pairs"]}
    actual = set()
    for target in TARGETS:
        offsets = document_pair_offsets(records[target])
        for pair in offsets[offsets.systematic_offset].itertuples(index=False):
            actual.add((target, pair.document_a, pair.document_b))
    assert expected == actual, "元数据标记文献对与默认诊断不一致"
    for pair in metadata["pairs"]:
        target, a, b = pair["target"], pair["document_a"], pair["document_b"]
        source = records[target]
        doc = document_values(source, merge_identical=False)
        wide = doc[doc.document.isin([a, b])].pivot(index="entry", columns="document", values="document_value")
        shared = wide[[a, b]].dropna()
        shared_union[target].update(shared.index)
        differences = shared[b] - shared[a]
        correlation_defined = bool(shared[a].nunique() > 1 and shared[b].nunique() > 1)
        rho = float(stats.spearmanr(shared[a], shared[b]).statistic) if correlation_defined else np.nan
        selected = g[(g.target == target) & g.entry.isin(shared.index)]
        base = {"target": target, "document_a": a, "document_b": b,
            "n_shared_entries": len(shared), "median_offset_b_minus_a": float(differences.median()),
            "spearman_rho": rho, "spearman_defined": correlation_defined,
            "spearman_note": "" if correlation_defined else "NA: a document has constant shared-entry values",
            "counted_prediction_model": "G",
            "metadata_source": metadata["retrieved_from"], "metadata_retrieved_at": metadata["retrieved_at"]}
        for suffix, document in [("a", a), ("b", b)]:
            info = metadata["documents"][document]
            base.update({f"{key}_{suffix}": info.get(key) for key in
                ("doc_type", "year", "patent_id", "journal", "title", "doi", "pubmed_id")})
        assert len(shared) == pair["shared"]
        assert abs(float(differences.median()) - pair["median_offset"]) < 1e-12
        rows.append({"row_type": "summary", **base, "scope": "shared entries in this pair",
            "n_entries": len(shared), "n_prediction_entries": selected.entry.nunique(),
            "n_predictions": len(selected), "n_agree": int(selected.agree.sum()),
            "agreement_rate": selected.agree.mean()})
        touched = set(source.loc[source.document_chembl_id.isin([a, b]), "molecule_identity_key"])
        touched_predictions = g[(g.target == target) & g.entry.isin(touched)]
        rows.append({"row_type": "summary", **base,
            "scope": "all entries touched by either document", "n_entries": len(touched),
            "n_prediction_entries": touched_predictions.entry.nunique(),
            "n_predictions": len(touched_predictions), "n_agree": int(touched_predictions.agree.sum()),
            "agreement_rate": touched_predictions.agree.mean()})
        grouped_raw = source[source.document_chembl_id.isin([a, b])].groupby(
            ["molecule_identity_key", "document_chembl_id"]).pchembl_value
        raw_values = grouped_raw.agg(lambda values: json.dumps(sorted(values.astype(float).tolist())))
        for entry, values in shared.iterrows():
            entry_predictions = selected[selected.entry == entry]
            rows.append({"row_type": "shared_entry", **base, "entry": entry,
                "value_a": values[a], "value_b": values[b], "offset_b_minus_a": values[b] - values[a],
                "class_a": activity_class(values[a], LOWER, UPPER),
                "class_b": activity_class(values[b], LOWER, UPPER),
                "retained_values_a": raw_values.loc[(entry, a)],
                "retained_values_b": raw_values.loc[(entry, b)],
                "n_predictions": len(entry_predictions), "n_agree": int(entry_predictions.agree.sum())})
    # 并集单独去重；同一条目若属于两对文献，不在这里重复计算。
    for target in TARGETS:
        for scope, entries in [("union of shared entries", shared_union[target]),
                               ("all entries touched by flagged documents", flagged_entries[target])]:
            selected = g[(g.target == target) & g.entry.isin(entries)]
            rows.append({"row_type": "union_summary", "target": target, "scope": scope,
                "n_entries": len(entries), "n_prediction_entries": selected.entry.nunique(),
                "counted_prediction_model": "G", "n_predictions": len(selected),
                "n_agree": int(selected.agree.sum()), "agreement_rate": selected.agree.mean()})
    return pd.DataFrame(rows)


def single_unit_thresholds(documents, mapped):
    rows, checks = [], []
    z = float(stats.norm.ppf(SUPPORT))
    for target in TARGETS:
        tau = fit_between_unit_sd(documents[target])["tau"]
        high, low = UPPER + z * tau, LOWER - z * tau
        medium_max = float(2 * stats.norm.cdf((UPPER - LOWER) / (2 * tau)) - 1)
        # 本任务三个数据集均满足此条件，因而中间类别无 0.8 支持区间。
        assert medium_max < SUPPORT, "出现中间类别可达 0.8 的 τ，需报告任务定义范围"
        selected = mapped[(mapped.target == target) & (mapped.n_evidence_units == 1)]
        target_checks = []
        for item in selected.itertuples(index=False):
            center = item.document_balanced_pchembl
            interval_supported = bool((item.activity_class == "active" and center >= high)
                or (item.activity_class == "inactive" and center <= low))
            frozen_supported = bool(item.class_confidence >= SUPPORT)
            p = float(_shown_probability([center], tau, CLASSES.index(item.activity_class))[0])
            target_checks.append({"target": target, "entry": item.molecule_identity_key,
                "displayed_class": item.activity_class, "displayed_median": item.activity_pchembl,
                "single_unit_center": center, "frozen_support": item.class_confidence,
                "recomputed_support": p, "interval_supported": interval_supported,
                "frozen_supported": frozen_supported, "status_matches": interval_supported == frozen_supported,
                "absolute_probability_error": abs(p - item.class_confidence)})
        checks.extend(target_checks)
        check = pd.DataFrame(target_checks)
        rows.append({"row_type": "summary", "target": target, "tau": tau, "z": z,
            "support_threshold": SUPPORT, "active_min_pchembl": high, "inactive_max_pchembl": low,
            "maximum_medium_support": medium_max, "medium_supported_interval": "none",
            "n_single_unit_entries": len(check), "n_supported": int(check.frozen_supported.sum()),
            "n_status_mismatches": int((~check.status_matches).sum()),
            "max_absolute_probability_error": check.absolute_probability_error.max(),
            "threshold_center": "document_balanced_pchembl (the single merged unit value)"})
    return pd.DataFrame(rows), pd.DataFrame(checks)


def run_descriptive(records, predictions, mapped, flagged_entries, out, qa):
    """计算并写出预定描述输出；输入预测不能有 population 复制行。"""
    out, qa = Path(out), Path(qa)
    mapped = mapped[(mapped.model == "G") & (mapped.evidence_policy == "merged")].copy()
    assert not mapped.duplicated(["target", "molecule_identity_key"]).any()
    documents = {target: document_values(records[target]) for target in TARGETS}
    tau, loo, draws = tau_influence(records, documents, flagged_entries)
    flips, flip_summary, all_counts = count_flips(records, documents, mapped)
    thresholds, threshold_checks = single_unit_thresholds(documents, mapped)
    outputs = {"tau_influence": tau, "pseudo_consensus": pseudo_consensus(documents),
        "count_flip_entries": flips, "count_flip_summary": flip_summary,
        "flagged_pairs_profile": flagged_pairs_profile(records, predictions, flagged_entries),
        "single_unit_thresholds": thresholds}
    quality = {"tau_leave_one_entry_all": loo, "tau_bootstrap_draws": draws,
        "count_support_all_mapped": all_counts, "single_unit_threshold_checks": threshold_checks,
        "single_unit_threshold_summary": thresholds.copy()}
    for directory, tables in [(out, outputs), (qa, quality)]:
        directory.mkdir(parents=True, exist_ok=True)
        for name, frame in tables.items():
            frame["sensitivity_note"] = NOTE
            frame.to_csv(directory / f"{name}.csv", index=False)
    # 先保存核对明细，再严格验证，便于失败时追溯实际差异。
    assert threshold_checks.status_matches.all(), "单单元支持度阈值与冻结映射不一致"
    assert threshold_checks.absolute_probability_error.max() <= 1e-12
    return {**outputs, **quality}
