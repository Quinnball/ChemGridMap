"""E1-A：按任务单固定规则审计三个靶点；冻结模块仅导入。"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

# 禁止导入时在冻结目录产生字节码。
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[3]
CGM = ROOT
REV = ROOT / "paper"
S = REV / "strengthening"
OUT = S / "output/E1A"
QA = S / "qa/E1A"
sys.path[:0] = [str(REV), str(ROOT / "src")]

import numpy as np
import pandas as pd
from scipy import stats
from reliability_activity import (entry_folds, gaussian_prediction, lodo_models,
                                  probability_scores, evaluate_predictions)
from chemgridmap.activity import activity_class
from chemgridmap.evidence import (CLASSES, LITERATURE_BETWEEN_UNIT_SD,
    document_values, fit_between_unit_sd, class_probabilities, document_pair_offsets)
from e1a_statistics import (auc_mann_whitney, cluster_weights, weighted_scores,
                            percentile_ci, capped_tau)

TARGETS = ("chembl205", "chembl204", "chembl240")
MODELS = ("B0", "B1", "G", "G-record", "G-fixed", "X")
KEY = ["target", "method", "fold", "entry", "held_out_unit"]
NOTE = "post hoc sensitivity analysis motivated by the flagged hERG pair"
BOOT_SEED, BOOT_N, FOLD_SEED, NULL_SEED = 20261028, 2000, 20261024, 20261030


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


def load_inputs():
    records = {t: pd.read_csv(CGM / f"paper/output/revision_v6/{t}/{t}_retained_records.csv")
               for t in TARGETS}
    frozen = pd.read_csv(CGM / "paper/output/reliability_v1/activity/lodo_predictions.csv")
    return records, frozen[frozen.evaluation_set == "merged"].copy()


def regression_check(records, frozen):
    """首先逐条重算；任何不一致都抛出异常，不执行后续分析。"""
    regenerated = []
    for target in TARGETS:
        predictions, _, _, _ = lodo_models(target, records[target])
        regenerated.append(predictions[predictions.method.isin(["B0", "G"])])
    actual = pd.concat(regenerated).set_index(KEY).sort_index()
    expected = frozen[frozen.method.isin(["B0", "G"])].set_index(KEY).sort_index()
    assert actual.index.is_unique and expected.index.is_unique
    assert actual.index.equals(expected.index), "冻结预测键不一致"
    errors = np.abs(actual[["p0", "p1", "p2"]].to_numpy() - expected[["p0", "p1", "p2"]].to_numpy())
    assert errors.max() <= 1e-12, f"概率回归失败：{errors.max()}"
    for name in ["observed", "displayed"]:
        assert np.array_equal(actual[name], expected[name]), f"{name} 回归失败"
    g = actual.reset_index().query("method == 'G'")
    score = probability_scores(g)
    assert score["n_predictions"] == 820 and score["n_entries"] == 289
    assert abs(score["ece"] - 0.0767537203890331) <= 1e-12
    assert abs(score["brier"] - 0.5218617813580732) <= 1e-12
    result = {"passed": True, "aligned_rows": len(actual), "max_probability_error": float(errors.max()),
              "G": score}
    write_json(QA / "regression.json", result)
    print("REGRESSION PASS", json.dumps(result), flush=True)
    return result


def permute_agreement(frozen):
    """R1：在靶点的全部 G 预测内置换一致指标，概率与类别保持不变。"""
    null = frozen[frozen.method == "G"].copy().reset_index(drop=True)
    null["agree"] = null.observed == null.displayed
    shuffle = np.random.default_rng(NULL_SEED)
    for target in TARGETS:
        selected = null.target == target
        null.loc[selected, "agree"] = shuffle.permutation(null.loc[selected, "agree"].to_numpy())
    return null


def negative_control(frozen):
    """负控仅用于 QA；无论结果如何，都不阻止后续分析。"""
    path, archive = QA / "negative_control.csv", QA / "negative_control_R0_misspecified.csv"
    if path.exists() and not archive.exists():
        old = pd.read_csv(path)
        assert old.permutation.str.contains("observed class labels").all(), "旧负控内容与 R0 不符"
        path.rename(archive)
    null = permute_agreement(frozen)
    rng = np.random.default_rng(BOOT_SEED)
    rows = []
    for target in ("all", *TARGETS):
        frame = null if target == "all" else null[null.target == target]
        draws = weighted_scores(frame, cluster_weights(frame, rng, BOOT_N))["auc"]
        low, high, valid = percentile_ci(draws)
        shown = frame[["p0", "p1", "p2"]].to_numpy()[np.arange(len(frame)), frame.displayed.to_numpy(int)]
        row = {"target": target, "auc": auc_mann_whitney(frame.agree, shown),
               "ci_low": low, "ci_high": high, "n_bootstrap_valid": valid,
               "contains_half": bool(low <= .5 <= high), "shuffle_seed": NULL_SEED,
               "bootstrap_seed": BOOT_SEED, "n_bootstrap": BOOT_N,
               "permutation": "binary agreement indicator permuted within each target",
               "sensitivity_note": NOTE}
        rows.append(row)
    result = pd.DataFrame(rows)
    result.to_csv(QA / "negative_control.csv", index=False)
    print("NEGATIVE CONTROL\n" + result.to_string(index=False), flush=True)
    return result


def evaluate_all(predictions):
    """每个评价总体单独抽样；同一抽样矩阵用于全部模型及成对差异。"""
    metrics, intervals, paired = [], [], []
    rng = np.random.default_rng(BOOT_SEED)
    align_key = ["target", "fold", "entry", "held_out_unit"]
    for population in ("all", "flag-clean"):
        pop = predictions if population == "all" else predictions[predictions.flag_clean]
        for target in ("all", *TARGETS):
            target_frame = pop if target == "all" else pop[pop.target == target]
            for stratum in ("all", "m=1", "m>=2"):
                frame = target_frame if stratum == "all" else target_frame[target_frame.stratum == stratum]
                base = {"population": population, "target": target, "stratum": stratum,
                        "sensitivity_note": NOTE}
                reference = frame[frame.method == "B0"].sort_values(align_key).reset_index(drop=True)
                assert len(reference), f"空评价组合：{base}"
                weights = cluster_weights(reference, rng, BOOT_N)
                point, draws = {}, {}
                for method in MODELS:
                    selected = frame[frame.method == method].sort_values(align_key).reset_index(drop=True)
                    assert selected[align_key].equals(reference[align_key])
                    score = probability_scores(selected)
                    score["n_entries"] = len(selected[["target", "entry"]].drop_duplicates())
                    extra = weighted_scores(selected, np.ones((1, len(selected))))
                    for name in ("auc", "agreement_rate", "mean_p_displayed", "calibration_bias"):
                        score[name] = float(extra[name][0])
                    metrics.append({**base, "method": method, **score})
                    point[method] = score
                    draws[method] = weighted_scores(selected, weights)
                    for metric, samples in draws[method].items():
                        low, high, valid = percentile_ci(samples)
                        intervals.append({**base, "method": method, "metric": metric,
                            "estimate": score[metric], "ci_low": low, "ci_high": high,
                            "n_bootstrap_valid": valid, "n_bootstrap": BOOT_N, "seed": BOOT_SEED})
                for other in (m for m in MODELS if m != "G"):
                    for metric in ("brier", "log_loss", "ece", "auc"):
                        low, high, valid = percentile_ci(draws["G"][metric] - draws[other][metric])
                        paired.append({**base, "reference": "G", "alternative": other,
                            "metric": metric, "difference": point["G"][metric] - point[other][metric],
                            "ci_low": low, "ci_high": high, "n_bootstrap_valid": valid,
                            "n_bootstrap": BOOT_N, "seed": BOOT_SEED})
                print(f"METRICS {population} {target} {stratum}: {len(reference)} predictions", flush=True)
    return pd.DataFrame(metrics), pd.DataFrame(intervals), pd.DataFrame(paired)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--regression-only", action="store_true")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    records, frozen = load_inputs()
    regression_check(records, frozen)
    if not args.regression_only:
        from e1a_models import build_predictions
        from e1a_descriptive import run_descriptive
        from e1a_reporting import write_results, write_manifest
        negative = negative_control(frozen)
        predictions, flagged_entries, fits, flags = build_predictions(records, frozen)
        fits.assign(sensitivity_note=NOTE).to_csv(QA / "fold_fits.csv", index=False)
        flags.assign(sensitivity_note=NOTE).to_csv(QA / "flag_clean_selection.csv", index=False)
        # population 为长表：flag-clean 是 all 的子集，不是额外的独立预测。
        marked = pd.concat([predictions.assign(population="all"),
                            predictions[predictions.flag_clean].assign(population="flag-clean")], ignore_index=True)
        marked["sensitivity_note"] = NOTE
        marked.to_csv(OUT / "predictions.csv", index=False)
        metrics, intervals, paired = evaluate_all(predictions)
        for name, frame in [("metrics", metrics), ("bootstrap_ci", intervals), ("paired_differences", paired)]:
            frame.to_csv(OUT / f"{name}.csv", index=False)
        mapped = pd.read_csv(CGM / "paper/output/reliability_v1/activity/mapped_entry_confidence.csv")
        descriptive = run_descriptive(records, predictions, mapped, flagged_entries, OUT, QA)
        write_results(metrics, intervals, paired, descriptive, flags)
        write_manifest(negative)
        print("E1-A COMPLETE; negative-control QA passed=" + str(bool(negative.contains_half.all())), flush=True)


if __name__ == "__main__":
    main()
