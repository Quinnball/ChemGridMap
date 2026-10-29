"""只读独立核对 E1-A 产物；展开整群样本复核合计 CI，不重拟合模型。"""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / "paper"), str(ROOT / "src")]

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from reliability_activity import probability_scores

S = ROOT / "paper/strengthening"
OUT, QA = S / "output/E1A", S / "qa/E1A"
MODELS = ("B0", "B1", "G", "G-record", "G-fixed", "X")
TARGETS = ("all", "chembl205", "chembl204", "chembl240")
POPULATIONS, STRATA = ("all", "flag-clean"), ("all", "m=1", "m>=2")
KEY = ["target", "fold", "entry", "held_out_unit"]
CELL = ["population", "target", "stratum", "method"]
NOTE = "post hoc sensitivity analysis motivated by the flagged hERG pair"
REQUIRED = (
    "predictions.csv", "metrics.csv", "bootstrap_ci.csv", "paired_differences.csv",
    "tau_influence.csv", "pseudo_consensus.csv", "count_flip_entries.csv",
    "count_flip_summary.csv", "flagged_pairs_profile.csv", "single_unit_thresholds.csv",
    "RESULTS_E1A.md", "manifest.json",
)


def same(actual, expected):
    np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-12, equal_nan=True)


def interval_summary(values):
    finite = np.asarray(values)[np.isfinite(values)]
    return [*np.percentile(finite, [2.5, 97.5]), len(finite)] if len(finite) else [np.nan, np.nan, 0]


def literal_auc(labels, scores):
    labels = np.asarray(labels, bool)
    n1, n0 = int(labels.sum()), int((~labels).sum())
    if not n1 or not n0:
        return np.nan
    ranks = rankdata(scores, method="average")
    return (ranks[labels].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def scores_from_arrays(p, observed, displayed):
    """对展开后的行直接取均值，不调用 e1a_statistics 的重数实现。"""
    confidence, predicted = p.max(axis=1), p.argmax(axis=1)
    correct, agree = predicted == observed, displayed == observed
    shown = p[np.arange(len(p)), displayed]
    bins = np.minimum((confidence * 10).astype(int), 9)
    high = confidence >= .8
    return {
        "accuracy": correct.mean(),
        "brier": ((p - np.eye(3)[observed]) ** 2).sum(axis=1).mean(),
        "log_loss": -np.log(np.clip(p[np.arange(len(p)), observed], 1e-6, 1)).mean(),
        "ece": sum(abs((correct - confidence)[bins == b].sum()) for b in range(10)) / len(p),
        "mean_confidence": confidence.mean(), "agreement_rate": agree.mean(),
        "mean_p_displayed": shown.mean(), "calibration_bias": shown.mean() - agree.mean(),
        "high_confidence_n": high.sum(),
        "high_confidence_mean_probability": confidence[high].mean() if high.any() else np.nan,
        "high_confidence_accuracy": correct[high].mean() if high.any() else np.nan,
        "auc": literal_auc(agree, shown),
    }


def literal_cluster_samples(frame, rng, collect):
    """按靶点抽取原有条目数；直接拼接被抽到的整条目原始行。"""
    samples = [[] for _ in range(2000)] if collect else None
    for target in sorted(frame.target.unique()):
        groups = [group.index.to_numpy() for _, group in
                  frame[frame.target == target].groupby("entry", sort=True)]
        draws = rng.integers(0, len(groups), size=(2000, len(groups)))
        if collect:
            for b, draw in enumerate(draws):
                samples[b].append(np.concatenate([groups[i] for i in draw]))
    return [np.concatenate(parts) for parts in samples] if collect else None


def verify_negative_control():
    frozen = pd.read_csv(ROOT / "paper/output/reliability_v1/activity/lodo_predictions.csv")
    g = frozen[(frozen.evaluation_set == "merged") & (frozen.method == "G")].reset_index(drop=True)
    g["agree"] = g.observed.eq(g.displayed)
    g["score"] = g[["p0", "p1", "p2"]].to_numpy()[np.arange(len(g)), g.displayed.to_numpy(int)]
    original = g.copy()
    shuffle = np.random.default_rng(20261030)
    for target in TARGETS[1:]:
        selected = g.target.eq(target)
        g.loc[selected, "agree"] = shuffle.permutation(g.loc[selected, "agree"].to_numpy())
    rng = np.random.default_rng(20261028)
    qa = pd.read_csv(QA / "negative_control.csv", float_precision="round_trip").set_index("target")
    assert qa.permutation.str.contains("binary agreement").all()
    for target in TARGETS:
        frame = (g if target == "all" else g[g.target == target]).reset_index(drop=True)
        labels, probabilities = frame.agree.to_numpy(), frame.score.to_numpy()
        samples = literal_cluster_samples(frame, rng, True)
        aucs = np.array([literal_auc(labels[indices], probabilities[indices]) for indices in samples])
        finite = aucs[np.isfinite(aucs)]
        same([literal_auc(labels, probabilities), *np.percentile(finite, [2.5, 97.5]), len(finite)],
             qa.loc[target, ["auc", "ci_low", "ci_high", "n_bootstrap_valid"]].to_numpy(float))
    # H(a,b)=1[a>b]+0.5*1[a=b]。靶点内贡献 n1_t*n0_t/2；跨靶点贡献 q_t*(1-q_u)*sum H。
    numerator = 0.
    for t, a in original.groupby("target"):
        n1, n0 = int(a.agree.sum()), int((~a.agree).sum())
        numerator += .5 * n1 * n0
        for u, b in original.groupby("target"):
            if t == u:
                continue
            differences = a.score.to_numpy()[:, None] - b.score.to_numpy()[None, :]
            numerator += a.agree.mean() * (1 - b.agree.mean()) * (
                (differences > 0).sum() + .5 * (differences == 0).sum())
    expected = numerator / (original.agree.sum() * (~original.agree).sum())
    print(f"PASS: 8000 literal negative-control bootstrap draws; pooled null expectation={expected:.17g}")


def main():
    assert all((OUT / name).is_file() for name in REQUIRED)
    # AUC 的同分判断需要保留输出浮点数的全部位数，默认解析会改变少量末位。
    tables = {name: pd.read_csv(OUT / name, float_precision="round_trip")
              for name in REQUIRED if name.endswith(".csv")}
    for name, frame in tables.items():
        assert frame.sensitivity_note.eq(NOTE).all(), name
    predictions, metrics = tables["predictions.csv"], tables["metrics.csv"]
    intervals, paired = tables["bootstrap_ci.csv"], tables["paired_differences.csv"]
    assert len(metrics) == 144 and len(intervals) == 1728 and len(paired) == 480
    expected_cells = set(itertools.product(POPULATIONS, TARGETS, STRATA, MODELS))
    assert set(metrics[CELL].itertuples(index=False, name=None)) == expected_cells
    assert not metrics.duplicated(CELL).any()
    assert not predictions.duplicated(["population", "method", *KEY]).any()
    all_rows = predictions[predictions.population == "all"].reset_index(drop=True)
    clean = predictions[predictions.population == "flag-clean"].reset_index(drop=True)
    expected_clean = all_rows[all_rows.flag_clean].assign(population="flag-clean").reset_index(drop=True)
    pd.testing.assert_frame_equal(clean, expected_clean)
    assert set(predictions.population) == set(POPULATIONS)
    assert all_rows.groupby("method").size().eq(820).all()
    assert clean.groupby("method").size().eq(457).all()
    assert predictions.stratum.eq(np.where(predictions.remaining_units == 1, "m=1", "m>=2")).all()
    p = predictions[["p0", "p1", "p2"]].to_numpy()
    same(p.sum(axis=1), 1.)
    same(predictions.p_displayed, p[np.arange(len(p)), predictions.displayed.to_numpy(int)])
    assert predictions.agree.eq(predictions.observed.eq(predictions.displayed)).all()
    frozen = pd.read_csv(ROOT / "paper/output/reliability_v1/activity/lodo_predictions.csv")
    frozen = frozen[(frozen.evaluation_set == "merged") & frozen.method.isin(["B0", "G", "X"])]
    observed = all_rows[all_rows.method.isin(["B0", "G", "X"])]
    order = ["method", *KEY]
    pd.testing.assert_frame_equal(observed[frozen.columns].sort_values(order).reset_index(drop=True),
                                  frozen.sort_values(order).reset_index(drop=True), rtol=0, atol=1e-12)
    scores = metrics.set_index(CELL)
    interval_index = intervals.set_index([*CELL, "metric"])
    pair_index = paired.set_index(["population", "target", "stratum", "alternative", "metric"])
    assert interval_index.index.is_unique and pair_index.index.is_unique
    assert intervals.n_bootstrap.eq(2000).all() and paired.n_bootstrap.eq(2000).all()
    assert intervals.seed.eq(20261028).all() and paired.seed.eq(20261028).all()
    for row in metrics.itertuples(index=False):
        frame = predictions[(predictions.population == row.population) & (predictions.method == row.method)]
        if row.target != "all":
            frame = frame[frame.target == row.target]
        if row.stratum != "all":
            frame = frame[frame.stratum == row.stratum]
        expected = probability_scores(frame)
        expected["n_entries"] = len(frame[["target", "entry"]].drop_duplicates())
        expected.update(scores_from_arrays(frame[["p0", "p1", "p2"]].to_numpy(),
                                          frame.observed.to_numpy(int), frame.displayed.to_numpy(int)))
        for metric, value in expected.items():
            same(getattr(row, metric), value)
            if metric not in ("n_predictions", "n_entries"):
                same(interval_index.loc[(row.population, row.target, row.stratum, row.method, metric), "estimate"], value)
    for row in paired.itertuples(index=False):
        cell = (row.population, row.target, row.stratum)
        same(row.difference, scores.loc[(*cell, "G"), row.metric] - scores.loc[(*cell, row.alternative), row.metric])
        assert row.reference == "G" and row.alternative != "G"
    print("PASS: 12 required files; 144 complete metric cells; 1728 CIs; 480 paired differences")
    print(f"PASS: population subset identity, no duplicate keys; all={len(all_rows)}, flag-clean={len(clean)} rows")
    print("PASS: all 144 point estimates match frozen definitions; B0/G/X frozen fields unchanged")

    rng = np.random.default_rng(20261028)
    verified_draws = 0
    for population, target, stratum in itertools.product(POPULATIONS, TARGETS, STRATA):
        frame = predictions[predictions.population == population]
        if target != "all":
            frame = frame[frame.target == target]
        if stratum != "all":
            frame = frame[frame.stratum == stratum]
        reference = frame[frame.method == "B0"].sort_values(KEY).reset_index(drop=True)
        samples = literal_cluster_samples(reference, rng, target == "all")
        model_frames = {method: frame[frame.method == method].sort_values(KEY).reset_index(drop=True)
                        for method in MODELS}
        for selected in model_frames.values():
            pd.testing.assert_frame_equal(selected[KEY], reference[KEY])
        if target != "all":
            continue
        distributions = {}
        for method, selected in model_frames.items():
            probabilities = selected[["p0", "p1", "p2"]].to_numpy()
            observed, displayed = selected.observed.to_numpy(int), selected.displayed.to_numpy(int)
            distributions[method] = pd.DataFrame([
                scores_from_arrays(probabilities[index], observed[index], displayed[index]) for index in samples])
            for metric, values in distributions[method].items():
                ci = interval_index.loc[(population, target, stratum, method, metric)]
                same(interval_summary(values),
                     ci[["ci_low", "ci_high", "n_bootstrap_valid"]].to_numpy(float))
        for alternative in (method for method in MODELS if method != "G"):
            for metric in ("brier", "log_loss", "ece", "auc"):
                differences = distributions["G"][metric] - distributions[alternative][metric]
                ci = pair_index.loc[(population, target, stratum, alternative, metric)]
                same(interval_summary(differences),
                     ci[["ci_low", "ci_high", "n_bootstrap_valid"]].to_numpy(float))
        verified_draws += len(samples) * len(MODELS)
        print(f"PASS: literal shared bootstrap, {population}/{target}/{stratum}", flush=True)
    print(f"PASS: {verified_draws} expanded model-bootstrap samples; all pooled CIs and paired CIs")
    verify_negative_control()

    manifest = json.loads((OUT / "manifest.json").read_text())
    for section in ("inputs_sha256", "source_sha256", "outputs_sha256", "qa_sha256"):
        for relative, expected in manifest[section].items():
            assert "strengthening/data/chembl37_api" not in relative
            assert hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() == expected, relative
    listed_outputs = {Path(relative).name for relative in manifest["outputs_sha256"]}
    assert set(REQUIRED) - {"manifest.json"} <= listed_outputs
    report = (OUT / "RESULTS_E1A.md").read_text()
    assert all(f"## 3.{i} " in report for i in range(10))
    print("PASS: all declared manifest hashes, output coverage and RESULTS sections 3.0-3.9")


if __name__ == "__main__":
    main()
