"""E1-A 的统计等价性、无折内泄漏与 R1 负控约束。"""
import numpy as np
import pandas as pd
import pytest
from scipy.stats import norm

from e1a_evidence_audit import document_values, fit_between_unit_sd, probability_scores
from e1a_statistics import auc_mann_whitney, cluster_weights, weighted_scores, capped_tau


def example_predictions():
    return pd.DataFrame({
        "target": ["a"] * 4 + ["b"] * 2,
        "entry": ["same", "same", "other", "third", "same", "same"],
        "observed": [0, 1, 2, 0, 0, 2], "displayed": [0, 1, 1, 0, 0, 2],
        "p0": [.8, .1, .1, .3, 1., .1], "p1": [.1, .8, .8, .4, 0., .1],
        "p2": [.1, .1, .1, .3, 0., .8],
    })


@pytest.mark.parametrize("labels,scores,expected", [
    ([0, 1, 0, 1], [.1, .2, .3, .4], .75),
    ([0, 1, 0, 1], [.1, .2, .2, .4], .875),
    ([0, 1, 0, 1], [.5, .5, .5, .5], .5),
])
def test_auc_hand_calculation_including_ties(labels, scores, expected):
    assert auc_mann_whitney(labels, scores) == expected


def test_auc_one_class_is_na():
    assert np.isnan(auc_mann_whitney([1, 1], [.1, .8]))


def test_weighted_metrics_match_literal_repetition_and_frozen_ece():
    frame = example_predictions()
    weights = np.array([[1, 1, 1, 1, 1, 1], [2, 2, 0, 1, 3, 3], [1, 0, 2, 3, 0, 0]])
    actual = weighted_scores(frame, weights)
    for i, counts in enumerate(weights):
        sampled = frame.iloc[np.repeat(np.arange(len(frame)), counts)]
        frozen = probability_scores(sampled)
        for name in ["accuracy", "brier", "log_loss", "ece", "mean_confidence",
                     "high_confidence_n", "high_confidence_mean_probability", "high_confidence_accuracy"]:
            assert actual[name][i] == pytest.approx(frozen[name], abs=1e-14)
        shown = sampled[["p0", "p1", "p2"]].to_numpy()[np.arange(len(sampled)), sampled.displayed.to_numpy()]
        labels = sampled.observed.to_numpy() == sampled.displayed.to_numpy()
        np.testing.assert_allclose(actual["auc"][i], auc_mann_whitney(labels, shown), equal_nan=True)
        assert actual["calibration_bias"][i] == pytest.approx(shown.mean() - labels.mean())


def test_bootstrap_reproducible_and_target_entry_clusters_preserved():
    frame = example_predictions()
    first = cluster_weights(frame, np.random.default_rng(20261028), 50)
    second = cluster_weights(frame, np.random.default_rng(20261028), 50)
    np.testing.assert_array_equal(first, second)
    np.testing.assert_array_equal(first[:, 0], first[:, 1])
    np.testing.assert_array_equal(first[:, 4], first[:, 5])
    np.testing.assert_array_equal(first[:, [0, 2, 3]].sum(axis=1), np.full(50, 3))
    np.testing.assert_array_equal(first[:, 4], np.ones(50))
    assert np.any(first[:, 0] != first[:, 4])  # 跨靶点同名条目不共用重数。


def test_g_record_tau_manual_example():
    # 五条目各三条记录 [0,0,3]：SSE=6、df=2，所以 pooled SD=sqrt(3)。
    records = pd.DataFrame([{"molecule_identity_key": f"e{e}", "activity_id": e * 3 + r,
        "document_chembl_id": "same_document", "pchembl_value": value}
        for e in range(5) for r, value in enumerate([0., 0., 3.])])
    units = document_values(records, document="activity_id", merge_identical=False)
    assert len(units) == 15
    assert fit_between_unit_sd(units)["tau"] == pytest.approx(np.sqrt(3))


def test_capped_df_formula():
    # a 的 df=1、s2=2；b 的 df=5、s2=3.5。封顶后权重是 1 与 4。
    doc = pd.DataFrame({"entry": ["a"] * 2 + ["b"] * 6,
                        "document_value": [0., 2., 0., 1., 2., 3., 4., 5.]})
    assert capped_tau(doc, cap=4) == pytest.approx(np.sqrt((2 + 4 * 3.5) / 5))


def test_b1_excludes_current_fold_and_other_targets():
    from e1a_models import b1_confidence

    # 训练条目的预测数不同，必须按全部预测求一致率，不能先算条目均值。
    b0 = pd.DataFrame({
        "target": ["a"] * 6 + ["b"] * 3,
        "entry": ["held", "held", "train1", "train1", "train1", "train2",
                  "other", "other", "other"],
        "fold": [0, 0, 1, 1, 1, 2, 1, 1, 2],
        "observed": [0, 1, 1, 0, 0, 1, 0, 0, 0],
        "displayed": [0] * 9,
    })
    assert b1_confidence(b0, "a", 0) == .5
    changed = b0.copy()
    changed.loc[(changed.target == "a") & (changed.fold == 0), "observed"] = 2
    changed.loc[changed.target == "b", "observed"] = 2
    assert b1_confidence(changed, "a", 0) == .5
    changed.loc[(changed.target == "a") & (changed.fold == 1), "observed"] = 0
    assert b1_confidence(changed, "a", 0) == .75


def test_g_record_excludes_every_record_of_merged_held_documents():
    from e1a_models import record_prediction

    # a、b 的中位数同为 3，合并留出时两篇文献的全部五条记录都必须排除。
    source = pd.DataFrame({
        "document_chembl_id": ["a", "a", "b", "b", "b", "c", "c", "d"],
        "pchembl_value": [2., 4., 1., 3., 5., 7., 9., 10.],
    })
    actual, remaining_n = record_prediction(source, "a;b", .5)
    center, sd = 26 / 3, .5 * np.sqrt(1 + 1 / 3)
    low, high_cdf = norm.cdf([(6 - center) / sd, (7 - center) / sd])
    np.testing.assert_allclose(actual, [low, high_cdf - low, 1 - high_cdf], atol=1e-15)
    assert remaining_n == 3
    changed = source.copy()
    changed.loc[changed.document_chembl_id.isin(["a", "b"]), "pchembl_value"] = 50.
    changed_p, changed_n = record_prediction(changed, "a;b", .5)
    np.testing.assert_array_equal(changed_p, actual)
    assert changed_n == remaining_n


def test_r1_permutation_preserves_target_marginals_and_all_probabilities():
    from e1a_evidence_audit import TARGETS, permute_agreement

    rows = []
    for target_i, target in enumerate(TARGETS):
        for i in range(24):
            rows.append({"target": target, "method": "G", "entry": f"e{i // 2}",
                         "observed": int(i % (target_i + 2) == 0), "displayed": 0,
                         "p0": .25 + i / 100, "p1": .5 - i / 100, "p2": .25})
    original = pd.DataFrame(rows)
    original = pd.concat([original, original.assign(method="B0")], ignore_index=True)
    before = original.copy(deep=True)
    null = permute_agreement(original)
    pd.testing.assert_frame_equal(original, before)
    expected_g = original[original.method == "G"].reset_index(drop=True)
    pd.testing.assert_frame_equal(null[expected_g.columns].reset_index(drop=True), expected_g)
    pd.testing.assert_frame_equal(null, permute_agreement(original))
    assert set(null.agree.unique()) == {False, True}
    for target in TARGETS:
        selected = null.target == target
        assert null.loc[selected, "agree"].sum() == (
            expected_g.loc[selected, "observed"] == expected_g.loc[selected, "displayed"]
        ).sum()
    assert np.any(null.agree.to_numpy() != (null.observed == null.displayed).to_numpy())


def test_agree_override_only_changes_binary_agreement_metrics():
    frame = example_predictions()
    weights = np.array([[1, 1, 1, 1, 1, 1], [2, 2, 0, 1, 3, 3]])
    original = weighted_scores(frame, weights)
    frame["agree"] = [False, True, True, False, True, False]
    actual = weighted_scores(frame, weights)
    for name in original.keys() - {"auc", "agreement_rate", "calibration_bias"}:
        np.testing.assert_allclose(actual[name], original[name], equal_nan=True)
    for i, counts in enumerate(weights):
        sampled = frame.iloc[np.repeat(np.arange(len(frame)), counts)]
        shown = sampled[["p0", "p1", "p2"]].to_numpy()[
            np.arange(len(sampled)), sampled.displayed.to_numpy()]
        labels = sampled.agree.to_numpy(bool)
        # 直接枚举正负配对作独立核验，避免与秩和实现共用计算路径。
        pair_scores = [float(a > b) + .5 * (a == b)
                       for a in shown[labels] for b in shown[~labels]]
        assert actual["auc"][i] == pytest.approx(np.mean(pair_scores))
        assert actual["agreement_rate"][i] == pytest.approx(labels.mean())
        assert actual["calibration_bias"][i] == pytest.approx(shown.mean() - labels.mean())
