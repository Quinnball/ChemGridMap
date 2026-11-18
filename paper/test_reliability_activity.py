"""Check offset recovery, held-out predictions, calibration and adoption rules."""
import numpy as np
import pandas as pd
import pytest

from reliability_activity import (
    adoption_decision, calibration_bins, corrected_values, excluded_prediction,
    fit_document_ridge, select_ridge_penalty, support_summary, unflagged_units,
)


def ridge_fixture():
    rows = []
    for i in range(30):
        for document, offset in (("A", -1.2), ("B", .3), ("C", .9)):
            rows.append({"entry": f"E{i:02d}", "document": document,
                         "documents": document, "document_value": 5. + i / 20 + offset})
    return pd.DataFrame(rows)


def test_ridge_recovers_known_offsets_with_unpenalized_entry_means():
    doc = ridge_fixture()
    fit = fit_document_ridge(doc, .1)
    assert fit["offsets"] == pytest.approx({"A": -1.2, "B": .3, "C": .9}, abs=.005)
    assert fit["entry_means"]["E00"] == pytest.approx(5.)
    assert fit["entry_means"]["E29"] == pytest.approx(6.45)
    assert sum(fit["offsets"].values()) == pytest.approx(0., abs=1e-10)


def test_internal_cv_selects_offset_recovery_and_accounts_for_every_unit():
    doc = ridge_fixture()
    fit, cv = select_ridge_penalty(doc)
    assert fit["penalty"] == .1
    assert cv.groupby("lambda").n_predictions.sum().eq(len(doc)).all()
    assert len(cv) == 5 * 7
    assert fit["variance_fit"]["tau"] < .01


def test_merged_offset_uses_representative_document_and_unseen_is_zero():
    doc = pd.DataFrame({"document": ["A", "unseen"], "documents": ["A;B", "unseen"],
                        "document_value": [6., 7.]})
    np.testing.assert_array_equal(corrected_values(doc, {"A": 1., "B": 3.}), [5., 7.])


def test_exclusion_checks_all_constituents_and_empty_rest_falls_back():
    doc = pd.DataFrame({"document": ["A"], "documents": ["A;B"], "document_value": [6.]})
    assert unflagged_units(doc, {"B"}).empty
    fallback = np.array([.2, .3, .5])
    actual, used, reverted = excluded_prediction(doc, {"B"}, {"tau": .5, "prior_df": np.inf}, fallback)
    np.testing.assert_array_equal(actual, fallback)
    assert used == 0 and reverted


def test_calibration_uses_displayed_probability_instead_of_maximum():
    frame = pd.DataFrame({"p0": [.1, .2, 0.], "p1": [.8, .3, 0.], "p2": [.1, .5, 1.],
                          "displayed": [0, 0, 2], "observed": [0, 1, 2]})
    bins = calibration_bins(frame)
    assert [row["n"] for row in bins] == [1, 1, 0, 0, 1]
    assert bins[0]["x"] == .1 and bins[0]["y"] == 1.
    assert bins[-1]["x"] == 1.
    assert not any(row["plot_eligible"] for row in bins)


def test_adoption_has_no_tolerance_and_never_selects_n():
    rows = [{"target": target, "method": method, "brier": .5, "ece": .1,
             "high_confidence_accuracy": .8 if method == "G" else .86}
            for target in ("all", "chembl205", "chembl204", "chembl240")
            for method in ("G", "O", "X", "N")]
    metrics = pd.DataFrame(rows)
    metrics.loc[(metrics.target == "all") & (metrics.method == "O"), "brier"] += 1e-12
    metrics.loc[(metrics.target == "all") & (metrics.method == "X"), "brier"] = .49
    decision, rules = adoption_decision(metrics)
    assert decision["default_model"] == "X"
    assert not rules.query("method == 'O' and rule == 'overall_brier'").iloc[0].passed
    assert decision["N_excluded_from_adoption"]


def test_support_summary_keeps_empty_strata_missing_and_threshold_boundaries():
    frame = pd.DataFrame({"target": ["chembl205"], "evidence_policy": ["merged"],
        "n_evidence_units": [1], "activity_class": ["active"], "class_confidence": [.8],
        "class_changed_by_document_balance": [0], "discordant_document_pairs": [0]})
    summary = support_summary(frame)
    empty = summary.query("scope == 'multi_unit'")
    assert empty.n.eq(0).all() and empty.supported_n.eq(0).all()
    assert empty.supported_fraction.isna().all()
    active = summary.query("target == 'chembl205' and scope == 'all' and display_class == 'active'")
    assert active.supported_n.tolist() == [1, 1, 0]
