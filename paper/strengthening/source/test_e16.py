import pandas as pd

import e16_identity_sensitivity as e16


def test_variant_kind_classification():
    assert e16.classify_variant_kind(
        "[H][C@](F)(Cl)Br", "F[C@@H](Cl)Br", ()
    ) == "explicit_H"
    assert e16.classify_variant_kind(
        "O=c1cccc[nH]1", "Oc1ccccn1", ()
    ) == "tautomer"
    assert e16.classify_variant_kind(
        "C[N+](=O)[O-]", "C[N+]([O-])=O", (-1, 1)
    ) == "resonance"


def test_combined_report_keeps_original_and_sensitivity_rows():
    columns = [
        "hypothesis", "endpoint", "population", "n_targets",
        "median_difference", "median_ci_low", "median_ci_high",
        "wilcoxon_statistic", "raw_p", "holm_p",
    ]
    original = pd.DataFrame([{
        **dict.fromkeys(columns, 0),
        "hypothesis": "H3",
        "endpoint": "IC50",
        "population": "all",
    }])
    updated = original.copy()
    updated["median_difference"] = 0.25
    combined = e16.combine_report_rows(original, updated, ("H3",))
    assert combined.scenario.tolist() == [
        "original", "with_reviewed_variants"
    ]
    assert combined.median_difference.tolist() == [0, 0.25]


def test_h1_comparison_uses_summary_columns():
    row = {
        "hypothesis": "H1", "endpoint": "IC50", "population": "all",
        "n_targets": 29, "median": 0.2, "q1": 0.1, "q3": 0.3,
        "pooled_estimate": 0.25, "pooled_ci_low": 0.2,
        "pooled_ci_high": 0.3,
    }
    combined = e16.combine_h1_rows(
        pd.DataFrame([row]), pd.DataFrame([row])
    )
    assert combined.scenario.tolist() == [
        "original", "with_reviewed_variants"
    ]
    assert combined["median"].tolist() == [0.2, 0.2]


def test_e16_model_set_matches_the_e2_hypothesis_inputs():
    assert e16.MODELS_E16 == (
        "B0", "G", "G-record", "G-fixed", "G-global", "X"
    )
