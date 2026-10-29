"""E2a 的 G-global、比较对象、flag-clean 与 E1-A 实现回归测试。"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from chemgridmap.evidence import document_values, fit_between_unit_sd
from e2_analysis import E1A_PREDICTIONS, FROZEN, ORIGINAL_TARGETS
from e2_models import (
    S,
    build_target_predictions,
    clean_training,
    comparator_scores,
    pooled_global_fit,
)
from e2b_analysis import (
    ORIGINAL_IC50_TARGETS,
    add_unit_years,
    apply_holm_corrections,
    hypothesis_tests_by_stratum,
    holm_adjust,
    load_original_ic50_records,
    one_sided_wilcoxon,
    time_split_predictions,
)

ROOT = Path(__file__).resolve().parents[3]
E2_OUTPUT = ROOT / "paper/strengthening/output/E2"


def two_target_docs():
    rows = []
    for target, base in (("a", 0.0), ("b", 6.0)):
        for entry in range(5):
            for unit in range(2):
                rows.append({
                    "target": target,
                    "entry": f"{target}{entry}",
                    "document_value": base + unit,
                })
    frame = pd.DataFrame(rows)
    return {
        target: group[["entry", "document_value"]].reset_index(drop=True)
        for target, group in frame.groupby("target")
    }


def test_g_global_excludes_current_target():
    docs = two_target_docs()
    first = pooled_global_fit(docs, {"a": set(), "b": set()}, "a", "all")
    changed = {key: value.copy() for key, value in docs.items()}
    changed["a"]["document_value"] = [100.0, -100.0] * 5
    second = pooled_global_fit(changed, {"a": set(), "b": set()}, "a", "all")

    donor = docs["b"].copy()
    donor["entry"] = "b::" + donor.entry
    expected = fit_between_unit_sd(donor)
    assert first["tau"] == pytest.approx(expected["tau"])
    assert second["tau"] == pytest.approx(expected["tau"])
    assert first["donor_targets"] == ["b"]


def test_comparators_exclude_held_document_and_single_record_edge():
    source = pd.DataFrame({
        "document_chembl_id": ["a", "a", "b", "b", "c", "d"],
        "pchembl_value": [1.0, 2.0, 3.0, 4.0, 5.0, 9.0],
    })
    original = comparator_scores(source, "a;b")
    changed = source.copy()
    changed.loc[changed.document_chembl_id.isin(["a", "b"]), "pchembl_value"] = 50.0
    assert comparator_scores(changed, "a;b") == original

    range_score, sd_score, n = comparator_scores(source, "a;b;c")
    assert range_score == 1.0
    assert np.isnan(sd_score)
    assert n == 1


def test_flag_clean_sigma_uses_only_entries_after_flagged_removal():
    doc = pd.DataFrame({
        "entry": ["flagged"] * 2 + [f"clean{i}" for i in range(5) for _ in range(2)],
        "document_value": [0.0, 20.0] + [
            value for i, value in enumerate((6.0, 6.2, 6.1, 6.3, 6.0, 6.1,
                                             6.2, 6.4, 6.1, 6.3))
        ],
    })
    clean = clean_training(doc, {"flagged"})
    assert "flagged" not in set(clean.entry)
    assert fit_between_unit_sd(clean)["tau"] == pytest.approx(
        fit_between_unit_sd(doc.iloc[2:])["tau"]
    )
    assert fit_between_unit_sd(doc)["tau"] != pytest.approx(
        fit_between_unit_sd(clean)["tau"]
    )


def test_e2_gaussian_implementation_regresses_e1a_exactly():
    if not E1A_PREDICTIONS.is_file():
        pytest.skip("Run paper/strengthening/source/e1a_evidence_audit.py first.")
    e1a = pd.read_csv(
        E1A_PREDICTIONS, usecols=[
            "target", "fold", "entry", "held_out_unit", "method",
            "population", "p0", "p1", "p2",
        ]
    )
    e1a = e1a.loc[(e1a.method == "G") & (e1a.population == "all")]
    errors = []
    for target in ORIGINAL_TARGETS:
        lower = target.lower()
        records = pd.read_csv(
            FROZEN / lower / f"{lower}_retained_records.csv", low_memory=False
        )
        actual, _ = build_target_predictions(
            records,
            target,
            "IC50",
            models=("B0", "G"),
        )
        actual = actual.loc[
            (actual.method == "G") & (actual.population == "all")
        ]
        expected = e1a.loc[e1a.target == target.lower()]
        keys = ["fold", "entry", "held_out_unit"]
        merged = actual.merge(
            expected,
            on=keys,
            suffixes=("_e2", "_e1a"),
            validate="one_to_one",
        )
        assert len(actual) == len(expected) == len(merged)
        for column in ("p0", "p1", "p2"):
            errors.append(
                float(
                    np.max(
                        np.abs(
                            merged[f"{column}_e2"] - merged[f"{column}_e1a"]
                        )
                    )
                )
            )
    assert max(errors) <= 1e-12


def temporal_records():
    rows = []
    for entry in range(5):
        base = 4.8 + entry / 10
        for document, year, value in (
            (f"past-a-{entry}", 2000, base),
            (f"past-b-{entry}", 2005, base + 0.3),
            (f"future-{entry}", 2020, 3.5 + entry / 10),
        ):
            rows.append({
                "molecule_identity_key": f"E{entry}",
                "document_chembl_id": document,
                "document_year": year,
                "pchembl_value": value,
            })
    return pd.DataFrame(rows)


def test_time_split_never_uses_future_for_sigma_or_display():
    records = temporal_records()
    before, _ = time_split_predictions("T", "IC50", records, 2015)
    changed_records = records.copy()
    changed_records.loc[
        changed_records.document_chembl_id.str.startswith("future"), "pchembl_value"
    ] = 9.5
    after, _ = time_split_predictions("T", "IC50", changed_records, 2015)

    keys = ["target", "endpoint", "cutoff_year", "entry", "future_unit"]
    compared = before.merge(
        after,
        on=keys,
        suffixes=("_before", "_after"),
        validate="one_to_one",
    )
    assert compared.sigma_past_before.equals(compared.sigma_past_after)
    assert compared.past_center_before.equals(compared.past_center_after)
    assert compared.displayed_before.equals(compared.displayed_after)
    assert (compared.displayed_before == 0).all()


def test_time_split_excludes_and_counts_missing_year_units():
    records = temporal_records()
    missing = records.iloc[[0]].copy()
    missing["molecule_identity_key"] = "MISSING"
    missing["document_chembl_id"] = "missing-year"
    missing["document_year"] = np.nan
    records = pd.concat([records, missing], ignore_index=True)
    doc = document_values(records, merge_identical=True)
    with_years, missing_count = add_unit_years(doc, records)
    assert missing_count == 1
    assert with_years.loc[with_years.entry == "MISSING", "year"].isna().all()

    pred, diagnostics = time_split_predictions("T", "IC50", records, 2015)
    assert diagnostics["missing_year_units"] == 1
    assert diagnostics["missing_year_entries"] == 1
    assert "MISSING" not in set(pred.entry)


def test_holm_adjustment_matches_hand_calculation():
    observed = holm_adjust([0.01, 0.04, 0.03])
    assert np.allclose(observed, [0.03, 0.06, 0.06])


def test_wilcoxon_alternative_direction_is_one_sided():
    _statistic, greater, _n = one_sided_wilcoxon(
        [1, 2, 3, 4, 5], alternative="greater"
    )
    _statistic, less, _n = one_sided_wilcoxon(
        [1, 2, 3, 4, 5], alternative="less"
    )
    assert greater == pytest.approx(0.03125)
    assert less > 0.95
    assert greater < less


def test_h8_has_a_separate_holm_family():
    rows = []
    for endpoint in ("IC50", "Ki"):
        for hypothesis, p_value in zip(
            ("H3", "H4", "H5", "H6", "H7"),
            (0.01, 0.02, 0.03, 0.04, 0.05),
        ):
            rows.append({
                "hypothesis": hypothesis,
                "endpoint": endpoint,
                "population": "all",
                "correction_family": "H3-H7",
                "raw_p": p_value,
            })
    for endpoint, p_value in (("IC50", 0.04), ("Ki", 0.03)):
        rows.append({
            "hypothesis": "H8",
            "endpoint": endpoint,
            "population": "all",
            "correction_family": "H8",
            "raw_p": p_value,
        })
    adjusted = apply_holm_corrections(pd.DataFrame(rows))
    h8 = adjusted.loc[adjusted.hypothesis == "H8"]
    assert np.allclose(h8.holm_p, holm_adjust([0.04, 0.03]))
    first_family = adjusted.loc[
        (adjusted.correction_family == "H3-H7")
        & (adjusted.endpoint == "IC50")
    ]
    assert np.allclose(
        first_family.holm_p,
        holm_adjust([0.01, 0.02, 0.03, 0.04, 0.05]),
    )


def test_time_split_includes_original_ic50_dataset():
    if not (S / "data/chembl37_api" / "activities").is_dir():
        pytest.skip("Download the selected ChEMBL 37 API records first.")
    records = load_original_ic50_records()
    assert set(records) == {
        (target, "IC50") for target in ORIGINAL_IC50_TARGETS
    }

    by_target = pd.read_csv(E2_OUTPUT / "time_split_by_target.csv")
    original = by_target.loc[
        (by_target.endpoint == "IC50")
        & by_target.target.isin(ORIGINAL_IC50_TARGETS)
    ]
    assert set(original.target) == set(ORIGINAL_IC50_TARGETS)
    assert set(original.cutoff_year) == {2010, 2015}
    assert original.eligible.all()
    assert (original.n_predictions > 0).all()


def test_original_ic50_targets_are_excluded_from_h1_to_h8():
    if not (E2_OUTPUT / "target_inclusion.csv").is_file():
        pytest.skip("Run paper/strengthening/source/e2_analysis.py first.")
    inclusion = pd.read_csv(E2_OUTPUT / "target_inclusion.csv")
    included_original_ic50 = inclusion.loc[
        (inclusion.endpoint == "IC50")
        & inclusion.target.isin(ORIGINAL_IC50_TARGETS)
        & inclusion.included
    ]
    assert included_original_ic50.empty

    tests = pd.read_csv(E2_OUTPUT / "hypothesis_tests.csv")
    ic50 = tests.loc[tests.endpoint == "IC50"]
    eligible = {
        target
        for value in ic50.eligible_targets.dropna().astype(str)
        for target in value.split(";")
        if target
    }
    assert eligible.isdisjoint(ORIGINAL_IC50_TARGETS)


def test_stratified_h3_matches_metrics_directly():
    if not (E2_OUTPUT / "metrics_by_target.csv").is_file():
        pytest.skip("Run paper/strengthening/source/e2_analysis.py first.")
    metrics = pd.read_csv(
        E2_OUTPUT / "metrics_by_target.csv",
        low_memory=False,
        float_precision="round_trip",
    )
    expected = hypothesis_tests_by_stratum(metrics)
    saved = pd.read_csv(
        E2_OUTPUT / "hypothesis_tests_by_stratum.csv",
        float_precision="round_trip",
    )
    assert saved.equals(expected)
    assert set(saved.family_note) == {
        "stratified; not part of the Holm families"
    }

    for row in saved.itertuples(index=False):
        layer = metrics.loc[
            (metrics.endpoint == row.endpoint)
            & (metrics.population == row.population)
            & (metrics.model == "G")
            & (metrics.stratum == row.stratum)
            & np.isfinite(metrics.auc)
        ].sort_values("target")
        values = layer.auc.to_numpy(float) - 0.5
        statistic, p_value, n_targets = one_sided_wilcoxon(
            values, alternative="greater"
        )
        assert n_targets == row.n_targets == len(layer)
        assert statistic == pytest.approx(row.h3_wilcoxon_statistic)
        assert p_value == pytest.approx(row.h3_raw_p)
        assert float(np.median(values)) == pytest.approx(
            row.h3_median_difference
        )
