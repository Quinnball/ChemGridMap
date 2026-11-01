"""E5 model formulas, regression, and fold-leakage tests."""
from __future__ import annotations

import gzip
from pathlib import Path
import pickle

import numpy as np
import pandas as pd
import pytest

from chemgridmap.evidence import document_values
from e2_models import FOLD_SEED
from e4c_analysis import near_duplicate_units
from e5_models import (
    b1c_probabilities,
    b1c_training_statistics,
    document_fold_map,
    document_fold_predictions,
    eb_prediction,
    empirical_bayes_prior,
    entry_fold_predictions,
    near_duplicate_predictions,
    pooled_sigma,
)
from reliability_activity import entry_folds

ROOT = Path(__file__).resolve().parents[3]
S = ROOT / "paper/strengthening"


def test_e5_g_matches_e2_row_by_row():
    if not (S / "qa/E2/curated_records.pkl.gz").is_file():
        pytest.skip("Run E2 before the archived row regression.")
    records = pickle.loads(
        gzip.decompress((S / "qa/E2/curated_records.pkl.gz").read_bytes())
    )[("CHEMBL1994", "IC50")]
    all_predictions = entry_fold_predictions(
        records, "CHEMBL1994", "IC50"
    )
    actual = all_predictions.loc[all_predictions.method == "G"]
    b1c = all_predictions.loc[all_predictions.method == "B1c"]
    assert len(b1c) == len(actual)
    assert not b1c.duplicated(
        ["fold", "entry", "held_out_unit"]
    ).any()
    expected = pd.read_csv(
        S / "output/E2/predictions.csv",
        low_memory=False,
        float_precision="round_trip",
    )
    expected = expected.loc[
        (expected.target == "CHEMBL1994")
        & (expected.endpoint == "IC50")
        & (expected.population == "all")
        & (expected.method == "G")
    ]
    keys = ["fold", "entry", "held_out_unit"]
    merged = actual.merge(
        expected, on=keys, suffixes=("_e5", "_e2"), validate="one_to_one"
    )
    assert len(actual) == len(expected) == len(merged)
    for column in ("p0", "p1", "p2"):
        assert np.max(
            np.abs(
                merged[f"{column}_e5"] - merged[f"{column}_e2"]
            )
        ) <= 1e-12


def test_eb_posterior_and_tau_floor_by_hand():
    training = pd.DataFrame(
        {
            "entry": ["a", "a", "b", "b"],
            "document_value": [6.0, 8.0, 7.0, 9.0],
        }
    )
    prior = empirical_bayes_prior(training, sigma=1.0)
    assert prior["tau2"] == pytest.approx(0.25**2)
    support, predictive, posterior_mean, posterior_var = eb_prediction(
        [6.0], 1.0, prior
    )
    precision = 1.0 / 0.25**2 + 1.0
    assert posterior_var == pytest.approx(1.0 / precision)
    assert posterior_mean == pytest.approx(
        (prior["mu0"] / prior["tau2"] + 6.0) / precision
    )
    assert support.sum() == pytest.approx(1.0)
    assert predictive.sum() == pytest.approx(1.0)


def test_b1c_uses_only_other_folds():
    frame = pd.DataFrame(
        [
            {
                "fold": fold,
                "method": "G",
                "displayed": 0,
                "observed": fold % 3,
                "m": 1,
            }
            for fold in range(10)
        ]
    )
    statistics = b1c_training_statistics(frame, held_fold=4)
    assert set(statistics.index) == set(range(10)) - {4}
    probabilities = b1c_probabilities(statistics, displayed=0, m=1)
    expected_counts = np.ones(3)
    for observed in [0, 1, 2, 0, 2, 0, 1, 2, 0]:
        expected_counts[observed] += 1
    assert np.allclose(probabilities, expected_counts / expected_counts.sum())
    assert 4 not in set(statistics.index)


def synthetic_records() -> pd.DataFrame:
    rows = []
    for entry_index in range(12):
        entry = f"E{entry_index:02d}"
        for unit_index, offset in enumerate((0.0, 0.4, 0.8)):
            rows.append(
                {
                    "molecule_identity_key": entry,
                    "document_chembl_id": (
                        f"D{entry_index:02d}{unit_index}"
                    ),
                    "pchembl_value": 5.0 + entry_index / 10 + offset,
                    "assay_type": "B" if unit_index < 2 else "F",
                }
            )
    return pd.DataFrame(rows)


def test_document_fold_sigma_excludes_held_documents():
    records = synthetic_records()
    predictions = document_fold_predictions(records, "T", "IC50")
    doc = document_values(records, merge_identical=True)
    doc["document_fold"] = document_fold_map(doc)
    for fold, group in predictions.groupby("fold"):
        expected = pooled_sigma(
            doc.loc[doc.document_fold != fold]
        )["tau"]
        assert np.allclose(group.sigma_doc, expected)


def test_near_duplicate_sigma_is_fold_out():
    records = synthetic_records()
    predictions = near_duplicate_predictions(records, "T", "IC50")
    original = document_values(records, merge_identical=True)
    assignments = entry_folds(original, seed=FOLD_SEED)
    units = near_duplicate_units(records)
    for fold, group in predictions.groupby("fold"):
        held_entries = {
            entry
            for entry, assigned in assignments.items()
            if assigned == fold
        }
        expected = pooled_sigma(
            units.loc[~units.entry.isin(held_entries)]
        )["tau"]
        assert np.allclose(group.sigma_doc, expected)
