from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from e15_leakfree import (
    E5,
    E5C,
    archived_document_g,
    archived_rule_frame,
    document_fold_predictions_e15,
    filtered_observations,
    load_e2_records,
    recompute_archived_flag_touched,
    strict_fold_masks,
    wilcoxon_details,
)
from e15_offset_audit import SwapAuditor, s2_status, s3_status


def test_pairs_containing_held_document_are_ignored():
    observations = pd.DataFrame(
        {
            "document_a": ["A", "A", "B"],
            "document_b": ["B", "C", "C"],
            "entry": ["e1", "e2", "e3"],
            "difference": [1.0, 2.0, 3.0],
        }
    )
    selected = filtered_observations(
        observations,
        held_documents={"A"},
    )
    assert selected[["document_a", "document_b"]].values.tolist() == [["B", "C"]]


def test_strict_document_fold_removes_merged_units_touching_held_fold():
    units = pd.DataFrame(
        {
            "documents": ["A", "B", "C;A"],
            "document_fold": [0, 1, 1],
        }
    )
    training, held, held_documents = strict_fold_masks(units, 0)
    assert held_documents == {"A"}
    assert held.tolist() == [True, False, False]
    assert training.tolist() == [False, True, False]


def _swap_records(target: str, swapped: bool) -> pd.DataFrame:
    parents = ["P1", "P2", "P3"]
    first = [5.0, 6.0, 7.0]
    second = [7.0, 8.0, 9.0]
    if swapped:
        first, second = second, first
    rows = []
    for document, values in (("D1", first), ("D2", second)):
        for parent, value in zip(parents, values):
            rows.append(
                {
                    "document_chembl_id": document,
                    "parent_molecule_chembl_id": parent,
                    "pchembl_value": value,
                    "standard_relation": "=",
                }
            )
    result = pd.DataFrame(rows)
    result["target_chembl_id"] = target
    return result


def test_s1_positive_and_negative():
    auditor = SwapAuditor(
        {
            ("T1", "IC50"): _swap_records("T1", False),
            ("T2", "IC50"): _swap_records("T2", True),
        }
    )
    assert auditor.evaluate("D1", "D2", "T1", "IC50")["is_swap"]

    auditor = SwapAuditor(
        {
            ("T1", "IC50"): _swap_records("T1", False),
            ("T2", "IC50"): _swap_records("T2", False),
        }
    )
    assert not auditor.evaluate("D1", "D2", "T1", "IC50")["is_swap"]


def test_s2_positive_and_negative():
    assert s2_status(np.array([2.9, 3.0, 3.1, 3.0]), 3.0)
    assert not s2_status(np.array([2.2, 2.5, 3.1, 3.8]), 2.8)


def test_s3_positive_and_negative():
    assert s3_status(np.array([1.1, 1.2, 1.3]), 1.2)
    assert not s3_status(np.array([0.5, 1.0, 1.5]), 1.0)


def test_wilcoxon_method_selection():
    exact = wilcoxon_details(np.arange(1, 11))
    approximate = wilcoxon_details(np.arange(1, 31))
    tied = wilcoxon_details(np.array([1, 1, 2, 3]))
    assert exact["wilcoxon_method"] == "exact"
    assert approximate["wilcoxon_method"] == "approx"
    assert tied["wilcoxon_method"] == "approx"


def _first_regression_dataset():
    if not (E5C / "predictions_e5c.csv").exists():
        pytest.skip("archived E5c outputs are not present")
    records_by_key = load_e2_records()
    archived = archived_rule_frame()
    first_key = sorted(records_by_key)[0]
    entry = archived.loc[
        (archived.split_design == "entry")
        & (archived.merge_basis == "same_value")
        & (archived.target == first_key[0])
        & (archived.endpoint == first_key[1])
    ]
    return first_key, records_by_key, archived, entry


def test_r1_archived_flags_with_data():
    key, records_by_key, _, entry = _first_regression_dataset()
    result = recompute_archived_flag_touched(
        {key: records_by_key[key]}, entry
    )
    assert result.matches.all()


def test_r4_document_replay_with_data():
    if not (E5 / "document_folds_predictions.csv").exists():
        pytest.skip("archived E5 document predictions are not present")
    key, records_by_key, archived, _ = _first_regression_dataset()
    replay = document_fold_predictions_e15(
        records_by_key[key],
        key[0],
        key[1],
        strict_leak_free=False,
    )
    expected = archived_document_g()
    expected = expected.loc[
        (expected.target == key[0]) & (expected.endpoint == key[1])
    ].reset_index(drop=True)
    columns = [
        "fold",
        "entry",
        "held_out_unit",
        "displayed",
        "observed",
        "p0",
        "p1",
        "p2",
        "p_displayed",
        "support",
        "agree",
    ]
    pd.testing.assert_frame_equal(
        expected[columns].sort_values(columns[:3]).reset_index(drop=True),
        replay[columns].sort_values(columns[:3]).reset_index(drop=True),
        check_exact=False,
        rtol=0,
        atol=1e-12,
    )
