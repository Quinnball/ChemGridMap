"""Tests for the 0.6.0 additions made in response to the v0.6 review."""
import numpy as np
import pandas as pd
import pytest

from chemgridmap.inspection import excluded_bound_records

pytest.importorskip("rdkit")
from chemgridmap.queries import placement_recall  # noqa: E402


def _excluded():
    return pd.DataFrame({
        "Molecule ChEMBL ID": ["CHEMBL1", "CHEMBL1", "CHEMBL2", "CHEMBL1"],
        "Standard Relation": ["'>'", "<=", ">", ">"],
        "Standard Value": [10000, 5, 30000, 1000],
        "Standard Units": ["nM", "nM", "nM", "ug.mL-1"],
        "Document ChEMBL ID": ["CHEMBL_D1", "CHEMBL_D2", "CHEMBL_D3", "CHEMBL_D4"],
        "exclusion_reason": ["non_exact_standard_relation", "non_exact_standard_relation",
                             "non_exact_standard_relation", "missing_units;non_exact_standard_relation"],
    })


def test_bound_records_are_listed_for_the_selected_entry_only():
    result = excluded_bound_records(_excluded(), "chembl1")
    assert result["count"] == 3
    first, second, third = result["records"]
    assert (first["relation"], first["p_relation"], first["document"]) == (">", "<", "CHEMBL_D1")
    assert first["p_bound"] == pytest.approx(5.0)
    assert (second["relation"], second["p_relation"]) == ("<=", "≥")
    assert second["p_bound"] == pytest.approx(-np.log10(5e-9))
    # Bounds in other units are listed but not converted to pActivity.
    assert third["units"] == "ug.mL-1" and third["p_bound"] is None


def test_bound_records_handle_missing_tables_and_limit():
    assert excluded_bound_records(None, "CHEMBL1") == {"count": 0, "records": []}
    assert excluded_bound_records(pd.DataFrame(), "CHEMBL1") == {"count": 0, "records": []}
    other = _excluded().assign(exclusion_reason="duplicate_activity_id")
    assert excluded_bound_records(other, "CHEMBL1")["count"] == 0
    limited = excluded_bound_records(_excluded(), "CHEMBL1", limit=1)
    assert limited["count"] == 3 and len(limited["records"]) == 1


def _reference():
    # Two series in opposite corners of a 10 x 10 grid.
    benzenes = ["c1ccccc1C", "c1ccccc1CC", "c1ccccc1CCC", "c1ccccc1CCCC", "c1ccccc1CCCCC"]
    pyridines = ["c1ccncc1N", "c1ccncc1NC", "c1ccncc1NCC", "c1ccncc1NCCC", "c1ccncc1NCCCC"]
    rows = [(s, 0, i) for i, s in enumerate(benzenes)] + [(s, 9, 9 - i) for i, s in enumerate(pyridines)]
    return pd.DataFrame(rows, columns=["canonical_smiles", "grid_row", "grid_col"])


def test_placement_recall_is_high_next_to_neighbours_and_low_far_away():
    reference = _reference()
    queries = pd.DataFrame({"canonical_smiles": ["c1ccccc1CCCCCC", "c1ccccc1CCCCCC"],
                            "grid_row": [1, 9], "grid_col": [2, 5]})
    result = placement_recall(reference, queries, k=5)
    near, far = result.placement_recall.tolist()
    assert near == pytest.approx(1.0)
    assert far < 0.5
    assert result.placement_recall_fingerprint_count.tolist() == [5, 5]
    assert result.placement_recall_shared.iloc[0] == 5


def test_placement_recall_keeps_ties_at_the_kth_similarity():
    reference = pd.DataFrame({"canonical_smiles": ["CCO", "CCO", "CCO", "c1ccccc1"],
                              "grid_row": [0, 0, 5, 9], "grid_col": [0, 1, 5, 9]})
    queries = pd.DataFrame({"canonical_smiles": ["CCO"], "grid_row": [0], "grid_col": [2]})
    result = placement_recall(reference, queries, k=2)
    # Three reference molecules share the query's fingerprint, so all three are counted.
    assert result.placement_recall_fingerprint_count.iloc[0] == 3
    assert result.placement_recall_shared.iloc[0] == 2
