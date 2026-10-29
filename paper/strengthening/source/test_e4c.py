"""E4c 的范围限定单元测试。"""
from __future__ import annotations

import numpy as np
import pandas as pd

from e4c_analysis import near_duplicate_units, shared_entry_stratum
from e2_models import (
    document_pair_observations,
    flags_from_observations,
    touch_flagged_entries,
)
from chemgridmap.evidence import document_values


def _records(rows):
    return pd.DataFrame(
        rows,
        columns=[
            "molecule_identity_key",
            "pchembl_value",
            "document_chembl_id",
        ],
    )


def test_near_duplicate_units_use_single_link_clusters_and_means():
    records = _records(
        [
            ("E1", 7.00, "D1"),
            ("E1", 7.05, "D2"),
            ("E1", 7.12, "D3"),
            ("E1", 7.31, "D4"),
            ("E2", 6.00, "D1"),
        ]
    )
    units = near_duplicate_units(records)
    entry = units.loc[units.entry == "E1"].sort_values("document_value")
    assert len(entry) == 2
    assert np.isclose(entry.document_value.iloc[0], np.mean([7.00, 7.05, 7.12]))
    assert set(entry.documents.iloc[0].split(";")) == {"D1", "D2", "D3"}
    assert np.isclose(entry.document_value.iloc[1], 7.31)


def test_shared_entry_strata_accept_float_counts_from_csv():
    assert shared_entry_stratum(3.0) == "3"
    assert shared_entry_stratum(4.0) == "4"
    assert shared_entry_stratum(5.0) == ">=5"


def test_training_fold_flags_do_not_use_test_entries():
    rows = []
    for entry in ("TRAIN1", "TRAIN2", "TRAIN3", "TEST1", "TEST2", "TEST3"):
        group = "TEST" if entry.startswith("TEST") else "TRAIN"
        for document, value in (("A_" + group, 5.0), ("B_" + group, 6.5)):
            rows.append((entry, value, document))
    records = _records(rows)
    training_entries = {"TRAIN1", "TRAIN2", "TRAIN3"}
    observations = document_pair_observations(records)
    assert len(flags_from_observations(observations)) == 4
    flags = flags_from_observations(observations, training_entries)
    assert flags == {"A_TRAIN", "B_TRAIN"}
    touched = touch_flagged_entries(
        document_values(records.loc[records.molecule_identity_key.isin(training_entries)]),
        flags,
    )
    assert touched == training_entries


def test_flagged_training_documents_exclude_touched_test_entries():
    rows = []
    for entry in ("TRAIN1", "TRAIN2", "TRAIN3"):
        for document, value in (("A_TRAIN", 5.0), ("B_TRAIN", 6.5)):
            rows.append((entry, value, document))
    rows.extend(
        [
            ("TEST_TOUCHED", 5.0, "A_TRAIN"),
            ("TEST_TOUCHED", 6.0, "C_TEST"),
            ("TEST_CLEAN", 5.0, "C_TEST"),
            ("TEST_CLEAN", 6.0, "D_TEST"),
        ]
    )
    records = _records(rows)
    training_entries = {"TRAIN1", "TRAIN2", "TRAIN3"}
    flags = flags_from_observations(
        document_pair_observations(records), training_entries
    )
    touched = touch_flagged_entries(document_values(records), flags)

    assert "TEST_TOUCHED" in touched
    assert "TEST_CLEAN" not in touched
