import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from chemgridmap.evidence import (
    DOCUMENT_OFFSET_COLUMNS,
    annotate_activity_confidence,
    document_pair_offsets,
)


def _records(rows):
    return pd.DataFrame(rows, columns=["molecule_identity_key", "document_chembl_id", "pchembl_value"])


def _molecules(records):
    return records.groupby("molecule_identity_key").pchembl_value.median().rename("activity_pchembl").reset_index()


def test_document_offsets_use_entry_document_medians_and_count_pairs():
    records = _records([
        ("E1", "A", 4), ("E1", "A", 6), ("E1", "B", 7), ("E1", "C", 5),
        ("E2", "A", 6), ("E2", "B", 8), ("E2", "C", 6),
        ("E3", "A", 7), ("E3", "B", 9), ("E4", "A", 6), ("E5", "D", 8),
    ])
    original = records.copy(deep=True)
    offsets = document_pair_offsets(records.rename(columns={
        "molecule_identity_key": "entry_id", "document_chembl_id": "source", "pchembl_value": "measurement",
    }), entry="entry_id", value="measurement", document="source")
    assert len(offsets) == 1
    row = offsets.iloc[0]
    assert (row.document_a, row.document_b, row.shared_entries) == ("A", "B", 3)
    assert row.median_offset == row.mean_offset == 2
    assert row.std_offset == 0 and row.sign_agreement == 1 and row.systematic_offset
    assert offsets.attrs["document_pair_possible_count"] == 6
    assert offsets.attrs["observed_document_pair_count"] == 3
    assert offsets.attrs["document_pair_count"] == 1

    annotated, fit = annotate_activity_confidence(_molecules(records), records)
    keyed = annotated.set_index("molecule_identity_key")
    assert keyed.loc[["E1", "E2", "E3"], "discordant_document_pairs"].tolist() == [1, 1, 1]
    assert keyed.loc[["E1", "E2", "E3"], "discordant_documents"].tolist() == ["A;B"] * 3
    assert keyed.loc[["E4", "E5"], "discordant_document_pairs"].tolist() == [0, 0]
    assert keyed.loc[["E4", "E5"], "discordant_documents"].tolist() == ["", ""]
    assert fit["document_pair_count"] == fit["discordant_document_pair_count"] == 1
    assert fit["document_offsets_top10"] == fit["document_offsets"] == offsets.to_dict("records")
    json.dumps(fit, allow_nan=False)
    pd.testing.assert_frame_equal(records, original)


@pytest.mark.parametrize("differences, agreement", [([-2, -2, -2, 2], 0.75), ([2, 2, 0, 0, 2], 0.6)])
def test_document_offsets_require_consistent_sign_and_respect_thresholds(differences, agreement):
    records = _records([(str(i), source, value)
                        for i, difference in enumerate(differences)
                        for source, value in [("A", 6), ("B", 6 + difference)]])
    result = document_pair_offsets(records).iloc[0]
    assert result.median_offset == np.median(differences)
    assert result.mean_offset == np.mean(differences)
    assert result.std_offset == pytest.approx(np.std(differences, ddof=1))
    assert result.sign_agreement == agreement
    assert not result.systematic_offset
    assert document_pair_offsets(records, sign_agreement=agreement).iloc[0].systematic_offset
    assert not document_pair_offsets(records, offset_threshold=2.01, sign_agreement=agreement).iloc[0].systematic_offset
    pd.testing.assert_frame_equal(document_pair_offsets(records), document_pair_offsets(records.iloc[::-1]))


def test_activity_evidence_unit_policy_can_keep_identical_documents():
    records = _records([("E1", "A", 6), ("E1", "B", 6), ("E1", "C", 8)])
    molecules = _molecules(records)
    merged, merged_fit = annotate_activity_confidence(molecules, records)
    separate, separate_fit = annotate_activity_confidence(molecules, records, merge_identical=False)
    assert merged.n_evidence_units.tolist() == [2]
    assert separate.n_evidence_units.tolist() == [3]
    assert merged.document_balanced_pchembl.iloc[0] == 7
    assert separate.document_balanced_pchembl.iloc[0] == pytest.approx(20 / 3)
    assert merged_fit["evidence_unit_policy"] == "document; identical values merged"
    assert separate_fit["evidence_unit_policy"] == "document; identical values kept separate"
    assert merged_fit["merge_identical_documents"] is True
    assert separate_fit["merge_identical_documents"] is False
    assert merged_fit["document_pair_count"] == separate_fit["document_pair_count"] == 0
    assert separate_fit["observed_document_pair_count"] == 3
    pd.testing.assert_frame_equal(merged[molecules.columns], molecules)
    pd.testing.assert_frame_equal(separate[molecules.columns], molecules)


def test_document_offsets_empty_result_has_export_columns():
    result = document_pair_offsets(_records([("E1", "A", 7), ("E2", "B", 8)]))
    assert result.empty and list(result.columns) == list(DOCUMENT_OFFSET_COLUMNS)
    assert result.attrs["document_pair_possible_count"] == 1
    assert result.attrs["observed_document_pair_count"] == result.attrs["document_pair_count"] == 0


def test_identical_unmerged_documents_give_finite_point_mass_probabilities():
    values = [6, 7, 6.5, 8, 5]
    records = _records([(str(i), document, value) for i, value in enumerate(values) for document in ("A", "B")])
    annotated, fit = annotate_activity_confidence(_molecules(records), records, merge_identical=False)
    assert fit["method"] == "pooled_between_unit_sd" and fit["tau"] == 0
    expected = np.eye(3)[[0, 2, 1, 2, 0]]
    for prefix in ("p_consensus_", "p_new_document_"):
        assert np.array_equal(annotated[[prefix + c for c in ("inactive", "medium", "active")]], expected)
    assert annotated.class_confidence.tolist() == [1] * 5
    assert annotated.class_status.tolist() == ["supported"] * 5
    json.dumps(fit, allow_nan=False)
    json.dumps(annotated.to_dict("records"), allow_nan=False)


def test_document_offsets_top_ten_are_ranked_by_shared_entries():
    records = _records([(str(i), f"D{document:02}", 4 + document)
                        for document in range(6) for i in range(3 + document)])
    _, fit = annotate_activity_confidence(_molecules(records), records)
    assert fit["discordant_document_pair_count"] == 15
    assert len(fit["document_offsets_top10"]) == 10
    counts = [pair["shared_entries"] for pair in fit["document_offsets"]]
    assert counts == sorted(counts, reverse=True)
    assert fit["document_offsets_top10"] == fit["document_offsets"][:10]


def test_herg_document_offset_matches_frozen_source_records():
    path = Path(__file__).resolve().parents[1] / "paper/output/revision_v6/chembl240/chembl240_retained_records.csv"
    if not path.exists():
        pytest.skip("The full frozen hERG source records are available in the research checkout.")
    records = pd.read_csv(path, usecols=["molecule_identity_key", "document_chembl_id", "pchembl_value"])
    offsets = document_pair_offsets(records)
    pair = offsets[(offsets.document_a == "CHEMBL5725903") & (offsets.document_b == "CHEMBL5727624")].iloc[0]
    assert pair.shared_entries == 30
    assert pair.median_offset == pytest.approx(3.5, abs=0.1)
    assert pair.sign_agreement == 1 and pair.systematic_offset
    relevant = records[records.document_chembl_id.isin([pair.document_a, pair.document_b])]
    annotated, fit = annotate_activity_confidence(_molecules(relevant), relevant)
    assert (annotated.discordant_document_pairs == 1).sum() == 30
    assert fit["discordant_document_pair_count"] == 1
