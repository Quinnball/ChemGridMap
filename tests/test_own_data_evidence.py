"""Own-data units and explicitly assumed evidence uncertainty."""
import json

import numpy as np
import pandas as pd
import pytest
from scipy.stats import norm

from chemgridmap.core import prepare_molecules
from chemgridmap.evidence import annotate_activity_confidence, prepare_generic_activity_records
from chemgridmap.inspection import inspect_entry


def _records():
    return pd.DataFrame({
        "compound": ["A", "A", "A", "B", "B", "C", "C"],
        "smiles": ["CCO"] * 3 + ["CCN"] * 2 + ["CCC"] * 2,
        "value": [5.2, 6.8, 7.6, 5.5, 7.2, 6.0, 7.8],
        "batch": ["D1", "D1", "D2", "D1", "D2", "D1", "D2"],
    })


def _prepare(frame, **kwargs):
    return prepare_generic_activity_records(
        frame, compound_column="compound", value_column="value",
        smiles_column="smiles", **kwargs,
    )


@pytest.mark.parametrize("unit,factor", [("nM", 1e-9), ("uM", 1e-6), ("M", 1), ("p", None)])
def test_record_unit_conversion_matches_pactivity(unit, factor):
    source = _records()
    expected, _, _, _ = _prepare(source, document_column="batch")
    original = source.value.copy()
    if factor:
        source["value"] = 10.0 ** -source.value / factor
    molecules, retained, _, report = _prepare(source, document_column="batch", value_unit=unit)
    np.testing.assert_allclose(retained.pchembl_value, original, atol=1e-12, rtol=0)
    np.testing.assert_allclose(molecules.activity_pchembl, expected.activity_pchembl, atol=1e-12, rtol=0)
    if factor:
        assert retained.raw_activity_value.tolist() == source.value.tolist()
        assert retained.raw_activity_unit.eq(unit).all()
        assert not any(report["excluded_reasons"].values())
        grid = molecules.assign(grid_row=np.arange(len(molecules)), grid_col=0)
        for identity in grid.molecule_identity_key:
            _, _, summary = inspect_entry(grid, retained, molecule_id=identity)
            assert summary["median_verified"]
    else:
        assert "raw_activity_value" not in retained
        assert "excluded_reasons" not in report


def test_per_record_unicode_units_and_exclusion_reasons():
    source = pd.DataFrame({
        "compound": ["A"] * 4 + ["B"] * 8,
        "smiles": ["CCO"] * 4 + ["CCN"] * 5 + ["invalid", "CCN", "CCN"],
        "value": [1, 1, 1000, 1e-6, 0, -2, "bad", float("inf"), 1, 1, 1, 1],
        "unit": ["µM", " μM ", "nM", "M", "uM", "uM", "uM", "uM", "NM", "uM", "uM", "uM"],
        "relation": ["="] * 11 + ["<"],
    })
    source.loc[10, "compound"] = ""
    molecules, retained, fit, report = _prepare(
        source, unit_column="unit", relation_column="relation", confidence=False,
    )
    assert len(molecules) == 1 and len(retained) == 4 and fit is None
    np.testing.assert_allclose(retained.pchembl_value, 6, atol=1e-12, rtol=0)
    assert retained.raw_activity_unit.tolist() == ["µM", " μM ", "nM", "M"]
    assert report["excluded_reasons"] == {
        "invalid_activity_value": 2, "nonpositive_activity_value": 2,
        "unrecognized_activity_unit": 1, "invalid_smiles": 1,
        "missing_compound_id": 1, "non_exact_relation": 1,
    }
    assert report["excluded_activity_records"] == sum(report["excluded_reasons"].values())
    assert report["excluded_non_exact_relation_records"] == 1


@pytest.mark.parametrize("unit,factor", [("nM", 1e-9), ("uM", 1e-6), ("M", 1), ("p", None)])
def test_molecule_table_unit_conversion_matches_pactivity(unit, factor):
    frame = pd.DataFrame({"smiles": ["CCO", "CCN", "CCC"], "value": [5.2, 6.7, 8.1]})
    expected = frame.value.copy()
    if factor:
        frame["value"] = 10.0 ** -frame.value / factor
    molecules, _ = prepare_molecules(frame, value_col="value", value_unit=unit)
    np.testing.assert_allclose(molecules.value, expected, atol=1e-12, rtol=0)
    assert molecules.activity_class.tolist() == ["inactive", "medium", "active"]


def test_molecule_units_preserve_source_rows_and_count_exclusions():
    frame = pd.DataFrame({"smiles": ["CCO", "CCN", "CCC", "CO", "CF"],
                          "value": [0, "bad", 1, 1, 1], "unit": ["nM", "nM", "µM", "bad", "μM"]})
    molecules, report = prepare_molecules(frame, value_col="value", unit_column="unit")
    assert molecules.source_row.tolist() == [2, 4]
    assert molecules.value.tolist() == [6, 6]
    assert report["input_rows"] == 5 and report["output_rows"] == 2
    assert report["excluded_reasons"] == {"invalid_activity_value": 1, "nonpositive_activity_value": 1,
                                          "unrecognized_activity_unit": 1, "invalid_smiles": 0}


def test_assumed_sigma_with_documents_matches_same_prefitted_sigma():
    molecules, records, _, _ = _prepare(_records(), document_column="batch", confidence=False)
    expected, _ = annotate_activity_confidence(
        molecules, records, document="batch", fit={"tau": 0.64, "prior_df": np.inf},
    )
    actual, fit = annotate_activity_confidence(molecules, records, document="batch", assumed_sigma=0.64)
    pd.testing.assert_frame_equal(actual[expected.columns], expected)
    assert fit["method"] == "assumed_sigma" and fit["assumed_sigma"] == 0.64
    assert fit["evidence_unit_source"] == "batch"
    assert fit["document_pair_count"] == 1 and fit["discordant_document_pair_count"] == 1
    assert fit["sigma_doc_degrees_of_freedom"] is None
    assert actual.sigma_doc.eq(0.64).all() and actual.sigma_doc_degrees_of_freedom.isna().all()
    assert fit["warnings"] == ["Class support uses an assumed σ of 0.64 log units; it was not estimated from these data."]
    json.dumps(fit, allow_nan=False)


def test_assumed_sigma_without_documents_uses_one_median_unit():
    molecules, records, fit, report = _prepare(_records().drop(columns="batch"), assumed_sigma=0.64)
    centers = molecules.activity_pchembl.to_numpy()
    for scale, prefix in ((0.64, "p_consensus_"), (0.64 * np.sqrt(2), "p_new_document_")):
        low, high = norm.cdf((6 - centers) / scale), norm.sf((7 - centers) / scale)
        expected = np.column_stack([low, 1 - low - high, high])
        np.testing.assert_allclose(molecules[[prefix + name for name in ("inactive", "medium", "active")]],
                                   expected, atol=1e-12, rtol=0)
        shown = np.array([{"inactive": 0, "medium": 1, "active": 2}[label] for label in molecules.activity_class])
        output = molecules.class_confidence if prefix == "p_consensus_" else molecules.class_reproducibility
        np.testing.assert_allclose(output, expected[np.arange(len(molecules)), shown], atol=1e-12, rtol=0)
    assert molecules.n_evidence_units.eq(1).all() and molecules.n_documents.eq(0).all()
    assert molecules.top_unit_documents.eq("").all() and not molecules.offset_pair_touched.any()
    assert fit["evidence_unit_source"] == "one unit per compound"
    assert fit["document_offset_detection"] == "not applicable"
    assert fit["document_pair_offsets"] == [] and report["activity_confidence"] == fit
    assert "document_chembl_id" not in records
    reused, reused_fit = annotate_activity_confidence(
        molecules, records, document=None, fit=fit,
    )
    pd.testing.assert_frame_equal(reused, molecules)
    assert reused_fit == fit


def test_assumed_sigma_is_not_refitted_when_excluding_offset_documents():
    molecules, _, fit, _ = _prepare(_records(), document_column="batch", assumed_sigma=0.2,
                                   evidence_model="exclude-flagged")
    assert fit["method"] == "assumed_sigma" and fit["tau"] == 0.2
    assert fit["evidence_unit_source"] == "batch"
    assert molecules.sigma_doc.eq(0.2).all()


@pytest.mark.parametrize("sigma", [0, -1, float("nan"), float("inf")])
def test_assumed_sigma_must_be_finite_and_positive(sigma):
    with pytest.raises(ValueError, match="finite number greater than 0"):
        _prepare(_records(), assumed_sigma=sigma)


def test_molecule_table_rejects_assumed_sigma_and_conflicting_units():
    with pytest.raises(ValueError, match="activity-record table"):
        prepare_molecules(_records(), value_col="value", assumed_sigma=0.64)
    with pytest.raises(ValueError, match="cannot be used together"):
        _prepare(_records(), value_unit="nM", unit_column="batch")


def test_default_units_and_sigma_do_not_add_output_fields():
    default = _prepare(_records(), document_column="batch")
    explicit = _prepare(_records(), document_column="batch", value_unit="p", unit_column=None, assumed_sigma=None)
    pd.testing.assert_frame_equal(default[0], explicit[0])
    pd.testing.assert_frame_equal(default[1], explicit[1])
    assert default[2:] == explicit[2:]
    assert "raw_activity_value" not in default[1] and "sigma_doc" not in default[0]
    assert "assumed_sigma" not in default[2] and "excluded_reasons" not in default[3]
    no_document = _prepare(_records(), confidence=False)
    assert no_document[2] is None and "class_confidence" not in no_document[0]
    assert "--assumed-sigma" in no_document[3]["warnings"][0]
