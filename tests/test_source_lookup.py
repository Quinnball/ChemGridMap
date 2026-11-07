"""Local structure lookup only enriches possible sources, never map inputs."""

from collections import Counter
from pathlib import Path
import shutil

import pandas as pd
import pytest

from chemgridmap.chembl import curate_chembl_activity_data
from chemgridmap.evidence import prepare_generic_activity_records
from chemgridmap import inspection
from chemgridmap.inspection import (
    SOURCE_LOOKUP_COLUMNS_ERROR, match_source_lookup, validate_source_lookup,
)


def measured(smiles, values=None, identities=None):
    identities = identities or [f"CMP-{i:04d}" for i in range(1, len(smiles) + 1)]
    molecules = pd.DataFrame({"molecule_identity_key": identities, "canonical_smiles": smiles})
    records = pd.DataFrame({"molecule_identity_key": identities,
                            "activity_pchembl_raw": values or [7.] * len(smiles)})
    return molecules, records


@pytest.mark.parametrize("columns", [{"Smiles": ["CCO"]}, {"Document ChEMBL ID": ["CHEMBL1"]}, {}])
def test_lookup_requires_smiles_and_document_columns(columns):
    with pytest.raises(ValueError) as error:
        validate_source_lookup(pd.DataFrame(columns))
    assert str(error.value) == SOURCE_LOOKUP_COLUMNS_ERROR


def test_chembl_column_aliases_are_detected():
    assert validate_source_lookup(pd.DataFrame(columns=["Smiles", "Document ChEMBL ID"]))["smiles"] == "Smiles"
    assert validate_source_lookup(pd.DataFrame(columns=["canonical_smiles", "document_chembl_id"]))["document_id"] == "document_chembl_id"


def test_salts_match_parents_but_stereoisomers_do_not():
    molecules, records = measured(["Oc1ccccc1", "N[C@@H](C)C(=O)O", "CCO"])
    lookup = pd.DataFrame({
        "smiles": ["Oc1ccccc1.Cl", "N[C@H](C)C(=O)O", "not-a-smiles", "CCO"],
        "document_chembl_id": ["CHEMBL1", "CHEMBL2", "CHEMBL3", "javascript:alert(1)"],
        "pchembl_value": [7., 7., 7., 7.],
    })
    result, summary = match_source_lookup(molecules, records, lookup)
    assert summary == {"found": 1, "total": 3, "same": 1, "target_count": 0}
    assert result["CMP-0001"][0]["same_values"] == ["7.00"]
    assert result["CMP-0002"] == result["CMP-0003"] == []
    assert result["CMP-0001"][0]["links"][-1]["url"].endswith("/CHEMBL1")


@pytest.mark.parametrize("identifier", ["chembl1", "CHEMBL1/evil", "CHEMBL1?x=1", "CHEMBL", "", None])
def test_invalid_document_ids_are_not_matched(identifier):
    result, summary = match_source_lookup(*measured(["CCO"]), pd.DataFrame({
        "smiles": ["CCO"], "document_chembl_id": [identifier],
    }))
    assert result == {"CMP-0001": []}
    assert summary["found"] == 0


def test_values_use_inclusive_tolerance_and_document_sorting(monkeypatch):
    monkeypatch.setattr(inspection, "load_documents", lambda: {
        "CHEMBL1": {"document_chembl_id": "CHEMBL1", "year": "2000", "doi": "javascript:alert(1)"},
        "CHEMBL2": {"document_chembl_id": "CHEMBL2", "year": "2002"},
        "CHEMBL3": {"document_chembl_id": "CHEMBL3", "year": "2002"},
    })
    molecules, records = measured(["CCO"], [6.995])
    lookup = pd.DataFrame({
        "smiles": ["CCO"] * 6,
        "document_chembl_id": ["CHEMBL1", "CHEMBL3", "CHEMBL2", "CHEMBL2", "CHEMBL2", "CHEMBL9"],
        "pchembl_value": [6.9899, 7., 6.99, 7., 7., 6.995],
    })
    result, summary = match_source_lookup(molecules, records, lookup)
    sources = result["CMP-0001"]
    assert [source["id"] for source in sources] == ["CHEMBL2", "CHEMBL3", "CHEMBL9", "CHEMBL1"]
    assert sources[0]["same_values"] == ["6.99", "7.00"]
    assert sources[0]["n_records"] == 3
    assert sources[-1]["same_values"] == []
    assert summary["same"] == 1
    assert all(link["url"].startswith("https://") for source in sources for link in source["links"])


def test_all_retained_values_are_used_without_target_or_activity_filtering():
    molecules, records = measured(["CCO"], [6.])
    records = pd.concat([records, pd.DataFrame({"molecule_identity_key": ["CMP-0001"],
                                               "activity_pchembl_raw": [8.]})], ignore_index=True)
    lookup = pd.DataFrame({
        "smiles": ["CCO"] * 3, "document_chembl_id": ["CHEMBL1", "CHEMBL2", "CHEMBL3"],
        "pchembl_value": [6., 7., 8.], "target_chembl_id": ["CHEMBL205", "CHEMBL999", ""],
        "standard_type": ["IC50", "Ki", "unknown"],
    })
    result, summary = match_source_lookup(molecules, records, lookup)
    assert len(result["CMP-0001"]) == 3
    assert {source["id"] for source in result["CMP-0001"] if source["same_values"]} == {"CHEMBL1", "CHEMBL3"}
    assert summary["target_count"] == 2


@pytest.mark.parametrize("lookup_values", [None, [None, "bad", float("inf")]])
def test_missing_or_invalid_lookup_pchembl_does_not_mark_same_value(lookup_values):
    lookup = pd.DataFrame({"smiles": ["CCO"] * 3, "document_chembl_id": ["CHEMBL1"] * 3})
    if lookup_values is not None:
        lookup["pchembl_value"] = lookup_values
    result, summary = match_source_lookup(*measured(["CCO"]), lookup)
    assert summary["found"] == 1 and summary["same"] == 0
    assert result["CMP-0001"][0]["n_records"] == 3
    assert result["CMP-0001"][0]["same_values"] == []


def test_raw_pchembl_is_preferred_and_legacy_pchembl_is_supported():
    molecules, records = measured(["CCO"], [7.])
    lookup = pd.DataFrame({"smiles": ["CCO"], "document_chembl_id": ["CHEMBL1"], "pchembl_value": [7.]})
    records["pchembl_value"] = 8.
    assert match_source_lookup(molecules, records, lookup)[1]["same"] == 1
    assert match_source_lookup(molecules, records.drop(columns="activity_pchembl_raw"), lookup)[1]["same"] == 0


def test_unique_smiles_are_normalized_once_across_both_inputs(monkeypatch):
    canonical_calls, parent_calls = Counter(), Counter()
    canonicalize, parent = inspection.canonicalize_smiles, inspection.parent_structure

    def counted_canonical(smiles):
        canonical_calls[smiles] += 1
        return canonicalize(smiles)

    def counted_parent(smiles):
        parent_calls[smiles] += 1
        return parent(smiles)

    monkeypatch.setattr(inspection, "canonicalize_smiles", counted_canonical)
    monkeypatch.setattr(inspection, "parent_structure", counted_parent)
    lookup = pd.DataFrame({"smiles": ["Oc1ccccc1", "Oc1ccccc1", "c1ccccc1O", "Oc1ccccc1.Cl", "invalid", "invalid"],
                           "document_chembl_id": ["CHEMBL1"] * 6})
    molecules, records = measured(["Oc1ccccc1", "Oc1ccccc1.Cl"])
    lookup_original, molecules_original, records_original = lookup.copy(deep=True), molecules.copy(deep=True), records.copy(deep=True)
    result, summary = match_source_lookup(molecules, records, lookup)
    assert summary["found"] == 2
    assert canonical_calls == {"Oc1ccccc1": 1, "c1ccccc1O": 1, "Oc1ccccc1.Cl": 1, "invalid": 1}
    assert parent_calls == {"Oc1ccccc1": 1, "Cl.Oc1ccccc1": 1}
    assert result["CMP-0001"][0]["n_records"] == 4
    pd.testing.assert_frame_equal(lookup, lookup_original)
    pd.testing.assert_frame_equal(molecules, molecules_original)
    pd.testing.assert_frame_equal(records, records_original)


@pytest.fixture(scope="module")
def caii_lookup(tmp_path_factory):
    source = Path(inspection.__file__).with_name("data") / "chembl205_ic50.csv"
    copied = tmp_path_factory.mktemp("source_lookup") / source.name
    shutil.copy2(source, copied)
    lookup = pd.read_csv(copied)
    curated = curate_chembl_activity_data(lookup, target_id="CHEMBL205")
    return lookup, curated


def test_bundled_caii_covers_retained_sources_and_additional_matching_documents(caii_lookup):
    lookup, curated = caii_lookup
    result, summary = match_source_lookup(curated.molecules, curated.retained_records, lookup)
    assert summary == {"found": 889, "total": 889, "same": 889, "target_count": 1}
    document_column = curated.report["detected_columns"]["document_id"]
    chembl20 = curated.retained_records[curated.retained_records.molecule_identity_key.eq("CHEMBL20")]
    assert chembl20[document_column].nunique() == 27
    # Lookup accepts every valid structure/document row, including rows excluded
    # from the map's activity curation: acetazolamide therefore has 49 sources.
    assert len(result["CHEMBL20"]) == 49
    for identity, group in curated.retained_records.groupby("molecule_identity_key"):
        assert set(group[document_column]) <= {source["id"] for source in result[identity]}


def test_bundled_caii_matching_does_not_depend_on_compound_ids(caii_lookup):
    lookup, curated = caii_lookup
    expected, expected_summary = match_source_lookup(curated.molecules, curated.retained_records, lookup)
    molecules, records = curated.molecules.copy(), curated.retained_records.copy()
    renamed = {key: f"CMP-{i:04d}" for i, key in enumerate(molecules.molecule_identity_key, 1)}
    for frame in (molecules, records):
        frame["molecule_identity_key"] = frame.molecule_identity_key.map(renamed)
    observed, summary = match_source_lookup(molecules, records, lookup)
    assert summary == expected_summary
    assert all(observed[renamed[key]] == sources for key, sources in expected.items())


def test_bundled_caii_one_median_value_per_molecule(caii_lookup):
    lookup, curated = caii_lookup
    medians = curated.retained_records.groupby("molecule_identity_key").activity_pchembl_raw.median()
    frame = curated.molecules[["molecule_identity_key", "canonical_smiles"]].copy()
    frame["IC50_nM"] = frame.molecule_identity_key.map(lambda key: 10 ** (9 - medians[key]))
    molecules, records, _, _ = prepare_generic_activity_records(
        frame, compound_column="molecule_identity_key", smiles_column="canonical_smiles",
        value_column="IC50_nM", value_unit="nM", assumed_sigma=.64,
    )
    result, summary = match_source_lookup(molecules, records, lookup)
    assert summary == {"found": 889, "total": 889, "same": 830, "target_count": 1}
    columns = validate_source_lookup(lookup)
    parent_by_smiles = {}
    for smiles in lookup[columns["smiles"]].dropna().unique():
        canonical = inspection.canonicalize_smiles(smiles)
        parent_by_smiles[smiles] = inspection.parent_structure(canonical)[0] if canonical else None
    groups = lookup.assign(parent=lookup[columns["smiles"]].map(parent_by_smiles)).groupby(
        ["parent", columns["document_id"]])
    for row in molecules.itertuples(index=False):
        parent = inspection.parent_structure(row.canonical_smiles)[0]
        for source in result[row.molecule_identity_key]:
            values = pd.to_numeric(groups.get_group((parent, source["id"]))[columns["pchembl"]], errors="coerce")
            expected = sorted({f"{value:.2f}" for value in values
                               if abs(value - medians[row.molecule_identity_key]) <= .005 + 1e-12}, key=float)
            assert source["same_values"] == expected
