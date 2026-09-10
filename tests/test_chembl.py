import unittest

import pandas as pd

from chemgridmap import curate_chembl_activity_data


def chembl_row(
    molecule_id,
    smiles,
    pchembl,
    *,
    target_id="CHEMBL205",
    standard_type="IC50",
    relation="=",
    units="nM",
    validity=None,
    potential_duplicate=0,
    parent_molecule_id=None,
):
    return {
        "Molecule ChEMBL ID": molecule_id,
        "Parent Molecule ChEMBL ID": parent_molecule_id,
        "Smiles": smiles,
        "pChEMBL Value": pchembl,
        "Target ChEMBL ID": target_id,
        "Target Name": "Carbonic anhydrase 2",
        "Standard Type": standard_type,
        "Standard Relation": relation,
        "Standard Units": units,
        "Data Validity Comment": validity,
        "Potential Duplicate": potential_duplicate,
        "Assay ChEMBL ID": "CHEMBL_A1",
        "Document ChEMBL ID": "CHEMBL_D1",
    }


class ChemblCurationTests(unittest.TestCase):
    def test_web_export_quoted_relations_preserve_exact_filter(self):
        data = pd.DataFrame([
            chembl_row("CHEMBL1", "CCO", 7.0, relation="'='"),
            chembl_row("CHEMBL2", "CCN", 6.0, relation="'>'"),
            chembl_row("CHEMBL3", "CCC", 6.5, relation="'<='"),
        ])
        result = curate_chembl_activity_data(data, target_id="CHEMBL205")
        self.assertEqual(len(result.retained_records), 1)
        self.assertEqual(result.retained_records.iloc[0]["Standard Relation"], "'='")
        self.assertEqual(len(result.excluded_records), 2)

    def test_standard_export_is_filtered_and_aggregated(self):
        data = pd.DataFrame(
            [
                chembl_row("CHEMBL1", "CCO", 5.8),
                chembl_row("CHEMBL1", "OCC", 6.0),
                chembl_row("CHEMBL2", "c1ccccc1", 7.5),
                chembl_row("CHEMBL3", "Cc1ccccc1", 6.5),
                chembl_row(
                    "CHEMBL4",
                    "Oc1ccccc1",
                    8.0,
                    potential_duplicate=1,
                ),
                chembl_row("CHEMBL5", "not-a-smiles", 5.0),
                chembl_row(
                    "CHEMBL6",
                    "n1ccccc1",
                    7.0,
                    target_id="CHEMBL999",
                ),
                chembl_row(
                    "CHEMBL7",
                    "Nc1ccccc1",
                    6.2,
                    validity="Potential author error",
                ),
                chembl_row(
                    "CHEMBL8",
                    "CC(=O)Oc1ccccc1C(=O)O",
                    6.9,
                    relation=">",
                ),
            ]
        )

        result = curate_chembl_activity_data(
            data,
            target_id="CHEMBL205",
        )

        self.assertEqual(result.report["retained_activity_records"], 4)
        self.assertEqual(result.report["excluded_activity_records"], 5)
        self.assertEqual(result.report["unique_molecules"], 3)
        self.assertEqual(result.report["repeated_molecules"], 1)
        self.assertEqual(result.report["repeated_unflagged_molecules"], 1)
        self.assertEqual(result.report["repeated_conflicted_molecules"], 0)

        ethanol = result.molecules[
            result.molecules["canonical_smiles"].eq("CCO")
        ].iloc[0]
        self.assertAlmostEqual(ethanol["activity_pchembl"], 5.9)
        self.assertEqual(ethanol["activity_class"], "inactive")
        self.assertEqual(ethanol["n_records"], 2)
        self.assertEqual(ethanol["repeat_status"], "repeated_unflagged")

        reasons = result.excluded_records["exclusion_reason"].tolist()
        self.assertTrue(
            any("chembl_potential_duplicate" in reason for reason in reasons)
        )
        self.assertTrue(any("invalid_or_missing_smiles" in reason for reason in reasons))
        self.assertTrue(any("target_mismatch" in reason for reason in reasons))
        self.assertTrue(any("data_validity_flag" in reason for reason in reasons))
        self.assertTrue(
            any("non_exact_standard_relation" in reason for reason in reasons)
        )

    def test_missing_optional_audit_columns_produce_warnings(self):
        data = pd.DataFrame(
            {
                "Compound": ["CHEMBL1", "CHEMBL2"],
                "smiles": ["CCO", "c1ccccc1"],
                "pchembl": [5.5, 7.5],
            }
        )
        result = curate_chembl_activity_data(
            data,
            target_id="CHEMBL205",
        )
        self.assertEqual(result.report["unique_molecules"], 2)
        self.assertGreaterEqual(len(result.report["warnings"]), 1)

    def test_multiple_targets_require_an_explicit_selection(self):
        data = pd.DataFrame(
            [
                chembl_row("CHEMBL1", "CCO", 5.5),
                chembl_row(
                    "CHEMBL2",
                    "c1ccccc1",
                    7.5,
                    target_id="CHEMBL999",
                ),
            ]
        )
        with self.assertRaisesRegex(ValueError, "Multiple Target ChEMBL IDs"):
            curate_chembl_activity_data(data)

    def test_parent_molecule_identity_groups_salt_records(self):
        data = pd.DataFrame(
            [
                chembl_row(
                    "CHEMBL_SALT1",
                    "CCN.Cl",
                    6.2,
                    parent_molecule_id="CHEMBL_PARENT1",
                ),
                chembl_row(
                    "CHEMBL_SALT2",
                    "CCN",
                    6.4,
                    parent_molecule_id="CHEMBL_PARENT1",
                ),
            ]
        )
        parent_result = curate_chembl_activity_data(data)
        smiles_result = curate_chembl_activity_data(
            data,
            molecule_identity="canonical-smiles",
        )
        self.assertEqual(parent_result.report["unique_molecules"], 1)
        self.assertEqual(smiles_result.report["unique_molecules"], 2)
        self.assertEqual(
            parent_result.molecules.loc[0, "molecule_identity_key"],
            "CHEMBL_PARENT1",
        )
        self.assertAlmostEqual(
            parent_result.molecules.loc[0, "activity_pchembl"],
            6.3,
        )
        self.assertEqual(parent_result.molecules.loc[0, "canonical_smiles"], "CCN")
        self.assertIn("CCN.Cl", parent_result.molecules.loc[0, "record_canonical_smiles"])

    def test_unflagged_does_not_mean_identical_display_labels(self):
        data = pd.DataFrame([chembl_row("CHEMBL1", "CCN", 6.74),
                             chembl_row("CHEMBL1", "CCN", 7.68)])
        result = curate_chembl_activity_data(data)
        row = result.molecules.iloc[0]
        self.assertEqual(row.repeat_status, "repeated_unflagged")
        self.assertEqual(row.has_class_boundary_crossing, 1)
        self.assertEqual(row.crosses_clean_boundary, 0)
        self.assertEqual(result.report["unflagged_class_boundary_crossing_molecules"], 1)

    def test_conflict_keeps_median_instead_of_forcing_medium(self):
        data = pd.DataFrame([chembl_row("CHEMBL1", "CCN", x) for x in [5., 7.2, 8.]])
        row = curate_chembl_activity_data(data).molecules.iloc[0]
        self.assertEqual(row.repeat_status, "repeated_conflicted")
        self.assertEqual(row.activity_class, "active")
        self.assertAlmostEqual(row.pchembl_median, 7.2)

    def test_missing_parent_id_is_recovered_only_unambiguously(self):
        data = pd.DataFrame([chembl_row("CHEMBL1", "CCN.Cl", 6.4, parent_molecule_id="P1"),
                             chembl_row("CHEMBL2", "CCN", 6.6)])
        result = curate_chembl_activity_data(data)
        self.assertEqual(len(result.molecules), 1)
        self.assertEqual(result.molecules.iloc[0].canonical_smiles, "CCN")
        self.assertIn("unambiguous_structure_parent_lookup", result.molecules.iloc[0].molecule_identity_source)

    def test_disagreeing_parent_structures_fail_instead_of_guessing(self):
        data = pd.DataFrame([chembl_row("CHEMBL1", "CCN", 6., parent_molecule_id="P1"),
                             chembl_row("CHEMBL2", "c1ccccc1", 7., parent_molecule_id="P1")])
        with self.assertRaisesRegex(ValueError, "multiple normalized parent structures"):
            curate_chembl_activity_data(data)

    def test_stereochemistry_is_retained(self):
        data = pd.DataFrame([chembl_row("CHEMBL1", "N[C@H](C)Cc1ccccc1.Cl", 6.),
                             chembl_row("CHEMBL2", "N[C@@H](C)Cc1ccccc1.Cl", 7.)])
        result = curate_chembl_activity_data(data)
        self.assertEqual(len(result.molecules), 2)
        self.assertTrue(all("@" in x for x in result.molecules.canonical_smiles))


if __name__ == "__main__":
    unittest.main()
