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
):
    return {
        "Molecule ChEMBL ID": molecule_id,
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
        self.assertEqual(result.report["repeated_consistent_molecules"], 1)
        self.assertEqual(result.report["repeated_conflicted_molecules"], 0)

        ethanol = result.molecules[
            result.molecules["canonical_smiles"].eq("CCO")
        ].iloc[0]
        self.assertAlmostEqual(ethanol["activity_pchembl"], 5.9)
        self.assertEqual(ethanol["activity_class"], "inactive")
        self.assertEqual(ethanol["n_records"], 2)
        self.assertEqual(ethanol["repeat_status"], "repeated_consistent")

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


if __name__ == "__main__":
    unittest.main()
