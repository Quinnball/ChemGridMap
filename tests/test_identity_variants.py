import tempfile
import unittest
from pathlib import Path

import pandas as pd

from chemgridmap import __version__, curate_chembl_activity_data
from chemgridmap.chembl import save_chembl_curation
from chemgridmap.provenance import environment_versions
from chemgridmap.structures import parent_standard_identity


def chembl_row(parent_id, smiles, pchembl=6.5):
    return {
        "Molecule ChEMBL ID": f"{parent_id}_record",
        "Parent Molecule ChEMBL ID": parent_id,
        "Smiles": smiles,
        "pChEMBL Value": pchembl,
        "Target ChEMBL ID": "CHEMBL_TEST",
        "Standard Type": "IC50",
        "Standard Relation": "=",
        "Standard Units": "nM",
        "Potential Duplicate": 0,
    }


class IdentityVariantTests(unittest.TestCase):
    def assert_same_key_variant(self, parent_id, smiles_a, smiles_b):
        self.assertEqual(
            parent_standard_identity(smiles_a)[1],
            parent_standard_identity(smiles_b)[1],
        )
        result = curate_chembl_activity_data(
            pd.DataFrame([
                chembl_row(parent_id, smiles_a, 6.4),
                chembl_row(parent_id, smiles_b, 6.6),
            ]),
            target_id="CHEMBL_TEST",
        )
        self.assertEqual(len(result.molecules), 1)
        self.assertEqual(len(result.identity_conflicts), 0)
        variants = result.structure_variants.loc[
            result.structure_variants.molecule_identity_key == parent_id
        ]
        self.assertEqual(len(variants), 2)
        self.assertEqual(variants.standard_inchikey.nunique(), 1)

    def test_reviewed_parent_variants_are_not_conflicts(self):
        cases = {
            "CHEMBL429910": (
                "[2H][C@](F)(Cl)Br",
                "F[C@@H](Cl)Br",
            ),
            "CHEMBL589305": (
                "[2H][C@](F)(Cl)Br",
                "F[C@@H](Cl)Br",
            ),
            "CHEMBL191083": (
                "CN(C)c1ccc2nc3ccc(=[N+](C)C)cc-3sc2c1",
                "CN(C)c1ccc2nc3ccc(N(C)C)cc3[s+]c2c1",
            ),
            "CHEMBL4467360": (
                "[2H][C@](F)(Cl)Br",
                "F[C@@H](Cl)Br",
            ),
            "CHEMBL24828": (
                "O=c1cccc[nH]1",
                "Oc1ccccn1",
            ),
            "CHEMBL498770": (
                "c1ccc2[nH]ncc2c1",
                "c1ccc2n[nH]cc2c1",
            ),
            "CHEMBL8565": (
                "O=c1cc[nH]c(=O)[nH]1",
                "O=c1nc(O)cc[nH]1",
            ),
            "CHEMBL506871": (
                "O=c1cccc[nH]1",
                "Oc1ccccn1",
            ),
        }
        for parent_id, smiles in cases.items():
            with self.subTest(parent_id=parent_id):
                self.assert_same_key_variant(parent_id, *smiles)

    def test_different_constitutions_remain_conflicts(self):
        data = pd.DataFrame([
            chembl_row("CHEMBL_CONFLICT", "CCO", 6.0),
            chembl_row("CHEMBL_CONFLICT", "c1ccccc1", 7.0),
        ])
        keys = {
            parent_standard_identity("CCO")[1],
            parent_standard_identity("c1ccccc1")[1],
        }
        with self.assertRaises(ValueError) as context:
            curate_chembl_activity_data(data, identity_conflicts="error")
        message = str(context.exception)
        self.assertIn("CHEMBL_CONFLICT", message)
        for key in keys:
            self.assertIn(key, message)

    def test_different_stereochemistry_remains_a_conflict(self):
        data = pd.DataFrame([
            chembl_row("CHEMBL_STEREO", "F[C@H](Cl)Br", 6.0),
            chembl_row("CHEMBL_STEREO", "F[C@@H](Cl)Br", 7.0),
        ])
        with self.assertRaises(ValueError):
            curate_chembl_activity_data(data, identity_conflicts="error")

    def test_variant_csv_is_written_only_when_present(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            variant = curate_chembl_activity_data(
                pd.DataFrame([
                    chembl_row("CHEMBL_VARIANT", "O=c1cccc[nH]1"),
                    chembl_row("CHEMBL_VARIANT", "Oc1ccccn1"),
                ]),
                target_id="CHEMBL_TEST",
            )
            paths = save_chembl_curation(variant, output, "variant")
            self.assertIn("structure_variants", paths)
            self.assertTrue(Path(paths["structure_variants"]).is_file())

            unchanged = curate_chembl_activity_data(
                pd.DataFrame([chembl_row("CHEMBL_ONE", "CCO")]),
                target_id="CHEMBL_TEST",
            )
            paths = save_chembl_curation(unchanged, output, "unchanged")
            self.assertNotIn("structure_variants", paths)
            self.assertFalse(
                (output / "unchanged_structure_variants.csv").exists()
            )

    def test_manifest_uses_package_version(self):
        self.assertEqual(
            environment_versions()["packages"]["chemgridmap"],
            __version__,
        )


if __name__ == "__main__":
    unittest.main()
