import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from chemgridmap.cli import main


class CliTests(unittest.TestCase):
    def test_cli_smoke_run(self):
        with tempfile.TemporaryDirectory() as directory:
            input_path = Path(directory, "molecules.csv")
            pd.DataFrame(
                {
                    "smiles": [
                        "c1ccccc1",
                        "Cc1ccccc1",
                        "Oc1ccccc1",
                        "Nc1ccccc1",
                        "n1ccccc1",
                        "CC(=O)Oc1ccccc1C(=O)O",
                    ],
                    "value": [78.1, 92.1, 94.1, 93.1, 79.1, 180.2],
                    "group": [
                        "small",
                        "small",
                        "small",
                        "small",
                        "small",
                        "medium",
                    ],
                }
            ).to_csv(input_path, index=False)
            return_code = main(
                [
                    str(input_path),
                    "--value-col",
                    "value",
                    "--label-col",
                    "group",
                    "--projection",
                    "pca",
                    "--grid-padding",
                    "2",
                    "--output-dir",
                    directory,
                    "--name",
                    "test_map",
                ]
            )
            self.assertEqual(return_code, 0)
            self.assertTrue(
                Path(directory, "test_map_grid_coordinates.csv").is_file()
            )

    def test_chembl_csv_to_map_end_to_end(self):
        with tempfile.TemporaryDirectory() as directory:
            input_path = Path(directory, "chembl_export.csv")
            rows = []
            structures = [
                ("CHEMBL1", "CCO", 5.8),
                ("CHEMBL1", "OCC", 6.0),
                ("CHEMBL2", "c1ccccc1", 7.5),
                ("CHEMBL3", "Cc1ccccc1", 6.5),
                ("CHEMBL4", "Oc1ccccc1", 8.0),
                ("CHEMBL5", "Nc1ccccc1", 5.5),
                ("CHEMBL6", "n1ccccc1", 7.2),
                ("CHEMBL7", "CC(=O)Oc1ccccc1C(=O)O", 6.8),
            ]
            for molecule_id, smiles, pchembl in structures:
                rows.append(
                    {
                        "Molecule ChEMBL ID": molecule_id,
                        "Smiles": smiles,
                        "pChEMBL Value": pchembl,
                        "Target ChEMBL ID": "CHEMBL205",
                        "Standard Type": "IC50",
                        "Standard Relation": "=",
                        "Standard Units": "nM",
                        "Data Validity Comment": None,
                        "Potential Duplicate": 0,
                    }
                )
            rows.append(
                {
                    "Molecule ChEMBL ID": "CHEMBL8",
                    "Smiles": "CCN",
                    "pChEMBL Value": 6.1,
                    "Target ChEMBL ID": "CHEMBL205",
                    "Standard Type": "IC50",
                    "Standard Relation": "=",
                    "Standard Units": "nM",
                    "Data Validity Comment": None,
                    "Potential Duplicate": 1,
                }
            )
            pd.DataFrame(rows).to_csv(input_path, index=False)

            return_code = main(
                [
                    str(input_path),
                    "--input-format",
                    "chembl",
                    "--target-id",
                    "CHEMBL205",
                    "--projection",
                    "pca",
                    "--grid-padding",
                    "2",
                    "--output-dir",
                    directory,
                    "--name",
                    "chembl_map",
                ]
            )

            self.assertEqual(return_code, 0)
            expected = [
                "chembl_map_molecule_level.csv",
                "chembl_map_retained_records.csv",
                "chembl_map_excluded_records.csv",
                "chembl_map_curation_report.json",
                "chembl_map_grid_coordinates.csv",
                "chembl_map_metrics.csv",
                "chembl_map.svg",
                "chembl_map.png",
                "chembl_map.pdf",
            ]
            for filename in expected:
                self.assertTrue(Path(directory, filename).is_file(), filename)

            with Path(directory, "chembl_map_curation_report.json").open() as handle:
                report = json.load(handle)
            self.assertEqual(report["input_rows"], 9)
            self.assertEqual(report["retained_activity_records"], 8)
            self.assertEqual(report["unique_molecules"], 7)
            self.assertEqual(
                report["exclusion_counts"]["chembl_potential_duplicate"],
                1,
            )


if __name__ == "__main__":
    unittest.main()
