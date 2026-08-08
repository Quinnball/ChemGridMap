import unittest
from pathlib import Path

import pandas as pd


class ExampleDatasetTests(unittest.TestCase):
    def test_project_example_is_complete(self):
        path = (
            Path(__file__).resolve().parents[1]
            / "examples"
            / "chembl205_embedding_umap.csv"
        )
        data = pd.read_csv(path)
        embedding_columns = [
            column for column in data.columns if column.startswith("emb_")
        ]

        self.assertEqual(len(data), 770)
        self.assertEqual(data["canonical_smiles"].nunique(), 770)
        self.assertEqual(len(embedding_columns), 128)
        self.assertEqual(data["activity_pchembl"].isna().sum(), 0)
        self.assertEqual(
            data["activity_class"].value_counts().to_dict(),
            {"active": 475, "inactive": 165, "medium": 130},
        )


if __name__ == "__main__":
    unittest.main()

