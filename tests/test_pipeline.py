import tempfile
import unittest
from pathlib import Path

import pandas as pd

from chemgridmap import build_grid_map


MOLECULES = pd.DataFrame(
    {
        "name": [
            "benzene",
            "toluene",
            "phenol",
            "aniline",
            "pyridine",
            "aspirin",
        ],
        "smiles": [
            "c1ccccc1",
            "Cc1ccccc1",
            "Oc1ccccc1",
            "Nc1ccccc1",
            "n1ccccc1",
            "CC(=O)Oc1ccccc1C(=O)O",
        ],
        "value": [78.1, 92.1, 94.1, 93.1, 79.1, 180.2],
        "group": ["small", "small", "small", "small", "small", "medium"],
        "x": [0.0, 0.1, 0.2, 0.25, 0.7, 1.0],
        "y": [0.0, 0.1, 0.15, 0.2, 0.8, 1.0],
    }
)


class PipelineTests(unittest.TestCase):
    def test_existing_coordinates_route_and_exports(self):
        with tempfile.TemporaryDirectory() as directory:
            result = build_grid_map(
                MOLECULES,
                output_dir=directory,
                name="coordinates",
                value_col="value",
                label_col="group",
                x_col="x",
                y_col="y",
                grid_padding=2,
                k=2,
            )
            self.assertEqual(len(result.data), len(MOLECULES))
            self.assertEqual(
                result.data[["grid_col", "grid_row"]].drop_duplicates().shape[0],
                len(MOLECULES),
            )
            for suffix in ["svg", "png", "pdf"]:
                self.assertTrue(Path(result.output_files[suffix]).is_file())

    def test_morgan_pca_route(self):
        result = build_grid_map(
            MOLECULES,
            value_col="value",
            label_col="group",
            representation="morgan",
            projection="pca",
            grid_padding=2,
            k=2,
        )
        self.assertEqual(result.projection.shape, (len(MOLECULES), 2))
        self.assertIn(
            "projection_to_grid_trustworthiness", result.metrics.columns
        )
        self.assertIn("projection_distance_metric", result.metrics.columns)
        self.assertEqual(
            result.metrics.loc[0, "representation_distance_metric"],
            "jaccard",
        )
        self.assertIn(
            "projection_grid_tie_aware_recall", result.metrics.columns
        )
        self.assertIn(
            "projection_grid_tie_aware_precision", result.metrics.columns
        )
        self.assertIn("grid_knn_value_abs_diff", result.metrics.columns)
        self.assertEqual(
            result.metrics.loc[0, "coordinate_scaling"],
            "isotropic",
        )
        self.assertEqual(
            result.metrics.loc[0, "grid_sizing_strategy"],
            "explicit_padding",
        )


if __name__ == "__main__":
    unittest.main()
