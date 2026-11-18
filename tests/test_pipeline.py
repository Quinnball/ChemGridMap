import tempfile
import unittest
from inspect import signature
from pathlib import Path

import pandas as pd
import pytest

from chemgridmap import build_grid_map, pipeline


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
    def test_library_defaults_match_the_0_5_public_defaults(self):
        parameters = signature(build_grid_map).parameters
        self.assertEqual(parameters["projection"].default, "umap")
        self.assertEqual(parameters["grid_refinement"].default, "neighbors")
        self.assertEqual(parameters["core_rule"].default, "coverage")
        self.assertEqual(parameters["fmcs_timeout"].default, 1)
        self.assertEqual(parameters["evidence_model"].default, "gaussian")

    def test_library_generic_records_compute_support_by_default(self):
        rows = []
        for index, smiles in enumerate(["CCO", "c1ccccc1", "Cc1ccccc1", "Oc1ccccc1"]):
            for document, offset in [("D1", 0.0), ("D2", 0.4)]:
                rows.append({
                    "compound": f"C{index}",
                    "smiles": smiles,
                    "p_activity": 5.5 + index * 0.3 + offset,
                    "document": document,
                })
        result = build_grid_map(
            pd.DataFrame(rows),
            projection="pca",
            grouping="none",
            compound_col="compound",
            value_col="p_activity",
            document_col="document",
            grid_padding=2,
            output_formats=["svg"],
            render_detail="overview",
        )
        self.assertIn("class_confidence", result.data)
        self.assertIn("document_balanced_pchembl", result.data)

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
        self.assertIn("neighbor_kth_similarity", result.data)
        self.assertTrue(result.data["neighbor_kth_similarity"].between(0, 1).all())
        self.assertEqual(
            result.metrics.loc[0, "coordinate_scaling"],
            "isotropic",
        )
        self.assertEqual(
            result.metrics.loc[0, "grid_sizing_strategy"],
            "explicit_padding",
        )


def test_grouping_limit_is_checked_before_projection(monkeypatch):
    projected = []
    monkeypatch.setattr(pipeline, "project_2d", lambda *args, **kwargs: projected.append(True))
    large = pd.DataFrame({"smiles": ["CCO"] * 10_001, "value": [6.0] * 10_001})
    with pytest.raises(ValueError, match=r"10,001.*--grouping none"):
        build_grid_map(large, value_col="value", grouping="butina", duplicate_policy="keep")
    assert projected == []


if __name__ == "__main__":
    unittest.main()
