import unittest

import numpy as np
import pandas as pd

from chemgridmap.core import (
    _tie_inclusive_neighbor_indices,
    assign_to_grid,
    canonicalize_smiles,
    choose_assignment_method,
    make_candidate_grid,
    mean_knn_absolute_difference,
    normalize_coordinates,
    prepare_molecules,
)


class CoreTests(unittest.TestCase):
    def test_canonicalize_smiles(self):
        self.assertEqual(canonicalize_smiles("C(C)O"), "CCO")
        self.assertIsNone(canonicalize_smiles("not-a-smiles"))

    def test_repeated_structures_raise_by_default(self):
        data = pd.DataFrame({"smiles": ["CCO", "OCC"]})
        with self.assertRaisesRegex(ValueError, "repeated canonical structures"):
            prepare_molecules(data)

    def test_explicit_first_policy_is_reported(self):
        data = pd.DataFrame({"smiles": ["CCO", "OCC", "c1ccccc1"]})
        prepared, report = prepare_molecules(
            data, duplicate_policy="first"
        )
        self.assertEqual(len(prepared), 2)
        self.assertEqual(report["duplicate_structures"], 1)

    def test_grid_cells_are_unique(self):
        points = np.array(
            [[0.0, 0.0], [0.2, 0.1], [0.8, 0.9], [1.0, 1.0]],
            dtype=float,
        )
        assigned, _, _ = assign_to_grid(points, padding=2)
        self.assertEqual(len(np.unique(assigned, axis=0)), len(points))

    def test_isotropic_scaling_preserves_axis_ratio(self):
        points = np.array([[0.0, 0.0], [10.0, 0.0], [0.0, 1.0]])
        isotropic = normalize_coordinates(points)
        independent = normalize_coordinates(points, mode="independent")
        self.assertAlmostEqual(np.ptp(isotropic[:, 0]), 1.0)
        self.assertAlmostEqual(np.ptp(isotropic[:, 1]), 0.1)
        self.assertAlmostEqual(np.ptp(independent[:, 1]), 1.0)

    def test_target_occupancy_controls_grid_size(self):
        grid, rows, cols = make_candidate_grid(892, target_occupancy=0.40)
        self.assertEqual(rows, 48)
        self.assertEqual(cols, 48)
        self.assertEqual(len(grid), rows * cols)
        self.assertLessEqual(892 / len(grid), 0.40)

    def test_invalid_target_occupancy_raises(self):
        with self.assertRaisesRegex(ValueError, "target_occupancy"):
            make_candidate_grid(10, target_occupancy=0.0)

    def test_zero_padding_still_provides_enough_cells(self):
        points = np.column_stack(
            [np.linspace(0.0, 1.0, 892), np.linspace(1.0, 0.0, 892)]
        )
        assigned, rows, cols = assign_to_grid(points, padding=0)
        self.assertGreaterEqual(rows * cols, len(points))
        self.assertEqual(len(np.unique(assigned, axis=0)), len(points))

    def test_sparse_grid_cells_are_unique_for_coincident_points(self):
        points = np.zeros((64, 2), dtype=float)
        assigned, _, _ = assign_to_grid(
            points,
            padding=4,
            method="sparse",
            sparse_neighbors=4,
        )
        self.assertEqual(len(np.unique(assigned, axis=0)), len(points))

    def test_auto_assignment_switches_by_pair_count(self):
        self.assertEqual(
            choose_assignment_method(20, padding=2, dense_max_pairs=10_000),
            "dense",
        )
        self.assertEqual(
            choose_assignment_method(200, padding=20, dense_max_pairs=10_000),
            "sparse",
        )

    def test_tie_inclusive_neighbors_keep_complete_grid_shell(self):
        points = np.array(
            [
                [0.0, 0.0],
                [1.0, 0.0],
                [-1.0, 0.0],
                [0.0, 1.0],
                [0.0, -1.0],
            ]
        )
        neighborhoods = _tie_inclusive_neighbor_indices(points, k=1)
        self.assertEqual(set(neighborhoods[0]), {1, 2, 3, 4})

    def test_jaccard_representation_metric_is_supported(self):
        fingerprints = np.array(
            [[1, 0, 0], [1, 1, 0], [0, 0, 1]],
            dtype=bool,
        )
        values = np.array([1.0, 2.0, 5.0])
        score = mean_knn_absolute_difference(
            fingerprints,
            values,
            k=1,
            metric="jaccard",
        )
        self.assertTrue(np.isfinite(score))


if __name__ == "__main__":
    unittest.main()
