import unittest

import numpy as np
import pandas as pd

from chemgridmap.core import (
    assign_to_grid,
    canonicalize_smiles,
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


if __name__ == "__main__":
    unittest.main()

