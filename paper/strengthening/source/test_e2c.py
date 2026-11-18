"""Focused checks for the E2c implementation."""
from __future__ import annotations

from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[3]
S = ROOT / "paper/strengthening"
sys.path[:0] = [
    str(S / "source"),
    str(ROOT / "paper"),
    str(ROOT / "src"),
]

import numpy as np
import pandas as pd
import pytest

from chemgridmap.reliability import grid_neighbor_sets, representation_neighbor_sets
from e1b_neighbors import neighborhood_metrics, spearman_rho
from e2c_analysis import layout_metrics, subset_positions
from reliability_analysis import OUTPUT, layout_points, load_target


def test_subset_selection_depends_only_on_identifiers():
    identifiers = np.array([f"CHEMBL{index}" for index in range(1300)])
    first = subset_positions(identifiers, limit=1200, seed=20261102)
    second = subset_positions(identifiers.copy(), limit=1200, seed=20261102)
    assert len(first) == 1200
    assert np.array_equal(first, second)
    assert np.array_equal(first, np.sort(first))
    # The helper's input type is identifiers only, so activity cannot affect selection.
    assert first.dtype.kind in {"i", "u"}


def test_seed42_chembl205_refined_rho_matches_e1b():
    if not (S / "output/E1B/neighbor_sar_summary.csv").is_file():
        pytest.skip("Run paper/strengthening/source/e1b_neighbors.py first.")
    target = "chembl205"
    data, fingerprints = load_target(target)
    frame = pd.read_csv(OUTPUT / "layouts" / target / "coordinates_seed_42.csv")
    assert frame.molecule_identity_key.tolist() == data.molecule_identity_key.tolist()
    points = layout_points(frame, "neighbor_refined")
    observed = layout_metrics(
        data.activity_pchembl.to_numpy(float),
        fingerprints,
        points,
    )
    expected_metrics = neighborhood_metrics(
        data.activity_pchembl.to_numpy(float),
        representation_neighbor_sets(fingerprints, k=10),
        grid_neighbor_sets(points, k=10),
    )
    observed_rho = spearman_rho(observed.recall, observed.excess_error)
    expected_rho = spearman_rho(expected_metrics.recall, expected_metrics.excess_error)
    summary = pd.read_csv(S / "output/E1B/neighbor_sar_summary.csv")
    recorded_rho = float(
        summary.loc[
            (summary.target == target) & (summary.layout == "neighbor_refined"),
            "rho_seed42",
        ].iloc[0]
    )
    assert abs(observed_rho - expected_rho) <= 1e-12
    assert abs(observed_rho - recorded_rho) <= 1e-12
