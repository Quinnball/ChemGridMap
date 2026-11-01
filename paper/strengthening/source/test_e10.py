from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

os.environ.setdefault("NUMBA_NUM_THREADS", "1")

from e6_core import tie_weights
from e10_analysis import changed_e2c_datasets, component_summary
from chemgridmap.projection import project_2d


def test_disconnected_projection_repeats_with_same_seed():
    matrix = np.array([
        [0, 0], [0, 1], [0, 2], [0, 3],
        [100, 0], [100, 1], [100, 2], [100, 3],
    ], dtype=float)
    first = project_2d(
        matrix, method="umap", random_state=42, umap_metric="euclidean"
    )
    second = project_2d(
        matrix, method="umap", random_state=42, umap_metric="euclidean"
    )
    assert np.array_equal(first, second)


def test_fixed_k_weights_sum_to_k_with_ties():
    cells = np.array([
        [0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [-1.0, 0.0]
    ])
    probabilities, _ = tie_weights(cells, k=2)
    np.testing.assert_allclose(
        probabilities.sum(axis=1), 2.0, rtol=0.0, atol=1e-9
    )


def test_umap_graph_components_detects_disconnected_and_connected_inputs():
    disconnected = np.array([
        [1, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0],
        [1, 1, 0, 0, 0, 0],
        [1, 1, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0, 0],
        [0, 0, 0, 1, 1, 0],
        [0, 0, 0, 1, 1, 0],
    ], dtype=np.float32)
    count, _ = component_summary(disconnected)
    assert count > 1

    connected = np.zeros((24, 8), dtype=np.float32)
    for index in range(len(connected)):
        connected[index, index % 8] = 1
        connected[index, (index + 1) % 8] = 1
    count, _ = component_summary(connected)
    assert count == 1


def test_disconnected_graphs_match_changed_e2c_datasets():
    root = Path(__file__).resolve().parents[3]
    output = root / "paper/strengthening/output"
    graph_path = output / "E10/graph_components.csv"
    old_path = output / "as_executed_pre05/E2c/layout_by_target.csv"
    new_path = output / "E2c/layout_by_target.csv"
    if not all(path.is_file() for path in (graph_path, old_path, new_path)):
        return

    graph = pd.read_csv(graph_path)
    disconnected = {
        (str(target), str(endpoint))
        for target, endpoint in graph.loc[
            graph.components.gt(1), ["target", "endpoint"]
        ].itertuples(index=False, name=None)
    }
    old = pd.read_csv(old_path, float_precision="round_trip")
    new = pd.read_csv(new_path, float_precision="round_trip")
    changed, _ = changed_e2c_datasets(old, new)
    assert changed == disconnected
