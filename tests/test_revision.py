import json
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd

from chemgridmap.core import (assign_to_grid, label_purity,
                             mean_knn_absolute_difference, tie_aware_knn_agreement)
from chemgridmap.cli import main


def test_tie_inclusive_annotations_are_invariant_to_row_order():
    xy = np.array([(x, y) for x in range(5) for y in range(5)], dtype=float)
    labels = np.arange(len(xy)) % 3
    values = np.arange(len(xy)) / 10
    expected = (label_purity(xy, labels), mean_knn_absolute_difference(xy, values))
    for seed in range(8):
        order = np.random.RandomState(seed).permutation(len(xy))
        assert np.allclose(expected, (label_purity(xy[order], labels[order]),
                                    mean_knn_absolute_difference(xy[order], values[order])))
        assert tie_aware_knn_agreement(xy[order], xy[order])["recall"] == 1


def test_sparse_manifest_reports_fallback_graph():
    xy = np.random.RandomState(3).normal(size=(40, 2))
    info = {}
    grid, _, _ = assign_to_grid(xy, method="sparse", sparse_neighbors=1, diagnostics=info,
                               sparse_adaptive=False)
    assert len(np.unique(grid, axis=0)) == len(xy)
    assert info["candidate_edges"] == 40 + info["fallback_edges_added"]
    assert info["max_candidates_per_molecule"] <= 2
    assert info["fallback_edges_used"] <= info["fallback_edges_added"]


def test_semicolon_web_export_runs_with_parent_fallback_and_manifest(tmp_path):
    data = pd.DataFrame({"Molecule ChEMBL ID": ["C1", "C2", "C3", "C4"],
                         "Smiles": ["CCN.Cl", "CCN", "c1ccccc1", "Cc1ccccc1"],
                         "pChEMBL Value": [6.8, 7.4, 5.4, 8.1],
                         "Standard Relation": ["="] * 4,
                         "Standard Units": ["nM"] * 4})
    source = tmp_path / "web.csv"
    data.to_csv(source, sep=";", index=False, encoding="utf-8-sig")
    out = tmp_path / "out"
    assert main([str(source), "--input-format", "chembl", "--output-dir", str(out),
                 "--output-formats", "svg"]) == 0
    molecules = pd.read_csv(out / "grid_map_molecule_level.csv")
    assert len(molecules) == 3
    assert "CCN" in molecules.canonical_smiles.values
    manifest = json.loads((out / "grid_map_run_manifest.json").read_text())
    assert manifest["source"]["input_file_sha256"]
    assert manifest["assignment"]["cost_dtype"] == "float64"
    root = ET.parse(out / "grid_map.svg").getroot()
    assert not root.findall(".//{http://www.w3.org/2000/svg}image")
    assert any("bond-" in x.get("class", "") for x in root.iter())
    tiles = [x for x in root.iter() if x.get("data-molecule-id") is not None]
    assert len(tiles) == len(molecules)


def test_adaptive_sparse_expands_and_checks_its_cost():
    xy = np.random.default_rng(9).normal(size=(80, 2))
    info = {}
    assigned, _, _ = assign_to_grid(xy, method="sparse", sparse_neighbors=2, diagnostics=info)
    dense, _, _ = assign_to_grid(xy, method="dense")
    assert info["converged"] and info["iterations"] > 1
    assert info["fallback_edges_used"] == 0
    from chemgridmap.core import normalize_coordinates
    points = normalize_coordinates(xy)
    assert np.isclose(np.square(points - assigned).sum(), np.square(points - dense).sum())
    costs = [step["squared_cost"] for step in info["candidate_history"]]
    assert np.all(np.diff(costs) <= 1e-10)
