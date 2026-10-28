"""Independent checks for adjacent contrasts and CLI consistency validation."""
import json

import numpy as np
import pandas as pd
import pytest

from reliability_contrasts import (eight_connected_pairs, tie_aware_neighbor_pairs,
    contrast_statistics, compare_figure1_entries, valid_cli_cache)
from reliability_analysis import sha256


def test_eight_connected_occupied_pairs_match_all_pair_distances():
    cells = np.array([[2, 1], [0, 0], [1, 1], [2, 0], [5, 5], [1, 0]])
    expected = {(i, j) for i in range(len(cells)) for j in range(i + 1, len(cells))
                if np.abs(cells[i] - cells[j]).max() == 1}
    assert set(map(tuple, eight_connected_pairs(cells))) == expected
    assert eight_connected_pairs([[0, 0], [2, 2]]).shape == (0, 2)
    with pytest.raises(ValueError, match="unique"):
        eight_connected_pairs([[0, 0], [0, 0]])


def test_tie_aware_union_includes_shell_ties_and_one_way_edges():
    cells = np.array([[0, 0], [1, 0], [-1, 0], [0, 1], [0, -3]])
    # The centre has three equidistant nearest neighbours; the remote cell
    # contributes a directed nearest-neighbour edge absent in the reverse set.
    assert set(map(tuple, tie_aware_neighbor_pairs(cells, k=1))) == {
        (0, 1), (0, 2), (0, 3), (0, 4)}


def test_contrast_counts_and_expectations_use_consensus_probabilities():
    stats = contrast_statistics([[0, 1], [0, 2], [1, 2]], ["active", "inactive", "active"],
                                [[.1, 0, .9], [.9, 0, .1], [.2, 0, .8]])
    assert stats["neighbor_pairs"] == 3
    assert stats["displayed_contrasts"] == 2
    assert stats["supported_contrasts_p_ge_0.8"] == 1
    assert stats["expected_true_contrasts_among_displayed"] == pytest.approx(.82 + .74)
    assert stats["expected_hidden_contrasts_same_colour"] == pytest.approx(.26)
    assert contrast_statistics([], ["active"], [[0, 0, 1]])["neighbor_pairs"] == 0


def figure1_fixture():
    cli = pd.DataFrame({"molecule_identity_key": ["A", "B"], "grid_col": [1, 3], "grid_row": [0, 2],
        "class_confidence": [.9, .7], "p_consensus_inactive": [0, .7],
        "p_consensus_medium": [.1, .2], "p_consensus_active": [.9, .1],
        "class_status": ["supported", "uncertain"], "activity_class": ["active", "inactive"]})
    layout = cli[["molecule_identity_key", "grid_col", "grid_row"]].rename(
        columns={"grid_col": "refined_grid_col", "grid_row": "refined_grid_row"})
    return cli, layout, cli.drop(columns=["grid_col", "grid_row"])


def test_figure1_comparison_aligns_identity_and_detects_coordinates_and_support():
    cli, layout, annotations = figure1_fixture()
    _, checks = compare_figure1_entries(cli, layout.iloc[::-1], annotations.iloc[::-1])
    assert checks["passed"]
    cli.loc[0, "grid_col"] += 1
    cli.loc[1, "class_confidence"] += .01
    comparison, checks = compare_figure1_entries(cli, layout, annotations)
    assert not checks["passed"]
    assert checks["grid_mismatched_entries"] == 1
    assert checks["class_support_mismatched_entries"] == 1
    assert not comparison.coordinates_match.iloc[0]
    assert not comparison.class_support_match.iloc[1]


def test_figure1_comparison_rejects_missing_or_duplicate_entries():
    cli, layout, annotations = figure1_fixture()
    with pytest.raises(ValueError, match="identity sets"):
        compare_figure1_entries(cli.iloc[:1], layout, annotations)
    with pytest.raises(ValueError, match="duplicate"):
        compare_figure1_entries(pd.concat([cli, cli]), layout, annotations)


def test_cli_cache_requires_signature_and_output_hashes(tmp_path):
    export = tmp_path / "map.csv"
    export.write_text("verified export\n")
    (tmp_path / "cli_cache.json").write_text(json.dumps({"signature": "verified", "returncode": 0,
                                                        "outputs_sha256": {export.name: sha256(export)}}))
    assert valid_cli_cache(tmp_path, "verified")
    assert not valid_cli_cache(tmp_path, "changed-parameters")
    export.write_text("changed export\n")
    assert not valid_cli_cache(tmp_path, "verified")
