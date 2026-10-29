"""E1-B：并列邻居、手算指标、重抽样和冻结数据的一致性。"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import spearmanr

from e1b_cores import size_stratum, stratify_cores
from e1b_neighbors import (
    neighborhood_metrics, bootstrap_spearman, paired_mean_ci,
    representation_neighbor_sets, grid_neighbor_sets,
    recall_bins, check_frozen_layout_metrics, load_target, layout_points,
)

ROOT = Path(__file__).resolve().parents[3]
OUTPUT = Path(__file__).resolve().parents[1] / "output/E1B"
FROZEN = ROOT / "paper/output/reliability_v1"


def test_fingerprint_neighbors_keep_kth_ties_and_exclude_self():
    bits = np.array([[1, 1, 0, 0], [1, 0, 1, 0],
                     [0, 1, 0, 1], [0, 0, 1, 1]])
    neighbors = representation_neighbor_sets(bits, k=1)
    assert set(neighbors[0]) == {1, 2}
    assert all(i not in group for i, group in enumerate(neighbors))


def test_grid_neighbors_keep_ties_and_distinct_duplicate_positions():
    points = np.array([[0, 0], [1, 0], [-1, 0], [0, 2]])
    assert set(grid_neighbor_sets(points, k=1)[0]) == {1, 2}
    duplicates = grid_neighbor_sets(np.array([[0, 0], [0, 0], [1, 0], [2, 0]]), k=1)
    assert set(duplicates[0]) == {1}
    assert set(duplicates[1]) == {0}
    assert all(i not in group for i, group in enumerate(duplicates))


def test_recall_and_activity_error_hand_calculation():
    reference = [np.array(x) for x in ([1, 2], [0], [1, 3], [2])]
    displayed = [np.array(x) for x in ([1, 3], [0, 2], [0, 1], [1])]
    result = neighborhood_metrics(np.array([0., 1., 4., 6.]), reference, displayed)
    np.testing.assert_allclose(result["recall"], [.5, 1., .5, 0.])
    np.testing.assert_allclose(result["e_fp"], [2.5, 1., 2.5, 2.])
    np.testing.assert_allclose(result["e_grid"], [3.5, 2., 3.5, 5.])
    np.testing.assert_allclose(result["excess_error"], [1., 1., 1., 3.])
    assert np.mean(result["e_grid"]) == 3.5  # 分子等权，而非把所有边混在一起平均。


@pytest.mark.parametrize("size, expected", [(2, "2"), (3, "3-4"), (4, "3-4"),
    (5, "5-9"), (9, "5-9"), (10, ">=10"), (42, ">=10")])
def test_prespecified_core_size_boundaries(size, expected):
    assert size_stratum(size) == expected


def test_no_core_is_missing_and_half_is_not_strict_majority():
    table = pd.DataFrame({"target": ["chembl205"] * 3, "size": [2] * 3,
        "rule": ["most_frequent"] * 3, "has_core": [True, False, True],
        "strict_majority": [False, False, True], "target_met": [False, False, True],
        "coverage_fraction": [.5, np.nan, 1.], "core_heavy_atoms": [6, 0, 10],
        "is_benzene": [True, False, False]})
    result = stratify_cores(table)
    row = result.loc[result.target.eq("overall") & result.size_stratum.eq("2")
                     & result.rule.eq("most_frequent")].iloc[0]
    assert (row.groups, row.no_core_groups, row.strict_majority_groups, row.target_met_groups) == (3, 1, 1, 1)
    assert row.mean_coverage == .75
    assert row.mean_core_heavy_atoms == 8
    assert row.benzene_cores == 1


def test_recall_bin_boundaries_and_empty_bin_means():
    values = [0., .199999, .2, .4, .6, .8, 1.]
    metrics = pd.DataFrame({"recall": values, "e_grid": np.arange(7.),
                            "e_fp": np.zeros(7), "excess_error": np.arange(7.)})
    bins = recall_bins(metrics)
    assert [row["n"] for row in bins] == [2, 1, 1, 1, 2]
    assert bins[-1]["mean_e_grid"] == 5.5
    empty = recall_bins(metrics.iloc[[0]])[1]
    assert empty["n"] == 0 and np.isnan(empty["mean_e_grid"]) and np.isnan(empty["mean_e_fp"])


def test_spearman_bootstrap_reranks_each_resample_and_reproduces():
    x = np.array([0., 0., .2, .4, .4, .6, .8, 1.])
    y = np.array([1., 4., 3., 4., 2., -1., 0., -2.])
    result = bootstrap_spearman(x, y, n_boot=64, seed=20261031)
    assert result == bootstrap_spearman(x, y, n_boot=64, seed=20261031)
    rng = np.random.default_rng(20261031)
    expected = []
    for _ in range(64):
        index = rng.integers(0, len(x), len(x))
        expected.append(spearmanr(x[index], y[index]).statistic)
    finite = np.array(expected)[np.isfinite(expected)]
    np.testing.assert_allclose([result["rho_ci_low"], result["rho_ci_high"]],
                               np.quantile(finite, [.025, .975]), atol=1e-14, rtol=0)
    assert result["bootstrap_valid_count"] == len(finite)


def test_paired_mean_ci_resamples_seed_differences():
    differences = np.array([-.1, -.2, .01, -.04, -.3, .03, -.07, -.12, -.09, -.01])
    result = paired_mean_ci(differences, n_boot=1000, seed=20261031)
    rng = np.random.default_rng(20261031)
    means = differences[rng.integers(0, 10, (1000, 10))].mean(axis=1)
    assert result["paired_difference_mean"] == pytest.approx(differences.mean())
    np.testing.assert_allclose([result["paired_ci_low"], result["paired_ci_high"]],
                               np.quantile(means, [.025, .975]), atol=1e-14, rtol=0)
    assert result["n_seeds"] == 10


@pytest.mark.parametrize("target", ["chembl205", "chembl204", "chembl240"])
def test_recomputed_seed42_activity_error_matches_frozen(target):
    if not (ROOT / "paper/output/revision_v6" / target / "map_molecules.csv").is_file():
        pytest.skip("Restore the frozen revision_v6 paper data first.")
    data, fingerprints = load_target(target)
    reference = representation_neighbor_sets(fingerprints, k=10)
    coordinates = pd.read_csv(FROZEN / "layouts" / target / "coordinates_seed_42.csv")
    frozen = pd.read_csv(FROZEN / "neighborhoods/layout_metrics.csv")
    frozen = frozen.loc[frozen.target.eq(target) & frozen.seed.eq(42)].set_index("layout")
    for layout in ("projection", "minimum_displacement", "neighbor_refined"):
        neighbors = grid_neighbor_sets(layout_points(coordinates, layout), k=10)
        metrics = neighborhood_metrics(data.activity_pchembl, reference, neighbors)
        assert metrics.e_grid.mean() == pytest.approx(frozen.loc[layout, "neighbor_abs_dpchembl"], abs=1e-12, rel=0)
        assert metrics.e_fp.mean() == pytest.approx(frozen.loc[layout, "fingerprint_neighbor_abs_dpchembl"], abs=1e-12, rel=0)


def test_all_90_exported_metrics_match_frozen_and_detect_corruption():
    if not (OUTPUT / "layout_activity_error.csv").is_file():
        pytest.skip("Run paper/strengthening/source/e1b_neighbors.py first.")
    exported = pd.read_csv(OUTPUT / "layout_activity_error.csv")
    actual = exported.loc[exported.row_type.eq("seed_layout")]
    frozen = pd.read_csv(FROZEN / "neighborhoods/layout_metrics.csv")
    check = check_frozen_layout_metrics(actual, frozen)
    assert check["passed"] and check["rows"] == 90
    corrupted = actual.copy()
    corrupted.loc[corrupted.index[0], "mean_e_grid"] += .01
    with pytest.raises(ValueError, match="neighbor_abs_dpchembl"):
        check_frozen_layout_metrics(corrupted, frozen)


def test_output_population_and_bootstrap_completeness():
    required = [
        OUTPUT / "neighbor_sar_summary.csv",
        OUTPUT / "neighbor_sar_bins.csv",
        OUTPUT / "layout_activity_error.csv",
        OUTPUT / "fmcs_budget_by_group.csv",
    ]
    if not all(path.is_file() for path in required):
        pytest.skip("Run the E1-B analysis scripts first.")
    summary = pd.read_csv(OUTPUT / "neighbor_sar_summary.csv")
    assert len(summary) == 27
    assert summary.bootstrap_valid_count.eq(1000).all()
    assert summary.loc[summary.layout_family.eq("main"), "n_seeds"].eq(10).all()
    assert summary.loc[summary.layout_family.eq("baseline"), "n_seeds"].eq(1).all()
    bins = pd.read_csv(OUTPUT / "neighbor_sar_bins.csv")
    assert len(bins) == 540
    counts = bins.groupby(["target", "seed", "layout"]).n.sum()
    errors = pd.read_csv(OUTPUT / "layout_activity_error.csv")
    expected_rows = errors.loc[errors.row_type.eq("seed_layout")].copy()
    # 汇总行的空 seed 使 CSV 读回为 float；核对整数值后统一标识类型。
    assert expected_rows.seed.eq(expected_rows.seed.round()).all()
    expected_rows["seed"] = expected_rows.seed.astype(int)
    expected = expected_rows.set_index(["target", "seed", "layout"]).n
    pd.testing.assert_series_equal(counts.sort_index().astype(int), expected.sort_index().astype(int), check_names=False)
    raw = pd.read_csv(OUTPUT / "fmcs_budget_by_group.csv")
    assert len(raw) == 1920
    assert raw.groupby("timeout_seconds").size().to_dict() == {1: 480, 2: 480, 5: 480, 10: 480}
    assert not raw.duplicated(["target", "group", "timeout_seconds"]).any()
    assert np.isfinite(raw.seconds).all()
