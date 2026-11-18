"""Independent checks for the paper analysis definitions and case exports."""
import json

import numpy as np
import pandas as pd
import pytest
from scipy.spatial.distance import cdist
from sklearn.manifold import trustworthiness

from reliability_analysis import compare_frozen, jaccard_distances, mean_jaccard, evaluate_neighborhoods
from reliability_cases import neighbor_case_table, compare_v2_neighbors


def test_binary_jaccard_matches_scipy_including_empty_fingerprints():
    bits = np.random.default_rng(7).integers(0, 2, (24, 20))
    bits[:2] = 0
    np.testing.assert_allclose(jaccard_distances(bits), cdist(bits, bits, metric="jaccard"), atol=1e-15)


def test_precomputed_trustworthiness_matches_original_metric():
    bits = np.random.default_rng(21).integers(0, 2, (35, 64))
    xy = np.random.default_rng(42).normal(size=(35, 2))
    actual, _, _ = evaluate_neighborhoods(bits, xy)
    expected = trustworthiness(bits.astype(bool), xy, n_neighbors=10, metric="jaccard")
    assert actual["trustworthiness_k10"] == pytest.approx(expected, abs=1e-14)


def test_seed_jaccard_averages_per_entry_not_pooled_counts():
    assert mean_jaccard([{1}, {0, 2, 3}], [{1}, {0, 4, 5}]) == pytest.approx(.6)


def frozen_fixture():
    return pd.DataFrame({"molecule_identity_key": ["A", "B"], "grid_col": [0, 1],
        "grid_row": [1, 0], "projection_x": [0., 1.], "projection_y": [1., 0.],
        "projection_x_raw": [0., 5.], "projection_y_raw": [5., 0.]})


def test_reproduction_detects_raw_coordinate_and_cell_changes():
    original = frozen_fixture()
    assert compare_frozen(original, original)["passed"]
    changed = original.copy()
    changed.loc[0, "projection_x_raw"] += 2e-5
    assert not compare_frozen(changed, original)["passed"]
    changed = original.copy()
    changed.loc[0, "grid_col"] += 1
    assert compare_frozen(changed, original)["grid_mismatched_entries"] == 1
    with pytest.raises(ValueError, match="identity order"):
        compare_frozen(original.iloc[::-1], original)


def case_fixture():
    data = pd.DataFrame({"molecule_identity_key": ["Q", "S", "F", "G"]})
    fp = np.array([[1, 0, 0, 0], [1, 1, 0, 0], [1, 0, 1, 0], [0, 0, 0, 1]])
    cells = np.array([[0, 0], [1, 0], [0, 3], [0, 1]])
    return neighbor_case_table(data, fp, cells, query_id="Q", layout="minimum_displacement", k=1)


def test_case_markers_include_ties_and_exclude_query():
    table, summary = case_fixture()
    assert dict(zip(table.molecule_identity_key, table.figure_marker)) == {"Q": "Q", "S": "S", "F": "F", "G": "G"}
    assert summary["fingerprint_neighbor_count"] == summary["grid_neighbor_count"] == 2
    assert summary["shared_neighbor_count"] == 1
    assert summary["neighbor_recall"] == summary["neighbor_precision"] == .5
    assert not table.loc[table.is_query, ["representation_neighbor", "grid_neighbor"]].any().any()


def test_case_v2_comparison_detects_coordinate_change():
    table, _ = case_fixture()
    old = table.loc[~table.is_query].rename(columns={"representation_distance_tanimoto": "representation_distance"}).copy()
    assert compare_v2_neighbors(table, old)["passed"]
    old.loc[old.index[0], "grid_col"] += 1
    assert not compare_v2_neighbors(table, old)["passed"]


@pytest.fixture
def baseline_checkpoints(tmp_path, monkeypatch):
    import reliability_baselines as baselines

    data = pd.DataFrame({"molecule_identity_key": ["A", "B", "C", "D"]})
    features = np.array([[0., 0.], [1., 0.], [0., 2.], [2., 3.]])
    high = cdist(features, features)
    cells = np.array([[0, 0], [1, 0], [0, 1], [1, 1]])
    layout = data.assign(projection_x=[0., .3, .1, 1.], projection_y=[0., .1, .7, 1.])
    layout_path = tmp_path / "layouts/example/coordinates_seed_42.csv"
    layout_path.parent.mkdir(parents=True)
    layout.to_csv(layout_path, index=False)
    out = tmp_path / "baselines"
    (out / "checkpoints").mkdir(parents=True)
    rows = []
    for occupancy in baselines.OCCUPANCIES:
        group = []
        for index, method in enumerate(baselines.METHODS):
            points = np.roll(cells, index, axis=0)
            scores = baselines.dpq_scores(high, points)
            group.append({"target": "example", "target_occupancy": occupancy,
                          "seed": 42, "method": method, "grid_side": 4 if occupancy == .4 else 2,
                          "dpq16_tanimoto": scores[16], "dpq2_tanimoto": scores[2],
                          "total_seconds": .12345678901234567 + index})
        csv_path, npz_path, metadata_path = baselines._checkpoint_paths(tmp_path, "example", occupancy, 42)
        pd.DataFrame(group).to_csv(csv_path, index=False)
        np.savez_compressed(npz_path, molecule_identity_key=data.molecule_identity_key.to_numpy(str),
                            **{method: np.roll(cells, i, axis=0) for i, method in enumerate(baselines.METHODS)})
        baselines.dump_json(metadata_path, {"cache_key": "old-source-hash",
                            "csv_sha256": baselines.sha256(csv_path),
                            "coordinates_sha256": baselines.sha256(npz_path)})
        rows.extend(group)
    frame = pd.DataFrame(rows)
    frame.to_csv(out / "baseline_metrics.csv", index=False)
    frame.groupby(["target", "target_occupancy", "method"])[["dpq16_tanimoto", "dpq2_tanimoto"]].agg(
        ["mean", "min", "max"]).to_csv(out / "summary.csv")
    (out / "runtime_seed42.csv").write_text("original timing bytes\n")
    provenance = {"code": {"source_sha256": {"old.py": "old-hash"}}, "parameters": {"phase": "all"}}
    baselines.dump_json(out / "manifest.json", {
        "inputs_sha256": {str(layout_path): baselines.sha256(layout_path)}, **provenance})
    monkeypatch.setattr(baselines, "_prepared", lambda target: (data, None, None, high, high, None))

    def forbid_layout_work(*args, **kwargs):
        raise AssertionError("metrics-only must not construct layouts or refresh checkpoints")

    for name in ("project_2d", "assign_to_grid", "refine_grid_layout", "las", "flas",
                 "random_grid_layout", "_measure_method", "_cache_key", "_save_group"):
        monkeypatch.setattr(baselines, name, forbid_layout_work)
    monkeypatch.setattr(baselines, "write_manifest", lambda directory, inputs, parameters: {
        "code": {"source_sha256": {"new.py": "new-hash"}}, "parameters": parameters})
    return tmp_path, baselines, provenance


def test_baseline_metrics_only_preserves_existing_results_and_checkpoints(baseline_checkpoints):
    root, baselines, provenance = baseline_checkpoints
    out = root / "baselines"
    old = pd.read_csv(out / "baseline_metrics.csv", float_precision="round_trip")
    old_summary = pd.read_csv(out / "summary.csv", header=[0, 1], index_col=[0, 1, 2],
                              float_precision="round_trip")
    protected = {path: path.read_bytes() for path in (out / "checkpoints").iterdir()}
    protected[out / "runtime_seed42.csv"] = (out / "runtime_seed42.csv").read_bytes()
    result = baselines.run_baselines(root, seeds=(42,), targets=("example",), phase="metrics-only")
    written = pd.read_csv(out / "baseline_metrics.csv", float_precision="round_trip")
    pd.testing.assert_frame_equal(result[old.columns], old, check_exact=True)
    pd.testing.assert_frame_equal(written[old.columns], old, check_exact=True)
    assert list(written.columns) == [*old.columns, *baselines.ADDED_METRICS]
    assert np.isfinite(written[list(baselines.ADDED_METRICS)]).all().all()
    summary = pd.read_csv(out / "summary.csv", header=[0, 1], index_col=[0, 1, 2],
                         float_precision="round_trip")
    pd.testing.assert_frame_equal(summary[old_summary.columns], old_summary, check_exact=True)
    assert all(path.read_bytes() == contents for path, contents in protected.items())
    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["layout_generation_provenance"] == provenance
    audit = json.loads((out / "metrics_only_verification.json").read_text())
    assert audit["passed"] and audit["rows"] == 12
    assert audit["checkpoint_dpq16_max_absolute_difference"] <= 1e-12
    assert not audit["layouts_recomputed"] and not audit["timings_recomputed"]


def test_baseline_metrics_only_rejects_wrong_dpq_before_writing(baseline_checkpoints):
    root, baselines, _ = baseline_checkpoints
    out = root / "baselines"
    csv_path, _, metadata_path = baselines._checkpoint_paths(root, "example", .4, 42)
    checkpoint = pd.read_csv(csv_path, float_precision="round_trip")
    checkpoint.loc[0, "dpq16_tanimoto"] += .01
    checkpoint.to_csv(csv_path, index=False)
    metadata = json.loads(metadata_path.read_text())
    metadata["csv_sha256"] = baselines.sha256(csv_path)
    baselines.dump_json(metadata_path, metadata)
    metrics = pd.read_csv(out / "baseline_metrics.csv", float_precision="round_trip")
    metrics.loc[0, "dpq16_tanimoto"] = checkpoint.loc[0, "dpq16_tanimoto"]
    metrics.to_csv(out / "baseline_metrics.csv", index=False)
    before = {path: path.read_bytes() for path in out.rglob("*") if path.is_file()}
    with pytest.raises(ValueError, match="Checkpoint DPQ mismatch"):
        baselines.run_baselines(root, seeds=(42,), targets=("example",), phase="metrics-only")
    after = {path: path.read_bytes() for path in out.rglob("*") if path.is_file()}
    assert after == before


# Keep the task's fixed test command inclusive of the T2b helper checks.
from test_reliability_activity import (
    test_ridge_recovers_known_offsets_with_unpenalized_entry_means,
    test_internal_cv_selects_offset_recovery_and_accounts_for_every_unit,
    test_merged_offset_uses_representative_document_and_unseen_is_zero,
    test_exclusion_checks_all_constituents_and_empty_rest_falls_back,
    test_calibration_uses_displayed_probability_instead_of_maximum,
    test_adoption_has_no_tolerance_and_never_selects_n,
    test_support_summary_keeps_empty_strata_missing_and_threshold_boundaries,
)
from test_reliability_cores import (
    test_missing_core_is_not_zero_coverage_and_majority_is_strict,
    test_fmcs_smiles_uses_actual_charge_from_member_and_keeps_query_authoritative,
    test_example_ranking_prefers_larger_cores_then_coverage_gain,
    test_fmcs_audit_detects_coverage_members_and_query_changes,
    test_timed_most_frequent_rule_matches_software_ties_and_missing_core,
    test_core_timer_summary_reconciles_current_history_and_timeout_denominators,
    test_full_fmcs_summary_separates_partial_timeout_and_completed_differences,
    test_timing_runs_every_group_and_preserves_old_outputs_with_timeout_difference,
)
from test_reliability_contrasts import (
    test_eight_connected_occupied_pairs_match_all_pair_distances,
    test_tie_aware_union_includes_shell_ties_and_one_way_edges,
    test_contrast_counts_and_expectations_use_consensus_probabilities,
    test_figure1_comparison_aligns_identity_and_detects_coordinates_and_support,
    test_figure1_comparison_rejects_missing_or_duplicate_entries,
    test_cli_cache_requires_signature_and_output_hashes,
)


def test_t2b_generated_section_preserves_frozen_prefix_and_is_rerunnable(tmp_path):
    from reliability_t2b import append_section
    path = tmp_path / "RESULTS.md"
    original = "# Accepted T2a\n\noriginal table 1.234\n\n"
    path.write_bytes(original.encode("utf-8"))
    append_section(path, "<!-- T2b -->", "new results")
    first = path.read_bytes()
    assert first.startswith(original.encode())
    append_section(path, "<!-- T2b -->", "new results")
    assert path.read_bytes() == first
    append_section(path, "<!-- T2b -->", "revised results")
    assert path.read_text(encoding="utf-8") == original + "<!-- T2b -->\n\nrevised results\n"
