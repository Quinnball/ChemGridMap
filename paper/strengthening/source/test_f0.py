"""Generate and verify the R1 F0 Prism-ready tables."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from f0_build import (
    DATA_OUT,
    MANIFEST_OUT,
    PLAN_OUT,
    WORKSPACE,
    atr_entry_frame,
    build_all,
    load_sources,
    sha256,
    source_paths,
)


def read(name: str) -> pd.DataFrame:
    return pd.read_csv(DATA_OUT / name, float_precision="round_trip")


def assert_no_absolute_paths() -> None:
    for path in sorted(DATA_OUT.glob("*.csv")):
        text = path.read_text(encoding="utf-8")
        assert "/Users/" not in text and str(WORKSPACE) not in text, path.name


def assert_manifest_sources_and_counts(manifest: dict) -> None:
    paths = source_paths()
    for relative, metadata in manifest["sources"].items():
        path = WORKSPACE / relative
        assert path.is_file(), relative
        assert metadata["sha256"] == sha256(path), relative

    panel_map = {panel["panel"]: panel for panel in manifest["panels"]}
    assert set(panel_map) == {
        "3a", "3b", "3c", "4a", "4b", "4c", "4d", "5a", "5b", "5d", "5e"
    }

    for panel_name, panel in panel_map.items():
        entries_by_endpoint = panel["series"] if isinstance(panel["series"], dict) else None
        endpoints = entries_by_endpoint.keys() if entries_by_endpoint else [None]
        for endpoint in endpoints:
            if endpoint is None:
                filename = panel["files"][0]
                entries = panel["series"]
            else:
                filename = f"{panel_name}_{endpoint}.csv"
                if filename not in panel["files"]:
                    continue
                entries = entries_by_endpoint[endpoint]
            frame = read(filename)
            for entry in entries:
                selected = frame
                if "filter" in entry:
                    item = entry["filter"]
                    selected = selected.loc[selected[item["column"]].eq(item["value"])]
                if entry.get("count_mode", "row") == "cells":
                    actual = int(
                        selected.loc[:, entry["columns"]].notna().to_numpy().sum()
                    )
                else:
                    actual = int(
                        selected.loc[:, entry["columns"]].notna().any(axis=1).sum()
                    )
                assert actual == entry["expected_count"], (
                    panel_name,
                    endpoint,
                    entry["name"],
                    actual,
                    entry["expected_count"],
                )

    for item in manifest["posthoc"]:
        frame = read(item["file"])
        for entry in item["series"]:
            actual = int(
                frame.loc[:, entry["columns"]].notna().any(axis=1).sum()
            )
            assert actual == entry["expected_count"], (item["file"], entry)


def test_f0() -> None:
    missing = [path for path in source_paths().values() if not path.is_file()]
    if missing:
        pytest.skip("Run the E1/E2/E4 analyses first.")

    manifest = build_all()
    data = load_sources(source_paths())

    expected_files = {
        "3a_B0_disagreement.csv",
        "3b_IC50.csv",
        "3b_Ki.csv",
        "3c_IC50.csv",
        "3c_Ki.csv",
        "4a_calibration_IC50.csv",
        "4a_calibration_Ki.csv",
        "4b_hERG_patent_pair.csv",
        "4c_IC50.csv",
        "4c_Ki.csv",
        "4d_IC50.csv",
        "4d_Ki.csv",
        "5a_recall.csv",
        "5b_layout_contrasts.csv",
        "5d_layout_correlation.csv",
        "5e_core_coverage_by_group_size.csv",
        "posthoc_sigma_gain.csv",
        "posthoc_atr_patent_family.csv",
    }
    assert {path.name for path in DATA_OUT.glob("*.csv")} == expected_files
    assert MANIFEST_OUT.is_file() and PLAN_OUT.is_file()
    assert_no_absolute_paths()
    assert_manifest_sources_and_counts(manifest)

    archive = DATA_OUT / "archive_R0"
    assert len(list(archive.glob("*.csv"))) == 15
    assert (archive / "data_manifest.json").is_file()
    assert (PLAN_OUT.parent / "archive_R0" / "gui_plan.md").is_file()

    # 3a remains the endpoint-wide target table from E2.
    metrics = data["e2_metrics_by_target"]
    b0 = metrics.loc[
        metrics.population.eq("all")
        & metrics.model.eq("B0")
        & metrics.stratum.eq("all")
    ].copy()
    b0["disagreement"] = 1 - b0.agreement_rate
    expected_3a = b0.pivot(index="target", columns="endpoint", values="disagreement")
    actual_3a = read("3a_B0_disagreement.csv").set_index("target_id")
    assert set(actual_3a.index) == set(expected_3a.index)
    np.testing.assert_allclose(
        actual_3a.to_numpy(float),
        expected_3a.reindex(index=actual_3a.index, columns=["IC50", "Ki"]).to_numpy(float),
        equal_nan=True,
    )
    assert len(actual_3a) == 52

    # 3b: G comes from the forest table and D/B-range from E4 post hoc.
    forest = data["e2_forest"]
    comparisons = data["e4_comparators"]
    for endpoint, expected_rows in (("IC50", 29), ("Ki", 23)):
        actual = read(f"3b_{endpoint}.csv")
        assert len(actual) == expected_rows
        assert actual.target_rank.tolist() == list(range(1, expected_rows + 1))
        assert actual.g_auc.is_monotonic_decreasing
        source = forest.loc[
            forest.population.eq("all") & forest.endpoint.eq(endpoint)
        ].set_index("target_id")
        np.testing.assert_allclose(actual.g_auc, source.loc[actual.target_id, "g_auc"], atol=1e-15)
        source_e4 = comparisons.loc[
            comparisons.population.eq("all") & comparisons.endpoint.eq(endpoint)
        ]
        np.testing.assert_allclose(
            actual.g_ci_low, source.loc[actual.target_id, "g_auc_ci_low"], atol=1e-15
        )
        np.testing.assert_allclose(
            actual.g_ci_high, source.loc[actual.target_id, "g_auc_ci_high"], atol=1e-15
        )
        for comparator, columns in {
            "D": ("d_auc", "d_auc_ci_low", "d_auc_ci_high"),
            "B-sd": ("b_range_auc", "b_range_auc_ci_low", "b_range_auc_ci_high"),
        }.items():
            expected = source_e4.loc[
                source_e4.comparator.eq(comparator), ["target", "auc", "auc_ci_low", "auc_ci_high"]
            ].set_index("target")
            np.testing.assert_allclose(
                actual.loc[:, columns[0]],
                expected.loc[actual.target_id, "auc"],
                atol=1e-15,
            )
            np.testing.assert_allclose(
                actual.loc[:, columns[1]],
                expected.loc[actual.target_id, "auc_ci_low"],
                atol=1e-15,
            )
            np.testing.assert_allclose(
                actual.loc[:, columns[2]],
                expected.loc[actual.target_id, "auc_ci_high"],
                atol=1e-15,
            )

    # 3c: alternative minus G; X exists only for flagged datasets.
    census = data["e2_flagged_pairs_census"]
    flagged = set(
        census.loc[
            census.row_type.eq("target_summary") & census.n_flagged_pairs.gt(0),
            ["target", "endpoint"],
        ].itertuples(index=False, name=None)
    )
    alternatives = {
        "g_record_minus_g": "G-record",
        "g_global_minus_g": "G-global",
        "g_fixed_minus_g": "G-fixed",
        "x_minus_g": "X",
    }
    assert sum(endpoint == "IC50" for _, endpoint in flagged) == 15
    assert sum(endpoint == "Ki" for _, endpoint in flagged) == 11
    for endpoint in ("IC50", "Ki"):
        actual = read(f"3c_{endpoint}.csv")
        targets = actual.loc[actual.row_type.eq("target")].set_index("target_id")
        loss = metrics.loc[
            metrics.population.eq("all")
            & metrics.endpoint.eq(endpoint)
            & metrics.stratum.eq("all")
            & metrics.model.isin(["G", *alternatives.values()])
        ].pivot(index="target", columns="model", values="log_loss")
        for output_column, model in alternatives.items():
            expected = loss[model] - loss["G"]
            if model == "X":
                expected = expected.loc[
                    [target for target, item_endpoint in flagged if item_endpoint == endpoint]
                ]
                assert targets.loc[expected.index, output_column].notna().all()
                assert int(targets[output_column].notna().sum()) == len(expected)
                np.testing.assert_allclose(
                    targets.loc[expected.index, output_column], expected, atol=1e-15
                )
            else:
                np.testing.assert_allclose(
                    targets.loc[expected.index, output_column], expected, atol=1e-15
                )
        summary = actual.loc[actual.row_type.eq("summary")].iloc[0]
        for output_column in alternatives:
            values = targets[output_column].dropna()
            assert np.isclose(summary[f"median_{output_column}"], values.median(), atol=1e-15)
            assert summary[f"ci_low_{output_column}"] <= summary[f"median_{output_column}"]
            assert summary[f"median_{output_column}"] <= summary[f"ci_high_{output_column}"]

    # 4a: eligible E4 bins are plotted directly; B0 is one point at probability 1.
    calibration = data["e4_calibration"]
    for endpoint in ("IC50", "Ki"):
        actual = read(f"4a_calibration_{endpoint}.csv")
        source = calibration.loc[
            calibration.endpoint.eq(endpoint)
            & calibration.eligible.astype(bool)
        ]
        g = source.loc[source.model.eq("G")].sort_values("bin")
        b0_source = source.loc[source.model.eq("B0")]
        assert len(b0_source) == 1
        np.testing.assert_allclose(actual.G_probability, g.mean_probability, atol=1e-15)
        np.testing.assert_allclose(
            actual.G_observed_agreement, g.observed_agreement, atol=1e-15
        )
        assert int(actual.G_n_predictions.notna().sum()) == len(g)
        assert actual.loc[actual.B0_probability.notna(), "B0_probability"].iloc[0] == 1.0
        assert np.isclose(
            actual.loc[
                actual.B0_observed_agreement.notna(), "B0_observed_agreement"
            ].iloc[0],
            b0_source.observed_agreement.iloc[0],
            atol=1e-15,
        )

    # 4b: unchanged 30-compound patent pair and median offset.
    retained = data["cgm_hERG_retained"]
    docs = {"CHEMBL5725903", "CHEMBL5727624"}
    units = (
        retained.loc[retained.document_chembl_id.isin(docs)]
        .groupby(["molecule_identity_key", "document_chembl_id"], as_index=False)
        .pchembl_value.median()
    )
    first = units.loc[units.document_chembl_id.eq("CHEMBL5725903")].set_index("molecule_identity_key")
    second = units.loc[units.document_chembl_id.eq("CHEMBL5727624")].set_index("molecule_identity_key")
    shared = first.index.intersection(second.index)
    table_4b = read("4b_hERG_patent_pair.csv")
    assert len(shared) == len(table_4b) == 30
    assert set(shared) == set(table_4b.compound_id)
    assert abs(float(table_4b.offset_b_minus_a.median()) - 3.47) <= 0.005

    # 4c/4d: latest-unit-year split at 2015, endpoints only.
    time_split = data["e4_time_split"]
    expected_4c = time_split.loc[
        time_split.record_type.eq("dataset_comparison")
        & time_split.variant.eq("latest_unit_year")
        & time_split.cutoff_year.eq(2015)
        & time_split.endpoint.isin(["IC50", "Ki"])
    ]
    table_4c = pd.concat(
        [read("4c_IC50.csv"), read("4c_Ki.csv")], ignore_index=True
    )
    assert len(expected_4c) == len(table_4c) == 33
    assert table_4c.supported_agreement.gt(table_4c.uncertain_agreement).sum() == 31
    for endpoint in ("IC50", "Ki"):
        actual = table_4c.loc[table_4c.index.isin([])]
        actual = read(f"4c_{endpoint}.csv").set_index("target_id")
        expected = expected_4c.loc[expected_4c.endpoint.eq(endpoint)].set_index("target")
        np.testing.assert_allclose(
            actual.supported_agreement,
            expected.loc[actual.index, "supported_agreement_rate"],
            atol=1e-15,
        )
        np.testing.assert_allclose(
            actual.uncertain_agreement,
            expected.loc[actual.index, "uncertain_agreement_rate"],
            atol=1e-15,
        )
    expected_4d = time_split.loc[
        time_split.record_type.eq("m_past_stratum")
        & time_split.variant.eq("latest_unit_year")
        & time_split.cutoff_year.eq(2015)
        & time_split.endpoint.isin(["IC50", "Ki"])
    ]
    table_4d = pd.concat(
        [read("4d_IC50.csv"), read("4d_Ki.csv")], ignore_index=True
    )
    assert len(expected_4d) == 8 and len(table_4d) == 4
    assert {"m=1", "m>=2"} == set(table_4d.m_past)

    # 5a remains ten seeds for each of three layouts and targets.
    recall_source = data["cgm_figure_3a"]
    recall = read("5a_recall.csv")
    assert recall.shape == (3, 31)
    for row in recall.itertuples(index=False):
        label = row.Row_title
        for layout, offset in (
            ("projection", 1),
            ("minimum_displacement", 11),
            ("neighbor_refined", 21),
        ):
            actual_mean = np.mean(row[offset:offset + 10])
            expected_mean = recall_source.loc[
                recall_source.target_label.eq(label) & recall_source.layout.eq(layout),
                "mean_morgan_neighbor_recall_k10",
            ].mean()
            assert np.isclose(actual_mean, expected_mean, atol=1e-12)

    # 5b: other layout minus refined, including the reversed minimum-displacement row.
    baseline = data["e4_layout_baselines"]
    expected_contrasts = baseline.loc[
        baseline.row_type.eq("paired_t_interval")
        & baseline.metric.eq("mean_e_grid")
    ].copy()
    table_5b = read("5b_layout_contrasts.csv")
    assert len(table_5b) == 12
    for row in table_5b.itertuples(index=False):
        source = expected_contrasts.loc[
            expected_contrasts.target.eq(row.target_id)
            & expected_contrasts.comparison.eq(row.comparison)
        ].iloc[0]
        expected = float(source.paired_difference_mean)
        if row.comparison == "neighbor_refined_minus_minimum_displacement":
            expected = -expected
        assert np.isclose(row.other_minus_refined_e_grid, expected, atol=1e-15)
    herg_las = table_5b.loc[
        table_5b.target_label.eq("hERG")
        & table_5b.comparison.eq("las_random_minus_neighbor_refined"),
        "other_minus_refined_e_grid",
    ].iloc[0]
    assert round(float(herg_las), 3) == -0.024
    herg_minimum = table_5b.loc[
        table_5b.target_label.eq("hERG")
        & table_5b.comparison.eq("neighbor_refined_minus_minimum_displacement"),
        "other_minus_refined_e_grid",
    ].iloc[0]
    assert round(float(herg_minimum), 3) == 0.005

    # 5d: 77 included sampled datasets and three seed-42 reference targets.
    layout_source = data["e2c_layout_by_target"]
    sampled_source = layout_source.loc[
        layout_source.included.astype(bool)
        & layout_source.layout.eq("neighbor_refined")
    ]
    table_5d = read("5d_layout_correlation.csv")
    sampled = table_5d.loc[table_5d.row_type.eq("sampled_dataset")]
    references = table_5d.loc[table_5d.row_type.eq("detailed_reference")]
    assert len(sampled) == len(sampled_source) == 77
    assert sampled.endpoint.value_counts().to_dict() == {"IC50": 48, "Ki": 29}
    assert len(references) == 3
    expected_sampled = sampled_source.assign(target_id=sampled_source.target.str.upper())
    actual_sampled = sampled.set_index(["endpoint", "target_id"]).sort_index()
    expected_sampled = expected_sampled.set_index(["endpoint", "target_id"]).sort_index()
    np.testing.assert_allclose(
        actual_sampled.rho, expected_sampled.rho.loc[actual_sampled.index], atol=1e-15
    )
    rs = data["e1b_neighbor_sar_summary"].loc[
        lambda frame: frame.layout.eq("neighbor_refined")
    ]
    assert set(references.target_label) == set(rs.target_label)
    np.testing.assert_allclose(
        references.rho.sort_values().to_numpy(),
        rs.rho_seed42.sort_values().to_numpy(),
        atol=1e-15,
    )

    # 5e: old core-coverage panel, with group counts and mean core size.
    core = data["e1b_core_by_group_size"]
    core_overall = core.loc[core.target.eq("overall")]
    table_5e = read("5e_core_coverage_by_group_size.csv")
    assert len(table_5e) == 12
    for rule in ("most_frequent", "coverage", "fmcs"):
        selected = table_5e.loc[table_5e.rule.eq(rule)]
        assert selected.groups.sum() == 480
        recomputed = (
            selected.coverage_ge_0_8
            + selected.coverage_gt_0_5_lt_0_8
            + selected.coverage_le_0_5
            + selected.no_core
        )
        assert recomputed.equals(selected.groups)
        for stratum, row in selected.set_index("size_stratum").iterrows():
            source = core_overall.loc[
                core_overall.rule.eq(rule) & core_overall.size_stratum.eq(stratum)
            ]
            assert int(row.coverage_ge_0_8) == int(source.target_met_groups.sum())
            assert int(row.no_core) == int(source.no_core_groups.sum())
            assert np.isclose(
                row.mean_core_heavy_atoms, source.mean_core_heavy_atoms.mean()
            )

    # Existing post hoc outputs remain unchanged and independently checkable.
    sigma = data["e2_sigma_doc"]
    sigma_targets = sigma.loc[
        sigma.row_type.eq("target") & sigma.measurement_type.eq("all"),
        ["target", "endpoint", "sigma_doc"],
    ]
    all_forest = forest.loc[forest.population.eq("all")][[
        "target_id", "endpoint", "g_minus_g_fixed_log_loss"
    ]].merge(
        sigma_targets,
        left_on=["target_id", "endpoint"],
        right_on=["target", "endpoint"],
    )
    all_forest["distance"] = (all_forest.sigma_doc - 0.64).abs()
    all_forest["gain"] = -all_forest.g_minus_g_fixed_log_loss
    for endpoint, expected in (("IC50", 0.57), ("Ki", 0.60)):
        rho = all_forest.loc[
            all_forest.endpoint.eq(endpoint), ["distance", "gain"]
        ].corr(method="spearman").iloc[0, 1]
        assert round(float(rho), 2) == expected

    entries, pair, shared = atr_entry_frame(data)
    multi = entries.loc[entries.n_distinct_values.ge(2)]
    assert len(entries) == 1011
    assert len(multi) == 654
    assert round(float(multi.range_le_0_1.mean()), 3) == 0.732
    atr = read("posthoc_atr_patent_family.csv")
    summary = atr.loc[atr.row_type.eq("summary")].iloc[0]
    assert int(summary.n_multi_value_entries) == 654
    assert np.isclose(summary.fraction_multi_value_range_le_0_1, multi.range_le_0_1.mean())
    doc_pair = atr.loc[atr.row_type.eq("document_pair")].iloc[0]
    assert {doc_pair.document_a, doc_pair.document_b} == set(pair)
    assert int(doc_pair.n_shared_entries) == int(shared)

    print("F0 R1 PASS")
    print("tables=18 panels=11 reused=1 posthoc=2")
    print("3a=52; 3b=29/23; 3c=29/23 with X=15/11")
    print("4a=G+1 eligible bins per endpoint; 4b=30; 4c=33 (31 supported>uncertain); 4d=4 rows (2 per endpoint)")
    print("5a=3x10; 5b=12; 5d=77+3; 5e=12")
    print(f"5b hERG LAS-random other-minus-refined={herg_las:.15g}")
    print(f"hERG patent median offset={table_4b.offset_b_minus_a.median():.15g}")
    print("PASS: hashes, row/series counts, values, and R1 targets verified")


if __name__ == "__main__":
    test_f0()
