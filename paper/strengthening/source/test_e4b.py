"""Tests for E4b temporal-split sensitivities."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from e2b_analysis import API
from e4b_analysis import (
    BASELINE,
    INCLUSION,
    RECORDS,
    add_unit_years,
    build_posthoc,
    build_temporal_frame,
    dataset_comparisons,
    load_records,
)

ORIGINAL_MANIFEST = API / "download_original_ic50_manifest.json"


def require_e4b_data():
    if not all(path.is_file() for path in (BASELINE, RECORDS, INCLUSION)):
        pytest.skip("Run the E2 analyses first.")
    if not ORIGINAL_MANIFEST.is_file():
        pytest.skip("Run the E2 analyses and download the selected ChEMBL 37 records first.")
    manifest = json.loads(ORIGINAL_MANIFEST.read_text())
    if not all((API / "activities" / item["file"]).is_file() for item in manifest["files"]):
        pytest.skip("Download the selected ChEMBL 37 records first.")


def test_latest_year_rebuild_matches_accepted_temporal_predictions():
    require_e4b_data()
    inclusion = pd.read_csv(
        INCLUSION,
        float_precision="round_trip",
    )
    rebuilt = build_temporal_frame(load_records(), inclusion, "latest")
    accepted = pd.read_csv(BASELINE, float_precision="round_trip")
    keys = ["target", "endpoint", "cutoff_year", "entry", "future_unit"]
    left = accepted.sort_values(keys).reset_index(drop=True)
    right = rebuilt.sort_values(keys).reset_index(drop=True)
    pd.testing.assert_frame_equal(
        left[keys],
        right[keys],
        check_dtype=False,
    )
    for column in (
        "m_past",
        "sigma_past",
        "past_center",
        "past_support",
        "displayed",
        "observed",
        "p0",
        "p1",
        "p2",
        "p_displayed",
    ):
        np.testing.assert_allclose(
            left[column].astype(float),
            right[column].astype(float),
            rtol=0,
            atol=1e-14,
        )
    assert left.display_status.equals(right.display_status)
    assert left.agree.equals(right.agree)


def test_earliest_year_uses_minimum_component_document_year():
    documents = pd.DataFrame({
        "entry": ["E1"],
        "document": ["D1"],
        "document_value": [6.0],
        "n_records": [2],
        "documents": ["D1;D2"],
    })
    records = pd.DataFrame({
        "document_chembl_id": ["D1", "D2"],
        "document_year": [2020, 2001],
    })
    earliest = add_unit_years(documents, records, "earliest")
    latest = add_unit_years(documents, records, "latest")
    assert earliest.year.tolist() == [2001.0]
    assert latest.year.tolist() == [2020.0]


def test_dataset_paired_summary_requires_five_predictions_per_group():
    rows = []
    for target in ("A", "B", "C"):
        for status, agree in (
            (["supported"] * 5, [True, True, True, True, False]),
            (["uncertain"] * 5, [True, False, False, False, False]),
        ):
            for item_status, item_agree in zip(status, agree):
                rows.append({
                    "target": target,
                    "endpoint": "IC50",
                    "cutoff_year": 2015,
                    "display_status": item_status,
                    "agree": item_agree,
                })
    rows.append({
        "target": "D",
        "endpoint": "IC50",
        "cutoff_year": 2015,
        "display_status": "supported",
        "agree": True,
    })
    frame = pd.DataFrame(rows)
    output = dataset_comparisons(frame, "test")
    datasets = [
        row
        for row in output
        if row["record_type"] == "dataset_comparison"
        and row["endpoint"] == "IC50"
    ]
    summaries = [
        row
        for row in output
        if row["record_type"] == "paired_summary"
        and row["endpoint"] == "IC50"
    ]
    assert len(datasets) == 3
    assert summaries[0]["n_datasets"] == 3
    assert summaries[0]["positive_difference_datasets"] == 3
    assert summaries[0]["wilcoxon_p_one_sided"] < 0.5


def test_posthoc_output_has_required_sections_and_post_hoc_role():
    require_e4b_data()
    output = build_posthoc()
    assert set(output.analysis_section) == {1, 2, 3, 4}
    assert set(output.variant) == {
        "latest_unit_year",
        "earliest_unit_year",
        "drop_near_duplicate_future_units",
    }
    assert (output.analysis_role == "post hoc").all()
    summary = output.loc[
        (output.record_type == "pooled_summary")
        & (output.sample_filter == "exclude_CHEMBL5113_Ki")
        & (output.endpoint == "ALL")
    ]
    assert len(summary) == 18
    assert summary.n_predictions.gt(0).all()
    near_duplicate = output.loc[
        output.variant == "drop_near_duplicate_future_units"
    ]
    latest = output.loc[output.variant == "latest_unit_year"]
    n_latest = latest.loc[
        (latest.analysis_section == 1)
        & (latest.record_type == "paired_summary")
        & (latest.endpoint == "ALL"),
        "n_predictions",
    ].iloc[0]
    n_dedup = near_duplicate.loc[
        (near_duplicate.analysis_section == 1)
        & (near_duplicate.record_type == "paired_summary")
        & (near_duplicate.endpoint == "ALL"),
        "n_predictions",
    ].iloc[0]
    assert n_dedup < n_latest


def test_inputs_are_read_only_with_respect_to_e4b():
    require_e4b_data()
    before = {
        path: path.stat().st_mtime_ns
        for path in (BASELINE, RECORDS)
    }
    _ = build_posthoc()
    after = {
        path: path.stat().st_mtime_ns
        for path in (BASELINE, RECORDS)
    }
    assert before == after
