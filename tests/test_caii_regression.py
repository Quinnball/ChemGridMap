"""Public regression checks on the bundled CA II export."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from chemgridmap.app import AppSession, DEMO
from chemgridmap.columns import describe


@pytest.fixture(scope="module")
def caii_runs(tmp_path_factory):
    runs = {}
    for projection in ("umap", "pca"):
        session = AppSession(tmp_path_factory.mktemp(projection))
        session.load(DEMO.read_bytes(), DEMO.name)
        session.start(
            {
                "target_id": "CHEMBL205",
                "projection": projection,
                "output_formats": ["svg"],
                "grouping": "butina",
                "render_detail": "overview",
            },
            background=False,
        )
        assert session.job["state"] == "done", session.job
        runs[projection] = session
    return runs


def test_caii_default_reliability_values(caii_runs):
    data = caii_runs["umap"].result.data
    fit = caii_runs["umap"].curation.report["activity_confidence"]
    assert len(data) == 889
    assert data.class_status.eq("supported").sum() == 506
    assert ((data.class_status.eq("supported")) & data.activity_class.eq("active")).sum() == 390
    assert ((data.class_status.eq("supported")) & data.activity_class.eq("inactive")).sum() == 116
    assert data.offset_pair_touched.sum() == 0
    assert fit["sigma_doc"] == pytest.approx(0.796, abs=5e-4)
    assert fit["sigma_doc_degrees_of_freedom"] == 64
    assert fit["multi_unit_entry_count"] == 29
    assert np.nanmean(data.neighbor_recall) == pytest.approx(0.585, abs=0.02)
    assert np.nanmean(data.neighbor_recall_fixed_k) == pytest.approx(0.536, abs=0.02)


def test_caii_pca_grid_coordinates_match_expectation(caii_runs):
    expected_path = Path(__file__).with_name("data") / "caii_pca_seed42_grid.csv"
    expected = pd.read_csv(expected_path, encoding="utf-8-sig")
    actual = pd.read_csv(
        caii_runs["pca"].result.output_files["coordinates"], encoding="utf-8-sig"
    )
    pd.testing.assert_frame_equal(
        actual[["molecule_identity_key", "grid_row", "grid_col"]],
        expected,
        check_dtype=False,
    )


def test_caii_output_columns_are_documented(caii_runs):
    files = caii_runs["umap"].files
    mapping = {
        "coordinate": "grid_coordinates.csv",
        "molecule": "molecule_level.csv",
        "metrics": "metrics.csv",
        "structure": "structure_members.csv",
        "offsets": "document_offsets.csv",
        "retained": "retained_records.csv",
        "excluded": "excluded_records.csv",
    }
    for key, label in mapping.items():
        candidates = [
            path for name, path in files.items()
            if name.startswith("map_") and name.endswith(f"_{key}.csv")
        ]
        if label == "grid_coordinates.csv":
            candidates = [files["map_grid_coordinates.csv"]]
        elif label == "molecule_level.csv":
            candidates = [files["map_molecule_level.csv"]]
        elif label == "metrics.csv":
            candidates = [files["map_metrics.csv"]]
        elif label == "structure_members.csv":
            candidates = [files["map_structure_members.csv"]]
        elif label == "document_offsets.csv":
            candidates = [files["map_document_offsets.csv"]]
        elif label == "retained_records.csv":
            candidates = [files["map_retained_records.csv"]]
        elif label == "excluded_records.csv":
            candidates = [files["map_excluded_records.csv"]]
        assert candidates, label
        describe(pd.read_csv(candidates[0], nrows=0).columns, label)
