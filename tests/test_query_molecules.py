"""Unmeasured structures must not alter measured evidence or group definitions."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

from chemgridmap.chembl import curate_chembl_activity_data
from chemgridmap.core import label_purity, mean_knn_absolute_difference, prepare_molecules
from chemgridmap.evidence import annotate_activity_confidence
from chemgridmap.pipeline import build_grid_map
from chemgridmap.queries import assign_query_groups, prepare_query_molecules


def small_reference():
    return pd.DataFrame({
        "molecule_id": [f"M{i}" for i in range(6)],
        "smiles": ["CCC", "CCO", "CCN", "CCCC", "CCCl", "CCBr"],
        "value": [5.0, 6.0, 7.0, 8.0, 5.5, 6.5],
    })


def test_check_records_every_source_row_and_does_not_copy_activity():
    reference, _ = prepare_molecules(small_reference(), value_col="value")
    query = pd.DataFrame({"SMILES": ["CCF", "not a smiles", "FCC", "CCC"],
                          "activity_pchembl": [99] * 4, "class_confidence": [1] * 4})
    mapped, check, report = prepare_query_molecules(query, reference)
    assert check.status.tolist() == ["mapped", "invalid_smiles", "duplicate_query", "in_reference"]
    assert check.source_row.tolist() == [0, 1, 2, 3]
    assert check.query_id.tolist() == ["Q00001", "Q00002", "Q00003", "Q00004"]
    assert check.iloc[2].first_query_id == "Q00001"
    assert check.iloc[3].reference_id == "M0"
    assert {s: report[s] for s in check.status} == {s: 1 for s in check.status}
    assert mapped.activity_class.tolist() == ["not_measured"]
    assert mapped.measurement_status.tolist() == ["not_measured"]
    assert "activity_pchembl" not in mapped and "class_confidence" not in mapped


def test_query_identifiers_and_missing_smiles_fail_clearly():
    reference, _ = prepare_molecules(small_reference(), value_col="value")
    with pytest.raises(ValueError, match="duplicate query IDs"):
        prepare_query_molecules(pd.DataFrame({"smiles": ["CCF", "CCCF"], "id": ["q", "q"]}),
                                reference, id_column="id")
    with pytest.raises(ValueError, match="SMILES column"):
        prepare_query_molecules(pd.DataFrame({"wrong": ["CCF"]}), reference)


def test_reference_parent_policy_applies_to_queries():
    reference = pd.DataFrame({"molecule_id": ["amine", "ethanol"],
                              "canonical_smiles": ["CCN", "CCO"],
                              "structure_policy": ["parent", "parent"]})
    mapped, check, report = prepare_query_molecules(pd.DataFrame({"smiles": ["CC[NH3+].[Cl-]"]}), reference)
    assert mapped.empty
    assert check.status.tolist() == ["in_reference"]
    assert check.canonical_smiles.tolist() == ["CCN"]
    assert report["structure_policy"] == "parent"


@pytest.mark.parametrize("representation", ["morgan", "descriptors"])
def test_measured_neighbors_match_brute_force_and_include_fifth_place_ties(representation):
    reference = small_reference()
    result = build_grid_map(reference, value_col="value", projection="pca",
                            representation=representation, grid_refinement="none", grouping="none",
                            radius=1, n_bits=256,
                            query_data=pd.DataFrame({"smiles": ["[Na+]", "CCF"]}))
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=1, fpSize=256)
    fingerprints = [generator.GetFingerprint(Chem.MolFromSmiles(s)) for s in reference.smiles]
    for query in result.data[result.data.measurement_status.eq("not_measured")].itertuples():
        scores = np.asarray(DataStructs.BulkTanimotoSimilarity(
            generator.GetFingerprint(Chem.MolFromSmiles(query.canonical_smiles)), fingerprints))
        cutoff = sorted(scores, reverse=True)[4]
        expected = {reference.molecule_id.iloc[i]: score for i, score in enumerate(scores) if score >= cutoff}
        actual = result.query_neighbors[result.query_neighbors.query_id.eq(query.query_id)]
        assert dict(zip(actual.measured_id, actual.tanimoto)) == pytest.approx(expected, abs=1e-15)
        assert query.nearest_measured_similarity == pytest.approx(scores.max())
        assert query.nearest_measured_id in reference.molecule_id[scores == scores.max()].tolist()
    tied = result.query_neighbors[result.query_neighbors.query_id.eq("Q00001")]
    assert len(tied) == 6 and tied["rank"].eq(1).all()
    assert result.query_neighbors.class_confidence.isna().all()


@pytest.fixture(scope="module")
def caii_reference_and_queries():
    source = Path(__file__).resolve().parents[1] / "src" / "chemgridmap" / "data" / "chembl205_ic50.csv"
    curated = curate_chembl_activity_data(pd.read_csv(source), target_id="CHEMBL205")
    annotated, fit = annotate_activity_confidence(
        curated.molecules, curated.retained_records, value="activity_pchembl_raw",
        document=curated.report["detected_columns"]["document_id"])
    query = pd.DataFrame({"smiles": ["[Xe]", "[Ne]", "[Xe]", "invalid", annotated.canonical_smiles.iloc[0]],
                          "class_confidence": [1.0] * 5, "activity_class": ["active"] * 5,
                          "activity_pchembl": [99.0] * 5})
    kwargs = dict(smiles_col="canonical_smiles", value_col="activity_pchembl", label_col="activity_class",
                  projection="pca", grid_refinement="none", core_rule="coverage", random_state=42)
    return (build_grid_map(annotated, **kwargs),
            build_grid_map(annotated, query_data=query, **kwargs), fit)


def test_caii_queries_preserve_every_reference_evidence_and_group_value(caii_reference_and_queries):
    baseline, with_queries, fit = caii_reference_and_queries
    measured = with_queries.data.iloc[:len(baseline.data)]
    varying = {"projection_x_raw", "projection_y_raw", "projection_x", "projection_y",
               "grid_x", "grid_y", "grid_col", "grid_row"}
    preserved = [c for c in baseline.data if c not in varying and "neighbor" not in c]
    pd.testing.assert_frame_equal(baseline.data[preserved], measured[preserved], check_dtype=False)
    for before, after in zip(baseline.structure_groups, with_queries.structure_groups):
        assert before == {key: value for key, value in after.items() if key != "n_query_members"}
    assert baseline.grouping_protocol == with_queries.grouping_protocol
    assert with_queries.query_check.status.tolist() == ["mapped", "mapped", "duplicate_query", "invalid_smiles", "in_reference"]
    unmeasured = with_queries.data.iloc[len(baseline.data):]
    assert unmeasured.class_confidence.isna().all()
    assert unmeasured.activity_pchembl.isna().all()
    assert unmeasured.class_status.isna().all()
    assert unmeasured.activity_class.eq("not_measured").all()


def test_annotation_metrics_use_only_measured_rows(caii_reference_and_queries):
    baseline, result, _ = caii_reference_and_queries
    n = len(baseline.data)
    expected_purity = label_purity(result.grid[:n], baseline.data.activity_class.to_numpy(), k=10)
    expected_difference = mean_knn_absolute_difference(result.grid[:n], baseline.data.activity_pchembl.to_numpy(), k=10)
    assert result.metrics.iloc[0].grid_label_purity == pytest.approx(expected_purity)
    assert result.metrics.iloc[0].grid_knn_value_abs_diff == pytest.approx(expected_difference)
    assert result.metrics.iloc[0].n_samples == n + 2
    assert result.metrics.iloc[0].n_measured == n
    assert result.query_report["metric_populations"]["activity_and_class"]["n_molecules"] == n


@pytest.mark.parametrize("group_ids", [("G002", "G001"), ("G1000", "G999")])
def test_query_group_ties_select_smallest_group_and_check_core(group_ids):
    queries = pd.DataFrame({"molecule_identity_key": ["query:q"], "canonical_smiles": ["Fc1ccccc1"]})
    groups = [dict(structure_group=group_ids[0], representative_smiles="c1ccccc1", color="#000000",
                   core_smiles="c1ccccc1"),
              dict(structure_group=group_ids[1], representative_smiles="c1ccccc1", color="#ffffff",
                   core_smiles="c1ccccc1")]
    mapped, summaries = assign_query_groups(queries, groups, grouping="butina", similarity=0.1,
                                            radius=2, n_bits=2048)
    assert mapped.structure_group.tolist() == [group_ids[1]]
    assert mapped.contains_group_core.tolist() == [True]
    assert mapped.query_group_center_similarity.iloc[0] >= 0.1
    assert [group["n_query_members"] for group in summaries] == [0, 1]
    assert all("n_query_members" not in group for group in groups)


def test_query_outputs_and_external_groups(tmp_path):
    reference = small_reference().assign(group=["A", "A", "A", "B", "B", "B"])
    result = build_grid_map(reference, value_col="value", projection="pca", grid_refinement="none",
                            grouping="external", cluster_col="group", cluster_source="test groups",
                            query_data=pd.DataFrame({"smiles": ["CCF"]}),
                            output_dir=tmp_path, output_formats=["svg"], render_detail="overview")
    query = result.data.iloc[-1]
    assert query.structure_group == ""
    assert all(group["n_query_members"] == 0 for group in result.structure_groups)
    manifest = json.loads(Path(result.output_files["run_manifest"]).read_text(encoding="utf-8"))
    assert manifest["query_molecules"]["mapped"] == 1
    assert manifest["metric_populations"]["structure_and_layout"]["n_molecules"] == 7
    assert any("External grouping" in warning for warning in manifest["warnings"])
    check = pd.read_csv(result.output_files["query_check"])
    assert check.status.tolist() == ["mapped"]
    members = pd.read_csv(result.output_files["structure_members"])
    assert not members.iloc[-1].core_member


def test_queries_need_measured_references():
    with pytest.raises(ValueError, match="at least two measured"):
        build_grid_map(pd.DataFrame({"smiles": ["CCC", "CCO"]}), projection="pca",
                       query_data=pd.DataFrame({"smiles": ["CCF"]}))


@pytest.mark.parametrize("missing", [None, np.nan, "", "  ", "not_measured", "<NA>"])
def test_missing_reference_labels_are_not_measured(missing):
    with pytest.raises(ValueError, match="at least two measured"):
        build_grid_map(pd.DataFrame({"smiles": ["CCC", "CCO"], "label": [missing, missing]}),
                       label_col="label", projection="pca",
                       query_data=pd.DataFrame({"smiles": ["CCF"]}))


def test_reference_rows_without_measurement_belong_in_query_table():
    with pytest.raises(ValueError, match="must be supplied in query_data"):
        build_grid_map(pd.DataFrame({"smiles": ["CCC", "CCO", "CCN"], "value": [5.0, 6.0, np.nan]}),
                       value_col="value", projection="pca",
                       query_data=pd.DataFrame({"smiles": ["CCF"]}))
