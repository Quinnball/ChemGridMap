from pathlib import Path
import shutil
import subprocess

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem

from chemgridmap import __all__
from chemgridmap.cores import mcs_core, select_core
from chemgridmap.evidence import (annotate_activity_confidence, class_probabilities, document_values,
                                  fit_between_unit_sd, predictive_parameters)
from chemgridmap.refinement import neighbor_objective, refine_grid_layout, symmetric_neighbor_lists
from chemgridmap.reliability import (grid_neighbor_weights, neighbor_recall,
                                     neighbor_recall_fixed_k, representation_neighbor_sets,
                                     representation_neighbor_weights)
from chemgridmap.representations import morgan_fingerprints
from chemgridmap.projection import project_2d
from chemgridmap.provenance import environment_versions, git_provenance


def test_tanimoto_neighbor_sets_keep_ties_and_exclude_self():
    bits = np.array([[1, 1, 0, 0], [1, 1, 0, 0], [1, 0, 1, 0], [0, 0, 1, 1]], dtype=bool)
    sets = representation_neighbor_sets(bits, k=1, metric="jaccard")
    assert list(sets[0]) == [1]
    assert 2 not in sets[2] or len(sets[2]) >= 1
    assert all(i not in row for i, row in enumerate(sets))
    tied = representation_neighbor_sets(np.array([[0, 0], [1, 0], [-1, 0], [0, 5.0]]), k=1, metric="euclidean")
    assert sorted(tied[0]) == [1, 2]
    assert neighbor_recall([np.array([1, 2])], [np.array([2, 3])])[0] == 0.5
    neighbors, weights, kth_distance = representation_neighbor_weights(
        bits, k=1, metric="jaccard", return_kth=True
    )
    assert len(neighbors) == len(weights) == len(kth_distance) == len(bits)
    assert np.all((0 <= kth_distance) & (kth_distance <= 1))


def test_environment_versions_record_umap_numba_stack():
    packages = environment_versions()["packages"]
    assert {"numba", "pynndescent", "llvmlite"} <= set(packages)


def test_refinement_keeps_unique_cells_and_lowers_neighbor_objective():
    rng = np.random.default_rng(0)
    side, n = 12, 40
    cells = np.array([(c, r) for r in range(side) for c in range(side)])[rng.choice(side * side, n, replace=False)]
    # Neighbour pairs placed far apart should be pulled together.
    neighbors = [[(i + 1) % n, (i - 1) % n] for i in range(n)]
    adjacency = symmetric_neighbor_lists(neighbors)
    before = neighbor_objective(cells, adjacency, 2.0)
    diagnostics = {}
    refined = refine_grid_layout(cells, cells.astype(float), neighbors, side, sigmas=(2.0,), diagnostics=diagnostics)
    assert len({tuple(row) for row in refined}) == n
    assert refined.min() >= 0 and refined.max() < side
    assert neighbor_objective(refined, adjacency, 2.0) < before
    assert diagnostics["moved_molecules"] > 0
    again = refine_grid_layout(cells, cells.astype(float), neighbors, side, sigmas=(2.0,))
    assert np.array_equal(refined, again)


def test_refinement_rejects_duplicate_cells():
    with pytest.raises(ValueError):
        refine_grid_layout(np.array([[0, 0], [0, 0]]), np.zeros((2, 2)), [[1], [0]], 4)


def _records():
    rows = []
    for doc, value in [("D1", 8.6), ("D2", 8.4), ("D3", 8.5)]:
        rows.append({"molecule_identity_key": "A", "document_chembl_id": doc, "pchembl_value": value})
    for doc, value in [("D1", 6.9), ("D4", 7.1)]:
        rows.append({"molecule_identity_key": "B", "document_chembl_id": doc, "pchembl_value": value})
    for k in range(6):  # a re-reported value and many condition variants in one document
        rows.append({"molecule_identity_key": "C", "document_chembl_id": "D5", "pchembl_value": 5.0 + 0.1 * k})
    rows.append({"molecule_identity_key": "C", "document_chembl_id": "D6", "pchembl_value": 7.8})
    rows.append({"molecule_identity_key": "C", "document_chembl_id": "D7", "pchembl_value": 7.8})
    for e in ("E", "F", "G"):
        rows += [{"molecule_identity_key": e, "document_chembl_id": "D8", "pchembl_value": 6.0},
                 {"molecule_identity_key": e, "document_chembl_id": "D9", "pchembl_value": 7.0}]
    return pd.DataFrame(rows)


def test_document_units_merge_identical_reports():
    doc = document_values(_records())
    c = doc[doc.entry == "C"]
    assert len(c) == 2  # D6 and D7 report the same value and form one unit
    assert set(c.documents) == {"D5", "D6;D7"}


def test_class_probabilities_sum_to_one_and_follow_thresholds():
    p = class_probabilities(np.array([5.0, 6.5, 8.0]), np.array([0.5, 0.5, 0.5]), np.array([np.inf] * 3), 6.0, 7.0)
    assert np.allclose(p.sum(1), 1)
    assert p[0].argmax() == 0 and p[1].argmax() == 1 and p[2].argmax() == 2
    loc, scale, df = predictive_parameters([7.0], [0.0], [1.0], 0.8, np.inf)
    assert scale[0] == pytest.approx(0.8 * np.sqrt(2))


def test_activity_confidence_keeps_display_and_reports_support():
    records = _records()
    molecules = records.groupby("molecule_identity_key").pchembl_value.median().rename("activity_pchembl").reset_index()
    out, fit = annotate_activity_confidence(molecules, records, model="gaussian")
    assert fit["method"] == "pooled_between_unit_sd" and fit["n_entries"] == 6
    assert list(out.activity_pchembl) == list(molecules.activity_pchembl)
    a = out.set_index("molecule_identity_key")
    assert a.loc["A", "class_status"] == "supported"  # three documents agree well above 7
    assert a.loc["B", "class_status"] == "uncertain"  # straddles the upper threshold
    assert a.loc["C", "n_evidence_units"] == 2 and a.loc["C", "n_documents"] == 3
    assert a.loc["C", "document_balanced_pchembl"] == pytest.approx((5.25 + 7.8) / 2)
    assert 0 <= a.class_confidence.min() and a.class_confidence.max() <= 1


def test_offset_pair_rule_keeps_probabilities_and_changes_only_display_status():
    rows = []
    for index in range(5):
        rows.extend([
            {"molecule_identity_key": f"TOUCHED{index}", "document_chembl_id": "D1",
             "pchembl_value": 5.5},
            {"molecule_identity_key": f"TOUCHED{index}", "document_chembl_id": "D2",
             "pchembl_value": 7.5},
        ])
    rows.append({"molecule_identity_key": "SINGLE_DOCUMENT", "document_chembl_id": "D1",
                 "pchembl_value": 8.0})
    rows.append({"molecule_identity_key": "CLEAN", "document_chembl_id": "D3",
                 "pchembl_value": 8.0})
    records = pd.DataFrame(rows)
    molecules = (
        records.groupby("molecule_identity_key").pchembl_value.median()
        .rename("activity_pchembl").reset_index()
    )
    out, fit = annotate_activity_confidence(molecules, records)
    indexed = out.set_index("molecule_identity_key")
    touched = indexed.loc[[f"TOUCHED{i}" for i in range(5)]]
    assert touched.offset_pair_touched.all()
    assert touched.status_reason.eq("offset_pair").all()
    assert touched.class_status.eq("uncertain").all()
    assert fit["offset_pair_touched_entry_count"] == 6
    assert indexed.loc["SINGLE_DOCUMENT", "discordant_document_pairs"] == 0
    assert indexed.loc["SINGLE_DOCUMENT", "offset_pair_touched"]
    assert indexed.loc["SINGLE_DOCUMENT", "status_reason"] == "offset_pair"
    assert indexed.loc["CLEAN", "status_reason"] in {"supported", "low_support"}
    assert {"class_confidence", "class_reproducibility"} <= set(out)
    assert fit["sigma_doc"] == fit["tau"]
    assert fit["document_offsets"][0]["touched_entry_count"] == 6


def test_offset_display_none_keeps_offsets_but_uses_support_only():
    records = _records()
    molecules = records.groupby("molecule_identity_key").pchembl_value.median().rename("activity_pchembl").reset_index()
    default, _ = annotate_activity_confidence(molecules, records)
    unchanged, fit = annotate_activity_confidence(
        molecules, records, offset_display="none"
    )
    indexed = unchanged.set_index("molecule_identity_key")
    expected = np.where(
        default.set_index("molecule_identity_key")["class_confidence"] >= 0.8,
        "supported", "uncertain"
    )
    assert indexed.loc[default.set_index("molecule_identity_key").index, "class_status"].tolist() == expected.tolist()
    assert set(indexed.status_reason) <= {"supported", "low_support"}
    assert indexed.offset_pair_touched.any()
    assert fit["offset_display"] == "none"


def test_offset_min_shared_filters_pairs_and_support_is_monotone():
    rows = []
    for index in range(3):
        rows.extend([
            {"molecule_identity_key": f"E{index}", "document_chembl_id": "D1", "pchembl_value": 5.5},
            {"molecule_identity_key": f"E{index}", "document_chembl_id": "D2", "pchembl_value": 7.5},
        ])
    records = pd.DataFrame(rows)
    molecules = records.groupby("molecule_identity_key").pchembl_value.median().rename("activity_pchembl").reset_index()
    loose, _ = annotate_activity_confidence(molecules, records, offset_min_shared=3)
    strict, _ = annotate_activity_confidence(molecules, records, offset_min_shared=5)
    assert loose.offset_pair_touched.all()
    assert not strict.offset_pair_touched.any()
    low, _ = annotate_activity_confidence(molecules, records, confidence_threshold=0.8)
    high, _ = annotate_activity_confidence(molecules, records, confidence_threshold=0.9)
    assert set(high.loc[high.class_status.eq("supported"), "molecule_identity_key"]) <= set(
        low.loc[low.class_status.eq("supported"), "molecule_identity_key"]
    )


@pytest.mark.parametrize("kwargs,message", [
    ({"offset_min_shared": 1}, "at least 2"),
    ({"offset_threshold": 0}, "greater than 0"),
    ({"offset_sign_agreement": 0.5}, r"\(0.5, 1\]"),
])
def test_offset_parameter_validation(kwargs, message):
    records = _records()
    molecules = records.groupby("molecule_identity_key").pchembl_value.median().rename("activity_pchembl").reset_index()
    with pytest.raises(ValueError, match=message):
        annotate_activity_confidence(molecules, records, **kwargs)


def test_between_unit_sd_falls_back_to_literature_value_for_small_inputs():
    doc = document_values(_records().query("molecule_identity_key in ['A', 'B']"))
    fit = fit_between_unit_sd(doc)
    assert fit["method"] == "literature_fallback" and fit["tau"] == pytest.approx(0.64)


def test_coverage_aware_core_prefers_largest_core_meeting_target():
    smiles = ["CC(=O)Oc1ccc2nc(S(N)(=O)=O)sc2c1", "NS(=O)(=O)c1nc2ccc(OCC)cc2s1",
              "NS(=O)(=O)c1nc2ccc(OC(=O)c3ccccc3)cc2s1", "CC(=O)Oc1ccc2sc(S(N)(=O)=O)cc2c1"]
    mols = [Chem.MolFromSmiles(s) for s in smiles]
    choice = select_core(mols, min_coverage=0.75)
    assert choice["core_smiles"] == "c1ccc2scnc2c1" and choice["coverage_count"] == 3 and choice["target_met"]
    strict = select_core(mols, min_coverage=1.0)
    assert strict["coverage_count"] == 4 and strict["core_heavy_atoms"] == 6  # only benzene covers all


def test_fmcs_core_exports_valid_attachment_smiles():
    mols = [Chem.MolFromSmiles(smiles) for smiles in [
        "c1ccccc1CC", "c1ccccc1CO", "c1ccccc1CN",
    ]]
    result = mcs_core(mols, min_coverage=1.0, timeout=10)
    assert result["core_smiles_attached"]
    attached = Chem.MolFromSmiles(result["core_smiles_attached"])
    assert attached is not None
    assert any(atom.GetAtomicNum() == 0 for atom in attached.GetAtoms())


def test_fmcs_core_search_is_nonchiral_but_containment_uses_chirality():
    mols = [Chem.MolFromSmiles(value) for value in ["C[C@H](O)F", "C[C@@H](O)F"]]
    result = mcs_core(mols, min_coverage=1.0, timeout=10)
    pattern = Chem.MolFromSmarts(result["core_smarts"])
    assert result["core_heavy_atoms"] == 4
    assert result["coverage_count"] == 2
    assert all(atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED for atom in pattern.GetAtoms())
    assert all(mol.HasSubstructMatch(pattern, useChirality=True) for mol in mols)


def _temporary_git_package(tmp_path, *, track):
    repo = tmp_path / "repo"
    package = repo / "src" / "sample_package"
    package.mkdir(parents=True)
    (package / "provenance.py").write_text("# test module\n", encoding="utf-8")
    (package / "module.py").write_text("VALUE = 1\n", encoding="utf-8")
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    subprocess.run(["git", "-C", str(repo), "config", "user.email", "test@example.invalid"], check=True)
    subprocess.run(["git", "-C", str(repo), "config", "user.name", "Test"], check=True)
    if track:
        subprocess.run(["git", "-C", str(repo), "add", "src/sample_package"], check=True)
        subprocess.run(["git", "-C", str(repo), "commit", "-q", "-m", "initial"], check=True)
    return repo, package


def test_git_provenance_ignores_an_untracked_module_in_a_repository(tmp_path):
    if shutil.which("git") is None:
        pytest.skip("git is not available")
    _, package = _temporary_git_package(tmp_path, track=False)
    assert git_provenance(package) == {"git_commit": None, "git_dirty": None}


def test_git_provenance_tracks_clean_and_modified_package_only(tmp_path):
    if shutil.which("git") is None:
        pytest.skip("git is not available")
    _, package = _temporary_git_package(tmp_path, track=True)
    clean = git_provenance(package)
    assert len(clean["git_commit"]) == 40
    assert clean["git_dirty"] is False

    (package / "module.py").write_text("VALUE = 2\n", encoding="utf-8")
    assert git_provenance(package)["git_dirty"] is True


def test_fixed_k_weights_sum_to_k_and_expected_recall_uses_tie_probabilities():
    reference = [{1: 1.0, 2: 0.5, 3: 0.5}]
    displayed = [{2: 0.5, 3: 0.5, 4: 1.0}]
    assert neighbor_recall_fixed_k(reference, displayed, k=2)[0] == pytest.approx(0.25)
    sets, weights = representation_neighbor_weights(
        np.array([[0, 0], [1, 0], [-1, 0], [0, 5.0]]), k=1, metric="euclidean"
    )
    assert all(sum(row.values()) == pytest.approx(1.0) for row in weights)
    assert sorted(sets[0]) == [1, 2]


def test_caii_fixed_k_recall_matches_e6_r1():
    root = Path(__file__).resolve().parents[1]
    coordinates_path = (
        root / "paper/output/reliability_v1/layouts/chembl205/coordinates_seed_42.csv"
    )
    if not coordinates_path.exists():
        pytest.skip("The frozen CA II coordinates are available in the research checkout.")
    coordinates = pd.read_csv(coordinates_path)
    fingerprints = morgan_fingerprints(coordinates.canonical_smiles.tolist())
    reference, reference_weights = representation_neighbor_weights(
        fingerprints, k=10, metric="jaccard"
    )
    displayed, displayed_weights = grid_neighbor_weights(
        coordinates[["refined_grid_col", "refined_grid_row"]].to_numpy(float), k=10
    )
    retained = neighbor_recall(reference, displayed)
    fixed = neighbor_recall_fixed_k(reference_weights, displayed_weights, k=10)
    assert np.nanmean(retained) == pytest.approx(0.5845, abs=1e-3)
    assert np.nanmean(fixed) == pytest.approx(0.5359, abs=1e-3)


def test_umap_disconnected_graph_is_repeatable_with_same_seed():
    matrix = np.array([
        [0, 0], [0, 1], [0, 2], [0, 3],
        [100, 0], [100, 1], [100, 2], [100, 3],
    ], dtype=float)
    first = project_2d(matrix, method="umap", random_state=42, umap_metric="euclidean")
    second = project_2d(matrix, method="umap", random_state=42, umap_metric="euclidean")
    assert np.array_equal(first, second)


def test_public_api_exports_evidence_neighbors_and_cores():
    expected = {
        "document_values", "fit_between_unit_sd", "annotate_activity_confidence",
        "document_pair_offsets", "neighbor_recall", "neighbor_recall_fixed_k",
        "select_core", "mcs_core", "structural_groups",
    }
    assert expected <= set(__all__)


def test_activity_confidence_report_is_strict_json():
    import json
    records = _records()
    molecules = records.groupby("molecule_identity_key").pchembl_value.median().rename("activity_pchembl").reset_index()
    out, fit = annotate_activity_confidence(molecules, records)
    json.dumps(fit, allow_nan=False)
    again, _ = annotate_activity_confidence(molecules, records, fit=fit)  # a reported fit can be reused
    assert np.allclose(again.class_confidence, out.class_confidence)


def test_exclude_flagged_evidence_model_reports_flagged_documents_and_fallback():
    rows = []
    for entry in ("A", "B", "C"):
        rows.extend([
            {"molecule_identity_key": entry, "document_chembl_id": "D1", "pchembl_value": 5.0},
            {"molecule_identity_key": entry, "document_chembl_id": "D2", "pchembl_value": 6.5},
        ])
    rows.extend([
        {"molecule_identity_key": "FALLBACK", "document_chembl_id": "D2", "pchembl_value": 7.0},
        {"molecule_identity_key": "UNIQUE", "document_chembl_id": "D3", "pchembl_value": 6.0},
    ])
    records = pd.DataFrame(rows)
    molecules = records.groupby("molecule_identity_key").pchembl_value.median().rename(
        "activity_pchembl"
    ).reset_index()
    out, fit = annotate_activity_confidence(
        molecules, records, model="exclude-flagged"
    )
    assert fit["model"] == "exclude_offset_documents"
    assert fit["flagged_documents"] == ["D1", "D2"]
    assert fit["excluded_entry_count"] == 4
    assert out.set_index("molecule_identity_key").loc["FALLBACK", "n_evidence_units"] == 1


def test_exclude_flagged_empty_entry_reuses_global_gaussian_support():
    rows = []
    for index in range(5):
        rows.extend([
            {"molecule_identity_key": f"EXCLUDED{index}", "document_chembl_id": "D1",
             "pchembl_value": 5.0},
            {"molecule_identity_key": f"EXCLUDED{index}", "document_chembl_id": "D2",
             "pchembl_value": 6.5},
        ])
    rows.extend([
        {"molecule_identity_key": "FALLBACK", "document_chembl_id": "D1", "pchembl_value": 5.0},
        {"molecule_identity_key": "FALLBACK", "document_chembl_id": "D2", "pchembl_value": 7.5},
    ])
    for index in range(5):
        rows.extend([
            {"molecule_identity_key": f"KEPT{index}", "document_chembl_id": "D3",
             "pchembl_value": 6.0},
            {"molecule_identity_key": f"KEPT{index}", "document_chembl_id": "D4",
             "pchembl_value": 6.2},
        ])
    records = pd.DataFrame(rows)
    molecules = (
        records.groupby("molecule_identity_key").pchembl_value.median()
        .rename("activity_pchembl").reset_index()
    )
    gaussian, _ = annotate_activity_confidence(molecules, records)
    excluded, fit = annotate_activity_confidence(molecules, records, model="exclude-flagged")
    support = [
        "class_confidence", "class_reproducibility",
        "p_consensus_inactive", "p_consensus_medium", "p_consensus_active",
        "p_new_document_inactive", "p_new_document_medium", "p_new_document_active",
    ]
    fallback = "FALLBACK"
    actual = excluded.set_index("molecule_identity_key").loc[fallback, support]
    expected = gaussian.set_index("molecule_identity_key").loc[fallback, support]
    assert fit["excluded_entry_count"] == 6
    pd.testing.assert_series_equal(actual, expected, check_names=False)
