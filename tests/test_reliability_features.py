import numpy as np
import pandas as pd
import pytest
from rdkit import Chem

from chemgridmap.cores import select_core
from chemgridmap.evidence import (annotate_activity_confidence, class_probabilities, document_values,
                                  fit_between_unit_sd, predictive_parameters)
from chemgridmap.refinement import neighbor_objective, refine_grid_layout, symmetric_neighbor_lists
from chemgridmap.reliability import neighbor_recall, representation_neighbor_sets


def test_tanimoto_neighbor_sets_keep_ties_and_exclude_self():
    bits = np.array([[1, 1, 0, 0], [1, 1, 0, 0], [1, 0, 1, 0], [0, 0, 1, 1]], dtype=bool)
    sets = representation_neighbor_sets(bits, k=1, metric="jaccard")
    assert list(sets[0]) == [1]
    assert 2 not in sets[2] or len(sets[2]) >= 1
    assert all(i not in row for i, row in enumerate(sets))
    tied = representation_neighbor_sets(np.array([[0, 0], [1, 0], [-1, 0], [0, 5.0]]), k=1, metric="euclidean")
    assert sorted(tied[0]) == [1, 2]
    assert neighbor_recall([np.array([1, 2])], [np.array([2, 3])])[0] == 0.5


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


def test_activity_confidence_report_is_strict_json():
    import json
    records = _records()
    molecules = records.groupby("molecule_identity_key").pchembl_value.median().rename("activity_pchembl").reset_index()
    out, fit = annotate_activity_confidence(molecules, records)
    json.dumps(fit, allow_nan=False)
    again, _ = annotate_activity_confidence(molecules, records, fit=fit)  # a reported fit can be reused
    assert np.allclose(again.class_confidence, out.class_confidence)
