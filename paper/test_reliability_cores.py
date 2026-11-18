"""Independent checks for core evidence and example selection."""
import json

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem

import reliability_cores
from reliability_cores import (compare_fmcs, describe_core, most_frequent_core, run_cores_timing,
                               select_examples, summarize_core_timing, summarize_cores,
                               summarize_fmcs_recomputation)


def test_missing_core_is_not_zero_coverage_and_majority_is_strict():
    mols = [Chem.MolFromSmiles(s) for s in ("c1ccccc1", "CC")]
    missing, matched = describe_core("", mols, ["A", "B"])
    assert not missing["has_core"] and np.isnan(missing["coverage_fraction"])
    assert matched == [False, False]
    benzene, matched = describe_core("c1ccccc1", mols, ["A", "B"])
    assert matched == [True, False]
    assert benzene["coverage_fraction"] == .5 and not benzene["strict_majority"]
    assert benzene["is_benzene"] and benzene["is_single_ring"]
    assert benzene["core_heavy_atoms"] == 6 and benzene["core_ring_count"] == 1
    summary = summarize_cores(pd.DataFrame([missing, benzene]).assign(rule="test"), ["rule"]).iloc[0]
    assert summary.no_core == 1 and summary.mean_coverage == .5 and summary.median_core_heavy_atoms == 6


def test_fmcs_smiles_uses_actual_charge_from_member_and_keeps_query_authoritative():
    mols = [Chem.MolFromSmiles("C[N+](=O)[O-]")]
    core = "[#6]-[#7](=[#8])-[#8]"
    description, matched = describe_core(core, mols, ["nitro"], is_smarts=True)
    assert matched == [True] and description["core_smarts"] == core
    assert description["core_smiles"] == "C[N+](=O)[O-]"
    assert description["core_smiles_source_id"] == "nitro"
    assert description["core_heavy_atoms"] == 4


def test_example_ranking_prefers_larger_cores_then_coverage_gain():
    rows = []
    for name, heavy, current_count, aware_count in [("small", 6, 1, 10), ("large", 9, 4, 8),
                                                   ("largest_gain", 10, 1, 9), ("already_majority", 12, 6, 10)]:
        for rule, count in [("most_frequent", current_count), ("coverage", aware_count)]:
            rows.append({"target": "test", "target_label": "Test", "group": name, "size": 10,
                         "rule": rule, "core_smiles": "core", "core_heavy_atoms": heavy,
                         "core_ring_count": 1, "coverage_count": count, "coverage_fraction": count / 10,
                         "strict_majority": count > 5, "target_met": count >= 8,
                         "uncovered_member_ids": "A|B"})
    actual = select_examples(pd.DataFrame(rows))
    assert actual.group.tolist() == ["largest_gain", "large", "small"]
    assert actual.candidate_rank.tolist() == [1, 2, 3]


def test_fmcs_audit_detects_coverage_members_and_query_changes():
    mols = [Chem.MolFromSmiles(s) for s in ("c1ccccc1", "CC")]
    cached = pd.Series({"mcs_core": "[#6]", "mcs_heavy_atoms": 1, "mcs_coverage": 1.,
                        "mcs_timed_out": False, "mcs_seconds": .1})
    recomputed = {"core_smarts": "[#6]", "timed_out": False, "seconds": .1}
    assert compare_fmcs(cached, recomputed, mols, ["A", "B"])["passed"]
    recomputed["core_smarts"] = "[#6]-[#6]"
    result = compare_fmcs(cached, recomputed, mols, ["A", "B"])
    assert not result["passed"] and not result["same_member_matches"] and not result["same_heavy_atoms"]
    assert not result["same_coverage_count"]
    assert result["recomputed_coverage_fraction"] == pytest.approx(.5)


def test_timed_most_frequent_rule_matches_software_ties_and_missing_core():
    mols = [Chem.MolFromSmiles(s) for s in ("c1ccccc1", "c1ccc2ccccc2c1", "CC")]
    result = most_frequent_core(mols)
    assert result["core_smiles"] == "c1ccc2ccccc2c1"
    assert result["coverage_count"] == 1 and result["seconds"] >= 0
    assert most_frequent_core([mols[-1]])["core_smiles"] == ""
    data = pd.DataFrame({"molecule_identity_key": ["A", "B", "C"],
                         "canonical_smiles": [Chem.MolToSmiles(m) for m in mols], "group": "one"})
    _, groups, _ = reliability_cores.structural_groups(data, cluster_col="group", cluster_source="test")
    assert groups[0]["core_smiles"] == result["core_smiles"]
    assert groups[0]["core_coverage_count"] == result["coverage_count"]


def test_core_timer_summary_reconciles_current_history_and_timeout_denominators():
    table = pd.DataFrame({"target": ["a", "a", "b"], "target_label": ["A", "A", "B"],
        "most_frequent_seconds": [1., 2., 4.], "coverage_seconds": [2., 4., 8.],
        "fmcs_seconds": [3., 6., 12.], "fmcs_timed_out": [False, False, True],
        "history_aware_seconds": [1., 2., 3.], "history_mcs_seconds": [4., 5., 6.],
        "history_mcs_timed_out": [False, True, True]})
    summary = summarize_core_timing(table).set_index(["target", "rule"])
    assert len(summary) == 15
    current = summary.loc[("overall", "most_frequent")]
    assert current.groups == 3 and current.total_seconds == 7
    assert current.median_seconds == 2 and current.p90_seconds == pytest.approx(3.6)
    assert current.max_seconds == 4
    assert summary.loc[("overall", "fmcs"), "timeout_groups"] == 1
    assert summary.loc[("overall", "historical_fmcs"), "timeout_groups"] == 2
    assert "historical_most_frequent" not in summary.index.get_level_values("rule")


def test_full_fmcs_summary_separates_partial_timeout_and_completed_differences():
    audit = pd.DataFrame({"cached_timed_out": [True, False, False],
                          "recomputed_timed_out": [True, False, False],
                          "same_smarts": [False, False, True], "passed": [False, False, True]})
    result = summarize_fmcs_recomputation(audit)
    assert result["different_groups"] == 2 and result["exact_matches"] == 1
    assert result["differences_with_timeout"] == result["differences_without_timeout"] == 1
    assert result["difference_counts"] == {"same_smarts": 2}
    assert not result["non_timeout_reproduction_passed"]


def test_timing_runs_every_group_and_preserves_old_outputs_with_timeout_difference(tmp_path, monkeypatch):
    input_root, output_root = tmp_path / "input", tmp_path / "output"
    (input_root / "toy").mkdir(parents=True)
    (output_root / "cores").mkdir(parents=True)
    sentinel = output_root / "cores/core_checks.json"
    sentinel.write_text("old output stays byte-identical\n")
    data = pd.DataFrame({"molecule_identity_key": ["A", "B", "C", "D"],
                         "canonical_smiles": ["c1ccccc1", "c1ccccc1", "C1CCCCC1", "C1CCCCC1"]})
    data.to_csv(input_root / "toy/map_molecules.csv", index=False)
    frame, groups, _ = reliability_cores.structural_groups(data)
    old_rows = []
    for group in groups:
        members = frame.loc[frame.structure_group.eq(group["structure_group"])]
        mols = [Chem.MolFromSmiles(s) for s in members.canonical_smiles]
        aware, fmcs = reliability_cores.select_core(mols), reliability_cores.mcs_core(mols)
        old_rows.append({"target": "Toy", "group": group["structure_group"], "size": len(mols),
                         "current_core": group["core_smiles"], "aware_core": aware["core_smiles"],
                         "aware_seconds": aware["seconds"], "mcs_core": fmcs["core_smarts"],
                         "mcs_heavy_atoms": fmcs["core_heavy_atoms"], "mcs_coverage": fmcs["coverage_fraction"],
                         "mcs_seconds": fmcs["seconds"], "mcs_timed_out": False})
    cache = tmp_path / "cache.csv"
    pd.DataFrame(old_rows).to_csv(cache, index=False)
    monkeypatch.setattr(reliability_cores, "TARGETS", {"toy": "Toy"})
    monkeypatch.setattr(reliability_cores, "EXPECTED_GROUPS", {"toy": 2})
    original_mcs, calls = reliability_cores.mcs_core, []

    def simulate_timeout(mols, *, min_coverage, timeout):
        calls.append((min_coverage, timeout))
        result = original_mcs(mols, min_coverage=min_coverage, timeout=timeout)
        if len(calls) == 1:
            result.update(core_smarts="[#6]", core_heavy_atoms=1, timed_out=True)
        return result

    monkeypatch.setattr(reliability_cores, "mcs_core", simulate_timeout)
    result = run_cores_timing(output_root, input_root=input_root, cache_path=cache)
    assert calls == [(.8, 10), (.8, 10)]
    assert sentinel.read_text() == "old output stays byte-identical\n"
    checks = json.loads((output_root / "cores/core_timing_checks.json").read_text())
    assert checks == result["summary"] and checks["passed"] and checks["groups_total"] == 2
    assert checks["fmcs_recomputation"]["differences_with_timeout"] == 1
    assert checks["fmcs_recomputation"]["differences_without_timeout"] == 0
    audit = pd.read_csv(output_root / "cores/fmcs_full_recomputation.csv")
    timing = pd.read_csv(output_root / "cores/core_rule_timing.csv")
    assert len(audit) == len(timing) == 2 and audit.group.is_unique
    assert audit.recomputed_seconds.tolist() == timing.fmcs_seconds.tolist()
    assert audit.recomputed_coverage_count.tolist() == timing.fmcs_coverage_count.tolist()
    metadata = json.loads((output_root / "cores/core_timing_metadata.json").read_text())
    assert metadata["workers"] == 1 and metadata["platform"] and metadata["cpu_model"]
