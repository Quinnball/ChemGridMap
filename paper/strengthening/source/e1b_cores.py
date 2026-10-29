"""E1-B: 冻结核心规则分层，以及顺序单进程 FMCS 时间预算扫描。"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

sys.dont_write_bytecode = True

import numpy as np
import pandas as pd
from rdkit import Chem, rdBase

ROOT = Path(__file__).resolve().parents[3]
CGM = ROOT
OUTPUT = Path(__file__).resolve().parents[1] / "output/E1B"
sys.path.insert(0, str(CGM / "src"))
sys.path.insert(0, str(CGM / "paper"))
from chemgridmap.comparison import structural_groups
from chemgridmap.cores import mcs_core
from reliability_cores import describe_core

TARGETS = {"chembl205": "CA II", "chembl204": "Thrombin", "chembl240": "hERG"}
EXPECTED = {"chembl205": 125, "chembl204": 152, "chembl240": 203}
STRATA = ("2", "3-4", "5-9", ">=10")


def size_stratum(size):
    """分组规则在任务单中预先指定；单分子组不属于分析总体。"""
    if size < 2:
        raise ValueError("Only nonsingleton groups are eligible.")
    return "2" if size == 2 else "3-4" if size <= 4 else "5-9" if size <= 9 else ">=10"


def stratify_cores(table):
    """No core 是分类状态，覆盖率与核心大小均排除缺失核心。"""
    rows = []
    groups = [(target, label, table.loc[table.target.eq(target)]) for target, label in TARGETS.items()]
    groups.append(("overall", "Overall", table))
    for target, label, population in groups:
        for stratum in STRATA:
            for rule in ("most_frequent", "coverage", "fmcs"):
                subset = population.loc[population["size"].map(size_stratum).eq(stratum) & population.rule.eq(rule)]
                present = subset.loc[subset.has_core]
                rows.append({"target": target, "target_label": label, "size_stratum": stratum,
                             "rule": rule, "groups": len(subset), "no_core_groups": int((~subset.has_core).sum()),
                             "strict_majority_groups": int(subset.strict_majority.sum()),
                             "target_met_groups": int(subset.target_met.sum()),
                             "mean_coverage": present.coverage_fraction.mean(),
                             "mean_core_heavy_atoms": present.core_heavy_atoms.mean(),
                             "benzene_cores": int(subset.is_benzene.sum())})
    return pd.DataFrame(rows)


def summarize_budgets(table):
    rows = []
    populations = [(target, label, table.loc[table.target.eq(target)]) for target, label in TARGETS.items()]
    populations.append(("overall", "Overall", table))
    for target, label, population in populations:
        for timeout, subset in population.groupby("timeout_seconds", sort=True):
            rows.append({"target": target, "target_label": label, "timeout_seconds": int(timeout),
                         "groups": len(subset), "target_met_groups": int(subset.target_met.sum()),
                         "no_core_groups": int((~subset.has_core).sum()),
                         "mean_coverage": subset.loc[subset.has_core, "coverage_fraction"].mean(),
                         "timeout_groups": int(subset.timed_out.sum()),
                         "total_seconds": float(subset.seconds.sum()),
                         "median_seconds": float(subset.seconds.median()),
                         "source": subset.source.iloc[0]})
    return pd.DataFrame(rows)


def member_ids(value):
    return str(value).split("|") if pd.notna(value) and str(value) else []


def prepare_groups(comparison, timing, audit):
    """重建固定分组，核对每组完整成员以及三条冻结规则的匹配证据。"""
    prepared, evidence, protocols = [], [], {}
    keys = ["target", "group"]
    frozen = comparison.set_index(keys + ["rule"])
    timer = timing.set_index(keys)
    old = audit.set_index(keys)
    if len(comparison) != 1440 or len(timing) != 480 or len(audit) != 480:
        raise ValueError("Expected 480 groups and 1440 frozen rule rows.")
    if not (frozen.index.is_unique and timer.index.is_unique and old.index.is_unique):
        raise ValueError("Duplicate frozen group keys.")
    for target, label in TARGETS.items():
        data = pd.read_csv(CGM / f"paper/output/revision_v6/{target}/map_molecules.csv", low_memory=False)
        data = data.sort_values("molecule_identity_key").reset_index(drop=True)
        frame, groups, protocols[target] = structural_groups(data, similarity=.6, core_rule="most_frequent")
        groups = [group for group in groups if group["size"] >= 2]
        if len(groups) != EXPECTED[target]:
            raise ValueError(f"Unexpected group count: {target} {len(groups)}")
        for group in groups:
            group_id = group["structure_group"]
            members = frame.loc[frame.structure_group.eq(group_id)]
            ids = members.molecule_identity_key.astype(str).tolist()
            mols = [Chem.MolFromSmiles(smiles) for smiles in members.canonical_smiles]
            assert all(mol is not None for mol in mols)
            cache_fmcs = frozen.loc[(target, group_id, "fmcs")]
            prior = old.loc[(target, group_id)]
            prior_ids = member_ids(prior.cached_matched_member_ids) + member_ids(cache_fmcs.uncovered_member_ids)
            assert len(prior_ids) == len(ids) == len(set(prior_ids))
            assert set(prior_ids) == set(ids)
            assert int(timer.loc[(target, group_id), "size"]) == len(ids) == int(prior["size"])
            for rule in ("most_frequent", "coverage", "fmcs"):
                previous = frozen.loc[(target, group_id, rule)]
                assert int(previous["size"]) == len(ids)
                core_field = "core_smarts" if rule == "fmcs" else "core_smiles"
                core = previous[core_field] if pd.notna(previous[core_field]) else ""
                description, matches = describe_core(core, mols, ids, is_smarts=rule == "fmcs")
                for field in ("has_core", "core_status", "core_heavy_atoms", "coverage_count", "strict_majority", "target_met", "is_benzene"):
                    assert description[field] == previous[field], (target, group_id, rule, field)
                assert description["uncovered_member_ids"] == "|".join(member_ids(previous.uncovered_member_ids))
                assert np.isclose(description["coverage_fraction"], previous.coverage_fraction, atol=1e-14, rtol=0, equal_nan=True)
            common = {"target": target, "target_label": label, "group": group_id, "size": len(ids)}
            prepared.append((common, mols))
            evidence.append({**common, "size_stratum": size_stratum(len(ids)),
                             "member_ids": "|".join(ids),
                             "member_ids_sha256": hashlib.sha256("\n".join(ids).encode()).hexdigest(),
                             "same_frozen_member_ids": True, "same_frozen_size": True,
                             "all_frozen_core_evidence_reproduced": True})
        print(f"Prepared and verified {target}: {len(groups)} groups", flush=True)
    rebuilt = {(common["target"], common["group"]) for common, _ in prepared}
    assert rebuilt == set(timer.index) == set(old.index) == set(comparison[keys].itertuples(index=False, name=None))
    return prepared, pd.DataFrame(evidence), protocols


def run(output=OUTPUT):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    frozen = CGM / "paper/output/reliability_v1/cores"
    input_paths = [frozen / name for name in ("group_core_comparison.csv", "core_rule_timing.csv", "fmcs_full_recomputation.csv", "core_timing_metadata.json")]
    input_paths += [CGM / f"paper/output/revision_v6/{target}/map_molecules.csv" for target in TARGETS]
    input_paths += [CGM / "paper/reliability_cores.py", CGM / "src/chemgridmap/cores.py", CGM / "src/chemgridmap/comparison.py"]
    comparison, timing, audit = (pd.read_csv(path) for path in input_paths[:3])
    frozen_metadata = json.loads(input_paths[3].read_text())
    table = stratify_cores(comparison)
    table.to_csv(output / "core_by_group_size.csv", index=False, na_rep="NA")
    prepared, membership, protocols = prepare_groups(comparison, timing, audit)
    membership.to_csv(output / "core_group_members.csv", index=False)
    cpu_model = platform.processor()
    physical_memory = None
    if platform.system() == "Darwin":
        hardware = json.loads(subprocess.check_output(["system_profiler", "SPHardwareDataType", "-json"], text=True))
        overview = hardware["SPHardwareDataType"][0]
        cpu_model = overview.get("chip_type") or overview.get("cpu_type") or cpu_model
        physical_memory = overview.get("physical_memory")
    metadata = {"started_at_utc": datetime.now(timezone.utc).isoformat(), "platform": platform.platform(),
                "machine": platform.machine(), "cpu_model": cpu_model, "logical_cpu_count": os.cpu_count(),
                "physical_memory": physical_memory,
                "python_version": platform.python_version(), "rdkit_version": rdBase.rdkitVersion,
                "process_id": os.getpid(), "workers": 1, "execution": "single process, sequential budgets and groups",
                "load_average_start": list(os.getloadavg()),
                "thread_environment": {key: os.environ.get(key) for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMBA_NUM_THREADS")},
                "timer": "mcs_core internal time.perf_counter: FindMCS, SMARTS parsing and member matching",
                "timing_excludes": "Group construction, SMILES parsing, validation, file I/O and export",
                "frozen_10s_metadata": frozen_metadata,
                "load_note": "Other E1-B compute jobs paused during this scan; unrelated system or background work was not controlled."}
    rows = []
    for timeout in (1, 2, 5):
        print(f"START timeout={timeout} load={os.getloadavg()}", flush=True)
        for index, (common, mols) in enumerate(prepared, 1):
            result = mcs_core(mols, min_coverage=.8, timeout=timeout)
            present = bool(result["core_smarts"] and result["core_heavy_atoms"])
            coverage = result["coverage_count"] / len(mols) if present else np.nan
            rows.append({**common, "timeout_seconds": timeout, "core_smarts": result["core_smarts"],
                         "core_heavy_atoms": result["core_heavy_atoms"] if present else np.nan,
                         "has_core": present, "core_status": "core" if present else "No core",
                         "coverage_count": result["coverage_count"], "coverage_fraction": coverage,
                         "target_met": bool(present and coverage >= .8 - 1e-12),
                         "timed_out": result["timed_out"], "seconds": result["seconds"],
                         "source": "E1B_sequential_recomputation"})
            if index % 25 == 0 or index == 480:
                print(f"FMCS timeout={timeout}: {index}/480 groups", flush=True)
        # 每个预算结束保存，方便明确识别已完成的原始计算。
        pd.DataFrame(rows).to_csv(output / "fmcs_budget_by_group.csv", index=False, na_rep="NA")
        print(f"END timeout={timeout} total_seconds={sum(row['seconds'] for row in rows if row['timeout_seconds'] == timeout):.6f}", flush=True)
    metadata["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
    metadata["load_average_finish"] = list(os.getloadavg())
    # 10 s 的覆盖和计时来自同一次冻结全量重算；不将历史缓存覆盖与新计时拼接。
    merged = timing.merge(audit, on=["target", "target_label", "group", "size"], validate="one_to_one")
    assert np.allclose(merged.fmcs_seconds, merged.recomputed_seconds, rtol=0, atol=1e-12)
    assert (merged.fmcs_timed_out == merged.recomputed_timed_out).all()
    assert (merged.fmcs_coverage_count == merged.recomputed_coverage_count).all()
    for row in merged.itertuples(index=False):
        present = bool(pd.notna(row.recomputed_smarts) and row.recomputed_smarts and row.recomputed_heavy_atoms)
        coverage = row.recomputed_coverage_fraction if present else np.nan
        rows.append({"target": row.target, "target_label": row.target_label, "group": row.group, "size": row.size,
                     "timeout_seconds": 10, "core_smarts": row.recomputed_smarts if present else "",
                     "core_heavy_atoms": row.recomputed_heavy_atoms if present else np.nan,
                     "has_core": present, "core_status": "core" if present else "No core",
                     "coverage_count": row.recomputed_coverage_count, "coverage_fraction": coverage,
                     "target_met": bool(present and coverage >= .8 - 1e-12), "timed_out": row.fmcs_timed_out,
                     "seconds": row.fmcs_seconds, "source": "frozen_10s_full_recomputation"})
    raw = pd.DataFrame(rows)
    assert len(raw) == 1920 and not raw.duplicated(["target", "group", "timeout_seconds"]).any()
    assert np.isfinite(raw.seconds).all() and raw.seconds.ge(0).all()
    raw.to_csv(output / "fmcs_budget_by_group.csv", index=False, na_rep="NA")
    summary = summarize_budgets(raw)
    summary.to_csv(output / "fmcs_budget_sweep.csv", index=False, na_rep="NA")
    completed = timing.loc[~timing.fmcs_timed_out]
    large = table.loc[table.target.eq("overall") & table.size_stratum.eq(">=10")]
    result = {"passed": True, "groups_total": len(prepared), "groups_per_target": EXPECTED,
              "frozen_core_rows_reproduced": len(comparison), "complete_member_identity_checks": len(membership),
              "large_group_strata": large.to_dict("records"),
              "frozen_10s_comparison_timing_coverage_counts_identical": bool((merged.cached_coverage_count == merged.fmcs_coverage_count).all()),
              "frozen_10s_completed": {"groups": len(completed), "excluded_timeout_groups": int(timing.fmcs_timed_out.sum()),
                                       "fmcs_total_seconds": float(completed.fmcs_seconds.sum()),
                                       "coverage_total_seconds_same_groups": float(completed.coverage_seconds.sum()),
                                       "fmcs_divided_by_coverage": float(completed.fmcs_seconds.sum() / completed.coverage_seconds.sum())},
              "overall_budget_sweep": summary.loc[summary.target.eq("overall")].to_dict("records")}
    parameters = {"grouping": protocols, "coverage_target": .8, "timeouts_seconds_recomputed": [1, 2, 5],
                  "timeout_seconds_frozen": 10, "workers": 1, "group_order": "chembl205, chembl204, chembl240; frozen group order",
                  "member_order": "molecule_identity_key ascending", "seed": None,
                  "complete_rings_only": True, "ring_matches_ring_only": True,
                  "atom_compare": "CompareElements", "bond_compare": "CompareOrder",
                  "size_strata": list(STRATA), "no_core_policy": "Categorical No core; missing coverage and size; excluded from means",
                  "frozen_10s_source": "fmcs_full_recomputation.csv recomputed fields paired and checked against core_rule_timing.csv",
                  "timing_comparison_limit": "Budgets 1,2,5 were measured in this run; 10 s is frozen and not contemporaneous."}
    provenance = {"inputs": [str(path.relative_to(ROOT)) for path in input_paths], "parameters": parameters,
                  "metadata": metadata, "summary": result}
    (output / "core_fmcs_provenance.json").write_text(json.dumps(provenance, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2, allow_nan=False), flush=True)
    return provenance


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    run(parser.parse_args().output)
