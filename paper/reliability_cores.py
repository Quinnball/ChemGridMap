"""Compare group-core rules and retain explicit matching evidence for T2b."""
from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import subprocess
import time

import numpy as np
import pandas as pd
from rdkit import Chem, rdBase
from rdkit.Chem.Scaffolds import MurckoScaffold

from chemgridmap.comparison import structural_groups
from chemgridmap.cores import mcs_core, select_core


ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "paper/output/revision_v6"
CACHE = ROOT.parent / "project1/JCIM_review_innovation_20261024/results/core_selection/group_core_comparison.csv"
HISTORY_ROOT = ROOT.parent / "project1/JCIM_review_innovation_20261024"
TARGETS = {"chembl205": "CA II", "chembl204": "Thrombin", "chembl240": "hERG"}
EXPECTED_GROUPS = {"chembl205": 125, "chembl204": 152, "chembl240": 203}
THRESHOLDS = (0.6, 0.7, 0.8, 0.9, 1.0)
SAMPLE_SEED = 20261024


def describe_core(core, mols, ids, *, is_smarts=False, coverage_target=0.8):
    """Measure the query; use a matched fragment only as its SMILES illustration."""
    pattern = (Chem.MolFromSmarts(core) if is_smarts else Chem.MolFromSmiles(core)) if core else None
    if core and pattern is None:
        raise ValueError("A core could not be parsed.")
    matches = [m.GetSubstructMatch(pattern, useChirality=not is_smarts) if pattern is not None else ()
               for m in mols]
    matched = [bool(match) for match in matches]
    count, source_id = int(sum(matched)), ""
    smiles = core if not is_smarts else ""
    if pattern is not None:
        Chem.GetSymmSSSR(pattern)
        rings, heavy = int(pattern.GetRingInfo().NumRings()), int(pattern.GetNumHeavyAtoms())
        if is_smarts and count:
            first = matched.index(True)
            mol, atom_map = Chem.Mol(mols[first]), matches[first]
            # Kekulize before cutting so partial aromatic fragments remain valid SMILES.
            Chem.Kekulize(mol, clearAromaticFlags=True)
            bonds = [mol.GetBondBetweenAtoms(atom_map[b.GetBeginAtomIdx()], atom_map[b.GetEndAtomIdx()]).GetIdx()
                     for b in pattern.GetBonds()]
            smiles = Chem.MolFragmentToSmiles(mol, atomsToUse=list(atom_map), bondsToUse=bonds,
                                             isomericSmiles=False, canonical=True)
            illustration = Chem.MolFromSmiles(smiles)
            if illustration is None:
                raise ValueError("An FMCS matched fragment is not a valid SMILES.")
            smiles = Chem.MolToSmiles(illustration)
            source_id = ids[first]
        benzene = bool(heavy == 6 and rings == 1 and smiles and
                       Chem.MolToSmiles(Chem.MolFromSmiles(smiles)) == "c1ccccc1")
    else:
        rings, heavy, benzene = 0, 0, False
    result = {"core_smiles": smiles, "core_smarts": core if is_smarts else "",
              "core_smiles_kind": "matched_member_fragment" if is_smarts else "exact_core",
              "core_smiles_source_id": source_id, "has_core": pattern is not None,
              "core_status": "core" if pattern is not None else "No core",
              "core_heavy_atoms": heavy, "core_ring_count": rings,
              "coverage_count": count, "coverage_fraction": count / len(mols) if pattern is not None else np.nan,
              "strict_majority": bool(pattern is not None and count > len(mols) / 2),
              "target_met": bool(pattern is not None and count / len(mols) >= coverage_target - 1e-12),
              "coverage_target": coverage_target, "is_benzene": benzene, "is_single_ring": rings == 1,
              "uncovered_member_ids": "|".join(key for key, yes in zip(ids, matched) if not yes)}
    return result, matched


def summarize_cores(table, by):
    """Summarize coverage and size; missing cores are not zero-coverage points."""
    rows = []
    for key, group in table.groupby(by, sort=False, dropna=False):
        key = key if isinstance(key, tuple) else (key,)
        present = group.loc[group.has_core]
        rows.append({**dict(zip(by, key)), "groups": len(group), "no_core": int((~group.has_core).sum()),
                     "strict_majority_groups": int(group.strict_majority.sum()),
                     "target_met_groups": int(group.target_met.sum()),
                     "mean_coverage": float(present.coverage_fraction.mean()) if len(present) else np.nan,
                     "median_coverage": float(present.coverage_fraction.median()) if len(present) else np.nan,
                     "median_core_heavy_atoms": float(present.core_heavy_atoms.median()) if len(present) else np.nan,
                     "median_core_ring_count": float(present.core_ring_count.median()) if len(present) else np.nan,
                     "benzene_cores": int(group.is_benzene.sum()),
                     "single_ring_cores": int(group.is_single_ring.sum())})
    return pd.DataFrame(rows)


def select_examples(table, limit=5):
    """Rank eligible improvements, preferring cores with at least nine heavy atoms."""
    index = ["target", "target_label", "group", "size"]
    fields = ["core_smiles", "core_heavy_atoms", "core_ring_count", "coverage_count", "coverage_fraction",
              "strict_majority", "target_met", "uncovered_member_ids"]
    current = table.loc[table.rule.eq("most_frequent"), index + fields].set_index(index).add_prefix("most_frequent_")
    aware = table.loc[table.rule.eq("coverage"), index + fields].set_index(index).add_prefix("coverage_")
    joined = current.join(aware, validate="one_to_one").reset_index()
    joined = joined.loc[~joined.most_frequent_strict_majority & joined.coverage_target_met].copy()
    joined["coverage_gain"] = joined.coverage_coverage_fraction - joined.most_frequent_coverage_fraction
    joined["preferred_core_at_least_9_heavy_atoms"] = joined.coverage_core_heavy_atoms.ge(9)
    joined = joined.sort_values(["preferred_core_at_least_9_heavy_atoms", "coverage_gain", "size", "target", "group"],
                                ascending=[False, False, False, True, True]).head(limit)
    joined.insert(0, "candidate_rank", range(1, len(joined) + 1))
    return joined


def compare_fmcs(cached, recalculated, mols, ids):
    """Keep strict SMARTS, graph-size and member-match comparisons separate."""
    old, old_matches = describe_core(str(cached.mcs_core), mols, ids, is_smarts=True)
    new, new_matches = describe_core(recalculated["core_smarts"], mols, ids, is_smarts=True)
    checks = {"same_smarts": str(cached.mcs_core) == recalculated["core_smarts"],
              "same_heavy_atoms": int(cached.mcs_heavy_atoms) == new["core_heavy_atoms"],
              "same_ring_count": old["core_ring_count"] == new["core_ring_count"],
              "same_member_matches": old_matches == new_matches,
              "same_coverage_count": old["coverage_count"] == new["coverage_count"],
              "same_coverage": bool(np.isclose(float(cached.mcs_coverage), new["coverage_fraction"], atol=1e-14, rtol=0)),
              "same_timeout_status": bool(cached.mcs_timed_out) == bool(recalculated["timed_out"])}
    return {"cached_smarts": str(cached.mcs_core), "recomputed_smarts": recalculated["core_smarts"],
            "cached_heavy_atoms": int(cached.mcs_heavy_atoms), "recomputed_heavy_atoms": new["core_heavy_atoms"],
            "cached_coverage_count": old["coverage_count"], "recomputed_coverage_count": new["coverage_count"],
            "cached_coverage_fraction": float(cached.mcs_coverage), "recomputed_coverage_fraction": new["coverage_fraction"],
            "cached_ring_count": old["core_ring_count"], "recomputed_ring_count": new["core_ring_count"],
            "cached_timed_out": bool(cached.mcs_timed_out), "recomputed_timed_out": bool(recalculated["timed_out"]),
            "cached_seconds": float(cached.mcs_seconds), "recomputed_seconds": float(recalculated["seconds"]),
            "cached_matched_member_ids": "|".join(key for key, yes in zip(ids, old_matches) if yes),
            "recomputed_matched_member_ids": "|".join(key for key, yes in zip(ids, new_matches) if yes),
            **checks, "passed": all(checks.values())}


def run_cores(output_root, *, input_root=FROZEN, cache_path=CACHE):
    """Write the 480-group comparison and return shared-manifest provenance."""
    output_root, input_root, cache_path = Path(output_root), Path(input_root), Path(cache_path)
    folder = output_root / "cores"
    folder.mkdir(parents=True, exist_ok=True)
    cache = pd.read_csv(cache_path).fillna("")
    # Normalize the frozen cache's legacy label in memory; never rewrite the cache.
    cache["target"] = cache.target.replace({"D3": "Thrombin"})
    cached = cache.set_index(["target", "group"])
    sample_indices = np.random.default_rng(SAMPLE_SEED).choice(len(cache), 20, replace=False)
    sampled = {(str(cache.iloc[i].target), str(cache.iloc[i].group)): rank
               for rank, i in enumerate(sample_indices, 1)}
    inputs = [cache_path]
    rows, sensitivity, members_out, audit, reproduction = [], [], [], [], []
    protocols, counts = {}, {}
    for target, label in TARGETS.items():
        path = input_root / target / "map_molecules.csv"
        inputs.append(path)
        data = pd.read_csv(path, low_memory=False).sort_values("molecule_identity_key").reset_index(drop=True)
        frame, groups, protocols[target] = structural_groups(data, similarity=0.6, core_rule="most_frequent")
        groups = [g for g in groups if g["size"] >= 2]
        counts[target] = len(groups)
        if len(groups) != EXPECTED_GROUPS[target]:
            raise ValueError(f"Unexpected nonsingleton group count for {target}: {len(groups)}.")
        for group in groups:
            group_id = group["structure_group"]
            members = frame.loc[frame.structure_group.eq(group_id)].copy()
            ids = members.molecule_identity_key.astype(str).tolist()
            mols = [Chem.MolFromSmiles(s) for s in members.canonical_smiles]
            common = {"target": target, "target_label": label, "group": group_id, "size": len(mols)}
            previous = cached.loc[(label, group_id)]
            rule_matches = {}
            current, rule_matches["most_frequent"] = describe_core(group["core_smiles"], mols, ids)
            rows.append({**common, "rule": "most_frequent", **current, "fmcs_timed_out": False})
            aware = None
            for threshold in THRESHOLDS:
                choice = select_core(mols, min_coverage=threshold)
                description, matched = describe_core(choice["core_smiles"], mols, ids, coverage_target=threshold)
                sensitivity.append({**common, **description, "candidate_count": int(choice["candidates"])})
                if threshold == 0.8:
                    aware = description
                    rule_matches["coverage"] = matched
            rows.append({**common, "rule": "coverage", **aware, "fmcs_timed_out": False})
            fmcs, rule_matches["fmcs"] = describe_core(str(previous.mcs_core), mols, ids, is_smarts=True)
            rows.append({**common, "rule": "fmcs", **fmcs, "fmcs_timed_out": bool(previous.mcs_timed_out)})
            checks = {"same_size": int(previous["size"]) == len(mols),
                      "same_most_frequent_core": previous.current_core == current["core_smiles"],
                      "same_coverage_core": previous.aware_core == aware["core_smiles"],
                      "same_fmcs_heavy_atoms": int(previous.mcs_heavy_atoms) == fmcs["core_heavy_atoms"]}
            for rule, old_field, description in (("most_frequent", "current_coverage", current),
                                                ("coverage", "aware_coverage", aware),
                                                ("fmcs", "mcs_coverage", fmcs)):
                old_coverage = previous[old_field]
                checks[f"same_{rule}_coverage"] = bool(old_coverage == "" and not description["has_core"] or
                    old_coverage != "" and np.isclose(float(old_coverage), description["coverage_fraction"], atol=1e-14, rtol=0))
            reproduction.append({**common, **checks, "passed": all(checks.values())})
            if (label, group_id) in sampled:
                result = mcs_core(mols, min_coverage=0.8, timeout=10)
                audit.append({"sample_rank": sampled[(label, group_id)], "sample_seed": SAMPLE_SEED,
                              **common, **compare_fmcs(previous, result, mols, ids)})
            if target == "chembl205" and group_id in {"G001", "G002"}:
                members["target"] = target
                members["target_label"] = label
                members["group"] = group_id
                for rule, matched in rule_matches.items():
                    members[f"matches_{rule}_core"] = matched
                members["matches_exact_scaffold"] = members.scaffold_smiles.eq(current["core_smiles"])
                members["most_frequent_core_smiles"] = current["core_smiles"]
                members["coverage_core_smiles"] = aware["core_smiles"]
                fields = ["target", "target_label", "group", "molecule_identity_key", "canonical_smiles",
                          "scaffold_smiles", "activity_pchembl", "activity_class", "matches_most_frequent_core",
                          "matches_coverage_core", "matches_fmcs_core", "matches_exact_scaffold",
                          "most_frequent_core_smiles", "coverage_core_smiles"]
                members_out.append(members[fields])
        print(f"Cores: {label}, {counts[target]} nonsingleton groups", flush=True)
    table, sens = pd.DataFrame(rows), pd.DataFrame(sensitivity)
    matches = pd.concat(members_out, ignore_index=True)
    fmcs_audit = pd.DataFrame(audit).sort_values("sample_rank")
    reproduction = pd.DataFrame(reproduction)
    tradeoff = pd.concat([summarize_cores(table, ["target", "target_label", "rule"]),
                         summarize_cores(table.assign(target="overall", target_label="Overall"),
                                         ["target", "target_label", "rule"])], ignore_index=True)
    sens_summary = summarize_cores(sens, ["target", "target_label", "coverage_target"])
    examples = select_examples(table)
    files = {"group_core_comparison.csv": table, "core_target_sensitivity.csv": sens,
             "core_target_sensitivity_summary.csv": sens_summary, "core_size_tradeoff.csv": tradeoff,
             "G001_G002_members.csv": matches, "fmcs_recomputation.csv": fmcs_audit,
             "working_copy_reproduction.csv": reproduction}
    for name, frame in files.items():
        frame.to_csv(folder / name, index=False)
    (output_root / "cases").mkdir(parents=True, exist_ok=True)
    examples.to_csv(output_root / "cases/core_rule_examples.csv", index=False)
    g002 = matches.loc[matches.group.eq("G002")]
    g001 = table.loc[table.target.eq("chembl205") & table.group.eq("G001")]
    case_check = {"size": len(g002), "coverage_core_smiles": sorted(g002.coverage_core_smiles.unique()),
                  "coverage_count": int(g002.matches_coverage_core.sum()),
                  "exact_scaffold_count": int(g002.matches_exact_scaffold.sum()),
                  "uncovered_ids": g002.loc[~g002.matches_coverage_core, "molecule_identity_key"].tolist()}
    case_check["passed"] = (case_check["size"] == 28 and case_check["coverage_count"] == 27 and
                            case_check["exact_scaffold_count"] == 17 and
                            case_check["coverage_core_smiles"] == ["c1ccc2scnc2c1"] and
                            case_check["uncovered_ids"] == ["CHEMBL85636"])
    result = {"nonsingleton_groups": counts, "groups_total": sum(counts.values()), "rule_rows": len(table),
              "working_copy_reproduction_passed": bool(reproduction.passed.all()),
              "fmcs_sample_seed": SAMPLE_SEED, "fmcs_recomputed_groups": len(fmcs_audit),
              "fmcs_recomputation_passed": bool(fmcs_audit.passed.all()),
              "fmcs_cached_timeout_groups": int(table.loc[table.rule.eq("fmcs"), "fmcs_timed_out"].sum()),
              "G002": case_check,
              "G001": g001[["rule", "core_smiles", "core_heavy_atoms", "coverage_count", "size"]].to_dict("records"),
              "candidate_groups": examples[["candidate_rank", "target", "group"]].to_dict("records"),
              "overall_tradeoff": tradeoff.loc[tradeoff.target.eq("overall")].to_dict("records")}
    result["passed"] = bool(result["working_copy_reproduction_passed"] and result["fmcs_recomputation_passed"]
                            and case_check["passed"] and len(fmcs_audit) == 20 and len(table) == 1440)
    (folder / "core_checks.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    if not result["passed"]:
        raise ValueError("Core reproduction failed; inspect cores/core_checks.json and the audit tables.")
    parameters = {"grouping": protocols, "coverage_targets": THRESHOLDS, "fmcs_threshold": 0.8,
                  "fmcs_complete_rings_only": True, "fmcs_ring_matches_ring_only": True,
                  "fmcs_timeout_seconds": 10, "fmcs_sample_seed": SAMPLE_SEED, "fmcs_sample_size": 20,
                  "fmcs_sample_order": "Existing cache row order; numpy default_rng choice without replacement.",
                  "fmcs_smiles": "Matched-member fragment for display; SMARTS is authoritative for matching.",
                  "no_core_policy": "No core is categorical; coverage is missing; size medians exclude missing cores.",
                  "example_order": "Core heavy atoms >= 9 first; coverage gain descending; size descending; target/group."}
    return {"inputs": inputs, "parameters": parameters, "summary": result}


def most_frequent_core(mols):
    """Time the software's exact-scaffold rule from prepared molecules, including matching."""
    started = time.perf_counter()
    scaffolds = [MurckoScaffold.MurckoScaffoldSmiles(mol=m, includeChirality=True) for m in mols]
    counts = Counter(s for s in scaffolds if s)
    ranked = sorted(counts, key=lambda s: (-counts[s], -Chem.MolFromSmiles(s).GetNumHeavyAtoms(), s))
    core = ranked[0] if ranked else ""
    pattern = Chem.MolFromSmiles(core) if core else None
    matched = [bool(pattern is not None and m.HasSubstructMatch(pattern, useChirality=True)) for m in mols]
    return {"core_smiles": core, "coverage_count": sum(matched),
            "coverage_fraction": sum(matched) / len(mols), "seconds": time.perf_counter() - started}


def summarize_core_timing(table):
    """Retain current and historical timer populations, without inventing a historical baseline."""
    rules = {"most_frequent": ("most_frequent_seconds", None),
             "coverage": ("coverage_seconds", None), "fmcs": ("fmcs_seconds", "fmcs_timed_out"),
             "historical_coverage": ("history_aware_seconds", None),
             "historical_fmcs": ("history_mcs_seconds", "history_mcs_timed_out")}
    rows = []
    populations = list(table.groupby(["target", "target_label"], sort=False)) + [(("overall", "Overall"), table)]
    for (target, label), group in populations:
        for rule, (seconds, timeout) in rules.items():
            values = group[seconds]
            rows.append({"target": target, "target_label": label, "rule": rule, "groups": len(group),
                         "total_seconds": float(values.sum()), "median_seconds": float(values.median()),
                         "p90_seconds": float(values.quantile(.9)), "max_seconds": float(values.max()),
                         "timeout_groups": int(group[timeout].sum()) if timeout else 0})
    return pd.DataFrame(rows)


def summarize_fmcs_recomputation(audit):
    """A timeout may change the partial query; report that separately from completed searches."""
    timeout = audit.cached_timed_out | audit.recomputed_timed_out
    fields = [name for name in audit if name.startswith("same_")]
    return {"groups": len(audit), "exact_matches": int(audit.passed.sum()),
            "different_groups": int((~audit.passed).sum()),
            "timeout_either_run_groups": int(timeout.sum()),
            "differences_with_timeout": int((~audit.passed & timeout).sum()),
            "differences_without_timeout": int((~audit.passed & ~timeout).sum()),
            "cached_timeout_groups": int(audit.cached_timed_out.sum()),
            "recomputed_timeout_groups": int(audit.recomputed_timed_out.sum()),
            "difference_counts": {name: int((~audit[name]).sum()) for name in fields},
            "non_timeout_reproduction_passed": bool(audit.loc[~timeout, "passed"].all())}


def core_timing_metadata():
    cpu_model = platform.processor()
    if platform.system() == "Darwin":
        hardware = json.loads(subprocess.check_output(["system_profiler", "SPHardwareDataType", "-json"], text=True))
        overview = hardware["SPHardwareDataType"][0]
        cpu_model = overview.get("chip_type") or overview.get("cpu_type") or cpu_model
    elif platform.system() == "Linux" and Path("/proc/cpuinfo").exists():
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                cpu_model = line.split(":", 1)[1].strip()
                break
    return {"started_at_utc": datetime.now(timezone.utc).isoformat(),
            "platform": platform.platform(), "cpu_model": cpu_model or platform.machine(),
            "machine": platform.machine(), "logical_cpu_count": os.cpu_count(),
            "python_version": platform.python_version(), "rdkit_version": rdBase.rdkitVersion,
            "process_id": os.getpid(), "workers": 1, "execution": "single process, sequential groups and rules",
            "thread_environment": {key: os.environ.get(key) for key in
                                   ("NUMBA_NUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")},
            "clock": "time.perf_counter; one measured invocation per rule and group",
            "excluded": "File I/O, SMILES parsing, fingerprint/group construction, audit, and export.",
            "included": {
                "most_frequent": "Murcko scaffold generation, frequency/tie ranking, query parsing, and member matching.",
                "coverage": "Software select_core: scaffold/network candidates, candidate parsing/matching, ranking and selection.",
                "fmcs": "Software mcs_core: FindMCS, returned SMARTS parsing, and member matching."},
            "audit_seconds": "Independent cached/current SMARTS, graph-size and member-match audit; excluded from rule times.",
            "preparation": "Groups and member molecules are prepared before timing; each rule recomputes its own core.",
            "rule_order": ["most_frequent", "coverage", "fmcs"],
            "member_order": "molecule_identity_key ascending, matching the historical cache input order",
            "historical_timers": {
                "aware_seconds": "Same select_core internal timer boundary as current coverage_seconds.",
                "mcs_seconds": "Same mcs_core internal timer boundary as current fmcs_seconds.",
                "most_frequent": "Not recorded; no historical time is inferred.",
                "machine_and_load": "Historical machine and system load were not recorded; timings are contextual comparisons."}}


def run_cores_timing(output_root, *, input_root=FROZEN, cache_path=CACHE):
    """Time all 480 groups sequentially and audit full FMCS recomputation without replacing old outputs."""
    output_root, input_root, cache_path = Path(output_root), Path(input_root), Path(cache_path)
    folder = output_root / "cores"
    folder.mkdir(parents=True, exist_ok=True)
    cache = pd.read_csv(cache_path).fillna("")
    # Normalize the frozen cache's legacy label in memory; never rewrite the cache.
    cache["target"] = cache.target.replace({"D3": "Thrombin"})
    cached = cache.set_index(["target", "group"])
    if not cached.index.is_unique:
        raise ValueError("The historical FMCS cache contains duplicate target/group keys.")
    metadata = core_timing_metadata()
    inputs = [cache_path, HISTORY_ROOT / "software/chemgridmap/src/chemgridmap/cores.py",
              HISTORY_ROOT / "analysis/evaluate_core_selection.py"]
    rows, audit_rows, protocols, counts = [], [], {}, {}
    for target, label in TARGETS.items():
        path = input_root / target / "map_molecules.csv"
        inputs.append(path)
        data = pd.read_csv(path, low_memory=False).sort_values("molecule_identity_key").reset_index(drop=True)
        frame, groups, protocols[target] = structural_groups(data, similarity=.6, core_rule="most_frequent")
        groups = [g for g in groups if g["size"] >= 2]
        counts[target] = len(groups)
        if len(groups) != EXPECTED_GROUPS[target]:
            raise ValueError(f"Unexpected nonsingleton group count for {target}: {len(groups)}.")
        for index, group in enumerate(groups, 1):
            group_id = group["structure_group"]
            members = frame.loc[frame.structure_group.eq(group_id)]
            ids = members.molecule_identity_key.astype(str).tolist()
            mols = [Chem.MolFromSmiles(s) for s in members.canonical_smiles]
            previous = cached.loc[(label, group_id)]
            common = {"target": target, "target_label": label, "group": group_id, "size": len(mols)}
            current = most_frequent_core(mols)
            aware = select_core(mols, min_coverage=.8)
            fmcs = mcs_core(mols, min_coverage=.8, timeout=10)
            started = time.perf_counter()
            audit = compare_fmcs(previous, fmcs, mols, ids)
            audit_seconds = time.perf_counter() - started
            audit_rows.append({**common, **audit})
            rows.append({**common, "most_frequent_seconds": current["seconds"],
                         "coverage_seconds": aware["seconds"], "fmcs_seconds": fmcs["seconds"],
                         "fmcs_timed_out": fmcs["timed_out"], "audit_seconds": audit_seconds,
                         "history_aware_seconds": float(previous.aware_seconds),
                         "history_mcs_seconds": float(previous.mcs_seconds),
                         "history_mcs_timed_out": bool(previous.mcs_timed_out),
                         "most_frequent_core": current["core_smiles"], "coverage_core": aware["core_smiles"],
                         "most_frequent_coverage_count": current["coverage_count"],
                         "coverage_coverage_count": aware["coverage_count"], "fmcs_coverage_count": fmcs["coverage_count"],
                         "same_cached_size": int(previous["size"]) == len(mols),
                         "same_cached_most_frequent_core": previous.current_core == current["core_smiles"],
                         "same_cached_coverage_core": previous.aware_core == aware["core_smiles"]})
            if index % 25 == 0 or index == len(groups):
                print(f"Core timing: {label}, {index}/{len(groups)} groups; FMCS cache audit recorded", flush=True)
    table, audit = pd.DataFrame(rows), pd.DataFrame(audit_rows)
    summary = summarize_core_timing(table)
    fmcs_checks = summarize_fmcs_recomputation(audit)
    seconds = table[[name for name in table if name.endswith("_seconds")]]
    expected_keys = set(cached.index)
    actual_keys = set(zip(table.target_label, table.group))
    checks = {"groups_total": len(table), "nonsingleton_groups": counts,
              "all_cached_groups_recomputed": expected_keys == actual_keys,
              "unique_target_group": not table.duplicated(["target", "group"]).any(),
              "finite_nonnegative_seconds": bool(np.isfinite(seconds).all().all() and seconds.ge(0).all().all()),
              "same_cached_size": bool(table.same_cached_size.all()),
              "same_cached_most_frequent_core": bool(table.same_cached_most_frequent_core.all()),
              "same_cached_coverage_core": bool(table.same_cached_coverage_core.all()),
              "fmcs_recomputation": fmcs_checks,
              "non_timeout_fmcs_reproduction_passed": fmcs_checks["non_timeout_reproduction_passed"],
              "overall_timing": summary.loc[summary.target.eq("overall")].to_dict("records")}
    checks["passed"] = bool(len(table) == sum(EXPECTED_GROUPS.values()) and all(checks[name] for name in
        ("all_cached_groups_recomputed", "unique_target_group", "finite_nonnegative_seconds", "same_cached_size",
         "same_cached_most_frequent_core", "same_cached_coverage_core", "non_timeout_fmcs_reproduction_passed")))
    metadata["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
    parameters = {"grouping": protocols, "coverage_target": .8, "fmcs_threshold": .8,
                  "fmcs_complete_rings_only": True, "fmcs_ring_matches_ring_only": True,
                  "fmcs_atom_compare": "CompareElements", "fmcs_bond_compare": "CompareOrder",
                  "fmcs_timeout_seconds": 10, "workers": 1, "full_recomputation_groups": len(table),
                  "timeout_difference_policy": "Report differences whenever either run timed out; do not fail them.",
                  "metadata": "core_timing_metadata.json"}
    metadata["parameters"] = parameters
    for name, values in {"core_rule_timing.csv": table, "fmcs_full_recomputation.csv": audit,
                         "core_rule_timing_summary.csv": summary}.items():
        values.to_csv(folder / name, index=False)
    for name, values in {"core_timing_metadata.json": metadata, "core_timing_checks.json": checks}.items():
        (folder / name).write_text(json.dumps(values, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    if not checks["passed"]:
        raise ValueError("Core timing validation failed; inspect cores/core_timing_checks.json and full recomputation.")
    return {"inputs": inputs, "parameters": parameters, "summary": checks}
