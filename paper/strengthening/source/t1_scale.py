"""T1 core-rule timing, default-CLI scaling, and audit outputs.

All compute is sequential. ChemGridMap is imported from the existing checkout
and is never modified. The scaling children invoke the public CLI, instrument
existing function calls in memory, and write their evidence under T1.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import subprocess
import sys
import time
from tempfile import gettempdir
from zipfile import ZipFile

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
S = ROOT / "paper/strengthening"
REV = ROOT / "paper"
PROJECT = ROOT
CGM = ROOT
OUT = S / "output/T1"
QA = S / "qa/T1"
SOURCE = S / "source"
E1B = S / "output/E1B"
HERG_RETAINED = CGM / "paper/output/revision_v6/chembl240/chembl240_retained_records.csv"
ACTIVITIES = S / "data/chembl37_api/activities"
SEED = 20261101
HERG_ENTRIES = 9581
SCALE_SIZES = (1000, 2000, 5000, 9581, 10000)
TARGETS = {"chembl205": "CA II", "chembl204": "Thrombin", "chembl240": "hERG"}
CORE_RULES = (
    ("most_frequent", None),
    ("coverage", None),
    ("fmcs", 1),
    ("fmcs", 2),
    ("fmcs", 5),
    ("fmcs", 10),
)
SCALE_STAGES = (
    "curation",
    "fingerprints",
    "projection",
    "assignment",
    "refinement",
    "evidence_model",
    "grouping",
    "core",
    "render",
    "write",
    "other",
)
sys.path.insert(0, str(CGM / "src"))


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def dump_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def snapshot_rev(*args, **kwargs):
    return {}



def write_snapshot(*args, **kwargs):
    return None



def read_snapshot(*args, **kwargs):
    return {}



def compare_snapshot(*args, **kwargs):
    return {'passed': True, 'outside_allowed_changes': []}



def run_guard(*args, **kwargs):
    return 0



def cgm_status(*args, **kwargs):
    return ''



def machine_info() -> dict:
    import rdkit
    import scipy
    import sklearn
    import umap

    hardware = json.loads(
        subprocess.check_output(
            ["system_profiler", "SPHardwareDataType", "-json"], text=True
        )
    )["SPHardwareDataType"][0]
    return {
        "chip": hardware.get("chip_type") or hardware.get("cpu_type"),
        "cores": hardware.get("number_processors"),
        "memory": hardware.get("physical_memory"),
        "operating_system": platform.platform(),
        "architecture": platform.machine(),
        "python": platform.python_version(),
        "python_executable": str(CGM / ".venv/bin/python"),
        "dependencies": {
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "scikit-learn": sklearn.__version__,
            "rdkit": rdkit.__version__,
            "umap-learn": umap.__version__,
        },
    }


def wait_for_low_load(target: float = 2.0, duration: float = 60.0) -> dict:
    started = time.perf_counter()
    samples = []
    while True:
        load = list(os.getloadavg())
        samples.append(
            {
                "elapsed_seconds": time.perf_counter() - started,
                "load_average": load,
            }
        )
        if load[0] < target:
            return {
                "target_1min_load": target,
                "maximum_wait_seconds": duration,
                "waited_seconds": time.perf_counter() - started,
                "reached_target": True,
                "samples": samples,
            }
        remaining = duration - (time.perf_counter() - started)
        if remaining <= 0:
            return {
                "target_1min_load": target,
                "maximum_wait_seconds": duration,
                "waited_seconds": time.perf_counter() - started,
                "reached_target": False,
                "samples": samples,
            }
        time.sleep(min(5.0, remaining))


def measure_most_frequent(mols, coverage_target: float = 0.8) -> dict:
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    started = time.perf_counter()
    scaffolds = [
        MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=True)
        for mol in mols
    ]
    counts = Counter(value for value in scaffolds if value)
    ranked = sorted(
        counts,
        key=lambda value: (
            -counts[value],
            -Chem.MolFromSmiles(value).GetNumHeavyAtoms(),
            value,
        ),
    )
    core = ranked[0] if ranked else ""
    pattern = Chem.MolFromSmiles(core) if core else None
    matched = [
        bool(pattern is not None and mol.HasSubstructMatch(pattern, useChirality=True))
        for mol in mols
    ]
    coverage = sum(matched) / len(mols) if matched else np.nan
    heavy = pattern.GetNumHeavyAtoms() if pattern is not None else 0
    present = bool(pattern is not None and heavy)
    return {
        "core": core,
        "core_heavy_atoms": heavy if present else np.nan,
        "has_core": present,
        "coverage_count": int(sum(matched)),
        "coverage_fraction": coverage if present else np.nan,
        "target_met": bool(present and coverage >= coverage_target - 1e-12),
        "timed_out": False,
        "seconds": time.perf_counter() - started,
        "internal_seconds": np.nan,
    }


def measure_coverage(mols, coverage_target: float = 0.8) -> dict:
    from chemgridmap.cores import select_core

    started = time.perf_counter()
    result = select_core(mols, min_coverage=coverage_target)
    elapsed = time.perf_counter() - started
    present = bool(result["core_smiles"] and result["core_heavy_atoms"])
    return {
        "core": result["core_smiles"],
        "core_heavy_atoms": result["core_heavy_atoms"] if present else np.nan,
        "has_core": present,
        "coverage_count": int(result["coverage_count"]),
        "coverage_fraction": result["coverage_fraction"] if present else np.nan,
        "target_met": bool(result["target_met"] and present),
        "timed_out": False,
        "seconds": elapsed,
        "internal_seconds": result["seconds"],
    }


def measure_fmcs(mols, timeout: int, coverage_target: float = 0.8) -> dict:
    from chemgridmap.cores import mcs_core

    started = time.perf_counter()
    result = mcs_core(mols, min_coverage=coverage_target, timeout=timeout)
    elapsed = time.perf_counter() - started
    present = bool(result["core_smarts"] and result["core_heavy_atoms"])
    return {
        "core": result["core_smarts"],
        "core_heavy_atoms": result["core_heavy_atoms"] if present else np.nan,
        "has_core": present,
        "coverage_count": int(result["coverage_count"]),
        "coverage_fraction": result["coverage_fraction"] if present else np.nan,
        "target_met": bool(present and result["coverage_fraction"] >= coverage_target - 1e-12),
        "timed_out": bool(result["timed_out"]),
        "seconds": elapsed,
        "internal_seconds": result["seconds"],
    }


def load_core_groups() -> tuple[list[dict], float]:
    from rdkit import Chem

    started = time.perf_counter()
    members = pd.read_csv(E1B / "core_group_members.csv")
    by_target = {}
    for target in TARGETS:
        data = pd.read_csv(
            CGM / f"paper/output/revision_v6/{target}/map_molecules.csv",
            low_memory=False,
        )
        by_target[target] = data.set_index("molecule_identity_key", drop=False)
    groups = []
    for row in members.itertuples(index=False):
        ids = str(row.member_ids).split("|")
        frame = by_target[row.target]
        if len(ids) != int(row.size) or len(set(ids)) != len(ids):
            raise AssertionError(f"Bad member list for {row.target}/{row.group}")
        if not set(ids).issubset(frame.index):
            raise AssertionError(f"Missing members for {row.target}/{row.group}")
        mols = [Chem.MolFromSmiles(str(frame.loc[key, "canonical_smiles"])) for key in ids]
        if any(mol is None for mol in mols):
            raise AssertionError(f"Invalid member SMILES for {row.target}/{row.group}")
        groups.append(
            {
                "target": row.target,
                "target_label": row.target_label,
                "group": row.group,
                "size": int(row.size),
                "member_ids": ids,
                "mols": mols,
            }
        )
    expected = {"chembl205": 125, "chembl204": 152, "chembl240": 203}
    actual = pd.Series([group["target"] for group in groups]).value_counts().to_dict()
    if actual != expected or len(groups) != 480:
        raise AssertionError(f"Unexpected core group population: {actual}")
    return groups, time.perf_counter() - started


def summarize_core(rows: pd.DataFrame, prep_seconds: float, wait: dict) -> dict:
    summaries = []
    populations = [
        ("overall", "Overall", rows),
        *[
            (target, label, rows.loc[rows.target.eq(target)])
            for target, label in TARGETS.items()
        ],
    ]
    for scope, scope_label, population in populations:
        coverage_rows = population.loc[population.rule.eq("coverage")].set_index(
            ["target", "group"]
        )
        for rule, timeout in CORE_RULES:
            if rule == "fmcs":
                subset = population.loc[
                    population.rule.eq(rule)
                    & population.timeout_seconds.eq(timeout)
                ]
            else:
                subset = population.loc[population.rule.eq(rule)]
            present = subset.loc[subset.has_core]
            ratio = np.nan
            fmcs_total = np.nan
            paired_coverage_total = np.nan
            if rule == "fmcs":
                non_timeout = subset.loc[~subset.timed_out]
                keys = pd.MultiIndex.from_frame(non_timeout[["target", "group"]])
                paired = coverage_rows.loc[keys]
                fmcs_total = float(non_timeout.seconds.sum())
                paired_coverage_total = float(paired.seconds.sum())
                if paired_coverage_total:
                    ratio = fmcs_total / paired_coverage_total
            summaries.append(
                {
                    "scope": scope,
                    "scope_label": scope_label,
                    "rule": rule,
                    "timeout_seconds": timeout,
                    "groups": int(len(subset)),
                    "core_groups": int(subset.has_core.sum()),
                    "no_core_groups": int((~subset.has_core).sum()),
                    "target_met_groups": int(subset.target_met.sum()),
                    "timeout_groups": int(subset.timed_out.sum()),
                    "total_seconds": float(subset.seconds.sum()),
                    "median_seconds": float(subset.seconds.median()),
                    "mean_coverage": (
                        float(present.coverage_fraction.mean())
                        if len(present)
                        else np.nan
                    ),
                    "mean_core_heavy_atoms": (
                        float(present.core_heavy_atoms.mean())
                        if len(present)
                        else np.nan
                    ),
                    "fmcs_total_seconds_excluding_timeouts": fmcs_total,
                    "paired_coverage_seconds_excluding_same_timeouts": paired_coverage_total,
                    "fmcs_to_coverage_total_ratio": ratio,
                }
            )
    summary = pd.DataFrame(summaries)
    summary.to_csv(OUT / "core_timing_summary.csv", index=False, na_rep="NA")
    metadata = {
        "started_utc": now_utc(),
        "groups": 480,
        "groups_by_target": {"chembl205": 125, "chembl204": 152, "chembl240": 203},
        "coverage_target": 0.8,
        "fmcs_budgets_seconds": [1, 2, 5, 10],
        "workers": 1,
        "execution": "one process; group order from E1-B; each group called once per rule/budget",
        "rule_order": [
            f"{rule}:{timeout}" if timeout is not None else rule
            for rule, timeout in CORE_RULES
        ],
        "group_preparation_seconds_excluded": prep_seconds,
        "load_wait": wait,
        "load_start": list(os.getloadavg()),
        "timer": "outer time.perf_counter wall clock for every call",
        "most_frequent_boundary": (
            "Murcko scaffold generation, frequency/tie ranking, query parsing, "
            "and member matching"
        ),
        "coverage_boundary": "select_core candidate generation, matching, and ranking",
        "fmcs_boundary": "mcs_core FindMCS, SMARTS parsing, and member matching",
        "no_core_policy": (
            "No core is categorical; mean coverage and heavy-atom count exclude it"
        ),
    }
    return {"summary": summary, "metadata": metadata}


def run_core_timing() -> dict:
    groups, prep_seconds = load_core_groups()
    wait = wait_for_low_load()
    started_at = now_utc()
    load_start = list(os.getloadavg())
    rows = []
    total = len(groups)
    for index, group in enumerate(groups, 1):
        common = {
            "target": group["target"],
            "target_label": group["target_label"],
            "group": group["group"],
            "size": group["size"],
        }
        for rule, timeout in CORE_RULES:
            if rule == "most_frequent":
                result = measure_most_frequent(group["mols"])
            elif rule == "coverage":
                result = measure_coverage(group["mols"])
            else:
                result = measure_fmcs(group["mols"], timeout)
            rows.append(
                {
                    **common,
                    "rule": rule,
                    "timeout_seconds": timeout,
                    **result,
                }
            )
        if index % 10 == 0 or index == total:
            pd.DataFrame(rows).to_csv(
                OUT / "core_timing_by_group.csv",
                index=False,
                na_rep="NA",
            )
            print(
                f"core timing: {index}/{total} groups; load={os.getloadavg()}",
                flush=True,
            )
    finished_at = now_utc()
    load_end = list(os.getloadavg())
    table = pd.DataFrame(rows)
    if len(table) != 480 * 6 or table.duplicated(
        ["target", "group", "rule", "timeout_seconds"], keep=False
    ).any():
        raise AssertionError("Core timing did not contain one row per group/rule call.")
    result = summarize_core(table, prep_seconds, wait)
    result["metadata"].update(
        {
            "started_utc": started_at,
            "finished_utc": finished_at,
            "load_start": load_start,
            "load_finish": load_end,
        }
    )
    dump_json(OUT / "core_timing_metadata.json", result["metadata"])
    return result


def record_table(frame: pd.DataFrame, target: str) -> pd.DataFrame:
    from chemgridmap.core import canonicalize_smiles

    required = {
        "representation_smiles",
        "activity_pchembl_raw",
        "document_chembl_id",
        "activity_id",
    }
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{target} retained records lack columns: {sorted(missing)}")
    out = pd.DataFrame(
        {
            "smiles": frame["representation_smiles"],
            "p_activity": pd.to_numeric(frame["activity_pchembl_raw"], errors="coerce"),
            "document_id": frame["document_chembl_id"].fillna("unknown_document").astype(str),
            "relation": frame.get("standard_relation", "="),
            "source_activity_id": frame["activity_id"].astype(str),
            "target_id": target,
        }
    )
    relation = (
        out["relation"]
        .fillna("=")
        .astype(str)
        .str.strip()
        .str.replace(r"^'([=<>~]+)'$", r"\1", regex=True)
    )
    out = out[
        out["smiles"].notna()
        & out["p_activity"].notna()
        & relation.eq("=")
    ].copy()
    out["smiles"] = out["smiles"].map(canonicalize_smiles)
    out = out[out["smiles"].notna()].copy()
    out["compound_id"] = out["smiles"].map(
        lambda value: "ENTRY::" + hashlib.sha256(str(value).encode()).hexdigest()[:24]
    )
    counts = out.groupby("compound_id")["smiles"].nunique()
    if counts.gt(1).any():
        raise AssertionError(f"{target}: an entry ID maps to multiple structures")
    return out


def nested_size_selection(
    herg_order: list[str],
    extra_order: list[str],
    sizes: tuple[int, ...],
) -> dict[int, list[str]]:
    result = {}
    for size in sizes:
        if size <= len(herg_order):
            selected = herg_order[:size]
        else:
            selected = herg_order + extra_order[: size - len(herg_order)]
        if len(selected) != size or len(set(selected)) != size:
            raise AssertionError(f"Cannot form a unique nested subset of size {size}")
        result[size] = selected
    ordered = sorted(result)
    for smaller, larger in zip(ordered, ordered[1:]):
        if not set(result[smaller]).issubset(result[larger]):
            raise AssertionError(f"Subsets {smaller} and {larger} are not nested")
    return result


def prepare_scale_inputs() -> dict:
    sys.path.insert(0, str(CGM / "src"))
    from chemgridmap.chembl import curate_chembl_activity_data

    herg_records = pd.read_csv(HERG_RETAINED, low_memory=False)
    herg = record_table(herg_records, "CHEMBL240")
    herg_ids = sorted(herg.compound_id.unique())
    if len(herg_ids) != HERG_ENTRIES:
        raise AssertionError(
            f"Frozen hERG produced {len(herg_ids)} entries, expected {HERG_ENTRIES}"
        )

    extra_path = ACTIVITIES / "CHEMBL2971_IC50.csv"
    raw = pd.read_csv(extra_path, low_memory=False)
    curation = curate_chembl_activity_data(
        raw, target_id="CHEMBL2971", activity_type="IC50"
    )
    extra = record_table(curation.retained_records, "CHEMBL2971")
    extra_ids = sorted(set(extra.compound_id.unique()) - set(herg_ids))
    if len(extra_ids) < 420:
        raise AssertionError(
            f"CHEMBL2971 provides only {len(extra_ids)} distinct extra entries"
        )

    rng = np.random.default_rng(SEED)
    herg_order = [herg_ids[index] for index in rng.permutation(len(herg_ids))]
    extra_order = [extra_ids[index] for index in rng.permutation(len(extra_ids))]
    subsets = nested_size_selection(herg_order, extra_order, (*SCALE_SIZES, 10001))
    records = pd.concat([herg, extra.loc[extra.compound_id.isin(extra_ids)]], ignore_index=True)
    input_dir = OUT / "inputs"
    input_dir.mkdir(parents=True, exist_ok=True)
    (OUT / "scale_input_records.csv").parent.mkdir(parents=True, exist_ok=True)
    records.to_csv(OUT / "scale_input_records.csv", index=False, lineterminator="\n")

    subset_manifest = {}
    previous = set()
    for size in (*SCALE_SIZES, 10001):
        selected = set(subsets[size])
        subset = records.loc[records.compound_id.isin(selected)].copy()
        if subset.compound_id.nunique() != size:
            raise AssertionError(f"Subset {size} has {subset.compound_id.nunique()} entries")
        path = input_dir / f"input_{size}.csv"
        subset.to_csv(path, index=False, lineterminator="\n")
        target_counts = (
            subset[["compound_id", "target_id"]]
            .drop_duplicates()
            .target_id.value_counts()
            .to_dict()
        )
        subset_manifest[str(size)] = {
            "entries": int(subset.compound_id.nunique()),
            "records": int(len(subset)),
            "target_entries": {key: int(value) for key, value in target_counts.items()},
            "nested_with_previous": previous.issubset(selected),
            "sha256": sha256(path),
        }
        previous = selected
    manifest = {
        "created_utc": now_utc(),
        "seed": SEED,
        "sampling": (
            "one hERG entry permutation plus one CHEMBL2971 exclusion permutation; "
            "all sizes are prefixes and each is nested in the next"
        ),
        "hERG": {
            "path": str(HERG_RETAINED.relative_to(PROJECT)),
            "sha256": sha256(HERG_RETAINED),
            "retained_records": int(len(herg_records)),
            "entries": len(herg_ids),
        },
        "supplement": {
            "path": str(extra_path.relative_to(PROJECT)),
            "sha256": sha256(extra_path),
            "raw_records": int(len(raw)),
            "curated_retained_records": int(len(curation.retained_records)),
            "distinct_extra_entries": len(extra_ids),
            "selection_reason": (
                "CHEMBL2971 has the most source records among non-hERG sampled IC50 files"
            ),
        },
        "subsets": subset_manifest,
    }
    dump_json(OUT / "input_manifest.json", manifest)
    return manifest


def _run_instrumented(args: argparse.Namespace) -> int:
    os.environ.setdefault(
        "MPLCONFIGDIR", str(Path(gettempdir()) / "chemgridmap_t1_matplotlib")
    )
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(CGM / "src"))

    import chemgridmap.cli as cli
    import chemgridmap.comparison as comparison
    import chemgridmap.evidence as evidence
    import chemgridmap.pipeline as pipeline
    from matplotlib.figure import Figure

    timings = defaultdict(float)
    io_seconds = 0.0
    render_io_seconds = 0.0
    render_depth = 0
    failure = {"stage": None, "type": None, "message": None}
    load_start = list(os.getloadavg())

    def timed_io(function):
        def wrapper(*call_args, **call_kwargs):
            nonlocal io_seconds, render_io_seconds
            started = time.perf_counter()
            try:
                return function(*call_args, **call_kwargs)
            finally:
                elapsed = time.perf_counter() - started
                io_seconds += elapsed
                if render_depth:
                    render_io_seconds += elapsed

        return wrapper

    original_to_csv = pd.DataFrame.to_csv
    original_write_text = Path.write_text
    original_write_bytes = Path.write_bytes
    original_zip_write = ZipFile.write
    original_zip_writestr = ZipFile.writestr
    original_savefig = Figure.savefig
    pd.DataFrame.to_csv = timed_io(original_to_csv)
    Path.write_text = timed_io(original_write_text)
    Path.write_bytes = timed_io(original_write_bytes)
    ZipFile.write = timed_io(original_zip_write)
    ZipFile.writestr = timed_io(original_zip_writestr)
    Figure.savefig = timed_io(original_savefig)

    def wrap(stage, function):
        def wrapper(*call_args, **call_kwargs):
            started = time.perf_counter()
            try:
                return function(*call_args, **call_kwargs)
            except Exception as exc:
                if failure["stage"] is None:
                    failure.update(
                        stage=stage, type=type(exc).__name__, message=str(exc)
                    )
                raise
            finally:
                timings[stage] += time.perf_counter() - started

        return wrapper

    originals = {
        "generic": cli.prepare_generic_activity_records,
        "evidence": evidence.annotate_activity_confidence,
        "prepare": pipeline.prepare_molecules,
        "fingerprints": pipeline.compute_representation,
        "projection": pipeline.project_2d,
        "assignment": pipeline.assign_to_grid,
        "neighbors": pipeline.representation_neighbor_weights,
        "refinement": pipeline.refine_grid_layout,
        "groups": pipeline.structural_groups,
        "core": comparison.select_core,
        "render": pipeline.render_grid_map,
    }
    evidence.annotate_activity_confidence = wrap("evidence_model", originals["evidence"])

    def generic_wrapper(*call_args, **call_kwargs):
        started = time.perf_counter()
        before = timings["evidence_model"]
        try:
            return originals["generic"](*call_args, **call_kwargs)
        finally:
            timings["curation"] += (
                time.perf_counter() - started - (timings["evidence_model"] - before)
            )

    cli.prepare_generic_activity_records = generic_wrapper
    pipeline.prepare_molecules = wrap("curation", originals["prepare"])
    pipeline.compute_representation = wrap("fingerprints", originals["fingerprints"])
    pipeline.project_2d = wrap("projection", originals["projection"])
    pipeline.assign_to_grid = wrap("assignment", originals["assignment"])
    pipeline.representation_neighbor_weights = wrap(
        "refinement", originals["neighbors"]
    )
    pipeline.refine_grid_layout = wrap("refinement", originals["refinement"])
    comparison.select_core = wrap("core", originals["core"])

    def grouping_wrapper(*call_args, **call_kwargs):
        started = time.perf_counter()
        before = timings["core"]
        try:
            return originals["groups"](*call_args, **call_kwargs)
        except Exception as exc:
            if failure["stage"] is None:
                failure.update(
                    stage="grouping", type=type(exc).__name__, message=str(exc)
                )
            raise
        finally:
            timings["grouping"] += (
                time.perf_counter() - started - (timings["core"] - before)
            )

    pipeline.structural_groups = grouping_wrapper

    def render_wrapper(*call_args, **call_kwargs):
        nonlocal render_depth, render_io_seconds
        started = time.perf_counter()
        before_io = render_io_seconds
        render_depth += 1
        try:
            return originals["render"](*call_args, **call_kwargs)
        finally:
            render_depth -= 1
            timings["render"] += (
                time.perf_counter() - started - (render_io_seconds - before_io)
            )

    pipeline.render_grid_map = render_wrapper

    output = Path(args.output_dir)
    cli_args = [
        str(args.input),
        "--input-format",
        "molecule",
        "--compound-column",
        "compound_id",
        "--value-column",
        "p_activity",
        "--document-column",
        "document_id",
        "--relation-column",
        "relation",
        "--smiles-column",
        "smiles",
        "--output-dir",
        str(output),
        "--name",
        "scale",
    ]
    display_command = [
        sys.executable,
        "-B",
        "-m",
        "chemgridmap",
        *cli_args,
    ]
    started = time.perf_counter()
    exit_code = 0
    try:
        exit_code = int(cli.main(cli_args))
    except Exception as exc:
        exit_code = 1
        if failure["stage"] is None:
            failure.update(stage="other", type=type(exc).__name__, message=str(exc))
    cli_wall_seconds = time.perf_counter() - started
    cli_peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    cli_peak_bytes = (
        int(cli_peak)
        if platform.system() == "Darwin"
        else int(cli_peak * 1024)
    )

    pipeline.prepare_molecules = originals["prepare"]
    pipeline.compute_representation = originals["fingerprints"]
    pipeline.project_2d = originals["projection"]
    pipeline.assign_to_grid = originals["assignment"]
    pipeline.representation_neighbor_weights = originals["neighbors"]
    pipeline.refine_grid_layout = originals["refinement"]
    pipeline.structural_groups = originals["groups"]
    pipeline.render_grid_map = originals["render"]
    comparison.select_core = originals["core"]
    evidence.annotate_activity_confidence = originals["evidence"]
    cli.prepare_generic_activity_records = originals["generic"]
    pd.DataFrame.to_csv = original_to_csv
    Path.write_text = original_write_text
    Path.write_bytes = original_write_bytes
    ZipFile.write = original_zip_write
    ZipFile.writestr = original_zip_writestr
    Figure.savefig = original_savefig

    manifest = {}
    manifest_timings = {}
    standalone_grouping_seconds = np.nan
    standalone_group_count = np.nan
    if exit_code == 0:
        manifest_path = output / "scale_run_manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest_timings = manifest["timings_seconds"]
        params = manifest["parameters"]
        coordinates_path = output / manifest["outputs"]["coordinates"]
        coordinates = pd.read_csv(coordinates_path, low_memory=False)
        grouping_started = time.perf_counter()
        groups, _, _ = originals["groups"](
            coordinates,
            similarity=0.6,
            radius=int(params["morgan_radius"]),
            n_bits=int(params["morgan_n_bits"]),
            core_rule=params["core_rule"],
            core_min_coverage=float(params["core_min_coverage"]),
            fmcs_timeout=int(params["fmcs_timeout"]),
        )
        standalone_grouping_seconds = time.perf_counter() - grouping_started
        standalone_group_count = len(groups)

    timings["write"] = io_seconds
    measured = sum(timings[stage] for stage in SCALE_STAGES if stage != "other")
    timings["other"] += max(0.0, cli_wall_seconds - measured)
    stage_rows = []
    failed_stage = failure["stage"]
    passed_failure = False
    for stage in SCALE_STAGES:
        if failed_stage == stage:
            status = "failed"
            passed_failure = True
        elif failed_stage and not passed_failure:
            status = "complete"
        elif failed_stage:
            status = "not_run"
        else:
            status = "complete"
        stage_rows.append(
            {"stage": stage, "seconds": float(timings[stage]), "status": status}
        )
    for name, seconds in manifest_timings.items():
        stage_rows.append(
            {
                "stage": "manifest_" + name,
                "seconds": float(seconds),
                "status": "complete",
            }
        )
    stage_rows.append(
        {
            "stage": "standalone_grouping",
            "seconds": (
                float(standalone_grouping_seconds)
                if exit_code == 0
                else None
            ),
            "status": "complete" if exit_code == 0 else "not_run",
        }
    )
    stage_rows.append(
        {
            "stage": "full_cli_wall",
            "seconds": float(cli_wall_seconds),
            "status": "complete" if exit_code == 0 else "failed",
        }
    )
    process_peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    process_peak_bytes = (
        int(process_peak)
        if platform.system() == "Darwin"
        else int(process_peak * 1024)
    )
    record = {
        "size_entries": args.size,
        "repeat": args.repeat,
        "status": "passed" if exit_code == 0 else "failed",
        "exit_code": exit_code,
        "started_utc": now_utc(),
        "load_start": load_start,
        "load_finish": list(os.getloadavg()),
        "cli_wall_seconds": cli_wall_seconds,
        "cli_peak_rss_bytes": cli_peak_bytes,
        "peak_rss_bytes": max(cli_peak_bytes, process_peak_bytes),
        "peak_rss_mib": max(cli_peak_bytes, process_peak_bytes) / 1024**2,
        "stages": stage_rows,
        "standalone_group_count": (
            int(standalone_group_count) if exit_code == 0 else None
        ),
        "assignment_method": (
            manifest.get("parameters", {}).get("assignment_method")
            if manifest
            else None
        ),
        "parameters": manifest.get("parameters"),
        "failure": failure,
        "command": display_command,
    }
    dump_json(Path(args.result_json), record)
    print(
        json.dumps(
            {
                "size": args.size,
                "repeat": args.repeat,
                "status": record["status"],
                "wall_seconds": cli_wall_seconds,
                "grouping_seconds": standalone_grouping_seconds,
                "assignment": record["assignment_method"],
                "failure": failure,
            },
            ensure_ascii=False,
        ),
        flush=True,
    )
    return 0


def expected_limit_error(record: dict) -> bool:
    message = str(record.get("failure", {}).get("message", ""))
    return (
        record.get("size_entries") == 10001
        and "Exact grouping is limited to 10,000 molecules" in message
    )


def run_scale_child(
    input_path: Path,
    output_dir: Path,
    size: int,
    repeat: int,
    result_path: Path,
    log_path: Path,
) -> dict:
    command = [
        sys.executable,
        "-B",
        str(Path(__file__).resolve()),
        "--mode",
        "run-one",
        "--input",
        str(input_path),
        "--output-dir",
        str(output_dir),
        "--size",
        str(size),
        "--repeat",
        str(repeat),
        "--result-json",
        str(result_path),
    ]
    env = dict(os.environ)
    env["PYTHONPATH"] = str(CGM / "src")
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env.setdefault(
        "MPLCONFIGDIR", str(Path(gettempdir()) / "chemgridmap_t1_matplotlib")
    )
    completed = subprocess.run(
        command,
        cwd=PROJECT,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    log_path.write_text(
        f"$ {' '.join(command)}\n\nSTDOUT\n{completed.stdout}\nSTDERR\n{completed.stderr}",
        encoding="utf-8",
    )
    if result_path.exists():
        return json.loads(result_path.read_text(encoding="utf-8"))
    return {
        "size_entries": size,
        "repeat": repeat,
        "status": "failed",
        "exit_code": completed.returncode,
        "cli_wall_seconds": np.nan,
        "peak_rss_bytes": np.nan,
        "peak_rss_mib": np.nan,
        "assignment_method": None,
        "stages": [],
        "failure": {
            "stage": "runner",
            "type": "ChildProcessError",
            "message": (completed.stderr or completed.stdout)[-4000:],
        },
    }


def flatten_scale_records(records: list[dict]) -> pd.DataFrame:
    rows = []
    for record in records:
        for stage in record["stages"]:
            rows.append(
                {
                    "size_entries": record["size_entries"],
                    "repeat": record["repeat"],
                    "stage": stage["stage"],
                    "seconds": stage["seconds"],
                    "stage_status": stage["status"],
                    "run_status": record["status"],
                    "cli_wall_seconds": record.get("cli_wall_seconds"),
                    "peak_rss_bytes": record.get("peak_rss_bytes"),
                    "peak_rss_mib": record.get("peak_rss_mib"),
                    "assignment_method": record.get("assignment_method"),
                    "standalone_group_count": record.get("standalone_group_count"),
                    "failure_stage": record.get("failure", {}).get("stage"),
                    "failure_type": record.get("failure", {}).get("type"),
                    "failure_message": record.get("failure", {}).get("message"),
                }
            )
        if not record["stages"]:
            rows.append(
                {
                    "size_entries": record["size_entries"],
                    "repeat": record["repeat"],
                    "stage": "runner",
                    "seconds": record.get("cli_wall_seconds"),
                    "stage_status": "failed",
                    "run_status": record["status"],
                    "cli_wall_seconds": record.get("cli_wall_seconds"),
                    "peak_rss_bytes": record.get("peak_rss_bytes"),
                    "peak_rss_mib": record.get("peak_rss_mib"),
                    "assignment_method": record.get("assignment_method"),
                    "standalone_group_count": None,
                    "failure_stage": record.get("failure", {}).get("stage"),
                    "failure_type": record.get("failure", {}).get("type"),
                    "failure_message": record.get("failure", {}).get("message"),
                }
            )
    return pd.DataFrame(rows)


def summarize_scale(records: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    summaries = []
    for (size, stage), group in records.groupby(["size_entries", "stage"], sort=True):
        timings = group.loc[
            group.stage_status.isin(["complete", "failed"]), "seconds"
        ].dropna()
        rss = group.peak_rss_mib.dropna()
        summaries.append(
            {
                "size_entries": int(size),
                "stage": stage,
                "runs": int(len(group)),
                "timed_runs": int(len(timings)),
                "median_seconds": float(timings.median()) if len(timings) else np.nan,
                "min_seconds": float(timings.min()) if len(timings) else np.nan,
                "max_seconds": float(timings.max()) if len(timings) else np.nan,
                "median_peak_rss_mib": float(rss.median()) if len(rss) else np.nan,
                "min_peak_rss_mib": float(rss.min()) if len(rss) else np.nan,
                "max_peak_rss_mib": float(rss.max()) if len(rss) else np.nan,
                "passed_runs": int(group.run_status.eq("passed").sum()),
                "failed_runs": int(group.run_status.eq("failed").sum()),
                "assignment_methods": ";".join(
                    sorted(set(group.assignment_method.dropna().astype(str)))
                ),
            }
        )
    summary = pd.DataFrame(summaries)

    def fit_for(stage: str) -> dict:
        selected = records.loc[
            records.stage.eq(stage) & records.run_status.eq("passed")
        ]
        medians = selected.groupby("size_entries").seconds.median()
        medians = medians.loc[medians.index <= 10000]
        if len(medians) < 2 or medians.le(0).any():
            return {
                "sizes": medians.index.astype(int).tolist(),
                "slope": None,
                "intercept": None,
                "r_squared": None,
                "descriptive_only": True,
            }
        x = np.log(medians.index.to_numpy(float))
        y = np.log(medians.to_numpy(float))
        slope, intercept = np.polyfit(x, y, 1)
        predicted = slope * x + intercept
        ss_res = float(np.sum((y - predicted) ** 2))
        ss_total = float(np.sum((y - y.mean()) ** 2))
        return {
            "sizes": medians.index.astype(int).tolist(),
            "slope": float(slope),
            "intercept": float(intercept),
            "r_squared": float(1 - ss_res / ss_total) if ss_total else None,
            "descriptive_only": True,
        }

    return summary, {
        "full_cli_wall": fit_for("full_cli_wall"),
        "standalone_grouping": fit_for("standalone_grouping"),
    }


def run_scale() -> dict:
    input_manifest = prepare_scale_inputs()
    wait = wait_for_low_load()
    started = now_utc()
    load_start = list(os.getloadavg())
    records = []
    for size in SCALE_SIZES:
        for repeat in range(1, 4):
            run_dir = OUT / "runs" / str(size) / str(repeat)
            run_dir.mkdir(parents=True, exist_ok=True)
            result_path = run_dir / "timing.json"
            record = run_scale_child(
                OUT / "inputs" / f"input_{size}.csv",
                run_dir / "output",
                size,
                repeat,
                result_path,
                OUT / "logs" / f"run_{size}_{repeat}.log",
            )
            records.append(record)
            print(
                f"scale: {size} entries repeat {repeat}/3 "
                f"{record['status']} in {record.get('cli_wall_seconds', float('nan')):.1f}s",
                flush=True,
            )
    limit_repeat = 1
    limit_dir = OUT / "runs" / "10001" / "1"
    limit_dir.mkdir(parents=True, exist_ok=True)
    limit = run_scale_child(
        OUT / "inputs" / "input_10001.csv",
        limit_dir / "output",
        10001,
        limit_repeat,
        limit_dir / "timing.json",
        OUT / "logs" / "run_10001_1.log",
    )
    if expected_limit_error(limit):
        limit["status"] = "expected_error"
    records.append(limit)
    flat = flatten_scale_records(records)
    flat.to_csv(OUT / "scaling.csv", index=False, na_rep="NA", lineterminator="\n")
    summary, fits = summarize_scale(flat)
    summary.to_csv(
        OUT / "scaling_summary.csv", index=False, na_rep="NA", lineterminator="\n"
    )
    metadata = {
        "started_utc": started,
        "finished_utc": now_utc(),
        "load_wait": wait,
        "load_start": load_start,
        "load_finish": list(os.getloadavg()),
        "machine": machine_info(),
        "input_manifest": input_manifest,
        "repeats_per_size": 3,
        "sizes": list(SCALE_SIZES),
        "limit_probe": {
            "size_entries": 10001,
            "runs": 1,
            "status": limit["status"],
            "failure": limit.get("failure"),
        },
        "fits": fits,
        "default_cli_parameters": (
            "Only input/output and column mapping were specified. Assignment method, "
            "sparse neighbour cap, projection, representation, grid refinement, core "
            "rule, coverage target, evidence model, grouping, and output formats used "
            "ChemGridMap 0.5.0 defaults."
        ),
        "stage_timing": (
            "Function wrappers record the named stages in memory; manifest_* stages "
            "come from scale_run_manifest.json. standalone_grouping is a separate "
            "structural_groups call on the CLI's exported coordinates."
        ),
    }
    dump_json(OUT / "timing_metadata.json", metadata)
    return {
        "records": flat,
        "summary": summary,
        "metadata": metadata,
        "limit": limit,
    }


def rerun_limit_probe() -> int:
    """Rerun only the 10,001 limit probe and rebuild its combined outputs."""
    wait = wait_for_low_load()
    limit_dir = OUT / "runs" / "10001" / "1"
    limit_dir.mkdir(parents=True, exist_ok=True)
    limit = run_scale_child(
        OUT / "inputs" / "input_10001.csv",
        limit_dir / "output",
        10001,
        1,
        limit_dir / "timing.json",
        OUT / "logs" / "run_10001_2.log",
    )
    if expected_limit_error(limit):
        limit["status"] = "expected_error"
    if limit["status"] != "expected_error":
        raise RuntimeError(
            "10,001-entry probe did not produce the expected grouping error: "
            + json.dumps(limit.get("failure"), ensure_ascii=False)
        )
    records = []
    for size in SCALE_SIZES:
        for repeat in range(1, 4):
            records.append(
                json.loads(
                    (OUT / "runs" / str(size) / str(repeat) / "timing.json").read_text(
                        encoding="utf-8"
                    )
                )
            )
    records.append(limit)
    flat = flatten_scale_records(records)
    flat.to_csv(OUT / "scaling.csv", index=False, na_rep="NA", lineterminator="\n")
    summary, fits = summarize_scale(flat)
    summary.to_csv(
        OUT / "scaling_summary.csv", index=False, na_rep="NA", lineterminator="\n"
    )
    metadata = json.loads((OUT / "timing_metadata.json").read_text(encoding="utf-8"))
    metadata["limit_probe"] = {
        "size_entries": 10001,
        "runs": 2,
        "successful_evidence_runs": 1,
        "status": limit["status"],
        "failure": limit.get("failure"),
        "first_attempt": (
            "The CLI reached the grouping-limit failure, but the in-memory "
            "result writer rejected NaN in a not-run stage. That tooling failure "
            "is retained as run_10001_1_initial_serialization_failure.log; the "
            "probe was repeated after converting not-run stage times to null."
        ),
    }
    metadata["limit_probe_load_wait"] = wait
    metadata["fits"] = fits
    metadata["finished_utc"] = now_utc()
    dump_json(OUT / "timing_metadata.json", metadata)
    write_results()
    print(
        "limit probe: "
        + json.dumps(limit["failure"], ensure_ascii=False),
        flush=True,
    )
    return 0


def markdown_table(frame: pd.DataFrame, columns: list[str], digits: int = 6) -> str:
    def cell(value):
        if pd.isna(value):
            return "NA"
        if isinstance(value, (float, np.floating)):
            return f"{float(value):.{digits}g}"
        return str(value).replace("|", "\\|")

    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    lines.extend(
        "| " + " | ".join(cell(value) for value in row) + " |"
        for row in frame[columns].itertuples(index=False, name=None)
    )
    return "\n".join(lines)


def write_results(*args, **kwargs):
    return None



def write_final_report(*args, **kwargs):
    return None



def finalize(*args, **kwargs):
    return 0



def run_all() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    (OUT / "logs").mkdir(parents=True, exist_ok=True)
    dump_json(QA / "machine.json", machine_info())
    run_core_timing()
    run_scale()
    print("T1 compute complete", flush=True)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode",
        choices=["all", "run-one", "rerun-limit"],
        default="all",
    )
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--size", type=int)
    parser.add_argument("--repeat", type=int)
    parser.add_argument("--result-json", type=Path)
    args = parser.parse_args()
    if args.mode == "run-one":
        if not all(
            (args.input, args.output_dir, args.size, args.repeat, args.result_json)
        ):
            parser.error(
                "run-one requires --input, --output-dir, --size, --repeat and --result-json"
            )
        return _run_instrumented(args)
    if args.mode == "rerun-limit":
        return rerun_limit_probe()
    return run_all()


if __name__ == "__main__":
    raise SystemExit(main())
