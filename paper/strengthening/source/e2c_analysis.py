"""E2c: E1-B neighbourhood results on the 80 sampled target-endpoints."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[3]
REV = ROOT / "paper"
S = REV / "strengthening"
API = S / "data/chembl37_api"
ACTIVITIES = API / "activities"
OUT = S / "output/E2c"
QA = S / "qa/E2c"
FROZEN = ROOT / "paper/output/revision_v6"
sys.path[:0] = [
    str(S / "source"),
    str(ROOT / "paper"),
    str(ROOT / "src"),
]

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

from chemgridmap.core import assign_to_grid, normalize_coordinates
from chemgridmap.projection import project_2d
from chemgridmap.refinement import REFINEMENT_DEFAULTS, refine_grid_layout
from chemgridmap.reliability import grid_neighbor_sets, representation_neighbor_sets
from chemgridmap.representations import morgan_fingerprints
from e1b_neighbors import bootstrap_spearman, neighborhood_metrics, spearman_rho
from e2_models import curate_activity_table
from grid_baselines import random_grid_layout

SEALED_PLAN = S / "plans/E2c_task_draft_sealed.md"
SEALED_SHA256 = "95d67fd71c384fbb39c31170786b75acedd4898033e89d2fc636d4dcc461cbb0"
MAP_SEED = 42
SUBSET_SEED = 20261102
SUBSET_LIMIT = 1200
MIN_ENTRIES = 200
K = 10
MOLECULE_BOOTSTRAP_N = 1000
MOLECULE_BOOTSTRAP_SEED = 20261031
TARGET_BOOTSTRAP_N = 2000
TARGET_BOOTSTRAP_SEED = 20261028
LAYOUTS = ("projection", "minimum_displacement", "neighbor_refined", "random")
TESTED_LAYOUTS = ("neighbor_refined", "minimum_displacement")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_safe(value):
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def subset_positions(entry_ids, limit: int = SUBSET_LIMIT,
                     seed: int = SUBSET_SEED) -> np.ndarray:
    """Choose rows from identifiers alone; activity is not an argument."""
    entry_ids = np.asarray(entry_ids, dtype=str)
    if len(entry_ids) <= limit:
        return np.arange(len(entry_ids))
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(len(entry_ids), size=limit, replace=False))


def summarize_entries(records: pd.DataFrame) -> pd.DataFrame:
    """Create one parent entry per identity with the displayed pChEMBL median."""
    records = records.copy()
    records["molecule_identity_key"] = records.molecule_identity_key.astype(str)
    # Identity is derived from the normalized parent structure, not the salt form
    # retained in a source record's canonical_smiles column.
    unique_smiles = records.groupby("molecule_identity_key").representation_smiles.nunique()
    if (unique_smiles > 1).any():
        bad = unique_smiles[unique_smiles > 1].index.tolist()[:5]
        raise ValueError(f"One parent identity has multiple canonical SMILES: {bad}")
    entries = (
        records.groupby("molecule_identity_key", sort=True)
        .agg(
            canonical_smiles=("representation_smiles", "first"),
            activity_pchembl=("pchembl_value", "median"),
            retained_records=("pchembl_value", "size"),
        )
        .reset_index()
        .sort_values("molecule_identity_key")
        .reset_index(drop=True)
    )
    if not np.isfinite(entries.activity_pchembl.to_numpy(float)).all():
        raise ValueError("A displayed pChEMBL median is not finite.")
    return entries


def layout_metrics(activity, fingerprints, points, reference=None):
    """Apply the accepted E1-B definitions to one target and layout."""
    reference = reference if reference is not None else representation_neighbor_sets(
        fingerprints, k=K
    )
    return neighborhood_metrics(
        activity, reference, grid_neighbor_sets(np.asarray(points, dtype=float), k=K)
    )


def evaluate_layout(activity, fingerprints, points, layout: str,
                    reference=None) -> dict:
    metrics = layout_metrics(activity, fingerprints, points, reference)
    row = {
        "layout": layout,
        "n": int(len(metrics)),
        "mean_recall": float(metrics.recall.mean()),
        "mean_e_grid": float(metrics.e_grid.mean()),
        "mean_e_fp": float(metrics.e_fp.mean()),
        "rho": spearman_rho(metrics.recall, metrics.excess_error),
        "rho_ci_low": np.nan,
        "rho_ci_high": np.nan,
        "bootstrap_valid_count": np.nan,
    }
    if layout in TESTED_LAYOUTS:
        row.update(bootstrap_spearman(
            metrics.recall,
            metrics.excess_error,
            n_boot=MOLECULE_BOOTSTRAP_N,
            seed=MOLECULE_BOOTSTRAP_SEED,
        ))
    return row


def checkpoint_path(target: str, endpoint: str) -> Path:
    return QA / "targets" / f"{target}_{endpoint}.csv"


def pipeline_hash() -> str:
    return sha256(Path(__file__))


def process_target(item) -> pd.DataFrame:
    target, endpoint = str(item.target_chembl_id), str(item.standard_type)
    path = ACTIVITIES / f"{target}_{endpoint}.csv"
    cp = checkpoint_path(target, endpoint)
    if cp.exists():
        cached = pd.read_csv(cp)
        if (
            len(cached) == len(LAYOUTS)
            and cached.target.eq(target).all()
            and cached.endpoint.eq(endpoint).all()
            and cached.pipeline_sha256.eq(pipeline_hash()).all()
        ):
            print(f"E2c {target} {endpoint}: checkpoint reused", flush=True)
            return cached

    base = {
        "target": target,
        "endpoint": endpoint,
        "role": "sampled",
        "api_records": np.nan,
        "retained_records": np.nan,
        "entries": np.nan,
        "subset_n": np.nan,
        "subset_seed": SUBSET_SEED,
        "included": False,
        "reason": "",
        "input_file": str(path.relative_to(ROOT)),
        "pipeline_sha256": pipeline_hash(),
    }
    try:
        raw = pd.read_csv(path, low_memory=False)
        records = curate_activity_table(raw, target, endpoint).retained_records
        entries = summarize_entries(records)
        base.update({
            "api_records": int(len(raw)),
            "retained_records": int(len(records)),
            "entries": int(len(entries)),
        })
        if len(entries) < MIN_ENTRIES:
            base["reason"] = f"fewer than {MIN_ENTRIES} curated parent entries"
            return pd.DataFrame([base])

        positions = subset_positions(entries.molecule_identity_key.to_numpy())
        selected = entries.iloc[positions].reset_index(drop=True)
        base.update({
            "subset_n": int(len(selected)),
            "included": True,
            "reason": "included",
        })
        fingerprints = morgan_fingerprints(selected.canonical_smiles.tolist())
        projection_raw = project_2d(
            fingerprints,
            method="umap",
            random_state=MAP_SEED,
            umap_metric="jaccard",
        )
        projection = normalize_coordinates(projection_raw, mode="isotropic")
        minimum, side, _ = assign_to_grid(
            projection_raw,
            padding=None,
            target_occupancy=0.40,
            coordinate_scaling="isotropic",
            method="dense",
        )
        cells = np.rint(minimum * (side - 1)).astype(int)
        reference = representation_neighbor_sets(fingerprints, k=K)
        refined = refine_grid_layout(
            cells,
            projection * (side - 1),
            reference,
            side,
            **REFINEMENT_DEFAULTS,
            random_state=MAP_SEED,
        )
        random = random_grid_layout(len(selected), (side, side), seed=MAP_SEED)
        points = {
            "projection": projection,
            "minimum_displacement": cells,
            "neighbor_refined": refined,
            "random": random,
        }
        rows = []
        for layout in LAYOUTS:
            row = evaluate_layout(
                selected.activity_pchembl.to_numpy(float),
                fingerprints,
                points[layout],
                layout,
                reference,
            )
            rows.append({**base, **row, "grid_side": side})
        result = pd.DataFrame(rows)
        random_e_grid = float(result.loc[result.layout.eq("random"), "mean_e_grid"].iloc[0])
        result["lambda"] = 1.0 - result.mean_e_fp / random_e_grid
        result.to_csv(cp, index=False)
        assert np.isfinite(result.mean_recall).all()
        assert np.isfinite(result.mean_e_grid).all()
        assert np.isfinite(result.mean_e_fp).all()
        print(
            f"E2c {target} {endpoint}: {len(selected)} entries, side {side}, "
            f"rho refined={result.loc[result.layout.eq('neighbor_refined'), 'rho'].iloc[0]:.6f}",
            flush=True,
        )
        return result
    except Exception as error:
        base["reason"] = f"pipeline failure: {type(error).__name__}: {error}"
        print(f"E2c {target} {endpoint}: FAILED {base['reason']}", flush=True)
        return pd.DataFrame([base])


def bootstrap_median(values, rng) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    draws = np.empty(TARGET_BOOTSTRAP_N)
    for index in range(TARGET_BOOTSTRAP_N):
        draws[index] = np.median(values[rng.integers(0, len(values), len(values))])
    return tuple(float(value) for value in np.quantile(draws, [0.025, 0.975]))


def holm_adjust(p_values) -> np.ndarray:
    p_values = np.asarray(p_values, dtype=float)
    order = np.argsort(p_values, kind="mergesort")
    adjusted = np.empty(len(p_values))
    running = 0.0
    for rank, index in enumerate(order):
        running = max(running, (len(p_values) - rank) * p_values[index])
        adjusted[index] = min(running, 1.0)
    return adjusted


def run_tests(layouts: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for endpoint in ("IC50", "Ki"):
        population = layouts.loc[
            layouts.endpoint.eq(endpoint) & layouts.included.astype(bool)
        ]
        refined = population.loc[
            population.layout.eq("neighbor_refined"),
            "rho",
        ].dropna().to_numpy(float)
        if len(refined) < 2:
            raise ValueError(f"L1 {endpoint} has fewer than two target values.")
        rng = np.random.default_rng(TARGET_BOOTSTRAP_SEED)
        low, high = bootstrap_median(refined, rng)
        test = wilcoxon(refined, alternative="less", zero_method="wilcox", method="auto")
        rows.append({
            "test_id": "L1",
            "endpoint": endpoint,
            "contrast": "median rho_refined < 0",
            "n_targets": len(refined),
            "statistic": float(test.statistic),
            "estimate": float(np.median(refined)),
            "ci_low": low,
            "ci_high": high,
            "raw_p": float(test.pvalue),
            "alternative": "less",
        })

        pivot = population.pivot(
            index=["target", "endpoint"],
            columns="layout",
            values="mean_e_grid",
        )
        differences = (
            pivot["neighbor_refined"] - pivot["minimum_displacement"]
        ).dropna().to_numpy(float)
        if len(differences) < 2:
            raise ValueError(f"L2 {endpoint} has fewer than two target pairs.")
        rng = np.random.default_rng(TARGET_BOOTSTRAP_SEED)
        low, high = bootstrap_median(differences, rng)
        test = wilcoxon(
            differences,
            alternative="less",
            zero_method="wilcox",
            method="auto",
        )
        rows.append({
            "test_id": "L2",
            "endpoint": endpoint,
            "contrast": "median e_grid(neighbor_refined) - e_grid(minimum_displacement) < 0",
            "n_targets": len(differences),
            "statistic": float(test.statistic),
            "estimate": float(np.median(differences)),
            "ci_low": low,
            "ci_high": high,
            "raw_p": float(test.pvalue),
            "alternative": "less",
        })
    result = pd.DataFrame(rows)
    result["holm_p"] = holm_adjust(result.raw_p)
    result["bootstrap_n"] = TARGET_BOOTSTRAP_N
    result["bootstrap_seed"] = TARGET_BOOTSTRAP_SEED
    result["ci_method"] = "target bootstrap percentile 95%"
    return result


def descriptive_summary(layouts: pd.DataFrame) -> dict:
    included = layouts.loc[layouts.included.astype(bool)]
    result = {}
    for endpoint in ("IC50", "Ki"):
        population = included.loc[
            (included.endpoint == endpoint)
            & (included.layout == "neighbor_refined")
        ].copy()
        result[endpoint] = {
            "n_targets": int(len(population)),
            "rho_ci_below_zero": int((population.rho_ci_high < 0).sum()),
            "rho_ci_includes_zero": int(
                ((population.rho_ci_low <= 0) & (population.rho_ci_high >= 0)).sum()
            ),
            "rho_ci_above_zero": int((population.rho_ci_low > 0).sum()),
            "rho_lambda_spearman": float(
                spearmanr(population.rho, population["lambda"]).statistic
            ),
        }
    refined = included.loc[included.layout.eq("neighbor_refined")]
    result["all"] = {
        "n_targets": int(len(refined)),
        "rho_lambda_spearman": float(
            spearmanr(refined.rho, refined["lambda"]).statistic
        ),
        "note": "Exploratory descriptive correlation; no inferential test performed.",
    }
    return result


def markdown_table(frame: pd.DataFrame) -> str:
    def cell(value):
        if isinstance(value, (float, np.floating)):
            return f"{value:.6g}"
        return str(value).replace("|", "\\|")
    lines = [
        "| " + " | ".join(frame.columns) + " |",
        "| " + " | ".join("---" for _ in frame.columns) + " |",
    ]
    lines.extend(
        "| " + " | ".join(cell(value) for value in row) + " |"
        for row in frame.itertuples(index=False, name=None)
    )
    return "\n".join(lines)


def write_results(layout_targets, layout_by_target, tests, descriptive) -> None:
    included = layout_targets.loc[layout_targets.included.astype(bool)]
    excluded = layout_targets.loc[~layout_targets.included.astype(bool)]
    table = tests[[
        "test_id", "endpoint", "n_targets", "estimate", "ci_low", "ci_high",
        "statistic", "raw_p", "holm_p",
    ]]
    lines = [
        "# E2c Results",
        "",
        "## Scope and inclusion",
        "",
        f"The sampled set contained {len(layout_targets)} target-endpoints. "
        f"{len(included)} met the planned minimum of {MIN_ENTRIES} curated parent entries; "
        f"{len(excluded)} did not. Endpoint counts were "
        f"{int((included.endpoint == 'IC50').sum())} IC50 and "
        f"{int((included.endpoint == 'Ki').sum())} Ki.",
        "",
        "The displayed activity was the median retained-record pChEMBL value within each "
        "parent identity. Targets above 1,200 entries used a deterministic random subset "
        "selected from entry identifiers with seed 20261102.",
        "",
        "## Planned hypothesis tests",
        "",
        "Wilcoxon signed-rank tests were one-sided. The four raw p-values were Holm-adjusted "
        "together. Intervals are 95% percentile intervals from 2,000 target-level bootstrap "
        "replicates with seed 20261028.",
        "",
        markdown_table(table),
        "",
        "## Descriptive results",
        "",
    ]
    for endpoint in ("IC50", "Ki"):
        item = descriptive[endpoint]
        lines.extend([
            f"- {endpoint} refined rho CI categories (n={item['n_targets']}): "
            f"{item['rho_ci_below_zero']} entirely below zero, "
            f"{item['rho_ci_includes_zero']} including zero, "
            f"{item['rho_ci_above_zero']} entirely above zero.",
            f"- {endpoint} Spearman correlation between rho_refined and lambda: "
            f"{item['rho_lambda_spearman']:.6g}.",
        ])
    lines.extend([
        f"- Across endpoints, Spearman correlation between rho_refined and lambda: "
        f"{descriptive['all']['rho_lambda_spearman']:.6g} (n={descriptive['all']['n_targets']}).",
        "- The rho-lambda association is exploratory and its p-value was not tested or reported.",
        "",
    ])
    if not excluded.empty:
        lines.extend([
            "## Excluded target-endpoints",
            "",
            markdown_table(excluded[["target", "endpoint", "entries", "reason"]]),
            "",
        ])
    (OUT / "RESULTS_E2c.md").write_text("\n".join(lines), encoding="utf-8")


def input_paths(selected: pd.DataFrame) -> list[Path]:
    paths = [SEALED_PLAN, API / "selected_targets.csv"]
    paths.extend(
        ACTIVITIES / f"{row.target_chembl_id}_{row.standard_type}.csv"
        for row in selected.itertuples(index=False)
    )
    return paths


def source_paths() -> list[Path]:
    sources = [
        Path(__file__),
        S / "source/test_e2c.py",
        S / "source/e1b_neighbors.py",
        S / "source/e2_models.py",
        S / "source/e1a_models.py",
        S / "source/e1a_statistics.py",
        ROOT / "paper/reliability_activity.py",
        ROOT / "paper/grid_baselines.py",
        ROOT / "src/chemgridmap/activity.py",
        ROOT / "src/chemgridmap/chembl.py",
        ROOT / "src/chemgridmap/evidence.py",
        ROOT / "src/chemgridmap/projection.py",
        ROOT / "src/chemgridmap/core.py",
        ROOT / "src/chemgridmap/refinement.py",
        ROOT / "src/chemgridmap/reliability.py",
        ROOT / "src/chemgridmap/representations.py",
    ]
    return [path for path in sources if path.exists()]


def write_manifest(inputs, sources, started: str, ended: str) -> dict:
    outputs = [
        path for path in sorted(OUT.iterdir())
        if path.is_file() and path.name != "manifest.json"
    ]
    manifest = {
        "task": "E2c",
        "sealed_plan_sha256": sha256(SEALED_PLAN),
        "started_at": started,
        "ended_at": ended,
        "seeds": {
            "map": MAP_SEED,
            "subset": SUBSET_SEED,
            "molecule_bootstrap": MOLECULE_BOOTSTRAP_SEED,
            "target_bootstrap": TARGET_BOOTSTRAP_SEED,
        },
        "parameters": {
            "minimum_entries": MIN_ENTRIES,
            "subset_limit": SUBSET_LIMIT,
            "projection": {
                "method": "UMAP",
                "metric": "Jaccard",
                "n_neighbors": 15,
                "min_dist": 0.1,
                "n_components": 2,
                "random_state": MAP_SEED,
            },
            "minimum_displacement": {
                "occupancy": 0.40,
                "coordinate_scaling": "isotropic",
                "assignment_method": "dense",
            },
            "neighbor_refined": {
                "parameters": json_safe(REFINEMENT_DEFAULTS),
                "random_state": MAP_SEED,
            },
            "random": {"seed": MAP_SEED},
            "k": K,
            "molecule_bootstrap_n": MOLECULE_BOOTSTRAP_N,
            "target_bootstrap_n": TARGET_BOOTSTRAP_N,
            "tests": [
                "one-sided Wilcoxon signed-rank L1 by endpoint",
                "one-sided paired Wilcoxon signed-rank L2 by endpoint",
                "Holm adjustment across four p-values",
            ],
        },
        "inputs_sha256": {
            str(path.relative_to(ROOT)): sha256(path) for path in inputs
        },
        "source_sha256": {
            str(path.relative_to(ROOT)): sha256(path) for path in sources
        },
        "outputs_sha256": {
            str(path.relative_to(ROOT)): sha256(path) for path in outputs
        },
        "software_versions": {
            "python": platform.python_version(),
            **{
                name: importlib.metadata.version(name)
                for name in ("numpy", "pandas", "scipy", "scikit-learn", "umap-learn")
            },
        },
        "manifest_self_hash": "omitted to avoid a self-referential hash",
    }
    (OUT / "manifest.json").write_text(
        json.dumps(json_safe(manifest), indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return manifest


def run() -> None:
    started = datetime.now(timezone.utc).isoformat()
    if sha256(SEALED_PLAN) != SEALED_SHA256:
        raise RuntimeError("The sealed E2c plan hash does not match the activation.")
    OUT.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    (QA / "targets").mkdir(exist_ok=True)
    (QA / "layouts").mkdir(exist_ok=True)
    selected = pd.read_csv(API / "selected_targets.csv")
    selected = selected.loc[selected.role.eq("sampled")].reset_index(drop=True)
    if len(selected) != 80:
        raise ValueError(f"Expected 80 sampled target-endpoints, found {len(selected)}")
    if selected.standard_type.value_counts().to_dict() != {"IC50": 50, "Ki": 30}:
        raise ValueError("Sampled endpoint counts do not match the sealed plan.")

    frames = []
    for index, item in enumerate(selected.itertuples(index=False), 1):
        print(f"START {index}/80 {item.target_chembl_id} {item.standard_type}", flush=True)
        frames.append(process_target(item))
        combined = pd.concat(frames, ignore_index=True)
        combined.to_csv(OUT / "layout_by_target.csv", index=False)

    layouts = pd.concat(frames, ignore_index=True)
    included = layouts.loc[layouts.included.astype(bool)].copy()
    if len(included) == 0:
        raise ValueError("No target-endpoint met the inclusion rule.")
    layout_targets = (
        layouts.sort_values(["target", "endpoint", "layout"])
        .drop_duplicates(["target", "endpoint"])
        .copy()
    )
    refined = included.loc[included.layout.eq("neighbor_refined")].set_index(
        ["target", "endpoint"]
    )
    layout_targets = layout_targets.set_index(["target", "endpoint"])
    layout_targets["rho_refined"] = refined.rho
    layout_targets["rho_ci_low"] = refined.rho_ci_low
    layout_targets["rho_ci_high"] = refined.rho_ci_high
    layout_targets["lambda"] = refined["lambda"]
    layout_targets = layout_targets.reset_index()
    layout_targets.to_csv(OUT / "layout_targets.csv", index=False)
    layouts.to_csv(OUT / "layout_by_target.csv", index=False)

    tests = run_tests(layouts)
    tests.to_csv(OUT / "layout_tests.csv", index=False)
    descriptive = descriptive_summary(layouts)
    write_results(layout_targets, layouts, tests, descriptive)
    ended = datetime.now(timezone.utc).isoformat()
    manifest = write_manifest(
        input_paths(selected),
        source_paths(),
        started,
        ended,
    )
    print(
        f"COMPLETE included={int(layout_targets.included.astype(bool).sum())} "
        f"outputs={len(manifest['outputs_sha256'])}",
        flush=True,
    )
    print(tests.to_string(index=False), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    run()
