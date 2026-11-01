"""Reproduce core coverage under three stereochemistry settings."""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import sys

sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parents[3]
S = ROOT / "paper/strengthening"
OUT = S / "output/E12c"

sys.path[:0] = [str(S / "source"), str(ROOT / "src")]
import chemgridmap

if not Path(chemgridmap.__file__).resolve().is_relative_to(ROOT / "src"):
    raise RuntimeError(f"chemgridmap imported outside repository: {chemgridmap.__file__}")

import numpy as np
import pandas as pd
from rdkit import Chem
from manifest_utils import write_manifest

SETTINGS = ("original", "unified_use_chirality", "no_use_chirality")
RULES = ("coverage", "most_frequent", "fmcs")
SIZE_ORDER = ("2", "3-4", "5-9", ">=10")


def core_pattern(core: str, rule: str):
    if not isinstance(core, str) or not core:
        return None
    pattern = Chem.MolFromSmarts(core) if rule == "fmcs" else Chem.MolFromSmiles(core)
    if pattern is None:
        raise ValueError(f"Could not parse {rule} core: {core}")
    return pattern


def coverage_bucket(value: float) -> str:
    if not np.isfinite(value):
        return "no_core"
    if value >= 0.8 - 1e-12:
        return "coverage_ge_0_8"
    if value > 0.5 + 1e-12:
        return "coverage_gt_0_5_lt_0_8"
    return "coverage_le_0_5"


def load_smiles() -> dict[str, dict[str, str]]:
    result = {}
    for target in ("chembl204", "chembl205", "chembl240"):
        path = S / f"data/chembl37_api/activities/{target.upper()}_IC50.csv"
        frame = pd.read_csv(path, low_memory=False)
        frame = frame.loc[
            frame.pchembl_value.notna(),
            ["parent_molecule_chembl_id", "canonical_smiles"],
        ].dropna()
        grouped = frame.groupby("parent_molecule_chembl_id").canonical_smiles.agg(
            lambda values: sorted(set(values))
        )
        result[target] = {str(key): values[0] for key, values in grouped.items()}
    return result


def load_core_rows() -> pd.DataFrame:
    members = pd.read_csv(S / "output/E1B/core_group_members.csv")
    members = members[["target", "group", "size", "size_stratum", "member_ids"]]
    timing = pd.read_csv(S / "output/T1/core_timing_by_group.csv")
    selected = timing.loc[
        ((timing.rule != "fmcs") & timing.timeout_seconds.isna())
        | ((timing.rule == "fmcs") & (timing.timeout_seconds == 10))
    ].copy()
    if len(selected) != 480 * 3:
        raise RuntimeError(f"Expected 1,440 core rows; found {len(selected)}")
    result = selected.merge(
        members,
        on=["target", "group", "size"],
        how="left",
        validate="many_to_one",
    )
    if result.member_ids.isna().any():
        raise RuntimeError("Missing group membership for a core row")
    return result


def group_smiles(member_ids: str, target: str, smiles: dict[str, dict[str, str]]) -> list:
    molecules = []
    for key in str(member_ids).split("|"):
        value = key[7:] if key.startswith("SMILES:") else smiles[target].get(key)
        if value is None:
            raise KeyError(f"No SMILES for {target}/{key}")
        molecule = Chem.MolFromSmiles(value)
        if molecule is None:
            raise ValueError(f"Invalid SMILES for {target}/{key}")
        molecules.append(molecule)
    return molecules


def compute_settings() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    rows = load_core_rows()
    smiles = load_smiles()
    molecule_cache = {
        (row.target, row.group): group_smiles(row.member_ids, row.target, smiles)
        for row in rows[["target", "group", "member_ids"]].drop_duplicates().itertuples(
            index=False
        )
    }

    detail_rows = []
    for row in rows.itertuples(index=False):
        pattern = core_pattern(row.core, row.rule)
        molecules = molecule_cache[(row.target, row.group)]
        for setting in SETTINGS:
            use_chirality = (
                row.rule != "fmcs"
                if setting == "original"
                else setting == "unified_use_chirality"
            )
            count = (
                sum(
                    bool(molecule.HasSubstructMatch(pattern, useChirality=use_chirality))
                    for molecule in molecules
                )
                if pattern is not None
                else 0
            )
            detail_rows.append(
                {
                    "setting": setting,
                    "target": row.target,
                    "group": row.group,
                    "size": int(row.size),
                    "size_stratum": row.size_stratum,
                    "rule": row.rule,
                    "core": row.core,
                    "has_core": pattern is not None,
                    "use_chirality": use_chirality,
                    "coverage_count": int(count),
                    "coverage_fraction": count / len(molecules) if pattern is not None else np.nan,
                    "core_heavy_atoms": (
                        int(row.core_heavy_atoms) if pattern is not None else np.nan
                    ),
                }
            )
    detail = pd.DataFrame(detail_rows)

    grouped_rows = []
    for (setting, rule, stratum), frame in detail.groupby(
        ["setting", "rule", "size_stratum"], sort=True, observed=True
    ):
        buckets = frame.coverage_fraction.map(coverage_bucket)
        grouped_rows.append(
            {
                "setting": setting,
                "rule": rule,
                "size_stratum": stratum,
                "coverage_ge_0_8": int((buckets == "coverage_ge_0_8").sum()),
                "coverage_gt_0_5_lt_0_8": int(
                    (buckets == "coverage_gt_0_5_lt_0_8").sum()
                ),
                "coverage_le_0_5": int((buckets == "coverage_le_0_5").sum()),
                "no_ring_scaffold_groups": int((buckets == "no_core").sum()),
                "groups": len(frame),
                "mean_core_heavy_atoms": float(
                    frame.loc[frame.has_core, "core_heavy_atoms"].mean()
                ),
            }
        )
    grouped = pd.DataFrame(grouped_rows)
    grouped["size_stratum"] = pd.Categorical(
        grouped.size_stratum, categories=SIZE_ORDER, ordered=True
    )
    grouped = grouped.sort_values(["setting", "rule", "size_stratum"]).reset_index(drop=True)
    grouped["analysis_role"] = "post hoc"

    e12b = pd.read_csv(S / "output/E12b/core_chirality.csv")
    comparison_columns = [
        "setting",
        "rule",
        "size_stratum",
        "coverage_ge_0_8",
        "coverage_gt_0_5_lt_0_8",
        "coverage_le_0_5",
        "no_ring_scaffold_groups",
        "groups",
    ]
    reproduced = grouped.loc[
        grouped.setting.isin(["original", "unified_use_chirality"]),
        comparison_columns,
    ].reset_index(drop=True)
    expected = e12b.loc[:, comparison_columns].copy()
    expected["size_stratum"] = expected["size_stratum"].astype(str)
    reproduced = reproduced.copy()
    reproduced["size_stratum"] = reproduced["size_stratum"].astype(str)
    expected = expected.sort_values(["setting", "rule", "size_stratum"]).reset_index(drop=True)
    reproduced = reproduced.sort_values(
        ["setting", "rule", "size_stratum"]
    ).reset_index(drop=True)
    if not reproduced.equals(expected):
        difference = reproduced.compare(expected)
        raise RuntimeError(f"E12b core-setting regression failed:\n{difference}")

    original = detail.loc[
        detail.setting == "original",
        ["target", "group", "rule", "coverage_count"],
    ].set_index(["target", "group", "rule"])
    no_chirality = detail.loc[
        detail.setting == "no_use_chirality",
        ["target", "group", "rule", "coverage_count"],
    ].set_index(["target", "group", "rule"])
    changed = original.loc[original.coverage_count != no_chirality.coverage_count]
    changed_by_rule = changed.index.get_level_values("rule").value_counts().to_dict()
    changed_groups = set(changed.index.droplevel("rule"))

    majority = (
        grouped.assign(
            strict_majority=grouped.coverage_ge_0_8
            + grouped.coverage_gt_0_5_lt_0_8
        )
        .groupby(["setting", "rule"])
        .strict_majority.sum()
    )
    checks = {
        "e12b_two_setting_regression": True,
        "no_chirality_strict_majority": {
            rule: int(majority.loc[("no_use_chirality", rule)]) for rule in RULES
        },
        "changed_group_rule_pairs": int(len(changed)),
        "changed_groups": int(len(changed_groups)),
        "changed_pairs_by_rule": {
            rule: int(changed_by_rule.get(rule, 0)) for rule in RULES
        },
        "rows": len(grouped),
        "detail_rows": len(detail),
    }
    if (
        checks["no_chirality_strict_majority"]
        != {"coverage": 470, "most_frequent": 272, "fmcs": 480}
        or checks["changed_group_rule_pairs"] != 9
        or checks["changed_groups"] != 9
        or checks["changed_pairs_by_rule"]
        != {"coverage": 0, "most_frequent": 9, "fmcs": 0}
    ):
        raise RuntimeError(f"Unexpected no-stereochemistry checks: {checks}")
    return grouped, detail, checks


def main() -> int:
    started = datetime.now(timezone.utc).isoformat()
    OUT.mkdir(parents=True, exist_ok=True)

    grouped, _, checks = compute_settings()
    output_path = OUT / "core_stereo_settings.csv"
    grouped.to_csv(output_path, index=False, float_format="%.17g", na_rep="NA")

    inputs = [
        S / "output/E1B/core_group_members.csv",
        S / "output/T1/core_timing_by_group.csv",
        S / "output/E12b/core_chirality.csv",
        *[
            S / f"data/chembl37_api/activities/{target}_IC50.csv"
            for target in ("CHEMBL204", "CHEMBL205", "CHEMBL240")
        ],
    ]
    write_manifest(
        OUT,
        stage="E12c",
        inputs=inputs,
        sources=[
            Path(__file__),
            ROOT / "src/chemgridmap/__init__.py",
            ROOT / "src/chemgridmap/structures.py",
        ],
        parameters={"analysis_role": "post hoc", "checks": checks},
        started_at=started,
        finished_at=datetime.now(timezone.utc).isoformat(),
    )
    counts = ", ".join(f"{rule}={count}" for rule, count in checks["no_chirality_strict_majority"].items())
    print(
        f"E12C COMPLETE rows={checks['rows']} "
        f"changed_group_rule_pairs={checks['changed_group_rule_pairs']} "
        f"strict_majority_without_chirality {counts}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
