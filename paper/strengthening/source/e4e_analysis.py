"""E4e: compare ChemGridMap 0.5.0 with evaluated reliability outputs."""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
import pandas as pd

from chemgridmap.evidence import (
    annotate_activity_confidence,
    document_values,
    fit_between_unit_sd,
)
from chemgridmap.reliability import (
    grid_neighbor_sets,
    neighbor_recall,
    representation_neighbor_sets,
)
from chemgridmap.representations import morgan_fingerprints


SOURCE = Path(__file__).resolve().parent
STRENGTHENING = SOURCE.parent
REV = STRENGTHENING.parent
WORKSPACE = REV.parents[1]
CGM = WORKSPACE / "chemgridmap"

OUT = STRENGTHENING / "output" / "E4"
QA = STRENGTHENING / "qa" / "E4"
RESULT = OUT / "regression_release_vs_evaluated.csv"

FROZEN_MAP = CGM / "paper/output/reliability_v1/activity/mapped_entry_confidence.csv"
FROZEN_LAYOUT = CGM / "paper/output/reliability_v1/layouts/chembl205/coordinates_seed_42.csv"
FROZEN_SIGMA = STRENGTHENING / "output/E2/sigma_doc_by_target.csv"
E2_INCLUSION = STRENGTHENING / "output/E2/target_inclusion.csv"
E2_RECORDS = STRENGTHENING / "qa/E2/curated_records.pkl.gz"
CAII_RAW = CGM / "validation_data/chembl37_chembl205/chembl37_chembl205_ic50_raw.csv"
REVISION = CGM / "paper/output/revision_v6"
PYTHON = CGM / ".venv/bin/python"

ANALYSIS_ROLE = "post hoc"
SEED = 42
TOLERANCE = 1e-12
SUPPORT_COLUMNS = [
    "document_balanced_pchembl",
    "consensus_scale",
    "p_consensus_inactive",
    "p_consensus_medium",
    "p_consensus_active",
    "p_new_document_inactive",
    "p_new_document_medium",
    "p_new_document_active",
    "class_confidence",
    "class_reproducibility",
]
SUPPORT_LABELS = [
    "document_balanced_class",
    "class_changed_by_document_balance",
    "class_status",
]


def _row(
    check: str,
    scope: str,
    n_items: int,
    max_error: float,
    n_differences: int,
    fields: str,
    reference: str,
    explanation: str,
    differing_items: list[str] | None = None,
) -> dict:
    return {
        "analysis_role": ANALYSIS_ROLE,
        "check": check,
        "scope": scope,
        "n_items": int(n_items),
        "max_absolute_error": float(max_error),
        "n_differences": int(n_differences),
        "fields_checked": fields,
        "reference": reference,
        "passed": bool(n_differences == 0),
        "differing_items": ";".join(differing_items or []),
        "explanation": explanation,
    }


def mapped_support_regression() -> list[dict]:
    """Recompute the three frozen map-entry Gaussian support annotations."""
    frozen = pd.read_csv(FROZEN_MAP, float_precision="round_trip")
    rows = []
    targets = {"chembl204": "Thrombin", "chembl205": "CA II", "chembl240": "hERG"}
    for target, label in targets.items():
        molecules = pd.read_csv(REVISION / target / "map_molecules.csv", low_memory=False)
        records = pd.read_csv(
            REVISION / target / f"{target}_retained_records.csv",
            low_memory=False,
        )
        current, _ = annotate_activity_confidence(molecules, records)
        expected = frozen.loc[frozen.target.eq(target)].set_index(
            "molecule_identity_key"
        )
        actual = current.set_index("molecule_identity_key")
        if not expected.index.is_unique or not actual.index.is_unique:
            raise ValueError(f"{target}: mapped identities are not unique")
        if set(actual.index) != set(expected.index):
            raise ValueError(f"{target}: current and frozen mapped identities differ")
        expected = expected.loc[actual.index]

        errors = {
            field: np.abs(
                actual[field].to_numpy(float) - expected[field].to_numpy(float)
            )
            for field in SUPPORT_COLUMNS
        }
        numeric_difference = np.zeros(len(actual), dtype=bool)
        for values in errors.values():
            numeric_difference |= values > TOLERANCE
        label_difference = np.zeros(len(actual), dtype=bool)
        for field in SUPPORT_LABELS:
            label_difference |= (
                actual[field].astype(str).to_numpy()
                != expected[field].astype(str).to_numpy()
            )
        different = numeric_difference | label_difference
        max_field = max(errors, key=lambda field: float(errors[field].max()))
        explanation = (
            "No difference after round-trip CSV parsing."
            if not different.any()
            else "Current 0.5.0 support differs from the frozen evaluation; see differing_items."
        )
        rows.append(
            _row(
                check="mapped_class_support",
                scope=f"{target} ({label})",
                n_items=len(actual),
                max_error=float(errors[max_field].max()),
                n_differences=int(different.sum()),
                fields=f"{max_field}; {', '.join(SUPPORT_LABELS)}",
                reference=FROZEN_MAP.name,
                explanation=explanation,
                differing_items=actual.index[different].astype(str).tolist(),
            )
        )
    return rows


def sigma_doc_regression() -> list[dict]:
    """Recompute full-data pooled document SD for every included E2 dataset."""
    frozen = pd.read_csv(FROZEN_SIGMA, float_precision="round_trip")
    frozen = frozen.loc[
        frozen.row_type.eq("target")
        & frozen.measurement_type.eq("all")
    ].set_index(["target", "endpoint"])
    if not frozen.index.is_unique:
        raise ValueError("Frozen sigma target/endpoint keys are not unique")
    inclusion = pd.read_csv(E2_INCLUSION)
    keys = list(
        inclusion.loc[inclusion.included.astype(bool), ["target", "endpoint"]]
        .itertuples(index=False, name=None)
    )
    if len(keys) != 52:
        raise ValueError(f"Expected 52 included E2 datasets, found {len(keys)}")
    if set(keys) != set(frozen.index):
        raise ValueError("E2 inclusion and frozen sigma table dataset sets differ")
    cached = pd.read_pickle(E2_RECORDS)
    rows = []
    for target, endpoint in keys:
        doc = document_values(cached[(target, endpoint)], merge_identical=True)
        current = float(fit_between_unit_sd(doc)["tau"])
        expected = float(frozen.loc[(target, endpoint), "sigma_doc"])
        error = abs(current - expected)
        rows.append(
            _row(
                check="sigma_doc_full_data",
                scope=f"{target}:{endpoint}",
                n_items=len(doc),
                max_error=error,
                n_differences=int(error > TOLERANCE),
                fields="sigma_doc",
                reference=FROZEN_SIGMA.name,
                explanation=(
                    "No difference."
                    if error <= TOLERANCE
                    else f"Current={current:.17g}; frozen={expected:.17g}."
                ),
                differing_items=[] if error <= TOLERANCE else [f"{target}:{endpoint}"],
            )
        )
    return rows


def _max_item_difference(current: pd.Series, frozen: pd.Series) -> tuple[float, list[str]]:
    errors = np.abs(current.astype(float).to_numpy() - frozen.astype(float).to_numpy())
    differing = current.index[errors > TOLERANCE].astype(str).tolist()
    return float(errors.max()), differing


def cli_layout_regression() -> list[dict]:
    """Run the default CA II CLI path and compare it with the frozen seed-42 layout."""
    environment = dict(os.environ)
    environment.update(
        {
            "PYTHONDONTWRITEBYTECODE": "1",
            "NUMBA_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "OMP_NUM_THREADS": "1",
            "VECLIB_MAXIMUM_THREADS": "1",
        }
    )
    with tempfile.TemporaryDirectory(prefix="e4e_caii_") as temporary:
        directory = Path(temporary)
        environment["MPLCONFIGDIR"] = str(directory / "mpl")
        command = [
            str(PYTHON),
            "-B",
            "-m",
            "chemgridmap",
            str(CAII_RAW),
            "--input-format",
            "chembl",
            "--target-id",
            "CHEMBL205",
            "--random-state",
            str(SEED),
            "--output-dir",
            str(directory),
            "--name",
            "e4e_caii",
        ]
        completed = subprocess.run(
            command,
            cwd=CGM,
            env=environment,
            capture_output=True,
            text=True,
        )
        if completed.returncode:
            raise RuntimeError(
                "CA II CLI failed:\n"
                f"stdout:\n{completed.stdout}\n"
                f"stderr:\n{completed.stderr}"
            )
        cli = pd.read_csv(directory / "e4e_caii_grid_coordinates.csv")
        frozen = pd.read_csv(FROZEN_LAYOUT)

    key = "molecule_identity_key"
    cli = cli.set_index(key)
    frozen = frozen.set_index(key)
    if not cli.index.is_unique or not frozen.index.is_unique:
        raise ValueError("CA II CLI or frozen layout identities are not unique")
    if set(cli.index) != set(frozen.index):
        raise ValueError("CA II CLI and frozen layout identities differ")
    frozen = frozen.loc[cli.index]

    rows = []
    for current_field, frozen_field in (
        ("grid_col", "refined_grid_col"),
        ("grid_row", "refined_grid_row"),
        ("projection_x", "projection_x"),
        ("projection_y", "projection_y"),
    ):
        error, differing = _max_item_difference(
            cli[current_field], frozen[frozen_field]
        )
        rows.append(
            _row(
                check="cli_layout_coordinate",
                scope=f"chembl205 seed={SEED}",
                n_items=len(cli),
                max_error=error,
                n_differences=len(differing),
                fields=f"cli.{current_field} vs frozen.{frozen_field}",
                reference=FROZEN_LAYOUT.name,
                explanation=(
                    "No difference."
                    if not differing
                    else "Default 0.5.0 layout differs from the frozen seed-42 layout."
                ),
                differing_items=differing,
            )
        )

    fingerprints = morgan_fingerprints(cli.canonical_smiles.tolist())
    reference = representation_neighbor_sets(fingerprints, k=10)
    cli_recall = neighbor_recall(
        reference,
        grid_neighbor_sets(cli[["grid_col", "grid_row"]].to_numpy(float), k=10),
    )
    frozen_recall = neighbor_recall(
        reference,
        grid_neighbor_sets(
            frozen[["refined_grid_col", "refined_grid_row"]].to_numpy(float),
            k=10,
        ),
    )
    recall_error, differing = _max_item_difference(
        pd.Series(cli_recall, index=cli.index),
        pd.Series(frozen_recall, index=cli.index),
    )
    rows.append(
        _row(
            check="cli_layout_recall",
            scope=f"chembl205 seed={SEED}",
            n_items=len(cli),
            max_error=recall_error,
            n_differences=len(differing),
            fields="per-molecule k=10 representation-neighbor recall",
            reference=FROZEN_LAYOUT.name,
            explanation=(
                f"Mean recall: current={np.mean(cli_recall):.12g}; "
                f"frozen={np.mean(frozen_recall):.12g}. "
                + (
                    "No per-molecule difference."
                    if not differing
                    else "Differences arise from the layout comparison."
                )
            ),
            differing_items=differing,
        )
    )
    return rows


def run_regression() -> pd.DataFrame:
    OUT.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    rows = mapped_support_regression() + sigma_doc_regression() + cli_layout_regression()
    result = pd.DataFrame(rows)
    result.to_csv(RESULT, index=False, float_format="%.17g", lineterminator="\n")
    return result


if __name__ == "__main__":
    table = run_regression()
    print(
        f"wrote {RESULT}: {len(table)} checks, "
        f"{int((table.n_differences > 0).sum())} with differences, "
        f"max error={table.max_absolute_error.max():.17g}"
    )
