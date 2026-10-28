"""Recover the measurements behind a molecular grid cell without remapping it."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

from .chembl import detect_chembl_columns
from .plotting import render_svg
from .activity import RULE_KEYS, RULE_COLUMNS, activity_annotations
from .provenance import file_digests


def annotation_rules(row, parameters=None):
    present = [column in row.index for column in RULE_COLUMNS]
    if any(present) and not all(present):
        raise ValueError("Incomplete saved annotation rules; supply an intact coordinate table.")
    embedded = {key: float(row[column]) for key, column in zip(RULE_KEYS, RULE_COLUMNS)} if all(present) else None
    if parameters is None:
        return embedded
    if not set(RULE_KEYS).issubset(parameters):
        raise ValueError("Curation parameters lack activity thresholds or the conflict range.")
    explicit = {key: float(parameters[key]) for key in RULE_KEYS}
    if embedded is not None and explicit != embedded:
        raise ValueError("Saved annotation rules differ from the supplied curation parameters.")
    return explicit


def verify_annotations(row, evidence, values, parameters):
    rules = annotation_rules(row, parameters)
    calculated = {
        "pchembl_min": float(values.min()), "pchembl_max": float(values.max()),
        "pchembl_range": float(values.max() - values.min()), "pchembl_median": float(values.median()),
        "pchembl_std": float(values.std(ddof=1)),
    }
    if rules is not None:
        calculated.update(activity_annotations(values.to_numpy(), **rules))
    verified = []
    for key, value in calculated.items():
        if key not in row.index:
            continue
        try:
            agrees = str(row[key]) == value if isinstance(value, str) else bool(
                np.isclose(float(row[key]), value, rtol=0, atol=1e-10, equal_nan=True))
        except (ValueError, TypeError):
            agrees = False
        if not agrees:
            raise ValueError(f"Recalculated annotation {key} does not match the map.")
        verified.append(key)
    # Identity links also bind the displayed structure to the retained representation.
    if "representation_smiles" in evidence:
        expected = sorted(evidence.representation_smiles.astype(str).unique())[0]
        for key in ("canonical_smiles", "representation_smiles", "smiles"):
            if key in row.index and str(row[key]) != expected:
                raise ValueError(f"Displayed structure {key} does not match the retained records.")
            if key in row.index:
                verified.append(key)
    columns = detect_chembl_columns(evidence)
    if "n_assays" in row.index and columns["assay_id"]:
        if float(row.n_assays) != evidence[columns["assay_id"]].nunique():
            raise ValueError("Source assay count does not match the retained records.")
        verified.append("n_assays")
    for key in ("molecule_id", "parent_molecule_id", "assay_id", "document_id", "activity_id", "target_id"):
        field = columns[key]
        if key in row.index and field:
            expected = {str(value).strip() for value in evidence[field].dropna() if str(value).strip()}
            actual = set() if pd.isna(row[key]) else {value.strip() for value in str(row[key]).split(";") if value.strip()}
            if actual != expected:
                raise ValueError(f"Source identifier {key} does not match the retained records.")
            verified.append(key)
    required_flags = {"activity_class", "is_conflicted", "has_class_boundary_crossing"}
    quality_verified = rules is not None and required_flags.issubset(verified)
    return calculated, rules, verified, quality_verified


def inspect_entry(coordinates, records, *, molecule_id=None, cell=None, parameters=None):
    """Select one entry by exact ID or zero-based (row, column), then verify its records."""
    if (molecule_id is None) == (cell is None):
        raise ValueError("Select either molecule_id or cell=(row, column).")
    if cell is not None:
        selected = coordinates[coordinates.grid_row.eq(cell[0]) & coordinates.grid_col.eq(cell[1])]
    else:
        wanted = str(molecule_id).strip()
        is_chembl_id = re.fullmatch(r"CHEMBL\d+", wanted, flags=re.IGNORECASE) is not None
        if is_chembl_id:
            wanted = wanted.upper()
        matches = pd.Series(False, index=coordinates.index)
        for field in ("molecule_identity_key", "molecule_id", "parent_molecule_id"):
            if field in coordinates:
                # SMILES keys are case-sensitive: C (aliphatic) and c (aromatic) differ.
                matches |= coordinates[field].fillna("").astype(str).map(
                    lambda value: wanted in {v.strip().upper() if is_chembl_id else v.strip()
                                             for v in value.split(";")}
                )
        selected = coordinates[matches]
    if len(selected) != 1:
        raise ValueError(f"Selection matched {len(selected)} entries; one occupied cell is required.")
    row = selected.iloc[0]
    required = {"source_rows", "n_records", "activity_pchembl", "molecule_identity_key"}
    if not required.issubset(selected.columns):
        raise ValueError("This map lacks ChEMBL record provenance; use ChEMBL input mode first.")
    source_rows = [int(v) for v in str(row.source_rows).split(";")]
    if "source_row" not in records or records.source_row.duplicated().any():
        raise ValueError("Use the retained-record table from the same run, with unique source_row values.")
    evidence = records[records.source_row.isin(source_rows)].copy().sort_values("source_row")
    if len(evidence) != len(source_rows) or len(evidence) != int(row.n_records):
        raise ValueError("Record counts do not match this cell; check the retained-record file.")
    if "molecule_identity_key" not in evidence or not evidence.molecule_identity_key.eq(row.molecule_identity_key).all():
        raise ValueError("The records belong to a different molecular identity.")
    identity_records = records[records.molecule_identity_key.eq(row.molecule_identity_key)]
    if set(identity_records.source_row) != set(source_rows):
        raise ValueError("Source links omit records belonging to this molecular identity.")
    values = pd.to_numeric(evidence.activity_pchembl_raw, errors="raise")
    median = float(values.median())
    if not np.isfinite(values).all() or not np.isclose(median, float(row.activity_pchembl), rtol=0, atol=1e-10):
        raise ValueError("Recalculated activity does not match the map annotation.")
    columns = detect_chembl_columns(evidence)
    calculated, rules, verified, quality_verified = verify_annotations(row, evidence, values, parameters)
    summary = {
        "molecule_identity_key": str(row.molecule_identity_key),
        "grid_row": int(row.grid_row), "grid_col": int(row.grid_col),
        "coordinate_convention": "Zero-based row and column from the grid-coordinate CSV; rows increase upwards.",
        "activity_class": calculated.get("activity_class"), "median_pchembl": median,
        "minimum_pchembl": float(values.min()), "maximum_pchembl": float(values.max()),
        "n_records": len(evidence),
        "is_conflicted": bool(calculated["is_conflicted"]) if rules is not None else None,
        "has_class_boundary_crossing": bool(calculated["has_class_boundary_crossing"]) if rules is not None else None,
        "record_count_verified": True, "median_verified": True,
        "quality_annotations_verified": quality_verified, "annotation_rules": rules,
        "verified_fields": verified,
        "warnings": [] if quality_verified else ["Quality annotations are not fully verified. Legacy maps need the original curation parameters; no default thresholds were assumed."],
        "interpretation": "Median annotation summarizes retained records; it does not harmonize assay conditions.",
    }
    for name, field in (("n_activity_ids", "activity_id"), ("n_assays", "assay_id"), ("n_documents", "document_id")):
        summary[name] = int(evidence[columns[field]].nunique()) if columns[field] else None
    return selected.copy(), evidence, summary


def load_audit_context(coordinates_path, records_path, *, manifest_path=None, report_path=None):
    """Check saved file digests before interpreting them as outputs of the same run."""
    explicit_manifest = manifest_path is not None
    manifest_path = Path(manifest_path) if manifest_path else coordinates_path.with_name(
        coordinates_path.name.replace("_grid_coordinates.csv", "_run_manifest.json"))
    parameters, verified, warnings = None, False, []
    if manifest_path != coordinates_path and manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        parameters = manifest.get("source", {}).get("curation_parameters")
        expected = {"coordinates": manifest.get("output_digests", {}).get("coordinates"),
                    "retained_records": manifest.get("source", {}).get("curation_outputs", {}).get("retained_records")}
        actual = file_digests({"coordinates": coordinates_path, "retained_records": records_path})
        for key, digest in expected.items():
            if digest and actual[key]["sha256"] != digest["sha256"]:
                raise ValueError(f"{key} digest differs from the run manifest; do not mix or modify run outputs.")
        verified = all(expected.values())
    elif explicit_manifest:
        raise ValueError("The supplied run manifest does not exist.")
    if report_path is not None:
        report = json.loads(Path(report_path).read_text())
        supplied = report["parameters"]
        if parameters is not None and any(parameters.get(key) != supplied.get(key) for key in RULE_KEYS):
            raise ValueError("Curation report and run manifest use different annotation rules.")
        parameters = supplied
    if not verified:
        warnings.append("Run binding is not verified: original output digests are unavailable. Current file hashes alone do not establish a shared run.")
    return parameters, verified, warnings


def main(argv=None):
    parser = argparse.ArgumentParser(prog="chemgridmap inspect", description=__doc__)
    parser.add_argument("coordinates", type=Path, help="Grid-coordinate CSV from a ChEMBL-to-map run.")
    parser.add_argument("--records", required=True, type=Path, help="Retained-record CSV from the same run.")
    parser.add_argument("--manifest", type=Path, help="Run manifest; automatically found beside standard coordinate exports.")
    parser.add_argument("--curation-report", type=Path, help="Original curation report for legacy maps without saved annotation rules.")
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--molecule-id", help="Exact molecule, parent, or molecular identity key.")
    selection.add_argument("--cell", nargs=2, type=int, metavar=("ROW", "COL"), help="Zero-based row and column in the grid CSV.")
    parser.add_argument("--output-dir", type=Path, default=Path("record_audit"))
    args = parser.parse_args(argv)
    parameters, run_verified, warnings = load_audit_context(
        args.coordinates, args.records, manifest_path=args.manifest, report_path=args.curation_report)
    selected, records, summary = inspect_entry(
        pd.read_csv(args.coordinates), pd.read_csv(args.records),
        molecule_id=args.molecule_id, cell=args.cell, parameters=parameters,
    )
    summary["run_binding_verified"] = run_verified
    summary["warnings"].extend(warnings)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.output_dir / "selected_molecule.csv", index=False)
    records.to_csv(args.output_dir / "source_records.csv", index=False)
    render_svg(selected, args.output_dir / "selected_molecule.svg", tile_size=320, molecule_margin=12)
    summary["inputs"] = {name: {"file": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                         for name, path in (("coordinates", args.coordinates), ("retained_records", args.records))}
    (args.output_dir / "audit_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(f"Evidence written to {args.output_dir}")
    return 0
