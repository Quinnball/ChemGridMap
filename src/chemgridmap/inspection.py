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


def inspect_entry(coordinates, records, *, molecule_id=None, cell=None):
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
    values = pd.to_numeric(evidence.activity_pchembl_raw, errors="raise")
    median = float(values.median())
    if not np.isfinite(values).all() or not np.isclose(median, float(row.activity_pchembl), rtol=0, atol=1e-10):
        raise ValueError("Recalculated activity does not match the map annotation.")
    columns = detect_chembl_columns(evidence)
    summary = {
        "molecule_identity_key": str(row.molecule_identity_key),
        "grid_row": int(row.grid_row), "grid_col": int(row.grid_col),
        "coordinate_convention": "Zero-based row and column from the grid-coordinate CSV; rows increase upwards.",
        "activity_class": str(row.activity_class), "median_pchembl": median,
        "minimum_pchembl": float(values.min()), "maximum_pchembl": float(values.max()),
        "n_records": len(evidence), "is_conflicted": bool(row.is_conflicted),
        "has_class_boundary_crossing": bool(row.has_class_boundary_crossing),
        "record_count_verified": True, "median_verified": True,
        "interpretation": "Median annotation summarizes retained records; it does not harmonize assay conditions.",
    }
    for name, field in (("n_activity_ids", "activity_id"), ("n_assays", "assay_id"), ("n_documents", "document_id")):
        summary[name] = int(evidence[columns[field]].nunique()) if columns[field] else None
    return selected.copy(), evidence, summary


def main(argv=None):
    parser = argparse.ArgumentParser(prog="chemgridmap inspect", description=__doc__)
    parser.add_argument("coordinates", type=Path, help="Grid-coordinate CSV from a ChEMBL-to-map run.")
    parser.add_argument("--records", required=True, type=Path, help="Retained-record CSV from the same run.")
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--molecule-id", help="Exact molecule, parent, or molecular identity key.")
    selection.add_argument("--cell", nargs=2, type=int, metavar=("ROW", "COL"), help="Zero-based row and column in the grid CSV.")
    parser.add_argument("--output-dir", type=Path, default=Path("record_audit"))
    args = parser.parse_args(argv)
    selected, records, summary = inspect_entry(
        pd.read_csv(args.coordinates), pd.read_csv(args.records),
        molecule_id=args.molecule_id, cell=args.cell,
    )
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
