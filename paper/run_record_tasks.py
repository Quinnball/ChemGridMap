"""Compare source-record retrieval with an explicit metadata-preserving pandas baseline.

Reuse saved layouts. This is an executable correctness comparison, not a user study.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".cache/matplotlib"))

import numpy as np
import pandas as pd

from chemgridmap.inspection import inspect_entry

ROOT = Path(__file__).resolve().parents[1]


def pandas_retrieve(coordinates, records, *, identity=None, cell=None):
    """Ordinary table filtering and joins, independent of ChemGridMap retrieval."""
    if identity is not None:
        chosen = coordinates.loc[coordinates.molecule_identity_key.eq(identity)]
    else:
        chosen = coordinates.loc[coordinates.grid_row.eq(cell[0]) & coordinates.grid_col.eq(cell[1])]
    if len(chosen) != 1:
        raise ValueError("Expected one molecular entry.")
    row = chosen.iloc[0]
    evidence = records.loc[records.molecule_identity_key.eq(row.molecule_identity_key)].sort_values("source_row")
    expected_rows = sorted(map(int, str(row.source_rows).split(";")))
    if evidence.source_row.tolist() != expected_rows or len(evidence) != row.n_records:
        raise ValueError("Source linkage does not match.")
    values = evidence.activity_pchembl_raw.astype(float)
    median = values.median()
    if not np.isclose(median, row.activity_pchembl, atol=1e-10, rtol=0):
        raise ValueError("Median does not match.")
    conflict = bool(values.max() - values.min() > 1 or (values.min() <= 6 and values.max() >= 7))
    return evidence, float(median), conflict


def compare_target(source, target):
    folder = source / target
    coordinates = pd.read_csv(folder / "coordinates_seed_42.csv")
    records = pd.read_csv(folder / f"{target}_retained_records.csv")
    parameters = json.loads((folder / f"{target}_curation_report.json").read_text())["parameters"]
    results = []
    for row in coordinates.itertuples():
        cell = (row.grid_row, row.grid_col)
        baseline_id, median, conflict = pandas_retrieve(coordinates, records, identity=row.molecule_identity_key)
        baseline_cell, _, _ = pandas_retrieve(coordinates, records, cell=cell)
        _, tool_id, summary = inspect_entry(coordinates, records, molecule_id=row.molecule_identity_key, parameters=parameters)
        _, tool_cell, _ = inspect_entry(coordinates, records, cell=cell, parameters=parameters)
        expected = sorted(records.loc[records.molecule_identity_key.eq(row.molecule_identity_key), "source_row"])
        passed = all(frame.source_row.tolist() == expected for frame in (baseline_id, baseline_cell, tool_id, tool_cell))
        assert passed and np.isclose(summary["median_pchembl"], median, rtol=0, atol=1e-10)
        assert summary["is_conflicted"] == conflict == bool(row.is_conflicted)
        results.append({"target": target.upper(), "identity": row.molecule_identity_key,
                        "n_records": len(expected), "is_conflicted": conflict,
                        "pandas_id_pass": passed, "pandas_cell_pass": passed,
                        "tool_id_pass": passed, "tool_cell_pass": passed,
                        "median_match": True, "conflict_match": True})
    return pd.DataFrame(results)


def reject_bad_inputs(coordinates, records):
    selected = coordinates.loc[coordinates.molecule_identity_key.eq("CHEMBL20")].copy()
    first_source = int(str(selected.iloc[0].source_rows).split(";")[0])
    wrong_median = coordinates.copy()
    wrong_median.loc[wrong_median.molecule_identity_key.eq("CHEMBL20"), "activity_pchembl"] += 0.5
    wrong_identity = records.copy()
    wrong_identity.loc[wrong_identity.source_row.eq(first_source), "molecule_identity_key"] = "WRONG_ID"
    cases = [("missing_source_row", coordinates, records[records.source_row.ne(first_source)]),
             ("changed_median", wrong_median, records),
             ("changed_identity", coordinates, wrong_identity)]
    outputs = []
    for name, coords, evidence in cases:
        for route in ("pandas_with_checks", "chemgridmap"):
            try:
                if route == "chemgridmap":
                    inspect_entry(coords, evidence, molecule_id="CHEMBL20")
                else:
                    pandas_retrieve(coords, evidence, identity="CHEMBL20")
            except ValueError as error:
                outputs.append({"case": name, "route": route, "rejected": True, "message": str(error)})
            else:
                raise AssertionError(f"{route} did not reject {name}")
    return outputs


def assay_context(source):
    all_records = pd.read_csv(source / "chembl205/acetazolamide_retained_records.csv")
    # Post hoc case follow-up: use the most represented document, not the largest effect.
    document_id = all_records.document_chembl_id.value_counts().index[0]
    assert document_id == "CHEMBL1141531"
    document_path = ROOT / "validation_data/assay_context" / f"{document_id}.json"
    metadata = json.loads(document_path.read_text())
    assert metadata["doi"] == "10.1016/j.bmcl.2008.02.008"
    subset = all_records[all_records.document_chembl_id.eq(document_id)].copy()
    names = {"DTT": "dithiothreitol", "2-ME": "mercaptoethanol",
             "TCEP": "phosphine", "Threitol": "threitol"}
    subset["additive"] = "Not specified"
    subset["concentration_mM"] = np.nan
    for index, row in subset.iterrows():
        description = str(row.assay_description).lower()
        for label, keyword in names.items():
            # Test DTT first because its full name contains 'threitol'.
            if keyword in description:
                match = re.search(r"(\d+(?:\.\d+)?)\s*mm", description)
                if match is None:
                    raise ValueError(f"Missing concentration: {description}")
                subset.loc[index, ["additive", "concentration_mM"]] = [label, float(match[1])]
                break
    specified = subset[subset.additive.ne("Not specified")]
    assert len(subset) == 25 and len(specified) == 24
    assert specified.groupby("additive").size().eq(6).all()
    assert not specified.duplicated(["additive", "concentration_mM"]).any()
    subset = subset.sort_values(["additive", "concentration_mM"])
    subset.to_csv(source / "acetazolamide_condition_series.csv", index=False)
    endpoints = specified.groupby("additive").apply(
        lambda g: pd.Series({"n": len(g),
                             "pchembl_low_concentration": g.loc[g.concentration_mM.idxmin(), "pchembl_value"],
                             "pchembl_high_concentration": g.loc[g.concentration_mM.idxmax(), "pchembl_value"]}),
    ).reset_index()
    endpoints.to_csv(source / "acetazolamide_condition_summary.csv", index=False)
    molecules = pd.read_csv(source / "chembl205/chembl205_molecule_level.csv")
    control = molecules[(molecules.n_records >= 3) & molecules.is_conflicted.eq(0)
                        & molecules.has_class_boundary_crossing.eq(0)].sort_values(
                            ["n_records", "molecule_identity_key"], ascending=[False, True]).iloc[0]
    records = pd.read_csv(source / "chembl205/chembl205_retained_records.csv")
    evidence = records[records.molecule_identity_key.eq(control.molecule_identity_key)]
    evidence.to_csv(source / "unflagged_control_records.csv", index=False)
    return {"selection": "Post hoc follow-up of the most represented document within the previously selected acetazolamide case.",
            "source_document": document_id, "doi": metadata["doi"], "pubmed_id": metadata["pubmed_id"],
            "metadata_url": f"https://www.ebi.ac.uk/chembl/api/data/document/{document_id}.json",
            "metadata_sha256": hashlib.sha256(document_path.read_bytes()).hexdigest(),
            "n_all_records": len(all_records), "n_document_records": len(subset),
            "n_explicit_condition_records": len(specified), "n_additives": 4,
            "reference_pchembl": float(subset.loc[subset.additive.eq("Not specified"), "pchembl_value"].iloc[0]),
            "evidence_scope": "Concentrations and values from archived ChEMBL assay descriptions and records; source DOI, metadata and abstract checked. Full article tables not independently transcribed.",
            "interpretation": "These concentration settings are not technical replicates under a single assay condition. No records were deleted, corrected or causally reclassified.",
            "control": {"selection": "Most records among groups with n>=3, no conflict flag, and no display-boundary crossing; ID breaks ties.",
                        "identity": str(control.molecule_identity_key), "n_records": len(evidence),
                        "n_assays": int(evidence.assay_chembl_id.nunique()),
                        "n_documents": int(evidence.document_chembl_id.nunique()),
                        "minimum": float(control.pchembl_min), "maximum": float(control.pchembl_max),
                        "median": float(control.activity_pchembl), "interpretation": "Unflagged, not proven assay comparability."}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=ROOT / "paper/output/revision_v6")
    args = parser.parse_args()
    source = args.source
    outputs = pd.concat([compare_target(source, target) for target in ("chembl205", "chembl204", "chembl240")])
    outputs.to_csv(source / "record_task_per_molecule.csv", index=False)
    summary = outputs.groupby("target").agg(n_molecules=("identity", "size"), n_records=("n_records", "sum"),
                n_conflicted=("is_conflicted", "sum"), pandas_id_pass=("pandas_id_pass", "sum"),
                pandas_cell_pass=("pandas_cell_pass", "sum"), tool_id_pass=("tool_id_pass", "sum"),
                tool_cell_pass=("tool_cell_pass", "sum"), median_match=("median_match", "sum"), conflict_match=("conflict_match", "sum"))
    summary.to_csv(source / "record_task_summary.csv")
    bad = reject_bad_inputs(pd.read_csv(source / "chembl205/coordinates_seed_42.csv"),
                           pd.read_csv(source / "chembl205/chembl205_retained_records.csv"))
    pd.DataFrame(bad).to_csv(source / "record_task_negative_controls.csv", index=False)
    result = {"passed": True, "timestamp_utc": datetime.now(timezone.utc).isoformat(),
              "design": "Same saved seed-42 layouts and full metadata for both routes; all mapped molecules tested, no outcome-based selection.",
              "n_molecules": len(outputs), "n_source_records": int(outputs.n_records.sum()),
              "n_id_cell_queries_per_route": 2 * len(outputs), "negative_controls_per_route": 3,
              "scope": "Task correctness and packaged functionality, not user efficiency, assay harmonization or superiority to mature visualization tools.",
              "case": assay_context(source)}
    (source / "record_tasks.json").write_text(json.dumps(result, indent=2) + "\n")
    print(summary.to_string())
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
