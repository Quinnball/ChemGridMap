"""Validate record retrieval and fixed-molecule neighborhoods using saved layouts."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import itertools
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".cache/matplotlib"))

import numpy as np
import pandas as pd

from chemgridmap.core import _tie_inclusive_neighbor_indices
from chemgridmap.inspection import inspect_entry
from chemgridmap.plotting import render_svg
from chemgridmap.representations import morgan_fingerprints

ROOT = Path(__file__).resolve().parents[1]


def dump(path, data):
    path.write_text(json.dumps(data, indent=2, default=str) + "\n", encoding="utf-8")


def agreement(reference, candidate):
    pairs = [(set(a), set(b)) for a, b in zip(reference, candidate)]
    return {"recall": float(np.mean([len(a & b) / len(a) for a, b in pairs])),
            "precision": float(np.mean([len(a & b) / len(b) for a, b in pairs]))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=ROOT / "paper/output/revision_v4")
    parser.add_argument("--output", type=Path, default=ROOT / "paper/output/revision_v5")
    args = parser.parse_args()
    if args.source.resolve() == args.output.resolve():
        raise ValueError("Keep the earlier validated layouts unchanged.")
    shutil.copytree(args.source, args.output, dirs_exist_ok=True, ignore=shutil.ignore_patterns("figures"))
    out = args.output
    full_path = []
    for target in ("chembl205", "chembl204", "chembl240"):
        case = out / target
        data = pd.read_csv(case / "coordinates_seed_42.csv")
        fp = morgan_fingerprints(data.canonical_smiles.tolist())
        upstream = _tie_inclusive_neighbor_indices(fp, 10, metric="jaccard")
        saved = json.loads((case / "upstream_fidelity.json").read_text())
        for destination, columns in (("projection", ["projection_x", "projection_y"]),
                                     ("grid", ["grid_x", "grid_y"])):
            downstream = _tie_inclusive_neighbor_indices(data[columns].to_numpy(), 10)
            full_path.append({"target": target.upper(), "destination": destination,
                              "seed": 42, "n": len(data), "k": 10,
                              **agreement(upstream, downstream),
                              "trustworthiness": saved[f"fingerprint_to_{destination}_trustworthiness"]})
    pd.DataFrame(full_path).to_csv(out / "full_path_fidelity.csv", index=False)

    case = out / "chembl205"
    main_coords = pd.read_csv(case / "coordinates_seed_42.csv")
    records = pd.read_csv(case / "chembl205_retained_records.csv")
    chosen, evidence, summary = inspect_entry(main_coords, records, molecule_id="CHEMBL20")
    cell = (int(chosen.iloc[0].grid_row), int(chosen.iloc[0].grid_col))
    _, same_evidence, _ = inspect_entry(main_coords, records, cell=cell)
    assert evidence.source_row.tolist() == same_evidence.source_row.tolist()
    assert summary["n_records"] == 54 and np.isclose(summary["median_pchembl"], 7.185)
    for selector, selection in (("id", ["--molecule-id", "CHEMBL20"]),
                                ("cell", ["--cell", str(cell[0]), str(cell[1])])):
        command = [sys.executable, "-m", "chemgridmap", "inspect", str(case / "coordinates_seed_42.csv"),
                   "--records", str(case / "chembl205_retained_records.csv"), *selection,
                   "--output-dir", str(out / f"inspection_{selector}")]
        result = subprocess.run(command, text=True, capture_output=True, check=True)
        (out / f"inspection_{selector}" / "command.log").write_text(result.stdout + result.stderr)
    render_svg(main_coords, case / "molecular_grid.svg")

    windows = pd.read_csv(out / "window_members.csv")
    anchors = windows[windows.window.isin(["Active-rich", "Mixed"])][["window", "molecule_identity_key"]].copy()
    anchors = pd.concat([anchors, pd.DataFrame([{"window": "Audited molecule", "molecule_identity_key": summary["molecule_identity_key"]}])], ignore_index=True)
    assert len(anchors) == 19 and anchors.molecule_identity_key.nunique() == 19
    neighbors = {}
    for seed in (7, 21, 42, 84, 126):
        layout = pd.read_csv(case / f"coordinates_seed_{seed}.csv").set_index("molecule_identity_key")
        for space, cols in (("projection", ["projection_x", "projection_y"]), ("grid", ["grid_x", "grid_y"])):
            sets = _tie_inclusive_neighbor_indices(layout[cols].to_numpy(), 10)
            neighbors[seed, space] = {key: set(layout.index[sets[layout.index.get_loc(key)]])
                                     for key in anchors.molecule_identity_key}
    stability = []
    for anchor in anchors.itertuples():
        for a, b in itertools.combinations((7, 21, 42, 84, 126), 2):
            row = {"window": anchor.window, "molecule_identity_key": anchor.molecule_identity_key,
                   "seed_a": a, "seed_b": b}
            for space in ("projection", "grid"):
                aa, bb = neighbors[a, space][anchor.molecule_identity_key], neighbors[b, space][anchor.molecule_identity_key]
                row[f"{space}_neighbor_jaccard"] = len(aa & bb) / len(aa | bb)
            stability.append(row)
    stability = pd.DataFrame(stability)
    stability.to_csv(out / "fixed_molecule_seed_stability.csv", index=False)
    stability.groupby("window")[["projection_neighbor_jaccard", "grid_neighbor_jaccard"]].mean().to_csv(out / "fixed_molecule_seed_summary.csv")
    pd.DataFrame([
        {"output": "Median-only color annotation", "displayed_information": "Active; median pChEMBL 7.185",
         "interpretation": "Insufficient to assess measurement consistency"},
        {"output": "Linked source-record export", "displayed_information": "54 records; range 4.73-8.47; 54 assays; 27 documents",
         "interpretation": "Active median with flagged measurement heterogeneity; assay review needed"},
    ]).to_csv(out / "inspection_task_comparison.csv", index=False)
    baseline = json.loads((out / "integrity_audit.json").read_text())
    dump(out / "application_validation.json", {
        "passed": True, "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "baseline_audit_passed": baseline["passed"], "saved_layout_source": str(args.source),
        "layout_policy": "Reuse all previous layouts without outcome-dependent reselection or rerunning UMAP.",
        "retrieval": summary, "id_and_cell_retrieve_identical_source_rows": True,
        "fixed_anchor_count": len(anchors), "anchor_seed_pairs": len(stability),
        "fixed_anchor_policy": "Track molecular identity across seeds, not the same grid coordinates.",
        "source_coordinate_sha256": hashlib.sha256((case / "coordinates_seed_42.csv").read_bytes()).hexdigest(),
        "scope": "Functional information-retrieval validation; not a user-efficiency study or competitor benchmark.",
    })
    print(pd.DataFrame(full_path).to_string(index=False))
    print(stability.groupby("window")[["projection_neighbor_jaccard", "grid_neighbor_jaccard"]].mean())
    print("Application validation passed.")


if __name__ == "__main__":
    main()
