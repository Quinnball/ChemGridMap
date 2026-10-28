"""Reproducible neighborhood examples and CA II case positions for T2a."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from chemgridmap.comparison import structural_groups
from chemgridmap.reliability import grid_neighbor_sets, representation_neighbor_sets
from chemgridmap.representations import morgan_fingerprints


ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "paper/output/revision_v6"
V2_SOURCE = ROOT.parent / "project1/JCIM_molecular_figure_revision_20261014_v2/figures/source_data"
MARKERS = {"query": "Q", "shared": "S", "grid_only": "G", "representation_only": "F"}


def neighbor_case_table(data, fingerprints, cells, *, query_id="CHEMBL81688", layout, k=10):
    """Return the query plus the tie-inclusive union of its two neighbor sets."""
    data = data.reset_index(drop=True)
    matches = np.flatnonzero(data.molecule_identity_key.astype(str).eq(query_id))
    if len(matches) != 1:
        raise ValueError("The case query must identify exactly one mapped entry.")
    position = int(matches[0])
    cells = np.asarray(cells, dtype=int)
    fingerprint_set = set(representation_neighbor_sets(fingerprints, k=k)[position])
    grid_set = set(grid_neighbor_sets(cells, k=k)[position])
    shared = fingerprint_set & grid_set
    distance = cdist(np.asarray(fingerprints)[[position]], fingerprints, metric="jaccard")[0]
    grid_distance = np.linalg.norm(cells - cells[position], axis=1)
    selected = [position] + sorted(fingerprint_set | grid_set,
                                  key=lambda i: (distance[i], str(data.molecule_identity_key.iloc[i])))
    fields = [name for name in ("molecule_identity_key", "molecule_id", "canonical_smiles") if name in data]
    table = data.iloc[selected][fields].copy().reset_index(drop=True)
    table.insert(0, "layout", layout)
    table.insert(0, "seed", 42)
    table.insert(0, "target", "chembl205")
    table["query_id"] = query_id
    table["is_query"] = [i == position for i in selected]
    table["grid_col"] = cells[selected, 0]
    table["grid_row"] = cells[selected, 1]
    table["representation_neighbor"] = [i in fingerprint_set for i in selected]
    table["grid_neighbor"] = [i in grid_set for i in selected]
    table["shared_neighbor"] = [i in shared for i in selected]
    table["neighbor_relation"] = ["query" if i == position else "shared" if i in shared else
                                  "representation_only" if i in fingerprint_set else "grid_only"
                                  for i in selected]
    table["figure_marker"] = table.neighbor_relation.map(MARKERS)
    table["representation_distance_tanimoto"] = distance[selected]
    table["grid_distance_cells"] = grid_distance[selected]
    summary = {"target": "chembl205", "seed": 42, "query_id": query_id, "layout": layout,
               "requested_k": int(k), "fingerprint_neighbor_count": len(fingerprint_set),
               "grid_neighbor_count": len(grid_set), "shared_neighbor_count": len(shared),
               "fingerprint_only_count": len(fingerprint_set - grid_set),
               "grid_only_count": len(grid_set - fingerprint_set),
               "neighbor_recall": len(shared) / len(fingerprint_set),
               "neighbor_precision": len(shared) / len(grid_set),
               "query_grid_col": int(cells[position, 0]), "query_grid_row": int(cells[position, 1]),
               "fingerprint_neighbors": sorted(data.molecule_identity_key.iloc[list(fingerprint_set)].tolist()),
               "grid_neighbors": sorted(data.molecule_identity_key.iloc[list(grid_set)].tolist()),
               "shared_neighbors": sorted(data.molecule_identity_key.iloc[list(shared)].tolist())}
    return table, summary


def compare_v2_neighbors(current, previous):
    """Verify the frozen v2 example before interpreting the refined example."""
    actual = current.loc[~current.is_query].set_index("molecule_identity_key").sort_index()
    expected = previous.set_index("molecule_identity_key").sort_index()
    same_members = actual.index.equals(expected.index)
    fields = ("neighbor_relation", "representation_neighbor", "grid_neighbor", "grid_col", "grid_row")
    equal_fields = {name: bool(same_members and actual[name].eq(expected[name]).all()) for name in fields}
    same_distance = bool(same_members and np.allclose(actual.representation_distance_tanimoto,
                         expected.representation_distance, atol=1e-12, rtol=0))
    return {"passed": same_members and all(equal_fields.values()) and same_distance,
            "same_neighbor_members": same_members, "same_fields": equal_fields,
            "same_tanimoto_distances": same_distance, "neighbor_union_count": len(actual),
            "v2_marker_definitions": MARKERS}


def annotated_case_positions(data, layout, previous_g002):
    """Calculate coverage-rule groups and verify that G002 retains its v2 members."""
    grouped, groups, protocol = structural_groups(data, core_rule="coverage")
    fields = ["molecule_identity_key", "structure_group", "scaffold_smiles"]
    positions = grouped[fields].merge(
        layout[["molecule_identity_key", "refined_grid_col", "refined_grid_row"]],
        on="molecule_identity_key", how="left", validate="one_to_one")
    if positions[["refined_grid_col", "refined_grid_row"]].isna().any().any():
        raise ValueError("Every CA II mapped entry must have a refined seed-42 position.")
    positions.insert(0, "seed", 42)
    positions.insert(0, "target", "chembl205")
    positions["is_G002"] = positions.structure_group.eq("G002")
    positions["is_CHEMBL20"] = positions.molecule_identity_key.eq("CHEMBL20")
    actual = set(positions.loc[positions.is_G002, "molecule_identity_key"])
    expected = set(previous_g002.comparison_id)
    group = next(row for row in groups if row["structure_group"] == "G002")
    check = {"passed": actual == expected and int(positions.is_CHEMBL20.sum()) == 1,
             "group_id": "G002", "group_member_count": len(actual),
             "members_match_v2": actual == expected, "missing_from_current": sorted(expected - actual),
             "new_in_current": sorted(actual - expected), "core_smiles": group["core_smiles"],
             "core_coverage_count": group["core_coverage_count"],
             "core_exact_scaffold_count": group["core_exact_scaffold_count"],
             "representative_id": group["representative_id"],
             "CHEMBL20_count": int(positions.is_CHEMBL20.sum())}
    return positions, check, protocol


def run_cases(output_root, *, input_root=FROZEN, layout_root=None, v2_root=V2_SOURCE):
    """Write case tables only; return provenance for the shared manifest writer."""
    output_root, input_root, v2_root = Path(output_root), Path(input_root), Path(v2_root)
    layout_root = Path(layout_root) if layout_root is not None else output_root / "layouts"
    inputs = [input_root / "chembl205/map_molecules.csv",
              layout_root / "chembl205/coordinates_seed_42.csv",
              v2_root / "Figure_3_example_neighbors.csv",
              v2_root / "Figure_1_all_28_member_categories.csv"]
    data, layout, old_neighbors, old_g002 = [pd.read_csv(path, low_memory=False) for path in inputs]
    data = data.sort_values("molecule_identity_key").reset_index(drop=True)
    layout = layout.set_index("molecule_identity_key").loc[data.molecule_identity_key].reset_index()
    fp = morgan_fingerprints(data.canonical_smiles.tolist())
    tables, summaries = [], []
    for name, columns in (("minimum_displacement", ["grid_col", "grid_row"]),
                          ("neighbor_refined", ["refined_grid_col", "refined_grid_row"])):
        table, summary = neighbor_case_table(data, fp, layout[columns].to_numpy(), layout=name)
        tables.append(table)
        summaries.append(summary)
    v2_check = compare_v2_neighbors(tables[0], old_neighbors)
    positions, group_check, grouping = annotated_case_positions(data, layout, old_g002)
    if not v2_check["passed"] or not group_check["passed"]:
        raise ValueError(f"Case source mismatch: v2={v2_check}; G002={group_check}")
    folder = output_root / "cases"
    folder.mkdir(parents=True, exist_ok=True)
    pd.concat(tables, ignore_index=True).to_csv(folder / "CHEMBL81688_neighbors_seed42.csv", index=False)
    positions.to_csv(folder / "chembl205_annotations_seed42.csv", index=False)
    positions.loc[positions.is_G002 | positions.is_CHEMBL20].to_csv(
        folder / "G002_CHEMBL20_positions_seed42.csv", index=False)
    result = {"CHEMBL81688": summaries, "v2_neighbor_check": v2_check, "G002_check": group_check}
    (folder / "case_summary.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return {"inputs": inputs, "parameters": {"seed": 42, "target": "chembl205", "k": 10,
            "tie_policy": "distance <= kth_distance * (1 + 1e-9) + 1e-12; self excluded",
            "markers": MARKERS, "grouping": grouping}, "summary": result}
