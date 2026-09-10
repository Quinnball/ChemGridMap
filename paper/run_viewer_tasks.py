"""Controlled viewer-component comparison, not a user study or CLI benchmark.

Use the unmodified MolCompassView callbacks at the pinned source commit. An
adapter supplies saved coordinates to avoid confounding the comparison with its
parametric projection. Both routes receive the same map table and record file.
"""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import importlib.metadata
import io
import json
import os
from pathlib import Path
import platform
import sys

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".cache/matplotlib"))

import numpy as np
import pandas as pd

from chemgridmap.inspection import main as inspect_cli
from run_record_tasks import pandas_retrieve

ROOT = Path(__file__).resolve().parents[1]
COMMIT = "27c19943a256d4791984944293999b389eb42e17"
SOURCE_BLOBS = {
    "actions.py": "90ad9d31c5a7c7bb184e702b2af1dfc6819a449d",
    "__init__.py": "c57d6361af7a9cefe4d3f1bb102fe48535b1b020",
}


def choose_cases(coordinates):
    """Three prespecified strata; record count then exact identity breaks ties."""
    ordered = coordinates.sort_values(["n_records", "molecule_identity_key"], ascending=[False, True])
    groups = {
        "flagged_repeat": ordered[(ordered.n_records > 1) & ordered.is_conflicted.eq(1)],
        "unflagged_repeat": ordered[(ordered.n_records > 1) & ordered.is_conflicted.eq(0)],
        "singleton": ordered[ordered.n_records.eq(1)],
    }
    return [(name, group.iloc[0]) for name, group in groups.items() if len(group)]


def components(node):
    """Walk the actual Dash callback response, including molecule SVG and table."""
    if hasattr(node, "to_plotly_json"):
        node = node.to_plotly_json()
    if isinstance(node, dict):
        yield node
        for value in node.values():
            yield from components(value)
    elif isinstance(node, (list, tuple)):
        for value in node:
            yield from components(value)


def git_blob_sha(path):
    data = path.read_bytes()
    return hashlib.sha1(f"blob {len(data)}\0".encode() + data).hexdigest()


def check_source(folder):
    for name, expected in SOURCE_BLOBS.items():
        if git_blob_sha(folder / "molcompview" / name) != expected:
            raise ValueError(f"Comparator source differs from pinned commit: {name}")


def viewer(coordinates):
    from dash import Dash, html
    from molcompview import DatasetState, __smiles_name__, __x_name__, __y_name__
    from molcompview.actions import ColumnType, generate_figure_from_data, init_callbacks

    # Retain every input column. The extra aliases are internal viewer conventions.
    data = coordinates.copy()
    data[__smiles_name__] = data.representation_smiles
    data[__x_name__], data[__y_name__] = data.projection_x, data.projection_y
    kinds = {"activity_pchembl": ColumnType.NUMERICAL}
    figure = generate_figure_from_data(data, "activity_pchembl", kinds,
                                     [data.activity_pchembl.min(), data.activity_pchembl.max()])
    app = Dash(__name__)
    app.layout = html.Div()
    init_callbacks(app, data, kinds, DatasetState.PROPERTY)
    callback = next(v["callback"].__wrapped__ for k, v in app.callback_map.items()
                    if "molcompass-graph-tooltip.children" in k)
    return figure, callback


def run(source, output, comparator):
    check_source(comparator)
    sys.path.insert(0, str(comparator))
    output.mkdir(parents=True, exist_ok=True)
    protocol = {
        "comparator": "MolCompassView visualization callbacks (not its complete CLI)",
        "comparator_commit": COMMIT,
        "comparator_url": f"https://github.com/sergsb/molcompview/tree/{COMMIT}",
        "source_blobs": SOURCE_BLOBS,
        "design": "Same saved seed-42 map tables and retained records; no projection recomputation or viewer source edits.",
        "selection": "Per target: greatest record count among flagged repeats, unflagged repeats, and singletons; exact identity breaks ties.",
        "tasks": ["Recover structure and supplied annotation from a selected point/cell",
                  "Retrieve all contributing source rows and recalculate median",
                  "Export a structure, source-record CSV and checked audit summary"],
        "scope": "Programmatic component task trace. No human task timing, error rate, click-count ranking, or exhaustive feature comparison.",
        "baseline_extension": "The same explicit pandas join/check baseline is supplied for the viewer-to-record bridge; both workflows may be scripted.",
        "source_file_hashes": {},
    }
    for target in ("chembl205", "chembl204", "chembl240"):
        for name in ("coordinates_seed_42.csv", f"{target}_retained_records.csv"):
            path = source / target / name
            protocol["source_file_hashes"][f"{target}/{name}"] = hashlib.sha256(path.read_bytes()).hexdigest()
    # Save the deterministic protocol before executing the selected tasks.
    (output / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n", encoding="utf-8")
    results = []
    for target in ("chembl205", "chembl204", "chembl240"):
        coordinate_path = source / target / "coordinates_seed_42.csv"
        record_path = source / target / f"{target}_retained_records.csv"
        coordinates, records = pd.read_csv(coordinate_path), pd.read_csv(record_path)
        figure, hover = viewer(coordinates)
        for stratum, row in choose_cases(coordinates):
            folder = output / target / stratum
            folder.mkdir(parents=True, exist_ok=True)
            _, _, card = hover({"points": [{"x": row.projection_x, "y": row.projection_y,
                                            "bbox": {"x0": 0, "x1": 1, "y0": 0, "y1": 1}}]}, [])
            nodes = list(components(card))
            images = [n["props"]["src"] for n in nodes if n.get("type") == "Img"]
            properties = {}
            for node in nodes:
                if node.get("type") == "Tr":
                    cells = node["props"]["children"]
                    if len(cells) == 2 and all(c.to_plotly_json()["type"] == "Td" for c in cells):
                        properties[str(cells[0].children)] = str(cells[1].children)
            assert len(images) == 1 and images[0].startswith("data:image/svg+xml;base64,")
            assert properties["molecule_identity_key"] == row.molecule_identity_key
            assert np.isclose(float(properties["activity_pchembl"]), row.activity_pchembl, atol=0.0005, rtol=0)
            assert properties["source_rows"] == str(row.source_rows)
            (folder / "viewer_properties.json").write_text(json.dumps(properties, indent=2) + "\n", encoding="utf-8")
            import base64
            (folder / "viewer_molecule.svg").write_bytes(base64.b64decode(images[0].split(",", 1)[1]))
            baseline, median, conflict = pandas_retrieve(coordinates, records, identity=row.molecule_identity_key)
            baseline.to_csv(folder / "viewer_plus_pandas_records.csv", index=False)
            with contextlib.redirect_stdout(io.StringIO()):
                inspect_cli([str(coordinate_path), "--records", str(record_path), "--cell",
                             str(row.grid_row), str(row.grid_col), "--output-dir", str(folder / "chemgridmap")])
            tool = pd.read_csv(folder / "chemgridmap/source_records.csv")
            summary = json.loads((folder / "chemgridmap/audit_summary.json").read_text(encoding="utf-8"))
            pd.testing.assert_frame_equal(tool, baseline.reset_index(drop=True), check_dtype=False)
            assert np.isclose(summary["median_pchembl"], median, rtol=0, atol=1e-10)
            assert summary["is_conflicted"] == conflict
            results.append({"target": target.upper(), "stratum": stratum, "identity": row.molecule_identity_key,
                            "n_records": len(tool), "median": median,
                            "viewer_structure": True, "viewer_annotation": True, "viewer_source_links": True,
                            "viewer_plus_pandas_record_match": True, "chemgridmap_cli_record_match": True,
                            "median_match": True, "flag_match": True})
        # A saved interactive plot verifies all input points, not only case crops.
        figure.write_html(output / f"{target}_viewer_scatter.html", include_plotlyjs="cdn")
    frame = pd.DataFrame(results)
    frame.to_csv(output / "task_results.csv", index=False)
    result = {"passed": True, "n_cases": len(frame), "n_records": int(frame.n_records.sum()),
              "python": platform.python_version(), "platform": platform.platform(),
              "packages": {name: importlib.metadata.version(name) for name in
                           ("dash", "dash-bootstrap-components", "plotly", "pandas", "rdkit", "chemgridmap")},
              "interpretation": "Native viewer callbacks recover structures, annotations and supplied row links. An explicit external join/check bridge recovers identical records. ChemGridMap supplies this bridge and evidence export through its inspection command.",
              "not_established": "Usability superiority, human time saved, competitor inability, stock MolCompass CLI behavior, or novel biological findings."}
    (output / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(frame[["target", "stratum", "n_records", "median"]].to_string(index=False))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=ROOT / "paper/output/revision_v6")
    parser.add_argument("--output", type=Path, default=ROOT / "paper/output/revision_v7/viewer_tasks")
    parser.add_argument("--comparator", type=Path, required=True, help="Extracted pinned MolCompassView source directory.")
    args = parser.parse_args()
    run(args.source, args.output, args.comparator)
