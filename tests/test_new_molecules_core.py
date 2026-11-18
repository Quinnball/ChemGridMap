"""Frozen-reference insertion: evidence, exact placement, and isolated artifacts."""
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import xml.etree.ElementTree as ET
import zipfile

import numpy as np
import pandas as pd
import pytest

from chemgridmap import build_grid_map, insert_new_molecules
from chemgridmap.queries import (QUERY_CHECK_COLUMNS, QUERY_NEIGHBOR_COLUMNS,
    place_new_molecules, summarize_neighbor_evidence)


@pytest.fixture
def reference(tmp_path):
    return build_grid_map(pd.DataFrame({"molecule_id": ["B", "A", "C", "D", "E", "F"],
        "smiles": ["CCC", "CCO", "CCN", "CCCC", "CCCl", "CCBr"], "value": [5., 6., 7., 8., 5.5, 6.5]}),
        value_col="value", projection="pca", grouping="none", grid_refinement="none",
        output_dir=tmp_path, output_formats=["svg"], render_detail="overview")


def test_six_neighbor_categories_threshold_support_and_measured_comparison():
    queries = pd.DataFrame({"query_id": list("abcdef"), "measurement_status": ["measured"] * 6,
                            "activity_class": ["active", "active", "active", "inactive", "inactive", "active"]})
    neighbors = pd.DataFrame([
        ("a", .59, "active", "supported"), ("b", .8, "active", "uncertain"),
        ("c", .8, "active", "supported"), ("d", .8, "medium", "supported"),
        ("e", .8, "inactive", "supported"), ("f", .8, "active", "supported"),
        ("f", .8, "inactive", "supported"),
    ], columns=["query_id", "tanimoto", "activity_class", "class_status"])
    result = summarize_neighbor_evidence(queries, neighbors)
    assert result.neighbor_evidence.tolist() == ["none", "insufficient", "active", "medium", "inactive", "mixed"]
    assert result.neighbor_comparison.tolist() == ["", "", "agrees", "disagrees", "agrees", ""]
    assert result.similar_reference_count.tolist() == [0, 1, 1, 1, 1, 2]
    assert result.supported_reference_count.tolist() == [0, 0, 1, 1, 1, 2]
    changed = summarize_neighbor_evidence(queries, neighbors, similarity=.59, support_assessed=False)
    assert changed.neighbor_evidence.tolist()[:2] == ["active", "active"]
    assert not changed.reference_support_assessed.any()
    assert changed.supported_reference_count.eq(0).all()
    assert changed.supported_reference_classes.eq("").all()
    assert changed.evidence_reference_count.tolist() == [1, 1, 1, 1, 1, 2]


def test_placement_exhaustive_distance_order_ties_and_immutable_reference():
    reference = pd.DataFrame({"molecule_id": ["B", "A"], "canonical_smiles": ["CCC", "CCO"],
                              "grid_row": [1, 2], "grid_col": [1, 2]})
    original = reference.copy(deep=True)
    queries = pd.DataFrame({"query_id": ["z", "b", "a"], "nearest_measured_id": ["A", "B", "B"],
                           "nearest_measured_similarity": [.5, .8, .8]})
    result = place_new_molecules(reference, queries, grid_rows=3, grid_cols=3)
    repeated = place_new_molecules(reference, queries, grid_rows=3, grid_cols=3)
    pd.testing.assert_frame_equal(result, repeated)
    pd.testing.assert_frame_equal(reference, original)
    occupied = set(zip(reference.grid_row, reference.grid_col))
    available = {(r, c) for r in range(3) for c in range(3)} - occupied
    lookup = reference.set_index("molecule_id")
    for identity in ["a", "b", "z"]:
        row = result[result.query_id.eq(identity)].iloc[0]
        target = lookup.loc[row.nearest_measured_id]
        expected = sorted(available, key=lambda p: ((p[0]-target.grid_row)**2+(p[1]-target.grid_col)**2,p))[0]
        assert (row.grid_row, row.grid_col) == expected
        assert row.placement_distance == pytest.approx(np.hypot(expected[0]-target.grid_row, expected[1]-target.grid_col))
        available.remove(expected)
    assert tuple(result[result.query_id.eq("a")][["grid_row", "grid_col"]].iloc[0]) == (0, 1)
    with pytest.raises(ValueError, match="empty cells"):
        place_new_molecules(reference, pd.concat([queries] * 3, ignore_index=True), grid_rows=3, grid_cols=3)


def test_outputs_raw_input_suffix_manifest_zip_and_reference_bytes(reference, tmp_path):
    before = {p.name: p.read_bytes() for p in tmp_path.iterdir() if p.is_file()}
    raw = b"id\tsmiles\r\nnew\tCCF\r\nmatch\tCCC\r\ndup\tFCC\r\nbad\tinvalid\r\n"
    frame = pd.read_csv(__import__("io").BytesIO(raw), sep="\t")
    one = insert_new_molecules(reference, frame, id_col="id", output_dir=tmp_path,
                               input_bytes=raw, input_filename="input.tsv")
    two = insert_new_molecules(reference, frame, id_col="id", output_dir=tmp_path,
                               input_bytes=raw, input_filename="input.tsv")
    assert {p.name: p.read_bytes() for p in tmp_path.iterdir() if p.is_file()} == before
    assert Path(one.output_files["new_molecules_input"]).parent.name == "input"
    assert Path(two.output_files["new_molecules_input"]).parent.name == "input-2"
    assert Path(one.output_files["new_molecules_input"]).read_bytes() == raw
    assert list(one.check) == QUERY_CHECK_COLUMNS
    assert list(one.neighbors) == QUERY_NEIGHBOR_COLUMNS
    assert one.check.status.tolist() == ["mapped", "in_reference", "duplicate_query", "invalid_smiles"]
    pd.testing.assert_frame_equal(one.data, two.data)
    manifest = json.loads(Path(one.output_files["new_molecules_manifest"]).read_text(encoding="utf-8"))
    assert manifest["input_sha256"] == hashlib.sha256(raw).hexdigest()
    assert manifest["reference_manifest_sha256"] == hashlib.sha256(Path(reference.output_files["run_manifest"]).read_bytes()).hexdigest()
    assert len(list(Path(one.output_files["new_molecules_input"]).parent.iterdir())) == 7
    with zipfile.ZipFile(one.output_files["new_molecules_results_zip"]) as bundle:
        assert len(bundle.namelist()) == 6
        for filename in bundle.namelist():
            assert bundle.read(filename) == Path(one.output_files["new_molecules_input"]).with_name(filename).read_bytes()
    json.dumps(one.report, allow_nan=False)
    assert one.report["insertion_seconds"] >= manifest["timings_seconds"]["through_svg_and_tables"]


def test_measured_units_aggregation_matches_and_exclusions(reference):
    query = pd.DataFrame({"id": ["new"]*4 + ["control"]*2,
        "smiles": ["CCF"]*4 + ["CCC"]*2, "value": [100, .1, .1, 30, 1000, 100],
        "unit": ["nM", "uM", "µM", "nM", "nM", "nM"], "relation": ["=", "=", "=", "<", "=", "="],
        "batch": ["R1", "R2", "R3", "R3", "R1", "R2"], "label": ["candidate"]*4+["control"]*2})
    result = insert_new_molecules(reference, query, id_col="id", value_col="value", unit_column="unit",
                                  relation_column="relation", document_column="batch", label_column="label")
    row = result.data.iloc[0]
    assert row.query_id == "new" and row.activity_pchembl == pytest.approx(7.)
    assert row.n_records == 3 and row.pchembl_range == pytest.approx(0)
    assert row.activity_class == "active" and row.label == "candidate"
    assert result.report["excluded_reasons"]["non_exact_relation"] == 1
    assert len(result.report["excluded_records"]) == 1
    assert result.report["excluded_record_columns"] == {
        "compound_id": "id", "evidence_unit": "batch", "value": "value",
        "unit": "unit", "relation": "relation"}
    assert result.report["excluded_record_defaults"] == {"unit": "p", "relation": "="}
    match = result.report["reference_matches"][0]
    assert match["reference_id"] == "B"
    assert match["median_difference"] == pytest.approx(1.5)
    assert result.report["reference_match_summary"]["median"] == pytest.approx(1.5)
    json.dumps(result.report, allow_nan=False)


@pytest.mark.parametrize("excluded_count", [0, 57])
def test_excluded_export_keeps_all_original_columns_rows_and_values(reference, tmp_path, excluded_count):
    query = pd.DataFrame({"id": [" new "] * (excluded_count + 1),
        "smiles": ["CCF"] * (excluded_count + 1),
        "amount": ["00100"] + ["<00030"] * excluded_count,
        "unit": ["nM"] * (excluded_count + 1),
        "original_row": [f"original-{i}" for i in range(excluded_count + 1)],
        "custom note": ["保留,原文"] * (excluded_count + 1)})
    original = query.copy(deep=True)
    result = insert_new_molecules(reference, query, id_col="id", value_col="amount",
                                  unit_column="unit", output_dir=tmp_path)
    pd.testing.assert_frame_equal(query, original)
    path = Path(result.output_files["new_molecules_excluded_records"])
    exported = pd.read_csv(path, dtype=str, keep_default_na=False)
    expected = original.iloc[1:].reset_index(drop=True)
    expected["exclusion_reason"] = "invalid_activity_value"
    pd.testing.assert_frame_equal(exported, expected)
    assert list(exported) == list(original) + ["exclusion_reason"]
    assert len(exported) == result.report["excluded_activity_records"] == excluded_count
    assert result.report["excluded_record_columns"]["relation"] is None
    with zipfile.ZipFile(result.output_files["new_molecules_results_zip"]) as archive:
        assert archive.read(path.name) == path.read_bytes()
    manifest = json.loads(Path(result.output_files["new_molecules_manifest"]).read_text(encoding="utf-8"))
    assert any(item["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
               for item in manifest["output_digests"].values())


@pytest.mark.parametrize("detail", ["overview", "full"])
def test_inserted_svg_html_parse_and_render_preserve_reference_and_new_cells(reference, tmp_path, detail):
    # 实际 HTML 解析检出 ns0:g；随后渲染该 HTML 解析结果，检查每个格子的像素范围。
    etree = pytest.importorskip("lxml.etree")
    html_parser = pytest.importorskip("lxml.html")
    pymupdf = pytest.importorskip("pymupdf")
    from chemgridmap.plotting import render_grid_map

    reference.output_files["svg"] = render_grid_map(reference.data, tmp_path / "reference_for_html",
        detail=detail, formats=["svg"])["svg"]
    result = insert_new_molecules(reference, pd.DataFrame({"id": ["NEW"], "smiles": ["CCF"]}),
                                  id_col="id", output_dir=tmp_path)
    raw = Path(result.output_files["map_with_new_molecules"]).read_bytes()
    assert b"ns0:" not in raw
    container = html_parser.fragment_fromstring(raw.decode(), create_parent="div")
    svg = container[0]
    assert svg.tag == "svg"
    cells = svg.xpath(".//*[@data-molecule-id]")
    assert len(cells) == len(reference.data) + len(result.data)
    for cell in cells:
        assert cell.tag == "g"
        assert cell.find("rect") is not None
        assert all(":" not in node.tag and "}" not in node.tag for node in cell.iter()
                   if isinstance(node.tag, str))
    parsed = etree.tostring(svg)
    # HTML 序列化仍保留默认 xmlns；独立打开时所有图形元素属于 SVG namespace。
    for document in (raw, parsed):
        xml = ET.fromstring(document)
        assert all(node.tag.startswith("{http://www.w3.org/2000/svg}") for node in xml.iter())

    def render(document):
        with pymupdf.open(stream=document, filetype="svg") as drawing:
            pixels = drawing[0].get_pixmap(alpha=False)
            return np.frombuffer(pixels.samples, dtype=np.uint8).reshape(pixels.height, pixels.width, 3)

    standalone, inline = render(raw), render(parsed)
    np.testing.assert_array_equal(inline, standalone)
    rows = reference.insertion_parameters["grid_rows"]
    for frame in (reference.data, result.data):
        for molecule in frame.itertuples():
            x, y = int(molecule.grid_col) * 100, (rows - 1 - int(molecule.grid_row)) * 100
            content = np.argwhere(np.any(inline[y+2:y+98, x+2:x+98] != 255, axis=2))
            assert len(content) > 50
            height, width = np.ptp(content, axis=0) + 1
            assert width > 10 and height > 10


@pytest.mark.parametrize("representation", ["descriptors", "embedding", "user_coordinates"])
def test_non_morgan_reference_disabled(reference, representation):
    reference.insertion_parameters["representation"] = representation
    with pytest.raises(ValueError, match="Morgan"):
        insert_new_molecules(reference, pd.DataFrame({"smiles": ["CCF"]}))


def test_limit_is_checked_before_fingerprints(reference, monkeypatch):
    import chemgridmap.queries as module
    original = module.prepare_query_molecules
    def excessive(*args, **kwargs):
        data, check, report = original(*args, **kwargs)
        return pd.concat([data] * 1001, ignore_index=True), check, report
    monkeypatch.setattr(module, "prepare_query_molecules", excessive)
    with pytest.raises(ValueError, match="1,000"):
        insert_new_molecules(reference, pd.DataFrame({"smiles": ["CCF"]}))


def test_reference_nearest_tie_uses_smallest_identifier(reference):
    result = insert_new_molecules(reference, pd.DataFrame({"smiles": ["[Xe]"]}))
    assert result.data.iloc[0].nearest_reference_id == "A"
    assert result.neighbors.measured_id.tolist() == ["A", "B", "C", "D", "E", "F"]
    assert result.neighbors["rank"].eq(1).all()


def test_insertion_does_not_repeat_reference_computation(reference, monkeypatch):
    import chemgridmap.pipeline as pipeline
    def forbidden(*args, **kwargs):
        pytest.fail("Insertion repeated reference computation")
    for name in ("project_2d", "compute_representation", "assign_to_grid", "refine_grid_layout",
                 "reliability_columns", "structural_groups", "compute_metrics", "render_grid_map"):
        monkeypatch.setattr(pipeline, name, forbidden)
    result = insert_new_molecules(reference, pd.DataFrame({"smiles": ["CCF"]}))
    assert len(result.data) == 1


@pytest.mark.parametrize("parameter,value", [("radius", 3), ("n_bits", 256), ("lower_threshold", 5.),
                                            ("upper_threshold", 8.), ("group_similarity", .8),
                                            ("grouping", "butina"), ("structure_policy", "parent")])
def test_new_molecule_settings_cannot_drift_from_reference(reference, parameter, value):
    with pytest.raises(ValueError, match="must match the reference map"):
        insert_new_molecules(reference, pd.DataFrame({"smiles": ["CCF"]}), **{parameter: value})
