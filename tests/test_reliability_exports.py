"""Reliability colors preserve cells and travel with browser exports."""

import io
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile

import pandas as pd
import pytest
import numpy as np
from PIL import Image

from chemgridmap.app import AppSession
from chemgridmap import pipeline
from chemgridmap.plotting import render_grid_map


@pytest.fixture
def colored_grid():
    return pd.DataFrame({
        "canonical_smiles": ["CCO", "CCN", "CCC"],
        "activity_class": ["active", "medium", "inactive"],
        "grid_col": [0, 2, 1], "grid_row": [0, 1, 2],
        "class_status": ["supported", "uncertain", "uncertain"],
        "neighbor_recall": [.1, .5, .9],
        "structure_group": ["G001", "G002", "G001"],
        "structure_group_color": ["#0072B2", "#D55E00", "#0072B2"],
    })


def svg_cells(path):
    root = ET.parse(path).getroot()
    cells = [node for node in root.iter() if "data-grid-row" in node.attrib]
    return root, [(node.attrib["data-grid-row"], node.attrib["data-grid-col"],
                   next(child for child in node if child.tag.endswith("rect")).attrib) for node in cells]


@pytest.mark.parametrize("mode,expected", [
    ("activity", ["#2A9D8F", "#E6A23C", "#C44E52"]),
    ("class_support", ["#2A9D8F", "#f6dfb8", "#ebc0c1"]),
    ("neighbor_recall", ["#cde1f5", "#7599bf", "#1d5689"]),
    ("structure_group", ["#0072B2", "#D55E00", "#0072B2"]),
])
def test_color_exports_preserve_geometry_and_have_legends(tmp_path, colored_grid, mode, expected):
    files = render_grid_map(colored_grid, tmp_path / mode, color_by=mode,
                            detail="overview", formats=["svg", "png"])
    root, cells = svg_cells(files["svg"])
    assert [(row, col, rect["x"], rect["y"], rect["width"], rect["height"])
            for row, col, rect in cells] == [
                ("0", "0", "0", "200", "100", "100"),
                ("1", "2", "200", "100", "100", "100"),
                ("2", "1", "100", "0", "100", "100")]
    assert [rect["fill"] for _, _, rect in cells] == expected
    legends = [node for node in root.iter() if node.get("data-color-legend") == mode]
    assert bool(legends) == (mode != "activity")
    assert Path(files["png"]).read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert Path(files["png"]).stat().st_size > 1000


def test_activity_svg_default_matches_and_uses_opaque_tiles(tmp_path, colored_grid):
    data = colored_grid[["canonical_smiles", "activity_class", "grid_col", "grid_row"]]
    implicit = render_grid_map(data, tmp_path / "implicit", detail="overview", formats=["svg", "png"])
    explicit = render_grid_map(data, tmp_path / "explicit", detail="overview", formats=["svg", "png"], color_by="activity")
    for extension in ("svg", "png"):
        assert Path(implicit[extension]).read_bytes() == Path(explicit[extension]).read_bytes()
    assert 'fill-opacity="1"' in Path(explicit["svg"]).read_text()


def test_unknown_color_mode_is_rejected(tmp_path, colored_grid):
    with pytest.raises(ValueError, match="Unknown map color mode"):
        render_grid_map(colored_grid, tmp_path / "bad", color_by="unknown")


@pytest.mark.parametrize("layer_formats,main_formats,expected", [
    (None, ["svg"], ["svg"]),
    (None, ["png"], ["svg", "png"]),
    (["svg"], ["svg", "png"], ["svg"]),
    (["png"], ["svg"], ["png"]),
])
def test_additional_layer_formats_are_independent_of_main_map(
        tmp_path, colored_grid, layer_formats, main_formats, expected):
    result = pipeline.build_grid_map(
        colored_grid, output_dir=tmp_path, name="map", smiles_col="canonical_smiles",
        grouping="none", render_detail="overview", output_formats=main_formats,
        map_color_modes=["neighbor_recall"], map_color_formats=layer_formats)
    assert {path.name for path in tmp_path.glob("map_neighbor_recall.*")} == {
        f"map_neighbor_recall.{extension}" for extension in expected}
    assert {path.name for path in tmp_path.glob("map.*")} == {
        f"map.{extension}" for extension in main_formats}
    manifest = json.loads(Path(result.output_files["run_manifest"]).read_text())
    assert manifest["parameters"]["map_color_formats"] == layer_formats


@pytest.mark.parametrize("extension", ["pdf", "jpeg"])
def test_additional_layers_reject_unsupported_formats(colored_grid, extension):
    with pytest.raises(ValueError, match="map_color_formats only supports 'svg' and 'png'"):
        pipeline.build_grid_map(colored_grid, map_color_formats=[extension])


@pytest.mark.parametrize("mode,column", [("class_support", "class_status"),
                                        ("neighbor_recall", "neighbor_recall")])
def test_reliability_color_requires_the_estimated_values(tmp_path, colored_grid, mode, column):
    with pytest.raises(ValueError, match=f"requires the {column} column"):
        render_grid_map(colored_grid.drop(columns=column), tmp_path / mode, color_by=mode)


def test_support_legend_uses_the_configured_threshold(tmp_path, colored_grid):
    colored_grid.attrs["activity_confidence_fit"] = {"confidence_threshold": .9}
    files = render_grid_map(colored_grid, tmp_path / "support", color_by="class_support", formats=["svg"])
    assert "class support below 0.9" in Path(files["svg"]).read_text()


@pytest.mark.parametrize("threshold", [.8, .9])
def test_pipeline_passes_support_threshold_to_static_layers(tmp_path, colored_grid, threshold):
    colored_grid.attrs["activity_confidence_fit"] = {
        "confidence_threshold": threshold,
        "offset_display": "uncertain",
        "merge_identical_documents": True,
        "document_offsets": [],
        "document_pair_count": 0,
        "observed_document_pair_count": 0,
        "document_pair_possible_count": 0,
    }
    files = pipeline.build_grid_map(
        colored_grid, output_dir=tmp_path, name="support", smiles_col="canonical_smiles",
        x_col="grid_col", y_col="grid_row", grid_padding=1, grouping="none",
        output_formats=["svg"], map_color_modes=["class_support"],
        map_color_formats=["svg"],
    ).output_files
    assert f"class support below {threshold:g}" in Path(files["class_support_svg"]).read_text()


def test_offset_touched_support_cell_has_distinct_marker_and_legend(tmp_path, colored_grid):
    touched = colored_grid.copy()
    touched["offset_pair_touched"] = [False, True, False]
    files = render_grid_map(touched, tmp_path / "support", color_by="class_support", formats=["svg"])
    svg = Path(files["svg"]).read_text()
    assert 'data-offset-pair-marker="true"' in svg
    assert "◣ offset document pair" in svg


def test_png_molecule_positions_are_identical_across_modes(tmp_path, colored_grid):
    images = {}
    for mode in ("activity", "class_support", "neighbor_recall", "structure_group"):
        files = render_grid_map(colored_grid, tmp_path / mode, color_by=mode,
                                detail="full", formats=["png"])
        images[mode] = np.asarray(Image.open(files["png"]))[:, :, :3]
    height, width, _ = images["activity"].shape
    # Black molecular bonds are independent of tile color and reveal any map movement.
    bonds = images["activity"].max(axis=2) < 10
    assert bonds.sum() > 1000
    for mode, image in images.items():
        assert image.shape[1] == width
        assert np.array_equal(bonds, image[:height].max(axis=2) < 10), mode


@pytest.mark.parametrize("grouping,merge_identical", [("butina", True), ("none", False)])
def test_browser_exports_offsets_colors_and_evidence_policy(tmp_path, monkeypatch, grouping, merge_identical):
    rendered_modes = []
    def record_render(*args, **kwargs):
        rendered_modes.append(kwargs.get("color_by", "activity"))
        return render_grid_map(*args, **kwargs)
    monkeypatch.setattr(pipeline, "render_grid_map", record_render)
    rows = []
    for i, smiles in enumerate(["CCO", "CCN", "CCC"]):
        for document, value in [("DOC_A", 5 + i * .2), ("DOC_B", 8.5 + i * .2)]:
            rows.append({"Smiles": smiles, "pChEMBL Value": value,
                         "Molecule ChEMBL ID": f"CHEMBL{i + 1}", "Document ChEMBL ID": document})
    rows.append({**rows[0], "Document ChEMBL ID": "DOC_COPY"})
    app = AppSession(tmp_path)
    app.load(pd.DataFrame(rows).to_csv(index=False).encode(), "offsets.csv")
    app.start({"output_formats": ["svg", "png", "pdf"], "grouping": grouping,
               "merge_identical_documents": merge_identical}, background=False)
    assert app.job["state"] == "done", app.job
    build_files = set(app.files) - {"chemgridmap_results.zip"}
    detail, _ = app.inspect(molecule_id="CHEMBL1")
    assert detail["molecule"]["n_evidence_units"] == (2 if merge_identical else 3)
    assert detail["molecule"]["discordant_document_pairs"] == 1
    assert detail["summary"]["offset_pair_touched"] is True
    assert detail["summary"]["status_reason"] == "offset_pair"
    assert any("Offset document pair: DOC_A vs DOC_B, median difference 3.5 log units, n shared 3; "
                   "the pair touches 3 entries. This entry has records from DOC_A, DOC_B, "
                   "so it is shown as uncertain." in warning
                   for warning in detail["summary"]["warnings"])
    modes = ["neighbor_recall", "class_support"]
    if grouping != "none":
        modes.append("structure_group")
    assert rendered_modes == ["activity", *modes]
    assert not any(name.startswith("map_activity.") for name in app.files)
    layer_stems = {f"map_{mode}" for mode in ("activity", "neighbor_recall", "class_support", "structure_group")}
    assert {path.name for path in Path(app.job["output_dir"]).iterdir() if path.stem in layer_stems} == {
        f"map_{mode}.svg" for mode in modes}
    manifest = json.loads(app.files["map_run_manifest.json"].read_text())
    assert detail["summary"]["sigma_doc"] == manifest["activity_confidence"]["sigma_doc"]
    assert detail["summary"]["sigma_doc_degrees_of_freedom"] == manifest["activity_confidence"]["sigma_doc_degrees_of_freedom"]
    assert detail["summary"]["multi_unit_entry_count"] == manifest["activity_confidence"]["multi_unit_entry_count"]
    assert manifest["parameters"]["map_color_modes"] == modes
    assert manifest["parameters"]["map_color_formats"] == ["svg"]
    with zipfile.ZipFile(app.files["chemgridmap_results.zip"]) as archive:
        assert set(archive.namelist()) == build_files
        for extension in ("svg", "png", "pdf"):
            assert archive.read(f"map.{extension}") == app.files[f"map.{extension}"].read_bytes()
        for mode in modes:
            name = f"map_{mode}.svg"
            assert name in app.job["files"]
            assert archive.read(name) == app.files[name].read_bytes()
        assert "class support below 0.8" in archive.read("map_class_support.svg").decode()
        offsets = pd.read_csv(io.BytesIO(archive.read("map_document_offsets.csv")))
        assert offsets.shared_entries.tolist() == [3]
        assert offsets.median_offset.tolist() == [3.5]
        options = json.loads(archive.read("app_options.json"))
        assert options["merge_identical_documents"] == merge_identical
