import json
import zipfile

import pandas as pd
import pytest

from chemgridmap.cli import build_parser, main


@pytest.mark.parametrize("merge_identical", [True, False])
def test_cli_evidence_policy_offsets_layers_and_archive(tmp_path, merge_identical):
    rows = []
    for i, smiles in enumerate(["CCO", "c1ccccc1", "Cc1ccccc1", "Oc1ccccc1", "Nc1ccccc1", "n1ccccc1"]):
        for document, offset in [("DOC_A", 0), ("DOC_B", 0), ("DOC_C", 2)]:
            rows.append({"molecule_chembl_id": f"CHEMBL{i + 1}", "canonical_smiles": smiles,
                         "pchembl_value": 5.0 + i * .2 + offset, "document_chembl_id": document,
                         "target_chembl_id": "CHEMBL205", "standard_type": "IC50",
                         "standard_relation": "=", "standard_units": "nM",
                         "potential_duplicate": 0})
    source = tmp_path / "input.csv"
    pd.DataFrame(rows).to_csv(source, index=False)
    modes = ["activity", "class_support", "neighbor_recall", "structure_group"]
    args = [str(source), "--input-format", "chembl", "--projection", "pca",
            "--output-dir", str(tmp_path), "--name", "check", "--render-detail", "overview",
            "--output-formats", "svg", "png", "--map-color-modes", *modes]
    if not merge_identical:
        args.append("--no-merge-identical-documents")
    assert main(args) == 0
    data = pd.read_csv(tmp_path / "check_grid_coordinates.csv")
    assert data.n_evidence_units.eq(2 if merge_identical else 3).all()
    assert data.discordant_document_pairs.eq(2).all()
    offsets = pd.read_csv(tmp_path / "check_document_offsets.csv")
    assert len(offsets) == 2
    assert offsets.shared_entries.eq(6).all()
    assert offsets.median_offset.eq(2).all()
    assert offsets.total_evaluated_pairs.eq(3).all()
    manifest = json.loads((tmp_path / "check_run_manifest.json").read_text())
    assert manifest["software_version"] == "0.5.0.dev0"
    params = manifest["parameters"]
    assert params["merge_identical_documents"] is merge_identical
    assert params["activity_confidence"] is True
    assert params["grid_refinement"] == "neighbors"
    assert params["core_rule"] == "coverage"
    assert params["map_color_modes"] == modes
    assert params["map_color_formats"] is None
    assert manifest["activity_confidence"]["discordant_document_pair_count"] == 2
    assert (tmp_path / "check.svg").read_bytes() == (tmp_path / "check_activity.svg").read_bytes()
    assert (tmp_path / "check.png").read_bytes() == (tmp_path / "check_activity.png").read_bytes()
    with zipfile.ZipFile(tmp_path / "check_results.zip") as archive:
        assert "check_document_offsets.csv" in archive.namelist()
        for mode in modes:
            for suffix in ("svg", "png"):
                name = f"check_{mode}.{suffix}"
                assert archive.read(name) == (tmp_path / name).read_bytes()


def test_cli_reliability_help_and_legacy_switches():
    parser = build_parser()
    help_text = parser.format_help()
    for option in ("--no-merge-identical-documents", "--map-color-modes", "--grid-refinement",
                   "--core-rule", "--core-min-coverage", "--no-activity-confidence"):
        assert option in help_text
    normalized_help = " ".join(help_text.split())
    assert "Coverage target for --core-rule coverage; range: 0.5–1 (default: 0.8)." in normalized_help
    args = parser.parse_args(["input.csv", "--grid-refinement", "none", "--core-rule", "most_frequent",
                              "--no-activity-confidence"])
    assert args.grid_refinement == "none" and args.core_rule == "most_frequent"
    assert args.core_min_coverage == 0.8
    assert args.no_activity_confidence


@pytest.mark.parametrize("coverage", ["0.3", "1.2"])
def test_cli_rejects_core_min_coverage_outside_range(coverage, capsys):
    with pytest.raises(SystemExit) as error:
        build_parser().parse_args(["input.csv", "--core-min-coverage", coverage])
    assert error.value.code == 2
    message = capsys.readouterr().err
    assert "--core-min-coverage" in message
    assert "must be between 0.5 and 1 (inclusive)" in message


@pytest.mark.parametrize("coverage", ["0.5", "1.0"])
def test_cli_accepts_core_min_coverage_boundaries(coverage):
    args = build_parser().parse_args(["input.csv", "--core-min-coverage", coverage])
    assert args.core_min_coverage == float(coverage)
