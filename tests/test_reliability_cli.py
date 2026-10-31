import json
import zipfile

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem

from chemgridmap.cli import build_parser, main
from chemgridmap.comparison import structural_groups
from chemgridmap.cores import _attached_candidate, _core_smiles_with_attachments, mcs_core

INDOLE_SULFONAMIDES = [
    "CC(=O)Oc1ccc2cc(S(N)(=O)=O)[nH]c2c1",
    "COc1ccc2cc(S(N)(=O)=O)[nH]c2c1",
    "COc1ccc2[nH]c(S(N)(=O)=O)cc2c1",
]

KEKULE_FALLBACK_CASES = [
    (
        [
            "Cc1nc(N(C)C(=O)Cc2ccc(-c3cc(F)ccc3F)cc2)sc1S(N)(=O)=O",
            "Cc1nc(N(C)C(=O)Cc2ccc(-c3ccccn3)cc2)sc1S(N)(=O)=O",
        ],
        22,
    ),
    (
        ["N#Cc1ccc(=O)c2ccccc2c1", "N#Cc1ccc2ccccc2c(=O)c1"],
        13,
    ),
]


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
    assert data.offset_pair_touched.all()
    assert data.status_reason.eq("offset_pair").all()
    assert data.wide_record_range.equals(data.is_conflicted)
    offsets = pd.read_csv(tmp_path / "check_document_offsets.csv")
    assert len(offsets) == 2
    assert offsets.shared_entries.eq(6).all()
    assert offsets.median_offset.eq(2).all()
    assert offsets.total_evaluated_pairs.eq(3).all()
    manifest = json.loads((tmp_path / "check_run_manifest.json").read_text())
    assert manifest["software_version"] == "0.5.0"
    assert len(manifest["software_source_sha256"]) == 64
    assert manifest["software_source_file_count"] > 0
    assert manifest["git_commit"] is None or len(manifest["git_commit"]) == 40
    assert manifest["git_dirty"] is None or isinstance(manifest["git_dirty"], bool)
    params = manifest["parameters"]
    assert params["merge_identical_documents"] is merge_identical
    assert params["activity_confidence"] is True
    assert params["grid_refinement"] == "neighbors"
    assert params["core_rule"] == "coverage"
    assert params["map_color_modes"] == modes
    assert params["map_color_formats"] is None
    assert manifest["activity_confidence"]["discordant_document_pair_count"] == 2
    assert manifest["activity_confidence"]["offset_pair_touched_entry_count"] == 6
    assert manifest["activity_confidence"]["sigma_doc"] == manifest["activity_confidence"]["tau"]
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
                   "--core-rule", "--core-min-coverage", "--fmcs-timeout",
                   "--no-activity-confidence", "--evidence-model", "--compound-column",
                   "--document-column", "--relation-column", "--smiles-column", "--value-column"):
        assert option in help_text
    normalized_help = " ".join(help_text.split())
    assert "Coverage target for coverage/FMCS rules; range: 0.5–1 (default: 0.8)." in normalized_help
    defaults = parser.parse_args(["input.csv"])
    assert defaults.projection == "umap"
    assert defaults.grid_refinement == "neighbors"
    assert defaults.core_rule == "coverage"
    assert defaults.evidence_model == "gaussian"
    assert defaults.fmcs_timeout == 1
    args = parser.parse_args(["input.csv", "--grid-refinement", "none", "--core-rule", "most_frequent",
                              "--no-activity-confidence"])
    assert args.grid_refinement == "none" and args.core_rule == "most_frequent"
    assert args.core_min_coverage == 0.8
    assert args.no_activity_confidence


def test_generic_activity_records_match_chembl_annotation_numerically(tmp_path):
    rows = []
    for index, smiles in enumerate(
        ["CCO", "c1ccccc1", "Cc1ccccc1", "Oc1ccccc1", "Nc1ccccc1", "n1ccccc1"]
    ):
        for document, offset in [("DOC_A", 0.0), ("DOC_B", 0.25), ("DOC_C", 1.5)]:
            rows.append({
                "molecule_chembl_id": f"CHEMBL{index + 1}",
                "canonical_smiles": smiles,
                "pchembl_value": 5.0 + index * 0.2 + offset,
                "document_chembl_id": document,
                "target_chembl_id": "CHEMBL205",
                "standard_type": "IC50",
                "standard_relation": "=",
                "standard_units": "nM",
                "potential_duplicate": 0,
            })
    source = tmp_path / "chembl.csv"
    pd.DataFrame(rows).to_csv(source, index=False)
    common = [
        "--projection", "pca", "--grid-padding", "2", "--render-detail", "overview",
        "--output-formats", "svg", "--output-dir", str(tmp_path),
    ]
    assert main([
        str(source), "--input-format", "chembl", "--name", "chembl_path", *common
    ]) == 0
    retained = pd.read_csv(tmp_path / "chembl_path_retained_records.csv")
    generic = retained[[
        "molecule_identity_key", "canonical_smiles", "activity_pchembl_raw",
        "document_chembl_id", "standard_relation",
    ]].rename(columns={
        "molecule_identity_key": "compound_id",
        "canonical_smiles": "structure",
        "activity_pchembl_raw": "p_activity",
        "document_chembl_id": "document_id",
        "standard_relation": "relation",
    })
    generic_path = tmp_path / "generic.csv"
    generic.to_csv(generic_path, index=False)
    assert main([
        str(generic_path), "--input-format", "molecule", "--name", "generic_path",
        "--compound-column", "compound_id", "--smiles-column", "structure",
        "--value-column", "p_activity", "--document-column", "document_id",
        "--relation-column", "relation", *common,
    ]) == 0
    chembl = pd.read_csv(tmp_path / "chembl_path_grid_coordinates.csv").set_index("molecule_identity_key")
    generic_map = pd.read_csv(tmp_path / "generic_path_grid_coordinates.csv").set_index("molecule_identity_key")
    columns = [
        "activity_pchembl", "document_balanced_pchembl", "class_confidence",
        "class_reproducibility", "n_evidence_units", "n_documents",
    ]
    assert list(generic_map.index) == list(chembl.index)
    np.testing.assert_allclose(
        generic_map[columns].to_numpy(float),
        chembl[columns].to_numpy(float),
        rtol=0.0,
        atol=1e-12,
    )


def test_cli_fmcs_core_end_to_end_records_timeout(tmp_path):
    source = tmp_path / "molecules.csv"
    pd.DataFrame({
        "smiles": ["c1ccccc1", "Cc1ccccc1", "Oc1ccccc1", "Nc1ccccc1", "n1ccccc1", "CC(=O)O"],
        "value": [7.1, 7.0, 6.8, 6.7, 6.9, 5.0],
    }).to_csv(source, index=False)
    assert main([
        str(source), "--projection", "pca", "--grid-padding", "2",
        "--core-rule", "fmcs", "--fmcs-timeout", "1", "--output-dir", str(tmp_path),
        "--name", "fmcs", "--output-formats", "svg", "--render-detail", "overview",
    ]) == 0
    groups = json.loads((tmp_path / "fmcs_structure_groups.json").read_text())["groups"]
    assert all(group["core_rule"] == "fmcs" and isinstance(group["fmcs_timed_out"], bool) for group in groups)
    members = pd.read_csv(tmp_path / "fmcs_structure_members.csv")
    group = max(groups, key=lambda item: item["size"])
    mols = [
        Chem.MolFromSmiles(value)
        for value in members.loc[members.structure_group.eq(group["structure_group"]), "canonical_smiles"]
    ]
    direct = mcs_core(mols, min_coverage=0.8, timeout=1)
    assert group["core_smarts"] == direct["core_smarts"]
    assert group["core_coverage_count"] == direct["coverage_count"]
    assert group["fmcs_timed_out"] == direct["timed_out"]


def test_mcs_core_writes_valid_smiles_from_matched_aromatic_member():
    mols = [Chem.MolFromSmiles(smiles) for smiles in INDOLE_SULFONAMIDES]
    result = mcs_core(mols, min_coverage=0.8, timeout=10)
    core = Chem.MolFromSmiles(result["core_smiles"])
    pattern = Chem.MolFromSmarts(result["core_smarts"])
    assert core is not None
    assert core.GetNumHeavyAtoms() == 13
    assert result["core_heavy_atoms"] == 13
    assert result["coverage_count"] == 3
    assert result["core_smiles_source"] == "matched_member"
    assert result["core_smiles_attached_source"] == "matched_member"
    assert all(mol.HasSubstructMatch(pattern) for mol in mols)


@pytest.mark.parametrize(("smiles", "expected_heavy_atoms"), KEKULE_FALLBACK_CASES)
def test_mcs_core_uses_kekule_fallback_for_invalid_aromatic_fragments(
    smiles, expected_heavy_atoms
):
    mols = [Chem.MolFromSmiles(value) for value in smiles]
    result = mcs_core(mols, min_coverage=0.8, timeout=10)
    core = Chem.MolFromSmiles(result["core_smiles"])
    pattern = Chem.MolFromSmarts(result["core_smarts"])
    assert core is not None
    assert core.GetNumHeavyAtoms() == expected_heavy_atoms
    assert result["core_heavy_atoms"] == expected_heavy_atoms
    assert result["core_smiles_source"] == "matched_member_kekule"
    assert all(mol.HasSubstructMatch(pattern) for mol in mols)


def test_structural_groups_uses_fmcs_heavy_atom_count():
    data = pd.DataFrame({
        "molecule_identity_key": [f"MOL{i + 1}" for i in range(3)],
        "canonical_smiles": INDOLE_SULFONAMIDES,
    })
    _, groups, _ = structural_groups(
        data, similarity=0.5, core_rule="fmcs", fmcs_timeout=10
    )
    assert len(groups) == 1
    assert groups[0]["size"] == 3
    assert groups[0]["core_heavy_atoms"] == 13
    assert Chem.MolFromSmiles(groups[0]["core_smiles"]) is not None


def test_attached_indole_core_uses_kekule_carbazole_without_rdkit_stderr(capfd):
    indole = Chem.MolFromSmiles("c1ccc2[nH]ccc2c1")
    carbazole = Chem.MolFromSmiles("c1ccc2c(c1)c1ccccc1[nH]2")
    attached, source = _core_smiles_with_attachments(indole, [carbazole])
    parsed = Chem.MolFromSmiles(attached)
    captured = capfd.readouterr()
    assert source == "matched_member_kekule"
    assert parsed is not None
    assert any(atom.GetAtomicNum() == 0 for atom in parsed.GetAtoms())
    assert "non-ring atom" not in captured.err


def test_attached_core_adds_missing_member_bond_once_at_each_end():
    pattern = Chem.MolFromSmiles("CCCC")
    member = Chem.MolFromSmiles("C1CCC1")
    attached = _attached_candidate(
        pattern, member, member.GetSubstructMatch(pattern), kekule=False
    )
    parsed = Chem.MolFromSmiles(attached)
    _, source = _core_smiles_with_attachments(pattern, [member])
    assert attached
    assert parsed is not None
    assert sum(atom.GetAtomicNum() == 0 for atom in parsed.GetAtoms()) == 2
    assert source == "matched_member"


def test_cli_fmcs_core_accepts_aromatic_heterocycles(tmp_path):
    source = tmp_path / "indoles.csv"
    pd.DataFrame({
        "smiles": INDOLE_SULFONAMIDES,
        "value": [7.1, 7.0, 6.8],
    }).to_csv(source, index=False)
    assert main([
        str(source), "--projection", "pca", "--grid-padding", "2",
        "--group-similarity", "0.5", "--core-rule", "fmcs",
        "--fmcs-timeout", "10", "--output-dir", str(tmp_path),
        "--name", "fmcs_indoles", "--output-formats", "svg",
        "--render-detail", "overview",
    ]) == 0
    groups = json.loads((tmp_path / "fmcs_indoles_structure_groups.json").read_text())["groups"]
    assert len(groups) == 1
    assert groups[0]["core_heavy_atoms"] == 13
    assert Chem.MolFromSmiles(groups[0]["core_smiles"]) is not None


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
