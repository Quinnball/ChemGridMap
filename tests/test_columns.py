"""Generated schema and API references stay synchronized with the package."""

from pathlib import Path

import pandas as pd

from chemgridmap import curate_chembl_activity_data
from chemgridmap.chembl import save_chembl_curation
from chemgridmap.cli import main
from chemgridmap.columns import api_markdown, describe, output_columns_markdown


ROOT = Path(__file__).resolve().parents[1]


def test_generated_schema_documents_match():
    offset_display_note = (
        "\n## Offset display default\n\n"
        "`offset_display` defaults to `none`: `offset_pair_touched` and\n"
        "`document_offsets.csv` are still written, and touched entries retain the\n"
        "support-layer marker, but `class_status` and `status_reason` follow support\n"
        "alone. Use `offset_display=\"uncertain\"` or `--offset-display uncertain` to\n"
        "show touched entries as uncertain.\n"
    )
    audit_export_note = (
        "\n## Audit export subsets\n\n"
        "`selected_molecule.csv` contains one row from `grid_coordinates.csv` with exactly the same columns.\n\n"
        "`source_records.csv` contains filtered rows from `retained_records.csv` with exactly the same columns.\n\n"
        "`comparison_members.csv` contains up to 12 rows from `grid_coordinates.csv` for the selected\n"
        "structure group, with exactly the same columns.\n"
    )
    assert (ROOT / "docs" / "OUTPUT_COLUMNS.md").read_text(encoding="utf-8") == (
        output_columns_markdown() + offset_display_note + audit_export_note
    )
    assert (ROOT / "docs" / "API.md").read_text(encoding="utf-8") == api_markdown()


def chembl_row(parent_id, smiles, pchembl=6.5):
    return {
        "Molecule ChEMBL ID": f"{parent_id}_record",
        "Parent Molecule ChEMBL ID": parent_id,
        "Smiles": smiles,
        "pChEMBL Value": pchembl,
        "Target ChEMBL ID": "CHEMBL_TEST",
        "Standard Type": "IC50",
        "Standard Relation": "=",
        "Standard Units": "nM",
        "Potential Duplicate": 0,
    }


def test_structure_variant_export_columns_match_documentation(tmp_path):
    curated = curate_chembl_activity_data(
        pd.DataFrame([
            chembl_row("CHEMBL_VARIANT", "O=c1cccc[nH]1"),
            chembl_row("CHEMBL_VARIANT", "Oc1ccccn1"),
        ]),
        target_id="CHEMBL_TEST",
    )
    paths = save_chembl_curation(curated, tmp_path, "variant")
    variants = pd.read_csv(paths["structure_variants"], encoding="utf-8-sig")

    assert tuple(variants.columns) == (
        "molecule_identity_key",
        "smiles",
        "normalized_smiles",
        "structure_record_count",
        "record_count",
        "standard_inchikey",
        "selected",
    )
    assert len(variants) == 2
    assert variants.standard_inchikey.nunique() == 1
    assert variants.selected.sum() == 1
    describe(variants.columns, "structure_variants.csv")
    assert "## `structure_variants.csv`" in output_columns_markdown(
        {"structure_variants.csv": variants}
    )


def test_inspect_exports_are_documented_subsets(tmp_path):
    raw = pd.read_csv(
        ROOT / "src" / "chemgridmap" / "data" / "chembl205_ic50.csv",
        nrows=1,
    )
    curated = curate_chembl_activity_data(
        raw
    )
    grid = curated.molecules.copy()
    grid["grid_row"] = 2
    grid["grid_col"] = 3
    identity = grid.iloc[0].molecule_identity_key
    coordinates = tmp_path / "grid.csv"
    records = tmp_path / "records.csv"
    grid.to_csv(coordinates, index=False)
    curated.retained_records.to_csv(records, index=False)

    assert main([
        "inspect",
        str(coordinates),
        "--records",
        str(records),
        "--molecule-id",
        str(identity),
        "--output-dir",
        str(tmp_path / "audit"),
    ]) == 0

    for filename, source, source_filename in (
        ("selected_molecule.csv", coordinates, "grid_coordinates.csv"),
        ("source_records.csv", records, "retained_records.csv"),
    ):
        exported = pd.read_csv(tmp_path / "audit" / filename, encoding="utf-8-sig")
        parent = pd.read_csv(source, encoding="utf-8-sig")
        assert tuple(exported.columns) == tuple(parent.columns)
        describe(exported.columns, source_filename)
