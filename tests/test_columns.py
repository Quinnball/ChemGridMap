"""Generated schema and API references stay synchronized with the package."""

from pathlib import Path

from chemgridmap.columns import api_markdown, output_columns_markdown


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
    assert (ROOT / "docs" / "OUTPUT_COLUMNS.md").read_text(encoding="utf-8") == (
        output_columns_markdown() + offset_display_note
    )
    assert (ROOT / "docs" / "API.md").read_text(encoding="utf-8") == api_markdown()
