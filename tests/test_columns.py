"""Generated schema and API references stay synchronized with the package."""

from pathlib import Path

from chemgridmap.columns import api_markdown, output_columns_markdown


ROOT = Path(__file__).resolve().parents[1]


def test_generated_schema_documents_match():
    assert (ROOT / "docs" / "OUTPUT_COLUMNS.md").read_text(encoding="utf-8") == output_columns_markdown()
    assert (ROOT / "docs" / "API.md").read_text(encoding="utf-8") == api_markdown()
