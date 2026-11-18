"""Optional input controls must stay separate from the default import path."""

from html.parser import HTMLParser

from chemgridmap.app import ASSETS


class Workspace(HTMLParser):
    def __init__(self):
        super().__init__()
        self.elements = {}
        self.options = {}
        self.select = None

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if "id" in attrs:
            self.elements[attrs["id"]] = attrs
        if tag == "select":
            self.select = attrs["id"]
            self.options[self.select] = []
        elif tag == "option" and self.select:
            self.options[self.select].append(attrs.get("value"))

    def handle_endtag(self, tag):
        if tag == "select":
            self.select = None


def workspace():
    html = Workspace()
    html.feed((ASSETS / "index.html").read_text(encoding="utf-8"))
    return html


def test_own_data_controls_offer_record_and_molecule_inputs():
    html = workspace()
    assert html.options["input-format"] == ["chembl", "generic", "molecule"]
    assert html.options["value-unit"][0] == "p"
    assert len(html.options["value-unit"]) == 6
    assert "unit-column" in html.elements
    assert "checked" not in html.elements["use-assumed-sigma"]
    assert html.elements["assumed-sigma"]["value"] == "0.64"
    assert "disabled" in html.elements["assumed-sigma"]
    assert "assumed-sigma-unavailable" in html.elements


def test_query_upload_is_optional_and_does_not_submit_the_map():
    html = workspace()
    upload = html.elements["query-csv"]
    assert upload["type"] == "file"
    assert {".csv", ".tsv"} <= set(upload["accept"].split(","))
    for control in ("query-csv", "query-id-col"):
        assert "required" not in html.elements[control]
    assert "required" in html.elements["query-smiles-col"]
    assert html.elements["remove-query-input"]["type"] == "button"
    assert "hidden" in html.elements["query-input-info"]
    assert html.elements["molecule-id"]["placeholder"] == "Find a molecule ID"


def test_query_results_have_distinct_filter_and_inspector_sections():
    html = workspace()
    assert html.options["measurement-filter"] == ["all", "reference", "new"]
    assert "hidden" in html.elements["measurement-filter-control"]
    assert "hidden" in html.elements["measurement-badge"]
    assert "hidden" in html.elements["query-selection"]
    for section in ("query-summary", "query-neighbor-table", "measured-selection"):
        assert section in html.elements
