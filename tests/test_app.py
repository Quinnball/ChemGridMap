import hashlib
import io
import json
from html.parser import HTMLParser
import re
import threading
import urllib.error
import urllib.request
import zipfile

import pandas as pd
import pytest

from chemgridmap.app import ASSETS, AppSession, DEMO, LocalServer, point_rows, validate_options
from chemgridmap import build_grid_map, curate_chembl_activity_data


RAW = b'Smiles,pChEMBL Value,Molecule ChEMBL ID\nCCO,5.2,CHEMBL20\nCCO,7.2,CHEMBL20\nCCO,8,CHEMBL20\nCCN,5,CHEMBL200\nCCC,6.4,CHEMBL300\nCCCC,8,CHEMBL400\n'


def test_app_uses_the_same_curation_and_grid_as_library(tmp_path):
    app = AppSession(tmp_path)
    info = app.load(RAW, "example.csv")
    assert info["rows"] == 6
    app.start({"output_formats": ["svg"]}, background=False)
    assert app.job["state"] == "done", app.job
    curated = curate_chembl_activity_data(pd.read_csv(io.BytesIO(RAW)))
    # This input has no document identifiers, so class support is skipped with a warning.
    assert any("no document identifiers" in text for text in app.curation.report["warnings"])
    expected = build_grid_map(curated.molecules, smiles_col="canonical_smiles", value_col="activity_pchembl",
                              label_col="activity_class", projection="umap", random_state=42, grid_refinement="neighbors",
                              core_rule="coverage", core_min_coverage=0.8)
    pd.testing.assert_frame_equal(app.result.data, expected.data)
    for column in ("neighbor_recall", "neighbor_shared_count", "representation_neighbor_count"):
        assert column in app.result.data
    detail, evidence = app.inspect(molecule_id="CHEMBL20")
    assert detail["summary"]["n_records"] == 3
    assert detail["summary"]["median_pchembl"] == 7.2
    assert detail["summary"]["is_conflicted"]
    assert len(evidence) == 3 and '<image' not in detail['svg']
    with zipfile.ZipFile(app.files["chemgridmap_results.zip"]) as archive:
        assert archive.read("input.csv") == RAW
        assert "map_retained_records.csv" in archive.namelist()
        assert json.loads(archive.read("app_options.json"))["projection"] == "umap"
        params = json.loads(archive.read("map_run_manifest.json"))["parameters"]
        assert params["map_color_modes"] == ["neighbor_recall", "structure_group"]
        assert params["map_color_formats"] == ["svg"]
        assert params["confidence_threshold"] == 0.8
        assert params["offset_display"] == "uncertain"
        assert params["offset_min_shared"] == 3
        assert params["offset_threshold"] == 1.0
        assert params["offset_sign_agreement"] == 0.8
        assert "map_neighbor_recall.svg" in archive.namelist()
        assert "map_structure_group.svg" in archive.namelist()
        assert not any(name.startswith(("map_activity.", "map_class_support.")) for name in archive.namelist())


def test_changing_input_discards_stale_map_but_not_saved_files(tmp_path):
    app = AppSession(tmp_path)
    app.load(RAW, "first.csv")
    app.start({"output_formats": ["svg"]}, background=False)
    file = app.files['map.svg']
    app.load(RAW, "second.csv")
    assert app.result is None and not app.files and file.exists()
    with pytest.raises(ValueError, match="Build a map"):
        app.inspect(row=0, col=0)


def test_running_job_cannot_replace_input(tmp_path):
    app = AppSession(tmp_path)
    app.load(RAW, "test.csv")
    app.update(state="running")
    with pytest.raises(ValueError, match="running"):
        app.load(RAW, "other.csv")


def test_curated_molecule_mode_has_no_invented_record_provenance(tmp_path):
    app = AppSession(tmp_path)
    app.load(b'smiles,value\nCCO,5\nCCN,8\nCCC,6.5\n', "molecules.csv")
    app.start({"input_format":"molecule", "value_col":"value", "output_formats":["svg"]}, background=False)
    assert app.job["state"] == "done", app.job
    row = app.result.data.iloc[0]
    detail, records = app.inspect(row=int(row.grid_row), col=int(row.grid_col))
    assert records.empty and not detail["summary"]["median_verified"]


def test_generic_activity_records_use_shared_support_code(tmp_path):
    app = AppSession(tmp_path)
    app.load(
        b'compound,smiles,p_activity,document\n'
        b'C1,CCO,5.8,D1\nC1,OCC,6.2,D2\n'
        b'C2,c1ccccc1,7.0,D1\nC2,c1ccccc1,7.4,D2\n'
        b'C3,Cc1ccccc1,6.4,D1\nC3,Cc1ccccc1,6.8,D2\n'
        b'C4,Oc1ccccc1,5.9,D1\nC4,Oc1ccccc1,6.1,D2\n',
        "generic.csv",
    )
    app.start({
        "input_format": "molecule",
        "compound_column": "compound",
        "smiles_col": "smiles",
        "value_col": "p_activity",
        "document_column": "document",
        "grouping": "none",
        "output_formats": ["svg"],
    }, background=False)
    assert app.job["state"] == "done", app.job
    assert {"class_confidence", "document_balanced_pchembl", "n_evidence_units"} <= set(app.result.data)
    assert app.job["options"]["evidence_model"] == "gaussian"


def test_browser_excludes_identity_conflicts(tmp_path):
    rows = pd.DataFrame([
        {"Molecule ChEMBL ID": "CHEMBL1", "Parent Molecule ChEMBL ID": "CHEMBL9",
         "Smiles": "CCO", "pChEMBL Value": 6.0},
        {"Molecule ChEMBL ID": "CHEMBL2", "Parent Molecule ChEMBL ID": "CHEMBL9",
         "Smiles": "c1ccccc1", "pChEMBL Value": 6.5},
        {"Molecule ChEMBL ID": "CHEMBL3", "Parent Molecule ChEMBL ID": "CHEMBL10",
         "Smiles": "CCC", "pChEMBL Value": 7.0},
        {"Molecule ChEMBL ID": "CHEMBL4", "Parent Molecule ChEMBL ID": "CHEMBL11",
         "Smiles": "CCN", "pChEMBL Value": 5.5},
    ])
    app = AppSession(tmp_path)
    app.load(rows.to_csv(index=False).encode(), "conflicts.csv")
    app.start({
        "projection": "pca", "grouping": "none", "identity_conflicts": "exclude",
        "output_formats": ["svg"],
    }, background=False)
    assert app.job["state"] == "done", app.job
    assert "map_identity_conflicts.csv" in app.files
    conflicts = pd.read_csv(app.files["map_identity_conflicts.csv"])
    assert conflicts.molecule_identity_key.tolist() == ["CHEMBL9"]
    options = json.loads(app.files["app_options.json"].read_text(encoding="utf-8"))
    assert options["identity_conflicts"] == "exclude"


def test_error_is_visible_and_does_not_leave_old_outputs(tmp_path):
    app = AppSession(tmp_path)
    app.load(b'smiles\nCCO\n', 'missing.csv')
    app.start({}, background=False)
    assert app.job['state'] == 'error' and 'pChEMBL' in app.job['message']
    assert app.result is None and not app.files


@pytest.mark.parametrize("raw,columns", [(b'smiles\nCCO\nCCC\n', ['smiles']),
    (b'Smiles;pChEMBL Value\nCCO;5.4\nCCC;7.2\n', ['Smiles', 'pChEMBL Value'])])
def test_csv_delimiters_and_single_column(tmp_path, raw, columns):
    info = AppSession(tmp_path).load(raw, 'input.csv')
    assert info['columns'] == columns and info['rows'] == 2


@pytest.mark.parametrize("options", [{"grid_occupancy":0}, {"lower_threshold":float('nan')},
                                    {"random_state":-1}, {"projection":"fake"}, {"output_dir":"../outside"},
                                    {"support_level":0.5}, {"offset_min_shared":1},
                                    {"offset_min_difference":0}, {"offset_sign_agreement":0.5},
                                    {"offset_display":"fake"}, {"identity_conflicts":"fake"}])
def test_invalid_options_rejected(options):
    with pytest.raises(ValueError):
        validate_options(options)


def test_packaged_example_matches_archived_source():
    metadata = json.loads(DEMO.with_name('chembl205_source.json').read_text(encoding='utf-8'))
    assert hashlib.sha256(DEMO.read_bytes()).hexdigest() == metadata['sha256']
    assert len(pd.read_csv(DEMO)) == 2281


@pytest.fixture
def local_server(tmp_path):
    app = AppSession(tmp_path)
    try:
        server = LocalServer(app)
    except PermissionError as exc:
        pytest.skip(f"Local HTTP sockets are unavailable in this sandbox: {exc}")
    worker = threading.Thread(target=server.serve_forever, daemon=True)
    worker.start()
    yield server
    server.shutdown()
    server.server_close()
    worker.join(timeout=5)


def request(server, path, *, token=True, headers=None, data=None):
    merged = {"X-ChemGridMap-Token":server.session.token} if token else {}
    merged.update(headers or {})
    return urllib.request.urlopen(urllib.request.Request(server.origin+path, data=data, headers=merged), timeout=10)


def test_local_api_requires_session_and_same_origin(local_server):
    for headers, token in [({},False), ({'Origin':'https://example.com'},True), ({'Host':'example.com'},True)]:
        with pytest.raises(urllib.error.HTTPError) as error:
            request(local_server,'/api/session',headers=headers,token=token)
        assert error.value.code == 403
    assert json.load(request(local_server,'/api/session'))['demo_available']
    with request(local_server,'/',token=False) as response:
        assert b'ChemGridMap' in response.read()
        assert 'frame-ancestors' in response.headers['Content-Security-Policy']


def test_http_input_session_restore_and_file_allowlist(local_server):
    payload = json.load(request(local_server,'/api/input',data=RAW,headers={'X-Filename':'sample.csv'}))
    assert payload['rows'] == 6
    restored = json.load(request(local_server,'/api/session'))
    assert restored['input']['filename'] == 'sample.csv'
    with pytest.raises(urllib.error.HTTPError):
        request(local_server,'/api/file?name=../../pyproject.toml')


def test_points_include_is_conflicted_without_class_support(tmp_path):
    app = AppSession(tmp_path)
    app.load(RAW, "sample.csv")
    app.start({
        "projection": "pca",
        "grouping": "none",
        "output_formats": ["svg"],
    }, background=False)
    assert app.job["state"] == "done", app.job
    points = point_rows(app.result.data)
    assert points and all("is_conflicted" in point for point in points)
    assert not any("class_status" in point for point in points)
    script = (ASSETS / "app.js").read_text(encoding="utf-8")
    assert "point.wide_record_range ?? point.is_conflicted" in script


def test_workspace_dom_matches_javascript_controls():
    class Elements(HTMLParser):
        def __init__(self):
            super().__init__()
            self.ids = []
            self.attributes = {}

        def handle_starttag(self, tag, attrs):
            attrs = dict(attrs)
            if 'id' in attrs:
                self.ids.append(attrs['id'])
                self.attributes[attrs['id']] = attrs

    html = Elements()
    html.feed((ASSETS / 'index.html').read_text(encoding='utf-8'))
    assert len(html.ids) == len(set(html.ids)), 'Duplicate IDs break UI navigation.'
    controls = set(re.findall(r"\$\('([^']+)'\)", (ASSETS / 'app.js').read_text(encoding='utf-8')))
    assert controls <= set(html.ids)
    for panel in ('page-configure', 'page-explore', 'inspector'):
        assert 'hidden' in html.attributes[panel]
    assert 'open' not in html.attributes['advanced']
    script = (ASSETS / 'app.js').read_text(encoding='utf-8')
    assert "return '> 99%'" in script and "return '< 1%'" in script
    assert "offset document pair" in script


def test_number_input_defaults_satisfy_html_constraints():
    class NumberInputs(HTMLParser):
        def __init__(self):
            super().__init__()
            self.inputs = []

        def handle_starttag(self, tag, attrs):
            attrs = dict(attrs)
            if tag == 'input' and attrs.get('type') == 'number':
                self.inputs.append(attrs)

    parser = NumberInputs()
    parser.feed((ASSETS / 'index.html').read_text(encoding='utf-8'))
    assert parser.inputs

    failures = []
    for attrs in parser.inputs:
        input_id = attrs.get('id', '<unnamed>')
        value = float(attrs['value'])
        minimum = float(attrs.get('min', 0))
        if value < minimum:
            failures.append(f'{input_id}: value {value} is below min {minimum}')
        if 'max' in attrs and value > float(attrs['max']):
            failures.append(f'{input_id}: value {value} exceeds max {attrs["max"]}')

        step = attrs.get('step', '1')
        if step != 'any':
            ratio = (value - minimum) / float(step)
            if abs(ratio - round(ratio)) > 1e-9:
                failures.append(
                    f'{input_id}: value {value} is not on step {step} from min {minimum}'
                )
    assert not failures, '\n'.join(failures)


def test_map_outlines_stay_fixed_in_screen_space_and_offset_legend_is_conditional():
    script = (ASSETS / 'app.js').read_text(encoding='utf-8')
    stylesheet = (ASSETS / 'app.css').read_text(encoding='utf-8')
    inverse = 'var(--zoom-inverse, 1)'

    assert "$('map-stage').style.setProperty('--zoom-inverse', String(1 / transform.s))" in script
    assert "points.some(point => point.offset_pair_touched)" in script
    for width in ('3.5', '2.5', '2'):
        assert f'calc({width} * {inverse})' in stylesheet
    for dash in ('6', '3', '1.5', '2.5'):
        assert f'calc({dash} * {inverse})' in stylesheet
    assert 'stroke-linecap: round' in stylesheet
    for selector in (
        '#map-stage g.selected-cell>rect',
        '#map-stage circle.selected-cell',
        '#map-stage circle.neighbor-shared',
        '#map-stage circle.neighbor-grid_only',
        '#map-stage circle.neighbor-representation_only',
    ):
        assert selector in stylesheet
    for color in ('#0072B2', '#D55E00', '#7A5195', '#193e91'):
        assert color in stylesheet


def test_inspector_uses_neutral_note_for_wide_record_range():
    script = (ASSETS / 'app.js').read_text(encoding='utf-8')
    stylesheet = (ASSETS / 'app.css').read_text(encoding='utf-8')

    assert "marker.className = !summary.quality_annotations_verified && summary.median_verified ? 'flag' : summary.wide_record_range ? 'note' : 'micro'" in script
    assert "warning.className = 'flag'" in script
    note_rule = re.search(r'\.note\s*\{([^}]*)\}', stylesheet)
    assert note_rule is not None
    declarations = note_rule.group(1)
    for declaration in ('color: var(--ink)', 'background: var(--pale)', 'padding: 9px', 'font-size: 11px', 'border-radius: 4px'):
        assert declaration in declarations
    assert '--red' not in declarations


def test_inspector_neighbor_terms_follow_representation():
    html = (ASSETS / 'index.html').read_text(encoding='utf-8')
    script = (ASSETS / 'app.js').read_text(encoding='utf-8')

    assert '<summary id="neighbor-title">Representation vs grid neighbors</summary>' in html
    assert "job.options.representation === 'morgan'" in script
    assert "${fingerprintNeighbors ? 'Fingerprint' : 'Representation'} vs grid neighbors" in script
    assert "${neighbors.representation_count} ${fingerprintNeighbors ? 'fingerprint' : 'representation'} neighbors" in script
    assert '${neighbors.representation_count} representation neighbors' not in script


def test_inspector_kth_neighbor_line_is_compact():
    script = (ASSETS / 'app.js').read_text(encoding='utf-8')

    assert 'function ordinal(n)' in script
    assert 'const lastTwo = n % 100;' in script
    assert "lastTwo >= 11 && lastTwo <= 13 ? 'th'" in script
    assert "({1: 'st', 2: 'nd', 3: 'rd'})[n % 10] || 'th'" in script
    assert '[`${ordinal(neighbors.effective_k)} fingerprint neighbor`, `Tanimoto ${summary.neighbor_kth_similarity.toFixed(2)}`]' in script
    assert '[`${ordinal(neighbors.effective_k)} representation neighbor`, `distance ${summary.neighbor_kth_distance.toFixed(2)}`]' in script
    assert 'summary.neighbor_kth_label' not in script

    suffixes = dict(re.findall(r"([123]): '([a-z]{2})'", script))
    examples = {1: '1st', 2: '2nd', 3: '3rd', 4: '4th', 11: '11th', 12: '12th',
                13: '13th', 21: '21st', 22: '22nd', 23: '23rd', 101: '101st'}
    for number, expected in examples.items():
        suffix = 'th' if 11 <= number % 100 <= 13 else suffixes.get(str(number % 10), 'th')
        assert f'{number}{suffix}' == expected


@pytest.mark.parametrize('filename', ['app.js', 'app.css', 'example-map.svg', 'example-molecule.svg'])
def test_workspace_assets_are_served_locally(local_server, filename):
    with request(local_server, '/' + filename, token=False) as response:
        assert response.read() == (ASSETS / filename).read_bytes()
