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

from chemgridmap.app import ASSETS, AppSession, DEMO, LocalServer, validate_options
from chemgridmap import build_grid_map, curate_chembl_activity_data


RAW = b'Smiles,pChEMBL Value,Molecule ChEMBL ID\nCCO,5.2,CHEMBL20\nCCO,7.2,CHEMBL20\nCCO,8,CHEMBL20\nCCN,5,CHEMBL200\nCCC,6.4,CHEMBL300\nCCCC,8,CHEMBL400\n'


def test_app_uses_the_same_curation_and_grid_as_library(tmp_path):
    app = AppSession(tmp_path)
    info = app.load(RAW, "example.csv")
    assert info["rows"] == 6
    app.start({"output_formats": ["svg"]}, background=False)
    assert app.job["state"] == "done", app.job
    curated = curate_chembl_activity_data(pd.read_csv(io.BytesIO(RAW)))
    expected = build_grid_map(curated.molecules, smiles_col="canonical_smiles", value_col="activity_pchembl",
                              label_col="activity_class", random_state=42)
    pd.testing.assert_frame_equal(app.result.data, expected.data)
    detail, evidence = app.inspect(molecule_id="CHEMBL20")
    assert detail["summary"]["n_records"] == 3
    assert detail["summary"]["median_pchembl"] == 7.2
    assert detail["summary"]["is_conflicted"]
    assert len(evidence) == 3 and '<image' not in detail['svg']
    with zipfile.ZipFile(app.files["chemgridmap_results.zip"]) as archive:
        assert archive.read("input.csv") == RAW
        assert "map_retained_records.csv" in archive.namelist()
        assert json.loads(archive.read("app_options.json"))["projection"] == "pca"


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
                                    {"random_state":-1}, {"projection":"fake"}, {"output_dir":"../outside"}])
def test_invalid_options_rejected(options):
    with pytest.raises(ValueError):
        validate_options(options)


def test_packaged_example_matches_archived_source():
    metadata = json.loads(DEMO.with_name('chembl205_source.json').read_text())
    assert hashlib.sha256(DEMO.read_bytes()).hexdigest() == metadata['sha256']
    assert len(pd.read_csv(DEMO)) == 2281


@pytest.fixture
def local_server(tmp_path):
    app = AppSession(tmp_path)
    server = LocalServer(app)
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
    html.feed((ASSETS / 'index.html').read_text())
    assert len(html.ids) == len(set(html.ids)), 'Duplicate IDs break UI navigation.'
    controls = set(re.findall(r"\$\('([^']+)'\)", (ASSETS / 'app.js').read_text()))
    assert controls <= set(html.ids)
    for panel in ('page-configure', 'page-explore', 'inspector'):
        assert 'hidden' in html.attributes[panel]
    assert 'open' not in html.attributes['advanced']


@pytest.mark.parametrize('filename', ['app.js', 'app.css', 'example-map.svg', 'example-molecule.svg'])
def test_workspace_assets_are_served_locally(local_server, filename):
    with request(local_server, '/' + filename, token=False) as response:
        assert response.read() == (ASSETS / filename).read_bytes()
