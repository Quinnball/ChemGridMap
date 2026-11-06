from html.parser import HTMLParser
import json
import re
from types import SimpleNamespace

import pytest

from chemgridmap import app as app_module
from chemgridmap.app import ASSETS, AppSession, LocalHandler, _group_drawing


RAW = (b'smiles,value,cluster\nCc1ccccc1,6,ring\nCCc1ccccc1,7,ring\n'
       b'Oc1ccccc1,8,ring\nCCO,5,acyclic\nCCN,6,acyclic\n')


@pytest.fixture(scope='module')
def group_app(tmp_path_factory):
    app = AppSession(tmp_path_factory.mktemp('group-highlight'))
    app.load(RAW, 'groups.csv')
    app.start({'input_format': 'molecule', 'value_col': 'value', 'projection': 'pca',
               'grouping': 'external', 'cluster_col': 'cluster', 'cluster_source': 'test input',
               'output_formats': ['svg']}, background=False)
    assert app.job['state'] == 'done', app.job
    return app


def get(app, path, *, token=True, headers=None):
    # Exercise the real routing and authorization without opening a local socket.
    handler = LocalHandler.__new__(LocalHandler)
    handler.server = SimpleNamespace(session=app, origin='http://127.0.0.1:8123')
    handler.path = path
    handler.headers = {'Host': '127.0.0.1:8123'}
    if token:
        handler.headers['X-ChemGridMap-Token'] = app.token
    handler.headers.update(headers or {})
    response = {}

    def send(payload, status=200):
        response.update(status=status, body=json.loads(json.dumps(payload, allow_nan=False)))

    handler.send = send
    handler.do_GET()
    return response['status'], response['body']


def output_files(app):
    return sorted(str(path.relative_to(app.output_root)) for path in app.output_root.rglob('*') if path.is_file())


def test_group_endpoint_matches_groups_and_comparison_without_writing_files(group_app):
    status, listing = get(group_app, '/api/groups')
    assert status == 200
    summary = next(group for group in listing['groups'] if group['core_smiles'])
    group_id = summary['structure_group']
    before = output_files(group_app)
    registered = dict(group_app.files)
    status, detail = get(group_app, '/api/group?group=' + group_id)
    assert status == 200 and set(detail) == {'summary', 'core_svg'}
    assert detail['summary'] == summary
    assert detail['core_svg'].startswith('<svg ')
    assert 'viewBox="0 0 180 180"' in detail['core_svg']
    assert output_files(group_app) == before and group_app.files == registered
    status, comparison = get(group_app, '/api/compare?group=' + group_id)
    assert status == 200
    assert detail['core_svg'] == comparison['drawings']['core_smiles']
    before = output_files(group_app)
    assert get(group_app, '/api/group?group=' + group_id) == (200, detail)
    assert output_files(group_app) == before


def test_group_endpoint_without_ring_scaffold(group_app):
    _, listing = get(group_app, '/api/groups')
    summary = next(group for group in listing['groups'] if not group['core_smiles'])
    before = output_files(group_app)
    status, detail = get(group_app, '/api/group?group=' + summary['structure_group'])
    assert status == 200
    assert detail == {'summary': summary, 'core_svg': ''}
    assert summary['size'] == 2 and summary['core_coverage_count'] == 0
    assert summary['core_exact_scaffold_count'] == 0
    assert output_files(group_app) == before


@pytest.mark.parametrize('group_id', ['', 'unknown'])
def test_group_endpoint_unknown_group_matches_comparison_error(group_app, group_id):
    response = get(group_app, '/api/group?group=' + group_id)
    assert response == get(group_app, '/api/compare?group=' + group_id)
    assert response == (400, {'error': 'Unknown structure group.'})


def test_group_endpoint_requires_a_map(tmp_path):
    app = AppSession(tmp_path)
    assert get(app, '/api/group?group=G001') == (400, {'error': 'Build a map first.'})


@pytest.mark.parametrize('token,headers', [
    (False, {}), (True, {'Origin': 'https://example.com'}), (True, {'Host': 'example.com'}),
])
def test_group_endpoint_uses_existing_api_authorization(group_app, token, headers):
    response = get(group_app, '/api/group?group=G001', token=token, headers=headers)
    assert response[0] == 403
    assert response == get(group_app, '/api/compare?group=G001', token=token, headers=headers)


def test_group_detail_is_cached_by_group(group_app, monkeypatch):
    groups = group_app.result.structure_groups
    details = {group['structure_group']: group_app.group_detail(group['structure_group']) for group in groups}

    def unexpected_drawing(*args):
        pytest.fail('Cached group details must not render another drawing.')

    monkeypatch.setattr(app_module, '_group_drawing', unexpected_drawing)
    for group_id, detail in details.items():
        assert group_app.group_detail(group_id) is detail


@pytest.mark.parametrize('action', ['load', 'start'])
def test_group_detail_cache_clears_for_new_input_or_map(tmp_path, monkeypatch, action):
    app = AppSession(tmp_path)
    app.load(RAW, 'groups.csv')
    app.group_detail_cache['G001'] = {'summary': {}}
    app.group_cache['G001'] = {'summary': {}}
    if action == 'load':
        app.load(RAW, 'replacement.csv')
    else:
        monkeypatch.setattr(app, 'run', lambda options: None)
        app.start({}, background=False)
    assert not app.group_detail_cache and not app.group_cache
    assert app.result is None


@pytest.mark.parametrize('attached,expected', [('*c1ccccc1', '*c1ccccc1'), ('', 'c1ccccc1'), (None, 'c1ccccc1')])
def test_shared_group_drawing_prefers_attached_core(monkeypatch, attached, expected):
    calls = []

    def draw(smiles, size):
        calls.append((smiles, size))
        return '<path/>'

    monkeypatch.setattr(app_module, '_svg_molecule_content', draw)
    summary = {'core_smiles': 'c1ccccc1', 'core_smiles_attached': attached}
    assert _group_drawing(summary, 'core_smiles') == (
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 180 180"><path/></svg>')
    assert calls == [(expected, 180)]


def test_group_highlight_controls_have_ids_and_default_to_none():
    class Elements(HTMLParser):
        def __init__(self):
            super().__init__()
            self.elements = {}

        def handle_starttag(self, tag, attrs):
            attrs = dict(attrs)
            if 'id' in attrs:
                self.elements[attrs['id']] = (tag, attrs)

    html = (ASSETS / 'index.html').read_text(encoding='utf-8')
    parser = Elements()
    parser.feed(html)
    for control in ('group-core-card', 'group-core-svg', 'group-core-title',
                    'group-core-coverage', 'group-core-rule', 'selected-group-highlight'):
        assert control in parser.elements
    assert 'hidden' in parser.elements['group-core-card'][1]
    assert 'disabled' in parser.elements['compare-group'][1]
    assert re.search(r'<select id="group-select">\s*<option value=""[^>]*>None</option>', html)


def test_inspect_does_not_change_highlight_selection():
    script = (ASSETS / 'app.js').read_text(encoding='utf-8')
    inspect = re.search(r'async function inspect\(query\) \{(.*?)\n\}', script, re.S)
    assert inspect is not None
    assert 'group-select' not in inspect.group(1)
    assert 'selected-group-highlight' in inspect.group(1)


def test_faded_cells_keep_inspection_outlines_in_a_separate_layer():
    script = (ASSETS / 'app.js').read_text(encoding='utf-8')
    stylesheet = (ASSETS / 'app.css').read_text(encoding='utf-8')
    # Parent opacity also fades child strokes; inspection outlines need an opaque sibling layer.
    assert "querySelector('[data-map-outlines]')?.remove()" in script
    assert "outlines.setAttribute('pointer-events', 'none')" in script
    assert "mapCells[0].parentNode.append(outlines)" in script
    assert "['class', 'data-grid-row', 'data-grid-col']" in script
    assert "outline.setAttribute('fill', 'none')" in script
    assert '#map-stage .group-faded { opacity: .2; }' in stylesheet
    for relation, color in [('shared', '#0072B2'), ('grid_only', '#D55E00'),
                            ('representation_only', '#7A5195'), ('selected', '#193e91')]:
        rule = re.search(r'\[data-outline=' + relation + r'\]\s*\{([^}]+)\}', stylesheet)
        assert rule is not None
        assert 'stroke: ' + color in rule.group(1)
        assert 'var(--zoom-inverse, 1)' in rule.group(1)
        assert 'stroke-opacity: 1' in rule.group(1)
