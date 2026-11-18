import hashlib
import json
import threading
from pathlib import Path

import pandas as pd
import pytest
from rdkit import Chem

from chemgridmap import cores, app as app_module
from chemgridmap.app import AppSession, validate_options
from chemgridmap.cli import build_parser
from chemgridmap.comparison import structural_groups


SMILES = ['CCNC(=O)c1ccc(O)cc1', 'CCCNC(=O)c1ccc(O)cc1', 'CCCCNC(=O)c1ccc(O)cc1']
RAW = ('smiles,name,value\n' + '\n'.join(f'{s},M{i},{5+i}' for i,s in enumerate(SMILES))).encode()
OPTIONS = {'input_format': 'molecule', 'value_col': 'value', 'projection': 'pca',
           'grid_refinement': 'none', 'group_similarity': .05,
           'output_formats': ['svg'], 'render_detail': 'overview'}


def make_session(tmp_path, rule='coverage'):
    session = AppSession(tmp_path)
    session.load(RAW, 'series.csv')
    session.start({**OPTIONS, 'core_rule': rule}, background=False)
    assert session.job['state'] == 'done', session.job
    return session


def hashes(directory):
    return {p.relative_to(directory).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in directory.rglob('*') if p.is_file() and 'session' not in p.relative_to(directory).parts}


def test_brics_candidates_are_ring_containing_dummy_free_and_deterministic():
    mols = [Chem.MolFromSmiles(s) for s in SMILES]
    first, guarded = cores.brics_candidate_cores(mols)
    second, _ = cores.brics_candidate_cores(list(reversed(mols)))
    assert first == second == sorted(set(first))
    assert guarded == 0
    for smiles in first:
        mol = Chem.MolFromSmiles(smiles)
        assert mol.GetRingInfo().NumRings() > 0
        assert all(atom.GetAtomicNum() != 0 for atom in mol.GetAtoms())
    assert cores.brics_candidate_cores([Chem.MolFromSmiles('CCO')])[0] == []


def test_brics_guard_applies_only_above_ten_bonds(monkeypatch):
    actual = cores.BRICS.BRICSDecompose
    calls = []
    def capture(mol, **options):
        calls.append(options)
        return actual(mol, **options)
    monkeypatch.setattr(cores.BRICS, 'BRICSDecompose', capture)
    monkeypatch.setattr(cores.BRICS, 'FindBRICSBonds', lambda mol: iter(range(mol.GetIntProp('bonds'))))
    mols = [Chem.MolFromSmiles(s) for s in SMILES[:2]]
    for mol, count in zip(mols, (10, 11)):
        mol.SetIntProp('bonds', count)
    _, guarded = cores.brics_candidate_cores(mols)
    assert guarded == 1
    assert [call['singlePass'] for call in calls] == [False, True]
    assert all(call['returnMols'] and call['keepNonLeafNodes'] for call in calls)


def test_brics_core_enlarges_shared_amide_series_and_retains_coverage():
    mols = [Chem.MolFromSmiles(s) for s in SMILES]
    old = cores.select_core(mols)
    new = cores.select_core(mols, include_brics=True)
    assert new['core_heavy_atoms'] > old['core_heavy_atoms']
    assert new['coverage_count'] == len(mols)
    assert new['target_met'] and new['core_rule'] == 'coverage_brics'
    assert new['brics_candidates'] == new['candidates'] - old['candidates'] > 0
    assert Chem.MolFromSmiles(new['core_smiles_attached']) is not None


@pytest.mark.parametrize('rule', ['coverage', 'coverage_brics', 'coverage_recap', 'most_frequent'])
def test_build_and_switch_use_identical_group_summaries(tmp_path, rule):
    session = make_session(tmp_path)
    switched = session.set_core_rule(rule)['groups']
    built = make_session(tmp_path / ('built_' + rule), rule=rule)
    assert switched == built.result.structure_groups
    assert validate_options({'core_rule': rule})['core_rule'] == rule
    assert build_parser().parse_args(['--example', 'caii', '--core-rule', rule]).core_rule == rule


def test_switch_preserves_groups_layout_and_all_previous_files_and_caches(tmp_path, monkeypatch):
    session = make_session(tmp_path)
    before_data = session.result.data.copy(deep=True)
    before_groups = list(session.result.structure_groups)
    directory = Path(session.job['output_dir'])
    before = hashes(directory)
    payload = session.set_core_rule('coverage_brics')
    assert not payload['cached']
    pd.testing.assert_frame_equal(session.result.data, before_data)
    assert session.result.structure_groups == before_groups
    after = hashes(directory)
    assert {name: after[name] for name in before} == before
    assert set(after) - set(before) == {'core_rules/map_cores_coverage_brics.csv',
                                       'core_rules/map_cores_coverage_brics.json'}
    metadata = json.loads((directory/'core_rules/map_cores_coverage_brics.json').read_text(encoding='utf-8'))
    assert metadata['run_manifest_sha256'] == before['map_run_manifest.json']
    assert set(metadata['structure_sha256']) == {g['structure_group'] for g in before_groups}
    table = pd.read_csv(directory/'core_rules/map_cores_coverage_brics.csv')
    assert {'candidates', 'seconds', 'brics_candidates', 'brics_single_pass_members'} <= set(table)
    def unexpected(*args, **kwargs):
        pytest.fail('Cached core rules must not be recomputed.')
    monkeypatch.setattr(session, '_calculate_group_core', unexpected)
    assert session.set_core_rule('coverage')['cached']
    assert session.set_core_rule('coverage_brics')['cached']
    assert hashes(directory) == after


def test_group_comparison_and_aligned_zip_are_rule_specific_and_cached(tmp_path, monkeypatch):
    session = make_session(tmp_path)
    group_id = session.groups()[0]['structure_group']
    original = session.compare_group(group_id)
    archive = session.files[original['archive']].read_bytes()
    before = hashes(Path(session.job['output_dir']))
    session.set_core_rule('coverage_brics')
    alternative = session.compare_group(group_id)
    assert alternative['summary']['core_rule'] == 'coverage_brics'
    assert alternative['archive'].startswith('core_rules/')
    assert alternative['archive'] != original['archive']
    comparison = session.group_cores(group_id)
    assert {row['core_rule'] for row in comparison['cores']} == {'coverage', 'coverage_brics', 'coverage_recap', 'most_frequent', 'fmcs'}
    assert comparison['fmcs_timeout'] == 10
    assert all('seconds' in row and 'core_svg' in row for row in comparison['cores'])
    table = pd.read_csv(session.files[comparison['file']])
    assert len(table) == 6 and 'fmcs_timed_out' in table
    assert table.core_rule.iloc[-1] == 'generic_murcko_supplementary'
    after = hashes(Path(session.job['output_dir']))
    assert {name: after[name] for name in before} == before
    assert all(name.startswith('core_rules/') for name in set(after)-set(before))
    monkeypatch.setattr(session, '_calculate_group_core', lambda *args, **kwargs: pytest.fail('not cached'))
    assert session.group_cores(group_id) == comparison
    session.set_core_rule('coverage')
    assert session.compare_group(group_id) == original
    assert session.files[original['archive']].read_bytes() == archive
    assert session.group_cores(group_id)['rule'] == 'coverage'


@pytest.mark.parametrize('measured', [False, True])
def test_new_molecule_uses_current_rule_without_rewriting_exports(tmp_path, measured):
    session = make_session(tmp_path)
    session.load_query(b'smiles,name,value\nOc1ccccc1,Q1,8\n', 'queries.csv')
    session.insert_new({'id_col': 'name', **({'value_col': 'value', 'value_unit': 'p'} if measured else {})})
    initial, _ = session.inspect(molecule_id='Q1')
    assert initial['molecule']['contains_group_core']
    assert initial['molecule']['measurement_status'] == ('measured' if measured else 'not_measured')
    assert session.groups()[0]['n_new_members'] == 1
    before = hashes(Path(session.job['output_dir']))
    session.set_core_rule('coverage_brics')
    detail, _ = session.inspect(molecule_id='Q1')
    assert detail['molecule']['contains_group_core'] is False
    assert detail['molecule']['contains_group_core_rule'] == 'coverage_brics'
    assert detail['group_core']['core_rule'] == 'coverage_brics'
    after = hashes(Path(session.job['output_dir']))
    assert {name: after[name] for name in before} == before


def test_core_computation_does_not_hold_session_lock(tmp_path, monkeypatch):
    session = make_session(tmp_path)
    original = app_module.summarize_group_cores
    entered, release = threading.Event(), threading.Event()
    def paused(*args, **kwargs):
        entered.set()
        assert release.wait(5)
        return original(*args, **kwargs)
    monkeypatch.setattr(app_module, 'summarize_group_cores', paused)
    errors = []
    def calculate():
        try:
            session.set_core_rule('coverage_brics')
        except Exception as error:
            errors.append(error)
    worker = threading.Thread(target=calculate)
    worker.start()
    try:
        assert entered.wait(5)
        assert session.lock.acquire(timeout=.5)
        session.lock.release()
        assert len(session.points()) == len(SMILES)
        assert session.snapshot()['state'] == 'done'
        from test_app_u8_session import api
        progress = api(session, '/api/core-progress')
        assert progress['status'] == 200 and progress['data'] == {
            'running': True, 'rule': 'coverage_brics', 'done': 0, 'total': len(session.result.structure_groups)}
        assert session.group_detail(session.groups()[0]['structure_group'])['summary']['core_rule'] == 'coverage'
    finally:
        release.set()
        worker.join(10)
    assert not errors and not worker.is_alive()
    assert session.core_progress['running'] is False
    assert session.core_progress['done'] == session.core_progress['total']


def test_invalid_rules_and_fmcs_build_rule(tmp_path):
    session = make_session(tmp_path)
    for rule in ('no-such-rule', 'fmcs', None):
        with pytest.raises(ValueError, match='core rule'):
            session.set_core_rule(rule)
    with pytest.raises(ValueError, match='Unknown structure group'):
        session.group_cores('missing')
    fmcs = make_session(tmp_path/'fmcs', rule='fmcs')
    expected = fmcs.result.structure_groups
    fmcs.set_core_rule('coverage_brics')
    assert fmcs.set_core_rule('fmcs')['groups'] == expected
    with pytest.raises(ValueError, match='Build a map'):
        AppSession(tmp_path/'empty').set_core_rule('coverage')


def test_core_http_contract_errors_and_comparison(tmp_path):
    from test_app_u8_session import api
    session = make_session(tmp_path)
    unknown = api(session, '/api/core-rule', b'{"rule":"unknown"}', method='POST')
    assert unknown['status'] == 400 and 'core rule' in unknown['data']['error']
    unsupported_fmcs = api(session, '/api/core-rule', b'{"rule":"fmcs"}', method='POST')
    assert unsupported_fmcs['status'] == 400
    switched = api(session, '/api/core-rule', b'{"rule":"coverage_brics"}', method='POST')
    assert switched['status'] == 200
    assert switched['data']['rule'] == 'coverage_brics'
    group_id = switched['data']['groups'][0]['structure_group']
    groups = api(session, '/api/groups')['data']
    assert groups['rule'] == 'coverage_brics' and groups['build_rule'] == 'coverage'
    comparison = api(session, '/api/group-cores?group=' + group_id)
    assert comparison['status'] == 200
    assert len(comparison['data']['cores']) == 5
    assert comparison['data']['generic']['size'] == len(SMILES)
    assert comparison['data']['fmcs_timeout'] == 10
    assert api(session, '/api/group-cores?group=unknown')['status'] == 400
    from urllib.parse import quote
    exported = api(session, '/api/file?name=' + quote(comparison['data']['file']))
    assert exported['status'] == 200 and b'fmcs_timed_out' in exported['data']


@pytest.mark.parametrize('rule', ['coverage', 'coverage_brics', 'most_frequent'])
def test_external_group_switch_reuses_build_coverage_target(tmp_path, rule):
    session = AppSession(tmp_path)
    raw = ('smiles,name,value,series\n' + '\n'.join(
        f'{s},M{i},{5+i},provided' for i,s in enumerate(SMILES))).encode()
    session.load(raw, 'external.csv')
    session.start({**OPTIONS, 'grouping': 'external', 'cluster_col': 'series',
                   'cluster_source': 'supplied series', 'core_min_coverage': .6}, background=False)
    assert session.job['state'] == 'done', session.job
    switched = session.set_core_rule(rule)['groups']
    _, expected, _ = structural_groups(session.result.data, cluster_col='series',
        cluster_source='supplied series', core_rule=rule, core_min_coverage=.6)
    assert switched == expected
    if rule != 'most_frequent':
        assert switched[0]['core_coverage_target'] == .6
