"""Protect scientific values and frozen provenance during the T2c label migration."""
import json
from pathlib import Path

import pytest

from relabel_chembl204 import (assert_csv_preserved, assert_json_preserved, csv_rows,
                             latest_source_record, relabel_csv, relabel_json_value, run, verify)


def test_csv_preserves_numeric_spelling_free_text_and_rule_flags():
    before = (b'target_label,prism_column,rule,adoption_d3_ece,value,assay_description\r\n'
              b'D3,D3 | coverage,d3_ece,False,1.2300e-05,"BL21 D3 strain, unchanged"\r\n')
    after, audit = relabel_csv(before)
    assert csv_rows(after) == [
        ['target_label', 'prism_column', 'rule', 'adoption_thrombin_ece', 'value', 'assay_description'],
        ['Thrombin', 'Thrombin | coverage', 'thrombin_ece', 'False', '1.2300e-05', 'BL21 D3 strain, unchanged']]
    assert audit['changed_cells'] == 3
    assert audit['changed_headers'] == 1
    with pytest.raises(AssertionError, match='Non-label cell'):
        assert_csv_preserved(before, after.replace(b'1.2300e-05', b'0.0000123'))
    assert relabel_csv(after)[0] == after


def test_json_limits_edits_to_label_fields_and_rule_identifiers():
    before = {'target_label': 'D3', 'description': 'D3', 'value': 1.234,
              'd3_ece': {'rule': 'd3_ece', 'passed': False},
              'core_prism_columns': ['D3 | coverage']}
    edits = []
    after = relabel_json_value(before, edits=edits)
    assert after['target_label'] == 'Thrombin'
    assert after['description'] == 'D3'
    assert after['core_prism_columns'] == before['core_prism_columns']
    assert after['thrombin_ece'] == {'rule': 'thrombin_ece', 'passed': False}
    assert_json_preserved(before, after, edits)
    after['value'] += .01
    with pytest.raises(AssertionError, match='Non-label JSON'):
        assert_json_preserved(before, after, edits)


def test_saved_audit_replays_and_second_run_is_read_only(tmp_path):
    activity = tmp_path / 'activity'
    activity.mkdir()
    table = activity / 'scores.csv'
    table.write_text('target_label,value\nD3,0.848201\n')
    sub = activity / 'manifest.json'
    sub.write_text(json.dumps({'outputs_sha256': {str(table): 'old'},
                               'parameters': {'historical': 'D3'}, 'code': {'run': 'old'}}))
    historical = {'inputs_sha256': {}, 'reports_sha256': {}, 'code': {'source_sha256': {}}}
    (tmp_path / 'manifest.json').write_text(json.dumps({'t2b': historical, 'outputs_sha256': {}}))
    (tmp_path / 'RESULTS.md').write_text('D3 results\n\nTarget-label issue requiring coordinated correction: old\n')
    # Historical provenance in this fixture is deliberate and falls outside the residual scan.
    old_sub = json.loads(sub.read_text())
    old_sub['parameters']['historical'] = 'original label retained'
    sub.write_text(json.dumps(old_sub))
    result = run(tmp_path)
    assert result['passed'] and result['csv_cells_changed'] == 1
    assert json.loads(sub.read_text())['parameters'] == old_sub['parameters']
    before = {p: p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}
    assert run(tmp_path) == result
    assert all(p.read_bytes() == data for p, data in before.items())
    table.write_text('target_label,value\nThrombin,0.9\n')
    with pytest.raises(AssertionError, match='hash mismatch'):
        verify(tmp_path)


def test_latest_source_record_uses_post_relabel_sections():
    manifest = {
        'label_correction': {'code': {'source_sha256': {'old.py': 'old'}}},
        'portability_fix': {'modified_files': {}},
        'software_checks_050': {'code': {'source_sha256': {'new.py': 'new'}}},
    }
    assert latest_source_record(manifest) is manifest['software_checks_050']
    del manifest['software_checks_050']
    manifest['portability_fix']['code'] = {'source_sha256': {'portable.py': 'portable'}}
    assert latest_source_record(manifest) is manifest['portability_fix']


def test_verify_allows_only_registered_post_relabel_sections(tmp_path):
    activity = tmp_path / 'activity'
    activity.mkdir()
    (activity / 'scores.csv').write_text('target_label,value\nD3,0.848201\n')
    (tmp_path / 'manifest.json').write_text(json.dumps({
        't2b': {'inputs_sha256': {}, 'reports_sha256': {}, 'code': {'source_sha256': {}}},
        'outputs_sha256': {},
    }))
    (tmp_path / 'RESULTS.md').write_text('D3 results\n')
    run(tmp_path)
    path = tmp_path / 'manifest.json'
    manifest = json.loads(path.read_text())
    manifest['portability_fix'] = {'code': {'source_sha256': {}}}
    manifest['software_checks_050'] = {'code': {'source_sha256': {}}}
    path.write_text(json.dumps(manifest))
    assert verify(tmp_path)['passed']

    manifest['unregistered_section'] = {}
    path.write_text(json.dumps(manifest))
    with pytest.raises(AssertionError):
        verify(tmp_path)
