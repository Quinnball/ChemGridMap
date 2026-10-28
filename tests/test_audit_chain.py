import json

import numpy as np
import pandas as pd
import pytest

from chemgridmap import curate_chembl_activity_data, inspect_entry
from chemgridmap.activity import activity_annotations, RULE_COLUMNS
from chemgridmap.app import AppSession
from chemgridmap.inspection import load_audit_context, main


RAW = b'Smiles,pChEMBL Value,Molecule ChEMBL ID\nCCO,5.2,CHEMBL20\nCCO,7.2,CHEMBL20\nCCO,8,CHEMBL20\nCCN,5,CHEMBL200\nCCC,6.4,CHEMBL300\n'


@pytest.fixture
def case():
    import io
    curated = curate_chembl_activity_data(pd.read_csv(io.BytesIO(RAW)))
    grid = curated.molecules.copy()
    grid['grid_row'], grid['grid_col'] = range(len(grid)), 0
    return grid, curated.retained_records


@pytest.mark.parametrize('values,expected', [
    ([6.0], ('inactive', 0, 0)), ([7.0], ('active', 0, 0)),
    ([6.5], ('medium', 0, 0)), ([5.0, 6.0], ('inactive', 0, 0)),
    ([6.0, 7.0], ('medium', 1, 1)), ([6.5, 7.0], ('medium', 0, 1)),
    ([5.0, 6.01], ('inactive', 1, 1)), ([5.0, 8.0, 9.0], ('active', 1, 1)),
])
def test_independently_specified_boundary_cases(values, expected):
    summary = activity_annotations(values, lower_threshold=6, upper_threshold=7, conflict_range_threshold=1)
    assert (summary['activity_class'], summary['is_conflicted'], summary['has_class_boundary_crossing']) == expected


@pytest.mark.parametrize('column,bad', [('is_conflicted', 0), ('has_class_boundary_crossing', 0),
    ('activity_class', 'inactive'), ('pchembl_range', 0), ('pchembl_min', 4), ('pchembl_std', 0),
    ('repeat_status', 'singleton'), ('representation_smiles', 'CCCC'), ('molecule_id', 'CHEMBL999')])
def test_changed_annotations_are_rejected(case, column, bad):
    grid, records = case
    index = grid.index[grid.molecule_id.eq('CHEMBL20')][0]
    identity = grid.loc[index, 'molecule_identity_key']
    grid.loc[index, column] = bad
    with pytest.raises(ValueError, match='annotation|structure|identifier'):
        inspect_entry(grid, records, molecule_id=identity)


def test_rules_are_not_silently_assumed_for_legacy_maps(case):
    grid, records = case
    legacy = grid.drop(columns=list(RULE_COLUMNS))
    summary = inspect_entry(legacy, records, molecule_id='CHEMBL20')[2]
    assert summary['median_verified'] and not summary['quality_annotations_verified']
    assert summary['is_conflicted'] is None and summary['warnings']
    rules = dict(lower_threshold=6, upper_threshold=7, conflict_range_threshold=1)
    assert inspect_entry(legacy, records, molecule_id='CHEMBL20', parameters=rules)[2]['quality_annotations_verified']


def test_custom_rules_survive_csv_roundtrip(tmp_path):
    raw = pd.DataFrame({'Smiles':['CCO', 'CCO'], 'pChEMBL Value':[5.5, 6.0]})
    result = curate_chembl_activity_data(raw, lower_threshold=5, upper_threshold=5.5, conflict_range_threshold=.2)
    grid = result.molecules.assign(grid_row=0, grid_col=0)
    path = tmp_path/'grid.csv'
    grid.to_csv(path, index=False)
    grid = pd.read_csv(path)
    summary = inspect_entry(grid, result.retained_records, cell=(0, 0))[2]
    assert summary['activity_class'] == 'active' and summary['is_conflicted']
    assert not summary['has_class_boundary_crossing'] and summary['quality_annotations_verified']
    with pytest.raises(ValueError, match='rules differ'):
        inspect_entry(grid, result.retained_records, cell=(0, 0),
                      parameters=dict(lower_threshold=6, upper_threshold=7, conflict_range_threshold=1))


def test_omitted_identity_record_is_detected_even_if_counts_are_changed(case):
    grid, records = case
    index = grid.index[grid.molecule_id.eq('CHEMBL20')][0]
    grid.loc[index, 'source_rows'], grid.loc[index, 'n_records'] = '0;1', 2
    with pytest.raises(ValueError, match='omit records'):
        inspect_entry(grid, records, molecule_id='CHEMBL20')


def test_manifest_binds_original_files_and_all_quality_flags(tmp_path):
    app = AppSession(tmp_path)
    app.load(RAW, 'raw.csv')
    app.start({'output_formats':['svg']}, background=False)
    assert app.job['state'] == 'done', app.job
    grid, records = app.files['map_grid_coordinates.csv'], app.files['map_retained_records.csv']
    assert load_audit_context(grid, records)[1]
    output = tmp_path/'inspect'
    assert main([str(grid), '--records', str(records), '--molecule-id', 'CHEMBL20', '--output-dir', str(output)]) == 0
    summary = json.loads((output/'audit_summary.json').read_text())
    assert summary['run_binding_verified'] and summary['quality_annotations_verified']
    other = tmp_path/'other_records.csv'
    other.write_bytes(records.read_bytes().replace(b'5.2', b'5.3'))
    with pytest.raises(ValueError, match='digest differs'):
        load_audit_context(grid, other)
    with pytest.raises(ValueError, match='does not exist'):
        load_audit_context(grid, records, manifest_path=tmp_path/'missing_run_manifest.json')


@pytest.mark.parametrize('value', [float('nan'), float('inf')])
def test_nonfinite_annotation_rules_are_rejected(value):
    with pytest.raises(ValueError, match='finite'):
        activity_annotations([5.0], lower_threshold=6, upper_threshold=7, conflict_range_threshold=value)
