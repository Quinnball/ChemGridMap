"""Fixed-scope claim checks. Archived geometry is retained, not tuned to pass.

Writes a separate current evidence directory. No user-study outcomes are generated.
"""
from __future__ import annotations

import argparse
import contextlib
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR', str(ROOT / '.cache/matplotlib'))

import numpy as np
import pandas as pd
from rdkit import RDLogger

from chemgridmap import __version__, curate_chembl_activity_data, inspect_entry
from chemgridmap.cli import main as cli
from chemgridmap.chembl import save_chembl_curation
from chemgridmap.neighborhoods import neighborhood_evidence
from chemgridmap.provenance import environment_versions, file_digests
from chemgridmap.representations import morgan_fingerprints
from run_record_tasks import pandas_retrieve
from verify_source_context import verify as verify_primary_source


def raw_file(target):
    folder = 'chembl37_chembl240_large_scale' if target == 'chembl240' else 'chembl37_' + target
    return ROOT / 'validation_data' / folder / f'chembl37_{target}_ic50_raw.csv'


def compare_curated(original, current):
    fields = ['canonical_smiles', 'activity_pchembl', 'n_records', 'source_rows',
              'activity_class', 'is_conflicted', 'has_class_boundary_crossing', 'repeat_status']
    a = original.set_index('molecule_identity_key').sort_index()[fields].copy()
    b = current.set_index('molecule_identity_key').sort_index()[fields].copy()
    a.source_rows, b.source_rows = a.source_rows.astype(str), b.source_rows.astype(str)
    pd.testing.assert_frame_equal(a, b, check_dtype=False, atol=1e-10, rtol=0)


def audit_target(source, output, target):
    previous = source / target
    report = json.loads((previous/f'{target}_curation_report.json').read_text())
    parameters = report['parameters']
    raw = pd.read_csv(raw_file(target))
    refreshed = curate_chembl_activity_data(raw, target_id=target.upper(),
        activity_type=parameters['activity_type'], lower_threshold=parameters['lower_threshold'],
        upper_threshold=parameters['upper_threshold'], conflict_range_threshold=parameters['conflict_range_threshold'],
        molecule_identity=parameters['molecule_identity_requested'], structure_policy=parameters['structure_policy'],
        exclude_potential_duplicates=parameters['exclude_potential_duplicates'], assay_types=parameters.get('assay_types'))
    compare_curated(pd.read_csv(previous/f'{target}_molecule_level.csv'), refreshed.molecules)
    archived_records = pd.read_csv(previous/f'{target}_retained_records.csv')
    serialized_records = pd.read_csv(io.StringIO(refreshed.retained_records.to_csv(index=False)))
    pd.testing.assert_frame_equal(archived_records, serialized_records, check_dtype=False, check_like=True)
    paths = save_chembl_curation(refreshed, output/target, target)
    coords = pd.read_csv(previous/'coordinates_seed_42.csv')
    enriched = refreshed.molecules.set_index('molecule_identity_key').loc[coords.molecule_identity_key]
    for key in ('audit_lower_threshold', 'audit_upper_threshold', 'audit_conflict_range_threshold'):
        coords[key] = enriched[key].to_numpy()
    coords.to_csv(output/target/'coordinates_seed_42.csv', index=False)
    results = []
    for row in coords.itertuples():
        _, evidence, summary = inspect_entry(coords, refreshed.retained_records, molecule_id=row.molecule_identity_key,
                                            parameters=parameters)
        _, by_cell, other = inspect_entry(coords, refreshed.retained_records, cell=(row.grid_row, row.grid_col),
                                          parameters=parameters)
        baseline, median, conflict = pandas_retrieve(coords, refreshed.retained_records, identity=row.molecule_identity_key)
        pd.testing.assert_frame_equal(evidence, baseline)
        pd.testing.assert_frame_equal(evidence, by_cell)
        assert summary['quality_annotations_verified'] and other['quality_annotations_verified']
        assert summary['median_pchembl'] == median and summary['is_conflicted'] == conflict
        results.append({'target':target.upper(), 'identity':row.molecule_identity_key, 'n_records':len(evidence),
                        'id_verified':True, 'cell_verified':True, 'quality_verified':True,
                        'baseline_agrees':True, 'verified_fields':';'.join(summary['verified_fields'])})
    matrices = morgan_fingerprints(coords.canonical_smiles.tolist())
    selected = sorted(range(len(coords)), key=lambda i: hashlib.sha256(str(coords.iloc[i].molecule_identity_key).encode()).hexdigest())[:5]
    if target == 'chembl205':
        selected = sorted(set(selected) | set(np.flatnonzero(coords.molecule_identity_key.eq('CHEMBL20'))))
    neighborhood_rows = []
    for index in selected:
        table, summary = neighborhood_evidence(coords, matrices, index)
        identity = str(coords.iloc[index].molecule_identity_key)
        table.to_csv(output/target/f'neighbors_{index}.csv', index=False)
        neighborhood_rows.append({'target':target.upper(), 'identity':identity, **summary})
    source_rows = set(refreshed.retained_records.source_row) | set(refreshed.excluded_records.source_row)
    assert source_rows == set(range(len(raw)))
    assert not set(refreshed.retained_records.source_row) & set(refreshed.excluded_records.source_row)
    return pd.DataFrame(results), neighborhood_rows, file_digests(paths)


def altered_controls(source, output):
    folder = source/'chembl205'
    coords = pd.read_csv(folder/'coordinates_seed_42.csv')
    records = pd.read_csv(folder/'chembl205_retained_records.csv')
    parameters = json.loads((folder/'chembl205_curation_report.json').read_text())['parameters']
    index = coords.index[coords.molecule_identity_key.eq('CHEMBL20')][0]
    cases = []
    for column, value in [('is_conflicted',0), ('has_class_boundary_crossing',0), ('activity_class','inactive'),
                           ('pchembl_range',0), ('activity_pchembl',0), ('representation_smiles','CC'),
                           ('document_id','CHEMBL0'), ('n_records',1)]:
        changed = coords.copy()
        changed.loc[index,column] = value
        cases.append((column, changed, records))
    first_row = int(str(coords.loc[index,'source_rows']).split(';')[0])
    cases.append(('missing_record', coords, records.loc[records.source_row.ne(first_row)]))
    wrong = records.copy()
    wrong.loc[wrong.source_row.eq(first_row),'molecule_identity_key'] = 'WRONG'
    cases.append(('record_identity', coords, wrong))
    checks = []
    for name, table, evidence in cases:
        try:
            inspect_entry(table, evidence, molecule_id='CHEMBL20', parameters=parameters)
        except ValueError as error:
            checks.append({'case':name, 'rejected':True, 'message':str(error)})
        else:
            raise AssertionError(f'Altered {name} was not detected')
    pd.DataFrame(checks).to_csv(output/'altered_input_controls.csv',index=False)
    return checks


def web_route(source, output):
    raw = source/'web_export_map/chembl205_web_raw.csv'
    folder = output/'web_csv_run'
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        cli([str(raw), '--input-format','chembl','--target-id','CHEMBL205','--projection','umap',
             '--output-dir',str(folder),'--name','web','--output-formats','svg'])
        cli(['inspect',str(folder/'web_grid_coordinates.csv'),'--records',str(folder/'web_retained_records.csv'),
             '--molecule-id','CHEMBL20','--output-dir',str(folder/'audit')])
    (folder/'commands.log').write_text(buffer.getvalue())
    table = pd.read_csv(folder/'web_grid_coordinates.csv')
    reference_path = source/'web_export_map/chembl205_web_grid_coordinates.csv'
    has_web_coordinates = reference_path.exists()
    reference = pd.read_csv(reference_path if has_web_coordinates else source/'chembl205/coordinates_seed_42.csv')
    left, right = table.set_index('canonical_smiles').sort_index(), reference.set_index('canonical_smiles').sort_index()
    assert left.index.equals(right.index)
    assert np.allclose(left.activity_pchembl, right.activity_pchembl, rtol=0,atol=1e-10)
    geometry_equal = bool(np.array_equal(left[['grid_row','grid_col']],right[['grid_row','grid_col']])) if has_web_coordinates else None
    summary = json.loads((folder/'audit/audit_summary.json').read_text())
    assert summary['run_binding_verified'] and summary['quality_annotations_verified']
    return {'molecules':len(table), 'raw_csv_sha256':hashlib.sha256(raw.read_bytes()).hexdigest(),
            'annotations_match_archive':True, 'grid_cells_match_archive':geometry_equal,
            'run_binding_verified':True, 'quality_annotations_verified':True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=ROOT/'paper/output/revision_v6')
    parser.add_argument('--output',type=Path,default=ROOT/'paper/output/current')
    parser.add_argument('--skip-web',action='store_true',help='Audit archived data only; report the direct CSV run as not executed.')
    args = parser.parse_args()
    output = args.output.resolve()
    if output == args.source.resolve():
        raise ValueError('Do not overwrite the archived evidence directory.')
    output.mkdir(parents=True,exist_ok=True)
    RDLogger.DisableLog('rdApp.error')
    protocol = {'version':2, 'targets':['chembl205','chembl204','chembl240'],
        'audit_population':'All seed-42 mapped molecules, without outcome selection.',
        'neighbor_examples':'Five smallest SHA256(identity) per target, plus the previously reported CHEMBL20.',
        'evaluation_scope':'Computational correctness, source-record inspection, and layout diagnostics.',
        'original_article_tables':'Condition-keyed reconciliation with the bundled CA II table transcription.',
        'no_new_geometry_optimization':True}
    (output/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    frames, neighborhoods, source_files = [], [], {}
    for target in protocol['targets']:
        frame, local, digests = audit_target(args.source,output,target)
        frames.append(frame); neighborhoods.extend(local)
        source_files[target] = {'raw':file_digests({'raw':raw_file(target)})['raw'],'curated':digests}
        print(f'{target}: {len(frame)} mapped entries verified',flush=True)
    records = pd.concat(frames,ignore_index=True)
    records.to_csv(output/'audited_molecules.csv',index=False)
    pd.DataFrame(neighborhoods).to_csv(output/'selected_neighborhoods.csv',index=False)
    controls = altered_controls(args.source,output)
    web = None if args.skip_web else web_route(args.source,output)
    primary = verify_primary_source(args.source, output)
    result = {'timestamp_utc':datetime.now(timezone.utc).isoformat(), 'software_version':__version__,
              'environment':environment_versions(), 'n_molecules':len(records),
              'n_source_records':int(records.n_records.sum()), 'n_id_cell_queries':2*len(records),
              'all_annotations_verified':bool(records.quality_verified.all()), 'altered_controls':len(controls),
              'direct_web_csv':web, 'sources':source_files,
              'primary_source_reconciliation':primary,
              'claims':{'source_record_traceability':'supported on archived test population',
                        'explicit_layout_diagnostics':'Representation and grid neighborhoods checked separately, including seed sensitivity.',
                        'primary_article_numeric_validation':f"{primary['n_primary_table_matches']} records match the transcribed source endpoint, units, and values."},
              'author_approval':'pending',
              'submission_ready':False,
              'code_digests':file_digests({str(p.relative_to(ROOT)):p for p in sorted((ROOT/'src').rglob('*.py'))})}
    (output/'claim_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('n_molecules','n_source_records','n_id_cell_queries','altered_controls','direct_web_csv','claims')},indent=2))


if __name__ == '__main__':
    main()
