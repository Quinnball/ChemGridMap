"""Compare archived ChEMBL records with the cited CA II table transcription."""
import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PRIMARY = ROOT / 'paper/data/innocenti_2008_ca2.csv'
PROVENANCE = ROOT / 'paper/data/innocenti_2008_source.json'


def condition(description):
    if description == 'Inhibition of human carbonic anhydrase 2':
        return 'none', 0.0
    match = re.fullmatch(
        r'Inhibition of human carbonic anhydrase 2 in presence of ([\d.]+) mM (.+)',
        description,
    )
    names = {'dithiothreitol': 'DTT', 'beta-mercaptoethanol': 'MET',
             'tris-(carboxyethyl)-phosphine': 'TCP', 'threitol': 'THT'}
    if not match or match[2] not in names:
        raise ValueError(f'Unmapped primary-source condition: {description}')
    return names[match[2]], float(match[1])


def compare_primary_tables(records, primary):
    records = records.copy()
    keys = records.assay_description.map(condition).tolist()
    records[['additive', 'concentration_mM']] = pd.DataFrame(keys, index=records.index)
    # Match the experimental condition, not the activity value being checked.
    checked = records.merge(primary, on=['additive', 'concentration_mM'],
                            how='outer', validate='one_to_one', indicator=True)
    if not checked['_merge'].eq('both').all():
        raise ValueError('Primary-table and database conditions are incomplete or mismatched')
    checked = checked.drop(columns='_merge')
    checked['primary_endpoint'] = 'IC50'
    checked['primary_units'] = 'nM'
    checked['primary_value_matches'] = np.isclose(
        checked.standard_value, checked.primary_ic50_nM, rtol=0, atol=1e-10)
    checked['primary_endpoint_matches'] = checked.standard_type.eq('IC50')
    checked['primary_units_match'] = checked.standard_units.eq('nM')
    checked['primary_table_verified'] = checked[
        ['primary_value_matches', 'primary_endpoint_matches', 'primary_units_match']
    ].all(axis=1)
    return checked.sort_values('source_row')


def verify(source, output):
    path = source / 'chembl205/acetazolamide_retained_records.csv'
    retained = pd.read_csv(path)
    records = retained[retained.document_chembl_id.eq('CHEMBL1141531')].copy()
    assert len(records) == 25
    assert records.standard_type.eq('IC50').all() and records.standard_units.eq('nM').all()
    assert records.standard_relation.eq('=').all() and records.standard_value.gt(0).all()
    records['recalculated_from_exported_nM'] = 9 - np.log10(records.standard_value)
    records['absolute_log_difference'] = abs(records.recalculated_from_exported_nM - records.pchembl_value)
    records['exceeds_half_last_pchembl_decimal'] = records.absolute_log_difference.gt(.005)
    records = compare_primary_tables(records, pd.read_csv(PRIMARY, keep_default_na=False))
    provenance = json.loads(PROVENANCE.read_text())
    output.mkdir(parents=True, exist_ok=True)
    columns = ['source_row','activity_id','assay_chembl_id','document_chembl_id','assay_description',
               'standard_type','standard_units','standard_value','pchembl_value','recalculated_from_exported_nM',
               'absolute_log_difference','exceeds_half_last_pchembl_decimal','primary_table_verified',
               'primary_table_page','primary_table_number','primary_endpoint','primary_units',
               'additive','concentration_mM','primary_ic50_nM','primary_value_matches',
               'primary_endpoint_matches','primary_units_match']
    records[columns].to_csv(output/'source_context_checks.csv',index=False)
    result = {'doi':'10.1016/j.bmcl.2008.02.008','pubmed_id':'18295485','n_records':len(records),
              'archived_records_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
              'maximum_absolute_log_difference':float(records.absolute_log_difference.max()),
              'n_exceeding_half_last_pchembl_decimal':int(records.exceeds_half_last_pchembl_decimal.sum()),
              'pchembl_values_modified':False,
              'arithmetic_scope':'Comparison to the exported standard_value, not reconstruction of original experimental precision. Differences are recorded, not silently corrected.',
              'full_text_tables_verified':bool(records.primary_table_verified.all()),
              'n_primary_table_matches':int(records.primary_table_verified.sum()),
              'primary_transcription_sha256':hashlib.sha256(PRIMARY.read_bytes()).hexdigest(),
              'primary_source':provenance,
              'supported':f'{int(records.primary_table_verified.sum())} of {len(records)} archived CA II records match Tables 2-5, pp. 1900-1901. Matching includes endpoint, units and numerical value after condition-keyed alignment. This is numerical source reconciliation, not experimental replication.',
              'pending':'Author review/signoff of the AI-assisted transcription; the original paper does not report pChEMBL, and the cause of the three small database log-value discrepancies remains unresolved.',
              'primary_sources':['https://pubmed.ncbi.nlm.nih.gov/18295485/',
                                 'https://doi.org/10.1016/j.bmcl.2008.02.008']}
    (output/'source_context_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    if not result['full_text_tables_verified']:
        raise ValueError('Primary-source discrepancies recorded; do not claim complete agreement')
    return result


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=ROOT/'paper/output/revision_v6')
    parser.add_argument('--output',type=Path,default=ROOT/'paper/output/current')
    args=parser.parse_args()
    verify(args.source,args.output)
