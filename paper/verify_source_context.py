"""Separate database arithmetic checks from still-unverified primary-paper tables."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


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
    records['primary_table_verified'] = False
    records['primary_table_page'] = ''
    records['primary_endpoint'] = ''
    output.mkdir(parents=True, exist_ok=True)
    columns = ['source_row','activity_id','assay_chembl_id','document_chembl_id','assay_description',
               'standard_type','standard_units','standard_value','pchembl_value','recalculated_from_exported_nM',
               'absolute_log_difference','exceeds_half_last_pchembl_decimal','primary_table_verified',
               'primary_table_page','primary_endpoint']
    records[columns].to_csv(output/'source_context_checks.csv',index=False)
    result = {'doi':'10.1016/j.bmcl.2008.02.008','pubmed_id':'18295485','n_records':len(records),
              'archived_records_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
              'maximum_absolute_log_difference':float(records.absolute_log_difference.max()),
              'n_exceeding_half_last_pchembl_decimal':int(records.exceeds_half_last_pchembl_decimal.sum()),
              'pchembl_values_modified':False,
              'arithmetic_scope':'Comparison to the exported standard_value, not reconstruction of original experimental precision. Differences are recorded, not silently corrected.',
              'full_text_tables_verified':False,
              'supported':'Database source links, described concentration groups, and comparison to the published abstract.',
              'pending':'Original table numbers, endpoint definition (including IC50 versus Ki), units and additive concentrations must be checked from full text before claiming independent primary-source numerical verification.',
              'primary_sources':['https://pubmed.ncbi.nlm.nih.gov/18295485/',
                                 'https://doi.org/10.1016/j.bmcl.2008.02.008']}
    (output/'source_context_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=ROOT/'paper/output/revision_v6')
    parser.add_argument('--output',type=Path,default=ROOT/'paper/output/current')
    args=parser.parse_args()
    verify(args.source,args.output)
