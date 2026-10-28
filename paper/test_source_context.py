"""Regression checks for condition-keyed primary-source reconciliation."""
import pandas as pd
import pytest

from verify_source_context import PRIMARY, compare_primary_tables


def sample_records():
    return pd.DataFrame([{
        'source_row': 436,
        'assay_description': 'Inhibition of human carbonic anhydrase 2 in presence of 0.01 mM dithiothreitol',
        'standard_value': 68.0, 'standard_type': 'IC50', 'standard_units': 'nM',
    }])


def primary_row():
    frame = pd.read_csv(PRIMARY, keep_default_na=False)
    return frame[frame.additive.eq('DTT') & frame.concentration_mM.eq(.01)]


def test_primary_match():
    result = compare_primary_tables(sample_records(), primary_row())
    assert result.primary_table_verified.all()
    assert result.primary_ic50_nM.iloc[0] == 68


@pytest.mark.parametrize('column,value', [('standard_value', 67), ('standard_type', 'Ki'),
                                         ('standard_units', 'uM')])
def test_primary_mismatch_is_not_verified(column, value):
    records = sample_records()
    records.loc[0, column] = value
    assert not compare_primary_tables(records, primary_row()).primary_table_verified.any()


def test_missing_condition_is_not_silently_dropped():
    records = sample_records()
    records.loc[0, 'assay_description'] = records.loc[0, 'assay_description'].replace('0.01', '0.05')
    with pytest.raises(ValueError, match='conditions are incomplete'):
        compare_primary_tables(records, primary_row())
