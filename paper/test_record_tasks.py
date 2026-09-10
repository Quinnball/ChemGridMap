import pandas as pd
import pytest

from run_record_tasks import pandas_retrieve


def test_metadata_preserving_baseline_recovers_all_measurements():
    coordinates = pd.DataFrame([dict(molecule_identity_key="CHEMBL20", grid_row=2, grid_col=4,
                                     source_rows="0;1;2", n_records=3, activity_pchembl=7.0)])
    records = pd.DataFrame(dict(molecule_identity_key=["CHEMBL20"] * 3,
                                source_row=[0, 1, 2], activity_pchembl_raw=[5.5, 7.0, 8.0]))
    selected, median, flag = pandas_retrieve(coordinates, records, cell=(2, 4))
    assert selected.source_row.tolist() == [0, 1, 2]
    assert median == 7.0 and flag
    with pytest.raises(ValueError, match="linkage"):
        pandas_retrieve(coordinates, records.iloc[:2], identity="CHEMBL20")
