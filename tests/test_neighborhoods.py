import numpy as np
import pandas as pd

from chemgridmap.neighborhoods import neighborhood_evidence


def test_ties_are_included_and_selected_molecule_is_excluded():
    grid = pd.DataFrame({'grid_x':[0,1,-1,0], 'grid_y':[0,0,0,2], 'molecule_id':['A','B','C','D']})
    representation = np.array([[0.], [1.], [3.], [-1.]])
    table, report = neighborhood_evidence(grid, representation, 0, k=1, metric='euclidean')
    assert report['representation_count'] == report['grid_count'] == 2
    assert report['shared_count'] == 1 and report['recall'] == .5
    assert set(table.molecule_id) == {'B','C','D'}
    assert table.set_index('molecule_id').loc['B', ['representation_neighbor','grid_neighbor']].all()


def test_small_datasets_and_identical_fingerprints():
    grid = pd.DataFrame({'grid_x':[0,1,2], 'grid_y':[0,0,0]})
    values = np.array([[1,0],[1,0],[0,1]])
    table, report = neighborhood_evidence(grid, values, 0, k=10)
    assert report['effective_k'] == 2 and report['recall'] == 1
    assert table.representation_distance.tolist() == [0., 1.]
