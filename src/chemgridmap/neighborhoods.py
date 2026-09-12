"""Separate layout neighbors from neighbors in the representation used for a map."""
from __future__ import annotations

import numpy as np
from scipy.spatial.distance import cdist


def neighborhood_evidence(data, representation, position, *, k=10, metric="jaccard"):
    if len(data) != len(representation) or len(data) < 2 or k < 1:
        raise ValueError("Neighborhood inspection needs aligned arrays, two molecules, and k >= 1.")
    representation = np.asarray(representation)
    distances = cdist(representation[[position]], representation, metric=metric)[0]
    grid = data[['grid_x', 'grid_y']].to_numpy(float)
    grid_distances = np.linalg.norm(grid - grid[position], axis=1)
    neighbors = []
    effective_k = min(k, len(data) - 1)
    for vector in (distances, grid_distances):
        vector[position] = np.inf
        radius = np.partition(vector, effective_k - 1)[effective_k - 1]
        neighbors.append(set(np.flatnonzero(vector <= radius + 1e-12)))
    original, displayed = neighbors
    overlap = original & displayed
    positions = sorted(original | displayed, key=lambda index: (distances[index], index))
    fields = [name for name in ('molecule_identity_key', 'molecule_id', 'canonical_smiles',
              'activity_pchembl', 'activity_class', 'grid_row', 'grid_col') if name in data]
    table = data.iloc[positions][fields].copy().reset_index(drop=True)
    table['representation_neighbor'] = [index in original for index in positions]
    table['grid_neighbor'] = [index in displayed for index in positions]
    table['representation_distance'] = distances[positions]
    table['grid_distance'] = grid_distances[positions]
    summary = {'requested_k': int(k), 'effective_k': effective_k, 'representation_metric': metric,
               'tie_policy': 'all_neighbors_within_kth_distance_plus_1e-12',
               'representation_count': len(original), 'grid_count': len(displayed), 'shared_count': len(overlap),
               'recall': len(overlap) / len(original), 'precision': len(overlap) / len(displayed),
               'jaccard': len(overlap) / len(original | displayed),
               'interpretation': 'Layout adjacency is not a chemical or biological similarity claim. Counts may exceed k because all distance ties are retained.'}
    return table, summary
