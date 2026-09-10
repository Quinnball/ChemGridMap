"""Molecule preparation, grid assignment and map-quality metrics."""

from __future__ import annotations

from typing import Dict, Optional, Tuple
import warnings

import numpy as np
import pandas as pd
from rdkit import Chem
from scipy.optimize import linear_sum_assignment
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import min_weight_full_bipartite_matching
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from sklearn.manifold import trustworthiness
from sklearn.neighbors import NearestNeighbors

def canonicalize_smiles(smiles: object) -> Optional[str]:
    """Return an RDKit canonical SMILES, or None for an invalid structure."""
    if pd.isna(smiles) or not str(smiles).strip():
        return None
    mol = Chem.MolFromSmiles(str(smiles).strip())
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def prepare_molecules(
    data: pd.DataFrame,
    smiles_col: str = "smiles",
    value_col: Optional[str] = None,
    label_col: Optional[str] = None,
    lower_threshold: float = 6.0,
    upper_threshold: float = 7.0,
    duplicate_policy: str = "error",
) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """Validate structures and prepare a molecule-level table.

    Duplicate canonical structures are not silently aggregated. The default
    policy raises an error so that repeated measurements can be curated before
    map construction.
    """
    if smiles_col not in data.columns:
        raise ValueError("SMILES column {!r} was not found.".format(smiles_col))
    if value_col and value_col not in data.columns:
        raise ValueError("Value column {!r} was not found.".format(value_col))
    if label_col and label_col not in data.columns:
        raise ValueError("Label column {!r} was not found.".format(label_col))
    if duplicate_policy not in {"error", "first", "keep"}:
        raise ValueError("duplicate_policy must be 'error', 'first' or 'keep'.")
    if lower_threshold >= upper_threshold:
        raise ValueError("lower_threshold must be smaller than upper_threshold.")

    rows = []
    invalid_count = 0
    for source_index, (_, row) in enumerate(data.iterrows()):
        canonical = canonicalize_smiles(row[smiles_col])
        if canonical is None:
            invalid_count += 1
            continue
        item = row.to_dict()
        item.setdefault("source_row", source_index)
        item["canonical_smiles"] = canonical
        rows.append(item)

    prepared = pd.DataFrame(rows)
    if prepared.empty:
        raise ValueError("No valid SMILES were found.")

    duplicate_mask = prepared.duplicated("canonical_smiles", keep=False)
    duplicate_rows = int(duplicate_mask.sum())
    duplicate_structures = int(prepared.loc[duplicate_mask, "canonical_smiles"].nunique())
    if duplicate_structures and duplicate_policy == "error":
        raise ValueError(
            "Found {} repeated canonical structures ({} rows). "
            "Curate repeated measurements first, or set duplicate_policy='first' "
            "or 'keep' explicitly.".format(duplicate_structures, duplicate_rows)
        )
    if duplicate_policy == "first":
        prepared = prepared.drop_duplicates("canonical_smiles", keep="first")

    if label_col:
        prepared["activity_class"] = prepared[label_col].astype(str).str.strip().str.lower()
    elif "activity_class" in prepared.columns:
        prepared["activity_class"] = (
            prepared["activity_class"].astype(str).str.strip().str.lower()
        )
    elif value_col:
        values = pd.to_numeric(prepared[value_col], errors="coerce")
        labels = np.full(len(prepared), "unlabelled", dtype=object)
        labels[values <= lower_threshold] = "inactive"
        labels[values >= upper_threshold] = "active"
        middle = values.gt(lower_threshold) & values.lt(upper_threshold)
        labels[middle] = "medium"
        prepared["activity_class"] = labels
    else:
        prepared["activity_class"] = "unlabelled"

    prepared = prepared.reset_index(drop=True)
    report = {
        "input_rows": int(len(data)),
        "valid_rows": int(len(rows)),
        "invalid_smiles": invalid_count,
        "duplicate_rows": duplicate_rows,
        "duplicate_structures": duplicate_structures,
        "output_rows": int(len(prepared)),
    }
    return prepared, report


def normalize_coordinates(
    points: np.ndarray,
    mode: str = "isotropic",
) -> np.ndarray:
    """Translate and scale 2D coordinates into a unit-square frame.

    Isotropic scaling uses one common factor for both axes, so distances and
    angles are preserved up to translation and uniform scale. Independent
    scaling is retained for reproducing maps generated before version 0.3.
    """
    array = np.asarray(points, dtype=np.float64)
    if array.ndim != 2 or array.shape[1] != 2:
        raise ValueError("Expected a two-column coordinate matrix.")
    if not np.isfinite(array).all():
        raise ValueError("Coordinates contain missing or non-finite values.")
    normalized_mode = str(mode).strip().lower()
    if normalized_mode not in {"isotropic", "independent"}:
        raise ValueError(
            "coordinate scaling must be 'isotropic' or 'independent'."
        )
    low = array.min(axis=0)
    high = array.max(axis=0)
    span = high - low
    if normalized_mode == "independent":
        scale = np.where(span == 0, 1.0, span)
        return (array - low) / scale

    common_scale = float(np.max(span))
    if common_scale == 0:
        return np.full_like(array, 0.5, dtype=np.float64)
    center = (low + high) / 2.0
    return (array - center) / common_scale + 0.5


def make_candidate_grid(
    n_points: int,
    padding: Optional[int] = None,
    target_occupancy: float = 0.40,
) -> Tuple[np.ndarray, int, int]:
    """Create a square lattice from target occupancy or explicit padding."""
    if n_points < 1:
        raise ValueError("At least one point is required.")
    if padding is not None:
        if int(padding) < 0:
            raise ValueError("padding must be zero or greater.")
        side = int(np.ceil(np.sqrt(n_points))) + int(padding)
    else:
        occupancy = float(target_occupancy)
        if not 0.0 < occupancy <= 1.0:
            raise ValueError("target_occupancy must be greater than 0 and at most 1.")
        side = int(np.ceil(np.sqrt(float(n_points) / occupancy)))
    side = max(2, side)
    axis = np.linspace(0.0, 1.0, side)
    grid = np.dstack(np.meshgrid(axis, axis)).reshape(-1, 2)
    return grid.astype(np.float64), side, side


def choose_assignment_method(
    n_points: int,
    padding: Optional[int] = None,
    target_occupancy: float = 0.40,
    method: str = "auto",
    dense_max_pairs: int = 20_000_000,
) -> str:
    """Choose dense exact matching or scalable sparse matching."""
    normalized = str(method).strip().lower()
    if normalized not in {"auto", "dense", "sparse"}:
        raise ValueError("assignment method must be 'auto', 'dense' or 'sparse'.")
    if normalized != "auto":
        return normalized
    grid, _, _ = make_candidate_grid(
        n_points,
        padding=padding,
        target_occupancy=target_occupancy,
    )
    n_pairs = int(n_points) * int(len(grid))
    return "dense" if n_pairs <= int(dense_max_pairs) else "sparse"


def _morton_codes(points: np.ndarray, bits: int = 16) -> np.ndarray:
    """Return deterministic 2D Morton codes for normalized coordinates."""
    maximum = (1 << bits) - 1
    integer = np.rint(np.clip(points, 0.0, 1.0) * maximum).astype(np.uint32)

    def spread(value: np.ndarray) -> np.ndarray:
        result = value.astype(np.uint32)
        result = (result | (result << 8)) & np.uint32(0x00FF00FF)
        result = (result | (result << 4)) & np.uint32(0x0F0F0F0F)
        result = (result | (result << 2)) & np.uint32(0x33333333)
        result = (result | (result << 1)) & np.uint32(0x55555555)
        return result

    return spread(integer[:, 0]) | (spread(integer[:, 1]) << np.uint32(1))


def _fallback_matching(points: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """Build a locality-aware unique matching used to complete a sparse graph."""
    point_codes = _morton_codes(points)
    grid_codes = _morton_codes(grid)
    point_order = np.argsort(point_codes, kind="mergesort")
    grid_order = np.argsort(grid_codes, kind="mergesort")
    sorted_grid_codes = grid_codes[grid_order]
    desired = np.searchsorted(
        sorted_grid_codes,
        point_codes[point_order],
        side="left",
    )

    selected_positions = np.empty(len(points), dtype=int)
    previous = -1
    n_grid = len(grid)
    n_points = len(points)
    for rank, preferred in enumerate(desired):
        upper = n_grid - (n_points - rank)
        position = min(max(int(preferred), previous + 1), upper)
        selected_positions[rank] = position
        previous = position

    selected = np.empty(n_points, dtype=int)
    selected[point_order] = grid_order[selected_positions]
    return selected


def _assign_sparse(
    points: np.ndarray,
    grid: np.ndarray,
    neighbors: int = 32,
    diagnostics: Optional[Dict] = None,
) -> np.ndarray:
    """Solve a sparse minimum-weight matching with a guaranteed full fallback."""
    n_points = len(points)
    n_grid = len(grid)
    candidate_count = min(max(1, int(neighbors)), n_grid)
    distances, indices = cKDTree(grid).query(points, k=candidate_count)
    if candidate_count == 1:
        distances = distances[:, None]
        indices = indices[:, None]

    rows = np.repeat(np.arange(n_points, dtype=int), candidate_count)
    columns = np.asarray(indices, dtype=int).reshape(-1)
    costs = np.square(np.asarray(distances, dtype=np.float64).reshape(-1)) + 1e-12

    fallback = _fallback_matching(points, grid)
    needs_fallback = np.all(indices != fallback[:, None], axis=1)
    if np.any(needs_fallback):
        fallback_rows = np.flatnonzero(needs_fallback)
        fallback_columns = fallback[fallback_rows]
        fallback_costs = (
            np.sum(
                np.square(points[fallback_rows] - grid[fallback_columns]),
                axis=1,
            )
            + 1e-12
        )
        rows = np.concatenate([rows, fallback_rows])
        columns = np.concatenate([columns, fallback_columns])
        costs = np.concatenate([costs, fallback_costs])

    graph = coo_matrix(
        (costs, (rows, columns)),
        shape=(n_points, n_grid),
        dtype=np.float64,
    ).tocsr()
    row_indices, column_indices = min_weight_full_bipartite_matching(graph)
    selected = np.empty(n_points, dtype=int)
    selected[row_indices] = column_indices
    if diagnostics is not None:
        diagnostics.update({
            "solver": "scipy.min_weight_full_bipartite_matching",
            "nearest_candidates": candidate_count,
            "fallback_edges_added": int(needs_fallback.sum()),
            "fallback_edges_used": int((needs_fallback & (selected == fallback)).sum()),
            "candidate_edges": int(graph.nnz),
            "max_candidates_per_molecule": int(np.diff(graph.indptr).max()),
            "cost_dtype": "float64",
        })
    return selected


def assign_to_grid(
    points_2d: np.ndarray,
    padding: Optional[int] = None,
    target_occupancy: float = 0.40,
    coordinate_scaling: str = "isotropic",
    method: str = "auto",
    dense_max_pairs: int = 20_000_000,
    sparse_neighbors: int = 128,
    diagnostics: Optional[Dict] = None,
    sparse_adaptive: bool = True,
    sparse_max_neighbors: int = 1024,
) -> Tuple[np.ndarray, int, int]:
    """Assign each projected point to a unique dense- or sparse-matched cell."""
    points = normalize_coordinates(points_2d, mode=coordinate_scaling)
    grid, rows, cols = make_candidate_grid(
        len(points),
        padding=padding,
        target_occupancy=target_occupancy,
    )
    selected_method = choose_assignment_method(
        len(points),
        padding=padding,
        target_occupancy=target_occupancy,
        method=method,
        dense_max_pairs=dense_max_pairs,
    )

    if selected_method == "sparse":
        candidate_count = min(max(1, int(sparse_neighbors)), len(grid))
        limit = min(max(candidate_count, int(sparse_max_neighbors)), len(grid))
        history = []
        previous_cost = None
        while True:
            info = {}
            selected = _assign_sparse(points, grid, neighbors=candidate_count, diagnostics=info)
            cost = float(np.square(points - grid[selected]).sum())
            relative_change = (None if previous_cost is None else
                               abs(previous_cost - cost) / max(previous_cost, 1e-15))
            converged = candidate_count == len(grid) or (
                relative_change is not None and relative_change <= 1e-3
                and info["fallback_edges_used"] == 0)
            history.append({"candidates": candidate_count, "squared_cost": cost,
                            "fallback_edges_used": info["fallback_edges_used"]})
            if not sparse_adaptive or converged or candidate_count >= limit:
                break
            previous_cost = cost
            candidate_count = min(candidate_count * 2, limit)
        info.update({"adaptive": bool(sparse_adaptive), "initial_candidates": int(sparse_neighbors),
                     "converged": bool(converged), "iterations": len(history),
                     "candidate_history": history, "convergence_tolerance": 1e-3})
        if diagnostics is not None:
            diagnostics.update(info)
        if sparse_adaptive and not converged:
            warnings.warn("Sparse assignment reached its candidate cap without convergence; "
                          "inspect the run manifest or increase sparse_max_neighbors.", RuntimeWarning)
    else:
        cost = cdist(grid, points, metric="sqeuclidean").astype(np.float64)
        row_indices, column_indices = linear_sum_assignment(cost)
        selected = np.empty(len(points), dtype=int)
        selected[column_indices] = row_indices
        if diagnostics is not None:
            diagnostics.update({"solver": "scipy.linear_sum_assignment",
                                "cost_dtype": "float64",
                                "candidate_edges": int(cost.size),
                                "fallback_edges_added": 0,
                                "fallback_edges_used": 0})

    if len(np.unique(selected)) != len(points):
        raise RuntimeError("Grid assignment did not produce unique cells.")
    return grid[selected], rows, cols


def _neighbor_indices(
    points: np.ndarray,
    k: int,
    batch_size: int = 256,
    metric: str = "euclidean",
) -> np.ndarray:
    points = np.asarray(points)
    if metric == "jaccard":
        points = points.astype(bool, copy=False)
    else:
        points = points.astype(np.float64, copy=False)
    if len(points) < 2:
        return np.empty((len(points), 0), dtype=int)
    neighbors = min(max(1, int(k)), len(points) - 1)
    model = NearestNeighbors(
        n_neighbors=neighbors + 1,
        metric=metric,
        algorithm="brute" if metric == "jaccard" else "auto",
    ).fit(points)
    output = np.empty((len(points), neighbors), dtype=int)
    step = max(1, int(batch_size))
    for start in range(0, len(points), step):
        stop = min(len(points), start + step)
        candidates = model.kneighbors(
            points[start:stop],
            return_distance=False,
        )
        for local_index, candidate_row in enumerate(candidates):
            source_index = start + local_index
            without_self = candidate_row[candidate_row != source_index]
            output[source_index] = without_self[:neighbors]
    return output


def mean_knn_absolute_difference(
    points: np.ndarray,
    values: np.ndarray,
    k: int = 10,
    metric: str = "euclidean",
    tie_inclusive: bool = True,
) -> float:
    """Mean per-molecule property difference over complete kth-distance shells."""
    if tie_inclusive:
        indices = _tie_inclusive_neighbor_indices(points, k, metric=metric)
        values = np.asarray(values, dtype=np.float64)
        means = [np.nanmean(np.abs(values[i] - values[row]))
                 for i, row in enumerate(indices) if len(row)]
        return float(np.nanmean(means)) if means else float("nan")
    indices = _neighbor_indices(points, k, metric=metric)
    if indices.shape[1] == 0:
        return float("nan")
    values = np.asarray(values, dtype=np.float64)
    differences = np.abs(values[:, None] - values[indices])
    return float(np.nanmean(differences))


def label_purity(points: np.ndarray, labels: np.ndarray, k: int = 10,
                 tie_inclusive: bool = True) -> float:
    """Mean same-label fraction; include complete distance ties by default."""
    if tie_inclusive:
        indices = _tie_inclusive_neighbor_indices(points, k)
        labels = np.asarray(labels).astype(str)
        means = [np.mean(labels[i] == labels[row])
                 for i, row in enumerate(indices) if len(row)]
        return float(np.mean(means)) if means else float("nan")
    indices = _neighbor_indices(points, k)
    if indices.shape[1] == 0:
        return float("nan")
    labels = np.asarray(labels).astype(str)
    return float(np.mean((labels[:, None] == labels[indices]).mean(axis=1)))


def knn_overlap(reference: np.ndarray, mapped: np.ndarray, k: int = 10) -> float:
    """Mean exact k-nearest-neighbor overlap between two coordinate spaces."""
    reference_indices = _neighbor_indices(reference, k)
    mapped_indices = _neighbor_indices(mapped, k)
    if reference_indices.shape[1] == 0:
        return float("nan")
    overlaps = []
    for before, after in zip(reference_indices, mapped_indices):
        overlaps.append(len(set(before).intersection(after)) / float(len(before)))
    return float(np.mean(overlaps))


def _tie_inclusive_neighbor_indices(
    points: np.ndarray,
    k: int,
    relative_tolerance: float = 1e-9,
    absolute_tolerance: float = 1e-12,
    metric: str = "euclidean",
) -> list[np.ndarray]:
    """Return all neighbors at or within the kth-neighbor distance.

    A square lattice contains many exactly equidistant neighbors. Including the
    complete boundary shell avoids making fidelity depend on an arbitrary order
    among tied grid distances.
    """
    points = np.asarray(points, dtype=np.float64)
    if len(points) < 2:
        return [np.empty(0, dtype=int) for _ in range(len(points))]

    if metric != "euclidean":
        # Chunked distances bound memory for fingerprint-space tie handling.
        neighborhoods = []
        neighbors = min(max(1, int(k)), len(points) - 1)
        for start in range(0, len(points), 128):
            distances = cdist(points[start:start + 128], points, metric=metric)
            for local, row in enumerate(distances):
                row[start + local] = np.inf
                radius = np.partition(row, neighbors - 1)[neighbors - 1]
                neighborhoods.append(np.flatnonzero(
                    row <= radius * (1 + relative_tolerance) + absolute_tolerance
                ))
        return neighborhoods

    neighbors = min(max(1, int(k)), len(points) - 1)
    model = NearestNeighbors(n_neighbors=neighbors + 1).fit(points)
    distances, indices = model.kneighbors(points, return_distance=True)
    kth_distances = np.empty(len(points), dtype=np.float64)
    for source_index, (distance_row, index_row) in enumerate(zip(distances, indices)):
        without_self = distance_row[index_row != source_index]
        kth_distances[source_index] = without_self[neighbors - 1]

    radii = (
        kth_distances * (1.0 + float(relative_tolerance))
        + float(absolute_tolerance)
    )
    tree = cKDTree(points)
    neighborhoods = tree.query_ball_point(points, radii)
    return [
        np.asarray(
            sorted(index for index in row if index != source_index),
            dtype=int,
        )
        for source_index, row in enumerate(neighborhoods)
    ]


def tie_aware_knn_agreement(
    reference: np.ndarray,
    mapped: np.ndarray,
    k: int = 10,
) -> Dict[str, float]:
    """Compare tie-inclusive kth-neighbor shells in two coordinate spaces."""
    reference_sets = _tie_inclusive_neighbor_indices(reference, k)
    mapped_sets = _tie_inclusive_neighbor_indices(mapped, k)
    if not reference_sets or not reference_sets[0].size:
        return {
            "recall": float("nan"),
            "precision": float("nan"),
            "reference_size": float("nan"),
            "mapped_size": float("nan"),
        }

    recalls = []
    precisions = []
    for reference_row, mapped_row in zip(reference_sets, mapped_sets):
        overlap = len(set(reference_row).intersection(mapped_row))
        recalls.append(overlap / float(len(reference_row)))
        precisions.append(overlap / float(len(mapped_row)))

    return {
        "recall": float(np.mean(recalls)),
        "precision": float(np.mean(precisions)),
        "reference_size": float(np.mean([len(row) for row in reference_sets])),
        "mapped_size": float(np.mean([len(row) for row in mapped_sets])),
    }


def compute_metrics(
    representation: np.ndarray,
    projection: np.ndarray,
    grid: np.ndarray,
    values: Optional[np.ndarray] = None,
    labels: Optional[np.ndarray] = None,
    k: int = 10,
    trustworthiness_sample_size: int = 3000,
    random_state: int = 42,
    representation_metric: str = "euclidean",
    coordinate_scaling: str = "isotropic",
) -> pd.DataFrame:
    """Compute projection-to-grid fidelity and optional annotation metrics."""
    projection_01 = normalize_coordinates(projection, mode=coordinate_scaling)
    grid = np.asarray(grid, dtype=np.float64)
    displacement = np.linalg.norm(projection_01 - grid, axis=1)
    n_samples = len(grid)
    effective_k = min(max(1, int(k)), max(1, n_samples - 1))

    if n_samples >= 3:
        sample_size = min(
            n_samples,
            max(3, int(trustworthiness_sample_size)),
        )
        if sample_size < n_samples:
            sample_indices = np.sort(
                np.random.RandomState(random_state).choice(
                    n_samples,
                    size=sample_size,
                    replace=False,
                )
            )
        else:
            sample_indices = np.arange(n_samples)
        trust_k = min(effective_k, max(1, (sample_size - 1) // 2))
        trust = float(
            trustworthiness(
                projection_01[sample_indices],
                grid[sample_indices],
                n_neighbors=trust_k,
            )
        )
    else:
        sample_size = n_samples
        trust_k = 0
        trust = float("nan")

    tie_aware = tie_aware_knn_agreement(
        projection_01,
        grid,
        k=effective_k,
    )
    row = {
        "n_samples": n_samples,
        "k": effective_k,
        "trustworthiness_k": trust_k,
        "trustworthiness_sample_size": sample_size,
        "trustworthiness_scope": "all_rows" if sample_size == n_samples else "induced_subsample",
        "annotation_neighborhood_policy": "tie_inclusive_kth_distance_shell",
        "projection_to_grid_trustworthiness": trust,
        "projection_grid_knn_overlap": knn_overlap(projection_01, grid, k=effective_k),
        "projection_grid_tie_aware_recall": tie_aware["recall"],
        "projection_grid_tie_aware_precision": tie_aware["precision"],
        "projection_tie_neighborhood_size": tie_aware["reference_size"],
        "grid_tie_neighborhood_size": tie_aware["mapped_size"],
        "mean_grid_displacement": float(np.mean(displacement)),
        "mean_squared_grid_displacement": float(np.mean(np.square(displacement))),
        "p95_grid_displacement": float(np.quantile(displacement, 0.95)),
    }

    if values is not None:
        row.update(
            {
                "representation_knn_value_abs_diff": mean_knn_absolute_difference(
                    representation,
                    values,
                    k=effective_k,
                    metric=representation_metric,
                ),
                "projection_knn_value_abs_diff": mean_knn_absolute_difference(
                    projection, values, k=effective_k
                ),
                "grid_knn_value_abs_diff": mean_knn_absolute_difference(
                    grid, values, k=effective_k
                ),
            }
        )
    if labels is not None:
        row.update(
            {
                "projection_label_purity": label_purity(
                    projection, labels, k=effective_k
                ),
                "grid_label_purity": label_purity(grid, labels, k=effective_k),
            }
        )
    return pd.DataFrame([row])
