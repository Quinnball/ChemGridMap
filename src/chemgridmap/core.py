"""Molecule preparation, grid assignment and map-quality metrics."""

from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
from rdkit import Chem
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
from sklearn.manifold import trustworthiness
from sklearn.neighbors import NearestNeighbors

try:
    from lap import lapjv
except ImportError:  # Optional speed-up for larger maps.
    lapjv = None


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
        item["source_row"] = source_index
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


def normalize_coordinates(points: np.ndarray) -> np.ndarray:
    """Scale each coordinate dimension independently to [0, 1]."""
    array = np.asarray(points, dtype=np.float64)
    if array.ndim != 2 or array.shape[1] != 2:
        raise ValueError("Expected a two-column coordinate matrix.")
    if not np.isfinite(array).all():
        raise ValueError("Coordinates contain missing or non-finite values.")
    low = array.min(axis=0)
    high = array.max(axis=0)
    scale = np.where((high - low) == 0, 1.0, high - low)
    return (array - low) / scale


def make_candidate_grid(n_points: int, padding: int = 20) -> Tuple[np.ndarray, int, int]:
    """Create a square candidate lattice with more cells than molecules."""
    if n_points < 1:
        raise ValueError("At least one point is required.")
    side = max(2, int(np.sqrt(n_points)) + int(padding))
    axis = np.linspace(0.0, 1.0, side)
    grid = np.dstack(np.meshgrid(axis, axis)).reshape(-1, 2)
    return grid.astype(np.float64), side, side


def assign_to_grid(
    points_2d: np.ndarray,
    padding: int = 20,
) -> Tuple[np.ndarray, int, int]:
    """Assign each projected point to a unique lattice cell at minimum cost."""
    points = normalize_coordinates(points_2d)
    grid, rows, cols = make_candidate_grid(len(points), padding=padding)
    cost = cdist(grid, points, metric="sqeuclidean").astype(np.float64)

    if lapjv is not None:
        scaled = cost * (100000.0 / max(float(cost.max()), 1e-12))
        _, _, column_to_row = lapjv(scaled.astype(np.float32), extend_cost=True)
        selected = np.asarray(column_to_row[: len(points)], dtype=int)
    else:
        row_indices, column_indices = linear_sum_assignment(cost)
        selected = np.empty(len(points), dtype=int)
        selected[column_indices] = row_indices

    if len(np.unique(selected)) != len(points):
        raise RuntimeError("Grid assignment did not produce unique cells.")
    return grid[selected], rows, cols


def _neighbor_indices(points: np.ndarray, k: int) -> np.ndarray:
    points = np.asarray(points, dtype=np.float64)
    if len(points) < 2:
        return np.empty((len(points), 0), dtype=int)
    neighbors = min(max(1, int(k)), len(points) - 1)
    model = NearestNeighbors(n_neighbors=neighbors + 1).fit(points)
    return model.kneighbors(points, return_distance=False)[:, 1:]


def mean_knn_absolute_difference(
    points: np.ndarray,
    values: np.ndarray,
    k: int = 10,
) -> float:
    """Mean absolute property difference to k nearest neighbors."""
    indices = _neighbor_indices(points, k)
    if indices.shape[1] == 0:
        return float("nan")
    values = np.asarray(values, dtype=np.float64)
    differences = np.abs(values[:, None] - values[indices])
    return float(np.nanmean(differences))


def label_purity(points: np.ndarray, labels: np.ndarray, k: int = 10) -> float:
    """Mean fraction of k nearest neighbors with the same label."""
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


def compute_metrics(
    representation: np.ndarray,
    projection: np.ndarray,
    grid: np.ndarray,
    values: Optional[np.ndarray] = None,
    labels: Optional[np.ndarray] = None,
    k: int = 10,
) -> pd.DataFrame:
    """Compute projection-to-grid fidelity and optional annotation metrics."""
    projection_01 = normalize_coordinates(projection)
    grid = np.asarray(grid, dtype=np.float64)
    displacement = np.linalg.norm(projection_01 - grid, axis=1)
    n_samples = len(grid)
    effective_k = min(max(1, int(k)), max(1, n_samples - 1))

    if n_samples >= 3:
        trust_k = min(effective_k, max(1, (n_samples - 1) // 2))
        trust = float(
            trustworthiness(projection_01, grid, n_neighbors=trust_k)
        )
    else:
        trust_k = 0
        trust = float("nan")

    row = {
        "n_samples": n_samples,
        "k": effective_k,
        "trustworthiness_k": trust_k,
        "projection_to_grid_trustworthiness": trust,
        "projection_grid_knn_overlap": knn_overlap(projection_01, grid, k=effective_k),
        "mean_grid_displacement": float(np.mean(displacement)),
        "p95_grid_displacement": float(np.quantile(displacement, 0.95)),
    }

    if values is not None:
        row.update(
            {
                "representation_knn_value_abs_diff": mean_knn_absolute_difference(
                    representation, values, k=effective_k
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
