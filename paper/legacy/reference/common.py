"""Shared loaders and neighbourhood metrics for the review analyses.

All analyses reuse the frozen revision_v6 layouts behind the manuscript
(CA II seeds 7/21/42/84/126; D3 and hERG seeds 7/21/42) without refitting UMAP.
"""
from __future__ import annotations

import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

# Historical snapshot: the original absolute workspace path is replaced by a
# placeholder; this file is kept for provenance and is not run by the pipeline.
ROOT = Path("<research workspace>")
WORK = ROOT / "project1/JCIM_review_innovation_20261024"
sys.path.insert(0, str(WORK / "software/chemgridmap/src"))

from chemgridmap.core import _tie_inclusive_neighbor_indices  # noqa: E402
from chemgridmap.representations import morgan_fingerprints  # noqa: E402

LAYOUTS = ROOT / "chemgridmap/paper/output/revision_v6"
RESULTS = WORK / "results"
K = 10
TARGETS = {
    "chembl205": ("CA II", (7, 21, 42, 84, 126)),
    "chembl204": ("D3", (7, 21, 42)),
    "chembl240": ("hERG", (7, 21, 42)),
}


@lru_cache(maxsize=None)
def load_layout(target: str, seed: int) -> pd.DataFrame:
    data = pd.read_csv(LAYOUTS / target / f"coordinates_seed_{seed}.csv", low_memory=False)
    return data.sort_values("molecule_identity_key").reset_index(drop=True)


@lru_cache(maxsize=None)
def load_records(target: str) -> pd.DataFrame:
    """Retained source records restricted to the mapped entries of the paper."""
    records = pd.read_csv(LAYOUTS / target / f"{target}_retained_records.csv", low_memory=False)
    mapped = set(load_layout(target, 42).molecule_identity_key)
    return records[records.molecule_identity_key.isin(mapped)].reset_index(drop=True)


@lru_cache(maxsize=None)
def fingerprints(target: str) -> np.ndarray:
    return morgan_fingerprints(load_layout(target, 42).canonical_smiles.tolist()).astype(bool)


@lru_cache(maxsize=None)
def fingerprint_neighbors(target: str, k: int = K) -> tuple:
    return tuple(_tie_inclusive_neighbor_indices(fingerprints(target), k, metric="jaccard"))


def space_neighbors(points: np.ndarray, k: int = K) -> list:
    return _tie_inclusive_neighbor_indices(np.asarray(points, dtype=float), k)


def mean_recall(reference, displayed) -> float:
    return float(np.mean([len(set(r) & set(d)) / len(r) for r, d in zip(reference, displayed)]))


def per_molecule_recall(reference, displayed) -> np.ndarray:
    return np.array([len(set(r) & set(d)) / len(r) for r, d in zip(reference, displayed)])


def mean_jaccard(first, second) -> float:
    return float(np.mean([len(set(a) & set(b)) / len(set(a) | set(b)) for a, b in zip(first, second)]))


def grid_xy(layout: pd.DataFrame) -> np.ndarray:
    return layout[["grid_x", "grid_y"]].to_numpy(float)


def projection_xy(layout: pd.DataFrame) -> np.ndarray:
    return layout[["projection_x_raw", "projection_y_raw"]].to_numpy(float)
