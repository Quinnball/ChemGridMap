"""Fingerprint-neighbour refinement of one-molecule-per-cell grid layouts.

Minimum-displacement assignment keeps molecules close to their projected
positions but ignores which molecules are neighbours in the fingerprint space.
This refinement starts from that assignment and applies greedy local moves:
a molecule may move to an empty cell or swap with another molecule when the
exchange brings fingerprint neighbours closer on the lattice. A Gaussian
saturation bounds the influence of distant neighbours, and a weak quadratic
anchor keeps every molecule near its projected position so that the global
arrangement of the projection is retained.
"""
from __future__ import annotations

from typing import Dict, Optional, Sequence

import numpy as np

REFINEMENT_DEFAULTS = {
    "sigmas": (3.0, 2.0, 1.5),
    "anchor": 0.003,
    "radius": 3,
    "jump_radius": 1,
    "max_sweeps": 40,
}


def symmetric_neighbor_lists(neighbor_lists: Sequence[Sequence[int]]) -> list:
    """Return undirected adjacency lists from directed neighbour lists."""
    sets = [set(int(j) for j in row) for row in neighbor_lists]
    for i, row in enumerate(neighbor_lists):
        for j in row:
            sets[int(j)].add(i)
    for i, row in enumerate(sets):
        row.discard(i)
    return [np.fromiter(sorted(row), dtype=int, count=len(row)) for row in sets]


def _disk(radius: int) -> np.ndarray:
    radius = int(radius)
    return np.array([(a, b) for a in range(-radius, radius + 1) for b in range(-radius, radius + 1)
                     if a * a + b * b <= radius * radius], dtype=int).reshape(-1, 2)


def neighbor_objective(cells: np.ndarray, adjacency: Sequence[np.ndarray], sigma: float) -> float:
    """Sum of saturated squared lattice distances over undirected neighbour edges."""
    cells = np.asarray(cells, dtype=float)
    total = 0.0
    for i, row in enumerate(adjacency):
        row = row[row > i]
        if row.size:
            total += float((1.0 - np.exp(-((cells[row] - cells[i]) ** 2).sum(1) / (2.0 * sigma ** 2))).sum())
    return total


def refine_grid_layout(
    cells: np.ndarray,
    projection: np.ndarray,
    neighbor_lists: Sequence[Sequence[int]],
    side: int,
    *,
    sigmas: Sequence[float] = REFINEMENT_DEFAULTS["sigmas"],
    anchor: float = REFINEMENT_DEFAULTS["anchor"],
    radius: int = REFINEMENT_DEFAULTS["radius"],
    jump_radius: int = REFINEMENT_DEFAULTS["jump_radius"],
    max_sweeps: int = REFINEMENT_DEFAULTS["max_sweeps"],
    random_state: int = 0,
    diagnostics: Optional[Dict] = None,
) -> np.ndarray:
    """Refine integer lattice cells so that fingerprint neighbours sit closer together.

    ``cells`` and ``projection`` are (n, 2) arrays in lattice units (column, row);
    ``projection`` is the normalized projection scaled to the same lattice.
    ``neighbor_lists`` are directed representation-space neighbours; they are
    symmetrized before optimization. Returns refined integer cells (one molecule
    per cell). A negative ``jump_radius`` disables neighbour-centroid moves.
    """
    x = np.asarray(cells, dtype=int).copy()
    p = np.asarray(projection, dtype=float)
    n = len(x)
    if x.shape != (n, 2) or p.shape != (n, 2) or len(neighbor_lists) != n:
        raise ValueError("cells, projection and neighbour lists must describe the same molecules.")
    if side < 2 or x.min() < 0 or x.max() >= side:
        raise ValueError("cells must lie on the side x side lattice.")
    if len({tuple(row) for row in x}) != n:
        raise ValueError("Initial cells must be unique.")
    if anchor < 0 or radius < 1 or not len(sigmas) or min(sigmas) <= 0:
        raise ValueError("Refinement needs a non-negative anchor, radius >= 1 and positive sigmas.")
    adjacency = symmetric_neighbor_lists(neighbor_lists)
    occupant = -np.ones((side, side), dtype=int)
    occupant[x[:, 0], x[:, 1]] = np.arange(n)
    local, jump = _disk(radius), _disk(max(jump_radius, 0))
    local = local[(local != 0).any(1)]
    rng = np.random.default_rng(random_state)
    start = {sigma: neighbor_objective(x, adjacency, sigma) for sigma in sigmas}
    history = []
    for sigma in sigmas:
        scale = 2.0 * float(sigma) ** 2
        for sweep in range(int(max_sweeps)):
            moves = 0
            for i in rng.permutation(n):
                xi, ni = x[i], adjacency[i]
                xn = x[ni]
                parts = [xi + local]
                if jump_radius >= 0 and ni.size:
                    parts.append(np.rint(np.median(xn, axis=0)).astype(int) + jump)
                cand = np.unique(np.vstack(parts), axis=0)
                cand = cand[(cand >= 0).all(1) & (cand < side).all(1) & (cand != xi).any(1)]
                if not len(cand):
                    continue
                old_i = 1.0 - np.exp(-((xi - xn) ** 2).sum(1) / scale)
                new_i = 1.0 - np.exp(-((cand[:, None, :] - xn[None]) ** 2).sum(-1) / scale)
                delta = (new_i.sum(1) - old_i.sum()
                         + anchor * (((cand - p[i]) ** 2).sum(1) - ((xi - p[i]) ** 2).sum()))
                occupants = occupant[cand[:, 0], cand[:, 1]]
                for idx in np.flatnonzero(occupants >= 0):
                    j = occupants[idx]
                    nj = adjacency[j]
                    nj = nj[nj != i]
                    xnj = x[nj]
                    delta[idx] += (
                        (1.0 - np.exp(-((xi - xnj) ** 2).sum(1) / scale)).sum()
                        - (1.0 - np.exp(-((cand[idx] - xnj) ** 2).sum(1) / scale)).sum()
                        + anchor * (((xi - p[j]) ** 2).sum() - ((cand[idx] - p[j]) ** 2).sum()))
                    shared = np.flatnonzero(ni == j)
                    if shared.size:  # a swapped neighbour pair keeps its distance
                        delta[idx] += old_i[shared[0]] - new_i[idx, shared[0]]
                best = int(np.argmin(delta))
                if delta[best] < -1e-9:
                    target = cand[best]
                    j = occupant[target[0], target[1]]
                    occupant[xi[0], xi[1]] = j
                    if j >= 0:
                        x[j] = xi.copy()
                    occupant[target[0], target[1]] = i
                    x[i] = target
                    moves += 1
            history.append({"sigma": float(sigma), "sweep": sweep + 1, "moves": moves})
            if moves == 0:
                break
    if len({tuple(row) for row in x}) != n:
        raise RuntimeError("Refinement produced duplicate cells.")
    if diagnostics is not None:
        diagnostics.update({
            "method": "fingerprint_neighbor_local_swaps",
            "sigmas": [float(s) for s in sigmas], "anchor": float(anchor), "radius": int(radius),
            "jump_radius": int(jump_radius), "random_state": int(random_state),
            "max_sweeps": int(max_sweeps),
            "undirected_edges": int(sum(len(row) for row in adjacency) // 2),
            "objective_start": {str(k): v for k, v in start.items()},
            "objective_end": {str(s): neighbor_objective(x, adjacency, s) for s in sigmas},
            "sweeps": history, "moved_molecules": int((x != np.asarray(cells)).any(1).sum()),
        })
    return x
