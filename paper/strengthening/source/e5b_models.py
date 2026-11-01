"""Grid posteriors with NPMLE priors for the E5b post hoc analysis."""
from __future__ import annotations

import ctypes
import hashlib
import os
from pathlib import Path
import subprocess
import sys
from tempfile import gettempdir

sys.dont_write_bytecode = True

import numpy as np
from scipy import ndimage, stats

MU_GRID = np.arange(1.0, 13.0 + 0.01, 0.02, dtype=float)
SMOOTHING_SD = 0.15
EM_TOLERANCE = 1e-8
EM_MAX_ITER = 3000
STUDENT_DFS = (3, 4, 8)
ERROR_LABELS = ("gaussian", "t3", "t4", "t8")
NEW_METHODS = ("G-t", "G-NP", "G-NP-t", "G-NP-t3", "G-NP-t8")

_NATIVE_C_SOURCE = r"""
#include <math.h>
#include <stdlib.h>
#include <Accelerate/Accelerate.h>

void em_npmle(const float *ll, int n, int k, float *out_g, int max_iter,
              float tol, int *out_iter, float *out_delta) {
    float *g = malloc(k * sizeof(float));
    float *work = malloc(k * sizeof(float));
    float *ng = malloc(k * sizeof(float));
    float *lg = malloc(k * sizeof(float));
    for (int j = 0; j < k; j++) g[j] = 1.0f / k;
    for (int it = 1; it <= max_iter; it++) {
        for (int j = 0; j < k; j++)
            lg[j] = logf(g[j] > 0 ? g[j] : 1e-30f);
        vDSP_vclr(ng, 1, k);
        for (int e = 0; e < n; e++) {
            const float *row = ll + (size_t)e * k;
            vDSP_vadd(row, 1, lg, 1, work, 1, k);
            float maximum;
            vDSP_maxv(work, 1, &maximum, k);
            float negative = -maximum;
            vDSP_vsadd(work, 1, &negative, work, 1, k);
            int count = k;
            vvexpf(work, work, &count);
            float denominator;
            vDSP_sve(work, 1, &denominator, k);
            float inverse = 1.0f / denominator;
            vDSP_vsma(work, 1, &inverse, ng, 1, ng, 1, k);
        }
        float inverse_n = 1.0f / (float)n;
        vDSP_vsmul(ng, 1, &inverse_n, ng, 1, k);
        float *difference = work;
        vDSP_vsub(g, 1, ng, 1, difference, 1, k);
        float delta;
        vDSP_maxmgv(difference, 1, &delta, k);
        for (int j = 0; j < k; j++) g[j] = ng[j];
        *out_iter = it;
        *out_delta = delta;
        if (delta < tol) break;
    }
    for (int j = 0; j < k; j++) out_g[j] = g[j];
    free(g);
    free(work);
    free(ng);
    free(lg);
}
"""
_NATIVE_LIBRARY: ctypes.CDLL | None = None
_NATIVE_READY = False


def ensure_native_em_kernel() -> bool:
    """Build the local Accelerate EM kernel; NumPy remains the fallback."""
    global _NATIVE_LIBRARY, _NATIVE_READY
    if _NATIVE_READY:
        return _NATIVE_LIBRARY is not None
    _NATIVE_READY = True
    if sys.platform != "darwin":
        return False
    digest = hashlib.sha256(_NATIVE_C_SOURCE.encode()).hexdigest()[:12]
    library_path = (
        Path(gettempdir()) / f"e5b_em_{digest}_{os.getpid()}.dylib"
    )
    if not library_path.exists():
        command = [
            "clang",
            "-O3",
            "-march=native",
            "-ffast-math",
            "-dynamiclib",
            "-fPIC",
            "-x",
            "c",
            "-framework",
            "Accelerate",
            "-o",
            str(library_path),
            "-",
        ]
        result = subprocess.run(
            command,
            input=_NATIVE_C_SOURCE,
            text=True,
            capture_output=True,
            check=False,
        )
        if result.returncode != 0:
            print(f"Accelerate EM kernel unavailable: {result.stderr}", flush=True)
            return False
    library = ctypes.CDLL(str(library_path))
    function = library.em_npmle
    function.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.c_float,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_float),
    ]
    function.restype = None
    _NATIVE_LIBRARY = library
    return True


def student_t_scale(sigma: float, df: int) -> float:
    """Scale a Student-t density so its variance is sigma squared."""
    return float(sigma) * np.sqrt((int(df) - 2.0) / int(df))


def _density_parameters(sigma: float, df: float) -> tuple[float, float]:
    if np.isinf(df):
        return float(sigma), np.inf
    return student_t_scale(sigma, int(df)), float(df)


def _log_density(values, center, sigma: float, df: float) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    center = np.asarray(center, dtype=float)
    scale, degrees = _density_parameters(sigma, df)
    if np.isinf(degrees):
        return stats.norm.logpdf(values, loc=center, scale=scale)
    return stats.t.logpdf(values, degrees, loc=center, scale=scale)


def entry_likelihood(values_by_entry, sigma: float, df: float) -> np.ndarray:
    """Log likelihood on the fixed grid for each training entry."""
    rows = np.empty((len(values_by_entry), len(MU_GRID)), dtype=float)
    for index, values in enumerate(values_by_entry):
        values = np.asarray(values, dtype=float)
        rows[index] = _log_density(
            values[:, None], MU_GRID[None, :], sigma, df
        ).sum(axis=0)
    return rows


def em_npmle(log_likelihood: np.ndarray) -> tuple[np.ndarray, int, float]:
    """Uniform-start EM for the fixed-grid nonparametric prior."""
    log_likelihood = np.asarray(log_likelihood, dtype=float)
    n_entries, n_grid = log_likelihood.shape
    if ensure_native_em_kernel():
        matrix = np.ascontiguousarray(log_likelihood, dtype=np.float32)
        prior = np.empty(n_grid, dtype=np.float32)
        iterations = ctypes.c_int()
        delta = ctypes.c_float()
        _NATIVE_LIBRARY.em_npmle(
            matrix.ctypes.data,
            n_entries,
            n_grid,
            prior.ctypes.data,
            EM_MAX_ITER,
            EM_TOLERANCE,
            ctypes.byref(iterations),
            ctypes.byref(delta),
        )
        return prior.astype(float), iterations.value, float(delta.value)
    prior = np.full(n_grid, 1.0 / n_grid)
    last_delta = np.inf
    iterations = 0
    for iterations in range(1, EM_MAX_ITER + 1):
        with np.errstate(divide="ignore"):
            joint = log_likelihood + np.log(prior)[None, :]
        maximum = joint.max(axis=1, keepdims=True)
        responsibility = np.exp(joint - maximum)
        responsibility /= responsibility.sum(axis=1, keepdims=True)
        updated = responsibility.mean(axis=0)
        last_delta = float(np.max(np.abs(updated - prior)))
        prior = updated
        if last_delta < EM_TOLERANCE:
            break
    return prior, iterations, last_delta


def smooth_prior(prior: np.ndarray) -> np.ndarray:
    """Gaussian-kernel smooth on the 0.02 pChEMBL grid, then normalize."""
    smoothed = ndimage.gaussian_filter1d(
        np.asarray(prior, dtype=float),
        sigma=SMOOTHING_SD / 0.02,
        mode="constant",
        truncate=5.0,
    )
    return smoothed / smoothed.sum()


def fit_npmle(
    values_by_entry,
    sigma: float,
    df: float,
) -> dict[str, float | int | np.ndarray]:
    """Fit and smooth one target-fold NPMLE on training entries only."""
    likelihood = entry_likelihood(values_by_entry, sigma, df)
    raw, iterations, delta = em_npmle(likelihood)
    prior = smooth_prior(raw)
    cumulative = np.cumsum(prior)
    quantiles = np.interp(
        [0.05, 0.25, 0.50, 0.75, 0.95],
        cumulative - 0.5 * prior,
        MU_GRID,
    )
    return {
        "prior": prior,
        "iterations": iterations,
        "last_delta": delta,
        "metric_weights": likelihood,
        **{
            f"q{name}": float(value)
            for name, value in zip(
                ("05", "25", "50", "75", "95"), quantiles
            )
        },
    }


def fit_error_priors(
    values_by_entry,
    sigma: float,
) -> dict[str, dict]:
    """Fit Gaussian and t3/t4/t8 NPMLE priors from one training fold."""
    return {
        label: fit_npmle(values_by_entry, sigma, df)
        for label, df in zip(ERROR_LABELS, (np.inf, 3.0, 4.0, 8.0))
    }


def flat_prior() -> np.ndarray:
    return np.full(len(MU_GRID), 1.0 / len(MU_GRID))


def _grid_class_probabilities(sigma: float, df: float) -> np.ndarray:
    scale, degrees = _density_parameters(sigma, df)
    low = stats.norm.cdf(
        (6.0 - MU_GRID) / scale
    )
    high = stats.norm.sf(
        (7.0 - MU_GRID) / scale
    )
    if not np.isinf(degrees):
        low = stats.t.cdf((6.0 - MU_GRID) / scale, degrees)
        high = stats.t.sf((7.0 - MU_GRID) / scale, degrees)
    medium = np.clip(1.0 - low - high, 0.0, 1.0)
    return np.column_stack([low, medium, high])


def _integrated_support(posterior: np.ndarray) -> np.ndarray:
    """Integrate posterior density over the three displayed intervals."""
    lower_index = int(np.argmin(np.abs(MU_GRID - 6.0)))
    upper_index = int(np.argmin(np.abs(MU_GRID - 7.0)))
    parts = np.array(
        [
            np.trapezoid(
                posterior[: lower_index + 1],
                MU_GRID[: lower_index + 1],
            ),
            np.trapezoid(
                posterior[lower_index : upper_index + 1],
                MU_GRID[lower_index : upper_index + 1],
            ),
            np.trapezoid(
                posterior[upper_index:],
                MU_GRID[upper_index:],
            ),
        ],
        dtype=float,
    )
    return parts / parts.sum()


def grid_prediction(
    values,
    sigma: float,
    df: float,
    prior: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, float, np.ndarray]:
    """Return support, predictive classes, posterior mean and prior."""
    values = np.asarray(values, dtype=float)
    if not len(values):
        raise ValueError("Grid posterior needs at least one remaining value.")
    if prior is None:
        prior = flat_prior()
    prior = np.asarray(prior, dtype=float)
    with np.errstate(divide="ignore"):
        joint = (
            _log_density(
                values[:, None], MU_GRID[None, :], sigma, df
            ).sum(axis=0)
            + np.log(prior)
        )
    posterior = np.exp(joint - joint.max())
    posterior /= posterior.sum()
    support = _integrated_support(posterior)
    predictive = posterior @ _grid_class_probabilities(sigma, df)
    posterior_mean = float(posterior @ MU_GRID)
    return support, predictive, posterior_mean, prior


def entry_values(frame, excluded_entries=()) -> list[np.ndarray]:
    """Return entry arrays after excluding specified entries."""
    excluded = set(excluded_entries)
    selected = frame.loc[~frame.entry.isin(excluded)]
    return [
        group.document_value.to_numpy(float)
        for _, group in selected.groupby("entry", sort=True)
    ]
