"""Offset-aware support models for the E5c post hoc analysis."""
from __future__ import annotations

import ctypes
from itertools import combinations
import hashlib
import os
from pathlib import Path
import subprocess
import sys
from tempfile import gettempdir

sys.dont_write_bytecode = True

import numpy as np
from scipy import ndimage, stats

from e5b_models import MU_GRID, SMOOTHING_SD, student_t_scale

ANALYSIS_ROLE = "post hoc"
STUDENT_DF = 4
EM_RELATIVE_TOLERANCE = 1e-9
EM_MAX_ITER = 5000
CLASS_BOUNDARIES = (6.0, 7.0)


def weighted_gaussian_prediction(
    values,
    unit_variances,
    base_sigma: float,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """G-F weighted merge with one variance per remaining unit."""
    values = np.asarray(values, dtype=float)
    unit_variances = np.asarray(unit_variances, dtype=float)
    if not len(values) or len(values) != len(unit_variances):
        raise ValueError("G-F needs one positive variance per unit.")
    if np.any(unit_variances <= 0):
        raise ValueError("Unit variances must be positive.")
    weights = 1.0 / unit_variances
    posterior_mean = float(np.sum(weights * values) / weights.sum())
    posterior_variance = float(1.0 / weights.sum())
    support = _normal_interval_probabilities(
        posterior_mean, posterior_variance
    )
    predictive = _normal_interval_probabilities(
        posterior_mean,
        posterior_variance + float(base_sigma) ** 2,
    )
    return (
        support,
        predictive,
        posterior_mean,
        posterior_variance,
    )


def _normal_interval_probabilities(
    mean: float, variance: float
) -> np.ndarray:
    scale = float(np.sqrt(variance))
    low = float(stats.norm.cdf((CLASS_BOUNDARIES[0] - mean) / scale))
    high = float(stats.norm.sf((CLASS_BOUNDARIES[1] - mean) / scale))
    return np.array([low, max(0.0, 1.0 - low - high), high])


def _t_log_density(
    values, centers, sigmas, df: int = STUDENT_DF
) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    centers = np.asarray(centers, dtype=float)
    scales = np.asarray(
        [student_t_scale(float(value), int(df)) for value in sigmas],
        dtype=float,
    )
    return stats.t.logpdf(
        values[:, None],
        int(df),
        loc=centers,
        scale=scales[:, None],
    )


def weighted_entry_likelihood(
    values_by_entry,
    sigmas_by_entry,
    df: int = STUDENT_DF,
) -> np.ndarray:
    """Fixed-grid likelihood for training entries with unit-specific scales."""
    if len(values_by_entry) != len(sigmas_by_entry):
        raise ValueError("Each training entry needs a matching sigma array.")
    rows = np.empty((len(values_by_entry), len(MU_GRID)), dtype=float)
    for index, (values, sigmas) in enumerate(
        zip(values_by_entry, sigmas_by_entry)
    ):
        rows[index] = _t_log_density(
            values, MU_GRID, sigmas, df
        ).sum(axis=0)
    return rows


_NATIVE_C_SOURCE = r"""
#include <math.h>
#include <stdlib.h>
#include <pthread.h>

typedef struct {
    const double *ll;
    const double *log_g;
    double *next;
    double *work;
    double *total_ll;
    int n;
    int k;
    int thread_id;
    int nthreads;
} worker_state;

typedef struct {
    worker_state *state;
    int thread_id;
} worker_arg;

void *em_worker(void *raw) {
    worker_arg *arg = (worker_arg *)raw;
    worker_state *state = arg->state;
    int thread_id = arg->thread_id;
    double *next = state->next + (size_t)thread_id * state->k;
    double *work = state->work + (size_t)thread_id * state->k;
    for (int j = 0; j < state->k; j++) next[j] = 0.0;
    double total_ll = 0.0;
    for (int e = thread_id; e < state->n; e += state->nthreads) {
        const double *row = state->ll + (size_t)e * state->k;
        double maximum = -INFINITY;
        for (int j = 0; j < state->k; j++) {
            work[j] = row[j] + state->log_g[j];
            if (work[j] > maximum) maximum = work[j];
        }
        double denominator = 0.0;
        for (int j = 0; j < state->k; j++) {
            work[j] = exp(work[j] - maximum);
            denominator += work[j];
        }
        total_ll += maximum + log(denominator);
        double inverse = 1.0 / denominator;
        for (int j = 0; j < state->k; j++)
            next[j] += work[j] * inverse;
    }
    state->total_ll[thread_id] = total_ll;
    return NULL;
}

void em_npmle_relative(
    const double *ll, int n, int k, const double *initial_g,
    double *out_g, int max_iter,
    int nthreads, double tol, int *out_iter, double *out_relative_change,
    double *out_delta) {
    double *g = malloc(k * sizeof(double));
    double *next = malloc(k * sizeof(double));
    double *log_g = malloc(k * sizeof(double));
    double *all_next = calloc((size_t)nthreads * k, sizeof(double));
    double *all_work = malloc((size_t)nthreads * k * sizeof(double));
    double *all_total = malloc((size_t)nthreads * sizeof(double));
    pthread_t *threads = malloc((size_t)nthreads * sizeof(pthread_t));
    worker_arg *args = malloc((size_t)nthreads * sizeof(worker_arg));
    worker_state state;
    state.ll = ll;
    state.log_g = log_g;
    state.next = all_next;
    state.work = all_work;
    state.total_ll = all_total;
    state.n = n;
    state.k = k;
    state.nthreads = nthreads;
    for (int j = 0; j < k; j++) g[j] = initial_g[j];
    double old_ll = 0.0;
    for (int j = 0; j < k; j++)
        log_g[j] = log(g[j] > 0 ? g[j] : 1e-300);
    for (int e = 0; e < n; e++) {
        const double *row = ll + (size_t)e * k;
        double maximum = -INFINITY;
        for (int j = 0; j < k; j++) {
            double value = row[j] + log_g[j];
            if (value > maximum) maximum = value;
        }
        double denominator = 0.0;
        for (int j = 0; j < k; j++)
            denominator += exp(row[j] + log_g[j] - maximum);
        old_ll += maximum + log(denominator);
    }
    old_ll /= (double)n;
    double relative = INFINITY;
    double delta = INFINITY;
    for (int it = 1; it <= max_iter; it++) {
        for (int j = 0; j < k; j++)
            log_g[j] = log(g[j] > 0 ? g[j] : 1e-300);
        for (int t = 0; t < nthreads; t++) {
            args[t].state = &state;
            args[t].thread_id = t;
            pthread_create(&threads[t], NULL, em_worker, &args[t]);
        }
        for (int t = 0; t < nthreads; t++) pthread_join(threads[t], NULL);
        double total_ll = 0.0;
        for (int j = 0; j < k; j++) next[j] = 0.0;
        for (int t = 0; t < nthreads; t++) {
            total_ll += all_total[t];
            for (int j = 0; j < k; j++)
                next[j] += all_next[(size_t)t * k + j];
        }
        double inverse_n = 1.0 / (double)n;
        delta = 0.0;
        for (int j = 0; j < k; j++) {
            next[j] *= inverse_n;
            double difference = fabs(next[j] - g[j]);
            if (difference > delta) delta = difference;
        }
        total_ll *= inverse_n;
        double denominator = fmax(fabs(old_ll), 1e-12);
        relative = fabs(total_ll - old_ll) / denominator;
        old_ll = total_ll;
        for (int j = 0; j < k; j++) g[j] = next[j];
        *out_iter = it;
        *out_relative_change = relative;
        *out_delta = delta;
        if (isfinite(relative) && relative < tol) break;
    }
    for (int j = 0; j < k; j++) out_g[j] = g[j];
    free(g);
    free(next);
    free(log_g);
    free(all_next);
    free(all_work);
    free(all_total);
    free(threads);
    free(args);
}
"""
_NATIVE_LIBRARY: ctypes.CDLL | None = None
_NATIVE_READY = False


def ensure_native_em_kernel() -> bool:
    """Build the local double-precision EM kernel."""
    global _NATIVE_LIBRARY, _NATIVE_READY
    if _NATIVE_READY:
        return _NATIVE_LIBRARY is not None
    _NATIVE_READY = True
    if sys.platform != "darwin":
        return False
    digest = hashlib.sha256(_NATIVE_C_SOURCE.encode()).hexdigest()[:12]
    library_path = (
        Path(gettempdir()) / f"e5c_em_{digest}_{os.getpid()}.dylib"
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
            print(
                f"Native E5c EM kernel unavailable: {result.stderr}",
                flush=True,
            )
            return False
    library = ctypes.CDLL(str(library_path))
    function = library.em_npmle_relative
    function.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_void_p,
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_double,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
    ]
    function.restype = None
    _NATIVE_LIBRARY = library
    return True


def _numpy_em(
    log_likelihood: np.ndarray,
    initial_prior: np.ndarray | None = None,
) -> tuple[np.ndarray, int, float, float]:
    if initial_prior is None:
        prior = np.full(
            log_likelihood.shape[1],
            1.0 / log_likelihood.shape[1],
            dtype=float,
        )
    else:
        prior = np.asarray(initial_prior, dtype=float).copy()
        prior /= prior.sum()
    initial_joint = log_likelihood + np.log(prior)[None, :]
    initial_maximum = initial_joint.max(axis=1)
    old_ll = float(
        np.mean(
            initial_maximum
            + np.log(
                np.exp(initial_joint - initial_maximum[:, None]).sum(axis=1)
            )
        )
    )
    relative = np.inf
    delta = np.inf
    for iteration in range(1, EM_MAX_ITER + 1):
        with np.errstate(divide="ignore"):
            joint = log_likelihood + np.log(prior)[None, :]
        maximum = joint.max(axis=1, keepdims=True)
        exponential = np.exp(joint - maximum)
        denominator = exponential.sum(axis=1, keepdims=True)
        responsibility = exponential / denominator
        updated = responsibility.mean(axis=0)
        delta = float(np.max(np.abs(updated - prior)))
        log_likelihood_value = float(
            np.mean(maximum[:, 0] + np.log(denominator[:, 0]))
        )
        relative = abs(log_likelihood_value - old_ll) / max(
            abs(old_ll), 1e-12
        )
        old_ll = log_likelihood_value
        prior = updated
        if np.isfinite(relative) and relative < EM_RELATIVE_TOLERANCE:
            break
    return prior, iteration, float(relative), delta


def em_npmle_relative(
    log_likelihood: np.ndarray,
    initial_prior: np.ndarray | None = None,
) -> tuple[np.ndarray, int, float, float]:
    """EM with the task's relative log-likelihood criterion."""
    log_likelihood = np.asarray(log_likelihood, dtype=float)
    if initial_prior is None:
        initial_prior = np.full(
            log_likelihood.shape[1],
            1.0 / log_likelihood.shape[1],
            dtype=float,
        )
    else:
        initial_prior = np.asarray(initial_prior, dtype=float).copy()
        initial_prior /= initial_prior.sum()
    if ensure_native_em_kernel():
        matrix = np.ascontiguousarray(log_likelihood, dtype=np.float64)
        initial = np.ascontiguousarray(
            initial_prior, dtype=np.float64
        )
        prior = np.empty(log_likelihood.shape[1], dtype=np.float64)
        iterations = ctypes.c_int()
        relative = ctypes.c_double()
        delta = ctypes.c_double()
        nthreads = max(
            1, min(4, int(os.environ.get("E5C_EM_THREADS", "2")))
        )
        _NATIVE_LIBRARY.em_npmle_relative(
            matrix.ctypes.data,
            log_likelihood.shape[0],
            log_likelihood.shape[1],
            initial.ctypes.data,
            prior.ctypes.data,
            EM_MAX_ITER,
            nthreads,
            EM_RELATIVE_TOLERANCE,
            ctypes.byref(iterations),
            ctypes.byref(relative),
            ctypes.byref(delta),
        )
        return (
            prior.astype(float),
            iterations.value,
            float(relative.value),
            float(delta.value),
        )
    return _numpy_em(log_likelihood, initial_prior)


def smooth_prior(prior: np.ndarray) -> np.ndarray:
    smoothed = ndimage.gaussian_filter1d(
        np.asarray(prior, dtype=float),
        sigma=SMOOTHING_SD / 0.02,
        mode="constant",
        truncate=5.0,
    )
    return smoothed / smoothed.sum()


def fit_npmle_weighted(
    values_by_entry,
    sigmas_by_entry,
    df: int = STUDENT_DF,
    initial_prior: np.ndarray | None = None,
) -> dict:
    """Fit one smoothed weighted t NPMLE from training entries."""
    likelihood = weighted_entry_likelihood(
        values_by_entry, sigmas_by_entry, df
    )
    raw, iterations, relative, delta = em_npmle_relative(
        likelihood, initial_prior
    )
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
        "last_relative_change": relative,
        "last_g_delta": delta,
        "converged": bool(relative < EM_RELATIVE_TOLERANCE),
        **{
            f"q{name}": float(value)
            for name, value in zip(
                ("05", "25", "50", "75", "95"), quantiles
            )
        },
    }


def _integrated_posterior(posterior: np.ndarray) -> np.ndarray:
    lower_index = int(np.argmin(np.abs(MU_GRID - CLASS_BOUNDARIES[0])))
    upper_index = int(np.argmin(np.abs(MU_GRID - CLASS_BOUNDARIES[1])))
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
            np.trapezoid(posterior[upper_index:], MU_GRID[upper_index:]),
        ],
        dtype=float,
    )
    return parts / parts.sum()


def weighted_grid_prediction(
    values,
    sigmas,
    base_sigma: float,
    prior: np.ndarray,
    df: int = STUDENT_DF,
) -> tuple[np.ndarray, np.ndarray, float, np.ndarray]:
    """G-NP-t-F posterior with one evidence scale per remaining unit."""
    values = np.asarray(values, dtype=float)
    sigmas = np.asarray(sigmas, dtype=float)
    if not len(values) or len(values) != len(sigmas):
        raise ValueError("Grid prediction needs one sigma per unit.")
    prior = np.asarray(prior, dtype=float)
    with np.errstate(divide="ignore"):
        joint = (
            _t_log_density(values, MU_GRID, sigmas, df).sum(axis=0)
            + np.log(prior)
        )
    posterior = np.exp(joint - joint.max())
    posterior /= posterior.sum()
    support = _integrated_posterior(posterior)
    predictive = posterior @ _t_interval_probabilities(base_sigma, df)
    posterior_mean = float(posterior @ MU_GRID)
    return support, predictive, posterior_mean, prior


def _t_interval_probabilities(sigma: float, df: int) -> np.ndarray:
    scale = student_t_scale(float(sigma), int(df))
    low = stats.t.cdf((CLASS_BOUNDARIES[0] - MU_GRID) / scale, int(df))
    high = stats.t.sf((CLASS_BOUNDARIES[1] - MU_GRID) / scale, int(df))
    medium = np.clip(1.0 - low - high, 0.0, 1.0)
    return np.column_stack([low, medium, high])


def delta_document_map(
    observations,
    training_entries: set[str] | None = None,
    expand_documents: bool = False,
) -> dict[str, float]:
    """Find flagged document pairs and retain each document's largest |median|."""
    selected = observations
    if training_entries is not None:
        selected = selected.loc[selected.entry.isin(training_entries)]
    delta: dict[str, float] = {}
    if selected.empty:
        return delta
    for _, group in selected.groupby(
        ["document_a", "document_b"], sort=True
    ):
        if len(group) < 3:
            continue
        median = float(np.median(group.difference))
        sign = float(np.sign(median))
        agreement = float(
            np.mean(group.difference.to_numpy(float) * sign > 0)
        )
        if abs(median) >= 1.0 and agreement >= 0.8:
            document_a = str(group.document_a.iloc[0])
            document_b = str(group.document_b.iloc[0])
            value = abs(median)
            documents_a = (
                document_a.split(";") if expand_documents else [document_a]
            )
            documents_b = (
                document_b.split(";") if expand_documents else [document_b]
            )
            for document in documents_a:
                delta[document] = max(delta.get(document, 0.0), value)
            for document in documents_b:
                delta[document] = max(delta.get(document, 0.0), value)
    return delta


def pair_observations_from_units(units) -> "pandas.DataFrame":
    """Recompute pair observations after a design-specific unit merge/filter."""
    import pandas as pd

    rows = []
    for entry, group in units.groupby("entry", sort=True):
        records = [
            (
                str(row.documents),
                set(str(row.documents).split(";")),
                float(row.document_value),
            )
            for row in group.itertuples(index=False)
        ]
        for (documents_a, _set_a, value_a), (
            documents_b,
            _set_b,
            value_b,
        ) in combinations(records, 2):
            rows.append(
                {
                    "document_a": documents_a,
                    "document_b": documents_b,
                    "entry": str(entry),
                    "difference": value_b - value_a,
                }
            )
    return pd.DataFrame(
        rows,
        columns=["document_a", "document_b", "entry", "difference"],
    )


def delta_document_map_from_units(units) -> dict[str, float]:
    """Flag merged units and propagate the pair offset to their documents."""
    observations = pair_observations_from_units(units)
    return delta_document_map(
        observations, expand_documents=True
    )


def unit_deltas(
    units,
    delta_by_document: dict[str, float],
    fixed_delta: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return each unit's extra SD and whether a training flag touches it."""
    values = []
    touched = []
    for row in units.itertuples(index=False):
        documents = str(row.documents).split(";")
        matched = [
            delta_by_document[document]
            for document in documents
            if document in delta_by_document
        ]
        is_touched = bool(matched)
        touched.append(is_touched)
        if not is_touched:
            values.append(0.0)
        elif fixed_delta is None:
            values.append(float(max(matched)))
        else:
            values.append(float(fixed_delta))
    return (
        np.asarray(values, dtype=float),
        np.asarray(touched, dtype=bool),
    )


def unit_sigmas(delta, base_sigma: float) -> np.ndarray:
    delta = np.asarray(delta, dtype=float)
    return np.sqrt(float(base_sigma) ** 2 + delta**2)


def entry_values_and_sigmas(units, base_sigma: float) -> tuple[list, list]:
    """Group marked training units by entry for weighted NPMLE fitting."""
    values_by_entry = []
    sigmas_by_entry = []
    for _, group in units.groupby("entry", sort=True):
        values_by_entry.append(
            group.document_value.to_numpy(float)
        )
        sigmas_by_entry.append(
            unit_sigmas(group.delta.to_numpy(float), base_sigma)
        )
    return values_by_entry, sigmas_by_entry
