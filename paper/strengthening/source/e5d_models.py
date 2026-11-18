"""Pure helpers for the E5d post hoc analysis."""
from __future__ import annotations

import numpy as np
import pandas as pd

from e7_analysis import conditional_agreement, two_stage_weights

SUPPORT_THRESHOLD = 0.8
BOOT_REPETITIONS = 2000
BOOT_SEED = 20261105
BIAS_BOOT_REPETITIONS = 500
BIAS_BOOT_SEED = 20261028


def add_rule_states(frame: pd.DataFrame) -> pd.DataFrame:
    """Add the unchanged G state and rule P's mutually exclusive display states."""
    out = frame.copy()
    support = out.support.to_numpy(float)
    touched = out.flag_touched.to_numpy(bool)
    out["g_status"] = np.where(
        support >= SUPPORT_THRESHOLD, "supported", "uncertain"
    )
    out["p_status"] = np.select(
        [
            touched,
            support < SUPPORT_THRESHOLD,
        ],
        [
            "uncertain_touched",
            "uncertain_low_support",
        ],
        default="supported",
    )
    return out


def b_straddle(values) -> int:
    """One unless the remaining units include both sides of the 6/7 class gap."""
    values = np.asarray(values, dtype=float)
    if not len(values):
        return 1
    return int(not (np.any(values <= 6.0) and np.any(values >= 7.0)))


def state_consistency(
    frame: pd.DataFrame,
    mask,
    weights: np.ndarray,
) -> dict[str, float | int]:
    """Point agreement and corrected two-stage percentile interval for a state."""
    mask = np.asarray(mask, dtype=bool)
    n = int(mask.sum())
    if not n:
        return {
            "n": 0,
            "agreement_rate": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "valid_repetitions": 0,
        }
    rate = float(frame.loc[mask, "agree"].astype(bool).mean())
    draws = conditional_agreement(frame, weights, mask)
    draws = draws[np.isfinite(draws)]
    if not len(draws):
        low = high = np.nan
    else:
        low, high = np.percentile(draws, [2.5, 97.5])
    return {
        "n": n,
        "agreement_rate": rate,
        "ci_low": float(low),
        "ci_high": float(high),
        "valid_repetitions": int(len(draws)),
    }


def corrected_two_stage_weights(
    frame: pd.DataFrame,
    repetitions: int = BOOT_REPETITIONS,
    seed: int = BOOT_SEED,
) -> np.ndarray:
    """Call the E7-R1 implementation directly rather than duplicating it."""
    return two_stage_weights(frame, seed=seed, repetitions=repetitions)


def bootstrap_bias_variance(
    frame: pd.DataFrame,
    display_index: int,
    repetitions: int = BIAS_BOOT_REPETITIONS,
    seed: int = BIAS_BOOT_SEED,
) -> tuple[float, int]:
    """Within-dataset entry bootstrap variance of a displayed-class bias."""
    selected = frame.loc[frame.displayed == display_index].copy()
    if selected.empty:
        return np.nan, 0
    selected["residual"] = (
        selected.p_displayed.to_numpy(float)
        - selected.agree.to_numpy(float)
    )
    moments = selected.groupby("entry", sort=True).residual.agg(
        ["sum", "count"]
    )
    n_entries = len(moments)
    if n_entries < 2:
        return 0.0, int(n_entries)
    rng = np.random.default_rng(seed)
    sums = moments["sum"].to_numpy(float)
    counts = moments["count"].to_numpy(float)
    draws = np.empty(repetitions, dtype=float)
    for start in range(0, repetitions, 100):
        stop = min(start + 100, repetitions)
        chosen = rng.integers(
            0, n_entries, size=(stop - start, n_entries)
        )
        draws[start:stop] = (
            sums[chosen].sum(axis=1)
            / counts[chosen].sum(axis=1)
        )
    return float(np.var(draws, ddof=1)), int(n_entries)
