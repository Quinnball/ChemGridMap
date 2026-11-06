"""Activity concentration units shared by record and molecule preparation."""
from __future__ import annotations

import numpy as np
import pandas as pd


UNIT_FACTORS = {"M": 1.0, "mM": 1e-3, "uM": 1e-6, "µM": 1e-6,
                "μM": 1e-6, "nM": 1e-9, "pM": 1e-12}


def convert_activity_values(values: pd.Series, *, value_unit: str = "p",
                            units: pd.Series | None = None):
    """Return pActivity values and one exclusion reason per invalid row.

    A per-row unit is a concentration unit, case-sensitive after stripping
    whitespace. The scalar ``p`` option preserves existing pActivity values.
    """
    if value_unit not in {"p", "M", "mM", "uM", "nM", "pM"}:
        raise ValueError("value_unit must be one of p, M, mM, uM, nM, pM.")
    if units is not None and value_unit != "p":
        raise ValueError("unit_column and value_unit cannot be used together.")
    numeric = pd.to_numeric(values, errors="coerce").astype(float)
    reasons = pd.Series("", index=values.index, dtype=object)
    reasons.loc[~np.isfinite(numeric)] = "invalid_activity_value"
    if units is None and value_unit == "p":
        return numeric.where(reasons.eq("")), reasons
    reasons.loc[reasons.eq("") & numeric.le(0)] = "nonpositive_activity_value"
    factors = (pd.Series(UNIT_FACTORS[value_unit], index=values.index)
               if units is None else units.fillna("").astype(str).str.strip().map(UNIT_FACTORS))
    reasons.loc[reasons.eq("") & factors.isna()] = "unrecognized_activity_unit"
    valid = reasons.eq("")
    converted = pd.Series(np.nan, index=values.index, dtype=float)
    # Separating logarithms also avoids underflow for very small concentrations.
    converted.loc[valid] = -np.log10(numeric.loc[valid]) - np.log10(factors.loc[valid])
    return converted, reasons


def validate_assumed_sigma(value):
    """Validate the user-specified between-unit SD, if supplied."""
    if value is not None and (not np.isfinite(value) or value <= 0):
        raise ValueError("assumed_sigma must be a finite number greater than 0.")
