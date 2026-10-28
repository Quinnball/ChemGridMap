"""Activity annotations shared by curation and evidence inspection."""
from __future__ import annotations

import numpy as np


RULE_KEYS = ("lower_threshold", "upper_threshold", "conflict_range_threshold")
RULE_COLUMNS = tuple("audit_" + key for key in RULE_KEYS)


def activity_class(value, lower_threshold, upper_threshold):
    if value <= lower_threshold:
        return "inactive"
    if value >= upper_threshold:
        return "active"
    return "medium"


def activity_annotations(values, *, lower_threshold, upper_threshold, conflict_range_threshold):
    rules = (lower_threshold, upper_threshold, conflict_range_threshold)
    if not np.isfinite(rules).all() or lower_threshold >= upper_threshold or conflict_range_threshold < 0:
        raise ValueError("Activity rules must be finite, ordered, and have a non-negative conflict range.")
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or not len(values) or not np.isfinite(values).all():
        raise ValueError("Activity measurements must be a non-empty finite vector.")
    median, low, high = float(np.median(values)), float(values.min()), float(values.max())
    inactive, active = low <= lower_threshold, high >= upper_threshold
    crossing = inactive and active
    large = high - low > conflict_range_threshold
    conflicted = crossing or large
    classes = {activity_class(value, lower_threshold, upper_threshold) for value in values}
    return {
        "activity_pchembl": median, "activity_class": activity_class(median, lower_threshold, upper_threshold),
        "n_records": len(values), "pchembl_min": low, "pchembl_max": high,
        "pchembl_median": median, "pchembl_std": float(np.std(values, ddof=1)) if len(values) > 1 else np.nan,
        "pchembl_range": high - low, "has_inactive_measurement": int(inactive),
        "has_active_measurement": int(active), "crosses_clean_boundary": int(crossing),
        "is_large_variation": int(large), "is_conflicted": int(conflicted),
        "has_class_boundary_crossing": int(len(classes) > 1), "n_measurement_classes": len(classes),
        "repeat_status": "singleton" if len(values) == 1 else "repeated_conflicted" if conflicted else "repeated_unflagged",
    }
