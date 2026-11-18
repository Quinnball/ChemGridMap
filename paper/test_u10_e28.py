"""Unit tests for the E28 aggregation (synthetic tables; no ChEMBL data needed)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import u10_neighbor_validation as v  # noqa: E402


def table(rows):
    return pd.DataFrame(rows, columns=["measured_class", "neighbor_evidence", "tanimoto", "nearest_reference_class"])


def test_dataset_metrics_follow_the_plan_definitions():
    rows = [
        # consistent: 3 active-evidence compounds (2 agree), 1 inactive-evidence compound (agrees)
        ("active", "active", 0.8, "active"), ("active", "active", 0.7, "active"), ("inactive", "active", 0.65, "inactive"),
        ("inactive", "inactive", 0.9, "medium"),
        # unsettled: nearest reference right once out of two
        ("medium", "insufficient", 0.7, "medium"), ("active", "mixed", 0.62, "inactive"),
        # none: nearest reference below 0.6, right once out of two
        ("active", "none", 0.3, "active"), ("inactive", "none", 0.2, "active"),
        # unmeasured compounds are ignored
        ("not_measured", "active", 0.9, "active"),
    ]
    m = v.dataset_metrics(table(rows))
    assert m["n"] == 8 and m["n_consistent"] == 4 and m["n_none"] == 2
    assert m["consistent_agreement"] == pytest.approx(3 / 4)
    # measured within C: active 2/4, inactive 2/4; evidence within C: active 3/4, inactive 1/4
    assert m["consistent_chance"] == pytest.approx(0.5 * 0.75 + 0.5 * 0.25)
    assert m["nearest_consistent_agreement"] == pytest.approx(3 / 4)
    assert m["nearest_unsettled_agreement"] == pytest.approx(1 / 2)
    assert m["nearest_none_agreement"] == pytest.approx(1 / 2)
    assert m["n2_difference"] == pytest.approx(0.25)
    assert m["n3_difference"] == pytest.approx(0.25)


def test_holm_is_monotone_and_capped():
    adjusted = v.holm(np.array([0.01, 0.04, 0.03]))
    assert adjusted.tolist() == pytest.approx([0.03, 0.06, 0.06])
    assert v.holm(np.array([0.6, 0.9])).max() <= 1.0


def test_wilson_interval_contains_the_estimate():
    low, high = v.wilson(8, 10)
    assert low < 0.8 < high
    assert np.isnan(v.wilson(0, 0)[0])
