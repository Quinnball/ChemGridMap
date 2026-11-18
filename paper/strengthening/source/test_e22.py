import pandas as pd
import pytest

from e22_class_reliability import (
    BIAS_CATEGORIES,
    BIAS_HOLM_FAMILY,
    E5C,
    OUT,
    add_support_status,
    bin_composition,
    class_mask,
    reliability_tests,
    validate_bin_composition,
    validate_expected,
)


def synthetic_predictions() -> pd.DataFrame:
    rows = []
    for endpoint in ("IC50", "Ki"):
        for target in ("T1", "T2"):
            for displayed, support, agreement in (
                (0, 0.9, 1),
                (0, 0.7, 0),
                (1, 0.9, 1),
                (1, 0.7, 0),
                (2, 0.9, 1),
                (2, 0.7, 0),
            ):
                for copy in range(5):
                    rows.append(
                        {
                            "target": target,
                            "endpoint": endpoint,
                            "entry": f"{target}-{displayed}-{copy}",
                            "displayed": displayed,
                            "support": support,
                            "p_displayed": 0.9 if support >= 0.8 else 0.4,
                            "agree": bool(agreement),
                        }
                    )
    frame = pd.DataFrame(rows)
    return add_support_status(frame, "E5c")


def test_support_status_and_class_split() -> None:
    frame = synthetic_predictions()
    assert set(frame.e22_status) == {"supported", "uncertain"}
    assert class_mask(frame, "inactive").sum() == 40
    assert class_mask(frame, "active").sum() == 40
    assert class_mask(frame, "inactive+active").sum() == 80
    assert class_mask(frame, "all").sum() == 120


def test_holm_family_is_six_tests_per_design() -> None:
    assert BIAS_CATEGORIES == ("inactive", "medium", "active")
    assert "all" not in BIAS_CATEGORIES
    assert BIAS_HOLM_FAMILY == "inactive, medium, active x IC50, Ki"
    frame = synthetic_predictions()
    designs = {name: frame for name in ("entry", "document", "near_duplicate", "time_2015", "E11")}
    tests = reliability_tests(designs)
    assert len(tests) == 30
    for _, group in tests.groupby("design"):
        assert len(group) == 6
        assert group.holm_two_sided_p.notna().all()
        assert (group.holm_two_sided_p >= group.two_sided_p).all()


def test_bin_composition_keeps_empty_bins_and_largest_dataset() -> None:
    composition = bin_composition(synthetic_predictions())
    assert len(composition) == 30
    assert composition.groupby(["endpoint", "display_class"]).size().eq(5).all()

    highest = composition.loc[
        composition.endpoint.eq("Ki")
        & composition.display_class.eq("medium")
        & composition.probability_bin.eq(4)
    ].iloc[0]
    assert highest.n_predictions == 10
    assert highest.n_datasets == 2
    assert highest.largest_dataset == "T1"
    assert highest.largest_dataset_fraction == 0.5
    assert highest.agreement_rate == 1.0

    empty = composition.loc[
        composition.endpoint.eq("Ki")
        & composition.display_class.eq("medium")
        & composition.probability_bin.eq(0)
    ].iloc[0]
    assert empty.n_predictions == 0
    assert empty.n_datasets == 0
    assert empty.largest_dataset == ""
    assert pd.isna(empty.largest_dataset_fraction)
    assert pd.isna(empty.agreement_rate)


def test_archived_e22_expected_values() -> None:
    required = (
        E5C / "predictions_e5c.csv",
        OUT / "class_reliability_by_design.csv",
        OUT / "class_reliability_tests.csv",
        OUT / "class_bias_dataset_counts.csv",
        OUT / "e11_class_calibration.csv",
        OUT / "e11_v2_default_display.csv",
        OUT / "e11_document_overlap_analyzed.csv",
        OUT / "fig4a_bin_composition.csv",
    )
    if not all(path.exists() for path in required):
        pytest.skip("E22 archived outputs are not available")
    validate_expected(
        pd.read_csv(OUT / "class_reliability_by_design.csv"),
        pd.read_csv(OUT / "class_reliability_tests.csv"),
        pd.read_csv(OUT / "class_bias_dataset_counts.csv"),
        pd.read_csv(OUT / "e11_class_calibration.csv"),
        pd.read_csv(OUT / "e11_v2_default_display.csv"),
        pd.read_csv(OUT / "e11_document_overlap_analyzed.csv"),
    )
    validate_bin_composition(pd.read_csv(OUT / "fig4a_bin_composition.csv"))
