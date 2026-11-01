"""Build final data tables and source inventories for Figures 3-5."""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


SOURCE_DIR = Path(__file__).resolve().parent
FIG_DIR = SOURCE_DIR.parent
ROOT = Path(__file__).resolve().parents[3]
S_OUTPUT = ROOT / "paper/strengthening/output"
PRISM_DATA = S_OUTPUT / "F0"
OUTPUT = FIG_DIR / "data_final"
TARGET_SYMBOLS = FIG_DIR / "data" / "target_symbols.csv"
DISPLAY_CLASSES = ("inactive", "medium", "active")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, float_precision="round_trip")


def write_csv(frame: pd.DataFrame, name: str) -> Path:
    path = OUTPUT / name
    frame.to_csv(path, index=False, lineterminator="\n")
    return path


def copy_unchanged(name: str, records: list[dict[str, object]]) -> Path:
    source = PRISM_DATA / name
    target = OUTPUT / name
    shutil.copyfile(source, target)
    source_hash = sha256(source)
    target_hash = sha256(target)
    if source_hash != target_hash:
        raise AssertionError(f"{name} changed hash after copying")
    frame = read_csv(target)
    records.append(
        {
            "file": name,
            "sha256": target_hash,
            "columns": list(frame.columns),
            "source_files": [str(source.relative_to(ROOT))],
            "source_hashes": [source_hash],
            "derivation": "byte-for-byte copy from accepted Prism data",
            "unchanged_copy": True,
        }
    )
    return target


def attach_symbols(frame: pd.DataFrame) -> pd.DataFrame:
    symbols = read_csv(TARGET_SYMBOLS)[["target_id", "symbol", "development_target"]]
    merged = frame.merge(symbols, on="target_id", how="left", validate="many_to_one")
    if merged["symbol"].isna().any():
        missing = merged.loc[merged["symbol"].isna(), "target_id"].unique().tolist()
        raise AssertionError(f"Missing target gene symbols: {missing}")
    return merged


def wilson_interval(successes: int, n: int) -> tuple[float, float]:
    if n == 0:
        return np.nan, np.nan
    z = float(stats.norm.ppf(0.975))
    rate = successes / n
    denominator = 1.0 + z**2 / n
    center = (rate + z**2 / (2 * n)) / denominator
    half = z * np.sqrt(rate * (1 - rate) / n + z**2 / (4 * n**2)) / denominator
    return float(center - half), float(center + half)


def calibration_curve(predictions: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for endpoint in ("IC50", "Ki"):
        endpoint_frame = predictions.loc[predictions["endpoint"].eq(endpoint)].copy()
        for population, mask in (
            ("untouched", ~endpoint_frame["flag_touched"].astype(bool)),
            ("touched", endpoint_frame["flag_touched"].astype(bool)),
        ):
            selected = endpoint_frame.loc[mask].copy()
            selected["probability_bin"] = np.clip(
                np.floor(selected["p_displayed"] * 5 + 1e-12).astype(int),
                0,
                4,
            )
            for bin_index in range(5):
                group = selected.loc[selected["probability_bin"].eq(bin_index)]
                n = len(group)
                low, high = wilson_interval(int(group["agree"].sum()), n)
                rows.append(
                    {
                        "analysis_role": "post hoc",
                        "design": "entry",
                        "endpoint": endpoint,
                        "population": population,
                        "display_class": "all",
                        "probability_bin": bin_index,
                        "n": n,
                        "mean_p": float(group["p_displayed"].mean()) if n else np.nan,
                        "agreement_rate": float(group["agree"].mean()) if n else np.nan,
                        "wilson_ci_low": low,
                        "wilson_ci_high": high,
                        "omit_in_figure": bool(n < 8),
                    }
                )
    return pd.DataFrame(rows)


def calibration_by_display_class(predictions: pd.DataFrame) -> pd.DataFrame:
    """Bin class reproducibility separately for each displayed class."""
    rows = []
    for endpoint in ("IC50", "Ki"):
        endpoint_frame = predictions.loc[predictions["endpoint"].eq(endpoint)].copy()
        endpoint_frame["probability_bin"] = np.clip(
            np.floor(endpoint_frame["p_displayed"] * 5 + 1e-12).astype(int),
            0,
            4,
        )
        for display_index, display_class in enumerate(DISPLAY_CLASSES):
            selected = endpoint_frame.loc[
                endpoint_frame["displayed"].eq(display_index)
            ]
            for bin_index in range(5):
                group = selected.loc[selected["probability_bin"].eq(bin_index)]
                n = len(group)
                low, high = wilson_interval(int(group["agree"].sum()), n)
                rows.append(
                    {
                        "analysis_role": "post hoc",
                        "design": "entry",
                        "endpoint": endpoint,
                        "display_class": display_class,
                        "probability_bin": bin_index,
                        "n": n,
                        "mean_p": float(group["p_displayed"].mean()) if n else np.nan,
                        "agreement_rate": float(group["agree"].mean()) if n else np.nan,
                        "wilson_ci_low": low,
                        "wilson_ci_high": high,
                        "omit_in_figure": bool(n < 8),
                    }
                )
    return pd.DataFrame(rows)


def build_3b(
    endpoint: str,
    forest: pd.DataFrame,
    fair: pd.DataFrame,
    records: list[dict[str, object]],
) -> Path:
    g = forest.loc[
        forest["endpoint"].eq(endpoint) & forest["population"].eq("all"),
        [
            "target_id",
            "g_auc",
            "g_auc_ci_low",
            "g_auc_ci_high",
            "n_predictions",
        ],
    ].rename(
        columns={
            "g_auc_ci_low": "g_ci_low",
            "g_auc_ci_high": "g_ci_high",
            "n_predictions": "g_n_predictions",
        }
    )
    d = fair.loc[
        fair["endpoint"].eq(endpoint)
        & fair["score"].eq("D")
        & fair["subset"].eq("E4 stored E2 all"),
        ["target", "auc", "n_predictions"],
    ].rename(
        columns={
            "target": "target_id",
            "auc": "d_auc",
            "n_predictions": "d_n_predictions",
        }
    )
    b2 = fair.loc[
        fair["endpoint"].eq(endpoint)
        & fair["score"].eq("B-range-2")
        & fair["subset"].eq("remaining_records>=2"),
        ["target", "auc", "n_predictions"],
    ].rename(
        columns={
            "target": "target_id",
            "auc": "b_range_2_auc",
            "n_predictions": "b_range_2_n_predictions",
        }
    )
    table = g.merge(d, on="target_id", how="left", validate="one_to_one")
    table = table.merge(b2, on="target_id", how="left", validate="one_to_one")
    table = attach_symbols(table)
    table = table.sort_values("g_auc", ascending=False, kind="stable").reset_index(drop=True)
    table.insert(0, "target_rank", np.arange(1, len(table) + 1))
    expected = 29 if endpoint == "IC50" else 23
    if len(table) != expected:
        raise AssertionError(f"3b {endpoint} has rows {len(table)} instead of {expected}")
    if table["g_auc"].isna().any() or table["d_auc"].isna().any():
        raise AssertionError(f"3b {endpoint} has a missing G or D AUC")
    if endpoint == "Ki":
        missing_b2 = table.loc[table["b_range_2_auc"].isna(), "target_id"].tolist()
        if missing_b2 != ["CHEMBL4005"]:
            raise AssertionError(f"Unexpected missing B-range-2 rows in 3b Ki: {missing_b2}")
    columns = [
        "target_rank",
        "target_id",
        "symbol",
        "development_target",
        "g_auc",
        "g_ci_low",
        "g_ci_high",
        "d_auc",
        "d_n_predictions",
        "b_range_2_auc",
        "b_range_2_n_predictions",
        "g_n_predictions",
    ]
    table = table[columns]
    name = f"3b_{endpoint}.csv"
    path = write_csv(table, name)
    sources = [
        S_OUTPUT / "E2" / "forest_plot_data.csv",
        S_OUTPUT / "E5d" / "comparators_fair.csv",
    ]
    records.append(
        {
            "file": name,
            "sha256": sha256(path),
            "columns": columns,
            "source_files": [str(path.relative_to(ROOT)) for path in sources],
            "source_hashes": [sha256(path) for path in sources],
            "derivation": (
                "G and 95% CI from E2 forest_plot_data.csv; D and B-range-2 "
                "from E5d comparators_fair.csv; rows sorted by G AUC"
            ),
            "unchanged_copy": False,
            "assertions": {
                "rows": len(table),
                "target_rank": list(range(1, len(table) + 1)),
                "g_sorted_descending": bool(table["g_auc"].is_monotonic_decreasing),
            },
        }
    )
    return path


def build_3c(
    endpoint: str,
    forest: pd.DataFrame,
    hypotheses: pd.DataFrame,
    old_ranks: pd.DataFrame,
    records: list[dict[str, object]],
) -> Path:
    old_ranks = old_ranks.loc[old_ranks["row_type"].eq("target")].copy()
    source_columns = {
        "g_minus_g_record_log_loss": "g_minus_record_log_loss",
        "g_minus_g_global_log_loss": "g_minus_global_log_loss",
        "g_minus_g_fixed_log_loss": "g_minus_fixed_log_loss",
        "g_minus_x_log_loss": "g_minus_x_log_loss",
    }
    values = forest.loc[
        forest["endpoint"].eq(endpoint) & forest["population"].eq("all"),
        ["target_id", *source_columns],
    ].rename(columns=source_columns)
    values = old_ranks[["target_id", "target_rank", "x_minus_g"]].merge(
        values,
        on="target_id",
        how="left",
        validate="one_to_one",
    )
    values.loc[values["x_minus_g"].isna(), "g_minus_x_log_loss"] = np.nan
    values = values.drop(columns="x_minus_g")
    old_columns = {
        "g_minus_record_log_loss": "g_record_minus_g",
        "g_minus_global_log_loss": "g_global_minus_g",
        "g_minus_fixed_log_loss": "g_fixed_minus_g",
        "g_minus_x_log_loss": "x_minus_g",
    }
    for new_column, old_column in old_columns.items():
        mask = values[new_column].notna() & old_ranks[old_column].notna()
        if not np.allclose(
            values.loc[mask, new_column],
            -old_ranks.loc[mask, old_column],
            rtol=0.0,
            atol=1e-15,
        ):
            raise AssertionError(f"3c {endpoint} {new_column} is not the negated original value")
    if values[list(source_columns.values())[:3]].isna().any().any():
        raise AssertionError(f"3c {endpoint} has targets missing a G-minus-comparator difference")
    values.insert(0, "row_type", "target")
    values = values.sort_values("target_rank").reset_index(drop=True)

    summary: dict[str, object] = {
        "row_type": "summary",
        "target_id": "median",
        "target_rank": np.nan,
    }
    hypothesis_map = {
        "H4": "g_minus_record_log_loss",
        "H5": "g_minus_global_log_loss",
        "H6": "g_minus_fixed_log_loss",
        "H7": "g_minus_x_log_loss",
    }
    for hypothesis_id, column in hypothesis_map.items():
        row = hypotheses.loc[
            hypotheses["population"].eq("all")
            & hypotheses["endpoint"].eq(endpoint)
            & hypotheses["hypothesis"].eq(hypothesis_id)
        ].iloc[0]
        summary[f"median_{column}"] = float(row["median_difference"])
        summary[f"ci_low_{column}"] = float(row["median_ci_low"])
        summary[f"ci_high_{column}"] = float(row["median_ci_high"])
    table = pd.concat([values, pd.DataFrame([summary])], ignore_index=True)
    columns = [
        "row_type",
        "target_id",
        "target_rank",
        "g_minus_record_log_loss",
        "g_minus_global_log_loss",
        "g_minus_fixed_log_loss",
        "g_minus_x_log_loss",
    ]
    for stem in (
        "g_minus_record_log_loss",
        "g_minus_global_log_loss",
        "g_minus_fixed_log_loss",
        "g_minus_x_log_loss",
    ):
        columns.extend([f"median_{stem}", f"ci_low_{stem}", f"ci_high_{stem}"])
    table = table.reindex(columns=columns)
    name = f"3c_{endpoint}.csv"
    path = write_csv(table, name)
    sources = [
        S_OUTPUT / "E2" / "forest_plot_data.csv",
        S_OUTPUT / "E2" / "hypothesis_tests.csv",
        PRISM_DATA / f"3c_{endpoint}.csv",
    ]
    records.append(
        {
            "file": name,
            "sha256": sha256(path),
            "columns": columns,
            "source_files": [str(path.relative_to(ROOT)) for path in sources],
            "source_hashes": [sha256(path) for path in sources],
            "derivation": (
                "G − comparator values from E2 forest_plot_data.csv; medians and "
                "95% CIs from E2 hypothesis_tests.csv; accepted target order retained"
            ),
            "unchanged_copy": False,
            "assertions": {
                "target_rows": int((table["row_type"] == "target").sum()),
                "X_rows": int(table["g_minus_x_log_loss"].notna().sum()),
            },
        }
    )
    return path


def build_4d(
    endpoint: str,
    source: pd.DataFrame,
    records: list[dict[str, object]],
) -> Path:
    selected = source.loc[
        source["design"].eq("time_2015")
        & source["endpoint"].eq(endpoint)
        & source["earlier_units_group"].isin(["m_eq_1", "m_ge_2"])
        & source["state"].isin(["P_supported", "uncertain_low_support", "touched_all"])
    ].copy()
    if len(selected) != 6:
        raise AssertionError(f"4d {endpoint} has rows {len(selected)} instead of 6")
    selected = selected.sort_values(
        ["earlier_units_group", "state"],
        key=lambda values: values.map(
            {"m_eq_1": 0, "m_ge_2": 1, "P_supported": 0,
             "uncertain_low_support": 1, "touched_all": 2}
        ),
        kind="stable",
    ).reset_index(drop=True)
    name = f"4d_{endpoint}.csv"
    path = write_csv(selected, name)
    source_path = S_OUTPUT / "E12a" / "fig4d_three_groups.csv"
    records.append(
        {
            "file": name,
            "sha256": sha256(path),
            "columns": list(selected.columns),
            "source_files": [str(source_path.relative_to(ROOT))],
            "source_hashes": [sha256(source_path)],
            "derivation": (
                "E12a 2015 time split; groups m=1 and m>=2; "
                "states P_supported, uncertain_low_support and touched_all"
            ),
            "unchanged_copy": False,
            "assertions": {
                "rows": len(selected),
                "omitted_bars_n_lt_8": selected.loc[
                    selected["n_predictions"] < 8,
                    ["earlier_units_group", "state", "n_predictions"],
                ].to_dict("records"),
            },
        }
    )
    return path


def main() -> int:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, object]] = []
    symbol_frame = read_csv(TARGET_SYMBOLS)
    if symbol_frame["target_id"].duplicated().any():
        raise AssertionError("target_symbols.csv contains duplicate target IDs")
    target_symbols_output = OUTPUT / TARGET_SYMBOLS.name
    shutil.copyfile(TARGET_SYMBOLS, target_symbols_output)
    if sha256(TARGET_SYMBOLS) != sha256(target_symbols_output):
        raise AssertionError("target_symbols.csv changed hash after copying")
    records.append(
        {
            "file": target_symbols_output.name,
            "sha256": sha256(target_symbols_output),
            "columns": list(symbol_frame.columns),
            "source_files": [str(TARGET_SYMBOLS.relative_to(ROOT))],
            "source_hashes": [sha256(TARGET_SYMBOLS)],
            "derivation": "target symbols supplied for the 52 evaluation datasets",
            "unchanged_copy": True,
        }
    )

    for name in (
        "3a_B0_disagreement.csv",
        "4b_hERG_patent_pair.csv",
        "4c_IC50.csv",
        "4c_Ki.csv",
        "5a_recall.csv",
        "5b_layout_contrasts.csv",
        "5d_layout_correlation.csv",
        "5e_core_coverage_by_group_size.csv",
    ):
        copy_unchanged(name, records)

    forest = read_csv(S_OUTPUT / "E2" / "forest_plot_data.csv")
    hypotheses = read_csv(S_OUTPUT / "E2" / "hypothesis_tests.csv")
    fair = read_csv(S_OUTPUT / "E5d" / "comparators_fair.csv")
    for endpoint in ("IC50", "Ki"):
        build_3b(endpoint, forest, fair, records)
        old_ranks = read_csv(PRISM_DATA / f"3c_{endpoint}.csv")
        build_3c(endpoint, forest, hypotheses, old_ranks, records)

    predictions_path = S_OUTPUT / "E5c" / "predictions_e5c.csv"
    predictions = pd.read_csv(
        predictions_path,
        usecols=[
            "split_design",
            "merge_basis",
            "target",
            "endpoint",
            "displayed",
            "p_displayed",
            "agree",
            "method",
            "flag_touched",
        ],
        float_precision="round_trip",
    )
    entry_predictions = predictions.loc[
        predictions["split_design"].eq("entry")
        & predictions["merge_basis"].eq("same_value")
        & predictions["method"].eq("G")
    ].reset_index(drop=True)
    legacy_calibration = OUTPUT / "4a_calibration_curves.csv"
    if legacy_calibration.exists():
        legacy_calibration.unlink()
    calibration = calibration_by_display_class(entry_predictions)
    calibration_name = "4a_calibration_by_class.csv"
    calibration_output = write_csv(calibration, calibration_name)
    records.append(
        {
            "file": calibration_name,
            "sha256": sha256(calibration_output),
            "columns": list(calibration.columns),
            "source_files": [str(predictions_path.relative_to(ROOT))],
            "source_hashes": [sha256(predictions_path)],
            "derivation": (
                "E5c predictions; entry design, model G, five equal-width bins, "
                "separately for each displayed class; Wilson 95% intervals; "
                "bins with fewer than 8 predictions omitted"
            ),
            "unchanged_copy": False,
            "assertions": {
                "rows": len(calibration),
                "omitted_rows": int(calibration["omit_in_figure"].sum()),
                "display_classes": list(calibration["display_class"].unique()),
                "recomputed_from_e5c": True,
            },
        }
    )

    without_orexin = entry_predictions.loc[
        ~entry_predictions["target"].eq("CHEMBL5113")
    ].reset_index(drop=True)
    si_4a = calibration_curve(without_orexin)
    si_name = "si_4a_without_CHEMBL5113.csv"
    si_output = write_csv(si_4a, si_name)
    records.append(
        {
            "file": si_name,
            "sha256": sha256(si_output),
            "columns": list(si_4a.columns),
            "source_files": [str(predictions_path.relative_to(ROOT))],
            "source_hashes": [sha256(predictions_path)],
            "derivation": (
                "same 4a calibration function as the main table, "
                "after removing target CHEMBL5113"
            ),
            "unchanged_copy": False,
            "assertions": {
                "rows": len(si_4a),
                "omitted_rows": int(si_4a["omit_in_figure"].sum()),
            },
        }
    )

    fig4d_path = S_OUTPUT / "E12a" / "fig4d_three_groups.csv"
    fig4d = read_csv(fig4d_path)
    for endpoint in ("IC50", "Ki"):
        build_4d(endpoint, fig4d, records)

    fixed_detail_path = S_OUTPUT / "E6" / "recall_fixedk_detailed.csv"
    fixed_detail = read_csv(fixed_detail_path)
    fixed = (
        fixed_detail.groupby(["target_label", "layout"], as_index=False, sort=False)
        .agg(
            retained_recall_mean=("retained_recall", "mean"),
            fixed_recall_mean=("fixed_recall", "mean"),
            n_seeds=("seed", "nunique"),
        )
    )
    baseline = read_csv(PRISM_DATA / "5a_recall.csv")
    layouts = ["projection", "minimum_displacement", "neighbor_refined"]
    for row in baseline.itertuples(index=False):
        for layout in layouts:
            values = np.array(
                [getattr(row, f"{layout}__r{seed:02d}") for seed in range(1, 11)],
                dtype=float,
            )
            expected = float(values.mean())
            actual = float(
                fixed.loc[
                    fixed["target_label"].eq(row.Row_title)
                    & fixed["layout"].eq(layout),
                    "retained_recall_mean",
                ].iloc[0]
            )
            if abs(actual - expected) > 1e-9:
                raise AssertionError(
                    f"5a {row.Row_title} {layout} retained mean mismatch: "
                    f"{actual} vs {expected}"
                )
    if not (fixed["n_seeds"] == 10).all():
        raise AssertionError("5a fixed-k summary does not contain ten seeds per combination")
    fixed_name = "5a_fixedk_means.csv"
    fixed_output = write_csv(fixed, fixed_name)
    records.append(
        {
            "file": fixed_name,
            "sha256": sha256(fixed_output),
            "columns": list(fixed.columns),
            "source_files": [str(fixed_detail_path.relative_to(ROOT))],
            "source_hashes": [sha256(fixed_detail_path)],
            "derivation": "mean of ten seed-level fixed_recall values in E6 recall_fixedk_detailed.csv",
            "unchanged_copy": False,
            "assertions": {
                "rows": len(fixed),
                "seeds_per_combination": sorted(fixed["n_seeds"].unique().tolist()),
                "retained_means_equal_5a_within_1e-9": True,
            },
        }
    )

    manifest = {
        "schema_version": 1,
        "generated_by": str(Path(__file__).resolve().relative_to(ROOT)),
        "output_directory": str(OUTPUT.relative_to(ROOT)),
        "record_count": len(records),
        "files": records,
    }
    (OUTPUT / "data_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        f"data_final files={len(records)} "
        f"unchanged_copies={sum(bool(item['unchanged_copy']) for item in records)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
