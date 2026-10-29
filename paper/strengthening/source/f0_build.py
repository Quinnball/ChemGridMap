"""Build Prism-ready tables for manuscript Figures 3-5.

This module only reads frozen analysis outputs. It does not import or execute
ChemGridMap and it does not access Prism.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import pickle
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd


SOURCE = Path(__file__).resolve().parent
STRENGTHENING = SOURCE.parent
REV = STRENGTHENING.parent
WORKSPACE = REV.parents[1]
CGM = WORKSPACE / "chemgridmap"

S_OUT = STRENGTHENING / "output"
E2_OUT = S_OUT / "E2"
E1A_OUT = S_OUT / "E1A"
E1B_OUT = S_OUT / "E1B"
E2C_OUT = S_OUT / "E2c"
E4_OUT = S_OUT / "E4"
CGM_OUT = CGM / "paper" / "output" / "reliability_v1"

DATA_OUT = REV / "prism_jcheminf" / "data"
PLAN_OUT = REV / "prism_jcheminf" / "source" / "gui_plan.md"
MANIFEST_OUT = DATA_OUT / "data_manifest.json"

ENDPOINTS = ("IC50", "Ki")
TARGET_COLORS = {"CA II": "#1D5689", "Thrombin": "#006354", "hERG": "#624839"}
ENDPOINT_COLORS = {"IC50": "#D55E00", "Ki": "#56B4E9"}
DISPLAY_STATUS_COLORS = {"supported": "#1D5689", "uncertain": "#995500"}
COMPARISON_COLORS = {
    "minimum_displacement": "#995500",
    "las_random": "#172536",
    "las_projection": "#624839",
    "flas_random": "#CC79A7",
}
LAYOUT_COLORS = {
    "projection": "#1D5689",
    "minimum_displacement": "#995500",
    "neighbor_refined": "#006354",
    "las_projection": "#624839",
    "las_random": "#172536",
    "flas_random": "#CC79A7",
    "random": "#666666",
}

ATR_DOCUMENTS = {
    "CHEMBL5727267": ("US-10392391-B2", 2019),
    "CHEMBL5727989": ("US-10787452-B2", 2020),
    "CHEMBL5728744": ("US-11117900-B2", 2021),
    "CHEMBL5729180": ("US-11370798-B2", 2022),
    "CHEMBL5727412": ("US-10479784-B2", 2019),
    "CHEMBL5728359": ("US-10961232-B2", 2021),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rel(path: Path) -> str:
    return path.relative_to(WORKSPACE).as_posix()


def source_paths() -> dict[str, Path]:
    return {
        "e2_forest": E2_OUT / "forest_plot_data.csv",
        "e2_metrics_by_target": E2_OUT / "metrics_by_target.csv",
        "e2_flagged_pairs_census": E2_OUT / "flagged_pairs_census.csv",
        "e2_sigma_doc": E2_OUT / "sigma_doc_by_target.csv",
        "e2_curated_cache": STRENGTHENING / "qa" / "E2" / "curated_records.pkl.gz",
        "e2_target_inclusion": E2_OUT / "target_inclusion.csv",
        "e1a_flagged_pairs": E1A_OUT / "flagged_pairs_profile.csv",
        "e1b_core_by_group_size": E1B_OUT / "core_by_group_size.csv",
        "e1b_neighbor_sar_summary": E1B_OUT / "neighbor_sar_summary.csv",
        "e2c_layout_by_target": E2C_OUT / "layout_by_target.csv",
        "e4_comparators": E4_OUT / "comparators_posthoc.csv",
        "e4_calibration": E4_OUT / "calibration_curve_displayed.csv",
        "e4_time_split": E4_OUT / "time_split_posthoc.csv",
        "e4_layout_baselines": E4_OUT / "layout_baselines_ten_seed.csv",
        "cgm_figure_3a": CGM_OUT / "figure_data" / "Figure_3A_neighbor_recall.csv",
        "cgm_figure_3c": CGM_OUT / "figure_data" / "Figure_3C_CHEMBL81688_neighbors.csv",
        "cgm_hERG_retained": (
            CGM_OUT / "software_checks_050" / "chembl240"
            / "chembl240_retained_records.csv"
        ),
    }


def load_sources(paths: dict[str, Path]) -> dict[str, pd.DataFrame | dict]:
    data: dict[str, pd.DataFrame | dict] = {}
    for key, path in paths.items():
        if path.suffix == ".gz":
            data[key] = pickle.loads(gzip.decompress(path.read_bytes()))
        elif path.suffix == ".csv":
            data[key] = pd.read_csv(path)
        else:
            raise ValueError(f"Unsupported source type: {path}")
    return data


def write_csv(frame: pd.DataFrame, name: str) -> Path:
    path = DATA_OUT / name
    frame.to_csv(path, index=False, float_format="%.15g", lineterminator="\n")
    return path


def target_rank(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.sort_values(["g_auc", "target_id"], ascending=[False, True]).copy()
    out.insert(0, "target_rank", np.arange(1, len(out) + 1))
    return out


def build_3a(data: dict) -> pd.DataFrame:
    metrics = data["e2_metrics_by_target"]
    b0 = metrics.loc[
        metrics.population.eq("all")
        & metrics.model.eq("B0")
        & metrics.stratum.eq("all"),
        ["target", "endpoint", "agreement_rate"],
    ].copy()
    b0["disagreement_rate"] = 1.0 - b0.agreement_rate
    out = (
        b0.pivot(index="target", columns="endpoint", values="disagreement_rate")
        .reindex(columns=list(ENDPOINTS))
        .rename_axis("target_id")
        .reset_index()
    )
    return out


def build_3b(data: dict) -> dict[str, pd.DataFrame]:
    forest = data["e2_forest"]
    comparators = data["e4_comparators"]
    out = {}
    for endpoint in ENDPOINTS:
        selected = forest.loc[
            forest.population.eq("all") & forest.endpoint.eq(endpoint),
            [
                "target_id",
                "n_predictions",
                "g_auc",
                "g_auc_ci_low",
                "g_auc_ci_high",
            ],
        ].rename(columns={"n_predictions": "g_n_predictions"}).copy()
        posthoc = comparators.loc[
            comparators.endpoint.eq(endpoint)
            & comparators.population.eq("all")
            & comparators.comparator.isin(["G", "D", "B-sd"]),
            [
                "target",
                "comparator",
                "n_predictions",
                "auc",
                "auc_ci_low",
                "auc_ci_high",
            ],
        ].rename(columns={"target": "target_id"})
        checks = posthoc.pivot(index="target_id", columns="comparator", values="auc")
        g_check = selected.set_index("target_id").g_auc
        if not np.allclose(
            checks["G"].reindex(g_check.index),
            g_check,
            rtol=0,
            atol=1e-15,
        ):
            raise AssertionError("forest and E4 G AUC disagree")
        d = posthoc.loc[posthoc.comparator.eq("D")].drop(columns="comparator")
        d = d.rename(
            columns={
                "n_predictions": "d_n_predictions",
                "auc": "d_auc",
                "auc_ci_low": "d_auc_ci_low",
                "auc_ci_high": "d_auc_ci_high",
            }
        )
        b_range = posthoc.loc[posthoc.comparator.eq("B-sd")].drop(
            columns="comparator"
        )
        b_range = b_range.rename(
            columns={
                "n_predictions": "b_range_n_predictions",
                "auc": "b_range_auc",
                "auc_ci_low": "b_range_auc_ci_low",
                "auc_ci_high": "b_range_auc_ci_high",
            }
        )
        selected = selected.merge(d, on="target_id", validate="one_to_one")
        selected = selected.merge(b_range, on="target_id", validate="one_to_one")
        selected = target_rank(selected)
        out[endpoint] = selected[
            [
                "target_rank",
                "target_id",
                "g_auc",
                "g_auc_ci_low",
                "g_auc_ci_high",
                "d_auc",
                "d_auc_ci_low",
                "d_auc_ci_high",
                "b_range_auc",
                "b_range_auc_ci_low",
                "b_range_auc_ci_high",
                "g_n_predictions",
                "d_n_predictions",
                "b_range_n_predictions",
            ]
        ].rename(columns={"g_auc_ci_low": "g_ci_low", "g_auc_ci_high": "g_ci_high"})
    return out


def bootstrap_median_ci(values: np.ndarray, seed: int = 20261028, repetitions: int = 2000):
    rng = np.random.default_rng(seed)
    draws = np.empty(repetitions)
    for index in range(repetitions):
        draws[index] = np.median(rng.choice(values, size=len(values), replace=True))
    low, high = np.percentile(draws, [2.5, 97.5])
    return float(np.median(values)), float(low), float(high)


def build_3c(data: dict, b_tables: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    metrics = data["e2_metrics_by_target"]
    alternatives = {
        "g_record_minus_g": "G-record",
        "g_global_minus_g": "G-global",
        "g_fixed_minus_g": "G-fixed",
        "x_minus_g": "X",
    }
    flagged = data["e2_flagged_pairs_census"]
    flagged_targets = set(
        flagged.loc[
            flagged.row_type.eq("target_summary") & flagged.n_flagged_pairs.gt(0),
            ["target", "endpoint"],
        ].itertuples(index=False, name=None)
    )
    out = {}
    for endpoint in ENDPOINTS:
        source = metrics.loc[
            metrics.population.eq("all")
            & metrics.endpoint.eq(endpoint)
            & metrics.stratum.eq("all")
            & metrics.model.isin(["G", *alternatives.values()]),
            ["target", "model", "log_loss"],
        ]
        losses = source.pivot(index="target", columns="model", values="log_loss")
        values = pd.DataFrame({"target_id": losses.index})
        for name, model in alternatives.items():
            values[name] = losses[model].to_numpy() - losses["G"].to_numpy()
        values.loc[
            ~values.target_id.isin(
                {target for target, item_endpoint in flagged_targets if item_endpoint == endpoint}
            ),
            "x_minus_g",
        ] = np.nan
        ranks = b_tables[endpoint][["target_id", "target_rank"]]
        values = values.merge(ranks, on="target_id", validate="one_to_one")
        values = values.sort_values("target_rank").copy()
        values.insert(0, "row_type", "target")

        summary = {"row_type": "summary", "target_id": "median", "target_rank": np.nan}
        for name in alternatives:
            median, low, high = bootstrap_median_ci(
                values[name].dropna().to_numpy(float)
            )
            summary[f"median_{name}"] = median
            summary[f"ci_low_{name}"] = low
            summary[f"ci_high_{name}"] = high
        combined = pd.concat([values, pd.DataFrame([summary])], ignore_index=True)
        out[endpoint] = combined
    return out


def build_4a(data: dict) -> dict[str, pd.DataFrame]:
    calibration = data["e4_calibration"]
    out = {}
    for endpoint in ENDPOINTS:
        selected = calibration.loc[
            calibration.endpoint.eq(endpoint) & calibration.eligible.astype(bool)
        ].copy()
        g = selected.loc[selected.model.eq("G")].sort_values("bin").reset_index(drop=True)
        b0 = selected.loc[selected.model.eq("B0")].reset_index(drop=True)
        if len(b0) != 1 or not np.isclose(float(b0.iloc[0].mean_probability), 1.0):
            raise AssertionError("B0 must contribute one calibration point at probability 1")
        rows = max(len(g), len(b0))
        table = pd.DataFrame(
            {
                "bin": g.bin.to_numpy(),
                "G_probability": g.mean_probability.to_numpy(float),
                "G_observed_agreement": g.observed_agreement.to_numpy(float),
                "G_n_predictions": g.n_predictions.to_numpy(int),
                "B0_probability": np.full(rows, np.nan),
                "B0_observed_agreement": np.full(rows, np.nan),
                "B0_n_predictions": np.full(rows, np.nan),
            }
        )
        table.loc[0, "B0_probability"] = float(b0.iloc[0].mean_probability)
        table.loc[0, "B0_observed_agreement"] = float(b0.iloc[0].observed_agreement)
        table.loc[0, "B0_n_predictions"] = int(b0.iloc[0].n_predictions)
        out[endpoint] = table
    return out


def build_4b(data: dict) -> pd.DataFrame:
    records = data["cgm_hERG_retained"]
    flagged = data["e1a_flagged_pairs"]
    docs = ("CHEMBL5725903", "CHEMBL5727624")
    selected = records.loc[records.document_chembl_id.isin(docs)].copy()
    units = (
        selected.groupby(
            ["molecule_identity_key", "document_chembl_id"], as_index=False
        )
        .pchembl_value.median()
    )
    first = units.loc[units.document_chembl_id.eq(docs[0])].set_index("molecule_identity_key")
    second = units.loc[units.document_chembl_id.eq(docs[1])].set_index("molecule_identity_key")
    shared = sorted(first.index.intersection(second.index))
    out = pd.DataFrame(
        {
            "compound_id": shared,
            "us_10131663_b2": first.loc[shared, "pchembl_value"].to_numpy(float),
            "us_10590130_b2": second.loc[shared, "pchembl_value"].to_numpy(float),
        }
    )
    out.insert(0, "compound_rank", np.arange(1, len(out) + 1))
    out["offset_b_minus_a"] = out.us_10590130_b2 - out.us_10131663_b2

    profile = flagged.loc[
        flagged.row_type.eq("shared_entry")
        & flagged.document_a.eq(docs[0])
        & flagged.document_b.eq(docs[1])
    ]
    if len(out) != 30 or len(profile) != 30:
        raise AssertionError("The frozen hERG pair must contain 30 shared entries")
    if not np.isclose(out.offset_b_minus_a.median(), 3.47, atol=0.005):
        raise AssertionError("The hERG pair median offset must be 3.47")
    return out


def wilson_interval(rate: float, n: float) -> tuple[float, float]:
    """95% Wilson interval for a binomial rate."""
    z = 1.959963984540054
    n_float = float(n)
    denominator = 1.0 + z * z / n_float
    centre = (rate + z * z / (2.0 * n_float)) / denominator
    half_width = (
        z
        * np.sqrt(rate * (1.0 - rate) / n_float + z * z / (4.0 * n_float**2))
        / denominator
    )
    return float(max(0.0, centre - half_width)), float(min(1.0, centre + half_width))


def build_4c(data: dict) -> dict[str, pd.DataFrame]:
    split = data["e4_time_split"]
    selected = split.loc[
        split.record_type.eq("dataset_comparison")
        & split.variant.eq("latest_unit_year")
        & split.cutoff_year.eq(2015)
        & split.endpoint.isin(ENDPOINTS)
    ].copy()
    out = {}
    for endpoint in ENDPOINTS:
        frame = selected.loc[selected.endpoint.eq(endpoint)].copy()
        frame["supported_interval"] = [
            wilson_interval(rate, n)
            for rate, n in zip(frame.supported_agreement_rate, frame.supported_n)
        ]
        frame["uncertain_interval"] = [
            wilson_interval(rate, n)
            for rate, n in zip(frame.uncertain_agreement_rate, frame.uncertain_n)
        ]
        frame["supported_ci_low"] = [item[0] for item in frame.supported_interval]
        frame["supported_ci_high"] = [item[1] for item in frame.supported_interval]
        frame["uncertain_ci_low"] = [item[0] for item in frame.uncertain_interval]
        frame["uncertain_ci_high"] = [item[1] for item in frame.uncertain_interval]
        frame["difference_supported_minus_uncertain"] = (
            frame.supported_agreement_rate - frame.uncertain_agreement_rate
        )
        frame = frame.sort_values(
            ["supported_agreement_rate", "uncertain_agreement_rate", "target"],
            ascending=[False, False, True],
        ).reset_index(drop=True)
        frame.insert(0, "dataset_rank", np.arange(1, len(frame) + 1))
        out[endpoint] = frame[
            [
                "dataset_rank",
                "target",
                "supported_agreement_rate",
                "supported_ci_low",
                "supported_ci_high",
                "uncertain_agreement_rate",
                "uncertain_ci_low",
                "uncertain_ci_high",
                "difference_supported_minus_uncertain",
                "supported_n",
                "uncertain_n",
            ]
        ].rename(
            columns={
                "target": "target_id",
                "supported_agreement_rate": "supported_agreement",
                "uncertain_agreement_rate": "uncertain_agreement",
                "supported_n": "supported_n_predictions",
                "uncertain_n": "uncertain_n_predictions",
            }
        )
    return out


def build_4d(data: dict) -> dict[str, pd.DataFrame]:
    split = data["e4_time_split"]
    selected = split.loc[
        split.record_type.eq("m_past_stratum")
        & split.variant.eq("latest_unit_year")
        & split.cutoff_year.eq(2015)
        & split.endpoint.isin(ENDPOINTS)
    ].copy()
    out = {}
    for endpoint in ENDPOINTS:
        rows = []
        frame = selected.loc[selected.endpoint.eq(endpoint)]
        for m_past, supported_m in (("m=1", True), ("m>=2", False)):
            stratum = frame.loc[frame.m_past.eq(m_past)].set_index("display_status")
            supported = stratum.loc["supported"]
            uncertain = stratum.loc["uncertain"]
            rows.append(
                {
                    "Row_title": f"{m_past} | supported/uncertain",
                    "m_past": m_past,
                    "supported_agreement": float(supported.supported_agreement_rate),
                    "uncertain_agreement": float(uncertain.uncertain_agreement_rate),
                    "mean_predicted_probability": float(
                        supported.mean_p_displayed
                    ),
                    "supported_n_predictions": int(supported.supported_n),
                    "uncertain_n_predictions": int(uncertain.uncertain_n),
                }
            )
        out[endpoint] = pd.DataFrame(rows)
    return out


def grouped_replicates(
    frame: pd.DataFrame, value: str, layouts: list[str]
) -> pd.DataFrame:
    pieces = []
    for layout in layouts:
        table = frame.loc[frame.layout.eq(layout)].pivot(
            index="target_label", columns="seed", values=value
        )
        table = table.reindex(columns=sorted(table.columns))
        table.columns = [
            f"{layout}__r{index:02d}" for index in range(1, len(table.columns) + 1)
        ]
        pieces.append(table)
    order = ["CA II", "Thrombin", "hERG"]
    return pd.concat(pieces, axis=1).reindex(order).rename_axis("Row_title").reset_index()


def build_5a(data: dict) -> pd.DataFrame:
    recall = data["cgm_figure_3a"].copy()
    recall = recall.loc[
        recall.layout.isin(["projection", "minimum_displacement", "neighbor_refined"])
    ]
    if not recall.groupby(["target_label", "layout"]).seed.nunique().eq(10).all():
        raise AssertionError("Figure 3A source must contain ten seeds per layout")
    return grouped_replicates(
        recall,
        "mean_morgan_neighbor_recall_k10",
        ["projection", "minimum_displacement", "neighbor_refined"],
    )


def build_5b(data: dict) -> pd.DataFrame:
    source = data["e4_layout_baselines"]
    comparisons = source.loc[
        source.row_type.eq("paired_t_interval")
        & source.metric.eq("mean_e_grid")
        & source.comparison.isin(
            [
                "neighbor_refined_minus_minimum_displacement",
                "las_random_minus_neighbor_refined",
                "las_projection_minus_neighbor_refined",
                "flas_random_minus_neighbor_refined",
            ]
        )
    ].copy()
    comparison_order = [
        "neighbor_refined_minus_minimum_displacement",
        "las_random_minus_neighbor_refined",
        "las_projection_minus_neighbor_refined",
        "flas_random_minus_neighbor_refined",
    ]
    seed_rows = source.loc[
        source.row_type.eq("seed_layout")
        & source.layout.isin(
            ["neighbor_refined", "minimum_displacement", "las_random", "las_projection", "flas_random"]
        )
    ]
    refined = (
        seed_rows.loc[seed_rows.layout.eq("neighbor_refined")]
        .groupby(["target", "target_label"], as_index=False)
        .mean_e_grid.mean()
        .rename(columns={"mean_e_grid": "refined_mean_e_grid"})
    )
    labels = {
        "neighbor_refined_minus_minimum_displacement": (
            "minimum_displacement",
            "Minimum-displacement grid",
        ),
        "las_random_minus_neighbor_refined": (
            "las_random",
            "LAS (random initialization)",
        ),
        "las_projection_minus_neighbor_refined": (
            "las_projection",
            "LAS (projection initialization)",
        ),
        "flas_random_minus_neighbor_refined": (
            "flas_random",
            "FLAS (random initialization)",
        ),
    }
    rows = []
    for target, target_label in (
        ("chembl205", "CA II"),
        ("chembl204", "Thrombin"),
        ("chembl240", "hERG"),
    ):
        for comparison in comparison_order:
            record = comparisons.loc[
                comparisons.target.eq(target)
                & comparisons.comparison.eq(comparison)
            ].iloc[0]
            other_layout, other_label = labels[comparison]
            difference = float(record.paired_difference_mean)
            ci_low = float(record.t_ci_low)
            ci_high = float(record.t_ci_high)
            if comparison == "neighbor_refined_minus_minimum_displacement":
                # Minimum displacement stores refined - other; display other - refined.
                difference = -difference
                ci_low, ci_high = -ci_high, -ci_low
            rows.append(
                {
                    "target_id": target,
                    "target_label": target_label,
                    "comparison": comparison,
                    "other_layout": other_layout,
                    "other_layout_label": other_label,
                    "refined_mean_e_grid": float(
                        refined.loc[
                            refined.target.eq(target), "refined_mean_e_grid"
                        ].iloc[0]
                    ),
                    "other_minus_refined_e_grid": difference,
                    "ci_low": ci_low,
                    "ci_high": ci_high,
                    "n_seeds": int(record.n_seeds),
                    "ci_level": float(record.ci_level),
                    "ci_method": record.ci_method,
                }
            )
    return pd.DataFrame(rows)


def build_5d(data: dict) -> pd.DataFrame:
    layouts = data["e2c_layout_by_target"]
    selected = layouts.loc[
        layouts.included.astype(bool) & layouts.layout.eq("neighbor_refined")
    ].copy()
    rows = []
    for endpoint in ENDPOINTS:
        frame = selected.loc[selected.endpoint.eq(endpoint)].sort_values(
            ["rho", "target"]
        ).reset_index(drop=True)
        for rank, record in enumerate(frame.itertuples(index=False), start=1):
            rows.append(
                {
                    "row_type": "sampled_dataset",
                    "endpoint": endpoint,
                    "dataset_rank": rank,
                    "target_id": record.target.upper(),
                    "target_label": "",
                    "n_entries": int(record.n),
                    "rho": float(record.rho),
                    "rho_ci_low": float(record.rho_ci_low),
                    "rho_ci_high": float(record.rho_ci_high),
                    "bootstrap_valid_count": int(record.bootstrap_valid_count),
                }
            )
    reference = data["e1b_neighbor_sar_summary"].loc[
        lambda frame: frame.layout.eq("neighbor_refined")
    ]
    reference_start = int(len(selected) + 1)
    for index, record in enumerate(
        reference.itertuples(index=False), start=reference_start
    ):
        rows.append(
            {
                "row_type": "detailed_reference",
                "endpoint": "reference",
                "dataset_rank": index,
                "target_id": record.target.upper(),
                "target_label": record.target_label,
                "n_entries": int(record.n),
                "rho": float(record.rho_seed42),
                "rho_ci_low": float(record.rho_ci_low),
                "rho_ci_high": float(record.rho_ci_high),
                "bootstrap_valid_count": int(record.bootstrap_valid_count),
            }
        )
    out = pd.DataFrame(rows)
    if len(out.loc[out.row_type.eq("sampled_dataset")]) != 77:
        raise AssertionError("5d must contain 77 sampled datasets")
    return out


def build_5e(data: dict) -> pd.DataFrame:
    source = data["e1b_core_by_group_size"]
    source = source.loc[source.target.eq("overall")].copy()
    rows = []
    for stratum in ("2", "3-4", "5-9", ">=10"):
        for rule in ("most_frequent", "coverage", "fmcs"):
            selected = source.loc[
                source.rule.eq(rule) & source.size_stratum.eq(stratum)
            ]
            groups = int(selected.groups.sum())
            no_core = int(selected.no_core_groups.sum())
            target_met = int(selected.target_met_groups.sum())
            strict = int(selected.strict_majority_groups.sum())
            rows.append(
                {
                    "Row_title": f"{rule} | {stratum}",
                    "rule": rule,
                    "size_stratum": stratum,
                    "coverage_ge_0_8": target_met,
                    "coverage_gt_0_5_lt_0_8": strict - target_met,
                    "coverage_le_0_5": groups - no_core - strict,
                    "no_core": no_core,
                    "groups": groups,
                    "mean_core_heavy_atoms": float(
                        selected.mean_core_heavy_atoms.mean()
                    ),
                }
            )
    out = pd.DataFrame(rows)
    counts = out[
        [
            "coverage_ge_0_8",
            "coverage_gt_0_5_lt_0_8",
            "coverage_le_0_5",
            "no_core",
        ]
    ].sum(axis=1)
    if not counts.eq(out.groups).all():
        raise AssertionError("Coverage bins do not sum to the group count")
    if not out.groupby("rule", sort=False)["groups"].sum().eq(480).all():
        raise AssertionError("Each core rule must contain 480 groups")
    return out


def build_posthoc_sigma(data: dict) -> pd.DataFrame:
    forest = data["e2_forest"]
    sigma = data["e2_sigma_doc"]
    all_rows = forest.loc[forest.population.eq("all")].copy()
    sig = sigma.loc[
        sigma.row_type.eq("target") & sigma.measurement_type.eq("all"),
        ["target", "endpoint", "sigma_doc"],
    ].rename(columns={"target": "target_id"})
    values = all_rows[[
        "target_id", "endpoint", "g_minus_g_fixed_log_loss"
    ]].merge(sig, on=["target_id", "endpoint"], validate="one_to_one")
    values["abs_sigma_minus_0_64"] = (values.sigma_doc - 0.64).abs()
    values["log_loss_g_fixed_minus_g"] = -values.g_minus_g_fixed_log_loss
    values.insert(0, "row_type", "target")
    values = values.drop(columns="g_minus_g_fixed_log_loss")

    summaries = []
    for endpoint in ENDPOINTS:
        selected = values.loc[values.endpoint.eq(endpoint)]
        rho = selected[["abs_sigma_minus_0_64", "log_loss_g_fixed_minus_g"]].corr(
            method="spearman"
        ).iloc[0, 1]
        summaries.append(
            {
                "row_type": "spearman",
                "target_id": "endpoint_summary",
                "endpoint": endpoint,
                "n_targets": int(len(selected)),
                "spearman_rho": float(rho),
            }
        )
    return pd.concat([values, pd.DataFrame(summaries)], ignore_index=True)


def atr_entry_frame(data: dict) -> tuple[pd.DataFrame, tuple[str, str], int]:
    cached = data["e2_curated_cache"]
    records = cached[("CHEMBL5024", "Ki")]
    units = (
        records.groupby(
            ["molecule_identity_key", "document_chembl_id"], as_index=False
        )
        .pchembl_value.median()
    )
    rows = []
    for entry, group in units.groupby("molecule_identity_key", sort=True):
        group = group.sort_values("document_chembl_id")
        document_values = ";".join(
            f"{row.document_chembl_id}={row.pchembl_value:.15g}"
            for row in group.itertuples(index=False)
        )
        distinct = np.sort(group.pchembl_value.round(10).unique())
        rows.append(
            {
                "entry_id": entry,
                "document_values": document_values,
                "n_documents": int(len(group)),
                "n_distinct_values": int(len(distinct)),
                "value_min": float(group.pchembl_value.min()),
                "value_max": float(group.pchembl_value.max()),
                "value_range": float(group.pchembl_value.max() - group.pchembl_value.min()),
                "range_le_0_1": bool(len(distinct) >= 2 and np.ptp(distinct) <= 0.1),
            }
        )
    entries = pd.DataFrame(rows)

    memberships = {
        document: set(frame.molecule_identity_key)
        for document, frame in units.groupby("document_chembl_id")
    }
    pair_counts = [
        (len(memberships[a] & memberships[b]), a, b)
        for a, b in combinations(sorted(memberships), 2)
    ]
    shared, doc_a, doc_b = max(pair_counts, key=lambda item: (item[0], item[1], item[2]))
    return entries, (doc_a, doc_b), shared


def build_posthoc_atr(data: dict) -> pd.DataFrame:
    entries, pair, shared = atr_entry_frame(data)
    multi = entries.loc[entries.n_distinct_values.ge(2)]
    rows = entries.assign(row_type="entry").to_dict("records")

    summary = {
        "row_type": "summary",
        "n_entries": int(len(entries)),
        "n_multi_value_entries": int(len(multi)),
        "n_multi_value_range_le_0_1": int(multi.range_le_0_1.sum()),
        "fraction_multi_value_range_le_0_1": float(multi.range_le_0_1.mean()),
    }
    rows.append(summary)
    doc_a, doc_b = pair
    rows.append(
        {
            "row_type": "document_pair",
            "document_a": doc_a,
            "patent_id_a": ATR_DOCUMENTS.get(doc_a, (None, None))[0],
            "year_a": ATR_DOCUMENTS.get(doc_a, (None, None))[1],
            "document_b": doc_b,
            "patent_id_b": ATR_DOCUMENTS.get(doc_b, (None, None))[0],
            "year_b": ATR_DOCUMENTS.get(doc_b, (None, None))[1],
            "n_shared_entries": int(shared),
        }
    )
    for document, (patent, year) in ATR_DOCUMENTS.items():
        rows.append(
            {
                "row_type": "document_metadata",
                "document_a": document,
                "patent_id_a": patent,
                "year_a": year,
            }
        )
    return pd.DataFrame(rows)


def series(
    name: str,
    columns: list[str],
    frame: pd.DataFrame,
    filter_values: dict | None = None,
    count_mode: str = "row",
) -> dict:
    selected = frame
    if filter_values:
        selected = selected.loc[
            selected[filter_values["column"]].eq(filter_values["value"])
        ]
    if count_mode == "cells":
        count = int(selected.loc[:, columns].notna().to_numpy().sum())
    elif count_mode == "row":
        count = int(selected.loc[:, columns].notna().any(axis=1).sum())
    else:
        raise ValueError(f"Unknown count mode: {count_mode}")
    output = {
        "name": name,
        "columns": columns,
        "expected_count": count,
        "count_mode": count_mode,
    }
    if filter_values:
        output["filter"] = filter_values
    return output


def source_refs(paths: dict[str, Path], keys: list[str]) -> list[str]:
    return [rel(paths[key]) for key in keys]


def make_manifest(
    paths: dict[str, Path],
    panels: list[dict],
    posthoc: list[dict],
    atr_entries: int,
    atr_multi: int,
    atr_fraction: float,
    sigma_rho: dict[str, float],
    r1_validation: dict,
) -> dict:
    return {
        "schema_version": 1,
        "generated_by": rel(SOURCE / "f0_build.py"),
        "task": "F0 Prism-ready data for manuscript Figures 3-5",
        "data_directory": rel(DATA_OUT),
        "sources": {
            rel(path): {"sha256": sha256(path)}
            for path in paths.values()
        },
        "panels": panels,
        "posthoc": posthoc,
        "reused_panels": [
            {
                "panel": "5c",
                "source": rel(paths["cgm_figure_3c"]),
                "source_sha256": sha256(paths["cgm_figure_3c"]),
                "rows": int(len(pd.read_csv(paths["cgm_figure_3c"]))),
                "instruction": "Reuse the frozen Figure_3C panel unchanged; no new CSV.",
            }
        ],
        "validation_targets": {
            "hERG_shared_compounds": 30,
            "hERG_median_offset": 3.47,
            "core_groups_per_rule": 480,
            "posthoc_atr_entries_with_ge_2_distinct_document_values": atr_entries,
            "posthoc_atr_fraction_range_le_0_1": atr_fraction,
            "posthoc_sigma_spearman": sigma_rho,
            "expected_fraction_range_le_0_1": 0.732,
            **r1_validation,
        },
        "prism": {
            "page_font": "Arial",
            "axis_and_tick_pt": 8,
            "line_pt": 0.75,
            "symbol_pt": 3.0,
            "background": "white",
            "target_colors": TARGET_COLORS,
            "endpoint_colors": ENDPOINT_COLORS,
            "display_status_colors": DISPLAY_STATUS_COLORS,
            "comparison_colors": COMPARISON_COLORS,
            "layout_colors": LAYOUT_COLORS,
            "hide_legends": True,
            "legend_owner": "F stage layout or final multi-panel assembly",
            "note": "No Prism document was opened, modified, or saved by F0.",
        },
    }


def gui_plan() -> str:
    return f"""# F0 Prism 操作清单：稿件图 3–5

仅在 F 阶段由人工/Codex 界面操作读取本清单。F0 生成了 CSV 与清单，但未打开、复制、修改或保存任何 Prism 工程。

## 通用设置

1. 每个 CSV 作为独立数据表导入；表型、图型、行标题列和系列列见逐图说明。不要把表中空值改为 0。
2. 页面为白底、Arial。轴标题与刻度 8 pt；数据线、误差线、柱边框、均值线、参考线 0.75 pt；散点符号 3 pt。若最终导出后小于 6 pt 或 0.5 pt，再在新图副本中放大，不改数据。
3. 图标题、图例、系列标签和面板字母均隐藏；在最终多面板布局中统一添加。不得增加 SD、SEM 或显著性标记。
4. 端点颜色固定为 IC50 `{ENDPOINT_COLORS['IC50']}`、Ki `{ENDPOINT_COLORS['Ki']}`。显示状态颜色固定为 supported `{DISPLAY_STATUS_COLORS['supported']}`、uncertain `{DISPLAY_STATUS_COLORS['uncertain']}`。三靶点颜色保持 CA II `{TARGET_COLORS['CA II']}`、Thrombin `{TARGET_COLORS['Thrombin']}`、hERG `{TARGET_COLORS['hERG']}`。
5. 每种布局颜色固定为 projection `{LAYOUT_COLORS['projection']}`、minimum_displacement `{LAYOUT_COLORS['minimum_displacement']}`、neighbor_refined `{LAYOUT_COLORS['neighbor_refined']}`、las_projection `{LAYOUT_COLORS['las_projection']}`、las_random `{LAYOUT_COLORS['las_random']}`、flas_random `{LAYOUT_COLORS['flas_random']}`、random `{LAYOUT_COLORS['random']}`。
6. R0 旧表及旧清单已在 `data/archive_R0/` 和 `source/archive_R0/` 留档；本清单只使用根目录下新表。

## 图 3a — `3a_B0_disagreement.csv`

- 表型：Grouped；`target_id` 为行标题；`IC50` 与 `Ki` 为两个系列；每格 1 个值。
- 图型：grouped bars。页面 3.3 × 2.3 in；柱从 0 起，边框 0.75 pt；不显示柱值；不另加均值线或误差线；隐藏图例。
- 排序：按 CSV 行顺序。若某端点没有该数据集，保留空格，不画零柱。
- 轴：X 为空分类轴，标签过密时横排目标 ID 并缩小绘图区，不删行；Y「B0 disagreement rate」[0, 1]。
- 系列：IC50 `{ENDPOINT_COLORS['IC50']}`，预期点数 {int(pd.read_csv(DATA_OUT / '3a_B0_disagreement.csv').IC50.notna().sum())}；Ki `{ENDPOINT_COLORS['Ki']}`，预期点数 {int(pd.read_csv(DATA_OUT / '3a_B0_disagreement.csv').Ki.notna().sum())}。
- 参考线：无。

## 图 3b — `3b_IC50.csv`、`3b_Ki.csv`

- 表型：XY，每格 1 个值；每个端点单独导入，保持 `target_rank` 升序，即以 G AUC 从高到低排列。
- 图型：horizontal forest plot。每端点页面 2.7 × 4.4 in。G 为端点色实心圆点；D 为 `#666666` 空心圆点；B-range 为 `#666666` 空心方形点。三者均带水平 95% CI，误差线 0.75 pt，无连接线；隐藏图例。
- 轴：X「AUC」[0, 1]；Y 按 `target_rank` 排列 target ID。垂直参考线 x = 0.5，`#666666`，dashed，不显示标签和端点符号。
- IC50 系列：G `g_auc` + `g_ci_low/g_ci_high`，D `d_auc` + `d_auc_ci_low/d_auc_ci_high`，B-range `b_range_auc` + CI；各 {len(pd.read_csv(DATA_OUT / '3b_IC50.csv'))} 点。
- Ki 系列：同上；各 {len(pd.read_csv(DATA_OUT / '3b_Ki.csv'))} 点。`*_n_predictions` 只作数据表信息，不画图。
- 参考线：仅 x = 0.5。

## 图 3c — `3c_IC50.csv`、`3c_Ki.csv`

- 表型：XY；`row_type=target` 的行画每个数据集，`row_type=summary` 行的 `median_*` 与 `ci_low_*`/`ci_high_*` 只画中位数及 bootstrap 95% CI。
- 图型：horizontal forest。每端点页面 2.7 × 4.4 in；目标点和区间 3 pt/0.75 pt；无连接线、柱值、SD 或 SEM；隐藏图例。
- 数值方向：`g_record_minus_g`、`g_global_minus_g`、`g_fixed_minus_g`、`x_minus_g` 均为“替代模型 − G”；正值表示 G 更好。
- X「Alternative minus G log loss」[−0.30, 0.30]；Y 按 `target_rank` 排列。
- 系列颜色：G-record `#995500`；G-global `#CC79A7`；G-fixed `#D55E00`；X `{LAYOUT_COLORS['las_projection']}`。
- IC50：每个目标系列 {int((pd.read_csv(DATA_OUT / '3c_IC50.csv').row_type == 'target').sum())} 点，X 仅 15 点；Ki：每个目标系列 {int((pd.read_csv(DATA_OUT / '3c_Ki.csv').row_type == 'target').sum())} 点，X 仅 11 点。无标记对的数据集保持 X 空值。
- 参考线：x = 0，`#666666`，dashed；不显示标签和端点符号。

## 图 4a — `4a_calibration_IC50.csv`、`4a_calibration_Ki.csv`

- 表型：XY；`G_probability`/`G_observed_agreement` 与 `B0_probability`/`B0_observed_agreement` 各为一组 X/Y；每个端点单独成图。仅含 n ≥ 8 的显示类别概率区间。
- 图型：G 用圆点加实线；B0 在概率 1 处用同色空心圆，仅一个点且不连线。页面 2.3 × 2.15 in；隐藏图例。
- IC50 颜色 `{ENDPOINT_COLORS['IC50']}`；Ki 颜色 `{ENDPOINT_COLORS['Ki']}`。
- 轴：X「Predicted probability of the displayed class」[0, 1]；Y「Observed agreement」[0, 1]。
- 参考线：y = x，`#666666`，dashed；不显示标签和端点符号。
- 预期点：IC50 G {int(pd.read_csv(DATA_OUT / '4a_calibration_IC50.csv').G_probability.notna().sum())} 点 + B0 1 点；Ki G {int(pd.read_csv(DATA_OUT / '4a_calibration_Ki.csv').G_probability.notna().sum())} 点 + B0 1 点。

## 图 4b — `4b_hERG_patent_pair.csv`

- 表型：XY，每格 1 个值；X=`us_10131663_b2`，Y=`us_10590130_b2`；30 个共享化合物。
- 图型：scatter。页面 2.5 × 2.5 in；点为 `{TARGET_COLORS['CA II']}`；圆点 3 pt、filled；无连接线；隐藏图例。
- 轴：X「US-10131663-B2 pChEMBL」[4.4, 9.3]；Y「US-10590130-B2 pChEMBL」[4.4, 9.6]。
- 参考线：y = x，`#666666`，dashed；y = x + 3.47，`{ENDPOINT_COLORS['IC50']}`，dotted；均不显示标签、范围标签或端点符号。

## 图 4c — `4c_IC50.csv`、`4c_Ki.csv`

- 表型：XY，每个数据集一行；X=`supported_agreement`/`uncertain_agreement`，Y=`dataset_rank`；配对双侧误差线仅使用 Wilson `*_ci_low`/`*_ci_high`。
- 图型：horizontal paired dot-and-line plot，supported 与 uncertain 两点之间用细灰线连接。每个端点页面 2.8 × 4.2 in；点 3 pt、误差线 0.75 pt；隐藏图例。
- 轴：X「Agreement with later documents」[0, 1]；Y 按 `dataset_rank` 排列，目标 ID 可另行标注。supported `{DISPLAY_STATUS_COLORS['supported']}`，uncertain `{DISPLAY_STATUS_COLORS['uncertain']}`。
- 数据范围：IC50 {len(pd.read_csv(DATA_OUT / '4c_IC50.csv'))} 个数据集；Ki {len(pd.read_csv(DATA_OUT / '4c_Ki.csv'))} 个；合计 33 个。全部为 2015 年、`latest_unit_year`，不使用 `ALL` 行。
- 参考线：无。

## 图 4d — `4d_IC50.csv`、`4d_Ki.csv`

- 表型：Grouped；`Row_title` 为 `m=1 | supported/uncertain` 和 `m>=2 | supported/uncertain`；每格 1 个值。
- 图型：grouped bars + open points。supported/uncertain 一致率为两根柱；`mean_predicted_probability` 用 `#172536` 空心菱形点叠加。每端点页面 2.8 × 2.3 in；显示 2 位小数柱值；隐藏图例。
- 轴：X 为空分类轴，按 CSV 行顺序；Y「Agreement / mean predicted probability」[0, 1]。
- 颜色：supported `{DISPLAY_STATUS_COLORS['supported']}`；uncertain `{DISPLAY_STATUS_COLORS['uncertain']}`。
- 参考线：无。

## 图 5a — `5a_recall.csv`

- 表型：Grouped，每格 10 个 seed 值；`Row_title` 为三靶点；每个布局使用 10 个 `__rNN` 列。
- 图型：bars_and_points。页面 3.3 × 2.55 in；每靶点 3 根柱；柱显示均值，点显示原始 seed；柱边框和均值线 0.75 pt；符号 3 pt；不显示柱值；不另加误差线；隐藏图例。
- 轴：X 为空分类轴；Y「Fingerprint-neighbor recall (k = 10)」[0, 1]。
- 系列与颜色：projection、minimum_displacement、neighbor_refined 按通用布局色；各 30 点（3 靶点 × 10 seed）。无参考线。

## 图 5b — `5b_layout_contrasts.csv`

- 表型：XY，12 行；`target_label` + `other_layout_label` 为行标签，X=`other_minus_refined_e_grid`，使用 `ci_low`/`ci_high` 的 95% t 区间，`n_seeds` 为 10。
- 图型：horizontal forest。页面 3.3 × 2.9 in；点 3 pt、区间 0.75 pt、无连接线；隐藏图例。
- 方向：所有值统一为“其他布局 − 精炼网格”；源文件中的“精炼 − 最小位移”已取反。正值表示精炼网格的活性邻居更相近。
- 轴：X「Other layout minus refined e_grid」[−0.05, 0.05]；Y 按 CSV 顺序。
- 颜色：minimum_displacement `{COMPARISON_COLORS['minimum_displacement']}`；LAS random `{COMPARISON_COLORS['las_random']}`；LAS projection `{COMPARISON_COLORS['las_projection']}`；FLAS `{COMPARISON_COLORS['flas_random']}`。
- 参考线：x = 0，`#666666`，dashed；不显示标签和端点符号。

## 图 5c — 冻结图 3C

- 直接复用 `chemgridmap/paper/output/reliability_v1/figure_data/Figure_3C_CHEMBL81688_neighbors.csv` 生成的现有面板；不新建数据表，不改数值、标签或布局。
- 保留现有字体、线宽和点样式；最终只需缩放到图 5 的面板槽位并检查可读性。隐藏 P 图标题和图例不影响原有 Q/S/G/F 标记。

## 图 5d — `5d_layout_correlation.csv`

- 表型：XY；`dataset_rank` 为 Y，`rho` 为 X，`rho_ci_low`/`rho_ci_high` 为非对称 bootstrap 95% CI。
- 图型：horizontal forest。页面 3.3 × 4.5 in。`sampled_dataset` 为实心圆点，共 77 个；`detailed_reference` 为 seed 42 的空心菱形参考点，共 3 个，位于第 78–80 行。点 3 pt、区间 0.75 pt；隐藏图例。
- 排序：先 IC50（48 个）后 Ki（29 个）；每个端点内按 ρ 从低到高。三个详细靶点单独列在末尾，不混入 77 个总体数据集。
- 轴：X「Spearman ρ(recall, e_grid − e_fp)」[−0.30, 0.20]；Y 按 `dataset_rank`。sampled 点用端点色：IC50 `{ENDPOINT_COLORS['IC50']}`、Ki `{ENDPOINT_COLORS['Ki']}`；参考点分别用 CA II `{TARGET_COLORS['CA II']}`、Thrombin `{TARGET_COLORS['Thrombin']}`、hERG `{TARGET_COLORS['hERG']}`。
- 参考线：x = 0，`#666666`，dashed；不显示标签和端点符号。

## 图 5e — `5e_core_coverage_by_group_size.csv`

- 表型：Grouped；`Row_title` 为 rule | size stratum，共 12 行；四个覆盖率分类列为堆叠系列，`mean_core_heavy_atoms` 为另一组点。
- 图型：stacked bars + open diamond points。页面 3.3 × 2.7 in；柱从 0 起、0.75 pt 边框；堆叠柱左轴，平均重原子数右轴。隐藏图例。
- 轴：X 为空分类轴，顺序为 2、3–4、5–9、≥10，每档 most_frequent→coverage→fmcs；左 Y「Number of groups」[0, 220]；右 Y「Mean core heavy atoms」[0, 25]。
- 堆叠颜色：≥0.8 `#006354`；0.5<c<0.8 `#1D5689`；≤0.5 `#666666`；No core `#D9D9D9`。平均重原子数点用 `#172536`。
- 数据约束：每条规则四档组数之和必须为 480；No core 不并入 ≤0.5。无参考线。旧 FMCS 时间预算图已移出本主图。
"""


def build_all() -> dict:
    paths = source_paths()
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing F0 sources: {missing}")
    DATA_OUT.mkdir(parents=True, exist_ok=True)
    PLAN_OUT.parent.mkdir(parents=True, exist_ok=True)
    data = load_sources(paths)

    table_3a = build_3a(data)
    table_3b = build_3b(data)
    table_3c = build_3c(data, table_3b)
    table_4a = build_4a(data)
    table_4b = build_4b(data)
    table_4c = build_4c(data)
    table_4d = build_4d(data)
    table_5a = build_5a(data)
    table_5b = build_5b(data)
    table_5d = build_5d(data)
    table_5e = build_5e(data)
    table_sigma = build_posthoc_sigma(data)
    table_atr = build_posthoc_atr(data)

    write_csv(table_3a, "3a_B0_disagreement.csv")
    for endpoint, frame in table_3b.items():
        write_csv(frame, f"3b_{endpoint}.csv")
    for endpoint, frame in table_3c.items():
        write_csv(frame, f"3c_{endpoint}.csv")
    for endpoint, frame in table_4a.items():
        write_csv(frame, f"4a_calibration_{endpoint}.csv")
    write_csv(table_4b, "4b_hERG_patent_pair.csv")
    for endpoint, frame in table_4c.items():
        write_csv(frame, f"4c_{endpoint}.csv")
    for endpoint, frame in table_4d.items():
        write_csv(frame, f"4d_{endpoint}.csv")
    write_csv(table_5a, "5a_recall.csv")
    write_csv(table_5b, "5b_layout_contrasts.csv")
    write_csv(table_5d, "5d_layout_correlation.csv")
    write_csv(table_5e, "5e_core_coverage_by_group_size.csv")
    write_csv(table_sigma, "posthoc_sigma_gain.csv")
    write_csv(table_atr, "posthoc_atr_patent_family.csv")
    PLAN_OUT.write_text(gui_plan(), encoding="utf-8")

    panels = {
        "3a": {
            "panel": "3a",
            "files": ["3a_B0_disagreement.csv"],
            "table_type": "Grouped",
            "row_title_column": "target_id",
            "rows": len(table_3a),
            "source_files": source_refs(paths, ["e2_metrics_by_target"]),
            "series": [
                series("IC50", ["IC50"], table_3a),
                series("Ki", ["Ki"], table_3a),
            ],
        },
        "3b": {
            "panel": "3b",
            "files": ["3b_IC50.csv", "3b_Ki.csv"],
            "table_type": "XY",
            "rows": {endpoint: len(frame) for endpoint, frame in table_3b.items()},
            "source_files": source_refs(paths, ["e2_forest", "e4_comparators"]),
            "series": {
                endpoint: [
                    series("G AUC", ["g_auc"], frame),
                    series("G 95% CI", ["g_ci_low", "g_ci_high"], frame),
                    series("D AUC", ["d_auc"], frame),
                    series("D 95% CI", ["d_auc_ci_low", "d_auc_ci_high"], frame),
                    series("B-range AUC", ["b_range_auc"], frame),
                    series(
                        "B-range 95% CI",
                        ["b_range_auc_ci_low", "b_range_auc_ci_high"],
                        frame,
                    ),
                ]
                for endpoint, frame in table_3b.items()
            },
        },
        "3c": {
            "panel": "3c",
            "files": ["3c_IC50.csv", "3c_Ki.csv"],
            "table_type": "XY",
            "rows": {endpoint: len(frame) for endpoint, frame in table_3c.items()},
            "source_files": source_refs(
                paths, ["e2_metrics_by_target", "e2_flagged_pairs_census"]
            ),
            "series": {
                endpoint: [
                    series("G-record minus G", ["g_record_minus_g"], frame),
                    series("G-global minus G", ["g_global_minus_g"], frame),
                    series("G-fixed minus G", ["g_fixed_minus_g"], frame),
                    series("X minus G", ["x_minus_g"], frame),
                ]
                for endpoint, frame in table_3c.items()
            },
        },
        "4a": {
            "panel": "4a",
            "files": ["4a_calibration_IC50.csv", "4a_calibration_Ki.csv"],
            "table_type": "XY",
            "rows": {endpoint: len(frame) for endpoint, frame in table_4a.items()},
            "source_files": source_refs(paths, ["e4_calibration"]),
            "series": {
                endpoint: [
                    series(
                        "G",
                        ["G_probability", "G_observed_agreement"],
                        frame,
                    ),
                    series(
                        "B0",
                        ["B0_probability", "B0_observed_agreement"],
                        frame,
                    ),
                ]
                for endpoint, frame in table_4a.items()
            },
        },
        "4b": {
            "panel": "4b",
            "files": ["4b_hERG_patent_pair.csv"],
            "table_type": "XY",
            "rows": len(table_4b),
            "source_files": source_refs(
                paths, ["cgm_hERG_retained", "e1a_flagged_pairs"]
            ),
            "series": [
                series(
                    "US-10590130-B2 vs US-10131663-B2",
                    ["us_10131663_b2", "us_10590130_b2"],
                    table_4b,
                )
            ],
            "references": [
                {"equation": "y=x", "style": "dashed"},
                {"equation": "y=x+3.47", "style": "dotted"},
            ],
        },
        "4c": {
            "panel": "4c",
            "files": ["4c_IC50.csv", "4c_Ki.csv"],
            "table_type": "XY",
            "rows": {endpoint: len(frame) for endpoint, frame in table_4c.items()},
            "source_files": source_refs(paths, ["e4_time_split"]),
            "series": {
                endpoint: [
                    series("Supported agreement", ["supported_agreement"], frame),
                    series(
                        "Supported 95% Wilson CI",
                        ["supported_ci_low", "supported_ci_high"],
                        frame,
                    ),
                    series("Uncertain agreement", ["uncertain_agreement"], frame),
                    series(
                        "Uncertain 95% Wilson CI",
                        ["uncertain_ci_low", "uncertain_ci_high"],
                        frame,
                    ),
                ]
                for endpoint, frame in table_4c.items()
            },
        },
        "4d": {
            "panel": "4d",
            "files": ["4d_IC50.csv", "4d_Ki.csv"],
            "table_type": "Grouped",
            "row_title_column": "Row_title",
            "rows": {endpoint: len(frame) for endpoint, frame in table_4d.items()},
            "source_files": source_refs(paths, ["e4_time_split"]),
            "series": {
                endpoint: [
                    series("Supported agreement", ["supported_agreement"], frame),
                    series("Uncertain agreement", ["uncertain_agreement"], frame),
                    series(
                        "Mean predicted probability",
                        ["mean_predicted_probability"],
                        frame,
                    ),
                ]
                for endpoint, frame in table_4d.items()
            },
        },
        "5a": {
            "panel": "5a",
            "files": ["5a_recall.csv"],
            "table_type": "Grouped",
            "row_title_column": "Row_title",
            "rows": len(table_5a),
            "source_files": source_refs(paths, ["cgm_figure_3a"]),
            "series": [],
        },
        "5b": {
            "panel": "5b",
            "files": ["5b_layout_contrasts.csv"],
            "table_type": "XY",
            "rows": len(table_5b),
            "source_files": source_refs(paths, ["e4_layout_baselines"]),
            "series": [
                series(
                    name,
                    ["other_minus_refined_e_grid"],
                    table_5b,
                    {"column": "comparison", "value": comparison},
                )
                for comparison, name in (
                    ("neighbor_refined_minus_minimum_displacement", "Minimum-displacement - refined"),
                    ("las_random_minus_neighbor_refined", "LAS random - refined"),
                    ("las_projection_minus_neighbor_refined", "LAS projection - refined"),
                    ("flas_random_minus_neighbor_refined", "FLAS - refined"),
                )
            ],
        },
        "5d": {
            "panel": "5d",
            "files": ["5d_layout_correlation.csv"],
            "table_type": "XY",
            "rows": {
                "sampled_datasets": int(table_5d.row_type.eq("sampled_dataset").sum()),
                "detailed_references": int(
                    table_5d.row_type.eq("detailed_reference").sum()
                ),
            },
            "source_files": source_refs(
                paths, ["e2c_layout_by_target", "e1b_neighbor_sar_summary"]
            ),
            "series": [
                series(
                    "Sampled IC50",
                    ["rho", "rho_ci_low", "rho_ci_high"],
                    table_5d,
                    {"column": "endpoint", "value": "IC50"},
                ),
                series(
                    "Sampled Ki",
                    ["rho", "rho_ci_low", "rho_ci_high"],
                    table_5d,
                    {"column": "endpoint", "value": "Ki"},
                ),
                series(
                    "Detailed targets",
                    ["rho", "rho_ci_low", "rho_ci_high"],
                    table_5d,
                    {"column": "row_type", "value": "detailed_reference"},
                ),
            ],
        },
        "5e": {
            "panel": "5e",
            "files": ["5e_core_coverage_by_group_size.csv"],
            "table_type": "Grouped",
            "row_title_column": "Row_title",
            "rows": len(table_5e),
            "source_files": source_refs(paths, ["e1b_core_by_group_size"]),
            "series": [
                series("Coverage >= 0.8", ["coverage_ge_0_8"], table_5e),
                series(
                    "0.5 < coverage < 0.8",
                    ["coverage_gt_0_5_lt_0_8"],
                    table_5e,
                ),
                series("Coverage <= 0.5", ["coverage_le_0_5"], table_5e),
                series("No core", ["no_core"], table_5e),
                series(
                    "Mean core heavy atoms",
                    ["mean_core_heavy_atoms"],
                    table_5e,
                ),
            ],
        },
    }

    # Replicate-column series are generated explicitly from the source order.
    recall_columns = [column for column in table_5a if "__r" in column]
    panels["5a"]["series"] = [
        series(
            "projection",
            [column for column in recall_columns if recall_columns.index(column) < 10],
            table_5a,
            count_mode="cells",
        ),
        series(
            "minimum_displacement",
            [
                column
                for column in recall_columns
                if 10 <= recall_columns.index(column) < 20
            ],
            table_5a,
            count_mode="cells",
        ),
        series(
            "neighbor_refined",
            [
                column
                for column in recall_columns
                if recall_columns.index(column) >= 20
            ],
            table_5a,
            count_mode="cells",
        ),
    ]
    panels = list(panels.values())

    multi_atr = table_atr.loc[table_atr.row_type.eq("entry") & table_atr.n_distinct_values.ge(2)]
    atr_fraction = float(multi_atr.range_le_0_1.mean())
    sigma_rho = (
        table_sigma.loc[table_sigma.row_type.eq("spearman")]
        .set_index("endpoint")
        .spearman_rho.astype(float)
        .to_dict()
    )
    small_time_split = pd.concat(table_4c.values(), ignore_index=True)
    sampled_5d = table_5d.loc[table_5d.row_type.eq("sampled_dataset")]
    herg_las = table_5b.loc[
        table_5b.target_label.eq("hERG")
        & table_5b.comparison.eq("las_random_minus_neighbor_refined"),
        "other_minus_refined_e_grid",
    ].iloc[0]
    r1_validation = {
        "time_split_2015_dataset_rows": int(len(small_time_split)),
        "time_split_2015_supported_gt_uncertain": int(
            small_time_split.supported_agreement.gt(
                small_time_split.uncertain_agreement
            ).sum()
        ),
        "five_b_hERG_las_random_other_minus_refined": float(herg_las),
        "five_d_sampled_datasets": int(len(sampled_5d)),
        "five_d_sampled_ic50": int(sampled_5d.endpoint.eq("IC50").sum()),
        "five_d_sampled_ki": int(sampled_5d.endpoint.eq("Ki").sum()),
        "three_c_flagged_ic50": int(table_3c["IC50"].x_minus_g.notna().sum()),
        "three_c_flagged_ki": int(table_3c["Ki"].x_minus_g.notna().sum()),
        "four_d_rows": int(sum(len(frame) for frame in table_4d.values())),
        "five_e_rows": int(len(table_5e)),
    }
    posthoc = [
        {
            "name": "posthoc_sigma_gain",
            "file": "posthoc_sigma_gain.csv",
            "rows": len(table_sigma),
            "source_files": source_refs(paths, ["e2_forest", "e2_sigma_doc"]),
            "series": [
                series("|sigma_doc - 0.64|", ["abs_sigma_minus_0_64"], table_sigma),
                series(
                    "log loss(G-fixed) - log loss(G)",
                    ["log_loss_g_fixed_minus_g"],
                    table_sigma,
                ),
            ],
        },
        {
            "name": "posthoc_atr_patent_family",
            "file": "posthoc_atr_patent_family.csv",
            "rows": len(table_atr),
            "source_files": source_refs(paths, ["e2_curated_cache", "e2_target_inclusion"]),
            "series": [
                series(
                    "ATR entry values",
                    ["document_values"],
                    table_atr.loc[table_atr.row_type.eq("entry")],
                )
            ],
        },
    ]
    manifest = make_manifest(
        paths,
        panels,
        posthoc,
        int(len(multi_atr)),
        int(multi_atr.range_le_0_1.sum()),
        atr_fraction,
        sigma_rho,
        r1_validation,
    )
    MANIFEST_OUT.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return manifest


if __name__ == "__main__":
    build_all()
