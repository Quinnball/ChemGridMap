"""E5 post hoc evidence-model variants and fold constructors."""
from __future__ import annotations

from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[3]
S = ROOT / "paper/strengthening"
sys.path[:0] = [
    str(S / "source"),
    str(ROOT / "paper"),
    str(ROOT / "src"),
]

import numpy as np
import pandas as pd

from chemgridmap.activity import activity_class
from chemgridmap.evidence import (
    CLASSES,
    class_probabilities,
    document_values,
    fit_between_unit_sd,
)
from e2_models import FOLD_SEED
from e4c_analysis import near_duplicate_units
from reliability_activity import entry_folds, gaussian_prediction

CLASS_BOUNDARIES = (6.0, 7.0)
DOCUMENT_FOLD_SEED = 20261102
SUPPORT_THRESHOLD = 0.8
ANALYSIS_ROLE = "post hoc"
TAU_FLOOR = 0.25
METHODS = ("G", "G-EB", "G-HS", "G-EB-HS", "B1c")


def classify(value: float) -> int:
    return CLASSES.index(activity_class(float(value), *CLASS_BOUNDARIES))


def hs_class(value: float) -> str:
    """The E5 heterogeneity strata use <6, 6-7 and >=7."""
    if value < 6.0:
        return "inactive"
    if value < 7.0:
        return "medium"
    return "active"


def pooled_sigma(doc: pd.DataFrame) -> dict[str, float | int | str]:
    return fit_between_unit_sd(doc)


def heterogeneity_sigmas(
    doc: pd.DataFrame, overall_sigma: float
) -> dict[str, dict[str, float | int | str]]:
    """Pooled within-entry SD by the category of each training-entry mean."""
    moments = doc.groupby("entry").document_value.agg(["size", "mean", "var"])
    moments = moments.loc[moments["size"] >= 2].copy()
    moments["stratum"] = moments["mean"].map(hs_class)
    result = {}
    for stratum in CLASSES:
        selected = moments.loc[moments.stratum == stratum]
        if len(selected) < 5:
            result[stratum] = {
                "sigma": float(overall_sigma),
                "n_entries": int(len(selected)),
                "method": "overall_fallback_lt_5_entries",
            }
            continue
        dfs = selected["size"].to_numpy(float) - 1.0
        variance = float(
            np.sum(selected["var"].to_numpy(float) * dfs) / np.sum(dfs)
        )
        result[stratum] = {
            "sigma": float(np.sqrt(variance)),
            "n_entries": int(len(selected)),
            "method": "class_specific_pooled",
        }
    return result


def empirical_bayes_prior(
    doc: pd.DataFrame, sigma: float
) -> dict[str, float | int | str]:
    """Fit the target-and-fold normal prior using equally weighted entry means."""
    moments = doc.groupby("entry").document_value.agg(["size", "mean"])
    if len(moments) < 1:
        raise ValueError("At least one training entry is needed for G-EB.")
    entry_means = moments["mean"].to_numpy(float)
    unit_counts = moments["size"].to_numpy(float)
    mu0 = float(entry_means.mean())
    between = (
        float(entry_means.var(ddof=1)) if len(entry_means) >= 2 else 0.0
    )
    sampling = float(np.mean(sigma**2 / unit_counts))
    tau2 = max(between - sampling, TAU_FLOOR**2)
    return {
        "mu0": mu0,
        "tau2": tau2,
        "tau": float(np.sqrt(tau2)),
        "raw_between_variance": between,
        "mean_sampling_variance": sampling,
        "n_entries": int(len(moments)),
        "n_units": int(len(doc)),
        "method": (
            "normal_empirical_bayes"
            if len(moments) >= 2
            else "single_entry_tau_floor_fallback"
        ),
    }


def eb_prediction(
    values, sigma: float, prior: dict
) -> tuple[np.ndarray, np.ndarray, float, float]:
    values = np.asarray(values, dtype=float)
    m = len(values)
    if not m:
        raise ValueError("Empirical Bayes prediction needs at least one value.")
    ybar = float(values.mean())
    precision = 1.0 / float(prior["tau2"]) + m / sigma**2
    posterior_var = 1.0 / precision
    posterior_mean = (
        float(prior["mu0"]) / float(prior["tau2"])
        + m * ybar / sigma**2
    ) / precision
    support = class_probabilities(
        [posterior_mean],
        [np.sqrt(posterior_var)],
        np.inf,
        *CLASS_BOUNDARIES,
    )[0]
    predictive = class_probabilities(
        [posterior_mean],
        [np.sqrt(posterior_var + sigma**2)],
        np.inf,
        *CLASS_BOUNDARIES,
    )[0]
    return support, predictive, float(posterior_mean), float(posterior_var)


def hs_prediction(
    values, sigmas: dict[str, dict]
) -> tuple[np.ndarray, np.ndarray, float]:
    values = np.asarray(values, dtype=float)
    m = len(values)
    ybar = float(values.mean())
    sigma = float(sigmas[hs_class(ybar)]["sigma"])
    support = class_probabilities(
        [ybar], [sigma / np.sqrt(m)], np.inf, *CLASS_BOUNDARIES
    )[0]
    predictive = class_probabilities(
        [ybar],
        [sigma * np.sqrt(1.0 + 1.0 / m)],
        np.inf,
        *CLASS_BOUNDARIES,
    )[0]
    return support, predictive, sigma


def eb_hs_prediction(
    values, sigmas: dict[str, dict], prior: dict
) -> tuple[np.ndarray, np.ndarray, float, float, float]:
    values = np.asarray(values, dtype=float)
    sigma = float(sigmas[hs_class(float(values.mean()))]["sigma"])
    support, predictive, posterior_mean, posterior_var = eb_prediction(
        values, sigma, prior
    )
    return support, predictive, posterior_mean, posterior_var, sigma


def b1c_probabilities(
    statistics: pd.DataFrame, displayed: int, m: int
) -> np.ndarray:
    """Laplace-smoothed observed classes from fitted out-of-fold G predictions."""
    m_layer = 1 if int(m) == 1 else 2
    selected = statistics.loc[
        (statistics.displayed == int(displayed))
        & (statistics.m_layer == m_layer)
    ]
    counts = np.ones(3, dtype=float)
    counts += np.bincount(
        selected.observed.to_numpy(int), minlength=3
    ).astype(float)
    return counts / counts.sum()


def b1c_training_statistics(
    frame: pd.DataFrame, held_fold: int
) -> pd.DataFrame:
    """Use only G predictions from the other nine entry folds."""
    selected = frame.loc[
        (frame.method == "G") & (frame.fold != int(held_fold)),
        ["displayed", "observed", "m"],
    ].copy()
    selected["m_layer"] = np.where(selected.m == 1, 1, 2)
    return selected


def document_fold_map(doc: pd.DataFrame) -> pd.Series:
    primary = doc.documents.map(
        lambda value: min(str(value).split(";"))
    )
    documents = np.array(sorted(primary.unique()))
    permuted = np.random.default_rng(DOCUMENT_FOLD_SEED).permutation(documents)
    folds = {document: index % 10 for index, document in enumerate(permuted)}
    return primary.map(folds)


def _project(row: dict) -> dict:
    support = np.asarray(row.pop("_support"), dtype=float)
    probabilities = np.asarray(row.pop("_p"), dtype=float)
    shown = int(row["displayed"])
    row["p0"], row["p1"], row["p2"] = probabilities
    row["p_displayed"] = float(probabilities[shown])
    row["support"] = float(support[shown])
    row["display_status"] = (
        "supported" if row["support"] >= SUPPORT_THRESHOLD else "uncertain"
    )
    row["agree"] = bool(int(row["observed"]) == shown)
    row["analysis_role"] = ANALYSIS_ROLE
    return row


def _remaining_raw(
    source: pd.DataFrame, held_documents: set[str]
) -> tuple[np.ndarray, list[str]]:
    selected = source.loc[
        ~source.document_chembl_id.astype(str).isin(held_documents)
    ]
    values = selected.pchembl_value.to_numpy(float)
    assay_types = selected.assay_type.fillna("missing").astype(str).tolist()
    return values, assay_types


def entry_fold_predictions(
    records: pd.DataFrame, target: str, endpoint: str
) -> pd.DataFrame:
    """Build E2-compatible held-unit predictions for all E5 models."""
    doc = document_values(records, merge_identical=True)
    assignments = entry_folds(doc, seed=FOLD_SEED)
    units_by_entry = dict(tuple(doc.groupby("entry", sort=True)))
    source_by_entry = dict(
        tuple(records.groupby("molecule_identity_key", sort=True))
    )
    rows = []
    for fold in range(10):
        test_entries = {
            entry for entry, assigned in assignments.items() if assigned == fold
        }
        train_doc = doc.loc[~doc.entry.isin(test_entries)]
        sigma_fit = pooled_sigma(train_doc)
        sigma = float(sigma_fit["tau"])
        prior = empirical_bayes_prior(train_doc, sigma)
        sigmas = heterogeneity_sigmas(train_doc, sigma)
        fold_rows = []
        for entry in sorted(test_entries):
            units = units_by_entry[entry]
            source = source_by_entry[entry]
            for held in units.itertuples(index=False):
                held_documents = set(str(held.documents).split(";"))
                remaining_units = units.loc[
                    ~units.documents.map(
                        lambda value: bool(
                            set(str(value).split(";")) & held_documents
                        )
                    )
                ]
                remaining_raw, assay_types = _remaining_raw(
                    source, held_documents
                )
                if remaining_units.empty or not len(remaining_raw):
                    raise ValueError("Held unit leaves no remaining evidence.")
                values = remaining_units.document_value.to_numpy(float)
                displayed = classify(np.median(remaining_raw))
                observed = classify(float(held.document_value))
                base = {
                    "split_design": "entry",
                    "merge_basis": "same_value",
                    "target": target,
                    "endpoint": endpoint,
                    "fold": int(fold),
                    "entry": entry,
                    "held_out_unit": str(held.documents),
                    "m": int(len(values)),
                    "displayed": displayed,
                    "observed": observed,
                    "held_out_value": float(held.document_value),
                    "sigma_doc": sigma,
                    "sigma_method": sigma_fit["method"],
                    "prior_mu0": float(prior["mu0"]),
                    "prior_tau2": float(prior["tau2"]),
                    "prior_tau": float(prior["tau"]),
                    "remaining_center": float(values.mean()),
                    "remaining_values": ";".join(
                        f"{value:.10g}" for value in values
                    ),
                    "remaining_assay_types": ";".join(assay_types),
                }

                g_predictive = gaussian_prediction(values, sigma_fit)
                g_support = class_probabilities(
                    [values.mean()],
                    [sigma / np.sqrt(len(values))],
                    np.inf,
                    *CLASS_BOUNDARIES,
                )[0]
                fold_rows.append(
                    {
                        **base,
                        "method": "G",
                        "_support": g_support,
                        "_p": g_predictive,
                    }
                )

                eb_support, eb_predictive, posterior_mean, posterior_var = (
                    eb_prediction(values, sigma, prior)
                )
                fold_rows.append(
                    {
                        **base,
                        "method": "G-EB",
                        "posterior_mean": posterior_mean,
                        "posterior_variance": posterior_var,
                        "_support": eb_support,
                        "_p": eb_predictive,
                    }
                )

                hs_support, hs_predictive, selected_sigma = hs_prediction(
                    values, sigmas
                )
                fold_rows.append(
                    {
                        **base,
                        "method": "G-HS",
                        "sigma_selected": selected_sigma,
                        "_support": hs_support,
                        "_p": hs_predictive,
                    }
                )

                (
                    ebhs_support,
                    ebhs_predictive,
                    ebhs_mean,
                    ebhs_variance,
                    ebhs_sigma,
                ) = eb_hs_prediction(values, sigmas, prior)
                fold_rows.append(
                    {
                        **base,
                        "method": "G-EB-HS",
                        "posterior_mean": ebhs_mean,
                        "posterior_variance": ebhs_variance,
                        "sigma_selected": ebhs_sigma,
                        "_support": ebhs_support,
                        "_p": ebhs_predictive,
                    }
                )
        rows.extend(fold_rows)

    frame = pd.DataFrame(rows)
    b1c_rows = []
    for row in frame.loc[frame.method == "G"].to_dict("records"):
        statistics = b1c_training_statistics(
            frame, int(row["fold"])
        )
        probabilities = b1c_probabilities(
            statistics, int(row["displayed"]), int(row["m"])
        )
        b1c_rows.append(
            {
                **{key: value for key, value in row.items() if key != "method"},
                "method": "B1c",
                "_support": probabilities,
                "_p": probabilities,
            }
        )
    frame = pd.concat([frame, pd.DataFrame(b1c_rows)], ignore_index=True)
    output = pd.DataFrame([_project(row) for row in frame.to_dict("records")])
    output["m_layer"] = np.where(output.m == 1, 1, 2)
    return output.sort_values(
        ["method", "fold", "entry", "held_out_unit"]
    ).reset_index(drop=True)


def document_fold_predictions(
    records: pd.DataFrame, target: str, endpoint: str
) -> pd.DataFrame:
    """Build document-grouped folds; all fitted statistics exclude the test fold."""
    doc = document_values(records, merge_identical=True)
    doc["document_fold"] = document_fold_map(doc)
    units_by_entry = dict(tuple(doc.groupby("entry", sort=True)))
    source_by_entry = dict(
        tuple(records.groupby("molecule_identity_key", sort=True))
    )
    rows = []
    for document_fold in range(10):
        train = doc.loc[doc.document_fold != document_fold]
        test = doc.loc[doc.document_fold == document_fold]
        sigma_fit = pooled_sigma(train)
        sigma = float(sigma_fit["tau"])
        prior = empirical_bayes_prior(train, sigma)
        sigmas = heterogeneity_sigmas(train, sigma)

        # B1c uses only leave-one-unit-out G predictions constructed inside train.
        b1c_statistics = []
        for entry, group in train.groupby("entry", sort=True):
            if len(group) < 2:
                continue
            source = source_by_entry[entry]
            train_documents = set(
                ";".join(group.documents.astype(str)).split(";")
            )
            train_source = source.loc[
                source.document_chembl_id.astype(str).isin(train_documents)
            ]
            for held in group.itertuples(index=False):
                held_documents = set(str(held.documents).split(";"))
                rest = group.loc[
                    ~group.documents.map(
                        lambda value: bool(
                            set(str(value).split(";")) & held_documents
                        )
                    )
                ]
                raw, _ = _remaining_raw(train_source, held_documents)
                if rest.empty or not len(raw):
                    continue
                b1c_statistics.append(
                    {
                        "displayed": classify(np.median(raw)),
                        "observed": classify(float(held.document_value)),
                        "m_layer": 1 if len(rest) == 1 else 2,
                    }
                )
        b1c_statistics = pd.DataFrame(b1c_statistics)

        fold_rows = []
        for held in test.itertuples(index=False):
            entry = str(held.entry)
            group = units_by_entry[entry]
            held_documents = set(str(held.documents).split(";"))
            remaining_units = group.loc[
                ~group.documents.map(
                    lambda value: bool(
                        set(str(value).split(";")) & held_documents
                    )
                )
            ].loc[lambda frame: frame.document_fold != document_fold]
            source = source_by_entry[entry]
            remaining_documents = set(
                ";".join(remaining_units.documents.astype(str)).split(";")
            )
            train_source = source.loc[
                source.document_chembl_id.astype(str).isin(
                    remaining_documents
                )
            ]
            remaining_raw, assay_types = _remaining_raw(
                train_source, held_documents
            )
            if remaining_units.empty or not len(remaining_raw):
                continue
            values = remaining_units.document_value.to_numpy(float)
            displayed = classify(np.median(remaining_raw))
            observed = classify(float(held.document_value))
            base = {
                "split_design": "document",
                "merge_basis": "same_value",
                "target": target,
                "endpoint": endpoint,
                "fold": int(document_fold),
                "entry": entry,
                "held_out_unit": str(held.documents),
                "m": int(len(values)),
                "displayed": displayed,
                "observed": observed,
                "held_out_value": float(held.document_value),
                "sigma_doc": sigma,
                "sigma_method": sigma_fit["method"],
                "prior_mu0": float(prior["mu0"]),
                "prior_tau2": float(prior["tau2"]),
                "prior_tau": float(prior["tau"]),
                "remaining_center": float(values.mean()),
                "remaining_values": ";".join(
                    f"{value:.10g}" for value in values
                ),
                "remaining_assay_types": ";".join(assay_types),
            }
            g_predictive = gaussian_prediction(values, sigma_fit)
            g_support = class_probabilities(
                [values.mean()],
                [sigma / np.sqrt(len(values))],
                np.inf,
                *CLASS_BOUNDARIES,
            )[0]
            fold_rows.append(
                {
                    **base,
                    "method": "G",
                    "_support": g_support,
                    "_p": g_predictive,
                }
            )
            eb_support, eb_predictive, eb_mean, eb_var = eb_prediction(
                values, sigma, prior
            )
            fold_rows.append(
                {
                    **base,
                    "method": "G-EB",
                    "posterior_mean": eb_mean,
                    "posterior_variance": eb_var,
                    "_support": eb_support,
                    "_p": eb_predictive,
                }
            )
            hs_support, hs_predictive, selected_sigma = hs_prediction(
                values, sigmas
            )
            fold_rows.append(
                {
                    **base,
                    "method": "G-HS",
                    "sigma_selected": selected_sigma,
                    "_support": hs_support,
                    "_p": hs_predictive,
                }
            )
            ebhs_support, ebhs_predictive, ebhs_mean, ebhs_var, ebhs_sigma = (
                eb_hs_prediction(values, sigmas, prior)
            )
            fold_rows.append(
                {
                    **base,
                    "method": "G-EB-HS",
                    "posterior_mean": ebhs_mean,
                    "posterior_variance": ebhs_var,
                    "sigma_selected": ebhs_sigma,
                    "_support": ebhs_support,
                    "_p": ebhs_predictive,
                }
            )
            b1c = b1c_probabilities(b1c_statistics, displayed, len(values))
            fold_rows.append(
                {
                    **base,
                    "method": "B1c",
                    "_support": b1c,
                    "_p": b1c,
                }
            )
        rows.extend(fold_rows)
    output = pd.DataFrame([_project(row) for row in rows])
    output["m_layer"] = np.where(output.m == 1, 1, 2)
    return output.sort_values(
        ["method", "fold", "entry", "held_out_unit"]
    ).reset_index(drop=True)


def near_duplicate_predictions(
    records: pd.DataFrame, target: str, endpoint: str
) -> pd.DataFrame:
    """Build 0.1-merged predictions with fold-out sigma and EB prior."""
    original = document_values(records, merge_identical=True)
    assignments = entry_folds(original, seed=FOLD_SEED)
    units = near_duplicate_units(records)
    units_by_entry = dict(tuple(units.groupby("entry", sort=True)))
    source_by_entry = dict(
        tuple(records.groupby("molecule_identity_key", sort=True))
    )
    rows = []
    for fold in range(10):
        test_entries = {
            entry for entry, assigned in assignments.items() if assigned == fold
        }
        train = units.loc[~units.entry.isin(test_entries)]
        sigma_fit = pooled_sigma(train)
        sigma = float(sigma_fit["tau"])
        prior = empirical_bayes_prior(train, sigma)
        for entry in sorted(test_entries):
            group = units_by_entry[entry]
            source = source_by_entry[entry]
            for held_index, held in group.iterrows():
                held_documents = set(str(held.documents).split(";"))
                remaining = group.drop(index=held_index)
                remaining_raw, assay_types = _remaining_raw(
                    source, held_documents
                )
                if remaining.empty or not len(remaining_raw):
                    continue
                values = remaining.document_value.to_numpy(float)
                displayed = classify(np.median(remaining_raw))
                observed = classify(float(held.document_value))
                base = {
                    "split_design": "entry",
                    "merge_basis": "near_duplicate_0.1",
                    "target": target,
                    "endpoint": endpoint,
                    "fold": int(fold),
                    "entry": entry,
                    "held_out_unit": str(held.documents),
                    "m": int(len(values)),
                    "displayed": displayed,
                    "observed": observed,
                    "held_out_value": float(held.document_value),
                    "sigma_doc": sigma,
                    "sigma_method": sigma_fit["method"],
                    "prior_mu0": float(prior["mu0"]),
                    "prior_tau2": float(prior["tau2"]),
                    "prior_tau": float(prior["tau"]),
                    "remaining_center": float(values.mean()),
                    "remaining_values": ";".join(
                        f"{value:.10g}" for value in values
                    ),
                    "remaining_assay_types": ";".join(assay_types),
                }
                g_predictive = gaussian_prediction(values, sigma_fit)
                g_support = class_probabilities(
                    [values.mean()],
                    [sigma / np.sqrt(len(values))],
                    np.inf,
                    *CLASS_BOUNDARIES,
                )[0]
                rows.append(
                    {
                        **base,
                        "method": "G",
                        "_support": g_support,
                        "_p": g_predictive,
                    }
                )
                eb_support, eb_predictive, eb_mean, eb_var = eb_prediction(
                    values, sigma, prior
                )
                rows.append(
                    {
                        **base,
                        "method": "G-EB",
                        "posterior_mean": eb_mean,
                        "posterior_variance": eb_var,
                        "_support": eb_support,
                        "_p": eb_predictive,
                    }
                )
    output = pd.DataFrame([_project(row) for row in rows])
    output["m_layer"] = np.where(output.m == 1, 1, 2)
    return output.sort_values(
        ["method", "fold", "entry", "held_out_unit"]
    ).reset_index(drop=True)
