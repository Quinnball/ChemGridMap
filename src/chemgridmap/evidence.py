"""Evidence-aware activity annotations.

A median over source records lets one publication with many condition
variants dominate an entry, and it hides how often independent publications
disagree about the displayed activity class. This module treats each source
document as an evidence unit (documents that report an identical value are
merged, because later publications often re-report an earlier measurement)
and models unit values within an entry as

    y_eu ~ N(mu_e, sigma^2),

where sigma is the pooled standard deviation between evidence units of the
same entry in the loaded dataset. Each entry receives a document-balanced
estimate (mean over units), the probability of each display class for its
consensus value mu_e (flat prior), and the predictive probability that a new
independent document falls in each class. The default Gaussian model uses
pooled between-unit variance. An optional ``exclude-flagged`` model removes
evidence units containing documents involved in systematic document-offset
pairs; entries without remaining units fall back to the Gaussian model.
Moderated-variance t models remain available to library callers for comparison.
"""
from __future__ import annotations

from collections import defaultdict
from itertools import combinations
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import optimize, stats

from .activity import activity_annotations, activity_class
from .core import canonicalize_smiles

CLASSES = ("inactive", "medium", "active")
VARIANCE_FLOOR = 1e-4  # pChEMBL values are reported to two decimals
# Fallback for datasets with too few repeated entries: Landrum & Riniker
# (J. Chem. Inf. Model. 2024, 64, 1560) report that 27% of IC50 pairs from
# different assays differ by more than one log unit; for Gaussian errors this
# corresponds to a per-measurement SD of 1 / (sqrt(2) * 1.103) = 0.64.
LITERATURE_BETWEEN_UNIT_SD = 0.64
DOCUMENT_OFFSET_COLUMNS = (
    "document_a", "document_b", "shared_entries", "median_offset", "mean_offset",
    "std_offset", "sign_agreement", "systematic_offset",
)


def document_values(records: pd.DataFrame, *, entry: str = "molecule_identity_key",
                    value: str = "pchembl_value", document: str = "document_chembl_id",
                    merge_identical: bool = True, decimals: int = 2) -> pd.DataFrame:
    """Collapse source records to evidence units: one median value per entry and document.

    With ``merge_identical``, documents reporting the same entry value (at the
    reported precision) form one unit.
    """
    missing = [name for name in (entry, value, document) if name not in records]
    if missing:
        raise ValueError(f"Records lack required columns: {missing}")
    frame = records[[entry, value, document]].copy()
    frame[value] = pd.to_numeric(frame[value], errors="coerce")
    frame = frame[np.isfinite(frame[value])]
    frame[document] = frame[document].fillna("unknown_document").astype(str)
    grouped = frame.groupby([entry, document], sort=True)[value]
    doc = pd.DataFrame({"document_value": grouped.median(), "n_records": grouped.size()}).reset_index().rename(
        columns={entry: "entry", document: "document"})
    doc["documents"] = doc["document"]
    if merge_identical:
        doc["rounded"] = doc["document_value"].round(decimals)
        merged = doc.groupby(["entry", "rounded"], sort=True)
        doc = pd.DataFrame({"document": merged["document"].first(), "document_value": merged["document_value"].median(),
                            "n_records": merged["n_records"].sum(),
                            "documents": merged["document"].agg(";".join)}).reset_index().drop(columns="rounded")
    return doc[["entry", "document", "document_value", "n_records", "documents"]]


def document_pair_offsets(records: pd.DataFrame, *, entry: str = "molecule_identity_key",
                          value: str = "pchembl_value", document: str = "document_chembl_id",
                          min_shared: int = 3, offset_threshold: float = 1.0,
                          sign_agreement: float = 0.8) -> pd.DataFrame:
    """Diagnose document offsets from the median value of each shared entry.

    Documents A/B are sorted by ID; offsets are B minus A. The standard
    deviation is the sample SD (zero for a single shared entry). Agreement is
    the fraction with the median offset's sign; zero differences do not agree.
    Only pairs with at least ``min_shared`` entries are returned. Counts of
    all possible, observed (at least one shared entry), and evaluated pairs
    are retained in ``attrs``.
    """
    if min_shared < 1:
        raise ValueError("min_shared must be at least 1.")
    if offset_threshold < 0 or not 0 <= sign_agreement <= 1:
        raise ValueError("offset_threshold must be nonnegative and sign_agreement must lie between 0 and 1.")
    doc = document_values(records, entry=entry, value=value, document=document, merge_identical=False)
    differences = defaultdict(list)
    for _, group in doc.groupby("entry", sort=False):
        values = group[["document", "document_value"]].itertuples(index=False, name=None)
        for (document_a, a), (document_b, b) in combinations(values, 2):
            differences[(document_a, document_b)].append(float(b - a))
    rows = []
    for (document_a, document_b), values in sorted(differences.items()):
        if len(values) < min_shared:
            continue
        values = np.asarray(values, dtype=float)
        median = float(np.median(values))
        agreement = float(np.mean(values * np.sign(median) > 0))
        rows.append({
            "document_a": document_a, "document_b": document_b, "shared_entries": len(values),
            "median_offset": median, "mean_offset": float(values.mean()),
            "std_offset": float(values.std(ddof=1)) if len(values) > 1 else 0.0,
            "sign_agreement": agreement,
            "systematic_offset": bool(abs(median) >= offset_threshold and agreement >= sign_agreement),
        })
    out = pd.DataFrame(rows, columns=DOCUMENT_OFFSET_COLUMNS)
    out = out.sort_values(["shared_entries", "document_a", "document_b"], ascending=[False, True, True]).reset_index(drop=True)
    n_documents = int(doc.document.nunique())
    out.attrs.update({"document_pair_possible_count": n_documents * (n_documents - 1) // 2,
                      "observed_document_pair_count": len(differences), "document_pair_count": len(out),
                      "min_shared": min_shared, "offset_threshold": offset_threshold,
                      "sign_agreement": sign_agreement})
    return out


def _entry_moments(doc: pd.DataFrame) -> pd.DataFrame:
    grouped = doc.groupby("entry")["document_value"]
    moments = pd.DataFrame({"n_units": grouped.size(), "mean": grouped.mean(),
                            "median": grouped.median(), "variance": grouped.var(ddof=1)})
    moments["variance"] = moments["variance"].fillna(0.0)
    return moments


def fit_between_unit_sd(doc: pd.DataFrame, *, min_entries: int = 5,
                        fallback_sd: Optional[float] = LITERATURE_BETWEEN_UNIT_SD) -> Dict[str, float]:
    """Pooled SD between evidence units of the same entry (Gaussian model)."""
    moments = _entry_moments(doc)
    multi = moments[moments.n_units >= 2]
    base = {"prior_df": float("inf"), "model": "gaussian", "n_entries": int(len(multi)),
            "n_units": int(multi.n_units.sum())}
    if len(multi) < min_entries:
        if fallback_sd is None:
            raise ValueError(f"At least {min_entries} entries with two or more evidence units are needed; "
                             f"found {len(multi)}.")
        return {**base, "tau": float(fallback_sd), "method": "literature_fallback"}
    dfs = multi.n_units.to_numpy(float) - 1.0
    pooled = float(np.sqrt(np.sum(multi.variance.to_numpy(float) * dfs) / np.sum(dfs)))
    return {**base, "tau": pooled, "method": "pooled_between_unit_sd", "degrees_of_freedom": int(dfs.sum())}


def fit_document_variance(doc: pd.DataFrame, *, min_entries: int = 5) -> Dict[str, float]:
    """Moderated t model: (tau, prior_df) by maximum likelihood over within-entry variances."""
    moments = _entry_moments(doc)
    multi = moments[moments.n_units >= 2]
    if len(multi) < min_entries:
        raise ValueError(f"At least {min_entries} entries with two or more evidence units are needed; found {len(multi)}.")
    variances = np.maximum(multi["variance"].to_numpy(float), VARIANCE_FLOOR)
    dfs = multi["n_units"].to_numpy(float) - 1.0

    def negative_log_likelihood(params):
        tau2, prior_df = np.exp(params)
        return -np.sum(stats.f.logpdf(variances / tau2, dfs, prior_df) - np.log(tau2))

    pooled = float(np.sum(variances * dfs) / np.sum(dfs))
    result = optimize.minimize(negative_log_likelihood, np.log([pooled, 4.0]), method="L-BFGS-B",
                               bounds=[(np.log(1e-3), np.log(25.0)), (np.log(0.5), np.log(1e3))])
    tau2, prior_df = np.exp(result.x)
    return {"tau": float(np.sqrt(tau2)), "prior_df": float(prior_df), "pooled_sd": float(np.sqrt(pooled)),
            "n_entries": int(len(multi)), "n_units": int(multi.n_units.sum()), "model": "moderated_t",
            "converged": bool(result.success), "method": "scaled_F_maximum_likelihood"}


def _held_out_units(doc: pd.DataFrame) -> Dict[str, np.ndarray]:
    """Leave-one-unit-out statistics for every entry with two or more evidence units."""
    observed, center, own, remaining = [], [], [], []
    for _, values in doc.groupby("entry")["document_value"]:
        values = values.to_numpy(float)
        if len(values) < 2:
            continue
        for k in range(len(values)):
            rest = np.delete(values, k)
            observed.append(values[k])
            center.append(rest.mean())
            own.append(max(rest.var(ddof=1), VARIANCE_FLOOR) if len(rest) > 1 else 0.0)
            remaining.append(len(rest))
    return {"observed": np.array(observed), "center": np.array(center), "own": np.array(own),
            "remaining": np.array(remaining, dtype=float)}


def predictive_parameters(center, own, remaining, tau: float, prior_df: float):
    """Location, scale and degrees of freedom of the prediction for a new unit.

    ``prior_df=inf`` gives the Gaussian model, in which the entry's own spread
    is not used; finite values give the moderated t model.
    """
    center, own, remaining = (np.asarray(v, dtype=float) for v in (center, own, remaining))
    if np.isinf(prior_df):
        return center, tau * np.sqrt(1.0 + 1.0 / remaining), np.full_like(remaining, np.inf)
    df = prior_df + remaining - 1.0
    s2 = (prior_df * tau ** 2 + (remaining - 1.0) * own) / df
    return center, np.sqrt(s2 * (1.0 + 1.0 / remaining)), df


def fit_predictive_variance(doc: pd.DataFrame, *, min_entries: int = 5) -> Dict[str, float]:
    """Moderated t model: (tau, prior_df) maximizing leave-one-unit-out predictive likelihood."""
    held = _held_out_units(doc)
    n_entries = int((doc.groupby("entry").size() >= 2).sum())
    if n_entries < min_entries:
        raise ValueError(f"At least {min_entries} entries with two or more evidence units are needed; found {n_entries}.")

    def negative_log_likelihood(params):
        tau, prior_df = np.exp(params)
        loc, scale, df = predictive_parameters(held["center"], held["own"], held["remaining"], tau, prior_df)
        return -np.sum(stats.t.logpdf(held["observed"], df, loc=loc, scale=scale))

    result = optimize.minimize(negative_log_likelihood, np.log([0.7, 4.0]), method="L-BFGS-B",
                               bounds=[(np.log(0.05), np.log(3.0)), (np.log(0.5), np.log(200.0))])
    tau, prior_df = np.exp(result.x)
    return {"tau": float(tau), "prior_df": float(prior_df), "n_entries": n_entries,
            "n_predictions": int(len(held["observed"])), "converged": bool(result.success),
            "model": "moderated_t", "method": "leave_one_unit_out_predictive_likelihood"}


def class_probabilities(center, scale, df, lower: float, upper: float) -> np.ndarray:
    """Probabilities of the inactive (<=lower), medium and active (>=upper) classes."""
    center, scale, df = (np.broadcast_to(np.asarray(v, dtype=float), np.shape(center)) for v in (center, scale, df))
    gaussian = np.isinf(df)
    finite_df = np.where(gaussian, 1.0, df)
    safe_scale = np.where(scale == 0, 1.0, scale)
    z_low, z_high = (lower - center) / safe_scale, (upper - center) / safe_scale
    low = np.where(gaussian, stats.norm.cdf(z_low), stats.t.cdf(z_low, finite_df))
    high = np.where(gaussian, stats.norm.sf(z_high), stats.t.sf(z_high, finite_df))
    probabilities = np.column_stack([low, np.clip(1.0 - low - high, 0.0, 1.0), high])
    # Zero between-document variance is a point mass, including at a class boundary.
    point_mass = np.column_stack([center <= lower, (center > lower) & (center < upper), center >= upper]).astype(float)
    return np.where(np.atleast_1d(scale == 0)[:, None], point_mass, probabilities)


def entry_posteriors(doc: pd.DataFrame, tau: float, prior_df: float) -> pd.DataFrame:
    """Document-balanced estimate with consensus and predictive scales per entry."""
    moments = _entry_moments(doc)
    m = moments["n_units"].to_numpy(float)
    if np.isinf(prior_df):
        s2, df = np.full_like(m, tau ** 2), np.full_like(m, np.inf)
    else:
        own = np.where(m > 1, np.maximum(moments["variance"].to_numpy(float), VARIANCE_FLOOR), 0.0)
        df = prior_df + m - 1.0
        s2 = (prior_df * tau ** 2 + (m - 1.0) * own) / df
    moments["posterior_df"] = df
    moments["posterior_scale"] = np.sqrt(s2 / m)
    moments["predictive_scale"] = np.sqrt(s2 * (1.0 + 1.0 / m))
    return moments


def annotate_activity_confidence(molecules: pd.DataFrame, records: pd.DataFrame, *, lower: float = 6.0,
                                 upper: float = 7.0, confidence_threshold: float = 0.8,
                                 model: str = "gaussian", fit: Optional[Dict[str, object]] = None,
                                 entry: str = "molecule_identity_key", value: str = "pchembl_value",
                                 document: str = "document_chembl_id", displayed: str = "activity_pchembl",
                                 merge_identical: bool = True):
    """Add document-balanced estimates and class confidence to entry-level annotations.

    Returns the annotated copy and the fit. Displayed classes, medians and
    source records are unchanged; the added columns describe how strongly the
    evidence units support the displayed class.
    """
    if not 0.5 < confidence_threshold < 1.0:
        raise ValueError("confidence_threshold must lie between 0.5 and 1.")
    if model not in {"gaussian", "exclude-flagged", "moderated_t"}:
        raise ValueError("model must be 'gaussian', 'exclude-flagged' or 'moderated_t'.")
    doc = document_values(records, entry=entry, value=value, document=document, merge_identical=merge_identical)
    offsets = document_pair_offsets(records, entry=entry, value=value, document=document)
    flagged = offsets[offsets.systematic_offset.astype(bool)].to_dict("records")
    if fit is None:
        fit = fit_predictive_variance(doc) if model == "moderated_t" else fit_between_unit_sd(doc)
    fit = dict(fit)
    if fit.get("prior_df") is None:  # JSON reports store the Gaussian model's infinite prior_df as null
        fit["prior_df"] = float("inf")
    gaussian_fit = dict(fit)

    def selected_posteriors(selected):
        posterior = entry_posteriors(selected, fit["tau"], fit["prior_df"])
        top = selected.sort_values(
            ["entry", "n_records"], ascending=[True, False]
        ).groupby("entry").first()
        posterior["top_unit_record_fraction"] = (
            top["n_records"] / selected.groupby("entry")["n_records"].sum()
        )
        posterior["top_unit_documents"] = top["documents"]
        return posterior

    if model == "exclude-flagged":
        flagged_documents = {
            str(pair[name])
            for pair in flagged
            for name in ("document_a", "document_b")
        }
        keep = doc["documents"].map(
            lambda text: not (set(str(text).split(";")) & flagged_documents)
        )
        remaining = doc.loc[keep]
        excluded_entries = sorted(set(doc["entry"]) - set(remaining["entry"]))
        fallback = doc[doc["entry"].isin(excluded_entries)]
        fit = fit_between_unit_sd(remaining) if not remaining.empty else fit_between_unit_sd(doc)
        if fit.get("prior_df") is None:
            fit["prior_df"] = float("inf")
        pieces = []
        if not remaining.empty:
            pieces.append(selected_posteriors(remaining))
        if not fallback.empty:
            fallback_post = entry_posteriors(
                fallback, gaussian_fit["tau"], gaussian_fit["prior_df"]
            )
            top = fallback.sort_values(
                ["entry", "n_records"], ascending=[True, False]
            ).groupby("entry").first()
            fallback_post["top_unit_record_fraction"] = (
                top["n_records"] / fallback.groupby("entry")["n_records"].sum()
            )
            fallback_post["top_unit_documents"] = top["documents"]
            pieces.append(fallback_post)
        post = pd.concat(pieces).sort_index()
    else:
        post = selected_posteriors(doc)

    out = molecules.copy()
    keyed = post.reindex(out[entry].to_numpy())
    if keyed["n_units"].isna().any():
        raise ValueError("Every entry needs at least one record with a finite value.")
    consensus = class_probabilities(keyed["mean"], keyed["posterior_scale"], keyed["posterior_df"], lower, upper)
    predictive = class_probabilities(keyed["mean"], keyed["predictive_scale"], keyed["posterior_df"], lower, upper)
    shown = np.array([activity_class(v, lower, upper) for v in out[displayed]])
    index = np.array([CLASSES.index(c) for c in shown])
    rows = np.arange(len(out))
    out["n_documents"] = records.groupby(entry)[document].nunique().reindex(out[entry]).fillna(0).astype(int).to_numpy()
    out["n_evidence_units"] = keyed["n_units"].astype(int).to_numpy()
    out["top_unit_documents"] = keyed["top_unit_documents"].to_numpy()
    out["top_unit_record_fraction"] = keyed["top_unit_record_fraction"].to_numpy()
    out["document_balanced_pchembl"] = keyed["mean"].to_numpy()
    out["document_balanced_class"] = [activity_class(v, lower, upper) for v in keyed["mean"]]
    out["class_changed_by_document_balance"] = (out["document_balanced_class"].to_numpy() != shown).astype(int)
    out["consensus_scale"] = keyed["posterior_scale"].to_numpy()
    for k, name in enumerate(CLASSES):
        out[f"p_consensus_{name}"] = consensus[:, k]
        out[f"p_new_document_{name}"] = predictive[:, k]
    out["class_confidence"] = consensus[rows, index]
    out["class_reproducibility"] = predictive[rows, index]
    out["class_status"] = np.where(out["class_confidence"] >= confidence_threshold, "supported", "uncertain")
    doc_entries = defaultdict(set)
    for item in doc.itertuples(index=False):
        for source_document in item.documents.split(";"):
            doc_entries[source_document].add(item.entry)
    discordant_pairs, discordant_documents = defaultdict(int), defaultdict(set)
    for pair in flagged:
        document_a, document_b = pair["document_a"], pair["document_b"]
        for key in doc_entries[document_a] & doc_entries[document_b]:
            discordant_pairs[key] += 1
            discordant_documents[key].update((document_a, document_b))
    out["discordant_document_pairs"] = [discordant_pairs[key] for key in out[entry]]
    out["discordant_documents"] = [";".join(sorted(discordant_documents[key])) for key in out[entry]]
    policy = "identical values merged" if merge_identical else "identical values kept separate"
    fit.update({"model": "exclude_offset_documents" if model == "exclude-flagged" else model,
                "confidence_threshold": confidence_threshold, "lower_threshold": lower,
                "upper_threshold": upper, "evidence_unit": f"{document} ({policy})",
                "evidence_unit_policy": f"document; {policy}", "merge_identical_documents": bool(merge_identical),
                "document_pair_count": offsets.attrs["document_pair_count"],
                "observed_document_pair_count": offsets.attrs["observed_document_pair_count"],
                "document_pair_possible_count": offsets.attrs["document_pair_possible_count"],
                "discordant_document_pair_count": len(flagged), "document_offsets_top10": flagged[:10],
                "document_offsets": flagged, "document_offset_min_shared": offsets.attrs["min_shared"],
                "document_offset_threshold": offsets.attrs["offset_threshold"],
                "document_offset_sign_agreement": offsets.attrs["sign_agreement"]})
    if model == "exclude-flagged":
        fit["flagged_documents"] = sorted(flagged_documents)
        fit["excluded_entry_count"] = len(excluded_entries)
    fit = {key: None if isinstance(val, float) and not np.isfinite(val) else val for key, val in fit.items()}
    out.attrs["activity_confidence_fit"] = fit
    return out, fit


def prepare_generic_activity_records(
    data: pd.DataFrame,
    *,
    compound_column: str,
    value_column: str,
    document_column: Optional[str],
    smiles_column: str,
    relation_column: Optional[str] = None,
    lower: float = 6.0,
    upper: float = 7.0,
    conflict_range_threshold: float = 1.0,
    confidence: bool = True,
    merge_identical: bool = True,
    evidence_model: str = "gaussian",
) -> Tuple[pd.DataFrame, pd.DataFrame, Optional[Dict[str, object]], Dict[str, object]]:
    """Curate a generic activity-record CSV into the molecule-level map schema.

    Each compound ID is one entry. Values must already be pActivity values
    (-log10 M). If a document column is supplied, class support is calculated
    by default with the same evidence functions used by the ChEMBL path.
    """
    required = [compound_column, value_column, smiles_column]
    if confidence:
        if not document_column:
            raise ValueError("document_column is required when activity confidence is enabled.")
        required.append(document_column)
    missing = [name for name in required if name not in data.columns]
    if missing:
        raise ValueError(f"Generic activity input lacks required columns: {missing}.")
    if relation_column and relation_column not in data.columns:
        raise ValueError(f"Relation column {relation_column!r} was not found.")
    if lower >= upper:
        raise ValueError("lower must be smaller than upper.")

    records = data.copy()
    records.insert(0, "source_row", np.arange(len(records), dtype=int))
    numeric = pd.to_numeric(records[value_column], errors="coerce")
    canonical = records[smiles_column].map(canonicalize_smiles)
    compounds = records[compound_column].fillna("").astype(str).str.strip()
    keep = np.isfinite(numeric).to_numpy() & canonical.notna().to_numpy() & compounds.ne("").to_numpy()
    excluded_non_exact = 0
    if relation_column:
        relation = records[relation_column].fillna("").astype(str).str.strip()
        relation = relation.str.replace(r"^'([=<>~]+)'$", r"\1", regex=True)
        exact = relation.eq("=").to_numpy()
        excluded_non_exact = int((keep & ~exact).sum())
        keep &= exact
    records = records.loc[keep].copy()
    if records.empty:
        raise ValueError("No finite non-exact relation, SMILES and compound records remained.")
    records["molecule_identity_key"] = compounds.loc[keep].to_numpy()
    records["pchembl_value"] = numeric.loc[keep].to_numpy(float)
    records["canonical_smiles"] = canonical.loc[keep].to_numpy()
    if confidence:
        records["document_chembl_id"] = (
            records[document_column].fillna("unknown_document").astype(str).str.strip()
        )

    structures = records.groupby("molecule_identity_key", sort=True)["canonical_smiles"].nunique()
    conflicting = structures[structures > 1]
    if not conflicting.empty:
        raise ValueError(
            "Each generic compound ID must map to one canonical structure; "
            f"conflicts: {conflicting.index.tolist()}."
        )

    rows = []
    for identity, group in records.groupby("molecule_identity_key", sort=True):
        annotations = activity_annotations(
            group["pchembl_value"].to_numpy(float),
            lower_threshold=lower,
            upper_threshold=upper,
            conflict_range_threshold=conflict_range_threshold,
        )
        smiles = str(group["canonical_smiles"].iloc[0])
        rows.append({
            "molecule_id": identity,
            "molecule_identity_key": identity,
            "molecule_identity_source": "generic_compound_id",
            "canonical_smiles": smiles,
            "smiles": smiles,
            "record_canonical_smiles": smiles,
            "n_record_structures": 1,
            "source_rows": ";".join(str(value) for value in group["source_row"].unique()),
            **annotations,
        })
    molecules = pd.DataFrame(rows).sort_values("molecule_identity_key").reset_index(drop=True)
    duplicate_structures = int(molecules["canonical_smiles"].duplicated(keep=False).sum())
    if duplicate_structures:
        raise ValueError(
            "Different generic compound IDs resolve to the same canonical structure. "
            "Merge them before mapping."
        )

    fit = None
    if confidence:
        molecules, fit = annotate_activity_confidence(
            molecules,
            records,
            lower=lower,
            upper=upper,
            model=evidence_model,
            merge_identical=merge_identical,
        )
    report = {
        "input_rows": int(len(data)),
        "retained_activity_records": int(len(records)),
        "excluded_activity_records": int(len(data) - len(records)),
        "excluded_non_exact_relation_records": excluded_non_exact,
        "unique_molecules": int(len(molecules)),
        "compound_column": compound_column,
        "value_column": value_column,
        "document_column": document_column,
        "smiles_column": smiles_column,
        "relation_column": relation_column,
        "aggregation": "median pActivity",
        "activity_confidence": fit,
    }
    return molecules, records, fit, report
