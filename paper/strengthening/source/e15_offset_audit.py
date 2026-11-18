"""E15 descriptive audit of systematic document offsets."""
from __future__ import annotations

from collections import Counter, defaultdict
import gzip
from pathlib import Path
import pickle
import sys

sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parents[3]
S = ROOT / "paper/strengthening"
SOURCE = S / "source"
sys.path[:0] = [str(SOURCE), str(ROOT / "paper"), str(ROOT / "src")]

import numpy as np
import pandas as pd
from scipy import stats

from chemgridmap.evidence import document_pair_offsets, document_values
from e2_models import curate_activity_table

S1_TOLERANCE = 0.05
MIN_S1_SHARED = 3
S1_AGREEMENT = 0.8
S2_TARGETS = (3.0, 6.0, 9.0)
S2_MEDIAN_TOLERANCE = 0.1
S2_FRACTION = 0.8
S2_HALF_WIDTH = 0.2
S3_IQR_LIMIT = 0.1
FEATURE_ORDER = ("S1", "S2", "S3", "S4", "S5", "S6", "S0")


def split_documents(value) -> set[str]:
    return {part for part in str(value).split(";") if part}


def pair_key(document_a: str, document_b: str) -> tuple[str, str]:
    return tuple(sorted((str(document_a), str(document_b))))


def exact_interval(successes: int, n: int) -> tuple[float, float]:
    if not n:
        return np.nan, np.nan
    interval = stats.binomtest(
        int(successes), int(n)
    ).proportion_ci(confidence_level=0.95, method="exact")
    return float(interval.low), float(interval.high)


def share_row(feature: str, successes: int, n: int) -> dict:
    low, high = exact_interval(successes, n)
    return {
        "feature": feature,
        "n": int(n),
        "successes": int(successes),
        "proportion": float(successes / n) if n else np.nan,
        "ci_low_exact": low,
        "ci_high_exact": high,
    }


def exact_values(records: pd.DataFrame) -> pd.DataFrame:
    required = [
        "document_chembl_id",
        "parent_molecule_chembl_id",
        "pchembl_value",
        "standard_relation",
    ]
    missing = [name for name in required if name not in records]
    if missing:
        raise ValueError(f"Records lack S1 columns: {missing}")
    frame = records.loc[
        records.standard_relation.astype(str) == "=",
        required,
    ].copy()
    frame["pchembl_value"] = pd.to_numeric(
        frame.pchembl_value, errors="coerce"
    )
    frame = frame.loc[
        np.isfinite(frame.pchembl_value)
        & frame.document_chembl_id.notna()
        & frame.parent_molecule_chembl_id.notna()
    ]
    return (
        frame.groupby(
            ["document_chembl_id", "parent_molecule_chembl_id"],
            sort=True,
        )
        .pchembl_value.median()
        .reset_index()
    )


def dominant_assay_types(records: pd.DataFrame) -> dict[tuple[str, str], str]:
    counts = defaultdict(Counter)
    for item in records[
        ["document_chembl_id", "molecule_identity_key", "assay_type"]
    ].itertuples(index=False):
        counts[(str(item.document_chembl_id), str(item.molecule_identity_key))][
            str(item.assay_type)
        ] += 1
    return {
        key: value.most_common(1)[0][0]
        for key, value in counts.items()
    }


def source_ids(records: pd.DataFrame) -> dict[str, str]:
    result = {}
    for document, group in records.groupby(
        "document_chembl_id", sort=True
    ):
        values = group.src_id.dropna().astype(str)
        result[str(document)] = values.iloc[0] if len(values) else ""
    return result


def pair_differences(
    records: pd.DataFrame,
    document_a: str,
    document_b: str,
) -> tuple[np.ndarray, list[str]]:
    doc = document_values(
        records,
        merge_identical=False,
    )
    pivot = doc.pivot(
        index="entry",
        columns="document",
        values="document_value",
    )
    shared = pivot.loc[
        pivot[document_a].notna() & pivot[document_b].notna()
    ]
    differences = (
        shared[document_b].to_numpy(float)
        - shared[document_a].to_numpy(float)
    )
    return differences, shared.index.astype(str).tolist()


def s2_status(differences: np.ndarray, median: float) -> bool:
    target_close = any(
        abs(abs(median) - target) <= S2_MEDIAN_TOLERANCE
        for target in S2_TARGETS
    )
    if not target_close:
        return False
    within = np.mean(np.abs(differences - median) <= S2_HALF_WIDTH)
    return bool(within >= S2_FRACTION)


def s3_status(differences: np.ndarray, median: float) -> bool:
    q1, q3 = np.percentile(differences, [25, 75])
    return bool(
        (q3 - q1) <= S3_IQR_LIMIT + 1e-12 and abs(median) >= 1.0
    )


class SwapAuditor:
    def __init__(
        self,
        records_by_key: dict[tuple[str, str], pd.DataFrame],
    ) -> None:
        self.values = {}
        self.parents = {}
        self.document_targets = defaultdict(lambda: defaultdict(set))
        for (target, endpoint), records in records_by_key.items():
            exact = exact_values(records)
            self.values[(target, endpoint)] = exact
            self.parents[(target, endpoint)] = exact.groupby(
                "document_chembl_id"
            ).parent_molecule_chembl_id.agg(set).to_dict()
            for document in exact.document_chembl_id.unique():
                self.document_targets[endpoint][str(document)].add(target)
        self._cache = {}

    def evaluate(
        self,
        document_a: str,
        document_b: str,
        target: str,
        endpoint: str,
    ) -> dict:
        key = (document_a, document_b, target, endpoint)
        if key in self._cache:
            return self._cache[key]
        candidates = (
            self.document_targets[endpoint][document_a]
            & self.document_targets[endpoint][document_b]
        ) - {target}
        best = None
        for target_b in sorted(candidates):
            base = self.values[(target, endpoint)]
            other = self.values[(target_b, endpoint)]
            base_a = base.loc[
                base.document_chembl_id == document_a
            ].set_index("parent_molecule_chembl_id").pchembl_value
            base_b = base.loc[
                base.document_chembl_id == document_b
            ].set_index("parent_molecule_chembl_id").pchembl_value
            other_a = other.loc[
                other.document_chembl_id == document_a
            ].set_index("parent_molecule_chembl_id").pchembl_value
            other_b = other.loc[
                other.document_chembl_id == document_b
            ].set_index("parent_molecule_chembl_id").pchembl_value
            parents = (
                base_a.index.intersection(base_b.index)
                .intersection(other_a.index)
                .intersection(other_b.index)
            )
            if len(parents) < MIN_S1_SHARED:
                continue
            swapped = (
                (
                    base_a.loc[parents] - other_b.loc[parents]
                ).abs()
                <= S1_TOLERANCE
            ) & (
                (
                    base_b.loc[parents] - other_a.loc[parents]
                ).abs()
                <= S1_TOLERANCE
            )
            count = int(swapped.sum())
            fraction = float(count / len(parents))
            current = {
                "is_swap": bool(fraction >= S1_AGREEMENT),
                "target_b": target_b,
                "shared_parents": int(len(parents)),
                "swapped_parents": count,
                "swap_fraction": fraction,
            }
            if best is None or (
                current["is_swap"],
                current["swapped_parents"],
                current["swap_fraction"],
            ) > (
                best["is_swap"],
                best["swapped_parents"],
                best["swap_fraction"],
            ):
                best = current
        if best is None:
            best = {
                "is_swap": False,
                "target_b": "",
                "shared_parents": 0,
                "swapped_parents": 0,
                "swap_fraction": np.nan,
            }
        self._cache[key] = best
        return best


def classify_feature(
    records: pd.DataFrame,
    row: pd.Series,
    swap_auditor: SwapAuditor,
    differences: np.ndarray,
    shared_entries: list[str],
    target: str,
    endpoint: str,
) -> dict:
    document_a = str(row.document_a)
    document_b = str(row.document_b)
    median = float(row.median_offset)
    swap = swap_auditor.evaluate(
        document_a,
        document_b,
        target,
        endpoint,
    )
    if swap["is_swap"]:
        return {
            "feature": "S1",
            "feature_label": "靶点互换",
            "evidence": (
                f"{swap['target_b']}: "
                f"{swap['swapped_parents']}/{swap['shared_parents']} "
                f"({swap['swap_fraction']:.1%}) simultaneously swapped"
            ),
            "s1_target": swap["target_b"],
            "s1_shared_parents": swap["shared_parents"],
            "s1_swapped_parents": swap["swapped_parents"],
            "s1_swap_fraction": swap["swap_fraction"],
        }
    if s2_status(differences, median):
        return {
            "feature": "S2",
            "feature_label": "单位换算",
            "evidence": (
                f"median={median:.3f}; "
                f"{np.mean(np.abs(differences - median) <= S2_HALF_WIDTH):.1%} "
                f"within median +/- {S2_HALF_WIDTH:g}"
            ),
            "s1_target": swap["target_b"],
            "s1_shared_parents": swap["shared_parents"],
            "s1_swapped_parents": swap["swapped_parents"],
            "s1_swap_fraction": swap["swap_fraction"],
        }
    q1, q3 = np.percentile(differences, [25, 75])
    if s3_status(differences, median):
        return {
            "feature": "S3",
            "feature_label": "恒定平移",
            "evidence": f"median={median:.3f}; IQR={q3 - q1:.3f}",
            "s1_target": swap["target_b"],
            "s1_shared_parents": swap["shared_parents"],
            "s1_swapped_parents": swap["swapped_parents"],
            "s1_swap_fraction": swap["swap_fraction"],
        }

    assay = dominant_assay_types(records)
    shared = set(shared_entries)
    assay_differences = [
        assay.get((document_a, entry), "")
        != assay.get((document_b, entry), "")
        for entry in shared
    ]
    if assay_differences and np.mean(assay_differences) > 0.5:
        return {
            "feature": "S4",
            "feature_label": "测定类型不同",
            "evidence": (
                f"{sum(assay_differences)}/{len(assay_differences)} "
                "shared entries have different dominant assay_type"
            ),
            "s1_target": swap["target_b"],
            "s1_shared_parents": swap["shared_parents"],
            "s1_swapped_parents": swap["swapped_parents"],
            "s1_swap_fraction": swap["swap_fraction"],
        }

    pair_records = records.loc[
        records.document_chembl_id.astype(str).isin((document_a, document_b))
        & records.molecule_identity_key.astype(str).isin(shared)
    ]
    validity = (
        pd.to_numeric(
            pair_records.potential_duplicate, errors="coerce"
        ).fillna(0)
        == 1
    ) | pair_records.data_validity_comment.fillna("").astype(str).ne("")
    if validity.any():
        count = int(validity.sum())
        return {
            "feature": "S5",
            "feature_label": "ChEMBL 已有标记",
            "evidence": f"{count}/{len(pair_records)} shared-entry records flagged",
            "s1_target": swap["target_b"],
            "s1_shared_parents": swap["shared_parents"],
            "s1_swapped_parents": swap["swapped_parents"],
            "s1_swap_fraction": swap["swap_fraction"],
        }

    sources = source_ids(records)
    if sources.get(document_a, "") != sources.get(document_b, ""):
        return {
            "feature": "S6",
            "feature_label": "来源不同",
            "evidence": (
                f"src_id {sources.get(document_a, '')} vs "
                f"{sources.get(document_b, '')}"
            ),
            "s1_target": swap["target_b"],
            "s1_shared_parents": swap["shared_parents"],
            "s1_swapped_parents": swap["swapped_parents"],
            "s1_swap_fraction": swap["swap_fraction"],
        }
    return {
        "feature": "S0",
        "feature_label": "未解释",
        "evidence": "",
        "s1_target": swap["target_b"],
        "s1_shared_parents": swap["shared_parents"],
        "s1_swapped_parents": swap["swapped_parents"],
        "s1_swap_fraction": swap["swap_fraction"],
    }


def touched_entries(records: pd.DataFrame, pair: pd.Series) -> set[str]:
    doc = document_values(records, merge_identical=False)
    selected = doc.loc[
        doc.document.isin((str(pair.document_a), str(pair.document_b)))
    ]
    return set(selected.entry.astype(str))


def audit_dataset(
    records: pd.DataFrame,
    *,
    dataset_kind: str,
    target: str,
    endpoint: str,
    swap_auditor: SwapAuditor,
) -> tuple[list[dict], list[dict]]:
    offsets = document_pair_offsets(records)
    flagged = offsets.loc[offsets.systematic_offset.astype(bool)]
    controls = offsets.loc[~offsets.systematic_offset.astype(bool)]
    feature_rows = []
    for row in flagged.itertuples(index=False):
        item = pd.Series(row._asdict())
        differences, shared = pair_differences(
            records, str(item.document_a), str(item.document_b)
        )
        touched = touched_entries(records, item)
        feature = classify_feature(
            records,
            item,
            swap_auditor,
            differences,
            shared,
            target,
            endpoint,
        )
        feature_rows.append(
            {
                "dataset_kind": dataset_kind,
                "target": target,
                "endpoint": endpoint,
                "document_a": str(item.document_a),
                "document_b": str(item.document_b),
                "shared_entries": int(item.shared_entries),
                "touched_entries": int(len(touched)),
                "touched_entry_ids": ";".join(sorted(touched)),
                "median_offset": float(item.median_offset),
                "mean_offset": float(item.mean_offset),
                "std_offset": float(item.std_offset),
                "sign_agreement": float(item.sign_agreement),
                "difference_iqr": float(
                    np.percentile(differences, 75)
                    - np.percentile(differences, 25)
                ),
                **feature,
            }
        )

    control_rows = []
    for row in controls.itertuples(index=False):
        item = pd.Series(row._asdict())
        swap = swap_auditor.evaluate(
            str(item.document_a),
            str(item.document_b),
            target,
            endpoint,
        )
        control_rows.append(
            {
                "dataset_kind": dataset_kind,
                "target": target,
                "endpoint": endpoint,
                "document_a": str(item.document_a),
                "document_b": str(item.document_b),
                "shared_entries": int(item.shared_entries),
                "touched_entries": len(touched_entries(records, item)),
                "s1": bool(swap["is_swap"]),
                "s1_target": swap["target_b"],
                "s1_shared_parents": swap["shared_parents"],
                "s1_swapped_parents": swap["swapped_parents"],
                "s1_swap_fraction": swap["swap_fraction"],
            }
        )
    return feature_rows, control_rows


def detailed_records() -> dict[tuple[str, str], pd.DataFrame]:
    result = {}
    for target in ("CHEMBL204", "CHEMBL205", "CHEMBL240"):
        endpoint = "IC50"
        path = S / f"data/chembl37_api/activities/{target}_{endpoint}.csv"
        raw = pd.read_csv(path, low_memory=False)
        curated = curate_activity_table(raw, target, endpoint)
        result[(target, endpoint)] = curated.retained_records
    return result


def run(
    *,
    records_by_key: dict[tuple[str, str], pd.DataFrame],
    e11_records: dict[tuple[str, str], pd.DataFrame],
    out_dir: Path,
) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    inclusion = pd.read_csv(S / "output/E11/target_inclusion.csv")
    included_e11 = {
        (str(row.target), str(row.endpoint))
        for row in inclusion.loc[inclusion.included].itertuples(index=False)
    }
    selected_e11 = {
        key: records for key, records in e11_records.items()
        if key in included_e11
    }
    detailed = detailed_records()

    all_e2_records = pickle.loads(
        gzip.decompress((S / "qa/E2/curated_records.pkl.gz").read_bytes())
    )
    background = {**all_e2_records, **e11_records, **detailed}
    swap_auditor = SwapAuditor(background)
    flagged_rows = []
    control_rows = []
    for kind, records_map in (
        ("52_evaluation", records_by_key),
        ("23_new", selected_e11),
        ("detailed", detailed),
    ):
        for (target, endpoint), records in sorted(records_map.items()):
            features, controls = audit_dataset(
                records,
                dataset_kind=kind,
                target=target,
                endpoint=endpoint,
                swap_auditor=swap_auditor,
            )
            flagged_rows.extend(features)
            control_rows.extend(controls)
    flagged = pd.DataFrame(flagged_rows)
    controls = pd.DataFrame(control_rows)
    if flagged.empty:
        flagged = pd.DataFrame(
            columns=[
                "dataset_kind",
                "target",
                "endpoint",
                "document_a",
                "document_b",
                "shared_entries",
                "touched_entries",
                "median_offset",
                "feature",
                "feature_label",
                "evidence",
            ]
        )
    flagged["document_pair_key"] = [
        ";".join(pair_key(a, b))
        for a, b in zip(flagged.document_a, flagged.document_b)
    ]
    flagged["feature_rank"] = flagged.feature.map(
        {name: index for index, name in enumerate(FEATURE_ORDER)}
    )
    flagged = flagged.sort_values(
        [
            "touched_entries",
            "shared_entries",
            "dataset_kind",
            "target",
            "document_a",
            "document_b",
        ],
        ascending=[False, False, True, True, True, True],
    ).reset_index(drop=True)
    flagged.to_csv(out_dir / "q3_flagged_pairs.csv", index=False)

    summary_rows = []
    summary_rows.append(
        {
            "metric": "flagged_pair_events",
            "n": int(len(flagged)),
            "unique_pairs": int(flagged.document_pair_key.nunique()),
            "touched_entries_sum": int(flagged.touched_entries.sum())
            if len(flagged)
            else 0,
        }
    )
    for feature in FEATURE_ORDER:
        selected = flagged.loc[flagged.feature == feature]
        pair_row = share_row(feature, len(selected), len(flagged))
        pair_row.update(
            {
                "row_type": "pair_proportion",
                "unique_pairs": int(selected.document_pair_key.nunique()),
                "touched_entries_sum": int(selected.touched_entries.sum()),
            }
        )
        summary_rows.append(pair_row)
        touched_total = int(flagged.touched_entries.sum()) if len(flagged) else 0
        touched_success = int(selected.touched_entries.sum())
        entry_row = share_row(feature, touched_success, touched_total)
        entry_row.update(
            {
                "row_type": "touched_entry_proportion",
                "unique_pairs": int(selected.document_pair_key.nunique()),
                "touched_entries_sum": touched_success,
            }
        )
        summary_rows.append(entry_row)

    data_error = flagged.loc[flagged.feature.isin(("S1", "S2", "S3"))]
    pair_fraction = len(data_error) / len(flagged) if len(flagged) else np.nan
    touched_total = int(flagged.touched_entries.sum()) if len(flagged) else 0
    entry_fraction = (
        data_error.touched_entries.sum() / touched_total
        if touched_total
        else np.nan
    )
    unique_error = data_error.drop_duplicates("document_pair_key")
    unique_datasets = set(
        zip(unique_error.target.astype(str), unique_error.endpoint.astype(str))
    )
    if pair_fraction >= 0.25 or entry_fraction >= 0.25:
        q3_decision = "strong"
    elif len(unique_error) >= 3 and len(unique_datasets) >= 2:
        q3_decision = "weak"
    else:
        q3_decision = "none"
    summary_rows.extend(
        [
            {
                "row_type": "data_error",
                "metric": "S1_S2_S3_pair_events",
                "n": int(len(flagged)),
                "successes": int(len(data_error)),
                "proportion": pair_fraction,
            },
            {
                "row_type": "data_error",
                "metric": "S1_S2_S3_touched_entries",
                "n": touched_total,
                "successes": int(data_error.touched_entries.sum()),
                "proportion": entry_fraction,
            },
            {
                "row_type": "data_error",
                "metric": "S1_S2_S3_unique_pairs",
                "n": int(len(unique_error)),
                "successes": int(len(unique_datasets)),
                "proportion": np.nan,
            },
            {
                "row_type": "decision",
                "metric": "Q3",
                "decision": q3_decision,
            },
        ]
    )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(out_dir / "q3_summary.csv", index=False)
    flagged.head(15).to_csv(out_dir / "q3_top15.csv", index=False)

    if controls.empty:
        controls = pd.DataFrame(
            columns=[
                "dataset_kind",
                "target",
                "endpoint",
                "document_a",
                "document_b",
                "s1",
                "touched_entries",
            ]
        )
    control_rows_out = [
        {
            **share_row(
                "S1_control_pairs",
                int(controls.s1.sum()),
                len(controls),
            ),
            "row_type": "pair_proportion",
            "touched_entries_sum": int(
                controls.loc[controls.s1, "touched_entries"].sum()
            ),
        },
        {
            **share_row(
                "S1_control_touched_entries",
                int(controls.loc[controls.s1, "touched_entries"].sum()),
                int(controls.touched_entries.sum()),
            ),
            "row_type": "touched_entry_proportion",
            "touched_entries_sum": int(
                controls.loc[controls.s1, "touched_entries"].sum()
            ),
        },
    ]
    pd.DataFrame(control_rows_out).to_csv(
        out_dir / "q3_control_s1.csv", index=False
    )
    return {
        "decision": q3_decision,
        "flagged_pair_events": int(len(flagged)),
        "flagged_unique_pairs": int(flagged.document_pair_key.nunique()),
        "flagged_touched_entries": touched_total,
        "S1_S2_S3_pair_events": int(len(data_error)),
        "S1_S2_S3_pair_fraction": float(pair_fraction),
        "S1_S2_S3_touched_entries": int(data_error.touched_entries.sum()),
        "S1_S2_S3_touched_fraction": float(entry_fraction),
        "S1_S2_S3_unique_pairs": int(len(unique_error)),
        "S1_S2_S3_datasets": int(len(unique_datasets)),
        "control_pairs": int(len(controls)),
        "control_S1_pairs": int(controls.s1.sum()),
        "control_S1_touched_entries": int(
            controls.loc[controls.s1, "touched_entries"].sum()
        ),
        "control_touched_entries": int(controls.touched_entries.sum()),
    }


if __name__ == "__main__":
    raise SystemExit(
        "Run e15_leakfree.py; it supplies the curated records to this audit."
    )
