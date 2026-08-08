"""ChEMBL activity-table import and molecule-level curation."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import re
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from .core import canonicalize_smiles


COLUMN_ALIASES = {
    "smiles": [
        "smiles",
        "canonical smiles",
        "canonical_smiles",
        "molecule smiles",
    ],
    "pchembl": [
        "pchembl",
        "pchembl value",
        "pchembl_value",
    ],
    "molecule_id": [
        "molecule chembl id",
        "molecule_chembl_id",
        "compound chembl id",
        "compound_chembl_id",
        "compound",
    ],
    "target_id": [
        "target chembl id",
        "target_chembl_id",
    ],
    "target_name": [
        "target name",
        "target_name",
        "target pref name",
    ],
    "standard_type": [
        "standard type",
        "standard_type",
    ],
    "standard_relation": [
        "standard relation",
        "standard_relation",
    ],
    "standard_units": [
        "standard units",
        "standard_units",
    ],
    "data_validity_comment": [
        "data validity comment",
        "data_validity_comment",
    ],
    "potential_duplicate": [
        "potential duplicate",
        "potential_duplicate",
    ],
    "assay_id": [
        "assay chembl id",
        "assay_chembl_id",
    ],
    "document_id": [
        "document chembl id",
        "document_chembl_id",
    ],
    "activity_id": [
        "activity id",
        "activity_id",
    ],
}


@dataclass
class ChemblCurationResult:
    """Tables and audit information produced from a ChEMBL activity export."""

    molecules: pd.DataFrame
    retained_records: pd.DataFrame
    excluded_records: pd.DataFrame
    report: Dict[str, object]


def _normalize_column_name(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).strip().lower())


def detect_chembl_columns(data: pd.DataFrame) -> Dict[str, Optional[str]]:
    """Match common ChEMBL web-export and API column names."""
    normalized = {}
    for column in data.columns:
        normalized.setdefault(_normalize_column_name(column), str(column))

    detected: Dict[str, Optional[str]] = {}
    for field, aliases in COLUMN_ALIASES.items():
        detected[field] = None
        for alias in aliases:
            match = normalized.get(_normalize_column_name(alias))
            if match is not None:
                detected[field] = match
                break
    return detected


def _truthy_flag(value: object) -> bool:
    if pd.isna(value):
        return False
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(value) != 0.0
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def _joined_unique(values: pd.Series) -> str:
    cleaned = {
        str(value).strip()
        for value in values
        if not pd.isna(value) and str(value).strip()
    }
    return ";".join(sorted(cleaned))


def _activity_class(
    value: float,
    lower_threshold: float,
    upper_threshold: float,
) -> str:
    if value <= lower_threshold:
        return "inactive"
    if value >= upper_threshold:
        return "active"
    return "medium"


def curate_chembl_activity_data(
    data: pd.DataFrame,
    target_id: Optional[str] = None,
    activity_type: str = "IC50",
    lower_threshold: float = 6.0,
    upper_threshold: float = 7.0,
    conflict_range_threshold: float = 1.0,
    exclude_potential_duplicates: bool = True,
) -> ChemblCurationResult:
    """Convert a ChEMBL activity export into an audited molecule-level table.

    Non-null pChEMBL values are used as supplied by ChEMBL. When the related
    metadata columns are present, the importer also enforces exact relations,
    nM standard units, accepted validity flags, the requested target, and the
    requested activity type.
    """
    if lower_threshold >= upper_threshold:
        raise ValueError("lower_threshold must be smaller than upper_threshold.")
    if conflict_range_threshold < 0:
        raise ValueError("conflict_range_threshold must be non-negative.")

    columns = detect_chembl_columns(data)
    if columns["smiles"] is None:
        raise ValueError(
            "No SMILES column was found. Expected a ChEMBL column such as "
            "'Smiles' or 'Canonical Smiles'."
        )
    if columns["pchembl"] is None:
        raise ValueError(
            "No pChEMBL column was found. Export activities with the "
            "'pChEMBL Value' field included."
        )

    warnings: List[str] = []
    effective_target_id = target_id
    if columns["target_id"] is not None and target_id is None:
        observed_targets = sorted(
            {
                str(value).strip().upper()
                for value in data[columns["target_id"]]
                if not pd.isna(value) and str(value).strip()
            }
        )
        if len(observed_targets) > 1:
            raise ValueError(
                "Multiple Target ChEMBL IDs were found: {}. "
                "Select one with --target-id.".format(
                    ", ".join(observed_targets)
                )
            )
        if len(observed_targets) == 1:
            effective_target_id = observed_targets[0]
    if columns["target_id"] is None:
        warnings.append(
            "Target ChEMBL ID column was absent; target identity could not be "
            "verified from the file."
        )
    if columns["standard_type"] is None:
        warnings.append(
            "Standard Type column was absent; the requested activity type "
            "could not be verified from the file."
        )
    if columns["standard_relation"] is None:
        warnings.append(
            "Standard Relation column was absent; exact-relation filtering "
            "could not be repeated from the file."
        )
    if columns["standard_units"] is None:
        warnings.append(
            "Standard Units column was absent; nM units could not be verified "
            "from the file."
        )
    if columns["data_validity_comment"] is None:
        warnings.append(
            "Data Validity Comment column was absent; ChEMBL validity flags "
            "could not be audited from the file."
        )
    if exclude_potential_duplicates and columns["potential_duplicate"] is None:
        warnings.append(
            "Potential Duplicate column was absent; duplicate-citation flags "
            "could not be removed."
        )

    records = data.copy()
    records.insert(0, "source_row", np.arange(len(records), dtype=int))
    exclusion_reasons: List[List[str]] = [[] for _ in range(len(records))]

    def flag(mask: pd.Series, reason: str) -> None:
        for index in np.flatnonzero(mask.to_numpy(dtype=bool)):
            exclusion_reasons[int(index)].append(reason)

    pchembl = pd.to_numeric(records[columns["pchembl"]], errors="coerce")
    flag(~np.isfinite(pchembl), "missing_or_non_numeric_pchembl")

    if effective_target_id and columns["target_id"] is not None:
        observed = records[columns["target_id"]].astype(str).str.strip().str.upper()
        flag(
            observed.ne(str(effective_target_id).strip().upper()),
            "target_mismatch",
        )

    if columns["standard_type"] is not None:
        observed = records[columns["standard_type"]].astype(str).str.strip().str.upper()
        flag(observed.ne(str(activity_type).strip().upper()), "activity_type_mismatch")

    if columns["standard_relation"] is not None:
        observed = records[columns["standard_relation"]].astype(str).str.strip()
        flag(observed.ne("="), "non_exact_standard_relation")

    if columns["standard_units"] is not None:
        observed = records[columns["standard_units"]].astype(str).str.strip().str.lower()
        flag(observed.ne("nm"), "non_nm_standard_units")

    if columns["data_validity_comment"] is not None:
        validity = records[columns["data_validity_comment"]]
        normalized = validity.fillna("").astype(str).str.strip().str.lower()
        accepted = normalized.isin({"", "nan", "none", "manually validated"})
        flag(~accepted, "data_validity_flag")

    if exclude_potential_duplicates and columns["potential_duplicate"] is not None:
        duplicate_flags = records[columns["potential_duplicate"]].map(_truthy_flag)
        flag(duplicate_flags, "chembl_potential_duplicate")

    canonical = records[columns["smiles"]].map(canonicalize_smiles)
    flag(canonical.isna(), "invalid_or_missing_smiles")

    records["activity_pchembl_raw"] = pchembl
    records["canonical_smiles"] = canonical
    records["exclusion_reason"] = [
        ";".join(sorted(set(reasons))) for reasons in exclusion_reasons
    ]

    excluded = records[records["exclusion_reason"].ne("")].copy()
    retained = records[records["exclusion_reason"].eq("")].copy()
    retained = retained.reset_index(drop=True)
    if retained.empty:
        raise ValueError("No ChEMBL activity records remained after curation.")

    grouped = retained.groupby("canonical_smiles", sort=True)
    molecule_rows = []
    for canonical_smiles, group in grouped:
        values = group["activity_pchembl_raw"].to_numpy(dtype=float)
        median = float(np.median(values))
        minimum = float(np.min(values))
        maximum = float(np.max(values))
        record_range = maximum - minimum
        has_inactive = bool(np.any(values <= lower_threshold))
        has_active = bool(np.any(values >= upper_threshold))
        crosses_clean_boundary = has_inactive and has_active
        is_large_variation = record_range > conflict_range_threshold
        is_conflicted = crosses_clean_boundary or is_large_variation
        n_records = int(len(group))

        row = {
            "canonical_smiles": canonical_smiles,
            "smiles": canonical_smiles,
            "activity_pchembl": median,
            "activity_class": _activity_class(
                median,
                lower_threshold=lower_threshold,
                upper_threshold=upper_threshold,
            ),
            "n_records": n_records,
            "pchembl_min": minimum,
            "pchembl_max": maximum,
            "pchembl_median": median,
            "pchembl_std": (
                float(np.std(values, ddof=1)) if len(values) > 1 else 0.0
            ),
            "pchembl_range": record_range,
            "has_inactive_measurement": int(has_inactive),
            "has_active_measurement": int(has_active),
            "crosses_clean_boundary": int(crosses_clean_boundary),
            "is_large_variation": int(is_large_variation),
            "is_conflicted": int(is_conflicted),
            "repeat_status": (
                "singleton"
                if n_records == 1
                else (
                    "repeated_conflicted"
                    if is_conflicted
                    else "repeated_consistent"
                )
            ),
            "source_rows": _joined_unique(group["source_row"]),
        }
        for field in [
            "molecule_id",
            "target_id",
            "target_name",
            "assay_id",
            "document_id",
            "activity_id",
        ]:
            column = columns[field]
            if column is not None:
                row[field] = _joined_unique(group[column])
        molecule_rows.append(row)

    molecules = pd.DataFrame(molecule_rows).sort_values("canonical_smiles")
    molecules = molecules.reset_index(drop=True)

    exclusion_counts: Dict[str, int] = {}
    for reasons in exclusion_reasons:
        for reason in set(reasons):
            exclusion_counts[reason] = exclusion_counts.get(reason, 0) + 1

    repeat_counts = molecules["repeat_status"].value_counts().to_dict()
    class_counts = molecules["activity_class"].value_counts().to_dict()
    report: Dict[str, object] = {
        "input_rows": int(len(records)),
        "retained_activity_records": int(len(retained)),
        "excluded_activity_records": int(len(excluded)),
        "unique_molecules": int(len(molecules)),
        "repeated_molecules": int((molecules["n_records"] > 1).sum()),
        "repeated_consistent_molecules": int(
            repeat_counts.get("repeated_consistent", 0)
        ),
        "repeated_conflicted_molecules": int(
            repeat_counts.get("repeated_conflicted", 0)
        ),
        "activity_class_counts": {
            str(key): int(value) for key, value in class_counts.items()
        },
        "exclusion_counts": {
            str(key): int(value) for key, value in sorted(exclusion_counts.items())
        },
        "detected_columns": columns,
        "parameters": {
            "target_id": effective_target_id,
            "activity_type": activity_type,
            "lower_threshold": float(lower_threshold),
            "upper_threshold": float(upper_threshold),
            "conflict_range_threshold": float(conflict_range_threshold),
            "exclude_potential_duplicates": bool(exclude_potential_duplicates),
            "molecule_identity": "RDKit canonical isomeric SMILES",
            "aggregation": "median pChEMBL",
            "salt_stripping": False,
        },
        "warnings": warnings,
    }
    return ChemblCurationResult(
        molecules=molecules,
        retained_records=retained,
        excluded_records=excluded.reset_index(drop=True),
        report=report,
    )


def save_chembl_curation(
    result: ChemblCurationResult,
    output_dir: Path,
    name: str,
) -> Dict[str, str]:
    """Write curated records, molecule-level data, exclusions, and an audit."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "molecule_level": output_dir / f"{name}_molecule_level.csv",
        "retained_records": output_dir / f"{name}_retained_records.csv",
        "excluded_records": output_dir / f"{name}_excluded_records.csv",
        "curation_report": output_dir / f"{name}_curation_report.json",
    }
    result.molecules.to_csv(paths["molecule_level"], index=False)
    result.retained_records.to_csv(paths["retained_records"], index=False)
    result.excluded_records.to_csv(paths["excluded_records"], index=False)
    with paths["curation_report"].open("w", encoding="utf-8") as handle:
        json.dump(result.report, handle, indent=2, ensure_ascii=True)
    return {key: str(path) for key, path in paths.items()}
