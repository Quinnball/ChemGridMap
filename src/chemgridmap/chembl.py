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
from .structures import parent_structure
from .activity import RULE_KEYS, activity_annotations


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
    "parent_molecule_id": [
        "parent molecule chembl id",
        "parent_molecule_chembl_id",
        "parent compound chembl id",
        "parent_compound_chembl_id",
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
    "assay_type": ["assay type", "assay_type"],
    "assay_description": ["assay description", "assay_description"],
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


def curate_chembl_activity_data(
    data: pd.DataFrame,
    target_id: Optional[str] = None,
    activity_type: str = "IC50",
    lower_threshold: float = 6.0,
    upper_threshold: float = 7.0,
    conflict_range_threshold: float = 1.0,
    exclude_potential_duplicates: bool = True,
    molecule_identity: str = "auto",
    structure_policy: str = "auto",
    assay_types: Optional[List[str]] = None,
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
    requested_identity = str(molecule_identity).strip().lower()
    if requested_identity not in {"auto", "parent-id", "canonical-smiles"}:
        raise ValueError(
            "molecule_identity must be 'auto', 'parent-id' or 'canonical-smiles'."
        )

    if structure_policy not in {"auto", "parent", "as-recorded"}:
        raise ValueError("structure_policy must be auto, parent or as-recorded.")
    effective_structure_policy = structure_policy
    if structure_policy == "auto":
        effective_structure_policy = (
            "as-recorded" if requested_identity == "canonical-smiles" else "parent"
        )

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
    if requested_identity == "parent-id" and columns["parent_molecule_id"] is None:
        raise ValueError(
            "Parent molecule identity was requested, but no "
            "parent_molecule_chembl_id column was found."
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

    if assay_types:
        if columns["assay_type"] is None:
            raise ValueError("Assay-type filtering requires an Assay Type column.")
        allowed = {str(value).strip().upper() for value in assay_types}
        observed = records[columns["assay_type"]].fillna("").astype(str).str.upper()
        flag(~observed.isin(allowed), "assay_type_mismatch")

    if columns["standard_relation"] is not None:
        # The web CSV wraps operators in apostrophes to avoid spreadsheet formulas.
        observed = records[columns["standard_relation"]].astype(str).str.strip().str.replace(
            r"^'([=<>~]+)'$", r"\1", regex=True
        )
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

    records["original_smiles"] = records[columns["smiles"]]
    records["activity_pchembl_raw"] = pchembl
    records["canonical_smiles"] = canonical
    records["record_canonical_smiles"] = canonical
    records["exclusion_reason"] = [
        ";".join(sorted(set(reasons))) for reasons in exclusion_reasons
    ]

    excluded = records[records["exclusion_reason"].ne("")].copy()
    retained = records[records["exclusion_reason"].eq("")].copy()
    retained = retained.reset_index(drop=True)
    if retained.empty:
        raise ValueError("No ChEMBL activity records remained after curation.")

    if effective_structure_policy == "parent":
        structures = retained["canonical_smiles"].map(parent_structure)
        retained["representation_smiles"] = structures.map(lambda item: item[0])
        retained["parent_extraction_excluded"] = structures.map(lambda item: int(item[1]))
    else:
        retained["representation_smiles"] = retained["canonical_smiles"]
        retained["parent_extraction_excluded"] = 0
    retained["structure_changed"] = retained["representation_smiles"].ne(
        retained["record_canonical_smiles"]
    ).astype(int)

    parent_column = columns["parent_molecule_id"]
    use_parent_identity = requested_identity == "parent-id"
    if requested_identity == "auto" and parent_column is not None:
        parent_values = retained[parent_column].fillna("").astype(str).str.strip()
        use_parent_identity = bool(parent_values.ne("").any())
    if use_parent_identity:
        parent_values = retained[parent_column].fillna("").astype(str).str.strip()
        parent_values = parent_values.str.upper()
        # A missing ID can be recovered only when the same representation has
        # one unambiguous parent elsewhere in this export.
        known = pd.DataFrame({"structure": retained["representation_smiles"],
                              "parent": parent_values})
        known = known[known["parent"].ne("")].groupby("structure")["parent"].unique()
        unique_parents = {key: value[0] for key, value in known.items() if len(value) == 1}
        recovered = retained["representation_smiles"].map(unique_parents).fillna("")
        inferred = parent_values.eq("") & recovered.ne("")
        parent_values = parent_values.mask(inferred, recovered)
        has_parent = parent_values.ne("")
        retained["molecule_identity_key"] = np.where(
            has_parent,
            parent_values,
            "SMILES:" + retained["representation_smiles"].astype(str),
        )
        retained["molecule_identity_source"] = np.where(
            has_parent,
            "parent_molecule_chembl_id",
            "canonical_smiles_fallback",
        )
        retained.loc[inferred, "molecule_identity_source"] = "unambiguous_structure_parent_lookup"
        effective_identity = "parent_molecule_chembl_id_with_smiles_fallback"
        identity_fallback_records = int((~has_parent).sum())
    else:
        retained["molecule_identity_key"] = (
            "SMILES:" + retained["representation_smiles"].astype(str)
        )
        retained["molecule_identity_source"] = "canonical_smiles"
        effective_identity = "canonical_representation_smiles"
        identity_fallback_records = 0
        if requested_identity == "auto" and parent_column is None:
            warnings.append(
                "Parent molecule ID column was absent; molecule identity fell "
                "back to canonical representation SMILES under the selected structure policy."
            )

    grouped = retained.groupby("molecule_identity_key", sort=True)
    molecule_rows = []
    for identity_key, group in grouped:
        structures = sorted(group["representation_smiles"].astype(str).unique())
        if len(structures) > 1 and effective_structure_policy == "parent":
            raise ValueError(
                "Identity {} maps to multiple normalized parent structures. "
                "Review the records, or explicitly use canonical-smiles identity.".format(identity_key)
            )
        canonical_smiles = structures[0]
        values = group["activity_pchembl_raw"].to_numpy(dtype=float)
        rules = dict(zip(RULE_KEYS, (lower_threshold, upper_threshold, conflict_range_threshold)))
        annotations = activity_annotations(values, **rules)

        row = {
            "molecule_identity_key": identity_key,
            "molecule_identity_source": _joined_unique(
                group["molecule_identity_source"]
            ),
            "canonical_smiles": canonical_smiles,
            "smiles": canonical_smiles,
            "representation_smiles": canonical_smiles,
            "record_canonical_smiles": _joined_unique(group["record_canonical_smiles"]),
            "structure_policy": effective_structure_policy,
            "structure_changed": int(group["structure_changed"].any()),
            "parent_extraction_excluded": int(group["parent_extraction_excluded"].any()),
            "n_record_structures": int(group["record_canonical_smiles"].nunique()),
            **annotations,
            **{"audit_" + key: float(value) for key, value in rules.items()},
            "source_rows": _joined_unique(group["source_row"]),
        }
        for field in [
            "molecule_id",
            "parent_molecule_id",
            "target_id",
            "target_name",
            "assay_id",
            "document_id",
            "activity_id",
            "assay_type",
        ]:
            column = columns[field]
            if column is not None:
                row[field] = _joined_unique(group[column])
        row["n_assays"] = (int(group[columns["assay_id"]].nunique())
                           if columns["assay_id"] else np.nan)
        molecule_rows.append(row)

    molecules = pd.DataFrame(molecule_rows).sort_values("molecule_identity_key")
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
        "repeated_unflagged_molecules": int(
            repeat_counts.get("repeated_unflagged", 0)
        ),
        "repeated_conflicted_molecules": int(
            repeat_counts.get("repeated_conflicted", 0)
        ),
        "class_boundary_crossing_molecules": int(molecules["has_class_boundary_crossing"].sum()),
        "unflagged_class_boundary_crossing_molecules": int((
            molecules["has_class_boundary_crossing"].eq(1) & molecules["is_conflicted"].eq(0)
        ).sum()),
        "structure_changed_records": int(retained["structure_changed"].sum()),
        "structure_changed_molecules": int(molecules["structure_changed"].sum()),
        "parent_extraction_excluded_records": int(retained["parent_extraction_excluded"].sum()),
        "retained_assay_type_counts": (
            {str(key): int(value) for key, value in retained[columns["assay_type"]].fillna("missing").value_counts().items()}
            if columns["assay_type"] else {}
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
            "molecule_identity_requested": requested_identity,
            "molecule_identity": effective_identity,
            "identity_fallback_records": identity_fallback_records,
            "aggregation": "median pChEMBL",
            "structure_policy": effective_structure_policy,
            "structure_standardizer": (
                "ChEMBL Structure Pipeline Standardizer + GetParent"
                if effective_structure_policy == "parent" else "none"
            ),
            "assay_types": assay_types,
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
