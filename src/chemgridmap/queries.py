"""Unmeasured structures and measured fingerprint neighbours, without predictions."""

from __future__ import annotations

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

from .comparison import molecule_keys
from .core import canonicalize_smiles


QUERY_CHECK_COLUMNS = [
    "source_row", "query_id", "input_smiles", "canonical_smiles", "status",
    "first_query_id", "reference_id",
]
QUERY_NEIGHBOR_COLUMNS = [
    "query_id", "rank", "measured_id", "tanimoto", "activity_class",
    "activity_pchembl", "class_confidence", "class_status", "n_evidence_units",
]


def prepare_query_molecules(
    query_data: pd.DataFrame,
    reference_data: pd.DataFrame,
    *,
    smiles_column: str | None = None,
    id_column: str | None = None,
    structure_policy: str = "canonical",
):
    """Check unmeasured molecules against the prepared measured reference table.

    Parameters
    ----------
    query_data : pandas.DataFrame
        Unmeasured structures, in source order. Activity annotations are ignored.
    reference_data : pandas.DataFrame
        Prepared measured molecules with ``canonical_smiles``. Pass the actual
        reference structure policy, including ``parent`` for ChEMBL parents.
    smiles_column : str, optional
        Query SMILES column. Otherwise recognize SMILES, canonical_smiles,
        canonical smiles, or representation_smiles, without case sensitivity.
    id_column : str, optional
        Unique, nonempty query identifiers. Otherwise generate Q00001, Q00002,
        and so on, counting every input row, including excluded rows.
    structure_policy : {"canonical", "as-recorded", "parent"}, default "canonical"
        Apply canonical isomeric SMILES or the same ChEMBL parent standardizer
        used for reference molecules. Identity comparisons use these SMILES.

    Returns
    -------
    molecules : pandas.DataFrame
        Valid, first-occurrence structures absent from the reference, with
        ``measurement_status`` and ``activity_class`` both ``not_measured``.
    check : pandas.DataFrame
        One row per input with zero-based source_row, IDs, input and canonical
        SMILES, status, first_query_id for duplicates, and matching reference_id.
    report : dict
        Input and four status counts, selected columns and structure policy.

    Notes
    -----
    This function never infers activity or class from structural similarity.
    ``build_grid_map(query_data=...)`` calls this same public preparation path.
    """
    if structure_policy == "canonical" and "structure_policy" in reference_data:
        policies = reference_data.structure_policy.dropna().unique()
        if len(policies) == 1:
            structure_policy = str(policies[0])
    if structure_policy not in {"canonical", "as-recorded", "parent"}:
        raise ValueError("query structure_policy must be canonical, as-recorded, or parent.")
    if smiles_column is None:
        available = {str(c).strip().lower(): c for c in query_data.columns}
        smiles_column = next((available[c] for c in (
            "smiles", "canonical_smiles", "canonical smiles", "representation_smiles"
        ) if c in available), None)
    if smiles_column is None or smiles_column not in query_data:
        raise ValueError("Query SMILES column was not found; supply query_smiles_column.")
    if id_column is not None and id_column not in query_data:
        raise ValueError(f"Query ID column {id_column!r} was not found.")
    if id_column is None:
        identifiers = [f"Q{i + 1:05d}" for i in range(len(query_data))]
    else:
        ids = query_data[id_column]
        if ids.isna().any() or ids.astype(str).str.strip().eq("").any():
            raise ValueError("Query IDs must be nonempty.")
        identifiers = ids.astype(str).str.strip().tolist()
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("Query IDs must be unique; duplicate query IDs were found.")

    reference_by_smiles = {}
    for smiles, identifier in zip(reference_data.canonical_smiles, molecule_keys(reference_data)):
        reference_by_smiles.setdefault(str(smiles), identifier)
    seen, checks, mapped = {}, [], []
    for source_row, (query_id, input_smiles) in enumerate(zip(identifiers, query_data[smiles_column])):
        canonical = canonicalize_smiles(input_smiles)
        if canonical is not None and structure_policy == "parent":
            from .structures import parent_structure
            canonical, _ = parent_structure(canonical)
        row = dict(source_row=source_row, query_id=query_id, input_smiles=input_smiles,
                   canonical_smiles=canonical, status="mapped", first_query_id="", reference_id="")
        if canonical is None:
            row["status"] = "invalid_smiles"
        elif canonical in seen:
            row.update(status="duplicate_query", first_query_id=seen[canonical],
                       reference_id=reference_by_smiles.get(canonical, ""))
        else:
            seen[canonical] = query_id
            if canonical in reference_by_smiles:
                row.update(status="in_reference", reference_id=reference_by_smiles[canonical])
            else:
                mapped.append(dict(source_row=source_row, query_id=query_id, molecule_id=query_id,
                                   molecule_identity_key=f"query:{query_id}", name=query_id,
                                   smiles=canonical, canonical_smiles=canonical,
                                   representation_smiles=canonical,
                                   measurement_status="not_measured", activity_class="not_measured"))
        checks.append(row)
    columns = ["source_row", "query_id", "molecule_id", "molecule_identity_key", "name", "smiles",
               "canonical_smiles", "representation_smiles", "measurement_status", "activity_class"]
    check = pd.DataFrame(checks, columns=QUERY_CHECK_COLUMNS)
    report = {"input_rows": len(query_data), **{
        status: int(check.status.eq(status).sum())
        for status in ("mapped", "invalid_smiles", "duplicate_query", "in_reference")
    }, "smiles_column": str(smiles_column), "id_column": id_column,
        "structure_policy": structure_policy, "warnings": []}
    return pd.DataFrame(mapped, columns=columns), check, report


def annotate_query_neighbors(reference, queries, *, radius=2, n_bits=2048, value_col=None):
    """Return query rows and the five most similar measured neighbours, including ties."""
    out = queries.copy()
    out["nearest_measured_id"] = ""
    out["nearest_measured_similarity"] = np.nan
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=int(radius), fpSize=int(n_bits))
    reference_fps = [generator.GetFingerprint(Chem.MolFromSmiles(s)) for s in reference.canonical_smiles]
    reference_ids = molecule_keys(reference)
    rows = []
    for index, query in out.iterrows():
        fingerprint = generator.GetFingerprint(Chem.MolFromSmiles(query.canonical_smiles))
        similarity = np.asarray(DataStructs.BulkTanimotoSimilarity(fingerprint, reference_fps))
        order = sorted(range(len(reference)), key=lambda i: (-similarity[i], reference_ids[i]))
        cutoff = similarity[order[min(5, len(order)) - 1]]
        out.at[index, "nearest_measured_id"] = reference_ids[order[0]]
        out.at[index, "nearest_measured_similarity"] = float(similarity[order[0]])
        rank, previous = 0, None
        for position, ref_index in enumerate(order, start=1):
            score = float(similarity[ref_index])
            if score < cutoff - 1e-12:
                break
            if previous is None or score != previous:
                rank = position
            previous = score
            measured = reference.iloc[ref_index]
            rows.append(dict(query_id=query.query_id, rank=rank, measured_id=reference_ids[ref_index],
                             tanimoto=score, activity_class=measured.get("activity_class"),
                             activity_pchembl=measured.get("activity_pchembl", measured.get(value_col)),
                             class_confidence=measured.get("class_confidence"),
                             class_status=measured.get("class_status"),
                             n_evidence_units=measured.get("n_evidence_units")))
    return out, pd.DataFrame(rows, columns=QUERY_NEIGHBOR_COLUMNS)


def assign_query_groups(queries, groups, *, grouping, similarity, radius, n_bits):
    """Attach queries to fixed reference centers without changing reference cores."""
    out = queries.copy()
    out["comparison_id"] = out.molecule_identity_key
    out["structure_group"] = ""
    out["structure_group_color"] = ""
    out["scaffold_smiles"] = [MurckoScaffold.MurckoScaffoldSmiles(
        mol=Chem.MolFromSmiles(s), includeChirality=True) for s in out.canonical_smiles]
    out["query_group_center_similarity"] = np.nan
    out["contains_group_core"] = pd.Series([None] * len(out), index=out.index, dtype=object)
    summaries = [{**group, "n_query_members": 0} for group in groups]
    if grouping != "butina" or not summaries or out.empty:
        return out, summaries
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=int(radius), fpSize=int(n_bits))
    centers = sorted(summaries, key=lambda group: int(group["structure_group"][1:]))
    fingerprints = [generator.GetFingerprint(Chem.MolFromSmiles(group["representative_smiles"]))
                    for group in centers]
    for index, query in out.iterrows():
        molecule = Chem.MolFromSmiles(query.canonical_smiles)
        scores = np.asarray(DataStructs.BulkTanimotoSimilarity(generator.GetFingerprint(molecule), fingerprints))
        best = int(np.argmax(scores))
        if scores[best] < similarity:
            continue
        group = centers[best]
        out.at[index, "structure_group"] = group["structure_group"]
        out.at[index, "structure_group_color"] = group["color"]
        out.at[index, "query_group_center_similarity"] = float(scores[best])
        pattern = Chem.MolFromSmarts(group["core_smarts"]) if group.get("core_smarts") else (
            Chem.MolFromSmiles(group["core_smiles"]) if group.get("core_smiles") else None)
        out.at[index, "contains_group_core"] = bool(
            pattern is not None and molecule.HasSubstructMatch(pattern, useChirality=True))
        group["n_query_members"] += 1
    return out, summaries
