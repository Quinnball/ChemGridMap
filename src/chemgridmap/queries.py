"""Unmeasured structures and measured fingerprint neighbours, without predictions."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import hashlib
import html
import json
import time
import zipfile
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

from .chembl import detect_chembl_columns
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


def annotate_query_neighbors(reference, queries, *, radius=2, n_bits=2048, value_col=None,
                             reference_fingerprints=None):
    """Return query rows and the five most similar measured neighbours, including ties."""
    out = queries.copy()
    out["nearest_measured_id"] = ""
    out["nearest_measured_similarity"] = np.nan
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=int(radius), fpSize=int(n_bits))
    reference_fps = reference_fingerprints if reference_fingerprints is not None else [
        generator.GetFingerprint(Chem.MolFromSmiles(s)) for s in reference.canonical_smiles]
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


# 新分子只占用参照图的空格；这里不调用投影、网格优化或参照证据计算。
NEIGHBOR_EVIDENCE_CATEGORIES = ("none", "insufficient", "active", "medium", "inactive", "mixed")


@dataclass
class NewMoleculesResult:
    data: pd.DataFrame
    check: pd.DataFrame
    neighbors: pd.DataFrame
    report: dict
    output_files: dict = field(default_factory=dict)
    records: pd.DataFrame = field(default_factory=pd.DataFrame)
    structure_groups: list = field(default_factory=list)
    all_molecules: pd.DataFrame = field(default_factory=pd.DataFrame)


def summarize_neighbor_evidence(queries, neighbors, *, similarity=0.6, support_assessed=True):
    """Summarize only the existing tie-inclusive five-reference neighbor table."""
    out = queries.copy()
    evidence = []
    for _, row in out.iterrows():
        close = neighbors.loc[neighbors.query_id.eq(row.query_id) & neighbors.tanimoto.ge(similarity)]
        supported = close[close.class_status.eq("supported")] if support_assessed else close
        classes = sorted(set(supported.activity_class.dropna()) & {"active", "medium", "inactive"})
        category = ("none" if close.empty else "insufficient" if not classes
                    else classes[0] if len(classes) == 1 else "mixed")
        own = row.get("activity_class")
        comparison = ("agrees" if own == category else "disagrees") if (
            row.get("measurement_status") == "measured" and category in {"active", "medium", "inactive"}
        ) else ""
        evidence.append({"neighbor_evidence": category, "similar_reference_count": len(close),
                         "supported_reference_count": len(supported) if support_assessed else 0,
                         "supported_reference_classes": ";".join(classes) if support_assessed else "",
                         "evidence_reference_count": len(supported),
                         "evidence_reference_classes": ";".join(classes),
                         "reference_support_assessed": bool(support_assessed),
                         "neighbor_comparison": comparison})
    return pd.concat([out.reset_index(drop=True), pd.DataFrame(evidence, columns=[
        "neighbor_evidence", "similar_reference_count", "supported_reference_count",
        "supported_reference_classes", "evidence_reference_count", "evidence_reference_classes",
        "reference_support_assessed", "neighbor_comparison"])], axis=1)


def _insertion_warnings(queries, query_data, reference_target_id=None):
    """只检查已放置分子的结构差异；靶点检查覆盖整个新分子输入。"""
    placed = len(queries)
    none = int(queries.neighbor_evidence.eq("none").sum())
    none_share = none / placed if placed else 0.0
    median = float(queries.tanimoto.median()) if placed else None
    warnings = []
    if none_share >= 0.9 and median is not None and median < 0.35:
        warnings.append({"code": "structural_mismatch", "placed": placed,
                         "none": none, "median": median})
    target_column = detect_chembl_columns(query_data)["target_id"]
    if reference_target_id and target_column:
        reference = str(reference_target_id).strip().upper()
        targets = sorted({str(value).strip().upper() for value in query_data[target_column].dropna()
                          if str(value).strip() and str(value).strip().upper() != reference})
        if targets:
            warnings.append({"code": "target_mismatch", "targets": targets, "reference": reference})
    return {"warnings": warnings, "none_share": none_share,
            "median_nearest_tanimoto": median}


def place_new_molecules(reference, queries, *, grid_rows, grid_cols):
    """Use Euclidean cell distance, then row/column; never move a reference."""
    out = queries.copy()
    occupied = set(zip(reference.grid_row.astype(int), reference.grid_col.astype(int)))
    available = {(row, col) for row in range(grid_rows) for col in range(grid_cols)} - occupied
    if len(out) > len(available):
        raise ValueError("Not enough empty cells in the reference grid.")
    lookup = dict(zip(molecule_keys(reference), zip(reference.grid_row, reference.grid_col)))
    for column in ("grid_row", "grid_col", "grid_x", "grid_y", "placement_distance"):
        out[column] = pd.Series(index=out.index, dtype=float)
    for index in sorted(out.index, key=lambda i: (-float(out.at[i, "nearest_measured_similarity"]),
                                                 str(out.at[i, "query_id"]))):
        target_row, target_col = lookup[str(out.at[index, "nearest_measured_id"])]
        row, col = min(available, key=lambda cell: (
            (cell[0] - target_row) ** 2 + (cell[1] - target_col) ** 2, cell[0], cell[1]))
        available.remove((row, col))
        out.loc[index, ["grid_row", "grid_col", "grid_x", "grid_y", "placement_distance"]] = [
            row, col, col / max(1, grid_cols - 1), row / max(1, grid_rows - 1),
            float(np.hypot(row - target_row, col - target_col))]
    out["grid_row"] = out.grid_row.astype(int)
    out["grid_col"] = out.grid_col.astype(int)
    out["nearest_reference_id"] = out.nearest_measured_id
    out["tanimoto"] = out.nearest_measured_similarity
    return out


def _new_molecules_svg(reference_result, queries, path, *, rows, cols, lower, upper, show_feedback=False):
    """Reuse existing reference SVG groups; depict only newly inserted structures."""
    from .plotting import DEFAULT_COLORS, _svg_molecule_content
    size = 100
    width, map_height = max(cols * size, 1000), rows * size
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{map_height+100}" '
             f'viewBox="0 0 {width} {map_height+100}">', '<rect width="100%" height="100%" fill="white"/>']
    reference = reference_result.data
    source = reference_result.output_files.get("svg")
    if source and Path(source).is_file():
        root = ET.fromstring(Path(source).read_text(encoding="utf-8"))
        dx = int(reference.grid_col.min()) * size
        dy = (rows - 1 - int(reference.grid_row.max())) * size
        parts.append(f'<g data-reference-map="true" transform="translate({dx},{dy})">')
        for node in root:
            if node.get("data-molecule-id") is not None:
                # HTML 的 innerHTML 不识别 ns0:g；根 SVG 已声明默认命名空间。
                # 只改当前复制子树，避免全局 namespace 注册影响其他输出。
                for child in node.iter():
                    child.tag = child.tag.removeprefix("{http://www.w3.org/2000/svg}")
                parts.append(ET.tostring(node, encoding="unicode"))
        parts.append('</g>')
    else:
        # 无 SVG 的内存/仅 PNG 参照图采用同一已固定格子的概览导出。
        for key, (_, row) in zip(molecule_keys(reference), reference.iterrows()):
            x, y = int(row.grid_col) * size, (rows - 1 - int(row.grid_row)) * size
            color = DEFAULT_COLORS.get(str(row.activity_class), DEFAULT_COLORS["unlabelled"])
            parts.append(f'<g data-molecule-id="{html.escape(key, quote=True)}"><rect x="{x}" y="{y}" '
                         f'width="100" height="100" fill="{color}" stroke="white"/></g>')
    for _, row in queries.iterrows():
        x, y = int(row.grid_col) * size, (rows - 1 - int(row.grid_row)) * size
        measured = row.measurement_status == "measured"
        color = DEFAULT_COLORS.get(str(row.activity_class), "#FFFFFF") if measured else "#FFFFFF"
        edge = "#173e34" if measured else "#8A99A6"
        parts.append(f'<g data-molecule-id="{html.escape(str(row.molecule_identity_key), quote=True)}" '
                     f'data-new-molecule="true" data-grid-row="{int(row.grid_row)}" data-grid-col="{int(row.grid_col)}">'
                     f'<title>{html.escape(str(row.query_id))}</title><rect x="{x+1}" y="{y+1}" width="98" height="98" '
                     f'fill="{color}" stroke="{edge}" stroke-width="1.5" vector-effect="non-scaling-stroke"/>')
        drawing = _svg_molecule_content(str(row.canonical_smiles), 74)
        if drawing:
            parts.append(f'<g transform="translate({x+13},{y+13})"><rect width="74" height="74" fill="white"/>{drawing}</g>')
        if measured and row.get("wide_record_range", 0) == 1:
            parts.append(f'<path data-conflict-marker="true" d="M{x+70},{y} L{x+100},{y} L{x+100},{y+30} Z" fill="#172536"/>')
        if show_feedback and row.get("neighbor_comparison") == "disagrees":
            marker_size = 98 * .22
            left, top = x + 1, y + 1
            parts.append(f'<path data-neighbor-mismatch="true" '
                         f'd="M{left+marker_size/2},{top+2} L{left+marker_size},{top+marker_size/2} '
                         f'L{left+marker_size/2},{top+marker_size} L{left+2},{top+marker_size/2} Z" '
                         'fill="#6f1d6b" stroke="#fff" stroke-width="1"/>')
        parts.append('</g>')
    parts.append('<g data-color-legend="activity" font-family="Arial, sans-serif" font-size="14">')
    entries = [(DEFAULT_COLORS["inactive"], f"Inactive ≤ {lower:g}", ""),
               (DEFAULT_COLORS["medium"], f"Medium ({lower:g}, {upper:g})", ""),
               (DEFAULT_COLORS["active"], f"Active ≥ {upper:g}", ""),
               ("#FFFFFF", "New molecules: not measured", "#8A99A6"),
               ("#FFFFFF", "New molecules: measured", "#173e34")]
    if show_feedback:
        entries.append(("#6f1d6b", "Measured result disagrees with neighbor evidence", "#6f1d6b"))
    for i, (color, label, edge) in enumerate(entries):
        x, y = 12 + (i % 3) * 300, map_height + 20 + (i // 3) * 25
        if show_feedback and i == len(entries) - 1:
            parts.append(f'<path data-neighbor-mismatch-legend="true" '
                         f'd="M{x+7},{y} L{x+14},{y+7} L{x+7},{y+14} L{x},{y+7} Z" fill="{color}"/>'
                         f'<text x="{x+20}" y="{y+12}">{html.escape(label)}</text>')
        else:
            parts.append(f'<rect x="{x}" y="{y}" width="14" height="14" fill="{color}" stroke="{edge or color}" '
                         f'stroke-width="1.5" vector-effect="non-scaling-stroke"/><text x="{x+20}" y="{y+12}">{html.escape(label)}</text>')
    parts.append('</g></svg>')
    Path(path).write_text("\n".join(parts), encoding="utf-8", newline="\n")


def insert_new_molecules(
    reference_result, query_data, *, output_dir=None, input_bytes=None,
    input_filename="new_molecules.csv", smiles_col=None, id_col=None, value_col=None,
    value_unit="p", unit_column=None, relation_column=None, document_column=None,
    label_column=None, structure_policy=None, radius=None, n_bits=None,
    lower_threshold=None, upper_threshold=None, conflict_range_threshold=None,
    grouping=None, group_similarity=None, reference_target_id=None,
):
    """Insert one batch into a frozen Morgan reference map, with independent files.

    ``data`` contains placed molecules; ``all_molecules`` contains every checked
    compound, including duplicates, invalid structures and reference matches.
    Reference data, diagnostics and output bytes are never modified.
    """
    from .evidence import prepare_generic_activity_records, document_pair_offsets
    from .provenance import source_tree_digest, file_digests
    from .inspection import offset_pair_details
    started = time.perf_counter()
    settings = getattr(reference_result, "insertion_parameters", {})
    reference = reference_result.data
    if reference_target_id is None and "target_id" in reference:
        targets = reference.target_id.dropna().astype(str).str.strip().str.upper()
        targets = targets[targets.ne("")].unique()
        if len(targets) == 1:
            reference_target_id = targets[0]
    representation = settings.get("representation", reference.attrs.get("representation_type", "morgan"))
    if representation != "morgan":
        raise ValueError("New molecules require a Morgan fingerprint reference map.")
    defaults = dict(radius=2, n_bits=2048, lower_threshold=6., upper_threshold=7.,
                    conflict_range_threshold=1., grouping="butina", group_similarity=.6,
                    structure_policy="canonical")
    supplied = dict(radius=radius, n_bits=n_bits, lower_threshold=lower_threshold,
                    upper_threshold=upper_threshold, conflict_range_threshold=conflict_range_threshold,
                    grouping=grouping, group_similarity=group_similarity, structure_policy=structure_policy)
    for key, value in supplied.items():
        if key != "structure_policy" and value is not None and key in settings and value != settings[key]:
            raise ValueError(f"New molecule settings must match the reference map: {key}.")
    parameters = {key: settings.get(key, fallback) if supplied[key] is None else supplied[key]
                  for key, fallback in defaults.items()}
    policies = reference.structure_policy.dropna().unique() if "structure_policy" in reference else []
    reference_policy = str(policies[0]) if len(policies) == 1 else settings.get("structure_policy", "canonical")
    if parameters["structure_policy"] not in {"canonical", reference_policy}:
        raise ValueError("New molecule settings must match the reference map: structure_policy.")
    parameters["structure_policy"] = reference_policy
    radius, n_bits = int(parameters["radius"]), int(parameters["n_bits"])
    lower_threshold, upper_threshold = parameters["lower_threshold"], parameters["upper_threshold"]
    observed = reference.activity_class.notna() & ~reference.activity_class.astype(str).str.lower().isin(
        {"unlabelled", "not_measured", "", "nan", "none", "<na>"})
    if int(observed.sum()) < 2:
        raise ValueError("Query mapping requires at least two measured reference molecules.")
    if not observed.all():
        raise ValueError("Reference rows without measured activity must be supplied in query_data instead.")
    source = query_data.copy().reset_index(drop=True)
    if smiles_col is None:
        available = {str(c).strip().lower(): c for c in source.columns}
        smiles_col = next((available[c] for c in ("smiles", "canonical_smiles", "canonical smiles", "representation_smiles")
                           if c in available), None)
    if smiles_col is None or smiles_col not in source:
        raise ValueError("New molecule SMILES column was not found.")
    if id_col is None:
        id_col = "_new_molecule_id"
        source[id_col] = [f"N{i+1:05d}" for i in range(len(source))]
    if id_col not in source:
        raise ValueError("New molecule ID column was not found.")
    source[id_col] = source[id_col].fillna("").astype(str).str.strip()
    if source[id_col].eq("").any():
        raise ValueError("New molecule IDs must be nonempty.")
    compounds = source.drop_duplicates(id_col, keep="first").copy() if value_col else source.copy()
    # 同一测定 ID 的多条记录必须描述同一个结构；结构重复在化合物层面检查。
    if value_col:
        canonical = source[smiles_col].map(canonicalize_smiles)
        if parameters["structure_policy"] == "parent":
            from .structures import parent_structure
            canonical = canonical.map(lambda s: parent_structure(s)[0] if s is not None else None)
        if canonical.groupby(source[id_col]).nunique().gt(1).any():
            raise ValueError("Each new compound ID must map to one canonical structure.")
    queries, check, report = prepare_query_molecules(
        compounds, reference, smiles_column=smiles_col, id_column=id_col,
        structure_policy=parameters["structure_policy"])
    original_rows = dict(enumerate(compounds.index))
    check["source_row"] = check.source_row.map(original_rows)
    queries["source_row"] = queries.source_row.map(original_rows)
    parameters["structure_policy"] = report["structure_policy"]
    if len(queries) > 1000:
        raise ValueError("At most 1,000 checked new molecules are supported.")
    report.update(input_rows=len(source), input_compounds=len(compounds), measured=bool(value_col))
    all_molecules = check.copy()
    all_molecules["label"] = compounds[label_column].fillna("").astype(str).to_numpy() if label_column else ""
    all_molecules["is_new"] = True
    all_molecules["measurement_status"] = "measured" if value_col else "not_measured"
    all_molecules["activity_class"] = "not_measured"
    records = pd.DataFrame()
    if value_col:
        from .units import convert_activity_values
        _, reasons = convert_activity_values(source[value_col], value_unit=value_unit,
            units=source[unit_column] if unit_column else None)
        reasons.loc[reasons.eq("") & source[smiles_col].map(canonicalize_smiles).isna()] = "invalid_smiles"
        if relation_column:
            relation = source[relation_column].fillna("").astype(str).str.strip().str.replace(
                r"^'([=<>~]+)'$", r"\1", regex=True)
            reasons.loc[reasons.eq("") & relation.ne("=")] = "non_exact_relation"
        report["excluded_reasons"] = {reason: int(reasons.eq(reason).sum()) for reason in (
            "invalid_activity_value", "nonpositive_activity_value", "unrecognized_activity_unit",
            "invalid_smiles", "missing_compound_id", "non_exact_relation")}
        report["excluded_activity_records"] = int(reasons.ne("").sum())
        excluded = source.loc[reasons.ne("")].copy()
        excluded["source_row"] = excluded.index
        excluded["exclusion_reason"] = reasons[reasons.ne("")]
        report["excluded_records"] = json.loads(excluded.to_json(orient="records", double_precision=15))
        report["excluded_record_columns"] = {
            "compound_id": id_col, "evidence_unit": document_column, "value": value_col,
            "unit": unit_column, "relation": relation_column}
        report["excluded_record_defaults"] = {"unit": value_unit, "relation": "="}
        eligible = check[check.status.isin(["mapped", "in_reference"])]
        selected = source[source[id_col].isin(eligible.query_id)].copy()
        selected["_new_source_row"] = selected.index
        canonical_by_id = dict(zip(eligible.query_id, eligible.canonical_smiles))
        selected[smiles_col] = selected[id_col].map(canonical_by_id)
        if not selected.empty and reasons.loc[selected.index].eq("").any():
            aggregated, records, fit, activity_report = prepare_generic_activity_records(
                selected, compound_column=id_col, smiles_column=smiles_col, value_column=value_col,
                document_column=document_column, relation_column=relation_column,
                lower=lower_threshold, upper=upper_threshold,
                conflict_range_threshold=parameters["conflict_range_threshold"],
                value_unit=value_unit, unit_column=unit_column, confidence=bool(document_column),
                offset_display="none", offset_min_shared=3, offset_threshold=1., offset_sign_agreement=.8)
            report["activity_records"] = activity_report
            records["source_row"] = records.pop("_new_source_row")
            source_rows = records.groupby("molecule_identity_key")["source_row"].agg(
                lambda values: ";".join(str(value) for value in values))
            aggregated["source_rows"] = aggregated.molecule_identity_key.map(source_rows)
            if "raw_activity_value" not in records:
                records["raw_activity_value"] = records[value_col]
                records["raw_activity_unit"] = value_unit
            records["query_id"] = records[id_col].astype(str)
            records["record_relation"] = records[relation_column] if relation_column else "="
            records["evidence_unit"] = records[document_column].fillna("unknown_document") if document_column else ""
            if "n_evidence_units" not in aggregated:
                aggregated["n_evidence_units"] = 1
            annotations = aggregated.drop(columns=[c for c in ("canonical_smiles", "smiles", "activity_class")
                                                   if c in aggregated]).rename(columns={"molecule_id": "query_id"})
            all_molecules = all_molecules.merge(annotations, on="query_id", how="left", sort=False)
            class_by_id = dict(zip(aggregated.molecule_id, aggregated.activity_class))
            all_molecules["activity_class"] = all_molecules.query_id.map(class_by_id).fillna("not_measured")
            if document_column:
                offsets = document_pair_offsets(records)
                report["document_pair_offsets"] = offsets.to_dict(orient="records")
                flagged = offsets[offsets.systematic_offset].to_dict(orient="records")
                report["offset_pair_details"] = offset_pair_details(
                    records, "document_chembl_id", {"document_offsets": flagged, "offset_display": "none"}, bool(flagged))[0]
        missing = all_molecules.status.eq("mapped") & all_molecules.activity_class.eq("not_measured")
        if missing.any():
            # 保留有效结构，并明确说明该化合物没有可用的精确测定。
            all_molecules.loc[missing, "measurement_status"] = "not_measured"
            report["no_retained_activity_compounds"] = all_molecules.loc[missing, "query_id"].tolist()
    # 将聚合与标签按稳定 ID 挂回有效结构，不让输入中的 activity/evidence 混入未测结果。
    metadata = all_molecules.drop(columns=[c for c in ("source_row", "input_smiles", "canonical_smiles", "status",
                                                        "first_query_id", "reference_id", "molecule_identity_key") if c in all_molecules])
    queries = queries.drop(columns=[c for c in metadata if c != "query_id" and c in queries]).merge(
        metadata, on="query_id", how="left", sort=False)
    reference_fps = getattr(reference_result, "_insertion_fingerprints", None)
    # Morgan 表征已在建图时计算，这里仅转换其存储形式，不重新计算参照指纹。
    if reference_fps is None and np.shape(reference_result.representation) == (len(reference), n_bits):
        reference_fps = [DataStructs.CreateFromBitString("".join("1" if bit else "0" for bit in row))
                         for row in reference_result.representation]
    queries, neighbors = annotate_query_neighbors(reference, queries, radius=radius, n_bits=n_bits,
        value_col=settings.get("value_col"), reference_fingerprints=reference_fps)
    support_assessed = "class_status" in reference and reference.class_status.notna().any()
    queries = summarize_neighbor_evidence(queries, neighbors,
        similarity=parameters["group_similarity"], support_assessed=support_assessed)
    queries, groups = assign_query_groups(queries, reference_result.structure_groups,
        grouping=parameters["grouping"], similarity=parameters["group_similarity"], radius=radius, n_bits=n_bits)
    rows = int(settings.get("grid_rows", reference_result.metrics.iloc[0].grid_side))
    cols = int(settings.get("grid_cols", rows))
    queries = place_new_molecules(reference, queries, grid_rows=rows, grid_cols=cols)
    for group in groups:
        group["n_new_members"] = group["n_query_members"]
    report["reference_support_assessed"] = bool(support_assessed)
    report["neighbor_evidence_counts"] = {key: int(queries.neighbor_evidence.eq(key).sum()) for key in NEIGHBOR_EVIDENCE_CATEGORIES}
    report.update(_insertion_warnings(queries, query_data, reference_target_id))
    report["parameters"] = parameters
    report["projection_note"] = "New molecules are displayed only in the grid map."
    for column in queries:
        if column not in all_molecules:
            all_molecules[column] = all_molecules.query_id.map(queries.set_index("query_id")[column])
    all_molecules["user_median"] = all_molecules.get("activity_pchembl", np.nan)
    ref_value = settings.get("value_col", "activity_pchembl")
    reference_values = reference.get("activity_pchembl", reference.get(ref_value, pd.Series(np.nan, index=reference.index)))
    ref_lookup = dict(zip(molecule_keys(reference), reference_values))
    all_molecules["reference_median"] = all_molecules.reference_id.map(ref_lookup)
    all_molecules["median_difference"] = all_molecules.user_median - all_molecules.reference_median
    matches = all_molecules[all_molecules.status.eq("in_reference")]
    report["reference_matches"] = json.loads(matches[["query_id", "reference_id", "user_median", "reference_median", "median_difference"]].to_json(orient="records", double_precision=15))
    differences = matches.median_difference.dropna()
    report["reference_match_summary"] = {"count": len(differences), "median": float(differences.median()) if len(differences) else None,
                                         "min": float(differences.min()) if len(differences) else None,
                                         "max": float(differences.max()) if len(differences) else None}
    files = {}
    if output_dir is not None:
        base = Path(output_dir) / "new_molecules" / (Path(input_filename).stem or "new_molecules")
        destination, suffix = base, 2
        while destination.exists():
            destination = base.with_name(base.name + f"-{suffix}")
            suffix += 1
        destination.mkdir(parents=True)
        raw = bytes(input_bytes) if input_bytes is not None else query_data.to_csv(index=False, lineterminator="\n").encode("utf-8")
        path = destination / "new_molecules_input.csv"
        path.write_bytes(raw)
        files["new_molecules_input"] = str(path)
        for key, table in (("new_molecules_check", check), ("new_molecules", all_molecules), ("new_molecules_neighbors", neighbors)):
            path = destination / f"{key}.csv"
            table.to_csv(path, index=False, encoding="utf-8", lineterminator="\n")
            files[key] = str(path)
        if value_col:
            # 保留输入原列与原值，不导出标准化 ID 或内部 source_row。
            excluded_input = query_data.reset_index(drop=True).loc[reasons.ne("")].copy()
            excluded_input["exclusion_reason"] = reasons.loc[reasons.ne("")]
            path = destination / "new_molecules_excluded_records.csv"
            excluded_input.to_csv(path, index=False, encoding="utf-8", lineterminator="\n")
            files["new_molecules_excluded_records"] = str(path)
        svg = destination / "map_with_new_molecules.svg"
        _new_molecules_svg(reference_result, queries, svg, rows=rows, cols=cols,
                           lower=lower_threshold, upper=upper_threshold)
        files["map_with_new_molecules"] = str(svg)
        reference_manifest = reference_result.output_files.get("run_manifest")
        manifest = {"input_filename": Path(input_filename).name, "input_sha256": hashlib.sha256(raw).hexdigest(),
                    "reference_run_directory": str(Path(output_dir).resolve()),
                    "reference_manifest_sha256": hashlib.sha256(Path(reference_manifest).read_bytes()).hexdigest() if reference_manifest else None,
                    "parameters": {**parameters, "smiles_col": smiles_col, "id_col": id_col, "value_col": value_col,
                                   "value_unit": value_unit, "unit_column": unit_column, "relation_column": relation_column,
                                   "document_column": document_column, "label_column": label_column},
                    "report": report, "software_source_sha256": source_tree_digest()["source_sha256"],
                    **{key: report[key] for key in ("warnings", "none_share", "median_nearest_tanimoto")},
                    "output_digests": file_digests(files),
                    "timings_seconds": {"through_svg_and_tables": time.perf_counter()-started},
                    "timing_note": "The returned report.insertion_seconds includes manifest and ZIP writing."}
        path = destination / "new_molecules_manifest.json"
        # pandas JSON conversion replaces numerical missing values with JSON null.
        manifest = json.loads(pd.Series(manifest).to_json(double_precision=15))
        path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8", newline="\n")
        files["new_molecules_manifest"] = str(path)
        archive = destination / "new_molecules_results.zip"
        with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as bundle:
            for filename in files.values():
                bundle.write(filename, Path(filename).name)
        files["new_molecules_results_zip"] = str(archive)
    report["insertion_seconds"] = time.perf_counter() - started
    report = json.loads(pd.Series(report).to_json(double_precision=15))
    return NewMoleculesResult(queries, check, neighbors, report, files, records, groups, all_molecules)
