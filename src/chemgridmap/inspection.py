"""Recover the measurements behind a molecular grid cell without remapping it."""
from __future__ import annotations

import argparse
import csv
from functools import lru_cache
import gzip
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import quote, urlsplit

import numpy as np
import pandas as pd

from .chembl import detect_chembl_columns
from .core import canonicalize_smiles
from .structures import parent_structure
from .plotting import render_svg
from .activity import RULE_KEYS, RULE_COLUMNS, activity_annotations
from .provenance import file_digests


@lru_cache(maxsize=1)
def load_documents():
    """Read the bundled ChEMBL 37 document table, keyed by document ID.

    Values retain the CSV strings, including empty strings. The cached table is
    shared by callers and must be treated as read-only. No network is used.
    """
    path = Path(__file__).with_name("data") / "chembl37_documents.csv.gz"
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return {row["document_chembl_id"]: row for row in csv.DictReader(handle)}


def _reference_text(value):
    if value is None or pd.isna(value):
        return ""
    if isinstance(value, (float, np.floating)) and float(value).is_integer():
        return str(int(value))
    return str(value)


def normalize_doi(value):
    """Remove supported DOI prefixes and surrounding whitespace, preserving case."""
    text = _reference_text(value).strip()
    return re.sub(r"^(?:doi:\s*|https?://(?:dx\.)?doi\.org/)", "", text,
                  flags=re.IGNORECASE).strip()


def reference_link(value):
    """Resolve a DOI, PMID or HTTP(S) reference into a link, or return ``None``.

    The returned ``kind``, ``label`` and ``url`` are plain strings. Renderers
    must escape labels and open links only after an explicit user click.
    """
    text = _reference_text(value).strip()
    doi = normalize_doi(text)
    if re.fullmatch(r"10\.\d{4,9}/.+", doi, flags=re.DOTALL):
        return {"kind": "doi", "label": "DOI " + doi,
                "url": "https://doi.org/" + quote(doi, safe="/")}
    pmid = (re.fullmatch(r"PMID:\s*(\d+)", text, flags=re.IGNORECASE)
            or re.fullmatch(r"(\d{1,9})", text))
    if pmid:
        identifier = pmid.group(1)
        return {"kind": "pubmed", "label": "PubMed " + identifier,
                "url": "https://pubmed.ncbi.nlm.nih.gov/" + identifier + "/"}
    try:
        parsed = urlsplit(text)
    except ValueError:
        return None
    if parsed.scheme.lower() in {"http", "https"} and parsed.netloc:
        return {"kind": "url", "label": text, "url": text}
    return None


def _document_links(document):
    links = []
    if document.get("doi"):
        link = reference_link(document["doi"])
        if link:
            links.append(link)
    if document.get("pubmed_id"):
        link = reference_link("PMID:" + document["pubmed_id"])
        if link:
            links.append(link)
    if document.get("patent_id"):
        patent = document["patent_id"]
        links.append({"kind": "patent", "label": "Patent " + patent,
                      "url": "https://patents.google.com/patent/" +
                             quote(patent.replace("-", ""), safe="")})
    identifier = document["document_chembl_id"]
    links.append({"kind": "chembl", "label": "ChEMBL " + identifier,
                  "url": "https://www.ebi.ac.uk/chembl/explore/document/" + identifier})
    return links


def build_sources(records, document_column, reference_column=None):
    """Group retained records into document/reference sources without changing them.

    Evidence-unit IDs take precedence for grouping. Without an evidence column,
    nonempty reference values form the groups. References remain visible even
    when their protocol cannot become a link. Results sort by year then ID,
    with unknown years last; all records and distinct references are retained.
    """
    group_column = document_column or reference_column
    if not group_column or group_column not in records:
        return []
    documents = load_documents()
    groups = records.groupby(records[group_column].map(_reference_text).str.strip(), sort=True)
    sources = []
    for identifier, group in groups:
        if not identifier:
            continue
        document = documents.get(identifier, {}) if document_column else {}
        links = _document_links(document) if document else []
        references = []
        seen_urls = {link["url"] for link in links}
        if reference_column and reference_column in group:
            values = dict.fromkeys(group[reference_column].map(_reference_text))
            for value in values:
                if not value.strip():
                    continue
                link = reference_link(value)
                references.append({"value": value, **(link or {
                    "kind": "text", "label": value, "url": ""})})
                if link and link["url"] not in seen_urls:
                    links.append(link)
                    seen_urls.add(link["url"])
        sources.append({"id": identifier, "doc_type": document.get("doc_type", ""),
                        "journal": document.get("journal", ""),
                        "year": document.get("year", ""), "n_records": len(group),
                        "links": links, "references": references})
    return sorted(sources, key=lambda source: (
        int(source["year"]) if source["year"].isdigit() else float("inf"), source["id"]))


SOURCE_LOOKUP_COLUMNS_ERROR = (
    "The source lookup file needs SMILES and ChEMBL document ID columns, "
    "as in a ChEMBL activity export."
)


def validate_source_lookup(frame):
    """Detect supported ChEMBL columns and require structure/document fields."""
    columns = detect_chembl_columns(frame)
    if not columns["smiles"] or not columns["document_id"]:
        raise ValueError(SOURCE_LOOKUP_COLUMNS_ERROR)
    return columns


def match_source_lookup(molecules, retained_records, lookup_frame):
    """Return possible documents per molecule and counts for the session audit.

    Parent normalization is shared across both inputs and evaluated once per
    distinct SMILES. Neither molecular IDs, targets nor activity types filter
    the match. Values only mark documents; they never select or reject sources.
    """
    columns = validate_source_lookup(lookup_frame)
    smiles_cache = {}
    parent_cache = {}

    def normalized_parent(value):
        raw = "" if pd.isna(value) else str(value).strip()
        if raw not in smiles_cache:
            canonical = canonicalize_smiles(raw)
            if canonical and canonical not in parent_cache:
                try:
                    parent_cache[canonical] = parent_structure(canonical)[0]
                except (ValueError, RuntimeError):
                    parent_cache[canonical] = None
            smiles_cache[raw] = parent_cache.get(canonical)
        return smiles_cache[raw]

    document_ids = lookup_frame[columns["document_id"]].fillna("").astype(str).str.strip()
    valid_ids = document_ids.str.fullmatch(r"CHEMBL\d+")
    lookup = pd.DataFrame({
        "document_id": document_ids[valid_ids],
        "parent": lookup_frame.loc[valid_ids, columns["smiles"]].map(normalized_parent),
    })
    lookup["pchembl"] = (pd.to_numeric(lookup_frame.loc[valid_ids, columns["pchembl"]], errors="coerce")
                         if columns["pchembl"] else np.nan)
    lookup = lookup.dropna(subset=["parent"])

    # Resolve the bundled metadata once per document, then reuse the safe links.
    documents = {source["id"]: source for source in build_sources(lookup, "document_id")}
    for identifier, source in documents.items():
        if not any(link["kind"] == "chembl" for link in source["links"]):
            source["links"].append({
                "kind": "chembl", "label": "ChEMBL " + identifier,
                "url": "https://www.ebi.ac.uk/chembl/explore/document/" + identifier,
            })
    by_parent = {}
    for (parent, identifier), group in lookup.groupby(["parent", "document_id"], sort=False):
        values = group["pchembl"].to_numpy(float)
        by_parent.setdefault(parent, []).append((
            identifier, len(group), np.unique(values[np.isfinite(values)]),
        ))

    value_column = ("activity_pchembl_raw" if "activity_pchembl_raw" in retained_records
                    else "pchembl_value")
    measured_values = {
        str(identity): pd.to_numeric(group[value_column], errors="coerce").to_numpy(float)
        for identity, group in retained_records.groupby("molecule_identity_key", sort=False)
    }
    by_molecule_id = {}
    for row in molecules.itertuples(index=False):
        identity = str(row.molecule_identity_key)
        parent = normalized_parent(row.canonical_smiles)
        values = measured_values.get(identity, np.array([], dtype=float))
        values = values[np.isfinite(values)]
        sources = []
        for identifier, count, reported in by_parent.get(parent, []):
            # A tiny float allowance keeps the inclusive 0.005 boundary exact.
            same = reported[(np.abs(reported[:, None] - values) <= .005 + 1e-12).any(axis=1)]
            source = {**documents[identifier], "n_records": count,
                      "same_values": sorted({f"{value:.2f}" for value in same}, key=float)}
            sources.append(source)
        by_molecule_id[identity] = sorted(sources, key=lambda source: (
            not bool(source["same_values"]),
            int(source["year"]) if source["year"].isdigit() else float("inf"),
            source["id"],
        ))

    targets = (lookup_frame[columns["target_id"]].dropna().astype(str).str.strip()
               if columns["target_id"] else pd.Series(dtype=str))
    summary = {
        "found": sum(bool(sources) for sources in by_molecule_id.values()),
        "total": len(molecules),
        "same": sum(any(source["same_values"] for source in sources)
                    for sources in by_molecule_id.values()),
        "target_count": int(targets[targets.ne("")].nunique()),
    }
    return by_molecule_id, summary



def document_coverage(sources):
    """Count linked sources and mutually exclusive ChEMBL coverage categories.

    Patents take precedence over DOI/PubMed sources to avoid double-counting the
    few ChEMBL patents with publication identifiers. ``linked`` also includes
    safe user-supplied URLs; ``unlinked`` preserves evidence units without links.
    """
    coverage = {"total": len(sources), "doi_or_pubmed": 0, "patent": 0,
                "chembl_only": 0, "linked": 0, "unlinked": 0}
    for source in sources:
        kinds = {link["kind"] for link in source["links"]}
        coverage["linked" if kinds else "unlinked"] += 1
        if "patent" in kinds:
            coverage["patent"] += 1
        elif kinds & {"doi", "pubmed"}:
            coverage["doi_or_pubmed"] += 1
        elif "chembl" in kinds:
            coverage["chembl_only"] += 1
    return coverage


def enrich_source_records(records, document_column=None):
    """Append four ChEMBL document fields to a copy of inspector source records.

    Existing columns and their order are preserved. These fields are only for
    inspector downloads, never for map calculation or normal run output files.
    An input field with the same name remains in place, yielding a repeated CSV heading.
    """
    additions = pd.DataFrame(index=records.index)
    if document_column is None and "document_chembl_id" in records:
        document_column = "document_chembl_id"
    documents = load_documents()
    identifiers = (records[document_column].map(_reference_text).str.strip()
                   if document_column and document_column in records
                   else pd.Series("", index=records.index))
    selected = [documents.get(identifier, {}) for identifier in identifiers]
    for output, source in (("document_doi", "doi"), ("document_pubmed_id", "pubmed_id"),
                           ("document_patent_id", "patent_id")):
        additions[output] = [document.get(source, "") for document in selected]
    additions["document_url"] = [
        "https://www.ebi.ac.uk/chembl/explore/document/" + document["document_chembl_id"]
        if document else "" for document in selected]
    # CSV permits repeated headings: even existing same-named input fields stay untouched.
    return pd.concat([records, additions], axis=1)


def annotation_rules(row, parameters=None):
    present = [column in row.index for column in RULE_COLUMNS]
    if any(present) and not all(present):
        raise ValueError("Incomplete saved annotation rules; supply an intact coordinate table.")
    embedded = {key: float(row[column]) for key, column in zip(RULE_KEYS, RULE_COLUMNS)} if all(present) else None
    if parameters is None:
        return embedded
    if not set(RULE_KEYS).issubset(parameters):
        raise ValueError("Curation parameters lack activity thresholds or the conflict range.")
    explicit = {key: float(parameters[key]) for key in RULE_KEYS}
    if embedded is not None and explicit != embedded:
        raise ValueError("Saved annotation rules differ from the supplied curation parameters.")
    return explicit


def verify_annotations(row, evidence, values, parameters):
    rules = annotation_rules(row, parameters)
    calculated = {
        "pchembl_min": float(values.min()), "pchembl_max": float(values.max()),
        "pchembl_range": float(values.max() - values.min()), "pchembl_median": float(values.median()),
        "pchembl_std": float(values.std(ddof=1)),
    }
    if rules is not None:
        calculated.update(activity_annotations(values.to_numpy(), **rules))
    verified = []
    for key, value in calculated.items():
        if key not in row.index:
            continue
        try:
            agrees = str(row[key]) == value if isinstance(value, str) else bool(
                np.isclose(float(row[key]), value, rtol=0, atol=1e-10, equal_nan=True))
        except (ValueError, TypeError):
            agrees = False
        if not agrees:
            raise ValueError(f"Recalculated annotation {key} does not match the map.")
        verified.append(key)
    # Identity links also bind the displayed structure to the retained representation.
    if "representation_smiles" in evidence:
        expected = sorted(evidence.representation_smiles.astype(str).unique())[0]
        for key in ("canonical_smiles", "representation_smiles", "smiles"):
            if key in row.index and str(row[key]) != expected:
                raise ValueError(f"Displayed structure {key} does not match the retained records.")
            if key in row.index:
                verified.append(key)
    columns = detect_chembl_columns(evidence)
    if "n_assays" in row.index and columns["assay_id"]:
        if float(row.n_assays) != evidence[columns["assay_id"]].nunique():
            raise ValueError("Source assay count does not match the retained records.")
        verified.append("n_assays")
    for key in ("molecule_id", "parent_molecule_id", "assay_id", "document_id", "activity_id", "target_id"):
        field = columns[key]
        if key in row.index and field:
            expected = {str(value).strip() for value in evidence[field].dropna() if str(value).strip()}
            actual = set() if pd.isna(row[key]) else {value.strip() for value in str(row[key]).split(";") if value.strip()}
            if actual != expected:
                raise ValueError(f"Source identifier {key} does not match the retained records.")
            verified.append(key)
    required_flags = {"activity_class", "is_conflicted", "has_class_boundary_crossing"}
    quality_verified = rules is not None and required_flags.issubset(verified)
    return calculated, rules, verified, quality_verified


def offset_pair_details(evidence, document_column, confidence, touched):
    """Describe flagged pairs involving any evidence unit listed for this entry."""
    documents = set()
    if document_column:
        documents = {
            str(value) for value in evidence[document_column].dropna().astype(str)
            if str(value)
        }
    matched_offsets = []
    offset_warnings = []
    if touched:
        offset_display = str(confidence.get("offset_display") or "uncertain")
        for pair in confidence.get("document_offsets", []):
            document_a, document_b = str(pair["document_a"]), str(pair["document_b"])
            present = [document for document in (document_a, document_b) if document in documents]
            if not present:
                continue
            matched_offsets.append(pair)
            touched_count = pair.get("touched_entry_count")
            touches = (
                f"the pair touches {int(touched_count)} entries"
                if touched_count is not None else "the pair touches an unrecorded number of entries"
            )
            record_text = ", ".join(present)
            if offset_display == "uncertain":
                display_rule = (
                    f"This entry has records from {record_text}, so it is shown as uncertain."
                )
            else:
                display_rule = (
                    f"This entry has records from {record_text}; "
                    "the display is unchanged because the offset rule is off."
                )
            offset_warnings.append(
                "Offset document pair: {} vs {}, median difference {:.3g} log units, "
                "n shared {}; {}. {}".format(
                    document_a, document_b, float(pair["median_offset"]),
                    int(pair["shared_entries"]), touches, display_rule,
                )
            )
        if not matched_offsets:
            offset_warnings.append("Offset document pair details are unavailable for this map.")
    return matched_offsets, offset_warnings


def inspect_entry(coordinates, records, *, molecule_id=None, cell=None, parameters=None,
                  activity_confidence=None, k=None):
    """Select one entry by exact ID or zero-based (row, column), then verify its records.

    Parameters
    ----------
    coordinates : pandas.DataFrame
        Grid-coordinate table from one completed run.
    records : pandas.DataFrame
        Retained source-record table from the same run.
    molecule_id : str, optional
        Exact molecule identity or source identifier to select.
    cell : tuple of int, optional
        Zero-based ``(row, column)`` cell to select.
    parameters : dict, optional
        Curation annotation thresholds for verification.
    activity_confidence : dict, optional
        Activity-confidence fit from the run manifest or curation report.
    k : int, optional
        Neighbour count used to label the kth-neighbour diagnostics.

    Returns
    -------
    selected : pandas.DataFrame
        One-row coordinate subset for the selected entry.
    evidence : pandas.DataFrame
        Source records linked to the selected map row.
    summary : dict
        Verified annotations, confidence diagnostics, offset pairs and warnings.
    """
    if (molecule_id is None) == (cell is None):
        raise ValueError("Select either molecule_id or cell=(row, column).")
    if cell is not None:
        selected = coordinates[coordinates.grid_row.eq(cell[0]) & coordinates.grid_col.eq(cell[1])]
    else:
        wanted = str(molecule_id).strip()
        is_chembl_id = re.fullmatch(r"CHEMBL\d+", wanted, flags=re.IGNORECASE) is not None
        if is_chembl_id:
            wanted = wanted.upper()
        matches = pd.Series(False, index=coordinates.index)
        for field in ("molecule_identity_key", "molecule_id", "parent_molecule_id"):
            if field in coordinates:
                # SMILES keys are case-sensitive: C (aliphatic) and c (aromatic) differ.
                matches |= coordinates[field].fillna("").astype(str).map(
                    lambda value: wanted in {v.strip().upper() if is_chembl_id else v.strip()
                                             for v in value.split(";")}
                )
        selected = coordinates[matches]
    if len(selected) != 1:
        raise ValueError(f"Selection matched {len(selected)} entries; one occupied cell is required.")
    row = selected.iloc[0]
    required = {"source_rows", "n_records", "activity_pchembl", "molecule_identity_key"}
    if not required.issubset(selected.columns):
        raise ValueError("This map lacks ChEMBL record provenance; use ChEMBL input mode first.")
    source_rows = [int(v) for v in str(row.source_rows).split(";")]
    if "source_row" not in records or records.source_row.duplicated().any():
        raise ValueError("Use the retained-record table from the same run, with unique source_row values.")
    evidence = records[records.source_row.isin(source_rows)].copy().sort_values("source_row")
    if len(evidence) != len(source_rows) or len(evidence) != int(row.n_records):
        raise ValueError("Record counts do not match this cell; check the retained-record file.")
    if "molecule_identity_key" not in evidence or not evidence.molecule_identity_key.eq(row.molecule_identity_key).all():
        raise ValueError("The records belong to a different molecular identity.")
    identity_records = records[records.molecule_identity_key.eq(row.molecule_identity_key)]
    if set(identity_records.source_row) != set(source_rows):
        raise ValueError("Source links omit records belonging to this molecular identity.")
    values = pd.to_numeric(evidence.activity_pchembl_raw, errors="raise")
    median = float(values.median())
    if not np.isfinite(values).all() or not np.isclose(median, float(row.activity_pchembl), rtol=0, atol=1e-10):
        raise ValueError("Recalculated activity does not match the map annotation.")
    columns = detect_chembl_columns(evidence)
    calculated, rules, verified, quality_verified = verify_annotations(row, evidence, values, parameters)
    if "offset_pair_touched" not in row.index or pd.isna(row["offset_pair_touched"]):
        touched = None
        status_reason = None
    else:
        touched = bool(row["offset_pair_touched"])
        status_reason = None if pd.isna(row.get("status_reason")) else row.get("status_reason")
    confidence = activity_confidence or {}
    document_column = columns["document_id"]
    matched_offsets, offset_warnings = offset_pair_details(
        evidence, document_column, confidence, touched)
    summary = {
        "molecule_identity_key": str(row.molecule_identity_key),
        "grid_row": int(row.grid_row), "grid_col": int(row.grid_col),
        "coordinate_convention": "Zero-based row and column from the grid-coordinate CSV; rows increase upwards.",
        "activity_class": calculated.get("activity_class"), "median_pchembl": median,
        "minimum_pchembl": float(values.min()), "maximum_pchembl": float(values.max()),
        "n_records": len(evidence),
        "is_conflicted": bool(calculated["is_conflicted"]) if rules is not None else None,
        "wide_record_range": bool(calculated["is_conflicted"]) if rules is not None else None,
        "has_class_boundary_crossing": bool(calculated["has_class_boundary_crossing"]) if rules is not None else None,
        "offset_pair_touched": touched,
        "status_reason": status_reason,
        "document_offsets": matched_offsets,
        "confidence_threshold": confidence.get("confidence_threshold"),
        "offset_display": confidence.get("offset_display"),
        "offset_min_shared": confidence.get("offset_min_shared"),
        "offset_threshold": confidence.get("offset_threshold"),
        "offset_sign_agreement": confidence.get("offset_sign_agreement"),
        "sigma_doc": confidence.get("sigma_doc", confidence.get("tau")),
        "sigma_doc_degrees_of_freedom": confidence.get(
            "sigma_doc_degrees_of_freedom", confidence.get("degrees_of_freedom")
        ),
        "multi_unit_entry_count": confidence.get(
            "multi_unit_entry_count", confidence.get("n_entries")
        ),
        "record_count_verified": True, "median_verified": True,
        "quality_annotations_verified": quality_verified, "annotation_rules": rules,
        "verified_fields": verified,
        "warnings": offset_warnings + (
            [] if quality_verified else [
                "Quality annotations are not fully verified. Legacy maps need the original curation parameters; no default thresholds were assumed."
            ]
        ),
        "interpretation": (
            "Median annotation summarizes retained records; it does not harmonize assay conditions. "
            "Class support accounts for agreement among documents."
        ),
    }
    for name, field in (("n_activity_ids", "activity_id"), ("n_assays", "assay_id"), ("n_documents", "document_id")):
        summary[name] = int(evidence[columns[field]].nunique()) if columns[field] else None
    neighbor_k = f"{int(k)}th" if k else "kth"
    if "neighbor_kth_similarity" in row.index and pd.notna(row["neighbor_kth_similarity"]):
        summary["neighbor_kth_similarity"] = float(row["neighbor_kth_similarity"])
        summary["neighbor_kth_label"] = (
            f"{neighbor_k} fingerprint neighbor: Tanimoto "
            f"{float(row['neighbor_kth_similarity']):.3f}"
        )
    elif "neighbor_kth_distance" in row.index and pd.notna(row["neighbor_kth_distance"]):
        summary["neighbor_kth_distance"] = float(row["neighbor_kth_distance"])
        summary["neighbor_kth_label"] = (
            f"{neighbor_k} representation neighbor: Euclidean distance "
            f"{float(row['neighbor_kth_distance']):.3f}"
        )
    return selected.copy(), evidence, summary


def load_audit_context(coordinates_path, records_path, *, manifest_path=None, report_path=None,
                       return_manifest=False):
    """Check saved file digests before interpreting them as outputs of the same run."""
    explicit_manifest = manifest_path is not None
    manifest_path = Path(manifest_path) if manifest_path else coordinates_path.with_name(
        coordinates_path.name.replace("_grid_coordinates.csv", "_run_manifest.json"))
    parameters, verified, warnings = None, False, []
    run_manifest = {}
    if manifest_path != coordinates_path and manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        run_manifest = manifest
        parameters = manifest.get("source", {}).get("curation_parameters")
        expected = {"coordinates": manifest.get("output_digests", {}).get("coordinates"),
                    "retained_records": manifest.get("source", {}).get("curation_outputs", {}).get("retained_records")}
        actual = file_digests({"coordinates": coordinates_path, "retained_records": records_path})
        for key, digest in expected.items():
            if digest and actual[key]["sha256"] != digest["sha256"]:
                raise ValueError(f"{key} digest differs from the run manifest; do not mix or modify run outputs.")
        verified = all(expected.values())
    elif explicit_manifest:
        raise ValueError("The supplied run manifest does not exist.")
    if report_path is not None:
        report = json.loads(Path(report_path).read_text(encoding="utf-8"))
        supplied = report["parameters"]
        if parameters is not None and any(parameters.get(key) != supplied.get(key) for key in RULE_KEYS):
            raise ValueError("Curation report and run manifest use different annotation rules.")
        parameters = supplied
    if not verified:
        warnings.append("Run binding is not verified: original output digests are unavailable. Current file hashes alone do not establish a shared run.")
    if return_manifest:
        return parameters, verified, warnings, run_manifest
    return parameters, verified, warnings


def main(argv=None):
    parser = argparse.ArgumentParser(prog="chemgridmap inspect", description=__doc__)
    parser.add_argument("coordinates", type=Path, help="Grid-coordinate CSV from a ChEMBL-to-map run.")
    parser.add_argument("--records", required=True, type=Path, help="Retained-record CSV from the same run.")
    parser.add_argument("--manifest", type=Path, help="Run manifest; automatically found beside standard coordinate exports.")
    parser.add_argument("--curation-report", type=Path, help="Original curation report for legacy maps without saved annotation rules.")
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--molecule-id", help="Exact molecule, parent, or molecular identity key.")
    selection.add_argument("--cell", nargs=2, type=int, metavar=("ROW", "COL"), help="Zero-based row and column in the grid CSV.")
    parser.add_argument("--output-dir", type=Path, default=Path("record_audit"))
    args = parser.parse_args(argv)
    parameters, run_verified, warnings, manifest = load_audit_context(
        args.coordinates, args.records, manifest_path=args.manifest,
        report_path=args.curation_report, return_manifest=True)
    selected, records, summary = inspect_entry(
        pd.read_csv(args.coordinates, encoding="utf-8-sig"),
        pd.read_csv(args.records, encoding="utf-8-sig"),
        molecule_id=args.molecule_id, cell=args.cell, parameters=parameters,
        activity_confidence=manifest.get("activity_confidence"),
        k=manifest.get("parameters", {}).get("k"),
    )
    summary["run_binding_verified"] = run_verified
    summary["warnings"].extend(warnings)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.output_dir / "selected_molecule.csv", index=False,
                    encoding="utf-8", lineterminator="\n")
    records.to_csv(args.output_dir / "source_records.csv", index=False,
                   encoding="utf-8", lineterminator="\n")
    render_svg(selected, args.output_dir / "selected_molecule.svg", tile_size=320, molecule_margin=12)
    summary["inputs"] = {name: {"file": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                         for name, path in (("coordinates", args.coordinates), ("retained_records", args.records))}
    (args.output_dir / "audit_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(summary, indent=2))
    print(f"Evidence written to {args.output_dir}")
    return 0
