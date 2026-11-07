"""Local browser interface. Scientific operations use the same API as the CLI."""
from __future__ import annotations

import argparse
import csv
import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import importlib.util
import io
import json
import math
import mimetypes
from pathlib import Path
import re
import secrets
import threading
import time
from urllib.parse import parse_qs, unquote, urlsplit
import webbrowser
import zipfile

import pandas as pd

from . import __version__
from .chembl import curate_chembl_activity_data, detect_chembl_columns, save_chembl_curation
from .evidence import annotate_activity_confidence, prepare_generic_activity_records
from .inspection import (build_sources, document_coverage, enrich_source_records,
                         inspect_entry, load_audit_context, offset_pair_details,
                         match_source_lookup, validate_source_lookup)
from .pipeline import build_grid_map
from .plotting import render_svg
from .provenance import file_digests
from .neighborhoods import neighborhood_evidence
from .comparison import group_comparison, molecule_keys
from .plotting import _svg_molecule_content
from .units import validate_assumed_sigma

ASSETS = Path(__file__).with_name("web")
DEMO = Path(__file__).with_name("data") / "chembl205_ic50.csv"
MAX_UPLOAD = 50 * 1024 * 1024


def json_rows(frame):
    return json.loads(frame.to_json(orient="records"))


POINT_COLUMNS = (
    "projection_x", "projection_y", "grid_x", "grid_y", "grid_row", "grid_col",
    "activity_class", "molecule_id", "molecule_identity_key", "structure_group",
    "structure_group_color", "comparison_id", "neighbor_recall",
    "neighbor_recall_fixed_k", "class_confidence", "class_status", "status_reason",
    "offset_pair_touched", "wide_record_range", "is_conflicted",
    "neighbor_kth_similarity", "neighbor_kth_distance", "n_evidence_units",
)


def point_rows(data):
    columns = [name for name in POINT_COLUMNS if name in data]
    if "measurement_status" in data:
        columns.append("measurement_status")
    return json_rows(data[columns])


def _read_input(content, *, dtype=None, keep_default_na=True):
    if not content or len(content) > MAX_UPLOAD:
        raise ValueError("Choose a non-empty CSV smaller than 50 MB.")
    try:
        sample = content.decode("utf-8-sig", errors="strict")[:65536]
        try:
            delimiter = csv.Sniffer().sniff(sample, delimiters=",;\t|").delimiter
        except csv.Error:
            delimiter = ","
        frame = pd.read_csv(io.BytesIO(content), sep=delimiter, encoding="utf-8-sig", dtype=dtype,
                            keep_default_na=keep_default_na)
    except (UnicodeError, pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
        raise ValueError("Could not read this UTF-8 CSV. Check its delimiter and encoding.") from exc
    if frame.empty:
        raise ValueError("The CSV has no data rows.")
    return frame


def _select_molecule(data, row, col, molecule_id):
    if molecule_id is None:
        chosen = data[data.grid_row.eq(row) & data.grid_col.eq(col)]
    else:
        wanted = str(molecule_id).strip()
        chembl_id = re.fullmatch(r"CHEMBL\d+", wanted, flags=re.IGNORECASE) is not None
        if chembl_id:
            wanted = wanted.upper()
        matches = pd.Series(False, index=data.index)
        for field in ("molecule_identity_key", "molecule_id", "parent_molecule_id", "query_id", "name",
                      "canonical_smiles"):
            if field in data:
                matches |= data[field].fillna("").astype(str).map(
                    lambda value: wanted in {v.strip().upper() if chembl_id else v.strip()
                                             for v in value.split(";")})
        chosen = data[matches]
    if len(chosen) != 1:
        raise ValueError(f"Selection matched {len(chosen)} entries; one occupied cell is required.")
    return chosen


def _query_selection(check, molecule_id):
    """Resolve uploaded IDs using the existing structure check, without rematching."""
    if molecule_id is None or check.empty:
        return molecule_id, None
    matches = check[check.query_id.eq(str(molecule_id).strip())]
    if matches.empty:
        return molecule_id, None
    match = matches.iloc[0][["query_id", "status", "reference_id", "first_query_id"]].to_dict()
    if match["status"] == "mapped":
        return molecule_id, None
    if match["status"] == "invalid_smiles":
        raise ValueError(f"Unmeasured molecule {match['query_id']} has invalid SMILES.")
    target = match["reference_id"] or "query:" + (match["first_query_id"] or match["query_id"])
    return target, match


def _group_drawing(summary, field):
    molecule = summary.get("core_smiles_attached") if field == 'core_smiles' else summary[field]
    molecule = molecule or summary[field]
    return ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 180 180">'
            + (_svg_molecule_content(molecule, 180) or '') + '</svg>') if molecule else ''


class AppSession:
    """One local workspace; a single worker avoids concurrent chemistry/render jobs."""

    def __init__(self, output_root):
        self.output_root = Path(output_root).expanduser().resolve()
        self.token = secrets.token_urlsafe(32)
        self.lock = threading.RLock()
        self.frame = None
        self.raw = b""
        self.filename = ""
        self.query_frame = None
        self.query_raw = b""
        self.query_filename = ""
        self.source_frame = None
        self.source_filename = ""
        self.possible_sources = None
        self.result = None
        self.curation = None
        self.generic_records = None
        self.activity_confidence = None
        self.files = {}
        self.group_cache = {}
        self.group_detail_cache = {}
        self.job = {"state": "empty", "message": "Choose a CSV or open the example."}

    def load(self, content, filename):
        with self.lock:
            if self.job["state"] == "running":
                raise ValueError("A map is running. Wait for it to finish before replacing the input.")
            frame = _read_input(content)
            self.frame, self.raw = frame, content
            self.filename = Path(filename.replace("\\", "/")).name
            self.query_frame, self.query_raw, self.query_filename = None, b"", ""
            self.source_frame, self.source_filename, self.possible_sources = None, "", None
            self.result, self.curation, self.files = None, None, {}
            self.generic_records = None
            self.activity_confidence = None
            self.group_cache = {}
            self.group_detail_cache = {}
            self.job = {"state": "ready", "message": "Input loaded. Review the options, then build the map."}
            return self.input_info()

    def load_query(self, content, filename):
        with self.lock:
            if self.job["state"] == "running":
                raise ValueError("A map is running. Wait for it to finish before replacing the input.")
            if self.frame is None:
                raise ValueError("Load the measured data first.")
            # Query IDs remain text, including leading zeroes.
            frame = _read_input(content, dtype=str)
            self.query_frame, self.query_raw = frame, content
            self.query_filename = Path(filename.replace("\\", "/")).name
            self._discard_map()
            return self.query_input_info()

    def clear_query(self):
        with self.lock:
            if self.job["state"] == "running":
                raise ValueError("A map is running. Wait for it to finish before replacing the input.")
            self.query_frame, self.query_raw, self.query_filename = None, b"", ""
            self._discard_map()
            return {"query_input": None}

    def load_source(self, content, filename):
        with self.lock:
            if self.job["state"] == "running":
                raise ValueError("A map is running. Wait for it to finish before replacing the input.")
            if self.frame is None:
                raise ValueError("Load the measured data first.")
            frame = _read_input(content, dtype=str)
            validate_source_lookup(frame)
            self.source_frame = frame
            self.source_filename = Path(filename.replace("\\", "/")).name
            self._discard_map()
            return self.source_input_info()

    def clear_source(self):
        with self.lock:
            if self.job["state"] == "running":
                raise ValueError("A map is running. Wait for it to finish before replacing the input.")
            self.source_frame, self.source_filename = None, ""
            self._discard_map()
            return {"source_input": None}

    def source_input_info(self):
        if self.source_frame is None:
            return None
        return {"filename": self.source_filename, "rows": len(self.source_frame),
                "columns": self.source_frame.columns.tolist()}

    def _discard_map(self):
        self.result, self.curation, self.files = None, None, {}
        self.generic_records = None
        self.activity_confidence = None
        self.possible_sources = None
        self.group_cache, self.group_detail_cache = {}, {}
        self.job = {"state": "ready" if self.frame is not None else "empty",
                    "message": "Input updated. Review the options, then build the map."}

    def query_input_info(self):
        if self.query_frame is None:
            return None
        columns = self.query_frame.columns.tolist()
        available = {str(column).strip().lower(): column for column in columns}
        smiles = next((available[name] for name in (
            "smiles", "canonical_smiles", "canonical smiles", "representation_smiles"
        ) if name in available), None)
        return {"filename": self.query_filename, "rows": len(self.query_frame), "columns": columns,
                "detected": {"smiles": smiles}, "sha256": hashlib.sha256(self.query_raw).hexdigest()}

    def input_info(self):
        if self.frame is None:
            return None
        frame = self.frame
        columns = detect_chembl_columns(frame)
        targets = sorted(frame[columns["target_id"]].dropna().astype(str).unique().tolist()) if columns["target_id"] else []
        return {"filename": self.filename, "rows": len(frame), "columns": frame.columns.tolist(),
                "preview": json_rows(frame.head(5)), "detected": columns, "targets": targets,
                "sha256": hashlib.sha256(self.raw).hexdigest()}

    def start(self, options, background=True):
        with self.lock:
            if self.frame is None:
                raise ValueError("Load a CSV first.")
            if self.job["state"] == "running":
                raise ValueError("A map is already running.")
            options = validate_options(options)
            for key in ("query_smiles_col", "query_id_col"):
                column = options.get(key)
                if column and (self.query_frame is None or column not in self.query_frame):
                    raise ValueError(f"{key} requires a column in the unmeasured-molecule input.")
            if options.get("unit_column") and options["unit_column"] not in self.frame:
                raise ValueError("The selected unit column was not found.")
            if options.get("reference_column") and options["reference_column"] not in self.frame:
                raise ValueError("The selected reference column was not found.")
            self.result, self.curation, self.files = None, None, {}
            self.generic_records = None
            self.activity_confidence = None
            self.possible_sources = None
            self.group_cache = {}
            self.group_detail_cache = {}
            self.job = {"state": "running", "message": "Preparing input records...", "started": time.time()}
        if background:
            threading.Thread(target=self.run, args=(options,), daemon=True).start()
        else:
            self.run(options)

    def run(self, options):
        directory = self.output_root / (time.strftime("%Y%m%d-%H%M%S") + "-" + secrets.token_hex(3))
        try:
            directory.mkdir(parents=True)
            (directory / "input.csv").write_bytes(self.raw)
            (directory / "app_options.json").write_text(
                json.dumps(options, indent=2), encoding="utf-8", newline="\n")
            source = {"input_file": self.filename, "input_file_sha256": hashlib.sha256(self.raw).hexdigest(),
                      "interface": "local-browser", "chemgridmap_version": __version__}
            if self.query_frame is not None:
                (directory / "query_input.csv").write_bytes(self.query_raw)
                source.update(query_file=self.query_filename,
                              query_file_sha256=hashlib.sha256(self.query_raw).hexdigest())
            source["reliability_parameters"] = {
                "activity_confidence": options["activity_confidence"],
                "evidence_model": options["evidence_model"],
                "merge_identical_documents": options["merge_identical_documents"],
                "confidence_threshold": options["support_level"],
                "offset_display": options["offset_display"],
                "offset_min_shared": options["offset_min_shared"],
                "offset_threshold": options["offset_min_difference"],
                "offset_sign_agreement": options["offset_sign_agreement"],
            }
            data = self.frame.copy()
            curation, curation_files = None, {}
            generic_report = None
            generic_records = None
            assumed_sigma = options.get("assumed_sigma")
            value_unit, unit_column = options.get("value_unit", "p"), options.get("unit_column")
            query_structure_policy = "canonical"
            if options["input_format"] == "chembl":
                self.update(message="Validating structures and auditing repeated measurements...")
                curation = curate_chembl_activity_data(
                    data, target_id=options["target_id"], activity_type=options["activity_type"],
                    lower_threshold=options["lower_threshold"], upper_threshold=options["upper_threshold"],
                    conflict_range_threshold=options["conflict_range_threshold"],
                    structure_policy=options["structure_policy"], molecule_identity=options["molecule_identity"],
                    identity_conflicts=options["identity_conflicts"],
                )
                curation.report.update(source)
                source["curation_parameters"] = curation.report["parameters"]
                query_structure_policy = curation.report["parameters"]["structure_policy"]
                document_column = options["document_column"] or curation.report["detected_columns"].get("document_id")
                if options["activity_confidence"] and (document_column or assumed_sigma is not None):
                    self.update(message="Estimating how well source documents support each activity class...")
                    try:
                        curation.molecules, fit = annotate_activity_confidence(
                            curation.molecules, curation.retained_records, lower=options["lower_threshold"],
                            upper=options["upper_threshold"], value="activity_pchembl_raw", document=document_column,
                            merge_identical=options["merge_identical_documents"],
                            model=options["evidence_model"],
                            confidence_threshold=options["support_level"],
                            offset_display=options["offset_display"],
                            offset_min_shared=options["offset_min_shared"],
                            offset_threshold=options["offset_min_difference"],
                            offset_sign_agreement=options["offset_sign_agreement"],
                            assumed_sigma=assumed_sigma)
                        curation.report["activity_confidence"] = fit
                        source["activity_confidence"] = fit
                        self.activity_confidence = fit
                        if assumed_sigma is not None:
                            curation.report["warnings"].extend(fit.get("warnings", []))
                    except ValueError as exc:
                        if assumed_sigma is not None or options["document_column"]:
                            raise
                        curation.report["warnings"].append(f"Activity class support was not estimated: {exc}")
                elif options["activity_confidence"]:
                    curation.report["warnings"].append(
                        "Activity class support was not estimated: the input has no document identifiers. "
                        "Use an assumed σ to calculate support with an explicitly assumed log-unit SD.")
                curation_files = save_chembl_curation(curation, directory, "map")
                source["curation_outputs"] = file_digests(curation_files)
                data = curation.molecules
                smiles_col, value_col, label_col = "canonical_smiles", "activity_pchembl", "activity_class"
            else:
                smiles_col, value_col, label_col = options["smiles_col"], options["value_col"], options["label_col"]
                if options["compound_column"]:
                    if not value_col:
                        raise ValueError("A pActivity column is required for generic activity records.")
                    confidence = bool((options["document_column"] or assumed_sigma is not None)
                                      and options["activity_confidence"])
                    data, generic_records, fit, report = prepare_generic_activity_records(
                        data,
                        compound_column=options["compound_column"],
                        value_column=value_col,
                        document_column=options["document_column"],
                        smiles_column=smiles_col,
                        relation_column=options["relation_column"],
                        lower=options["lower_threshold"],
                        upper=options["upper_threshold"],
                        conflict_range_threshold=options["conflict_range_threshold"],
                        confidence=confidence,
                        confidence_threshold=options["support_level"],
                        merge_identical=options["merge_identical_documents"],
                        evidence_model=options["evidence_model"],
                        offset_display=options["offset_display"],
                        offset_min_shared=options["offset_min_shared"],
                        offset_threshold=options["offset_min_difference"],
                        offset_sign_agreement=options["offset_sign_agreement"],
                        value_unit=value_unit, unit_column=unit_column, assumed_sigma=assumed_sigma,
                    )
                    source["generic_activity"] = report
                    generic_report = report
                    if fit is not None:
                        source["activity_confidence"] = fit
                        self.activity_confidence = fit
                    smiles_col, value_col, label_col = (
                        "canonical_smiles", "activity_pchembl", "activity_class"
                    )
            self.update(message="Calculating the projection, assigning grid cells and rendering. UMAP's first run may take longer.")
            color_modes = ["neighbor_recall"]
            if "class_status" in data:
                color_modes.append("class_support")
            if options["grouping"] != "none":
                color_modes.append("structure_group")
            result = build_grid_map(
                data, output_dir=directory, name="map", smiles_col=smiles_col,
                value_col=value_col, label_col=label_col, representation=options["representation"],
                projection=options["projection"], lower_threshold=options["lower_threshold"],
                upper_threshold=options["upper_threshold"], random_state=options["random_state"],
                grid_occupancy=options["grid_occupancy"], render_detail=options["render_detail"],
                output_formats=options["output_formats"], provenance=source,
                grouping=options["grouping"], group_similarity=options["group_similarity"],
                cluster_col=options["cluster_col"], cluster_source=options["cluster_source"],
                grid_refinement=options["grid_refinement"], core_rule=options["core_rule"],
                core_min_coverage=options["core_min_coverage"], fmcs_timeout=options["fmcs_timeout"],
                confidence_threshold=options["support_level"],
                offset_display=options["offset_display"],
                offset_min_shared=options["offset_min_shared"],
                offset_threshold=options["offset_min_difference"],
                offset_sign_agreement=options["offset_sign_agreement"],
                map_color_modes=color_modes, map_color_formats=["svg"],
                value_unit=value_unit if options["input_format"] == "molecule" and not options["compound_column"] else "p",
                unit_column=unit_column if options["input_format"] == "molecule" and not options["compound_column"] else None,
                query_data=self.query_frame, query_smiles_col=options.get("query_smiles_col"),
                query_id_col=options.get("query_id_col"), query_structure_policy=query_structure_policy,
            )
            result.output_files.update(curation_files)
            files = {Path(path).name: Path(path) for path in result.output_files.values()}
            files.update({name: directory / name for name in ("input.csv", "app_options.json")})
            if self.query_frame is not None:
                files["query_input.csv"] = directory / "query_input.csv"
            archive = directory / "chemgridmap_results.zip"
            with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as bundle:
                for name, path in files.items():
                    bundle.write(path, name)
            files[archive.name] = archive
            report = curation.report if curation else result.preprocessing
            if generic_report is not None:
                # The grid preparation counts molecules; the browser summary counts source records.
                report = {**report, "generic_activity": generic_report,
                          "activity_confidence": self.activity_confidence,
                          "warnings": generic_report["warnings"],
                          "input_rows": generic_report["input_rows"],
                          "retained_activity_records": generic_report["retained_activity_records"],
                          "repeated_conflicted_molecules": int(data["wide_record_range"].eq(1).sum())}
            if result.query_report:
                report = {**report, "query_molecules": result.query_report}
            confidence_report = report.get("activity_confidence")
            if confidence_report is not None:
                # Include unflagged evaluated pairs in the browser audit only.
                # The fitted model, saved reports and manifest remain untouched.
                pairs = confidence_report.get("document_pair_offsets", [])
                top_pairs = sorted(pairs, key=lambda pair: (
                    not pair["systematic_offset"], -abs(pair["median_offset"]),
                    str(pair["document_a"]), str(pair["document_b"])))[:10]
                report = {**report, "activity_confidence": {
                    **confidence_report, "document_offsets_top10": top_pairs}}
            # Literature annotations belong to the browser audit, never the map outputs or fit.
            if curation is not None:
                document_column = curation.report["detected_columns"].get("document_id")
                sources = build_sources(curation.retained_records, document_column)
                report = {**report, "document_coverage": {
                    "kind": "chembl", **document_coverage(sources)}}
            elif generic_records is not None and options.get("reference_column"):
                reference_records = self._reference_records(generic_records, options)
                sources = build_sources(reference_records, options["document_column"],
                                        options["reference_column"])
                report = {**report, "document_coverage": {
                    "kind": "reference", **document_coverage(sources)}}
            if result.query_report:
                report = {**report, "query_reference_matches": self._query_reference_matches(
                    result, curation, generic_records, options)}
            possible_sources = None
            if generic_records is not None and self.source_frame is not None:
                self.update(message="Matching structures to the source lookup file...")
                possible_sources, lookup_summary = match_source_lookup(
                    data, generic_records, self.source_frame)
                report = {**report, "source_lookup": {
                    "file": self.source_filename, **lookup_summary}}
            with self.lock:
                self.result, self.curation, self.files = result, curation, files
                self.generic_records = generic_records
                self.possible_sources = possible_sources
                self.job.update(state="done", message="Map ready. Select a cell to inspect its evidence.",
                                output_dir=str(directory), molecules=len(result.data),
                                report=report, metrics=json_rows(result.metrics)[0],
                                options=options, files=list(files), elapsed=time.time() - self.job["started"])
        except Exception as exc:
            self.update(state="error", message=str(exc), output_dir=str(directory))

    def update(self, **values):
        with self.lock:
            self.job.update(values)

    def group_detail(self, group_id):
        with self.lock:
            if self.result is None:
                raise ValueError("Build a map first.")
            if group_id in self.group_detail_cache:
                return self.group_detail_cache[group_id]
            summary = next((g for g in self.result.structure_groups if g['structure_group'] == group_id), None)
            if summary is None:
                raise ValueError("Unknown structure group.")
            payload = {'summary': summary, 'core_svg': _group_drawing(summary, 'core_smiles')}
            self.group_detail_cache[group_id] = payload
            return payload

    def compare_group(self, group_id):
        with self.lock:
            if self.result is None:
                raise ValueError("Build a map first.")
            if group_id in self.group_cache:
                return self.group_cache[group_id]
            summary = next((g for g in self.result.structure_groups if g['structure_group'] == group_id), None)
            if summary is None:
                raise ValueError("Unknown structure group.")
            selected, svg, atom_audit = group_comparison(self.result.data, summary)
            members = self.result.data[self.result.data.structure_group.eq(group_id)]
            drawings = {field: _group_drawing(summary, field) for field in ('core_smiles', 'representative_smiles')}
            payload = {'summary': summary, 'svg': svg, 'drawings': drawings,
                       'members': json_rows(members), 'selected_count': len(selected),
                       'selection_rule': 'Representative first; then scaffold, canonical SMILES and identifier. At most 12 cards; all members are exported.',
                       'atom_matches': atom_audit}
            files = {'comparison.svg': svg.encode(), 'members.csv': members.to_csv(
                         index=False, lineterminator="\n").encode(),
                     'selected_members.csv': selected.to_csv(
                         index=False, lineterminator="\n").encode(),
                     'atom_matches.json': json.dumps(atom_audit, indent=2).encode(),
                     'group_summary.json': json.dumps({'summary': summary, 'protocol': self.result.grouping_protocol}, indent=2).encode()}
            digests = {name: hashlib.sha256(content).hexdigest() for name, content in files.items()}
            files['file_digests.json'] = json.dumps(digests, indent=2).encode()
            name = 'comparison_' + hashlib.sha256(group_id.encode()).hexdigest()[:12] + '.zip'
            path = Path(self.job['output_dir']) / name
            with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as archive:
                for filename, content in files.items():
                    archive.writestr(filename, content)
            self.files[name] = path
            payload['archive'] = name
            self.group_cache[group_id] = payload
            return payload

    def snapshot(self):
        with self.lock:
            return dict(self.job)

    def _reference_records(self, records, options):
        """Recover selected reference strings without altering calculation records."""
        column = options.get("reference_column")
        if not column:
            return records
        references = _read_input(self.raw, dtype=str, keep_default_na=False)[column]
        records = records.copy()
        records[column] = records.source_row.map(references)
        return records

    def _query_reference_matches(self, result, curation, generic_records, options):
        """Count retained documents for the browser only; never create audit exports."""
        checks = result.query_check[result.query_check.status.eq("in_reference")]
        if checks.empty:
            return []
        records = curation.retained_records if curation is not None else generic_records
        document_column = (curation.report["detected_columns"].get("document_id")
                           if curation is not None else options["document_column"])
        if records is not None:
            records = self._reference_records(records, options)
            if document_column:
                columns = list(dict.fromkeys(["molecule_identity_key", document_column,
                                              *([options["reference_column"]]
                                                if options.get("reference_column") else [])]))
                records = records[columns]
        # Query references already contain exact keys from the prepared map.
        molecules = dict(zip(molecule_keys(result.data), json_rows(result.data)))
        retained_groups = (dict(tuple(records.groupby("molecule_identity_key", sort=False)))
                           if records is not None and document_column else {})
        source_cache = {}
        matches = []
        for check in checks.itertuples(index=False):
            molecule = molecules[check.reference_id]
            sources = []
            # Reference-only text must not invent document identifiers or counts.
            if records is not None and document_column:
                retained = retained_groups[check.reference_id]
                # Counts depend on distinct document/reference pairs, not repeated rows.
                # Many compounds share a paper; resolve those same links only once.
                reference_column = options.get("reference_column")
                signature = (frozenset(retained[[document_column, reference_column]].itertuples(
                    index=False, name=None)) if reference_column else frozenset(retained[document_column]))
                if signature not in source_cache:
                    source_cache[signature] = build_sources(
                        retained, document_column, options.get("reference_column"))
                sources = source_cache[signature]
            matches.append({"query_id": check.query_id, "reference_id": check.reference_id,
                            "molecule_id": molecule.get("molecule_id") or check.reference_id,
                            "parent_molecule_id": molecule.get("parent_molecule_id") or "",
                            "source_count": len(sources),
                            "linked_source_count": sum(bool(source["links"]) for source in sources)})
        return matches

    def _generic_record_audit(self, chosen, summary):
        row = chosen.iloc[0]
        records = self.generic_records
        source_rows = [int(value) for value in str(row.source_rows).split(";")]
        evidence = records[records.molecule_identity_key.eq(row.molecule_identity_key)].copy()
        evidence = evidence.sort_values("source_row")
        if (len(evidence) != int(row.n_records) or evidence.source_row.duplicated().any()
                or set(evidence.source_row) != set(source_rows)):
            raise ValueError("Record counts do not match this cell; check the retained-record file.")
        # Legacy pActivity-only input lacks this display alias in retained records.
        evidence["activity_pchembl_raw"] = pd.to_numeric(evidence.pchembl_value, errors="raise")
        values = evidence.activity_pchembl_raw
        median = float(values.median())
        if not values.map(math.isfinite).all() or not math.isclose(
                median, float(row.activity_pchembl), rel_tol=0, abs_tol=1e-9):
            raise ValueError("Recalculated activity does not match the map annotation.")
        options = self.job["options"]
        touched = bool(row.get("offset_pair_touched", False))
        document_column = ("document_chembl_id" if "document_chembl_id" in evidence
                           else options["document_column"])
        pairs, warnings = offset_pair_details(
            evidence, document_column, self.activity_confidence or {}, touched)
        summary.update(
            input_format="generic_activity", record_count_verified=True, median_verified=True,
            n_records=len(evidence), median_pchembl=median, document_offsets=pairs, warnings=warnings,
            interpretation=("Generic activity records: the listed rows are the retained records "
                            "for this compound; values were converted to pActivity."))
        priority = ["source_row", options["compound_column"], options["document_column"],
                    options["value_col"], options.get("unit_column"), options["relation_column"],
                    "activity_pchembl_raw"]
        hidden = {"molecule_identity_key", "pchembl_value", "canonical_smiles", "document_chembl_id"}
        if options.get("unit_column"):
            hidden.update(("raw_activity_value", "raw_activity_unit"))
        # Explicit source-column selections are not internal duplicate columns.
        columns = list(dict.fromkeys(name for name in priority if name in evidence))
        columns.extend(name for name in evidence if name not in hidden and name not in columns)
        return evidence[columns]

    def inspect(self, row=None, col=None, molecule_id=None):
        with self.lock:
            if self.result is None:
                raise ValueError("Build a map first.")
            data, curation = self.result.data, self.curation
            molecule_id, query_match = _query_selection(self.result.query_check, molecule_id)
            if "measurement_status" in data:
                selected = _select_molecule(data, row, col, molecule_id)
                if selected.iloc[0].measurement_status == "not_measured":
                    query = selected.iloc[0]
                    neighbors = self.result.query_neighbors[
                        self.result.query_neighbors.query_id.eq(query.query_id)]
                    svg = ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 320 320">'
                           + (_svg_molecule_content(query.canonical_smiles, 320) or '') + '</svg>')
                    detail = {"molecule": json_rows(selected)[0], "query_neighbors": json_rows(neighbors),
                              "group_similarity": self.job["options"]["group_similarity"], "svg": svg}
                    if query_match is not None:
                        detail["query_match"] = query_match
                    return detail, pd.DataFrame()
            if curation is not None:
                _, run_bound, binding_warnings, manifest = load_audit_context(
                    Path(self.result.output_files['coordinates']),
                    Path(self.result.output_files['retained_records']),
                    return_manifest=True)
                chosen, records, summary = inspect_entry(
                    data, curation.retained_records, molecule_id=molecule_id,
                    cell=(row, col) if molecule_id is None else None,
                    parameters=curation.report["parameters"],
                    activity_confidence=manifest.get("activity_confidence"),
                    k=int(self.result.metrics.iloc[0].get("k", 10)))
                summary['run_binding_verified'] = run_bound
                summary['warnings'].extend(binding_warnings)
            else:
                chosen = _select_molecule(data, row, col, molecule_id)
                records = pd.DataFrame()
                summary = {"interpretation": "Molecule-table input: no assay-record audit is available.",
                           "record_count_verified": False, "median_verified": False,
                           "quality_annotations_verified": False, "run_binding_verified": False,
                           "sigma_doc": (self.activity_confidence or {}).get("sigma_doc"),
                           "sigma_doc_degrees_of_freedom": (
                               self.activity_confidence or {}
                           ).get("sigma_doc_degrees_of_freedom"),
                           "multi_unit_entry_count": (
                               self.activity_confidence or {}
                           ).get("multi_unit_entry_count")}
                if self.generic_records is not None:
                    records = self._generic_record_audit(chosen, summary)
            options = self.job["options"]
            document_column = (curation.report["detected_columns"].get("document_id")
                               if curation is not None else options["document_column"])
            source_records = self._reference_records(records, options)
            sources = build_sources(source_records, document_column, options.get("reference_column"))
            # Keep the on-screen record schema intact; the four appended fields are downloads.
            display_records = records
            if not records.empty:
                records = enrich_source_records(records, document_column)
            fit = self.activity_confidence or {}
            if fit.get("method") == "assumed_sigma":
                summary.update(method="assumed_sigma", assumed_sigma=fit["assumed_sigma"],
                               evidence_unit_source=fit["evidence_unit_source"])
                summary.pop("sigma_doc_degrees_of_freedom", None)
            selected = chosen.iloc[0]
            row, col = int(selected.grid_row), int(selected.grid_col)
            neighbors = data[data.grid_row.sub(row).abs().le(1) & data.grid_col.sub(col).abs().le(1)]
            directory = Path(self.job["output_dir"])
            # Fixed per-cell files make reinspection reproducible and avoid re-running the projection.
            image = render_svg(chosen, directory / f"cell_{row}_{col}.svg", tile_size=320, molecule_margin=12)
            neighborhood = render_svg(neighbors, directory / f"neighborhood_{row}_{col}.svg", tile_size=200)
            position = data.index.get_loc(chosen.index[0])
            neighbor_table, neighbor_summary = neighborhood_evidence(
                data, self.result.representation, position,
                k=int(self.result.metrics.iloc[0].get('k', 10)),
                metric=str(self.result.metrics.iloc[0].representation_distance_metric))
            neighbor_summary['representation'] = str(self.result.metrics.iloc[0].representation_type)
            summary['neighborhood'] = neighbor_summary
            files = {
                'selected_molecule.svg': image.read_bytes(),
                'source_records.csv': records.to_csv(
                    index=False, lineterminator="\n").encode('utf-8'),
                'selected_molecule.csv': chosen.to_csv(
                    index=False, lineterminator="\n").encode('utf-8'),
                'neighbors.csv': neighbor_table.to_csv(
                    index=False, lineterminator="\n").encode('utf-8'),
            }
            summary['evidence_digests'] = {name: hashlib.sha256(content).hexdigest() for name, content in files.items()}
            files['audit_summary.json'] = json.dumps(summary, indent=2, allow_nan=False).encode('utf-8')
            archive = directory / f'audit_{row}_{col}.zip'
            with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED) as bundle:
                for name, content in files.items():
                    bundle.writestr(name, content)
            self.files[archive.name] = archive
            detail = {"molecule": json_rows(chosen)[0], "summary": summary,
                    "sources": sources,
                    "records": json_rows(display_records if self.generic_records is not None else display_records.head(100)),
                    "total_records": len(records),
                    "neighbors": json_rows(neighbor_table), "audit_file": archive.name,
                    "svg": image.read_text(encoding="utf-8"),
                    "neighborhood_svg": neighborhood.read_text(encoding="utf-8")}
            if self.generic_records is not None:
                detail["record_columns"] = display_records.columns.tolist()
            if self.possible_sources is not None and not any(source["links"] for source in sources):
                detail["possible_sources"] = self.possible_sources.get(
                    str(selected.molecule_identity_key), [])
                detail["source_lookup_file"] = self.source_filename
            if query_match is not None:
                detail["query_match"] = query_match
            return detail, records


def validate_options(options):
    defaults = {"input_format": "chembl", "target_id": None, "activity_type": "IC50",
                "smiles_col": "smiles", "value_col": None, "label_col": None,
                "representation": "morgan", "projection": "umap", "lower_threshold": 6.0,
                "upper_threshold": 7.0, "conflict_range_threshold": 1.0, "random_state": 42,
                "grid_occupancy": .4, "render_detail": "auto", "structure_policy": "auto",
                "grouping": "butina", "group_similarity": .6, "cluster_col": None, "cluster_source": None,
                "molecule_identity": "auto", "output_formats": ["svg", "png", "pdf"],
                "grid_refinement": "neighbors", "core_rule": "coverage", "core_min_coverage": .8,
                "fmcs_timeout": 1, "compound_column": None, "document_column": None,
                "relation_column": None, "activity_confidence": True,
                "evidence_model": "gaussian", "merge_identical_documents": True,
                "support_level": .8, "offset_display": "none",
                "offset_min_shared": 3, "offset_min_difference": 1.0,
                "offset_sign_agreement": .8, "identity_conflicts": "error"}
    optional = {"value_unit", "unit_column", "assumed_sigma", "query_smiles_col", "query_id_col",
                "reference_column"}
    if not isinstance(options, dict) or set(options) - set(defaults) - optional:
        raise ValueError("Unknown map option.")
    values = {**defaults, **options}
    # An unused reference selector must preserve the BASE app_options.json byte for byte.
    if values.get("reference_column") is None or values.get("reference_column") == "":
        values.pop("reference_column", None)
    if values.get("value_unit", "p") not in {"p", "M", "mM", "uM", "nM", "pM"}:
        raise ValueError("Invalid value_unit.")
    for key in ("unit_column", "query_smiles_col", "query_id_col", "reference_column"):
        if key in values and values[key] is not None:
            if not isinstance(values[key], str) or not values[key].strip():
                raise ValueError(f"{key} must be a non-empty column name.")
    if values.get("reference_column") and (values["input_format"] != "molecule" or not values["compound_column"]):
        raise ValueError("Reference columns require generic activity records.")
    if values.get("unit_column") and "value_unit" in options:
        raise ValueError("Unit and unit column cannot be used together.")
    if values["input_format"] == "chembl" and (values.get("value_unit", "p") != "p" or values.get("unit_column")):
        raise ValueError("ChEMBL input already contains pChEMBL values; units apply to generic input.")
    if values.get("assumed_sigma") is not None:
        values["assumed_sigma"] = float(values["assumed_sigma"])
        validate_assumed_sigma(values["assumed_sigma"])
        if values["input_format"] == "molecule" and not values["compound_column"]:
            raise ValueError("Assumed σ requires activity records; select Generic activity records and a compound column.")
        if not values["activity_confidence"]:
            raise ValueError("Assumed σ requires activity class support.")
    choices = {"input_format": {"chembl", "molecule"}, "representation": {"morgan", "descriptors"},
               "projection": {"pca", "tsne", "umap"}, "render_detail": {"auto", "full", "overview"},
               "structure_policy": {"auto", "parent", "as-recorded"},
               "molecule_identity": {"auto", "parent-id", "canonical-smiles"},
               "grid_refinement": {"none", "neighbors"},
               "core_rule": {"most_frequent", "coverage", "fmcs"},
               "evidence_model": {"gaussian", "exclude-flagged"},
               "offset_display": {"uncertain", "none"},
               "identity_conflicts": {"error", "exclude"}}
    choices['grouping'] = {'butina', 'external', 'none'}
    for key, allowed in choices.items():
        if values[key] not in allowed:
            raise ValueError(f"Invalid {key}.")
    for key in ("lower_threshold", "upper_threshold", "conflict_range_threshold", "grid_occupancy"):
        values[key] = float(values[key])
    if not 0.1 <= values["grid_occupancy"] <= 1:
        raise ValueError("Grid occupancy must be between 10% and 100%.")
    if not (-100 < values["lower_threshold"] < values["upper_threshold"] < 100):
        raise ValueError("Activity thresholds must be finite, with the lower threshold below the upper.")
    if not 0 <= values["conflict_range_threshold"] < 100:
        raise ValueError("Wide record range threshold must be finite and non-negative.")
    values["support_level"] = float(values["support_level"])
    if not 0.5 < values["support_level"] < 1:
        raise ValueError("Support level must be between 0.5 and 1 (exclusive).")
    min_shared = values["offset_min_shared"]
    if isinstance(min_shared, bool) or int(min_shared) != min_shared or int(min_shared) < 2:
        raise ValueError("Offset minimum shared entries must be an integer of at least 2.")
    values["offset_min_shared"] = int(min_shared)
    values["offset_min_difference"] = float(values["offset_min_difference"])
    if not values["offset_min_difference"] > 0:
        raise ValueError("Offset minimum difference must be greater than 0.")
    values["offset_sign_agreement"] = float(values["offset_sign_agreement"])
    if not 0.5 < values["offset_sign_agreement"] <= 1:
        raise ValueError("Offset sign agreement must lie in (0.5, 1].")
    values["random_state"] = int(values["random_state"])
    values['group_similarity'] = float(values['group_similarity'])
    if not 0 < values['group_similarity'] <= 1:
        raise ValueError("Group similarity must be in (0, 1].")
    values["core_min_coverage"] = float(values["core_min_coverage"])
    if not 0.5 <= values["core_min_coverage"] <= 1:
        raise ValueError("Core coverage target must be between 50% and 100%.")
    values["activity_confidence"] = bool(values["activity_confidence"])
    values["merge_identical_documents"] = bool(values["merge_identical_documents"])
    values["fmcs_timeout"] = int(values["fmcs_timeout"])
    if values["fmcs_timeout"] < 1:
        raise ValueError("FMCS timeout must be a positive integer.")
    if values["compound_column"] and not values["value_col"]:
        raise ValueError("Generic activity records require a pActivity value column.")
    if values['grouping'] == 'external' and (values['input_format'] != 'molecule' or not values['cluster_col'] or not values['cluster_source']):
        raise ValueError("External groups require molecule-table input, a group column, and its source.")
    if not 0 <= values["random_state"] < 2**32:
        raise ValueError("The seed must be between 0 and 4294967295.")
    if not isinstance(values["output_formats"], list) or not set(values["output_formats"]).issubset({"svg", "png", "pdf"}):
        raise ValueError("Choose SVG, PNG and/or PDF output.")
    values["output_formats"] = list(dict.fromkeys(["svg"] + values["output_formats"]))
    return values


class LocalServer(ThreadingHTTPServer):
    daemon_threads = True

    def __init__(self, session, port=0):
        self.session = session
        super().__init__(("127.0.0.1", port), LocalHandler)
        self.origin = f"http://127.0.0.1:{self.server_port}"


class LocalHandler(BaseHTTPRequestHandler):
    def log_message(self, *args):
        # Do not log CSV contents, record lookups, or the session credential.
        pass

    def send(self, data, content_type="application/json", status=200, filename=None):
        if not isinstance(data, bytes):
            data = json.dumps(data, allow_nan=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Referrer-Policy", "no-referrer")
        self.send_header("Content-Security-Policy", "default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; img-src 'self' data: blob:; object-src 'none'; frame-ancestors 'none'")
        if filename:
            self.send_header("Content-Disposition", f'attachment; filename="{filename}"')
        self.end_headers()
        self.wfile.write(data)

    def authorized(self, api=True):
        expected_host = urlsplit(self.server.origin).netloc
        if self.headers.get("Host") != expected_host:
            self.send({"error": "Unexpected host."}, status=403)
            return False
        origin = self.headers.get("Origin")
        if origin and origin != self.server.origin:
            self.send({"error": "Only this local interface can access the workspace."}, status=403)
            return False
        if api and not secrets.compare_digest(self.headers.get("X-ChemGridMap-Token", ""), self.server.session.token):
            self.send({"error": "Session expired. Reopen the app from its launcher."}, status=403)
            return False
        return True

    def do_GET(self):
        parsed = urlsplit(self.path)
        path = parsed.path
        if not self.authorized(api=path.startswith("/api/")):
            return
        session = self.server.session
        try:
            if path in ("/", "/app.js", "/i18n.json", "/app.css", "/example-map.svg", "/example-molecule.svg"):
                asset = ASSETS / ("index.html" if path == "/" else path[1:])
                self.send(asset.read_bytes(), mimetypes.guess_type(asset.name)[0] or "text/plain")
            elif path == "/api/session":
                payload = {"version": __version__, "umap": importlib.util.find_spec("umap") is not None,
                           "output_root": str(session.output_root), "demo_available": DEMO.exists(),
                           "input": session.input_info(), "job": session.snapshot()}
                if session.query_frame is not None:
                    payload["query_input"] = session.query_input_info()
                if session.source_frame is not None:
                    payload["source_input"] = session.source_input_info()
                self.send(payload)
            elif path == "/api/job":
                self.send(session.snapshot())
            elif path == "/api/points":
                if session.result is None:
                    raise ValueError("Build a map first.")
                self.send(point_rows(session.result.data))
            elif path == "/api/groups":
                if session.result is None:
                    raise ValueError("Build a map first.")
                self.send({'groups': session.result.structure_groups, 'protocol': session.result.grouping_protocol})
            elif path == "/api/group":
                self.send(session.group_detail(parse_qs(parsed.query).get('group', [''])[0]))
            elif path == "/api/compare":
                self.send(session.compare_group(parse_qs(parsed.query).get('group', [''])[0]))
            elif path == "/api/file":
                name = parse_qs(parsed.query).get("name", [""])[0]
                file = session.files.get(name)
                if file is None:
                    raise ValueError("This output file is not available.")
                self.send(file.read_bytes(), mimetypes.guess_type(name)[0] or "application/octet-stream", filename=name)
            elif path in ("/api/inspect", "/api/evidence.csv"):
                query = parse_qs(parsed.query)
                selection = {"molecule_id": query["id"][0]} if "id" in query else {
                    "row": int(query["row"][0]), "col": int(query["col"][0])}
                detail, records = session.inspect(**selection)
                if path == "/api/evidence.csv":
                    self.send(records.to_csv(index=False, lineterminator="\n").encode("utf-8-sig"),
                              "text/csv", filename="source_records.csv")
                else:
                    self.send(detail)
            else:
                self.send({"error": "Not found."}, status=404)
        except (ValueError, KeyError, OSError) as exc:
            self.send({"error": str(exc)}, status=400)

    def do_POST(self):
        if not self.authorized():
            return
        session = self.server.session
        try:
            length = int(self.headers.get("Content-Length", "0"))
            if not 0 <= length <= MAX_UPLOAD:
                raise ValueError("Maximum upload size is 50 MB.")
            content = self.rfile.read(length)
            path = urlsplit(self.path).path
            if path == "/api/input":
                self.send(session.load(content, unquote(self.headers.get("X-Filename", "input.csv"))))
            elif path == "/api/query-input":
                self.send(session.load_query(content, unquote(self.headers.get("X-Filename", "queries.csv"))))
            elif path == "/api/query-input/clear":
                self.send(session.clear_query())
            elif path == "/api/source-input":
                self.send(session.load_source(content, unquote(self.headers.get("X-Filename", "sources.csv"))))
            elif path == "/api/source-input/clear":
                self.send(session.clear_source())
            elif path == "/api/demo":
                self.send(session.load(DEMO.read_bytes(), DEMO.name))
            elif path == "/api/build":
                session.start(json.loads(content))
                self.send(session.snapshot(), status=202)
            elif path == "/api/shutdown":
                if session.snapshot()["state"] == "running":
                    raise ValueError("Wait for the current map to finish before closing the app.")
                self.send({"message": "Local app stopped. Your saved results remain on disk."})
                threading.Thread(target=self.server.shutdown, daemon=True).start()
            else:
                self.send({"error": "Not found."}, status=404)
        except (ValueError, OSError) as exc:
            self.send({"error": str(exc)}, status=400)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Open ChemGridMap's local browser interface.")
    parser.add_argument("--output-dir", type=Path, default=Path.home() / "ChemGridMap" / "runs")
    parser.add_argument("--port", type=int, default=0, help="Local port; default selects an available port.")
    parser.add_argument("--no-browser", action="store_true")
    parser.add_argument("--smoke-test", action="store_true", help="Run the bundled ChEMBL example without a browser.")
    parser.add_argument("--smoke-projection", choices=["pca", "umap"], default="umap")
    args = parser.parse_args(argv)
    session = AppSession(args.output_dir)
    if args.smoke_test:
        session.load(DEMO.read_bytes(), DEMO.name)
        formats = ["svg", "png", "pdf"] if args.smoke_projection == "pca" else ["svg"]
        session.start({"projection": args.smoke_projection, "target_id": "CHEMBL205", "output_formats": formats}, background=False)
        if session.job["state"] != "done":
            raise RuntimeError(session.job["message"])
        detail, _ = session.inspect(molecule_id="CHEMBL20")
        assert session.job["molecules"] == 889 and detail["summary"]["n_records"] == 54
        assert detail["summary"]["median_verified"] and abs(detail["summary"]["median_pchembl"] - 7.185) < 1e-8
        assert all(session.files[f"map.{ext}"].stat().st_size > 0 for ext in formats)
        print(json.dumps({"molecules": 889, "acetazolamide_records": 54,
                          "median_verified": True, "verified_formats": formats}))
        return 0
    server = LocalServer(session, args.port)
    # The credential stays in the URL fragment, never in HTTP requests or access logs.
    url = server.origin + "/#" + session.token
    print(f"ChemGridMap {__version__} | Local files stay on this computer.\n{url}\nUse Quit in the interface or Ctrl+C to stop.", flush=True)
    if not args.no_browser:
        threading.Timer(.5, webbrowser.open, args=(url,)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
