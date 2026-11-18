"""Local browser interface. Scientific operations use the same API as the CLI."""
from __future__ import annotations

import argparse
import csv
import copy
from dataclasses import fields, is_dataclass
from datetime import datetime, timezone
import gzip
import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import importlib.util
import io
import json
import math
import mimetypes
import pickle
import platform
from pathlib import Path
import re
import secrets
import signal
import threading
import time
from urllib.parse import parse_qs, unquote, urlsplit
import webbrowser
import zipfile

import numpy as np
import pandas as pd
from rdkit import Chem, rdBase

from . import __version__
from .chembl import curate_chembl_activity_data, detect_chembl_columns, save_chembl_curation
from .evidence import annotate_activity_confidence, prepare_generic_activity_records
from .inspection import (build_sources, build_source_index, document_coverage, enrich_source_records, excluded_bound_records,
                         inspect_entry, load_audit_context, offset_pair_details,
                         match_source_lookup, validate_source_lookup)
from .pipeline import build_grid_map
from .queries import (insert_new_molecules, prepare_query_molecules, place_new_molecules,
                      _new_molecules_svg, NEIGHBOR_EVIDENCE_CATEGORIES)
from .plotting import render_svg
from .provenance import file_digests, source_tree_digest
from .neighborhoods import neighborhood_evidence
from .comparison import (group_comparison, molecule_keys, summarize_group_core, summarize_group_cores,
                         browser_core_jobs, shutdown_core_workers, group_core_members)
from .cores import generic_scaffold_summary
from .plotting import _svg_molecule_content
from .units import validate_assumed_sigma, convert_activity_values

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
    if "is_new_molecule" in data or "is_new" in data:
        columns.extend(name for name in (
            "is_new_molecule", "is_new", "query_id", "label", "neighbor_evidence", "neighbor_evidence_class",
            "nearest_reference_id", "nearest_measured_similarity", "tanimoto", "placement_distance",
            "activity_pchembl", "neighbor_comparison", "similar_reference_count",
            "supported_reference_count", "supported_reference_classes", "reference_support_assessed",
        ) if name in data and name not in columns)
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
        return "query:" + match["query_id"], None
    if match["status"] == "invalid_smiles":
        raise ValueError(f"New molecule {match['query_id']} has invalid SMILES.")
    target = match["reference_id"] or "query:" + (match["first_query_id"] or match["query_id"])
    return target, match


def _group_drawing(summary, field):
    molecule = summary.get("core_smiles_attached") if field == 'core_smiles' else summary[field]
    molecule = molecule or summary[field]
    return ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 180 180">'
            + (_svg_molecule_content(molecule, 180) or '') + '</svg>') if molecule else ''


# Manual overlays share the existing session module so the source-file inventory stays fixed.
MANUAL_COLUMNS = ["query_id", "input_smiles", "canonical_smiles", "value", "unit", "measured_class",
                  "class_source", "batch", "label", "note", "entered_at"]
FEEDBACK_COLUMNS = ["query_id", "canonical_smiles", "label", "result_source", "measured_value",
                    "measured_unit", "activity_pchembl", "measured_class", "class_source",
                    "neighbor_evidence", "neighbor_comparison", "nearest_reference_id",
                    "nearest_reference_similarity", "similar_reference_count", "evidence_reference_classes",
                    "batch", "note", "recorded_at"]
CLASSES = {"active", "medium", "inactive"}


def _rows(frame):
    return json.loads(frame.to_json(orient="records", double_precision=15))


def normalize_measurement(payload, *, lower, upper, optional=False):
    """Use the upload converter, accepting one exact positive value or one class."""
    if not payload:
        if optional:
            return None
        raise ValueError("Enter an exact value and unit, or choose a measured class.")
    mode = payload.get("mode", "value" if payload.get("value") not in (None, "") else "class")
    value, unit, pactivity = None, "", None
    if mode in {"value", "numeric"}:
        value = payload.get("value")
        numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
        if not np.isfinite(numeric) or numeric <= 0:
            raise ValueError("The measured value must be an exact number greater than 0.")
        unit = str(payload.get("unit", "")).strip()
        unit = {"p": "pActivity", "uM": "µM", "μM": "µM"}.get(unit, unit)
        canonical_unit = {"pActivity": "p", "µM": "uM", "μM": "uM"}.get(unit, unit)
        converted, reasons = convert_activity_values(pd.Series([numeric]), value_unit=canonical_unit)
        if reasons.iloc[0]:
            raise ValueError("The measured value or unit is invalid.")
        pactivity = float(converted.iloc[0])
        measured_class = "active" if pactivity >= upper else "inactive" if pactivity <= lower else "medium"
        value = float(numeric)
        class_source = "value"
    elif mode == "class":
        measured_class = payload.get("measured_class", payload.get("class"))
        if measured_class not in CLASSES:
            raise ValueError("Choose active, medium or inactive for the measured class.")
        class_source = "user_class"
    else:
        raise ValueError("Choose a numeric value or a measured class.")
    return {"measured_value": value, "measured_unit": unit, "activity_pchembl": pactivity,
            "measured_class": measured_class, "class_source": class_source,
            "batch": str(payload.get("batch") or ""), "note": str(payload.get("note") or "")}


def neighbor_comparison(measured_class, evidence):
    if evidence not in CLASSES:
        return "no_comparison"
    return "agrees" if measured_class == evidence else "disagrees"


class ManualMoleculesMixin:
    """App-session operations; the session lock covers each complete mutation."""

    def _reset_manual(self):
        self.manual_entries = []
        self.manual_entry_notes = {}
        self.measured_results = {}
        self.measurement_originals = {}

    def _manual_ready(self):
        if self.result is None:
            raise ValueError("Build a reference map first.")
        if self.job["state"] == "running":
            raise ValueError("A map is running. Wait for it to finish.")
        if self.job["options"]["representation"] != "morgan":
            raise ValueError("New molecules require a reference map arranged using Morgan fingerprints.")

    def _measurement(self, payload, *, optional=False):
        options = self.job["options"]
        return normalize_measurement(payload, lower=options["lower_threshold"],
                                     upper=options["upper_threshold"], optional=optional)

    def preview_new(self, payload):
        with self.lock:
            self._manual_ready()
            smiles = str(payload.get("smiles") or "").strip()
            _, check, _ = prepare_query_molecules(pd.DataFrame({"smiles": [smiles]}), self.result.data)
            row = check.iloc[0]
            status = {"mapped": "ok", "invalid_smiles": "invalid"}.get(row.status, row.status)
            canonical = row.canonical_smiles
            result = {"status": status, "canonical_smiles": canonical, "svg": ""}
            if canonical:
                result["svg"] = ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 240 240">'
                                 + (_svg_molecule_content(canonical, 240) or "") + "</svg>")
            if status == "in_reference":
                result["reference_id"] = row.reference_id
                reference_index = molecule_keys(self.result.data).index(row.reference_id)
                reference = self.result.data.iloc[reference_index]
                result.update(grid_row=int(reference.grid_row), grid_col=int(reference.grid_col))
            elif status == "ok" and self.new_result is not None:
                matches = self.new_result.data[self.new_result.data.canonical_smiles.eq(canonical)]
                if not matches.empty:
                    result.update(status="duplicate_new", query_id=str(matches.iloc[0].query_id))
            return result

    def measurement_rows(self):
        """Decorate API/export rows without changing any uploaded molecule fields."""
        if self.new_result is None:
            return pd.DataFrame()
        rows = _rows(self.new_result.all_molecules)
        self._ensure_core_state()
        groups = {group["structure_group"]: group for group in self.groups()}
        manual = {entry["query_id"]: entry for entry in self.manual_entries}
        records = self.new_result.records
        for row in rows:
            identity = row["query_id"]
            if row.get("activity_class") not in CLASSES:
                row["measurement_status"] = "not_measured"
            result = self.measured_results.get(identity)
            row.update(result_source="", measured_value=None, measured_unit="", measured_class="",
                       class_source="", batch="", note=manual.get(identity, {}).get("note", ""), recorded_at="")
            if result:
                row.update(result)
            elif row.get("measurement_status") == "measured":
                retained = (records[records.query_id.astype(str).eq(str(identity))]
                            if not records.empty and "query_id" in records else pd.DataFrame())
                value, unit, batch = row.get("activity_pchembl"), "pActivity", ""
                if len(retained) == 1:
                    value = retained.iloc[0].get("raw_activity_value", value)
                    unit = retained.iloc[0].get("raw_activity_unit", unit)
                if not retained.empty and "evidence_unit" in retained:
                    batch = ";".join(dict.fromkeys(retained.evidence_unit.fillna("").astype(str)))
                row.update(result_source="upload", measured_value=value, measured_unit=unit,
                           measured_class=row.get("activity_class"), class_source="value", batch=batch)
            elif identity in manual:
                row["result_source"] = "manual"
            if row.get("measurement_status") == "measured":
                row["neighbor_comparison"] = neighbor_comparison(row["measured_class"], row.get("neighbor_evidence"))
            row["nearest_reference_similarity"] = row.get("tanimoto")
            group = groups.get(row.get("structure_group"))
            if group and self.current_core_rule != self.job["options"]["core_rule"]:
                pattern = (Chem.MolFromSmarts(group["core_smarts"]) if group.get("core_smarts")
                           else Chem.MolFromSmiles(group["core_smiles"]) if group.get("core_smiles") else None)
                row["contains_group_core"] = bool(pattern is not None and
                    Chem.MolFromSmiles(row["canonical_smiles"]).HasSubstructMatch(pattern, useChirality=True))
            row["contains_group_core_rule"] = self.current_core_rule
        return pd.DataFrame(rows)

    def feedback_counts(self):
        rows = self.measurement_rows()
        measured = rows[rows.measurement_status.eq("measured")] if not rows.empty else rows
        return {"recorded": len(measured), **{key: int(measured.neighbor_comparison.eq(key).sum())
                if not measured.empty else 0 for key in ("agrees", "disagrees", "no_comparison")}}

    def measured_result_info(self, query_id):
        rows = self.measurement_rows()
        row = _rows(rows[rows.query_id.eq(query_id)])[0]
        fields = ["result_source", "measured_value", "measured_unit", "measured_class", "activity_pchembl",
                  "class_source", "batch", "note", "recorded_at", "neighbor_evidence", "neighbor_comparison",
                  "measurement_status"]
        return {**{field: row.get(field) for field in fields}, "editable": row["result_source"] != "upload",
                "directory": str(Path(self.new_result.output_files["new_molecules_input"]).parent)}

    def _apply_measurement(self, query_id, measurement):
        """Only measured fields of this molecule change; evidence and coordinates are untouched."""
        for table in (self.new_result.data, self.new_result.all_molecules):
            mask = table.query_id.eq(query_id)
            values = {"activity_class": measurement["measured_class"], "measurement_status": "measured",
                      "activity_pchembl": measurement["activity_pchembl"],
                      "user_median": measurement["activity_pchembl"]}
            evidence = table.loc[mask].iloc[0].get("neighbor_evidence")
            values["neighbor_comparison"] = neighbor_comparison(measurement["measured_class"], evidence)
            for field, value in values.items():
                if field not in table:
                    table[field] = None
                table.loc[mask, field] = value

    def _replace_result_records(self, query_id, measurement):
        """Prepare one retained exact record with the same aggregation used by uploads."""
        records = self.new_result.records
        if not records.empty and "query_id" in records:
            records = records[~records.query_id.eq(query_id)].copy()
        if measurement["class_source"] == "user_class":
            self.new_result.records = records
            return
        from .evidence import prepare_generic_activity_records
        row = self.new_result.data[self.new_result.data.query_id.eq(query_id)].iloc[0]
        source = pd.DataFrame({"query_id": [query_id], "smiles": [row.canonical_smiles],
                               "value": [measurement["measured_value"]], "batch": [measurement["batch"]]})
        unit = {"pActivity": "p", "µM": "uM", "μM": "uM"}.get(measurement["measured_unit"], measurement["measured_unit"])
        settings = self.job["options"]
        aggregated, retained, _, _ = prepare_generic_activity_records(source, compound_column="query_id",
            smiles_column="smiles", value_column="value", value_unit=unit,
            document_column="batch" if measurement["batch"] else None,
            lower=settings["lower_threshold"], upper=settings["upper_threshold"],
            conflict_range_threshold=settings["conflict_range_threshold"], confidence=bool(measurement["batch"]),
            offset_display="none", offset_min_shared=3, offset_threshold=1., offset_sign_agreement=.8)
        ignored = {"molecule_id", "molecule_identity_key", "smiles", "canonical_smiles", "activity_class"}
        annotations = {key: value for key, value in aggregated.iloc[0].items() if key not in ignored}
        for table in (self.new_result.data, self.new_result.all_molecules):
            for field, value in annotations.items():
                if field not in table:
                    table[field] = None
                table.loc[table.query_id.eq(query_id), field] = value
        retained["query_id"] = query_id
        retained["source_row"] = row.source_row
        retained["raw_activity_value"] = measurement["measured_value"]
        retained["raw_activity_unit"] = measurement["measured_unit"]
        retained["record_relation"] = "="
        retained["evidence_unit"] = measurement["batch"]
        self.new_result.records = pd.concat([records, retained], ignore_index=True)

    def _restore_unmeasured(self, query_id):
        originals = self.measurement_originals.get(query_id)
        if originals:
            for table, original in zip((self.new_result.data, self.new_result.all_molecules), originals):
                mask = table.query_id.eq(query_id)
                for field in table:
                    value = original.get(field)
                    if value is None and pd.api.types.is_bool_dtype(table[field]):
                        value = False
                    table.loc[mask, field] = value

    def _remember_unmeasured(self, query_id):
        if query_id not in self.measurement_originals:
            self.measurement_originals[query_id] = [table.loc[table.query_id.eq(query_id)].iloc[0].to_dict()
                                                  for table in (self.new_result.data, self.new_result.all_molecules)]

    def add_new(self, payload):
        with self.lock:
            self._manual_ready()
            preview = self.preview_new(payload)
            if preview["status"] != "ok":
                return preview
            old = self.new_result
            ids = set(old.all_molecules.query_id.astype(str)) if old is not None else set()
            if old is not None and len(old.data) >= 1000:
                raise ValueError("At most 1,000 checked new molecules are supported.")
            identity = str(payload.get("query_id") or "").strip()
            if not identity:
                counter = 1
                while f"M{counter:03d}" in ids:
                    counter += 1
                identity = f"M{counter:03d}"
            if identity in ids:
                raise ValueError("New molecule IDs must be unique; this ID already exists.")
            measurement = self._measurement(payload.get("measurement"), optional=True)
            source = pd.DataFrame({"query_id": [identity], "smiles": [str(payload.get("smiles")).strip()],
                                   "label": [str(payload.get("label") or "")]})
            # Reuse the exact file-insertion evidence, grouping and structure pipeline.
            new = insert_new_molecules(self.result, source, id_col="query_id", label_column="label")
            unmeasured = [table.iloc[0].to_dict() for table in (new.data, new.all_molecules)]
            if measurement and measurement["class_source"] == "value":
                # Keep the upload's exact aggregation fields and retained-record schema.
                source["value"] = measurement["measured_value"]
                source["batch"] = measurement["batch"]
                unit = {"pActivity": "p", "µM": "uM", "μM": "uM"}.get(
                    measurement["measured_unit"], measurement["measured_unit"])
                new = insert_new_molecules(self.result, source, id_col="query_id", label_column="label",
                    value_col="value", value_unit=unit,
                    document_column="batch" if measurement["batch"] else None)
            if old is not None:
                occupied = pd.concat([self.result.data, old.data], ignore_index=True)
                settings = self.result.insertion_parameters
                new.data = place_new_molecules(occupied, new.data, grid_rows=settings["grid_rows"],
                                               grid_cols=settings["grid_cols"])
                for field in ("grid_row", "grid_col", "grid_x", "grid_y", "placement_distance"):
                    new.all_molecules[field] = new.data[field].to_numpy()
                source_row = int(old.check.source_row.max()) + 1 if len(old.check) else 0
                for table in (new.data, new.all_molecules, new.check):
                    table["source_row"] = source_row
                for name in ("data", "all_molecules", "check", "neighbors"):
                    setattr(new, name, pd.concat([getattr(old, name), getattr(new, name)], ignore_index=True))
                new.records = pd.concat([old.records, new.records], ignore_index=True)
                new.report, new.output_files = dict(old.report), dict(old.output_files)
            else:
                directory = Path(self.job["output_dir"]) / "new_molecules" / "manual"
                base, suffix = directory, 2
                while directory.exists():
                    directory = base.with_name(f"manual-{suffix}")
                    suffix += 1
                directory.mkdir(parents=True)
                path = directory / "new_molecules_input.csv"
                path.write_bytes(b"query_id,smiles,label\n")
                new.output_files = {"new_molecules_input": str(path)}
                self.query_options = {}
            self.new_result = new
            stamp = datetime.now(timezone.utc).isoformat()
            note = str(payload.get("note") or (measurement or {}).get("note") or "")
            self.manual_entries.append({"query_id": identity, "input_smiles": source.smiles.iloc[0],
                "canonical_smiles": preview["canonical_smiles"], "value": (measurement or {}).get("measured_value"),
                "unit": (measurement or {}).get("measured_unit", ""),
                "measured_class": (measurement or {}).get("measured_class", ""),
                "class_source": (measurement or {}).get("class_source", ""),
                "batch": (measurement or {}).get("batch", ""), "label": source.label.iloc[0],
                "note": note, "entered_at": stamp})
            self.manual_entry_notes[identity] = note
            if measurement:
                self._remember_unmeasured(identity)
                # Restore this entry to its actual unmeasured state on removal, retaining its final cell.
                for baseline, original in zip(unmeasured, self.measurement_originals[identity]):
                    for field in ("source_row", "grid_row", "grid_col", "grid_x", "grid_y", "placement_distance"):
                        baseline[field] = original.get(field)
                self.measurement_originals[identity] = unmeasured
                measurement.update(result_source="manual", recorded_at=stamp, note=note)
                self.measured_results[identity] = measurement
                self._apply_measurement(identity, measurement)
            self._write_manual_outputs()
            return self._manual_response(identity)

    def _editable_query(self, query_id):
        self._manual_ready()
        if self.new_result is None or query_id not in set(self.new_result.data.query_id):
            raise ValueError("Choose a mapped new molecule.")
        if self.measured_result_info(query_id)["result_source"] == "upload":
            raise ValueError("Measured results from the uploaded file are read-only.")

    def record_result(self, payload):
        with self.lock:
            identity = str(payload.get("query_id") or "")
            self._editable_query(identity)
            measurement = self._measurement(payload.get("measurement"))
            self._remember_unmeasured(identity)
            self._restore_unmeasured(identity)
            manual = next((entry for entry in self.manual_entries if entry["query_id"] == identity), None)
            measurement.update(result_source=self.measured_results.get(identity, {}).get("result_source", "recorded"),
                               recorded_at=datetime.now(timezone.utc).isoformat())
            self.measured_results[identity] = measurement
            self._replace_result_records(identity, measurement)
            self._apply_measurement(identity, measurement)
            if manual:
                manual.update(value=measurement["measured_value"], unit=measurement["measured_unit"],
                              **{key: measurement[key] for key in ("measured_class", "class_source", "batch", "note")})
            self._write_manual_outputs()
            return self._manual_response(identity)

    def clear_result(self, payload):
        with self.lock:
            identity = str(payload.get("query_id") or "")
            self._editable_query(identity)
            self._restore_unmeasured(identity)
            self.measurement_originals.pop(identity, None)
            self.measured_results.pop(identity, None)
            if not self.new_result.records.empty and "query_id" in self.new_result.records:
                self.new_result.records = self.new_result.records[~self.new_result.records.query_id.eq(identity)].copy()
            manual = next((entry for entry in self.manual_entries if entry["query_id"] == identity), None)
            if manual:
                manual.update(value=None, unit="", measured_class="", class_source="", batch="",
                              note=self.manual_entry_notes.get(identity, ""))
            self._write_manual_outputs()
            return self._manual_response(identity)

    def _manual_response(self, identity):
        row = _rows(self.measurement_rows().loc[lambda frame: frame.query_id.eq(identity)])[0]
        point = next(point for point in self.points() if point.get("query_id") == identity)
        return {"status": "ok", "query_id": identity, "row": row, "point": point, "job": self.snapshot()}

    def _write_manual_outputs(self):
        """Only the overlay directory is rewritten; original input bytes are retained."""
        result = self.new_result
        directory = Path(result.output_files["new_molecules_input"]).parent
        rows = self.measurement_rows()
        feedback = rows[rows.measurement_status.eq("measured")].reindex(columns=FEEDBACK_COLUMNS)
        exported = result.all_molecules.copy()
        for field in rows:
            if field not in exported:
                exported[field] = rows[field].to_numpy()
        exported["contains_group_core_rule"] = self.job["options"]["core_rule"]
        tables = {"new_molecules": exported, "new_molecules_check": result.check,
                  "new_molecules_neighbors": result.neighbors, "neighbor_evidence_feedback": feedback}
        if self.manual_entries:
            tables["manual_entries"] = pd.DataFrame(self.manual_entries, columns=MANUAL_COLUMNS)
        for name, table in tables.items():
            path = directory / f"{name}.csv"
            table.to_csv(path, index=False, encoding="utf-8", lineterminator="\n")
            result.output_files[name] = str(path)
        counts = self.feedback_counts()
        result.report.update(manual_count=len(self.manual_entries), feedback_counts=counts,
                             feedback_file="neighbor_evidence_feedback.csv", mapped=len(result.data),
                             input_compounds=len(result.all_molecules),
                             neighbor_evidence_counts={key: int(result.data.neighbor_evidence.eq(key).sum())
                                                       for key in NEIGHBOR_EVIDENCE_CATEGORIES})
        for group in result.structure_groups:
            count = int(result.data.structure_group.eq(group["structure_group"]).sum())
            group.update(n_query_members=count, n_new_members=count)
        settings = self.result.insertion_parameters
        svg = directory / "map_with_new_molecules.svg"
        _new_molecules_svg(self.result, result.data, svg, rows=settings["grid_rows"], cols=settings["grid_cols"],
                           lower=settings["lower_threshold"], upper=settings["upper_threshold"], show_feedback=True)
        result.output_files["map_with_new_molecules"] = str(svg)
        manifest_path = directory / "new_molecules_manifest.json"
        if manifest_path.exists():
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        else:
            reference_manifest = self.result.output_files.get("run_manifest")
            manifest = {"input_filename": "manual", "reference_run_directory": self.job["output_dir"],
                        "reference_manifest_sha256": hashlib.sha256(Path(reference_manifest).read_bytes()).hexdigest()
                        if reference_manifest else None, "parameters": dict(settings)}
        raw = Path(result.output_files["new_molecules_input"]).read_bytes()
        manifest.update(input_sha256=hashlib.sha256(raw).hexdigest(), report=result.report,
                        software_source_sha256=source_tree_digest()["source_sha256"],
                        feedback={"filename": "neighbor_evidence_feedback.csv", **counts},
                        updated_at=datetime.now(timezone.utc).isoformat())
        if self.manual_entries:
            path = Path(result.output_files["manual_entries"])
            manifest["manual_entries"] = {"filename": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                                          "rows": len(self.manual_entries)}
        manifest["output_digests"] = file_digests({key: path for key, path in result.output_files.items()
            if key not in {"new_molecules_manifest", "new_molecules_results_zip"}})
        manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False),
                                 encoding="utf-8", newline="\n")
        result.output_files["new_molecules_manifest"] = str(manifest_path)
        archive = directory / "new_molecules_results.zip"
        with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as bundle:
            for key, path in result.output_files.items():
                if key != "new_molecules_results_zip":
                    bundle.write(path, Path(path).name)
        result.output_files["new_molecules_results_zip"] = str(archive)
        self.files = dict(self.reference_files)
        self.files.update({Path(path).name: Path(path) for path in result.output_files.values()})
        self.job.update(report={**self.job["report"], "new_molecules": result.report, "query_molecules": result.report},
                        files=list(self.files), molecules=len(self.result.data) + len(result.data))
        self.group_detail_cache = {}
        self.source_index_cache.pop("new", None)
        self.save_session(overlay=True)

    def feedback_csv(self):
        with self.lock:
            if self.new_result is None:
                raise ValueError("Insert new molecules first.")
            # Reading an existing feedback export must not rewrite its manifest/ZIP.
            if "neighbor_evidence_feedback" not in self.new_result.output_files:
                self._write_manual_outputs()
            return Path(self.new_result.output_files["neighbor_evidence_feedback"]).read_bytes()


# Snapshots are private local caches; none of these files enter exported manifests or ZIPs.
class SessionSnapshotError(ValueError):
    def __init__(self, code, details="", files=()):
        self.code, self.details, self.files = code, str(details), list(files)
        super().__init__(f"{code}: {details}. Please regenerate the map.")

    def payload(self):
        return {"error": str(self), "code": self.code, "details": self.details, "files": self.files}


class SessionUnpickler(pickle.Unpickler):
    """Allow only the constructors used by our tables/arrays/result containers.

    A module-prefix allowlist would also allow numpy.load or pandas.read_pickle.
    Exact globals deliberately exclude such I/O and executable callables.
    """
    allowed = {
        "builtins": {"set", "frozenset", "slice", "complex", "range", "bytearray"},
        "collections": {"OrderedDict", "defaultdict"},
        "datetime": {"datetime", "date", "time", "timedelta", "timezone"},
        "numpy": {"ndarray", "dtype"},
        "numpy.core.multiarray": {"_reconstruct", "scalar"},
        "numpy._core.multiarray": {"_reconstruct", "scalar"},
        "numpy.core.numeric": {"_frombuffer"},
        "numpy._core.numeric": {"_frombuffer"},
        "pandas": {"DataFrame", "Series", "Index", "RangeIndex", "StringDtype"},
        "pandas.arrays": {"StringArray"},
        "pandas.core.frame": {"DataFrame"},
        "pandas.core.series": {"Series"},
        "pandas.core.internals.managers": {"BlockManager", "SingleBlockManager"},
        "pandas._libs.internals": {"_unpickle_block"},
        "pandas._libs.arrays": {"__pyx_unpickle_NDArrayBacked"},
        "pandas.core.indexes.base": {"_new_Index", "Index"},
        "pandas.core.indexes.range": {"RangeIndex"},
        "pandas.core.indexes.multi": {"MultiIndex"},
        "pandas.core.arrays.string_": {"StringArray", "StringDtype"},
        "pandas.core.arrays.boolean": {"BooleanArray", "BooleanDtype"},
        "pandas.core.arrays.integer": {"IntegerArray", "Int64Dtype"},
        "pandas._libs.missing": {"NAType"},
        "pandas._libs.tslibs.nattype": {"__nat_unpickle"},
        "pandas._libs.tslibs.timestamps": {"_unpickle_timestamp"},
        "chemgridmap.pipeline": {"GridMapResult"},
        "chemgridmap.chembl": {"ChemblCurationResult"},
        "chemgridmap.queries": {"NewMoleculesResult"},
    }

    def find_class(self, module, name):
        if name not in self.allowed.get(module, set()):
            raise SessionSnapshotError("session_unsafe", f"Disallowed snapshot global: {module}.{name}")
        return super().find_class(module, name)


def _session_path(root, relative):
    """Reject absolute paths, traversal and symlinks leaving the selected run."""
    path = Path(str(relative))
    if path.is_absolute() or ".." in path.parts:
        raise SessionSnapshotError("session_invalid", "Invalid relative snapshot path", [str(relative)])
    resolved = (root / path).resolve()
    if not resolved.is_relative_to(root.resolve()):
        raise SessionSnapshotError("session_invalid", "Snapshot path leaves the run", [str(relative)])
    return resolved


def _session_paths(value, root, *, restore=False):
    """Relocate only explicit paths, leaving DataFrame values and arrays untouched."""
    if restore and isinstance(value, dict) and set(value) == {"__session_path__", "kind"}:
        path = _session_path(root, value["__session_path__"])
        return path if value["kind"] == "path" else str(path)
    if not restore and isinstance(value, (str, Path)):
        if isinstance(value, Path) or str(value).startswith(str(root) + "/") or str(value) == str(root):
            relative = Path(value).resolve().relative_to(root.resolve()).as_posix()
            return {"__session_path__": relative, "kind": "path" if isinstance(value, Path) else "str"}
    if is_dataclass(value) and not isinstance(value, type):
        result = copy.copy(value)
        for field in fields(value):
            setattr(result, field.name, _session_paths(getattr(value, field.name), root, restore=restore))
        return result
    if isinstance(value, dict):
        return {key: _session_paths(item, root, restore=restore) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        items = [_session_paths(item, root, restore=restore) for item in value]
        return tuple(items) if isinstance(value, tuple) else items
    return value


def _session_versions():
    return {"python": platform.python_version(), "pandas": pd.__version__,
            "numpy": np.__version__, "rdkit": rdBase.rdkitVersion}


def _session_sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class SavedSessionMixin:
    _reference_fields = ("filename", "result", "curation", "generic_records", "activity_confidence",
                         "possible_sources", "source_frame", "source_filename")
    _runtime_fields = ("files", "reference_files", "reference_report", "job", "core_rule_cache",
                       "group_core_cache", "current_core_rule", "group_cache")
    _overlay_fields = ("new_result", "query_options", "query_frame", "query_raw", "query_filename",
                       "manual_entries", "manual_entry_notes", "measured_results", "measurement_originals")

    def _write_session_blob(self, directory, name, state):
        root = Path(self.job["output_dir"]).resolve()
        encoded = pickle.dumps(_session_paths(state, root), protocol=5)
        path = directory / name
        temporary = path.with_suffix(path.suffix + ".tmp")
        temporary.write_bytes(gzip.compress(encoded, compresslevel=1, mtime=0))
        temporary.replace(path)
        return path

    def _session_manifest(self, directory, blobs, artifacts, **extra):
        root = Path(self.job["output_dir"]).resolve()
        run_manifest = Path(self.result.output_files["run_manifest"])
        manifest = {"format_version": 1, "source_tree_digest": source_tree_digest(),
                    "versions": _session_versions(),
                    "created_at": datetime.fromtimestamp(self.job["started"], timezone.utc).isoformat(),
                    "run_manifest": {"file": run_manifest.relative_to(root).as_posix(),
                                     "sha256": _session_sha(run_manifest)},
                    "files": {path.relative_to(root).as_posix(): _session_sha(path) for path in blobs},
                    "artifacts": {Path(path).relative_to(root).as_posix(): _session_sha(path)
                                  for path in artifacts}, **extra}
        path = directory / "session_manifest.json"
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
        temporary.replace(path)
        return path

    def save_session(self, *, initial=False, overlay=False):
        """Persist completed state without touching any existing scientific export."""
        if self.result is None or self.job.get("state") != "done":
            return
        started = time.perf_counter()
        self._ensure_core_state()
        root = Path(self.job["output_dir"]).resolve()
        directory = root / "session"
        directory.mkdir(exist_ok=True)
        reference_path = directory / "reference.pkl.gz"
        if initial or not reference_path.exists():
            self._write_session_blob(directory, reference_path.name,
                                     {key: getattr(self, key) for key in self._reference_fields})
        active = None
        overlay_manifest = None
        if self.new_result is not None:
            overlay_root = Path(self.new_result.output_files["new_molecules_input"]).parent
            active = overlay_root.relative_to(root).as_posix()
            overlay_directory = overlay_root / "session"
            overlay_manifest = overlay_directory / "session_manifest.json"
            if overlay or not overlay_manifest.exists():
                overlay_directory.mkdir(exist_ok=True)
                blob = self._write_session_blob(overlay_directory, "overlay.pkl.gz",
                                               {key: getattr(self, key) for key in self._overlay_fields})
                self._session_manifest(overlay_directory, [blob], self.new_result.output_files.values())
        state = {key: getattr(self, key) for key in self._runtime_fields}
        state["active_set"] = active
        runtime_path = self._write_session_blob(directory, "state.pkl.gz", state)
        blobs = [reference_path, runtime_path]
        if overlay_manifest:
            blobs.append(overlay_manifest)
        # On-demand inspection archives are regenerated by inspect(), so they are
        # not part of the durable map/cache integrity contract.
        artifacts = [path for path in self.reference_files.values()
                     if not Path(path).name.startswith("audit_")]
        info = self.reference_map_info()
        self._session_manifest(directory, blobs, artifacts, active_set=active,
            summary={"filename": self.filename, "reference_map": {**info, "label": info["description"]},
                     "molecules": len(self.result.data), "core_rule": self.job["options"]["core_rule"],
                     "new_molecules": len(self.new_result.data) if self.new_result is not None else 0})
        self.session_write_seconds = time.perf_counter() - started

    @staticmethod
    def _validate_session_manifest(root, directory):
        try:
            manifest = json.loads((directory / "session_manifest.json").read_text(encoding="utf-8"))
            if (manifest["format_version"] != 1 or manifest["versions"] != _session_versions()
                    or manifest["source_tree_digest"] != source_tree_digest()):
                raise SessionSnapshotError("session_version", "Software source or dependency versions differ")
            binding = manifest["run_manifest"]
            manifest_path = _session_path(root, binding["file"])
            expected = {binding["file"]: binding["sha256"], **manifest["files"], **manifest["artifacts"]}
            changed = [name for name, digest in expected.items()
                       if not _session_path(root, name).is_file() or _session_sha(_session_path(root, name)) != digest]
            if changed:
                raise SessionSnapshotError("session_modified", "Files were changed or deleted", changed)
            run_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            output_digests = dict(run_manifest["output_digests"])
            output_digests.update(run_manifest.get("source", {}).get("curation_outputs", {}))
            changed = [item["file"] for item in output_digests.values()
                       if not _session_path(root, item["file"]).is_file()
                       or _session_sha(_session_path(root, item["file"])) != item["sha256"]]
            if changed:
                raise SessionSnapshotError("session_modified", "Run outputs were changed or deleted", changed)
            return manifest
        except SessionSnapshotError:
            raise
        except FileNotFoundError as exc:
            raise SessionSnapshotError("session_missing", "No complete session snapshot") from exc
        except (OSError, ValueError, KeyError, TypeError) as exc:
            raise SessionSnapshotError("session_invalid", "Invalid session manifest") from exc

    def _validate_saved_run(self, root):
        manifest = self._validate_session_manifest(root, root / "session")
        summary_fields = {"filename", "reference_map", "molecules", "core_rule", "new_molecules"}
        if (not isinstance(manifest.get("created_at"), str)
                or not isinstance(manifest.get("summary"), dict)
                or not summary_fields <= set(manifest["summary"])):
            raise SessionSnapshotError("session_invalid", "Invalid previous-result metadata")
        required = {"session/reference.pkl.gz", "session/state.pkl.gz"}
        if manifest.get("active_set"):
            directory = _session_path(root, manifest["active_set"]) / "session"
            required.add(directory.relative_to(root).as_posix() + "/session_manifest.json")
            overlay = self._validate_session_manifest(root, directory)
            if directory.relative_to(root).as_posix() + "/overlay.pkl.gz" not in overlay["files"]:
                raise SessionSnapshotError("session_invalid", "Missing overlay checksum")
        if not required <= set(manifest["files"]):
            raise SessionSnapshotError("session_invalid", "Missing snapshot checksum")
        return manifest

    def previous_runs(self):
        runs, unavailable = [], {"missing_snapshot": 0, "incompatible": 0, "modified": 0}
        if self.output_root.exists():
            for root in self.output_root.iterdir():
                if not root.is_dir() or root.is_symlink():
                    continue
                try:
                    manifest = self._validate_saved_run(root)
                    runs.append({"run_name": root.name, "created_at": manifest["created_at"],
                                 **manifest["summary"]})
                except SessionSnapshotError as exc:
                    key = ("missing_snapshot" if exc.code == "session_missing" else
                           "incompatible" if exc.code == "session_version" else "modified")
                    unavailable[key] += 1
        runs.sort(key=lambda row: (row["created_at"], row["run_name"]), reverse=True)
        self._listed_run_names = {row["run_name"] for row in runs[:20]}
        return {"runs": runs[:20], "unavailable": unavailable}

    @staticmethod
    def _read_session_blob(root, relative):
        try:
            with gzip.open(_session_path(root, relative), "rb") as stream:
                state = SessionUnpickler(stream).load()
            if not isinstance(state, dict):
                raise ValueError("Expected a session state mapping")
            return _session_paths(state, root, restore=True)
        except SessionSnapshotError:
            raise
        except Exception as exc:
            raise SessionSnapshotError("session_invalid", "Cannot decode the session snapshot") from exc

    def open_run(self, run_name):
        # Build/open/cache operations cannot replace one another mid-restore.
        with self.core_compute_lock, self.lock:
            if self.job.get("state") == "running":
                raise SessionSnapshotError("map_running", "Wait for the running map")
            if (not isinstance(run_name, str) or run_name in {"", ".", ".."}
                    or Path(run_name).name != run_name or "/" in run_name or "\\" in run_name):
                raise SessionSnapshotError("run_not_listed", "Choose a run from the previous-results list")
            root = self.output_root / run_name
            if not root.is_dir() or root.is_symlink():
                raise SessionSnapshotError("run_not_listed", "Choose a run from the previous-results list")
            # Validate before checking membership so modified files receive a useful reason.
            manifest = self._validate_saved_run(root)
            if getattr(self, "_listed_run_names", None) is None:
                self.previous_runs()
            if run_name not in self._listed_run_names:
                raise SessionSnapshotError("run_not_listed", "Choose a run from the previous-results list")
            restored = self._read_session_blob(root, "session/reference.pkl.gz")
            runtime = self._read_session_blob(root, "session/state.pkl.gz")
            restored.update(runtime)
            if "active_set" not in restored or restored.pop("active_set") != manifest.get("active_set"):
                raise SessionSnapshotError("session_invalid", "Active new-molecule set differs")
            overlay_state = {"new_result": None, "query_options": None, "query_frame": None,
                             "query_raw": b"", "query_filename": "", "manual_entries": [],
                             "manual_entry_notes": {}, "measured_results": {}, "measurement_originals": {}}
            if manifest.get("active_set"):
                overlay_state = self._read_session_blob(root, manifest["active_set"] + "/session/overlay.pkl.gz")
            restored.update(overlay_state)
            expected = set(self._reference_fields + self._runtime_fields + self._overlay_fields)
            if (set(restored) != expected or not isinstance(restored.get("job"), dict)
                    or restored["job"].get("state") != "done"):
                raise SessionSnapshotError("session_invalid", "Incomplete completed session")
            try:
                from .pipeline import GridMapResult
                if not isinstance(restored["result"], GridMapResult):
                    raise ValueError("Missing reference result")
                candidate = AppSession(self.output_root)
                candidate.__dict__.update(restored)
                candidate.raw = (root / "input.csv").read_bytes()
                candidate.frame = _read_input(candidate.raw)
                candidate.core_result = candidate.result
                # Exercise the same view required by the browser before committing any live state.
                candidate.session_info()
                summary = {"run_name": run_name,
                           "new_set": Path(manifest["active_set"]).name if manifest.get("active_set") else None,
                           "total": len(candidate.new_result.data) if candidate.new_result is not None else 0,
                           "manual": len(candidate.manual_entries), "recorded": len(candidate.measured_results)}
                response = {**candidate.snapshot(), "restored": summary}
            except (ValueError, TypeError, AttributeError, KeyError, OSError) as exc:
                raise SessionSnapshotError("session_invalid", "Incomplete session state") from exc
            # Preserve this server's token and locks; restoration itself writes no files.
            restored.update(raw=candidate.raw, frame=candidate.frame, core_result=candidate.result,
                            pending_query=None, group_detail_cache={}, source_index_cache={})
            self.__dict__.update(restored)
            return response


class AppSession(SavedSessionMixin, ManualMoleculesMixin):
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
        self.query_options = None
        self.new_result = None
        self._reset_manual()
        self.pending_query = None
        self.reference_files = {}
        self.reference_report = {}
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
        self.core_compute_lock = threading.Lock()
        self.core_progress = {"running": False, "rule": None, "done": 0, "total": 0}
        self.core_result = None
        self.core_rule_cache = {}
        self.group_core_cache = {}
        self.current_core_rule = None
        self.source_index_cache = {}
        self.job = {"state": "empty", "message": "Choose a CSV or open the example."}

    def load(self, content, filename):
        with self.lock:
            if self.job["state"] == "running":
                raise ValueError("A map is running. Wait for it to finish before replacing the input.")
            frame = _read_input(content)
            self.frame, self.raw = frame, content
            self.filename = Path(filename.replace("\\", "/")).name
            self.query_frame, self.query_raw, self.query_filename = None, b"", ""
            self.query_options, self.new_result = None, None
            self._reset_manual()
            self.pending_query = None
            self.reference_files, self.reference_report = {}, {}
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
            filename = Path(filename.replace("\\", "/")).name
            if self.result is not None:
                self.pending_query = (frame, content, filename)
            else:
                self.query_frame, self.query_raw, self.query_filename = frame, content, filename
                self.query_options = None
            return self.query_input_info()

    def clear_query(self):
        with self.lock:
            if self.job["state"] == "running":
                raise ValueError("A map is running. Wait for it to finish before replacing the input.")
            self.query_frame, self.query_raw, self.query_filename = None, b"", ""
            self.query_options, self.new_result = None, None
            self._reset_manual()
            self.pending_query = None
            if self.result is not None:
                self.files = dict(self.reference_files)
                self.job.update(report=self.reference_report, files=list(self.files),
                                molecules=len(self.result.data))
            self.group_detail_cache = {}
            self.save_session()
            return {"query_input": None}

    def new_molecules_info(self):
        """Current overlay survives browser refresh; original run files stay separate."""
        if self.new_result is None:
            return None
        directory = str(Path(self.new_result.output_files["new_molecules_input"]).parent)
        return {"rows": json_rows(self.measurement_rows()),
                "report": {**self.new_result.report, "feedback_counts": self.feedback_counts(),
                           "manual_count": len(self.manual_entries)}, "options": self.query_options,
                "input": self.query_input_info(pending=False), "directory": directory,
                "replacement": {"total": len(self.new_result.data), "manual": len(self.manual_entries),
                                "recorded": len(self.measured_results), "directory": directory},
                "files": [Path(path).name for path in self.new_result.output_files.values()]}

    def points(self):
        if self.result is None:
            raise ValueError("Build a map first.")
        if self.new_result is None:
            return point_rows(self.result.data)
        reference = self.result.data.assign(measurement_status="measured", is_new_molecule=False)
        new = self.new_result.data.assign(is_new_molecule=True)
        return point_rows(pd.concat([reference, new], ignore_index=True))

    def groups(self):
        self._ensure_core_state()
        groups = self.core_rule_cache[self.current_core_rule]
        if self.new_result is None:
            return groups
        counts = {g["structure_group"]: {key: g[key] for key in ("n_query_members", "n_new_members")
                                          if key in g}
                  for g in self.new_result.structure_groups}
        return [{**group, **counts.get(group["structure_group"], {})} for group in groups]

    def _ensure_core_state(self):
        if self.result is None:
            raise ValueError("Build a map first.")
        if self.core_result is not self.result:
            self.core_result = self.result
            self.current_core_rule = self.job["options"]["core_rule"]
            self.core_rule_cache = {self.current_core_rule: self.result.structure_groups}
            self.group_core_cache = {}

    @staticmethod
    def _core_members(result, summary):
        members = result.data[result.data.structure_group.eq(summary["structure_group"])].copy()
        members = members.sort_values(["canonical_smiles", "comparison_id"])
        center = members.comparison_id.eq(summary["representative_id"])
        # Both Butina and external groups use center first, then canonical input order.
        return pd.concat([members[center], members[~center]])

    @staticmethod
    def _calculate_group_core(result, summary, rule, target, timeout=1):
        members = AppSession._core_members(result, summary)
        mols = [Chem.MolFromSmiles(smiles) for smiles in members.canonical_smiles]
        core, metrics = summarize_group_core(
            mols, members.comparison_id.tolist(), core_rule=rule,
            core_min_coverage=target, fmcs_timeout=timeout, return_metrics=True)
        identity = {key: summary[key] for key in (
            "structure_group", "size", "color", "representative_id", "representative_smiles",
            "representative_role")}
        return {**identity, **core, "minimum_center_similarity": summary["minimum_center_similarity"]}, metrics

    @staticmethod
    def _generic_scaffold(result, summary):
        """Supplementary generic Murcko scaffold of a fixed group, with a drawing."""
        members = AppSession._core_members(result, summary)
        generic = generic_scaffold_summary([Chem.MolFromSmiles(s) for s in members.canonical_smiles])
        smiles = generic["generic_scaffold"]
        generic["generic_svg"] = (
            '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 180 180">'
            + (_svg_molecule_content(smiles, 180) or '') + '</svg>') if smiles else ''
        return generic

    def _register_core_file(self, path):
        name = path.relative_to(Path(self.job["output_dir"])).as_posix()
        self.files[name] = path
        self.reference_files[name] = path
        self.job["files"] = list(self.files)
        return name

    @staticmethod
    def _core_table_rows(groups, metrics):
        fields = ("structure_group", "size", "core_rule", "core_smiles", "core_smiles_attached",
                  "core_heavy_atoms", "core_coverage_count", "core_coverage_fraction",
                  "core_coverage_target", "core_target_met")
        rows = []
        for group, timing in zip(groups, metrics):
            row = {key: group.get(key) for key in fields}
            row.update(candidates=timing["candidates"], seconds=timing["seconds"])
            for key in ("brics_candidates", "brics_single_pass_members",
                        "recap_candidates", "recap_skipped_members"):
                if key in group:
                    row[key] = group[key]
            rows.append(row)
        return rows

    def set_core_rule(self, rule):
        """Recompute only fixed-group cores; never mutate an original map or its files."""
        with self.core_compute_lock:
            try:
                return self._set_core_rule(rule)
            finally:
                self.core_progress["running"] = False

    def _set_core_rule(self, rule):
        # The public entry point owns core_compute_lock throughout this operation.
        with self.lock:
            self._ensure_core_state()
            built_rule = self.job["options"]["core_rule"]
            allowed = {"most_frequent", "coverage", "coverage_brics", "coverage_recap"}
            if built_rule == "fmcs":
                allowed.add("fmcs")
            if not isinstance(rule, str) or rule not in allowed:
                raise ValueError("Unknown or unavailable core rule.")
            if not self.result.structure_groups:
                raise ValueError("This map has no structure groups.")
            cached = rule in self.core_rule_cache
            result, options = self.result, dict(self.job["options"])
        if not cached:
            # Chemistry runs outside the session lock, so point/inspect requests remain available.
            started = time.perf_counter()
            self.core_progress = {"running": True, "rule": rule, "done": 0,
                                  "total": len(result.structure_groups)}
            try:
                members_by_group = group_core_members(result.data, result.structure_groups)
                calculated = summarize_group_cores(
                    [(members.canonical_smiles.tolist(), members.comparison_id.tolist())
                     for members in members_by_group], core_rule=rule,
                    core_min_coverage=options["core_min_coverage"],
                    fmcs_timeout=options["fmcs_timeout"], jobs=browser_core_jobs(),
                    progress=lambda done, total: self.core_progress.update(done=done, total=total))
                identity_fields = ("structure_group", "size", "color", "representative_id",
                                   "representative_smiles", "representative_role", "minimum_center_similarity")
                summaries = [{**{key: group[key] for key in identity_fields}, **core}
                             for group, (core, _) in zip(result.structure_groups, calculated)]
                metrics = [metrics for _, metrics in calculated]
            except Exception:
                self.core_progress["running"] = False
                raise
            elapsed = time.perf_counter() - started
            hashes = {}
            for group, members in zip(result.structure_groups, members_by_group):
                members = members.sort_values(["canonical_smiles", "comparison_id"])
                hashes[group["structure_group"]] = hashlib.sha256(
                    members.canonical_smiles.to_csv(index=False, lineterminator="\n").encode("utf-8")).hexdigest()
            with self.lock:
                if self.result is not result:
                    raise ValueError("The map changed during core calculation. Try again.")
                directory = Path(self.job["output_dir"]) / "core_rules"
                directory.mkdir(exist_ok=True)
                stem = Path(result.output_files["structure_groups"]).name.removesuffix("_structure_groups.json")
                path = directory / f"{stem}_cores_{rule}.csv"
                pd.DataFrame(self._core_table_rows(summaries, metrics)).to_csv(
                    path, index=False, encoding="utf-8", lineterminator="\n")
                self._register_core_file(path)
                path = directory / f"{stem}_cores_{rule}.json"
                metadata = {"core_rule": rule, "parameters": {"core_min_coverage": options["core_min_coverage"]},
                            "rdkit_version": rdBase.rdkitVersion, "structure_sha256": hashes,
                            "run_manifest_sha256": hashlib.sha256(
                                Path(result.output_files["run_manifest"]).read_bytes()).hexdigest(),
                            "seconds": elapsed}
                path.write_text(json.dumps(metadata, indent=2, allow_nan=False), encoding="utf-8", newline="\n")
                self._register_core_file(path)
                self.core_rule_cache[rule] = summaries
        with self.lock:
            if self.result is not result:
                raise ValueError("The map changed during core calculation. Try again.")
            self.current_core_rule = rule
            self.group_detail_cache = {}
            self.save_session()
            self.core_progress = {"running": False, "rule": rule,
                                  "done": len(result.structure_groups), "total": len(result.structure_groups)}
            return {"rule": rule, "build_rule": built_rule, "groups": self.groups(),
                    "files": list(self.files), "cached": cached}

    def group_cores(self, group_id):
        """Compare five core rules for one fixed group; FMCS has a ten-second budget.

        The generic Murcko scaffold is added as supplementary information, not as a rule.
        """
        with self.core_compute_lock:
            with self.lock:
                self._ensure_core_state()
                if group_id in self.group_core_cache:
                    return {**self.group_core_cache[group_id], "rule": self.current_core_rule}
                summary = next((g for g in self.result.structure_groups
                                if g["structure_group"] == group_id), None)
                if summary is None:
                    raise ValueError("Unknown structure group.")
                result = self.result
                target = self.job["options"]["core_min_coverage"]
            rows, cores = [], []
            for rule in ("most_frequent", "coverage", "coverage_brics", "coverage_recap", "fmcs"):
                group, metrics = self._calculate_group_core(result, summary, rule, target, timeout=10)
                row = self._core_table_rows([group], [metrics])[0]
                row["fmcs_timed_out"] = group["fmcs_timed_out"]
                rows.append(row)
                cores.append({**group, **metrics, "core_svg": _group_drawing(group, "core_smiles")})
            generic = self._generic_scaffold(result, summary)
            rows.append({"structure_group": summary["structure_group"], "size": generic["size"],
                         "core_rule": "generic_murcko_supplementary",
                         "core_smiles": generic["generic_scaffold"], "core_smiles_attached": "",
                         "core_heavy_atoms": generic["generic_scaffold_heavy_atoms"],
                         "core_coverage_count": generic["generic_scaffold_count"],
                         "core_coverage_fraction": generic["generic_scaffold_fraction"],
                         "core_coverage_target": None, "core_target_met": None})
            with self.lock:
                if self.result is not result:
                    raise ValueError("The map changed during core calculation. Try again.")
                directory = Path(self.job["output_dir"]) / "core_rules"
                directory.mkdir(exist_ok=True)
                name = "group_" + hashlib.sha256(group_id.encode()).hexdigest()[:12] + "_core_comparison.csv"
                path = directory / name
                pd.DataFrame(rows).to_csv(path, index=False, encoding="utf-8", lineterminator="\n")
                payload = {"group": group_id, "cores": cores, "generic": generic,
                           "file": self._register_core_file(path), "fmcs_timeout": 10}
                self.group_core_cache[group_id] = payload
                self.save_session()
                return {**payload, "rule": self.current_core_rule}

    def insert_new(self, options=None):
        with self.lock:
            if self.result is None:
                raise ValueError("Build a reference map first.")
            if self.job["state"] == "running":
                raise ValueError("A map is running. Wait for it to finish before inserting new molecules.")
            if self.job["options"]["representation"] != "morgan":
                raise ValueError("New molecules require a reference map arranged using Morgan fingerprints.")
            frame, raw, filename = self.pending_query or (self.query_frame, self.query_raw, self.query_filename)
            if frame is None:
                raise ValueError("Upload new molecules first.")
            aliases = {"activity_col": "value_col", "unit": "value_unit", "unit_col": "unit_column",
                       "relation_col": "relation_column", "document_col": "document_column",
                       "label_col": "label_column"}
            allowed = {"smiles_col", "id_col", "value_col", "value_unit", "unit_column",
                       "relation_column", "document_column", "label_column"}
            chosen = {aliases.get(key, key): value for key, value in (options or {}).items()}
            if set(chosen) - allowed:
                raise ValueError("Unknown new-molecule option.")
            chosen = {key: value for key, value in chosen.items() if value not in (None, "")}
            for key, value in chosen.items():
                if key != "value_unit" and (not isinstance(value, str) or value not in frame):
                    raise ValueError(f"The selected {key} column was not found.")
            build_options = self.job["options"]
            policy = (self.curation.report["parameters"]["structure_policy"]
                      if self.curation is not None else "canonical")
            inserted = insert_new_molecules(
                self.result, frame, output_dir=Path(self.job["output_dir"]),
                input_bytes=raw, input_filename=filename,
                structure_policy=policy, lower_threshold=build_options["lower_threshold"],
                upper_threshold=build_options["upper_threshold"],
                conflict_range_threshold=build_options["conflict_range_threshold"],
                grouping=build_options["grouping"], group_similarity=build_options["group_similarity"],
                reference_target_id=self.reference_map_info()["target_id"],
                **chosen)
            self.new_result, self.query_options = inserted, chosen
            self._reset_manual()
            self.query_frame, self.query_raw, self.query_filename = frame, raw, filename
            self.pending_query = None
            self.files = dict(self.reference_files)
            self.files.update({Path(path).name: Path(path) for path in inserted.output_files.values()})
            matches = self._query_reference_matches(self.result, self.curation, self.generic_records,
                                                    build_options, checks=inserted.check)
            if chosen.get("value_col"):
                rows = {row["query_id"]: row for row in json_rows(inserted.all_molecules)}
                for match in matches:
                    own = rows[match["query_id"]]
                    own_value = own.get("user_median")
                    ref_value = own.get("reference_median")
                    match.update(own_median=own_value, reference_median=ref_value,
                                 median_difference=(own_value - ref_value
                                                    if own_value is not None and ref_value is not None else None),
                                 own_n_records=own.get("n_records"))
            differences = [row["median_difference"] for row in matches
                           if row.get("median_difference") is not None]
            report = {**self.reference_report, "query_molecules": inserted.report,
                      "new_molecules": inserted.report, "query_reference_matches": matches}
            if differences:
                report["query_reference_differences"] = {
                    "median": float(pd.Series(differences).median()),
                    "min": min(differences), "max": max(differences), "count": len(differences)}
            self.group_detail_cache = {}
            self.job.update(report=report, files=list(self.files), molecules=len(self.result.data) + len(inserted.data))
            self.save_session(overlay=True)
            return self.snapshot()

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
        self.new_result = None
        self._reset_manual()
        self.reference_files, self.reference_report = {}, {}
        self.generic_records = None
        self.activity_confidence = None
        self.possible_sources = None
        self.group_cache, self.group_detail_cache = {}, {}
        self.job = {"state": "ready" if self.frame is not None else "empty",
                    "message": "Input updated. Review the options, then build the map."}

    def query_input_info(self, *, pending=True):
        frame, raw, filename = (self.pending_query if pending and self.pending_query else
                                (self.query_frame, self.query_raw, self.query_filename))
        if frame is None:
            return None
        columns = frame.columns.tolist()
        available = {str(column).strip().lower(): column for column in columns}
        smiles = next((available[name] for name in (
            "smiles", "canonical_smiles", "canonical smiles", "representation_smiles"
        ) if name in available), None)
        return {"filename": filename, "rows": len(frame), "columns": columns,
                "detected": {"smiles": smiles}, "sha256": hashlib.sha256(raw).hexdigest()}

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
            # Legacy Python callers may still supply selectors at build time. They
            # are overlay settings and must not enter the reference run's files.
            legacy_query = {key.removeprefix("query_"): options.pop(key)
                            for key in ("query_smiles_col", "query_id_col") if options.get(key)}
            for key in ("query_smiles_col", "query_id_col"):
                options.pop(key, None)
            if legacy_query:
                self.query_options = legacy_query
            if options.get("unit_column") and options["unit_column"] not in self.frame:
                raise ValueError("The selected unit column was not found.")
            if options.get("reference_column") and options["reference_column"] not in self.frame:
                raise ValueError("The selected reference column was not found.")
            # Rebuilding restores the committed group, not a still-unconfirmed upload.
            self.pending_query = None
            self.result, self.curation, self.files = None, None, {}
            self.new_result = None
            self._reset_manual()
            self.reference_files, self.reference_report = {}, {}
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
                jobs=browser_core_jobs(),
                confidence_threshold=options["support_level"],
                offset_display=options["offset_display"],
                offset_min_shared=options["offset_min_shared"],
                offset_threshold=options["offset_min_difference"],
                offset_sign_agreement=options["offset_sign_agreement"],
                map_color_modes=color_modes, map_color_formats=["svg"],
                value_unit=value_unit if options["input_format"] == "molecule" and not options["compound_column"] else "p",
                unit_column=unit_column if options["input_format"] == "molecule" and not options["compound_column"] else None,
            )
            # Generic record aggregation happens above the pipeline; retain its
            # actual span threshold for later insertion without changing exports.
            result.insertion_parameters["conflict_range_threshold"] = options["conflict_range_threshold"]
            result.output_files.update(curation_files)
            files = {Path(path).name: Path(path) for path in result.output_files.values()}
            files.update({name: directory / name for name in ("input.csv", "app_options.json")})
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
            possible_sources = None
            if generic_records is not None and self.source_frame is not None:
                self.update(message="Matching structures to the source lookup file...")
                possible_sources, lookup_summary = match_source_lookup(
                    data, generic_records, self.source_frame)
                report = {**report, "source_lookup": {
                    "file": self.source_filename, **lookup_summary}}
            with self.lock:
                self.result, self.curation, self.files = result, curation, files
                self.reference_files, self.reference_report = dict(files), report
                self.generic_records = generic_records
                self.possible_sources = possible_sources
                self.job.update(state="done", message="Map ready. Select a cell to inspect its evidence.",
                                output_dir=str(directory), molecules=len(result.data),
                                report=report, metrics=json_rows(result.metrics)[0],
                                options=options, files=list(files), elapsed=time.time() - self.job["started"])
                if self.query_frame is not None:
                    self.insert_new(self.query_options)
                else:
                    self.save_session(initial=True)
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
            summary = next((g for g in self.groups() if g['structure_group'] == group_id), None)
            if summary is None:
                raise ValueError("Unknown structure group.")
            payload = {'summary': summary, 'core_svg': _group_drawing(summary, 'core_smiles'),
                       'generic': self._generic_scaffold(self.result, summary)}
            self.group_detail_cache[group_id] = payload
            return payload

    def compare_group(self, group_id):
        with self.lock:
            if self.result is None:
                raise ValueError("Build a map first.")
            self._ensure_core_state()
            cache_key = (self.current_core_rule, group_id)
            if cache_key in self.group_cache:
                return self.group_cache[cache_key]
            summary = next((g for g in self.core_rule_cache[self.current_core_rule]
                            if g['structure_group'] == group_id), None)
            if summary is None:
                raise ValueError("Unknown structure group.")
            selected, svg, atom_audit = group_comparison(self.result.data, summary)
            members = self.result.data[self.result.data.structure_group.eq(group_id)]
            drawings = {field: _group_drawing(summary, field) for field in ('core_smiles', 'representative_smiles')}
            payload = {'summary': summary, 'svg': svg, 'drawings': drawings,
                       'members': json_rows(members), 'selected_count': len(selected),
                       'selection_rule': 'Representative first; then scaffold, canonical SMILES and identifier. At most 12 cards; all members are exported.',
                       'atom_matches': atom_audit}
            protocol = self.result.grouping_protocol
            is_build_rule = self.current_core_rule == self.job["options"]["core_rule"]
            if not is_build_rule:
                protocol = {**protocol, "core_rule": self.current_core_rule,
                            "core_min_coverage": summary["core_coverage_target"]}
            files = {'comparison.svg': svg.encode(), 'members.csv': members.to_csv(
                         index=False, lineterminator="\n").encode(),
                     'selected_members.csv': selected.to_csv(
                         index=False, lineterminator="\n").encode(),
                     'atom_matches.json': json.dumps(atom_audit, indent=2).encode(),
                     'group_summary.json': json.dumps({'summary': summary, 'protocol': protocol}, indent=2).encode()}
            digests = {name: hashlib.sha256(content).hexdigest() for name, content in files.items()}
            files['file_digests.json'] = json.dumps(digests, indent=2).encode()
            name = 'comparison_' + hashlib.sha256(group_id.encode()).hexdigest()[:12] + '.zip'
            if not is_build_rule:
                name = 'core_rules/' + name.removesuffix('.zip') + '_' + self.current_core_rule + '.zip'
            path = Path(self.job['output_dir']) / name
            path.parent.mkdir(exist_ok=True)
            with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as archive:
                for filename, content in files.items():
                    archive.writestr(filename, content)
            self.files[name] = path
            self.reference_files[name] = path
            payload['archive'] = name
            self.group_cache[cache_key] = payload
            self.save_session()
            return payload

    def reference_map_info(self):
        """参照说明仅用于界面，不改建图参数、参照报告或任何导出文件。"""
        options = self.job.get("options", {})
        target_id = None
        if options.get("input_format") == "chembl":
            parameters = self.curation.report["parameters"] if self.curation is not None else options
            target_id = parameters.get("target_id")
            target_id = str(target_id).strip().upper() if target_id else None
            parts = [target_id] if target_id else []
            columns = detect_chembl_columns(self.frame)
            if target_id and columns["target_id"] and columns["target_name"]:
                selected = self.frame[columns["target_id"]].fillna("").astype(str).str.strip().str.upper().eq(target_id)
                names = sorted({str(value).strip() for value in self.frame.loc[selected, columns["target_name"]].dropna()
                                if str(value).strip()})
                if names:
                    parts.append("; ".join(names))
            if options.get("activity_type"):
                parts.append(options["activity_type"])
        else:
            parts = [self.filename]
            activity_column = options.get("value_col") or options.get("label_col")
            if activity_column:
                parts.append(activity_column)
        return {"description": " · ".join(parts), "count": len(self.result.data), "target_id": target_id}

    def snapshot(self):
        with self.lock:
            snapshot = dict(self.job)
            if self.result is not None:
                snapshot["reference_map"] = self.reference_map_info()
            return snapshot

    def session_info(self):
        payload = {"version": __version__, "umap": importlib.util.find_spec("umap") is not None,
                   "output_root": str(self.output_root), "demo_available": DEMO.exists(),
                   "input": self.input_info(), "job": self.snapshot()}
        if self.query_input_info() is not None:
            payload["query_input"] = self.query_input_info()
        if self.new_result is not None:
            payload["new_molecules"] = self.new_molecules_info()
        if self.source_frame is not None:
            payload["source_input"] = self.source_input_info()
        return payload

    def _uploaded_new_source_records(self):
        """Use retained measurements, or uploaded batches for unmeasured molecules."""
        options = self.query_options or {}
        column = options.get("document_column")
        if self.new_result is None or not column:
            return pd.DataFrame(), column
        if options.get("value_col"):
            records = self.new_result.records.copy()
            if not records.empty:
                records["molecule_identity_key"] = "query:" + records.query_id.astype(str)
            return records, column
        checks = self.new_result.check
        records = self.query_frame.copy().reset_index(drop=True)
        identities = {}
        for row in checks.itertuples(index=False):
            if row.status == "mapped":
                identities[row.source_row] = "query:" + str(row.query_id)
            elif row.status == "duplicate_query" and not row.reference_id:
                identities[row.source_row] = "query:" + str(row.first_query_id)
        records["molecule_identity_key"] = records.index.map(identities)
        return records.dropna(subset=["molecule_identity_key"]), column

    def _new_source_records(self):
        records, column = self._uploaded_new_source_records()
        if not records.empty and self.measured_results:
            # Current manual evidence is added once below, using the upload's source column.
            identities = {"query:" + identity for identity in self.measured_results}
            records = records[~records.molecule_identity_key.isin(identities)].copy()
        manual = []
        if self.new_result is not None:
            for identity, result in self.measured_results.items():
                if result.get("batch"):
                    row = self.new_result.data[self.new_result.data.query_id.eq(identity)].iloc[0]
                    manual.append({"molecule_identity_key": "query:" + identity,
                                   "canonical_smiles": row.canonical_smiles,
                                   (column or "_manual_evidence_unit"): result["batch"]})
        if manual:
            return pd.concat([records, pd.DataFrame(manual)], ignore_index=True), column or "_manual_evidence_unit"
        return records, column

    def _source_index(self, scope="reference"):
        if self.result is None:
            raise ValueError("Build a map first.")
        if scope not in {"reference", "new"}:
            raise ValueError("Unknown source scope.")
        owner = self.result if scope == "reference" else self.new_result
        cached = self.source_index_cache.get(scope)
        if cached is not None and cached[0] is owner:
            return cached[1]
        if scope == "new":
            records, document = self._new_source_records()
            molecules = owner.data if owner is not None else pd.DataFrame()
            reference = None
        else:
            options = self.job["options"]
            records = self.curation.retained_records if self.curation is not None else self.generic_records
            document = (self.curation.report["detected_columns"].get("document_id")
                        if self.curation is not None else options["document_column"])
            reference = options.get("reference_column")
            if records is not None:
                records = self._reference_records(records, options)
            molecules = self.result.data
        index = build_source_index(molecules, records, document, reference)
        self.source_index_cache[scope] = (owner, index)
        return index

    def source_members(self, sources, scope="reference"):
        """Return a source union without rebuilding or writing any map files."""
        with self.lock:
            index = self._source_index(scope)
            identifiers = list(dict.fromkeys(str(source).strip() for source in sources))
            if not identifiers:
                raise ValueError("Choose at least one source.")
            for identifier in identifiers:
                if identifier not in index:
                    raise ValueError(f"Unknown source: {identifier}.")
            members = {member["molecule_identity_key"]: member for identifier in identifiers
                       for member in index[identifier]["members"]}
            return {"sources": [{key: value for key, value in index[identifier].items() if key != "members"}
                                for identifier in identifiers],
                    "map_molecule_count": len(members), "scope": scope,
                    "members": [members[key] for key in sorted(members)]}

    def _reference_records(self, records, options):
        """Recover selected reference strings without altering calculation records."""
        column = options.get("reference_column")
        if not column:
            return records
        references = _read_input(self.raw, dtype=str, keep_default_na=False)[column]
        records = records.copy()
        records[column] = records.source_row.map(references)
        return records

    def _query_reference_matches(self, result, curation, generic_records, options, checks=None):
        """Count retained documents for the browser only; never create audit exports."""
        checks = result.query_check if checks is None else checks
        checks = checks[checks.status.eq("in_reference")]
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
            checks = self.new_result.check if self.new_result is not None else self.result.query_check
            molecule_id, query_match = _query_selection(checks, molecule_id)
            if (self.new_result is not None and not self.new_result.data.empty
                    and not (query_match and query_match["reference_id"])):
                combined = pd.concat([data, self.new_result.data], ignore_index=True)
                selected = _select_molecule(combined, row, col, molecule_id)
                if selected.index[0] >= len(data):
                    # Recover the unmerged schema, so reference-only fields do not
                    # become artificial null values in the new-molecule inspector.
                    selected = self.new_result.data.iloc[[int(selected.index[0]) - len(data)]]
                    query = selected.iloc[0]
                    neighbors = self.new_result.neighbors[
                        self.new_result.neighbors.query_id.eq(query.query_id)]
                    records = self.new_result.records
                    if not records.empty:
                        record_key = next((column for column in (
                            "query_id", (self.query_options or {}).get("id_col"), "molecule_identity_key"
                        ) if column and column in records), None)
                        records = (records[records[record_key].astype(str).isin(
                            [str(query.query_id), "query:" + str(query.query_id)])]
                            if record_key else pd.DataFrame())
                    svg = ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 320 320">'
                           + (_svg_molecule_content(query.canonical_smiles, 320) or '') + '</svg>')
                    detail = {"molecule": json_rows(selected)[0], "query_neighbors": json_rows(neighbors),
                              "group_similarity": self.job["options"]["group_similarity"], "svg": svg,
                              "records": json_rows(records), "record_columns": records.columns.tolist(),
                              "total_records": len(records), "summary": {
                                  "reference_support_assessed": self.new_result.report["reference_support_assessed"]}}
                    identity = str(query.molecule_identity_key)
                    detail["new_molecule_sources"] = [
                        {key: value for key, value in source.items() if key != "members"}
                        for source in self._source_index("new").values()
                        if any(str(member["molecule_identity_key"]) == identity for member in source["members"])]
                    group = next((group for group in self.groups()
                                  if group["structure_group"] == query.structure_group), None)
                    if group is not None:
                        pattern = (Chem.MolFromSmarts(group["core_smarts"]) if group.get("core_smarts")
                                   else Chem.MolFromSmiles(group["core_smiles"]) if group["core_smiles"] else None)
                        detail["molecule"]["contains_group_core"] = bool(
                            pattern is not None and Chem.MolFromSmiles(query.canonical_smiles).HasSubstructMatch(
                                pattern, useChirality=True))
                    detail["molecule"]["contains_group_core_rule"] = self.current_core_rule
                    detail["group_core"] = group
                    detail["measured_result"] = self.measured_result_info(str(query.query_id))
                    if query_match is not None:
                        detail["query_match"] = query_match
                    return detail, records
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
                summary['excluded_bounds'] = excluded_bound_records(
                    curation.excluded_records, summary['molecule_identity_key'])
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
            index = self._source_index()
            sources = [{**source, "map_molecule_count": index.get(source["id"], {}).get("map_molecule_count", 0)}
                       for source in sources]
            source_union_count = len({member["molecule_identity_key"] for source in sources
                                      for member in index.get(source["id"], {}).get("members", [])})
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
            self.reference_files[archive.name] = archive
            detail = {"molecule": json_rows(chosen)[0], "summary": summary,
                    "sources": sources,
                    "source_union_molecule_count": source_union_count,
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
                measurement = next((match for match in self.job["report"].get("query_reference_matches", [])
                                    if match["query_id"] == query_match["query_id"]), None)
                if measurement and "own_median" in measurement:
                    detail["new_molecule_measurement"] = {key: measurement[key] for key in (
                        "own_median", "reference_median", "median_difference", "own_n_records")}
            if self.new_result is not None:
                detail["molecule"]["measurement_status"] = "measured"
            detail["group_core"] = next((group for group in self.groups()
                                         if group["structure_group"] == selected.get("structure_group")), None)
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
               "core_rule": {"most_frequent", "coverage", "coverage_brics", "coverage_recap", "fmcs"},
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
                self.send(session.session_info())
            elif path == "/api/core-progress":
                self.send(dict(session.core_progress))
            elif path == "/api/previous-runs":
                self.send(session.previous_runs())
            elif path == "/api/job":
                self.send(session.snapshot())
            elif path == "/api/points":
                self.send(session.points())
            elif path == "/api/new-molecules":
                self.send(session.new_molecules_info() or {"rows": [], "report": {}, "options": {},
                                                          "input": session.query_input_info(), "files": []})
            elif path == "/api/new-molecules/feedback":
                self.send(session.feedback_csv(), "text/csv", filename="neighbor_evidence_feedback.csv")
            elif path == "/api/groups":
                if session.result is None:
                    raise ValueError("Build a map first.")
                self.send({'groups': session.groups(), 'protocol': session.result.grouping_protocol,
                           'rule': session.current_core_rule, 'build_rule': session.job['options']['core_rule']})
            elif path == "/api/source-members":
                query = parse_qs(parsed.query)
                self.send(session.source_members(query.get("source", []), query.get("scope", ["reference"])[0]))
            elif path == "/api/group":
                self.send(session.group_detail(parse_qs(parsed.query).get('group', [''])[0]))
            elif path == "/api/group-cores":
                self.send(session.group_cores(parse_qs(parsed.query).get('group', [''])[0]))
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
        except SessionSnapshotError as exc:
            self.send(exc.payload(), status=400)
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
            elif path in ("/api/query-input", "/api/new-molecules/input"):
                if path == "/api/new-molecules/input":
                    if session.result is None:
                        raise ValueError("Build a reference map first.")
                    if session.job["options"]["representation"] != "morgan":
                        raise ValueError("New molecules require a reference map arranged using Morgan fingerprints.")
                self.send(session.load_query(content, unquote(self.headers.get("X-Filename", "queries.csv"))))
            elif path in ("/api/query-input/clear", "/api/new-molecules/clear"):
                self.send(session.clear_query())
            elif path == "/api/new-molecules/preview":
                self.send(session.preview_new(json.loads(content)))
            elif path == "/api/new-molecules/add":
                self.send(session.add_new(json.loads(content)))
            elif path == "/api/new-molecules/result":
                self.send(session.record_result(json.loads(content)))
            elif path == "/api/new-molecules/result/clear":
                self.send(session.clear_result(json.loads(content)))
            elif path in ("/api/new-molecules", "/api/new-molecules/insert"):
                self.send(session.insert_new(json.loads(content)))
            elif path == "/api/source-input":
                self.send(session.load_source(content, unquote(self.headers.get("X-Filename", "sources.csv"))))
            elif path == "/api/source-input/clear":
                self.send(session.clear_source())
            elif path == "/api/demo":
                self.send(session.load(DEMO.read_bytes(), DEMO.name))
            elif path == "/api/build":
                session.start(json.loads(content))
                self.send(session.snapshot(), status=202)
            elif path == "/api/core-rule":
                self.send(session.set_core_rule(json.loads(content).get("rule")))
            elif path == "/api/open-run":
                restored = session.open_run(json.loads(content).get("run_name"))["restored"]
                self.send({**session.session_info(), "restored": restored})
            elif path == "/api/shutdown":
                if session.snapshot()["state"] == "running":
                    raise ValueError("Wait for the current map to finish before closing the app.")
                with session.core_compute_lock:
                    shutdown_core_workers()
                self.send({"message": "Local app stopped. Your saved results remain on disk."})
                threading.Thread(target=self.server.shutdown, daemon=True).start()
            else:
                self.send({"error": "Not found."}, status=404)
        except SessionSnapshotError as exc:
            self.send(exc.payload(), status=400)
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
    def stop_app(signum, frame):
        raise KeyboardInterrupt

    previous_term = signal.signal(signal.SIGTERM, stop_app)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        shutdown_core_workers()
        server.server_close()
        signal.signal(signal.SIGTERM, previous_term)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
