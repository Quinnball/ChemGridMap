"""Local browser interface. Scientific operations use the same API as the CLI."""
from __future__ import annotations

import argparse
import csv
import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import importlib.util
import io
import json
import mimetypes
from pathlib import Path
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
from .inspection import inspect_entry, load_audit_context
from .pipeline import build_grid_map
from .plotting import render_svg
from .provenance import file_digests
from .neighborhoods import neighborhood_evidence
from .comparison import group_comparison
from .plotting import _svg_molecule_content

ASSETS = Path(__file__).with_name("web")
DEMO = Path(__file__).with_name("data") / "chembl205_ic50.csv"
MAX_UPLOAD = 50 * 1024 * 1024


def json_rows(frame):
    return json.loads(frame.to_json(orient="records"))


class AppSession:
    """One local workspace; a single worker avoids concurrent chemistry/render jobs."""

    def __init__(self, output_root):
        self.output_root = Path(output_root).expanduser().resolve()
        self.token = secrets.token_urlsafe(32)
        self.lock = threading.RLock()
        self.frame = None
        self.raw = b""
        self.filename = ""
        self.result = None
        self.curation = None
        self.files = {}
        self.group_cache = {}
        self.job = {"state": "empty", "message": "Choose a CSV or open the example."}

    def load(self, content, filename):
        with self.lock:
            if self.job["state"] == "running":
                raise ValueError("A map is running. Wait for it to finish before replacing the input.")
            if not content or len(content) > MAX_UPLOAD:
                raise ValueError("Choose a non-empty CSV smaller than 50 MB.")
            try:
                sample = content.decode("utf-8-sig", errors="strict")[:65536]
                try:
                    delimiter = csv.Sniffer().sniff(sample, delimiters=",;\t|").delimiter
                except csv.Error:
                    delimiter = ","
                frame = pd.read_csv(io.BytesIO(content), sep=delimiter, encoding="utf-8-sig")
            except (UnicodeError, pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
                raise ValueError("Could not read this UTF-8 CSV. Check its delimiter and encoding.") from exc
            if frame.empty:
                raise ValueError("The CSV has no data rows.")
            self.frame, self.raw = frame, content
            self.filename = Path(filename.replace("\\", "/")).name
            self.result, self.curation, self.files = None, None, {}
            self.group_cache = {}
            self.job = {"state": "ready", "message": "Input loaded. Review the options, then build the map."}
            return self.input_info()

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
            self.result, self.curation, self.files = None, None, {}
            self.group_cache = {}
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
            }
            data = self.frame.copy()
            curation, curation_files = None, {}
            if options["input_format"] == "chembl":
                self.update(message="Validating structures and auditing repeated measurements...")
                curation = curate_chembl_activity_data(
                    data, target_id=options["target_id"], activity_type=options["activity_type"],
                    lower_threshold=options["lower_threshold"], upper_threshold=options["upper_threshold"],
                    conflict_range_threshold=options["conflict_range_threshold"],
                    structure_policy=options["structure_policy"], molecule_identity=options["molecule_identity"],
                )
                curation.report.update(source)
                source["curation_parameters"] = curation.report["parameters"]
                document_column = curation.report["detected_columns"].get("document_id")
                if options["activity_confidence"] and document_column:
                    self.update(message="Estimating how well source documents support each activity class...")
                    try:
                        curation.molecules, fit = annotate_activity_confidence(
                            curation.molecules, curation.retained_records, lower=options["lower_threshold"],
                            upper=options["upper_threshold"], value="activity_pchembl_raw", document=document_column,
                            merge_identical=options["merge_identical_documents"],
                            model=options["evidence_model"])
                        curation.report["activity_confidence"] = fit
                        source["activity_confidence"] = fit
                    except ValueError as exc:
                        curation.report["warnings"].append(f"Activity class support was not estimated: {exc}")
                elif options["activity_confidence"]:
                    curation.report["warnings"].append(
                        "Activity class support was not estimated: the input has no document identifiers.")
                curation_files = save_chembl_curation(curation, directory, "map")
                source["curation_outputs"] = file_digests(curation_files)
                data = curation.molecules
                smiles_col, value_col, label_col = "canonical_smiles", "activity_pchembl", "activity_class"
            else:
                smiles_col, value_col, label_col = options["smiles_col"], options["value_col"], options["label_col"]
                if options["compound_column"]:
                    if not value_col:
                        raise ValueError("A pActivity column is required for generic activity records.")
                    confidence = bool(options["document_column"] and options["activity_confidence"])
                    data, _, fit, report = prepare_generic_activity_records(
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
                        merge_identical=options["merge_identical_documents"],
                        evidence_model=options["evidence_model"],
                    )
                    source["generic_activity"] = report
                    if fit is not None:
                        source["activity_confidence"] = fit
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
                map_color_modes=color_modes, map_color_formats=["svg"],
            )
            result.output_files.update(curation_files)
            files = {Path(path).name: Path(path) for path in result.output_files.values()}
            files.update({name: directory / name for name in ("input.csv", "app_options.json")})
            archive = directory / "chemgridmap_results.zip"
            with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as bundle:
                for name, path in files.items():
                    bundle.write(path, name)
            files[archive.name] = archive
            report = curation.report if curation else result.preprocessing
            with self.lock:
                self.result, self.curation, self.files = result, curation, files
                self.job.update(state="done", message="Map ready. Select a cell to inspect its evidence.",
                                output_dir=str(directory), molecules=len(result.data),
                                report=report, metrics=json_rows(result.metrics)[0],
                                options=options, files=list(files), elapsed=time.time() - self.job["started"])
        except Exception as exc:
            self.update(state="error", message=str(exc), output_dir=str(directory))

    def update(self, **values):
        with self.lock:
            self.job.update(values)

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
            drawings = {}
            for field in ('core_smiles', 'representative_smiles'):
                drawings[field] = ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 180 180">'
                                   + (_svg_molecule_content(summary[field], 180) or '') + '</svg>') if summary[field] else ''
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

    def inspect(self, row=None, col=None, molecule_id=None):
        with self.lock:
            if self.result is None:
                raise ValueError("Build a map first.")
            data, curation = self.result.data, self.curation
            if curation is not None:
                _, run_bound, binding_warnings = load_audit_context(
                    Path(self.result.output_files['coordinates']), Path(self.result.output_files['retained_records']))
                chosen, records, summary = inspect_entry(
                    data, curation.retained_records, molecule_id=molecule_id,
                    cell=(row, col) if molecule_id is None else None,
                    parameters=curation.report["parameters"])
                summary['run_binding_verified'] = run_bound
                summary['warnings'].extend(binding_warnings)
            else:
                if molecule_id is not None:
                    raise ValueError("For molecule-table input, select a cell on the map.")
                chosen = data[data.grid_row.eq(row) & data.grid_col.eq(col)]
                if len(chosen) != 1:
                    raise ValueError("Select one occupied cell.")
                records = pd.DataFrame()
                summary = {"interpretation": "Molecule-table input: no assay-record audit is available.",
                           "record_count_verified": False, "median_verified": False,
                           "quality_annotations_verified": False, "run_binding_verified": False}
            selected = chosen.iloc[0]
            if curation is not None and selected.get("discordant_document_pairs", 0):
                document_column = curation.report["detected_columns"].get("document_id")
                documents = set(records[document_column].dropna().astype(str))
                offsets = curation.report.get("activity_confidence", {}).get("document_offsets", [])
                summary["document_offsets"] = [pair for pair in offsets
                    if pair["document_a"] in documents and pair["document_b"] in documents]
                for pair in summary["document_offsets"]:
                    summary["warnings"].append(
                        "Source documents disagree systematically: {} vs {} "
                        "(median offset {:.3g} log units, {} shared entries).".format(
                            pair["document_a"], pair["document_b"], pair["median_offset"], pair["shared_entries"]))
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
            return {"molecule": json_rows(chosen)[0], "summary": summary,
                    "records": json_rows(records.head(100)), "total_records": len(records),
                    "neighbors": json_rows(neighbor_table), "audit_file": archive.name,
                    "svg": image.read_text(encoding="utf-8"),
                    "neighborhood_svg": neighborhood.read_text(encoding="utf-8")}, records


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
                "evidence_model": "gaussian", "merge_identical_documents": True}
    if not isinstance(options, dict) or set(options) - set(defaults):
        raise ValueError("Unknown map option.")
    values = {**defaults, **options}
    choices = {"input_format": {"chembl", "molecule"}, "representation": {"morgan", "descriptors"},
               "projection": {"pca", "tsne", "umap"}, "render_detail": {"auto", "full", "overview"},
               "structure_policy": {"auto", "parent", "as-recorded"},
               "molecule_identity": {"auto", "parent-id", "canonical-smiles"},
               "grid_refinement": {"none", "neighbors"},
               "core_rule": {"most_frequent", "coverage", "fmcs"},
               "evidence_model": {"gaussian", "exclude-flagged"}}
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
        raise ValueError("Conflict range must be finite and non-negative.")
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
            if path in ("/", "/app.js", "/app.css", "/example-map.svg", "/example-molecule.svg"):
                asset = ASSETS / ("index.html" if path == "/" else path[1:])
                self.send(asset.read_bytes(), mimetypes.guess_type(asset.name)[0] or "text/plain")
            elif path == "/api/session":
                self.send({"version": __version__, "umap": importlib.util.find_spec("umap") is not None,
                           "output_root": str(session.output_root), "demo_available": DEMO.exists(),
                           "input": session.input_info(), "job": session.snapshot()})
            elif path == "/api/job":
                self.send(session.snapshot())
            elif path == "/api/points":
                if session.result is None:
                    raise ValueError("Build a map first.")
                columns = [name for name in ("projection_x", "projection_y", "grid_x", "grid_y", "grid_row", "grid_col",
                           "activity_class", "molecule_id", "molecule_identity_key", "structure_group", "structure_group_color", "comparison_id",
                           "neighbor_recall", "class_confidence", "class_status", "n_evidence_units") if name in session.result.data]
                self.send(json_rows(session.result.data[columns]))
            elif path == "/api/groups":
                if session.result is None:
                    raise ValueError("Build a map first.")
                self.send({'groups': session.result.structure_groups, 'protocol': session.result.grouping_protocol})
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
