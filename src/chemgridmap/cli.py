"""Command-line entry point."""

from __future__ import annotations

import argparse
import hashlib
import sys
import zipfile
from pathlib import Path
from typing import Optional, Sequence

import pandas as pd

from .chembl import curate_chembl_activity_data, save_chembl_curation
from .evidence import annotate_activity_confidence, prepare_generic_activity_records
from .pipeline import build_grid_map


def _core_min_coverage(value: str) -> float:
    coverage = float(value)
    if not 0.5 <= coverage <= 1:
        raise argparse.ArgumentTypeError("Core coverage target must be between 0.5 and 1 (inclusive).")
    return coverage


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="chemgridmap",
        description=(
            "Convert molecular representations or existing 2D coordinates into "
            "a one-molecule-per-cell grid chemical map."
        ),
        epilog="To retrieve measurements behind an existing cell, run: chemgridmap inspect --help",
    )
    parser.add_argument("input", type=Path, help="Input CSV file.")
    parser.add_argument(
        "--input-format",
        choices=["molecule", "chembl"],
        default="molecule",
        help=(
            "Use 'chembl' for a raw ChEMBL activity CSV that should be "
            "filtered and aggregated before mapping."
        ),
    )
    parser.add_argument(
        "-o", "--output-dir", type=Path, default=Path("gridmap_output")
    )
    parser.add_argument("--name", default="grid_map", help="Output file prefix.")
    parser.add_argument("--smiles-col", "--smiles-column", dest="smiles_col", default="smiles")
    parser.add_argument("--value-col", "--value-column", dest="value_col")
    parser.add_argument("--label-col")
    parser.add_argument(
        "--compound-column",
        help="Compound identifier column for a generic activity-record CSV.",
    )
    parser.add_argument(
        "--document-column",
        help="Source-document column for a generic activity-record CSV.",
    )
    parser.add_argument(
        "--relation-column",
        help="Optional relation column for a generic activity-record CSV; only '=' is retained.",
    )
    parser.add_argument(
        "--representation",
        choices=["morgan", "descriptors", "embedding"],
        default="morgan",
    )
    parser.add_argument(
        "--projection", choices=["pca", "tsne", "umap"], default="umap"
    )
    parser.add_argument(
        "--umap-metric",
        choices=["auto", "euclidean", "jaccard", "cosine"],
        default="auto",
        help=(
            "Distance used by UMAP. Auto selects Jaccard for Morgan "
            "fingerprints and Euclidean distance otherwise."
        ),
    )
    parser.add_argument("--x-col", help="Existing 2D x-coordinate column.")
    parser.add_argument("--y-col", help="Existing 2D y-coordinate column.")
    parser.add_argument("--embedding-prefix", default="emb_")
    parser.add_argument("--radius", type=int, default=2)
    parser.add_argument("--n-bits", type=int, default=2048)
    parser.add_argument("--grouping", choices=["butina", "external", "none"], default="butina")
    parser.add_argument("--group-similarity", type=float, default=.6,
                        help="Morgan Tanimoto threshold to a Butina center; not an all-pairs bound.")
    parser.add_argument("--cluster-col", help="Existing group column in molecule-table input.")
    parser.add_argument("--cluster-source", help="Provenance of external group assignments.")
    parser.add_argument("--lower-threshold", type=float, default=6.0)
    parser.add_argument("--upper-threshold", type=float, default=7.0)
    parser.add_argument(
        "--target-id",
        help="Expected Target ChEMBL ID, for example CHEMBL205.",
    )
    parser.add_argument(
        "--activity-type",
        default="IC50",
        help="ChEMBL Standard Type retained in chembl input mode.",
    )
    parser.add_argument(
        "--conflict-range-threshold",
        type=float,
        default=1.0,
        help="pChEMBL range above which repeated measurements are marked as a wide record range.",
    )
    parser.add_argument(
        "--keep-potential-duplicates",
        action="store_true",
        help="Keep rows flagged by ChEMBL as potential duplicate citations.",
    )
    parser.add_argument(
        "--molecule-identity",
        choices=["auto", "parent-id", "canonical-smiles"],
        default="auto",
        help=(
            "Group ChEMBL records by parent molecule ID when available, or "
            "fall back to canonical SMILES."
        ),
    )
    parser.add_argument(
        "--duplicate-policy",
        choices=["error", "first", "keep"],
        default="error",
        help="How repeated canonical structures are handled.",
    )
    parser.add_argument("--structure-policy", choices=["auto", "parent", "as-recorded"],
                        default="auto", help="Normalize parent structures, or explicitly retain record forms.")
    parser.add_argument("--assay-types", nargs="+", help="Optional ChEMBL assay-type filter; not a reliability score.")
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("-k", type=int, default=10)
    parser.add_argument(
        "--trustworthiness-sample-size",
        type=int,
        default=3000,
        help="Maximum deterministic sample used for trustworthiness.",
    )
    parser.add_argument(
        "--grid-padding",
        type=int,
        help=(
            "Explicit cells added to each grid side. When omitted, grid size "
            "is determined by --grid-occupancy."
        ),
    )
    parser.add_argument(
        "--grid-occupancy",
        type=float,
        default=0.40,
        help="Target fraction of occupied candidate cells (default: 0.40).",
    )
    parser.add_argument(
        "--coordinate-scaling",
        choices=["isotropic", "independent"],
        default="isotropic",
        help=(
            "Scale both axes uniformly to preserve projection geometry, or "
            "use legacy independent axis scaling."
        ),
    )
    parser.add_argument(
        "--assignment-method",
        choices=["auto", "dense", "sparse"],
        default="auto",
        help="Use exact dense matching, scalable sparse matching, or auto-select.",
    )
    parser.add_argument(
        "--sparse-neighbors",
        type=int,
        default=128,
        help=(
            "Initial nearest candidate cells per molecule in adaptive sparse assignment "
            "(default: 128)."
        ),
    )
    parser.add_argument("--sparse-max-neighbors", type=int, default=1024,
                        help="Candidate cap for adaptive sparse matching (default: 1024).")
    parser.add_argument("--fixed-sparse-candidates", action="store_true",
                        help="Disable candidate expansion, for explicit calibration only.")
    parser.add_argument(
        "--render-detail",
        choices=["auto", "full", "overview"],
        default="auto",
        help=(
            "Draw molecule structures in every cell ('full') or color-only "
            "cells ('overview'). Auto uses overview above 2,000 molecules."
        ),
    )
    parser.add_argument(
        "--output-formats",
        nargs="+",
        choices=["svg", "png", "pdf"],
        default=["svg", "png", "pdf"],
        help="Map formats to write.",
    )
    parser.add_argument(
        "--grid-refinement", choices=["neighbors", "none"], default="neighbors",
        help="Refine the minimum-displacement grid so representation neighbours sit closer "
             "('none' reproduces 0.4.0 layouts).")
    parser.add_argument(
        "--core-rule", choices=["coverage", "most_frequent", "fmcs"], default="coverage",
        help="Group core: largest scaffold or ring fragment covering --core-min-coverage of members, "
             "the most frequent scaffold (0.4.0 rule), or FMCS.")
    parser.add_argument(
        "--fmcs-timeout", type=int, default=1,
        help="FMCS search budget in seconds (default: 1). Timeout may make results machine-dependent.")
    parser.add_argument(
        "--core-min-coverage", type=_core_min_coverage, default=0.8,
        help="Coverage target for coverage/FMCS rules; range: 0.5–1 (default: 0.8).")
    parser.add_argument(
        "--no-activity-confidence", action="store_true",
        help="Skip document-level class-support estimation.")
    parser.add_argument(
        "--evidence-model", choices=["gaussian", "exclude-flagged"], default="gaussian",
        help="Activity evidence model (default: gaussian).")
    parser.add_argument(
        "--no-merge-identical-documents", action="store_true",
        help="Keep documents with identical entry medians as separate evidence units (default: merge).")
    parser.add_argument(
        "--map-color-modes", nargs="+", default=[],
        choices=["activity", "class_support", "neighbor_recall", "structure_group"],
        help="Additional map layers to export as SVG and, when requested in --output-formats, PNG.")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "compare":
        from .comparison import main as compare_main
        return compare_main(argv[1:])
    if argv and argv[0] == "inspect":
        from .inspection import main as inspect_main
        return inspect_main(argv[1:])
    args = build_parser().parse_args(argv)
    data = pd.read_csv(args.input, sep=None, engine="python", encoding="utf-8-sig")
    source = {"input_file": args.input.name,
              "input_file_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
            "reliability_parameters": {"activity_confidence": not args.no_activity_confidence,
                                         "evidence_model": args.evidence_model,
                                         "merge_identical_documents": not args.no_merge_identical_documents}}
    curation_files = {}
    if args.input_format == "chembl":
        curation = curate_chembl_activity_data(
            data,
            target_id=args.target_id,
            activity_type=args.activity_type,
            lower_threshold=args.lower_threshold,
            upper_threshold=args.upper_threshold,
            conflict_range_threshold=args.conflict_range_threshold,
            exclude_potential_duplicates=not args.keep_potential_duplicates,
            molecule_identity=args.molecule_identity,
            structure_policy=args.structure_policy,
            assay_types=args.assay_types,
        )
        curation.report["input_file"] = str(args.input.resolve())
        curation.report["input_file_sha256"] = source["input_file_sha256"]
        source["curation_parameters"] = curation.report["parameters"]
        document_column = curation.report["detected_columns"].get("document_id")
        if not args.no_activity_confidence and document_column:
            try:
                curation.molecules, fit = annotate_activity_confidence(
                    curation.molecules, curation.retained_records, lower=args.lower_threshold,
                    upper=args.upper_threshold, value="activity_pchembl_raw", document=document_column,
                    merge_identical=not args.no_merge_identical_documents,
                    model=args.evidence_model)
                curation.report["activity_confidence"] = fit
                source["activity_confidence"] = fit
            except ValueError as exc:
                curation.report["warnings"].append(f"Activity class support was not estimated: {exc}")
        elif not args.no_activity_confidence:
            curation.report["warnings"].append(
                "Activity class support was not estimated: the input has no document identifiers.")
        curation_files = save_chembl_curation(
            curation,
            output_dir=args.output_dir,
            name=args.name,
        )
        data = curation.molecules
        from .provenance import file_digests
        source["curation_outputs"] = file_digests(curation_files)
        smiles_col = "canonical_smiles"
        value_col = "activity_pchembl"
        label_col = "activity_class"
        duplicate_policy = "error"
    else:
        curation = None
        smiles_col = args.smiles_col
        value_col = args.value_col
        label_col = args.label_col
        duplicate_policy = args.duplicate_policy
        if args.compound_column:
            if not value_col:
                raise ValueError("--value-column is required with --compound-column.")
            confidence = bool(args.document_column and not args.no_activity_confidence)
            data, records, fit, generic_report = prepare_generic_activity_records(
                data,
                compound_column=args.compound_column,
                value_column=value_col,
                document_column=args.document_column,
                smiles_column=smiles_col,
                relation_column=args.relation_column,
                lower=args.lower_threshold,
                upper=args.upper_threshold,
                conflict_range_threshold=args.conflict_range_threshold,
                confidence=confidence,
                merge_identical=not args.no_merge_identical_documents,
                evidence_model=args.evidence_model,
            )
            source["generic_activity"] = generic_report
            if fit is not None:
                source["activity_confidence"] = fit
            smiles_col = "canonical_smiles"
            value_col = "activity_pchembl"
            label_col = "activity_class"
            duplicate_policy = "error"

    result = build_grid_map(
        data,
        grouping=args.grouping, group_similarity=args.group_similarity,
        cluster_col=args.cluster_col, cluster_source=args.cluster_source,
        output_dir=args.output_dir,
        name=args.name,
        smiles_col=smiles_col,
        value_col=value_col,
        label_col=label_col,
        representation=args.representation,
        projection=args.projection,
        umap_metric=args.umap_metric,
        x_col=args.x_col,
        y_col=args.y_col,
        embedding_prefix=args.embedding_prefix,
        radius=args.radius,
        n_bits=args.n_bits,
        lower_threshold=args.lower_threshold,
        upper_threshold=args.upper_threshold,
        duplicate_policy=duplicate_policy,
        random_state=args.random_state,
        k=args.k,
        trustworthiness_sample_size=args.trustworthiness_sample_size,
        grid_padding=args.grid_padding,
        grid_occupancy=args.grid_occupancy,
        coordinate_scaling=args.coordinate_scaling,
        assignment_method=args.assignment_method,
        sparse_neighbors=args.sparse_neighbors,
        sparse_adaptive=not args.fixed_sparse_candidates,
        sparse_max_neighbors=args.sparse_max_neighbors,
        render_detail=args.render_detail,
        output_formats=args.output_formats,
        provenance=source,
        grid_refinement=args.grid_refinement,
        core_rule=args.core_rule,
        fmcs_timeout=args.fmcs_timeout,
        core_min_coverage=args.core_min_coverage,
        map_color_modes=args.map_color_modes,
    )
    result.output_files.update(curation_files)
    archive = args.output_dir / (args.name + "_results.zip")
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as bundle:
        for path in result.output_files.values():
            bundle.write(path, Path(path).name)
    result.output_files["results_zip"] = str(archive)

    report = result.preprocessing
    print(
        "Prepared {output_rows} molecules from {input_rows} input rows "
        "({invalid_smiles} invalid SMILES).".format(**report)
    )
    if curation is not None:
        audit = curation.report
        print(
            "ChEMBL curation retained {retained_activity_records} records and "
            "produced {unique_molecules} unique molecules "
            "({repeated_conflicted_molecules} with a wide record range).".format(
                **audit
            )
        )
        for warning in audit["warnings"]:
            print("Warning: {}".format(warning))
    print("Saved files:")
    for key, path in result.output_files.items():
        print("  {:12s} {}".format(key, path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
