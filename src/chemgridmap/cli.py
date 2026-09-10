"""Command-line entry point."""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path
from typing import Optional, Sequence

import pandas as pd

from .chembl import curate_chembl_activity_data, save_chembl_curation
from .pipeline import build_grid_map


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
    parser.add_argument("--smiles-col", default="smiles")
    parser.add_argument("--value-col")
    parser.add_argument("--label-col")
    parser.add_argument(
        "--representation",
        choices=["morgan", "descriptors", "embedding"],
        default="morgan",
    )
    parser.add_argument(
        "--projection", choices=["pca", "tsne", "umap"], default="pca"
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
        help="pChEMBL range above which repeated measurements are flagged.",
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
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "inspect":
        from .inspection import main as inspect_main
        return inspect_main(argv[1:])
    args = build_parser().parse_args(argv)
    data = pd.read_csv(args.input, sep=None, engine="python", encoding="utf-8-sig")
    source = {"input_file": args.input.name,
              "input_file_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest()}
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

    result = build_grid_map(
        data,
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
    )
    result.output_files.update(curation_files)

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
            "({repeated_conflicted_molecules} conflicted repeats).".format(
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
