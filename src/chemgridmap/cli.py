"""Command-line entry point."""

from __future__ import annotations

import argparse
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
        "--duplicate-policy",
        choices=["error", "first", "keep"],
        default="error",
        help="How repeated canonical structures are handled.",
    )
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("-k", type=int, default=10)
    parser.add_argument(
        "--trustworthiness-sample-size",
        type=int,
        default=3000,
        help="Maximum deterministic sample used for trustworthiness.",
    )
    parser.add_argument("--grid-padding", type=int, default=20)
    parser.add_argument(
        "--assignment-method",
        choices=["auto", "dense", "sparse"],
        default="auto",
        help="Use exact dense matching, scalable sparse matching, or auto-select.",
    )
    parser.add_argument(
        "--sparse-neighbors",
        type=int,
        default=32,
        help="Nearest candidate cells per molecule in sparse assignment.",
    )
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
    args = build_parser().parse_args(argv)
    data = pd.read_csv(args.input)
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
        )
        curation.report["input_file"] = str(args.input.resolve())
        curation_files = save_chembl_curation(
            curation,
            output_dir=args.output_dir,
            name=args.name,
        )
        data = curation.molecules
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
        assignment_method=args.assignment_method,
        sparse_neighbors=args.sparse_neighbors,
        render_detail=args.render_detail,
        output_formats=args.output_formats,
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
