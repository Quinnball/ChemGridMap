from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from chemgridmap import build_grid_map


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, default=Path("paper/output"))
    parser.add_argument("--smiles-col", default="smiles")
    parser.add_argument("--value-col", default="activity_pchembl")
    parser.add_argument("--label-col", default="activity_class")
    return parser.parse_args()


def main():
    args = parse_args()
    data = pd.read_csv(args.input)
    label_col = args.label_col if args.label_col in data.columns else None
    common = {
        "data": data,
        "smiles_col": args.smiles_col,
        "value_col": args.value_col,
        "label_col": label_col,
        "duplicate_policy": "error",
        "random_state": 42,
        "k": 10,
    }

    results = []
    projection_dir = args.output / "projection_comparison"
    for method in ["pca", "tsne", "umap"]:
        results.append(
            build_grid_map(
                output_dir=projection_dir,
                name="morgan_{}".format(method),
                representation="morgan",
                projection=method,
                **common
            )
        )

    representation_dir = args.output / "representation_comparison"
    for representation in ["morgan", "descriptors"]:
        results.append(
            build_grid_map(
                output_dir=representation_dir,
                name="{}_umap".format(representation),
                representation=representation,
                projection="umap",
                **common
            )
        )

    if any(str(column).startswith("emb_") for column in data.columns):
        results.append(
            build_grid_map(
                output_dir=representation_dir,
                name="embedding_umap",
                representation="embedding",
                projection="umap",
                **common
            )
        )

    args.output.mkdir(parents=True, exist_ok=True)
    metrics = pd.concat([result.metrics for result in results], ignore_index=True)
    metrics.to_csv(args.output / "comparison_metrics.csv", index=False)
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()

