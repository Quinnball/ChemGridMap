from pathlib import Path

import pandas as pd

from chemgridmap import build_grid_map


HERE = Path(__file__).resolve().parent


def main():
    molecules = pd.read_csv(HERE / "chembl205_embedding_umap.csv")
    result = build_grid_map(
        molecules,
        output_dir=HERE / "output",
        name="chembl205_embedding_umap",
        value_col="activity_pchembl",
        label_col="activity_class",
        x_col="projection_x",
        y_col="projection_y",
        grid_occupancy=0.40,
        coordinate_scaling="isotropic",
    )
    print(result.metrics.to_string(index=False))


if __name__ == "__main__":
    main()
