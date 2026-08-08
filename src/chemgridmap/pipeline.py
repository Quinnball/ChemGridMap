"""High-level grid chemical map workflow."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from .core import (
    assign_to_grid,
    compute_metrics,
    normalize_coordinates,
    prepare_molecules,
)
from .plotting import render_grid_map
from .projection import project_2d
from .representations import compute_representation


@dataclass
class GridMapResult:
    """Data and files produced by one grid-map run."""

    data: pd.DataFrame
    representation: np.ndarray
    projection: np.ndarray
    grid: np.ndarray
    metrics: pd.DataFrame
    preprocessing: Dict[str, int]
    feature_names: List[str] = field(default_factory=list)
    output_files: Dict[str, str] = field(default_factory=dict)


def build_grid_map(
    data: pd.DataFrame,
    output_dir: Optional[Path] = None,
    name: str = "grid_map",
    smiles_col: str = "smiles",
    value_col: Optional[str] = None,
    label_col: Optional[str] = None,
    representation: str = "morgan",
    projection: str = "pca",
    x_col: Optional[str] = None,
    y_col: Optional[str] = None,
    embedding_prefix: str = "emb_",
    radius: int = 2,
    n_bits: int = 2048,
    lower_threshold: float = 6.0,
    upper_threshold: float = 7.0,
    duplicate_policy: str = "error",
    random_state: int = 42,
    k: int = 10,
    grid_padding: int = 20,
    color_map: Optional[Dict[str, str]] = None,
) -> GridMapResult:
    """Build and optionally save a one-molecule-per-cell grid chemical map."""
    if bool(x_col) != bool(y_col):
        raise ValueError("x_col and y_col must be supplied together.")

    prepared, report = prepare_molecules(
        data,
        smiles_col=smiles_col,
        value_col=value_col,
        label_col=label_col,
        lower_threshold=lower_threshold,
        upper_threshold=upper_threshold,
        duplicate_policy=duplicate_policy,
    )
    if len(prepared) < 2:
        raise ValueError("At least two valid molecule rows are required.")

    if x_col and y_col:
        for column in [x_col, y_col]:
            if column not in prepared.columns:
                raise ValueError("Coordinate column {!r} was not found.".format(column))
        projected = (
            prepared[[x_col, y_col]]
            .apply(pd.to_numeric, errors="raise")
            .to_numpy(np.float64)
        )
        representation_matrix = projected.copy()
        feature_names = [x_col, y_col]
        representation_name = "user_coordinates"
        projection_name = "coordinates"
    else:
        representation_matrix, feature_names = compute_representation(
            prepared,
            representation=representation,
            radius=radius,
            n_bits=n_bits,
            embedding_prefix=embedding_prefix,
        )
        projected = project_2d(
            representation_matrix,
            method=projection,
            random_state=random_state,
        )
        representation_name = representation.lower()
        projection_name = projection.lower()

    grid, grid_rows, grid_cols = assign_to_grid(projected, padding=grid_padding)
    projected_01 = normalize_coordinates(projected)
    result_data = prepared.copy()
    result_data["projection_x"] = projected_01[:, 0]
    result_data["projection_y"] = projected_01[:, 1]
    result_data["grid_x"] = grid[:, 0]
    result_data["grid_y"] = grid[:, 1]
    result_data["grid_col"] = np.rint(grid[:, 0] * (grid_cols - 1)).astype(int)
    result_data["grid_row"] = np.rint(grid[:, 1] * (grid_rows - 1)).astype(int)

    values = None
    if value_col:
        numeric = pd.to_numeric(result_data[value_col], errors="coerce").to_numpy()
        if np.isfinite(numeric).any():
            values = numeric

    label_values = result_data["activity_class"].astype(str).to_numpy()
    if np.any(label_values == "unlabelled"):
        label_values = None

    metrics = compute_metrics(
        representation_matrix,
        projected,
        grid,
        values=values,
        labels=label_values,
        k=k,
    )
    metrics.insert(0, "projection_method", projection_name)
    metrics.insert(0, "representation_type", representation_name)
    metrics.insert(0, "map_name", name)

    output_files = {}
    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        coordinates_path = output_dir / "{}_grid_coordinates.csv".format(name)
        metrics_path = output_dir / "{}_metrics.csv".format(name)
        result_data.to_csv(coordinates_path, index=False)
        metrics.to_csv(metrics_path, index=False)
        output_files.update(
            {
                "coordinates": str(coordinates_path),
                "metrics": str(metrics_path),
            }
        )
        output_files.update(
            render_grid_map(
                result_data,
                output_dir / name,
                color_map=color_map,
            )
        )

    return GridMapResult(
        data=result_data,
        representation=representation_matrix,
        projection=projected,
        grid=grid,
        metrics=metrics,
        preprocessing=report,
        feature_names=feature_names,
        output_files=output_files,
    )

