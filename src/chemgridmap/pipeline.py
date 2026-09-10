"""High-level grid chemical map workflow."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import hashlib
import json
import time
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from .core import (
    assign_to_grid,
    choose_assignment_method,
    compute_metrics,
    normalize_coordinates,
    prepare_molecules,
)
from .plotting import render_grid_map
from .projection import project_2d
from .representations import compute_representation
from .provenance import environment_versions


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
    umap_metric: str = "auto",
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
    trustworthiness_sample_size: int = 3000,
    grid_padding: Optional[int] = None,
    grid_occupancy: float = 0.40,
    coordinate_scaling: str = "isotropic",
    assignment_method: str = "auto",
    dense_max_pairs: int = 20_000_000,
    sparse_neighbors: int = 128,
    render_detail: str = "auto",
    output_formats: Optional[List[str]] = None,
    color_map: Optional[Dict[str, str]] = None,
    provenance: Optional[Dict] = None,
    sparse_adaptive: bool = True,
    sparse_max_neighbors: int = 1024,
) -> GridMapResult:
    """Build and optionally save a one-molecule-per-cell grid chemical map."""
    started = time.perf_counter()
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
        representation_distance_metric = "euclidean"
        selected_umap_metric = "not_applicable"
    else:
        representation_matrix, feature_names = compute_representation(
            prepared,
            representation=representation,
            radius=radius,
            n_bits=n_bits,
            embedding_prefix=embedding_prefix,
        )
        selected_umap_metric = str(umap_metric).strip().lower()
        if selected_umap_metric == "auto":
            selected_umap_metric = (
                "jaccard" if representation.lower() == "morgan" else "euclidean"
            )
        projected = project_2d(
            representation_matrix,
            method=projection,
            random_state=random_state,
            umap_metric=selected_umap_metric,
        )
        representation_name = representation.lower()
        projection_name = projection.lower()
        representation_distance_metric = (
            "jaccard" if representation_name == "morgan" else "euclidean"
        )

    selected_assignment_method = choose_assignment_method(
        len(projected),
        padding=grid_padding,
        target_occupancy=grid_occupancy,
        method=assignment_method,
        dense_max_pairs=dense_max_pairs,
    )
    projection_seconds = time.perf_counter() - started
    assignment_started = time.perf_counter()
    assignment_diagnostics = {}
    grid, grid_rows, grid_cols = assign_to_grid(
        projected,
        padding=grid_padding,
        target_occupancy=grid_occupancy,
        coordinate_scaling=coordinate_scaling,
        method=selected_assignment_method,
        dense_max_pairs=dense_max_pairs,
        sparse_neighbors=sparse_neighbors,
        diagnostics=assignment_diagnostics,
        sparse_adaptive=sparse_adaptive,
        sparse_max_neighbors=sparse_max_neighbors,
    )
    assignment_seconds = time.perf_counter() - assignment_started
    projected_01 = normalize_coordinates(projected, mode=coordinate_scaling)
    result_data = prepared.copy()
    result_data["projection_x_raw"] = projected[:, 0]
    result_data["projection_y_raw"] = projected[:, 1]
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

    metric_started = time.perf_counter()
    metrics = compute_metrics(
        representation_matrix,
        projected,
        grid,
        values=values,
        labels=label_values,
        k=k,
        trustworthiness_sample_size=trustworthiness_sample_size,
        random_state=random_state,
        representation_metric=representation_distance_metric,
        coordinate_scaling=coordinate_scaling,
    )
    metric_seconds = time.perf_counter() - metric_started
    metrics.insert(0, "projection_method", projection_name)
    metrics.insert(0, "representation_type", representation_name)
    metrics.insert(0, "map_name", name)
    metrics["projection_distance_metric"] = (
        selected_umap_metric if projection_name == "umap" else "not_applicable"
    )
    metrics["representation_distance_metric"] = representation_distance_metric
    metrics["grid_assignment_method"] = selected_assignment_method
    metrics["grid_candidate_cells"] = int(grid_rows * grid_cols)
    metrics["grid_side"] = int(grid_rows)
    metrics["grid_sizing_strategy"] = (
        "target_occupancy" if grid_padding is None else "explicit_padding"
    )
    metrics["grid_target_occupancy"] = (
        float(grid_occupancy) if grid_padding is None else np.nan
    )
    metrics["grid_actual_occupancy"] = float(len(projected) / (grid_rows * grid_cols))
    metrics["grid_padding_cells"] = (
        np.nan if grid_padding is None else int(grid_padding)
    )
    metrics["coordinate_scaling"] = str(coordinate_scaling).strip().lower()
    for key, value in assignment_diagnostics.items():
        metrics["assignment_" + key] = json.dumps(value) if isinstance(value, (dict, list)) else value

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
        render_started = time.perf_counter()
        output_files.update(
            render_grid_map(
                result_data,
                output_dir / name,
                color_map=color_map,
                detail=render_detail,
                formats=output_formats,
            )
        )
        manifest = {
            "map_name": name,
            "input_table_sha256": hashlib.sha256(data.to_csv(index=False).encode()).hexdigest(),
            "source": provenance or {},
            "environment": environment_versions(),
            "parameters": {
                "representation": representation_name, "projection": projection_name,
                "umap_metric": selected_umap_metric, "umap_n_neighbors": min(15, max(2, len(data) - 1)),
                "umap_min_dist": 0.1, "random_state": random_state,
                "morgan_radius": radius, "morgan_n_bits": n_bits,
                "coordinate_scaling": coordinate_scaling,
                "grid_occupancy": grid_occupancy, "grid_padding": grid_padding,
                "assignment_method": selected_assignment_method, "sparse_neighbors": sparse_neighbors,
                "sparse_adaptive": sparse_adaptive, "sparse_max_neighbors": sparse_max_neighbors,
                "k": k, "trustworthiness_sample_size": trustworthiness_sample_size,
                "annotation_neighborhood_policy": "tie_inclusive_kth_distance_shell",
                "lower_threshold": lower_threshold, "upper_threshold": upper_threshold,
                "duplicate_policy": duplicate_policy, "render_detail": render_detail,
            },
            "assignment": assignment_diagnostics,
            "timings_seconds": {"preparation_representation_projection": projection_seconds,
                                "assignment": assignment_seconds, "metrics": metric_seconds,
                                "render_export": time.perf_counter() - render_started,
                                "total_map_workflow": time.perf_counter() - started},
            "outputs": {key: Path(value).name for key, value in output_files.items()},
        }
        path = output_dir / (name + "_run_manifest.json")
        path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        output_files["run_manifest"] = str(path)

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
