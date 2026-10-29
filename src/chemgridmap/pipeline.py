"""High-level grid chemical map workflow."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import hashlib
import json
import time
from typing import Dict, List, Optional, Sequence

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
from .evidence import DOCUMENT_OFFSET_COLUMNS, prepare_generic_activity_records
from .projection import project_2d
from .refinement import refine_grid_layout
from .reliability import grid_neighbor_sets, reliability_columns, representation_neighbor_sets
from .representations import compute_representation
from .provenance import environment_versions, file_digests
from .comparison import structural_groups, save_group_outputs


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
    structure_groups: List[Dict] = field(default_factory=list)
    grouping_protocol: Dict = field(default_factory=dict)


def build_grid_map(
    data: pd.DataFrame,
    output_dir: Optional[Path] = None,
    name: str = "grid_map",
    smiles_col: str = "smiles",
    value_col: Optional[str] = None,
    label_col: Optional[str] = None,
    representation: str = "morgan",
    projection: str = "umap",
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
    grouping: str = "butina",
    group_similarity: float = 0.6,
    cluster_col: Optional[str] = None,
    cluster_source: Optional[str] = None,
    grid_refinement: str = "neighbors",
    core_rule: str = "coverage",
    core_min_coverage: float = 0.8,
    fmcs_timeout: int = 1,
    compound_col: Optional[str] = None,
    document_col: Optional[str] = None,
    relation_col: Optional[str] = None,
    evidence_model: str = "gaussian",
    merge_identical_documents: bool = True,
    map_color_modes: Optional[List[str]] = None,
    map_color_formats: Optional[Sequence[str]] = None,
) -> GridMapResult:
    """Build and optionally save a one-molecule-per-cell grid chemical map.

    ``grid_refinement="neighbors"`` refines the minimum-displacement assignment
    so that representation-space neighbours sit closer on the lattice. Every
    run reports per-molecule neighbour reliability (``neighbor_recall``).
    Extra color layers use ``map_color_formats`` (SVG/PNG); ``None`` keeps SVG
    plus PNG when enabled by ``output_formats``. Main-map formats are unchanged.
    """
    started = time.perf_counter()
    if compound_col:
        if not value_col:
            raise ValueError("value_col is required with compound_col.")
        data, _, evidence_fit, generic_report = prepare_generic_activity_records(
            data,
            compound_column=compound_col,
            value_column=value_col,
            document_column=document_col,
            smiles_column=smiles_col,
            relation_column=relation_col,
            lower=lower_threshold,
            upper=upper_threshold,
            conflict_range_threshold=1.0,
            confidence=bool(document_col),
            merge_identical=merge_identical_documents,
            evidence_model=evidence_model,
        )
        smiles_col, value_col, label_col = (
            "canonical_smiles", "activity_pchembl", "activity_class"
        )
        provenance = dict(provenance or {})
        provenance["generic_activity"] = generic_report
        if evidence_fit is not None:
            provenance["activity_confidence"] = evidence_fit
    if bool(x_col) != bool(y_col):
        raise ValueError("x_col and y_col must be supplied together.")
    if grid_refinement not in {"none", "neighbors"}:
        raise ValueError("grid_refinement must be 'none' or 'neighbors'.")
    map_color_modes = list(dict.fromkeys(map_color_modes or []))
    if set(map_color_modes) - {"activity", "class_support", "neighbor_recall", "structure_group"}:
        raise ValueError("Unknown map color mode.")
    if map_color_formats is not None:
        map_color_formats = list(dict.fromkeys(map_color_formats))
        if set(map_color_formats) - {"svg", "png"}:
            raise ValueError("map_color_formats only supports 'svg' and 'png'.")
    confidence_fit = data.attrs.get("activity_confidence_fit") or (provenance or {}).get("activity_confidence")
    reliability_parameters = {"activity_confidence": confidence_fit is not None,
                              "merge_identical_documents": True}
    reliability_parameters.update((provenance or {}).get("reliability_parameters", {}))
    if confidence_fit is not None:
        reliability_parameters["merge_identical_documents"] = confidence_fit["merge_identical_documents"]

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

    projection_fallback = None
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
        effective_projection = projection.lower()
        if effective_projection == "umap" and len(representation_matrix) < 5:
            # UMAP needs a neighborhood large enough for its spectral initialization.
            effective_projection = "pca"
            projection_fallback = "fewer_than_5_molecules"
        projected = project_2d(
            representation_matrix,
            method=effective_projection,
            random_state=random_state,
            umap_metric=selected_umap_metric,
        )
        if effective_projection == "umap" and not np.isfinite(projected).all():
            effective_projection = "pca"
            projection_fallback = "non_finite_umap_coordinates"
            projected = project_2d(
                representation_matrix,
                method=effective_projection,
                random_state=random_state,
                umap_metric=selected_umap_metric,
            )
        representation_name = representation.lower()
        projection_name = effective_projection
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
    reliability_started = time.perf_counter()
    reference_neighbors = representation_neighbor_sets(
        representation_matrix, k=k,
        metric="jaccard" if representation_distance_metric == "jaccard" else "euclidean")
    refinement_diagnostics = {"method": "none"}
    if grid_refinement == "neighbors":
        cells = np.rint(grid * np.array([grid_cols - 1, grid_rows - 1])).astype(int)
        refined = refine_grid_layout(cells, projected_01 * np.array([grid_cols - 1, grid_rows - 1]),
                                     reference_neighbors, int(grid_rows), random_state=random_state,
                                     diagnostics=refinement_diagnostics)
        grid = refined / np.array([grid_cols - 1, grid_rows - 1], dtype=np.float64)
    reliability = reliability_columns(
        reference_neighbors,
        grid_neighbor_sets(np.rint(grid * np.array([grid_cols - 1, grid_rows - 1])), k=k))
    reliability_seconds = time.perf_counter() - reliability_started
    result_data = prepared.copy()
    result_data["projection_x_raw"] = projected[:, 0]
    result_data["projection_y_raw"] = projected[:, 1]
    result_data["projection_x"] = projected_01[:, 0]
    result_data["projection_y"] = projected_01[:, 1]
    result_data["grid_x"] = grid[:, 0]
    result_data["grid_y"] = grid[:, 1]
    result_data["grid_col"] = np.rint(grid[:, 0] * (grid_cols - 1)).astype(int)
    result_data["grid_row"] = np.rint(grid[:, 1] * (grid_rows - 1)).astype(int)
    for column, values in reliability.items():
        result_data[column] = values

    if grouping not in {"butina", "none", "external"}:
        raise ValueError("Grouping must be butina, external, or none.")
    groups, grouping_protocol = [], {"method": "none"}
    if grouping != "none":
        if grouping == "external" and not cluster_col:
            raise ValueError("External grouping requires cluster_col and cluster_source.")
        result_data, groups, grouping_protocol = structural_groups(
            result_data, similarity=group_similarity, radius=radius, n_bits=n_bits,
            cluster_col=cluster_col if grouping == "external" else None,
            cluster_source=cluster_source if grouping == "external" else None,
            core_rule=core_rule, core_min_coverage=core_min_coverage,
            fmcs_timeout=fmcs_timeout)

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
    metrics["projection_fallback"] = projection_fallback
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
    metrics["grid_refinement"] = grid_refinement
    metrics["representation_grid_tie_aware_recall"] = float(np.nanmean(reliability["neighbor_recall"]))
    for key, value in assignment_diagnostics.items():
        metrics["assignment_" + key] = json.dumps(value) if isinstance(value, (dict, list)) else value

    output_files = {}
    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        coordinates_path = output_dir / "{}_grid_coordinates.csv".format(name)
        metrics_path = output_dir / "{}_metrics.csv".format(name)
        result_data.to_csv(coordinates_path, index=False, encoding="utf-8", lineterminator="\n")
        metrics.to_csv(metrics_path, index=False, encoding="utf-8", lineterminator="\n")
        output_files.update(
            {
                "coordinates": str(coordinates_path),
                "metrics": str(metrics_path),
            }
        )
        if groups:
            output_files.update(save_group_outputs(result_data, groups, grouping_protocol, output_dir, name))
        if confidence_fit is not None:
            offsets = pd.DataFrame(confidence_fit["document_offsets"], columns=DOCUMENT_OFFSET_COLUMNS)
            offsets["total_evaluated_pairs"] = confidence_fit["document_pair_count"]
            offsets["total_observed_pairs"] = confidence_fit["observed_document_pair_count"]
            offsets["total_possible_pairs"] = confidence_fit["document_pair_possible_count"]
            offset_path = output_dir / (name + "_document_offsets.csv")
            offsets.to_csv(offset_path, index=False, encoding="utf-8", lineterminator="\n")
            output_files["document_offsets"] = str(offset_path)
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
        mode_formats = map_color_formats if map_color_formats is not None else (
            ["svg"] + (["png"] if output_formats is None or "png" in output_formats else []))
        for mode in map_color_modes:
            mode_files = render_grid_map(result_data, output_dir / (name + "_" + mode),
                                         color_map=color_map, detail=render_detail,
                                         formats=mode_formats, color_by=mode)
            output_files.update({mode + "_" + key: path for key, path in mode_files.items()})
        from . import __version__
        manifest = {
            "software_version": __version__,
            "map_name": name,
            "input_table_sha256": hashlib.sha256(
                data.to_csv(index=False, lineterminator="\n").encode("utf-8")
            ).hexdigest(),
            "source": provenance or {},
            "environment": environment_versions(),
            "parameters": {
                "representation": representation_name, "projection": projection_name,
                "projection_fallback": projection_fallback,
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
                "grid_refinement": grid_refinement, "core_rule": core_rule,
                "core_min_coverage": core_min_coverage, "fmcs_timeout": fmcs_timeout,
                "evidence_model": evidence_model,
                "merge_identical_documents": merge_identical_documents,
                "map_color_modes": map_color_modes,
                "map_color_formats": map_color_formats,
                **reliability_parameters,
            },
            "activity_confidence": confidence_fit,
            "assignment": assignment_diagnostics,
            "grid_refinement": refinement_diagnostics,
            "structure_grouping": grouping_protocol,
            "timings_seconds": {"preparation_representation_projection": projection_seconds,
                                "assignment": assignment_seconds,
                                "neighbor_reliability_and_refinement": reliability_seconds,
                                "metrics": metric_seconds,
                                "render_export": time.perf_counter() - render_started,
                                "total_map_workflow": time.perf_counter() - started},
            "outputs": {key: Path(value).name for key, value in output_files.items()},
            "output_digests": file_digests(output_files),
        }
        path = output_dir / (name + "_run_manifest.json")
        path.write_text(json.dumps(manifest, indent=2), encoding="utf-8", newline="\n")
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
        structure_groups=groups,
        grouping_protocol=grouping_protocol,
    )
