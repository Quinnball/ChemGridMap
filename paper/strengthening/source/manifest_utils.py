"""Write portable analysis manifests."""
from __future__ import annotations

from datetime import datetime
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
from typing import Iterable, Mapping

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
SOFTWARE_PACKAGES = (
    "numpy",
    "pandas",
    "scipy",
    "scikit-learn",
    "umap-learn",
    "pynndescent",
    "numba",
    "rdkit",
    "chemgridmap",
)


def _json_safe(value):
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _relative_path(path: Path) -> str:
    resolved = Path(path).resolve()
    return Path(os.path.relpath(resolved, ROOT.resolve())).as_posix()


def _paths(values: Iterable[Path | str]) -> list[Path]:
    if isinstance(values, (str, Path)):
        values = (values,)
    return [Path(value) for value in values]


def _file_hashes(
    paths: Iterable[Path | str],
    *,
    base: Path = ROOT,
) -> dict[str, str]:
    return {
        Path(os.path.relpath(Path(path).resolve(), base.resolve())).as_posix(): _sha256(path)
        for path in _paths(paths)
    }


def _software_versions() -> dict[str, str]:
    versions = {}
    for package in SOFTWARE_PACKAGES:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not installed"
    return {"python": platform.python_version(), **versions}


def write_manifest(
    out_dir: Path | str,
    *,
    stage: str,
    inputs: Iterable[Path | str],
    sources: Iterable[Path | str],
    parameters: Mapping,
    started_at: str,
    finished_at: str,
) -> dict:
    """Write one manifest whose scientific paths are repository-relative."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    output_paths = [
        path
        for path in sorted(out_dir.rglob("*"))
        if path.is_file() and path.name != "manifest.json"
    ]
    start = datetime.fromisoformat(started_at)
    finish = datetime.fromisoformat(finished_at)
    manifest_parameters = dict(parameters)
    analysis_role = str(manifest_parameters.pop("analysis_role", "post hoc"))
    manifest = {
        "stage": stage,
        "analysis_role": analysis_role,
        "started_at": started_at,
        "finished_at": finished_at,
        "elapsed_seconds": float((finish - start).total_seconds()),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "software_versions": _software_versions(),
        "parameters": manifest_parameters,
        "inputs_sha256": _file_hashes(inputs),
        "source_sha256": _file_hashes(sources),
        "outputs_sha256": _file_hashes(output_paths, base=out_dir),
    }
    (out_dir / "manifest.json").write_text(
        json.dumps(
            _json_safe(manifest),
            indent=2,
            ensure_ascii=False,
            allow_nan=False,
        )
        + "\n",
        encoding="utf-8",
    )
    return manifest
