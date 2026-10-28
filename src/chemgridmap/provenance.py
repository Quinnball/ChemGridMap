"""Small, portable manifests for reproducing a map or an experiment."""

from importlib import metadata
import hashlib
from pathlib import Path
import os
import platform
import sys


def file_digests(files):
    return {key: {"file": Path(path).name, "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest()}
            for key, path in files.items()}


def environment_versions():
    packages = {}
    for name in ("chemgridmap", "rdkit", "chembl_structure_pipeline", "numpy",
                 "pandas", "scipy", "scikit-learn", "umap-learn", "matplotlib"):
        try:
            packages[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            packages[name] = "not_installed"
    return {"python": sys.version, "platform": platform.platform(),
            "machine": platform.machine(), "processor": platform.processor(),
            "logical_cpu_count": os.cpu_count(), "packages": packages}
