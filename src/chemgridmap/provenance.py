"""Small, portable manifests for reproducing a map or an experiment."""

from importlib import metadata
import hashlib
from pathlib import Path
import os
import platform
import subprocess
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


def source_tree_digest(root=None):
    """SHA-256 over sorted source paths and their individual file hashes."""
    root = Path(root or Path(__file__).resolve().parent)
    files = sorted(
        path for path in root.rglob("*")
        if path.is_file() and path.suffix.lower() in {".py", ".js", ".css", ".html"}
    )
    combined = hashlib.sha256()
    for path in files:
        relative = path.relative_to(root).as_posix().encode("utf-8")
        digest = hashlib.sha256(path.read_bytes()).digest()
        combined.update(relative)
        combined.update(b"\0")
        combined.update(digest)
    return {"source_sha256": combined.hexdigest(), "source_file_count": len(files)}


def git_provenance(root=None):
    """Return the current commit and dirty state when Git metadata is available."""
    package_dir = Path(root or Path(__file__).resolve().parent)
    try:
        repo_root = Path(subprocess.run(
            ["git", "-C", str(package_dir), "rev-parse", "--show-toplevel"],
            check=True, capture_output=True, text=True, timeout=10,
        ).stdout.strip())
        subprocess.run(
            ["git", "-C", str(package_dir), "ls-files", "--error-unmatch", "provenance.py"],
            check=True, capture_output=True, text=True, timeout=10,
        )
        commit = subprocess.run(
            ["git", "-C", str(repo_root), "rev-parse", "HEAD"],
            check=True, capture_output=True, text=True, timeout=10,
        ).stdout.strip()
        package_path = package_dir.resolve().relative_to(repo_root.resolve()).as_posix()
        status = subprocess.run(
            ["git", "-C", str(repo_root), "status", "--porcelain", "--", package_path],
            check=True, capture_output=True, text=True, timeout=10,
        ).stdout
        return {"git_commit": commit, "git_dirty": bool(status.strip())}
    except (OSError, subprocess.CalledProcessError, subprocess.TimeoutExpired, ValueError):
        return {"git_commit": None, "git_dirty": None}
