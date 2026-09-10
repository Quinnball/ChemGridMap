"""Fail clearly before tests or paper experiments use the wrong interpreter."""

import importlib.metadata
import importlib.util
from pathlib import Path
import sys


def main():
    if sys.version_info < (3, 10):
        raise SystemExit("Use the project .venv (Python >= 3.10), not legacy Python.")
    missing = [name for name in ("pytest", "rdkit", "chembl_structure_pipeline",
                                "pandas", "sklearn", "umap")
               if importlib.util.find_spec(name) is None]
    if missing:
        raise SystemExit("Missing: {}. Run make setup first.".format(", ".join(missing)))
    print("Python:", sys.executable)
    print("Environment:", Path(sys.prefix))
    for name in ("pytest", "rdkit", "chembl_structure_pipeline", "numpy",
                 "scipy", "scikit-learn", "umap-learn"):
        print("{}=={}".format(name, importlib.metadata.version(name)))


if __name__ == "__main__":
    main()
