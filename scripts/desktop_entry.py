"""Entry point for the frozen local application."""
import multiprocessing
import os
from pathlib import Path
import tempfile

if __name__ == "__main__":
    multiprocessing.freeze_support()
    cache = Path(tempfile.gettempdir()) / "chemgridmap-cache"
    cache.mkdir(exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache / "matplotlib"))
    os.environ.setdefault("NUMBA_CACHE_DIR", str(cache / "numba"))
    from chemgridmap.app import main
    raise SystemExit(main())
