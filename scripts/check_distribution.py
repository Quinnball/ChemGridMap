"""Test the source distribution and wheel outside the editable checkout."""
from pathlib import Path
import json
import os
import subprocess
import sys
import tarfile
import tempfile
from chemgridmap import __version__

ROOT = Path(__file__).resolve().parents[1]


def run(args, cwd, env):
    result = subprocess.run(args, cwd=cwd, env=env, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(result.stdout + result.stderr)
    return result.stdout.strip()


def main():
    dist = ROOT / "dist"
    wheel = dist / f"chemgridmap-{__version__}-py3-none-any.whl"
    source = dist / f"chemgridmap-{__version__}.tar.gz"
    env = dict(os.environ, MPLCONFIGDIR=str(ROOT / ".cache/matplotlib"))
    with tempfile.TemporaryDirectory(prefix="chemgridmap-distribution-") as temp:
        directory = Path(temp)
        with tarfile.open(source) as archive:
            archive.extractall(directory, filter="data")
        checkout = directory / f"chemgridmap-{__version__}"
        required = ["Makefile", "scripts/check_environment.py", "requirements-paper.lock.txt",
                    "examples/chembl205_embedding_umap.csv", "paper/run_revision_validation.py"]
        assert all((checkout / name).is_file() for name in required)
        source_tests = run([sys.executable, "-m", "pytest", "-q"], checkout, env)
        installed = directory / "wheel_install"
        run([sys.executable, "-m", "pip", "install", "--no-deps", "--no-compile",
             "--target", str(installed), str(wheel)], directory, env)
        wheel_env = dict(env, PYTHONPATH=str(installed))
        import_path = run([sys.executable, "-c", "import chemgridmap; print(chemgridmap.__file__)"],
                          directory, wheel_env)
        assert Path(import_path).is_relative_to(installed)
        wheel_tests = run([sys.executable, "-m", "pytest", "-q", "-o", "pythonpath="],
                          checkout, wheel_env)
    result = {"passed": True, "python": sys.version, "source_distribution": source_tests,
              "wheel_distribution": wheel_tests,
              "wheel_import_outside_editable_checkout": True,
              "scope": "Both distributions tested with the existing locked dependency environment; not a remote cross-platform test."}
    (dist / "distribution_check.json").write_text(json.dumps(result, indent=2) + "\n")
    print(source_tests)
    print(wheel_tests)
    print("Source distribution and wheel checks passed.")


if __name__ == "__main__":
    main()
