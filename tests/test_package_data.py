"""Keep browser assets and bundled datasets available in installed packages."""

import ast
from fnmatch import fnmatchcase
from pathlib import Path
import re


def test_every_web_and_data_file_is_covered_by_package_data():
    root = Path(__file__).resolve().parents[1]
    configuration = (root / "pyproject.toml").read_text(encoding="utf-8")
    # This setting is a literal string list; use the stdlib on Python 3.10 too.
    match = re.search(
        r"(?ms)^\[tool\.setuptools\.package-data\]\s*chemgridmap\s*=\s*(\[.*?\])",
        configuration,
    )
    assert match is not None, "chemgridmap package-data must be configured"
    patterns = ast.literal_eval(match.group(1))
    package = root / "src" / "chemgridmap"
    resources = [
        path.relative_to(package).as_posix()
        for directory in ("web", "data")
        for path in (package / directory).rglob("*")
        if path.is_file()
    ]
    assert resources
    uncovered = [
        name for name in resources
        if not any(fnmatchcase(name, pattern) for pattern in patterns)
    ]
    assert not uncovered, f"Missing package-data rules: {uncovered}"
