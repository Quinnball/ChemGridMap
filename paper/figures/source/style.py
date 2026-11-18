"""Shared manuscript style for Figures 3-5."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
from matplotlib import font_manager


# Subscripts of 7 pt math text should reach 6 pt; 0.86 gives about 6.0 pt.
mpl._mathtext.SHRINK_FACTOR = 0.86

FONT_DIR = Path("/System/Library/Fonts/Supplemental")
FONT_FILES = [
    FONT_DIR / "Arial.ttf",
    FONT_DIR / "Arial Bold.ttf",
    FONT_DIR / "Arial Italic.ttf",
]
# Liberation Sans has the same metrics as Arial; it is used where Arial is not installed
# (for example on Linux CI), so that the layout and every QA measurement stay the same.
FALLBACK_DIR = Path("/usr/share/fonts/truetype/liberation")
FALLBACK_FILES = [
    FALLBACK_DIR / "LiberationSans-Regular.ttf",
    FALLBACK_DIR / "LiberationSans-Bold.ttf",
    FALLBACK_DIR / "LiberationSans-Italic.ttf",
]
FONT_NAME = "Arial"
_FONTS_REGISTERED = False


def apply_style() -> None:
    """Register Arial and configure the shared print style."""
    global _FONTS_REGISTERED, FONT_NAME
    if not _FONTS_REGISTERED:
        files = FONT_FILES
        if not all(path.is_file() for path in FONT_FILES):
            if not all(path.is_file() for path in FALLBACK_FILES):
                raise FileNotFoundError(f"Missing Arial font: {FONT_FILES[0]} (and no Liberation Sans fallback)")
            files, FONT_NAME = FALLBACK_FILES, "Liberation Sans"
        for path in files:
            font_manager.fontManager.addfont(path)
        _FONTS_REGISTERED = True

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [FONT_NAME],
            "mathtext.fontset": "custom",
            "mathtext.rm": FONT_NAME,
            "mathtext.it": f"{FONT_NAME}:italic",
            "mathtext.bf": f"{FONT_NAME}:bold",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.unicode_minus": True,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.dpi": 600,
            "axes.grid": False,
        }
    )


if __name__ == "__main__":
    apply_style()
