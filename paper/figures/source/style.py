"""F1 稿件图 3–5 的统一样式。"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
from matplotlib import font_manager


# 7 pt 数学文字的下标需达到 6 pt；0.86 使缩小后的下标约为 6.0 pt。
mpl._mathtext.SHRINK_FACTOR = 0.86

FONT_DIR = Path("/System/Library/Fonts/Supplemental")
FONT_FILES = [
    FONT_DIR / "Arial.ttf",
    FONT_DIR / "Arial Bold.ttf",
    FONT_DIR / "Arial Italic.ttf",
]
_FONTS_REGISTERED = False


def apply_style() -> None:
    """注册 Arial，并设置所有图共用的印刷样式。"""
    global _FONTS_REGISTERED
    if not _FONTS_REGISTERED:
        for path in FONT_FILES:
            if not path.is_file():
                raise FileNotFoundError(f"缺少 Arial 字体：{path}")
            font_manager.fontManager.addfont(path)
        _FONTS_REGISTERED = True

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial"],
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
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
