"""Build the v0.7 Figure 4 and the Additional file 1 hERG-pair figure (PDF and 600-dpi PNG)."""
from __future__ import annotations

import os
import sys
from pathlib import Path

SOURCE_DIR = Path(__file__).resolve().parent
FIG_DIR = SOURCE_DIR.parent
os.environ.setdefault("MPLCONFIGDIR", str(FIG_DIR / ".mplconfig"))
os.environ.setdefault("SOURCE_DATE_EPOCH", "1795046400")  # 2026-11-19T00:00:00Z
sys.dont_write_bytecode = True

import build_figures as bf  # noqa: E402
import fig4  # noqa: E402
import fig6  # noqa: E402
import fig_si_herg_pair  # noqa: E402


def main() -> int:
    out = FIG_DIR / "output_v07"
    out.mkdir(parents=True, exist_ok=True)
    (FIG_DIR / ".mplconfig").mkdir(parents=True, exist_ok=True)
    only = set(sys.argv[1:])
    builders = {"Figure_4": lambda: (fig4.build(), None), "Figure_A_hERG_pair": lambda: (fig_si_herg_pair.build(), None),
                "Figure_6": fig6.build}
    for name, builder in builders.items():
        if only and name not in only:
            continue
        fig, audit = builder()
        pdf, png = out / f"{name}.pdf", out / f"{name}.png"
        bf._save_pdf(fig, pdf)
        if name == "Figure_6":
            fig6.insert_images(pdf, audit)
        bf._crop_pdf_bottom(pdf)
        bf._render_pdf(pdf, png, 600)
        print(name, pdf.stat().st_size)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
