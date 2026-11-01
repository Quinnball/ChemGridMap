"""生成图 2–5，渲染 PNG，并执行完整 QA。"""

from __future__ import annotations

import os
import sys
import math
from pathlib import Path


SOURCE_DIR = Path(__file__).resolve().parent
FIG_DIR = SOURCE_DIR.parent
os.environ.setdefault("MPLCONFIGDIR", str(FIG_DIR / ".mplconfig"))
os.environ.setdefault("SOURCE_DATE_EPOCH", "1793318400")  # 2026-10-30T00:00:00Z
sys.dont_write_bytecode = True

import matplotlib.pyplot as plt
import pymupdf

import fig2
import fig3
import fig4
import fig5
import fig6
import molecule
import qa_figures


MM_PER_INCH = 25.4


def _save_pdf(fig, path: Path) -> None:
    fig.savefig(
        path,
        format="pdf",
        dpi=600,
        facecolor="white",
        edgecolor="white",
        bbox_inches=None,
        metadata={"Creator": "F1 Matplotlib figure build"},
    )


def _render_pdf(pdf_path: Path, png_path: Path, dpi: int) -> None:
    with pymupdf.open(pdf_path) as pdf:
        pixmap = pdf[0].get_pixmap(
            dpi=dpi,
            alpha=False,
            colorspace=pymupdf.csRGB,
        )
        pixmap.set_dpi(dpi, dpi)
        pixmap.save(png_path)


def _is_nonwhite(color) -> bool:
    return color is not None and any(float(channel) < 1.0 - 1e-6 for channel in color)


def _pdf_text_sequence(page) -> list[str]:
    return [
        span["text"]
        for block in page.get_text("dict")["blocks"]
        if block["type"] == 0
        for line in block["lines"]
        for span in line["spans"]
        if span["text"]
    ]


def _crop_pdf_bottom(pdf_path: Path) -> dict[str, object]:
    """按非白色文字与绘图的最低点裁去页面底部空白。"""
    temp_path = pdf_path.with_name(f".{pdf_path.name}.crop.tmp")
    with pymupdf.open(pdf_path) as source_doc:
        source = source_doc[0]
        content_rects = []
        for block in source.get_text("dict")["blocks"]:
            if block["type"] == 0:
                for line in block["lines"]:
                    for span in line["spans"]:
                        if span["text"].strip():
                            content_rects.append(pymupdf.Rect(span["bbox"]))
        for drawing in source.get_drawings():
            if _is_nonwhite(drawing.get("fill")) or _is_nonwhite(drawing.get("color")):
                content_rects.append(pymupdf.Rect(drawing["rect"]))

        if not content_rects:
            raise RuntimeError(f"{pdf_path.name} 未找到可见文字或绘图")

        content_bottom_mm = max(rect.y1 for rect in content_rects) * MM_PER_INCH / 72
        crop_height_mm = math.ceil((content_bottom_mm + 1.5) * 2) / 2
        source_width_mm = source.rect.width * MM_PER_INCH / 72
        source_height_mm = source.rect.height * MM_PER_INCH / 72
        if not 0 < crop_height_mm < source_height_mm:
            raise RuntimeError(
                f"{pdf_path.name} 裁切高度 {crop_height_mm} mm 不在原页面范围内"
            )

        crop_height_pt = crop_height_mm * 72 / MM_PER_INCH
        stripped_items = [rect for rect in content_rects if rect.y1 > crop_height_pt + 1e-6]
        text_before = _pdf_text_sequence(source)
        clip = pymupdf.Rect(0, 0, source.rect.width, crop_height_pt)
        cropped = pymupdf.open()
        target = cropped.new_page(width=170.0 * 72 / MM_PER_INCH, height=crop_height_pt)
        target.show_pdf_page(target.rect, source_doc, 0, clip=clip)
        cropped.set_metadata(source_doc.metadata)
        cropped.save(
            temp_path,
            garbage=4,
            deflate=True,
            no_new_id=True,
            reproducible=True,
        )
        cropped.close()

    os.replace(temp_path, pdf_path)
    return {
        "source_size_mm": [source_width_mm, source_height_mm],
        "content_bottom_mm": content_bottom_mm,
        "height_mm": crop_height_mm,
        "bottom_padding_mm": crop_height_mm - content_bottom_mm,
        "cropped_strip_item_count": len(stripped_items),
        "text_before": text_before,
    }


def main() -> int:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    output_dir = FIG_DIR / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    qa_dir = FIG_DIR / "qa"
    qa_dir.mkdir(parents=True, exist_ok=True)
    (FIG_DIR / ".mplconfig").mkdir(parents=True, exist_ok=True)

    figure_2, figure_2_placements, figure_2_assertions = fig2.build()
    figures = {
        2: figure_2,
        3: fig3.build(),
        4: fig4.build(),
    }
    figure_5, smiles = fig5.build()
    figures[5] = figure_5
    figure_6, figure_6_audit = fig6.build()
    figures[6] = figure_6

    pdf_paths = {number: output_dir / f"Figure_{number}.pdf" for number in figures}
    png_paths = {number: output_dir / f"Figure_{number}.png" for number in figures}
    preview_paths = {number: qa_dir / f"preview_Figure_{number}.png" for number in figures}

    for number, fig in figures.items():
        _save_pdf(fig, pdf_paths[number])
    fig6.insert_images(pdf_paths[6], figure_6_audit)

    # 页面顶部坐标：150 - (100.5 + 14)。
    molecule_measurements = {
        "Figure_5c": molecule.insert_molecule(
            pdf_paths[5],
            smiles,
            (114.0, 35.5, 27.0, 14.0),
        )
    }
    figure_2_measurements = [
        molecule.insert_fixed_bond_molecule(
            pdf_paths[2],
            placement["smiles"],
            center_x_mm=float(placement["x"]),
            center_y_mm=float(placement["center_y"]),
            canvas_height_mm=fig2.SOURCE_HEIGHT_MM,
            printed_bond_length_mm=2.4,
        )
        for placement in figure_2_placements
    ]
    molecule_measurements["Figure_2"] = figure_2_measurements

    crop_audits = {
        number: _crop_pdf_bottom(pdf_paths[number])
        for number in (2, 3, 4, 5)
    }
    with pymupdf.open(pdf_paths[6]) as pdf:
        figure_6_text = _pdf_text_sequence(pdf[0])
    crop_audits[6] = {
        "source_size_mm": [fig6.WIDTH_MM, fig6.HEIGHT_MM],
        "content_bottom_mm": fig6.HEIGHT_MM - 1.5,
        "height_mm": fig6.HEIGHT_MM,
        "bottom_padding_mm": 1.5,
        "cropped_strip_item_count": 0,
        "text_before": figure_6_text,
    }

    for number in figures:
        _render_pdf(pdf_paths[number], png_paths[number], 600)
        _render_pdf(pdf_paths[number], preview_paths[number], 200)

    supplementary_checks = {
        "assertions": [
            *figure_2_assertions,
            qa_figures.figure_2_supplementary(
                pdf_paths[2],
                figure_2_measurements,
            ),
        ],
        "passed": bool(figure_6_audit["passed"]),
        "figure_6": {
            "skipped": False,
            "reason": None,
            "build": figure_6_audit,
        },
    }
    report = qa_figures.run_qa(
        figures,
        pdf_paths,
        png_paths,
        preview_paths,
        molecule_measurements,
        crop_audits,
        supplementary_checks,
        allowed_raster_counts={2: 0, 3: 0, 4: 0, 5: 0, 6: 4},
    )
    for number, fig in figures.items():
        plt.close(fig)
    heights = {number: crop_audits[number]["height_mm"] for number in figures}
    print(
        f"F5 all_passed={report['all_passed']} crop_heights_mm={heights}"
    )
    return 0 if report["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
