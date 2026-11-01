"""图 6：地图、分子检查器、邻居与来源记录的四面板版式。"""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
from tempfile import TemporaryDirectory

import matplotlib.text as mtext
import pymupdf
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from PIL import Image, ImageDraw, ImageFont

import common as c


FIG_DIR = Path(__file__).resolve().parent.parent
SCREENSHOT_DIR = FIG_DIR / "data_fig6"
MAP_FILE = SCREENSHOT_DIR / "6a_map.png"
INSPECTOR_FILE = SCREENSHOT_DIR / "6b_inspector.png"
NEIGHBOR_FILE = SCREENSHOT_DIR / "6c_neighbors.png"
RECORDS_TABLE_FILE = SCREENSHOT_DIR / "6d_records_table.png"
CAPTURE_REPORT = SCREENSHOT_DIR / "capture_report.json"
Arial_BOLD = Path("/System/Library/Fonts/Supplemental/Arial Bold.ttf")
WIDTH_MM = 170.0
HEIGHT_MM = 220.0

ACTIVITY_COLORS = {
    "inactive": ("#C44E52", "#ebc0c1"),
    "medium": ("#E6A23C", "#f6dfb8"),
    "active": ("#2A9D8F", "#b8e0db"),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load_capture_report() -> dict[str, object]:
    report = json.loads(CAPTURE_REPORT.read_text(encoding="utf-8"))
    if report.get("status") != "passed":
        raise AssertionError(f"图 6 capture_report 状态不是 passed：{report.get('status')}")
    if not math.isclose(float(report.get("device_scale_factor", 0)), 3.0):
        raise AssertionError("图 6 capture_report 的 device_scale_factor 不是 3")
    expected = {
        item["file"]: item["sha256"]
        for item in report.get("screenshots", [])
        if isinstance(item, dict) and "file" in item and "sha256" in item
    }
    required_files = {
        MAP_FILE.name,
        INSPECTOR_FILE.name,
        NEIGHBOR_FILE.name,
        RECORDS_TABLE_FILE.name,
    }
    if not required_files <= set(expected):
        raise AssertionError(
            f"图 6 capture_report 缺少截图：{sorted(required_files - set(expected))}"
        )
    hashes = {}
    for name in sorted(required_files):
        path = SCREENSHOT_DIR / name
        actual = _sha256(path)
        if actual != expected[name]:
            raise AssertionError(f"{name} SHA-256 与 capture_report 不符")
        hashes[name] = actual
    return {"report": report, "hashes": hashes}


def _crop_from_css_rect(image: Image.Image, rect: dict[str, float], scale: float) -> Image.Image:
    left = round(float(rect["x"]) * scale)
    top = round(float(rect["y"]) * scale)
    right = round((float(rect["x"]) + float(rect["width"])) * scale)
    bottom = round((float(rect["y"]) + float(rect["height"])) * scale)
    return image.convert("RGB").crop((left, top, right, bottom)).copy()


def _draw_neighbor_labels(
    image: Image.Image,
    report: dict[str, object],
    scale: float,
    rendered_width_mm: float,
) -> list[dict[str, object]]:
    css_width = float(report["neighbor_clip_css"]["width"])
    font_px = round(
        6.5 * 25.4 / 72 * css_width * scale / rendered_width_mm
    )
    font = ImageFont.truetype(str(Arial_BOLD), font_px)
    draw = ImageDraw.Draw(image)
    records = []
    class_labels = {
        "neighbor-shared": "S",
        "neighbor-grid_only": "G",
        "neighbor-representation_only": "F",
    }
    for cell in report["neighbor_cells"]:
        letter = class_labels[cell["class"]]
        center = (
            (float(cell["x"]) + float(cell["width"]) / 2) * scale,
            (float(cell["y"]) + float(cell["height"]) / 2) * scale,
        )
        draw.text(
            center,
            letter,
            font=font,
            fill="#172536",
            stroke_width=3,
            stroke_fill="white",
            anchor="mm",
        )
        records.append({"letter": letter, "center_px": list(center)})
    for cell in report["selected_cell_in_clip_css"]:
        center = (
            (float(cell["x"]) + float(cell["width"]) / 2) * scale,
            (float(cell["y"]) + float(cell["height"]) / 2) * scale,
        )
        draw.text(
            center,
            "Q",
            font=font,
            fill="#172536",
            stroke_width=3,
            stroke_fill="white",
            anchor="mm",
        )
        records.append({"letter": "Q", "center_px": list(center)})
    return records


def _border_axis(fig, rect_mm):
    ax = c.add_axes_mm(fig, rect_mm, (WIDTH_MM, HEIGHT_MM))
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#DCE2E7")
        spine.set_linewidth(0.4)
    return ax


def _mm_rect(bbox, dpi: float, canvas_height_mm: float) -> tuple[float, float, float, float]:
    """Matplotlib 外框转成与版面矩形一致的左下原点毫米坐标。"""
    return (
        bbox.x0 * 25.4 / dpi,
        bbox.y0 * 25.4 / dpi,
        bbox.x1 * 25.4 / dpi,
        bbox.y1 * 25.4 / dpi,
    )


def _rect_intersection(a, b) -> bool:
    return min(a[2], b[2]) > max(a[0], b[0]) and min(a[3], b[3]) > max(a[1], b[1])


def _rect_gap(a, b) -> float:
    dx = max(a[0] - b[2], b[0] - a[2], 0.0)
    dy = max(a[1] - b[3], b[1] - a[3], 0.0)
    return math.hypot(dx, dy)


def _draw_map_legend(fig) -> None:
    entries = (
        ("Inactive (≤6)", "inactive", 6.5, 27.0),
        ("Medium (6, 7)", "medium", 38.0, 59.0),
        ("Active (≥7)", "active", 70.0, 84.0),
    )
    for label, key, label_x, swatch_x in entries:
        c.mark_legend_text(
            fig.text(
                label_x / WIDTH_MM,
                139.5 / HEIGHT_MM,
                label,
                ha="left",
                va="center",
                fontsize=6.5,
                color=c.INK,
            )
        )
        for offset, color in zip((0.0, 3.3), ACTIVITY_COLORS[key], strict=True):
            fig.add_artist(
                Rectangle(
                    (
                        (swatch_x + offset) / WIDTH_MM,
                        138.3 / HEIGHT_MM,
                    ),
                    2.5 / WIDTH_MM,
                    2.5 / HEIGHT_MM,
                    transform=fig.transFigure,
                    facecolor=color,
                    edgecolor="white",
                    linewidth=0.4,
                    clip_on=False,
                )
            )
    c.mark_legend_text(
        fig.text(
            6.5 / WIDTH_MM,
            135.8 / HEIGHT_MM,
            "Pale: uncertain (class support below 0.8, or an offset document pair)",
            ha="left",
            va="center",
            fontsize=6.5,
            color=c.INK,
        )
    )


def _draw_neighbor_legend(fig) -> None:
    c.mark_legend_text(
        fig.text(
            69.0 / WIDTH_MM,
            117.5 / HEIGHT_MM,
            "Neighbor outlines",
            ha="left",
            va="center",
            fontsize=6.5,
            color=c.INK,
        )
    )
    rows = (
        ("Shared", "#0072B2", "-", None),
        ("Grid only", "#D55E00", (0, (3, 2)), None),
        ("Fingerprint only", "#7A5195", (0, (1, 1.4)), None),
        ("Selected (Q)", "#193e91", "none", "s"),
    )
    for index, (label, color, linestyle, marker) in enumerate(rows):
        y = 113.0 - index * 5.5
        if marker is None:
            line = Line2D(
                (69.0 / WIDTH_MM, 78.0 / WIDTH_MM),
                (y / HEIGHT_MM, y / HEIGHT_MM),
                transform=fig.transFigure,
                color=color,
                linewidth=0.8,
                linestyle=linestyle,
                solid_capstyle="butt",
            )
        else:
            line = Line2D(
                (73.5 / WIDTH_MM,),
                (y / HEIGHT_MM,),
                transform=fig.transFigure,
                color=color,
                linestyle="none",
                marker=marker,
                markersize=3.2,
            )
        fig.add_artist(line)
        c.mark_legend_text(
            fig.text(
                80.0 / WIDTH_MM,
                y / HEIGHT_MM,
                label,
                ha="left",
                va="center",
                fontsize=6.5,
                color=c.INK,
            )
        )


def _layout_checks(fig, layout: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    image_rects = {
        panel: tuple(float(value) for value in item["rect_mm"])
        for panel, item in layout.items()
    }
    letter_boxes = {
        text._qa_panel_letter: _mm_rect(
            text.get_window_extent(renderer),
            fig.dpi,
            HEIGHT_MM,
        )
        for text in fig.findobj(mtext.Text)
        if hasattr(text, "_qa_panel_letter")
    }
    letter_image_intersections = [
        [letter, panel]
        for letter, box in letter_boxes.items()
        for panel, rect in image_rects.items()
        if _rect_intersection(box, rect)
    ]
    c_title = next(
        text for text in fig.findobj(mtext.Text) if hasattr(text, "_qa_fig6_c_title")
    )
    c_title_box = _mm_rect(c_title.get_window_extent(renderer), fig.dpi, HEIGHT_MM)
    c_title_intersections = [
        rect_name
        for rect_name, rect in {
            **{f"image_{panel}": value for panel, value in image_rects.items()},
            "letter_c": letter_boxes["c"],
        }.items()
        if _rect_intersection(c_title_box, rect)
    ]

    region_rects = {
        **{f"panel_{panel}": rect for panel, rect in image_rects.items()},
        "legend_a": (5.0, 134.0, 100.0, 142.0),
        "legend_c": (68.5, 94.0, 108.5, 120.0),
    }
    spacing_failures = []
    minimum_gap = math.inf
    for (first_name, first), (second_name, second) in itertools.combinations(
        region_rects.items(), 2
    ):
        gap = _rect_gap(first, second)
        minimum_gap = min(minimum_gap, gap)
        if gap < 2.0 - 1e-9:
            spacing_failures.append(
                {
                    "first": first_name,
                    "second": second_name,
                    "gap_mm": round(gap, 6),
                }
            )
    return [
        {
            "name": "panel letters outside all screenshot rectangles",
            "expected": [],
            "actual": letter_image_intersections,
            "passed": not letter_image_intersections,
        },
        {
            "name": "panel c title outside image and panel letter",
            "expected": [],
            "actual": c_title_intersections,
            "passed": not c_title_intersections,
        },
        {
            "name": "panels and legends separated by at least 2 mm",
            "expected": 2.0,
            "actual": {
                "minimum_gap_mm": round(minimum_gap, 6),
                "failures": spacing_failures,
            },
            "passed": not spacing_failures,
        },
    ]


def insert_images(pdf_path: Path, audit: dict[str, object]) -> None:
    """按审计后的毫米矩形嵌入裁好的截图。"""
    points_per_mm = 72 / 25.4

    def pdf_rect(rect_mm) -> pymupdf.Rect:
        x0, y0, x1, y1 = rect_mm
        return pymupdf.Rect(
            x0 * points_per_mm,
            (HEIGHT_MM - y1) * points_per_mm,
            x1 * points_per_mm,
            (HEIGHT_MM - y0) * points_per_mm,
        )

    with TemporaryDirectory(prefix="chemgridmap-fig6-f5r1-") as temporary:
        temporary_dir = Path(temporary)
        a_path = temporary_dir / "6a_map_viewport.png"
        b_path = temporary_dir / "6b_inspector_crop.png"
        c_path = temporary_dir / "6c_neighbors_labeled.png"
        a_path.write_bytes(Path(audit["prepared_files"]["a"]).read_bytes())
        b_path.write_bytes(Path(audit["prepared_files"]["b"]).read_bytes())
        c_path.write_bytes(Path(audit["prepared_files"]["c"]).read_bytes())
        source_files = {
            "a": a_path,
            "b": b_path,
            "c": c_path,
            "d": RECORDS_TABLE_FILE,
        }
        with pymupdf.open(pdf_path) as pdf:
            page = pdf[0]
            for panel in ("a", "b", "c", "d"):
                page.insert_image(
                    pdf_rect(audit["layout_mm"][panel]),
                    filename=str(source_files[panel]),
                )
            temporary_pdf = pdf_path.with_name(f".{pdf_path.name}.images.tmp")
            pdf.save(
                temporary_pdf,
                garbage=4,
                deflate=True,
                no_new_id=True,
                reproducible=True,
            )
        temporary_pdf.replace(pdf_path)

    actual_rects = {}
    with pymupdf.open(pdf_path) as pdf:
        page = pdf[0]
        for image in page.get_images(full=True):
            width_px = int(image[2])
            panel = next(
                (
                    name
                    for name, item in audit["panels"].items()
                    if int(item["pixel_width"]) == width_px
                ),
                None,
            )
            if panel is None:
                continue
            rects = page.get_image_rects(image[0])
            if len(rects) != 1:
                raise AssertionError(f"图 6{panel} 的 PDF 图像矩形数为 {len(rects)}")
            rect = rects[0]
            actual_rects[panel] = [
                float(rect.x0 * 25.4 / 72),
                float((page.rect.height - rect.y1) * 25.4 / 72),
                float(rect.x1 * 25.4 / 72),
                float((page.rect.height - rect.y0) * 25.4 / 72),
            ]
    rects_match = all(
        all(
            abs(actual - expected) <= 0.1
            for actual, expected in zip(
                actual_rects.get(panel, []),
                audit["layout_mm"][panel],
                strict=False,
            )
        )
        for panel in ("a", "b", "c", "d")
    )
    check = {
        "name": "PDF screenshot rectangles match the mm layout",
        "expected": {panel: audit["layout_mm"][panel] for panel in ("a", "b", "c", "d")},
        "actual": actual_rects,
        "passed": set(actual_rects) == {"a", "b", "c", "d"} and rects_match,
    }
    audit["checks"].append(check)
    audit["pdf_image_rects_mm"] = actual_rects
    audit["passed"] = bool(audit["passed"] and check["passed"])


def _write_plotted_values(layout: dict[str, dict[str, object]]) -> None:
    path = c.QA_DIR / "Figure_6_plotted_values.csv"
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=[
                "figure",
                "panel",
                "series",
                "row_identifier",
                "variable",
                "source_file",
                "value",
            ],
        )
        writer.writeheader()
        for panel in ("a", "b", "c", "d"):
            item = layout[panel]
            for variable in (
                "pixel_width",
                "pixel_height",
                "rendered_width_mm",
                "rendered_height_mm",
                "effective_dpi",
            ):
                writer.writerow(
                    {
                        "figure": 6,
                        "panel": panel,
                        "series": "screenshot",
                        "row_identifier": item["file"],
                        "variable": variable,
                        "source_file": item["file"],
                        "value": float(item[variable]),
                    }
                )


def build() -> tuple[object, dict[str, object]]:
    c.init()
    prepared = _load_capture_report()
    report: dict[str, object] = prepared["report"]
    hashes: dict[str, str] = prepared["hashes"]
    scale = float(report["device_scale_factor"])

    with Image.open(MAP_FILE) as image:
        map_image = _crop_from_css_rect(
            image,
            report["map_viewport_in_clip_css"],
            scale,
        )
    with Image.open(INSPECTOR_FILE) as image:
        inspector_top_px = math.floor(
            (float(report["heading_top_css_px"]) - 8.0) * scale
        )
        inspector_bottom_px = math.ceil(
            (float(report["evidence_bottom_css_px"]) + 8.0) * scale
        )
        inspector_image = image.convert("RGB").crop(
            (
                0,
                inspector_top_px,
                image.width,
                inspector_bottom_px,
            )
        ).copy()
    with Image.open(NEIGHBOR_FILE) as image:
        neighbor_image = image.convert("RGB").copy()
    with Image.open(RECORDS_TABLE_FILE) as image:
        records_image = image.convert("RGB").copy()

    prepared_dir = c.QA_DIR / "F6_fig6_crops"
    prepared_dir.mkdir(parents=True, exist_ok=True)
    a_prepared = prepared_dir / "6a_map_viewport.png"
    b_prepared = prepared_dir / "6b_inspector_crop.png"
    c_prepared = prepared_dir / "6c_neighbors_labeled.png"
    map_image.save(a_prepared)
    inspector_image.save(b_prepared)

    map_width_mm = 98.0
    map_height_mm = map_width_mm * map_image.height / map_image.width
    inspector_width_mm = 58.0
    inspector_height_mm = (
        inspector_width_mm * inspector_image.height / inspector_image.width
    )
    neighbor_size_mm = 60.0
    neighbor_labels = _draw_neighbor_labels(
        neighbor_image,
        report,
        scale,
        neighbor_size_mm,
    )
    neighbor_image.save(c_prepared)
    records_width_mm = 117.0
    records_height_mm = (
        records_width_mm * records_image.height / records_image.width
    )
    inspector_y = 50.9
    layout = {
        "a": {
            "file": MAP_FILE.name,
            "rect_mm": (
                5.5,
                146.0,
                5.5 + map_width_mm,
                146.0 + map_height_mm,
            ),
            "pixel_width": map_image.width,
            "pixel_height": map_image.height,
        },
        "b": {
            "file": INSPECTOR_FILE.name,
            "rect_mm": (
                111.5,
                inspector_y,
                111.5 + inspector_width_mm,
                inspector_y + inspector_height_mm,
            ),
            "pixel_width": inspector_image.width,
            "pixel_height": inspector_image.height,
        },
        "c": {
            "file": NEIGHBOR_FILE.name,
            "rect_mm": (
                5.5,
                63.0,
                5.5 + neighbor_size_mm,
                63.0 + neighbor_size_mm,
            ),
            "pixel_width": neighbor_image.width,
            "pixel_height": neighbor_image.height,
        },
        "d": {
            "file": RECORDS_TABLE_FILE.name,
            "rect_mm": (5.5, 2.0, 5.5 + records_width_mm, 2.0 + records_height_mm),
            "pixel_width": records_image.width,
            "pixel_height": records_image.height,
        },
    }
    for item in layout.values():
        x0, y0, x1, y1 = item["rect_mm"]
        item["rendered_width_mm"] = x1 - x0
        item["rendered_height_mm"] = y1 - y0
        item["effective_dpi"] = item["pixel_width"] / (item["rendered_width_mm"] / 25.4)

    fig = c.mm_figure(WIDTH_MM, HEIGHT_MM)
    # b、d 截图自带完整边界；a、c 不再叠加会穿入文字的外框。
    for panel in ("b", "d"):
        _border_axis(fig, layout[panel]["rect_mm"])

    c.panel_letter(fig, "a", 1.5, 219.0, (WIDTH_MM, HEIGHT_MM))
    c.panel_letter(fig, "b", 107.0, 219.0, (WIDTH_MM, HEIGHT_MM))
    c.panel_letter(fig, "c", 1.5, 132.0, (WIDTH_MM, HEIGHT_MM))
    c.panel_letter(fig, "d", 1.5, 45.0, (WIDTH_MM, HEIGHT_MM))
    c_title = c.mark_panel_title(
        fig.text(
            9.0 / WIDTH_MM,
            125.5 / HEIGHT_MM,
            "CHEMBL18 and its neighbors",
            fontsize=6.5,
            ha="left",
            va="bottom",
            color=c.INK,
        )
    )
    c_title._qa_fig6_c_title = True  # type: ignore[attr-defined]
    c.mark_panel_title(
        fig.text(
            5.5 / WIDTH_MM,
            38.7 / HEIGHT_MM,
            "Source activity records",
            fontsize=7,
            ha="left",
            va="bottom",
            color=c.INK,
        )
    )
    _draw_map_legend(fig)
    _draw_neighbor_legend(fig)

    _write_plotted_values(layout)
    checks = [
        {
            "name": "capture_report status and hashes",
            "expected": "passed and four matching SHA-256 values",
            "actual": {"status": report["status"], "hashes": hashes},
            "passed": report["status"] == "passed" and len(hashes) == 4,
        },
        {
            "name": "panel a uses only map_viewport_in_clip_css",
            "expected": report["map_viewport_in_clip_css"],
            "actual": {
                "pixel_width": map_image.width,
                "pixel_height": map_image.height,
            },
            "passed": map_image.width
            == round(float(report["map_viewport_in_clip_css"]["width"]) * scale)
            and map_image.height
            == round(float(report["map_viewport_in_clip_css"]["height"]) * scale),
        },
        {
            "name": "panel effective dpi",
            "expected": ">=300",
            "actual": {
                name: round(float(item["effective_dpi"]), 2)
                for name, item in layout.items()
            },
            "passed": all(float(item["effective_dpi"]) >= 300 for item in layout.values()),
        },
        {
            "name": "inspector crop uses heading_top_css_px - 8 through evidence_bottom_css_px + 8",
            "expected": {
                "top_px": inspector_top_px,
                "bottom_px": inspector_bottom_px,
                "height_px": inspector_bottom_px - inspector_top_px,
            },
            "actual": {
                "top_px": inspector_top_px,
                "bottom_px": inspector_bottom_px,
                "height_px": inspector_image.height,
            },
            "passed": inspector_image.height
            == inspector_bottom_px - inspector_top_px,
        },
        {
            "name": "four panels a-d",
            "expected": ["a", "b", "c", "d"],
            "actual": sorted(layout, key=lambda value: "abcd".index(value)),
            "passed": set(layout) == {"a", "b", "c", "d"},
        },
    ]
    css_widths = {
        "b": inspector_image.width / scale,
        "d": float(report["records_table_clip_css"]["width"]),
    }
    effective_font_pt = {
        "b": 11.0 * layout["b"]["rendered_width_mm"] / css_widths["b"] * 72 / 25.4,
        "c": 6.5,
        "d": 11.0 * layout["d"]["rendered_width_mm"] / css_widths["d"] * 72 / 25.4,
    }
    checks.extend(
        [
            {
                "name": "panel b uses all 300 CSS px at >=58 mm",
                "expected": {"css_width": 300.0, "rendered_width_mm": ">=58"},
                "actual": {
                    "css_width": round(float(css_widths["b"]), 6),
                    "rendered_width_mm": layout["b"]["rendered_width_mm"],
                },
                "passed": math.isclose(float(css_widths["b"]), 300.0)
                and float(layout["b"]["rendered_width_mm"]) >= 58.0,
            },
            {
                "name": "panel d uses 6d_records_table.png at least 117 mm",
                "expected": {
                    "file": RECORDS_TABLE_FILE.name,
                    "rendered_width_mm": ">=117",
                },
                "actual": {
                    "file": layout["d"]["file"],
                    "rendered_width_mm": layout["d"]["rendered_width_mm"],
                },
                "passed": layout["d"]["file"] == RECORDS_TABLE_FILE.name
                and float(layout["d"]["rendered_width_mm"]) >= 117.0,
            },
            {
                "name": "neighbor cell letters S/G/F and selected Q",
                "expected": {"S": 5, "G": 7, "F": 6, "Q": 1},
                "actual": {
                    letter: sum(item["letter"] == letter for item in neighbor_labels)
                    for letter in ("S", "G", "F", "Q")
                },
                "passed": len(neighbor_labels) == 19,
            },
            {
                "name": "screenshot UI effective font size",
                "expected": "b, c and d >= 6 pt",
                "actual": {
                    panel: round(float(value), 6)
                    for panel, value in effective_font_pt.items()
                },
                "passed": all(value >= 6.0 - 1e-9 for value in effective_font_pt.values()),
            },
        ]
    )
    checks.extend(_layout_checks(fig, layout))
    audit = {
        "passed": all(item["passed"] for item in checks),
        "checks": checks,
        "source_size_mm": [WIDTH_MM, HEIGHT_MM],
        "figure_size_mm": [WIDTH_MM, HEIGHT_MM],
        "layout_mm": {
            panel: list(item["rect_mm"])
            for panel, item in layout.items()
        },
        "panels": layout,
        "screenshot_hashes": hashes,
        "capture_report": report,
        "map_legend": "software legend removed; matplotlib recreation below panel a",
        "inspector_crop_height_px": inspector_image.height,
        "inspector_crop_px": {
            "top": inspector_top_px,
            "bottom": inspector_bottom_px,
            "height": inspector_image.height,
        },
        "inspector_crop_planned_css_px": {
            "top": float(report["heading_top_css_px"]) - 8.0,
            "bottom": float(report["evidence_bottom_css_px"]) + 8.0,
        },
        "prepared_files": {
            "a": str(a_prepared),
            "b": str(b_prepared),
            "c": str(c_prepared),
        },
        "neighbor_labels": neighbor_labels,
        "effective_font_pt": effective_font_pt,
    }
    return fig, audit
