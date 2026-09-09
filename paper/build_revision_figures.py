"""Compose real-data manuscript panels as native SVG, PDF, and 600-dpi PNG."""
from __future__ import annotations

import argparse
import html
import io
import json
import os
from pathlib import Path
import xml.etree.ElementTree as ET

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".cache/matplotlib"))
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymupdf

from chemgridmap.plotting import DEFAULT_COLORS, _svg_molecule_content
from figure_geometry import window_bounds

ROOT = Path(__file__).resolve().parents[1]
NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", NS)
ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
COLORS = DEFAULT_COLORS
INK = "#172536"
BLUE = "#2766A5"
TEAL = "#188978"
ORANGE = "#CA7C23"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.labelsize": 9, "axes.titlesize": 10,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "svg.fonttype": "path", "pdf.fonttype": 42,
                     "axes.edgecolor": "#56616B", "axes.linewidth": 0.6})


def text(x, y, value, size=17, bold=False, color=INK, anchor="start"):
    return (f'<text x="{x}" y="{y}" font-family="Arial,Helvetica,sans-serif" font-size="{size}" '
            f'font-weight="{700 if bold else 400}" fill="{color}" text-anchor="{anchor}">'
            f'{html.escape(str(value))}</text>')


def rect(x, y, w, h, color="white", stroke="none", opacity=1, sw=1):
    return (f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{color}" '
            f'fill-opacity="{opacity}" stroke="{stroke}" stroke-width="{sw}"/>')


def document(content, width=1000, height=650):
    inches = 7.0
    return (f'<svg xmlns="{NS}" width="{inches}in" height="{height / width * inches}in" '
            f'viewBox="0 0 {width} {height}">' + rect(0, 0, width, height) + content + '</svg>')


def export(svg, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.with_suffix(".svg").write_text(svg, encoding="utf-8")
    source = pymupdf.open(stream=svg.encode(), filetype="svg")
    pdf = pymupdf.open(stream=source.convert_to_pdf(), filetype="pdf")
    pdf.save(path.with_suffix(".pdf"))
    pdf[0].get_pixmap(dpi=600, alpha=False).save(path.with_suffix(".png"))
    # Molecules must remain vector paths, not embedded bitmap surrogates.
    assert not ET.fromstring(svg).findall(f".//{{{NS}}}image")
    assert not pdf[0].get_images()


def embed_svg(svg, x, y, w, h):
    root = ET.fromstring(svg)
    _, _, vw, vh = [float(v) for v in root.get("viewBox").split()]
    scale = min(w/vw, h/vh)
    group = ET.Element(f"{{{NS}}}g", {"transform": f"translate({x+(w-vw*scale)/2},{y+(h-vh*scale)/2}) scale({scale})"})
    group.extend(list(root))
    return ET.tostring(group, encoding="unicode")


def chart_svg(fig):
    stream = io.StringIO()
    fig.savefig(stream, format="svg")
    plt.close(fig)
    return stream.getvalue()


def mol_tile(row, x, y, size, with_id=False, mark_conflict=True):
    color = COLORS[str(row.activity_class)]
    result = rect(x, y, size, size, color, color, 0.12, 0.8)
    draw_height = size - 16 if with_id else size*0.92
    content = _svg_molecule_content(row.canonical_smiles, 220)
    if content:
        s = draw_height / 220
        result += f'<g transform="translate({x + (size - draw_height) / 2},{y + size*0.035}) scale({s})">{content}</g>'
    if with_id:
        identity = str(row.get("parent_molecule_id", row.get("molecule_id", ""))).split(";")[0]
        result += text(x + size / 2, y + size - 4, identity, 10, anchor="middle")
    if mark_conflict and row.get("is_conflicted", 0):
        result += f'<path d="M{x+size-9},{y+1} L{x+size-1},{y+1} L{x+size-1},{y+9} Z" fill="{INK}"/>'
    return result


def map_panel(data, kind="map"):
    side = int(max(data.grid_row.max(), data.grid_col.max())) + 1
    out = ""
    for row in data.itertuples():
        px, py = row.projection_x * 300, (1 - row.projection_y) * 300
        gx, gy = row.grid_x * 300, (1 - row.grid_y) * 300
        color = COLORS[row.activity_class]
        if kind == "projection":
            out += f'<circle cx="{px}" cy="{py}" r="1.9" fill="{color}"/>'
        elif kind == "arrows":
            dx, dy = gx-px, gy-py
            length = max(np.hypot(dx, dy), 1e-9)
            ux, uy = dx/length, dy/length
            ax, ay = gx - 2.5*ux, gy - 2.5*uy
            out += f'<path d="M{px},{py} L{gx},{gy} M{ax-uy},{ay+ux} L{gx},{gy} L{ax+uy},{ay-ux}" stroke="{BLUE}" stroke-width="0.48" opacity="0.6" fill="none"/>'
            out += f'<circle cx="{px}" cy="{py}" r="0.55" fill="{INK}"/>'
        else:
            size = 300/(side-1)*0.91
            if kind == "overview":
                out += rect(gx-size/2, gy-size/2, size, size, color, opacity=0.72)
                if row.is_conflicted:
                    out += f'<circle cx="{gx}" cy="{gy}" r="2" fill="none" stroke="{INK}" stroke-width="0.8"/>'
            else:
                out += mol_tile(pd.Series(row._asdict()), gx-size/2, gy-size/2, size, mark_conflict=False)
    return out


def legend(x, y):
    result = ""
    for offset, label, desc in [(0, "inactive", "Inactive <=6"), (205, "medium", "Medium 6-7"), (420, "active", "Active >=7")]:
        result += rect(x+offset, y-12, 16, 12, COLORS[label], opacity=0.72)
        result += text(x+offset+24, y, desc, 15)
    return result


def select_windows(data, outdir):
    index = {(int(r.grid_col), int(r.grid_row)): i for i, r in data.iterrows()}
    candidates = []
    for col, row in sorted(index):
        cells = [(col+dx, row+dy) for dx in range(3) for dy in range(3)]
        if not all(cell in index for cell in cells):
            continue
        ids = [index[cell] for cell in cells]
        members = data.iloc[ids]
        counts = members.activity_class.value_counts()
        p = counts.values/9
        candidates.append({"col": col, "row": row, "active": int(counts.get("active", 0)),
                           "inactive": int(counts.get("inactive", 0)), "entropy": float(-(p*np.log(p)).sum()),
                           "mean_pchembl": float(members.activity_pchembl.mean()), "indices": ids})
    table = pd.DataFrame(candidates)
    if table.empty:
        raise ValueError("This map has no completely occupied 3x3 windows.")
    selected = []
    used = set()
    orders = [("Active-rich", ["active", "mean_pchembl"], [False, False]),
              ("Mixed", ["entropy", "col", "row"], [False, True, True]),
              ("Inactive-rich", ["inactive", "mean_pchembl"], [False, True])]
    for name, by, ascending in orders:
        for _, candidate in table.sort_values(by, ascending=ascending, kind="stable").iterrows():
            if not used.intersection(candidate["indices"]):
                selected.append({"name": name, **candidate.to_dict()})
                used.update(candidate["indices"])
                break
    table.to_csv(outdir / "window_candidates.csv", index=False)
    all_members = []
    for win in selected:
        frame = data.iloc[win["indices"]].copy()
        frame["window"] = win["name"]
        all_members.append(frame)
    pd.concat(all_members).to_csv(outdir / "window_members.csv", index=False)
    (outdir / "selected_windows.json").write_text(json.dumps(selected, indent=2))
    return selected


def window_svg(data, win, size=260, with_ids=True):
    step = size/3
    result = ""
    for _, row in data.iloc[win["indices"]].iterrows():
        x = (row.grid_col-win["col"])*step
        y = (win["row"]+2-row.grid_row)*step
        result += mol_tile(row, x, y, step-3, with_id=with_ids)
    return result


def figure1(data, out):
    result = ""
    for x, letter, title, kind in [(25, "A", "Molecular projection", "projection"),
                                   (350, "B", "Projection-to-grid vectors", "arrows"),
                                   (675, "C", "Molecule-resolved grid", "map")]:
        result += text(x, 28, f"{letter}  {title}", 18, True)
        panel = map_panel(data, kind)
        result += f'<g transform="translate({x},75)">{panel}</g>'
        export(document(panel, 300, 300), out / "components" / f"Figure_1{letter}")
    result += legend(230, 420)
    result += text(25, 468, "D  Auditable CSV-to-map workflow", 18, True)
    stages = [("ChEMBL CSV", ["Target / endpoint", "Record identifiers"]),
              ("Record audit", ["Exclusion reasons", "Parent structures"]),
              ("Molecule table", ["Median / range", "Conflict + source links"]),
              ("Projection + grid", ["Isotropic scaling", "Unique-cell assignment"]),
              ("Inspect + validate", ["Vector maps / CSV", "Metrics / run manifest"])]
    for i, (title, lines) in enumerate(stages):
        x = 25+i*195
        result += rect(x, 494, 174, 86, "#F5F7F9", "#BDC8D3", sw=1)
        result += text(x+87, 518, title, 16, True, anchor="middle")
        for j, line in enumerate(lines):
            result += text(x+87, 541+j*20, line, 14, anchor="middle")
        if i < 4:
            result += f'<path d="M{x+177},538 H{x+190} l-4,-3 m4,3 l-4,3" fill="none" stroke="{INK}" stroke-width="1.4"/>'
    result += text(500, 618, "889 molecules | 48 x 48 candidate lattice | 38.6% occupied", 15, anchor="middle")
    export(document(result, height=645), out / "Figure_1_workflow")


def figure2(data, windows, source, out):
    result = text(25, 28, "A  Map with audit markers", 18, True)
    result += f'<g transform="translate(30,65) scale(1.06)">{map_panel(data, "overview")}</g>'
    audit_row = data[data.molecule_identity_key.eq("CHEMBL20")].iloc[0]
    dx, dy = 30+audit_row.grid_x*318, 65+(1-audit_row.grid_y)*318
    result += f'<circle cx="{dx}" cy="{dy}" r="7" fill="none" stroke="{INK}" stroke-width="1.7"/>'
    result += text(dx+8, dy-6, "D", 16, True)
    for i, win in enumerate(windows[:2]):
        col = BLUE if i == 0 else ORANGE
        gx, gy, width, height = window_bounds(data, win)
        result += rect(gx, gy, width, height, "none", col, sw=2)
        result += text(gx, gy-5, "B" if i == 0 else "C", 16, True, col)
        x = 402+i*300
        result += text(x, 28, ("B  " if i == 0 else "C  ")+win["name"]+" window", 18, True)
        result += f'<g transform="translate({x},68)">{window_svg(data, win)}</g>'
        result += text(x, 352, f"9 cells; mean pChEMBL {win['mean_pchembl']:.2f}", 15)
        export(document(window_svg(data, win), 260, 260), out / "components" / f"window_{i+1}")
    result += text(30, 400, "Circles / corner marks: conflict flag", 14)
    result += text(402, 380, "Fill denotes the median-derived activity class.", 14)
    result += text(402, 400, "Selected windows are illustrative, not discovered biological domains.", 13)
    result += '<path d="M25,424 H975" stroke="#CCD4DC"/>'
    result += text(25, 454, "D  Cell-to-record audit: acetazolamide (CHEMBL20)", 19, True)
    row = data[data.molecule_identity_key.eq("CHEMBL20")].iloc[0]
    result += mol_tile(row, 33, 492, 165)
    result += text(30, 683, "Color: active median", 14, True)
    result += text(30, 705, "Marker: inspect the records", 13)
    result += text(245, 480, "54 records / 54 assays / 27 documents; median 7.185; range 4.73-8.47", 14)
    if (source / "acetazolamide_condition_series.csv").exists():
        result += condition_plot(source)
    else:
        result += description_plot(source)
    result += rect(25, 741, 950, 49, "#F0F5F7", "#D2DDE3")
    result += text(45, 762, "Retrieve by molecule ID or grid cell", 15, True)
    result += text(45, 781, "chemgridmap inspect coordinates.csv --records retained_records.csv --molecule-id CHEMBL20", 13)
    result += text(500, 817, "Source-linked records support assay review; a median color alone does not establish measurement agreement.", 14, anchor="middle")
    export(document(result, height=840), out / "Figure_2_record_audit")


def condition_plot(source):
    records = pd.read_csv(source / "acetazolamide_condition_series.csv")
    fig, ax = plt.subplots(figsize=(6.3, 2.2))
    for label, color, marker in [("DTT", BLUE, "o"), ("2-ME", ORANGE, "s"),
                                  ("TCEP", TEAL, "^"), ("Threitol", "#806755", "D")]:
        group = records[records.additive.eq(label)].sort_values("concentration_mM")
        ax.plot(group.concentration_mM, group.pchembl_value, color=color, marker=marker,
                ms=3.5, lw=0.8, label=label)
    ax.axhline(8.05, color=INK, ls="--", lw=0.7)
    ax.text(0.011, 8.10, "Reference record: additive not specified", fontsize=7.5)
    ax.set(xscale="log", xlabel="Additive concentration (mM; log scale)", ylabel="pChEMBL",
           xlim=(0.008, 13), ylim=(6.2, 8.35))
    ax.legend(frameon=False, ncol=4, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.20))
    ax.grid(axis="y", color="#E2E7EC", linewidth=0.6)
    ax.tick_params(labelsize=8)
    fig.subplots_adjust(left=0.10, right=0.99, bottom=0.23, top=0.82)
    result = embed_svg(chart_svg(fig), 235, 500, 740, 223)
    result += text(605, 730, "25 records from one study: four additive series + one reference record", 13, anchor="middle")
    return result


def description_plot(source):
    records = pd.read_csv(source / "chembl205/acetazolamide_retained_records.csv")
    fig, ax = plt.subplots(figsize=(6.3, 2.2))
    labels = ["CO2 / stopped-flow mentioned", "Nitrophenyl substrate mentioned", "Other / insufficient description"]
    for i, name in enumerate(labels):
        group = records[records.description_group.eq(name)].sort_values("activity_id")
        offsets = np.linspace(-0.18, 0.18, len(group))
        yy = i+offsets
        xx = group.pchembl_value.to_numpy()
        cc = [COLORS["inactive" if v <= 6 else "active" if v >= 7 else "medium"] for v in xx]
        ax.scatter(xx, yy, c=cc, s=14, edgecolor="white", linewidth=0.25, zorder=3)
    ax.axvline(7.185, color=INK, ls="--", lw=0.8)
    ax.set_yticks(range(3), ["CO2 / stopped-flow\nmentioned (n=7)", "Nitrophenyl substrate\nmentioned (n=11)", "Other / insufficient\ndescription (n=36)"])
    ax.set(xlabel="pChEMBL of retained records", xlim=(4.4, 8.8), ylim=(2.5, -0.5))
    ax.tick_params(labelsize=8)
    ax.grid(axis="x", color="#E2E7EC", linewidth=0.6)
    fig.subplots_adjust(left=0.24, right=0.98, bottom=0.24, top=0.98)
    return embed_svg(chart_svg(fig), 220, 500, 760, 220)


def figure3(source, out):
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.2), layout="constrained")
    occupancy = pd.read_csv(source / "occupancy.csv")
    ax = axes[0, 0]
    for key, label, color, marker in [("projection_to_grid_trustworthiness", "Trustworthiness", BLUE, "o"),
                                      ("projection_grid_tie_aware_recall", "Tie-aware recall", TEAL, "s")]:
        ax.plot(occupancy.actual_occupancy*100, occupancy[key], marker=marker, ms=3.5, color=color, label=label)
    ax.axvline(38.585, color="#8995A2", lw=0.7, ls=":")
    ax.set(title="A  Occupancy sensitivity", xlabel="Actual grid occupancy (%)", ylabel="Projection-to-grid fidelity",
           ylim=(0, 1.03), xlim=(92, 15))
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    baseline = pd.read_csv(source / "baselines.csv")
    targets = ["CHEMBL205", "CHEMBL204", "CHEMBL240"]
    ax = axes[0, 1]
    for offset, method, color, marker in [(-0.13, "Greedy", "#87949E", "o"), (0.13, "Minimum-cost", TEAL, "s")]:
        means = []
        for i, target in enumerate(targets):
            values = baseline.query("target == @target and method == @method").groupby("seed").projection_grid_tie_aware_recall.mean()
            ax.scatter(i+offset+np.linspace(-0.04, 0.04, len(values)), values, c=color, s=17, marker=marker)
            means.append(values.mean())
        ax.plot(np.arange(3)+offset, means, color=color, marker=marker, lw=0, label=method)
    ax.set_xticks(range(3), ["CA II", "D3", "hERG"])
    ax.set(title="B  Assignment controls", ylabel="Tie-aware recall", ylim=(0, 1))
    ax.legend(frameon=False, fontsize=8)
    stability = pd.read_csv(source / "seed_stability.csv")
    ax = axes[1, 0]
    for offset, column, label, color, marker in [(-0.13, "projection_neighbor_jaccard", "Projection", BLUE, "o"),
                                               (0.13, "grid_neighbor_jaccard", "Grid", ORANGE, "s")]:
        for i, target in enumerate(targets):
            values = stability[stability.target.eq(target)][column]
            ax.scatter(i+offset+np.linspace(-0.045, 0.045, len(values)), values, s=12, c=color, marker=marker,
                       label=label if i == 0 else None)
            ax.plot([i+offset-0.07, i+offset+0.07], [values.mean()]*2, color=INK, lw=1)
    ax.set_xticks(range(3), ["CA II", "D3", "hERG"])
    ax.set(title="C  Seed-dependent neighborhoods", ylabel="Neighbor-set Jaccard", ylim=(0, 1))
    ax.legend(frameon=False, fontsize=8)
    fixed = pd.read_csv(source / "sparse_calibration.csv")
    adaptive = pd.read_csv(source / "adaptive_validation.csv")
    ax = axes[1, 1]
    for offset, table, label, color, marker in [(-0.13, fixed[fixed.nearest_candidates.eq(128)], "Fixed 128", "#87949E", "o"),
                                               (0.13, adaptive, "Adaptive", TEAL, "s")]:
        for i, target in enumerate(targets):
            values = table[table.target.eq(target)].dense_cost_ratio
            ax.scatter(i+offset+np.linspace(-0.04, 0.04, len(values)), values, s=17, c=color, marker=marker,
                       label=label if i == 0 else None)
    ax.axhline(1, color=INK, ls="--", lw=0.7)
    ax.set_xticks(range(3), ["CA II", "D3", "hERG"])
    ax.set(title="D  Sparse matching calibration", ylabel="Cost / dense optimum", ylim=(0.85, max(2.6, fixed.query('nearest_candidates == 128').dense_cost_ratio.max()*1.05)))
    ax.legend(frameon=False, fontsize=8)
    for ax in axes.flat:
        ax.grid(axis="y", color="#E2E7EC", linewidth=0.6)
        title = ax.get_title()
        ax.set_title("")
        ax.set_title(title, loc="left", fontsize=9.5, weight="bold")
    export(chart_svg(fig), out / "Figure_3_validation")


def supplement(source, out, data, windows):
    panel = description_plot(source)
    export(document(text(25, 30, "All 54 retained acetazolamide records", 19, True)
                    + f'<g transform="translate(-185,-435)">{panel}</g>', 800, 310),
           out / "Figure_S2_all_records")
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.7), layout="constrained")
    for ax, filename, title in [(axes[0], "runtime.csv", "A  Dense assignment"),
                                (axes[1], "large_scale_runtime.csv", "B  Adaptive sparse assignment")]:
        table = pd.read_csv(source / filename)
        for n, group in table.groupby("n"):
            ax.scatter([n]*len(group), group.seconds, color=BLUE, s=10, alpha=0.6)
        summary = table.groupby("n").seconds.median()
        ax.plot(summary.index, summary, color=TEAL, marker="o", ms=4)
        ax.set(title=title, xlabel="Number of molecules", ylabel="Assignment time (s)", ylim=(0, None))
        ax.grid(axis="y", color="#E2E7EC", linewidth=0.6)
    export(chart_svg(fig), out / "Figure_S1_runtime")
    win = windows[2]
    result = text(20, 25, "Inactive-rich illustrative window", 18, True)
    result += f'<g transform="translate(20,55)">{window_svg(data, win, 360)}</g>'
    result += text(200, 447, f"9 molecules; mean pChEMBL {win['mean_pchembl']:.2f}", 16, anchor="middle")
    export(document(result, 400, 465), out / "Figure_S3_inactive_window")


def toc(data, win, out):
    # Every plotted point is retained. Labels occupy their own space, not white masks.
    result = text(30, 25, "ChEMBL records", 18, True)
    for i in range(5):
        result += rect(30, 48+i*20, 128, 19, "#EEF2F6", "#C6D0DA", sw=0.7)
        result += text(40, 61+i*20, f"Record {i+1} / assay", 11)
    result += text(94, 188, "Audit + aggregate", 14, anchor="middle")
    result += '<path d="M165,107 H200 l-8,-5 m8,5 l-8,5" stroke="#172536" fill="none" stroke-width="2"/>'
    result += text(225, 25, "Grid chemical map", 18, True)
    result += f'<g transform="translate(239,52) scale(0.48)">{map_panel(data, "overview")}</g>'
    result += '<path d="M395,107 H430 l-8,-5 m8,5 l-8,5" stroke="#172536" fill="none" stroke-width="2"/>'
    result += text(450, 25, "Inspectable evidence", 18, True)
    result += f'<g transform="translate(467,45)">{window_svg(data, win, 145, with_ids=False)}</g>'
    export(document(result, 665, 215), out / "TOC")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=ROOT / "paper/output/revision_v4")
    parser.add_argument("--output", type=Path, default=ROOT / "paper/output/revision_v4/figures")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(args.source / "chembl205/coordinates_seed_42.csv")
    windows = select_windows(data, args.source)
    figure1(data, args.output)
    figure2(data, windows, args.source, args.output)
    figure3(args.source, args.output)
    if (args.source / "large_scale_runtime.csv").exists():
        supplement(args.source, args.output, data, windows)
    toc(data, windows[1], args.output)
    print(f"Vector figures and 600-dpi PNGs: {args.output}")


if __name__ == "__main__":
    main()
