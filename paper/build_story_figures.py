"""Task-led manuscript figures using frozen results and a real interface crop."""
import argparse
import base64
import hashlib
import io
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymupdf
from PIL import Image

import build_submission_figures as base
from build_revision_figures import document, text, rect, embed_svg, chart_svg
from chemgridmap.neighborhoods import neighborhood_evidence
from chemgridmap.representations import morgan_fingerprints

SOURCE, ARCHIVE = base.SOURCE, base.ARCHIVE
INK, BLUE, ORANGE, TEAL = base.INK, base.BLUE, base.ORANGE, base.TEAL
OUT = None
EXTRA_CHECKS = {}


def export(svg, name):
    path = OUT / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.with_suffix('.svg').write_text(svg)
    source = pymupdf.open(stream=svg.encode(), filetype='svg')
    pdf = pymupdf.open(stream=source.convert_to_pdf(), filetype='pdf')
    pdf.save(path.with_suffix('.pdf'))
    pdf[0].get_pixmap(dpi=600, alpha=False).save(path.with_suffix('.png'))


def arrow(x1, y1, x2, y2, color=base.MUTED):
    ux, uy = np.array([x2-x1, y2-y1]) / np.hypot(x2-x1, y2-y1)
    bx, by = x2-ux*6, y2-uy*6
    return (f'<path d="M{x1} {y1} L{x2} {y2}" stroke="{color}" stroke-width="1.3" fill="none"/>'
            f'<path d="M{x2} {y2} L{bx-uy*2.7} {by+ux*2.7} L{bx+uy*2.7} {by-ux*2.7} Z" fill="{color}"/>')


def figure1(data, groups, case):
    main = case['group']
    shown = [main] + [g for g in groups if g['structure_group'] != main['structure_group']][:2]
    body = '<defs><marker id="flow-arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="5" markerHeight="5" orient="auto"><path d="M0 0 L10 5 L0 10" fill="#536575"/></marker></defs>'
    body += base.panel(20, 28, 'A', 'Select a structural group')
    body += '<g transform="translate(54,54) scale(.80)">' + base.previous.group_map(data, {g['structure_group']: g['color'] for g in shown}) + '</g>'
    # External labels point to actual member cells rather than inferred regions.
    for g, (tx, ty) in zip(shown, [(343, 163), (304, 298), (29, 155)]):
        rows = data[data.structure_group.eq(g['structure_group'])]
        center = rows[['grid_x', 'grid_y']].median().to_numpy()
        near = rows.iloc[np.argmin(np.square(rows[['grid_x', 'grid_y']].to_numpy()-center).sum(axis=1))]
        px, py = 54 + near.grid_x*256, 54 + (1-near.grid_y)*256
        body += arrow(tx, ty+6, px, py, g['color'])
        body += rect(tx-23, ty-15, 49, 20, 'white')
        body += text(tx, ty, g['structure_group'], 13, True, g['color'], 'middle')
    body += text(25, 336, '889 molecules  |  314 structure-defined groups', 12, True)
    body += text(25, 355, 'Gray: other groups. Group IDs remain fixed across layouts.', 11, color=base.MUTED)
    body += text(25, 382, 'Actual G002 comparison window', 12, True)
    screenshot = SOURCE / 'browser/desktop_comparison.png'
    crop_box = (260, 98, 1178, 401)
    crop = Image.open(screenshot).crop(crop_box)
    stream = io.BytesIO(); crop.save(stream, format='PNG')
    uri = base64.b64encode(stream.getvalue()).decode()
    body += f'<image x="25" y="392" width="350" height="116" href="data:image/png;base64,{uri}"/>'
    body += rect(25, 392, 350, 116, 'none', base.LINE)
    EXTRA_CHECKS['interface_crop'] = {'source': str(screenshot), 'box': crop_box,
                                     'processing': 'Crop only; no content replacement.',
                                     'sha256': hashlib.sha256(screenshot.read_bytes()).hexdigest()}
    body += base.panel(410, 28, 'B', 'Check shared-core coverage')
    for i, g in enumerate(shown):
        y = 50 + i*127
        body += rect(410, y, 570, 115, '#F8FAFB', base.LINE) + rect(410, y, 4, 115, g['color'])
        body += base.previous.molecule(g['core_smiles'], 425, y+6, 136, 101)
        body += text(579, y+27, f"{g['structure_group']}  |  {g['size']} molecules", 17, True, g['color'])
        body += text(579, y+53, f"Core coverage  {g['core_coverage_count']}/{g['size']}", 15, True)
        body += text(801, y+53, f"Exact scaffold  {g['core_exact_scaffold_count']}/{g['size']}", 13, color=base.MUTED)
        body += text(579, y+78, f"Representative  {g['representative_id']}", 13)
        status = 'Majority coverage' if g['core_coverage_count'] > g['size']/2 else 'No majority coverage'
        body += text(579, y+101, status, 12, color=base.MUTED)
    body += arrow(694, 432, 694, 486)
    body += text(713, 458, 'Open G002 members', 14, True)
    body += text(713, 479, 'Compare at the same core orientation', 12, color=base.MUTED)
    body += base.panel(20, 544, 'C', 'Compare peripheral structures')
    body += text(980, 544, 'Orange: outside the core  |  strip: activity', 12, color=base.MUTED, anchor='end')
    members = pd.read_csv(SOURCE/'chembl205/case_members.csv').head(3)
    body += embed_svg(base.cards(members, main['core_smiles'], 'Figure_1'), 20, 563, 960, 224)
    body += text(500, 814, 'Select a group  >  inspect its core  >  compare members  >  follow the evidence', 14, True, anchor='middle')
    export(document(body, height=832), 'Figure_1_structural_comparison')


def figure2(data, case):
    captured = []
    original = base.export
    base.export = lambda svg, path: captured.append(svg)
    try:
        base.figure2(data, case)
    finally:
        base.export = original
    root = ET.fromstring(captured[0])
    changes = {
        'A structural exception within G002': 'Check the structural summary',
        'Selected benzothiazole core: 27 of 28 members': 'Task 1  |  G002: does the selected core cover every member?',
        'An annotation to inspect': 'Check an activity annotation',
        'Separate audit case, outside G002': 'Task 2  |  separate molecule, outside G002',
        'Assay conditions explain part of the variation': 'Follow the annotation to its experimental context',
        'The unmatched member remains visible.': 'Coverage flags the exception for member-level review.',
        '25 IC50 values match source tables': '25 records traced to one publication',
        'Median colors summarize observations;': 'Use the records to interpret the color;',
        'they do not harmonize assays.': 'compare potency within matched conditions.'}
    for node in root.iter(base.NS+'text'):
        if node.text in changes:
            node.text = changes[node.text]
            if node.text.startswith('Task 1'): node.set('font-size', '12')
            if node.text == 'Check an activity annotation': node.set('font-size', '16')
    export(ET.tostring(root, encoding='unicode'), 'Figure_2_structure_and_evidence')


def neighbor_panel(data, case):
    center_id = case['group']['representative_id']
    idx = int(np.flatnonzero(data.comparison_id.eq(center_id))[0])
    current, summary = neighborhood_evidence(data, morgan_fingerprints(data.canonical_smiles.tolist()), idx)
    archived = pd.read_csv(SOURCE/'chembl205/case_center_neighbors.csv')
    left = current.set_index('molecule_id').sort_index()
    right = archived.set_index('molecule_id').sort_index()
    assert left.neighbor_relation.equals(right.neighbor_relation)
    assert np.allclose(left.representation_distance, right.representation_distance)
    assert np.allclose(left.grid_distance, right.grid_distance)
    EXTRA_CHECKS['neighbor_example'] = {'center': center_id, **summary, 'archived_sets_reproduced': True}
    current.to_csv(OUT/'neighbor_example.csv', index=False)
    body = base.panel(516, 355, 'D', 'Inspect apparent neighbors')
    body += text(516, 382, f'Query: {center_id}  |  k = 10, including ties', 12, color=base.MUTED)
    coordinates = data.set_index('comparison_id')
    query = coordinates.loc[center_id]
    palette = {'shared': BLUE, 'grid_only': ORANGE, 'representation_only': TEAL}
    member_relation = dict(zip(current.molecule_id, current.neighbor_relation))
    selected = coordinates.loc[list(member_relation)+[center_id]]
    lo_c, hi_c = int(selected.grid_col.min())-1, int(selected.grid_col.max())+1
    lo_r, hi_r = int(selected.grid_row.min())-1, int(selected.grid_row.max())+1
    step = min(214/(hi_c-lo_c+1), 211/(hi_r-lo_r+1))
    for identity, row in coordinates.iterrows():
        if not(lo_c <= row.grid_col <= hi_c and lo_r <= row.grid_row <= hi_r): continue
        x, y = 520+(row.grid_col-lo_c)*step, 402+(hi_r-row.grid_row)*step
        relation = member_relation.get(identity)
        color = palette.get(relation, '#E7ECEF')
        body += rect(x, y, step*.84, step*.84, color, 'white', sw=.5)
        if identity == center_id:
            body += rect(x, y, step*.84, step*.84, INK)
            body += text(x+step*.42, y+step*.61, 'Q', step*.58, True, 'white', 'middle')
        elif relation:
            mark = {'shared':'S', 'grid_only':'G', 'representation_only':'F'}[relation]
            body += text(x+step*.42, y+step*.61, mark, step*.48, True, 'white', 'middle')
    counts = current.neighbor_relation.value_counts()
    for i, (key, label, mark) in enumerate([('shared','Shared neighbors','S'), ('grid_only','Grid only','G'), ('representation_only','Fingerprint only','F')]):
        y = 427+i*57
        body += rect(760, y-15, 22, 22, palette[key])
        body += text(771, y+1, mark, 13, True, 'white', 'middle')
        body += text(796, y+1, f'{int(counts[key])}  {label}', 14, True)
    body += text(759, 599, '12 grid / 11 fingerprint neighbors', 11, color=base.MUTED)
    body += text(516, 651, 'Open a marked cell to inspect its structure and identity.', 12, True)
    return body


def figure3(data, case):
    # Keep the full previous quantitative figure in SI, including every group.
    original = base.export
    base.export = lambda svg, path: export(svg, 'supplement/Full_layout_diagnostics')
    try:
        base.figure3()
    finally:
        base.export = original
    old = ET.parse(OUT/'supplement/Full_layout_diagnostics.svg').getroot()
    # Retain original A-C paths; relocate the entire D panel to SI without rescaling its data.
    for parent in old.iter():
        for child in list(parent):
            if child.get('id') == 'axes_4': parent.remove(child)
    inner = ''.join(ET.tostring(n, encoding='unicode') for n in old)
    vw, vh = map(float, old.get('viewBox').split()[2:])
    scale = 1000/vw
    body = f'<g transform="scale({scale})">{inner}</g>'
    body += neighbor_panel(data, case)
    export(document(body, height=681), 'Figure_3_structure_fidelity')


def main():
    global OUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    OUT = parser.parse_args().output_dir
    OUT.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(SOURCE/'chembl205/coordinates_seed_42.csv')
    groups = json.loads((SOURCE/'chembl205/map_structure_groups.json').read_text())['groups']
    case = json.loads((SOURCE/'chembl205/case.json').read_text())
    figure1(data, groups, case); figure2(data, case); figure3(data, case)
    inputs = [SOURCE/p for p in ['chembl205/coordinates_seed_42.csv', 'chembl205/case_members.csv',
              'chembl205/map_structure_groups.json', 'chembl205/case.json', 'chembl205/case_center_neighbors.csv',
              'group_macro_metrics.csv', 'group_layout_metrics.csv']]
    inputs += [ARCHIVE/p for p in ['acetazolamide_condition_series.csv', 'full_path_fidelity.csv', 'seed_stability.csv']]
    report = {'scope': 'Task-led presentation revision; no layout, membership, activity or metric changes.',
              'input_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
              'drawing_checks': base.CHECKS, **EXTRA_CHECKS}
    (OUT/'figure_provenance.json').write_text(json.dumps(report, indent=2))
    print(OUT)


if __name__ == '__main__':
    main()
