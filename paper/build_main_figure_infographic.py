"""A task-led Figure 1 using the Infographic Engine layout and frozen chemistry."""
from pathlib import Path
import argparse
import hashlib
import json
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import pymupdf
from PIL import Image
from rdkit import Chem

import build_submission_figures as base
from build_revision_figures import text, rect, embed_svg
from build_story_figures import arrow

SOURCE = base.SOURCE / 'chembl205'
NS = '{http://www.w3.org/2000/svg}'
INK, MUTED, LINE = base.INK, base.MUTED, base.LINE


def heading(x, y, letter, title):
    return text(x, y, letter, 23, True) + text(x+31, y, title, 21, True)


def line(x1, y1, x2, y2, color=LINE):
    return f'<path d="M{x1} {y1} L{x2} {y2}" stroke="{color}" stroke-width="1" fill="none"/>'


def overview(data, groups):
    colors = {group['structure_group']: group['color'] for group in groups}
    body = heading(20, 31, 'A', 'Select a group')
    # Cell centers retain the exact archived geometry; only the figure scale changes.
    side = int(max(data.grid_row.max(), data.grid_col.max())) + 1
    extent, ox, oy = 250, 49, 63
    step = extent / (side-1)
    for row in data.itertuples():
        x, y = ox + row.grid_x*extent, oy+(1-row.grid_y)*extent
        color = colors.get(row.structure_group, '#E2E7EA')
        body += f'<g data-cell="{row.comparison_id}" data-row="{row.grid_row}" data-col="{row.grid_col}">'
        body += rect(x-step*.44, y-step*.44, step*.88, step*.88, color)
        body += '</g>'
    callouts = []
    for group, (tx, ty) in zip(groups, [(351, 128), (342, 293), (27, 160)]):
        members = data[data.structure_group.eq(group['structure_group'])]
        center = members[['grid_x', 'grid_y']].median().to_numpy()
        row = members.iloc[np.argmin(np.square(members[['grid_x', 'grid_y']].to_numpy()-center).sum(axis=1))]
        px, py = ox+row.grid_x*extent, oy+(1-row.grid_y)*extent
        # Stop outside the text instead of hiding the leader beneath a white label.
        half_width = pymupdf.Font('hebo').text_length(group['structure_group'], fontsize=15)/2
        label_side = -1 if px < tx else 1
        ex, ey = tx+label_side*(half_width+6), ty-5.25
        assert ex < tx-half_width if label_side == -1 else ex > tx+half_width
        body += f'<g data-callout="{group["structure_group"]}">'
        body += line(px, py, ex, ey, group['color'])
        body += f'<circle cx="{px}" cy="{py}" r="2.5" fill="white" stroke="{group["color"]}" stroke-width="1.3"/>'
        body += text(tx, ty, group['structure_group'], 15, True, group['color'], 'middle')
        body += '</g>'
        callouts.append({'group': group['structure_group'], 'molecule_id': row.comparison_id,
                         'grid_row': int(row.grid_row), 'grid_col': int(row.grid_col), 'figure_center': [px, py],
                         'leader_end': [ex, ey], 'label_x_bounds': [tx-half_width, tx+half_width],
                         'label_gap': 6, 'label_background': False})
    body += text(20, 337, '889 molecules  |  314 structure groups', 15, True)
    body += line(20, 349, 392, 349)
    body += text(20, 370, 'Other groups', 14, True)
    body += text(233, 370, 'Core coverage', 14, anchor='middle')
    body += text(350, 370, 'Exact scaffold', 14, anchor='middle')
    for group, y in zip(groups[1:], [392, 414]):
        body += rect(21, y-10, 8, 8, group['color'])
        body += text(38, y, group['structure_group'], 14, True, group['color'])
        body += text(233, y, f"{group['core_coverage_count']}/{group['size']}", 14, anchor='middle')
        body += text(350, y, f"{group['core_exact_scaffold_count']}/{group['size']}", 14, anchor='middle')
    return body, callouts


def core_panel(group):
    color = group['color']
    body = heading(441, 31, 'B', 'Identify the shared core')
    body += text(459, 75, group['structure_group'], 22, True, color)
    body += text(537, 75, f"{group['size']} molecules", 18, True)
    body += base.previous.molecule(group['core_smiles'], 464, 90, 225, 181)
    body += text(576, 290, 'Benzothiazole core', 16, True, anchor='middle')
    body += text(725, 119, 'Core coverage', 17, True)
    body += text(725, 166, f"{group['core_coverage_count']}/{group['size']}", 38, True, color)
    body += text(725, 191, 'with the selected core', 14, color=MUTED)
    body += text(725, 231, 'Exact scaffold', 17, True)
    body += text(725, 278, f"{group['core_exact_scaffold_count']}/{group['size']}", 38, True)
    body += text(725, 303, 'with this Murcko scaffold', 14, color=MUTED)
    body += line(459, 323, 972, 323)
    body += text(459, 346, '1 unmatched member remains available for inspection', 14)
    body += arrow(683, 372, 683, 415, MUTED)
    body += text(705, 386, 'Compare members', 15, True)
    body += text(705, 407, 'Same core orientation and scale', 14, color=MUTED)
    return body


def comparison_panel(members, core):
    body = heading(20, 458, 'C', 'Compare molecular differences')
    body += text(978, 456, 'Orange: outside the core  |  strip: activity', 14, color=MUTED, anchor='end')
    svg = base.cards(members, core, 'Figure_1_infographic', height=246)
    node = ET.fromstring(svg)
    for element in node.iter(NS+'rect'):
        if element.get('rx'):
            element.set('stroke', 'none')
    for element in node.iter(NS+'text'):
        element.set('font-size', '17' if element.get('font-weight') == '600' else '14')
    body += embed_svg(ET.tostring(node, encoding='unicode'), 20, 477, 960, 246)
    for x in [340, 660]:
        body += line(x, 484, x, 720)
    return body


def export(body, output, name, bounds=(0, 0, 1000, 740)):
    x, y, width, height = bounds
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{7*width/1000}in" '
           f'height="{7*height/1000}in" viewBox="{x} {y} {width} {height}">'
           + rect(x, y, width, height) + body + '</svg>')
    target = output/name
    target.parent.mkdir(parents=True, exist_ok=True)
    target.with_suffix('.svg').write_text(svg)
    source = pymupdf.open(stream=svg.encode(), filetype='svg')
    pdf = pymupdf.open(stream=source.convert_to_pdf(), filetype='pdf')
    pdf.save(target.with_suffix('.pdf'))
    pdf[0].get_pixmap(dpi=600, alpha=False).save(target.with_suffix('.png'))
    return target.with_suffix('.svg')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    output = parser.parse_args().output_dir
    output.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(SOURCE/'coordinates_seed_42.csv')
    groups = json.loads((SOURCE/'map_structure_groups.json').read_text())['groups']
    main_group = json.loads((SOURCE/'case.json').read_text())['group']
    displayed = [main_group]+[g for g in groups if g['structure_group'] != main_group['structure_group']][:2]
    members = pd.read_csv(SOURCE/'case_members.csv').head(3)
    assert members.comparison_id.tolist() == ['CHEMBL81688', 'CHEMBL79340', 'CHEMBL309608']
    assert len(data) == 889 and len(groups) == 314
    checks = []
    for group in displayed:
        molecules = data[data.structure_group.eq(group['structure_group'])]
        core = Chem.MolFromSmiles(group['core_smiles'])
        coverage = sum(Chem.MolFromSmiles(s).HasSubstructMatch(core, useChirality=True) for s in molecules.canonical_smiles)
        exact = int(molecules.scaffold_smiles.eq(group['core_smiles']).sum())
        assert coverage == group['core_coverage_count'] and exact == group['core_exact_scaffold_count']
        checks.append({'group': group['structure_group'], 'size': len(molecules), 'coverage': coverage, 'exact': exact})
    a, callouts = overview(data, displayed)
    b = core_panel(main_group)
    c = comparison_panel(members, main_group['core_smiles'])
    body = a + line(416, 53, 416, 418) + b + line(20, 432, 980, 432) + c
    path = export(body, output, 'Figure_1_main_infographic')
    export(a, output, 'panels/A_group_overview', (0, 0, 410, 428))
    export(b, output, 'panels/B_shared_core', (432, 0, 568, 428))
    export(c, output, 'panels/C_molecular_comparison', (0, 438, 1000, 302))
    drawing_checks = base.CHECKS['Figure_1_infographic']
    assert all(check['structure_preserved'] and check['aligned'] for check in drawing_checks)
    coordinates = [np.asarray(check['core_draw_coordinates']) for check in drawing_checks]
    assert all(np.allclose(coordinates[0], value, atol=1e-5) for value in coordinates[1:])
    root = ET.parse(path).getroot()
    assert len([g for g in root.iter(NS+'g') if g.get('data-cell')]) == len(data)
    assert not list(root.iter(NS+'image')), 'Final scientific figure must contain only native vectors.'
    originals = ['coordinates_seed_42.csv', 'map_structure_groups.json', 'case.json', 'case_members.csv']
    provenance = {'template': 'infographic-engine', 'skill': 'gpt-image-2-style-library',
                  'scope': 'Figure 1 only. Image generation used for layout exploration, not scientific content.',
                  'input_sha256': {name: hashlib.sha256((SOURCE/name).read_bytes()).hexdigest() for name in originals},
                  'all_map_cells_preserved': len(data), 'group_counts_recomputed': checks,
                  'callout_targets': callouts, 'aligned_structures': drawing_checks,
                  'core_pixel_coordinates_equal': True, 'final_raster_embeds': 0}
    (output/'figure_checks.json').write_text(json.dumps(provenance, indent=2))
    image = Image.open(output/'Figure_1_main_infographic.png')
    image.thumbnail((1700, 1700))
    image.save(output/'Figure_1_preview.png')
    print(output/'Figure_1_main_infographic.png')


if __name__ == '__main__':
    main()
