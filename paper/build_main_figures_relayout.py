"""Recompose the three manuscript figures from frozen project evidence.

Run from the repository: .venv/bin/python paper/build_main_figures_relayout.py
Layout reference: Infographic Engine (cases 334, 1, 8); no generated data art.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymupdf
from PIL import Image
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Scaffolds import MurckoScaffold

import build_merged_figures as previous
from build_revision_figures import text, rect, document, export, embed_svg, chart_svg, COLORS

ROOT = Path(__file__).resolve().parents[1]
SOURCE, ARCHIVE, AUDIT = previous.SOURCE, previous.ARCHIVE, previous.AUDIT
PREVIOUS_OUT = ROOT.parent / 'project1/JCIM_main_figures_redesign_20260916'
DEFAULT_OUT = ROOT.parent / 'project1/JCIM_main_figures_redesign_20260916_review_v2'
CAPTIONS = ROOT / 'paper/main_figure_captions_review_v2.md'
OLD_FIGURES = ROOT.parent / 'project1/JCIM_final_polish/figures'
NS = '{http://www.w3.org/2000/svg}'
INK, BLUE, TEAL, ORANGE = previous.INK, previous.BLUE, previous.TEAL, previous.ORANGE
MUTED, LINE = previous.MUTED, previous.LINE
TARGETS = ['chembl205', 'chembl204', 'chembl240']
NAMES = ['CA II', 'D3', 'hERG']


def rule(x1, y1, x2, y2, color=LINE, width=1):
    return f'<path d="M{x1},{y1} L{x2},{y2}" fill="none" stroke="{color}" stroke-width="{width}"/>'


def arrow(x1, y1, x2, y2, color=BLUE):
    a = np.arctan2(y2-y1, x2-x1)
    p = [(x2-7*np.cos(a)+s*3.5*np.sin(a), y2-7*np.sin(a)-s*3.5*np.cos(a)) for s in [-1, 1]]
    return rule(x1, y1, x2, y2, color, 1.5) + f'<path d="M{x2},{y2} L{p[0][0]},{p[0][1]} L{p[1][0]},{p[1][1]} Z" fill="{color}"/>'


def heading(x, y, letter, label):
    return text(x, y, letter, 22, True) + text(x+29, y, label, 19, True)


def molecule(smiles, x, y, width, height):
    mol = Chem.MolFromSmiles(smiles)
    before = Chem.MolToSmiles(mol, isomericSmiles=True)
    rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opt = drawer.drawOptions()
    opt.clearBackground = False
    opt.bondLineWidth = 1.7
    opt.minFontSize = 16
    opt.maxFontSize = 18
    opt.padding = .08
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    assert Chem.MolToSmiles(mol, isomericSmiles=True) == before
    return embed_svg(drawer.GetDrawingText(), x, y, width, height)


def plot(width, height, margins=(.18, .97, .22, .92)):
    fig, ax = plt.subplots(figsize=(width*.007, height*.007))
    left, right, bottom, top = margins
    fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
    ax.grid(axis='y', color=LINE, lw=.5)
    ax.tick_params(length=2.5, width=.6, pad=3)
    return fig, ax


def place_plot(fig, x, y, width, height, name, out):
    svg = chart_svg(fig)
    (out/'panels'/f'{name}.svg').write_text(svg, encoding='utf-8')
    return embed_svg(svg, x, y, width, height)


def data_copy(frame, name, out):
    frame.to_csv(out/'source_data'/name, index=False)


def finish(body, height, name, out):
    # Nested Matplotlib SVGs reuse IDs such as axes_1 and clip-path IDs.
    # Prefix IDs per nested group to preserve independent clipping on every renderer.
    root = ET.fromstring(document(body, height=height))
    for i, group in enumerate(list(root)):
        ids = {n.attrib['id']: f'p{i}_{n.attrib["id"]}' for n in group.iter() if 'id' in n.attrib}
        for node in group.iter():
            for key, value in list(node.attrib.items()):
                if key == 'id':
                    node.set(key, ids[value])
                elif value.startswith('#') and value[1:] in ids:
                    node.set(key, '#'+ids[value[1:]])
                elif 'url(#' in value:
                    for old, new in ids.items():
                        value = value.replace(f'url(#{old})', f'url(#{new})')
                    node.set(key, value)
    svg = ET.tostring(root, encoding='unicode')
    export(svg, out/name)
    with Image.open(out/f'{name}.png') as im:
        im.thumbnail((1700, 1500))
        im.save(out/'previews'/f'{name}.png')
    with pymupdf.open(out/f'{name}.pdf') as pdf:
        pdf[0].get_pixmap(dpi=144, alpha=False).save(out/'previews'/f'{name}_print_scale.png')


def comparison_cards(members, absent, core):
    raw = ET.fromstring(previous.cards(members, core, 'relayout_matched', columns=3,
                                     width=312, height=282, compact=True))
    nodes = [g for g in raw if g.get('data-comparison-id')]
    scale = previous.base.CHECKS['relayout_matched'][0]['print_median_bond_length_px']
    other = ET.fromstring(previous.cards(absent, core, 'relayout_unmatched', columns=1,
                                        width=312, height=282, compact=True, bond_length=scale/1.5))
    nodes += [g for g in other if g.get('data-comparison-id')]
    checks = previous.base.CHECKS['relayout_matched'] + previous.base.CHECKS['relayout_unmatched']
    assert all(np.isclose(scale, a['print_median_bond_length_px'], rtol=.02) for a in checks)
    return nodes, checks


def figure1(data, case, out, checks):
    group = case['group']
    assert group['structure_group'] == 'G002'
    assert group['size'] == 28 and group['core_coverage_count'] == 27
    assert group['core_exact_scaffold_count'] == 17
    body = heading(20, 28, 'A', 'Locate the group')
    extent, ox, oy = 254, 33, 65
    step = extent / (max(data.grid_row.max(), data.grid_col.max()))
    selected = data[data.structure_group.eq('G002')]
    for r in data.itertuples():
        x, y = ox+r.grid_x*extent, oy+(1-r.grid_y)*extent
        body += f'<g data-cell="{r.comparison_id}" data-row="{r.grid_row}" data-col="{r.grid_col}">'
        body += rect(x-.44*step, y-.44*step, .88*step, .88*step,
                     BLUE if r.structure_group == 'G002' else '#E0E6EA')+'</g>'
    center = selected[['grid_x', 'grid_y']].median().to_numpy()
    chosen = selected.iloc[np.argmin(((selected[['grid_x','grid_y']].to_numpy()-center)**2).sum(axis=1))]
    x, y = ox+chosen.grid_x*extent, oy+(1-chosen.grid_y)*extent
    body += rule(x, y, 213, 112, BLUE, 1.5)
    body += f'<circle cx="{x}" cy="{y}" r="2.6" fill="white" stroke="{BLUE}" stroke-width="1.2"/>'
    body += text(221, 118, 'G002', 18, True, BLUE)
    body += text(20, 343, '889 entries | 314 groups', 16, True)
    body += text(20, 365, 'Blue: G002; gray: other groups', 14, color=MUTED)
    body += rule(20, 387, 311, 387)
    body += heading(20, 419, 'B', 'Selected scaffold')
    body += text(20, 449, 'G002 | 28 members', 17, True, BLUE)
    body += molecule(group['core_smiles'], 21, 470, 140, 128)
    body += text(91, 621, 'Benzothiazole', 15, True, anchor='middle')
    for y, value, label in [(496, '27 / 28', 'Core coverage'), (570, '17 / 28', 'Exact scaffold')]:
        body += text(180, y, value, 25, True, BLUE if y == 496 else INK)
        body += text(180, y+22, label, 14, color=MUTED)
    core = Chem.MolFromSmiles(group['core_smiles'])
    core_id = Chem.MolToSmiles(core, isomericSmiles=True)
    categories = []
    for row in selected.itertuples():
        mol = Chem.MolFromSmiles(row.canonical_smiles)
        present = mol.HasSubstructMatch(core, useChirality=True)
        exact = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol), isomericSmiles=True) == core_id
        category = 'same_scaffold' if exact else 'core_present_different_scaffold' if present else 'core_absent'
        categories.append(dict(comparison_id=row.comparison_id, core_present=present, exact_scaffold=exact, category=category))
    category_counts = pd.DataFrame(categories).category.value_counts().to_dict()
    assert category_counts == {'same_scaffold':17, 'core_present_different_scaffold':10, 'core_absent':1}
    body += text(20, 650, '28 members = 17 + 10 + 1', 15, True)
    for y, count, label in [(671,17,'Same scaffold'), (692,10,'Core present; different scaffold'), (713,1,'Core absent')]:
        body += text(20, y, count, 15, True, BLUE)
        body += text(51, y, label, 14, color=MUTED)
    body += rule(330, 45, 330, 710)
    body += heading(350, 28, 'C', 'Compare aligned members')
    body += text(350, 54, 'Orange: atoms and bonds outside the selected core', 15, color=MUTED)
    members = pd.read_csv(SOURCE/'chembl205/case_members.csv').head(3)
    absent = data[data.comparison_id.eq('CHEMBL85636')]
    nodes, audits = comparison_cards(members, absent, group['core_smiles'])
    roles = []
    all_rows = pd.concat([members, absent], ignore_index=True)
    for i, (node, row, pos) in enumerate(zip(nodes, all_rows.itertuples(), [(348,78),(674,78),(348,398),(674,398)])):
        x, y = pos
        mol = Chem.MolFromSmiles(row.canonical_smiles)
        matched = mol.HasSubstructMatch(core, useChirality=True)
        exact = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol), isomericSmiles=True) == core_id
        label = ('Same scaffold' if exact else 'Core present; different scaffold') if matched else 'Core absent'
        panel = f'C{i+1}' if i < 3 else 'D'
        body += rule(x+8, y, x+305, y, COLORS[row.activity_class], 2)
        if i < 3:
            body += text(x+8, y+27, panel, 16, True, BLUE)
        else:
            body += heading(x+8, y+27, 'D', 'Retain the exception')
            body += text(x+8, y+52, 'No forced core alignment', 14, color=MUTED)
        drawing = copy.deepcopy(next(g for g in node if g.tag == NS+'g'))
        drawing.set('transform', 'translate(5,31)')
        body += f'<g data-comparison-id="{row.comparison_id}" transform="translate({x},{y})">'
        body += ET.tostring(drawing, encoding='unicode')+'</g>'
        body += text(x+8, y+254, row.comparison_id, 16, True)
        body += text(x+8, y+276, f'pChEMBL {row.activity_pchembl:.2f} | {row.activity_class}', 15, color=MUTED)
        body += text(x+8, y+299, label, 15, True)
        roles.append(dict(panel=panel, comparison_id=row.comparison_id, core_matched=matched,
                          exact_scaffold=exact, pchembl=float(row.activity_pchembl)))
    body += text(350, 735, 'C1-C3: aligned core and bond scale. Teal strip: median pChEMBL >=7 (display class).', 13.5, color=MUTED)
    checks['figure1'] = dict(member_roles=roles, drawing_audit=audits, all_grid_entries=len(data),
                             highlighted_entries=len(selected), callout_id=chosen.comparison_id,
                             callout_label_has_background=False, group_categories=category_counts,
                             shared_bond_length_px=audits[0]['print_median_bond_length_px'])
    assert [r['exact_scaffold'] for r in roles] == [True, False, True, False]
    data_copy(pd.DataFrame(roles), 'Figure_1_member_roles.csv', out)
    data_copy(pd.DataFrame(categories), 'Figure_1_all_28_member_categories.csv', out)
    finish(body, 750, 'Figure_1_structural_comparison', out)


def figure2(data, out, checks):
    row = data[data.comparison_id.eq('CHEMBL20')].iloc[0]
    body = heading(20, 28, 'A', 'Select an annotation')
    sub = data[data.grid_col.between(row.grid_col-4,row.grid_col+4) & data.grid_row.between(row.grid_row-4,row.grid_row+4)]
    for r in sub.itertuples():
        x = 23+(r.grid_col-row.grid_col+4)*13
        y = 63+(row.grid_row+4-r.grid_row)*13
        body += f'<g data-cell="{r.comparison_id}" data-row="{r.grid_row}" data-col="{r.grid_col}">'
        body += rect(x, y, 11.5, 11.5, COLORS[r.activity_class], opacity=.72)
        if r.comparison_id == row.comparison_id:
            body += rect(x-1.7, y-1.7, 14.9, 14.9, 'none', BLUE, sw=1.5)
        body += '</g>'
    body += molecule(row.canonical_smiles, 157, 58, 178, 135)
    body += text(20, 222, 'CHEMBL20 | acetazolamide', 17, True)
    body += text(20, 249, 'Median pChEMBL 7.185 | active', 16, True)
    body += text(20, 274, '54 records | 27 documents', 15)
    body += text(20, 300, 'Range 4.73-8.47 | conflict flagged', 14, color=MUTED)
    body += heading(402, 28, 'B', 'Linked ChEMBL records')
    body += rule(370, 50, 370, 309)
    records = pd.read_csv(ARCHIVE/'chembl205/acetazolamide_retained_records.csv').sort_values(['pchembl_value','source_row'])
    records['primary_table_reconciled'] = records.document_chembl_id.eq('CHEMBL1141531')
    assert len(records) == 54 and records.primary_table_reconciled.sum() == 25
    assert np.isclose(records.pchembl_value.median(), 7.185)
    body += f'<circle cx="415" cy="57" r="4.2" fill="{BLUE}"/>'
    body += text(428, 62, '25 author-matched', 14, color=MUTED)
    body += f'<circle cx="616" cy="57" r="4.2" fill="#A9B2B8"/>'
    body += text(629, 62, '29 other records', 14, color=MUTED)
    fig, ax = plot(580, 237, (.15, .97, .25, .94))
    for status, color in [(False,'#A9B2B8'),(True,BLUE)]:
        mask = records.primary_table_reconciled.to_numpy() == status
        ax.scatter(np.arange(1,55)[mask], records.pchembl_value.to_numpy()[mask], c=color, s=15,
                   linewidth=.5, edgecolor='white', zorder=3)
    ax.axhline(7.185, c=INK, lw=.7, ls='--')
    ax.text(1, 7.32, 'Median 7.185', color=INK, fontsize=7.5)
    ax.set(xlabel='Records ordered by pChEMBL', ylabel='Database pChEMBL',
           ylim=(4.5,8.8), xlim=(-1,56), yticks=[5,6,7,8])
    body += place_plot(fig, 398, 73, 580, 237, 'Figure_2B_all_records', out)
    body += text(20, 334, 'Software: record access / CSV export', 14, color=MUTED)
    body += text(407, 330, '25 / 54 linked to one publication', 15, True, BLUE)
    body += text(20, 355, 'Blue outline marks CHEMBL20.', 14, color=MUTED)
    body += text(407, 351, 'Other 29: not individually source-checked here', 14, color=MUTED)
    body += rule(20, 375, 980, 375)
    body += rule(705, 325, 825, 325, BLUE, 1.5) + arrow(825, 325, 825, 389)
    body += heading(20, 410, 'C', 'Author-led source-table verification')
    body += text(20, 437, '24 additive conditions + 1 shared zero-additive reference', 15, color=MUTED)
    source = pd.read_csv(ROOT/'paper/data/innocenti_2008_ca2.csv')
    fig, ax = plot(802, 296, (.11,.985,.24,.88))
    for key, label, color, marker in [('DTT','DTT',BLUE,'o'),('MET','2-ME',ORANGE,'s'),
                                      ('TCP','TCEP',TEAL,'^'),('THT','Threitol','#866452','D')]:
        s = source[source.additive.eq(key)].sort_values('concentration_mM')
        assert len(s) == 6
        ax.plot(s.concentration_mM, s.primary_ic50_nM, color=color, marker=marker,
                ms=4, lw=1, label=label)
    ax.axhline(8.9, c=INK, lw=.7, ls='--', label='0 mM reference (8.9 nM)')
    ax.set(xscale='log', xlabel='Nonzero additive concentration (mM)', ylabel=r'IC$_{50}$ (nM)',
           xlim=(.008,13), ylim=(0,490))
    ax.legend(frameon=False, ncol=5, loc='lower left', bbox_to_anchor=(-.015,1.025),
              fontsize=7, columnspacing=.8, handlelength=1.5, handletextpad=.45, borderaxespad=0)
    body += place_plot(fig, 8, 452, 802, 296, 'Figure_2C_experimental_conditions', out)
    body += text(839, 518, 'DTT', 17, True, BLUE)
    body += text(839, 557, '49.6-fold', 26, True, BLUE)
    body += text(839, 583, 'relative to 0 mM', 14, color=MUTED)
    body += text(839, 607, '8.9 to 441 nM', 15)
    body += text(839, 630, '0 to 10 mM', 14, color=MUTED)
    body += rule(839, 647, 980, 647)
    body += text(839, 670, 'Not same-condition', 14, True)
    body += text(839, 691, 'replicate records.', 14, True)
    body += text(20, 769, 'Source: Innocenti et al. (2008), Tables 2-5. Lines connect reported values; no model fit.', 14, color=MUTED)
    data_copy(records, 'Figure_2_all_54_records.csv', out)
    data_copy(source, 'Figure_2_source_conditions.csv', out)
    checks['figure2'] = dict(records=54, matched=25, other=29, nonzero_conditions=24, shared_reference=1,
                            median=7.185, reference_nM=8.9, endpoint_nM=441, fold=441/8.9,
                            selected_id=row.comparison_id, selected_row=int(row.grid_row), selected_col=int(row.grid_col),
                            software_scope='Linked record access and CSV export',
                            source_table_scope='Author-led verification for this study; not automatic PDF extraction',
                            zero_reference_display='One horizontal line, not points on the log concentration axis')
    finish(body, 786, 'Figure_2_activity_evidence', out)


def figure3(data, case, out, checks):
    body = heading(18, 28, 'A', 'Selected-core coverage')
    body += heading(352, 28, 'B', 'Neighbor retention')
    body += heading(686, 28, 'C', 'Seed stability')
    body += text(18, 53, 'Counts: coverage > 0.5', 14, color=MUTED)
    body += text(352, 53, 'Mean recall across molecules', 14, color=MUTED)
    body += text(686, 53, 'Pairs share layout runs', 14, color=MUTED)
    coverage = pd.read_csv(AUDIT/'structures/all_group_coverage.csv')
    coverage = coverage[coverage['size'] > 1].copy()
    core_status = {}
    for target in TARGETS:
        groups = json.loads((SOURCE/target/'map_structure_groups.json').read_text())['groups']
        core_status.update({(target,g['structure_group']): bool(g['core_smiles']) for g in groups})
    coverage['has_ring_core'] = [core_status[(r.target,r.group)] for r in coverage.itertuples()]
    fig, ax = plot(331, 263, (.23,.96,.31,.90))
    # Core availability is categorical; do not place undefined coverage on zero.
    no_core = fig.add_axes([.23,.20,.73,.07], sharex=ax)
    no_core.set_facecolor('#F0F3F5')
    no_core.set_ylim(0,1)
    no_core.set_yticks([])
    for spine in no_core.spines.values():
        spine.set_visible(False)
    no_core.tick_params(axis='x', length=0, pad=5, labelsize=7)
    no_core.text(-.025,.5,'No core',transform=no_core.transAxes,ha='right',va='center',fontsize=6.5,color=MUTED)
    coverage_counts = []
    for i, (target, color) in enumerate(zip(TARGETS,[BLUE,TEAL,'#866452'])):
        vals = coverage[coverage.target.eq(target)].sort_values(['coverage','group'])
        jitter = np.random.default_rng(42+i).uniform(-.22,.22,len(vals))
        present = vals.has_ring_core.to_numpy()
        ax.scatter(i+jitter[present], vals.coverage.to_numpy()[present], s=7, color=color, alpha=.65, linewidth=0)
        if (~present).any():
            no_core.scatter(np.full((~present).sum(),i), np.full((~present).sum(),.5),
                            s=23, facecolor='none', edgecolor=color, marker='^', lw=.8)
        count = int((vals.coverage > .5).sum())
        ax.text(i, 1.09, f'{count}/{len(vals)}', fontsize=7.5, fontweight='bold', ha='center')
        coverage_counts.append([target,count,len(vals)])
    ax.axhline(.5, c=INK, lw=.7)
    ax.set(ylabel='Selected-core coverage', ylim=(-.03,1.19), xticks=range(3), xticklabels=NAMES,
           yticks=[0,.25,.5,.75,1])
    ax.tick_params(axis='x', labelbottom=False, length=0)
    no_core.set_xticks(range(3),NAMES)
    body += place_plot(fig, 0, 66, 331, 263, 'Figure_3A_core_coverage', out)
    body += text(18, 338, 'No core: undefined, not zero coverage', 14, color=MUTED)
    full = pd.read_csv(ARCHIVE/'full_path_fidelity.csv')
    fig, ax = plot(331, 263, (.23,.96,.20,.90))
    for offset, stage, color, mark in [(-.18,'projection',BLUE,'P'),(.18,'grid',ORANGE,'G')]:
        vals = [full[full.target.eq(t.upper()) & full.destination.eq(stage)].recall.iloc[0] for t in TARGETS]
        bars = ax.bar(np.arange(3)+offset, vals, .32, color=color, edgecolor='white', lw=.5)
        ax.bar_label(bars, fmt='%.3f', padding=3, fontsize=6.8)
        for i in range(3):
            ax.text(i+offset,.055,mark,ha='center',color='white',fontsize=6.5,fontweight='bold')
    ax.set(ylabel='Tie-aware neighbor recall', ylim=(0,1), xticks=range(3), xticklabels=NAMES,
           yticks=[0,.25,.5,.75,1])
    body += place_plot(fig, 334, 66, 331, 263, 'Figure_3B_neighbor_retention', out)
    body += rect(355, 328, 12, 12, BLUE) + text(374, 338, 'P: Projection', 14, color=MUTED)
    body += rect(510, 328, 12, 12, ORANGE) + text(529, 338, 'G: Grid', 14, color=MUTED)
    stability = pd.read_csv(ARCHIVE/'seed_stability.csv')
    seed_counts = {t:len(set(stability.loc[stability.target.eq(t.upper()),'seed_a']) |
                         set(stability.loc[stability.target.eq(t.upper()),'seed_b'])) for t in TARGETS}
    fig, ax = plot(331, 263, (.23,.96,.20,.90))
    for offset, space, color, mark in [(-.14,'projection',BLUE,'o'),(.14,'grid',ORANGE,'s')]:
        for i, target in enumerate(TARGETS):
            vals = stability[stability.target.eq(target.upper())][space+'_neighbor_jaccard']
            ax.scatter(i+offset+np.linspace(-.035,.035,len(vals)), vals, c=color, marker=mark, s=12)
            ax.plot([i+offset-.075,i+offset+.075],[vals.mean()]*2,c=INK,lw=1)
    ax.set(ylabel='Cross-seed neighbor Jaccard', ylim=(0,1), xticks=range(3),
           xticklabels=[f'{n}\n{seed_counts[t]} seeds' for t,n in zip(TARGETS,NAMES)],
           yticks=[0,.25,.5,.75,1])
    body += place_plot(fig, 668, 66, 331, 263, 'Figure_3C_seed_stability', out)
    body += f'<circle cx="692" cy="334" r="4.5" fill="{BLUE}"/>' + text(704,338,'Projection',14,color=MUTED)
    body += rect(824,329,9,9,ORANGE) + text(842,338,'Grid',14,color=MUTED)
    body += rule(20, 365, 980, 365)
    body += heading(20, 398, 'D', 'Inspect apparent neighbors')
    body += text(20, 424, 'Query CHEMBL81688 | k = 10, including ties', 15, color=MUTED)
    center_id = case['group']['representative_id']
    assert center_id == 'CHEMBL81688'
    current = pd.read_csv(SOURCE/'chembl205/case_center_neighbors.csv')
    expected = pd.read_csv(OLD_FIGURES/'neighbor_example.csv')
    pd.testing.assert_frame_equal(current.set_index('molecule_id').sort_index(), expected.set_index('molecule_id').sort_index(), check_dtype=False)
    coordinates = data.set_index('comparison_id')
    query = coordinates.loc[center_id]
    body += molecule(query.canonical_smiles, 18, 455, 245, 150)
    body += text(140, 635, 'Q  CHEMBL81688', 17, True, anchor='middle')
    body += arrow(269, 532, 323, 532)
    relations = dict(zip(current.molecule_id,current.neighbor_relation))
    palette = {'shared':BLUE,'grid_only':ORANGE,'representation_only':TEAL}
    selected = coordinates.loc[list(relations)+[center_id]]
    lo_c, hi_c = int(selected.grid_col.min())-1, int(selected.grid_col.max())+1
    lo_r, hi_r = int(selected.grid_row.min())-1, int(selected.grid_row.max())+1
    step = min(245/(hi_c-lo_c+1),223/(hi_r-lo_r+1))
    for identity, row in coordinates.iterrows():
        if not(lo_c <= row.grid_col <= hi_c and lo_r <= row.grid_row <= hi_r):
            continue
        x, y = 344+(row.grid_col-lo_c)*step, 443+(hi_r-row.grid_row)*step
        relation = relations.get(identity)
        color = INK if identity == center_id else palette.get(relation,'#E7ECEF')
        body += f'<g data-cell="{identity}" data-row="{row.grid_row}" data-col="{row.grid_col}">'
        body += rect(x,y,step*.84,step*.84,color)
        mark = 'Q' if identity == center_id else {'shared':'S','grid_only':'G','representation_only':'F'}.get(relation)
        if mark:
            body += text(x+step*.42,y+step*.63,mark,step*.51,True,'white','middle')
        body += '</g>'
    counts = current.neighbor_relation.value_counts()
    for y, key, mark, label in [(466,'shared','S','Shared neighbors'),(527,'grid_only','G','Grid only'),(588,'representation_only','F','Fingerprint only')]:
        body += rect(650,y-18,26,26,palette[key])
        body += text(663,y+1,mark,16,True,'white','middle')
        body += text(691,y+3,str(counts[key]),29,True,palette[key])
        body += text(728,y+2,label,17,True)
    body += text(650,635,'12 grid / 11 fingerprint neighbors',14,color=MUTED)
    body += text(20, 696, 'Spatial proximity is an inspection cue; check explicit fingerprint neighbors before interpreting similarity.', 15, True)
    coverage['plotted_coverage'] = coverage.coverage.where(coverage.has_ring_core)
    coverage['display_status'] = np.where(coverage.has_ring_core,'numeric_coverage','no_ring_core_category')
    data_copy(coverage,'Figure_3_core_coverage.csv',out)
    data_copy(full,'Figure_3_neighbor_retention.csv',out)
    data_copy(stability,'Figure_3_seed_stability.csv',out)
    data_copy(current,'Figure_3_example_neighbors.csv',out)
    checks['figure3'] = dict(coverage_counts=coverage_counts, no_ring_core_groups=int((~coverage.has_ring_core).sum()),
                             neighbor_counts=counts.to_dict(), neighbor_ids_equal_to_previous=True,
                             seed_pairs=stability.groupby('target').size().to_dict(), seed_runs=seed_counts,
                             seed_pairs_independent=False,
                             undefined_core_display='Separate categorical band; retained in group denominators',
                             recall_definition='Mean over molecules of |fingerprint neighbors intersect destination neighbors| / |fingerprint neighbors|; ties included at k=10')
    assert counts.to_dict() == {'grid_only':9,'representation_only':8,'shared':3}
    finish(body, 715, 'Figure_3_coverage_and_neighborhoods', out)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    out = args.output_dir
    for sub in ['', 'panels', 'previews', 'source_data']:
        (out/sub).mkdir(parents=True, exist_ok=True)
    sources = [SOURCE/'chembl205/coordinates_seed_42.csv', SOURCE/'chembl205/case.json',
               SOURCE/'chembl205/case_members.csv', SOURCE/'chembl205/case_center_neighbors.csv',
               ARCHIVE/'chembl205/acetazolamide_retained_records.csv', ROOT/'paper/data/innocenti_2008_ca2.csv',
               AUDIT/'structures/all_group_coverage.csv', ARCHIVE/'full_path_fidelity.csv', ARCHIVE/'seed_stability.csv']
    sources += [SOURCE/t/'map_structure_groups.json' for t in TARGETS]
    untouched = list(OLD_FIGURES.glob('Figure_[123]*')) + list(PREVIOUS_OUT.rglob('*'))
    untouched += list((OLD_FIGURES.parent/'manuscript_english').glob('*'))
    untouched += list((OLD_FIGURES.parent/'manuscript_chinese').glob('*'))
    untouched += list(OLD_FIGURES.parent.glob('*.docx')) + list(OLD_FIGURES.parent.glob('*.pdf'))
    preserved = {str(p):sha(p) for p in sources+untouched if p.is_file()}
    data = pd.read_csv(sources[0])
    case = json.loads(sources[1].read_text())
    checks = {'style_reference':'Infographic Engine; cases 334, 1, 8', 'generated_image_content':False,
              'width_in':7, 'png_dpi':600, 'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sources}}
    with plt.rc_context({'font.family':'DejaVu Sans', 'font.size':7.5, 'axes.labelsize':7.5,
                         'xtick.labelsize':7, 'ytick.labelsize':7, 'svg.fonttype':'path',
                         'pdf.fonttype':42, 'axes.spines.top':False, 'axes.spines.right':False,
                         'axes.edgecolor':MUTED, 'axes.linewidth':.6, 'axes.axisbelow':True,
                         'hatch.linewidth':.5}):
        figure1(data,case,out,checks)
        figure2(data,out,checks)
        figure3(data,case,out,checks)
    assert all(sha(Path(p)) == digest for p,digest in preserved.items())
    checks['original_data_and_figures_unchanged'] = True
    checks['manuscripts_and_previous_delivery_unchanged'] = True
    unchanged_tables = []
    for path in (PREVIOUS_OUT/'source_data').glob('*.csv'):
        old = pd.read_csv(path)
        current = pd.read_csv(out/'source_data'/path.name)
        pd.testing.assert_frame_equal(old,current[old.columns],check_dtype=False)
        unchanged_tables.append(path.name)
    checks['numeric_tables_equal_to_previous_delivery'] = unchanged_tables
    checks['captions_sha256'] = sha(CAPTIONS) if CAPTIONS.is_file() else None
    checks['captions_status'] = 'included' if CAPTIONS.is_file() else 'skipped: optional caption source not found'
    checks['exports'] = {}
    for path in out.glob('Figure_*.pdf'):
        with pymupdf.open(path) as pdf:
            assert len(pdf) == 1 and not pdf[0].get_images()
            page = pdf[0]
            spans = [s for b in page.get_text('dict')['blocks'] if 'lines' in b
                     for line in b['lines'] for s in line['spans']]
            assert all(page.rect.contains(pymupdf.Rect(s['bbox'])) for s in spans)
            for i, a in enumerate(spans):
                for b in spans[i+1:]:
                    overlap = pymupdf.Rect(a['bbox']) & pymupdf.Rect(b['bbox'])
                    assert overlap.is_empty or overlap.width <= 1 or overlap.height <= 1, (a['text'], b['text'])
            checks['exports'][path.stem] = {'pdf_width_pt':pdf[0].rect.width,'pdf_height_pt':pdf[0].rect.height,
                                           'raster_images_in_pdf':0, 'native_text_overflow_or_overlap':False,
                                           'minimum_native_text_pt':min(s['size'] for s in spans)}
        root = ET.parse(path.with_suffix('.svg')).getroot()
        ids = [n.get('id') for n in root.iter() if n.get('id')]
        assert len(ids) == len(set(ids))
        with Image.open(path.with_suffix('.png')) as im:
            checks['exports'][path.stem].update(png_pixels=im.size,png_dpi=im.info.get('dpi'))
    (out/'figure_checks.json').write_text(json.dumps(checks,indent=2,default=lambda x:x.item()),encoding='utf-8')
    if CAPTIONS.is_file():
        (out/'Figure_captions_EN_ZH.md').write_bytes(CAPTIONS.read_bytes())
    else:
        print(f'Warning: optional caption source not found; skipping {CAPTIONS.name}.')
    with zipfile.ZipFile(out/'ChemGridMap_three_main_figures.zip','w',zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(out.rglob('*')):
            if path.is_file() and path.suffix != '.zip':
                archive.write(path, path.relative_to(out))
        archive.write(Path(__file__), 'source/chemgridmap/paper/'+Path(__file__).name)
        if CAPTIONS.is_file():
            archive.write(CAPTIONS, 'source/chemgridmap/paper/'+CAPTIONS.name)
    print(out)
    print(json.dumps(checks['exports'],indent=2))


if __name__ == '__main__':
    main()
