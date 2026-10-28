"""Restyle the manuscript figures without changing data, layouts or membership."""
from pathlib import Path
import argparse
import hashlib
import json
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import build_structure_figures as previous
from build_revision_figures import (ROOT, text, rect, document, export, embed_svg,
                                    chart_svg, INK, BLUE, ORANGE, TEAL)
from chemgridmap.comparison import aligned_comparison

SOURCE, ARCHIVE = previous.SOURCE, previous.ARCHIVE
OUT = SOURCE / 'figures_submission'
MUTED, LINE = '#536575', '#D9E1E7'
NS = '{http://www.w3.org/2000/svg}'
CHECKS = {}
plt.rcParams.update({'font.size': 8, 'axes.labelsize': 8, 'axes.titlesize': 9,
                     'xtick.labelsize': 8, 'ytick.labelsize': 8,
                     'legend.fontsize': 7, 'axes.titlepad': 9, 'axes.axisbelow': True})


def panel(x, y, letter, title):
    return text(x, y, letter, 19, True) + text(x + 27, y, title, 17, True)


def cards(frame, core, name, columns=3, width=320, height=224):
    svg, checks = aligned_comparison(frame, core, columns=columns,
                                     cell_width=width, cell_height=height)
    assert all(item['structure_preserved'] for item in checks)
    CHECKS[name] = checks
    root = ET.fromstring(svg)
    # Restyle surrounding cards only; retain all atom paths and correspondences.
    for node in root.iter():
        if node.tag == NS + 'rect' and node.get('rx'):
            node.set('rx', '3'); node.set('stroke', LINE)
        if node.tag == NS + 'text':
            node.set('font-family', 'Arial,Helvetica,sans-serif')
            if node.get('font-size') == '13':
                node.set('font-weight', '600'); node.set('font-size', '14')
            else:
                node.set('font-size', '12')
    return ET.tostring(root, encoding='unicode')


def figure1(data, groups, case):
    main = case['group']
    display = [main] + [g for g in groups if g['structure_group'] != main['structure_group']][:2]
    colors = {g['structure_group']: g['color'] for g in display}
    body = panel(20, 28, 'A', 'Structure-defined groups')
    body += '<g transform="translate(37,52) scale(.90)">' + previous.group_map(data, colors) + '</g>'
    body += text(25, 369, '889 molecules  |  314 groups', 14, True)
    body += text(25, 390, 'Gray: other groups; positions do not define membership.', 12, color=MUTED)
    body += panel(410, 28, 'B', 'Scaffold coverage')
    for i, g in enumerate(display):
        y = 48 + i * 113
        body += rect(410, y, 570, 103, '#F8FAFB', LINE) + rect(410, y, 4, 103, g['color'])
        body += previous.molecule(g['core_smiles'], 426, y + 2, 137, 98)
        body += text(579, y + 25, f"{g['structure_group']}  |  {g['size']} molecules", 17, True, g['color'])
        body += text(579, y + 49, f"Core coverage  {g['core_coverage_count']}/{g['size']}", 15, True)
        body += text(801, y + 49, f"Exact scaffold  {g['core_exact_scaffold_count']}/{g['size']}", 13, color=MUTED)
        body += text(579, y + 73, f"Representative  {g['representative_id']}", 13)
        status = 'Majority coverage' if g['core_coverage_count'] > g['size'] / 2 else 'No majority coverage'
        body += text(579, y + 93, status, 12, color=MUTED)
    body += panel(20, 427, 'C', 'Aligned structural comparison')
    body += text(980, 427, 'G002  |  orange: outside the matched core', 12, color=MUTED, anchor='end')
    members = pd.read_csv(SOURCE / 'chembl205/case_members.csv').head(3)
    body += embed_svg(cards(members, main['core_smiles'], 'Figure_1'), 20, 442, 960, 224)
    for label, x in [('ChEMBL CSV', 98), ('Defined groups', 292), ('Aligned structures', 516), ('Linked evidence', 786)]:
        body += text(x, 699, label, 13, True, anchor='middle')
    for a, b in [(157, 212), (363, 424), (608, 698)]:
        body += f'<path d="M{a} 695 H{b} l-6 -4 m6 4 l-6 4" stroke="{MUTED}" stroke-width="1.2" fill="none"/>'
    export(document(body, height=715), OUT / 'Figure_1_structural_comparison')


def figure2(data, case):
    main = case['group']
    group = data[data.structure_group.eq(main['structure_group'])]
    absent = group[~group.comparison_id.isin(main['core_member_ids'])].sort_values('comparison_id').iloc[0]
    representative = data[data.comparison_id.eq(main['representative_id'])].iloc[0]
    body = panel(20, 28, 'A', 'A structural exception within G002')
    body += text(20, 54, 'Selected benzothiazole core: 27 of 28 members', 14, color=MUTED)
    svg = cards(pd.DataFrame([representative, absent]), main['core_smiles'], 'Figure_2', columns=2, width=295)
    body += embed_svg(svg, 20, 69, 590, 224)
    body += text(25, 320, 'Same fingerprint group; different ring systems.', 14, True)
    body += text(25, 343, 'The unmatched member remains visible.', 13, color=MUTED)
    body += panel(642, 28, 'B', 'An annotation to inspect')
    row = data[data.comparison_id.eq('CHEMBL20')].iloc[0]
    body += previous.molecule(row.canonical_smiles, 686, 50, 240, 139)
    for y, value, size, bold in [(220, 'Acetazolamide  |  CHEMBL20', 16, True),
                                (247, 'Median pChEMBL  7.185', 15, False),
                                (273, '54 records  |  27 documents', 14, False),
                                (299, 'Observed range  4.73-8.47', 14, False)]:
        body += text(642, y, value, size, bold)
    body += text(642, 333, 'Separate audit case, outside G002', 12, color=MUTED)
    body += f'<path d="M20 370 H980" stroke="{LINE}"/>'
    body += panel(20, 401, 'C', 'Assay conditions explain part of the variation')
    source = pd.read_csv(ARCHIVE / 'acetazolamide_condition_series.csv')
    fig, ax = plt.subplots(figsize=(4.35, 2.05))
    for label, color, marker in [('DTT', BLUE, 'o'), ('2-ME', ORANGE, 's'), ('TCEP', TEAL, '^'), ('Threitol', '#806755', 'D')]:
        subset = source[source.additive.eq(label)].sort_values('concentration_mM')
        ax.plot(subset.concentration_mM, subset.pchembl_value, c=color, marker=marker, ms=3.5, lw=1, label=label)
    ax.axhline(8.05, c=MUTED, lw=.7)
    ax.text(.01, 8.10, 'Zero-additive reference', fontsize=6.8, color=MUTED)
    ax.set(xscale='log', xlabel='Additive concentration (mM)', ylabel='pChEMBL', ylim=(6.2, 8.35), xlim=(.008, 13))
    ax.legend(frameon=False, ncol=4, fontsize=7, loc='upper center', bbox_to_anchor=(.5, 1.30))
    ax.grid(axis='y', color=LINE, lw=.5)
    fig.subplots_adjust(left=.115, right=.985, bottom=.23, top=.78)
    body += embed_svg(chart_svg(fig), 15, 427, 595, 268)
    body += rect(646, 453, 326, 111, '#F3F7F9')
    body += text(663, 480, '25 IC50 values match source tables', 14, True)
    body += text(663, 510, 'DTT: 0 to 10 mM', 14)
    body += text(663, 540, '8.9 to 441 nM  |  49.6-fold', 17, True, BLUE)
    body += text(650, 597, 'Different conditions,', 15, True)
    body += text(650, 620, 'not same-condition replicates.', 15, True)
    body += text(650, 656, 'Median colors summarize observations;', 12, color=MUTED)
    body += text(650, 676, 'they do not harmonize assays.', 12, color=MUTED)
    export(document(body, height=710), OUT / 'Figure_2_structure_and_evidence')


def figure3():
    fig, axes = plt.subplots(2, 2, figsize=(7, 4.65))
    fig.subplots_adjust(left=.085, right=.98, top=.92, bottom=.12, wspace=.31, hspace=.62)
    targets, names = ['CHEMBL205', 'CHEMBL204', 'CHEMBL240'], ['CA II', 'D3', 'hERG']
    full = pd.read_csv(ARCHIVE / 'full_path_fidelity.csv')
    ax = axes[0, 0]
    for offset, stage, label, color in [(-.18, 'projection', 'Fingerprint to projection', BLUE), (.18, 'grid', 'Fingerprint to grid', ORANGE)]:
        values = [full[(full.target.eq(t)) & (full.destination.eq(stage))].recall.iloc[0] for t in targets]
        bars = ax.bar(np.arange(3) + offset, values, width=.32, label=label, color=color)
        ax.bar_label(bars, fmt='%.3f', padding=3, fontsize=6.6, color=INK)
    ax.set_xticks(range(3), names); ax.set(ylabel='Neighbor recall', ylim=(0, 1))
    ax.legend(frameon=False, fontsize=6.5, loc='upper right', handlelength=1)
    ax.set_title('A  End-to-end retention', loc='left', fontweight='bold')
    stability = pd.read_csv(ARCHIVE / 'seed_stability.csv')
    macro = pd.read_csv(SOURCE / 'group_macro_metrics.csv')
    for ax, title, ylabel, is_macro in [(axes[0, 1], 'B  Seed sensitivity', 'Cross-seed neighbor Jaccard', False),
                                        (axes[1, 0], 'C  Group-level enrichment', 'Same-group neighbor fraction', True)]:
        for offset, space, label, color, mark in [(-.14, 'projection', 'Projection', BLUE, 'o'), (.14, 'grid', 'Grid', ORANGE, 's')]:
            for i, t in enumerate(targets):
                vals = (macro[(macro.target.eq(t.lower())) & (macro.space.eq(space))].same_group_fraction
                        if is_macro else stability[stability.target.eq(t)][space + '_neighbor_jaccard'])
                ax.scatter(i + offset + np.linspace(-.04, .04, len(vals)), vals, c=color, marker=mark, s=12,
                           label=label if i == 0 and not is_macro else None)
                ax.plot([i + offset - .075, i + offset + .075], [vals.mean()] * 2, c=INK, lw=.8)
        ax.set_xticks(range(3), names); ax.set(ylabel=ylabel, ylim=(-.01, .44) if is_macro else (0, 1))
        ax.set_title(title, loc='left', fontweight='bold')
        if is_macro:
            null = [macro[macro.target.eq(t.lower())].random_label_expectation.mean() for t in targets]
            ax.scatter(range(3), null, c=INK, marker='_', s=65, label='Random-label expectation', clip_on=False)
            ax.set_yticks([0, .1, .2, .3, .4])
            ax.legend(frameon=False, fontsize=6.5, loc='upper left', bbox_to_anchor=(-.03, -.22), handlelength=1)
        else:
            ax.legend(frameon=False, fontsize=7, loc='upper right', handlelength=1)
    table = pd.read_csv(SOURCE / 'group_layout_metrics.csv')
    avg = table[table['size'] > 1].groupby(['target', 'structure_group', 'space']).same_group_fraction.mean().unstack()
    ax = axes[1, 1]
    for t, label, color, marker in zip(targets, names, [BLUE, TEAL, '#99704D'], ['o', '^', 's']):
        sub = avg.loc[t.lower()]
        ax.scatter(sub.projection, sub.grid, s=10, c=color, marker=marker, alpha=.65, edgecolors='none', label=label)
    ax.plot([0, 1], [0, 1], c=MUTED, lw=.7, ls='--')
    ax.set(xlabel='Projection same-group fraction', ylabel='Grid same-group fraction', xlim=(-.025, 1.025), ylim=(-.025, 1.025))
    ax.set_title('D  Per-group retention', loc='left', fontweight='bold')
    ax.legend(frameon=False, fontsize=7, loc='upper left', handlelength=1)
    for ax in axes.flat:
        ax.grid(axis='y', c=LINE, lw=.5); ax.tick_params(width=.6, length=3)
    export(chart_svg(fig), OUT / 'Figure_3_structure_fidelity')


def main():
    global OUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=OUT)
    OUT = parser.parse_args().output_dir; OUT.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(SOURCE / 'chembl205/coordinates_seed_42.csv')
    groups = json.loads((SOURCE / 'chembl205/map_structure_groups.json').read_text())['groups']
    case = json.loads((SOURCE / 'chembl205/case.json').read_text())
    figure1(data, groups, case); figure2(data, case); figure3()
    previous.OUT = OUT
    previous.calibration(); previous.toc(data, case)
    nine = pd.read_csv(SOURCE / 'chembl205/case_members.csv')
    svg = cards(nine, case['group']['core_smiles'], 'Figure_S1', height=245)
    export(document(embed_svg(svg, 10, 10, 980, 750), height=775), OUT / 'supplement/Structure_case_nine_members')
    inputs = [SOURCE / f for f in ['chembl205/coordinates_seed_42.csv', 'chembl205/case_members.csv',
                                   'chembl205/map_structure_groups.json', 'chembl205/case.json',
                                   'group_macro_metrics.csv', 'group_layout_metrics.csv']]
    inputs += [ARCHIVE / f for f in ['acetazolamide_condition_series.csv', 'full_path_fidelity.csv',
                                     'seed_stability.csv', 'occupancy.csv', 'sparse_calibration.csv', 'adaptive_validation.csv']]
    manifest = {'scope': 'Presentation-only revision. No refitting, membership or numerical changes.',
                'inputs': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
                'drawing_checks': CHECKS}
    (OUT / 'figure_provenance.json').write_text(json.dumps(manifest, indent=2))
    print(OUT)


if __name__ == '__main__':
    main()
