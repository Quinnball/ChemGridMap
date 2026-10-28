"""Native-vector figures for structural comparison, with unchanged saved layouts."""
from pathlib import Path
import argparse
import json

import numpy as np
import pandas as pd
from rdkit import Chem
import matplotlib.pyplot as plt

from build_revision_figures import (ROOT, text, rect, document, export, embed_svg, chart_svg,
                                    _svg_molecule_content, INK, BLUE, ORANGE, TEAL, map_panel)
from chemgridmap.comparison import aligned_comparison

SOURCE = ROOT/'paper/output/structure_comparison'
ARCHIVE = ROOT/'paper/output/revision_v6'
OUT = ROOT/'paper/output/structure_comparison/figures'


def molecule(smiles, x, y, w, h):
    raw = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 240 240">' + (_svg_molecule_content(smiles, 240) or '') + '</svg>'
    return embed_svg(raw, x, y, w, h)


def group_map(data, highlighted):
    side = int(max(data.grid_row.max(), data.grid_col.max())) + 1
    step = 320 / (side - 1)
    body = ''
    for row in data.itertuples():
        color = highlighted.get(row.structure_group, '#D9DFE3')
        body += rect(row.grid_x*320-step*.44, (1-row.grid_y)*320-step*.44,
                     step*.88, step*.88, color, opacity=.9 if row.structure_group in highlighted else .65)
    return body


def figure1(data, groups, case):
    main = case['group']
    display = [main] + [g for g in groups if g['structure_group'] != main['structure_group']][:2]
    colors = {g['structure_group']:g['color'] for g in display}
    body = text(20, 27, 'A  Structure-defined groups on a fixed grid', 18, True)
    body += '<g transform="translate(25,63)">' + group_map(data, colors) + '</g>'
    body += text(25, 407, '889 molecules | 314 groups | three groups highlighted', 14)
    body += text(25, 428, 'Gray: other groups. Position does not define membership.', 14)
    body += text(430, 27, 'B  Scaffold coverage and representative IDs', 18, True)
    for i, g in enumerate(display):
        y = 46 + i*120
        body += rect(433, y, 545, 112, '#FAFBFC', '#D5DDE3')
        body += rect(433, y, 5, 112, g['color'])
        body += molecule(g['core_smiles'], 445, y+3, 130, 100)
        body += text(587, y+24, f"{g['structure_group']} | {g['size']} members", 18, True, g['color'])
        body += text(587, y+49, f"Core matches {g['core_coverage_count']}/{g['size']}; exact scaffold in {g['core_exact_scaffold_count']}", 15)
        body += text(587, y+73, f"Representative: {g['representative_id']}", 15)
        body += text(587, y+96, g['core_status'].capitalize(), 14)
    body += text(20, 469, 'C  Core-aligned comparison within G002', 18, True)
    body += text(20, 491, 'Three structure-selected members. Orange: outside the matched core; top strip: observed activity.', 14)
    members = pd.read_csv(SOURCE/'chembl205/case_members.csv').head(3)
    svg, _ = aligned_comparison(members, main['core_smiles'], columns=3, cell_width=320, cell_height=215)
    body += embed_svg(svg, 20, 503, 960, 215)
    body += text(500, 748, 'ChEMBL CSV  >  normalized structures  >  fingerprint groups + map  >  aligned comparison  >  evidence review', 14, anchor='middle')
    export(document(body, height=765), OUT/'Figure_1_structural_comparison')


def figure2(data, groups, case):
    main = case['group']
    group = data[data.structure_group.eq(main['structure_group'])]
    absent = group[~group.comparison_id.isin(main['core_member_ids'])].sort_values('comparison_id').iloc[0]
    representative = data[data.comparison_id.eq(main['representative_id'])].iloc[0]
    pair = pd.DataFrame([representative, absent])
    body = text(20, 28, 'A  Similarity membership is not a shared-core guarantee', 17, True)
    body += text(20, 56, 'G002: benzothiazole core matches 27 of 28 members.', 15)
    svg, _ = aligned_comparison(pair, main['core_smiles'], columns=2, cell_width=310, cell_height=260)
    body += embed_svg(svg, 20, 74, 590, 245)
    body += text(25, 345, 'Left: core present. Right: the selected core is absent.', 16, True)
    body += text(25, 370, 'The latter contains a benzothiophene ring system.', 15)
    body += text(25, 394, 'Both belong to the same fingerprint-defined group;', 15)
    body += text(25, 416, 'the coverage check prevents a false common-core label.', 15)
    body += text(650, 28, 'B  Inspect the observed annotation', 17, True)
    row = data[data.comparison_id.eq('CHEMBL20')].iloc[0]
    body += molecule(row.canonical_smiles, 670, 62, 220, 160)
    body += text(650, 248, 'Acetazolamide | CHEMBL20', 16, True)
    body += text(650, 277, 'Median pChEMBL 7.185 | active', 15)
    body += text(650, 305, '54 measurements | 27 documents', 15)
    body += text(650, 333, 'Range 4.73-8.47', 15)
    body += text(650, 370, 'A separate evidence-audit example:', 14)
    body += text(650, 392, 'not a member of the G002 comparison.', 14)
    body += '<path d="M20,442 H980" stroke="#D5DDE3"/>'
    body += text(20, 471, 'C  Original records reveal deliberately varied assay conditions', 18, True)
    source = pd.read_csv(ARCHIVE/'acetazolamide_condition_series.csv')
    fig, ax = plt.subplots(figsize=(4.45, 2.05))
    for label, color, marker in [('DTT',BLUE,'o'),('2-ME',ORANGE,'s'),('TCEP',TEAL,'^'),('Threitol','#806755','D')]:
        subset = source[source.additive.eq(label)].sort_values('concentration_mM')
        ax.plot(subset.concentration_mM, subset.pchembl_value, c=color, marker=marker, ms=3.5, lw=.8, label=label)
    ax.axhline(8.05, c=INK, ls='--', lw=.7)
    ax.set(xscale='log', xlabel='Additive concentration (mM)', ylabel='pChEMBL', ylim=(6.2,8.35), xlim=(.008,13))
    ax.legend(frameon=False, ncol=4, fontsize=7.5, loc='upper center', bbox_to_anchor=(.5,1.25))
    ax.grid(axis='y', color='#E3E8EC', lw=.5)
    fig.subplots_adjust(left=.12, right=.98, bottom=.24, top=.80)
    body += embed_svg(chart_svg(fig), 20, 495, 600, 265)
    body += text(650, 524, '25 values reconciled with source tables', 15, True)
    body += text(650, 554, 'DTT: 0 to 10 mM', 15)
    body += text(650, 580, 'IC50: 8.9 to 441 nM (~49.6-fold)', 15)
    body += text(650, 622, 'These are different assay conditions,', 15)
    body += text(650, 646, 'not same-condition technical replicates.', 15)
    body += text(650, 686, 'Records inform interpretation;', 15, True)
    body += text(650, 710, 'median colors do not harmonize assays.', 15)
    export(document(body, height=775), OUT/'Figure_2_structure_and_evidence')


def figure3():
    fig, axes = plt.subplots(2,2,figsize=(7,5.3), layout='constrained')
    targets = ['CHEMBL205','CHEMBL204','CHEMBL240']
    names = ['CA II','D3','hERG']
    full = pd.read_csv(ARCHIVE/'full_path_fidelity.csv')
    print('fidelity fields:', full.columns.tolist(), flush=True)
    ax = axes[0,0]
    for offset, stage, label, color in [(-.18,'projection','Fingerprint to projection',BLUE),(.18,'grid','Fingerprint to grid',ORANGE)]:
        values = [full[(full.target.eq(t)) & (full.destination.eq(stage))].recall.iloc[0] for t in targets]
        ax.bar(np.arange(3)+offset, values, width=.34, label=label, color=color)
    ax.set_xticks(range(3),names)
    ax.set(ylabel='Tie-aware neighbor recall', ylim=(0,1))
    ax.legend(frameon=False, fontsize=7.3, loc='upper right')
    ax.set_title('A  End-to-end neighbor retention',loc='left',fontweight='bold',fontsize=9)
    stability = pd.read_csv(ARCHIVE/'seed_stability.csv')
    ax=axes[0,1]
    for offset,col,label,color,mark in [(-.14,'projection_neighbor_jaccard','Projection',BLUE,'o'),(.14,'grid_neighbor_jaccard','Grid',ORANGE,'s')]:
        for i,t in enumerate(targets):
            vals=stability[stability.target.eq(t)][col]
            ax.scatter(i+offset+np.linspace(-.04,.04,len(vals)),vals,c=color,marker=mark,s=13,label=label if i==0 else None)
            ax.plot([i+offset-.07,i+offset+.07],[vals.mean()]*2,c=INK,lw=1)
    ax.set_xticks(range(3),names)
    ax.set(ylabel='Cross-seed neighbor Jaccard',ylim=(0,1))
    ax.set_title('B  Seed-sensitive neighborhoods',loc='left',fontweight='bold',fontsize=9)
    ax.legend(frameon=False,fontsize=8)
    macro=pd.read_csv(SOURCE/'group_macro_metrics.csv')
    ax=axes[1,0]
    for offset,space,label,color,mark in [(-.14,'projection','Projection',BLUE,'o'),(.14,'grid','Grid',ORANGE,'s')]:
        for i,t in enumerate(targets):
            vals=macro[(macro.target.eq(t.lower())) & (macro.space.eq(space))].same_group_fraction
            ax.scatter(i+offset+np.linspace(-.04,.04,len(vals)),vals,c=color,marker=mark,s=15,label=label if i==0 else None)
            ax.plot([i+offset-.07,i+offset+.07],[vals.mean()]*2,c=INK,lw=1)
    random=[macro[macro.target.eq(t.lower())].random_label_expectation.mean() for t in targets]
    ax.scatter(range(3),random,c=INK,marker='_',s=100,label='Random-label expectation')
    ax.set_xticks(range(3),names)
    ax.set(ylabel='Macro same-group neighbor fraction',ylim=(0,.48))
    ax.set_title('C  Group-level enrichment',loc='left',fontweight='bold',fontsize=9)
    ax.legend(frameon=False,fontsize=7.2,loc='upper right')
    table=pd.read_csv(SOURCE/'group_layout_metrics.csv')
    avg=table[table['size']>1].groupby(['target','structure_group','space']).same_group_fraction.mean().unstack()
    ax=axes[1,1]
    for t,label,color,marker in zip(targets,names,[BLUE,TEAL,'#99704D'],['o','^','s']):
        sub=avg.loc[t.lower()]
        ax.scatter(sub.projection,sub.grid,s=10,c=color,marker=marker,alpha=.6,label=label)
    ax.plot([0,1],[0,1],c=INK,lw=.8,ls='--')
    ax.set(xlabel='Projection same-group fraction',ylabel='Grid same-group fraction',xlim=(-.02,1.02),ylim=(-.02,1.02))
    ax.set_title('D  Per-group retention',loc='left',fontweight='bold',fontsize=9)
    ax.legend(frameon=False,fontsize=8)
    for ax in axes.flat:
        ax.grid(axis='y',c='#E3E8EC',lw=.5)
        ax.tick_params(labelsize=8)
    export(chart_svg(fig),OUT/'Figure_3_structure_fidelity')


def main():
    global OUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=OUT)
    OUT = parser.parse_args().output_dir
    OUT.mkdir(parents=True, exist_ok=True)
    data=pd.read_csv(SOURCE/'chembl205/coordinates_seed_42.csv')
    groups=json.loads((SOURCE/'chembl205/map_structure_groups.json').read_text())['groups']
    case=json.loads((SOURCE/'chembl205/case.json').read_text())
    figure1(data,groups,case)
    figure2(data,groups,case)
    figure3()
    calibration()
    toc(data, case)
    export((SOURCE/'chembl205/case_comparison.svg').read_text(),OUT/'supplement/Structure_case_nine_members')
    body = text(20,28,'All projection-to-grid displacement vectors',18,True)
    body += '<g transform="translate(40,60) scale(2.8)">' + map_panel(data,'arrows') + '</g>'
    export(document(body,height=940),OUT/'supplement/Projection_displacement_vectors')
    print(OUT)


def calibration():
    fig, axes = plt.subplots(1, 2, figsize=(7, 2.9), layout='constrained')
    occ = pd.read_csv(ARCHIVE/'occupancy.csv')
    x = occ['occupancy_actual'] if 'occupancy_actual' in occ else occ['actual_occupancy']
    axes[0].plot(x*100, occ.projection_grid_tie_aware_recall, '-o', color=BLUE, ms=4)
    axes[0].set(xlabel='Actual grid occupancy (%)', ylabel='Projection-to-grid recall', ylim=(0,1))
    axes[0].set_title('A  Space and neighbor retention', loc='left', fontsize=9, fontweight='bold')
    fixed = pd.read_csv(ARCHIVE/'sparse_calibration.csv')
    adaptive = pd.read_csv(ARCHIVE/'adaptive_validation.csv')
    candidate_col = 'candidates' if 'candidates' in fixed else 'nearest_candidates'
    fixed = fixed[fixed[candidate_col].eq(128)]
    for i, target in enumerate(['CHEMBL205','CHEMBL204','CHEMBL240']):
        for delta, table, color, marker in [(-.15,fixed,ORANGE,'s'),(.15,adaptive,TEAL,'o')]:
            values = table[table.target.eq(target)].dense_cost_ratio
            axes[1].scatter(np.full(len(values),i+delta), values, color=color, marker=marker, s=16,
                            label=('Fixed 128' if delta < 0 else 'Adaptive') if i==0 else None)
    axes[1].axhline(1, color=INK, lw=.7)
    axes[1].set(xticks=[0,1,2], xticklabels=['CA II','D3','hERG'], ylabel='Objective / dense objective', ylim=(.9,2.6))
    axes[1].set_title('B  Sparse approximation check',loc='left',fontsize=9,fontweight='bold')
    axes[1].legend(frameon=False,fontsize=8)
    for ax in axes:
        ax.grid(axis='y',color='#E3E8EC',lw=.5)
    export(chart_svg(fig), OUT/'supplement/Assignment_calibration')


def toc(data, case):
    members = pd.read_csv(SOURCE/'chembl205/case_members.csv').head(2)
    body = text(18,28,'Structure groups',20,True) + text(233,28,'Matched cores and differences',20,True)
    body += '<g transform="translate(28,68) scale(.50)">' + group_map(data,{case['group']['structure_group']:ORANGE}) + '</g>'
    svg, _ = aligned_comparison(members,case['group']['core_smiles'],columns=2,cell_width=240,cell_height=230)
    body += embed_svg(svg,228,43,440,240)
    body += '<path d="M197 156 H221 l-7 -5 m7 5 l-7 5" stroke="#172536" fill="none" stroke-width="2"/>'
    body += text(20,315,'Defined membership',19,True)
    body += text(245,315,'Checked coverage',19,True)
    body += text(465,315,'Linked evidence',19,True)
    svg = document(body,width=690,height=365).replace('width="7.0in"', 'width="3.25in"')
    svg = svg.replace(f'height="{365/690*7.0}in"',f'height="{365/690*3.25}in"')
    export(svg,OUT/'TOC')


if __name__ == '__main__':
    main()
