"""Task-led manuscript figures from frozen coordinates and checked source records."""
from pathlib import Path
import argparse
import hashlib
import json
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image
from rdkit import Chem, Geometry
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Scaffolds import MurckoScaffold

import build_submission_figures as base
import build_story_figures as story
from build_main_figure_infographic import line, heading
from build_revision_figures import text, rect, document, export, embed_svg, chart_svg, COLORS

ROOT = Path(__file__).resolve().parents[1]
SOURCE, ARCHIVE = base.SOURCE, base.ARCHIVE
OUT = ROOT.parent/'project1/JCIM_final_polish/figures'
AUDIT = ROOT/'paper/output/merge_20260915'
INK, BLUE, ORANGE, TEAL, MUTED, LINE = base.INK, base.BLUE, base.ORANGE, base.TEAL, base.MUTED, base.LINE
NS = '{http://www.w3.org/2000/svg}'
CHECKS = {}


def cards(frame, core, name, columns=3, width=240, height=310, bond_length=None, compact=False):
    raw = base.cards(frame, core, name, columns=columns, width=width, height=height)
    root = ET.fromstring(raw)
    # Print-only redraw: retain verified atom maps, but enlarge atom labels.
    mols=[]
    template=Chem.MolFromSmiles(core) if core else None
    if template is not None: rdDepictor.Compute2DCoords(template)
    for row,audit in zip(frame.itertuples(),base.CHECKS[name]):
        mol=Chem.MolFromSmiles(row.canonical_smiles)
        if audit['aligned']:
            rdDepictor.GenerateDepictionMatching2DStructure(mol,template,list(enumerate(audit['core_atom_indices'])))
        else: rdDepictor.Compute2DCoords(mol)
        assert Chem.MolToSmiles(Chem.MolFromMolBlock(Chem.MolToMolBlock(mol)),isomericSmiles=True)==audit['canonical_isomeric_smiles_before']
        mols.append(mol)
    positions=np.concatenate([m.GetConformer().GetPositions()[:,:2] for m in mols])
    margin=.35 if compact else .45
    lo,hi=positions.min(axis=0)-margin,positions.max(axis=0)+margin
    entries=[g for g in root if g.get('data-comparison-id')]
    print_points=[]
    for node,mol,audit in zip(entries,mols,base.CHECKS[name]):
        drawer=rdMolDraw2D.MolDraw2DSVG(width-10,height-62)
        opts=drawer.drawOptions();opts.clearBackground=False;opts.bondLineWidth=1.7
        opts.minFontSize=16;opts.maxFontSize=18;opts.highlightBondWidthMultiplier=6
        if compact:opts.padding=.015
        if bond_length is None:
            drawer.SetScale(width-10,height-62,Geometry.Point2D(*lo),Geometry.Point2D(*hi))
        else:
            # Cap the independent depiction at the same bond scale as the aligned cards.
            opts.fixedBondLength=bond_length
        atoms,bonds=audit['highlight_atom_indices'],audit['highlight_bond_indices']
        c=(.98,.73,.40)
        drawer.DrawMolecule(mol,highlightAtoms=atoms,highlightBonds=bonds,
                            highlightAtomColors={a:c for a in atoms},highlightBondColors={b:c for b in bonds})
        lengths=[(drawer.GetDrawCoords(b.GetBeginAtomIdx())-drawer.GetDrawCoords(b.GetEndAtomIdx())).Length() for b in mol.GetBonds()]
        audit['print_median_bond_length_px']=float(np.median(lengths))
        if audit['aligned']:print_points.append([[drawer.GetDrawCoords(a).x,drawer.GetDrawCoords(a).y] for a in audit['core_atom_indices']])
        drawer.FinishDrawing()
        drawing=next(g for g in node if g.tag==NS+'g')
        drawing[:]=list(ET.fromstring(drawer.GetDrawingText()))
    if print_points:
        assert all(np.allclose(print_points[0],p,atol=1e-5) for p in print_points)
    for node in root.iter(NS+'rect'):
        if node.get('rx'):
            node.set('stroke', 'none')
    for node in root.iter(NS+'text'):
        node.set('font-size', '17' if node.get('font-weight') == '600' else '15')
        if 'core absent' in (node.text or ''):
            node.text = node.text.split(' | core absent')[0]
    return ET.tostring(root, encoding='unicode')


def figure1(data, case):
    group = case['group']
    body = heading(20, 28, 'A', 'Locate the group')
    side = int(max(data.grid_row.max(), data.grid_col.max()))+1
    extent, ox, oy = 270, 40, 58
    step = extent/(side-1)
    for r in data.itertuples():
        x, y = ox+r.grid_x*extent, oy+(1-r.grid_y)*extent
        body += f'<g data-cell="{r.comparison_id}" data-row="{r.grid_row}" data-col="{r.grid_col}">'
        body += rect(x-step*.44, y-step*.44, step*.88, step*.88, BLUE if r.structure_group=='G002' else '#E0E6EA')+'</g>'
    sub = data[data.structure_group.eq('G002')]
    center = sub[['grid_x','grid_y']].median().to_numpy()
    r = sub.iloc[np.argmin(((sub[['grid_x','grid_y']].to_numpy()-center)**2).sum(axis=1))]
    x, y = ox+r.grid_x*extent, oy+(1-r.grid_y)*extent
    body += line(x, y, 319, 126, BLUE)+text(327, 132, 'G002', 18, True, BLUE)
    body += f'<circle cx="{x}" cy="{y}" r="3" fill="white" stroke="{BLUE}"/>'
    CHECKS['group_callout'] = dict(molecule=r.comparison_id,row=int(r.grid_row),col=int(r.grid_col),
                                   leader_end=319,label_start=327,white_label_background=False)
    body += text(20, 351, '889 entries | 314 fingerprint-defined groups', 16, True)
    body += text(20, 374, 'Blue: G002. Gray: other groups.', 16, color=MUTED)
    body += heading(452, 28, 'B', 'Check the selected scaffold')
    body += base.previous.molecule(group['core_smiles'], 461, 73, 230, 199)
    body += text(574, 297, 'Benzothiazole', 18, True, anchor='middle')
    body += text(705, 92, 'G002 | 28 members', 19, True, BLUE)
    body += text(705, 135, 'Core coverage', 18, True)
    body += text(705, 170, '27 / 28', 29, True, BLUE)
    body += text(705, 194, 'contain the selected core', 16, color=MUTED)
    body += text(705, 240, 'Exact scaffold frequency', 18, True)
    body += text(705, 275, '17 / 28', 29, True)
    body += text(705, 299, 'have this typed Murcko scaffold', 16, color=MUTED)
    body += text(455, 351, 'Membership, core coverage and scaffold identity', 16)
    body += text(455, 374, 'are different descriptions of the same group.', 16)
    body += line(427, 49, 427, 381)+line(20, 397, 980, 397)
    body += heading(20, 434, 'C', 'Compare aligned members')
    body += text(20, 460, 'Orange: atoms and bonds outside the selected core', 16, color=MUTED)
    members = pd.read_csv(SOURCE/'chembl205/case_members.csv').head(3)
    raw=ET.fromstring(cards(members,group['core_smiles'],'main_three',3,470,260,compact=True))
    entries=[g for g in raw if g.get('data-comparison-id')]
    roles=[]
    core=Chem.MolFromSmiles(group['core_smiles'])
    core_identity=Chem.MolToSmiles(core,isomericSmiles=True)
    for i,(node,row,(x,y)) in enumerate(zip(entries,members.itertuples(),[(20,476),(510,476),(20,771)]),1):
        mol=Chem.MolFromSmiles(row.canonical_smiles)
        scaffold=Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol),isomericSmiles=True)
        assert scaffold==row.scaffold_smiles
        matched=mol.HasSubstructMatch(core,useChirality=True)
        exact=scaffold==core_identity
        assert matched
        label='Core match; same scaffold' if exact else 'Core match; different scaffold'
        roles.append(dict(panel=f'C{i}',comparison_id=row.comparison_id,core_matched=matched,
                          exact_scaffold=exact,scaffold_smiles=scaffold,label=label))
        node.set('transform',f'translate({x},{y})')
        body+=ET.tostring(node,encoding='unicode')
        body+=text(x+12,y+32,f'C{i}',17,True,BLUE)
        body+=text(x+458,y+229,label,15,color=MUTED,anchor='end')
    body += heading(510, 754, 'D', 'Keep the unmatched member')
    absent = data[data.comparison_id.eq('CHEMBL85636')]
    assert not Chem.MolFromSmiles(absent.iloc[0].canonical_smiles).HasSubstructMatch(Chem.MolFromSmiles(group['core_smiles']))
    shared_length=base.CHECKS['main_three'][0]['print_median_bond_length_px']
    # RDKit fixedBondLength is the pixel-per-coordinate-unit scale (2D bonds are 1.5 units).
    body += embed_svg(cards(absent,group['core_smiles'],'unmatched',1,464,238,bond_length=shared_length/1.5,compact=True),510,791,464,238)
    unmatched_length=base.CHECKS['unmatched'][0]['print_median_bond_length_px']
    CHECKS['matched_unmatched_bond_scale']={'matched_px':shared_length,'unmatched_px':unmatched_length}
    assert np.isclose(shared_length,unmatched_length,rtol=.02),CHECKS['matched_unmatched_bond_scale']
    body += text(512, 780, 'Selected core absent; no forced core alignment', 16, color=MUTED)
    body += line(500, 480, 500, 726)
    body += rect(504,733,476,300,'none',LINE,sw=1.5)
    body += text(20, 1055, 'C: common core orientation and scale. Thin strips and values indicate median-based activity classes.', 16)
    absent_scaffold=Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(Chem.MolFromSmiles(absent.iloc[0].canonical_smiles)),isomericSmiles=True)
    roles.append(dict(panel='D',comparison_id='CHEMBL85636',core_matched=False,exact_scaffold=False,
                      scaffold_smiles=absent_scaffold,label='Selected core absent'))
    assert [r['exact_scaffold'] for r in roles]==[True,False,True,False]
    pd.DataFrame(roles).to_csv(OUT/'Figure_1_member_roles.csv',index=False)
    CHECKS['figure1_member_roles']=roles
    CHECKS['figure1_bond_enlargement_ratio']=shared_length/20.512163892445585
    assert CHECKS['figure1_bond_enlargement_ratio']>1.15
    export(document(body,height=1073),OUT/'Figure_1_structural_comparison')
    pts=[np.array(a['core_draw_coordinates']) for a in base.CHECKS['main_three']]
    assert all(np.allclose(pts[0],p,atol=1e-5) for p in pts)
    CHECKS['three_members_core_coordinates_equal']=True


def figure2(data):
    row=data[data.comparison_id.eq('CHEMBL20')].iloc[0]
    body=heading(20,28,'A','Select the activity annotation')
    sub=data[data.grid_col.between(row.grid_col-4,row.grid_col+4)&data.grid_row.between(row.grid_row-4,row.grid_row+4)]
    for r in sub.itertuples():
        x=22+(r.grid_col-row.grid_col+4)*17
        y=68+(row.grid_row+4-r.grid_row)*17
        body+=rect(x,y,15,15,COLORS[r.activity_class],opacity=.72)
        if r.comparison_id=='CHEMBL20':
            body+=rect(x-2,y-2,19,19,'none',BLUE,sw=2)
            CHECKS['acetazolamide_cell']={'row':int(r.grid_row),'col':int(r.grid_col),'x':x,'y':y}
    body+=base.previous.molecule(row.canonical_smiles,197,57,230,159)
    body+=text(198,240,'CHEMBL20 | acetazolamide',17,True)
    body+=text(20,274,'Median pChEMBL 7.185 | active',18,True)
    body+=text(20,300,'54 records | 27 documents',17)
    body+=text(20,326,'Conflict flag: yes | range 4.73-8.47',16,color=MUTED)
    body+=text(20,349,'Blue outline: selected cell on the saved grid',15,color=MUTED)
    body+=heading(490,28,'B','Expand all 54 source records')
    records=pd.read_csv(ARCHIVE/'chembl205/acetazolamide_retained_records.csv').sort_values(['pchembl_value','source_row'])
    records['primary_table_reconciled']=records.document_chembl_id.eq('CHEMBL1141531')
    assert len(records)==54 and records.primary_table_reconciled.sum()==25
    records.to_csv(OUT/'acetazolamide_all_54_records.csv',index=False)
    fig,ax=plt.subplots(figsize=(3.4,2.05))
    for status,color,label in [(False,'#A9B2B8','29 other records'),(True,BLUE,'25 source-matched records')]:
        selected=records.primary_table_reconciled.to_numpy()==status
        ax.scatter(np.arange(1,55)[selected],records.pchembl_value.to_numpy()[selected],c=color,s=18,label=label,linewidth=.3,edgecolor='white')
    ax.axhline(7.185,color=INK,lw=.7,ls='--')
    ax.text(1,7.265,'Median 7.185',fontsize=8,color=INK,ha='left',va='bottom')
    ax.set(xlabel='Records ordered by pChEMBL',ylabel='Database pChEMBL',ylim=(4.5,8.7),xlim=(-1,56))
    ax.legend(frameon=False,fontsize=8,loc='upper left',bbox_to_anchor=(.0,1.22),ncol=1)
    ax.grid(axis='y',color=LINE,lw=.5)
    fig.subplots_adjust(left=.16,right=.98,top=.76,bottom=.23)
    body+=embed_svg(chart_svg(fig),483,48,505,258)
    body+=text(494,324,'25 / 54 reconciled to one original publication',16,True,BLUE)
    body+=text(494,349,'Other 29: not individually source-checked here',15,color=MUTED)
    body+=line(20,369,980,369)
    body+=heading(20,403,'C','Read the matched experimental conditions')
    source=pd.read_csv(ROOT/'paper/data/innocenti_2008_ca2.csv')
    fig,ax=plt.subplots(figsize=(4.35,2.35))
    for key,label,color,marker in [('DTT','DTT',BLUE,'o'),('MET','2-ME',ORANGE,'s'),('TCP','TCEP',TEAL,'^'),('THT','Threitol','#866452','D')]:
        s=source[source.additive.eq(key)].sort_values('concentration_mM')
        assert len(s)==6
        ax.plot(s.concentration_mM,s.primary_ic50_nM,c=color,marker=marker,ms=4,lw=1,label=label)
    ax.axhline(8.9,c=INK,lw=.8,ls='--',label='0 mM reference (8.9 nM)')
    ax.set(xscale='log',xlabel='Nonzero additive concentration (mM)',ylabel='IC50 (nM)',xlim=(.008,13),ylim=(0,485))
    ax.legend(frameon=False,ncol=3,loc='upper left',bbox_to_anchor=(-.06,1.31),fontsize=8,columnspacing=.8)
    ax.grid(axis='y',color=LINE,lw=.5)
    fig.subplots_adjust(left=.13,right=.98,top=.74,bottom=.21)
    body+=embed_svg(chart_svg(fig),16,422,619,322)
    for y,t,size,bold in [(454,'24 conditions + 1 shared reference',17,True),
                          (482,'25 IC50 values match Tables 2-5',17,True),
                          (517,'DTT: 0 to 10 mM',17,False),(548,'8.9 to 441 nM | 49.6-fold',19,True),
                          (591,'Different assay conditions,',17,True),(617,'not 25 same-condition replicates.',17,True),
                          (660,'Innocenti et al. (2008), pp 1900-1901',15,False),
                          (689,'Lines connect reported values; no model fit.',15,False)]:
        body+=text(653,y,t,size,bold)
    CHECKS['records']={'all':54,'matched':25,'not_reconciled_here':29,'nonzero':24,'shared_zero':1,'reference_nM':8.9,'dtt_endpoint_nM':441,'fold':441/8.9,'database_unchanged':True}
    export(document(body,height=752),OUT/'Figure_2_activity_evidence')


def axis_panel(fig,ax,title,ylabel,ylim=(0,1)):
    ax.set_title(title,loc='left',fontweight='bold',fontsize=10)
    ax.set(ylabel=ylabel,ylim=ylim)
    ax.set_xticks(range(3),['CA II','D3','hERG'])
    ax.grid(axis='y',color=LINE,lw=.5)
    fig.subplots_adjust(left=.18,right=.97,bottom=.18,top=.80)
    return chart_svg(fig)


def figure3(data,case):
    coverage=pd.read_csv(AUDIT/'structures/all_group_coverage.csv')
    coverage=coverage[coverage['size']>1]
    targets=['chembl205','chembl204','chembl240']
    core_status={}
    for t in targets:
        groups=json.loads((SOURCE/t/'map_structure_groups.json').read_text())['groups']
        core_status.update({(t,g['structure_group']):bool(g['core_smiles']) for g in groups})
    coverage['has_ring_core']=[core_status[(r.target,r.group)] for r in coverage.itertuples()]
    assert coverage[~coverage.has_ring_core].coverage.eq(0).all()
    assert coverage[coverage.has_ring_core].coverage.gt(0).all()
    CHECKS['ring_core_absent_groups']=coverage.loc[~coverage.has_ring_core,['target','group','size']].to_dict('records')
    coverage.to_csv(OUT/'nonsingleton_coverage_status.csv',index=False)
    fig,ax=plt.subplots(figsize=(3.5,2.3))
    for i,(t,c) in enumerate(zip(targets,[BLUE,TEAL,'#866452'])):
        vals=coverage[coverage.target.eq(t)].sort_values(['coverage','group'])
        # Deterministic horizontal jitter separates points without changing coverage.
        jitter=np.random.default_rng(42+i).uniform(-.22,.22,len(vals))
        present=vals.has_ring_core.to_numpy()
        ax.scatter(i+jitter[present],vals.coverage.to_numpy()[present],s=9,c=c,alpha=.65,linewidth=0)
        ax.scatter(i+jitter[~present],vals.coverage.to_numpy()[~present],s=26,facecolor='none',edgecolor=c,marker='^',linewidth=.8)
        ax.text(i,1.09,f'{(vals.coverage>.5).sum()}/{len(vals)}',ha='center',fontsize=9,fontweight='bold')
    ax.axhline(.5,c=INK,lw=.75)
    ax.set(ylabel='Selected-core coverage',ylim=(-.05,1.16))
    ax.set_xticks(range(3),['CA II','D3','hERG'])
    ax.grid(axis='y',color=LINE,lw=.5)
    fig.subplots_adjust(left=.18,right=.97,bottom=.18,top=.74)
    fig.text(.18,.86,'A  Selected-core coverage',fontsize=10,fontweight='bold')
    fig.text(.18,.80,'Line: 0.5 cutoff | Triangle: no ring core',fontsize=7.6)
    body=embed_svg(chart_svg(fig),0,0,490,333)
    full=pd.read_csv(ARCHIVE/'full_path_fidelity.csv')
    fig,ax=plt.subplots(figsize=(3.5,2.3))
    for offset,stage,label,c in [(-.17,'projection','Fingerprint to projection',BLUE),(.17,'grid','Fingerprint to grid',ORANGE)]:
        vals=[full[(full.target.eq(t.upper()))&full.destination.eq(stage)].recall.iloc[0] for t in targets]
        bars=ax.bar(np.arange(3)+offset,vals,.31,label=label,color=c)
        ax.bar_label(bars,fmt='%.3f',padding=2,fontsize=8)
    ax.legend(frameon=False,fontsize=8,loc='upper right')
    body+=embed_svg(axis_panel(fig,ax,'B  End-to-end neighbor retention','Tie-aware neighbor recall'),500,0,490,333)
    stability=pd.read_csv(ARCHIVE/'seed_stability.csv')
    fig,ax=plt.subplots(figsize=(3.5,2.3))
    for offset,space,label,c,mark in [(-.14,'projection','Projection',BLUE,'o'),(.14,'grid','Grid',ORANGE,'s')]:
        for i,t in enumerate(targets):
            vals=stability[stability.target.eq(t.upper())][space+'_neighbor_jaccard']
            ax.scatter(i+offset+np.linspace(-.035,.035,len(vals)),vals,c=c,marker=mark,s=15,label=label if i==0 else None)
            ax.plot([i+offset-.07,i+offset+.07],[vals.mean()]*2,c=INK,lw=1)
    ax.legend(frameon=False,fontsize=8,loc='upper right')
    body+=embed_svg(axis_panel(fig,ax,'C  Seed-sensitive spatial neighbors','Cross-seed neighbor Jaccard'),0,310,490,333)
    story.OUT=OUT
    body+='<g transform="translate(0,13)">'+story.neighbor_panel(data,case).replace('font-size="17"','font-size="20"')+'</g>'
    CHECKS['neighbors']=story.EXTRA_CHECKS['neighbor_example']
    export(document(body,height=684),OUT/'Figure_3_coverage_and_neighborhoods')


def engineering():
    fig,axes=plt.subplots(2,2,figsize=(7,5.2),layout='constrained')
    baseline=pd.read_csv(ARCHIVE/'baselines.csv')
    ax=axes[0,0]
    for offset,method,color in [(-.24,'Minimum-cost',BLUE),(0,'Greedy',ORANGE),(.24,'Random','#9EADB6')]:
        means=baseline[baseline.method.eq(method)].groupby('target').projection_grid_tie_aware_recall.mean()
        ax.bar(np.arange(3)+offset,[means[t] for t in ['CHEMBL205','CHEMBL204','CHEMBL240']],width=.22,color=color,label=method)
    ax.set_xticks(range(3),['CA II','D3','hERG']);ax.set(ylabel='Projection-to-grid recall',ylim=(0,1))
    ax.legend(frameon=False,fontsize=7);ax.set_title('A  Assignment controls',loc='left',fontweight='bold')
    occupancy=pd.read_csv(ARCHIVE/'occupancy.csv')
    col='occupancy_actual' if 'occupancy_actual' in occupancy else 'actual_occupancy'
    axes[0,1].plot(occupancy[col]*100,occupancy.projection_grid_tie_aware_recall,'o-',c=BLUE,ms=4)
    axes[0,1].set(xlabel='Actual occupancy (%)',ylabel='Projection-to-grid recall',ylim=(0,1))
    axes[0,1].set_title('B  Occupancy sensitivity',loc='left',fontweight='bold')
    fixed=pd.read_csv(ARCHIVE/'sparse_calibration.csv');adaptive=pd.read_csv(ARCHIVE/'adaptive_validation.csv')
    col='candidates' if 'candidates' in fixed else 'nearest_candidates'
    fixed=fixed[fixed[col].eq(128)]
    for i,t in enumerate(['CHEMBL205','CHEMBL204','CHEMBL240']):
        for offset,tab,c,m,label in [(-.12,fixed,ORANGE,'s','Fixed 128'),(.12,adaptive,TEAL,'o','Adaptive')]:
            vals=tab[tab.target.eq(t)].dense_cost_ratio
            axes[1,0].scatter(i+offset+np.linspace(-.03,.03,len(vals)),vals,c=c,marker=m,s=13,label=label if i==0 else None)
    axes[1,0].axhline(1,c=INK,lw=.7);axes[1,0].set_xticks(range(3),['CA II','D3','hERG'])
    axes[1,0].set(ylabel='Sparse cost / dense cost');axes[1,0].legend(frameon=False,fontsize=7)
    axes[1,0].set_title('C  Sparse candidate calibration',loc='left',fontweight='bold')
    large=json.loads((ARCHIVE/'large_scale_extension.json').read_text())
    hist=large['assignment']['candidate_history']
    axes[1,1].plot([r['candidates'] for r in hist],[r['squared_cost'] for r in hist],'o-',c=TEAL,ms=4)
    axes[1,1].set(xscale='log',xlabel='Candidates per molecule',ylabel='Assignment squared cost')
    candidates=[r['candidates'] for r in hist]
    axes[1,1].set_xticks(candidates,[str(value) for value in candidates])
    axes[1,1].set_xticks([],minor=True)
    axes[1,1].set_title('D  Full hERG (9581 entries)',loc='left',fontweight='bold')
    for ax in axes.flat: ax.grid(axis='y',color=LINE,lw=.5)
    export(chart_svg(fig),OUT/'supplement/Engineering_calibration')


def toc(data, case):
    members=pd.read_csv(SOURCE/'chembl205/case_members.csv').head(2)
    body=text(18,28,'Structure groups',20,True)+text(233,28,'Core-aligned members',20,True)
    body+='<g transform="translate(28,68) scale(.50)">'+base.previous.group_map(data,{'G002':BLUE})+'</g>'
    raw,_=base.previous.aligned_comparison(members,case['group']['core_smiles'],columns=2,cell_width=240,cell_height=230)
    body+=embed_svg(raw,228,43,440,240)
    body+='<path d="M197 156 H221 l-7 -5 m7 5 l-7 5" stroke="#172536" fill="none" stroke-width="2"/>'
    body+=text(20,315,'27/28 core coverage',18,True)
    body+=text(245,315,'Orange: outside the selected core',18,True)
    raw=document(body,width=690,height=365).replace('width="7.0in"','width="3.25in"')
    raw=raw.replace(f'height="{365/690*7.0}in"',f'height="{365/690*3.25}in"')
    export(raw,OUT/'TOC')


def main():
    global OUT
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=OUT)
    parser.add_argument('--figure-one-only',action='store_true',help='Update Figure 1 without redrawing other figures or the TOC.')
    parser.add_argument('--engineering-only',action='store_true',help='Update the engineering supplement without redrawing the main figures or TOC.')
    args=parser.parse_args()
    OUT=args.output_dir
    OUT.mkdir(parents=True,exist_ok=True)
    if args.engineering_only:
        engineering()
        print(OUT)
        return
    data=pd.read_csv(SOURCE/'chembl205/coordinates_seed_42.csv')
    case=json.loads((SOURCE/'chembl205/case.json').read_text())
    if args.figure_one_only:
        checks_path=OUT/'figure_checks.json'
        previous=json.loads(checks_path.read_text()) if checks_path.exists() else {}
        figure1(data,case)
        previous.update(CHECKS)
        previous.setdefault('drawing_checks',{}).update(base.CHECKS)
        checks_path.write_text(json.dumps(previous,indent=2,default=lambda value:value.item()))
        preview=Image.open(OUT/'Figure_1_structural_comparison.png')
        preview.thumbnail((1500,1500));preview.save(OUT/'Figure_1_structural_comparison_preview.png')
        print(OUT)
        return
    figure1(data,case); figure2(data); figure3(data,case); engineering()
    # Preserve the previous full quantitative result, including group enrichment, in SI.
    old=base.OUT; base.OUT=OUT/'supplement'; base.figure3(); base.OUT=old
    toc(data,case)
    whole=data[data.structure_group.eq('G002')].sort_values(['scaffold_smiles','canonical_smiles','comparison_id'])
    whole.to_csv(OUT/'G002_all_28_members.csv',index=False)
    for page,start in enumerate(range(0,28,9),1):
        frame=whole.iloc[start:start+9]
        raw=cards(frame,case['group']['core_smiles'],f'all_members_{page}',3,320,280)
        export(raw,OUT/f'supplement/G002_members_{page}')
    CHECKS['drawing_checks']=base.CHECKS
    CHECKS['input_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
        for p in [SOURCE/'chembl205/coordinates_seed_42.csv',SOURCE/'chembl205/case_members.csv',AUDIT/'structures/all_group_coverage.csv',ROOT/'paper/data/innocenti_2008_ca2.csv']}
    (OUT/'figure_checks.json').write_text(json.dumps(CHECKS,indent=2,default=lambda value:value.item()))
    for p in OUT.glob('Figure_*.png'):
        if '_preview' in p.stem:
            continue
        image=Image.open(p); image.thumbnail((1500,1500)); image.save(p.with_name(p.stem+'_preview.png'))
    print(OUT)


if __name__=='__main__':
    main()
