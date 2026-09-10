"""Prepare an offline, two-view task preview. Produces no participant results."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd
from rdkit import Chem
from chemgridmap.inspection import inspect_entry
from chemgridmap.plotting import _svg_molecule_content
from chemgridmap.core import assign_to_grid, normalize_coordinates

ROOT = Path(__file__).resolve().parents[1]


def ordered(frame):
    return frame.assign(_order=frame.molecule_identity_key.map(lambda value: hashlib.sha256(str(value).encode()).hexdigest())).sort_values('_order')


def prepare(source, output):
    output.mkdir(parents=True, exist_ok=True)
    directory = source/'chembl205'
    coords = pd.read_csv(directory/'coordinates_seed_42.csv')
    records = pd.read_csv(directory/'chembl205_retained_records.csv')
    parameters = json.loads((directory/'chembl205_curation_report.json').read_text())['parameters']
    query = Chem.MolFromSmarts('[#7;a]')
    coords['motif'] = coords.canonical_smiles.map(lambda s: Chem.MolFromSmiles(s).HasSubstructMatch(query))
    positive, negative = ordered(coords[coords.motif]), ordered(coords[~coords.motif])
    panels, keys = [], {}
    for version in range(2):
        selected = pd.concat([positive.iloc[version*6:(version+1)*6], negative.iloc[version*6:(version+1)*6]])
        selected = ordered(selected).drop(columns='_order')
        assert len(selected) == 12
        identifier = 'structure_' + str(version+1)
        panels.append({'id':identifier,'kind':'structure','rows':selected})
        keys[identifier] = {str(row.molecule_identity_key): {'contains_aromatic_nitrogen':bool(row.motif), 'activity_class':row.activity_class}
                            for row in selected.itertuples()}
    # Two-record groups keep retrieval workload comparable; no selection by conflict outcome.
    repeats = ordered(coords[coords.n_records.eq(2)])
    assert len(repeats) >= 2
    for version in range(2):
        selected = repeats.iloc[[version]].drop(columns='_order')
        identity = selected.iloc[0].molecule_identity_key
        # Distractors are fixed without using their labels or layout quality.
        distractors = ordered(coords[~coords.molecule_identity_key.eq(identity)]).iloc[version*11:(version+1)*11]
        frame = ordered(pd.concat([selected,distractors])).drop(columns='_order')
        identifier = 'records_' + str(version+1)
        panels.append({'id':identifier,'kind':'records','target_id':identity,'rows':frame})
        summary = inspect_entry(coords,records,molecule_id=identity,parameters=parameters)[2]
        keys[identifier] = {key:summary[key] for key in ('n_records','median_pchembl','n_documents','is_conflicted')}
    payload = []
    for panel in panels:
        points = []
        projection = normalize_coordinates(panel['rows'][['projection_x', 'projection_y']].to_numpy())
        grid, grid_rows, _ = assign_to_grid(projection, target_occupancy=.4, method='dense')
        for index, row in enumerate(panel['rows'].itertuples()):
            chosen, evidence, summary = inspect_entry(coords,records,molecule_id=row.molecule_identity_key,parameters=parameters)
            fields = [key for key in ('activity_id','assay_chembl_id','document_chembl_id','activity_pchembl_raw','assay_description') if key in evidence]
            points.append({'id':str(row.molecule_identity_key),'label':row.activity_class,
                'median':float(row.activity_pchembl),'gx':float(grid[index,0]),'gy':float(grid[index,1]),
                'px':float(projection[index,0]),'py':float(projection[index,1]),
                'svg':_svg_molecule_content(row.canonical_smiles,240),'summary':summary,
                'records':json.loads(evidence[fields].to_json(orient='records'))})
        spacing = 1.0 / (grid_rows - 1)
        payload.append({key:value for key,value in panel.items() if key != 'rows'} | {'points':points,'grid_step':spacing})
    protocol = {'status':'prepared_not_run', 'participants_observed':0, 'target':'CHEMBL205',
        'comparison':'Two visualization modes of the same evidence interface; not a comparison of complete third-party applications.',
        'fixed_features':['same molecules and source metadata within each form','same inspector and record download','same zoom and pan controls','same activity colors'],
        'display_difference':'Native structure tiles versus clickable scatter points for each 12-molecule task form. Saved projection positions are isotropically rescaled; dense assignment uses the default 0.40 occupancy. This task layout is separate from the full-data scientific figures.',
        'structure_task':'Identify aromatic nitrogen and report the activity class for every displayed molecule; SMARTS [#7;a] defines the answer key.',
        'structure_sampling':'Two disjoint sets; six motif-positive and six motif-negative entries per set, ordered by SHA256(identity). Not prevalence estimation.',
        'record_task':'Report measurement count, median, document count, and conflict flag for a specified two-record molecule.',
        'record_sampling':'First two two-record groups by SHA256(identity); fixed distractors.',
        'assignment':'Four counterbalanced sequences; each form appears in both views across participants. No participant repeats a form.',
        'primary_outcome':'Participant-level correctness for the structure task; do not pool molecules as independent participants.',
        'secondary_outcomes':['elapsed time, with failures retained','record-task correctness','task abandonment and notes'],
        'analysis':'Exploratory paired descriptive differences by participant, stratified by task. No general superiority, equivalence, or statistical significance claim from this preview.',
        'before_collection':['Human review of molecular answer keys','Confirm institutional ethics and consent requirements','Fix recruitment count and stopping rule before viewing outcomes','Provide equal practice on both views using non-test examples'],
        'data_policy':'Anonymous study code only; files remain local; no names, email, credentials or automatic uploads.',
        'provenance':{name:hashlib.sha256((directory/name).read_bytes()).hexdigest() for name in ('coordinates_seed_42.csv','chembl205_retained_records.csv','chembl205_curation_report.json')}}
    encoded = json.dumps(payload,ensure_ascii=True).replace('</','<\\/')
    template = (Path(__file__).parent/'task_viewer.html').read_text()
    (output/'task_preview.html').write_text(template.replace('/* TASK_DATA */', 'const TASKS = ' + encoded + ';'),encoding='utf-8')
    (output/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    (output/'answer_key.json').write_text(json.dumps(keys,indent=2)+'\n')
    print(f'Prepared 4 task forms in {output}; zero human observations.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=ROOT/'paper/output/revision_v6')
    parser.add_argument('--output',type=Path,default=ROOT/'paper/output/current/task_materials')
    arguments = parser.parse_args()
    prepare(arguments.source,arguments.output)
