"""Score real task responses descriptively; never synthesize missing participants."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def score(observation, key):
    answers = observation['answers']
    if observation['kind'] == 'structure':
        indexed = {row['id']:row for row in answers}
        if len(indexed) != len(answers) or set(indexed) - set(key):
            raise ValueError('Duplicate or unknown molecular answers.')
        correct = []
        for identity, expected in key.items():
            supplied = indexed.get(identity,{})
            correct.extend([supplied.get('motif') == ('Yes' if expected['contains_aromatic_nitrogen'] else 'No'),
                            supplied.get('activity_class') == expected['activity_class']])
    else:
        correct = []
        for field, expected in key.items():
            if field == 'is_conflicted':
                correct.append(answers.get(field) == ('Yes' if expected else 'No'))
            else:
                try:
                    tolerance = .0005 if field == 'median_pchembl' else 0
                    correct.append(abs(float(answers.get(field,''))-expected) <= tolerance + 1e-12)
                except (TypeError, ValueError):
                    correct.append(False)
    return sum(correct)/len(correct)


def analyze(paths, key_path, output):
    key = json.loads(key_path.read_text())
    rows = []
    for path in paths:
        payload = json.loads(path.read_text())
        if payload.get('schema') != 1 or 'observations' not in payload:
            raise ValueError(f'Unsupported response schema in {path.name}')
        for observation in payload['observations']:
            if observation['view'] not in ('grid','scatter') or observation['status'] not in ('submitted','abandoned'):
                raise ValueError('Invalid view or task status.')
            seconds = float(observation['elapsed_seconds'])
            if not np.isfinite(seconds) or seconds < 0:
                raise ValueError('Elapsed time must be finite and non-negative.')
            rows.append({field:observation[field] for field in ('participant_code','sequence','form','kind','view','status','elapsed_seconds')} |
                        {'correctness':score(observation,key[observation['form']]), 'notes':observation.get('notes','')})
    if not rows:
        raise ValueError('No human observations supplied; no performance result can be generated.')
    table = pd.DataFrame(rows)
    if table.duplicated(['participant_code','form']).any():
        raise ValueError('A participant completed the same form more than once; review the raw responses.')
    output.mkdir(parents=True,exist_ok=True)
    table.to_csv(output/'task_scores.csv',index=False)
    # Participants, not molecules or clicks, are the unit of descriptive comparison.
    by_person = table.groupby(['participant_code','kind','view']).agg(correctness=('correctness','mean'),
        elapsed_seconds=('elapsed_seconds','mean'), n_abandoned=('status',lambda s:sum(s.eq('abandoned'))))
    by_person.to_csv(output/'participant_summary.csv')
    differences = []
    for (person,kind), group in by_person.groupby(level=[0,1]):
        sub = group.droplevel([0,1])
        if {'grid','scatter'} <= set(sub.index):
            differences.append({'participant_code':person,'kind':kind,
                'grid_minus_scatter_correctness':float(sub.loc['grid','correctness']-sub.loc['scatter','correctness']),
                'grid_minus_scatter_seconds':float(sub.loc['grid','elapsed_seconds']-sub.loc['scatter','elapsed_seconds'])})
    pd.DataFrame(differences).to_csv(output/'paired_descriptive_differences.csv',index=False)
    summary = {'n_participants':int(table.participant_code.nunique()), 'n_observations':len(table),
        'n_abandoned':int(table.status.eq('abandoned').sum()),'n_paired_task_comparisons':len(differences),
        'scope':'Exploratory descriptive responses only. No significance, superiority, equivalence, or external-application claim is computed.',
        'human_data_authenticity':'Must be confirmed by the study facilitator; a valid JSON file is not proof of human participation.'}
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('responses',nargs='+',type=Path)
    parser.add_argument('--answer-key',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    print(json.dumps(analyze(args.responses,args.answer_key,args.output),indent=2))
