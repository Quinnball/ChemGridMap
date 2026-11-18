"""Recheck frozen groups and report the denominator of every depiction check."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd
from rdkit import Chem
from chemgridmap.comparison import group_comparison, structural_groups
from chemgridmap.provenance import environment_versions
from chemgridmap.neighborhoods import neighborhood_evidence
from chemgridmap.representations import morgan_fingerprints

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'paper/output/structure_comparison'
OUT = ROOT / 'paper/output/merge_20260915/structures'


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    groups_checked, instances, populations = [], [], []
    hashes = {}
    for target in ['chembl205', 'chembl204', 'chembl240']:
        data_path = SOURCE / target / 'coordinates_seed_42.csv'
        data = pd.read_csv(data_path)
        saved = json.loads((SOURCE / target / 'map_structure_groups.json').read_text())['groups']
        current, groups, _ = structural_groups(data)
        assert current.set_index('comparison_id').structure_group.sort_index().equals(
            data.set_index('comparison_id').structure_group.sort_index())
        assert groups == saved
        hashes[str(data_path.relative_to(ROOT))] = hashlib.sha256(data_path.read_bytes()).hexdigest()
        checked = set()
        for group in groups:
            members = data[data.structure_group.eq(group['structure_group'])]
            core = Chem.MolFromSmiles(group['core_smiles']) if group['core_smiles'] else None
            coverage = sum(Chem.MolFromSmiles(s).HasSubstructMatch(core, useChirality=True)
                           for s in members.canonical_smiles) if core else 0
            exact = int(members.scaffold_smiles.fillna('').eq(group['core_smiles']).sum()) if core else 0
            assert coverage == group['core_coverage_count']
            assert exact == group['core_exact_scaffold_count']
            selected, _, audit = group_comparison(data, group)
            assert len(selected) == min(12, len(members))
            points = [np.asarray(a['core_draw_coordinates']) for a in audit if a['aligned']]
            assert all(np.allclose(points[0], p, atol=1e-5) for p in points)
            groups_checked.append(dict(target=target, group=group['structure_group'], size=len(members),
                coverage_count=coverage, coverage=coverage/len(members), exact_count=exact,
                selected=len(selected), not_displayed=len(members)-len(selected), verified=True))
            for a in audit:
                assert a['structure_preserved']
                assert a['comparison_id'] not in checked
                checked.add(a['comparison_id'])
                instances.append(dict(target=target, group=group['structure_group'], **a))
            (OUT / f'{target}_{group["structure_group"]}_atom_matches.json').write_text(json.dumps(audit, indent=2))
        trial = pd.DataFrame([a for a in instances if a['target'] == target])
        populations.append(dict(target=target, entries=len(data), group_summaries=len(groups),
            summaries_checked=len(groups), summaries_passed=len(groups), eligible_entries=len(data),
            depiction_instances=len(trial), unique_depicted_entries=len(checked), passed=int(trial.structure_preserved.sum()),
            core_matched=int(trial.core_matched.sum()), aligned=int(trial.aligned.sum()),
            unmatched=int((~trial.core_matched).sum()),
            stereo_fallback=int(trial.alignment_status.eq('stereo_safe_unaligned').sum()),
            not_displayed=len(data)-len(checked), reason='Default group preview limit of 12; full membership remains exportable'))
        print(populations[-1], flush=True)
    pd.DataFrame(groups_checked).to_csv(OUT / 'all_group_coverage.csv', index=False)
    pd.DataFrame(populations).to_csv(OUT / 'depiction_denominators.csv', index=False)
    pd.DataFrame(instances).drop(columns=['core_draw_coordinates']).to_csv(OUT / 'depiction_instances.csv', index=False)
    report = dict(environment=environment_versions(), input_sha256=hashes, per_target=populations,
        summaries_checked=len(groups_checked), depiction_instances=len(instances),
        all_passed=True, population='One default preview per group, at most 12 entries; case redraws are separate',
        code_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                     for p in sorted((ROOT/'src/chemgridmap').rglob('*.py'))})
    (OUT/'validation.json').write_text(json.dumps(report, indent=2))
    ca=pd.read_csv(SOURCE/'chembl205/coordinates_seed_42.csv')
    query='CHEMBL81688'
    checks=[]
    for seed in [7,21,42,84,126]:
        coords=ca if seed==42 else pd.read_csv(ROOT/f'paper/output/revision_v6/chembl205/coordinates_seed_{seed}.csv')
        coords=coords.set_index('molecule_identity_key').loc[ca.molecule_identity_key].reset_index()
        assert coords.canonical_smiles.tolist()==ca.canonical_smiles.tolist()
        coords['comparison_id']=ca.comparison_id.to_numpy()
        idx=int(np.flatnonzero(coords.comparison_id.eq(query))[0])
        neighbors,stats=neighborhood_evidence(coords,morgan_fingerprints(coords.canonical_smiles.tolist()),idx)
        neighbors.to_csv(OUT/f'query_neighbors_seed_{seed}.csv',index=False)
        counts=neighbors.neighbor_relation.value_counts()
        common=int(counts.get('shared',0));nf=common+int(counts.get('representation_only',0));ng=common+int(counts.get('grid_only',0))
        checks.append(dict(seed=seed,query=query,fingerprint_neighbors=nf,grid_neighbors=ng,shared=common,recall=common/nf,precision=common/ng,jaccard=common/(nf+ng-common)))
    pd.DataFrame(checks).to_csv(OUT/'fixed_query_seed_summary.csv',index=False)


if __name__ == '__main__':
    main()
