"""Frozen-layout structural comparisons. No activity-based case selection."""

from pathlib import Path
import hashlib
import json
import time

import numpy as np
import pandas as pd
from rdkit import Chem

from chemgridmap.comparison import (structural_groups, save_group_outputs, group_layout_metrics,
                                   aligned_comparison, group_comparison, fingerprint_distances)
from chemgridmap.neighborhoods import neighborhood_evidence
from chemgridmap.representations import morgan_fingerprints
from chemgridmap.provenance import environment_versions

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'paper/output/structure_comparison'
CURRENT = ROOT / 'paper/output/current'
ARCHIVE = ROOT / 'paper/output/revision_v6'
CASE_RULE = ('CA II groups with >=6 members, selected scaffold >=9 heavy atoms and >=80% substructure coverage; '
             'rank by covered-member count descending, core heavy atoms descending, then group ID. '
             'Within covered members: representative if covered, otherwise canonical-first; '
             'farthest-first Morgan Tanimoto diversity, ties by canonical SMILES then identifier; max 9.')


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False), encoding='utf-8')


def choose_case(data, groups):
    eligible = [g for g in groups if g['size'] >= 6 and g['core_coverage_fraction'] >= .8
                and Chem.MolFromSmiles(g['core_smiles']).GetNumHeavyAtoms() >= 9]
    if not eligible:
        raise ValueError('No group satisfies the prespecified structure-only case rule.')
    selected = sorted(eligible, key=lambda g: (-g['core_coverage_count'], -Chem.MolFromSmiles(g['core_smiles']).GetNumHeavyAtoms(), g['structure_group']))[0]
    candidates = data[data.comparison_id.isin(selected['core_member_ids'])].sort_values(['canonical_smiles', 'comparison_id']).reset_index(drop=True)
    _, distances = fingerprint_distances(candidates.canonical_smiles)
    rep = np.flatnonzero(candidates.comparison_id.eq(selected['representative_id']))
    chosen = [int(rep[0])] if len(rep) else [0]
    while len(chosen) < min(9, len(candidates)):
        remaining = [i for i in range(len(candidates)) if i not in chosen]
        chosen.append(max(remaining, key=lambda i: (float(distances[i, chosen].min()), -i)))
    return selected, candidates.iloc[chosen].copy()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    manifest = {'case_selection_rule': CASE_RULE, 'similarity': .6, 'threshold_sensitivity': [.5, .6, .7],
                'k': 10, 'tie_policy': 'all within kth distance plus 1e-12',
                'aggregation': 'equal weight per nonsingleton group; no biological replicate interpretation',
                'inputs': {}, 'environment': environment_versions()}
    # Record the rule before computing or viewing the case.
    dump(OUT/'protocol.json', manifest)
    metric_tables, sensitivity, datasets, checks = [], [], [], []
    for target, seeds in [('chembl205', [7, 21, 42, 84, 126]), ('chembl204', [7, 21, 42]), ('chembl240', [7, 21, 42])]:
        path = CURRENT/target/'coordinates_seed_42.csv'
        manifest['inputs'][str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
        base = pd.read_csv(path)
        data, groups, protocol = structural_groups(base)
        output = OUT/target
        save_group_outputs(data, groups, protocol, output)
        data.to_csv(output/'coordinates_seed_42.csv', index=False)
        for threshold in [.5, .6, .7]:
            if threshold == .6:
                trial, trial_groups = data, groups
            else:
                trial, trial_groups, _ = structural_groups(base, similarity=threshold)
            sizes = [g['size'] for g in trial_groups]
            sensitivity.append({'target': target, 'similarity': threshold, 'groups': len(sizes),
                                'singletons': sum(s == 1 for s in sizes), 'largest_group': max(sizes),
                                'molecules_in_nonsingleton_groups': sum(s for s in sizes if s > 1)})
        group_by_id = data.set_index('molecule_identity_key')[['structure_group', 'comparison_id']]
        for seed in seeds:
            source = (CURRENT if seed == 42 else ARCHIVE)/target/f'coordinates_seed_{seed}.csv'
            manifest['inputs'][str(source.relative_to(ROOT))] = hashlib.sha256(source.read_bytes()).hexdigest()
            coords = pd.read_csv(source)
            assert coords.molecule_identity_key.is_unique
            assert set(coords.molecule_identity_key) == set(data.molecule_identity_key)
            ref = data.set_index('molecule_identity_key').loc[coords.molecule_identity_key]
            assert ref.canonical_smiles.tolist() == coords.canonical_smiles.tolist()
            coords['structure_group'] = group_by_id.loc[coords.molecule_identity_key, 'structure_group'].to_numpy()
            for space, columns in [('projection', ['projection_x', 'projection_y']), ('grid', ['grid_x', 'grid_y'])]:
                table = group_layout_metrics(coords, coords[columns], space=space)
                table.insert(0, 'seed', seed)
                table.insert(0, 'target', target)
                metric_tables.append(table)
        for group in groups:
            selected, _, audit = group_comparison(data, group)
            dump(output/'atom_matches'/f"{group['structure_group']}.json", audit)
            coverage = sum(Chem.MolFromSmiles(str(s)).HasSubstructMatch(Chem.MolFromSmiles(group['core_smiles']), useChirality=True)
                           for s in data.loc[data.structure_group.eq(group['structure_group']), 'canonical_smiles']) if group['core_smiles'] else 0
            assert coverage == group['core_coverage_count']
            core_points = [np.array(a['core_draw_coordinates']) for a in audit if a['aligned']]
            alignment_error = max((float(np.max(np.abs(p-core_points[0]))) for p in core_points), default=0.)
            checks.append({'target': target, 'group': group['structure_group'], 'size': group['size'],
                           'coverage_verified': coverage == group['core_coverage_count'],
                           'selected_members': len(selected), 'matched_selected_members': len(core_points),
                           'stereo_safe_fallbacks': sum(a['alignment_status'] == 'stereo_safe_unaligned' for a in audit),
                           'all_depictions_preserve_structure': all(a['structure_preserved'] for a in audit),
                           'maximum_core_pixel_difference': alignment_error})
        datasets.append({'target': target, 'molecules': len(data), 'groups': len(groups),
                         'nonsingleton_groups': sum(g['size'] > 1 for g in groups),
                         'molecules_in_nonsingleton_groups': sum(g['size'] for g in groups if g['size'] > 1),
                         'majority_core_nonsingleton_groups': sum(g['size'] > 1 and g['core_coverage_fraction'] > .5 for g in groups)})
        if target == 'chembl205':
            case, members = choose_case(data, groups)
            svg, audit = aligned_comparison(members, case['core_smiles'])
            (output/'case_comparison.svg').write_text(svg, encoding='utf-8')
            members.to_csv(output/'case_members.csv', index=False)
            dump(output/'case_atom_matches.json', audit)
            dump(output/'case.json', {'group': case, 'selection_rule': CASE_RULE, 'shown_members': members.comparison_id.tolist()})
            # Inspect both the representative and the familiar evidence-audit molecule.
            fps = morgan_fingerprints(data.canonical_smiles.tolist())
            for label, identity in [('case_center', case['representative_id']), ('acetazolamide', 'CHEMBL20')]:
                idx = int(np.flatnonzero(data.comparison_id.eq(identity))[0])
                neighbors, summary = neighborhood_evidence(data, fps, idx)
                neighbors.to_csv(output/f'{label}_neighbors.csv', index=False)
                dump(output/f'{label}_neighbors.json', summary)
        print(target, datasets[-1], flush=True)
    metrics = pd.concat(metric_tables, ignore_index=True)
    metrics.to_csv(OUT/'group_layout_metrics.csv', index=False)
    macro = metrics[metrics['size'] > 1].groupby(['target', 'seed', 'space'])[['same_group_fraction', 'random_label_expectation', 'excess_fraction']].mean().reset_index()
    macro.to_csv(OUT/'group_macro_metrics.csv', index=False)
    pd.DataFrame(datasets).to_csv(OUT/'datasets.csv', index=False)
    pd.DataFrame(sensitivity).to_csv(OUT/'threshold_sensitivity.csv', index=False)
    pd.DataFrame(checks).to_csv(OUT/'structure_checks.csv', index=False)
    assert all(c['all_depictions_preserve_structure'] and c['coverage_verified'] and c['maximum_core_pixel_difference'] < 1e-5 for c in checks)
    manifest.update(seconds=time.perf_counter()-started, groups_checked=len(checks),
                    depictions_checked=sum(c['selected_members'] for c in checks),
                    aligned_depictions_checked=sum(c['matched_selected_members'] for c in checks),
                    stereo_safe_fallbacks=sum(c['stereo_safe_fallbacks'] for c in checks),
                    all_structure_checks_passed=True)
    dump(OUT/'validation.json', manifest)
    print(macro.groupby(['target', 'space']).mean(numeric_only=True).to_string(), flush=True)


if __name__ == '__main__':
    main()
