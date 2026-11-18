"""U10: RECAP fragments as a core-candidate source and the generic Murcko scaffold as supplementary information."""
from pathlib import Path

import pandas as pd
import pytest
from rdkit import Chem

from chemgridmap import cores
from chemgridmap.comparison import summarize_group_core

from test_u9a_core_rules import SMILES, hashes, make_session

# Sulfonamide aryl esters: RECAP cleaves the ester and amide bonds and keeps the
# sulfonyl-aryl ring fragment that typed Murcko scaffolds reduce to benzene.
ESTERS = ['CCOC(=O)c1ccc(S(N)(=O)=O)cc1', 'CCCOC(=O)c1ccc(S(N)(=O)=O)cc1',
          'CC(C)OC(=O)c1ccc(S(N)(=O)=O)cc1', 'CCCCOC(=O)c1ccc(S(N)(=O)=O)cc1']


def mols(smiles):
    return [Chem.MolFromSmiles(s) for s in smiles]


def test_recap_candidates_are_ring_containing_dummy_free_and_order_independent():
    first, skipped = cores.recap_candidate_cores(mols(ESTERS))
    second, _ = cores.recap_candidate_cores(list(reversed(mols(ESTERS))))
    assert first == second == sorted(set(first)) and first
    assert skipped == 0
    for smiles in first:
        mol = Chem.MolFromSmiles(smiles)
        assert mol.GetRingInfo().NumRings() > 0
        assert all(atom.GetAtomicNum() != 0 for atom in mol.GetAtoms())
    assert cores.recap_candidate_cores(mols(['CCOC(=O)CC']))[0] == []


def test_recap_skips_only_members_above_heavy_atom_limit(monkeypatch):
    monkeypatch.setattr(cores, 'RECAP_MAX_HEAVY_ATOMS', 12)
    small, large = mols(['CCOC(=O)c1ccccc1', 'CCCCCCCCOC(=O)c1ccc(S(N)(=O)=O)cc1'])
    assert small.GetNumHeavyAtoms() <= 12 < large.GetNumHeavyAtoms()
    _, skipped = cores.recap_candidate_cores([small, large])
    assert skipped == 1


def test_recap_rule_enlarges_core_and_reports_its_own_columns():
    group = mols(ESTERS)
    plain = cores.select_core(group)
    recap = cores.select_core(group, include_recap=True)
    assert recap['core_rule'] == 'coverage_recap' and recap['target_met']
    assert recap['core_heavy_atoms'] > plain['core_heavy_atoms']
    assert recap['coverage_count'] == len(group)
    assert recap['recap_candidates'] == recap['candidates'] - plain['candidates'] > 0
    assert recap['recap_skipped_members'] == 0
    assert 'brics_candidates' not in recap
    assert Chem.MolFromSmiles(recap['core_smiles_attached']) is not None
    with pytest.raises(ValueError, match='BRICS or RECAP'):
        cores.select_core(group, include_brics=True, include_recap=True)


def test_existing_rules_do_not_change_when_recap_exists():
    group = mols(ESTERS + SMILES)
    ids = [f'M{i}' for i in range(len(group))]
    for rule in ('coverage', 'coverage_brics', 'most_frequent'):
        summary = summarize_group_core(group, ids, core_rule=rule)
        assert not any(key.startswith('recap_') for key in summary)
    recap = summarize_group_core(group, ids, core_rule='coverage_recap')
    assert recap['core_rule'] == 'coverage_recap' and 'recap_candidates' in recap
    assert recap['core_coverage_target'] == .8 and recap['core_target_met'] in (True, False)


def test_generic_scaffold_summary_counts_shared_shapes():
    # Benzene and pyridine analogues share one generic scaffold but two typed scaffolds.
    group = mols(['CCc1ccccc1', 'CCc1ccncc1', 'CCc1ccccn1', 'CCC1CCCC1'])
    summary = cores.generic_scaffold_summary(group)
    assert summary['generic_scaffold'] == 'C1CCCCC1'
    assert summary['generic_scaffold_count'] == 3 and summary['size'] == 4
    assert summary['generic_scaffold_fraction'] == pytest.approx(.75)
    assert summary['distinct_generic_scaffolds'] == 2
    assert summary['generic_scaffold_heavy_atoms'] == 6
    empty = cores.generic_scaffold_summary(mols(['CCO', 'CCN']))
    assert empty['generic_scaffold'] == '' and empty['generic_scaffold_count'] == 0


def test_switch_to_recap_writes_only_new_core_files(tmp_path):
    session = make_session(tmp_path)
    directory = Path(session.job['output_dir'])
    before = hashes(directory)
    payload = session.set_core_rule('coverage_recap')
    assert payload['rule'] == 'coverage_recap' and not payload['cached']
    after = hashes(directory)
    assert {name: after[name] for name in before} == before
    assert set(after) - set(before) == {'core_rules/map_cores_coverage_recap.csv',
                                       'core_rules/map_cores_coverage_recap.json'}
    table = pd.read_csv(directory / 'core_rules/map_cores_coverage_recap.csv')
    assert {'recap_candidates', 'recap_skipped_members'} <= set(table)
    assert 'brics_candidates' not in table


def test_group_card_and_comparison_include_generic_scaffold(tmp_path):
    session = make_session(tmp_path)
    group_id = session.groups()[0]['structure_group']
    detail = session.group_detail(group_id)
    generic = detail['generic']
    assert generic['size'] == len(SMILES)
    assert generic['generic_scaffold'] == 'C1CCCCC1'
    assert generic['generic_svg'].startswith('<svg')
    comparison = session.group_cores(group_id)
    assert [core['core_rule'] for core in comparison['cores']] == [
        'most_frequent', 'coverage', 'coverage_brics', 'coverage_recap', 'fmcs']
    assert comparison['generic'] == generic
    table = pd.read_csv(session.files[comparison['file']])
    last = table.iloc[-1]
    assert last.core_rule == 'generic_murcko_supplementary'
    assert last.core_coverage_count == generic['generic_scaffold_count']
    assert pd.isna(last.core_target_met) and pd.isna(last.core_coverage_target)
