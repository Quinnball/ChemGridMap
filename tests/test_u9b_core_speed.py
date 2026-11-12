"""The faster evaluation order preserves the U9a exhaustive selection result."""
import math
from collections import Counter

import pytest
from rdkit import Chem

from chemgridmap import cores


def exhaustive_selection(mols, target, include_brics):
    """U9a's independent all-candidates/all-members ranking oracle."""
    candidates = cores.candidate_cores(mols)
    extra = {}
    rule = 'coverage_aware'
    if include_brics:
        fragments, guarded = cores.brics_candidate_cores(mols)
        extra = {'brics_candidates': len(set(fragments) - set(candidates)),
                 'brics_single_pass_members': guarded}
        candidates = sorted(set(candidates).union(fragments))
        rule = 'coverage_brics'
    scored = []
    for smiles in candidates:
        pattern = Chem.MolFromSmiles(smiles)
        matched = [bool(m.HasSubstructMatch(pattern, useChirality=True)) for m in mols]
        fraction = sum(matched) / len(mols)
        meets = fraction >= target - 1e-12
        heavy = pattern.GetNumHeavyAtoms()
        rank = (-heavy, -fraction) if meets else (-fraction, -heavy)
        scored.append((not meets, rank, smiles, matched, heavy))
    if not scored:
        return {'core_smiles': '', 'core_smiles_attached': '', 'core_smiles_attached_source': '',
                'core_rule': rule, 'coverage_target': target,
                'coverage_count': 0, 'coverage_fraction': 0.0, 'core_heavy_atoms': 0,
                'matched': [False] * len(mols), 'target_met': False, 'candidates': 0, **extra}
    missed, _, smiles, matched, heavy = min(scored, key=lambda item: (item[0], item[1], item[2]))
    attached, source = cores._core_smiles_with_attachments(Chem.MolFromSmiles(smiles), mols)
    return {'core_smiles': smiles, 'core_smiles_attached': attached,
            'core_smiles_attached_source': source, 'core_rule': rule, 'coverage_target': target,
            'coverage_count': int(sum(matched)), 'coverage_fraction': sum(matched) / len(mols),
            'core_heavy_atoms': int(heavy), 'matched': matched, 'target_met': not missed,
            'candidates': len(candidates), **extra}


def without_seconds(result):
    return {key: value for key, value in result.items() if key != 'seconds'}


@pytest.mark.parametrize('smiles', [
    [], ['CCO', 'CCCC'],
    ['CCNC(=O)c1ccc(O)cc1', 'CCCNC(=O)c1ccc(O)cc1', 'CCCCNC(=O)c1ccc(O)cc1'],
    ['c1ccccc1', 'Cc1ccccc1', 'c1ccncc1', 'C1CCCCC1'],
    ['C[C@H]1CC[C@@H](C)CC1', 'C[C@@H]1CC[C@@H](C)CC1', 'CC1CCCCC1'],
    ['c1ccc2ccccc2c1', 'c1ccc(-c2ccccc2)cc1', 'CC1CCCCC1', 'CCO'],
])
@pytest.mark.parametrize('include_brics', [False, True])
@pytest.mark.parametrize('target', [.01, .5, .75, .8, 1.0])
def test_matches_exhaustive_u9a_result(smiles, include_brics, target):
    mols = [Chem.MolFromSmiles(value) for value in smiles]
    expected = exhaustive_selection(mols, target, include_brics)
    actual = cores.select_core(mols, min_coverage=target, include_brics=include_brics)
    assert without_seconds(actual) == expected


def matrix_members(monkeypatch, matrix):
    calls = Counter()
    class Member:
        def __init__(self, index):
            self.index = index

        def HasSubstructMatch(self, pattern, *, useChirality):
            assert useChirality is True
            smiles = Chem.MolToSmiles(pattern)
            calls[smiles] += 1
            return matrix[smiles][self.index]

    monkeypatch.setattr(cores, 'candidate_cores', lambda mols: sorted(matrix))
    monkeypatch.setattr(cores, '_core_smiles_with_attachments', lambda *args: ('', ''))
    return [Member(index) for index in range(len(next(iter(matrix.values()))))], calls


def test_early_rejection_finishes_winning_size_layer_and_skips_smaller(monkeypatch):
    matrix = {
        'C': [True] * 10,
        'CC': [True] * 10,
        'CCC': [True] * 8 + [False] * 2,
        'CCN': [True] * 9 + [False],
        'CCO': [True] * 9 + [False],
        'CCCC': [False] * 10,
    }
    mols, calls = matrix_members(monkeypatch, matrix)
    actual = cores.select_core(mols)
    assert actual['core_smiles'] == 'CCN'  # Coverage first, canonical SMILES second.
    assert actual['matched'] == matrix['CCN']
    assert actual['candidates'] == len(matrix)
    assert calls == {'CCCC': 3, 'CCC': 10, 'CCN': 10, 'CCO': 10}
    assert without_seconds(actual) == exhaustive_selection(mols, .8, False)


def test_fallback_recomputes_all_candidates_to_full_coverage(monkeypatch):
    matrix = {'C': [False, True, True, True], 'CC': [False, True, True, True],
              'CN': [False, True, True, True], 'CCC': [False, False, True, True]}
    mols, calls = matrix_members(monkeypatch, matrix)
    actual = cores.select_core(mols, min_coverage=1)
    assert actual['core_smiles'] == 'CC'
    assert not actual['target_met']
    assert calls == {smiles: 5 for smiles in matrix}
    assert without_seconds(actual) == exhaustive_selection(mols, 1, False)


@pytest.mark.parametrize('target', [
    2 / 3, 2 / 3 + 1e-12,
    math.nextafter(2 / 3 + 1e-12, 0), math.nextafter(2 / 3 + 1e-12, 1),
])
def test_target_tolerance_uses_identical_float_comparison(monkeypatch, target):
    mols, _ = matrix_members(monkeypatch, {'C': [True, True, True], 'CC': [False, True, True]})
    expected = exhaustive_selection(mols, target, False)
    assert without_seconds(cores.select_core(mols, min_coverage=target)) == expected
