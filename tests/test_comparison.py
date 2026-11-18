import json
import zipfile

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem

from chemgridmap import build_grid_map
from chemgridmap.app import AppSession
from chemgridmap.comparison import (aligned_comparison, structural_groups, group_layout_metrics,
                                     group_comparison, main)


def test_browser_default_group_threshold_is_a_valid_number_input():
    from decimal import Decimal
    from html.parser import HTMLParser
    from pathlib import Path
    import chemgridmap
    found = []
    class Inputs(HTMLParser):
        def handle_starttag(self, tag, attrs):
            attrs = dict(attrs)
            if tag == 'input' and attrs.get('id') == 'group-similarity':
                found.append(attrs)
    Inputs().feed((Path(chemgridmap.__file__).parent/'web/index.html').read_text(encoding="utf-8"))
    item = found[0]
    assert (Decimal(item['value'])-Decimal(item['min'])) % Decimal(item['step']) == 0


def molecules():
    return pd.DataFrame({'canonical_smiles': ['Cc1ccccc1', 'CCc1ccccc1', 'Oc1ccccc1',
                                             'C[C@H](O)c1ccccc1', 'CCO', 'C1CCCCC1'],
                         'molecule_id': ['M1', 'M2', 'M3', 'M4', 'M5', 'M6'],
                         'activity_class': ['active'] * 6, 'activity_pchembl': np.arange(6) + 5.})


def test_groups_independent_of_input_order_and_activity_and_layout():
    data = molecules()
    first, groups, protocol = structural_groups(data, similarity=.4)
    shuffled = data.sample(frac=1, random_state=42).copy()
    shuffled['activity_pchembl'] = -10
    shuffled['projection_x'] = 1e6
    second, groups2, _ = structural_groups(shuffled, similarity=.4)
    pd.testing.assert_series_equal(first.set_index('molecule_id').structure_group.sort_index(),
                                   second.set_index('molecule_id').structure_group.sort_index())
    assert groups == groups2
    assert protocol['distance'] == '1 - binary Morgan Tanimoto'
    assert all(g['minimum_center_similarity'] >= .4 - 1e-12 for g in groups)


def test_external_groups_need_source_and_report_real_core_coverage():
    data = molecules().assign(cluster_id='source-A')
    with pytest.raises(ValueError, match='source'):
        structural_groups(data, cluster_col='cluster_id')
    frame, groups, protocol = structural_groups(data, cluster_col='cluster_id', cluster_source='supplied test set')
    assert groups[0]['core_coverage_count'] == 4
    assert groups[0]['core_exact_scaffold_count'] == 4
    assert groups[0]['size'] == 6
    assert groups[0]['acyclic_members'] == 1
    assert protocol['method'] == 'external'
    selected, _, audit = group_comparison(frame, groups[0])
    assert len(selected) == 6 and sum(a['aligned'] for a in audit) == 4


def test_core_alignment_highlighting_and_stereochemistry():
    data = molecules().iloc[:4]
    svg, audit = aligned_comparison(data, 'c1ccccc1')
    assert '<image' not in svg and len(audit) == 4
    reference = np.array(audit[0]['core_draw_coordinates'])
    for a, smiles in zip(audit, data.canonical_smiles):
        assert a['structure_preserved'] and a['molblock_roundtrip_preserves_structure']
        assert a['aligned'] and len(a['core_atom_indices']) == 6
        assert not set(a['core_atom_indices']) & set(a['highlight_atom_indices'])
        assert len(a['core_atom_indices']) + len(a['highlight_atom_indices']) == Chem.MolFromSmiles(smiles).GetNumAtoms()
        np.testing.assert_allclose(reference, a['core_draw_coordinates'], atol=1e-5)


def test_no_scaffold_is_not_invented():
    data = molecules().iloc[[4]]
    _, groups, _ = structural_groups(data)
    assert groups[0]['core_smiles'] == '' and groups[0]['core_coverage_count'] == 0
    _, audit = aligned_comparison(data, '')
    assert not audit[0]['aligned'] and not audit[0]['highlight_atom_indices']


def test_incompatible_ez_core_falls_back_without_changing_structure():
    from rdkit.Chem.Scaffolds import MurckoScaffold
    smiles = r'C/C(=C(\F)C(=O)Nc1ccc(-c2ccccc2S(N)(=O)=O)cc1Br)c1cccc(C(=N)N)c1'
    core = MurckoScaffold.MurckoScaffoldSmiles(smiles=smiles, includeChirality=True)
    _, audits = aligned_comparison(pd.DataFrame({'canonical_smiles': [smiles]}), core)
    assert audits[0]['structure_preserved'] and audits[0]['molblock_roundtrip_preserves_structure']
    assert audits[0]['core_matched']


def test_group_metrics_include_ties_and_singletons():
    data = pd.DataFrame({'structure_group': ['A', 'A', 'B', 'C']})
    coords = np.array([[0, 0], [1, 0], [0, 1], [-1, 0]])
    metrics = group_layout_metrics(data, coords, k=1).set_index('structure_group')
    assert metrics.loc['A', 'same_group_fraction'] == pytest.approx((1/3 + 1) / 2)
    assert metrics.loc['A', 'random_label_expectation'] == pytest.approx(1/3)
    assert pd.isna(metrics.loc['C', 'enrichment'])


def test_pipeline_grouping_does_not_change_coordinates(tmp_path):
    data = molecules()
    with_groups = build_grid_map(data, tmp_path, smiles_col='canonical_smiles', output_formats=['svg'])
    without = build_grid_map(data, smiles_col='canonical_smiles', grouping='none')
    np.testing.assert_array_equal(with_groups.grid, without.grid)
    np.testing.assert_array_equal(with_groups.projection, without.projection)
    assert 'structure_groups' in with_groups.output_files
    assert main([with_groups.output_files['coordinates'], '-o', str(tmp_path/'comparison')]) == 0
    assert (tmp_path/'comparison'/'atom_matches.json').exists()


def test_app_exports_core_comparison_evidence(tmp_path):
    app = AppSession(tmp_path)
    data = molecules().rename(columns={'canonical_smiles': 'smiles'})
    app.load(data.to_csv(index=False).encode(), 'structures.csv')
    app.start({'input_format': 'molecule', 'value_col': 'activity_pchembl', 'output_formats': ['svg']}, background=False)
    assert app.job['state'] == 'done', app.job
    group = app.result.structure_groups[0]
    detail = app.compare_group(group['structure_group'])
    assert detail['selected_count'] <= group['size']
    with zipfile.ZipFile(app.files[detail['archive']]) as archive:
        audit = json.loads(archive.read('atom_matches.json'))
        assert all(row['structure_preserved'] for row in audit)
        assert 'members.csv' in archive.namelist()
