"""BRICS skips invalid aromatic fragments quietly without changing candidates."""

import pytest
from rdkit import Chem

from chemgridmap import comparison, cores


NOISY_SMILES = ["CCn1cccc1", "CCn1ccnc1", "CCn1cccn1", "CCn1ccc2ccccc21"]


def test_real_brics_fragment_fails_sanitization_and_logging_is_scoped(capfd):
    molecule = Chem.MolFromSmiles(NOISY_SMILES[0])
    fragments = cores.BRICS.BRICSDecompose(
        molecule, keepNonLeafNodes=True, returnMols=True)
    fragment = next(fragment for fragment in fragments
                    if Chem.MolToSmiles(fragment) == "[9*]n1cccc1")
    clean = Chem.RWMol(fragment)
    dummy = next(atom.GetIdx() for atom in clean.GetAtoms() if atom.GetAtomicNum() == 0)
    clean.RemoveAtom(dummy)
    with pytest.raises(Chem.KekulizeException):
        Chem.SanitizeMol(clean.GetMol())
    assert "Can't kekulize mol" in capfd.readouterr().err

    assert cores.brics_candidate_cores([molecule]) == ([NOISY_SMILES[0]], 0)
    assert capfd.readouterr().err == ""

    # Leaving the candidate-generation scope restores normal RDKit diagnostics.
    with pytest.raises(Chem.KekulizeException):
        Chem.SanitizeMol(clean.GetMol())
    assert "Can't kekulize mol" in capfd.readouterr().err


@pytest.mark.parametrize("smiles", NOISY_SMILES)
def test_noisy_brics_candidates_match_base_without_stderr(smiles, capfd):
    # Exact candidates and guard count verified against BASE cc12cac.
    assert cores.brics_candidate_cores([Chem.MolFromSmiles(smiles)]) == ([smiles], 0)
    assert capfd.readouterr().err == ""


def test_serial_and_real_workers_keep_results_and_emit_no_brics_logs(monkeypatch, capfd):
    monkeypatch.setattr(comparison, "CORE_PARALLEL_MIN_GROUPS", 2)
    groups = [([smiles], [f"M{index}"]) for index, smiles in enumerate(NOISY_SMILES)]
    # Start fresh workers so their stderr belongs to this capture scope.
    comparison.shutdown_core_workers()
    try:
        serial = comparison.summarize_group_cores(groups, core_rule="coverage_brics")
        assert capfd.readouterr().err == ""
        parallel = comparison.summarize_group_cores(groups, core_rule="coverage_brics", jobs=2)
        assert comparison._core_executor is not None
        assert len(comparison._core_executor._processes) == 2
    finally:
        comparison.shutdown_core_workers()
    assert capfd.readouterr().err == ""
    for actual, expected in zip(parallel, serial):
        assert actual[0] == expected[0]
        assert {key: value for key, value in actual[1].items() if key != "seconds"} == {
            key: value for key, value in expected[1].items() if key != "seconds"}
