"""Coverage-aware selection of a representative core for a structure group.

The most frequent Bemis-Murcko scaffold can match only a minority of a
fingerprint group. Here the candidates are every member scaffold plus the
typed ring-system fragments of the RDKit scaffold network; the selected core
is the largest candidate (heavy atoms) that matches at least ``min_coverage``
of the members. If no candidate reaches that coverage, the candidate with the
highest coverage is used and the summary reports that the target was missed.
"""
from __future__ import annotations

import time
from typing import Dict, Sequence

from rdkit import Chem, rdBase
from rdkit.Chem import rdFMCS
from rdkit.Chem.Scaffolds import MurckoScaffold, rdScaffoldNetwork


def _network_parameters():
    params = rdScaffoldNetwork.ScaffoldNetworkParams()
    params.includeGenericScaffolds = False
    params.includeGenericBondScaffolds = False
    params.includeScaffoldsWithAttachments = False
    params.includeScaffoldsWithoutAttachments = True
    return params


def candidate_cores(mols: Sequence[Chem.Mol]) -> list:
    """Canonical SMILES of member scaffolds and their typed scaffold-network fragments."""
    found = set()
    for mol in mols:
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=True)
        if scaffold:
            found.add(scaffold)
    try:
        network = rdScaffoldNetwork.CreateScaffoldNetwork(list(mols), _network_parameters())
        nodes = list(network.nodes)
    except Exception:  # pragma: no cover - network generation is optional enrichment
        nodes = []
    for smiles in nodes:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None or not mol.GetRingInfo().NumRings():
            continue
        # The network also contains the input molecules; keep only ring scaffolds without side chains.
        canonical = Chem.MolToSmiles(mol)
        if MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=True) == canonical:
            found.add(canonical)
    return sorted(found)


def _coverage(pattern: Chem.Mol, mols: Sequence[Chem.Mol]) -> list:
    return [bool(m.HasSubstructMatch(pattern, useChirality=True)) for m in mols]


def _attached_candidate(pattern: Chem.Mol, member: Chem.Mol, match: tuple, *,
                        kekule: bool) -> str:
    source = Chem.Mol(member)
    if kekule:
        Chem.Kekulize(source, clearAromaticFlags=True)
    drawing = Chem.RWMol()
    for pattern_atom, member_index in zip(pattern.GetAtoms(), match):
        member_atom = source.GetAtomWithIdx(member_index)
        atom = Chem.Atom(pattern_atom.GetAtomicNum())
        atom.SetFormalCharge(member_atom.GetFormalCharge())
        if not kekule:
            atom.SetIsAromatic(any(
                bond.GetBondType() == Chem.BondType.AROMATIC
                for bond in pattern_atom.GetBonds()
            ))
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(member_atom.GetTotalNumHs())
        drawing.AddAtom(atom)
    for bond in pattern.GetBonds():
        if kekule:
            member_bond = source.GetBondBetweenAtoms(
                match[bond.GetBeginAtomIdx()], match[bond.GetEndAtomIdx()]
            )
            bond_type = member_bond.GetBondType()
        else:
            bond_type = bond.GetBondType()
        drawing.AddBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond_type)
        drawing.GetBondBetweenAtoms(
            bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        ).SetIsAromatic(bond_type == Chem.BondType.AROMATIC and not kekule)

    member_atoms = set(match)
    member_to_core = {member_index: core_index for core_index, member_index in enumerate(match)}
    handled_member_bonds = set()

    def add_dummy(core_index, bond_type):
        dummy_index = drawing.AddAtom(Chem.Atom(0))
        drawing.GetAtomWithIdx(dummy_index).SetNoImplicit(True)
        drawing.AddBond(core_index, dummy_index, bond_type)
        drawing.GetBondBetweenAtoms(core_index, dummy_index).SetIsAromatic(
            bond_type == Chem.BondType.AROMATIC and not kekule
        )

    for core_index, member_index in enumerate(match):
        for member_bond in source.GetAtomWithIdx(member_index).GetBonds():
            other_index = member_bond.GetOtherAtomIdx(member_index)
            bond_type = member_bond.GetBondType()
            if other_index not in member_atoms:
                add_dummy(core_index, bond_type)
                continue
            other_core_index = member_to_core[other_index]
            if drawing.GetBondBetweenAtoms(core_index, other_core_index) is None:
                if member_bond.GetIdx() in handled_member_bonds:
                    continue
                add_dummy(core_index, bond_type)
                add_dummy(other_core_index, bond_type)
                handled_member_bonds.add(member_bond.GetIdx())
    molecule = drawing.GetMol()
    Chem.SanitizeMol(molecule)
    return Chem.MolToSmiles(molecule, isomericSmiles=True)


def _core_smiles_with_attachments(pattern: Chem.Mol, mols: Sequence[Chem.Mol]) -> tuple[str, str]:
    """Encode a matched core and report its aromatic or Kekulé source."""
    with rdBase.BlockLogs():
        matches = [(member, member.GetSubstructMatch(pattern)) for member in mols]
        matches = [(member, match) for member, match in matches if match]
        for member, match in matches:
            try:
                return _attached_candidate(pattern, member, match, kekule=False), "matched_member"
            except Exception:
                continue
        for member, match in matches:
            try:
                return _attached_candidate(pattern, member, match, kekule=True), "matched_member_kekule"
            except Exception:
                continue
    return "", ""


def select_core(mols: Sequence[Chem.Mol], *, min_coverage: float = 0.8) -> Dict:
    """Largest candidate core matching at least ``min_coverage`` of the members.

    Candidates meeting the target are ranked by heavy atoms, then coverage;
    if none meets it, the highest-coverage candidate is used. Remaining ties
    are broken by canonical SMILES.
    """
    if not 0 < min_coverage <= 1:
        raise ValueError("min_coverage must lie in (0, 1].")
    started = time.perf_counter()
    candidates = candidate_cores(mols)
    scored = []
    for smiles in candidates:
        pattern = Chem.MolFromSmiles(smiles)
        matched = _coverage(pattern, mols)
        fraction = sum(matched) / len(mols)
        meets = fraction >= min_coverage - 1e-12
        heavy = pattern.GetNumHeavyAtoms()
        rank = (-heavy, -fraction) if meets else (-fraction, -heavy)
        scored.append((not meets, rank, smiles, matched, heavy))
    if not scored:
        return {"core_smiles": "", "core_smiles_attached": "", "core_smiles_attached_source": "",
                "core_rule": "coverage_aware", "coverage_target": min_coverage,
                "coverage_count": 0, "coverage_fraction": 0.0, "core_heavy_atoms": 0, "matched": [False] * len(mols),
                "target_met": False, "candidates": 0, "seconds": time.perf_counter() - started}
    missed, _, smiles, matched, heavy = min(scored, key=lambda item: (item[0], item[1], item[2]))
    pattern = Chem.MolFromSmiles(smiles)
    attached, attached_source = _core_smiles_with_attachments(pattern, mols)
    return {"core_smiles": smiles,
            "core_smiles_attached": attached,
            "core_smiles_attached_source": attached_source,
            "core_rule": "coverage_aware", "coverage_target": min_coverage,
            "coverage_count": int(sum(matched)), "coverage_fraction": sum(matched) / len(mols),
            "core_heavy_atoms": int(heavy), "matched": matched, "target_met": not missed,
            "candidates": len(candidates), "seconds": time.perf_counter() - started}


def mcs_core(mols: Sequence[Chem.Mol], *, min_coverage: float = 0.8, timeout: int = 10) -> Dict:
    """Maximum common substructure present in at least ``min_coverage`` of members (comparison method)."""
    started = time.perf_counter()
    if len(mols) < 2:
        return {"core_smiles": "", "core_smiles_attached": "", "core_smiles_attached_source": "",
                "core_smarts": "", "core_rule": "fmcs_threshold",
                "coverage_target": min_coverage, "coverage_count": 0, "coverage_fraction": 0.0,
                "core_heavy_atoms": 0, "timed_out": False, "core_smiles_source": "",
                "seconds": time.perf_counter() - started}
    result = rdFMCS.FindMCS(list(mols), threshold=min_coverage, ringMatchesRingOnly=True, completeRingsOnly=True,
                            atomCompare=rdFMCS.AtomCompare.CompareElements,
                            bondCompare=rdFMCS.BondCompare.CompareOrder,
                            timeout=int(timeout))
    pattern = Chem.MolFromSmarts(result.smartsString) if result.numAtoms else None
    matched = [
        bool(pattern is not None and m.HasSubstructMatch(pattern, useChirality=True))
        for m in mols
    ]
    core_smiles = ""
    core_smiles_source = ""
    for mol, is_match in zip(mols, matched):
        if not is_match:
            continue
        atom_indices = list(mol.GetSubstructMatch(pattern, useChirality=True))
        atom_set = set(atom_indices)
        bond_indices = [
            bond.GetIdx() for bond in mol.GetBonds()
            if bond.GetBeginAtomIdx() in atom_set and bond.GetEndAtomIdx() in atom_set
        ]
        candidate = Chem.MolFragmentToSmiles(
            mol, atomsToUse=atom_indices, bondsToUse=bond_indices, canonical=True
        )
        parsed = Chem.MolFromSmiles(candidate)
        if parsed is not None and parsed.GetNumHeavyAtoms() == result.numAtoms:
            core_smiles = candidate
            core_smiles_source = "matched_member"
            break
        try:
            kekule_mol = Chem.Mol(mol)
            Chem.Kekulize(kekule_mol, clearAromaticFlags=True)
        except Exception:
            break
        candidate = Chem.MolFragmentToSmiles(
            kekule_mol, atomsToUse=atom_indices, bondsToUse=bond_indices,
            canonical=True, kekuleSmiles=True,
        )
        parsed = Chem.MolFromSmiles(candidate)
        if parsed is not None and parsed.GetNumHeavyAtoms() == result.numAtoms:
            core_smiles = candidate
            core_smiles_source = "matched_member_kekule"
        break
    attached, attached_source = (
        _core_smiles_with_attachments(pattern, mols) if pattern is not None else ("", "")
    )
    return {"core_smiles": core_smiles,
            "core_smiles_attached": attached,
            "core_smiles_attached_source": attached_source,
            "core_smarts": result.smartsString, "core_rule": "fmcs_threshold", "coverage_target": min_coverage,
            "coverage_count": int(sum(matched)), "coverage_fraction": sum(matched) / len(mols),
            "core_heavy_atoms": int(result.numAtoms), "timed_out": bool(result.canceled),
            "core_smiles_source": core_smiles_source,
            "seconds": time.perf_counter() - started}
