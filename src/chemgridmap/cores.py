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
from rdkit.Chem import BRICS, Recap, rdFMCS
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


def brics_candidate_cores(mols: Sequence[Chem.Mol]) -> tuple[list, int]:
    """Return ring-containing BRICS fragments without dummy atoms and the guard count."""
    found, single_pass_members = set(), 0
    with rdBase.BlockLogs():
        for mol in mols:
            single_pass = len(list(BRICS.FindBRICSBonds(mol))) > 10
            single_pass_members += int(single_pass)
            fragments = BRICS.BRICSDecompose(
                mol, keepNonLeafNodes=True, returnMols=True, singlePass=single_pass)
            for fragment in fragments:
                clean = Chem.RWMol(fragment)
                for index in reversed([atom.GetIdx() for atom in clean.GetAtoms()
                                       if atom.GetAtomicNum() == 0]):
                    clean.RemoveAtom(index)
                clean = clean.GetMol()
                try:
                    Chem.SanitizeMol(clean)
                except Exception:
                    continue
                if clean.GetRingInfo().NumRings():
                    found.add(Chem.MolToSmiles(clean, isomericSmiles=True))
    return sorted(found), single_pass_members


def _ring_fragment_without_dummies(fragment: Chem.Mol) -> str:
    """Canonical SMILES of a fragment after removing dummy atoms, or "" if acyclic or invalid."""
    clean = Chem.RWMol(fragment)
    for index in reversed([atom.GetIdx() for atom in clean.GetAtoms() if atom.GetAtomicNum() == 0]):
        clean.RemoveAtom(index)
    clean = clean.GetMol()
    try:
        Chem.SanitizeMol(clean)
    except Exception:
        return ""
    if not clean.GetRingInfo().NumRings():
        return ""
    return Chem.MolToSmiles(clean, isomericSmiles=True)


RECAP_MAX_HEAVY_ATOMS = 80


def recap_candidate_cores(mols: Sequence[Chem.Mol]) -> tuple[list, int]:
    """Return ring-containing RECAP fragments without dummy atoms and the skipped-member count.

    Every node of each member's RECAP hierarchy is used (leaves and intermediate
    fragments), mirroring ``keepNonLeafNodes=True`` for BRICS. The full hierarchy
    can grow combinatorially, so members with more than ``RECAP_MAX_HEAVY_ATOMS``
    heavy atoms are not decomposed; they still count in coverage.
    """
    found, skipped = set(), 0
    with rdBase.BlockLogs():
        for mol in mols:
            if mol.GetNumHeavyAtoms() > RECAP_MAX_HEAVY_ATOMS:
                skipped += 1
                continue
            for node in Recap.RecapDecompose(mol).GetAllChildren().values():
                if node.mol is None:
                    continue
                smiles = _ring_fragment_without_dummies(node.mol)
                if smiles:
                    found.add(smiles)
    return sorted(found), skipped


def generic_scaffold_summary(mols: Sequence[Chem.Mol]) -> Dict:
    """Most frequent generic Murcko scaffold (all atoms carbon, all bonds single).

    Reported as supplementary information only: a generic scaffold shows a shared
    ring-and-linker shape but cannot be aligned to members as a typed core.
    Ties: member count descending, heavy atoms descending, then SMILES.
    """
    generic = []
    with rdBase.BlockLogs():
        for mol in mols:
            try:
                scaffold = MurckoScaffold.GetScaffoldForMol(mol)
                if not scaffold.GetNumAtoms():
                    generic.append("")
                    continue
                generic.append(Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(scaffold)))
            except Exception:
                generic.append("")
    counts: Dict[str, int] = {}
    for smiles in generic:
        if smiles:
            counts[smiles] = counts.get(smiles, 0) + 1
    if not counts:
        return {"generic_scaffold": "", "generic_scaffold_count": 0, "generic_scaffold_fraction": 0.0,
                "generic_scaffold_heavy_atoms": 0, "distinct_generic_scaffolds": 0, "size": len(mols)}
    best = min(counts, key=lambda s: (-counts[s], -Chem.MolFromSmiles(s).GetNumHeavyAtoms(), s))
    return {"generic_scaffold": best, "generic_scaffold_count": counts[best],
            "generic_scaffold_fraction": counts[best] / len(mols),
            "generic_scaffold_heavy_atoms": Chem.MolFromSmiles(best).GetNumHeavyAtoms(),
            "distinct_generic_scaffolds": len(counts), "size": len(mols)}


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


def select_core(mols: Sequence[Chem.Mol], *, min_coverage: float = 0.8,
                include_brics: bool = False, include_recap: bool = False) -> Dict:
    """Largest candidate core matching at least ``min_coverage`` of the members.

    Candidates meeting the target are ranked by heavy atoms, then coverage;
    if none meets it, the highest-coverage candidate is used. Remaining ties
    are broken by canonical SMILES.
    """
    if not 0 < min_coverage <= 1:
        raise ValueError("min_coverage must lie in (0, 1].")
    if include_brics and include_recap:
        raise ValueError("Choose BRICS or RECAP fragments, not both.")
    started = time.perf_counter()
    candidates = candidate_cores(mols)
    extra = {}
    rule = "coverage_aware"
    if include_brics:
        fragments, guarded = brics_candidate_cores(mols)
        extra = {"brics_candidates": len(set(fragments) - set(candidates)),
                 "brics_single_pass_members": guarded}
        candidates = sorted(set(candidates).union(fragments))
        rule = "coverage_brics"
    if include_recap:
        fragments, guarded = recap_candidate_cores(mols)
        extra = {"recap_candidates": len(set(fragments) - set(candidates)),
                 "recap_skipped_members": guarded}
        candidates = sorted(set(candidates).union(fragments))
        rule = "coverage_recap"
    patterns = [(smiles, Chem.MolFromSmiles(smiles)) for smiles in candidates]
    patterns.sort(key=lambda item: (-item[1].GetNumHeavyAtoms(), item[0]))
    scored = []
    winning_heavy = None
    for smiles, pattern in patterns:
        heavy = pattern.GetNumHeavyAtoms()
        if winning_heavy is not None and heavy < winning_heavy:
            break
        matched, missed = [], 0
        for mol in mols:
            matches = bool(mol.HasSubstructMatch(pattern, useChirality=True))
            matched.append(matches)
            missed += not matches
            # Use the original floating-point comparison, including its tolerance.
            # Even if all remaining members match, this candidate cannot meet it.
            if (len(mols) - missed) / len(mols) < min_coverage - 1e-12:
                break
        else:
            fraction = sum(matched) / len(mols)
            scored.append((False, (-heavy, -fraction), smiles, matched, heavy))
            winning_heavy = heavy
    if not scored:
        # No target-reaching candidate: preserve the full-coverage fallback ranking.
        for smiles, pattern in patterns:
            matched = _coverage(pattern, mols)
            fraction = sum(matched) / len(mols)
            heavy = pattern.GetNumHeavyAtoms()
            scored.append((True, (-fraction, -heavy), smiles, matched, heavy))
    if not scored:
        return {"core_smiles": "", "core_smiles_attached": "", "core_smiles_attached_source": "",
                "core_rule": rule, "coverage_target": min_coverage,
                "coverage_count": 0, "coverage_fraction": 0.0, "core_heavy_atoms": 0, "matched": [False] * len(mols),
                "target_met": False, "candidates": 0, "seconds": time.perf_counter() - started, **extra}
    missed, _, smiles, matched, heavy = min(scored, key=lambda item: (item[0], item[1], item[2]))
    pattern = Chem.MolFromSmiles(smiles)
    attached, attached_source = _core_smiles_with_attachments(pattern, mols)
    return {"core_smiles": smiles,
            "core_smiles_attached": attached,
            "core_smiles_attached_source": attached_source,
            "core_rule": rule, "coverage_target": min_coverage,
            "coverage_count": int(sum(matched)), "coverage_fraction": sum(matched) / len(mols),
            "core_heavy_atoms": int(heavy), "matched": matched, "target_met": not missed,
            "candidates": len(candidates), "seconds": time.perf_counter() - started, **extra}


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
