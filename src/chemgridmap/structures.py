"""Parent structures for representations, without discarding record provenance."""

from functools import lru_cache

from chembl_structure_pipeline import standardizer
from rdkit import Chem


@lru_cache(maxsize=20000)
def parent_structure(canonical_smiles: str):
    """Return ChEMBL-normalized parent SMILES and the pipeline exclusion flag.

    GetParent uses curated salt/solvent rules. It is not an unconditional
    largest-fragment operation and does not collapse stereoisomers or tautomers.
    """
    mol = Chem.MolFromSmiles(canonical_smiles)
    standardized = standardizer.standardize_mol(mol)
    parent, excluded = standardizer.get_parent_mol(standardized)
    Chem.SanitizeMol(parent)
    return Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True), bool(excluded)
