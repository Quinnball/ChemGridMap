"""Molecular representations supported by ChemGridMap."""

from __future__ import annotations

from typing import List, Tuple

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors
from sklearn.preprocessing import StandardScaler


DESCRIPTOR_NAMES = [
    "MolWt",
    "MolLogP",
    "TPSA",
    "NumHDonors",
    "NumHAcceptors",
    "NumRotatableBonds",
    "RingCount",
    "HeavyAtomCount",
    "FractionCSP3",
    "BalabanJ",
]

DESCRIPTOR_FUNCTIONS = [
    Descriptors.MolWt,
    Descriptors.MolLogP,
    Descriptors.TPSA,
    Descriptors.NumHDonors,
    Descriptors.NumHAcceptors,
    Descriptors.NumRotatableBonds,
    Descriptors.RingCount,
    Descriptors.HeavyAtomCount,
    Descriptors.FractionCSP3,
    Descriptors.BalabanJ,
]


def morgan_fingerprints(
    smiles: List[str],
    radius: int = 2,
    n_bits: int = 2048,
) -> np.ndarray:
    """Calculate radius-2, 2048-bit Morgan fingerprints by default."""
    rows = []
    for item in smiles:
        mol = Chem.MolFromSmiles(str(item))
        if mol is None:
            raise ValueError("Invalid canonical SMILES encountered: {!r}".format(item))
        fingerprint = AllChem.GetMorganFingerprintAsBitVect(
            mol, int(radius), nBits=int(n_bits)
        )
        array = np.zeros((int(n_bits),), dtype=np.float32)
        DataStructs.ConvertToNumpyArray(fingerprint, array)
        rows.append(array)
    return np.vstack(rows)


def rdkit_descriptors(smiles: List[str]) -> np.ndarray:
    """Calculate and standardize ten common RDKit descriptors."""
    rows = []
    for item in smiles:
        mol = Chem.MolFromSmiles(str(item))
        if mol is None:
            raise ValueError("Invalid canonical SMILES encountered: {!r}".format(item))
        rows.append([float(function(mol)) for function in DESCRIPTOR_FUNCTIONS])
    matrix = np.asarray(rows, dtype=np.float32)
    return StandardScaler().fit_transform(matrix)


def external_embedding(
    data: pd.DataFrame,
    prefix: str = "emb_",
) -> Tuple[np.ndarray, List[str]]:
    """Read an external embedding from columns sharing a prefix."""
    columns = sorted(
        [column for column in data.columns if str(column).startswith(prefix)]
    )
    if not columns:
        raise ValueError(
            "No embedding columns were found with prefix {!r}.".format(prefix)
        )
    matrix = data[columns].apply(pd.to_numeric, errors="raise").to_numpy(np.float32)
    return matrix, columns


def compute_representation(
    data: pd.DataFrame,
    representation: str = "morgan",
    radius: int = 2,
    n_bits: int = 2048,
    embedding_prefix: str = "emb_",
) -> Tuple[np.ndarray, List[str]]:
    """Build one of the supported molecular representation matrices."""
    method = representation.lower()
    smiles = data["canonical_smiles"].astype(str).tolist()
    if method == "morgan":
        matrix = morgan_fingerprints(smiles, radius=radius, n_bits=n_bits)
        feature_names = ["bit_{:04d}".format(index) for index in range(n_bits)]
    elif method == "descriptors":
        matrix = rdkit_descriptors(smiles)
        feature_names = list(DESCRIPTOR_NAMES)
    elif method == "embedding":
        matrix, feature_names = external_embedding(data, prefix=embedding_prefix)
    else:
        raise ValueError(
            "Unsupported representation {!r}. Choose morgan, descriptors or embedding.".format(
                representation
            )
        )
    return np.asarray(matrix, dtype=np.float32), feature_names

