from pathlib import Path

import pandas as pd

from examples.python_api import run


def test_python_api_example_runs_with_small_chEMBL_data(tmp_path):
    source = tmp_path / "small.csv"
    pd.DataFrame([
        {"Molecule ChEMBL ID": "CHEMBL1", "Smiles": "CCO",
         "pChEMBL Value": 5.8, "Document ChEMBL ID": "D1",
         "Target ChEMBL ID": "CHEMBL205", "Standard Type": "IC50",
         "Standard Relation": "="},
        {"Molecule ChEMBL ID": "CHEMBL1", "Smiles": "CCO",
         "pChEMBL Value": 6.2, "Document ChEMBL ID": "D2",
         "Target ChEMBL ID": "CHEMBL205", "Standard Type": "IC50",
         "Standard Relation": "="},
        {"Molecule ChEMBL ID": "CHEMBL2", "Smiles": "c1ccccc1",
         "pChEMBL Value": 7.2, "Document ChEMBL ID": "D1",
         "Target ChEMBL ID": "CHEMBL205", "Standard Type": "IC50",
         "Standard Relation": "="},
        {"Molecule ChEMBL ID": "CHEMBL3", "Smiles": "Cc1ccccc1",
         "pChEMBL Value": 6.5, "Document ChEMBL ID": "D2",
         "Target ChEMBL ID": "CHEMBL205", "Standard Type": "IC50",
         "Standard Relation": "="},
        {"Molecule ChEMBL ID": "CHEMBL4", "Smiles": "Oc1ccccc1",
         "pChEMBL Value": 8.0, "Document ChEMBL ID": "D3",
         "Target ChEMBL ID": "CHEMBL205", "Standard Type": "IC50",
         "Standard Relation": "="},
    ]).to_csv(source, index=False)
    output = tmp_path / "output"
    run(source, output, projection="pca")
    summary = pd.read_csv(output / "entry_reliability.csv")
    cores = pd.read_csv(output / "group_cores.csv")
    assert {"class_status", "neighbor_recall_fixed_k", "structure_group"} <= set(summary)
    assert "core_smiles_attached" in cores
    assert Path(output / "caii.svg").exists()
