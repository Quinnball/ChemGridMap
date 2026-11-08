"""Public interface for ChemGridMap."""

from .chembl import ChemblCurationResult, curate_chembl_activity_data
from .comparison import structural_groups
from .cores import mcs_core, select_core
from .core import assign_to_grid, canonicalize_smiles, normalize_coordinates, prepare_molecules
from .evidence import (
    annotate_activity_confidence,
    class_probabilities,
    document_pair_offsets,
    document_values,
    fit_between_unit_sd,
    prepare_generic_activity_records,
)
from .pipeline import GridMapResult, build_grid_map
from .queries import NewMoleculesResult, insert_new_molecules, prepare_query_molecules
from .inspection import inspect_entry
from .reliability import (
    grid_neighbor_sets,
    neighbor_recall,
    neighbor_recall_fixed_k,
    representation_neighbor_sets,
)

__all__ = [
    "ChemblCurationResult",
    "GridMapResult",
    "NewMoleculesResult",
    "insert_new_molecules",
    "annotate_activity_confidence",
    "assign_to_grid",
    "build_grid_map",
    "canonicalize_smiles",
    "class_probabilities",
    "curate_chembl_activity_data",
    "document_pair_offsets",
    "document_values",
    "fit_between_unit_sd",
    "prepare_generic_activity_records",
    "prepare_molecules",
    "prepare_query_molecules",
    "grid_neighbor_sets",
    "mcs_core",
    "neighbor_recall",
    "neighbor_recall_fixed_k",
    "normalize_coordinates",
    "inspect_entry",
    "representation_neighbor_sets",
    "select_core",
    "structural_groups",
]

__version__ = "0.5.0"
