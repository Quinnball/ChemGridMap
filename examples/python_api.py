"""Curate ChEMBL records, build a map, and save the reliability summary."""

from pathlib import Path

import pandas as pd

from chemgridmap import (
    annotate_activity_confidence,
    build_grid_map,
    curate_chembl_activity_data,
)


def run(input_path, output_dir, *, target_id="CHEMBL205",
        activity_type="IC50", projection="umap"):
    raw = pd.read_csv(input_path, encoding="utf-8-sig")
    curated = curate_chembl_activity_data(
        raw, target_id=target_id, activity_type=activity_type
    )
    document_column = curated.report["detected_columns"]["document_id"]
    molecules = curated.molecules
    if document_column:
        molecules, fit = annotate_activity_confidence(
            molecules,
            curated.retained_records,
            value="activity_pchembl_raw",
            document=document_column,
        )
        curated.report["activity_confidence"] = fit
    result = build_grid_map(
        molecules,
        output_dir=output_dir,
        name="caii",
        smiles_col="canonical_smiles",
        value_col="activity_pchembl",
        label_col="activity_class",
        projection=projection,
        map_color_modes=(
            ["class_support", "neighbor_recall"]
            if "class_status" in molecules else ["neighbor_recall"]
        ),
    )
    columns = [
        "molecule_identity_key", "activity_pchembl", "class_confidence",
        "class_status", "status_reason", "offset_pair_touched",
        "neighbor_recall", "neighbor_recall_fixed_k", "structure_group",
    ]
    result.data[columns].to_csv(
        Path(output_dir) / "entry_reliability.csv", index=False
    )
    core_columns = [
        "structure_group", "size", "core_smiles", "core_smiles_attached",
        "core_coverage_count", "core_coverage_fraction",
    ]
    pd.DataFrame(result.structure_groups)[core_columns].to_csv(
        Path(output_dir) / "group_cores.csv", index=False
    )
    return result


if __name__ == "__main__":
    run("chembl205_ic50.csv", "caii_map")
