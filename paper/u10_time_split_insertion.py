"""U10 (post hoc, exploratory): neighbor evidence for later compounds in a time split.

For a ChEMBL IC50 export, a reference map is built from the records of documents
published up to ``--cutoff`` (default 2015). Compounds whose first document is
later than the cutoff are inserted into this frozen map as new molecules with
their own later records (document = evidence unit). Their measured class is then
compared with the neighbor-evidence category that ChemGridMap shows for them.

Neighbor evidence depends only on fingerprints and on the class support of the
reference entries, not on the layout, so the reference map uses PCA without
refinement to save time. Usage::

    python paper/u10_time_split_insertion.py --raw validation_data/chembl37_chembl205/chembl37_chembl205_ic50_raw.csv \
        --target CHEMBL205 --output paper/output/u10_time_split/chembl205
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

from chemgridmap import annotate_activity_confidence, build_grid_map, curate_chembl_activity_data
from chemgridmap.queries import insert_new_molecules

DOCUMENTS = Path(__file__).resolve().parents[1] / "src/chemgridmap/data/chembl37_documents.csv.gz"
CLASSES = ("inactive", "medium", "active")
BATCH = 1000


def document_years():
    with gzip.open(DOCUMENTS, "rt", encoding="utf-8") as stream:
        table = pd.read_csv(stream, usecols=["document_chembl_id", "year"])
    return table.dropna().drop_duplicates("document_chembl_id").set_index("document_chembl_id").year.astype(int)


def wilson(k, n, z=1.959964):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return (centre - half, centre + half)


def run(raw_path, target, output, *, cutoff=2015, activity_type="IC50"):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(raw_path, low_memory=False, encoding="utf-8-sig")
    years = document_years()
    raw["u10_year"] = raw.document_chembl_id.map(years)
    if "document_year" in raw:
        raw["u10_year"] = raw.u10_year.fillna(pd.to_numeric(raw.document_year, errors="coerce"))
    full = curate_chembl_activity_data(raw, target_id=target, activity_type=activity_type)
    records = full.retained_records.copy()
    records["u10_year"] = records.document_chembl_id.map(years)
    if "document_year" in records:
        records["u10_year"] = records.u10_year.fillna(pd.to_numeric(records.document_year, errors="coerce"))
    undated = int(records.u10_year.isna().sum())
    records = records[records.u10_year.notna()]
    parent = records.parent_molecule_chembl_id.fillna(records.molecule_chembl_id).astype(str)
    first_year = records.groupby(parent).u10_year.min()
    early_parents = set(first_year[first_year <= cutoff].index)
    late_parents = set(first_year[first_year > cutoff].index)

    reference_raw = raw[raw.u10_year.notna() & (raw.u10_year <= cutoff)]
    curated = curate_chembl_activity_data(reference_raw, target_id=target, activity_type=activity_type)
    molecules, fit = annotate_activity_confidence(
        curated.molecules, curated.retained_records, value="activity_pchembl_raw",
        document="document_chembl_id", confidence_threshold=0.8)
    result = build_grid_map(molecules, output_dir=output / "reference", name="reference",
                            smiles_col="canonical_smiles", value_col="activity_pchembl",
                            label_col="activity_class", projection="pca", grid_refinement="none",
                            output_formats=("svg",), render_detail="overview",
                            map_color_modes=("class_support",))

    late = records[parent.loc[records.index].isin(late_parents)].copy()
    late["query_id"] = parent.loc[late.index].to_numpy()
    late = late[["query_id", "canonical_smiles", "activity_pchembl_raw", "document_chembl_id"]]
    # One structure per new compound: the most frequent curated structure of its records.
    structure = late.groupby("query_id").canonical_smiles.agg(lambda s: s.value_counts().index[0])
    late["canonical_smiles"] = late.query_id.map(structure)
    ids = sorted(late.query_id.unique())
    tables = []
    for start in range(0, len(ids), BATCH):
        chunk = late[late.query_id.isin(ids[start:start + BATCH])]
        inserted = insert_new_molecules(
            result, chunk.reset_index(drop=True), output_dir=output / f"new_{start // BATCH + 1}",
            input_filename=f"later_compounds_{start // BATCH + 1}.csv", smiles_col="canonical_smiles",
            id_col="query_id", value_col="activity_pchembl_raw", value_unit="p",
            document_column="document_chembl_id")
        tables.append(inserted.data)
    table = pd.concat(tables, ignore_index=True)
    from chemgridmap.comparison import molecule_keys
    reference_class = dict(zip(molecule_keys(result.data), result.data.activity_class.astype(str)))
    table["nearest_reference_class"] = table.nearest_reference_id.astype(str).map(reference_class)
    # A random nearest neighbour above the same similarity threshold, for the similarity-matched baseline.
    table["nearest_reference_similar"] = table.tanimoto >= 0.6
    table.to_csv(output / "later_compounds.csv", index=False, lineterminator="\n")
    return summarize(table, target=target, cutoff=cutoff, fit=fit, n_reference=len(molecules),
                     n_late=len(ids), undated=undated, output=output)


def summarize(table, *, target, cutoff, fit, n_reference, n_late, undated, output):
    evidence_col = "neighbor_evidence"
    measured_col = next(c for c in ("measured_class", "activity_class") if c in table)
    placed = table[table[evidence_col].notna() & table[measured_col].isin(CLASSES)]
    rows = []
    for category, part in placed.groupby(evidence_col):
        row = {"neighbor_evidence": category, "n": int(len(part))}
        for cls in CLASSES:
            row[f"measured_{cls}"] = int(part[measured_col].eq(cls).sum())
        if category in CLASSES:
            agree = int(part[measured_col].eq(category).sum())
            low, high = wilson(agree, len(part))
            row.update(agree=agree, agreement=agree / len(part), ci_low=low, ci_high=high)
        rows.append(row)
    by_category = pd.DataFrame(rows)
    by_category.to_csv(output / "agreement_by_category.csv", index=False, lineterminator="\n")
    consistent = placed[placed[evidence_col].isin(CLASSES)]
    # Baselines on the same compounds: nearest reference class, and the most common later class.
    nearest_col = "nearest_reference_class"
    summary = {
        "target": target, "cutoff": cutoff, "reference_entries": n_reference,
        "reference_sigma_doc": fit.get("sigma_doc") if isinstance(fit, dict) else None,
        "later_compounds": n_late, "placed_with_measured_class": int(len(placed)),
        "undated_records_dropped": undated,
        "categories": {str(k): int(v) for k, v in placed[evidence_col].value_counts().items()},
        "consistent_n": int(len(consistent)),
        "consistent_agreement": (float(consistent[measured_col].eq(consistent[evidence_col]).mean())
                                 if len(consistent) else None),
        "later_class_distribution": {c: int(placed[measured_col].eq(c).sum()) for c in CLASSES},
    }
    if nearest_col:
        known = placed[placed[nearest_col].isin(CLASSES)]
        summary["nearest_reference_agreement_all"] = float(known[measured_col].eq(known[nearest_col]).mean())
        sub = consistent[consistent[nearest_col].isin(CLASSES)]
        summary["nearest_reference_agreement_on_consistent"] = (
            float(sub[measured_col].eq(sub[nearest_col]).mean()) if len(sub) else None)
    similar = placed[placed.nearest_reference_similar & placed[nearest_col].isin(CLASSES)]
    summary["nearest_reference_agreement_similar"] = (
        float(similar[measured_col].eq(similar[nearest_col]).mean()) if len(similar) else None)
    summary["nearest_reference_similar_n"] = int(len(similar))
    if len(placed):
        summary["majority_class_rate"] = float(placed[measured_col].value_counts(normalize=True).iloc[0])
    (output / "summary.json").write_text(json.dumps(summary, indent=2, default=float), encoding="utf-8")
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cutoff", type=int, default=2015)
    args = parser.parse_args(argv)
    print(json.dumps(run(args.raw, args.target, args.output, cutoff=args.cutoff), indent=2, default=float))


if __name__ == "__main__":
    main()
