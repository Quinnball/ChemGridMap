"""E28: neighbor evidence for later compounds inserted into an earlier reference map.

Implements the plan ``E28_neighbor_evidence_plan.md`` (Additional file 2). For each
ChEMBL target-endpoint dataset:

1. Records are dated by their document year (bundled ChEMBL 37 document table,
   then the export's ``document_year``); undated records are dropped and counted.
2. Entries are curated from all dated records with default settings. An entry whose
   earliest dated record is at or before the cutoff is a reference entry; every other
   entry is a later compound.
3. The reference map is built from the records dated at or before the cutoff only
   (default curation, class support at 0.8, PCA without refinement: neighbor
   evidence uses fingerprints and class support, not the layout).
4. Later compounds are inserted with ``insert_new_molecules`` using their own records
   (document = evidence unit); their measured class is the class of their median.
5. Each later compound receives ChemGridMap's neighbor-evidence category.

Outputs per dataset (``compounds.csv``, ``dataset.json``) and, with ``aggregate``,
per-dataset metrics and the predeclared tests N1-N3. Usage::

    python paper/u10_neighbor_validation.py dataset --raw FILE --target CHEMBLxxx --endpoint IC50 --output DIR
    python paper/u10_neighbor_validation.py aggregate --root DIR --output DIR
"""
from __future__ import annotations

import argparse
import gzip
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from chemgridmap import annotate_activity_confidence, build_grid_map, curate_chembl_activity_data
from chemgridmap.comparison import molecule_keys
from chemgridmap.queries import insert_new_molecules

DOCUMENTS = Path(__file__).resolve().parents[1] / "src/chemgridmap/data/chembl37_documents.csv.gz"
CLASSES = ("inactive", "medium", "active")
CONSISTENT = set(CLASSES)
BATCH = 1000
MIN_REFERENCE, MIN_LATER = 100, 20
MIN_CONSISTENT, MIN_GROUP = 10, 5


def document_years():
    with gzip.open(DOCUMENTS, "rt", encoding="utf-8") as stream:
        table = pd.read_csv(stream, usecols=["document_chembl_id", "year"])
    return table.dropna().drop_duplicates("document_chembl_id").set_index("document_chembl_id").year.astype(int)


def _date(frame, years):
    dated = frame.document_chembl_id.map(years)
    if "document_year" in frame:
        dated = dated.fillna(pd.to_numeric(frame.document_year, errors="coerce"))
    return dated


def _curate(raw, target, endpoint):
    try:
        return curate_chembl_activity_data(raw, target_id=target, activity_type=endpoint), "error"
    except ValueError as error:
        if "identity" not in str(error).lower() and "InChIKey" not in str(error):
            raise
        return curate_chembl_activity_data(raw, target_id=target, activity_type=endpoint,
                                           identity_conflicts="exclude"), "exclude"


def run_dataset(raw_path, target, endpoint, output, *, cutoff=2015):
    started = time.perf_counter()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(raw_path, low_memory=False, encoding="utf-8-sig")
    years = document_years()
    raw["e27_year"] = _date(raw, years)
    info = {"target": target, "endpoint": endpoint, "cutoff": cutoff, "raw_records": int(len(raw)),
            "undated_raw_records": int(raw.e27_year.isna().sum())}
    raw = raw[raw.e27_year.notna()].copy()
    full, info["identity_conflicts_full"] = _curate(raw, target, endpoint)
    records = full.retained_records
    first = records.assign(e27_year=_date(records, years)).groupby("molecule_identity_key").e27_year.min()
    later_keys = set(first[first > cutoff].index)
    info["entries_all"] = int(len(full.molecules))
    reference_raw = raw[raw.e27_year <= cutoff]
    info["reference_raw_records"] = int(len(reference_raw))
    late_records = records[records.molecule_identity_key.isin(later_keys)].copy()
    info["later_compounds"] = int(len(later_keys))
    if len(reference_raw) == 0:
        info.update(included=False, reason="no reference records")
        (output / "dataset.json").write_text(json.dumps(info, indent=2), encoding="utf-8")
        return info
    curated, info["identity_conflicts_reference"] = _curate(reference_raw, target, endpoint)
    info["reference_entries"] = int(len(curated.molecules))
    if info["reference_entries"] < MIN_REFERENCE or info["later_compounds"] < MIN_LATER:
        info.update(included=False, reason=f"needs >= {MIN_REFERENCE} reference entries and >= {MIN_LATER} later compounds")
        (output / "dataset.json").write_text(json.dumps(info, indent=2), encoding="utf-8")
        return info
    molecules, fit = annotate_activity_confidence(
        curated.molecules, curated.retained_records, value="activity_pchembl_raw",
        document="document_chembl_id", confidence_threshold=0.8)
    info["reference_sigma_doc"] = fit.get("sigma_doc")
    info["reference_supported"] = int(molecules.class_status.eq("supported").sum()) if "class_status" in molecules else None
    map_options = dict(output_dir=output / "reference", name="reference",
                       smiles_col="canonical_smiles", value_col="activity_pchembl",
                       label_col="activity_class", projection="pca", grid_refinement="none",
                       output_formats=["svg"], render_detail="overview", grouping="none",
                       map_color_modes=["class_support"])
    try:
        result = build_grid_map(molecules, **map_options)
    except ValueError as error:
        if "repeated canonical structures" not in str(error):
            raise
        # Deviation D1 (run error): distinct curated entries (stereoisomers) share one
        # canonical structure, so the map is built with every entry kept.
        result = build_grid_map(molecules, duplicate_policy="keep", **map_options)
        info["reference_duplicate_policy"] = "keep"
    late = late_records[["molecule_identity_key", "canonical_smiles", "activity_pchembl_raw",
                         "document_chembl_id"]].rename(columns={"molecule_identity_key": "query_id"})
    structure = late.groupby("query_id").canonical_smiles.agg(lambda s: s.value_counts().sort_index().idxmax())
    late["canonical_smiles"] = late.query_id.map(structure)
    ids = sorted(late.query_id.unique())
    placed, checks = [], []
    # Deviation D2 (run error): a batch may not exceed the empty cells of the reference
    # grid. Neighbor evidence does not depend on the batch, so smaller batches change
    # nothing else; datasets with enough empty cells keep batches of BATCH.
    grid = getattr(result, "insertion_parameters", {})
    rows = int(grid.get("grid_rows", result.metrics.iloc[0].grid_side))
    cols = int(grid.get("grid_cols", rows))
    batch = max(1, min(BATCH, rows * cols - len(result.data)))
    if batch < BATCH:
        info["insertion_batch"] = int(batch)
    for start in range(0, len(ids), batch):
        chunk = late[late.query_id.isin(ids[start:start + batch])].reset_index(drop=True)
        inserted = insert_new_molecules(
            result, chunk, output_dir=output / f"later_{start // batch + 1}",
            input_filename=f"later_compounds_{start // batch + 1}.csv", smiles_col="canonical_smiles",
            id_col="query_id", value_col="activity_pchembl_raw", value_unit="p",
            document_column="document_chembl_id")
        placed.append(inserted.data)
        checks.append(inserted.all_molecules)
    table = pd.concat(placed, ignore_index=True)
    checked = pd.concat(checks, ignore_index=True)
    info["later_status_counts"] = {str(k): int(v) for k, v in checked.status.value_counts().items()} if "status" in checked else {}
    reference_class = dict(zip(molecule_keys(result.data), result.data.activity_class.astype(str)))
    table["nearest_reference_class"] = table.nearest_reference_id.astype(str).map(reference_class)
    keep = ["query_id", "canonical_smiles", "activity_class", "neighbor_evidence", "similar_reference_count",
            "supported_reference_count", "nearest_reference_id", "tanimoto", "nearest_reference_class"]
    table = table[[c for c in keep if c in table]].rename(columns={"activity_class": "measured_class"})
    table.insert(0, "endpoint", endpoint)
    table.insert(0, "target", target)
    table.to_csv(output / "compounds.csv", index=False, lineterminator="\n")
    info.update(included=True, placed=int(len(table)), seconds=time.perf_counter() - started)
    (output / "dataset.json").write_text(json.dumps(info, indent=2, default=float), encoding="utf-8")
    return info


def wilson(k, n, z=1.959964):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return (centre - half, centre + half)


def dataset_metrics(table):
    """Per-dataset quantities for N1-N3 (definitions in the plan, section 4)."""
    rows = table[table.measured_class.isin(CLASSES)]
    consistent = rows[rows.neighbor_evidence.isin(CONSISTENT)]
    out = {"n": int(len(rows)), "n_consistent": int(len(consistent))}
    for category in ("none", "insufficient", "mixed", *CLASSES):
        out[f"n_{category}"] = int(rows.neighbor_evidence.eq(category).sum())
    if len(consistent):
        observed = float(consistent.measured_class.eq(consistent.neighbor_evidence).mean())
        p_measured = consistent.measured_class.value_counts(normalize=True)
        p_category = consistent.neighbor_evidence.value_counts(normalize=True)
        expected = float(sum(p_measured.get(c, 0) * p_category.get(c, 0) for c in CLASSES))
        out.update(consistent_agreement=observed, consistent_chance=expected, n1_difference=observed - expected)
    similar = rows[rows.tanimoto >= 0.6 - 1e-12]
    nearest = similar[similar.nearest_reference_class.isin(CLASSES)]
    agree = nearest.measured_class.eq(nearest.nearest_reference_class)
    groups = {"consistent": nearest.neighbor_evidence.isin(CONSISTENT),
              "unsettled": nearest.neighbor_evidence.isin({"insufficient", "mixed"})}
    for name, mask in groups.items():
        out[f"nearest_{name}_n"] = int(mask.sum())
        out[f"nearest_{name}_agreement"] = float(agree[mask].mean()) if mask.any() else np.nan
    out["n2_difference"] = out["nearest_consistent_agreement"] - out["nearest_unsettled_agreement"]
    far = rows[rows.neighbor_evidence.eq("none") & rows.nearest_reference_class.isin(CLASSES)]
    out["nearest_none_n"] = int(len(far))
    out["nearest_none_agreement"] = float(far.measured_class.eq(far.nearest_reference_class).mean()) if len(far) else np.nan
    out["n3_difference"] = out["nearest_consistent_agreement"] - out["nearest_none_agreement"]
    return out


def holm(pvalues):
    order = np.argsort(pvalues)
    adjusted, running = np.empty(len(pvalues)), 0.0
    for rank, index in enumerate(order):
        running = max(running, min(1.0, (len(pvalues) - rank) * pvalues[index]))
        adjusted[index] = running
    return adjusted


def aggregate(root, output, *, primary_role="sampled"):
    root, output = Path(root), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    infos, tables = [], []
    for path in sorted(root.glob("*/dataset.json")):
        info = json.loads(path.read_text(encoding="utf-8"))
        info["dataset"] = path.parent.name
        infos.append(info)
        if info.get("included") and (path.parent / "compounds.csv").exists():
            tables.append(pd.read_csv(path.parent / "compounds.csv"))
    inclusion = pd.DataFrame(infos)
    inclusion.to_csv(output / "dataset_inclusion.csv", index=False, lineterminator="\n")
    compounds = pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()
    metrics = []
    for (target, endpoint), table in compounds.groupby(["target", "endpoint"], sort=True):
        metrics.append({"target": target, "endpoint": endpoint, **dataset_metrics(table)})
    metrics = pd.DataFrame(metrics)
    metrics.to_csv(output / "metrics_by_dataset.csv", index=False, lineterminator="\n")
    tests = []
    for endpoint, part in metrics.groupby("endpoint"):
        rows = []
        for name, column, mask in (
                ("N1", "n1_difference", part.n_consistent >= MIN_CONSISTENT),
                ("N2", "n2_difference", (part.nearest_consistent_n >= MIN_GROUP) & (part.nearest_unsettled_n >= MIN_GROUP)),
                ("N3", "n3_difference", (part.nearest_consistent_n >= MIN_GROUP) & (part.nearest_none_n >= MIN_GROUP))):
            values = part.loc[mask, column].dropna().to_numpy()
            if len(values) >= 5 and np.any(values != 0):
                p = float(stats.wilcoxon(values, alternative="greater").pvalue)
            else:
                p = np.nan
            rows.append({"endpoint": endpoint, "test": name, "datasets": int(len(values)),
                         "median_difference": float(np.median(values)) if len(values) else np.nan,
                         "positive": int((values > 0).sum()), "p_one_sided": p})
        valid = [i for i, row in enumerate(rows) if not np.isnan(row["p_one_sided"])]
        adjusted = holm(np.array([rows[i]["p_one_sided"] for i in valid])) if valid else []
        for i, value in zip(valid, adjusted):
            rows[i]["p_holm"] = float(value)
        tests.extend(rows)
    tests = pd.DataFrame(tests)
    tests.to_csv(output / "hypothesis_tests.csv", index=False, lineterminator="\n")
    pooled = []
    labelled = compounds[compounds.measured_class.isin(CLASSES)] if len(compounds) else compounds
    for (endpoint, category), part in labelled.groupby(["endpoint", "neighbor_evidence"]):
        row = {"endpoint": endpoint, "neighbor_evidence": category, "n": int(len(part)), "datasets": int(part.target.nunique())}
        for cls in CLASSES:
            row[f"measured_{cls}"] = int(part.measured_class.eq(cls).sum())
        if category in CONSISTENT:
            k = int(part.measured_class.eq(category).sum())
            row.update(agree=k, agreement=k / len(part), ci_low=wilson(k, len(part))[0], ci_high=wilson(k, len(part))[1])
        pooled.append(row)
    pd.DataFrame(pooled).to_csv(output / "pooled_by_category.csv", index=False, lineterminator="\n")
    summary = {"datasets": int(len(inclusion)), "included": int(inclusion.get("included", pd.Series(dtype=bool)).fillna(False).sum()),
               "tests": tests.to_dict(orient="records")}
    (output / "summary.json").write_text(json.dumps(summary, indent=2, default=float), encoding="utf-8")
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    one = sub.add_parser("dataset")
    one.add_argument("--raw", type=Path, required=True)
    one.add_argument("--target", required=True)
    one.add_argument("--endpoint", default="IC50")
    one.add_argument("--output", type=Path, required=True)
    one.add_argument("--cutoff", type=int, default=2015)
    many = sub.add_parser("aggregate")
    many.add_argument("--root", type=Path, required=True)
    many.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.command == "dataset":
        result = run_dataset(args.raw, args.target, args.endpoint, args.output, cutoff=args.cutoff)
    else:
        result = aggregate(args.root, args.output)
    print(json.dumps(result, indent=2, default=float))


if __name__ == "__main__":
    main()
