"""Post hoc E28 analyses (R11, R12, S18, S19): rebuild each reference table without the map.

Repeats u10_neighbor_validation.run_dataset up to annotate_activity_confidence (same code,
same curation, same class support) and writes the reference entries with their class
status and documents, and the later compounds with their documents. No layout is built,
because neighbor evidence does not depend on it.
"""
import json, os, sys, time
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # paper/, for u10_neighbor_validation
import u10_neighbor_validation as u
from chemgridmap import annotate_activity_confidence

def main(raw_path, target, endpoint, out, cutoff=2015):
    t0 = time.perf_counter()
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(raw_path, low_memory=False, encoding="utf-8-sig")
    years = u.document_years()
    raw["e27_year"] = u._date(raw, years)
    raw = raw[raw.e27_year.notna()].copy()
    full, _ = u._curate(raw, target, endpoint)
    records = full.retained_records.assign(e27_year=lambda r: u._date(r, years))
    first = records.groupby("molecule_identity_key").e27_year.min()
    later_keys = set(first[first > cutoff].index)
    reference_raw = raw[raw.e27_year <= cutoff]
    curated, _ = u._curate(reference_raw, target, endpoint)
    molecules, fit = annotate_activity_confidence(
        curated.molecules, curated.retained_records, value="activity_pchembl_raw",
        document="document_chembl_id", confidence_threshold=0.8)
    ref_docs = curated.retained_records.groupby("molecule_identity_key").document_chembl_id.agg(
        lambda s: ";".join(sorted(set(map(str, s)))))
    all_docs = records.groupby("molecule_identity_key").document_chembl_id.agg(
        lambda s: ";".join(sorted(set(map(str, s)))))
    keep = [c for c in ["molecule_identity_key", "canonical_smiles", "activity_class", "class_status",
                        "class_confidence", "class_reproducibility", "n_evidence_units", "activity_pchembl",
                        "document_balanced_pchembl"] if c in molecules]
    ref = molecules[keep].copy()
    ref["documents_reference"] = ref.molecule_identity_key.map(ref_docs)
    ref["documents_all"] = ref.molecule_identity_key.map(all_docs)
    ref.to_csv(out / "reference.csv", index=False)
    late = records[records.molecule_identity_key.isin(later_keys)]
    later = late.groupby("molecule_identity_key").document_chembl_id.agg(
        lambda s: ";".join(sorted(set(map(str, s))))).rename("documents").reset_index()
    later.rename(columns={"molecule_identity_key": "query_id"}).to_csv(out / "later_documents.csv", index=False)
    # assay types of the reference records, for the record
    info = {"target": target, "endpoint": endpoint, "reference_entries": int(len(ref)),
            "sigma_doc": fit.get("sigma_doc"), "seconds": time.perf_counter() - t0}
    (out / "info.json").write_text(json.dumps(info, indent=2, default=float))
    print(json.dumps(info, default=float))

if __name__ == "__main__":
    main(*sys.argv[1:5])
