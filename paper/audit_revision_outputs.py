"""Check the saved record-to-cell evidence chain without rerunning projections."""
from pathlib import Path
import argparse
import hashlib
import json
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd


def audit(root):
    checks = {}
    datasets = pd.read_csv(root / "datasets.csv")
    for item in datasets.itertuples():
        slug = item.target.lower()
        directory = root / slug
        retained = pd.read_csv(directory / f"{slug}_retained_records.csv")
        excluded = pd.read_csv(directory / f"{slug}_excluded_records.csv")
        molecules = pd.read_csv(directory / f"{slug}_molecule_level.csv")
        key = "molecule_identity_key"
        grouped = retained.groupby(key).activity_pchembl_raw
        medians = grouped.median().reindex(molecules[key]).to_numpy()
        counts = grouped.size().reindex(molecules[key]).to_numpy()
        checks[slug] = {
            "raw_partition": len(retained) + len(excluded) == item.raw_records,
            "retained_count": len(retained) == item.retained_records,
            "molecule_count": len(molecules) == item.curated_molecules,
            "identities_unique": bool(molecules[key].is_unique),
            "median_recalculation": bool(np.allclose(medians, molecules.activity_pchembl)),
            "record_count_recalculation": bool(np.array_equal(counts, molecules.n_records)),
            "partition_disjoint": not bool(set(retained.source_row) & set(excluded.source_row)),
            "excluded_reasons_present": bool(excluded.exclusion_reason.notna().all()),
        }
        layouts = list(directory.glob("coordinates_seed_*.csv"))
        checks[slug]["unique_grid_cells_all_seeds"] = all(
            not pd.read_csv(path).duplicated(["grid_col", "grid_row"]).any()
            for path in layouts
        ) and bool(layouts)
    main = pd.read_csv(root / "chembl205/coordinates_seed_42.csv").set_index("molecule_identity_key").sort_index()
    direct_path = root / "direct_csv_check/chembl205_grid_coordinates.csv"
    if direct_path.exists():
        direct = pd.read_csv(direct_path).set_index("molecule_identity_key").sort_index()
        checks["public_cli"] = {
            "same_identities": bool(main.index.equals(direct.index)),
            "same_projection": bool(np.array_equal(main[["projection_x_raw", "projection_y_raw"]], direct[["projection_x_raw", "projection_y_raw"]])),
            "same_grid": bool(np.array_equal(main[["grid_col", "grid_row"]], direct[["grid_col", "grid_row"]])),
        }
    adaptive = pd.read_csv(root / "adaptive_validation.csv")
    checks["adaptive_calibration"] = {
        "eleven_layouts": len(adaptive) == 11,
        "dense_cost_agreement": bool(np.allclose(adaptive.dense_cost_ratio, 1, rtol=1e-8)),
    }
    svg = ET.parse(root / "chembl205/molecular_grid.svg").getroot()
    checks["molecular_svg"] = {
        "native_paths": any(el.tag.endswith("}path") for el in svg.iter()),
        "no_bitmap_molecules": not any(el.tag.endswith("}image") for el in svg.iter()),
        "all_molecules_traceable": sum("data-molecule-id" in el.attrib for el in svg.iter()) == len(main),
    }
    passed = all(all(group.values()) for group in checks.values())
    result = {"passed": passed, "checks": checks,
              "scope": "Saved-output consistency audit, not independent biological validation.",
              "csv_sha256": {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
                             for path in sorted(root.rglob("*.csv"))}}
    (root / "integrity_audit.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"passed": passed, "checks": checks}, indent=2))
    if not passed:
        raise SystemExit("Saved-output integrity check failed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path(__file__).parent / "output/revision_v4")
    audit(parser.parse_args().source.resolve())
