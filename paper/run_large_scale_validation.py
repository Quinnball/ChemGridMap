"""Download and benchmark a large ChEMBL target dataset with ChemGridMap."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import resource
import sys
import time
from typing import Dict, List
from urllib.parse import urlencode, urljoin
from urllib.request import Request, urlopen

import pandas as pd
import matplotlib
import numpy as np
import rdkit
import scipy
import sklearn

from chemgridmap import build_grid_map, curate_chembl_activity_data
from chemgridmap.chembl import save_chembl_curation


API_ROOT = "https://www.ebi.ac.uk/chembl/api/data/"
ACTIVITY_FIELDS = [
    "activity_id",
    "assay_chembl_id",
    "canonical_smiles",
    "data_validity_comment",
    "document_chembl_id",
    "molecule_chembl_id",
    "pchembl_value",
    "potential_duplicate",
    "standard_relation",
    "standard_type",
    "standard_units",
    "standard_value",
    "target_chembl_id",
    "target_pref_name",
]


def fetch_json(url: str, attempts: int = 5) -> Dict[str, object]:
    """Fetch one ChEMBL JSON response with bounded retry."""
    request = Request(url, headers={"User-Agent": "ChemGridMap/large-validation"})
    for attempt in range(attempts):
        try:
            with urlopen(request, timeout=90) as response:
                return json.load(response)
        except Exception:
            if attempt + 1 == attempts:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError("Unreachable retry state.")


def download_activities(target_id: str) -> List[Dict[str, object]]:
    """Download IC50 records with non-null pChEMBL values for one target."""
    query = urlencode(
        {
            "target_chembl_id": target_id,
            "standard_type": "IC50",
            "pchembl_value__isnull": "false",
            "limit": 1000,
            "only": ",".join(ACTIVITY_FIELDS),
        }
    )
    next_url = f"{API_ROOT}activity.json?{query}"
    records: List[Dict[str, object]] = []
    while next_url:
        payload = fetch_json(next_url)
        records.extend(payload["activities"])
        next_path = payload["page_meta"]["next"]
        next_url = urljoin(API_ROOT, next_path) if next_path else ""
        print(f"Downloaded {len(records):,} activity records.", flush=True)
    return records


def peak_rss_mib() -> float:
    """Return peak resident memory in MiB on macOS or Linux."""
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    if sys.platform == "darwin":
        return value / (1024.0 * 1024.0)
    return value / 1024.0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-id", default="CHEMBL240")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("validation_data/chembl37_chembl240_large_scale"),
    )
    parser.add_argument("--skip-download", action="store_true")
    parser.add_argument("--skip-map", action="store_true")
    parser.add_argument(
        "--hardware-label",
        help="Optional human-readable hardware label stored with benchmark metadata.",
    )
    parser.add_argument(
        "--timing-repeats",
        type=int,
        default=1,
        help="Number of complete curation and map runs used for timing summaries.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    target_id = args.target_id.strip().upper()
    raw_path = output_dir / f"chembl37_{target_id.lower()}_ic50_raw.csv"

    started = time.perf_counter()
    status_path = output_dir / "chembl_status.json"
    target_path = output_dir / "target_metadata.json"
    if args.skip_download and status_path.exists() and target_path.exists():
        status = json.loads(status_path.read_text(encoding="utf-8"))
        target = json.loads(target_path.read_text(encoding="utf-8"))
    else:
        status = fetch_json(f"{API_ROOT}status.json")
        target = fetch_json(f"{API_ROOT}target/{target_id}.json")
        status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
        target_path.write_text(json.dumps(target, indent=2), encoding="utf-8")

    if args.skip_download:
        raw = pd.read_csv(raw_path)
    else:
        raw = pd.DataFrame(download_activities(target_id))
        raw.to_csv(raw_path, index=False)
    download_seconds = time.perf_counter() - started

    if args.timing_repeats < 1:
        raise ValueError("--timing-repeats must be at least 1.")

    curation_times = []
    curation = None
    for _ in range(args.timing_repeats):
        curation_started = time.perf_counter()
        curation = curate_chembl_activity_data(
            raw,
            target_id=target_id,
            activity_type="IC50",
        )
        curation_times.append(time.perf_counter() - curation_started)
    assert curation is not None
    curation.report["source_url"] = (
        f"{API_ROOT}activity.json?target_chembl_id={target_id}"
        "&standard_type=IC50&pchembl_value__isnull=false"
    )
    curation.report["chembl_version"] = status.get("chembl_db_version")
    curation.report["chembl_release_date"] = status.get("chembl_release_date")
    curation_files = save_chembl_curation(
        curation,
        output_dir=output_dir,
        name=f"chembl37_{target_id.lower()}_ic50",
    )
    curation_seconds = float(np.median(curation_times))

    benchmark = {
        "target_id": target_id,
        "target_name": target.get("pref_name"),
        "chembl_version": status.get("chembl_db_version"),
        "chembl_release_date": status.get("chembl_release_date"),
        "raw_activity_records": int(len(raw)),
        "curation": curation.report,
        "timing_seconds": {
            "download_and_metadata": download_seconds,
            "curation": curation_seconds,
            "curation_runs": curation_times,
            "curation_q25": float(np.quantile(curation_times, 0.25)),
            "curation_q75": float(np.quantile(curation_times, 0.75)),
        },
        "peak_rss_mib_after_curation": peak_rss_mib(),
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "logical_cpu_count": os.cpu_count(),
            "hardware_label": args.hardware_label,
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "scikit_learn": sklearn.__version__,
            "rdkit": rdkit.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "curation_files": curation_files,
    }

    if not args.skip_map:
        map_times = []
        result = None
        first_map_peak_rss_mib = None
        for repeat_index in range(args.timing_repeats):
            map_started = time.perf_counter()
            result = build_grid_map(
                curation.molecules,
                output_dir=output_dir,
                name=f"chembl37_{target_id.lower()}_morgan_pca_grid",
                smiles_col="canonical_smiles",
                value_col="activity_pchembl",
                label_col="activity_class",
                representation="morgan",
                projection="pca",
                assignment_method="auto",
                render_detail="auto",
                output_formats=["svg", "png", "pdf"],
                duplicate_policy="error",
            )
            map_times.append(time.perf_counter() - map_started)
            if repeat_index == 0:
                first_map_peak_rss_mib = peak_rss_mib()
        assert result is not None
        map_seconds = float(np.median(map_times))
        combined_runs = [
            curation_time + map_time
            for curation_time, map_time in zip(curation_times, map_times)
        ]
        benchmark["timing_seconds"].update(
            {
                "map_total": map_seconds,
                "map_total_runs": map_times,
                "map_total_q25": float(np.quantile(map_times, 0.25)),
                "map_total_q75": float(np.quantile(map_times, 0.75)),
                "curation_and_map": float(np.median(combined_runs)),
                "curation_and_map_runs": combined_runs,
                "curation_and_map_q25": float(np.quantile(combined_runs, 0.25)),
                "curation_and_map_q75": float(np.quantile(combined_runs, 0.75)),
            }
        )
        benchmark["peak_rss_mib_after_first_map"] = first_map_peak_rss_mib
        benchmark["peak_rss_mib_after_all_map_repeats"] = peak_rss_mib()
        benchmark["map_metrics"] = result.metrics.iloc[0].to_dict()
        benchmark["map_files"] = result.output_files

    benchmark_path = output_dir / "large_scale_benchmark.json"
    benchmark_path.write_text(
        json.dumps(benchmark, indent=2, default=float), encoding="utf-8"
    )
    print(json.dumps(benchmark, indent=2, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
