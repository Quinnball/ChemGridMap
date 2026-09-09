"""Download reproducible target-specific ChEMBL activity tables."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import hashlib
from pathlib import Path
import time
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd


API_ROOT = "https://www.ebi.ac.uk/chembl/api/data"


def fetch_json(url: str, attempts: int = 4) -> dict:
    last_error = None
    for attempt in range(attempts):
        try:
            request = Request(url, headers={"User-Agent": "ChemGridMap/0.3"})
            with urlopen(request, timeout=90) as response:
                return json.load(response)
        except Exception as error:  # Network retries are recorded in the manifest.
            last_error = error
            if attempt + 1 < attempts:
                time.sleep(2 ** attempt)
    raise RuntimeError("Failed to download {}: {}".format(url, last_error))


def download_target(
    target_id: str,
    output_root: Path,
    activity_type: str = "IC50",
    page_size: int = 1000,
) -> dict:
    target_id = str(target_id).strip().upper()
    status = fetch_json(API_ROOT + "/status.json")
    release = str(status["chembl_db_version"]).lower().replace("_", "")
    target_dir = output_root / "{}_{}".format(release, target_id.lower())
    target_dir.mkdir(parents=True, exist_ok=True)
    (target_dir / "chembl_status.json").write_text(json.dumps(status, indent=2))

    target_url = "{}/target/{}.json".format(API_ROOT, target_id)
    target_metadata = fetch_json(target_url)
    (target_dir / "target_metadata.json").write_text(
        json.dumps(target_metadata, indent=2),
        encoding="utf-8",
    )

    activities = []
    page_urls = []
    offset = 0
    total_count = None
    while total_count is None or offset < total_count:
        query = urlencode(
            {
                "target_chembl_id__exact": target_id,
                "standard_type__exact": activity_type,
                "limit": page_size,
                "offset": offset,
            }
        )
        url = "{}/activity.json?{}".format(API_ROOT, query)
        payload = fetch_json(url)
        page_path = target_dir / "activity_offset_{:04d}.json".format(offset)
        page_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        rows = payload.get("activities", [])
        activities.extend(rows)
        page_urls.append(url)
        total_count = int(payload.get("page_meta", {}).get("total_count", len(rows)))
        if not rows:
            raise RuntimeError("Empty page before the advertised record count was reached.")
        offset += len(rows)

    if fetch_json(API_ROOT + "/status.json")["chembl_db_version"] != status["chembl_db_version"]:
        raise RuntimeError("ChEMBL changed release during download; rerun the download.")
    if len(activities) != total_count:
        raise RuntimeError("Downloaded row count differs from the API total.")

    table = pd.json_normalize(activities, sep=".")
    csv_path = target_dir / "{}_{}_{}_raw.csv".format(
        release,
        target_id.lower(),
        activity_type.lower(),
    )
    table.to_csv(csv_path, index=False)
    manifest = {
        "downloaded_at_utc": datetime.now(timezone.utc).isoformat(),
        "api_root": API_ROOT,
        "chembl_status": status,
        "target_id": target_id,
        "target_pref_name": target_metadata.get("pref_name"),
        "activity_type": activity_type,
        "raw_activity_rows": int(len(table)),
        "reported_total_count": total_count,
        "page_size": int(page_size),
        "page_urls": page_urls,
        "raw_csv": str(csv_path),
        "raw_csv_sha256": hashlib.sha256(csv_path.read_bytes()).hexdigest(),
    }
    (target_dir / "download_manifest.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("targets", nargs="+", help="Target ChEMBL IDs.")
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("validation_data"),
    )
    parser.add_argument("--activity-type", default="IC50")
    args = parser.parse_args()

    for target in args.targets:
        manifest = download_target(
            target,
            output_root=args.output_root,
            activity_type=args.activity_type,
        )
        print(
            "{target_id}: {raw_activity_rows} {activity_type} records".format(
                **manifest
            )
        )


if __name__ == "__main__":
    main()
