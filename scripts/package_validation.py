"""Bundle public ChEMBL inputs and saved numeric evidence, without manuscript drafts."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=ROOT / "paper/output/revision_v6")
    parser.add_argument("--output", type=Path, default=ROOT / "validation_data/chemgridmap-0.3.0-validation.zip")
    args = parser.parse_args()
    source = args.source
    files = list(source.glob("*.csv")) + [source / name for name in (
        "record_tasks.json", "protocol.json", "selected_windows.json", "large_scale_extension.json")]
    for target in ("chembl205", "chembl204", "chembl240"):
        files.extend((source / target).glob("*.csv"))
        files.extend((source / target).glob("*curation_report.json"))
    files.extend((source / "figures").glob("Figure_*.svg"))
    web = source / "web_export_map/chembl205_web_raw.csv"
    files.append(web)
    members = {str(path.relative_to(ROOT)): path for path in files}
    for name in ("chembl37_chembl205", "chembl37_chembl204", "chembl37_chembl240_large_scale"):
        for path in (ROOT / "validation_data" / name).glob("*.csv"):
            if path.name.endswith("_ic50_raw.csv"):
                members[str(path.relative_to(ROOT))] = path
    for path in (ROOT / "validation_data/assay_context").glob("*.json"):
        members[str(path.relative_to(ROOT))] = path
    members["validation_data/README.md"] = ROOT / "validation_data/README.md"
    digests = {}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(args.output, "w", zipfile.ZIP_DEFLATED) as archive:
        for name, path in sorted(members.items()):
            archive.write(path, name)
            digests[name] = hashlib.sha256(path.read_bytes()).hexdigest()
        archive.writestr("SHA256.json", json.dumps(digests, indent=2) + "\n")
    print(json.dumps({"archive": str(args.output), "files": len(members),
                      "bytes": args.output.stat().st_size,
                      "sha256": hashlib.sha256(args.output.read_bytes()).hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
