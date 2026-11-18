"""Restore hash-checked manuscript inputs without replacing existing local files."""
import argparse
import hashlib
import json
from pathlib import Path
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = 'chemgridmap-0.3.0-validation.zip'
SHA256 = 'b8b51f5893fc035825a2f0169fdde31d5d805dd4a39d342615ffb4cfe08e1ba3'
URL = 'https://github.com/Quinnball/ChemGridMap/releases/download/v0.3.0/' + ARCHIVE


def restore(archive, root):
    if hashlib.sha256(archive.read_bytes()).hexdigest() != SHA256:
        raise ValueError('Unexpected validation archive digest; no files were extracted.')
    with zipfile.ZipFile(archive) as bundle:
        hashes = json.loads(bundle.read('SHA256.json'))
        files = []
        for name, digest in hashes.items():
            if name == 'validation_data/README.md':
                continue
            path = (root/name).resolve()
            if not path.is_relative_to(root.resolve()):
                raise ValueError('Archive path leaves the repository.')
            content = bundle.read(name)
            if hashlib.sha256(content).hexdigest() != digest:
                raise ValueError(f'Archive member failed its digest check: {name}')
            if path.exists() and hashlib.sha256(path.read_bytes()).hexdigest() != digest:
                raise ValueError(f'Existing file differs: {path}. Preserve it and restore into a clean checkout instead.')
            files.append((path,content))
        # Validate every member before adding any missing file.
        for path, content in files:
            if not path.exists():
                path.parent.mkdir(parents=True,exist_ok=True)
                path.write_bytes(content)
    print(f'Verified {len(files)} archived files; existing files were not replaced.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--download',action='store_true',help='Fetch the public v0.3.0 data archive if absent.')
    parser.add_argument('--root',type=Path,default=ROOT)
    arguments = parser.parse_args()
    archive = arguments.root/'validation_data'/ARCHIVE
    if not archive.exists():
        if not arguments.download:
            raise SystemExit(f'Download {URL} into {archive.parent}, or rerun with --download.')
        archive.parent.mkdir(parents=True,exist_ok=True)
        with urllib.request.urlopen(URL,timeout=120) as response:
            content = response.read()
        if hashlib.sha256(content).hexdigest() != SHA256:
            raise ValueError('Downloaded archive digest mismatch; not saved.')
        archive.write_bytes(content)
    restore(archive,arguments.root)
