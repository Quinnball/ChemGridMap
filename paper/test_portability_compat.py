"""Checks shared by the Python 3.10 and 3.12 CI environments."""

import hashlib

from relabel_chembl204 import file_hash


def test_file_hash_matches_standard_sha256(tmp_path):
    data = b"portable hash\x00with non-ASCII: \xe6\xb5\x8b\xe8\xaf\x95\n"
    path = tmp_path / "payload.bin"
    path.write_bytes(data)
    assert file_hash(path) == hashlib.sha256(data).hexdigest()
