from __future__ import annotations

import json
from pathlib import Path

from manifest_utils import write_manifest


def test_write_manifest_fields_and_relative_paths(tmp_path):
    root = Path(__file__).resolve().parents[3]
    input_path = root / "paper/strengthening/source/test_manifest_utils.py"
    source_path = root / "paper/strengthening/source/manifest_utils.py"
    out_dir = tmp_path / "output"
    out_dir.mkdir()
    (out_dir / "result.csv").write_text("value\n1\n", encoding="utf-8")

    manifest = write_manifest(
        out_dir,
        stage="unit-test",
        inputs=[input_path],
        sources=[source_path],
        parameters={"analysis_role": "test", "seed": 42},
        started_at="2026-11-01T00:00:00+00:00",
        finished_at="2026-11-01T00:00:02+00:00",
    )

    assert manifest["stage"] == "unit-test"
    assert manifest["analysis_role"] == "test"
    assert manifest["elapsed_seconds"] == 2.0
    assert set(manifest["software_versions"]) >= {
        "python",
        "numpy",
        "pandas",
        "scipy",
        "scikit-learn",
        "umap-learn",
        "pynndescent",
        "numba",
        "rdkit",
        "chemgridmap",
    }
    assert manifest["outputs_sha256"]
    text = (out_dir / "manifest.json").read_text(encoding="utf-8")
    assert str(root) not in text
    assert str(tmp_path) not in text
    assert "task_sha256" not in json.loads(text)
