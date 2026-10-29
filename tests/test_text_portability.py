"""Generated text exports must not depend on the host newline convention."""

from pathlib import Path

import pandas as pd

from chemgridmap.pipeline import build_grid_map


def test_generated_svg_csv_and_json_outputs_use_lf(tmp_path):
    data = pd.DataFrame({
        "canonical_smiles": ["CCO", "CCN", "CCC"],
        "x": [0.0, 1.0, 2.0],
        "y": [0.0, 2.0, 1.0],
        "activity_class": ["活性", "中等", "非活性"],
    })
    result = build_grid_map(
        data,
        output_dir=tmp_path,
        name="portable",
        smiles_col="canonical_smiles",
        x_col="x",
        y_col="y",
        grouping="none",
        render_detail="overview",
        output_formats=["svg"],
    )

    outputs = {
        Path(result.output_files["coordinates"]),
        Path(result.output_files["metrics"]),
        Path(result.output_files["run_manifest"]),
        Path(result.output_files["svg"]),
    }
    assert {path.suffix for path in outputs} == {".csv", ".json", ".svg"}
    for path in outputs:
        contents = path.read_bytes()
        assert b"\r\n" not in contents, path.name
        contents.decode("utf-8")
