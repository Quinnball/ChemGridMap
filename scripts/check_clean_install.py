"""Check an independently installed wheel using the real ChEMBL web CSV."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python", required=True, type=Path)
    parser.add_argument("--output", type=Path, default=ROOT / "paper/output/revision_v5")
    args = parser.parse_args()
    out = args.output.resolve()
    raw = out / "web_export_map/chembl205_web_raw.csv"
    env = dict(os.environ, MPLCONFIGDIR=str(ROOT / ".cache/matplotlib"), NUMBA_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    logs = out / "clean_environment"
    logs.mkdir(exist_ok=True)
    commands = []
    with tempfile.TemporaryDirectory(prefix="chemgridmap-installed-check-") as directory:
        def run(name, arguments):
            command = [str(args.python), *arguments]
            result = subprocess.run(command, cwd=directory, env=env, text=True, capture_output=True)
            (logs / (name + ".log")).write_text(result.stdout + result.stderr, encoding="utf-8")
            commands.append({"name": name, "command": command, "returncode": result.returncode})
            if result.returncode:
                raise RuntimeError(f"{name} failed; inspect {logs / (name + '.log')}")
            return result.stdout
        installation = json.loads(run("installed_package", ["-c",
            "import json,sys,chemgridmap; from pathlib import Path; "
            "print(json.dumps(dict(prefix=sys.prefix,base_prefix=sys.base_prefix,path=chemgridmap.__file__,"
            "version=chemgridmap.__version__,config=(Path(sys.prefix)/'pyvenv.cfg').read_text())))"]))
        assert installation["prefix"] != installation["base_prefix"]
        assert "include-system-site-packages = false" in installation["config"]
        assert Path(installation["path"]).is_relative_to(installation["prefix"])
        assert not Path(installation["path"]).is_relative_to(ROOT)
        run("pip_check", ["-m", "pip", "check"])
        run("pytest", ["-m", "pytest", "-q", "-o", "pythonpath=", str(ROOT / "tests")])
        run("pip_freeze", ["-m", "pip", "freeze"])
        mapped = logs / "web_map"
        run("web_csv_to_map", ["-m", "chemgridmap", str(raw), "--input-format", "chembl",
            "--target-id", "CHEMBL205", "--activity-type", "IC50", "--representation", "morgan",
            "--projection", "umap", "--output-dir", str(mapped), "--name", "web"])
        run("inspect", ["-m", "chemgridmap", "inspect", str(mapped / "web_grid_coordinates.csv"),
            "--records", str(mapped / "web_retained_records.csv"), "--molecule-id", "CHEMBL20",
            "--output-dir", str(logs / "record_audit")])

    fresh = pd.read_csv(logs / "web_map/web_grid_coordinates.csv").sort_values("molecule_identity_key")
    original = pd.read_csv(out / "web_export_map/chembl205_web_grid_coordinates.csv").sort_values("molecule_identity_key")
    assert fresh.molecule_identity_key.tolist() == original.molecule_identity_key.tolist()
    comparison = {"molecular_identities_match": True,
                  "raw_projection_exact_match": bool(np.array_equal(fresh[["projection_x_raw", "projection_y_raw"]], original[["projection_x_raw", "projection_y_raw"]])),
                  "grid_cells_exact_match": bool(np.array_equal(fresh[["grid_row", "grid_col"]], original[["grid_row", "grid_col"]]))}
    assert all(comparison.values())
    web = pd.read_csv(raw, sep=";")
    retained = pd.read_csv(logs / "web_map/web_retained_records.csv")
    archive = pd.read_csv(out / "chembl205/map_molecules.csv")
    left = fresh.set_index("canonical_smiles").activity_pchembl.sort_index()
    right = archive.set_index("canonical_smiles").activity_pchembl.sort_index()
    assert left.index.equals(right.index) and np.allclose(left, right)
    report = json.loads((logs / "web_map/web_curation_report.json").read_text())
    audit = json.loads((logs / "record_audit/audit_summary.json").read_text())
    result = {"passed": True, "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "installation": installation, "commands": commands, "repeat_comparison": comparison,
        "web_download": {"source": "ChEMBL Activities web interface CSV download",
            "query": "target_chembl_id:CHEMBL205 AND standard_type:IC50",
            "generated_url": "https://www.ebi.ac.uk/chembl/interface_api/delayed_jobs/outputs/DOWNLOAD-tob8Q5KByo6XoEv_6e-_t66enBuYgCQnTwKypmifomQ_eq_/DOWNLOAD-tob8Q5KByo6XoEv_6e-_t66enBuYgCQnTwKypmifomQ_eq_.zip",
            "sha256": hashlib.sha256(raw.read_bytes()).hexdigest(),
            "raw_records": len(web), "retained_records": len(retained), "molecules": len(fresh),
            "parent_identifiers_present": "Parent Molecule ChEMBL ID" in web,
            "activity_identifiers_present": "Activity ID" in web,
            "canonical_structures_and_medians_match_api_archive": True,
            "warnings": report["warnings"], "cell_audit": audit},
        "scope": "Fresh virtual environment, independently installed wheel and pinned dependencies on one macOS arm64 host; not a cross-platform claim."}
    (out / "clean_install_validation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
