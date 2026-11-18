"""Reliability revision entry point: independently rerunnable analysis stages."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import tempfile

# Fixed thread counts make measured runs comparable and keep UMAP deterministic.
for variable in ("NUMBA_NUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[variable] = "1"
os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "chemgridmap-t2a-mpl"))
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
import sys
sys.dont_write_bytecode = True

from reliability_analysis import OUTPUT, run_layouts, run_neighborhoods, run_stability, write_manifest

STEPS = ("layouts", "neighborhoods", "stability", "baselines", "scale", "activity", "cores", "cores-timing",
         "cases", "figure-data", "results")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps", choices=STEPS, nargs="+", default=list(STEPS))
    parser.add_argument("--workers", type=int, default=2,
                        help="Workers for non-timing baseline seeds; seed 42 always runs sequentially")
    parser.add_argument("--baseline-phase", choices=("all", "timing", "remaining", "metrics-only"),
                        default="all", help="metrics-only reads frozen checkpoints without layout construction or timing")
    parser.add_argument("--resume-baselines", action="store_true",
                        help="Reuse baseline checkpoints only when their input/code signature matches")
    args = parser.parse_args()
    for step in args.steps:
        print(f"=== {step} ===", flush=True)
        if step == "layouts":
            run_layouts(OUTPUT)
        elif step == "neighborhoods":
            run_neighborhoods(OUTPUT)
        elif step == "stability":
            run_stability(OUTPUT)
        elif step == "baselines":
            from reliability_baselines import run_baselines
            run_baselines(OUTPUT, workers=args.workers, phase=args.baseline_phase,
                          resume=args.resume_baselines)
        elif step == "scale":
            from reliability_scale import run_scale, SCALE_INPUTS
            result = run_scale(OUTPUT)
            write_manifest(OUTPUT / "scale", SCALE_INPUTS, {"seed": 42, "n": 9581,
                           "summary": "full_herg_summary.json", "runtime": "sequential"})
        elif step == "activity":
            from reliability_activity import run_activity
            run_activity(OUTPUT)
            from reliability_contrasts import run_contrasts
            run_contrasts(OUTPUT)
        elif step == "cores":
            from reliability_cores import run_cores
            run_cores(OUTPUT)
        elif step == "cores-timing":
            from reliability_cores import run_cores_timing
            from reliability_t2b import update_rework_outputs, refresh_rework_manifests
            run_cores_timing(OUTPUT)
            update_rework_outputs(OUTPUT)
            refresh_rework_manifests(OUTPUT)
        elif step == "cases":
            if not (OUTPUT / "cases/case_summary.json").exists():
                from reliability_cases import run_cases
                result = run_cases(OUTPUT)
                write_manifest(OUTPUT / "cases", result["inputs"], result["parameters"])
            if (OUTPUT / "activity/adoption.json").exists():
                from reliability_t2b import run_activity_case
                run_activity_case(OUTPUT)
        elif step == "figure-data":
            if (OUTPUT / "activity/adoption.json").exists():
                from reliability_t2b import run_t2b_figure_data
                run_t2b_figure_data(OUTPUT)
            else:
                from reliability_figure_data import run_figure_data
                run_figure_data(OUTPUT)
        elif step == "results":
            from reliability_t2b import run_results
            run_results(OUTPUT)
        print(f"completed: {step}", flush=True)


if __name__ == "__main__":
    main()
