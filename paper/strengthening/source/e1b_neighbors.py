"""E1-B 邻居 SAR：只读取冻结布局，按任务单的固定规则估计相关和误差。"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[3]
CGM = ROOT
S = ROOT / "paper/strengthening"
OUT = S / "output/E1B"
sys.path[:0] = [str(CGM / "paper"), str(CGM / "src")]

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from reliability_analysis import FROZEN, OUTPUT, TARGETS, SEEDS, LAYOUTS, load_target, layout_points
from chemgridmap.reliability import (
    ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE, grid_neighbor_sets,
    representation_neighbor_sets, neighbor_recall,
)

BASELINES = ("projection_lap", "projection_lap_refined", "las_random",
             "las_projection", "flas_random", "random")
K, BOOT_N, BOOT_SEED = 10, 1000, 20261031
BIN_EDGES = np.array([0., .2, .4, .6, .8, 1.])
BIN_LABELS = ("[0,0.2)", "[0.2,0.4)", "[0.4,0.6)", "[0.6,0.8)", "[0.8,1]")
FROZEN_METRICS = OUTPUT / "neighborhoods/layout_metrics.csv"
FROZEN_RECALL = OUTPUT / "neighborhoods/per_molecule_recall_seed42.csv"


def neighborhood_metrics(activity, reference, neighbors):
    """以每个分子为单位；每组邻居等权，保留并列，不补齐或截断邻居。"""
    activity = np.asarray(activity, dtype=float)
    if len(activity) != len(reference) or len(activity) != len(neighbors):
        raise ValueError("活动值与邻居列表长度不同")
    if not np.isfinite(activity).all():
        raise ValueError("显示活性包含缺失或非有限值")
    e_fp = np.array([np.abs(activity[group] - activity[i]).mean()
                     for i, group in enumerate(reference)])
    e_grid = np.array([np.abs(activity[group] - activity[i]).mean()
                       for i, group in enumerate(neighbors)])
    return pd.DataFrame({
        "recall": neighbor_recall(reference, neighbors),
        "e_grid": e_grid, "e_fp": e_fp, "excess_error": e_grid - e_fp,
        "fingerprint_neighbor_count": [len(group) for group in reference],
        "grid_neighbor_count": [len(group) for group in neighbors],
        "shared_neighbor_count": [len(np.intersect1d(a, b, assume_unique=True))
                                  for a, b in zip(reference, neighbors)],
    })


def spearman_rho(x, y):
    """平均秩的 Pearson 相关；不计算或解释显著性。"""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if len(x) != len(y) or len(x) < 2 or not (np.isfinite(x).all() and np.isfinite(y).all()):
        return float("nan")
    rx, ry = rankdata(x), rankdata(y)
    rx, ry = rx - rx.mean(), ry - ry.mean()
    denominator = np.sqrt(np.dot(rx, rx) * np.dot(ry, ry))
    return float(np.dot(rx, ry) / denominator) if denominator else float("nan")


def bootstrap_spearman(x, y, n_boot=BOOT_N, seed=BOOT_SEED):
    """分子成对重抽样，邻居保持冻结；每个靶点/布局重启同一个固定 RNG。"""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    rng = np.random.default_rng(seed)
    draws = np.empty(n_boot)
    for iteration in range(n_boot):
        indices = rng.integers(0, len(x), len(x))
        draws[iteration] = spearman_rho(x[indices], y[indices])
    valid = draws[np.isfinite(draws)]
    low, high = np.quantile(valid, [.025, .975]) if len(valid) else (np.nan, np.nan)
    return {"rho_ci_low": float(low), "rho_ci_high": float(high),
            "bootstrap_valid_count": int(len(valid)), "bootstrap_n": int(n_boot)}


def paired_mean_ci(differences, n_boot=BOOT_N, seed=BOOT_SEED):
    """成对 seed 差值的均值；按 seed 有放回抽样的百分位 95% CI。"""
    differences = np.asarray(differences, dtype=float)
    if not len(differences) or not np.isfinite(differences).all():
        raise ValueError("成对 seed 差值必须非空且有限")
    rng = np.random.default_rng(seed)
    draws = differences[rng.integers(0, len(differences), (n_boot, len(differences)))].mean(axis=1)
    low, high = np.quantile(draws, [.025, .975])
    return {"paired_difference_mean": float(differences.mean()),
            "paired_ci_low": float(low), "paired_ci_high": float(high),
            "bootstrap_n": int(n_boot), "n_seeds": int(len(differences))}


def recall_bins(metrics):
    """左闭右开，最后一箱含 1；空箱 n=0，均值保留 NA。"""
    recalls = metrics.recall.to_numpy(float)
    if not ((recalls >= 0) & (recalls <= 1)).all():
        raise ValueError("recall 必须属于 [0,1]")
    assignments = np.searchsorted(BIN_EDGES[1:-1], recalls, side="right")
    rows = []
    for index, label in enumerate(BIN_LABELS):
        group = metrics.loc[assignments == index]
        rows.append({"recall_bin": label, "n": len(group),
                     "mean_e_grid": group.e_grid.mean(), "mean_e_fp": group.e_fp.mean(),
                     "mean_excess_error": group.excess_error.mean()})
    return rows


def check_identity(expected, observed, description):
    """不依赖隐式行位置；身份唯一且严格同序后才使用冻结坐标。"""
    expected, observed = np.asarray(expected, dtype=str), np.asarray(observed, dtype=str)
    if len(np.unique(expected)) != len(expected) or not np.array_equal(expected, observed):
        raise ValueError(f"分子身份缺失、重复或排序不一致：{description}")


def check_frozen_layout_metrics(seed_metrics, frozen, tolerance=1e-12):
    """逐 target/seed/layout 核对所有主布局的 e_grid、e_fp 与平均 recall。"""
    key = ["target", "seed", "layout"]
    actual = seed_metrics.loc[seed_metrics.layout.isin(LAYOUTS)].set_index(key).sort_index()
    expected = frozen.set_index(key).sort_index()
    if not actual.index.is_unique or not expected.index.is_unique or not actual.index.equals(expected.index):
        raise ValueError("冻结布局指标的键不一致")
    result = {"passed": True, "rows": len(actual), "absolute_tolerance": tolerance}
    for new, old in (("mean_e_grid", "neighbor_abs_dpchembl"),
                     ("mean_e_fp", "fingerprint_neighbor_abs_dpchembl"),
                     ("mean_recall", "recall_k10")):
        difference = np.abs(actual[new].to_numpy(float) - expected[old].to_numpy(float))
        maximum = float(difference.max())
        if not np.isfinite(difference).all() or maximum > tolerance:
            raise ValueError(f"冻结 {old} 不一致，最大误差 {maximum}")
        result[f"{old}_max_absolute_difference"] = maximum
    return result


def input_paths():
    """显式白名单；不扫描或读取 strengthening/data。"""
    paths = [FROZEN_METRICS, FROZEN_RECALL, CGM / "paper/reliability_analysis.py",
             CGM / "src/chemgridmap/reliability.py", CGM / "src/chemgridmap/representations.py"]
    for target in TARGETS:
        paths.append(FROZEN / target / "map_molecules.csv")
        paths.extend(OUTPUT / "layouts" / target / f"coordinates_seed_{seed}.csv" for seed in SEEDS)
        paths.append(OUTPUT / "baselines/checkpoints" / f"{target}_occupancy0.4_seed42.npz")
    return paths


def run(output=OUT):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    seed_rows, molecule_frames, bin_rows = [], [], []
    bootstraps, identities, matching_baselines = {}, [], []

    def evaluate(target, seed, layout, family, points, data, reference):
        if points.shape != (len(data), 2) or not np.isfinite(points).all():
            raise ValueError(f"无效坐标：{target}/{seed}/{layout}")
        metrics = neighborhood_metrics(data.activity_pchembl, reference, grid_neighbor_sets(points, k=K))
        metadata = {"target": target, "target_label": TARGETS[target], "seed": seed,
                    "layout": layout, "layout_family": family, "occupancy": .4}
        seed_rows.append({**metadata, "row_type": "seed_layout", "n": len(data),
                          "mean_recall": metrics.recall.mean(), "mean_e_grid": metrics.e_grid.mean(),
                          "mean_e_fp": metrics.e_fp.mean(), "mean_excess_error": metrics.excess_error.mean(),
                          "rho": spearman_rho(metrics.recall, metrics.excess_error)})
        bin_rows.extend({**metadata, **row} for row in recall_bins(metrics))
        if seed == 42:
            molecules = metrics.assign(**metadata, molecule_identity_key=data.molecule_identity_key.to_numpy(),
                                       activity_pchembl=data.activity_pchembl.to_numpy(float))
            molecule_frames.append(molecules)
            bootstraps[target, layout] = bootstrap_spearman(metrics.recall, metrics.excess_error)

    for target in TARGETS:
        data, fingerprints = load_target(target)
        check_identity(data.molecule_identity_key, data.molecule_identity_key, target)
        if not np.allclose(data.activity_pchembl, data.pchembl_median, rtol=0., atol=1e-12):
            raise ValueError(f"显示活性不是冻结 pChEMBL 中位数：{target}")
        reference = representation_neighbor_sets(fingerprints, k=K)
        for seed in SEEDS:
            path = OUTPUT / "layouts" / target / f"coordinates_seed_{seed}.csv"
            frame = pd.read_csv(path)
            check_identity(data.molecule_identity_key, frame.molecule_identity_key, path.name)
            if not np.allclose(data.activity_pchembl, frame.activity_pchembl, rtol=0., atol=1e-12):
                raise ValueError(f"坐标文件活性与映射活性不一致：{path}")
            identities.append({"target": target, "seed": seed, "kind": "main", "n": len(data)})
            if seed == 42:
                seed42_frame = frame
            for layout in LAYOUTS:
                evaluate(target, seed, layout, "main", layout_points(frame, layout), data, reference)
        path = OUTPUT / "baselines/checkpoints" / f"{target}_occupancy0.4_seed42.npz"
        with np.load(path, allow_pickle=False) as checkpoint:
            check_identity(data.molecule_identity_key, checkpoint["molecule_identity_key"], path.name)
            identities.append({"target": target, "seed": 42, "kind": "baseline", "n": len(data)})
            for baseline, main in (("projection_lap", "minimum_displacement"),
                                   ("projection_lap_refined", "neighbor_refined")):
                if not np.array_equal(checkpoint[baseline], layout_points(seed42_frame, main)):
                    raise ValueError(f"主布局与 occupancy 0.4 基线坐标不一致：{target}/{baseline}")
                matching_baselines.append({"target": target, "baseline": baseline,
                                           "main": main, "seed": 42, "exact_match": True})
            for layout in BASELINES:
                evaluate(target, 42, layout, "baseline", checkpoint[layout], data, reference)
        print(f"NEIGHBOR SAR {target}: {len(data)} molecules; 30 main and 6 baseline layouts complete", flush=True)

    seed_metrics = pd.DataFrame(seed_rows)
    regression = check_frozen_layout_metrics(seed_metrics, pd.read_csv(FROZEN_METRICS))
    molecules = pd.concat(molecule_frames, ignore_index=True)
    key = ["target", "seed", "layout", "molecule_identity_key"]
    observed = molecules.loc[molecules.layout_family == "main"].set_index(key).sort_index()
    expected = pd.read_csv(FROZEN_RECALL).set_index(key).sort_index()
    if not observed.index.is_unique or not expected.index.is_unique or not observed.index.equals(expected.index):
        raise ValueError("冻结 seed42 单分子召回率身份不一致")
    recall_maximum = float(np.abs(observed.recall - expected.recall_k10).max())
    if recall_maximum > 1e-12:
        raise ValueError(f"冻结单分子 recall 不一致：{recall_maximum}")
    regression.update({"per_molecule_recall_rows": len(observed),
                       "per_molecule_recall_max_absolute_difference": recall_maximum})

    summaries, paired = [], []
    for (target, layout), group in seed_metrics.groupby(["target", "layout"], sort=False):
        seed42 = group.loc[group.seed == 42].iloc[0]
        summaries.append({"target": target, "target_label": TARGETS[target], "layout": layout,
                          "layout_family": seed42.layout_family, "occupancy": .4, "n": int(seed42.n),
                          "seed_for_ci": 42, "rho_seed42": float(seed42.rho),
                          **bootstraps[target, layout], "n_seeds": len(group),
                          "rho_seed_median": float(group.rho.median()), "rho_seed_min": float(group.rho.min()),
                          "rho_seed_max": float(group.rho.max())})
    for target in TARGETS:
        selected = seed_metrics.loc[(seed_metrics.target == target) & seed_metrics.layout.isin(LAYOUTS)]
        pivot = selected.pivot(index="seed", columns="layout", values="mean_e_grid").reindex(SEEDS)
        differences = pivot.neighbor_refined - pivot.minimum_displacement
        paired.append({"target": target, "target_label": TARGETS[target], "row_type": "paired_seed_difference",
                       "layout": "neighbor_refined-minus-minimum_displacement", "layout_family": "main",
                       "occupancy": .4, **paired_mean_ci(differences.to_numpy()),
                       "ci_method": "paired-seed bootstrap percentile 95%", "bootstrap_seed": BOOT_SEED})
    errors = pd.concat([seed_metrics, pd.DataFrame(paired)], ignore_index=True)
    tables = {"neighbor_sar_by_molecule_seed42.csv": molecules,
              "neighbor_sar_summary.csv": pd.DataFrame(summaries),
              "neighbor_sar_bins.csv": pd.DataFrame(bin_rows), "layout_activity_error.csv": errors}
    for filename, frame in tables.items():
        frame.to_csv(output / filename, index=False)
    checks = {
        "frozen_regression": regression, "identity_checks": identities,
        "main_baseline_coordinate_checks": matching_baselines,
        "inputs_sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                          for path in input_paths()},
        "parameters": {"targets": list(TARGETS), "main_layouts": list(LAYOUTS),
                       "main_seeds": list(SEEDS), "baseline_layouts": list(BASELINES),
                       "baseline_seed": 42, "occupancy": .4, "k": K,
                       "ties": {"relative_tolerance": RELATIVE_TOLERANCE, "absolute_tolerance": ABSOLUTE_TOLERANCE},
                       "fingerprint": "frozen load_target: Morgan radius 2, 2048 bits; binary Tanimoto",
                       "rho": "Spearman average ranks: recall versus e_grid - e_fp",
                       "bootstrap_n": BOOT_N, "bootstrap_seed": BOOT_SEED,
                       "rho_ci": "molecule-pair bootstrap, fixed neighborhoods, percentile 95%; seed 42 only",
                       "paired_ci": "mean of 10 seed differences; paired-seed bootstrap percentile 95%",
                       "rng": "numpy.random.default_rng, reinitialized per target/layout or target contrast",
                       "bins": list(BIN_LABELS), "significance_test": "none; effect estimation only"},
        "table_rows": {filename: len(frame) for filename, frame in tables.items()},
        "paired_differences": paired,
        "uncertainty_limit": "Molecule bootstrap resamples observed pairs with fixed overlapping neighborhoods; it does not account for neighborhood dependence. Seed contrasts describe layout-seed variability for these fixed targets.",
    }
    (output / "neighbor_checks.json").write_text(json.dumps(checks, ensure_ascii=False, indent=2,
                                                            allow_nan=False) + "\n")
    print("FROZEN REGRESSION PASS " + json.dumps(regression), flush=True)
    print("OUTPUT ROWS " + json.dumps(checks["table_rows"]), flush=True)
    print("PAIRED DIFFERENCES " + json.dumps(paired), flush=True)
    return checks


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUT)
    run(parser.parse_args().output)
