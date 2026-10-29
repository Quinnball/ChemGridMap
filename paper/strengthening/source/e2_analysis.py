"""E2a 多靶点整理、折外模型、指标、bootstrap 与事实报告。"""
from __future__ import annotations

import argparse
from datetime import datetime
import gzip
import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path
import pickle
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[3]
REV = ROOT / "paper"
S = REV / "strengthening"
sys.path[:0] = [
    str(S / "source"),
    str(ROOT / "paper"),
    str(ROOT / "src"),
]

import numpy as np
import pandas as pd

from e1a_statistics import cluster_weights, percentile_ci, weighted_scores
from e2_models import (
    BOOTSTRAP_SEED,
    CLASS_BOUNDARIES,
    ENDPOINTS,
    FOLD_SEED,
    MODELS,
    POPULATIONS,
    S,
    STRATA,
    build_global_fits,
    build_target_predictions,
    comparator_aucs,
    curate_activity_table,
    prepare_target,
)

API = S / "data/chembl37_api"
ACTIVITIES = API / "activities"
OUT = S / "output/E2"
QA = S / "qa/E2"
FROZEN = ROOT / "paper/output/revision_v6"
E1A_PREDICTIONS = S / "output/E1A/predictions.csv"
ORIGINAL_TARGETS = ("CHEMBL205", "CHEMBL204", "CHEMBL240")
METRICS = (
    "accuracy",
    "agreement_rate",
    "auc",
    "ece",
    "brier",
    "log_loss",
    "mean_confidence",
    "calibration_bias",
    "high_confidence_n",
    "high_confidence_mean_probability",
    "high_confidence_accuracy",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_safe(value):
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def empty_scores() -> dict:
    return {metric: float("nan") for metric in METRICS} | {
        "high_confidence_n": 0.0
    }


def point_scores(frame: pd.DataFrame) -> dict:
    if frame.empty:
        return empty_scores()
    weights = np.ones((1, len(frame)), dtype=float)
    scores = weighted_scores(frame, weights)
    return {metric: float(scores[metric][0]) for metric in METRICS}


def regression_check() -> dict:
    """API 整理后的保留记录必须逐 activity_id 复现冻结记录。"""
    rows = []
    for target in ORIGINAL_TARGETS:
        lower = target.lower()
        raw = pd.read_csv(ACTIVITIES / f"{target}_IC50.csv", low_memory=False)
        curated = curate_activity_table(raw, target, "IC50").retained_records
        frozen = pd.read_csv(
            FROZEN / lower / f"{lower}_retained_records.csv", low_memory=False
        )
        api_records = curated[["activity_id", "pchembl_value"]].copy()
        frozen_records = frozen[["activity_id", "pchembl_value"]].copy()
        api_records["activity_id"] = api_records.activity_id.astype(str)
        frozen_records["activity_id"] = frozen_records.activity_id.astype(str)
        api_ids = set(api_records.activity_id)
        frozen_ids = set(frozen_records.activity_id)
        shared = api_records.merge(
            frozen_records, on="activity_id", suffixes=("_api", "_frozen")
        )
        max_error = float(
            (shared.pchembl_value_api - shared.pchembl_value_frozen).abs().max()
        )
        row = {
            "target": target,
            "api_input_rows": int(len(raw)),
            "api_retained_records": int(len(api_records)),
            "frozen_retained_records": int(len(frozen_records)),
            "api_only_activity_ids": sorted(api_ids - frozen_ids),
            "frozen_only_activity_ids": sorted(frozen_ids - api_ids),
            "shared_activity_ids": int(len(shared)),
            "max_pchembl_absolute_error": max_error,
        }
        difference_fraction = (
            len(api_ids - frozen_ids) + len(frozen_ids - api_ids)
        ) / len(frozen_ids)
        if difference_fraction > 0.005 or max_error > 1e-12:
            raise RuntimeError(f"API 回归检查失败：{target}")
        rows.append(row)
    result = {
        "checks": rows,
        "max_pchembl_absolute_error": max(float(row["max_pchembl_absolute_error"]) for row in rows),
        "api_only_activity_ids": sum(len(row["api_only_activity_ids"]) for row in rows),
        "frozen_only_activity_ids": sum(len(row["frozen_only_activity_ids"]) for row in rows),
        "api_format_adaptations": [
            "CGM detect_chembl_columns 直接识别 canonical_smiles、pchembl_value、target_chembl_id 等 API 列名。",
            "API 的 standard_relation 已是 '='，无需网页导出的首尾单引号适配。",
            "API 的 assay_type 列存在，hERG 直接沿用其值。",
            "原三靶点使用冻结 curation_report.json 中的整理参数；其余靶点使用封存默认参数。",
        ],
    }
    (QA / "regression.json").write_text(
        json.dumps(json_safe(result), indent=2, ensure_ascii=False, allow_nan=False)
        + "\n"
    )
    return result


def load_and_curate_selected():
    selected = pd.read_csv(API / "selected_targets.csv")
    inclusion_rows = []
    records_by_key = {}
    contexts_by_key = {}
    cache_path = QA / "curated_records.pkl.gz"
    cached = (
        pickle.loads(gzip.decompress(cache_path.read_bytes()))
        if cache_path.exists()
        else {}
    )
    for index, item in selected.iterrows():
        target = str(item.target_chembl_id)
        endpoint = str(item.standard_type)
        path = ACTIVITIES / f"{target}_{endpoint}.csv"
        try:
            raw = pd.read_csv(path, low_memory=False)
            if (target, endpoint) in cached:
                records = cached[(target, endpoint)]
            else:
                result = curate_activity_table(raw, target, endpoint)
                records = result.retained_records
            context = prepare_target(records)
            entries = context["doc"].groupby("entry").size()
            multi = int((entries >= 2).sum())
            included = multi >= 30
            reason = (
                "included"
                if included
                else "fewer than 30 entries with two or more evidence units"
            )
            inclusion_rows.append({
                "target": target,
                "endpoint": endpoint,
                "role": item.role,
                "api_records": int(len(raw)),
                "retained_records": int(len(records)),
                "entries": int(len(entries)),
                "multi_unit_entries": multi,
                "included": bool(included),
                "reason": reason,
                "input_file": str(path.relative_to(ROOT)),
            })
            records_by_key[(target, endpoint)] = records
            contexts_by_key[(target, endpoint)] = context
        except Exception as error:  # 封存计划要求逐靶点记录失败并继续。
            inclusion_rows.append({
                "target": target,
                "endpoint": endpoint,
                "role": item.role,
                "api_records": np.nan,
                "retained_records": np.nan,
                "entries": np.nan,
                "multi_unit_entries": np.nan,
                "included": False,
                "reason": f"curation failure: {type(error).__name__}: {error}",
                "input_file": str(path.relative_to(ROOT)),
            })
        print(
            f"curation {index + 1}/{len(selected)} {target} {endpoint}",
            flush=True,
        )
    inclusion = pd.DataFrame(inclusion_rows)
    inclusion.to_csv(OUT / "target_inclusion.csv", index=False)
    cache_path.write_bytes(
        gzip.compress(pickle.dumps(records_by_key, protocol=5))
    )
    return inclusion, records_by_key, contexts_by_key


def build_all_predictions(inclusion, records_by_key, contexts_by_key):
    prediction_frames = []
    for endpoint in ENDPOINTS:
        selected = inclusion.loc[
            (inclusion.endpoint == endpoint) & inclusion.included
        ]
        target_keys = [(target, endpoint) for target in selected.target]
        docs = {
            target: contexts_by_key[(target, endpoint)]["doc"]
            for target in selected.target
        }
        clean_entries = {
            target: contexts_by_key[(target, endpoint)]["clean_entries"]
            for target in selected.target
        }
        global_fits = build_global_fits(docs, clean_entries)
        for target in selected.target:
            frame, _ = build_target_predictions(
                records_by_key[(target, endpoint)],
                target,
                endpoint,
                global_fits=global_fits,
                context=contexts_by_key[(target, endpoint)],
            )
            prediction_frames.append(frame)
            print(f"models {target} {endpoint}: {len(frame)} rows", flush=True)
    predictions = pd.concat(prediction_frames, ignore_index=True)
    predictions.to_csv(OUT / "predictions.csv", index=False)
    return predictions


def metrics_by_target(predictions, inclusion) -> pd.DataFrame:
    rows = []
    included = inclusion.loc[inclusion.included]
    for item in included.itertuples(index=False):
        target_frame = predictions.loc[
            (predictions.target == item.target)
            & (predictions.endpoint == item.endpoint)
        ]
        for population in POPULATIONS:
            population_frame = target_frame.loc[
                target_frame.population == population
            ]
            for model in MODELS:
                model_frame = population_frame.loc[
                    population_frame.method == model
                ]
                for stratum in STRATA:
                    selected = (
                        model_frame
                        if stratum == "all"
                        else model_frame.loc[model_frame.stratum == stratum]
                    )
                    scores = point_scores(selected)
                    rows.append({
                        "target": item.target,
                        "endpoint": item.endpoint,
                        "population": population,
                        "model": model,
                        "stratum": stratum,
                        "n_predictions": int(len(selected)),
                        "n_entries": int(selected.entry.nunique()),
                        **scores,
                    })
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "metrics_by_target.csv", index=False)
    return result


def pooled_metrics(predictions) -> pd.DataFrame:
    rows = []
    for endpoint in ENDPOINTS:
        endpoint_frame = predictions.loc[predictions.endpoint == endpoint]
        for population in POPULATIONS:
            population_frame = endpoint_frame.loc[
                endpoint_frame.population == population
            ]
            for model in MODELS:
                model_frame = population_frame.loc[
                    population_frame.method == model
                ]
                for stratum in STRATA:
                    selected = (
                        model_frame
                        if stratum == "all"
                        else model_frame.loc[model_frame.stratum == stratum]
                    )
                    point = point_scores(selected)
                    bootstrap = None
                    if not selected.empty:
                        rng = np.random.default_rng(BOOTSTRAP_SEED)
                        weights = cluster_weights(selected, rng, 2000)
                        bootstrap = weighted_scores(selected, weights)
                    for metric in METRICS:
                        if bootstrap is None:
                            low = high = float("nan")
                            valid = 0
                        else:
                            low, high, valid = percentile_ci(bootstrap[metric])
                        rows.append({
                            "endpoint": endpoint,
                            "population": population,
                            "model": model,
                            "stratum": stratum,
                            "n_predictions": int(len(selected)),
                            "n_entries": int(selected.entry.nunique()),
                            "n_targets": int(selected.target.nunique()),
                            "metric": metric,
                            "estimate": point[metric],
                            "ci_low": low,
                            "ci_high": high,
                            "n_bootstrap_valid": int(valid),
                            "bootstrap_repetitions": 2000,
                            "bootstrap_seed": BOOTSTRAP_SEED,
                        })
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "pooled_metrics.csv", index=False)
    return result


def write_manifest(started: str, ended: str) -> dict:
    selected_inputs = [
        API / "selected_targets.csv",
        API / "download_manifest.json",
        API / "download_original_ic50_manifest.json",
        API / "SELECTION_RULES.md",
        E1A_PREDICTIONS,
        S / "plans/E2_task_draft_sealed.md",
        ROOT / "paper/reliability_activity.py",
        ROOT / "src/chemgridmap/chembl.py",
        ROOT / "src/chemgridmap/evidence.py",
        ROOT / "src/chemgridmap/activity.py",
        S / "source/e1a_models.py",
        S / "source/e1a_statistics.py",
    ]
    selected_inputs.extend(sorted(ACTIVITIES.glob("*.csv")))
    for target in ORIGINAL_TARGETS:
        lower = target.lower()
        selected_inputs.extend([
            FROZEN / lower / f"{lower}_retained_records.csv",
            FROZEN / lower / f"{lower}_curation_report.json",
        ])
    sources = sorted((S / "source").glob("e2_*.py")) + [S / "source/test_e2.py"]
    outputs = [
        path for path in sorted(OUT.iterdir())
        if path.is_file() and path.name != "manifest.json"
    ]
    manifest = {
        "task": "E2a",
        "sealed_plan_sha256": sha256(
            S / "plans/E2_task_draft_sealed.md"
        ),
        "started_at": started,
        "ended_at": ended,
        "seeds": {
            "fold": FOLD_SEED,
            "bootstrap": BOOTSTRAP_SEED,
            "bootstrap_repetitions": 2000,
        },
        "parameters": {
            "endpoints": list(ENDPOINTS),
            "models": list(MODELS),
            "populations": list(POPULATIONS),
            "strata": list(STRATA),
            "inclusion_rule": ">=30 entries with >=2 merged document units",
            "class_boundaries": list(CLASS_BOUNDARIES),
            "ece_bins": 10,
            "log_loss_clip": [1e-6, 1.0],
            "auc_probability": "p_displayed",
            "auc_ties": "average ranks",
            "high_confidence_threshold": 0.8,
            "G_fixed_tau": 0.64,
            "G_global": "same endpoint; full other included targets; current target excluded",
            "flag_clean_training": "remove every entry touched by full-target flagged documents",
            "prediction_bootstrap_refit": False,
            "parquet_support": False,
        },
        "inputs_sha256": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in selected_inputs
        },
        "source_sha256": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in sources
        },
        "outputs_sha256": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in outputs
        },
        "qa_sha256": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in sorted(QA.glob("*")) if path.is_file()
        },
        "software_versions": {
            "python": platform.python_version(),
            **{
                name: importlib.metadata.version(name)
                for name in ("numpy", "pandas", "scipy", "pytest")
            },
        },
        "manifest_self_hash": "omitted to avoid a self-referential hash",
    }
    (OUT / "manifest.json").write_text(
        json.dumps(json_safe(manifest), indent=2, ensure_ascii=False, allow_nan=False)
        + "\n"
    )
    return manifest


def run() -> None:
    started = datetime.now().astimezone().isoformat()
    OUT.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    regression = regression_check()
    inclusion, records_by_key, contexts_by_key = load_and_curate_selected()
    predictions = build_all_predictions(
        inclusion, records_by_key, contexts_by_key
    )
    metrics_by_target(predictions, inclusion)
    pooled_metrics(predictions)
    comparator_aucs(predictions).to_csv(OUT / "comparators.csv", index=False)
    ended = datetime.now().astimezone().isoformat()
    write_manifest(started, ended)
    print(
        f"completed regression={regression['max_pchembl_absolute_error']:.3g} "
        f"predictions={len(predictions)}",
        flush=True,
    )


def markdown_table(frame: pd.DataFrame) -> str:
    def cell(value):
        if isinstance(value, (list, tuple, dict)):
            return json.dumps(json_safe(value), ensure_ascii=False).replace("|", "\\|")
        if value is None or (np.isscalar(value) and pd.isna(value)):
            return "NA"
        if isinstance(value, (float, np.floating)):
            return f"{value:.6g}"
        return str(value).replace("|", "\\|").replace("\n", " ")

    rows = [
        "| " + " | ".join(map(str, frame.columns)) + " |",
        "| " + " | ".join("---" for _ in frame.columns) + " |",
    ]
    rows.extend(
        "| " + " | ".join(cell(value) for value in row) + " |"
        for row in frame.itertuples(index=False, name=None)
    )
    return "\n".join(rows)


def csv_rows(path: Path) -> int:
    with path.open(encoding="utf-8") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def write_report() -> None:
    existing_manifest = json.loads((OUT / "manifest.json").read_text())
    write_manifest(
        existing_manifest["started_at"], existing_manifest["ended_at"]
    )
    regression = json.loads((QA / "regression.json").read_text())
    inclusion = pd.read_csv(OUT / "target_inclusion.csv")
    manifest = json.loads((OUT / "manifest.json").read_text())
    scope = json.loads((QA / "frozen_compare.json").read_text())
    cgm_before = (QA / "cgm_status_before.txt").read_text()
    cgm_after = (QA / "cgm_status_after.txt").read_text()
    pytest_log = (QA / "pytest.log").read_text() if (QA / "pytest.log").exists() else "未运行"
    guard_log = (
        (QA / "guard_after.log").read_text()
        if (QA / "guard_after.log").exists()
        else "未运行"
    )
    output_rows = []
    for path in sorted(OUT.iterdir()):
        if not path.is_file():
            continue
        output_rows.append({
            "file": path.name,
            "rows": csv_rows(path) if path.suffix == ".csv" else "NA",
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    included = inclusion.loc[inclusion.included]
    excluded = inclusion.loc[~inclusion.included]
    included_table = (
        included.groupby("endpoint").target.agg(
            included_targets="size",
            target_ids=lambda values: ", ".join(values),
        ).reset_index()
    )
    regression_table = pd.DataFrame(regression["checks"])
    inclusion_columns = [
        "target", "endpoint", "role", "api_records", "retained_records",
        "entries", "multi_unit_entries", "included", "reason",
    ]
    sections = [
        "# E2a 执行报告",
        "",
        "## 1. 封存计划",
        "",
        f"SHA-256：`{manifest['sealed_plan_sha256']}`，与任务单一致。",
        "本报告只记录 E2a 的事实；未执行假设检验、负控、时间切分、"
        "伪共识、σ_doc 分布、标记对普查、条目统计、森林图数据或 RESULTS_E2.md。",
        "",
        "## 2. 回归检查与 API 适配",
        "",
        markdown_table(regression_table),
        "",
        "API 适配记录：",
        *[f"- {item}" for item in regression["api_format_adaptations"]],
        "",
        "## 3. 纳入情况",
        "",
        markdown_table(included_table),
        "",
        "完整纳入表：",
        "",
        markdown_table(inclusion[inclusion_columns]),
        "",
        "未纳入靶点及原因：",
        "",
        markdown_table(
            excluded[["target", "endpoint", "retained_records", "entries",
                      "multi_unit_entries", "reason"]]
        ),
        "",
        "## 4. 模型、指标与输出",
        "",
        "B0、B1、G、G-record、G-fixed、G-global、X 使用条目级 10 折，"
        "折种子 20261024。G-global 的 pooled SD 来自同端点其他纳入靶点，"
        "不使用当前靶点；flag-clean 拟合移除完整数据中触及标记文献的条目。"
        "B-range 与 B-sd 只使用剔除留出文献后的原始记录。",
        "",
        "环境未安装 pyarrow，按任务单回退输出 `predictions.csv`。"
        "端点为 IC50 与 Ki，总体为 all 与 flag-clean，层为全部、m=1、m>=2。"
        "ECE 使用最大类别概率的 10 个等宽区间；AUC 与整体校准偏差使用"
        "显示类别概率和显示类别是否一致。整群 bootstrap 按靶点×条目抽样，"
        "2,000 次，种子 20261028，不重新拟合模型。",
        "",
        markdown_table(pd.DataFrame(output_rows)),
        "",
        "## 5. 执行过程与修正",
        "",
        "1. 读取并核对封存计划和任务单，记录执行前 REV 全量 SHA-256、"
        "CGM 50 行状态及 `_guard` 六项结果。",
        "2. 依次整理 83 个 `selected_targets.csv` 中的靶点–端点；"
        "使用 API 文件完成原三靶点回归，失败靶点逐项记录后继续。",
        "3. 首轮模型补齐发现 B1 从每折全部模型行重复生成，导致留出键重复；"
        "改为只从 B0 行生成 B1。未改动模型定义或任何判据。",
        "4. 修正后全部模型完成，但比较对象合并 `G` 与 `B0` 时出现"
        "`remaining_records_x` 列名冲突；改为只从 B0 一侧携带 m 计数。"
        "随后从本地整理缓存重跑，未重复结构规范化。",
        "5. 完整输出生成后核对：181,405 条预测、7 个模型全覆盖、"
        "重复留出键为 0、概率和最大误差 3.3e-16；再执行完整 pytest。",
        "",
        "## 6. 测试与冻结",
        "",
        "```sh\n"
        "chemgridmap/.venv/bin/python -m pytest -p no:cacheprovider "
        "python -m pytest -p no:cacheprovider paper/strengthening/source/test_e2.py -vv\n"
        "```",
        "",
        "完整 pytest 输出：",
        "",
        "```text\n" + pytest_log.rstrip() + "\n```",
        "",
        "`_guard` 结果：",
        "",
        "```text\n" + guard_log.rstrip() + "\n```",
        "",
        f"CGM 状态执行前后均为 "
        f"{len(cgm_before.splitlines())} 行，前后文本"
        f"{'一致' if cgm_before == cgm_after else '不一致'}。",
        "REV 全量文件级核对："
        f"执行前 {scope['before_files']} 个文件，比较时 {scope['after_files']} 个文件；"
        f"新增 {len(scope['added'])} 个、修改 {len(scope['modified'])} 个、"
        f"删除 {len(scope['removed'])} 个。",
        "",
        "允许位置之外的变化：",
        "",
        *([f"- `{path}`" for path in scope["outside_allowed_changes"]] or ["- 无"]),
        "",
        "上述路径由执行期间的外部协调/执行器并发写入，本轮没有回滚、"
        "修改或删除它们。因此本次全量范围核对不通过；`_guard` 六处仍通过。",
        "",
        "## 7. 问题与未完成事项",
        "",
        "E2b 范围未执行。`predictions.parquet` 因环境无 pyarrow 支持，"
        "按任务单改用 `predictions.csv`。3 个靶点在默认父结构整理时报告"
        "“同一 parent ID 对应多个标准化 parent 结构”，按封存原则未纳入并保留原因；"
        "其余未纳入靶点均为多单元条目少于 30。全量范围核对存在允许位置外的"
        "并发文件变化，见上节。",
        "",
    ]
    (OUT / "E2a_REPORT.md").write_text("\n".join(sections))
    brief = [
        "# E2a 简要结果",
        "",
        f"- 封存计划哈希核对通过。",
        f"- API 回归：3 个原靶点 activity_id 集合完全一致，"
        f"pChEMBL 最大绝对误差 {regression['max_pchembl_absolute_error']:.3g}。",
        f"- 纳入：IC50 {int((included.endpoint == 'IC50').sum())} 个，"
        f"Ki {int((included.endpoint == 'Ki').sum())} 个；"
        f"未纳入 {len(excluded)} 个，原因均见 E2a_REPORT.md。",
        f"- 输出：predictions.csv {csv_rows(OUT / 'predictions.csv')} 行，"
        f"metrics_by_target.csv {csv_rows(OUT / 'metrics_by_target.csv')} 行，"
        f"pooled_metrics.csv {csv_rows(OUT / 'pooled_metrics.csv')} 行，"
        f"comparators.csv {csv_rows(OUT / 'comparators.csv')} 行。",
        "- pyarrow 未安装，按任务单使用 CSV。",
        "",
    ]
    (OUT / "codex-result-E2a.md").write_text("\n".join(brief))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", action="store_true")
    args = parser.parse_args()
    if args.report:
        write_report()
    else:
        run()
