"""E1-A 的英文事实表与可复现清单；不改动冻结材料。"""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
REV = ROOT / "paper"
S = REV / "strengthening"
OUT, QA = S / "output/E1A", S / "qa/E1A"
NOTE = "post hoc sensitivity analysis motivated by the flagged hERG pair"


def markdown_table(frame):
    """不依赖额外排版库，NA 保持可见。"""
    def cell(value):
        if pd.isna(value):
            return "NA"
        if isinstance(value, (float, np.floating)):
            return f"{value:.6g}"
        return str(value).replace("|", "\\|").replace("\n", " ")
    rows = ["| " + " | ".join(map(str, frame.columns)) + " |",
            "| " + " | ".join(["---"] * len(frame.columns)) + " |"]
    rows.extend("| " + " | ".join(cell(v) for v in row) + " |"
                for row in frame.itertuples(index=False, name=None))
    return "\n".join(rows)


def metric_ci_table(metrics, intervals, target=None):
    frame = metrics if target is None else metrics[metrics.target == target]
    keys = ["population", "target", "stratum", "method"]
    result = frame[keys + ["n_predictions", "n_entries"]].copy()
    for metric in ("auc", "ece", "brier"):
        ci = intervals[intervals.metric == metric].copy()
        ci[metric + " [95% CI]"] = ci.apply(
            lambda r: f"{r.estimate:.6f} [{r.ci_low:.6f}, {r.ci_high:.6f}]", axis=1)
        result = result.merge(ci[keys + [metric + " [95% CI]"]], on=keys, validate="one_to_one")
    return result


def summary_frame(frame):
    if "row_type" in frame:
        frame = frame[frame.row_type.isin(["summary", "union_summary"])]
    return frame.drop(columns=[c for c in ("sensitivity_note", "row_type") if c in frame]).dropna(axis=1, how="all")


def write_results(metrics, intervals, paired, descriptive, flags):
    regression = json.loads((QA / "regression.json").read_text())
    text = ["# E1-A factual results", "", "All analyses follow docs/task.md, including rework R1.",
            "Sensitivity population: **" + NOTE + "**.",
            "", "## 3.0 Regression", "",
            f"B0/G: {regression['aligned_rows']} aligned predictions; maximum absolute probability error "
            f"{regression['max_probability_error']:.17g}. Observed and displayed classes match exactly.",
            markdown_table(pd.DataFrame([regression["G"]])),
            "", "## 3.1 Models", "",
            "B0, G and X retain the frozen merged predictions. B1 uses target-specific B0 agreement from entries "
            "outside the test fold. G-record fits dispersion over activity_id units and removes every record "
            "from every document in the held-out merged unit. G-fixed uses tau = 0.64. All six models share "
            "folds, held-out document units and displayed classes; each has 820 predictions from 289 target-entry clusters.",
            "", "## 3.2 Evaluation populations", "", NOTE + ".",
            markdown_table(flags.drop(columns=[c for c in ["sensitivity_note"] if c in flags])),
            "predictions.csv uses population rows: flag-clean repeats the retained subset of all. "
            "Its stratum is m=1 or m>=2, where m always counts remaining merged document units. "
            "Metrics additionally include the unsplit stratum all.",
            "", "## 3.3 Metrics", "",
            "Accuracy and ECE use the argmax class and maximum probability, respectively. "
            "AUC and overall calibration bias use the displayed-class probability and displayed/observed agreement. "
            "Brier is the unscaled sum over three classes; log loss clips probabilities at 1e-6. "
            "ECE uses 10 equal-width bins. High confidence means maximum class probability >= 0.8.",
            markdown_table(metric_ci_table(metrics, intervals, "all")),
            "", "## 3.4 Uncertainty", "",
            "Percentile 95% intervals use 2,000 entry-cluster bootstrap samples, seed 20261028. "
            "Each target is resampled separately in pooled estimates. Populations, targets and strata are "
            "evaluated separately. All models and each G-minus-alternative difference use the same draws "
            "within an evaluation cell; models are not refitted during this prediction bootstrap.",
            f"metrics.csv contains {len(metrics)} cells; bootstrap_ci.csv contains {len(intervals)} intervals; "
            f"paired_differences.csv contains {len(paired)} paired differences. "
            "Valid bootstrap counts are reported; single-class AUC replicates are NA. Negative-control results are stored in QA.",
            markdown_table(paired[(paired.target == "all") & (paired.stratum == "all")].drop(
                columns=["sensitivity_note", "n_bootstrap", "seed"])),
            "", "## 3.5 Tau influence", "",
            "Full-data and flag-clean fits use the frozen default evidence-unit merge. Leave-one-entry changes "
            "are removal-fit tau minus full-fit tau. Capped weights are min(n_units - 1, 4). "
            "Tau bootstrap resamples the full retained entry population, including single-unit entries.",
            markdown_table(summary_frame(descriptive["tau_influence"])),
            "The leave-one-entry rows show the five largest absolute changes per target; all removals are retained in QA/tau_leave_one_entry_all.csv.",
            "", "## 3.6 Pseudo-consensus", "",
            "Entries with at least five units contribute one observation per unit. "
            "Training tau comes from the entry's held-out fold. p_cons uses SD=tau; p_pred uses "
            "SD=tau*sqrt(1+1/(n_units-1)); the outcome is the class of the other units' mean. "
            "These are descriptive observations from small samples.",
            markdown_table(summary_frame(descriptive["pseudo_consensus"])),
            "", "## 3.7 Counting-based status changes", "",
            "Only the frozen G/merged mapped entries are evaluated. Document support is frozen class_confidence. "
            "Record support integrates the displayed class using the record mean and tau_rec/sqrt(n_records). "
            "The supported threshold is 0.8. Largest-document fraction counts original source documents, not merged units.",
            markdown_table(summary_frame(descriptive["count_flip_summary"])),
            "", "## 3.8 Flagged document pairs", "",
            "Offsets are document B minus document A, based on per-document medians for shared entries. "
            "Metadata is taken only from the supplied local JSON; unavailable metadata is NA. "
            "Prediction contributions count one shared evaluation design per entry, not six copies across models.",
            markdown_table(summary_frame(descriptive["flagged_pairs_profile"])),
            "", "## 3.9 Single-unit support thresholds", "",
            "z = Phi^-1(0.8). Active support >=0.8 requires unit center >=7+z*tau; inactive support >=0.8 "
            "requires unit center <=6-z*tau. The unit center is the document-balanced pChEMBL; "
            "the displayed class remains the frozen record-median class. Medium-class support is also checked.",
            markdown_table(summary_frame(descriptive["single_unit_thresholds"])), ""]
    (OUT / "RESULTS_E1A.md").write_text("\n\n".join(text))


def write_manifest(negative=None):
    if negative is None:
        negative = pd.read_csv(QA / "negative_control.csv")
    targets = ("chembl205", "chembl204", "chembl240")
    inputs = [ROOT / f"paper/output/revision_v6/{t}/{t}_retained_records.csv" for t in targets]
    inputs += [ROOT / "paper/output/reliability_v1/activity/lodo_predictions.csv",
               ROOT / "paper/output/reliability_v1/activity/mapped_entry_confidence.csv",
               ROOT / "paper/reliability_activity.py", ROOT / "src/chemgridmap/evidence.py",
               ROOT / "src/chemgridmap/activity.py", ROOT / "paper/reliability_analysis.py",
               S / "data/flagged_pairs_metadata.json", S / "SELECTION_RULES.md"]
    def hashes(paths):
        return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
    parameters = {
        "targets": targets, "models": ["B0", "B1", "G", "G-record", "G-fixed", "X"],
        "evaluation_set": "merged", "folds": 10, "fold_seed": 20261024,
        "bootstrap_seed": 20261028, "bootstrap_repetitions": 2000, "percentiles": [2.5, 97.5],
        "negative_control_seed": 20261030, "negative_control_permutation": "binary agree within target; p fixed",
        "negative_control_is_gate": False, "class_boundaries": [6.0, 7.0], "confidence_threshold": 0.8,
        "log_loss_clip": [1e-6, 1.0], "ece_bins": 10, "ece_probability": "maximum class probability",
        "auc_probability": "p_displayed", "auc_ties": "average ranks", "auc_single_class": "NA",
        "output_csv_readback": "pandas.read_csv(..., float_precision='round_trip') to preserve exact AUC ties",
        "merge_identical": True, "merge_decimals": 2, "G_fixed_tau": 0.64,
        "B1": "B0 agreement over all predictions of same-target entries outside test fold",
        "G_record_units": "activity_id; merge_identical=False; fit_between_unit_sd defaults",
        "G_record_prediction": "remaining record mean; tau_rec*sqrt(1+1/n_records); exclude all held-unit documents",
        "fit_between_unit_sd_min_entries": 5, "literature_fallback_tau": 0.64,
        "cap_degrees_of_freedom": 4, "pseudo_consensus_min_units": 5,
        "flag_min_shared": 3, "flag_offset_threshold": 1.0, "flag_sign_agreement": 0.8,
        "flag_clean_selection": "remove any entry with any constituent document flagged in full target records",
        "sensitivity_note": NOTE, "prediction_strata": ["m=1", "m>=2"], "metric_strata": ["all", "m=1", "m>=2"],
        "m_definition": "remaining merged document units for every model",
        "bootstrap_iteration": "population all,flag-clean; target all,chembl205,chembl204,chembl240; stratum all,m=1,m>=2",
        "bootstrap_target_order_within_pool": "sorted target ids", "bootstrap_entry_order": "sorted entry ids",
        "prediction_bootstrap_refit": False, "paired_bootstrap": "same entry multiplicities for every model",
        "tau_bootstrap_population": "all retained entries, including single-unit entries; duplicate draws count separately",
        "tau_bootstrap_iteration": "independent stream reset to seed 20261028; sequential targets chembl205,chembl204,chembl240; sorted entries; batches of 100",
        "prediction_population_encoding": "long table: all plus repeated flag-clean subset",
        "frozen_probabilities_comparison_tolerance": 1e-12,
    }
    manifest = {"completed_task": True, "analysis_completed": True,
        "acceptance_complete": bool(negative.contains_half.all()),
        "qa_negative_control_passed": bool(negative.contains_half.all()),
        "unresolved_acceptance_items": [] if negative.contains_half.all() else
            ["R1 negative-control CI excludes 0.5 for: " + ", ".join(negative.loc[~negative.contains_half, "target"])],
        "revision": "R1", "parameters": parameters,
        "inputs_sha256": hashes(inputs),
        "source_sha256": hashes(list((S / "source").glob("e1a_*.py")) + [S / "source/test_e1a.py"]),
        "outputs_sha256": hashes([p for p in OUT.iterdir() if p.is_file() and p.name != "manifest.json"]),
        "qa_sha256": hashes([QA / "regression.json", QA / "negative_control.csv", QA / "negative_control_R0_misspecified.csv"]),
        "manifest_self_hash": "omitted to avoid a self-referential hash",
        "software_versions": {"python": platform.python_version(), **{name: importlib.metadata.version(name)
            for name in ("numpy", "pandas", "scipy", "pytest")}},
        "excluded_without_access": "strengthening/data/chembl37_api"}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


def append_report():
    """保留 R0 原文，在末尾维护一节 R1；完整运行证据取自 QA 日志。"""
    load = lambda name: pd.read_csv(OUT / (name + ".csv"))
    metrics, intervals = load("metrics"), load("bootstrap_ci")
    tau, pseudo, pairs = load("tau_influence"), load("pseudo_consensus"), load("flagged_pairs_profile")
    negative = pd.read_csv(QA / "negative_control.csv")
    ledger = json.loads((QA / "commands_R1.json").read_text())
    scope_path = QA / "workspace_R1_changes.json"
    scope = json.loads(scope_path.read_text()) if scope_path.exists() else {}
    # 数值表按用途选列，CSV 保留全部元数据及逐条记录。
    tau_summary = tau[tau.row_type == "summary"][["target", "n_entries", "n_multi_unit_entries",
        "tau_full", "tau_flag_clean", "tau_cap", "tau_ci_low", "tau_ci_high"]]
    tau_top = tau[tau.row_type == "leave_one_entry"][["target", "entry", "n_units",
        "degrees_of_freedom_fraction", "tau_without_entry", "delta_tau", "absolute_delta_rank"]]
    pseudo_summary = pseudo[pseudo.row_type == "summary"][["target", "n_units", "n_entries",
        "mean_p_cons", "mean_p_pred", "n_agree", "agreement_rate"]]
    pair_summary = pairs[pairs.row_type.isin(["summary", "union_summary"])][["target", "document_a",
        "document_b", "scope", "n_entries", "n_prediction_entries", "n_shared_entries",
        "median_offset_b_minus_a", "spearman_rho", "n_predictions", "n_agree"]]
    sections = ["## 9. 返工 R1", "",
        "**本节为本轮最终状态：3.0–3.9 分析与第 4 节全部输出已完成；回归、12 项单元测试、"
        "独立数值复核、冻结 guard 与 CGM 状态检查通过。工作区范围比对发现并发范围外新增文件，未通过。"
        "R1 负控已按新定义执行，但 all 与 hERG 的 CI 不含 0.5，"
        "因此负控验收未通过。依 R1 第 2 条，此结果没有阻止其余分析。前述 R0 停止状态仅保留为历史记录。**",
        "### R1.1 修改与完成范围", "",
        "旧 qa/E1A/negative_control.csv 原字节重命名为 negative_control_R0_misspecified.csv；"
        "新负控仅在靶点内置换 agree=(observed==displayed)，保持原始类别与全部概率不变，seed=20261030。"
        "沿用回归与加权统计函数，新增六模型构建、全部指标与成对整群 bootstrap、3.5–3.9 描述分析、"
        "英文事实表、清单及独立验证辅助模块。所有修改均限于任务允许位置。",
        "全程离线，不安装依赖，不做 git 写操作，不操作 Prism/PowerPoint。"
        "使用已有 chemgridmap/.venv/bin/python；PYTHONDONTWRITEBYTECODE、MPLCONFIGDIR 与 pytest 的 "
        "no:cacheprovider 只控制缓存位置，不改变分析规则。",
        "### R1.2 口径与输出完整性", "",
        "六模型为 B0、B1、G、G-record、G-fixed、X。B0/G/X 原冻结预测逐值保留；B1 仅用同靶点折外条目全部 B0 预测；"
        "G-record 的拟合使用 activity_id 单元，预测时删除留出单元所有构成文献的全部记录。"
        "所有模型共用固定条目 folds、文献留出单元与剩余记录中位数显示类别。m 始终为剩余合并文献单元数。",
        "flag-clean 标注：**" + NOTE + "**。完整 hERG 数据标出 5 对、9 篇文献，涉及 265 个条目；"
        "其中 104 个有留出预测的条目被整体剔除，减少 363 次预测。CA II 与凝血酶均为 0 对。"
        "因此 all 每模型 820 次/289 条目；flag-clean 每模型 457 次/185 条目。"
        "该定义包含全部五对标记文献，不能与 E0 仅剔除首对文献的结果混用。",
        "predictions.csv 是总体长表：all 的 4920 行与 flag-clean 子集的 2742 行，共 7662 行；"
        "这不代表 7662 次独立预测。指标含 6 模型 × 4 靶点层级（含合计）× 2 总体 × 3 分层 = 144 格；"
        "bootstrap_ci.csv 含 1728 行，paired_differences.csv 含 480 行；第 4 节 12 个输出均齐全。"
        "准确率/ECE 采用 argmax/最大概率，AUC 与整体校准偏差采用 p_displayed/显示一致指标，二者没有混用。",
        "按条目整群重抽样 2000 次，seed=20261028；合计在各靶点内分别抽原条目数；"
        "每个总体×靶点×分层独立抽样，同格所有模型及 G 减替代模型差异共享重数。"
        "预测 bootstrap 不重新拟合模型。τ bootstrap 在完整条目总体（含单单元）抽样；重复抽中条目独立计入自由度。",
        "CSV 复核应使用 pandas.read_csv(..., float_precision='round_trip')。默认解析会改变少量浮点末位，"
        "从而改变 AUC 的严格同分关系；首次独立验证因此退出 1，修正验证器的输出读取方式后全部点估计及区间复核通过。"
        "没有舍入分数、改平局规则或改分析结果。一次并行复核的末尾也因报告模块仍在编辑、manifest 尚未刷新而退出 1；"
        "最终源文件稳定后刷新清单，再完整复核。",
        "### R1.3 指标表（点估计及百分位 95% CI）", "",
        markdown_table(metric_ci_table(metrics, intervals)),
        "### R1.4 τ 的影响", "", markdown_table(tau_summary),
        "各靶点 |Δτ| 最大的 5 个条目如下。Δτ=剔除后τ−完整τ；自由度占比=(n_units−1)/总自由度。"
        "全部 289 个多单元条目的留一结果位于 qa/E1A/tau_leave_one_entry_all.csv。",
        markdown_table(tau_top),
        "### R1.5 计数翻转与伪共识", "",
        markdown_table(summary_frame(load("count_flip_summary"))),
        "360 次状态翻转均为 uncertain → supported；完整明细含记录数、单元数、最大原始文献记录占比、"
        "两种支持度、显示类别与中位数。3289 个映射条目的按文献支持度另经完整 τ 重算，与冻结注释误差≤1e-12。",
        markdown_table(pseudo_summary),
        "伪共识只作描述：使用所在条目折外训练 τ，预测其余单元均值的类别，未据此作显著性或因果推断。",
        "### R1.6 文献对档案与单单元阈值", "", markdown_table(pair_summary),
        "CHEMBL5725903 / CHEMBL5727624 的 30 个共享条目贡献 60 次预测，0 次一致；"
        "两文献涉及的全部 67 个条目中有 32 个预测条目，贡献 64 次预测，0 次一致。"
        "两种口径分别报告。完整档案保存共享条目两侧原始记录值、文献中位数、类别、差值、Spearman、类型、年份、专利或期刊。"
        "CHEMBL1136167 / CHEMBL3085767 的一侧值恒定，Spearman 记 NA；缺失来源元数据保留 NA。",
        markdown_table(summary_frame(load("single_unit_thresholds"))),
        "单单元阈值按 document_balanced_pchembl（单个合并单元的中心）核对，显示类别仍为冻结记录中位数类别；"
        "全部 3208 个单单元映射条目状态差异为 0。三个 τ 下 medium 类理论最高支持度均低于 0.8。",
        "### R1.7 负控结果与未通过项", "",
        markdown_table(negative[["target", "auc", "ci_low", "ci_high", "n_bootstrap_valid", "contains_half"]]),
        "按要求只置换一次，未更换种子、筛选置乱结果或另行调整规则。独立将全部 8000 次整群抽样显式展开，"
        "用 scipy.stats.rankdata(method='average') 重算 AUC，和 QA 输出一致。",
        "靶点内二分类置换使每个单靶点 AUC 的零假设期望为 0.5，但不同靶点的一致率与概率分布不同，"
        "合计 AUC 的精确期望为 0.5130782469382973。计算式为："
        "E(AUC)=[Σ_t n1_t·n0_t/2 + Σ_(t≠u) q_t·(1−q_u)·Σ_(i∈t,j∈u)H(s_i,s_j)]/(N1·N0)，"
        "H(a,b)=1[a>b]+0.5·1[a=b]，q_t=n1_t/n_t，s=p_displayed。"
        "这一事实不能解释或消除本次 hERG 的未通过；本轮固定置换实测结果仍按原验收判据报告为失败。",
        "**未完成分析：无。未通过的验收项：负控 all 与 hERG 的 CI 不含 0.5；须由协调者决定如何处理该 QA 判据，"
        "本轮没有自行替换规则。工作区范围检查因外部并发产物未通过，S/data 全目录未变也不能成立。"
        "E1-A 完成执行不等于全部验收通过，也不代表 E2 或整篇论文补强已完成。**",
        "### R1.8 验证命令、退出码及原始输出", ""]
    for item in ledger:
        sections.extend([f"#### {item['name']}（退出码 {item['exit_code']}）", "",
            "```sh\n" + item["command"] + "\n```", "",
            "```text\n" + (QA / item["log"]).read_text().rstrip() + "\n```", ""])
    sections += ["### R1.9 冻结与写入范围", "",
        "初始与最终 verify_guard.sh 均六处 OK。CGM 状态前后逐字一致、45 行；"
        "指定冻结输入由 manifest SHA-256 核对。工作区比对使用文件大小与 mtime_ns，作为 guard/输入哈希的补充。"
        "本轮没有直接写入 S/data；指定 flagged_pairs_metadata.json 哈希未变，"
        "但快照发现并发新增 S/data/prior_art/doi_verified.tsv 和 notes.md，因此不声称 S/data 全目录未变。"
        "chembl37_api 在进入目录前排除，未读取或写入，"
        "因此不对 Opus 下载子目录的外部变化作未变声明。",
        "范围外还出现 E0_prior_art_summary.md、E2_task_draft_sealed.md 及其 sha256、"
        "download_chembl37_sample_parallel.py、下载器字节码、新版 v3.1 交接文件，以及 E1A_R1_run.log 增长。"
        "新 E0 摘要标明来源为独立调研代理和 Opus 复核；lsof 确认运行日志由现有 zsh/node/codex 的 stdout/stderr 打开。"
        "这些文件不在本执行器及子代理的直接写入清单内。ps 进程查询受沙箱限制返回 1，"
        "因此不对每个外部文件的具体写入进程作已核实声明；没有回滚、删除或修改这些并发产物。",
        "范围核对结果：", "```json\n" + json.dumps(scope, indent=2, ensure_ascii=False) + "\n```",
        "新增/修改/重命名文件详见上述清单及 qa/E1A/workspace_R1_changes.json；"
        "旧 negative_control.csv 的删除与新 R0 文件的出现是规定的重命名。"
        "报告追加保留 R0 全部原文；docs/codex-result.md 只保留本轮简要结论。", ""]
    report_path = OUT / "E1A_REPORT.md"
    prior = report_path.read_text().split("\n## 9. 返工 R1", 1)[0]
    report_path.write_text(prior + ("\n" if prior.endswith("\n") else "\n\n") + "\n\n".join(sections))
    (ROOT / "docs/codex-result.md").write_text(
        "# E1-A R1 简要结论\n\n"
        "已完成 R1 负控修正、3.0–3.9 全部分析及第 4 节 12 个输出；旧负控原样保留为 negative_control_R0_misspecified.csv。\n\n"
        "完整入口、回归、pytest（12 passed）、最终独立复核和冻结 guard 退出码均为 0；guard 六处 OK，CGM 状态保持 45 行。\n\n"
        "负控验收命令退出码 1：all 和 hERG 的 95% CI 不含 0.5。范围检查退出码 1：发现外部并发新增 E2、prior_art、下载脚本及交接文件；未直接改写这些文件。按 R1 完成其余分析，未改 seed 或规则。未完成分析无；上述验收项待协调者处理。\n\n"
        "详细报告：paper/strengthening/output/E1A/E1A_REPORT.md。\n")
