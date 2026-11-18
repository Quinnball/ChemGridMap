# 当前任务：E2c 邻居层在随机抽样靶点上的检验（事先声明）

- **状态：**草稿（Opus，2026-10-29）。在 E2b 验收后、E3-A 之前执行，或在用户的 DeepSeek 预算允许时执行。**本计划在任何人分析抽样靶点数据之前写定并封存**，封存时间以本文件的 SHA-256 记录为准。
- **执行：**DeepSeek（Codex CLI 的 DeepSeek profile），`--sandbox workspace-write`，**不联网**，不做 git 写操作。
- **工作目录与 Python：**与 E2a 相同。
- **依据：**
  - E1-B（`E1B_ACCEPTANCE.md`）中，事先假设「低召回 → 网格邻居的超额活性误差更大」只在三个靶点之一中成立；
  - 精炼对独立活性指标的改进很小。
  - 三个靶点不足以说明这一关系在多大范围内成立。本任务在随机抽样的靶点上检验 E1-B 的两项结论。
- **原则：**
  - 所有被选中的靶点都报告；
  - 某个靶点失败就照实记录并继续；
  - 判据需要改动时，停下并报告。

## 1. 输入（只读）

- `S/data/chembl37_api/`：`selected_targets.csv` 中 `role = sampled` 的 80 个靶点–端点（50 个 IC50，30 个 Ki）及其 `activities/<target>_<type>.csv`。
  - 3 个原靶点的 Ki 不参与，因为它们的布局分析已在 E1-B 中按 IC50 完成。
- 整理：与 E2 的回归检查相同的默认整理。条目为母体化合物，显示值为 pChEMBL 中位数。
- 代码（只导入）：`chemgridmap` 的投影、分配与精炼函数；`S/source/e1b_neighbors.py`，即 E1-B 已验收的邻居、召回与 e 的计算。

## 2. 规则

1. **纳入：**整理后条目数 ≥ 200。不满足的靶点列出原因。
2. **子集：**条目数 > 1,200 时，随机抽取 1,200 个条目，与活性无关，seed 20261102；否则用全部条目。
3. **布局**（每个靶点一个 seed：42）：
   - projection：UMAP，Jaccard 距离、15 个邻居、min_dist 0.1；
   - minimum_displacement：occupancy 0.40，等比缩放；
   - neighbor_refined：软件的默认精炼参数；
   - random：随机布局，seed 42，用于计算局部性指标。
   - **参数一律显式给出**，不依赖默认值。
4. **逐分子指标：**k = 10，保留并列，不含自身，定义与 E1-B 完全相同：recall、e_grid、e_fp，y 为显示用的 pChEMBL 中位数。
5. **每个靶点的统计量：**
   - ρ_refined = Spearman ρ(recall, e_grid − e_fp)，在 neighbor_refined 上计算（主指标）；minimum_displacement 上同样计算；
   - 分子 bootstrap 1,000 次，seed 20261031，得到 ρ_refined 的 95% CI；
   - 各布局的平均 recall 与平均 e_grid；
   - 局部性指标 λ = 1 − e_fp / e_random，其中 e_random 为 random 布局的平均 e_grid。

## 3. 事先声明的检验

每个端点分别检验，单位为靶点：

| 编号 | 内容 | 方法 |
|---|---|---|
| L1 | 召回率低的格子，网格邻居的超额活性误差更大：ρ_refined 的中位数 < 0 | 对各靶点的 ρ_refined 做 Wilcoxon 符号秩检验，单侧，与 0 比较 |
| L2 | 精炼降低网格邻居的活性误差：e_grid(neighbor_refined) < e_grid(minimum_displacement) | 配对 Wilcoxon 检验，单侧 |

- L1 与 L2 在两个端点上共 4 个检验，一起做 Holm 校正。报告原始与校正后的 p 值，以及中位数（或中位差）的 bootstrap CI：按靶点重抽样，2,000 次，seed 20261028。
- **描述性**（不检验）：
  - ρ_refined 的 CI 完全低于 0、包含 0、完全高于 0 的靶点数；
  - 各靶点 ρ_refined 与 λ 的 Spearman 相关。**这是探索性分析，须在结果中标明。**

## 4. 输出（`S/output/E2c/`）

- `layout_targets.csv`：纳入情况与子集大小；
- `layout_by_target.csv`：每个靶点 × 布局的平均 recall、e_grid、e_fp，以及 ρ 及其 CI；
- `layout_tests.csv`；
- `RESULTS_E2c.md`：英文，只写事实；
- `manifest.json`：输入与源码的哈希、输出的哈希、版本、种子与运行时间。

## 5. 验证与汇报

- **测试**（`S/source/test_e2c.py`）：
  - 在三个原靶点中的一个上，用本任务的流程与 E1-B 的参数重算 seed 42 的 ρ，与 `S/output/E1B/neighbor_sar_summary.csv` 一致，误差 ≤ 1e-12；
  - 子集抽取与活性无关：抽取时不读活性列。
- **冻结核对：**与 E2a 相同；`S/output/E2/` 不得改动。
- **报告：**`REV/codex_tasks/S/E2c_REPORT.md`，只写事实；简要版写入 `docs/codex-result.md`。
