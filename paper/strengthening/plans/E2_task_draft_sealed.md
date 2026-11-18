# 当前任务：E2 多靶点评价（ChEMBL 37 随机样本，IC50 与 Ki）

- **状态：**草稿（Opus，2026-10-28）。数据下载完成后启用。**本任务单中的假设、指标与判据在任何人查看结果之前写定**；写定时间以本文件的 SHA-256 记录为准。
- **执行、工作目录、Python：**与 E1-A 相同。**不联网。**
- **依据：**`REV/codex_tasks/S/E0_review_summary.md` §3 第 1、3、4、5 项；`E0_prior_art_summary.md` §5。
- **原则：**
  - 所有被选中的靶点都要报告，不得事后剔除；
  - 某个靶点或某一步失败，就照实记录并继续，不要停止整个任务；
  - 如果判据需要改动，停下并报告。

## 1. 输入（只读）

- `S/data/chembl37_api/`：
  - `selected_targets.csv`（50 个 IC50 靶点与 30 个 Ki 靶点，随机抽样；另含 3 个原靶点的 Ki）；
  - `activities/<target>_<type>.csv`；
  - `download_manifest.json`；
  - `download_original_ic50_manifest.json`（3 个原靶点的 IC50，用于回归检查与 hERG 的年份）；
  - `SELECTION_RULES.md`（含 Amendment 1）。
- 冻结的保留记录：`CGM/paper/output/revision_v6/<t>/<t>_retained_records.csv`。整理参数取自同一目录下的 `<t>_curation_report.json`。
- 代码（只导入）：
  - `chemgridmap.chembl.curate_chembl_activity_data`；
  - `chemgridmap.evidence`；
  - `CGM/paper/reliability_activity.py`；
  - E1-A 的统计模块 `S/source/e1a_statistics.py`（已验收的版本，复用其中的 AUC、ECE、bootstrap 与各模型实现）。

## 2. 修改范围

- **允许：**
  - 新建 `S/source/e2_*.py` 与 `S/source/test_e2.py`；
  - 写入 `S/output/E2/` 与 `S/qa/E2/`；
  - 报告写入 `REV/codex_tasks/S/E2_REPORT.md`。
- **冻结：**其余全部，包括 `S/data/`、`S/output/E1A/`、`CGM/`（45 行）。

## 3. 要做的事

### 3.0 回归检查

用 API 格式下载的 3 个原靶点 IC50，按冻结的整理参数整理后，保留记录的 `activity_id` 集合与 pChEMBL 值必须与冻结的保留记录完全一致。

- 如果不一致：列出差异记录及其原因（例如 API 字段的表示方式），在报告中说明，然后继续。
- 只有当差异超过 0.5% 的记录、或者原因不明时，才停止并报告。

### 3.1 纳入规则（与 `SELECTION_RULES.md` 相同）

- 一个靶点–端点在默认整理与默认合并之后，须有 ≥ 30 个单元数 ≥ 2 的条目，才纳入评价。
- 不满足的列出原因。
- IC50 与 Ki 分开分析，**不混合**。

### 3.2 模型

与 E1-A 相同：B0、B1、G、G-record、G-fixed、G-flag（各靶点按条目分折，seed 20261024）。另新增：

- **G-global（留一靶点的共同 σ_doc）**：对靶点 t，σ_global 是同一端点下**其他所有**纳入靶点的多单元条目合并计算的 pooled SD。合并方法与 `fit_between_unit_sd` 相同，自由度加权。其余与 G 相同。
- **G-flag 与 flag-clean 总体**：定义与 E1-A 相同，标记对用该靶点的完整数据确定。这一规则在 E2 中是事先声明的。

### 3.3 事先声明的假设与检验

每个端点分别检验。单位是靶点，n = 纳入的靶点数。

| 编号 | 内容 | 方法 |
|---|---|---|
| H1（问题的规模） | 中位数着色与留出文献不一致的比例 | 描述性，不检验。报告各靶点 1 − acc(B0) 的分布（中位数、IQR、范围），以及合计值（整群 bootstrap CI） |
| H2（校准） | G 的校准程度 | 描述性。报告各靶点 ECE 与整体校准偏差的分布；合计校准曲线（5 个等宽区间） |
| H3（区分能力） | G 的 AUC > 0.5 | 对各靶点的 AUC − 0.5 做 Wilcoxon 符号秩检验（单侧） |
| H4（按文献计数） | G 的 log loss 低于 G-record | 对各靶点的平均 log loss 差做配对 Wilcoxon 检验（单侧） |
| H5（靶点特异的 σ_doc） | G 的 log loss 低于 G-global | 同上 |
| H6（估计的 σ_doc） | G 的 log loss 低于 G-fixed | 同上 |

- H3–H6 在每个端点内用 Holm 法校正。报告原始 p 值与校正后的 p 值，以及中位差及其 bootstrap CI（按靶点重抽样，2,000 次，seed 20261028）。
- Brier 与 ECE 的同类比较作为次要结果，只描述，不检验。
- 每个假设都在 all 与 flag-clean 两个总体上报告。**主分析是 all 总体。**

### 3.4 只剩一个单元时的评价（地图上大多数格子的情形）

H2、H3 在 m = 1 与 m ≥ 2 的预测上分别报告。

### 3.5 时间切分（前瞻性检验，事先声明）

- **主截止年份为 2015，敏感性分析用 2010。**
- 单元的年份 = 其构成文献的最晚 `document_year`。年份缺失的单元不参与本项分析，并计数报告。
- 对每个靶点：
  - 「过去」单元的年份 ≤ 截止年，「将来」单元的年份 > 截止年；
  - σ_past 只用过去单元拟合（过去单元数 ≥ 2 的条目）；少于 5 个这样的条目时，该靶点不参与，并记录原因；
  - 符合条件的条目须同时有 ≥ 1 个过去单元与 ≥ 1 个将来单元；
  - 显示类别 = 过去记录中位数的类别；
  - 预测每个将来单元：中心 = 过去单元的均值，SD = σ_past·√(1 + 1/m_past)；
  - 一致 = 将来单元的类别等于显示类别。
- **在端点内跨靶点合并后**报告：
  - B0 的不一致率；
  - G 的 ECE、整体校准偏差与 AUC，均带整群 bootstrap CI（按靶点 × 条目抽样）；
  - 按 m_past = 1 与 ≥ 2 分层；
  - 「过去显示为 supported（consensus ≥ 0.8）」与「过去显示为 uncertain」两组条目在将来的一致率。这是地图颜色最直接的前瞻含义。
- 将来单元 ≥ 30 个的靶点，另外单独报告。
- 3 个原靶点的 IC50 使用 API 下载的数据（其中有年份）。

### 3.6 伪共识检验

- 单元数 ≥ 5 的条目，在端点内跨靶点合并。定义与 E1-A 的 3.6 相同。
- 报告 p_cons 与 p_pred 的校准曲线、ECE，以及 n。

### 3.7 σ_doc 的分布与文献对照

- 列出各靶点–端点的 σ_doc 及其 bootstrap CI（按条目抽样，1,000 次）、自由度、多单元条目数。
- 按端点给出分布。
- 写出文献参照值（单次测量尺度）：
  - Ki ≈ 0.47（Kalliokoski 2013，Table 3，ΔpKi 的 σ 0.66–0.68 ÷ √2）；
  - Ki 0.54（Kramer 2012）；
  - IC50 ≈ 0.57–0.61（Kalliokoski 2013）。
  - 只描述，不检验。
- 描述性附加：同时具备 B 型与 F 型测定、且两类各有 ≥ 30 个多单元条目的靶点，分别给出两类的 σ_doc。
- CA II 的 Ki：报告最大的文献簇（同一通讯作者或同一系列）对 σ_doc 的影响。做法是去掉贡献自由度最多的 5 篇文献后重新计算。依据见 `E0_prior_art_summary.md` §5。

### 3.8 被标记文献对的普查

对每个靶点，列出：

- 系统性偏移文献对的数量；
- 被标记文献所涉及条目与预测的占比；
- 偏移最大的 10 对（文献 ID、共享条目数、中位偏移、Spearman 相关）。

这些是潜在的数据错误，作为案例候选，只描述。

### 3.9 条目层面的统计（对应「地图显示什么」）

对每个靶点（全部条目，不限于映射子集），报告：

- 单单元条目的比例；
- 记录数大于单元数的条目比例；
- 支持度 ≥ 0.8 的比例；
- 按记录计数与按文献计数时，supported 与 uncertain 状态发生翻转的条目数，以及翻转方向。

## 4. 输出（`S/output/E2/`）

- `target_inclusion.csv`
- `predictions.parquet`（若无 parquet 支持则用 CSV）
- `metrics_by_target.csv`
- `hypothesis_tests.csv`
- `pooled_metrics.csv`
- `time_split.csv`
- `time_split_by_target.csv`
- `pseudo_consensus.csv`
- `sigma_doc_by_target.csv`
- `flagged_pairs_census.csv`
- `entry_level_by_target.csv`
- `forest_plot_data.csv`（每个靶点的 AUC、ECE 及其 CI，以及 G 与各替代模型的 log loss 差，供 Prism 使用）
- `RESULTS_E2.md`（英文，只写事实）
- `manifest.json`

## 5. 验证

1. 3.0 的回归检查结果写入报告。
2. pytest：`S/source/test_e2.py`，至少覆盖：
   - G-global 不使用靶点 t 自身的数据；
   - 时间切分不让将来的数据进入 σ_past 或显示类别；
   - Holm 校正；
   - Wilcoxon 的方向（单侧）；
   - 年份缺失时的处理。
3. 负控：把各靶点的二分类一致指标置换后，H3 的检验应不显著。
4. 冻结核对与 E1-A 相同。

## 6. 汇报

- `REV/codex_tasks/S/E2_REPORT.md`，内容包括：
  - 纳入情况；
  - 回归检查；
  - H1–H6 的结果表；
  - 时间切分；
  - 伪共识检验；
  - σ_doc；
  - 标记对；
  - 问题与未完成事项。
- `docs/codex-result.md` 只写简要结论。

---

## Amendment 1（2026-10-28 20:31；此时任何抽样靶点的活性数据都尚未分析）

依据：E1-A 的验收（`E1A_ACCEPTANCE.md` 第 2、4 节）。各项假设与判据都没有放宽，只修正 QA 的定义并补充明确性：

1. **负控**（替换第 5 节第 3 条）：在每个靶点内，把二分类的一致指标置换 200 次（seed 20261030），每次重新做 H3 的 Wilcoxon 检验。
   - **QA 通过的条件**：在 α = 0.05 下，拒绝率 ≤ 0.10。
   - 另外，对每个靶点给出 G 的 AUC 的置换 p 值（1,000 次置换），作为描述性结果。
   - 合计 AUC 的零分布不是 0.5，因为不同靶点的一致率不同，**不得**与 0.5 比较。
2. **G-flag 即模型 X**：实现为 `reliability_activity.excluded_prediction`。
   - 标记文献用当前折训练条目的记录、由 `flagged_documents` 确定；
   - 拟合 σ_doc 与预测时都排除带标记文献的单元；
   - 没有剩余单元时回退到 G。
   - 与 E1-A 相同，只是推广到多靶点。
3. **H7**（新增，与 H3–H6 一起做 Holm 校正）：X 的 log loss 低于 G，用配对 Wilcoxon 单侧检验，单位为靶点；只纳入至少有一对标记文献的靶点，并报告这些靶点的数量。
   - **决策规则（事先声明）**：如果 H7 在两个端点中都显著，0.5.0 的默认模型改为 X；否则保留 G，把 X 作为选项。
4. **flag-clean 总体中的公平比较**：在 flag-clean 总体中，所有模型的 σ_doc 都改用剔除标记条目后的训练数据重新拟合，以避免标记对的扭曲进入训练。在 all 总体中保持原样。

## Amendment 2（2026-10-28 21:53；数据已下载完毕，只核对了文件数、行数与哈希，任何抽样靶点的活性数据都尚未分析）

依据：J. Cheminform. 的 Software 类文章通常要求与现有软件直接比较（交接文档 v3.0 §8）。现有工具与数据集给出的可靠性标注是冲突标记或离散度统计：Isigkeit et al. 2022、CAPRICHO 按极差分级，Papyrus 给出 SD，ChemGridMap 自己的冲突三角也属于这一类。本修订只**增加**比较对象与描述性输出，H1–H7 及其判据、Holm 校正的范围都不变。

1. **代表现有做法的比较对象**（与 G 用同一批预测、同样的折与留出单元）：
   - **B-range（冲突标记）**：只用测试条目**剩余**的记录计算，与 B0 相同，被留出文献的记录不参与。剩余记录的极差 ≤ 1 个对数单位时得分为 1，否则为 0。剩余记录只有一条时，极差为 0，得分为 1。
   - **B-sd（离散度）**：得分 = −SD（剩余记录）。只在剩余记录 ≥ 2 条的预测上计算，并报告这一子集的 n；同一子集上也报告 G 的 AUC，以便比较。
   - 指标：以显示类别是否一致为结局的 AUC。并列按标准的秩方法处理。
   - **H8**（次要假设）：G 的 AUC 高于 B-range 的 AUC。
     - 方法：各靶点的 AUC 差，配对 Wilcoxon 单侧检验，单位为靶点，每个端点分别做。
     - H8 自成一族：在 {H8-IC50, H8-Ki} 内做 Holm 校正，**不并入** H3–H7 的校正。
     - 在 all 与 flag-clean 两个总体上报告，主分析为 all。
2. **按显示状态分层的一致率**（描述性，对应「饱和色意味着什么」）：
   - 对每个留出预测，用测试条目剩余的单元计算显示类别的**共识支持度**（SD = σ_doc/√m，与地图相同），并按 ≥ 0.8（显示为 supported）与 < 0.8（uncertain）分组。
   - 在端点内跨靶点合并，报告两组被留出文献的一致率及其整群 bootstrap CI（按靶点 × 条目抽样，2,000 次，seed 20261028），并按 m = 1 与 m ≥ 2 分层。
   - 这与 3.5 的时间切分中「过去 supported / uncertain」的比较相同，只是用于留一单元的评价。
3. **新增输出**（`S/output/E2/`）：
   - `comparators.csv`：各靶点、各总体中 B-range、B-sd 与 G 的 AUC，以及 n；
   - `agreement_by_display_status.csv`；
   - `forest_plot_data.csv` 增加 B-range 的 AUC 一列；
   - `hypothesis_tests.csv` 增加 H8 的两行，并注明它们单独校正。
4. **测试**：`test_e2.py` 须另外覆盖：
   - B-range 与 B-sd 不使用被留出文献的记录；
   - 剩余记录只有一条时 B-range 为 1，B-sd 不计算。
