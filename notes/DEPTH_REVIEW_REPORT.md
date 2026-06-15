# depth 命令多 agent 协同审阅报告

**对象**: `src/commands/depth.rs` (4007 行) + `notes/DEPTH_COMMAND.md`
**方法**: 7 轮迭代，共 14 位独立审阅员 / 验证员
- **R1** (4 并行)：设计 / 架构 / 正确性 / 性能
- **R2** (2 并行)：设计+架构验证 / 正确性验证
- **R3** (3 并行)：性能条目验证 / 测试覆盖审阅 / 深度正确性重审
- **R4** (3 并行)：R3 新发现验证 / 跨模块复用调查 / 边角情况扫描
- **R5** (1)：R4 边角发现验证 + 最后一轮扫描（产出 1 条新发现）
- **R6** (1)：R5 新发现验证 + 相邻模块扫描（R5 发现被推翻；产出 2 条相邻模块边角）
- **R7** (1)：R6 边角验证 + depth.rs 本体终检

**收敛判据**: 连续一轮在 depth.rs 本体无新发现即停。**R7 宣告收敛** —— depth.rs 本体终检未发现新问题；R6 抛出的 2 条均为 depth.rs 以外辅助模块的低严重度边角，其中一条在 R7 被推翻。

---

## 1. 结论速览

整体评价：**depth 模块在核心算法层是正确的**，两位独立审阅员各自独立核对后，第 1 轮提出的 3 条"高危正确性 bug"（sweep-line tie-break 反向、Phase1 TOCTOU、chunk 边界重复计数）经验证均不成立。真正留下的问题集中在**可维护性、命名/文档歧义、少量边界舍入不一致**与**大规模并发下的 I/O / 内存瓶颈**。

| 维度 | 等级 | 一句话 |
|------|------|--------|
| 算法正确性 | 良好 | 主路径无系统性 bug；留有 2 处小边界隐患 |
| 代码架构 | 待改进 | 单文件 4007 行、单函数 837 行是主要技术债 |
| 设计清晰度 | 中等 | `--use-BFS` / `--samples` / `--merge-tolerance` 语义存在易误解点 |
| 性能与并发 | 中等 | Writer 锁与 stats 内存在大机器 / 大基因组下会成瓶颈 |

---

## 2. 正确性（审阅员 C → 验证员 1）

| # | 发现 | 判定 | 证据 |
|---|------|------|------|
| C1 | sweep-line tie-break 把 start 排在 end 之前导致坐标相切点多计 | **❌ 驳回** | `depth.rs:1141, 3308` 与主循环 `2353-2454`：主循环先 emit 区间再更新 active，tie-break 方向正确 |
| C2 | Phase 1 chunk 并发存在 TOCTOU，导致相邻 chunk 重复发现中间序列 | **❌ 驳回** | `ConcurrentProcessedTracker` (`1348-1397`) 每个 seq_id 独立 `Mutex<IntervalSet>`；Phase 1 并发的是不同 hub 序列的 chunk |
| C3 | `get_unprocessed` 用全局 `total_length` 快速返回 | **⚠️ 潜伏 bug（当前不触发）** | `depth.rs:1369-1371`。R3 深挖：`compute_depth_global` 是唯一调用者，始终传入 `(seq_id, 0, seq_len)`，此时 `total_length ≥ seq_len` 等价于"整个序列已处理完"，语义正确。**一旦未来有调用者传入子区间就会静默漏计**。修复成本低，强烈建议先行修复以锁定接口语义 |
| C4 | hop 1+ 坐标投影未按 strand 分组，正反链合入同一外包围框 | **✅ 确认（中）** | `depth.rs:1728-1743, 2160-2178`。R3 反例：某中间序列 M 在 anchor 上同时有正链 [100,200)→[50,150) 与反链 [800,900)→[60,140) 两段 hop0；bounding box 合并为 query [100,900) → anchor [50,150)。M 中间间隙区 [500,600) 的 hop1 命中会被线性投影到 anchor [100,110)，而该坐标**根本未被 hop0 覆盖**。误差量级可达 anchor_len（数百 kb 级，对倒位+重复尤其严重） |
| C5 | 线性插值舍入不一致：`as i64` 截断 vs `.round() as i64` | **✅ 确认** | `depth.rs:1585-1593` vs `1104-1109`：两处语义相同但一个截断一个四舍五入，差值 ≤1 bp。建议复用 `map_target_to_query_linear` |
| C6 | 5 MB chunk 边界 BFS 跨段重复计数 | **❌ 驳回** | `depth.rs:3076-3095` 顺序切分；BFS 内部 `clipped_target_start/end` (`1934-1935`) 裁剪到 chunk 范围 |
| C7 | `is_self_alignment` 未过滤完全自比对 (`query_id == target_seq_id`) | **⚠️ 确认为可维护性问题** | `depth.rs:1155-1161` 条件语义混乱，但所有调用点（`1760, 2200, 2272`）都被上层 `query_sample_id == anchor_sample_id { continue }` 冗余拦截，运行时行为正确；**建议简化为 `query_sample_id == anchor_sample_id`** 以消除歧义 |

**可落地的正确性修复（优先级顺序）**
1. **C3** — 快速返回条件改为"该局部区间已被完全覆盖"检测。
2. **C4** — `seq_anchor_coverage` key 改为 `(query_id, is_reverse)`；或退化为 hop 1+ 时跳过投影。
3. **C5** — 统一舍入策略；抽出 `map_target_to_query_linear` 共用。
4. **C7** — 简化 `is_self_alignment`，去掉多余的 `query_id != target_seq_id` 子条件。

---

## 3. 设计（审阅员 A → 验证员 2）

| # | 发现 | 判定 | 严重度 |
|---|------|------|--------|
| A1 | `--samples` / `--samples-file` 在全局模式被静默忽略 | ✅ | 中 — 应在入口检测并 warning 或报错 |
| A2 | `--use-BFS` + `--transitive-dfs` 组合时实际走 CIGAR-**DFS**，命名误导 | ✅ | 低 — 建议重命名为 `--cigar-precise` 或 `--use-cigar` |
| A3 | `hub_threshold = ceil(max_degree/2)` 仅对接近星形拓扑有保证 | ⚠️ | 低 — 有 `--ref` 作逃生口，文档补充说明即可 |
| A4 | `--merge-tolerance` 合并时 depth 取 max，破坏"每位点准确计数"语义 | ✅ | 低-中 — 需在 CLI 帮助与文档中醒目警告；可考虑 mean/min 备选 |
| A5 | 全局模式下 `--transitive-dfs` 被忽略 | ❌ 驳回 | — `config.transitive_dfs` 实际通过 `is_transitive` 参与路径选择 (`depth.rs:2497`) |
| A6 | region vs global 两种输出格式完全不同 | ✅ | 低 — 文档中补一节"格式差异说明" |

---

## 4. 架构（审阅员 B → 验证员 2）

| # | 发现 | 判定 | 严重度 |
|---|------|------|--------|
| B1 | 单文件 4007 行，混合 6+ 职责（CLI / BFS / sweep-line / tracker / 输出 / 解析） | ✅ | 中 — 建议拆成 `depth/{bfs,sweep_line,tracker,output,parse}.rs`，每子模块 ≤ 500 行 |
| B2 | `compute_depth_global` **837 行**，内含 **130 行 `write_results` 捕获闭包**（`depth.rs:2483-3319, 2709-2840`） | ✅ | 中 — 闭包捕获 10+ 外部变量；提取为独立 `run_phase1` / `run_phase2` / `write_results` 函数 |
| B3 | `query_region_depth` **434 行**三路模式混杂（`depth.rs:3487-3920`） | ✅ | 中 — 拆成 `build_alignment_info_{cigar_bfs, raw_bfs, direct}`，主函数只分发 |
| B4 | tree cache 生命周期（`set_tree_cache_enabled` / `clear_sub_index_cache`）散落在 `compute_depth_global` 各处 | ⚠️ | 低 — 注释说明了 RSS 控制原因；可选用 RAII guard 封装 |
| B5 | 两套 sweep-line 事件结构并存：`CompactDepthEvent { sample_id: u16 }` (`1129`) vs `DepthEventMulti { sample: String }` (`3297`) | ✅ | 低-中 — `String` 字段在热路径上引入不必要堆分配；建议统一为泛型 `SweepEvent<SampleKey>` |
| B6 | 并行热路径中的裸 `unwrap()` (`2862, 2915, 3001, 3116`) | ⚠️ 降级 | 可忽略 — 绝大多数是对 `ProgressStyle::template(常量)` 的 unwrap，不会 panic |

**参考基线**：`impg.rs` 3483 行，`main.rs` 5082 行 —— `depth.rs` 行数并非离群值，**真正问题是内聚度，而非总长**。

---

## 5. 性能与并发（审阅员 D，未经第二轮验证，保留为待办）

| # | 问题 | 严重度 | 影响 | 建议 |
|---|------|--------|------|------|
| D1 | `BufWriter` 被 `parking_lot::Mutex` 包住作为写入点（`depth.rs:2611, 2628, 2889, 2960, 3133`） | **低-中（R3 降级）** | R3 确认 worker 先在线程本地 `Vec<u8>` 上完成全部计算，仅最后 `write_all(&buf)` 时短暂持锁做 memcpy，并非在计算期间持锁。写入本身是内存拷贝+缓冲，不是磁盘阻塞 | 仍可优化为 channel + I/O 线程，但非紧急 |
| D2 | sweep-line `active_alns[sid].retain(|x| x != idx)` 对每个 end 事件做 O(K) 扫描 (`depth.rs:2450`) | 中 | 高覆盖 hub (K 数百) 区段局部 O(K²) | `SmallVec<[usize;4]>` + `swap_remove`，或 `IndexSet` |
| D3 | `DepthStatsWithSamples.intervals` / `depth_intervals` 在 `--stats-combined` 下无上界，克隆完整 sample 名列表 | **高** | 300 样本 × 数百万 intervals → 数十 GB | 流式输出 stats；或只保留分布直方图 |
| D4 | `TRANSITIVE_CHUNK_SIZE = 5_000_000` 魔法常量（`depth.rs:2462`），针对 128 线程调优 | 中 | 小机器内存压力、大机器利用率不足 | 基于 `rayon::current_num_threads()` 与内存估算动态设置，或暴露为 CLI |
| D5 | ~~Phase 2 非 transitive 路径跳过 `clear_sub_index_cache`~~ | **❌ R3 驳回** | `depth.rs:3138-3151` 的注释明确说明这是有意设计：Phase 2 序列小且部分已处理，跨序列复用 sub-index 避免反复反序列化；内存上界 = 文件数 × sub_index_size，不随序列数增长 | 无需修改 |
| D6 | stdout 路径嵌套双 `BufWriter` (`depth.rs:2616-2618`) | 低 | 略增 syscall 开销 | stdout 时直接用 `std::io::stdout()` |
| D7 | `compute_alignment_degrees` 预扫描禁用 tree cache 后重复加载 sub-index (`depth.rs:1412-1439, 2583`) | 中 | 大基因组 per-file 索引初始化缓慢 | 按 sub-index 文件分组任务，同文件序列复用已加载 tree；或提供 `--skip-hub-detection` |

> D 系列多属"大机器 / 大基因组"场景下的瓶颈，**建议在实际 profile 后再动手**（先测后优），避免投机性优化。D1 经 R3 核查已降级，D3 仍为真实大数据 OOM 风险。

---

## 5a. R3 新增发现

### 新正确性 / 数据完整性条目

| # | 发现 | 严重度 | 证据 |
|---|------|------|------|
| **N1** | `compute_depth_global` 签名无 `sample_filter` 参数，`--samples` 仅在 region 模式生效、**全局模式下静默失效**（R4 修正：R3 的"完全无效"是夸大） | 中 | `depth.rs:2483-2496, 3494`（region 路径有效）；`main.rs:2814-2827` 全局路径调用点不传 |
| **N2** | multi_impg 场景下，**同名但长度不同**的序列在 `seqidx.rs:30-32` `or_insert` 中被静默丢弃后来的长度。`compact_lengths.get_length(seq_id)` 可能返回错误值，进而影响 Phase 2 的 `(0, seq_len)` 上界、`get_unprocessed` 范围甚至 sweep-line 区间 | **中** | `seqidx.rs:30-32`；消费点 `depth.rs:3009, 3016` |
| **N3** | `CigarOp::new` 中 `range_end - range_start` 从 i64 强转 i32 后再压缩到 29 位（上限 ~537 Mb）。大染色体的自对齐合成 CIGAR 在 `--use-BFS` 路径上可能溢出截断；下游 `adjust_len` 会读到错误长度 | 低 | `impg.rs:1971, 2162`；触发条件：单段 >537 Mb 比对 + `--use-BFS` |

### 新性能条目

| # | 发现 | 严重度 | 证据 |
|---|------|------|------|
| **N4** | 每个 chunk / 每次序列处理都 `Vec::new()` 新建 `buf: Vec<u8>`，高并发下产生大量短命堆分配 | 低 | `depth.rs:2885, 2956, 3129` |
| **N5** | Phase 1 transitive 路径在**每个 5 MB chunk 结束后**都 `clear_sub_index_cache()`，若相邻 chunk 属于同一 hub 序列则反复重装 sub-index；与 Phase 2 的复用策略形成非对称 | 中 | `depth.rs:2893` vs `3138-3151` |
| **N6** | `DepthStats.depth_intervals`（**非** `--stats-combined`）路径同样全量累积 `Vec<(String, i64, i64)>`，每次 push 都 `seq_name.to_string()`。D3 只报告了 combined 路径，普通 stats 路径被遗漏 | 中 | `depth.rs:108, 139` |

---

## 5c. R4-R7 补充发现

### R4 追加（边角正确性 / 防御性）

| # | 发现 | 严重度 | 证据 / 备注 |
|---|------|------|------|
| **E1** | 若 sample 数 > 65535，`i as u16` / `idx as u16` 截断导致 sample_id 混乱 | 低（理论存在 / HPRC 尺度无法触发，属于防御性问题） | `depth.rs:764, 917`；R5 确认 bitvec/SampleBitmap 受同一 `num_samples` 约束、同源 |
| **E3** | `active_alns[sample_id as usize]` 直接索引无 `.get()` 保护 | 低（正常路径不可越界，防御性差） | `depth.rs:2369, 2447, 2450`；`num_samples` 与 event sample_id 同源 |
| **E4** | `--ref` 与 `-r`/`-b` region 模式之间无 clap 互斥，region 模式下 `ref_sample` 被静默忽略 | 中 | `main.rs:2752-2761`（region 路径 ref_sample 参数位传 `None`） |
| **E5** | `--merge-tolerance nan` 被 clap 接受，`diff <= NaN` 恒为 false，效果是"完全不合并"而无警告 | 低 | `depth.rs:531, 1300` |
| ~~E2~~ | ~~TSV 注入（seq/sample 名含 \t/\n）~~ | **❌ R5 驳回** | PAF 解析阶段已按 `\t` 切列，含 `\t`/`\n` 的名字在入库时已被拒 |
| ~~E6~~ | ~~`pangenome_bases as f64 * chunk_frac` Tb 级累加精度损失~~ | **❌ R5 驳回** | 单区间乘 [0,1] 分数，远未触及 f64 ↔ i64 的 2^53 边界 |

### R5/R6/R7 追加

| # | 发现 | 严重度 | 证据 / 备注 |
|---|------|------|------|
| **X1** | stdout 路径同样存在 D6 报告的双层 `BufWriter` 嵌套，文件路径下更明显 | 低 | `depth.rs:2614-2618` |
| **X2** | `faidx.rs` 解析未剥离 UTF-8 BOM，Windows 工具生成的 FAI 第一条序列名会带 BOM 字节导致名称匹配失败 | 低 | `faidx.rs:129`；samtools/htslib 生成的 FAI 无 BOM，正常工作流不受影响 |
| ~~SR1~~ | ~~R5 声称 `SortedRanges::insert` 的 `start` 修改导致 BFS 漏探~~ | **❌ R6 推翻** | `start` 的"吸附"仅在 `min_distance > 0` 且接近阈值时触发，是有意的去碎片合并，返回空集合语义正确 |
| ~~F2~~ | ~~`extract_sample(separator="")` 退化使所有 sample 名变空~~ | **❌ R7 驳回** | clap 默认 `--separator="#"`，空字符串不可达 |

---

## 5b. 测试覆盖专项（R3 审阅）

**诊断：🔴 红灯**

| 项目 | 数量 |
|------|------|
| `depth.rs` 内部 `#[cfg(test)]` / `#[test]` | **0** |
| `tests/` 中调用 `depth` 子命令的集成测试 | **0** |
| 测试数据中 depth 专用 PAF 场景（跨样本 hub、5MB 边界、反向链 depth） | **0** |

现存 `tests/test_transitive_integrity.rs`（10 个测试）全部走 `query` 子命令，与 depth 的 sweep-line、两阶段 hub、BFS 模式切换、sample 过滤等核心逻辑**零交集**。

**关键未测路径**：`sweep_line_depth` / 两阶段 global depth / `ConcurrentProcessedTracker` / raw vs CIGAR BFS / `--transitive-dfs` / 5MB chunk 边界 / 反向链 depth / `--samples` / `--merge-tolerance` / `--ref` —— 无一覆盖。

**建议的最小回归测试集**（R3 审阅员提出，按优先级）：

1. `test_get_unprocessed_subinterval_semantics` — 构造 `total_length > (end-start)` 但 `[start,end)` 有未处理片段的场景，断言不漏计（锁定 C3 接口语义）。
2. `test_sweep_line_reverse_strand_depth` — 反向链 PAF → 全局 depth 数值断言。
3. `test_phase1_phase2_no_double_count` — 3 样本 hub 图，断言总行数 = 3，无区间重叠。
4. `test_is_self_alignment_same_sample_two_contigs` — 同 sample 不同 contig，断言深度不重复（锁定 C7 语义）。
5. `test_merge_tolerance_boundary` — 相邻 depth [10,10,9]，`--merge-tolerance` 0.05 不合并，0.15 合并为一行。
6. `test_samples_filter_global_mode` — 全局模式指定 `--samples A,B`，断言要么正确生效、要么显式 warning/error（锁定 N1 / A1）。
7. `test_hop1_inversion_repeat_projection` — 构造倒位+重复场景，断言 hop1+ 的 anchor 坐标不越出真实覆盖范围（锁定 C4）。

---

## 6. 建议的修复优先级

**P0（正确性 / 低成本 / 应立即修复）**
- **C3** 快速路径按区间判断（当前潜伏，接口一旦被子区间调用就漏计）
- **C7** 简化 `is_self_alignment` 条件
- **C5** 统一线性插值舍入
- **N1 / A1** 全局模式检测 `--samples` 并 warning 或报错
- **N2** multi_impg 同名序列长度冲突时报错或 warning，不得静默丢弃
- **E4** clap 中为 `--ref` 与 `-r`/`-b` 添加 `conflicts_with`，或在 region 路径 warning

**P1（正确性 / 中成本）**
- **C4** `seq_anchor_coverage` 按 `(seq_id, is_reverse)` 分组，或 hop1+ 投影时按分组分别插值
- **N3** `CigarOp::new` 对 >537 Mb 区间检测并拒绝 / 扩展位宽
- **测试红灯**：先补上述 7 条最小回归测试（P0 级别的"质量闸门"）

**P2（设计 / 用户体验）**
- A2 `--use-BFS` 重命名或文档明确说明与 `--transitive-dfs` 组合时实为 CIGAR-DFS
- A4 `--merge-tolerance` CLI 帮助与文档加粗警告 depth=max 语义
- A6 文档补一节"global vs region 输出格式差异说明"

**P3（架构重构，见效慢但收益持久）**
- B2 拆出 `compute_depth_global` 的 `write_results` 闭包与 Phase 函数
- B3 拆分 `query_region_depth` 三路
- B5 统一两套 sweep-line 事件结构（顺带消除 `String` 堆分配）
- B1 整体 `depth/` 子模块化

**P4（性能，待 profile 后动手）**
- **D3 / N6** stats 流式化 —— **真实 OOM 风险**，建议优先 profile
- N5 Phase 1 chunk 内 `clear_sub_index_cache` 改为"同 seq 最后一个 chunk 才清"
- D2 `active_alns[sid].retain` → `SmallVec` + `swap_remove`
- D4 `TRANSITIVE_CHUNK_SIZE` 基于 `rayon::current_num_threads()` 动态化
- D7 `compute_alignment_degrees` 按 sub-index 文件分组任务减少重复加载
- N4 worker buffer 线程本地复用
- D1 writer 锁 → channel + I/O 线程（R3 已降级，收益有限，最后动手）
- D6 stdout 路径去掉双 BufWriter

---

## 7. 已明确的"非问题"（诚实声明）

以下几项在早期轮次被质疑，后续轮次**独立核实后确认不是问题**，记录在此以避免重复质疑：

- **`--approximate` × depth**：R4 确认 depth 子命令的 clap 定义**根本不接受** `--approximate`（该 flag 仅存在于 `query` / `partition` 子命令），clap 会在参数解析阶段直接拒绝。这是有意设计，不是功能缺口。
- **depth 对其他子命令的影响面**：R4 确认 depth 模块是**架构上自洽**的，所有核心函数均为私有，仅通过 `main.rs` 被外部调用，partition/refine/similarity/lace/query 零引用。**修复 depth 不会影响其他子命令**。
- **TSV 注入（E2）**：PAF 解析器在入库阶段已按 `\t` 切列，含 `\t`/`\n` 的序列名根本进不来。
- **SortedRanges::insert 状态不一致（SR1）**：R5 的推理在 R6 被构造反例推翻，行为是有意的去碎片合并。
- **sample_index 非确定性**：`sorted_samples.sort()` 在从 `FxHashSet` 收集后显式排序，sample_id 分配完全确定。

## 8. 仍未完成的工作（诚实声明）

- **性能条目 profile 验证**：D2 / D3 / N5 / N6 的实际影响仍需在真实大基因组数据上测量；纸上分析已到位。
- **CLAUDE.md 已声明的测试要求未被满足**：文档要求"verify approximate mode results to normal mode for accuracy"、"test with both PAF and .1aln alignment formats"，但 depth 子命令的集成测试为零。

---

## 9. 迭代收敛说明

本报告经过 **7 轮迭代**（14 个独立 agent）。新发现增量轨迹：

| 轮次 | 类型 | depth.rs 本体新发现 | 相邻模块新发现 | 被推翻 |
|------|------|---|---|---|
| R1 | 4 并行初审 | ~22 条（C1-7, D1-7, B1-6, A1-6） | — | — |
| R2 | 2 并行验证 | 0 | 0 | C1, C2, C6, A5 |
| R3 | 3 并行深审 + 测试 | N1-N6（6 条）+ 测试红灯 | — | — |
| R4 | 3 并行边角扫描 | E1, E3, E4, E5（4 条防御/边角）+ u16、TSV、ref 互斥 | — | E2, E6 |
| R5 | 终检前扫描 | 0 | — | SR1 (后被 R6 推翻前即被 R5 提出) |
| R6 | 跨模块扫描 | 0 | faidx BOM、extract_sample 空 separator | SR1 |
| R7 | depth.rs 本体终检 | **0** | — | extract_sample 空 separator |

**收敛判据**: 连续一轮在 depth.rs 本体无新发现即停。**R5 起 depth.rs 本体新发现 = 0**。R6 只发现相邻模块边角（faidx BOM）；R7 对 depth.rs 本体做终检确认零新增、并推翻 R6 的一条。综合判定：**审阅已达稳定状态，继续迭代预期只能在非 depth 模块产出边角。**

**最重要的 7 条可立即动手的修复**（按收益/成本比排序）：
1. **C3** 快速路径（<10 行改动，锁死接口语义）
2. **C7** `is_self_alignment` 条件（<5 行）
3. **N1/A1** 全局模式 `--samples` warning（<20 行）
4. **N2** multi_impg 同名长度冲突 warning（<15 行）
5. **E4** clap `conflicts_with` 约束（<10 行）
6. **补 P1 级最小回归测试集**（5b 节列出的 7 条）
7. **C4** `seq_anchor_coverage` 按 `(seq, strand)` 分组（本报告中剩下的唯一"误差量级可达 anchor_len"级别的正确性问题）

---

*审阅流程：R1 (4) → R2 (2) → R3 (3) → R4 (3) → R5 (1) → R6 (1) → R7 (1)，共 7 轮 14 员。每轮独立核对前轮行号与结论，不信任前轮描述。最终收敛于 R7。*
