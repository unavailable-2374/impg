# 方案：impg depth --use-BFS 的 positions 坐标改为 CIGAR 精确

状态：**待审核**
日期：2026-05-19
涉及仓库：`/scratch/10779/shuocao2374/tool/impg`（分支 `fix/depth-review-p0-p1`）

---

## 1. 背景与问题

`impg depth` 输出 TSV，列为 `#id  length  depth  positions`。`positions` 列对每个窗口给出每个基因组的投影区间 `ASSEMBLY#0#CONTIG:start-end`。

实测发现：即使在 `--use-BFS` 模式下，`positions` 的坐标是**线性插值估算**的，不是 CIGAR 精确的。在一个真实窗口上测得起点偏移 17 bp。

### 1.1 代码层面的根因（已核实）

- `commands/depth.rs:1084` `map_target_to_query_linear` —— 纯线性插值；形参带 `cigar` 但 `#[allow(unused_variables)]` 忽略它。
- 三个 sweep 函数 `sweep_line_depth`(`depth.rs:3696`)、`sweep_line_depth_streaming`(`depth.rs:3836`)、`compute_region_sweep_compact`(`depth.rs:6466`)，写 `positions` 时都调用 `map_target_to_query_linear(&[], ...)`，传空 CIGAR。
- `commands/depth.rs:774` `CompactAlignmentInfo` —— sweep 操作的结构体，只存包围盒（query/target 起止 + is_reverse），**无 CIGAR 字段**。
- **关键**：所有 BFS/DFS/query 调用都传 `store_cigar=false`（已核实 `depth.rs:3484` 等）。所以 CIGAR 在 depth 路径里**根本没被计算**，`overlap.1` 是空 `Vec`。
- BFS 内部（`impg.rs:2893-2905`）：每个结果只携带**单跳**的 CIGAR，`overlap.2.metadata` 是**直接上一级父节点**，不是 anchor。
- BFS 的 frontier（`impg.rs:2959`）只携带裸 `(seq_id, start, end)`，且 `impg.rs:2971-2989` 会把不同父节点的区间合并 —— **接力链路被彻底丢弃，跑完后无法回溯**。

### 1.2 窗口的两种来源

depth 表每个窗口对每个基因组的坐标，分两类：

- **hop-0（直接比对）**：参考与目标基因组之间 PAF 里有直接比对记录。一条 CIGAR，照着走即可。
- **hop≥1（接力/传递比对）**：参考与目标基因组无直接比对，impg 靠接力发现同源：`参考 → 中间基因组B → 目标Z`。链路上有 ≥2 条 CIGAR，**没有一条直接连接 参考↔Z**。

线性插值的 17 bp 漂移主要来自 hop-0；hop≥1 当前是「双重近似」（`project_hop0_coords` 先线性把 hub 坐标投回 anchor，sweep 再线性投到 query）。

---

## 2. 目标

`--use-BFS` 模式下，`positions` 坐标全部由**真实 CIGAR 逐碱基走出**，hop-0 与 hop≥1 都精确到 bp。默认（非 `--use-BFS`）行为保持**逐字节不变**。

---

## 3. 设计总览

统一契约：**当 `store_cigar=true` 时，每个 BFS 结果都携带一条 anchor→query 的 CIGAR，且 `overlap.2` 用 anchor 坐标** —— 这样 hop-0 与 hop≥1 在下游完全同构，sweep 端无需区分跳数。

分两部分：

- **Part A**：CIGAR 合成（解决 hop≥1）—— 改 impg 核心搜索。
- **Part B**：depth sweep 用 CIGAR 投影（解决 hop-0 + 落地）—— 改 depth 输出路径。

---

## 4. Part A —— CIGAR 合成（`impg.rs`）

### 4.1 新增合成原语 `compose_cigars`

```
fn compose_cigars(
    ab: &[CigarOp], strand_ab: Strand,   // anchor(A) -> hub(B)
    bc: &[CigarOp], strand_bc: Strand,   // hub(B)    -> query(C)
) -> (Vec<CigarOp>, Strand)              // 返回 A->C CIGAR + 合成 strand
```

沿共享的 hub 坐标 B 同时走两条 CIGAR。AB 在其 query 侧消耗 B，BC 在其 target 侧消耗 B；合成结果消耗 A（来自 AB.target）与 C（来自 BC.query）。

操作对照表（已核对正确）：

| AB 操作 (A,B) | BC 操作 (B,C) | 合成进 A→C | 备注 |
|---|---|---|---|
| `=`/`X` (1,1) | `=`/`X` (1,1) | 两边都 `=` 则 `=`，否则 `X` | 标准情形 |
| `=`/`X` (1,1) | `D` (1,0) | `D`（A 有碱基，C 缺失） | |
| `=`/`X` (1,1) | `I` (0,1) | `I`（C 有碱基，无 A 对应） | BC 的 I 不消耗 B |
| `D` in AB (1,0) | 任意 | `D`（A 有碱基，C 缺失） | AB 的 D 不消耗 B，独立发出 |
| `I` in AB (0,1) | `=`/`X` (1,1) | `I`（C 有碱基，无 A 对应） | |
| `I` in AB (0,1) | `D` (1,0) | **不发出任何操作** | hub 碱基在两侧都被删 —— 有损情形（见 §10） |
| `I` in AB (0,1) | `I` (0,1) | BC 的 I 独立发出 `I`；AB 的 I 待后续 B 消耗操作配对 | 按确定顺序规则处理 |

- **顺序规则**：每步先 flush 所有前导 AB-`D`（仅 A），再 flush 所有前导 BC-`I`（仅 C），再处理消耗 B 的重叠部分 —— 保证输出规范、可测。
- **游程合并**：经累加器合并相邻同类操作，并按 `CIGAR_OP_MAX_LEN`（`impg.rs:79`）分块，避免逐碱基爆 op。
- **strand**：`strand_AC = strand_AB XOR strand_BC`；合成前把反向跳的 CIGAR 归一到 forward-target 朝向，strand 单独作为标量按 XOR 累乘。

### 4.2 在 BFS/DFS 内增量合成

post-hoc 重建链路不可能（§1.1：frontier 已丢父节点），合成必须在搜索过程中边走边做。

新增 BFS/DFS 私有 frontier 结构：

```
struct TransitiveRange {
    seq_id: u32,
    start: i64, end: i64,             // seq_id 上的坐标
    anchor_cigar: Arc<Vec<CigarOp>>,  // 已合成的 anchor -> seq_id CIGAR
    anchor_strand: Strand,            // anchor..seq_id 的累计 strand
    anchor_span: (i64, i64),          // 投影自的 anchor 区间
}
```

每跳流程：用 `project_target_range_through_alignment` 把携带的 anchor→hub CIGAR 裁到子命中区间 → `compose_cigars` 合成出 anchor→C CIGAR → 结果的 `overlap.2` 改为 **anchor 坐标**。hop-0 是恒等基例（identity ∘ pairwise = pairwise），与纯 hop-0 结果逐字节一致。

### 4.3 关闭 frontier 合并（仅 `store_cigar=true`）

`impg.rs:2971-2989` 的跨父节点合并会破坏 CIGAR 归属，必须在 `store_cigar=true` 时跳过：`if !store_cigar { ...原合并... }`。`visited_ranges` 仍去重，发现的样本集合不变，只是 frontier 略大、结果数/顺序在边缘可能微移（靠测试守住）。默认路径 `store_cigar=false` 走原合并 → **逐字节不变**。

---

## 5. Part B —— depth sweep 用 CIGAR 投影（`commands/depth.rs`）

- 三处 BFS/DFS/query 调用把 `store_cigar` 由 `false` 改 `true`（`depth.rs:3468, 3484, 6671, 6687, 6899`）。
- `CompactAlignmentInfo`（`depth.rs:774`）新增字段 `cigar_idx: u32`（`u32::MAX` = 无 CIGAR → 走线性）。新增构造器，~14 个非 BFS 调用点继续用旧构造器。
- 每个 anchor region 一张侧表 `Vec<CigarEntry>`（CIGAR ops + 真实 target/query span）；处理完即释放，内存按 region 计、不跨 region 累积。
- 因为每个结果都已是 anchor 坐标，CIGAR 路径里删除 `seq_anchor_coverage` pass-1、`project_hop0_coords`、`inverse_map_query_to_target` 的调用（这些函数保留，供非 `--use-BFS` 的 raw 路径继续用）。
- 三个 sweep 函数新增参数 `cigars: &[CigarEntry]`；逐窗口按 `cigar_idx` 分支：有 CIGAR 走 `project_target_range_coords_through_alignment`（`impg.rs:3386`，需改 `pub(crate)`），否则回退线性。
- **单调游标 `CigarCursor`**：每个活跃比对维护 `(op_idx, op_consumed, target_pos, query_pos)`，随 sweep 升序窗口只前进不回退 → 每窗口均摊 **O(1)**，而非每窗口 O(CIGAR 长度)。对合成的 hop≥1 CIGAR 同样适用。
- 窗口边界落在 `D` 操作内 → 投影函数返回 `None` → 该窗口回退线性（保证输出非空）。

---

## 6. 已定决策

| # | 决策点 | 选择 |
|---|---|---|
| 1 | hop≥1 是否也修 | **是**，hop≥1 也要 CIGAR 精确 |
| 2 | 单调游标优化 | **本次就做**，不推迟 |
| 3 | checkpoint schema 版本 | **`CKPT_SCHEMA_VERSION` 1→2**（`depth_checkpoint.rs:44`），格式不变，仅版本号，使旧 checkpoint 在 resume 时被拒，避免线性行与 CIGAR 行混在同一 TSV |
| 4 | 合成是否单独开关 | **不加新开关**。合成只要 `store_cigar=true` 就发生，`impg depth --use-BFS` 与 `impg query -x` 行为一致 |
| 5 | `impg query -x` 输出随之改变 | **接受**。其 PAF 从「逐跳 hub-target 记录」变为「合成 anchor-target 记录」（更正确的传递比对）；非 `--use-BFS`、非传递路径不受影响 |
| 6 | 超长合成 CIGAR 的处理 | **硬报错中止**。合成 CIGAR 超过配置阈值则整个任务报错退出，绝不偷偷近似（阈值需设得很宽，仅病态深链路才触发） |

---

## 7. 完整改动清单

| 文件 | 改动 |
|---|---|
| `src/impg.rs` | 新增 `compose_cigars`；新增 BFS/DFS 私有 `TransitiveRange`；每跳合成；`impg.rs:2971-2989` frontier 合并加 `if !store_cigar` 门控；hop≥1 的 `overlap.2` 改 anchor 坐标；超长合成 CIGAR 硬报错；`project_target_range_coords_through_alignment` 改 `pub(crate)`。`query`（非传递）不动。 |
| `src/commands/depth.rs` | `CompactAlignmentInfo` 加 `cigar_idx`；侧表 `Vec<CigarEntry>`；翻 `store_cigar=true`；CIGAR 路径删除 `project_hop0_coords`/`inverse_map_query_to_target`/pass-1 调用；新增 `CigarCursor`；三个 sweep 函数加 `cigars` 参数并分支。非 `--use-BFS` 路径不变。 |
| `src/commands/depth_checkpoint.rs` | `CKPT_SCHEMA_VERSION` 1→2。 |
| `src/main.rs` | 确认 `impg query -x` 的 PAF writer 能正确处理 anchor 坐标的 `overlap.2`（它本就泛型地写 `overlap.2.metadata`）。 |

**不变量**：默认（非 `--use-BFS`）的 depth 与非 BFS 的 `query` 输出**逐字节不变**，由 `store_cigar=false` 门控保证。

---

## 8. 测试计划

- `compose_cigars` 单元测试：覆盖全部操作对，含有损 `I∘D`、游程合并、`CIGAR_OP_MAX_LEN` 分块。
- 三序列 A→B→C 链路、两跳都有 indel：断言合成 CIGAR 与手算一致，且优于线性。
- 反向链路、多次翻转链路。
- 超长链路：断言硬报错触发。
- 窗口落在 `D` 内 → 线性回退。
- golden 回归：非 `--use-BFS` depth 改前改后逐字节一致。
- e2e：`depth --use-BFS` 的 hop≥1 坐标与手算真值一致；v1 checkpoint 在 v2 升级后被拒。
- 传递结果集合等价性检查：`--use-BFS` 与 raw 路径发现的样本集合一致。

---

## 9. 利弊分析

### 利

- `positions` 坐标精确到 bp，直接窗口与接力窗口都不再有几十 bp 漂移。
- 修好后可直接信赖 depth 表坐标做精确下游分析，不必回 PAF 核对。
- 顺带把 `impg query -x` 的传递输出修正成真正的 anchor↔query 整体比对。

### 弊

- **改动伤面大**：Part A 动 impg 核心搜索（`query_transitive_bfs/dfs`），新增算法与数据结构，bug 面增大。
- **变慢、变耗内存**：搜索全程携带、合成 CIGAR；深链路合成 CIGAR 长度随跳数累加增长；关闭 frontier 合并致搜索冗余探索。`--use-BFS` 的 depth 会比现在慢、占内存多。
- **`impg query -x` 输出格式变化**：任何解析其旧 PAF 的脚本会失效。
- **硬报错风险**：单个病态深链路可能让整个长任务中止（阈值需设宽）。
- **传递结果可能微移**：关闭 frontier 合并使探索范围变化，边缘上可能多/少出结果，依赖测试守住。

### 成本分布（重要）

**约 80% 的收益来自简单的 Part B/hop-0；约 80% 的成本与风险来自 Part A/hop≥1。** 多数 depth 窗口是直接比对（`filter/` 有约 67 万个 PAF，直接比对很密集），接力窗口是少数。

---

## 10. 两种实施路径（待审核选择）

### 方案 A —— 分两步（推荐）

1. **第一步**：只做 hop-0（Part B 中不依赖合成的部分 + 直接比对的 CIGAR 投影 + 游标）。改动小、风险低，直接窗口立即精确到 bp。
2. **第二步**：hop≥1 合成（Part A）作为独立后续，单独评估、单独测试，不拖累第一步上线。

优点：快速拿下大头收益，风险隔离。缺点：接力窗口要等第二步。

### 方案 B —— 一次全做

按本文 Part A + Part B + 游标完整实现，一次到位。

优点：一步到位，所有窗口同时精确。缺点：开发与测试周期更长，`--use-BFS` 运行更慢更耗内存，回归面更大。

---

## 11. 实施期需盯的风险与待确认项

1. **frontier 合并关闭**改变 `--use-BFS` 的 BFS 探索广度 —— 需做样本集合等价性检查，并确认无下游消费者依赖 `results` 顺序。
2. **`impg query -x` 输出变化**虽是有意为之，仍需写入 changelog 提醒用户。
3. **tracepoint vs PAF 的 CIGAR 朝向** —— 确认 `.1aln`/`.tpa` 重建的 CIGAR 与 PAF 的一样满足 `compose_cigars` 假定的 forward-target 契约。
4. **硬报错阈值**取值 —— 需设得足够宽，正常泛基因组深度绝不触发。
5. **有损 `I∘D` 情形**（hub 私有序列在 anchor 与 query 两侧都无对应而被丢弃）—— 语义上正确（从未触及 A 或 C），需确认符合「depth = anchor 碱基覆盖」的预期。
6. **多次 strand 翻转**的链路（A→B(rev)→C(rev) 合成为 forward）—— XOR 逻辑简单，但 `invert_cigar_ops` 的 I↔D 反转与合成顺序交互需人工核一个 3+ 跳用例。

---

## 12. 关键代码位置索引

- `impg.rs:74` `CigarOp` 定义；`:141` `target_delta`；`:149` `query_delta`；`:172` `invert_cigar_ops`
- `impg.rs:2125` `query`（非传递，不动）
- `impg.rs:2425` `query_transitive_dfs`；`:2697` `query_transitive_bfs`
- `impg.rs:2893-2905` 结果 push（单跳 CIGAR + 直接父节点）
- `impg.rs:2959` frontier push（裸三元组）；`:2971-2989` 跨父节点合并
- `impg.rs:3242` `project_target_range_through_alignment`；`:3386` `project_target_range_coords_through_alignment`
- `commands/depth.rs:774` `CompactAlignmentInfo`
- `commands/depth.rs:1084` `map_target_to_query_linear`；`:1192` `project_hop0_coords`；`:1121` `inverse_map_query_to_target`
- `commands/depth.rs:3441` `process_anchor_region_transitive_cigar`
- `commands/depth.rs:3696` `sweep_line_depth`；`:3836` `sweep_line_depth_streaming`；`:6466` `compute_region_sweep_compact`
- `commands/depth.rs:3765 / 3895 / 6517` 三处 `map_target_to_query_linear(&[], ...)` 调用
- `commands/depth_checkpoint.rs:44` `CKPT_SCHEMA_VERSION`
