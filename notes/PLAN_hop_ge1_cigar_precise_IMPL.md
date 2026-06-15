# 落地级实现方案：impg depth --use-BFS 让 hop≥1 也 CIGAR-precise（Part A）

状态：**已实施（2026-06-15）**
日期：2026-06-15
分支：`fix/depth-review-p0-p1`

## 实施记录（changelog）

- `src/impg.rs`：新增 `compose_cigars_shared_b`（纯函数 + `CigarRunEmitter`）、`compose_hop`、`TransitiveRange` / `BfsHit` / `NextCarry` 结构、`MAX_SYNTHESIZED_CIGAR_OPS` 上限（超限 `assert!` 硬中止）。`query_transitive_bfs` 与 `query_transitive_dfs` 在 `carry = store_cigar && !approximate_mode` 时携带 `B→A` CIGAR 并逐跳 slice→transpose→compose，输出 `overlap.2 = anchor` 坐标 + 合成 `A→C`；`carry` 时关闭跨父 frontier 合并。非 carry 路径逐字节不变（结构等价）。
- `src/commands/depth.rs`：`process_anchor_region_transitive_cigar` 现对 hop≥1 自动走 `is_hop0` 分支注册 `CigarEntry`；`project_hop0_coords` / `seq_anchor_coverage` 在该路径变为 vestigial（注释标注，保留供防御）。D-interior 窗口仍由 `project_window_for_sweep` 线性回退（已有注释）。
- `src/commands/depth_checkpoint.rs`：`CKPT_SCHEMA_VERSION` 2→3，拒绝 v2 checkpoint resume（hop≥1 行从线性变 CIGAR-precise）。
- **行为变化（有意）**：`impg query -x` 传递 PAF 输出从「逐跳 hub-target」变为「合成 anchor-target」记录（main.rs PAF writer 本就泛型，无需改）。
- 测试：`src/impg.rs` 11 个 compose 单元测试；`tests/transitive_cigar_compose.rs` 6 个集成测试（forward 2-hop 精确 CIGAR、非 carry 保持 hub、DFS 等价、reverse 方向、**reverse+indel 单碱基回归**、mixed-strand 3-hop 链 XOR 累计）。全套 `cargo test` 通过；e2e 手验 `depth --use-BFS` 与 raw 样本集一致、坐标由合成 CIGAR 走出。

### 真实数据验证（VGP/cmaes，单碱基级）

用 `GCA_003287225.2_vs_*` 星形 PAF（PanSN，eqx CIGAR）验证 `leaf1 → center → leaf2` 真实 2-hop 链。方法：独立用「两次受信单跳投影逐碱基串联」（`M2[M1[p]]`）作为 ground truth，对比 `query -x` 合成 CIGAR 走出的 `leaf1→leaf2` 映射。

- **修复一个真实 bug（风险 §8.1）**：原 `compose_hop` 用 `invert_cigar_ops(carried, anchor_strand)` 转置，当 `anchor_strand=Reverse` 时会**反转 op 顺序**，使 carried 的共享轴 B 方向与 pairwise 的 B 方向相反 → 合成时 B 轴错位、indel 错位（实测 reverse∘reverse 链漂移最多 ~6bp 累积）。单元测试只用 reverse-identity（无 indel）故未暴露。**修正**：转置只 swap I↔D 不反转（`invert_cigar_ops(.., Forward)`），保持 B 前向；当 `anchor_strand=Reverse` 时再把合成输出 op 列表反转以归一到 anchor 前向。
- 修复后：leaf2（reverse∘reverse → strand +）**1401/1401 碱基精确**、leaf3 **1382/1382 碱基精确**，0 mismatch。线性插值相对 ground truth 漂移最多 22bp（均值 ~8.7bp）——证明 CIGAR-precise 确实必要且单碱基正确。

### query 输出合并修复（`merge_adjusted_intervals` overlap-trim 复活，main.rs）

核查 `trim_cigar_prefix` 对 reverse 合成记录是否引入偏移时,发现 `merge_adjusted_intervals` 的 **overlap-trim 分支对正反向都是死代码**:重叠长度算成 `next - current`,真实重叠时恒为负,`>0` 守卫永远失败 → 重叠记录从不合并(只是分开输出,各自精确,无偏移)。**修复**:
- 重叠长度符号改为 `current - next`(forward)/ `current.first - next.last`、`next.last - current.first`(reverse)。
- reverse 分支几何原本整体镜像写反(trim 了 next 而非 current、check 比较端反了、拼接顺序反了)→ 改为:reverse 时 keep next、trim current 的 prefix、`next ++ trimmed_current`,check 交换实参比较 `prefix(current)==suffix(next)`。
- `trim_cigar_prefix` 由浮点比例切 partial op 改为**整数运算**,消除大 op 上 ±1 的舍入偏移风险。
- 重叠合并后补 `merge_consecutive_cigar_ops` 归一。
- 测试:`src/main.rs` 4 个单元测试(reverse 连续/重叠/partial-op-trim、forward 重叠),均以逐碱基 base-map 断言无偏移;重叠 seam 恰为 deletion 时 check 保守不合并(仍无偏移,仅不拼接)。全套 `cargo test` 通过;VGP 复验仍 **1401/1401** 单碱基精确(合并修复未影响合成)。
前置：Part B / hop-0（Method A）已上线（`CigarEntry`/`CigarCursor`/`process_anchor_region_transitive_cigar`/`sweep_line_depth`）。
本文档对应 `notes/PLAN_cigar_precise_depth_positions.md` 的 **方案 A 第二步 / Part A**，并**修正**其 §4.1 的合成向心框架。

---

## 1. Context（为什么做）

`impg depth --use-BFS` 的 `positions` 列输出每个基因组在每个窗口的投影区间。当前：
- **hop-0（直接比对）**：已用 CIGAR 精确（Method A，`overlap.2.metadata == anchor_seq_id` 时注册 `CigarEntry`，sweep 用 `CigarCursor` 精确投影）。
- **hop≥1（传递比对）**：仍是「双重线性近似」——`project_hop0_coords` 先线性把中间序列坐标投回 anchor，sweep 再线性投到 query。实测有数十 bp 漂移。

根因（已核实）：`src/impg.rs` 的 `query_transitive_bfs`(2697) / `query_transitive_dfs`(2425) 的 frontier 只携带 `(seq_id,start,end)`，**每跳丢弃父节点 CIGAR**，且 2971-2989 跨父节点合并 range。所以 anchor→query 的 CIGAR 必须**在搜索过程中边走边合成**，无法事后重建。

目标：`--use-BFS` 下 hop≥1 的 `positions` 坐标由真实合成 CIGAR 逐碱基走出，精确到 bp。默认（非 `--use-BFS`）路径**逐字节不变**（由 `store_cigar=false` 门控）。

已定决策（沿用 `notes/PLAN_..md` §6）：合成仅在 `store_cigar=true` 时发生，不加新开关 → **`impg query -x` 的传递 PAF 输出随之改变（从逐跳 hub-target 记录变为合成 anchor-target 记录），属有意为之**；checkpoint `CKPT_SCHEMA_VERSION` 1→2 以拒绝混用旧 TSV；超长合成 CIGAR 硬报错。

---

## 2. 核心设计：carry target=hub（已验证，修正原 §4.1）

**不变量：每个 frontier item 携带一条 CIGAR，方向为 `target = 当前 hub(B)`，`query = anchor(A)`，记作 `B→A`。**

理由（决定性、可从代码证明）：
- 仅有的切片原语 `project_target_range_through_alignment`(impg.rs:3242，返回切片 CIGAR) 与 `project_target_range_coords_through_alignment`(impg.rs:3386) **都只按 target 范围切**（`target_pos > last_target_pos` 即 break）。**没有 query-range 切片器。**
- 树查询返回的 pairwise CIGAR（`project_overlapping_interval`）天然是 `target=B(hub)`，`query=C`。
- 因此下一跳时，carried `B→A`（target=B）与新 pairwise `B→C`（target=B）**共享 target=B**，二者都能用现成的 target-range 原语切到同一 B 子区间。
- 若按原 §4.1 携带 `anchor→hub`（target=anchor），切到 B 子区间需按 *query* 切 → 需要不存在的 query-range 切片器，且 reverse-strand query 行走是本仓库最易错处。**故拒绝 §4.1 框架。**

**消费端要求 target=anchor**（已核实 depth.rs:3838-3851 + `CigarCursor::project` depth.rs:895 按 `target_delta` 从 `entry.target_start` 行走）。因此每跳合成产出**两条**方向相反的 CIGAR：
- **RESULT** `A→C`（target=anchor, query=C）→ 存入 `overlap.1`，被 `CigarEntry`/`CigarCursor` 消费。使 hop≥1 与 hop-0 在下游同构。
- **NEXT frontier** `C→A`（target=C, query=anchor）→ 携带进下一跳，保持「target=hub」不变量。
- 二者互为转置，用 `invert_cigar_ops` 互转，**不重复跑合成**。

---

## 3. 合成算法（共享 B，输入 `A→B` 与 `B→C`，输出 `A→C`）

「沿共享 target B 合成」无法直接单遍表达，因为合成本质是链接 *一方的 query = 另一方的 target*。共享轴 B 是两者的 *target*，故需先把一方转置使 B 变成它的 query：

**步骤：**
1. **转置** carried `B→A`（target=B）为 `A→B`（target=A, query=B），用 `invert_cigar_ops(ops, anchor_strand)`（impg.rs:172，swap I↔D 且 strand=Reverse 时反转 op 顺序）。strand 作标量单独保留。
2. **合成** `A→B`(query=B) ∘ `B→C`(target=B)，**沿 B 同步行走** → 输出 `A→C`。
3. **NEXT frontier** `C→A = invert_cigar_ops(A→C, strand_AC)`。

**为何不直接携带 `A→B`：** 那样下一跳 frontier 分裂时按 C 子区间再切，又需 query-range 切片器。携带 `B→A`（target=hub）使分裂坐标轴(C)恰是携带 CIGAR 的 target 轴 → 现成原语直接可切（见 §5.5）。每跳一次转置是正确权衡。

**行走与 case table（关键：按 `(consumes_A, consumes_C)` 二元组定 op，避免符号错误）：**

把每个行走步映射到 `(consumes_A, consumes_C)`：
- `(1,1)` → 两边都 `=` 则 `=`，否则 `X`
- `(1,0)` → `D`（A 有碱基，C 缺失）
- `(0,1)` → `I`（C 有碱基，A 缺失）
- `(0,0)` → **不发出任何 op**（B 碱基在两侧都被删，有损 I∘D 情形）

双游标 `i`(走 `A→B`)、`j`(走 `B→C`)，各维护当前 op 剩余长度。B-consuming op 驱动同步，取 `min(remaining)` 推进。

**顺序规则（保证规范、可测）：** 每个 B 位置先 flush 所有「消耗 A 不消耗 B」的 `A→B` op（发 `D`），再 flush 所有「消耗 C 不消耗 B」的 `B→C` op（发 `I`），再处理消耗 B 的重叠。

**strand：** `strand_AC = anchor_strand XOR strand_BC`。发出的 `A→C` op 无方向（纯 =/X/I/D 游程），方向由 `strand_AC` 标量随附（与 hop-0 一致，depth.rs:3846）。

**游程合并 + 分块：** 累加进 `last_op`，同类延长、异类 flush；**flush 必须用 `CigarOp::new_run(len, op)`（impg.rs:105）**，绝不直接 `CigarOp::new`（impg.rs:91 会 assert，长游程 panic）。

**anchor_span 重算：** 因 `(0,0)` 有损情形会缩短 A/C span，结果 `overlap.2` 的 anchor 坐标必须**走实际发出的 `A→C` 的 target deltas 重算**，不能用输入区间算术假定。

---

## 4. 复用的现有原语（无需新转置器）

- `invert_cigar_ops(ops, strand)`（impg.rs:172）= 转置（swap I↔D + strand=Reverse 时反转）。两处用：carried `B→A`→`A→B`（传 `anchor_strand`）；result `A→C`→`C→A`（传 `strand_AC`）。**传错 strand → reverse 链路静默损坏**，需专门 3-hop 混合 strand 测试。
- `project_target_range_through_alignment`（impg.rs:3242）切 carried CIGAR 到 B 子区间 / C 子区间。
- `CigarOp::new_run`（impg.rs:105）、`CIGAR_OP_MAX_LEN`（impg.rs:79）。

**新增：** `compose_cigars_shared_b(a_to_b, b_to_c) -> Vec<CigarOp>`（合成原语）；`TransitiveRange` frontier 结构。

---

## 5. impg.rs 改动清单

### 5.1 frontier 结构（BFS+DFS 共用）
```rust
struct TransitiveRange {
    seq_id: u32, start: i64, end: i64,
    hub_to_anchor_cigar: Arc<Vec<CigarOp>>, // B→A，target=hub
    anchor_strand: Strand,                   // A 相对 hub 的累计 strand
    anchor_id: u32,
    anchor_span: (i64, i64),                 // 投影自的 anchor 区间
}
```
用 `Arc` 让 frontier 分裂共享，再切时产生新 `Vec`（不原地改）。

### 5.2 base case / hop-0 identity（BFS init 2757-2761；DFS 2461-2485）
seed：`hub_to_anchor_cigar = =`-run over `[range_start,range_end]`，`anchor_strand=Forward`，`anchor_id=target_id`，`anchor_span=(filtered_start,filtered_end)`。
hop-0 = compose(identity, pairwise) = pairwise，与今日逐字节一致（写 golden test 断言）。

### 5.3 per-hop 合成（BFS par_iter 体 2784-2861；DFS into_par_iter 2520-2589）
拿到 pairwise `(query=C, cigar=B→C, target=B-sub, overlap=[os,oe])` 后：
1. 用 `project_target_range_through_alignment` 把 carried `B→A` 切到 `(os,oe)`（record=`(hub_start,hub_end,anchor_q_start,anchor_q_end,anchor_strand)`）。
2. 转置切片 → `A→B`。
3. compose `A→B ∘ B→C` → `A→C`(RESULT)；转置 → `C→A`(NEXT)。
4. `local_results` 新增携带：`A→C` cigar、`anchor_id`、重算的 anchor span、`C→A` cigar + `strand_AC`。

### 5.4 results push（BFS 2893-2905；DFS 2609-2611）
`overlap.1 = A→C`；`overlap.2 = Interval{anchor 坐标, anchor_id}`。
→ 使 depth.rs:3744 的 `is_hop0`(`target_interval.metadata == anchor_seq_id`) **对每个 hit 都为真**，hop≥1 走与 hop-0 同一分支拿 `CigarEntry`；`project_hop0_coords`/`seq_anchor_coverage`/`inverse_map_query_to_target` 在 `--use-BFS` 路径变为 dead，可绕过。
**reverse 归一化**：设 `overlap.0.first/last` 使 `first>last` ⟺ `strand_AC==Reverse`（否则 depth.rs:3789 `is_reverse` 推导错、`CigarEntry.strand` 错）。

### 5.5 frontier 分裂再切（BFS 2953-2961；DFS 对应处）
`ranges.insert(...)` 返回 C 上的非重叠子区间 `(new_start,new_end)`。携带的是 `C→A`（target=C）→ **用 target-range 原语把 `C→A` 切到每个 `(new_start,new_end)`** 再 push。（这是携带 target=hub 的收益。）

### 5.6 关闭跨父合并（BFS 2971-2989；DFS stack merge 2675-2690）
`if !store_cigar { ...merge... }`。`store_cigar=true` 时跳过（合并会孤立携带 CIGAR）。`visited_ranges.insert` 去重仍跑 → 发现样本集不变（需做 §7 等价性检查）。

### 5.7 超长合成 CIGAR 硬报错
合成结果 op 数 / 长度超阈值 → 返回 Err 整任务中止（阈值设宽，仅病态深链路触发）。

### 5.8 `project_target_range_coords_through_alignment` 改 `pub(crate)`（若 depth 端需要）

---

## 6. depth.rs 改动清单

- `process_anchor_region_transitive_cigar`（depth.rs:3670）：因每个 hit 的 `overlap.2` 已是 anchor 坐标 + 合成 CIGAR，hop≥1 自动走 `is_hop0` 分支注册 `CigarEntry`（depth.rs:3838）。
- pass-1 `seq_anchor_coverage` 构建、`project_hop0_coords`、hop≥1 的 `inverse_map_query_to_target` 调用在此路径变 dead → 可移除/绕过（函数保留供 raw 路径）。
- `store_cigar` 已是 `true`（depth.rs:3699/3715），无需改。
- `CigarCursor::project` 落在 D 内返回 `None` → `project_window_for_sweep`(depth.rs:1026) 回退线性。**显式决定**：合成 hop≥1 CIGAR 的 D-interior 窗口接受线性回退（与 hop-0 一致），并在注释写明，避免「静默漂移」。

## 6.1 depth_checkpoint.rs
`CKPT_SCHEMA_VERSION` 1→2（depth_checkpoint.rs:44），拒绝旧 checkpoint resume，避免线性行与 CIGAR 行混在同一 TSV。

## 6.2 main.rs
确认 `impg query -x` PAF writer 泛型写 `overlap.2.metadata`，能处理 anchor 坐标（本就泛型，预期无需改，需核对）。

---

## 7. 测试计划

1. `compose_cigars_shared_b` 单元测试：覆盖 `(consumes_A,consumes_C)` 全部 2×2、有损 `(0,0)`、游程合并、`CIGAR_OP_MAX_LEN` 分块。
2. 三序列 A→B→C 链路，两跳都有 indel：合成 CIGAR 与手算一致且优于线性。
3. **reverse 链路 / 多次翻转链路**（最高风险，§8.1）：混合 strand 3-hop，断言坐标方向正确。
4. 超长链路：断言硬报错触发。
5. 窗口落在 D 内 → 线性回退（验证显式行为）。
6. **golden 回归**：非 `--use-BFS` depth 改前改后逐字节一致；hop-0 identity-compose 与今日逐字节一致。
7. e2e：`depth --use-BFS` hop≥1 坐标 vs 手算真值；v1 checkpoint 升 v2 后被拒。
8. **传递样本集等价性**：`--use-BFS` 与 raw 路径发现的样本集合一致（守住关闭 frontier 合并的影响）。

验证命令：
```bash
cargo test                      # 单元+集成
cargo build --release
cargo run --release -- depth -a tests/test_data/<paf> --use-BFS -r <region>  # e2e 手验
```

---

## 8. 风险排序（实施期紧盯）

1. **reverse-strand 方向（compose+transpose 三处交互）** —— `anchor_strand XOR strand_BC`、`invert_cigar_ops` 反转、`query_interval.first>last` 编码 reverse。必须手验混合 strand 3-hop。
2. **`CIGAR_OP_MAX_LEN` panic** —— 长 `=` 累加超 `(1<<29)-1`。flush 一律 `new_run`，绝不 `new`。
3. **有损 `(0,0)` / I∘D** —— 缩短 span，anchor_span 必须走实际发出 op 重算。
4. **空/退化合成区间** —— 切片落在 D 内返回 `None` → 跳过该 hit，不 push 零长或 None-cigar（否则静默回退线性掩盖 bug）。
5. **D-interior 窗口线性回退** —— 消费端 `CigarCursor` 返回 None 回退线性，会重新引入漂移；显式决定+注释。
6. **min/max reverse 归一化** —— `overlap.0` 的 first/last 与 `strand_AC` 一致。
7. **Arc 共享 vs 分裂可变** —— 切片产生新 `Vec`(新 `Arc`)，不原地改。

---

## 9. 成本/收益（沿用 §9）

- 收益：hop≥1 坐标精确到 bp；顺带修正 `query -x` 传递输出为真正 anchor↔query 整体比对。
- 代价：动核心搜索（bug 面大）；`--use-BFS` 变慢/耗内存（全程携带+合成 CIGAR，深链路 CIGAR 随跳数增长，关闭 frontier 合并致冗余探索）；`query -x` 输出格式变化（需写 changelog）；病态深链路硬报错；传递结果边缘微移（靠测试守）。
- 约 80% 风险集中于此 Part A（核心搜索）。

---

## 10. 关键文件索引

- `src/impg.rs`：BFS 2697-3002，DFS 2425-2695，投影原语 3242-3480，`invert_cigar_ops` 172，`CigarOp`/`CIGAR_OP_MAX_LEN` 79-166。新增 `compose_cigars_shared_b` + `TransitiveRange`。
- `src/commands/depth.rs`：`process_anchor_region_transitive_cigar` 3670-3923，`CigarEntry`/`CigarCursor` 845-1004，`project_window_for_sweep` 1011-1033。
- `src/commands/depth_checkpoint.rs`：`CKPT_SCHEMA_VERSION` 44。
- `src/main.rs`：`query -x` PAF writer（核对）。
- `notes/PLAN_cigar_precise_depth_positions.md`：原设计契约（§4.1 框架已被本文档修正为 carry target=hub）。

---

## 11. 实施顺序建议（resume-friendly，应对 compute 时间有限）

1. 写 `compose_cigars_shared_b` + 单元测试（纯函数，无需大编译/数据，最先做、最易测）。
2. 写 `TransitiveRange` + BFS 集成（5.1-5.7），先只跑 BFS（depth 默认 BFS）。
3. depth.rs 绕过 dead 路径（§6）+ e2e 手验。
4. DFS 集成（镜像 BFS）。
5. checkpoint 版本 bump + changelog（`query -x` 输出变化）。
6. 全测试套件 + golden 回归。

每步可独立编译/测试，便于在有限 compute 窗口内分批推进、断点续作。
