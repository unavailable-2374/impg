# `impg depth` 断点续跑设计方案

> Status: 设计稿 (2026-05-03)
> Scope: `src/commands/depth.rs`（全局/ref-anchored/`--stats` 模式）
> 不在范围: `-r`/`-b` 区域查询模式（短任务，单次 BFS，没必要做 checkpoint）；`--stats --combined-output` 在最终阶段做全局排序合并，方案末尾单独说明

---

## 0. 设计目标

| 目标 | 衡量标准 |
|------|---------|
| 进程被 SIGKILL / OOM / 节点重启后，重新运行能从上次最后一次 checkpoint 之后继续 | 重跑总耗时 ≈ 剩余工作时间 + 1×checkpoint 加载（≤ 1 min） |
| **结果与一次跑完字节级一致**（在排序/归一化后） | TSV diff（按 anchor 排序后）为空；`--stats` summary 数字相同 |
| 不破坏现有"无 `--resume` 时性能" | 无 checkpoint flag 时性能损耗 < 2% |
| 状态文件原子更新，崩溃永远不会留下半成品 ckpt | 通过 `tmp + fsync + rename` 保证 |
| 状态文件大小可控 | < 1% 输出 TSV 大小（HPRC 580 样本上估算 < 200 MB） |

显式**非目标**：

- 不保证"任意一行 TSV"的恢复粒度；最小恢复粒度是一个 **batch**（5MB chunk × 数千～数万个）。
- 不保证 Phase 2 内部的输出顺序与一次跑完一致（原本就非确定性，DEPTH_COMMAND.md:316 已明示）。
- 不做"两个不同进程并发对同一 prefix 续跑"的合并，崩了就单进程重启。

---

## 1. 现状回顾（仅列与 checkpoint 相关的事实）

### 1.1 输出路径

`DepthWriter`（`depth.rs:4043-4106`）：
- 工作线程 `w.send(buf)` → 容量 256 的 bounded MPSC → `depth-writer` 线程 → `BufWriter<File>`（8 MB capacity）
- TSV header (`#id\tlength\tdepth\tpositions`) 由 writer 线程在循环开始前写出
- `finish()` 在末尾 drop tx + join，得到 IO 错误

**关键观察**：worker 把整块 4 MB 缓冲交给 writer，writer 顺序写 file。**我们可以注入一个 barrier**：往 channel 里塞一个特殊消息，让 writer 写完该消息前所有积压的 chunk 后，回报当前 `BufWriter` 的字节偏移。

### 1.2 状态对象

| 对象 | 定义位置 | 内容 |
|------|----------|------|
| `tracker: ConcurrentProcessedTracker` | `depth.rs:4295` | 每个 seq 一个 `Mutex<IntervalSet>`；记录"哪些 anchor 区间已被处理" |
| `global_used: Arc<ConcurrentProcessedTracker>` | `depth.rs:4296` | 同样结构；进入 `process_anchor_region_*` 用来防止 BFS 重复计入 |
| `row_counter: AtomicUsize` | `depth.rs:4346` | TSV `#id` 列的递增源 |
| `intervals_counter: AtomicUsize` | `depth.rs:4347` | 跨 chunk 全局区间编号 |
| `processed_count: AtomicUsize` | `depth.rs:4345` | 进度条用，非状态 |
| `stats_accumulator` / `stats_combined_acc` | `depth.rs:4276-4285` | `--stats` 模式下的全局直方图 |

`IntervalSet` 内部是 `BTreeMap<i64, i64>` (`depth.rs:884`)，自然可序列化。

### 1.3 工作单元划分

Phase 1（hub）：
- 非 transitive：`phase1_chunks: Vec<(u32, i64, i64)>`，5 MB 切分，rayon par_iter
- transitive raw：同样切 5 MB，先 `batch_depth_bfs(all phase1_chunks)`（单 batch！），再 par sweep
- transitive CIGAR：per-chunk par_iter

Phase 2（leaf）：
- 非 transitive：`phase2_chunks` → batch of `PHASE2_BATCH_SIZE=65_536` → batch_query + par sweep
- transitive raw：分成 `nontrans_chunks` + `transitive_chunks` 两池
  - `nontrans_chunks`: batch of `PHASE2_NONTRANS_BATCH_SIZE=65_536`
  - `transitive_chunks`: batch of `PHASE2_TRANS_BATCH_SIZE=8_192`
- transitive CIGAR：per-seq par_iter（不切 batch）

**问题**：Phase 1 transitive raw 是**一个**大 batch。崩溃如果落在 Phase 1 中段，没有任何中间状态可恢复。 → 方案里强制 Phase 1 也按 `PHASE1_BATCH_SIZE` 切批。

---

## 2. 总体架构

引入两个新文件（与 `<prefix>.depth.tsv` 同目录）：

```
<prefix>.depth.tsv          ← 主输出，append 写
<prefix>.depth.ckpt         ← 当前已提交的 checkpoint（原子更新）
<prefix>.depth.ckpt.tmp     ← 写入中的暂存（rename 后消失）
```

**单一不变量**：当 `<prefix>.depth.ckpt` 存在且解析成功时，

> `*.depth.tsv` 的前 `ckpt.tsv_byte_offset` 个字节 ≡ "处理完 batch_id ≤ ckpt.last_committed_batch_id 的所有 chunk 后的合法 TSV 输出"。

所有恢复逻辑都从这一条派生出来。

### 2.1 Batch ID 编号

把工作划分成**全局唯一**、**确定性**的 batch 序列：

```
batch_id 取值空间（u64 编码）：
  0..2^32 范围
  低 24 位：phase 内的 batch index
  位 24..28 (4 bit): pass_id（区分 Phase1/Phase2-nontrans/Phase2-trans/Phase2-fallback...）
  位 28..32 (4 bit): phase_id (1 或 2)

排序顺序：phase_id 升 → pass_id 升 → batch index 升
```

每个 phase / pass 在启动时把自己的 chunks 按**确定性**顺序排好（按 `(seq_id, start, end)` 排序），然后切成固定大小的 batch。**这一步必须确定性**，否则恢复时算出来的 batch_id 与上次不同。

### 2.2 Checkpoint 触发点

只有**两类**触发点：

1. **每个 batch 完成后**（worker 全部发送完毕，写线程也写完该 batch）。绝大多数 ckpt 在这里写出。
2. **Phase 1 → Phase 2 边界**（tracker 状态发生质变，单独打一次）。

不在每个 chunk 完成后 checkpoint —— 那是浪费。Batch size 已经是分钟级粒度。

### 2.3 Checkpoint 内容

```rust
#[derive(serde::Serialize, serde::Deserialize)]
struct DepthCheckpoint {
    /// schema 版本，破坏性变更时 +1
    schema_version: u32,

    /// 配置指纹：alignment 文件 mtime/size + DepthConfig + cli flags
    /// 不匹配则放弃恢复（参考 DegreesCache 做法）
    invalidation_hash: u64,

    /// 已提交完成的最大 batch_id（含），由它推出"哪些 batch 跳过"
    last_committed_batch_id: u64,

    /// 此 ckpt 提交时 *.depth.tsv 的字节长度。恢复时把 TSV truncate 到这个长度。
    tsv_byte_offset: u64,

    /// 行号/区间号，让恢复后 #id 连续
    row_counter: u64,
    intervals_counter: u64,

    /// 已完成的 batch 集合（防止跨 phase 顺序混乱时漏跳）
    /// 用排序 Vec 而非 HashSet，便于二分 + 紧凑序列化
    completed_batches: Vec<u64>,

    /// 主 tracker 与 global_used 的快照
    /// 仅持久化非空 IntervalSet，用 (seq_id, intervals) pair 列表
    tracker_snapshot: TrackerSnapshot,
    global_used_snapshot: TrackerSnapshot,

    /// --stats 模式下的累积直方图（仅在 --stats 时填）
    stats_snapshot: Option<DepthStatsSnapshot>,
    stats_combined_snapshot: Option<DepthStatsWithSamplesSnapshot>,
}

#[derive(serde::Serialize, serde::Deserialize)]
struct TrackerSnapshot {
    /// 仅含 intervals 非空的 seq_id
    /// 每条 intervals 用 delta-encoded varint 列表压缩
    entries: Vec<(u32, Vec<(i64, i64)>)>,
}
```

**大小估算**（HPRC 580 样本，CHM13 anchor）：
- tracker：每个 anchor seq 几十～几百个 interval；总条目 ~ 2k seq × 200 interval × 16 byte = 6.4 MB（未压缩）
- global_used：~10× tracker = 60 MB
- 加上其它字段 < 100 MB —— 完全可接受，每次落盘 1～2 秒

### 2.4 写入协议

```
batch 完成路径（编号 K）：
  worker_0..N: w.send(buf_i)         // 已有
  ---- batch barrier ----
  writer.barrier_and_offset() -> off // 新增：drain 已发送，flush BufWriter，返回当前 file.stream_position()
  ---- snapshot ----
  ckpt = build_checkpoint(K, off, tracker, global_used, counters, stats)
  write(ckpt, "<prefix>.depth.ckpt.tmp")
  ckpt_tmp.fsync()
  rename("<prefix>.depth.ckpt.tmp", "<prefix>.depth.ckpt")  // POSIX 原子
  // 注意：rename 之后才视为 checkpoint 已提交
```

`barrier_and_offset` 实现：

```rust
// DepthWriter 新增
enum WriterMsg {
    Chunk(Vec<u8>),
    /// 让 writer 把已收到的所有 Chunk 写完 + flush BufWriter，
    /// 通过 oneshot 把当前底层 File 的偏移发回来
    Barrier(oneshot::Sender<io::Result<u64>>),
}
```

写线程的 `BufWriter` 需要替换成 `BufWriter<File>` 并保留底层 `File` 的句柄；`flush()` 后用 `file.stream_position()?` 拿偏移。

### 2.5 启动 / 恢复路径

```
let resume_state: Option<ResumeState> = if config.resume {
    try_load_checkpoint(prefix)?
        .filter(|ck| ck.invalidation_hash == current_hash)
} else {
    None
};

if let Some(rs) = &resume_state {
    // 1. 截断 TSV 到 rs.tsv_byte_offset（丢弃 ckpt 之后的脏写）
    File::open(tsv_path)?.set_len(rs.tsv_byte_offset)?;
    // 2. 用 append=true 模式打开 DepthWriter，跳过 header
    writer = DepthWriter::open_existing(prefix, rs.tsv_byte_offset)?;
    // 3. 恢复 tracker / global_used / counters / stats
    tracker.restore_from(&rs.tracker_snapshot);
    global_used.restore_from(&rs.global_used_snapshot);
    row_counter.store(rs.row_counter, Ordering::Relaxed);
    intervals_counter.store(rs.intervals_counter, Ordering::Relaxed);
    if stats_mode { stats_accumulator.restore_from(&rs.stats_snapshot); ... }
} else {
    // 全新跑：清掉可能存在的旧 ckpt
    let _ = fs::remove_file(format!("{}.depth.ckpt", prefix));
}

// 各 phase / pass 在迭代 batch 时：
for (batch_id, batch) in enumerate_batches() {
    if let Some(rs) = &resume_state {
        if batch_id <= rs.last_committed_batch_id ||
           rs.completed_batches.binary_search(&batch_id).is_ok() {
            continue;  // 已完成，跳过
        }
    }
    run_batch(batch);
    commit_checkpoint(batch_id);
}
```

**为什么 `tsv_byte_offset` 是 ground truth 而非"重算行数"**：因为 worker 输出顺序非确定，ckpt 之后但崩溃之前可能已经写了**部分 batch 的部分行**到 TSV（writer 线程是单线程顺序写，但接收顺序按 `send` 到达顺序）。截断到上次 ckpt 时刻的偏移，这些"半成品行"被物理丢弃。下次重跑这些 batch 会重新生成对应行。**行号 `#id` 也从 ckpt 恢复**，所以行号连续不冲突。

### 2.6 Phase 1 改造

当前 Phase 1 transitive raw 是单个大 batch（`depth.rs:4566`）。改造：

```rust
const PHASE1_BATCH_SIZE: usize = 8_192;  // 与 PHASE2_TRANS_BATCH_SIZE 对齐

for (batch_idx, batch) in phase1_chunks.chunks(PHASE1_BATCH_SIZE).enumerate() {
    let batch_id = encode_batch_id(phase=1, pass=0, batch_idx);
    if should_skip(batch_id) { continue; }
    let chunk_coords: Vec<_> = batch.iter().map(|...|).collect();
    let all_hits = batch_depth_bfs(impg, &chunk_coords, ...);
    batch.par_iter().zip(all_hits.into_par_iter()).try_for_each(...)?;
    commit_checkpoint(batch_id);  // ← 新增
}
```

性能影响：`batch_depth_bfs` 的"每文件加载一次"摊销在 batch 内有效，把 Phase 1 拆成 ~150 个 batch（HPRC 580 样本，~1200 chunks / 8192 = 1 batch；CHM13 ~50 hub × 50 chunks = 2500 / 8192 = 1 batch）—— 实际上 HPRC 规模 Phase 1 chunks 数量未必能撑满一个 batch，**对小规模无负担**；只有当 chunks 真的超过 8192 时才会切，那时 ckpt 价值最高。

---

## 3. CLI 接口

```
impg depth ... \
    --resume \                       # 启用断点续跑（启动时尝试加载 ckpt）
    --checkpoint-interval-batches N  # 每 N 个 batch 落一次盘，默认 1
    --no-checkpoint                  # 显式禁用 ckpt 写入（即便 --resume）
```

| 行为 | `--resume` 缺省 | `--resume` 启用 |
|------|----------------|----------------|
| 启动时无 ckpt | 全新跑，**不写 ckpt**（向后兼容默认） | 全新跑，开始写 ckpt |
| 启动时有 ckpt（hash 匹配） | 报错 "found stale ckpt at X, refusing to overwrite. Use --resume or remove the file." | 加载 ckpt，从断点继续 |
| 启动时有 ckpt（hash 不匹配） | 同上 | 报错并退出，提示 "ckpt config mismatch, use --no-resume to start fresh" |
| 跑完成功 | 不动 | 删除 `.ckpt`（成功标志） |

`--checkpoint-interval-batches` 默认 1（每 batch 一次）。设置为 4 等可降低 ckpt 频率（更省 IO，崩溃丢更多）。

---

## 4. 关键正确性证明

### 4.1 为什么 `tsv_byte_offset` 截断后能恢复

一次完整运行的 TSV 字节流可分解为：

```
[header][batch_K0_chunks_in_some_order][batch_K1_...]...[batch_Kn_...]
```

在 ckpt 提交瞬间（rename 后），writer 线程已经把"所有 batch_id ≤ K 的 chunks"全部写到 BufWriter 并 `flush` 到 OS pagecache（且 `fsync` 到磁盘——见 4.2）。这之后任何 batch K+1 的 send 都还在 channel 里或 BufWriter buffer 里，崩溃时**最多**有部分 batch K+1 的字节被 OS 异步刷到 disk。

恢复时 `set_len(tsv_byte_offset)` 物理上抛弃这些字节，TSV 文件回到"恰好包含 batch ≤ K"的状态。重跑 batch K+1 后，TSV 末尾追加的内容与上次相同（chunks 集合相同；chunks 内部确定，但顺序可能不同——但本来就是非确定的，符合现状 invariant）。

### 4.2 fsync 顺序

```
1. writer.flush() → BufWriter 把所有字节交给 OS
2. file.sync_data() → OS 把这些字节落到磁盘介质
3. ckpt.fsync(); rename(tmp, dst); dir.fsync();
```

第 2 步必须在第 3 步之前。否则可能出现：ckpt 已落盘说"TSV 字节到 X"，但 TSV 实际只到 X-Δ —— 恢复时截断到 X 反而保留了不存在的字节（`set_len` 会用 0 填充 holes）。`sync_data` 比 `sync_all` 更便宜（不刷 metadata）。

如果性能太差，可降级为只 fsync ckpt——崩溃时若 TSV 末尾少几 MB 字节，截断后恢复行号/tracker，但 set_len 会扩展文件用 0 填充——这会**坏掉** TSV。所以 TSV 的 sync_data 不能省。**优化**：每个 ckpt 的 sync_data 实际只需要刷 ckpt 提交点之前的字节，可以用 `sync_file_range(SYNC_FILE_RANGE_WAIT_AFTER, 0, offset)` 减少全量落盘的开销（Linux 特定）。

### 4.3 同一 chunk 重做不会引入双计数

`global_used`（防 BFS 双计数）在 ckpt 中持久化。重做某 batch 的 chunk 时，**它对 `global_used` 的写入与上次相同**——`mark_processed` 是幂等的。`tracker.mark_processed_batch(&result.discovered_regions)` 同样幂等。因此 Phase 2 看到的"未处理区域"集合与原跑一致。

`row_counter` 和 `intervals_counter` 在 ckpt 中持久化的是**该 ckpt 之前的最终值**。重做 chunk 时它们从该值开始递增，重新分配的 `#id` 与原跑可能不同（因 chunk 内部并行顺序非确定），但在**同一次重启的 lifecycle 内不会冲突**。`#id` 本来就不保证跨 run 一致，可接受。

### 4.4 stats 模式

`--stats` 不写 TSV，写的是最终的 `summary.txt` + 多个 `.depthN.bed` + 可选 `combined.bed`。这些**只在所有 batch 完成后**生成。所以 stats 模式下：

- 不需要 `tsv_byte_offset`
- 把 `DepthStats` / `DepthStatsWithSamples` 累积器序列化进 ckpt
- 重启后从 ckpt 恢复累积器，跳过已完成 batch，继续累积，最后生成 stats 文件

`combined.bed` 在生成前会做全局排序，所以累积器的内部顺序无关，只要内容一致即可。

### 4.5 hash 不匹配的处理

`invalidation_hash` 包含：
- 所有 alignment 文件路径 + mtime + size
- `DepthConfig`（transitive、max_depth、min_transitive_len、use_cigar_bfs 等）
- ref_sample / ref_only / min_seq_length / window_size / merge_tolerance / stats / stats_combined
- `seq_included` 位图（影响 phase1/phase2 划分）
- `degrees` 哈希（degrees 决定 hub 划分；degrees cache 已有自己的 invalidation，这里再叠一层）

hash 不匹配 → 拒绝恢复并报错；用户必须显式 `--no-resume` 或删 ckpt。**不要静默回退到全新跑**——可能会和已存在的 TSV 拼接出乱码。

---

## 5. 实施路径（建议拆 PR）

### PR 1：Phase 1 batching 化（无 resume 行为，纯重构）
- 把 Phase 1 transitive raw 单 batch 改成 `chunks(PHASE1_BATCH_SIZE)` 循环
- Phase 1 非 transitive 已经是 chunk-level par_iter，不动
- 验证：HPRC 580 样本，跑前/跑后 TSV 一致；性能不退化

### PR 2：DepthWriter 加 Barrier
- 引入 `WriterMsg::{Chunk, Barrier}`
- 实现 `barrier_and_offset() -> io::Result<u64>`
- 单元测试：write 几次后 barrier，验证返回字节偏移；再 write，再 barrier

### PR 3：Checkpoint 数据结构 + 序列化
- 新建 `src/commands/depth_checkpoint.rs`
- `DepthCheckpoint` / `TrackerSnapshot` / `DepthStatsSnapshot`
- `IntervalSet::snapshot()` / `restore_from()`
- `ConcurrentProcessedTracker::snapshot()` / `restore_from()`
- 单元测试：序列化往返

### PR 4：CLI flags + 集成
- 在 `Args::Depth` 加 `--resume` / `--checkpoint-interval-batches` / `--no-checkpoint`
- 在 `compute_depth_global` 入口加载/初始化 ckpt
- 每个 batch 循环里加 `commit_checkpoint(batch_id)` 调用
- Phase 1→Phase 2 边界加一次额外 ckpt
- 集成测试：跑到一半 SIGKILL，重启 `--resume`，对比 TSV

### PR 5：fsync + 原子性
- 启用 `sync_data` + 目录 fsync
- 测试：用 `eatmydata` 等工具模拟掉电

### PR 6（可选）：增量 ckpt
- 当前每次 ckpt 写整个 tracker 快照（百 MB 量级）
- 优化为 base + delta 链：每 N 个 batch 一次完整 base，中间写 delta
- 仅在大规模实测发现 ckpt IO 成为瓶颈时再做

---

## 6. 验证清单

- [ ] **正确性**：HPRC 30 样本测试集，对比"一次跑完" vs "跑到 Phase 1 中间 SIGKILL → 续跑" 的 TSV，按 `(seq, start)` 排序后 `diff` 为空
- [ ] **正确性**：同上，Phase 2 中间 SIGKILL
- [ ] **正确性**：`--stats` 模式下，summary.txt 和 combined.bed 跨断点续跑后结果一致
- [ ] **正确性**：hash 不匹配（修改 alignment 文件 mtime）时正确拒绝
- [ ] **性能**：`--resume` 启用但全程不崩溃时，总耗时与无 `--resume` 差距 < 2%
- [ ] **性能**：ckpt 平均落盘耗时 < 2s（HPRC 580 样本）
- [ ] **健壮性**：ckpt 写到一半 SIGKILL（残留 `*.tmp`），下次启动正常恢复到上一个 ckpt
- [ ] **健壮性**：手动删除 `.ckpt`，`--resume` 不报错，从头开始跑
- [ ] **代码质量**：`cargo clippy -- -D warnings` 通过；`cargo test depth_checkpoint` 通过

---

## 7. 已知 trade-off / 不做的事

1. **不做按行恢复**：worker 输出非确定顺序，按行恢复需要给每行打 batch_id 标签，输出格式破坏向后兼容。Batch 粒度已经是分钟级，足够。
2. **不做并行 writer**：当前 writer 单线程是 ordering 保证的源头。改成多 writer 会让 `tsv_byte_offset` 失去意义。
3. **不在每个 chunk 后 ckpt**：会让小批量任务的 ckpt overhead 占比过高。可通过 `--checkpoint-interval-batches=N` 进一步降低频率。
4. **不支持跨节点恢复**（如把 ckpt 拷到另一台机器）：原则上没问题，但需要 alignment 文件路径相同；一致性 hash 会发现 mtime 差异并拒绝。当前不做特殊设计，留给用户用 `cp -p` 自行处理。
5. **`--ref-only` 与 ckpt**：`ref_only` 影响 should_output 判断而非 chunk 划分，可正常工作。
6. **`merge-tolerance` 跨 batch 边界**：当前 merge 逻辑在每个 anchor result 内部做（`split_intervals_by_window` + `windowed_intervals` 循环），不跨 batch，所以 ckpt 恢复后 merge 行为一致。

---

## 8. 文件改动一览（预估）

| 文件 | 修改类型 | 估算行数 |
|------|---------|---------|
| `src/commands/depth.rs` | 修改：Phase 1 batching；DepthWriter 加 Barrier；compute_depth_global 加 resume 路径；每 batch 循环加 commit_checkpoint | +300 |
| `src/commands/depth_checkpoint.rs`（新文件） | 新增：DepthCheckpoint 结构、序列化、原子写、加载、hash 计算 | +400 |
| `src/main.rs` | 加 3 个 CLI flag，传到 `compute_depth_global` | +20 |
| `notes/DEPTH_COMMAND.md` | 加一节"断点续跑" | +60 |
| `tests/depth_resume.rs`（新文件） | 集成测试（SIGKILL + 续跑） | +200 |

总改动 ≈ 1000 行（含测试）。

---

## 9. 一句话总结

**把每个 batch 当作 commit 单元，用 "writer barrier 拿字节偏移 + tracker 快照 + 原子 rename" 三件套提交 checkpoint；恢复时按字节截断 TSV、回放 tracker、跳过已完成 batch。** 现有 chunk-level 划分天然适合这个粒度，唯一硬改动是把 Phase 1 transitive raw 的"单 batch"拆成多 batch。
