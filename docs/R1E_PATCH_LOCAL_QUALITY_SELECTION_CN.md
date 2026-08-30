# R1：short-face 事务的 patch-local quality selection

日期：2026-08-30
基线 commit：`7c0be4ca40e803a6d6caef2a59d4eeab9e97a586`（R1-D checkpoint）
本轮 commit：见 `artifacts/r1/local-quality-selection-manifest.json` 的
`generator_commit`
验证工作树：`.claude/worktrees/sad-fermat-52670c`
构建：`build-r1`，`Release`，`/Applications/mesasdk/bin/c++`，单线程

## 0. 边界

本轮只做 R1 的下一步：把 short-face transaction 的候选筛选改成纯 patch-local。
没有迁移 sharp-tail，没有 rephase/resample，没有新的 Q1 修复算法，没有改动任何
质量阈值。narrow-gap 的 solver 网格、solver-quality JSON 和 quality-contract JSON
与基线逐字节相同，因此本轮不是质量提升，而是决策路径与成本结构的改造。

## 1. 调用链前后变化

改造前（R1-D）：

```text
repairSolverShortFaces2D
  -> evaluateSolverQuality2D(全网格)                     基线
  -> 对每个候选 pair：
       local topology delta
       -> evaluateTopologyPatchTransactionOracle2D
            -> buildGlobalTopology(全部 source cells)     每候选一次
       -> evaluateSolverQuality2D(全网格)                 每候选一次
       -> accept / 继续下一个候选
```

改造后：

```text
repairSolverShortFaces2D
  阶段 1（只做局部，不许 global）
    -> evaluateSolverQuality2D(全网格)                    基线，一次
    -> 对每个候选 pair：
         local topology delta
         -> buildPatchLocalScope2D(patch + 1-ring halo)
         -> evaluatePatchLocalQuality2D(base scope)
         -> evaluatePatchLocalQuality2D(candidate scope)
         -> patchLocalQualityNoWorse2D + local short-face 严格改善
         -> patchLocalRank2D / patchLocalRankBetter2D 更新唯一 winner
  阶段 2（只对 winner）
    -> evaluateTopologyPatchTransactionOracle2D
         -> buildGlobalTopology                            恰好一次
    -> evaluateSolverQuality2D(全网格)                     恰好一次
    -> 与局部判决不一致则 fail closed，不再尝试其他候选
```

关键差别：全局 topology 重建与全局 quality 从「每候选一次」变成「每事务一次」，
且只发生在已经选定的 winner 上。候选阶段的全局成本由 process 计数器实测为 0，
不是断言。

## 2. patch-local quality evaluation

新文件 `include/cartmesh2d/quality/PatchLocalQuality2D.hpp` 与
`src/quality/PatchLocalQuality2D.cpp`。

作用域定义：`buildPatchLocalScope2D` 取被替换的 patch cells 加它们的一圈
neighbour halo。stable identity 来自 `EdgeIncidenceStore2D`；physical boundary
标记直接取 `topology.edges[...].neighbour.has_value()`，不由几何反推。

覆盖的指标，全部来自 commit gate 使用的同一组：

- `face_length / local_h`（无量纲 termination 合同，含 hard count 与 severity）
- face weight
- volume ratio
- interior angle（最小值）
- non-orthogonality
- internal skewness、boundary skewness
- hydraulic aspect
- 附带 concavity、compactness、min face length，与全局报告字段一一对应

公式不另造。为此把 `evaluateSolverQuality2D` 的 per-entity 计算抽成三个
authoritative kernel，全局与局部共用：

- `evaluateSolverCellMetrics2D`
- `evaluateSolverInternalFaceMetrics2D`
- `evaluateSolverBoundaryFaceSkewness2D`

抽取后 narrow-gap 的 `hybrid.solver-quality.json` SHA-256 与抽取前完全相同
(`cb44b8e0…f83e59`)，证明这是纯重构而非公式改写。

正确性依据：一次 boundary-locked patch transaction 只改动 patch 内的 cell 几何，
以及至少一侧在 patch 内的 face。被排除的 cell 和 face 在 base scope 与 candidate
scope 中逐个相同，因此对每个单调聚合量（max 类与 min 类），局部不变差蕴含全局
不变差。这是把 global no-worse gate 换成 local no-worse gate 的唯一理由。

fail-closed 行为：如果某个 patch face 在作用域内只出现一次且未被声明为 physical
boundary，说明 halo 被截断，该 face 的 metric 无法计算——此时评估直接失败，而不是
跳过这条 face。测试
`tests/solver_export_test.cpp` 中删掉一个 halo cell 后断言评估失败。

一致性验证：作用域取整张网格时，局部评估的 10 项聚合量与 `issueCount` 必须与
`evaluateSolverQuality2D` 逐位相等；这是 CTest 里的硬断言，不是抽样比较。

排序（`patchLocalRank2D` / `patchLocalRankBetter2D`）：不同候选的 patch 不同，
作用域也不同，因此两个候选的**绝对**局部聚合量不可直接比较（一个作用域可能本来
就含有更差的邻居）。排序键以「候选相对自身 base scope 的改善量」为主：hard
short-face 数差、severity 差、solver issue 数差与 severity 差；只有改善量完全相同
时才退回绝对值，最后以 `(firstCellId, secondCellId)` 作为与枚举顺序无关的终局
tiebreak。

## 3. identity 约束

未新增第二套 truth source，未改动 stable ID / `StableEdgeKey2D` / tombstone
架构。`identityByPoint` 仍然只用于把「已存在的 exact vertex」迁移到 replacement
polygon 上，其查表结果全部来自 base topology 已有的 stable id。本轮 narrow-gap
路径的 replacement polygon 只有 merged union 与 simplified immutable 两种，其顶点
逐个由 `identityByPoint` 取到 base 已有 id，因此按构造不产生新 interior vertex。
一旦将来产生新点，`buildTopologyDelta2D` 仍然强制要求 producer 提供 typed
`PatchGenerated` key，否则 delta 无效
（`new replacement vertex lacks a unique typed patch identity`）。

## 4. 硬验收逐项

真实生产 run：

```bash
build-r1/cartmesh2d_hybrid_cli examples/h4_3/narrow_gap.xy build-r1/evidence/t1-narrow_gap 8 3 8 4 0.012 1.15 1.0
```

| 验收项 | 结果 | 证据 |
|---|---|---|
| candidate-loop global topology builds = 0 | PASS | `r1_candidate_global_topology_build_count = 0`，由 `globalTopologyBuildCount2D()` 差值实测 |
| candidate-loop full global quality evaluations = 0 | PASS | `r1_candidate_full_global_quality_evaluation_count = 0`，由 `solverQualityEvaluationCount2D()` 差值实测 |
| 每次 transaction 的 winner global oracle builds <= 1 | PASS | `r1_maximum_winner_global_oracle_builds_per_transaction = 1`，`r1_global_oracle_build_count = 4` 对应 4 次 accepted |
| local winner 与最终 authoritative global result 一致 | PASS | `r1_local_winner_matches_global_authority = true`；4 次 winner 全部通过 `legacyNoWorse` + `better` |
| patch 外 stable IDs 不变 | PASS | `r1_patch_outside_stable_ids_unchanged = true` |
| 完整 CTest | PASS | 75/75，56.16 s |
| OpenFOAM 真实 checkMesh | NOT_RUN | docker daemon 未运行（`unix:///Users/Zhuanz/.docker/run/docker.sock` 不存在）；独立 reader 不能替代 |

hard short-face = 0：`face_length_over_local_background_h` 的 hard 违规计数为 0
（12 条剩余违规全部是 `preferred` 级），`worst = 0.0146196 >= 0.01`。

`local_delta_matches_global_oracle = true` 同时保持。

## 5. narrow-gap 前后指标

对比同机、同 Release 构建的基线二进制（`/tmp/r1-base`，commit `7c0be4c`）。

| 项 | 基线 | 本轮 | 变化 |
|---|---:|---:|---|
| solver cells | 3185 | 3185 | 不变 |
| area error | -3.375077994860476e-14 | 同值 | 不变 |
| interface edges | 592 | 592 | 不变 |
| min face/local_h | 0.014619607843137943 | 同值 | 不变 |
| max non-orthogonality | 68.81332467755374 | 同值 | 不变 |
| min face weight | 0.051574510290106325 | 同值 | 不变 |
| min volume ratio | 0.011853263410370605 | 同值 | 不变 |
| min interior angle | 2.465828803420265 | 同值 | 不变 |
| max internal skewness | 3.442394377776235 | 同值 | 不变 |
| max cell aspect | 63.3863048155097 | 同值 | 不变 |
| Q1 hard 违规（face weight / volume ratio / angle / nonorth / skew / aspect） | 95 / 295 / 20 / 14 / 4 / 2 | 同值 | 不变 |
| Q1 status | FAIL | FAIL | 不变 |
| solver CM2D SHA-256 | `6ef16a5d…e0bf60` | 同值 | 不变 |

也就是说面积、patch、确定性、以及全部六项质量指标都没有回归——它们逐字节相同。

候选与成本：

| 项 | 基线 | 本轮 |
|---|---:|---:|
| candidate count | 10 | 12 |
| local-quality evaluations | 0 | 8 |
| candidate-loop global topology builds | 0 | 0 |
| candidate-loop full global quality evals | 未计数（结构上每个通过 oracle 的候选一次） | 0（实测） |
| winner global oracle builds | 4 | 4 |
| authoritative full-quality evaluations | 未计数 | 9 |
| accepted transactions | 4 | 4 |

基线二进制没有这两个 process 计数器，所以「基线的 candidate-loop full quality
evals」只能从代码结构读出（`evaluateTopologyPatchTransactionOracle2D` 成功后立即
`evaluateSolverQuality2D`，即每个通过 oracle 的候选一次），无法给出实测数字。
基线的 `candidate_global_topology_build_count` 之所以已经是 0，是因为 R1-D 把候选
的全局重建记在了 `globalOracleBuildCount` 名下而不是候选计数里——这正是本轮改用
实测计数器的原因：R1-D 的那个 0 是命名产物，本轮的 0 是测出来的。

候选数从 10 升到 12 是预期的：基线在第一个通过的候选上就 return，本轮必须枚举完
全部候选才能确定 rank 最优的唯一 winner。这多出的 2 个候选只花局部成本。

wall time（同机三次，单线程，不可复现值）：

| | run 1 | run 2 | run 3 |
|---|---:|---:|---:|
| 基线 | 45.71 s | 45.73 s | 46.01 s |
| 本轮 | 45.79 s | 45.80 s | 45.80 s |

端到端没有可测的加速，也没有回归。原因由本轮新增的 phase profile 直接给出：

- `r1_repair_seconds ≈ 2.8 s`，其中 winner 的 4 次全局 oracle 是主要部分；
- `build_global_topology_seconds ≈ 0.90 s`，`build_global_topology_call_count = 490`，
  `input_cell_total = 256914`——即 490 次调用平均 524 个 cell，绝大多数是
  legacy agglomeration/repartition 的小 patch 构建，不是全网格重建；
- `solver_profile_source_repair_seconds ≈ 0.95 s`；
- 三项相加远小于 45.8 s。

因此 narrow-gap 的 42–46 s 主耗时不在 short-face 修复，也不在 topology 重建。
`sample(1)` 采样把 356/368 个主线程样本落在 `buildConformalHybridMesh2D`，其中
279 落在 `buildSolverTopology2D`、192 落在 `partitionSourcePolygons`；而
`partitionSourcePolygons` 内部只有约 52+50+20 个样本在 `buildGlobalTopology`。
这与 R1 审计的假设 **[H]** 3 不同：真正的热点是 solver partition 阶段的
polygon/predicate 工作，而不是候选期的全局重建。这一项本轮只是测出来，没有修改。

对照 sharp-tail（`build-r1/evidence/r1-sharp_tail`）：`r1_candidates = 0`，
本轮的迁移门（要求失败面至少一侧是 immutable layer support）没有放行任何候选，
`min face/local_h = 0.00042`，Q1 仍为 FAIL。这正是「本轮不迁移 sharp-tail」的
可测证据，不是遗漏。

## 6. local-quality 与 global authoritative 是否 mismatch

没有 mismatch。

- 结构性一致：作用域取整张网格时，局部评估与 `evaluateSolverQuality2D` 的 10 项
  聚合量与 issue 计数逐位相等（CTest 断言）。
- 生产路径一致：narrow-gap 的 4 次 transaction 中，局部选出的 winner 全部通过了
  之后的 authoritative 全局 no-worse gate，`r1_local_winner_matches_global_authority`
  为 true，且最终网格与基线逐字节相同——基线是「每候选都做全局判决」的路径，
  两者选出了同一个结果。
- 若将来出现不一致，行为是 fail closed：`repairSolverShortFaces2D` 直接返回
  `patch-local winner disagrees with the authoritative global quality result`，
  上层 `buildConformalHybridMesh2D` 把它变成
  `SolverTopologyFailed`，不会静默换候选重试。

## 7. 是否具备迁移 sharp-tail 的条件

**GO**，但有一个前置条件。

支持 GO 的理由：

1. patch-local 评估已覆盖 sharp-tail 需要的全部指标，且与 authoritative 定义同源；
2. 候选选择已经与全局解耦，扩大候选集不再按候选数乘上全局重建成本；
3. fail-closed 通道已存在，迁移出错会显式失败而不是产生更差的网格；
4. 迁移门是 `SolverTopology2D.cpp` 中一处显式判断（要求失败面至少一侧
   immutable），移除它即可让 mutable/mutable patch 进入同一条局部路径。

前置条件：sharp-tail 的失败面是 mutable/mutable，patch 的 1-ring halo 会比
narrow-gap 大，而当前没有 k/N 的实测分布（审计中的 **[H]** 1 仍未关闭）。迁移的
第一步应该是先把 `PatchLocalScope2D` 的 cell 数写进 profile，测出 sharp-tail 的
patch closure 规模，再放开迁移门；如果 closure 接近全网格，局部路径的成本优势会
消失，那时需要先处理 closure 而不是先放开候选。

同时必须清楚：本轮没有改善任何质量指标。narrow-gap 与 sharp-tail 的 Q1 都仍是
FAIL，剩余 hard 违规集中在 volume ratio（295 / 152）与 face weight（95 / 133）。
这些不是 short-face 事务能解决的，属于后续的 construction-time quality 阶段。
