# R1 patch-local transaction 收口

日期：2026-08-30

## 结论

R1 的 narrow-gap 范围在本提交中收口。其目标不是让 Q1 质量合同转为 PASS，而是把 short-face repair 的候选决策路径改造成可证明、可测量的 patch-local transaction：候选循环不得进行全局 topology build 或 full-mesh solver-quality evaluation，唯一 winner 才允许一次 authoritative global oracle，并且 local 判决不得依赖一个可能在恢复未变区域后反转的排序条件。

sharp trailing edge 的 mutable/mutable 路径**没有迁移**，仍为 NO-GO。R1 收口后下一条质量主线应进入 construction-time quality；性能主线必须先补完整的顶层 phase attribution，不能把 `partitionSourcePolygons` 当作已经证实的 45 秒瓶颈。

机器可读证据：`artifacts/r1/closeout-manifest.json`。

## 1. 本轮实际修正

Claude 的 `c3f8650` 已完成主体架构：candidate 只做 patch-local topology delta 与 patch-local quality，选出一个 winner 后再进行 global oracle 与 authoritative full quality。该结构保留。

本轮修正一个未闭合的数学条件。旧 local short-face gate 使用：

```text
(hardCount, maximumSeverity, totalSeverity)
```

做局部字典序严格改善。仅凭这个条件，局部改善不一定在恢复未改变的 patch 外区域后仍是 global 改善。例如：

```text
base local      = (1, 5, 5)
candidate local = (1, 4, 6)
outside max severity = 10
```

局部看 candidate 因 maximumSeverity 从 5 降到 4 而更好；恢复 outside 后，两者 global maximum 都被 10 支配，此时 totalSeverity 反而从 `outsideTotal+5` 变成 `outsideTotal+6`，global 排序反转。

现在的充分条件是：

1. `candidateHard < baseHard`：直接接受 short-face strict improvement；或
2. hard count 不变时，必须同时满足 `candidateMaximum <= baseMaximum` 且 `candidateTotal < baseTotal`。

这样恢复任意未改变的 outside contribution 后，global short-face 字典序仍必然严格改善。最终 authoritative global gate 继续保留，作为 fail-closed 二次验证。

该反例以及 scope-relative ranking / stable cell-id tiebreak 已写成编译期回归断言，生产 comparator 与这些断言共用同一实现，不存在“测试一套、生产另一套”的重复逻辑。

## 2. CI 自身的旧断言修正

第一次独立 GitHub CI 在功能代码全部通过后，停在 Q1 workflow 的 superellipse 特例。原因不是本轮 R1 改动，而是 workflow 同时存在两个互相矛盾的期望：

- 前面的 Python 断言仍要求历史 Q1 superellipse micro-face hard failure 存在；
- 后面的 `generate_q1_baselines.py` 已明确使用 `--expect-superellipse-short-faces absent`，要求 Q2 已经消除该 hard short-face failure。

Q1 baseline 工具本身也明确标注：`present` 用于重现历史 Q1，`absent` 用于验证 Q2 superellipse fix。因此本轮把前面的陈旧断言改成与当前 Q2 基线一致：`face_length/local_h` worst 必须不低于 0.01，且不存在任何 `face_length_over_*` hard issue。没有修改任何质量阈值。

## 3. 独立 GitHub 验证

被验证代码 head：`6068f59da894162bb4547efa1c19f7d92a9f7e5f`

GitHub PR：#2

Actions run：`33315949353`，run #24，结论 **SUCCESS**。PR runner 实际 checkout 的 merge-test commit 为 `be72d66ef647874de52d88b436058443763c7cb4`，即本轮 head 合入 Claude base 后的测试结果。

硬结果：

- Release build：PASS
- 完整 CTest：**75/75 PASS**
- 五个 H4 acceptance case：PASS
- independent OpenFOAM reader：PASS
- Q1 typed dimensionless contract structural verification：PASS
- Q1 scale invariance：PASS
- **OpenFOAM v2606 `checkMesh`：五个 case 全部 `Mesh OK`**
- Q0 deterministic replay：PASS
- Q1 deterministic compact baseline replay：PASS
- acceptance artifact upload：PASS

GitHub artifact：`cartmesh2d-h4-validation`，artifact id `9733481868`，SHA-256 `c7c8ec33286a95226ab01fcb1c3a117df3ce6dca7753f92327c317c484f529a3`。

## 4. narrow-gap R1 硬验收

GitHub runner 生产 case：

```text
solver cells                                      3187
area error                          -3.375077994860476e-14
interface edges                                    592
R1 candidates                                       12
local candidates                                     4
local-quality evaluations                            8
candidate-loop global topology builds                0
candidate-loop full global quality evaluations       0
winner global oracle builds                          4
max oracle builds / transaction                      1
authoritative full quality evaluations               9
accepted transactions                                4
min(face/local_h) before           0.007843137254903904
min(face/local_h) after            0.014619607843136694
local winner matches global authority              true
patch-outside stable IDs unchanged                 true
local delta matches global oracle                  true
solver-quality issue count                           0
```

其中最关键的两个 0 来自 process counter 的实测差值，不是硬编码。

OpenFOAM v2606 对 narrow-gap 报告：

```text
cells: 3187
faces: 13501
internal faces: 6831
Max aspect ratio = 10.32388664 OK
Mesh non-orthogonality Max = 68.81332468
Max skewness = 2.878117895 OK
Mesh OK
```

因此 R1 的 solver-ready 验证现在已有真实 OpenFOAM 外部证据，不再只有自有 reader。

## 5. Q1 状态没有被掩盖

R1 收口并不等于网格质量已经完成。narrow-gap 的 Q1 仍为 **FAIL**，当前 hard issue 为：

| metric | hard count |
|---|---:|
| volume ratio | 295 |
| face weight | 95 |
| minimum interior angle | 20 |
| non-orthogonality | 14 |
| hydraulic aspect | 2 |

solver-quality hard safety 本身为 PASS，OpenFOAM `checkMesh` 为 PASS；Q1 是更严格的 dimensionless 产品质量合同。因此下一条质量主线应直接处理 transition / template / refinement construction 对 volume-ratio 与 face-weight 的形成机制，而不是继续扩大 short-face repair framework。

## 6. sharp-tail 暂不迁移

GitHub runner 上 sharp trailing edge：

```text
R1 candidates = 0
local candidates = 0
global oracle builds = 0
min(face/local_h) = 0.0004202383138623607
Q1 = FAIL
```

当前 mutable/mutable gate 仍然生效。这是有意保留，不是遗漏。迁移前至少需要：

1. 记录真实 patch / one-ring closure cell count 与 k/N；
2. shadow-mode 对候选做 local 与测试专用 global 对照；
3. 证明扩大候选域后不会因 closure 过大失去局部路径意义；
4. 再删除 mutable/mutable migration gate。

因此本轮对 sharp-tail 的正式结论是 **NO-GO**。

## 7. 性能结论收紧

GitHub runner 的 narrow-gap profile：

```text
R1 repair                       1.696 s
solver repair (including R1)    2.120 s
all buildGlobalTopology calls   0.425 s
```

而 acceptance command 的端到端运行仍约为数十秒。现有计时并不能解释绝大部分 wall time，所以本轮**撤回“partitionSourcePolygons 已被证明是主要热点”的说法**。sample stack 可以用于提出假设，但不能替代 phase accounting。

R1 收口不做性能微优化。后续性能任务的硬前置条件是：顶层 phase timer 的累计时间必须能够解释端到端 wall time，再根据实测热点优化。任何直接以 `partitionSourcePolygons` 为既定瓶颈的改动都不应进入主线。

## 8. provenance 边界

Claude 的 macOS/MesaSDK 本地 artifact 曾给出 narrow-gap solver cells = 3185；本轮 GitHub Ubuntu/GCC 独立 runner 给出 3187。两者的核心 R1 count、area/interface、hard short-face 清除和质量约束一致，但 solver topology 并非跨编译器逐字节相同。

因此本轮不宣称 cross-compiler byte identity。正式 R1 closeout authority 是：同一 GitHub runner 内的 deterministic replay、R1 measured counters、stable-ID/oracle gates、independent reader、以及 OpenFOAM v2606 `checkMesh`。

## 9. 收口判定

**R1 patch-local short-face transaction（narrow-gap scope）：CLOSED。**

不再继续扩展 transaction framework。后续分成两条独立工作：

- 质量：construction-time quality，优先 volume ratio / face weight；
- 性能：先完成顶层 phase attribution，再决定优化对象。

sharp-tail migration 保持关闭，直到其专门验证满足上述前置条件。
