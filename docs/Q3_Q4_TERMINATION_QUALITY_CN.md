# Q3/Q4：termination 局部质量修复的四个变体

日期：2026-09-02
范围：把 Q3-1、Q3-2、Q3-3、Q4-1 四轮合并成一份记录。Q3-2 与 Q3-3 此前**只有
`artifacts/` 里的 manifest，没有任何文档**；本文补上这一段，证据链因此闭合。

本文不修改任何算法，也不改动任何阈值。所有数字来自各轮当时的 manifest。

## 0. 一句话结论

四个变体都是 narrow-gap 专用的**事后修复**，全部默认关闭、需要显式 CLI flag。
它们把 narrow-gap 的 Q1 hard volume-ratio 从 295 降到 239、face-weight 从 95 降到 50，
但**没有**把任何案例转成 Q1 PASS，而且在 circle / superellipse 上完全 inert。

R1F 第 5 节与 R2 交接文档都已指出：下一条质量主线应该是 construction-time quality，
**不要继续扩大这个修复框架**。本文的作用是把已有结论固定下来，而不是邀请继续加变体。

## 1. 四个变体与 CLI flag

| 变体 | CLI flag | 事务上限 | 候选模式 |
|---|---|---:|---|
| Q3-1 | `--q3-termination-quality` | 32 | `Agglomeration` |
| Q3-2 | `--q3-termination-repartition` | 16 | `Repartition` |
| Q3-3 | `--q3-termination-grouped` | 8 | `GroupedRepartition` |
| Q4-1 | `--q4-termination-construction` | 8 | `ConstructionAgglomeration`（pre-commit） |

四个 flag 是**累积的**，不是正交的（`apps/cartmesh2d_hybrid_cli.cpp`）：
`--q3-termination-repartition` 同时打开 Q3-1；`--q3-termination-grouped` 同时打开
Q3-1 与 Q3-2。因此 Q3-2 与 Q3-3 从未被单独测量过，它们的数字都是**叠加在前一级之上的增量**。

四个模式共用同一个内核 `repairSolverTerminationQuality2D`
（`src/quality/SolverTopology2D.cpp`），通过 `TerminationQualityCandidateMode2D`
分派。模式特有代码集中在候选选择那一段。

## 2. narrow-gap 链条（唯一一条四级都有效的链）

基线（无 flag）：3185 solver cells，hard volume-ratio 295，face-weight 95。

| 阶段 | solver cells | hard VR | hard FW | 候选数 | 接受事务 | 上限触顶 | repair 秒 |
|---|---:|---:|---:|---:|---:|---|---:|
| baseline | 3185 | 295 | 95 | — | — | — | — |
| Q3-1 | 3153 | 265 | 65 | 256 | **32** | 是 | 23.56 |
| Q3-2 | 3153 | 241 | 50 | 79 | **16** | 是 | 11.35 |
| Q3-3 | 3153 | 239 | 50 | 27 | **2** | 否 | 1.44 |
| Q4-1 | 3149 | — | — | — | **4** | 否 | 2.20 |

Q3-2 之后的其余 hard 项在两轮里都没动：minimum angle 11、non-orthogonality 11、
skewness 4、hydraulic aspect 2。

两点必须直说：

- **Q3-3 的实测收益是 -2 个 volume-ratio、face-weight 零变化**，27 个候选里只接受 2 个。
  CI 对 Q3-3 的 face-weight 断言因此被放宽成 `<=`，而 Q3-1 与 Q3-2 用的是严格 `<`。
- Q3-1 与 Q3-2 都**触顶**（32/32、16/16），也就是说它们不是"改完了"，而是"预算用光了"。
  按 `bounded-search-needs-a-feasibility-check` 的教训，触顶不等于问题解决。

每一级的 `q3_*_cell_count_neutral`、`patch_outside_stable_ids_unchanged`、
`local_delta_matches_global_oracle`、`local_winner_matches_global_authority` 均为
true；`solver_quality_issue_count` 均为 0；`area_error` 为 `-3.375e-14`；
OpenFOAM v2606 `checkMesh` 每级均 `Mesh OK`；CTest 当时 75/75。

## 3. 其他几何：基本 inert

| 几何 | Q3-1 | Q3-2 | Q3-3 | Q4-1 |
|---|---|---|---|---|
| narrow_gap | VR 295→265 | 265→241 | 241→239 | 接受 4 |
| sharp_trailing_edge | VR 156→156 | 156→150 | 150→149 | **接受 0**（declined） |
| concave_l | VR 107→59 | — | 49→49 | **接受 0**（8 候选） |
| circle | inert | inert | inert | **0 候选** |
| superellipse | inert | inert | inert | **0 候选** |

circle 与 superellipse 的 `termination_cell_count = 0`，所以这一族修复对它们
**根本不适用**——它们的 hard face 全部在未被处理的 RemainderCut 家族里。
这就是 `q3-q4-only-apply-to-termination-geometries` 记录的事实。

Q3-1 是唯一在第二个几何上产生可观改善的变体（concave-L 107→59）。

## 4. 文档已记录的否决

- **Q4-1**：`不建议继续扩同一 convex-union template`——第 5 个局部 winner 未通过
  construction-context 权威校验，模板的安全产出已经耗尽。
- **R1F**：sharp-tail 的 mutable/mutable 迁移正式结论 **NO-GO**；R1 在 narrow-gap
  范围内 CLOSED，`不再继续扩展 transaction framework`。
- **R2 交接**：`不要继续扩大 Q3/Q4/R1 的 repair framework`。

更根本的一条在 R2 交接文档里：narrow_gap 的 4 层结构总厚 0.0599，而缝宽只有 0.08，
两侧对撞 0.12 > 0.08。**这才是 narrow_gap 需要事后修复的原因。** 如果 R2/W2 的
medial-axis 厚度上限落地，narrow_gap 在构造期就不再产生这些 hard face，这一族修复
也就失去了它唯一被证明有效的用例。

## 5. 删除这一族的代价（如果将来要删）

不是零代价，而且代价在证据里而不在代码里：CI 的 `q41_narrow_gap` 案例同时带
`--q4-termination-construction --q3-termination-grouped`，并且把 Q4-1 的 `q33_*` 键
与独立的 Q3-3 案例对照。删掉 Q3-3 会改变三个 CI 案例的产出网格，使
`artifacts/q3-3/manifest.json` 与 `artifacts/q4-1/manifest.json` 里记录的
`solver_cm2d` SHA-256 失效，workflow 的断言段需要重写。

因此正确的顺序是：**先做 W2 的 medial-axis 上限，确认 narrow_gap 不再产生这些
hard face，再把这一族连同它的基线一起退役。** 那时删除是水到渠成，而不是一次判断。

## 6. 复现命令

```bash
export DYLD_LIBRARY_PATH=/Applications/mesasdk/lib
./build/cartmesh2d_hybrid_cli examples/h4_3/narrow_gap.xy <out>/narrow_gap \
  8 3 8 4 0.012 1.15 1.0 <out>/narrow_gap-case 0.01 --q3-termination-grouped
```

把 flag 换成 `--q3-termination-quality` / `--q3-termination-repartition` 得到
Q3-1 / Q3-2；加 `--q4-termination-construction` 得到 Q4-1。

wall time 只写进 `<prefix>.hybrid.profile.json`（标记 `reproducible: false`），
不进入被逐字节比较的 `.hybrid.json`。
