# Q4-1 termination 构造期局部质量选择

日期：2026-08-31

## 结论

Q4-1 只处理 narrow-gap，并保留 Q3-1/2/3 作为 fallback。最终相对 Q3-3 将
volume-ratio hard count 从 239 降到 231，face-weight hard count 从 50 降到
46；hard short-face 保持 0。angle、non-orthogonality、skewness 和 hydraulic
aspect 的 hard count 分别保持 11、11、4、2。

机器可读证据为 `artifacts/q4-1/manifest.json`，网格预览为
`artifacts/q4-1/narrow-gap-after.png`。

## 1. 剩余 hard face 的构造来源

Q3-3 的 239 个 volume-ratio hard face 中，224 个来自 graded termination front
外侧的 remainder termination construction，15 个来自 graded transition strip；
50 个 face-weight hard face 全部来自前者。把 224 个 remainder-derived VR hard
按源 polygon 拆分后为 triangle 62、quad 111、其他 polygon 51。因此主要来源不是
ear-clipping 三角形，而是 termination remainder 的凸四边形/多边形与相邻 remainder
partition 形成的面积和质心跳变。

Q4-1 只选择这一构造点：对唯一最差 termination-adjacent hard face 的有界一环，
比较“保留两个 source cell”和“提交其严格凸 union”。没有新增对角线、没有移动 wall、
没有更改 pure Cut-cell fallback，也没有为 sharp-tail 启用该能力。

## 2. pre-commit candidate 机制

- 候选在 unified H4 topology commit 前产生；graded termination strip 作为 immutable
  context，保证 interface face 具有真实两侧 incidence。
- 每轮最多检查最差 face 两端的一环；整个 build 最多 8 轮。本例 4 轮接受后收敛。
- stable vertex identity、boundary lock、`TopologyDelta2D`、patch-local authoritative
  quality kernels 和 source lineage 均复用现有实现。
- candidate loop 的 global topology build 和 full global quality evaluation 均为 0；
  每个局部 winner 最多一次 construction-context global oracle。
- 临时 construction context 的 inner boundary 尚未接上 H4-1 immutable cells，因此
  不把它冒充最终物理边界；局部 hard gate 和 winner topology oracle仍为 authoritative，
  最终完整 H4 topology 继续执行 unchanged solver-quality gate。
- Q1 hard 阈值未改变。局部排序要求 VR/FW 各自不增加且合计严格下降，short-face、
  angle、non-orth、skewness、aspect 不恶化。

## 3. 指标

| metric | Q3-3 baseline | Q4-1 + Q3 fallback | change |
|---|---:|---:|---:|
| solver cells | 3153 | 3149 | -4 |
| volume-ratio hard | 239 | 231 | -8 |
| face-weight hard | 50 | 46 | -4 |
| hard short-face | 0 | 0 | unchanged |
| hard minimum angle | 11 | 11 | unchanged |
| hard non-orthogonality | 11 | 11 | unchanged |
| hard skewness | 4 | 4 | unchanged |
| hard hydraulic aspect | 2 | 2 | unchanged |
| minimum volume ratio | 0.01712166482 | 0.01712166482 | unchanged |
| minimum face weight | 0.05740815582 | 0.05740815582 | unchanged |
| maximum non-orthogonality | 68.81332468 | 68.81332468 | unchanged |
| maximum skewness | 3.442394378 | 3.442394378 | unchanged |

构造上下文中 VR hard 为 168→160，FW hard 为 60→56；4 个既有 context
short-face 保持 4。完整 solver topology 的 hard short-face 为 0。

执行数据：16 generated candidates、11 valid local candidates、22 local-quality
evaluations、5 winner oracle builds、4 accepted constructions；candidate-loop global
builds=0，candidate-loop full-quality evaluations=0，maximum winner oracle/build=1，
authoritative full-quality evaluations=10。本机 construction selection wall time 为
2.239792 s。

## 4. 验收

- area error `-3.375077994860476e-14`；owner/neighbour、boundary patch、interface、
  stable ID 和 source lineage invariants PASS；lineage full-scan oracle mismatch=0。
- deterministic replay 的 hybrid CM2D、solver CM2D、VTK、quality-contract JSON、
  solver-quality JSON 和 construction-quality JSON 逐字节一致。
- CTest 75/75 PASS。
- circle、superellipse、concave-L、sharp-tail、narrow-gap 五个 H4 case 均成功，
  solver quality PASS；独立 OpenFOAM reader 全部 PASS。
- OpenFOAM v2606 对五个 H4 case 和 Q4-1 narrow-gap 共 6 个 case 均为 `Mesh OK`。

## 5. 下一步判断

不建议继续扩同一 convex-union template。本轮在 4 次接受后，第 5 个局部 winner 未能
通过全 construction-context authority，说明安全收益已接近该模板边界。下一步应转向
另一类 H4 termination/transition 构造：优先处理不能形成严格凸 union 的
termination remainder 与 graded-strip interface partition，但仍保持 pre-commit、typed、
bounded，不回到 Q3 post-repair 模板堆叠。

## 6. 补记：declined Q4-1 不得代价整层（2026-08-31 后续修复）

初版 Q4-1 只在 narrow-gap 上验证过。在 sharp trailing edge 上打开
`--q4-termination-construction` 会把整个网格降级为 pure Cut-cell fallback：

```text
mesh_mode=pure_cutcell_fallback fallback_stage=hybrid_candidate
solver_cells=1888 boundary_layer_cell_count=0
hybrid_detail=solver quality remains invalid after constrained repair: issues=8
```

原因不在 Q4-1 的候选或 gate，而在 `buildAutomaticHybridWithConstruction2D` 的
候选族：它只用 requested policy 试 1.45/1.55/1.50 三个 termination growth ratio。
一个通过局部 authority、却在最终 unchanged solver-quality gate 失败的 Q4-1 提交，
因此耗尽全部 hybrid 候选，`buildRobustH4Mesh2D` 直接落到 pure Cut-cell。

Q4-1 是可选质量改进，不是产出 hybrid 网格的前提，所以修复是：当请求了 Q4-1 而
所有启用它的尝试都失败时，先用同一几何、关闭 Q4-1 重试一次，之后调用方才允许
考虑 pure Cut-cell。没有修改任何阈值、gate、壁面几何或 Q3 各遍。

修复后 sharp trailing edge：

| 指标 | 修复前 | 修复后 |
|---|---:|---:|
| mesh_mode | pure_cutcell_fallback | hybrid |
| solver cells | 1888 | 3391 |
| boundary-layer cells | 0 | 1044 |

新增 `q41_construction_selection_declined` 报告位。declined 路径的
`hybrid.cm2d`、`hybrid.solver.cm2d`、`hybrid.quality-contract.json` 与完全不加
flag 的构建逐字节相同，即 declined Q4-1 是真正的 no-op；narrow-gap 的
accepted 路径 solver `.cm2d` 与修复前同为
`e058ec5bca3c84f8b6985ba0752d09a95ee3b5dc9b56c33fd793db84465250cd`。

回归门：`tests/h4_local_termination_test.cpp` 断言 sharp taper 在 Q4-1 policy 下
仍为 Hybrid、上报 declined、且 boundary-layer cell 数与非 Q4 构建相同；CI 增加
`q41_sharp_trailing_edge` case 与其 OpenFOAM `checkMesh` 覆盖。CTest 75/75 PASS。

本机 `opencfd/openfoam-run:2606` 对新的 declined case 与原 accepted case 均为
`Mesh OK`：sharp-tail max aspect 13.912、max non-orth 69.3953455、max skewness
3.34005325；narrow-gap max aspect 11.268、max non-orth 68.81332468、max skewness
3.442394378。

机器可读证据：`artifacts/q4-1/fallback-safety-manifest.json`，其中还记录了
circle / superellipse / concave-L / sharp-tail 四例的 Q3/Q4 适用性实测。

## 7. 剩余未做项（按证据排序）

1. **RemainderCut 家族没有任何 construction-time 处理。** circle 与 superellipse
   的 `termination_cell_count=0`，Q3/Q4 全部 inert，但两者 Q1 仍为 FAIL，且全部
   hard face 都在 RemainderCut：circle volume-ratio 56 / face-weight 24，
   superellipse volume-ratio 10 / face-weight 16。这是覆盖面最广的空白。
2. **Q4-1 模板只对 narrow-gap 有收益。** concave-L 产生 8 个候选、0 接受；
   sharp-tail 0 接受。与 §5 的判断一致：需要新模板，而不是放宽同一个。
3. **BoundaryLayer 六项指标仍为 OBSERVED。** 没有 hard/preferred 数值，因此不参与
   总判级，Q1 的 PASS 定义目前不覆盖近壁网格质量。
4. **顶层 phase attribution 仍未补完。** narrow-gap Q4+Q3 端到端约 84 s，而
   `q3+q32+q33+q41` 计时合计约 37 s；R1F §7 的前置条件未满足，性能优化不应开工。
5. **sharp-tail 的 Q1 hard 数仍高**（volume-ratio 149、face-weight 136），
   Q3-2/Q3-3 只带来个位数改善。
