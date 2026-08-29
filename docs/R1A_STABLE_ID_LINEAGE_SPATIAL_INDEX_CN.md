# R1-A：稳定身份、source lineage 与 feature-aware 空间索引

日期：2026-08-30

## 1. 阶段结论

R1-A 已把稳定 vertex record、typed key、source reference、feature owner 和
multilevel spatial index 以 shadow mode 接入 Q2-A shared construction；同时把
source lineage 从 CutCell 贯穿到 global topology、solver partition、source merge 和
local repartition。

生产路径不再用 solver cell centroid 扫描全部 hybrid source polygon。完整扫描仍可通过
CLI 的 `--verify-source-lineage` 显式启用，并在任何 cell 映射不一致时 fail closed。

本阶段不改变 transition/cut polygon 坐标或 owner/neighbour topology，不宣称修复 Q1
短面、face weight、volume ratio 等既有失败；也未运行 OpenFOAM。

## 2. 数据模型与不变量

- `StableVertexId2D` 是 registry 内单调分配的 64-bit ID，不由后续浮点坐标反算；
- `StableVertexKey2D` 明确区分 legacy shadow、grid、source、wall-grid intersection、
  transition 和 patch-generated identity；
- `ConstructionVertexRecord2D` 保存原始/当前位置、local h、feature class/owner、
  source refs、decision/reason 和 creation revision；
- `FeatureVertexIndex2D` 按 support 和 local-h exponent 分桶，查询结果按 stable ID 排序；
- legacy canonicalization 仍是 R1-A 的坐标决策权威，新索引逐查询对照并在结果不同时失败；
- `CutCell2D::sourceLineage` 和 `TopologyCell2D::sourceLineage` 始终排序去重；partition
  子单元继承 lineage，merge/repartition 单元取 lineage 并集；
- production metadata lookup 只检查 lineage 中的 source；显式 oracle 才检查全部 source。

## 3. 定向回归

构建目录：`build-q2a`。

```text
cartmesh2d_stageQ2_intersection_registry_tests  PASS
cartmesh2d_stageQ2A_shared_construction_tests   PASS
cartmesh2d_stageS1_solver_export_tests          PASS
cartmesh2d_stageH4_2_conformal_hybrid_tests     PASS
```

新增断言覆盖：stable record 与 canonical vertex 一一对应、索引 shadow query、事件
source refs、clean solver lineage、source agglomeration lineage 并集、partition child
lineage 继承，以及 opt-in full-scan oracle。

随后运行当前 CMake 清单中的完整项目回归：`75/75 PASS`，总 wall time `53.74 s`。
该清单包含项目 e2e writer 与独立读取器，但本轮没有启动 Docker/OpenFOAM
`checkMesh`；`AGENTS.md` 中“73 项”已落后于当前实际 75 项。

## 4. narrow-gap 真实网格

生产命令：

```text
build-q2a/cartmesh2d_hybrid_cli examples/h4_3/narrow_gap.xy \
  build-q2a/evidence/r1a-narrow_gap 8 3 8 4 0.012 1.15 1.0
```

oracle 命令仅多加：

```text
--verify-source-lineage
```

结果：

| 项目 | 值 |
|---|---:|
| construction cells | 3244 |
| solver cells | 3189 |
| boundary-layer cells | 820 |
| interface edges / vertices | 592 / 592 |
| expected / actual fluid area | 10.24 / 10.239999999999966 |
| area error | -3.3750779948604759e-14 |
| source-lineage candidate checks | 3247 |
| full-scan oracle candidate checks | 5189542 |
| mismatched solver cells | 0 |
| production single-run wall time | 44.00 s |
| oracle single-run wall time | 43.06 s |

单次 timing 受系统负载影响，只作为现场记录；候选检查数是确定性的工作量证据。lineage
路径只检查 oracle 的约 `0.0626%` 候选，约减少 1598 倍。

四类 mesh artifact 与 Q2-A shared narrow-gap 基线逐字节相同：

| artifact | SHA-256 |
|---|---|
| construction CM2D | `04f012148ab0da40f421ec23563553275709545a72a217f997b1cbba9e8cf12b` |
| solver CM2D | `00d5f20abdaaf20df1aabc4599c089c3c89b28184e713e598314d487beae0505` |
| construction VTK | `3bd1c25e48f2f58fa8aaaef8191b2314ac07d170284c5f5c1af0e67f148a6bc5` |
| solver VTK | `b9f7f6c7b4dedbefdc2974bd37a54831e45edbf04503fc7eb689c4b224dbc3bb` |

独立 `check_hybrid_mesh2d.py` 读取 construction VTK：3244 cells、3947 vertices、
6896 internal edges、零 overlap、零 non-manifold，最小正 signed area
`2.2656250000003685e-07`，总面积 `10.24`，状态 `valid=true`。

证据入口：

- `build-q2a/evidence/r1a-narrow_gap.hybrid.solver.cm2d`
- `build-q2a/evidence/r1a-narrow_gap.hybrid.quality-contract.json`
- `build-q2a/evidence/r1a-narrow_gap.hybrid.json`
- `build-q2a/evidence/r1a-narrow_gap.hybrid.construction.json`
- `build-q2a/evidence/r1a-narrow_gap.independent.json`
- `build-q2a/evidence/r1a-narrow_gap.hybrid.solver.png`
- `build-q2a/evidence/r1a-narrow_gap-verify.hybrid.json`

## 5. 未完成与下一步

- narrow-gap 完整 Q1 仍为 `FAIL`；minimum face/local h 仍为
  `0.0078431372549030553 < 0.01`，没有隐藏该失败；
- R1-A 的 typed keys 目前对 legacy vertices 是 shadow key，R1-B 才让 exact grid/source
  key 和 feature compatibility 成为构造决策入口；
- half-edge-lite、edge-incidence patch、TopologyDelta/revision 和 patch-local quality
  属于 R1-D，尚未实现；
- 下一步 R1-B：统一 exact identity 与 proximity proposal decision API，补
  sharp/concave/gap-side compatibility matrix，并以 superellipse 与 sharp-tail 真实网格
  验证 feature 不被移动。
