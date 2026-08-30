# R1-D checkpoint：half-edge-lite 与 patch transaction oracle

## 1. 边界

这是 R1-D 的第一个可提交 checkpoint，不是 R1-D 完成。当前新增了确定性的
edge-incidence/half-edge-lite 数据模型、revisioned patch transaction、完整 patch
boundary lock、面积与 hard-quality 提交门，并在生产 hybrid build 对 construction 与
solver topology 都执行 incidence 审计。

当前事务仍使用 `buildGlobalTopology()` 作为候选拓扑 oracle。它证明事务语义与回滚门，
但还没有证明 patch-local 增量更新与全局结果等价；因此不得登记 R1-D 完成，也不得把
global rebuild profile 记为已经消除。

## 2. half-edge-lite / edge-incidence

`buildEdgeIncidenceStore2D()` 从已通过全局拓扑审计的网格建立：

- 每个 cell edge 生成一个有向 `HalfEdgeLite2D`；
- cell loop 保存 `previous/next`；
- internal edge 的两个 incidence 保存互反 `twin`；
- edge incidence 必须与 owner/neighbour 完全一致；
- internal twin 必须反向；
- stable edge key 使用 construction canonical handle；无 registry 的独立 fixture 使用
  确定性的 dense topology vertex ID；
- store 带 revision，供后续 transaction 检查 stale base。

任何 incidence 数不符、endpoint 不符、owner/neighbour 不符、twin 方向不符或 stable
edge key 重复均显式失败。

## 3. patch transaction oracle

`prepareTopologyPatchTransaction2D()` 对排序去重后的 cell patch 锁定：

- patch 与未选中 neighbour 之间的全部 interface edges；
- patch 上的全部 physical boundary edges 及其 patch classification；
- 原 patch 面积与 base revision。

`evaluateTopologyPatchTransactionOracle2D()` 只有在以下条件全部成立时提交：

1. base revision 未过期；
2. 调用者提供的权威 typed hard-quality gate 为 PASS；
3. replacement source identity 唯一且被一对一保留；
4. replacement patch 面积在集中 tolerance 内守恒；
5. 全局 topology oracle 有效；
6. 每条锁定 interface/physical boundary 在候选中原坐标、邻接语义与 patch 类型不变；
7. 候选 edge-incidence 审计有效。

接受时 revision 只增加 1；任何拒绝都返回未改变的 base topology、base incidence 与原
revision。测试覆盖 hard-quality rollback、boundary-lock rollback 与合法 2-cell -> 1-cell
area-preserving commit。

## 4. 验证

- 完整构建：PASS；
- CTest：75/75 PASS，51.90 s；
- topology fixture：coarse/fine shared face、next/previous/twin、稳定 edge keys、损坏
  owner/neighbour fail-closed、revision commit/rollback 全部通过；
- hybrid tests：construction 与 solver incidence 均通过，并满足
  `half_edges = edges + internal_edges`。

真实 narrow-gap：

```text
build-q2a/cartmesh2d_hybrid_cli examples/h4_3/narrow_gap.xy \
  build-q2a/evidence/r1d-incidence-narrow_gap 8 3 8 4 0.012 1.15 1.0
```

- construction：3244 cells，6896 internal edges，14088 half-edges，6896 twin pairs；
- solver：3189 cells，6841 internal faces，13978 half-edges，6841 twin pairs；
- area error：`-3.3750779948604759e-14`；
- independent hybrid reader：PASS，0 overlap，0 non-manifold，fluid area 10.24；
- solver CM2D SHA-256：
  `00d5f20abdaaf20df1aabc4599c089c3c89b28184e713e598314d487beae0505`，与 R1-C
  完全相同，说明本 checkpoint 未暗改网格；
- solver SVG：`build-q2a/evidence/r1d-incidence-narrow_gap.hybrid.solver.svg`；
- Q1 typed contract 仍为 FAIL。本 checkpoint 是 topology transaction 基础设施，不得冒充
  已修复 narrow-gap 的 construction quality。

Docker/OpenFOAM 本 checkpoint 未运行；独立 reader 不能替代真实 `checkMesh`。

## 5. 下一步

1. 以 locked boundary 为输入只重建 patch interior incidence，形成 `TopologyDelta2D`；
2. 与本文件的全局 oracle 逐 vertex/edge/cell/incidence 比较；
3. 记录 patch size、global oracle rebuild count、cache invalidation 与 revision；
4. 在等价门通过后，把 R1-C bounded rephase/resample 的真实 polygon rebuild 接入事务；
5. 候选必须比较 area、feature compatibility 与完整 typed Q1 hard metrics，不能只比较
   short face，也不能降低阈值或删除失败 cell。
