# Stage 2D-5 验证记录

## 状态

**PASS / CLOSED**

2D-5A small-cell 检测/候选分析与 2D-5B topology-safe agglomeration 已完成。用户已显式执行 `验证-5`；封口时确认 PR head 未发生代码漂移，仍为此前实际通过根工程 9/9 CTest 的同一实现提交。

## 2D-5A

已实现：

- configurable `areaFractionThreshold`
- Cut-cell alpha histogram
- deterministic `alpha < threshold` marking
- topology internal-edge neighbour discovery
- aggregate shared-interface length
- stable-neighbour preference
- deterministic tie-break
- explicit unresolved report

shifted-circle fixture：threshold `0.1` 下 small cells = `8`，minimum alpha ≈ `0.00120311`，unresolved = `0`。

## 2D-5B

已实现：

- `AgglomeratedCell2D` / `AgglomerationResult2D`
- small -> stable-target grouping
- 组内共享 directed edge fragment 消除
- 单闭环 exterior reconstruction
- 共线冗余顶点简化
- member-area / merged-area 一致性
- total fluid-area conservation
- global topology rebuild
- Stage 2D-4 topology audit 再验证
- unsafe small->small / self-target / 断链 / 分叉 / 多环显式失败

真实 shifted-circle 聚合：

```text
small cells detected = 8
merged small cells = 8
output cells = input cells - 8
total fluid-area error <= 1e-10

duplicateVertices = 0
duplicateEdges = 0
orphanInternalEdges = 0
nonManifoldEdges = 0
unclassifiedBoundaryEdges = 0
openCellLoops = 0
areaMismatches = 0
```

## 已执行根工程门禁

```text
cartmesh2d_stage0_geometry_tests .................. Passed
cartmesh2d_stage1_grid_tests ...................... Passed
cartmesh2d_stage2_quadtree_tests .................. Passed
cartmesh2d_stage3_cutcell_tests ................... Passed
cartmesh2d_stage3_quadtree_cutcell_audit_tests .... Passed
cartmesh2d_stage4_topology_tests .................. Passed
cartmesh2d_stage4_adaptive_topology_audit_tests ... Passed
cartmesh2d_stage5a_small_cell_tests ............... Passed
cartmesh2d_stage5b_agglomeration_tests ............ Passed

100% tests passed, 0 tests failed out of 9
```

## 隔离

二维实现仍限制在根二维入口/文档与 `cartmesh2d/**`；未修改三维 `include/cartmesh/**`、根 `src/**`、根 `apps/**`、根 `tests/**` 算法代码。

Stage 2D-5 正式关闭；后续进入 Stage 2D-6 quality/export/final acceptance。
