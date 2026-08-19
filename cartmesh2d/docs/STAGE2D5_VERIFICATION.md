# Stage 2D-5 验证记录

## 状态

**IN PROGRESS — 2D-5A small-cell 检测 / alpha 统计 / 最佳邻居候选已实现并通过当前回归；2D-5B topology-safe agglomeration 尚未开始，因此 Stage 2D-5 不能 CLOSED。**

当前严格停留在 small-cell stabilization；未实现 2D-6 quality/export 或 visualization。

## 阶段拆分

### 2D-5A — Detection / Candidate

已实现：

- `SmallCellPolicy2D`
- configurable `areaFractionThreshold`
- Cut-cell alpha histogram
- deterministic small-cell marking
- topology internal-edge neighbour discovery
- aggregate shared-interface length per neighbour
- deterministic best-neighbour scoring
- prefer non-small neighbour over small neighbour
- shared length / neighbour alpha / neighbour area / topology id tie-break
- explicit unresolved report

默认 small-cell 规则：

`alpha < threshold`

threshold equality 在 centralized tolerance 内视为 stable，不发生边界抖动。

### 2D-5B — Topology-safe agglomeration

尚未开始。

必须在真正修改 topology 前满足：

1. 以 2D-5A 的 candidate graph 为输入；
2. 聚合后的几何边界必须由原 cell boundary edges 合法重构；
3. 聚合不能制造重叠、洞、重复 edge 或 non-manifold edge；
4. 聚合前后总流体面积守恒；
5. 聚合后重新通过 Stage 2D-4 topology audit；
6. 无法安全聚合时必须保持 unresolved/explicit failure，不允许静默删除 tiny cell。

## 2D-5A 人工 sliver fixture

流体由一个 `alpha = 0.1` 的窄 Cut-cell 与一个完整邻居组成。

policy threshold = `0.2`：

- small-cell count = `1`
- unresolved = `0`
- best neighbour = stable full cell
- shared internal edge length = `1.0`

policy threshold = `0.1`：

- small-cell count = `0`

因此阈值行为是严格、确定的。

## deterministic tie-break fixture

中央 tiny cell 左右各有一个完全等价的稳定邻居。

候选的 shared length / alpha / area 全部相同，因此最终按 deterministic topology id 选择较小 ID。

重复运行选择不依赖容器遍历顺序。

## unresolved fixture

单独存在的 tiny fluid polygon 没有 internal owner-neighbour edge。

结果：

- small-cell count = `1`
- unresolved = `1`
- report contains `NoNeighbourCandidate`
- report is not valid

不会把“检测到了 small cell 但无法处理”伪装成成功。

## 自适应真实 fixture

使用：

`64-segment shifted circle -> Quadtree level 4 -> 2:1 -> CutCell2D -> global topology -> small-cell analysis`

圆心平移 `(0.07, 0.03)`，避免边界与 Cartesian 网格偶然过度对齐。

默认 threshold `0.1` 的当前确定性结果：

```text
fluid cells = 56
cut cells = 32
full cells = 24
small cells = 8
unresolved = 0
minimum small alpha ≈ 0.00120311
candidates targeting another small cell = 0
```

Cut-cell alpha histogram：

```text
<=0.01 : 3
<=0.05 : 3
<=0.10 : 2
<=0.25 : 5
<=0.50 : 6
<=0.75 : 5
<=1.00 : 8
```

说明当前分析器能在真实自适应 Cut-cell 网格中识别极小面积分数单元，并为这些 tiny cells 找到稳定邻居候选。

## 当前 CMake / CTest

从仓库根 CMake 开启二维并显式构建全部二维测试 target，当前本地门禁：

```text
cartmesh2d_stage0_geometry_tests .................. Passed
cartmesh2d_stage1_grid_tests ...................... Passed
cartmesh2d_stage2_quadtree_tests .................. Passed
cartmesh2d_stage3_cutcell_tests ................... Passed
cartmesh2d_stage3_quadtree_cutcell_audit_tests .... Passed
cartmesh2d_stage4_topology_tests .................. Passed
cartmesh2d_stage4_quadtree_topology_audit_tests ... Passed
cartmesh2d_stage5a_small_cell_tests ............... Passed

100% tests passed, 0 tests failed out of 8
```

## Stage 2D-5 尚未满足的关闭门槛

- [x] alpha histogram
- [x] threshold marking
- [x] best-neighbour candidate
- [x] deterministic candidate selection
- [x] unresolved explicit report
- [x] real adaptive tiny-cell fixture
- [ ] topology-safe cell agglomeration
- [ ] pre/post stabilization total area conservation
- [ ] post-agglomeration topology audit = PASS
- [ ] no negative area / duplicate edge / non-manifold topology after treatment
- [ ] final current-head root CMake / CTest closure gate

因此当前只能记为 **2D-5A PASS / Stage 2D-5 IN PROGRESS**。
