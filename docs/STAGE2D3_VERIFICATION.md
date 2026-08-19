# Stage 2D-3 验证记录

## 状态

**PASS / CLOSED — 真正 Cut-cell polygon、解析基准、病态拓扑拒绝、Quadtree 跨单元面积守恒以及根目录集成门禁全部通过。**

当前仍严格停止在 2D-3 局部 Cut-cell 几何层；未开始 2D-4 全局 vertex-edge-cell-neighbour topology、2D-5 small-cell、2D-6 export/quality 或可视化。

## 流体区域语义

当前第一版遵循项目总纲的单个外边界模型：

- `BoundaryLoop` 内部 = 保留流体区域；
- `Inside` leaf -> 完整矩形流体 polygon；
- `Outside` leaf -> 空；
- `Intersected` leaf -> 计算 `leaf AABB ∩ boundary interior` 的真实 polygon。

障碍物/孔洞语义留给后续多边界扩展，不在本阶段隐式反转 inside/outside。

## 已实现

- `CutCellKind::{Empty, Full, Cut, Unsupported}`
- `CutCellIssueCode` / `CutCellIssue2D`
- `CutCell2D`
- `buildCutCell(AABB2D, CellClass, BoundaryLoop)`
- `buildCutCell(QuadtreeLeaf2D, BoundaryLoop)`
- 真实 fluid polygon
- area / centroid / area fraction
- embedded-boundary fragments
- source leaf id/key 传播
- CCW normalization
- tolerance 下重复顶点清理
- 自交/退化输入拒绝
- tangent/corner 零面积接触 -> `Empty`
- multi-component cell intersection -> `Unsupported`

## 解析基准

### 斜线切单位格

单位 cell `[0,1] x [0,1]`，cell 内流体为三角形 `(0,0),(1,0),(0,1)`：

- area = `0.5`
- alpha = `0.5`
- centroid = `(1/3, 1/3)`
- embedded boundary length = `sqrt(2)`

全部通过，误差阈值 `1e-12`。

### 竖直 25% 切割

边界 `x=0.25`，保留左侧：

- area = `0.25`
- alpha = `0.25`
- centroid = `(0.125, 0.5)`

全部通过，误差阈值 `1e-12`。

## 永久回归：Quadtree + Cut-cell 面积守恒

新增 `cartmesh2d/tests/cutcell_quadtree_audit_test.cpp`。

fixture：64 边离散单位圆、domain `[-2,2]^2`、Quadtree boundary level 5、distance band level 4、2:1 balance。

独立审计结果：

```text
leaves = 256
cut = 60
full = 56
empty = 140
unsupported = 0
fluid_area = 3.13655
input_polygon_area = 3.13655
absolute_area_error = 1.33227e-15
invalid_polygon = 0
bad_area_fraction = 0
```

该测试验证：

- 所有有效非空 CutCell 的 `0 < alpha <= 1`；
- fluid polygon 通过 `BoundaryLoop::diagnose()`；
- convex circle fixture 无 `Unsupported`；
- 所有 leaf 流体面积之和与输入边界 polygon 面积一致到浮点舍入量级；
- 2:1 balance 后再构造 Cut-cell 仍保持面积守恒。

## 根目录最终集成门禁

从仓库根 CMake 开启二维并显式构建 2D-0/1/2/3 全部测试目标：

```text
cartmesh2d_stage0_geometry_tests ................. Passed
cartmesh2d_stage1_grid_tests ..................... Passed
cartmesh2d_stage2_quadtree_tests ................. Passed
cartmesh2d_stage3_cutcell_tests .................. Passed
cartmesh2d_stage3_quadtree_cutcell_audit_tests ... Passed

100% tests passed, 0 tests failed out of 5
```

同时 `CARTMESH_BUILD_2D=OFF` 根目录配置继续 PASS。

## 2D-3 验收对照

- [x] 相交 leaf 生成真实 fluid polygon
- [x] polygon area > tolerance
- [x] polygon 非自交
- [x] polygon orientation = CCW
- [x] `0 < area_fraction <= 1` 对非空有效流体 cell 成立
- [x] 解析直线切矩形面积正确
- [x] 解析质心正确
- [x] embedded-boundary fragment 提取
- [x] fragment 端点位于输入边界
- [x] tangent/corner 零面积不制造假 Cut-cell
- [x] invalid boundary 显式失败
- [x] unsupported multi-component 显式失败
- [x] Quadtree wrapper 保留 source id/key
- [x] Quadtree 跨 leaf 总流体面积守恒
- [x] 2D-0/1/2/3 根 CTest 全部通过
- [x] 三维构建默认路径保持隔离

## 分支隔离审计

相对 `main`，二维开发仍只修改：

- 根 `AGENTS.md` 的二维例外说明；
- 根 `CMakeLists.txt` 的可选二维入口；
- `cartmesh2d/**`。

未修改三维算法目录：

- `include/cartmesh/**`
- 根 `src/**`
- 根 `apps/**`
- 根 `tests/**`

## 明确未开始

- 2D-4 global vertex deduplication / Edge2D / owner-neighbour / boundary patch
- coarse-fine edge splitting
- 2D-5 small-cell stabilization/agglomeration
- 2D-6 quality/export
- 2D-V visualization

Stage 2D-3 到此正式关闭。
