# Stage 2D-3 验证记录

## 状态

**IN PROGRESS — 真正 Cut-cell polygon 实现、解析数值基准和病态拓扑拒绝逻辑已完成并推送；等待当前 GitHub 完整分支 CMake/CTest 门禁后正式 CLOSED。**

当前严格停留在 Cut-cell 局部几何层；未实现 2D-4 全局 vertex-edge-cell-neighbour topology、small-cell 聚合、求解器导出或可视化。

## 流体区域语义

当前第一版遵循项目总纲的“单个外边界”模型：

- `BoundaryLoop` 内部 = 保留的流体区域；
- `Inside` leaf -> 完整矩形流体 polygon；
- `Outside` leaf -> 空；
- `Intersected` leaf -> 计算 `leaf AABB ∩ boundary interior` 的真实 polygon。

外流场实体/孔洞语义属于后续多边界扩展，不在 2D-3 中偷偷反转 inside/outside。

## 已实现对象

- `CutCellKind::{Empty, Full, Cut, Unsupported}`
- `CutCellIssueCode`
- `CutCellIssue2D`
- `CutCell2D`
- `buildCutCell(AABB2D, CellClass, BoundaryLoop)`
- `buildCutCell(QuadtreeLeaf2D, BoundaryLoop)`

`CutCell2D` 实际保存：

- source leaf id/key
- background bounds
- real fluid `Polygon2D`
- area
- area fraction
- centroid
- embedded-boundary segment fragments
- explicit issues

## 核心算法

1. 先复用 2D-0 `BoundaryLoop::diagnose()` 拒绝非法、自交、退化边界。
2. 合法输入统一安全归一化为 CCW。
3. 对 `Intersected` leaf 使用四个 axis-aligned half-plane 顺序裁剪，构造 `boundary interior ∩ leaf AABB`。
4. 删除 tolerance 下连续重复顶点并统一输出 CCW polygon。
5. 计算真实 polygon area / centroid / area fraction。
6. 使用参数裁剪提取输入 boundary segment 位于 leaf 内的片段。
7. 位于背景 cell 外框上的片段不冒充 embedded boundary。
8. embedded fragments 按端点连通性审计；若一个 leaf 中出现多个不连通 embedded-boundary component，则显式 `Unsupported`。
9. 切线/角点仅接触导致零流体面积时输出 `Empty`，不制造零面积假 Cut-cell。
10. partial polygon 若缺失 embedded boundary、面积越界、退化或自交则显式失败。

## 解析基准

### 1. 斜线切单位格

单位 cell `[0,1] x [0,1]`，流体边界取三角形：

`(-1,-1) -> (2,-1) -> (-1,2)`

实际 cell 内流体区域为三角形：

`(0,0), (1,0), (0,1)`

解析值：

- area = `0.5`
- area fraction = `0.5`
- centroid = `(1/3, 1/3)`
- embedded boundary = `(1,0) -> (0,1)`
- embedded length = `sqrt(2)`

本地 C++20 解析 harness 已逐项对照通过。

### 2. 竖直切割 25%

边界在 cell 内为 `x = 0.25`，保留左侧流体。

解析值：

- area = `0.25`
- area fraction = `0.25`
- centroid = `(0.125, 0.5)`

本地解析 harness 已通过。

## 回归测试已写入

`cartmesh2d/tests/cutcell_test.cpp` 当前覆盖：

- diagonal analytic triangle area/centroid/alpha
- embedded fragment endpoint 位于原始 boundary
- embedded fragment 解析长度
- vertical 25% analytic cut
- cut polygon CCW orientation
- Inside -> real full rectangle polygon
- Outside -> true empty result
- tangent/corner zero-area contact -> Empty，不生成假 cut polygon
- clockwise valid input normalization
- bow-tie self-intersection rejection
- disconnected multi-component cell intersection explicit failure
- Quadtree leaf id/key 传播

## 多组件失败案例

测试包含一个单个合法 simple boundary（两个上部流体块通过 cell 外部走廊连接）。

该 boundary 与单位 cell 的交集在 cell 内分裂为两个流体 component。2D-3 第一版不声称支持 multi-component Cut-cell，因此必须返回：

- `CutCellKind::Unsupported`
- topology/geometry issue

严禁把两个 component 用伪造 bridge 连成一个 polygon。

## 当前门禁状态

已完成：

- [x] 真正 fluid polygon，而非整格删除/保留
- [x] analytic straight cut area
- [x] analytic centroid
- [x] `0 < alpha < 1` 对真实 Cut-cell 成立
- [x] polygon CCW
- [x] embedded boundary fragment
- [x] fragment 位于输入边界
- [x] tangent zero-area 不制造坏 cell
- [x] invalid boundary 显式失败
- [x] unsupported multi-component 显式失败
- [x] 本地 C++20 解析 harness 通过
- [ ] 当前 GitHub 完整 branch 的根 CMake build + 2D-0/1/2/3 CTest 最终门禁

## 明确未开始

- global vertex deduplication
- global Edge2D
- owner/neighbour
- boundary patch
- coarse-fine edge splitting
- small-cell stabilization/agglomeration
- solver export
- visualization

以上均属于 2D-4 或后续阶段。
