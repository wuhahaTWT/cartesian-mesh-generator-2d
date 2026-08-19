# Stage 2D-1 验证记录

## 状态

**IN PROGRESS — 核心实现已完成并推送，最终完整分支构建/CTest 门禁尚未封口。**

当前仍严格停留在 Uniform Cartesian + classification；未实现 Quadtree、2:1 balance、Cut-cell、全局拓扑或可视化。

## 已实现

- `Domain2D`
- `CartesianCell2D`
- `UniformCartesianGrid2D`
- `Nx × Ny` 构造
- spacing-based 构造（按 `ceil(domain / target spacing)`）
- row-major 确定性 cell ID：`id = j * nx + i`
- `CellClass::Outside / Inside / Intersected`
- `segmentIntersectsClosedAABB`
- `classifyCartesianCell`
- `classifyGrid`
- `ClassificationSummary`

## 分类规则

1. cell 使用闭 AABB。
2. 只要任一边界 segment 与 cell 闭 AABB 有非空交集，该 cell 为 `Intersected`。
3. 因此 tangent / corner touch / boundary-on-grid-line 都算 `Intersected`。
4. 只有在没有边界 segment 与 cell 相交时，才允许使用 cell center 的 point-in-polygon 结果判定完整 `Inside / Outside`。
5. `classifyGrid` 开始前必须通过 2D-0 `BoundaryLoop::diagnose()`；自交或退化边界直接拒绝。

## 当前测试矩阵

已加入：

- 4 × 3 网格 cell 数、dx/dy、中心坐标
- row-major ID 稳定性
- cell 总面积覆盖 domain
- 横/纵相邻 cell 无 gap/overlap
- spacing-based grid
- segment 穿过 AABB
- segment 沿 AABB edge
- segment 与 AABB 分离
- segment 只切触 AABB corner
- segment 完全位于 AABB 内
- rectangle 边界恰落在 grid line
- rectangle inside/outside/intersected 精确统计
- 重复运行 ID 与 classification 完全一致
- diagonal / tangent boundary
- 64-segment circle fixture 精确统计
- zero grid dimension 拒绝
- bow-tie self-intersecting boundary 在分类前拒绝

## 可复核统计

### Rectangle

Domain `[0,4] × [0,4]`，8 × 8 grid，rectangle `[1,3] × [1,3]`。
闭 AABB 相交规则下：

- Inside = 4
- Intersected = 32
- Outside = 28
- Total = 64

### 64-segment circle

Domain `[-2,2] × [-2,2]`，32 × 32 grid，单位圆用 64 条线段离散。
当前确定性结果：

- Inside = 164
- Intersected = 68
- Outside = 792
- Total = 1024

## 本地实现验证

新 2D-1 源码已在 C++20 / CMake 环境中按 2D-0 公共 API 编译并执行 `grid_test`，结果：

```text
cartmesh2d 2D-1 grid/classification tests passed
```

注意：本记录此时不把阶段声明为 CLOSED，因为仍需在 GitHub 当前完整分支 checkout 上执行一次完整 CMake/CTest 门禁，确认 2D-0 + 2D-1 与根 CMake 集成共同通过。

## 2D-1 验收对照（当前）

- [x] `Nx*Ny` cell 数正确
- [x] cell 总面积完整覆盖 domain
- [x] 相邻 cell 无 gap/overlap
- [x] rectangle fixture 分类统计可复核
- [x] circle fixture 分类统计可复核
- [x] tangent 有测试
- [x] boundary-on-grid-line 有测试
- [x] corner touch 有测试
- [x] 重复运行 cell ID 稳定
- [x] 重复运行 classification 稳定
- [x] 非法 boundary 不静默分类
- [ ] GitHub 当前完整分支 CMake build + CTest 最终门禁

## 尚未开始

- 2D-2 Quadtree
- 2:1 balance
- distance / curvature refinement
- 2D-3 Cut-cell
- 2D-4 topology
- 2D-5 small-cell
- 2D-6 export/quality
- 2D-V visualization
