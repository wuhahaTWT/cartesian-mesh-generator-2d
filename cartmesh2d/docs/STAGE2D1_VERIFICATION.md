# Stage 2D-1 验证记录

## 状态

**PASS / CLOSED — 2D-1 已完成并通过根目录 CMake + CTest 集成门禁。**

当前严格停留在 Uniform Cartesian + classification；未实现 Quadtree、2:1 balance、Cut-cell、全局拓扑或可视化。

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

## 测试矩阵

已通过：

- 4 × 3 网格 cell 数、dx/dy、中心坐标
- row-major ID 稳定性
- cell 完整覆盖 domain
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

32 × 32 grid，单位圆用 64 条线段离散。
当前确定性结果：

- Inside = 164
- Intersected = 68
- Outside = 792
- Total = 1024

## 最终集成门禁

验证环境：GNU C++ 14.2.0，C++20。

从仓库根目录执行：

```sh
cmake -S . -B build-verify1 \
  -DCMAKE_BUILD_TYPE=Release \
  -DCARTMESH_BUILD_2D=ON \
  -DCARTMESH_BUILD_TESTS=ON \
  -DCARTMESH_BUILD_BENCHMARKS=OFF

cmake --build build-verify1 \
  --target cartmesh2d_geometry_tests cartmesh2d_grid_tests

ctest --test-dir build-verify1 \
  -R 'cartmesh2d_stage(0|1)' \
  --output-on-failure
```

结果：

```text
cartmesh2d_geometry_tests  built successfully
cartmesh2d_grid_tests      built successfully

1/2 cartmesh2d_stage0_geometry_tests ... Passed
2/2 cartmesh2d_stage1_grid_tests ....... Passed

100% tests passed, 0 tests failed out of 2
```

同时验证默认关闭二维时根配置仍正常：

```sh
cmake -S . -B build-verify1-off \
  -DCMAKE_BUILD_TYPE=Release \
  -DCARTMESH_BUILD_2D=OFF \
  -DCARTMESH_BUILD_TESTS=OFF \
  -DCARTMESH_BUILD_BENCHMARKS=OFF
```

结果：PASS。

## 2D-1 验收对照

- [x] `Nx*Ny` cell 数正确
- [x] cell 完整覆盖 domain
- [x] 相邻 cell 无 gap/overlap
- [x] rectangle fixture 分类统计可复核
- [x] circle fixture 分类统计可复核
- [x] tangent 有测试
- [x] boundary-on-grid-line 有测试
- [x] corner touch 有测试
- [x] 重复运行 cell ID 稳定
- [x] 重复运行 classification 稳定
- [x] 非法 boundary 不静默分类
- [x] 根目录 CMake build + 2D-0/2D-1 CTest 最终门禁
- [x] `CARTMESH_BUILD_2D=OFF` 根配置正常
- [x] 未修改三维算法目录

## 隔离检查

与 `main` 比较，2D-1 相关新增仍位于 `cartmesh2d/**`；根目录仅保留 2D 子项目开关/接入。未修改：

- `include/cartmesh/**`
- 根 `src/**`
- 根 `apps/**`
- 根 `tests/**`

## 尚未开始

- 2D-2 Quadtree
- 2:1 balance
- distance / curvature refinement
- 2D-3 Cut-cell
- 2D-4 topology
- 2D-5 small-cell
- 2D-6 export/quality
- 2D-V visualization

**结论：Stage 2D-1 关闭。未经用户明确批准，不得自动进入 2D-2。**
