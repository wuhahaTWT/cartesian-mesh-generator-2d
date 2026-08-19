# Codex 下一步执行单：2D-0 已启动

> 本文件原用于二维架构基线确认后启动 Stage 2D-0。2026-08-19 用户已明确批准开始编码，2D-0 几何内核现已实现并进入验证状态。

## 当前阶段

只允许推进 `cartmesh2d` 的 **Stage 2D-0：二维几何内核**。

开始或继续修改前必须依次阅读：

1. 根目录 `AGENTS.md` 中“二维并行子项目例外”；
2. `cartmesh2d/AGENTS.md`；
3. `cartmesh2d/docs/PROJECT_BRIEF_CN.md`；
4. `cartmesh2d/docs/ARCHITECTURE_CN.md`；
5. `cartmesh2d/docs/STAGE_PLAN_CN.md`；
6. `cartmesh2d/docs/ACCEPTANCE_CN.md`；
7. `cartmesh2d/docs/STAGE2D0_VERIFICATION.md`。

## 修改范围

允许：

- `cartmesh2d/**`
- 顶层 CMake 中仅与 `cartmesh2d` 接入直接相关的必要小改动

禁止：

- 修改 `include/cartmesh/**`
- 修改根 `src/**`
- 修改根 `apps/**`
- 修改根 `tests/**`
- 修改三维 Stage 6/7 实现
- 做 GUI/绘图
- 实现 Quadtree、Cut-cell 或后续阶段

## 2D-0 当前已实现

- `Point2D`, `Vector2D`
- `Segment2D`
- `AABB2D`
- `Polygon2D`
- `BoundaryLoop`
- tolerance policy/config
- cross/dot/orientation
- signed area / area
- polygon centroid
- point-on-segment
- segment intersection：`None / Point / Overlap`
- point-in-polygon：`Inside / Outside / Boundary`
- BoundaryLoop diagnostics：
  - too few unique vertices
  - duplicate consecutive vertex
  - zero-length edge
  - self-intersection
  - CW/CCW orientation
  - safe normalization to project-standard CCW

## 当前测试

已建立并通过：

- rectangle
- triangle
- concave L polygon
- discretized circle
- crossing segments
- endpoint touching segments
- parallel non-intersecting segments
- collinear overlapping segments
- bow-tie self-intersecting loop（失败案例）
- duplicate consecutive point（失败案例）

测试验证数值结果，不依赖图片。

## 构建

二维子项目已成为真实可编译 target 和 CTest target，可独立执行：

```sh
cmake -S cartmesh2d -B cartmesh2d/build -DCMAKE_BUILD_TYPE=Release
cmake --build cartmesh2d/build
ctest --test-dir cartmesh2d/build --output-on-failure
```

## 停线要求

在用户再次明确批准前，不得自动开始 2D-1。后续若只做顶层 `CARTMESH_BUILD_2D` 构建接入，应视为 2D-0 的工程集成收尾，不得借机加入 Cartesian grid、Quadtree、Cut-cell 或可视化功能。
