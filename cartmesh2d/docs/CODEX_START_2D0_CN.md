# Codex 下一步执行单：只启动 2D-0

> 本文件用于二维架构基线确认后，用户明确要求开始编码时交给 Codex。

## 任务

只实现 `cartmesh2d` 的 **Stage 2D-0：二维几何内核**。

开始前必须依次阅读：

1. 根目录 `AGENTS.md` 中“二维并行子项目例外”；
2. `cartmesh2d/AGENTS.md`；
3. `cartmesh2d/docs/PROJECT_BRIEF_CN.md`；
4. `cartmesh2d/docs/ARCHITECTURE_CN.md`；
5. `cartmesh2d/docs/STAGE_PLAN_CN.md`；
6. `cartmesh2d/docs/ACCEPTANCE_CN.md`。

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

## 2D-0 必须实现

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
- segment intersection，明确区分 `None / Point / Overlap`
- point-in-polygon，明确区分 `Inside / Outside / Boundary`
- BoundaryLoop diagnostics：
  - too few unique vertices
  - duplicate consecutive vertex
  - zero-length edge
  - self-intersection
  - CW/CCW orientation
  - normalization to project-standard CCW when safe

## 测试要求

至少建立：

- rectangle
- triangle
- concave L polygon
- discretized circle
- crossing segments
- endpoint touching segments
- parallel non-intersecting segments
- collinear overlapping segments
- bow-tie self-intersecting loop（必须失败）
- duplicate consecutive point（必须失败或给明确 invalid diagnostic）

测试必须验证数值结果，而不是依赖图片。

## 构建要求

把 2D-0 变成真实可编译 target 和测试 target；保持 `CARTMESH_BUILD_2D=OFF` 时现有三维默认构建行为不变。

## 完成时必须提交的结果

1. 实际修改文件列表；
2. 2D-0 测试命令和结果；
3. 对照 `ACCEPTANCE_CN.md` 的逐项 PASS/FAIL；
4. 仍未实现的内容；
5. 明确停止在 2D-0，不得自动开始 2D-1。
