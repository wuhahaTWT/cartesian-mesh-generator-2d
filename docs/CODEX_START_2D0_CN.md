# Codex 阶段状态：2D-0 / 2D-1 已关闭

## 当前状态

- Stage 2D-0：PASS / CLOSED
- Stage 2D-1：PASS / CLOSED
- Stage 2D-2：尚未开始

开始或继续任何二维修改前必须阅读：

1. 根目录 `AGENTS.md` 中二维并行子项目例外；
2. `cartmesh2d/AGENTS.md`；
3. `cartmesh2d/docs/PROJECT_BRIEF_CN.md`；
4. `cartmesh2d/docs/ARCHITECTURE_CN.md`；
5. `cartmesh2d/docs/STAGE_PLAN_CN.md`；
6. `cartmesh2d/docs/ACCEPTANCE_CN.md`；
7. `cartmesh2d/docs/STAGE2D0_VERIFICATION.md`；
8. `cartmesh2d/docs/STAGE2D1_VERIFICATION.md`。

## 已关闭能力

2D-0 已提供原生二维几何内核、鲁棒线段相交、point-in-polygon 和 BoundaryLoop 诊断。

2D-1 已提供：

- `Domain2D`
- `CartesianCell2D`
- `UniformCartesianGrid2D`
- Nx × Ny / spacing-based 均匀 Cartesian 网格
- row-major 确定性 cell IDs
- `Outside / Inside / Intersected` 分类
- closed-AABB boundary intersection
- tangent / grid-line / corner-touch 明确规则
- invalid boundary 分类前拒绝
- rectangle / circle 可复核统计

2D-0 与 2D-1 已从仓库根目录通过 CMake + CTest 集成门禁。

## 停线要求

未经用户明确批准，不得自动开始 Stage 2D-2。

尤其禁止提前实现：

- Quadtree refinement
- 2:1 balance
- distance / curvature refinement
- Cut-cell polygon
- 全局 topology
- GUI / visualization

三维核心目录仍不得为二维功能修改。
