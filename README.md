# cartmesh2d — 原生二维自适应 Cartesian / Cut-cell 网格生成器

`cartmesh2d` 是主仓库中的**独立、封闭、原生二维**子项目。

它不是三维 `cartmesh` 的 `z=0` 模式，也不通过把 `Point3D`、Octree 或 3D Cut-cell 代码模板化来复用三维核心。二维项目首先追求一个小而完整、可验证、可作为结题成果的二维 CFD 网格生成核心；三维项目继续独立推进。

## 核心流水线

```text
2D closed boundary
    -> geometry validation
    -> Cartesian background grid
    -> boundary/cell intersection classification
    -> adaptive Quadtree refinement
    -> 2:1 balance
    -> inside / outside / intersected classification
    -> Cut-cell polygon construction
    -> small-cell handling
    -> cell-edge-neighbor topology
    -> mesh quality validation
    -> solver/standard export
    -> visualization (last)
```

## 阶段

- **2D-0**：二维几何内核
- **2D-1**：均匀 Cartesian 网格与分类
- **2D-2**：Quadtree 自适应与 2:1 平衡
- **2D-3**：真实二维 Cut-cell polygon
- **2D-4**：完整 cell-edge-neighbor 拓扑
- **2D-5**：小 Cut-cell 检测与稳定化/聚合
- **2D-6**：质量、导出、最终验收
- **2D-V**：可视化；只有 2D-6 核心验收后才进入

## 开发前必读

1. `AGENTS.md`
2. `docs/PROJECT_BRIEF_CN.md`
3. `docs/ARCHITECTURE_CN.md`
4. `docs/STAGE_PLAN_CN.md`
5. `docs/ACCEPTANCE_CN.md`

当前基线只建立项目边界、架构与验收标准。除非用户明确要求开始某个阶段，否则不得跨阶段实现算法。
