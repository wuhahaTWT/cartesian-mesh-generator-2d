# Codex 阶段状态：2D-0 ~ 2D-4 已关闭

## 当前状态

- Stage 2D-0：PASS / CLOSED
- Stage 2D-1：PASS / CLOSED
- Stage 2D-2：PASS / CLOSED
- Stage 2D-3：PASS / CLOSED
- Stage 2D-4：PASS / CLOSED
- Stage 2D-5：NOT STARTED

开始或继续任何二维修改前必须阅读：

1. 根目录 `AGENTS.md` 中二维并行子项目例外；
2. `cartmesh2d/AGENTS.md`；
3. `cartmesh2d/docs/PROJECT_BRIEF_CN.md`；
4. `cartmesh2d/docs/ARCHITECTURE_CN.md`；
5. `cartmesh2d/docs/STAGE_PLAN_CN.md`；
6. `cartmesh2d/docs/ACCEPTANCE_CN.md`；
7. `cartmesh2d/docs/STAGE2D0_VERIFICATION.md`；
8. `cartmesh2d/docs/STAGE2D1_VERIFICATION.md`；
9. `cartmesh2d/docs/STAGE2D2_VERIFICATION.md`；
10. `cartmesh2d/docs/STAGE2D3_VERIFICATION.md`；
11. `cartmesh2d/docs/STAGE2D4_VERIFICATION.md`。

## 已关闭能力

### 2D-0

原生二维几何内核、鲁棒线段相交、point-in-polygon 和 `BoundaryLoop` 诊断。

### 2D-1

均匀 Cartesian 网格、确定性 cell IDs、`Outside / Inside / Intersected` 分类及 tangent/grid-line/corner-touch 明确规则。

### 2D-2

原生 Quadtree 1->4 refinement、boundary/distance refinement、deterministic leaf key/ID、face-neighbor discovery 和 2:1 balance。

### 2D-3

真实 `CutCell2D` fluid polygon、area / centroid / area fraction、embedded boundary fragment、解析切割基准和病态 multi-component 显式拒绝。

### 2D-4

全局 `Vertex2D / Edge2D / TopologyCell2D`、owner/neighbour、boundary patch、coarse-fine hanging-node edge splitting、deterministic topology IDs 和完整 topology audit。

最终 2D-4 根目录门禁：

```text
2D-0 / 2D-1 / 2D-2 / 2D-3 / 2D-4
100% tests passed, 0 tests failed out of 7
```

并且 adaptive-circle topology 独立审计满足：

- duplicate vertex = 0
- duplicate edge = 0
- orphan internal edge = 0
- non-manifold edge = 0
- unclassified boundary edge = 0
- open cell loop = 0
- area mismatch = 0

## 停线要求

未经用户明确批准，不得自动开始 Stage 2D-5。

尤其禁止提前实现：

- small-cell detection / agglomeration
- quality/export
- GUI / visualization

三维核心目录仍不得为二维功能修改。
