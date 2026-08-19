# Codex 阶段状态：2D-0 ~ 2D-4 已关闭，2D-5 进行中

## 当前状态

- Stage 2D-0：PASS / CLOSED
- Stage 2D-1：PASS / CLOSED
- Stage 2D-2：PASS / CLOSED
- Stage 2D-3：PASS / CLOSED
- Stage 2D-4：PASS / CLOSED
- Stage 2D-5A：PASS（检测 / alpha histogram / best-neighbour candidate）
- Stage 2D-5B：NOT STARTED（topology-safe agglomeration）
- Stage 2D-5：IN PROGRESS
- Stage 2D-6：NOT STARTED

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
11. `cartmesh2d/docs/STAGE2D4_VERIFICATION.md`；
12. `cartmesh2d/docs/STAGE2D5_VERIFICATION.md`。

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

## 2D-5A 当前能力

新增 `SmallCell2D` 分析模块：

- configurable area-fraction threshold；
- Cut-cell alpha histogram；
- deterministic small-cell marking；
- 通过 Stage 2D-4 internal edges 寻找邻居；
- 按 stable target、shared interface length、target alpha、target area、topology id 确定最佳 candidate；
- 无候选 tiny cell 显式 `Unresolved`；
- shifted-circle 自适应真实 fixture 在 threshold=0.1 下稳定检测 8 个 small cells，minimum alpha 约 0.00120311，unresolved=0。

当前根目录二维回归：

```text
2D-0 / 2D-1 / 2D-2 / 2D-3 / 2D-4 / 2D-5A
100% tests passed, 0 tests failed out of 8
```

## 当前停线

Stage 2D-5 尚未 CLOSED。

下一允许工作仅为 **2D-5B topology-safe agglomeration**：

- 基于 2D-5A candidate graph 聚合；
- 聚合前后流体面积守恒；
- 重建合法 polygon / topology；
- Stage 2D-4 audit 重新 PASS；
- 不能产生 duplicate/non-manifold/open-loop；
- 无法安全聚合必须显式保留失败状态。

不得提前实现：

- 2D-6 quality/export；
- GUI / visualization。

三维核心目录仍不得为二维功能修改。
