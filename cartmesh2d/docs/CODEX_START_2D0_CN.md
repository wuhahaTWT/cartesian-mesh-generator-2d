# Codex 阶段状态：2D-0 ~ 2D-4 已关闭，2D-5 等待验证

## 当前状态

- Stage 2D-0：PASS / CLOSED
- Stage 2D-1：PASS / CLOSED
- Stage 2D-2：PASS / CLOSED
- Stage 2D-3：PASS / CLOSED
- Stage 2D-4：PASS / CLOSED
- Stage 2D-5A：PASS
- Stage 2D-5B：PASS
- Stage 2D-5：READY FOR VALIDATION
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
均匀 Cartesian 网格、确定性 cell IDs、`Outside / Inside / Intersected` 分类。

### 2D-2
原生 Quadtree 1->4 refinement、deterministic leaf key/ID、face-neighbor discovery 和 2:1 balance。

### 2D-3
真实 `CutCell2D` fluid polygon、area / centroid / area fraction、embedded boundary fragment 和病态输入显式拒绝。

### 2D-4
全局 `Vertex2D / Edge2D / TopologyCell2D`、owner/neighbour、boundary patch、coarse-fine hanging-node edge splitting 和完整 topology audit。

## 2D-5A 已通过

`SmallCell2D` 已提供：

- configurable alpha threshold；
- Cut-cell alpha histogram；
- deterministic small-cell marking；
- internal-edge neighbour discovery；
- stable-neighbour preference；
- deterministic best candidate；
- unresolved explicit failure。

shifted-circle 自适应 fixture 在 threshold=0.1 下：small=8，minimum alpha 约 0.00120311，unresolved=0。

## 2D-5B 已通过当前实现门禁

`Agglomeration2D` 已提供：

- `AgglomeratedCell2D`；
- topology-safe group merge；
- 组内共享 edge fragment 消除；
- 单闭环 exterior reconstruction；
- 共线冗余顶点简化；
- member-area / merged-area 一致性检查；
- total fluid-area conservation；
- post-agglomeration global topology rebuild；
- Stage 2D-4 topology audit 再验证；
- unsafe small->small、断链、多环、分叉、退化 polygon 显式失败。

当前根目录二维回归：

```text
2D-0 / 2D-1 / 2D-2 / 2D-3 / 2D-4 / 2D-5A / 2D-5B
100% tests passed, 0 tests failed out of 9
```

真实 shifted-circle 聚合：

- detected small cells = 8
- merged small cells = 8
- output cells = input cells - 8
- total area error <= 1e-10
- duplicate/orphan/non-manifold/unclassified/open-loop/area-mismatch = 0

## 当前停线

Stage 2D-5 尚未正式 CLOSED。

下一允许动作仅为用户显式 **`验证-5`** 后的 current-head 封口复跑与状态关闭。

在 `验证-5` 通过前不得提前实现：

- 2D-6 quality/export；
- GUI / visualization。

三维核心目录仍不得为二维功能修改。
