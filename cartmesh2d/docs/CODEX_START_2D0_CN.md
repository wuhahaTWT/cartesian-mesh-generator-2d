# Codex 阶段状态：2D-0 ~ 2D-6 已关闭

## 当前状态

- Stage 2D-0：PASS / CLOSED
- Stage 2D-1：PASS / CLOSED
- Stage 2D-2：PASS / CLOSED
- Stage 2D-3：PASS / CLOSED
- Stage 2D-4：PASS / CLOSED
- Stage 2D-5：PASS / CLOSED
- Stage 2D-6：PASS / CLOSED
- Stage 2D-V：NOT STARTED

开始或继续任何二维修改前必须阅读：

1. 根目录 `AGENTS.md` 中二维并行子项目例外；
2. `cartmesh2d/AGENTS.md`；
3. `cartmesh2d/docs/PROJECT_BRIEF_CN.md`；
4. `cartmesh2d/docs/ARCHITECTURE_CN.md`；
5. `cartmesh2d/docs/STAGE_PLAN_CN.md`；
6. `cartmesh2d/docs/ACCEPTANCE_CN.md`；
7. `cartmesh2d/docs/STAGE2D0_VERIFICATION.md` ~ `STAGE2D6_VERIFICATION.md`。

## 已关闭能力

二维核心产品链已完成：

```text
BoundaryLoop
-> Cartesian domain
-> Quadtree adaptive refinement
-> 2:1 balance
-> true Cut-cell polygon
-> global owner/neighbour topology
-> small-cell detection
-> topology-safe agglomeration
-> quality report
-> VTK / CM2D / JSON export
-> end-to-end CLI
-> independent CM2D read-back
```

## Stage 2D-6 最终验收

GitHub Actions `cartmesh2d-stage6` run #3 对 implementation head `15b383bee6292476e8348e28e1b0d7b9ce4d46ea` 完成真实根工程验证：

- root `CARTMESH_BUILD_2D=ON` configure PASS
- full native-2D build PASS
- Stage 0~6 CTest：15/15 PASS
- rectangle E2E PASS
- circle E2E PASS
- concave E2E PASS
- airfoil-like E2E PASS
- VTK / CM2D / JSON artifacts present and validated
- JSON `valid=true`
- all topology audit counters = 0
- `min_cell_area > 0`
- `min_edge_length > 0`
- root `CARTMESH_BUILD_2D=OFF` configure PASS

详见 `cartmesh2d/docs/STAGE2D6_VERIFICATION.md`。

## 下一允许阶段

下一阶段仅为 **2D-V visualization**。

2D-V 必须保持薄层：读取已经导出的网格/报告来展示 cell edges、Cut-cells、Quadtree level、boundary 与 quality/small-cell flags；不得复制或重写核心 meshing 算法。

三维核心目录仍不得为二维功能修改。