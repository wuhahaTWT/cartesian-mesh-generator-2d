# Codex 阶段状态：2D-0 ~ 2D-V 全部关闭

## 当前状态

- Stage 2D-0：PASS / CLOSED
- Stage 2D-1：PASS / CLOSED
- Stage 2D-2：PASS / CLOSED
- Stage 2D-3：PASS / CLOSED
- Stage 2D-4：PASS / CLOSED
- Stage 2D-5：PASS / CLOSED
- Stage 2D-6：PASS / CLOSED
- Stage 2D-V：PASS / CLOSED

开始或修改任何二维功能前必须阅读：

1. 根目录 `AGENTS.md` 中二维并行子项目例外；
2. `cartmesh2d/AGENTS.md`；
3. `cartmesh2d/docs/PROJECT_BRIEF_CN.md`；
4. `cartmesh2d/docs/ARCHITECTURE_CN.md`；
5. `cartmesh2d/docs/STAGE_PLAN_CN.md`；
6. `cartmesh2d/docs/ACCEPTANCE_CN.md`；
7. `cartmesh2d/docs/STAGE2D0_VERIFICATION.md` ~ `STAGE2D6_VERIFICATION.md`；
8. `cartmesh2d/docs/STAGE2DV_VERIFICATION.md`。

## 已关闭产品链

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
-> thin-layer visualization
```

## Stage 2D-V 最终状态

2D-V 为薄层后处理：

```text
CM2D + quality.json + viz.json
-> dependency-free Python renderer
-> standalone SVG
```

支持：

- final solver cell polygons / edges；
- adaptive level / area coloring；
- embedded / domain / unclassified boundary；
- Cut/boundary cell identification；
- source background-cell bounds；
- source small-cell exact overlay（pre-agglomeration）；
- cell id labels；
- quality/topology audit panel；
- invalid topology/quality visible banner。

CLI 输出 `<prefix>.viz.json`，只保存核心 pipeline 已计算的展示元数据，不产生新的分类或 meshing 判断。

最终 GitHub Actions run #16 (`32335057216`) 对 exact implementation head `dd2546a1ac2a95ea5a7baf07fbebcfb1a9a9c287` 验证成功：

- root 2D-ON configure PASS；
- full native-2D build PASS；
- Stage 0~6：15/15 CTest PASS；
- four E2E mesh exports PASS；
- viz sidecar validation PASS；
- renderer regression PASS；
- rectangle/circle/concave/airfoil-like 四个真实 SVG PASS；
- visualization artifact upload PASS；
- root 2D-OFF configure PASS。

详见 `cartmesh2d/docs/STAGE2DV_VERIFICATION.md`。

## 后续修改规则

二维 2D-0 ~ 2D-V 已全部封口。后续如需新增 GUI、更多输入格式、求解器接口或新网格能力，应新开阶段/任务，不应回写已关闭阶段的验收定义来伪装为原阶段工作。

三维核心目录仍不得因二维展示功能修改。
