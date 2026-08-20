# Codex 阶段状态：2D-0 ~ 2D-6 已关闭，2D-V 待最终验证

## 当前状态

- Stage 2D-0：PASS / CLOSED
- Stage 2D-1：PASS / CLOSED
- Stage 2D-2：PASS / CLOSED
- Stage 2D-3：PASS / CLOSED
- Stage 2D-4：PASS / CLOSED
- Stage 2D-5：PASS / CLOSED
- Stage 2D-6：PASS / CLOSED
- Stage 2D-V：READY FOR VALIDATION

开始或继续任何二维修改前必须阅读：

1. 根目录 `AGENTS.md` 中二维并行子项目例外；
2. `cartmesh2d/AGENTS.md`；
3. `cartmesh2d/docs/PROJECT_BRIEF_CN.md`；
4. `cartmesh2d/docs/ARCHITECTURE_CN.md`；
5. `cartmesh2d/docs/STAGE_PLAN_CN.md`；
6. `cartmesh2d/docs/ACCEPTANCE_CN.md`；
7. `cartmesh2d/docs/STAGE2D0_VERIFICATION.md` ~ `STAGE2D6_VERIFICATION.md`；
8. `cartmesh2d/docs/STAGE2DV_VERIFICATION.md`。

## 已关闭核心能力

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

## Stage 2D-V 当前实现

2D-V 保持薄层，不复制核心算法：

```text
CM2D + quality.json + viz.json
-> dependency-free Python renderer
-> standalone SVG
```

当前支持：

- final solver cell polygons / edges；
- adaptive level / area coloring；
- embedded / domain / unclassified boundary；
- Cut/boundary cell identification；
- source background-cell bounds；
- source small-cell exact overlay（pre-agglomeration）；
- cell id labels；
- quality/topology audit panel；
- invalid topology/quality visible banner。

CLI 额外输出 `<prefix>.viz.json`，仅保存已经由核心 pipeline 算出的 source-cell 展示元数据，不产生新的分类或 meshing 判断。

GitHub Actions run #14 (`32334947692`) 已实际通过：完整 C++ build、原 Stage 0~6 15/15 CTest、viz sidecar 校验、renderer regression、rectangle/circle/concave/airfoil-like 四个真实 SVG、artifact upload、2D-OFF configure 全部成功。

详见 `cartmesh2d/docs/STAGE2DV_VERIFICATION.md`。

## 当前停线

等待用户显式 `验证-v` 后对当前最终 head 做封口并决定 PASS / CLOSED。

在此之前：

- 不把 2D-V 宣称 CLOSED；
- 不新增 GUI；
- 不为展示效果修改核心 meshing 算法；
- 三维核心目录仍不得为二维功能修改。
