# Codex 阶段状态：2D-0 ~ 2D-5 已关闭，2D-6 进行中

## 当前状态

- Stage 2D-0：PASS / CLOSED
- Stage 2D-1：PASS / CLOSED
- Stage 2D-2：PASS / CLOSED
- Stage 2D-3：PASS / CLOSED
- Stage 2D-4：PASS / CLOSED
- Stage 2D-5：PASS / CLOSED
- Stage 2D-6：IN PROGRESS
- Stage 2D-V：NOT STARTED

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
12. `cartmesh2d/docs/STAGE2D5_VERIFICATION.md`；
13. `cartmesh2d/docs/STAGE2D6_VERIFICATION.md`。

## 已关闭能力

2D-0 ~ 2D-5 已完成：原生二维 geometry、Cartesian/grid classification、Quadtree + 2:1、true Cut-cell polygon、global owner/neighbour topology、small-cell detection 与 topology-safe agglomeration。

Stage 2D-5 最近一次实际根工程回归为 9/9 PASS；验证-5 时确认测试通过的 implementation head 未漂移后正式 CLOSED。

## 2D-6 当前实现

### 2D-6A Quality

已加入：

- geometry/topology quality evaluator；
- min area / min edge；
- max edge aspect ratio；
- centroid skewness；
- min Cut-cell alpha；
- level distribution；
- topology audit aggregation；
- deterministic JSON report。

### 2D-6B Export / Read-back

已加入：

- Legacy VTK polygon export；
- deterministic `CM2D v1` solver-topology format；
- CM2D independent reader；
- byte-determinism / owner-neighbour / cell-loop read-back tests。

### 2D-6C CLI / Acceptance

已加入 `cartmesh2d_cli`：

```text
boundary.xy
-> Quadtree
-> 2:1
-> Cut-cell
-> topology
-> small-cell stabilization
-> quality
-> VTK/CM2D/JSON
-> independent CM2D read-back
```

已注册四类 end-to-end fixture：

- rectangle
- circle
- concave
- airfoil-like

## 当前停线

2D-6 尚未 CLOSED。当前实现必须先完成 exact-head 根 CMake 编译与所有 Stage 0~6 CTest；四类 CLI fixture 必须实际 PASS，导出文件必须独立核对。

在 2D-6 正式关闭前禁止开始：

- Stage 2D-V visualization；
- GUI；
- 为展示效果复制/重写核心 meshing 算法。

三维核心目录仍不得为二维功能修改。
