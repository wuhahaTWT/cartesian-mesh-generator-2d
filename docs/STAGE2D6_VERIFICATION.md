# Stage 2D-6 验证记录

## 状态

**PASS / CLOSED — 2D-6A quality、2D-6B export/read-back、2D-6C end-to-end CLI 与四类 acceptance fixture 已在 GitHub Actions 对当前实现 head 完成真实根工程编译与验收。**

2D-V visualization 仍未开始。

## 2D-6A — Quality

已实现：

- `include/cartmesh2d/quality/Quality2D.hpp`
- `src/quality/Quality2D.cpp`
- `tests/quality_test.cpp`

质量报告包含：

- vertex / edge / cell count
- internal / boundary edge count
- source Cut-cell / Full-cell / small-cell count
- `minCellArea`
- `minEdgeLength`
- `maxEdgeAspectRatio = max(cell max-edge / min-edge)`
- `maxCentroidSkewness = |area-centroid - vertex-mean| / sqrt(area)`
- `minCutCellAreaFraction`
- Quadtree level distribution
- Stage 2D-4 topology audit 全部计数
- deterministic JSON report

## 2D-6B — Export / Read-back

已实现：

- Legacy VTK `UNSTRUCTURED_GRID`
- deterministic `CM2D v1` solver-topology format
- independent `readCm2dTopology(...)`

CM2D 保存 vertices、edges、owner、neighbour、boundary patch、cell vertex/edge loops、source ids/keys 和 topology audit。

门禁：

- invalid topology 拒绝导出
- 相同 topology 两次 CM2D 输出字节级一致
- read-back 恢复 vertex/edge/cell 数量
- owner/neighbour 与 cell loops 精确一致
- corrupt / unsupported CM2D 显式失败

## 2D-6C — End-to-end CLI

CLI：

```text
cartmesh2d_cli <boundary.xy> <output-prefix> [max-level=5] [padding-fraction=0.25] [small-alpha=0.10]
```

完整流水线：

```text
boundary file
-> BoundaryLoop diagnostics / CCW normalization
-> padded Cartesian domain
-> Quadtree boundary refinement
-> 2:1 balance
-> true CutCell2D
-> global topology
-> small-cell analysis
-> topology-safe agglomeration
-> quality report
-> VTK + CM2D + JSON export
-> independent CM2D read-back
-> exact topology comparison
```

任何 unsupported Cut-cell、unresolved small cell、agglomeration failure、topology audit failure、quality failure 或 read-back mismatch 均非零退出。

## Final acceptance fixtures

- `rectangle.xy`
- `circle.xy`
- `concave.xy`
- `airfoil_like.xy`

四条均作为真正 CLI end-to-end CTest 执行。

## Exact-head GitHub Actions 验收

Workflow：`.github/workflows/cartmesh2d-stage6.yml`

成功 run：

```text
workflow: cartmesh2d-stage6
run_number: 3
run_id: 32330927546
conclusion: success
validated implementation head: 15b383bee6292476e8348e28e1b0d7b9ce4d46ea
PR merge checkout: e855bfa480f6eddd604ca91f9a7a3b9adaf9bdf1
runner: Ubuntu 24.04 / GCC 13.3.0
```

### Root configure / build

`CARTMESH_BUILD_2D=ON` 根工程 configure 成功；Stage 0~6 全部二维测试 target 与 `cartmesh2d_cli` 编译成功。

仅有一个非致命 warning：`Topology2D.cpp` 中 `findVertexId(...)` 当前未使用；不影响验收结果，后续可清理。

### Full CTest gate

```text
cartmesh2d_stage0_geometry_tests .................. Passed
cartmesh2d_stage1_grid_tests ...................... Passed
cartmesh2d_stage2_quadtree_tests .................. Passed
cartmesh2d_stage3_cutcell_tests ................... Passed
cartmesh2d_stage3_quadtree_cutcell_audit_tests .... Passed
cartmesh2d_stage4_topology_tests .................. Passed
cartmesh2d_stage4_adaptive_topology_audit_tests ... Passed
cartmesh2d_stage5a_small_cell_tests ............... Passed
cartmesh2d_stage5b_agglomeration_tests ............ Passed
cartmesh2d_stage6a_quality_tests .................. Passed
cartmesh2d_stage6b_io_tests ....................... Passed
cartmesh2d_stage6c_rectangle_e2e .................. Passed
cartmesh2d_stage6c_circle_e2e ..................... Passed
cartmesh2d_stage6c_concave_e2e .................... Passed
cartmesh2d_stage6c_airfoil_like_e2e ............... Passed

100% tests passed, 0 tests failed out of 15
```

### Exported artifact gate

对 rectangle / circle / concave / airfoil_like 四组输出，CI 均要求并通过：

- `.vtk` 存在且非空
- `.cm2d` 存在且非空
- `.quality.json` 存在且非空
- CM2D header 为 `CM2D 1`
- JSON `valid == true`
- 所有 `topology_audit` 计数均为 0
- vertices / edges / cells 均 > 0
- `min_cell_area > 0`
- `min_edge_length > 0`

### 2D disabled gate

根工程再次以 `CARTMESH_BUILD_2D=OFF` configure，成功。

## 分支隔离最终审计

当前 branch 相对 `main`：

- `behind_by = 0`
- merge base 仍为 main baseline `8bec26d98eb8bd84033625ed2a41184c8cb223f1`
- 改动仅为 `.github/workflows/cartmesh2d-stage6.yml`、根二维入口/说明和 `cartmesh2d/**`
- 未修改三维算法目录 `include/cartmesh/**`、root `src/**`、root `apps/**`、root `tests/**`

## Stage 2D-6 closure gate

- [x] quality evaluator
- [x] deterministic JSON
- [x] VTK export
- [x] solver-topology CM2D export
- [x] CM2D independent read-back
- [x] end-to-end CLI
- [x] rectangle fixture PASS
- [x] circle fixture PASS
- [x] concave fixture PASS
- [x] airfoil-like fixture PASS
- [x] exact-head root CMake configure PASS
- [x] full 2D build PASS
- [x] Stage 0~6 CTest 15/15 PASS
- [x] exported artifact validation PASS
- [x] `CARTMESH_BUILD_2D=OFF` PASS
- [x] branch isolation PASS

因此 Stage 2D-6 正式状态：**PASS / CLOSED**。

二维核心产品链现已完成 2D-0 ~ 2D-6；下一阶段仅为 2D-V visualization。