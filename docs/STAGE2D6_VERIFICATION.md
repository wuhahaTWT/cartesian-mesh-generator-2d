# Stage 2D-6 验证记录

## 状态

**IN PROGRESS — 2D-6A quality report、2D-6B export/read-back、2D-6C end-to-end CLI 与四类 acceptance fixture 已实现到分支；尚未对当前最终 head 完成根工程编译/CTest，因此不能标记 PASS / CLOSED。**

2D-V visualization 仍未开始。

## 2D-6A — Quality

新增：

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

新增：

- `include/cartmesh2d/io/MeshIO2D.hpp`
- `src/io/MeshIO2D.cpp`
- `tests/io_test.cpp`

输出：

1. Legacy VTK `UNSTRUCTURED_GRID`：polygon cell 几何、geometry area、source id；
2. `CM2D v1`：确定性 solver-topology 文本格式，包含 vertices、edges、owner、neighbour、boundary patch、cell vertex/edge loops、source ids/keys、topology audit；
3. 独立 `readCm2dTopology(...)` parser。

门禁设计：

- invalid topology 拒绝导出；
- 同一 topology 两次 CM2D 输出要求字节完全一致；
- read-back 要恢复 vertex/edge/cell 数量；
- owner/neighbour 与 cell loops 必须与内存 topology 完全一致；
- corrupt / unsupported CM2D 显式失败。

## 2D-6C — End-to-end CLI

新增：

- `apps/cartmesh2d_cli.cpp`

CLI：

```text
cartmesh2d_cli <boundary.xy> <output-prefix> [max-level=5] [padding-fraction=0.25] [small-alpha=0.10]
```

边界文件为每行一个 `x y`；最后一点无需重复第一点。

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

## 最终 acceptance fixtures

新增：

- `examples/acceptance/rectangle.xy`
- `examples/acceptance/circle.xy`
- `examples/acceptance/concave.xy`
- `examples/acceptance/airfoil_like.xy`

CMake 已注册四条 CLI end-to-end CTest。每条测试自身完成 VTK/CM2D/JSON 输出与独立 read-back 核对。

## 当前待完成门禁

- [x] quality evaluator implementation
- [x] deterministic JSON report implementation
- [x] VTK export implementation
- [x] solver-topology CM2D export implementation
- [x] CM2D read-back implementation
- [x] end-to-end CLI implementation
- [x] rectangle fixture registered
- [x] circle fixture registered
- [x] concave fixture registered
- [x] airfoil-like fixture registered
- [ ] current-head root CMake configure
- [ ] current-head full 2D build
- [ ] Stage 0~6 unit/regression CTest
- [ ] four CLI end-to-end CTest actual PASS
- [ ] output file independent inspection/count verification
- [ ] branch isolation final audit
- [ ] final acceptance report / CLOSED marker

因此当前 Stage 2D-6 只能记为 **IN PROGRESS / implementation complete enough for first integration build**，不能提前宣称产品完成。
