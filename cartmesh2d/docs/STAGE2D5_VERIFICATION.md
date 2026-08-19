# Stage 2D-5 验证记录

## 状态

**IN PROGRESS — 2D-5A small-cell 检测/候选分析与 2D-5B topology-safe agglomeration 均已实现，并通过当前根工程回归；等待用户显式 `验证-5` 后对当前最终 head 做封口复跑并正式 CLOSED。**

当前严格停留在 small-cell stabilization；未实现 2D-6 quality/export 或 visualization。

## 2D-5A — Detection / Candidate

已实现：

- configurable `areaFractionThreshold`
- Cut-cell alpha histogram
- deterministic `alpha < threshold` marking
- topology internal-edge neighbour discovery
- aggregate shared-interface length
- stable-neighbour preference
- deterministic tie-break: shared length -> target alpha -> target area -> topology id
- explicit unresolved report

人工 sliver：`alpha=0.1`，threshold `0.2` 时 small=1 且找到稳定邻居；threshold `0.1` 时不标记 small。

真实 shifted-circle fixture（中心偏移 `(0.07,0.03)`）：

```text
fluid cells = 56
cut cells = 32
full cells = 24
small cells = 8
unresolved = 0
minimum small alpha ~= 0.00120311
small -> small candidate = 0
```

## 2D-5B — Topology-safe agglomeration

新增：

- `AgglomeratedCell2D`
- `AgglomerationResult2D`
- `AgglomerationIssueCode2D`
- `agglomerateSmallCells(...)`
- `src/stabilization/Agglomeration2D.cpp`
- `tests/agglomeration_test.cpp`

### 设计原则

聚合不把结果伪装成普通单-background-leaf `CutCell2D`。公开结果使用独立 `AgglomeratedCell2D` 保存：

- member topology cell ids
- member source ids
- real merged polygon
- area
- centroid

`CutCell2D` 只在内部作为已验证 2D-4 topology builder 的临时 geometry adapter；聚合后的公开语义不会错误继承原始 `areaFraction`。

### 几何聚合算法

1. 使用 2D-5A candidate graph；
2. 当前第一版只允许 small -> non-small target；small -> small 链式聚合显式失败；
3. 将同一 agglomeration group 中所有 Stage 2D-4 directed edge fragments 收集；
4. 组内共享 edge 必须恰好成对、方向相反，并从 exterior 中消去；
5. 剩余 directed edges 必须形成 one-in/one-out 的单闭环；
6. 分叉、断链、多 boundary loop、非流形 incident 立即失败；
7. polygon orientation 统一；
8. tolerance-aware 删除共线冗余顶点；
9. merged polygon area 必须等于成员 cell area 之和；
10. 用 merged polygons 重建全局 topology；
11. 重建 topology 必须再次通过完整 Stage 2D-4 audit。

因此 coarse-fine 接口已经由 2D-4 拆成 edge fragments 后，2D-5B 可以直接在规范拓扑上安全消除组内共享界面，而不会重新制造 T-junction。

## 2D-5B 人工 sliver fixture

两个相邻 fluid cells：

- small Cut-cell：area fraction `0.1`
- stable full cell：area `1.0`

threshold `0.2` 后：

```text
input cells = 2
merged small cells = 1
output cells = 1
area before = 1.1
area after = 1.1
area error <= 1e-12
rebuilt topology = PASS
```

共享内部边被真正移除；矩形 union 的共线冗余点被简化，最终 merged polygon 为 4 顶点矩形。

另有 no-small identity case，确保阈值下没有 tiny cell 时 cell count 与总面积保持不变。

非法 self-target candidate 会显式失败，不会生成假 union。

## 2D-5B 真实 shifted-circle 全链路

```text
64-segment shifted circle
-> Quadtree level 4
-> 2:1 balance
-> CutCell2D
-> global topology
-> 2D-5A
-> 2D-5B agglomeration
-> rebuilt global topology
```

当前确定性门禁：

```text
small cells detected = 8
unresolved = 0
merged small cells = 8
output cells = input cells - 8
total fluid-area error <= 1e-10
```

重建后的 `TopologyAudit2D`：

```text
duplicateVertices = 0
duplicateEdges = 0
orphanInternalEdges = 0
nonManifoldEdges = 0
unclassifiedBoundaryEdges = 0
openCellLoops = 0
areaMismatches = 0
```

每个 stabilized cell 都必须具有正面积、至少 3 个 polygon vertices，并且 merged polygon 通过 simple-boundary diagnostics。

## 当前根目录 CMake / CTest

按与 GitHub 推送相同的独立文件布局：

- `SmallCell2D.cpp`：2D-5A
- `Agglomeration2D.cpp`：2D-5B

从仓库根 CMake 开启二维并显式构建全部二维测试 target：

```text
cartmesh2d_stage0_geometry_tests .................. Passed
cartmesh2d_stage1_grid_tests ...................... Passed
cartmesh2d_stage2_quadtree_tests .................. Passed
cartmesh2d_stage3_cutcell_tests ................... Passed
cartmesh2d_stage3_quadtree_cutcell_audit_tests .... Passed
cartmesh2d_stage4_topology_tests .................. Passed
cartmesh2d_stage4_quadtree_topology_audit_tests ... Passed
cartmesh2d_stage5a_small_cell_tests ............... Passed
cartmesh2d_stage5b_agglomeration_tests ............ Passed

100% tests passed, 0 tests failed out of 9
```

一次并行 `make` 曾出现本地子进程 wait 异常；改为 `-j1` 后 exact layout 完整编译与 9/9 CTest PASS。该异常不是编译器错误或测试失败。

## 分支隔离

当前 branch 相对 `main` 仍只修改：

- 根 `AGENTS.md` 二维例外说明；
- 根 `CMakeLists.txt` 二维可选入口；
- `cartmesh2d/**`。

未修改三维算法目录：

- `include/cartmesh/**`
- root `src/**`
- root `apps/**`
- root `tests/**`

## Stage 2D-5 关闭门槛

- [x] alpha histogram
- [x] threshold marking
- [x] best-neighbour candidate
- [x] deterministic candidate selection
- [x] unresolved explicit report
- [x] real adaptive tiny-cell fixture
- [x] topology-safe cell agglomeration
- [x] pre/post stabilization total area conservation
- [x] post-agglomeration topology audit = PASS
- [x] no negative area / duplicate edge / non-manifold topology after treatment
- [x] current implementation root CMake / CTest = 9/9 PASS
- [ ] 用户显式 `验证-5` 后 current-head 封口复跑与 CLOSED 标记

因此当前状态为：**2D-5A PASS / 2D-5B PASS / Stage 2D-5 READY FOR VALIDATION**。
