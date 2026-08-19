# Stage 2D-4 验证记录

## 状态

**IN PROGRESS — 全局 topology 核心、coarse-fine T-junction 处理、真实 Cut-cell 全链路和根 CMake 集成门禁已完成；等待用户显式“验证-4”后对当前最终分支再跑一次封口门禁并正式 CLOSED。**

当前严格停留在 2D-4；未实现 small-cell 聚合、quality/export 或 visualization。

## 目标

把 2D-3 中“每个 leaf 自己正确”的流体 polygon 转换为求解器可用的全局二维拓扑：

- global `Vertex2D`
- global `Edge2D`
- `TopologyCell2D`
- owner / neighbour
- boundary patch
- coarse-fine edge splitting
- deterministic IDs / ordering
- topology audit

## 核心实现

新增：

- `cartmesh2d/include/cartmesh2d/topology/Topology2D.hpp`
- `cartmesh2d/src/topology/Topology2D.cpp`
- `cartmesh2d/tests/topology_test.cpp`
- `cartmesh2d/tests/topology_quadtree_audit_test.cpp`

### 数据模型

`Vertex2D`
- deterministic global id
- canonical point

`Edge2D`
- global id
- canonical `(v0,v1)`
- owner
- optional neighbour
- `BoundaryPatch2D::{None, EmbeddedBoundary, DomainBoundary, Unclassified}`

`TopologyCell2D`
- deterministic global cell id
- source CutCell id/key
- source geometry area
- ordered vertex loop
- ordered edge loop

`TopologyAudit2D`
- duplicate vertices
- duplicate edges
- orphan internal edges
- non-manifold edges
- unclassified boundary edges
- open cell loops
- area mismatches

## coarse-fine / hanging-node 规则

不能直接把每个 fluid polygon 的原始边作为全局 edge。

实现采用：

1. 收集所有非空流体 polygon 顶点；
2. tolerance 下做 deterministic 全局顶点去重；
3. 对每条 cell polygon 边，寻找所有落在该线段上的 global vertex；
4. 按线段参数排序；
5. 将粗边按 hanging-node 拆成 edge fragments；
6. 每个 fragment 用无向 canonical vertex pair 去重；
7. 聚合 incident cells。

因此粗 cell 的一条长 face 可以与两个细 cell 的短 face 正确匹配，而不会保留跨过 hanging node 的 T-junction 长边。

## incident / patch 规则

- 2 incident cells -> internal edge，必须 owner + neighbour，patch=None；
- 1 incident cell -> boundary edge，必须识别为 embedded boundary 或 domain boundary；
- >2 incident cells -> non-manifold，显式失败；
- boundary edge 无法分类 -> `Unclassified` + issue，显式失败。

## fixture 1：人工 coarse-fine 接口

Domain `[0,2] x [0,1]`：

- 左侧一个 `1 x 1` coarse cell；
- 右侧两个 `1 x 0.5` fine cells。

实际拓扑：

- vertices = `8`
- edges = `10`
- internal edge fragments = `3`
- outer boundary fragments = `7`

其中 `x=1` coarse-fine interface 被正确拆成两条长度 `0.5` 的 internal edge：

- coarse <-> lower fine
- coarse <-> upper fine

不存在跨过 `(1,0.5)` hanging node 的长度 `1.0` 长 edge。

原测试曾错误写成 outer boundary fragments = 6；打印完整 edge 列表后确认数学上应为 7（left 1 + bottom 2 + top 2 + right 2），已修正测试期望，算法无需修改。

## fixture 2：自适应圆形 Cut-cell 全链路

`64-segment circle -> Quadtree -> 2:1 balance -> CutCell2D -> global topology`

独立验收统计：

```text
vertices = 125
edges = 172
cells = 48
internal edges = 84
embedded-boundary edges = 88
domain-boundary edges = 0

duplicate vertex = 0
duplicate edge = 0
bad owner/neighbour = 0
bad cell loop = 0

TopologyAudit2D:
duplicateVertices = 0
duplicateEdges = 0
orphanInternalEdges = 0
nonManifoldEdges = 0
unclassifiedBoundaryEdges = 0
openCellLoops = 0
areaMismatches = 0

global fluid-area error = 8.88178e-16
```

除了 builder 自身 `TopologyAudit2D` 外，永久回归中又加入独立 duplicate-vertex、duplicate-edge、owner/neighbour 和 cell-loop 扫描，避免只相信构造器自身报告。

## 根目录集成门禁

从仓库根 CMake 开启二维，并显式构建全部二维测试 target：

```text
cartmesh2d_stage0_geometry_tests .................. Passed
cartmesh2d_stage1_grid_tests ...................... Passed
cartmesh2d_stage2_quadtree_tests .................. Passed
cartmesh2d_stage3_cutcell_tests ................... Passed
cartmesh2d_stage3_quadtree_cutcell_audit_tests .... Passed
cartmesh2d_stage4_topology_tests .................. Passed
cartmesh2d_stage4_quadtree_topology_audit_tests ... Passed

100% tests passed, 0 tests failed out of 7
```

说明 2D-4 不是 standalone-only；它已通过根工程集成路径。

## 分支隔离审计

当前 `agent/native-2d-baseline` 相对 `main` 的改动仍只包含：

- 根 `AGENTS.md` 的二维例外说明；
- 根 `CMakeLists.txt` 的二维可选入口；
- `cartmesh2d/**`。

未修改：

- `include/cartmesh/**`
- 根 `src/**`
- 根 `apps/**`
- 根 `tests/**`

三维 Stage 6 / Stage 7 算法没有被二维拓扑工作改动。

## 当前门禁状态

- [x] global vertex model
- [x] global edge model
- [x] owner / neighbour model
- [x] boundary patch model
- [x] coarse-fine edge splitting algorithm
- [x] Cut-cell embedded boundary integration
- [x] deterministic ordering strategy
- [x] coarse-fine T-junction fixture
- [x] adaptive Cut-cell topology fixture
- [x] root CMake build
- [x] 2D-0/1/2/3/4 全量 CTest
- [x] duplicate vertex/edge 独立审计
- [x] owner/neighbour 独立审计
- [x] cell-loop 独立审计
- [x] branch isolation audit
- [ ] 用户显式“验证-4”后的最终 current-head 复跑与 CLOSED 标记

## 明确未开始

- 2D-5 small-cell detection/agglomeration
- 2D-6 quality/export
- 2D-V visualization
