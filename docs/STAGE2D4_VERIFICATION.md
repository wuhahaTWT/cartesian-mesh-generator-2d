# Stage 2D-4 验证记录

## 状态

**IN PROGRESS — 全局 topology 核心实现和两类硬回归 fixture 已加入；等待完整分支编译/CTest 与最终独立拓扑审计后关闭。**

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

## 已加入 fixture 1：人工 coarse-fine 接口

Domain `[0,2] x [0,1]`：

- 左侧一个 `1 x 1` coarse cell；
- 右侧两个 `1 x 0.5` fine cells。

硬检查：

- hanging node `(1,0.5)` 成为 global vertex；
- coarse interface 被拆成两条长度 `0.5` 的 internal fragments；
- 两条 fragment 均有 owner + neighbour；
- 不存在长度 `1.0` 的未拆 coarse interface edge；
- 每个 cell vertex/edge loop 闭合；
- canonical edge 无重复；
- global vertex 无重复；
- topology area 与 source cell area 一致；
- duplicate source cell 显式拒绝。

## 已加入 fixture 2：自适应圆形 Cut-cell 全链路

`64-segment circle -> Quadtree -> 2:1 balance -> CutCell2D -> global topology`

目标检查：

- 所有非空 Cut-cell 都进入 topology；
- unsupported = 0；
- internal owner/neighbour edge > 0；
- embedded boundary edge > 0；
- 圆位于 domain 内，因此 domain-boundary edge = 0；
- duplicate/orphan/non-manifold/unclassified/open-loop/area-mismatch 全部 = 0；
- topology 总流体面积等于 2D-3 Cut-cell 总流体面积 / 输入 polygon 面积；
- 重复运行 vertex/edge/cell 顺序确定。

## 当前门禁状态

已完成代码与测试定义：

- [x] global vertex model
- [x] global edge model
- [x] owner / neighbour model
- [x] boundary patch model
- [x] coarse-fine edge splitting algorithm
- [x] Cut-cell embedded boundary integration path
- [x] deterministic ordering strategy
- [x] coarse-fine T-junction fixture
- [x] adaptive Cut-cell topology fixture
- [ ] 当前 GitHub 完整 branch CMake build
- [ ] 2D-0/1/2/3/4 全量 CTest
- [ ] 最终 duplicate vertex/edge 独立审计
- [ ] 最终 branch isolation audit

在以上最终门禁完成前，Stage 2D-4 不标记 CLOSED。

## 明确未开始

- 2D-5 small-cell detection/agglomeration
- 2D-6 quality/export
- 2D-V visualization
