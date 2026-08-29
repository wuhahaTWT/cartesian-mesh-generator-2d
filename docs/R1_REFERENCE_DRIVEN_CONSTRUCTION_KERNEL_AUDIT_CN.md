# R1：参考实现驱动的构造内核与局部拓扑架构审计

日期：2026-08-30
审计工作树：`/Users/Zhuanz/Desktop/cartesian-mesh-generator-2d`
参考源码根目录：`/Users/Zhuanz/Desktop/mesh-reference-sources`

## 0. 结论与证据等级

本报告只做事实核查、参考源码审计和架构设计，没有修改网格算法，没有合并
Q2-B，也没有运行完整 CTest 或五案例 OpenFOAM。

文中使用三种证据等级：

- **[F] 源码确认**：由当前 checkout、已有机器可读 artifact 或本轮定向运行确认；
- **[I] 推断**：由参考实现和当前调用链推导，尚未在 cartmesh2d 中实现；
- **[H] 假设**：必须通过后续小实验才能确认，不能作为完成声明。

推荐选择：**以路线 B 为主，吸收路线 C 的轻量数据模型**。保留 H1--H4 的上层
接口、物理语义、boundary-layer 规划、CLI/导出和验证合同；替换 transition/cut
construction kernel，并引入统一稳定 ID、source lineage、feature compatibility、
edge-incidence graph、patch-local transaction 和增量质量缓存。当前证据不支持立即
引入完整 DCEL。

理由：

1. **[F]** Q2-A 的共享交点身份和公共分割已经进入真实 solver 链路，10 万规则格
   隔离 topology benchmark 从 `0.982997 s` 降到 `0.449512 s`，应保留其逻辑身份
   思想；但 registry 的 legacy proximity path 仍线性扫描，source metadata 回绑仍
   是 cell × source 全扫描。
2. **[F]** Q2-B 证明了 narrow-gap 和 sharp-tail 存在面积守恒的局部可行拓扑，
   但每个候选仍重建全局 topology 并重算全局 solver quality。当前最终产物实际记录
   `40` 和 `1832` 个候选，不应演化成生产架构。
3. **[F]** Q2-B 只消除了两例的 dimensionless hard short-face；两例完整 Q1 仍为
   `FAIL`，主要剩余问题是 face weight、volume ratio、interior angle、non-orthogonality
   等。只扩大 snap 或继续事后组合搜索无法达到共同的精度/速度目标。
4. **[I]** 一个 half-edge-lite/edge-incidence patch model 已足够支持二维流形审计、
   patch 边界锁定和局部提交；完整 DCEL 的迁移成本与必要性尚无源码证据。

## 1. Git、阶段与工作区现场事实

### 1.1 Git

- **[F]** 当前分支：`codex/q2b-constrained-local-repair`。
- **[F]** 当前 HEAD：`2e71f45e215c1b56300f9556595438de204ec479`，提交说明
  `Complete Q2-B constrained local short-face repair`。
- **[F]** `main` 与 `origin/main` 都指向
  `eb7ab48880dd2274960a7d532dcdec3454739ec6`，即 Q2-A。
- **[F]** Q2-B 相对 Q2-A 修改 7 个文件，`690 insertions / 4 deletions`；主要新增
  位于 `src/quality/SolverTopology2D.cpp`。
- **[F]** 本轮开始与报告写入前，代码工作树没有用户未提交修改。Q2-B 未合并 main。

### 1.2 阶段边界

- **[F]** `AGENTS.md` 明确默认流体域为 `Domain2D - solid interior`，并要求真实
  fluid polygon、面积守恒、统一 owner/neighbour、确定性以及真实 OpenFOAM 验证。
- **[F]** H4-3 已形成局部降层、stepped envelope、termination buffer 和事务回退；
  pure Cut-cell 只允许作为最后 fallback。事实源为
  `docs/STAGE2DH4_3_LOCAL_TERMINATION_CN.md`。
- **[F]** Q0 区分 construction/solver/OpenFOAM 三类质量；Q1 是更严格的无量纲、
  分类型诊断合同。OpenFOAM `Mesh OK` 不等于 Q1 通过。
- **[F]** Q2-A 是 shared construction 基础设施，`q2_full_status` 仍为
  `PARTIAL_NOT_ACCEPTED`。Q2-B 是两个短面案例的 correctness prototype，不是完整
  Q1 修复，也未执行冻结后的五案例正式验收。

## 2. 当前实现的源码级调用链与瓶颈

### 2.1 交点与公共分割

当前真实路径：

```text
HybridMesh2D
  -> IntersectionRegistry2D::configureGrid/registerSegment/intersectGridLine
  -> CutCell2D 携带 canonicalVertexIds
  -> buildGlobalTopology
  -> SharedEdgePartition2D::partition
  -> owner/neighbour + boundary patch audit
```

- **[F]** `include/cartmesh2d/geometry/IntersectionRegistry2D.hpp:16-54` 已定义 source、
  feature、canonical vertex 和 provenance record，但 `CanonicalVertex2D` 只有一个
  `supportId`，没有可表达 narrow-gap 两侧 ownership 集合、snap 拒绝原因或 grid/
  wall/transition 的统一 typed key。
- **[F]** construction event 使用 `(supportId, GridLineIdentity2D)` 缓存；grid line
  是 axis + 最大层级整数坐标。见
  `src/geometry/IntersectionConstruction2D.cpp:131-201`。
- **[F]** arithmetic snap budget 是
  `64 * machine_epsilon * min(query_h, grid_h, endpoint_h, support_length)`，只处理舍入，
  不是 Q1 大小的几何移动。见同文件 `:152-185`。
- **[F]** legacy `IntersectionRegistry2D::canonicalize` 对 registry 全部 vertex 做线性
  扫描，复杂度为每次查询 `O(V)`；见
  `src/geometry/IntersectionRegistry2D.cpp:67-106`。
- **[F]** `SharedEdgePartition2D` 建立按 x/y 的有序索引，并缓存无向端点区间；但一般
  斜边仍对 bounding rectangle 内候选调用 `pointOnSegment`。见
  `src/topology/SharedEdgePartition2D.cpp:7-43`。
- **[F]** `buildGlobalTopology` 每次都复制、排序全部 cells，重建全部 vertices、edge
  uses、owner/neighbour、patch 和面积审计；见 `src/topology/Topology2D.cpp:242-505`。

结论：Q2-A 已消除重复求交与部分公共边重复分割，但还不是 arrangement，也没有
patch mutation API。`std::map<double,...>` 依赖精确坐标键，适合作为当前桥接，不适合
成为最终统一身份模型。

### 2.2 solver repair 与 Q2-B

当前路径：

```text
buildHybridMesh2D
  -> buildSolverTopology2D
  -> qualityMetadataForSolver
  -> repeat <= 32:
       repairSolverShortFaces2D
         -> scan all edges for worst short face
         -> enumerate pair / three-cell candidates
         -> rebuild all cells and buildGlobalTopology(candidate)
         -> evaluateSolverQuality2D(candidate over all cells and edges)
       commit one candidate
  -> global interface audit
  -> global solver quality
  -> global Q1 contract
```

- **[F]** 外层每次只提交一个 transaction，最多 32 次；见
  `src/hybrid/HybridMesh2D.cpp:1062-1104`。
- **[F]** Q2-B 每轮扫描所有 edges 形成 face score 与最坏短面；见
  `src/quality/SolverTopology2D.cpp:1374-1459`。
- **[F]** pair/triple 候选的局部 polygon union/split 本身是 patch-local 的，但
  `agglomerateCellPair`、`repartitionPair`、`replaceCellPatch` 会把全网格重新包装为
  CutCell 并调用 `buildGlobalTopology`。例如同文件 `:695-729`、`:732-766` 和
  `:1534-1631`。
- **[F]** 每个 topology-valid candidate 都执行完整 `evaluateSolverQuality2D`；见
  `:1534-1557`。该质量函数遍历全部 cells 和全部 edges，见
  `src/quality/SolverQuality2D.cpp:43-240`。
- **[F]** solver cell 到 hybrid source 的回绑在每个 cell centroid 上顺序扫描所有
  source polygon，复杂度最坏 `O(C*S*P)`；见
  `src/hybrid/HybridMesh2D.cpp:288-325`。这正是需要 source lineage 的直接证据。
- **[F]** topology rebuild 的主要容器是有序 `std::map`，单次通常约
  `O((V+E) log(V+E))`；全质量是 `O(V+E)`。因此 Q2-B 总成本近似
  `O(K * ((V+E)log(V+E)))`，`K` 为累计候选数，尚未计 polygon union/split。
- **[F]** 最终文件的现场定向读取显示 narrow-gap 为 40 candidates / 4 commits，
  sharp-tail 为 1832 candidates / 14 commits。仓库内旧汇总 JSON 仍写 6/198，已证实
  与当前最终生成物不一致；应在下一次证据冻结时重新生成，而不是继续引用旧数。

### 2.3 构造质量没有进入 template 决策

- **[F]** H4 termination growth 候选仍是固定 `{1.45, 1.55, 1.50}`，候选完整构建
  后才通过 solver-quality 门。
- **[F]** Q1 contract 在 solver topology 形成后评估；它能定位问题，但没有参与
  transition/template 的局部选择。
- **[I]** face weight、volume ratio、non-orthogonality 都能在一个候选 patch 的
  polygon centroid、area 与边 incidence 已知后局部计算，无需等到全局重建。

## 3. 修改前的可复现基线

### 3.1 当前五案例事实

下表的 circle/superellipse/concave-L 来自 Q2-A `artifacts/q2a/comparison.json` 和
对应 `build-q2a/evidence/shared-*.quality-contract.json`；narrow-gap/sharp-tail 的
最终列来自当前 Q2-B `build-q2b/final-targeted/`，并在本轮用只读定向脚本复核。
所有 `max/min` 都是 solver/Q1 内部公式，不冒充外部 OpenFOAM。

| case | solver cells | min face/local_h | max non-orth | min face weight | min volume ratio | area error | Q1 |
|---|---:|---:|---:|---:|---:|---:|---|
| circle (Q2-A) | 728 | 0.0459428 | 55.3968 | 0.0872106 | 0.0270379 | -1.07e-14 | FAIL |
| superellipse (Q2-A) | 795 | 0.0162229 | 63.4869 | 0.0780183 | 0.0191827 | 3.02e-14 | FAIL |
| concave L (Q2-A) | 5452 | 0.0144000 | 69.8254 | 0.0570773 | 0.0106997 | -6.75e-14 | FAIL |
| narrow gap (Q2-B) | 3185 | 0.0146196 | 68.8133 | 0.0515745 | 0.0118533 | -3.38e-14 | FAIL |
| sharp trailing edge (Q2-B) | 3379 | 0.0120053 | 69.3953 | 0.0513456 | 0.0108295 | 2.13e-14 | FAIL |

Q2-B 后仍存在的 hard issue：

| case | face weight | volume ratio | min angle | non-orth | skew | hydraulic aspect |
|---|---:|---:|---:|---:|---:|---:|
| narrow gap | 95 | 295 | 20 | 14 | 4 | 2 |
| sharp trailing edge | 133 | 152 | 23 | 11 | 2 | 2 |

因此 **[F]** “Q2-B 两例短面通过”与“完整 Q1 仍失败”必须同时保留。

### 3.2 性能、确定性与证据缺口

- **[F]** Q2-A 五案例历史端到端 shared wall time：circle `0.106 s`、superellipse
  `0.181 s`、concave-L `5.084 s`、narrow-gap `43.171 s`、sharp-tail `1.542 s`。
- **[F]** 本轮使用现有 Q2-B 二进制现场单次运行：narrow-gap `42.30 s`，sharp-tail
  `4.91 s`。这不是严谨 benchmark，但 sharp-tail 已约为 Q2-A 历史值的 `3.18x`，
  与 1832 次候选全局审计一致。
- **[F]** 当前 Q2-B 两例重复生成的 solver CM2D、quality、geometry、construction
  和 intersection 文件逐字节一致；boundary patch 计数和面积误差未变化。
- **[F]** `/usr/bin/time -lp` 在当前 sandbox 中因 `sysctl kern.clockrate` 权限未输出
  peak RSS。当前 hybrid JSON 也没有各阶段 profile 字段。
- **[H]** 五案例统一 peak RSS、construction/intersection/topology/quality 各阶段耗时
  仍未知；正式 R1 修改前必须补齐，不得从 H4-3 旧平台数据外推。

### 3.3 建议的基线记录格式

新增一个非 Git 大文件目录和一个可入 Git 的小 manifest：

```text
build-r1-baseline/<case>/...
artifacts/r1/baseline-manifest.json
```

manifest 每例至少记录：source commit、命令、输入 hash、参数、solver cell/face、
三项 face ratio、max non-orth、min face weight、min volume ratio、面积误差、boundary
patch counts、construction/solver hash、wall/user/sys、peak RSS、各阶段秒数、registry
query/hit、source-index query/hit、global/patch topology rebuild count、full/local quality
count。必须把“测不到”写成 `null + reason`，不能省略字段。

## 4. 参考源码逐项目审计

`mesh-reference-sources/papers/` 当前文件数为 0；本报告不声称审读了本地论文。
五个源码目录均无 `.git` 元数据，因此版本身份只采用目录名/源码自报版本，不能给出
上游 commit hash。

### 4.1 AMReX development / EB2

许可证：**[F] BSD-3-Clause 风格**；`amrex-development/LICENSE:1-31`。分发源码需
保留 copyright、条件和免责声明，二进制分发需在材料中复现；另保留 NOTICE。

核心位置与事实：

- `Src/EB/AMReX_EB2_GeometryShop.H:41-175`：在 Cartesian edge 上用 Brent root
  finder 求隐式面截距。double 路径含固定 `1e-12` 终止 tolerance，不能直接采用为
  cartmesh2d 的尺度合同。
- `Src/EB/AMReX_EB2_2D_C.cpp:198-350`：由 edge intercept 构造 x/y face area
  fraction 和 face centroid；一个 cell 多于两个 cuts 时默认 abort，可选 cover。
- 同文件 `:7-137`：从 face aperture 计算 volume fraction、cell/boundary centroid 和
  normal；这是构造期保存局部几何量、避免事后重复积分的可复用思想。
- 同文件 `:158-194`：`vfrac < small_volfrac` 时直接设 covered。该 small-cell 删除
  与本项目真实性规则冲突，不能复用。
- `Src/EB/AMReX_EB2_Graph.H:22-105`：Cell/Face/Neighbor/Edge/Vertex 构成轻量 EB
  graph；`VertexID = pair<IntVect,int>` 表达一个网格 cell 内的第几个 EB vertex。
- `Src/EB/AMReX_EB2_Level.H:230-250`：以 patch arrays 保存 flags、volfrac、centroid、
  boundary data、area fraction；天然适合 patch-local 更新和线性内存访问。

对问题的回答：

- 交点唯一性依赖 Cartesian edge array slot，而不是全局距离聚类；同一 grid edge
  只有一个 intercept slot。
- local scale 由 `dx` 和归一化 face/volume fractions 主导，但 root finder 仍含固定
  absolute tolerance。
- 不保存输入 polyline 的 feature/source lineage；sharp corner 与 narrow-gap 两侧
  ownership 不是 EB2 2D data model 的强项。
- refinement/transition 属于 AMR patch 层，不输出本项目所需的共形 polygon solver
  mesh。multi-cut 和 small-cell 的默认处置更不能直接采用。

可复用思想：**edge-slot identity、patch-array cache、构造时 area/centroid/aperture**。
不可复制部分：small-cell covering、fixed root tolerance、single-valued cut 限制。

### 4.2 p4est 2.8.7

许可证：**[F] GPL-2.0**；`p4est-2.8.7/COPYING`。只借鉴公开算法思想；若复制或链接
源码，必须评估 GPL 对整个分发的影响并保留版权/许可证。

核心位置与事实：

- `src/p4est.h:61-125`：quadrant 使用整数 `x,y,level`；边长由 level 的位移精确定义。
- `src/p4est_bits.h:57-150`：Morton identity、比较和 hash 基于逻辑坐标与 tree，
  不依赖浮点 proximity。
- `src/p4est_extended.h:445-494`：refine/coarsen 后显式执行 2:1 balance。
- `src/p4est_mesh.h:68-163`：平衡 forest 的 `quad_to_quad/quad_to_face/quad_to_half`
  直接表达邻接和 hanging neighbours，内存为 `O(N)`。
- `src/p4est_search.h:35-55,114-136,221-230`：top-down local tree search；源码明确
  建议在可能时维护 `O(1)` tree context，避免反复 binary search。
- `src/p4est.h:160-201`：forest 有 revision，refine/coarsen/partition/balance 的有效
  修改会递增，适合缓存失效版本号。

对问题的回答：p4est 不做 wall intersection、feature snap、constrained polygon 或
solver quality；它提供的是**精确背景 cell identity、2:1 propagation、邻接编码、
revisioned local search**。这些正适合 R1 的 grid/source key 和局部 refinement 队列，
但不能替代二维 cut construction。

### 4.3 cfMesh 1.2.0

许可证：**[F] GPL-3.0-or-later**；`cfMesh-1.2.0/README` License 段及每个源码头。

核心位置与事实：

- `meshLibrary/cartesian2DMesh/cartesian2DMeshGenerator/cartesian2DMeshGenerator.C:63-246`
  是真实 workflow：octree template -> surface topology cleanup -> projection -> patch
  assignment -> edge/corner remap -> surface optimize -> layer -> refine layer。
- 同文件 `:77-100` 会循环修复 irregular connections，并调用
  `checkNonMappableCellConnections(...).removeCells()`。这违反本项目“不删除坏 cell
  隐藏问题”的规则。
- 同文件 `:102-135` 将 mesh boundary 映射到 surface 后再 optimize/untangle；这是
  body-fitted morphing，不是保留真实 CutPolygon 的 cut-cell construction。
- `utilities/octrees/meshOctree/meshOctreeModifier/
  meshOctreeModifierEnsureCorrectRegularity.C:44-178`：从被标记 leaf 经真实 neighbour
  查询传播 regularity refinement，是局部 propagation 的直接实现位置。
- `utilities/octrees/meshOctree/meshOctreeAddressing/
  meshOctreeAddressingCreation.C:649-1106`：建立 face owner/neighbour、leaf faces 和
  transition/hanging addressing，避免事后几何猜邻接。
- `utilities/triSurfaceTools/triSurfaceDetectFeatureEdges/
  triSurfaceDetectFeatureEdgesFunctions.C:40-103`：按相邻 facet normal 与 angle
  tolerance 分类 feature edge，并保留非流形/边界类别 bit。

对问题的回答：feature 先分类再参与 patch/edge mapping；regularity 由 neighbour graph
传播；topology addressing 是独立缓存。它没有可直接移植的统一 intersection lineage，
且其删除/形变式修复与本项目硬约束不兼容。

### 4.4 Gmsh 4.15.2

许可证：**[F] GPL-2.0-or-later with stated exception**；
`gmsh-4.15.2-source/LICENSE.txt:1-16`。独立重写思想；复制核心代码或链接 Gmsh 会带来
GPL 分发义务，第三方 `contrib/` 还需逐项核对许可证。

核心位置与事实：

- `src/geo/MVertex.h:23-90`：mesh vertex 有创建后稳定的 `_num`，并持有 `GEntity*
  _ge`；ID 与临时 generation index 分离。
- 同文件 `:137-166`：edge vertex 额外保存曲线参数 `u` 与 local mesh size `lc`。
- `src/mesh/meshGFaceDelaunayInsertion.cpp:14,52-60,495-572`：segment test、incircle、
  orientation 使用 robust predicates，而非扩大 epsilon。
- 同文件 `:1847-2035`：Delaunay 建网后显式 recover constrained edges；恢复通过
  edge swaps 执行。
- `src/mesh/meshGFace.cpp:644-742`：约束边无法恢复时，不吞掉 feature，而是沿原
  `GEdge` 参数域插入中点、重建 line 并重新网格化相邻 faces。
- `src/mesh/BackgroundMeshTools.cpp:132-258`：mesh size 组合实体、曲率、point/field
  信息；`Field.cpp` 还支持多个 size/metric field 的 intersection。

对问题的回答：Gmsh 最有价值的是**stable mesh ID + geometry entity lineage + robust
predicate + constrained-edge recovery failure triggers source-edge resampling**。它的主
2D kernel 是非结构 Delaunay/BAMG，不保持 Cartesian/quadtree/H4 cell identity，不能
整体移植；但可指导 unsafe intersection 的 fallback 顺序。

### 4.5 OpenFOAM v2606 / snappyHexMesh / hexRef8

许可证：**[F] GPL-3.0-or-later**；`OpenFOAM-v2606-source/LICENSE.md:1-7`。

核心位置与事实：

- `src/dynamicMesh/polyTopoChange/polyTopoChange/hexRef8/hexRef8.H:75-96`：持久保存
  `cellLevel`、`pointLevel`、refinement history 以及 saved level maps。
- `hexRef8.C:1565-1651`：通过 owner/neighbour level 在 faces 上传播标记，直到满足
  2:1；这是基于邻接的 refinement closure，不是全域几何重分类。
- `hexRef8.H:119-276,495-502`：所有 point/face/cell change 先记录到
  `polyTopoChange`，随后由 map 更新 lineage；其 transaction/delta 思想值得借鉴，
  但它是三维 split-hex 实现，不能压成 z=0。
- `src/mesh/snappyHexMesh/meshRefinement/meshRefinementGapRefine.C:1229-1315`：用
  local cell size 和期望 gap cell 数计算 gap size/refinement level，分别保留 surface
  与 shell ownership。
- `src/mesh/snappyHexMesh/snappyHexMeshDriver/snapParameters/snapParameters.H:61-71`
  分离 snap tolerance、explicit feature snap、implicit feature snap。
- `src/mesh/snappyHexMesh/snappyHexMeshDriver/snappyLayerDriver.C:4000-4089`：layer
  iteration 在 move 前组合 quality constraints 并检查候选；但后期可切换 relaxed
  constraints。后者违反本项目“不降低阈值”，不能采用。
- `src/meshTools/cellQuality/cellQuality.C:42-166` 是 external solver quality 公式的
  独立事实源之一，但不能代替真实 `checkMesh`。

对问题的回答：OpenFOAM 提供**refinement lineage、face-local consistency closure、
delta topology、gap size 与 cell size 联动、feature snap 开关、构造过程质量门**。
不能复制三维 template、cell removal、relaxed quality 或将默认 checkMesh 当作 Q1。

## 5. 参考实现到本项目的源码映射表

| 参考实现/算法 | 源码位置 | 解决的问题 | 复杂度/内存 | 可复用思想 | 不能直接复用 | 许可证/署名 | 对 H4/Q 的影响 |
|---|---|---|---|---|---|---|---|
| p4est logical quadrant identity | `p4est.h:61-125`, `p4est_bits.h:57-150` | 浮点 cell identity 与确定性 | lookup/hash 期望 O(1) 或 ordered O(log N)，O(N) | typed dyadic `CellKey/GridVertexKey` | 无 wall/cut/quality | GPL-2.0；文档署名，独立实现 | 保留 H1/H2；替换散落浮点 grid key |
| p4est 2:1 balance + mesh adjacency | `p4est_extended.h:478-494`, `p4est_mesh.h:68-163` | refinement propagation/hanging neighbour | 通常与受影响 leaves/closure 近线性；O(N) | 邻接驱动 closure、half neighbour encoding、revision | 并行 forest API 不直接引入 | GPL-2.0 | 强化 H2/H4 termination，不改变物理边界 |
| AMReX edge-slot EB data | `AMReX_EB2_2D_C.cpp:198-350`, `AMReX_EBData.H:21-28` | 同一 grid edge 截距共享、局部矩量 | 每 patch O(cells+faces)，array O(N) | edge slot identity、构造时 aperture/centroid cache | fixed tol、multi-cut cover、small-cell delete | BSD-3-Clause + NOTICE | Q2-A identity 扩展；加速局部 Q metric |
| AMReX lightweight EB graph | `AMReX_EB2_Graph.H:22-105` | irregular cell 邻接 | O(V+E) | cell-local vertex id、typed face/neighbor | EB2 graph 不输出共形 polygon mesh | BSD-3-Clause | R1-D 的 edge-incidence 参考 |
| cfMesh regularity propagation | `meshOctreeModifierEnsureCorrectRegularity.C:44-178` | 局部 refinement closure | 约 O(k*d) 每轮，k 为 closure | neighbour queue/bitset 标记 | OpenFOAM/cfMesh data types、3D octree | GPL-3+ | H2/H4 局部 refinement fallback |
| cfMesh cached addressing | `meshOctreeAddressingCreation.C:649-1106` | owner/neighbour 与 transition addressing | O(N+E) storage | topology cache 与 geometry 分离 | body-fitted template、cell deletion | GPL-3+ | 替换 H4 事后全局猜测，不复用删除路径 |
| cfMesh feature classification | `triSurfaceDetectFeatureEdgesFunctions.C:40-103` | sharp/boundary/nonmanifold feature | O(surface edges) | feature first, movement second | 3D facet-normal 分类公式不直接用于 2D | GPL-3+ | R1-A feature class/ownership |
| Gmsh stable vertex + entity lineage | `MVertex.h:23-90,137-166` | ID、临时 index、source entity 混淆 | O(V) | stable ID、entity ref、curve parameter/local size | GEntity 指针模型过重 | GPL-2+ | 统一 grid/wall/transition identity 与 provenance |
| Gmsh robust constrained Delaunay | `meshGFaceDelaunayInsertion.cpp:52-60,1847-2035` | predicate robustness、feature edge 保留 | 常见 O(n log n)，recover 依 cavity/swaps | exact/adaptive predicate；约束必须恢复或失败 | 非 Cartesian kernel；GPL code | GPL-2+ | transition patch 的最后局部 partition 选项 |
| Gmsh unrecovered-edge resampling | `meshGFace.cpp:644-742` | 不安全约束不能静默丢失 | 只重网格相邻 faces，近似 patch-local | 沿 source parameter 插点并重试 | 全 face remesh 仍偏重 | GPL-2+ | R1-C wall resampling fallback |
| OpenFOAM hexRef8 consistency | `hexRef8.C:1565-1651` | refinement 2:1 closure | face sweeps至收敛；O(iter*E) | owner/neighbour level propagation | 3D split hex | GPL-3+ | R1-C refinement propagation |
| OpenFOAM polyTopoChange | `hexRef8.H:119-276,495-502` | 局部拓扑 delta 与映射 | O(changed patch + mapping) 目标 | add/modify/remove delta、commit map、history | 三维类型与复杂 framework | GPL-3+ | R1-D transaction 模型 |
| OpenFOAM gap refinement | `meshRefinementGapRefine.C:1229-1315` | narrow gap 被两侧错误吸附 | spatial query + local propagation | gap ownership、cells-across-gap、local cell size | snappy 3D surface modes | GPL-3+ | R1-A/C gap pair identity 与 refinement |
| OpenFOAM construction quality gate | `snappyLayerDriver.C:4000-4089` | layer movement 产生坏 cell | 每候选检查受影响 faces | quality 在构造决策中参与 | relaxed thresholds 禁止；不照搬 3D mover | GPL-3+ | R1-E quality-guided template |

## 6. 三条重构路线比较

| 项目 | A：现有 H4 + 构造期约束内核 | B：重做 transition/cut kernel，保留上层 | C：引入成熟参考数据模型核心 |
|---|---|---|---|
| 精度上限 | 中。能解决 near-coincident，但旧 source 回绑、全局 topology 与 template 后验质量仍限制 | 高。可在 polygon/edge 生成前阻止坏 topology，并保留 H4 物理语义 | 潜在最高，但取决于选择的数据模型；完整 DCEL 不等于自动正确 |
| 时间预期 | 短期最好；若 Q2-B 保留则最坏仍 `O(K*N log N)` | 目标构建 `O((N+I)log N)`，局部事务 `O(k log N + k log k)`，阶段末一次全局 O(N log N) | 初期最慢；成熟后可达 B，但迁移/调试常数大 |
| 数据结构 | 扩展 registry + index；旧 topology 不变 | typed stable IDs、source store、edge incidence、patch delta、local metric cache | half-edge/DCEL/arrangement 或外部框架式 model |
| H1--H4 影响 | 最小，但 H4 termination 仍被旧 rebuild 约束 | H1/H2 sizing/quadtree 保留；H4 layer planning 保留；transition/remainder construction 适配 | H3/H4 topology 几乎整体迁移，回归风险最大 |
| 渐进迁移 | 容易，双路径开关 | 可按 R1-A--G 逐步 shadow/compare | 困难，需长时间双表示与转换器 |
| 最大失败风险 | 做成“更聪明的 snap”，仍无法解决全 Q1 和速度 | patch 边界/lineage 不完整导致局部提交与全局不一致 | scope 爆炸、GPL 污染、确定性/接口漂移 |
| 建议 | 仅作为 R1-A/B 的过渡，不作为终局 | **推荐主路线** | 只吸收 half-edge-lite/edge-incidence、stable ID、delta 思想；暂不完整引入 |

路线 A 单独不足的源码证据是：`qualityMetadataForSolver` 的 `O(C*S)`、Q2-B 的 candidate
global rebuild、固定 termination candidates 都不由 registry 扩展自动消失。路线 C 暂不
选完整 DCEL 的原因是：当前二维 topology 的必要查询只有 vertex identity、edge
incidence、cell loop、patch boundary 与局部 replacement；尚无高频任意 arrangement
walk、复杂布尔层叠或多 region overlay 证明 DCEL 必要。

## 7. 推荐的新底层架构

### 7.1 稳定身份与 source lineage

建议核心类型：

```text
VertexId = monotonic stable id (never coordinate-derived after creation)
VertexKey = variant<
  GridVertexKey(tree, ix, iy, level-normalized),
  WallVertexKey(loop, sourceVertex),
  WallGridIntersectionKey(sourceSegment, GridEdgeKey),
  TransitionVertexKey(templatePatch, localVertex),
  PatchGeneratedKey(transaction, localOrdinal)>

SourceRef = (kind, objectId, subEntityId, parameterRange, side/region)
FeatureClass = Smooth | ConvexSharp | ConcaveSharp | GapSideA | GapSideB |
               DomainCorner | Grid | TransitionFixed | TransitionMutable

VertexRecord = {
  id, key, originalPosition, position, localH,
  featureClass, featureOwner, sourceRefs[],
  displacement, decision, decisionReason, creationRevision
}
```

关键不变量：

1. feature compatibility 先于距离；不同 `featureOwner`、gap pair 两侧或两个非 incident
   source 不得 snap，即使坐标很近。
2. exact logical identity 可直接合并；几何 proximity 只生成 proposal，不直接改点。
3. displacement budget 取 `min(local_h, incident edge length, source curvature/gap scale)`
   的无量纲 fraction，并受 feature clearance 约束。
4. 每个 proposal 必须记录 `accepted/rejected/refine/resample/rephase` 原因；零位移事件也
   有 lineage。
5. solver cell 直接保存 `SourceCellId[]`，任何 merge/split 都由 transaction 映射更新；
   禁止再用 centroid 对所有 source polygon 做回绑。

### 7.2 compatibility 决策顺序

```text
same exact VertexKey
  -> reuse
else feature/source compatible?
  -> no: preserve both; mark conflict
  -> yes: within dimensionless safe displacement and clearance?
       -> yes: canonicalize with deterministic priority
       -> no: choose local refine / grid-phase / source resample
```

确定性 priority 不能依赖注册时序：建议按
`physical feature > fixed layer > source-param vertex > exact grid vertex > mutable transition > smooth intersection`
排序，再按完整 typed key 排序。

### 7.3 空间索引

- grid vertex/grid edge：dyadic integer key 的 flat hash 或 ordered map；查询期望 `O(1)`，
  冻结输出时按 key 排序。
- source segments/features：静态 packed R-tree/BVH，build `O(S log S)`、query
  `O(log S + m)`；R1 初版可使用 deterministic uniform bins，但必须记录 worst-bin。
- mutable canonical vertices：以 finest-local-h bucket 的 multilevel spatial hash；每桶
  候选按 typed key 排序。禁止用一个全局 epsilon bucket。
- topology edge：`unordered_map<Undirected(VertexId,VertexId), EdgeRecord>`；输出前按
  ordered key compact。hash 只影响查找，不影响 ID/输出顺序。
- index update 只接受 transaction delta；失败 transaction 不污染主 index。每个 index
  带 revision，旧 handle 显式拒绝。

最坏情况：某一 spatial bucket/segment BVH leaf 过密时，不能静默退化全扫描；超过上限
触发局部再分桶/BVH rebuild，并把 `max_bucket_size`、`fallback_scan_count` 写入 profile。

### 7.4 half-edge-lite / edge-incidence topology

不引入完整 DCEL。最低结构：

```text
CellRecord { stableCellId, orderedVertexIds, edgeIds, sourceRefs, kind, localH }
EdgeRecord { stableEdgeId, v0, v1, owner, optional neighbour, patch, sourceRefs }
VertexIncidence { incidentEdges[], incidentCells[] }
Patch { cellSet, oneRingHalo, lockedBoundaryEdges, localVertex/Edge/Cell views }
TopologyDelta { add/update/remove records, lineageMap, affectedMetricEntities }
```

局部事务：

1. 从 conflict/bad metric edge 通过 incidence 取 cell set + 一圈 halo，`O(k*d)`；
2. 锁定 patch boundary 的有序 VertexId/EdgeId 与外侧 owner；
3. 在 shadow store 中构造候选，执行正面积、simple polygon、edge incidence 1/2、
   patch boundary equality、面积和 source/feature ownership 审计；
4. 只计算 patch cells、内部 edges 和跨 patch boundary edges 的质量；
5. 候选通过后原子提交 delta、递增 revision、更新 index/metric cache；
6. 每个 R1 子阶段结束执行一次 authoritative global rebuild/audit，验证增量结果与全局
   结果逐项相同。代码冻结前仍不以局部 audit 取代正式全局/OpenFOAM 验收。

目标复杂度：初始全局构建 `O((N+I) log(N+I))`、内存 `O(N+E+I+S)`；单次 patch
事务 `O(k log N + k log k)`，局部质量 `O(k+e_patch)`。这里的 `k` 必须由 profile
实际证明远小于 N；若 closure 扩展到全网格，必须显式报告而非仍称“局部”。

### 7.5 构造期 quality guidance

transition/template candidate 在生成 polygon 后、写入主 topology 前计算：

- 硬几何：positive area、simple loop、feature/source preservation、最短边三种 ratio；
- 局部 solver：non-orthogonality、skewness、face weight、volume ratio、min angle、
  hydraulic aspect；
- 选择 rank：先 hard violations 数，再 maximum/total normalized severity，再 cell count/
  template cost；不能先优化 cell count；
- neighbouring unbuilt cell 的 metric 用可证明上下界；无法给界时标为 `unknown` 并扩大
  patch，不能假定通过；
- unsafe intersection fallback 顺序：局部 quadtree refinement -> 合法 grid phase
  candidate -> 沿 source parameter wall/transition resampling -> constrained local polygon
  partition -> 明确失败。pure Cut-cell 仍只作为 H4 总 transaction 最后 fallback。

## 8. R1-A 到 R1-G 的迁移计划

每个子阶段独立 commit、真实 solver mesh、相关案例前后质量/时间、area/patch/hash、
无删 cell 证明。以下是计划，不是已完成范围。

### R1-A：lineage、feature classification、spatial index

- 新增 stable typed key、SourceRef/FeatureOwner、VertexRecord 和 index；
- 从 Q2-A registry 双写 shadow record，不改变坐标与 topology；
- solver source lineage 从已有 mapping 直接传播，同时保留旧 centroid scan 作 oracle；
- 验证两者逐 cell 相同后才能关闭 `O(C*S)` scan；
- 关键案例：narrow gap；预期 topology/hash 不变，source mapping time 显著下降。

### R1-B：构造期 canonicalization 与 compatibility

- 将 exact grid event 与 proximity proposal 统一到一个 decision API；
- sharp/concave/gap side compatibility matrix 单元测试；
- 把 Q2 partial superellipse 修复迁移为显式 source decision，不依赖散落 resampling；
- 关键案例：superellipse + sharp trailing edge；不允许移动原 wall feature。

### R1-C：不安全交点的 local refinement/resampling

- conflict 生成 typed request，而不是 throw 后由上层猜测；
- neighbour-driven 2:1 closure，记录 refined cell lineage；
- wall/transition resampling 只沿 source parameter，禁止跨 feature；
- grid phase 只作为 bounded deterministic candidates，必须比较面积/feature/quality；
- 关键案例：narrow gap + sharp trailing edge；保留当前失败最小坐标为回归。

### R1-D：patch-local topology 与增量质量

- 引入 edge-incidence store、patch boundary lock、TopologyDelta/revision；
- 局部 topology/quality 与全局 oracle 双算；
- profile global rebuild count、patch size、cache invalidation；
- 只有五案例逐项等价后，候选路径才停止 full rebuild/full quality。

### R1-E：quality-driven transition/template construction

- 固定 growth list 改为有界 template candidate generator；
- Q1 hard metrics 进入 commit 前 rank；
- 不降低合同；BoundaryLayer 仍保持独立 OBSERVED，除非另行批准阈值；
- 关键案例：concave L + sharp trailing edge。

### R1-F：隔离/删除 Q2-B 全局搜索

- 先保留 Q2-B 为 debug oracle，只在显式测试开关下运行；
- 对同一 failure patch 比较新 construction result 与 Q2-B 可行拓扑；
- 当两例及新最小回归均由构造期通过后，删除 production 调用和 32 次外层循环；
- 保留 Q2-B 文档、失败案例与历史 artifact provenance，不抹除历史。

### R1-G：冻结与正式验收

- 先由用户确认再运行完整 CTest；
- 五案例每例至少两次生成，比较稳定 IDs、report、solver hash；
- 五案例独立 reader + 真实 OpenFOAM v2606 `checkMesh -writeAllFields`；
- 报告 Q1 每项 hard count，目标不是只让 short-face 通过；
- 与本报告 baseline 同机同构建比较 wall/profile/peak RSS。若质量通过但时间或内存显著
  回归，不登记 R1 完成。

## 9. 现有代码的保留、替换、隔离与删除建议

### 保留

- `Domain2D/BoundaryRegion2D/Point2D/Segment2D/Polygon2D` 与外流物理语义；
- H1 sizing、H2 quadtree、H4 layer planning/local dropping/stepped envelope；
- Q0 provenance、Q1 typed contract、旧 solver hard safety；
- CLI、CM2D/VTK/OpenFOAM writer、独立 reader；
- Q2-A 的 dyadic grid-line identity、support event cache、active handle isolation；
- 所有面积、boundary classification、determinism 和 failure regressions。

### 替换

- `IntersectionRegistry2D::canonicalize` 的全 vertex scan -> feature-aware spatial index；
- `std::map<pair<double,double>>` 作为主身份 -> typed logical/source key；
- `qualityMetadataForSolver` centroid × all sources -> transaction lineage；
- candidate 全局 `buildGlobalTopology` -> patch-local delta + 阶段末 global oracle；
- candidate 全局 `evaluateSolverQuality2D` -> local metric cache + commit 后 affected edges；
- 固定 termination 候选后验试错 -> construction-time quality rank。

### 隔离

- legacy construction path：保留至 R1-D 作为 oracle，禁止默认生产使用；
- Q2-B repair：R1-F 前只作为 correctness oracle/debug；
- full global rebuild/audit：保留为 authoritative verification，不在逐候选热路径；
- constrained triangulation：仅 patch fallback，不把全网格改成 unstructured Delaunay。

### 最终删除条件

- 只有新路径在五案例、最小失败集、scale/determinism 和性能门全部通过后，才删除
  production Q2-B loop、centroid source scan 和旧 proximity canonicalization；
- 不删除历史文档、失败坐标、回归测试或 provenance；
- 不通过删除 solver cells、放宽 threshold、切 relaxed quality、吞掉 conflict 达标。

## 10. 第一批最小回归与性能实验

在任何大改前，先加/运行以下最小集合；开发期不运行完整 75 项与五案例 OpenFOAM。

### 单元/属性测试

1. 同一 wall segment × grid edge 从相邻 coarse/fine cell 请求，得到同 VertexId；
2. exact grid corner 与 arithmetic near-corner 的顺序无关；
3. 两条 narrow-gap side 距离小于 snap radius 仍得到不同 owner/VertexId；
4. sharp/concave feature 不被 smooth intersection 吸附；
5. source parameter resample 保持原 segment/curve、patch、面积与 feature endpoint；
6. patch replace 后 boundary edge sequence 完全相同，incidence 只能 1 或 2；
7. local metric 与 full Q1 对 affected entities 逐项相同；
8. transaction rollback 后 store/index/revision/hash 与前状态一致；
9. 整体 scale `1e-6/1/1e6` 后 decision、typed status、无量纲 worst 不变；
10. 随机注册顺序不影响 stable key、compact output 和 provenance（固定 seed）。

### 最小真实案例

- R1-A：narrow-gap，一次生成；证明 lineage 替代 scan 且网格 hash 不变；
- R1-B：superellipse + sharp-tail；证明 near-coincident 与 feature compatibility；
- R1-C/D：narrow-gap + sharp-tail；证明局部 fallback 和 patch transaction；
- R1-E：concave-L + sharp-tail；证明 quality-guided transition；
- 每个阶段只跑最相关 1--2 例，保留 solver mesh 和预览。

### microbench

- registry：10k/100k/1M vertices，均匀与最坏 cluster，记录 p50/p95 query、max bucket；
- source mapping：1k/10k/100k cells × sources，比较 lineage O(C) 与旧 O(C*S)；
- topology patch：固定全网 N=100k，patch k=4/16/64/256，比较 patch delta 与 full rebuild；
- quality：相同 k 测 local/full，验证结果等价；
- refinement closure：窄缝、尖角、长链三种 closure，记录 k/N 和迭代次数。

建议初始性能门 **[H]**：在质量不回归前提下，单 patch candidate 的 wall time 随 N
增长应接近 `log N` 而非线性；五案例端到端不得比 Q2-A 基线慢 20% 以上；sharp-tail
应从当前 Q2-B `4.91 s` 回到不高于 Q2-A `1.542 s` 的 1.2 倍。数值门需在 R1-A
基线工具补齐并由用户确认后固定。

## 11. 仍不知道、必须实验验证的事项

1. **[H]** 当前五案例和未来复杂输入中，patch closure 的 k/N 分布；若经常接近 N，
   patch-local 架构的速度收益会下降。
2. **[H]** sharp-tail 的 1832 candidates 中，时间分别花在 polygon union、global
   topology、global quality 的比例；当前 hybrid JSON 未输出 profile。
3. **[H]** narrow-gap 当前 42 s 的主耗时是否主要来自 H4 layer/termination 候选而非
   Q2-B；需要分阶段 timer 才能定量拆分。
4. **[H]** peak RSS 五案例基线；本轮系统权限阻止 `/usr/bin/time -lp` 提供该项。
5. **[H]** 2D source polyline 上最稳健的 adaptive predicate 组合；Gmsh 的 robust
   predicates 能指导，但不能直接复制 GPL 代码。
6. **[H]** local refinement 是否足以修复 immutable layer corner 邻近短面，还是必须
   grid phase 或 source resampling 联合；Q2 partial 的历史试验曾产生更坏 non-orth/
   face weight，不能假定 refinement 必然成功。
7. **[H]** 一个 patch 是否需要 constrained triangulation；先实现 simple polygon
   convex partition + incidence transaction，只有出现不能处理的最小失败例才升级。
8. **[H]** full DCEL 是否必要。目前没有证据；只有出现多 loop、多 region overlap、
   高频 face walk 且 half-edge-lite 无法保持不变量时再重新评估。
9. **[H]** grid phase 的候选数、可接受位移域和跨 patch 一致性；必须先写成有界、
   确定性的实验，不得成为全局连续优化。
10. **[H]** BoundaryLayer 六项 OBSERVED 的最终 hard/preferred 阈值；未获用户批准前
    不参与总 PASS，也不能借此放宽 ordinary contract。
11. **[H]** 新 topology store 的稳定 ID 在并行生成时如何保持调度无关；R1 初期应
    先固定单线程 deterministic commit order。
12. **[H]** GPL 参考实现思想与独立实现的法律边界不是技术测试能最终判定的；若未来
    要复制/链接 Gmsh/p4est/cfMesh/OpenFOAM 代码，需要单独许可证审查。

## 12. 决策请求

建议批准：**路线 B 主导 + 路线 C 的 half-edge-lite/edge-incidence、stable ID、delta
transaction 思想**，按 R1-A 开始；Q2-B 保持实验分支、不合并 main，直到 R1-F。

本报告不请求批准完整 DCEL、不请求复制 GPL 源码、不请求完整验收。若路线获批，
R1-A 的第一目标应是“lineage/index shadow mode，网格字节不变”，而不是立即改变
几何或声称 Q1 修复。
