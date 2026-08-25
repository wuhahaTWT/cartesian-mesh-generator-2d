# H2 开源复用研究与审计

日期：2026-08-25  
范围：Stage 2D-H2/H3；只研究大网格数据结构、空间索引、AMR/EB、尺寸场、solver topology transaction 与未来边界层参考。本文不是法律意见；本仓库未发现顶层许可证文件，任何外部代码的链接、复制或分发均须先完成许可证复核。

## 1. 复用等级

- Level A：直接依赖并链接；仅在接口、维护性和许可证均明确时采用。
- Level B：移植独立模块；保留原版权、许可证、commit、修改记录和第三方声明。
- Level C：阅读论文与源码后独立重写项目所需的精简算法；记录思想来源，不复制表达性源码。
- Level D：只借鉴架构思想或数据组织，不形成源码派生关系。

H2 实测过邻接 `sort+sweep` Level C 原型，但它在当前二维 lattice coordinate
分布下慢于 production 的 coordinate-bucket `std::map`，因此没有采用。稀疏 Cut-cell
存储、尺寸场组合、分阶段质量回退和 solution-driven AMR 仅为 Level D。H2 没有新增
外部依赖，也没有复制第三方源码；该二维结果不能直接外推到三维 octree。

## 2. 项目审计

### 2.1 p4est / p8est

- Project / repository：[p4est](https://github.com/cburstedde/p4est)，检查 `master` commit `5a891d2c69833eb05fd8299fc0f19243587d685a`。
- Paper：Burstedde, Wilcox, Ghattas, [p4est: Scalable Algorithms for Parallel Adaptive Mesh Refinement on Forests of Octrees](https://p4est.github.io/papers/BursteddeWilcoxGhattas11.pdf)。
- License：p4est 源文件为 GPL-2.0-or-later；其基础库 libsc 为 LGPL-2.1-or-later。直接链接或移植前均为 **LICENSE REVIEW REQUIRED**。
- 实读模块：`src/p4est.h`、`p4est_bits.h`、`p4est_search.h`、`p4est_balance.c`；关注 quadrant 整数坐标/level、连续 `sc_array_t`、Morton 顺序、祖先/重叠判断、递归搜索、balance 与 ghost/partition 接口。
- 与当前二维的相似处：整数 `(level, ix, iy)`、Morton key、连续 leaf vector、批量 refine、2:1 balance、确定性排序。
- 明显差距：当前每次邻接调用重建四棵 `std::map<coordinate, vector<face>>`；balance 多轮重复生成全邻接；leaf 同时保存可重建的 AABB 与分类；没有 coarsen、forest connectivity、ghost、partition、MPI。
- 2D 复用：H2 实现并 benchmark 了连续 `FaceRecord` 排序扫描 Level C 原型，但在
  486k leaves 下最好仍慢于 coordinate buckets，production 保留 `std::map`。三维需
  单独评估 p8est/Morton/octree 邻接，不能沿用二维结论。
- 3D 复用：p8est 对大规模八叉树、2:1、ghost、SFC partition 很有价值，优先做隔离后端原型，而不是把当前二维类型过度模板化。若产品许可证不兼容 GPL，则只可 Level D/C，需法务复核。
- 风险：p8est 解决 AMR 拓扑与并行分区，不自动提供 STL Cut-cell、多面体求解器拓扑或 prism layer。

### 2.2 AMReX AMR / Embedded Boundary

- Project / repository：[AMReX](https://github.com/AMReX-Codes/amrex)，检查 `development` commit `28ba7777ca50c76287aa119252e3f8275e4638f7`；可复现实验时优先固定正式 release。
- Papers：[JOSS framework paper](https://doi.org/10.21105/joss.01370)；[block-structured AMR paper](https://arxiv.org/abs/2009.12009)。
- License：[官方 LICENSE](https://raw.githubusercontent.com/AMReX-Codes/amrex/development/LICENSE) 为 BSD-3-Clause 风格；仍需保留 copyright/条件/免责声明。
- 实读模块与文档：[Embedded Boundary chapter](https://amrex-codes.github.io/amrex/docs_html/EB_Chapter.html)、`Src/EB/AMReX_MultiCutFab.H`、`AMReX_EBCellFlag.H`。`EBCellFlag` 区分 regular/covered/single-valued/multi-valued 并编码连通性；`MultiCutFab` 只为含 cut cell 的 Fab 分配昂贵数据。
- 2D 复用：如果 H2 benchmark 证明 CutCell 对每个 regular leaf 的动态成员占用成为热点，采用 `cellKind[] + regular-cell implicit geometry + sparse cut-record[]`，Level D；本阶段不预先改。
- 3D 复用：可评估 `CartMesh3D geometry frontend + AMReX AMR/EB backend`，Level A 原型。要先验证多值 EB、外部多面体导出、确定性 ID 与 OpenFOAM 接口，不能假定 AMReX 可直接替代网格生成器。
- 风险：block-structured hierarchy 与本项目逐 leaf solver topology 不同；GPU/MPI/Fortran 生态增加构建复杂度。

### 2.3 OpenFOAM snappyHexMesh

- Project / repository：[OpenFOAM-dev](https://github.com/OpenFOAM/OpenFOAM-dev)，检查 `master` commit `51ac9f2729626eb7838db77dc7ba71ba18e0bcb5`。
- Documentation：[snappyHexMesh user guide](https://www.openfoam.com/documentation/user-guide/4-mesh-generation-and-conversion/4.4-mesh-generation-with-the-snappyhexmesh-utility)、[layer controls](https://doc.openfoam.com/2212/tools/pre-processing/mesh/generation/snappyhexmesh/layers/)。
- License：[COPYING](https://github.com/OpenFOAM/OpenFOAM-dev/blob/master/COPYING) 与源码头为 GPL-3.0-or-later；复制或链接为 **LICENSE REVIEW REQUIRED**。
- 实读源码：`src/mesh/snappyHexMesh/snappyHexMeshDriver/snappyLayerDriver.C` 及其调用的 `addPatchCellLayer`、mesh mover、质量检查与 extrusion unmark 路径。
- 关键架构：castellation、snap、layer 是独立阶段；layer 通过表面法向/平滑、网格收缩、逐层挤出、质量验证、失败面取消或减薄并重试实现；`featureAngle`、`minThickness`、层数与 expansion 控制必须和 termination/corner 策略共同工作。
- 2D/3D 复用：未来二维 body-fitted layer 与三维 prism layer只借鉴阶段机、失败回退和质量门，Level D；H2 明确不实现 layer。
- 风险：其拓扑操作紧耦合 OpenFOAM polyMesh；直接移植既有许可证问题，也会带入大量框架依赖。

### 2.4 CGAL AABB Tree

- Project / repository：[CGAL](https://github.com/CGAL/cgal)，检查 HEAD `ff859cc761a68caf612822e2aa8ba6ca95deb0bc`；接口以 [AABB Tree reference](https://doc.cgal.org/latest/AABB_tree/group__PkgAABBTreeRef.html) 为准。
- License：[CGAL licensing](https://www.cgal.org/license.html) 指出 AABB Tree package 使用 GPL（另有商业许可），直接采用为 **LICENSE REVIEW REQUIRED**。
- 相关能力：静态 primitive array、层次包围盒、intersection/closest-point/distance query、遍历剪枝；部分距离查询可用额外搜索结构加速。
- 当前 `BoundarySegmentIndex2D`：已是连续 `segments/order/nodes` 的静态 BVH，按最长轴中心稳定二分，stack 查询，distance 用 lower-bound priority queue。它与成熟 AABB tree 的核心布局一致，但缺少 SAH/增量更新/统一 primitive trait，也会在每次查询分配 stack/queue。
- 2D 复用：继续维护轻量自有 BVH，Level D；只有 profile 指向该处才优化复用 scratch buffer 或 node layout。
- 3D 复用：segment primitive 可替换为 triangle primitive，保留 flat node/index/query-pruning 思路；若直接用 CGAL，须评估 GPL/商业许可、robust kernel 开销与确定性。

### 2.5 Gmsh size fields

- Project / repository：[Gmsh mirror](https://github.com/live-clones/gmsh)，检查 HEAD `827f12445e204afb85137a6eb843783d0495542d`；[官方手册](https://gmsh.info/doc/texinfo/gmsh.html)。
- License：[官方说明](https://gmsh.info/) 为 GPL-2.0-or-later，并提供 linking exception/商业选项；复制 `Field.cpp` 仍为 **LICENSE REVIEW REQUIRED**。
- 实读/核对：`src/mesh/Field.cpp` 的 field registry 与 Distance、Threshold、Box、Min/Max；手册中的 background mesh/field composition。
- H1 对应：`minimumLevel` 是全局 floor；distance band 对应离散化的 Distance+Threshold；refine box 对应 Box；多个场取 `max(targetLevel)` 等价于取最小目标物理尺寸。
- 2D/3D 路线：长期将 policy 的判定集中成只读 `targetLevel(x,y[,z], cellBounds)` 组合接口，Level D；H2 不增加更多 field 类型，也不把维度抽成复杂模板体系。
- 风险：level 是 dyadic 离散量，不能声称等价于 Gmsh 连续 size field；边界 AABB 距离和点距离语义也需明确。

### 2.6 Basilisk

- Project / source：[Basilisk](https://basilisk.fr/)，检查 2026-08-25 在线源码快照；站点没有稳定 Git commit，因此版本标识不足，**VERSION REVIEW REQUIRED**。
- Source：[tree-common.h / adapt_wavelet](https://basilisk.fr/src/grid/tree-common.h)；[tutorial](https://basilisk.fr/Tutorial)。
- License：[src/COPYING](https://basilisk.fr/src/COPYING) 为 GPL-3.0，复制/链接为 **LICENSE REVIEW REQUIRED**。
- 算法要点：对选定求解场做 restriction/prolongation，按用户容差估计表示误差；超阈值 refine，足够小则 coarsen，并维护依赖字段与边界更新。
- 2D/3D 复用：`OpenFOAM field -> error indicator -> target-level marks -> balance -> regenerate/remap` 仅作 Level D 设计参考。未来必须显式处理场读取、守恒映射、滞回、最小/最大 level、确定性和网格运动；H2 不实现 solution-driven AMR。

## 3. H2 实施决策与许可结论边界

1. P0 曾独立实现 flat record + deterministic sort + coordinate group sweep；benchmark
   证明它在当前二维分布下回退，因此未进入 production。最终只复用 balance 同一状态的
   邻接结果；没有复制 p4est/Gmsh/OpenFOAM/CGAL/Basilisk 源码。
2. P1 避免 Quadtree leaf/vector 的重复复制与排序，必须由 benchmark 证明后再改。
3. P2 topology CSR、P3 sparse cut-cell、P4 spatial index 只在前序 profile 显示为瓶颈时进入，避免为未来 3D 提前过度抽象。
4. AMReX 是本轮唯一具备宽松许可证、值得做未来直接依赖原型的候选；p8est 的算法/并行价值最高，但 GPL 兼容性必须先评审。
5. 本文不认定任何许可证与未来商业/闭源目标兼容；所有 Level A/B 都要在项目许可证明确后再次审核。

### 3.1 H3 补充：OpenFOAM 批量 topology change

- 参考：[OpenFOAM Foundation `polyTopoChange` API](https://cpp.openfoam.org/dev/polyTopoChange_8H_source.html)、
  [`refineMesh` 调用路径](https://cpp.openfoam.org/dev/refineMesh_8C_source.html) 和
  [`polyTopoChange::changeMesh`](https://cpp.openfoam.org/dev/polyTopoChange_8C_source.html)。
- 源码架构先把 point/face/cell 的 add、modify、remove actions 记录到 mesh changer，随后
  一次 `changeMesh()` 应用并返回 topology map，再由 mesh/object 执行映射更新。
- H3 复用等级：Level D/C。只借鉴“收集局部 actions -> 确定性批量 apply -> authoritative
  validation”的 transaction 思想；不复制 GPL-3.0-or-later 源码，不链接 OpenFOAM。
- 当前边界：H3 第一版仍通过一次全局 topology rebuild 提交批次，避免直接维护局部
  vertex/edge owner-neighbour 和 ID remap。只有 profile 证明 accepted-batch rebuild 仍是
  主瓶颈，才单独评估真正 incremental topology。
