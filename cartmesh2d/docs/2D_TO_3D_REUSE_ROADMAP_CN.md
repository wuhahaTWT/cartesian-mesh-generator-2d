# 2D 到 3D 的复用路线（H2）

## 目标边界

二维是算法试验场，不是把三维代码压成 `z=0`，也不应现在就把稳定的二维类型泛化为 `Dim` 模板。H2 只确定可迁移的概念、接口边界和证据门；不启动 CartMesh3D、boundary/prism layer、solution-driven AMR 或并行实现。

## 可迁移层次

| 二维验证对象 | 三维目标 | 推荐来源/等级 | 进入 3D 前的证据门 |
|---|---|---|---|
| `(level,ix,iy)` + Morton + flat leaves | quadrant/octant + SFC forest | p4est/p8est Level A 或 D/C | 许可证、1M octant、2:1、确定性、ghost/partition 原型 |
| flat face sort+sweep | 6-face octant adjacency | p4est思想，Level C | 2D 500k 的时间/RSS与邻接等价性；3D 内存估算 |
| `BoundarySegmentIndex2D` flat BVH | STL triangle BVH | CGAL思想 Level D；CGAL Level A需许可复核 | closest/intersection exactness、退化三角形、确定顺序、百万三角形 profile |
| full/cut/covered 分类 | 3D regular/cut/covered + multi-valued EB | AMReX Level A/D | sparse EB内存、multi-fragment、体积/面积闭合、AMR层间一致性 |
| H1 sizing policy | `targetLevel(x,y,z,bounds)` | Gmsh思想 Level D | 连续尺寸到 dyadic level 规则、composition、determinism |
| Cut polygon + global topology | Cut polyhedron + solver topology | 自有几何前端；后端另行选型 | 多三角共分割、共享面一致、正体积、真实 checkMesh |
| 小单元分析/合并 | 小体积分数 stabilization | 自有规则 + AMReX EB参考 | 守恒、连通、质量门，不允许 writer-side 掩盖 |
| 分阶段质量检查 | castellation/snap/layer pipeline | snappyHexMesh Level D | 各阶段可回滚、失败面显式、真实 OpenFOAM 质量 |
| 外部解场误差标记 | 求解驱动 AMR | Basilisk Level D | 守恒 remap、滞回、refine/coarsen、循环确定性 |

## 后端决策点

### 路线 A：p8est AMR 拓扑 + 自有 Cut-cell

优点是 octree、2:1、ghost、SFC partition 与 MPI 已成熟；自有几何/求解器拓扑仍可保持产品语义。缺点是 GPL 许可待评审，且 p8est 不解决 triangle cutting 和 polyhedral quality。只有许可与接口原型同时通过才进入 Level A，否则按论文/思想独立实现 Level C/D。

### 路线 B：AMReX AMR/EB backend

BSD-3-Clause 更适合直接依赖评估，并具备 block AMR、regular/cut/covered 与 sparse EB 数据。风险是 AMReX 的 patch/Fab 与本项目逐 leaf OpenFOAM mesh 的数据模型不同。先做只读转换原型：自有 boundary -> AMReX geometry -> 导回 cell flags/volume fractions -> 与自有面积/体积对照；不一开始就把产品架构绑死。

### 路线 C：保持轻量自有串行核心

适合中小规模、确定性和研究可控性。若 2D flat leaves/邻接/BVH 在 1M 仍合理，可保留为 reference backend；3D 规模、并行与 ghost 需求出现后，再用 A/B 后端替换 AMR 基础设施。reference backend 也用于交叉验证外部库。

## 数据布局演进原则

1. leaf/octant 核心保持 trivially movable 的整数坐标与 level；物理 bounds 可按需重建，是否删除缓存必须由 P1 benchmark 决定。
2. regular geometry 隐式保存；只有 cut/multi-fragment 单元进入 sparse geometry pool。H2 先测 `sizeof(QuadtreeLeaf2D)` 与 CutCell 内存占比，再决定 P3。
3. 邻接、拓扑与导出使用连续数组和稳定排序；CSR 只在 `vector`-per-cell 已被证明是内存热点时进入。
4. 2D 与 3D 共享的是行为契约和测试语义，不强制共享模板实现：分类、守恒、2:1、deterministic IDs、失败可见、外部检查。

## 未来阶段顺序与停线条件

1. H2：完成 100k 全流水线，目标 500k，记录 P0 前后时间/RSS；H1 correctness 不回退。
2. 3D backend spike：仅在用户批准后，分别测 p8est 与 AMReX 的最小原型；不接 Cut-cell 产品路径。
3. 3D geometry：triangle BVH、共同分割、shared-face coverage 与最小失败案例；任何负体积/非闭合立即停线。
4. solver mesh：独立 reader + OpenFOAM `checkMesh`；内部测试不能代替外部验收。
5. boundary/prism layer：最后单独立项，借鉴 snappy 的阶段机和失败回退；不得把简单法向挤出称为成熟 boundary layer。

P0 若无法在不改变邻接集合、2:1 与确定性输出的前提下显著改善 500k 级成本，则停止继续微优化，依据 profile 决定 P1/P2/P3，而不是顺序性重构全部容器。
