# CartMesh2D H4：共形 hybrid mesh solver-ready 收口

日期：2026-08-27  
状态：circle / superellipse 原生二维 hybrid mesh 已通过既有 solver-quality 门禁

## 1. 实现边界

本轮保留 H4-1 固定 boundary-layer strip 与 H4-2 outer-envelope 共形接口，不实现
local layer dropping、复杂 termination、Delaunay、overset 或三维功能，也未修改三维核心。

## 2. Solver-ready 路径

1. H4-1 的 wall、层数、layer quad 和 outer envelope 全部保持不变。
2. outer envelope 外增加三圈可修复的渐进 transition fan；第一圈与每条 envelope edge
   一一共享，后续圈逐级加密切向分辨率，避免长 layer face 直接连接许多细小 Cut-cell face。
3. transition 宽度使用确定性、事务式候选修复：每个候选都完整运行 remainder Cut-cell、
   small-cell analysis、agglomeration、统一拓扑、`buildSolverTopology2D` 和原 solver-quality
   阈值；仅返回首个真实 PASS 的候选，失败候选不会修改 H4-1 输入。
4. solver repair 新增 immutable/preserve 约束：固定 layer cell 不允许合并或重分区；
   transition/remainder 可以参与既有 quality agglomeration 和 local repartition。
5. `BoundaryLayer`、`Transition`、`RemainderCut`、`RemainderCartesian` 使用明确的 hybrid
   source 类型。layer/transition 没有伪造 quadtree key，质量报告只从真实 remainder
   Cut-cell 统计 cut 数量与 level distribution。
6. H4-1 strip 审计显式要求每个 quad 严格凸，并包含制造的正面积凹四边形回归。

## 3. 验收指标

| 样例 | layer | transition | remainder cut | Cartesian | base/solver cells | 接口边 | 面积误差 | max non-orth | min face weight | min volume ratio |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| shifted 32-segment circle，4 层 | 128 | 128 | 188 | 304 | 700 / 728 | 32 | `-1.07e-14` | `55.3968°` | `0.0872106` | `0.0270379` |
| `superellipse_24`，3 层 | 72 | 96 | 200 | 336 | 680 / 703 | 24 | `5.15e-14` | `69.8923°` | `0.0937449` | `0.0104554` |

两例均为：conformal interface PASS、area conservation PASS、topology PASS、
solver-quality PASS。阈值仍为 non-orthogonality `70°`、face weight `0.05`、volume ratio
`0.01`，没有放宽或隐藏失败。

独立 hybrid reader 验证两例均无 overlap、non-manifold edge 或非正面积；OpenFOAM
独立 reader 验证 solver topology 写出的 case 闭合、owner/neighbour 有效且最小体积为正：

- circle：728 cells，3212 faces，最小体积 `4.14567e-05`；
- superellipse：703 cells，3100 faces，最小体积 `2.93687e-05`。

完整构建与回归：`98/98 PASS`；H4 专项：`19/19 PASS`。circle 的 base CM2D、solver
CM2D 和 JSON 重复输出均字节一致。

## 4. H3 scalability

Release、单进程/单线程、Apple M1 MacBook Air 8 GB，NACA0012 约 100k 路径：

| 指标 | 修改前 | 修改后 | 变化 |
|---|---:|---:|---:|
| leaves | 101734 | 101734 | 0 |
| solver cells / faces | 102218 / 204678 | 102218 / 204678 | 0 |
| solver topology | 6.143073 s | 6.049187 s | -1.5% |
| total internal | 10.167817 s | 10.028113 s | -1.4% |
| external wall clock | 10.21 s | 10.07 s | -1.4% |
| peak RSS | 381008 KiB | 373328 KiB | -2.0% |

任意方向 edge splitting 使用窄坐标范围候选搜索；本轮未观察到 H3 scalability 退化。

## 5. 产物与 OpenFOAM 状态

`cartmesh2d/artifacts/h4_2/` 保存 circle 与 superellipse 的：

- 分类可视化 `.hybrid.vtk`、base `.hybrid.cm2d`、JSON/quality 报告；
- solver-ready `.hybrid.solver.vtk` 与 `.hybrid.solver.cm2d`；
- 独立检查 JSON 与 `.svg`/`.png` 可视化。

CLI 已用 `solverTopology` 写出两套真实 OpenFOAM case，独立 reader PASS。本机
`checkMesh` 不在 PATH，因此真实 `checkMesh` 状态为 **UNAVAILABLE / NOT RUN**，没有伪报。

## 6. 已知边界

- 当前仍仅支持 H4-1 已支持的闭合外流、固定层数 strip；
- nested wall、local dropping 和复杂 termination 继续 fail closed；
- 已验收 circle/superellipse 不存在已知非共形接口或 solver-quality 失败。
