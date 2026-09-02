# R2 参考项目审计

日期：2026-09-01
范围：R2（加密鲁棒性）为止实际读过的外部与内部参考。格式沿用
`docs/OPEN_SOURCE_REUSE_RESEARCH_CN.md` 的复用等级。本文不是法律意见。

## 0. 复用等级与本轮决定

沿用既有定义：Level A 直接依赖并链接；Level B 移植独立模块并保留版权与修改记录；
Level C 读论文与源码后独立重写项目所需的精简算法，记录思想来源，不复制表达性源码；
Level D 只借鉴架构思想。

**R2 全程按 Level C/D 执行**：不复制表达性源码，不链接任何外部库，不产生派生关系。
本轮**没有新增任何外部依赖**，也没有从任何参考树复制源码。

## 1. 本地参考树的许可证（逐一在本地核对）

| 项目 | 本地路径 | 许可证文件 | 结论 |
|---|---|---|---|
| OpenFOAM v2606 | `mesh-reference-sources/OpenFOAM-v2606-source` | `LICENSE.md` | GPL-3.0-or-later。**只读不抄** |
| cfMesh 1.2.0 | `mesh-reference-sources/cfMesh-1.2.0` | 无顶层文件；源码头声明 GPL v3+ | GPL-3.0-or-later。**只读不抄** |
| gmsh 4.15.2 | `mesh-reference-sources/gmsh-4.15.2-source` | `LICENSE.txt` | GPL-2.0-or-later，带 Netgen/METIS 等链接例外。**只读不抄** |
| p4est 2.8.7 | `mesh-reference-sources/p4est-2.8.7` | `COPYING` | GPL-2.0（v2 正文）。**只读不抄** |
| AMReX (development) | `mesh-reference-sources/amrex-development` | `LICENSE` | BSD-3-Clause 风格（LBNL）。若要 Level B，须保留版权、条件与免责声明 |

把 GPL 源码复制进本仓库会使本仓库成为 GPL 派生作品；本轮明确不这样做。AMReX 是
唯一许可宽松、值得将来做 Level B 的候选，但 R2 也没有用到它的代码。

## 2. 已消费的借鉴（W0）

### 2.1 OpenFOAM 质量指标定义与阈值 — Level D（交叉确认）

- 实读：`src/OpenFOAM/meshes/polyMesh/polyMeshCheck/polyMeshTools.C:175`
  （`faceWeights`：`min(dNei,dOwn)/(dNei+dOwn+VSMALL)`，其中
  `dOwn = |Sf·(Cf-Cown)|`）与 `:232`（`volRatio`：
  `min(Vown,Vnei)/(max(Vown,Vnei)+VSMALL)`）。
- 实读：`etc/caseDicts/meshQualityDict` 全文，取到 `maxNonOrtho 65`、
  `maxInternalSkewness 4`、`maxBoundarySkewness 20`、`minFaceWeight 0.05`、
  `minVolRatio 0.01`、`minDeterminant 0.001`、`minTwist 0.02`。
- 本轮用途：**只做交叉确认**。本仓库 `SolverQualityPolicy2D` 的
  `minFaceWeight 0.05` / `minVolumeRatio 0.01` / `maxInternalSkewness 4` 与之一致，
  `maxNonOrthogonalityDeg 70` 略宽，`maxBoundarySkewness 4` 则比 OpenFOAM 的 `20`
  **严格 5 倍**，Q1 contract 更严（65° / 3 / 0.10 / 0.05）。这个差别有实际后果：
  W0 记录的 circle level 10 失败是 boundary skewness `8.5955`，OpenFOAM 会接受它，
  是本仓库更严格的门拒绝的。相关判断见 `docs/R2_REFINEMENT_ROBUSTNESS_CN.md` 第 4.1 节。
  `tools/verification/refinement_ladder.py` 里的三个硬限是**镜像**自本仓库的
  C++ header，注释里同时标注了 OpenFOAM 的对应值；C++ 侧仍是唯一权威。
- 没有复制任何公式实现：本仓库的
  `evaluateSolverInternalFaceMetrics2D` 早于本轮存在，本轮未改动。

### 2.2 三维仓库的同现象归因 — Level D（仅设计参考）

- 实读：`cartesian-mesh-generator/STAGE6_8_CFD_QUALITY_PLAN_CN.md` 的「当前基线」
  与 6.8.1 节。
- 取到两条结论：(a) full cube 876 个、fixed icosahedron 744 个 non-orthogonality
  **全部位于 Cut ↔ Layer 接口**，采用的解法是 fragment-matched layer column
  （沿原始 STL 三角面重心插值细分，几何不移动）；(b)「把 layer facet 外侧全部
  Cut-cell 粗暴聚并成一个大 transition cell」的探针**已被否决**
  （non-orthogonality 恶化到约 111° / 129°）。
- 边界：`AGENTS.md` §2 禁止二维引入三维头文件或链接三维 library。本文只记录
  设计结论，W2 将在二维独立实现，不复用三维代码。

## 3. 已读、留给 W1–W4 的借鉴

| 需求 | 实读位置与深度 | 取什么 | 等级 |
|---|---|---|---|
| 按局部长度比例塌缩退化边 | `OpenFOAM/src/dynamicMesh/polyTopoChange/polyTopoChange/edgeCollapser.{C,H}`、`applications/utilities/mesh/advanced/collapseEdges`（**仅确认存在，尚未逐行读**） | 「退化判据用局部边长比例而非绝对值」这一原则。**W1 已消费**：`IntersectionRegistryPolicy2D::gridCornerWeldFractionOfLocalH` 是 h 的比例而不是绝对 epsilon。同时 W1 实测确认了 OpenFOAM 把 snap 与 collapse 分成两个阶段的理由——在本项目里，逐 face 塌缩会让共享 spur 的两个面不一致，逐 leaf 焦合会让邻居 leaf 不一致，只有全局 canonical vertex 空间可行 | C/D |
| 相对层厚与最小厚度语义 | `OpenFOAM/src/mesh/snappyHexMesh/snappyHexMeshDriver/layerParameters/layerParameters.H`（读成员与访问器清单：`relativeSizes_`、`firstLayerThickness_`、`finalLayerThickness_`、`expansionRatio_`、`minThickness_`、`featureAngle_`、`nGrow_`、`maxFaceThicknessRatio_`、`nBufferCellsNoExtrude_`、`nLayerIter_`、`nRelaxedIter_`） | 参数语义：厚度可相对局部尺寸给定；低于 `minThickness` 就不挤出；失败后有 relaxed 迭代 | D |
| 窄缝层厚上限 | `OpenFOAM/src/mesh/snappyHexMesh/externalDisplacementMeshMover/medialAxisMeshMover.C:1378`（`maxThicknessToMedialRatio`）、`:151`（`minMedialAxisAngle`）、`:182`（`nMedialAxisIter`）（**读参数取用点，未逐行读求解过程**） | 「总层厚不超过到 medial axis 距离的固定比例」这一思想。二维将用「到最近非相邻 wall segment 的距离」作代理，且**不移动网格** | C/D |
| 质量失败后的退层重试阶段机 | `OpenFOAM/src/mesh/snappyHexMesh/snappyHexMeshDriver/snappyLayerDriver.C`（**仅确认存在**；H2 审计曾读过其 extrusion unmark 路径） | 「挤出 → 质量检查 → unmark/减薄 → 重试」的分阶段回退结构 | D |
| 层可独立于输入表面离散度再细分 | `cfMesh-1.2.0/meshLibrary/utilities/boundaryLayers/refineBoundaryLayers/`（**仅列出文件与规模**：`refineBoundaryLayersCells.C` 1921 行、`refineBoundaryLayersFaces.C` 1291 行等） | 存在性证据：成熟工具会在生成后再细分层，切向与法向都可以 | D |
| 退化拓扑清理项分类 | `cfMesh-1.2.0/meshLibrary/utilities/smoothers/topology/topologicalCleaner/` 及同级 `checkNonMappableCellConnections`、`checkCellConnectionsOverFaces`、`checkIrregularSurfaceConnections`、`checkBoundaryFacesSharingTwoEdges`（**仅列出目录**） | 清理项的**分类**。本轮不引入 smoothing / untangling（会移动顶点，与「壁面点必须落在原输入 segment 上」冲突） | D |
| 二维 Cartesian mesher 的存在 | `cfMesh-1.2.0/meshLibrary/cartesian2DMesh/cartesian2DMeshGenerator`（**仅列出目录**） | 仅作对照，未取用 | D |
| 小体积分数阈值概念 | `amrex-development/Src/EB/AMReX_EB2_2D_C.cpp:164,187`（`small_volfrac`：`vfrac < small_volfrac` 时把单元置为 covered）、`Src/EB/AMReX_EB_Redistribution.H` 头部说明 | 「用体积分数阈值定义小单元」这一概念。**其解法不适用**：AMReX 直接 cover 掉小单元再靠 flux redistribution 补守恒，那是求解器侧策略；本产品导出 polyMesh 给通用求解器，只能**合并不能删除**（`AGENTS.md` §4.10） | D |

## 4. 本轮不动的参考

- `p4est`：H2 已实测本二维分布下 coordinate-bucket 优于 sort+sweep，且 W0 的阶段
  归因显示当前热点在 solver topology 阶段，不在邻接构造。本轮不重新评估。
- `gmsh`：size field 的对应关系已在 H2 审计登记，R2 没有新增 field 类型。

## 5. 后续义务

W1–W4 每消费一条上表条目，必须回到本文把「仅确认存在 / 仅列出目录」升级为实读
记录，并写明最终实现与参考的差异。任何要升到 Level A/B 的候选都必须在本仓库
许可证明确之后再次评审。
