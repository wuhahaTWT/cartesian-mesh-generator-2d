# H4 前置源码研究：二维边界层与未来三维复用路线

> 状态：研究结论，不是 H4 实施许可，也不是生产代码设计冻结。  
> 基线：`agent/native-2d-robustness`，`96616c68718eb2fd825544831a1cac6ec563e03c`。  
> 研究日期：2026-08-27。  
> 当前工程：`/Users/Zhuanz/Desktop/cartesian-mesh-generator-2d`。  
> 第三方只读源码实际位置：`/Users/Zhuanz/Desktop/mesh-reference-sources`。

## 0. 范围、事实边界与许可边界

本轮只做“读源码—理解算法—对照 CartMesh2D—形成 H4 技术路线”。没有修改生产代码，没有开始 H4，没有把第三方代码复制进仓库，也没有实现 p4est/p8est、AMReX、OpenFOAM 或三维功能。

用户提示中把参考源码描述为放在“网络生成器”的 `sources/` 下；现场检查后，它们实际位于独立目录 `mesh-reference-sources/`。本报告按实际只读路径研究，不混入旧三维工程 `/Users/Zhuanz/Desktop/网络生成器`。

本地 `papers/` 目录当前为空，因此论文部分由公开的一手论文、项目官方文档与上游源码补齐；本报告会把“本地源码证据”和“在线论文证据”分开陈述。

| 参考项目 | 本地快照 | 版本/提交可追溯性 | 许可判断 | 本项目允许的复用方式 |
|---|---|---|---|---|
| cfMesh | `cfMesh-1.2.0` | 目录名为 1.2.0；无 `.git` 提交可查 | 源文件头为 GPL v3 或更高版本 | 只借鉴算法思想；不得复制实现 |
| OpenFOAM | `OpenFOAM-v2606-source` | v2606 源码包；无 `.git` 提交可查 | GPL | 只借鉴调用链、事务式质量控制思想 |
| p4est/p8est | `p4est-2.8.7` | 2.8.7；无 `.git` 提交可查 | GPL v2 或更高版本；未发现通用链接例外 | 当前不能作为“无条件可直接链接”的后端；须先决定产品许可并做法律审查 |
| AMReX | `amrex-development` | development 快照；无 `.git` 提交可查 | 三条款式宽松许可证 | 可借鉴数据布局；直接依赖仍须版本锁定、归属声明和架构评估 |
| Gmsh | `gmsh-4.15.2-source` | 4.15.2；无 `.git` 提交可查 | GPL v2 或更高版本，例外仅覆盖列明依赖 | 只借鉴 sizing field 抽象，不复制实现 |

CartMesh2D 当前未找到明确的根许可证文件。在许可证未定之前，即使技术上适配 p8est，也不应把 GPL 库列为默认“可直接复用”组件。许可证兼容性是未来三维后端选型的硬门槛，不是发布前再补的文书问题。

## 1. cfMesh 边界层调用链

### 1.1 入口与阶段顺序

三维 Cartesian 生成器入口：

`meshLibrary/cartesianMesh/cartesianMeshGenerator/cartesianMeshGenerator.C`

主流程为：

1. `templateGeneration/createCartesianMesh`：生成并细化 Cartesian 背景网格；
2. `surfaceTopology/surfacePreparation`：准备表面拓扑；
3. `surfaceProjection/mapMeshToSurface`：把网格边界映射到输入表面；
4. `patchAssignment/extractPatches`：分配 patch；
5. `edgeExtraction/mapEdgesAndCorners` 与表面优化：映射特征边角；
6. `boundaryLayerGeneration/generateBoundaryLayers`：调用 `boundaryLayers::addLayerForAllPatches()`；
7. `meshOptimisation/optimiseFinalMesh`：体网格、表面和边界层优化与解缠；
8. `boundaryLayerRefinement/refBoundaryLayers`：按层数和增长率细分已有包裹层；
9. 重编号并写出边界。

二维生成器 `meshLibrary/cartesian2DMesh/cartesian2DMeshGenerator/cartesian2DMeshGenerator.C` 的顺序相同，但调用 `activate2DMode()`，排除 empty patch 后仍生成一个薄三维网格。它不是 CartMesh2D 这种原生二维多边形网格，因此只能迁移思想，不能把其“2D 模式”直接当作本项目实现。

### 1.2 第一层不是逐层前沿推进，而是“先造包裹单元”

核心目录：

- `meshLibrary/utilities/boundaryLayers/boundaryLayers/`
- `boundaryLayersCreateVertices.C`
- `boundaryLayersFacesAndCells.C`
- `boundaryLayerCells.C`
- `boundaryLayersCheckTopologyOfBndFaces.C`

`boundaryLayers::addLayerForAllPatches()` 先 `createNewVertices()`，再 `createLayerCells()`。其关键思想是：

- 使用 `meshSurfaceEngine.pointNormals()` 得到面积加权点法向；
- patch 交界处把法向投影到相邻 patch 的切平面，角点则沿特征交线处理；
- 初始推进尺度取局部邻点距离量级，并受边界边距离约束；
- 新旧顶点坐标交换，使原表面点继续留在真实边界，内侧点形成一层包裹单元；
- 对共享边、角和层终止位置建立专门的面与单元拓扑。

这是一条“拓扑先行”的路线：先获得一层合法、可识别的 layer cell，再沿其“hair edge”细分，而不是每一层都重新推进整个前沿。

### 1.3 多层细分、法向与厚度控制

主要文件：

- `meshLibrary/utilities/boundaryLayers/refineBoundaryLayers/refineBoundaryLayers.C`
- `refineBoundaryLayersFunctions.C`
- `refineBoundaryLayersFaces.C`
- `refineBoundaryLayersCells.C`
- `meshLibrary/utilities/optimisation/boundaryLayerOptimisation/`

`refineBoundaryLayers` 读取全局或逐 patch 的 `nLayers`、`thicknessRatio`、`maxFirstLayerThickness`、排除 patch 和是否允许不连续等设置，随后执行：

`analyseLayers → generateNewVertices → generateNewFaces → generateNewCells`

它识别从壁面到层外包络的 hair/split edges，并在已有 hair edge 上按几何级数放置新点。若首层厚度为 \(h_1\)、增长率为 \(r\)、层数为 \(n\)，总厚度是：

\[
H = h_1\frac{r^n-1}{r-1},\quad r\ne1
\]

优化阶段会重新计算和光顺 hair 法向，根据曲率、局部特征尺度、射线与面相交、自相交风险限制厚度，并对负体积或纠缠区域做平滑/解缠。cfMesh 用户指南也明确说明：生成过程本身不保证永远没有负体积，后续质量优化是必要防线，而不是可以省略的美化步骤。

### 1.4 凹角、窄缝与终止

`findPatchesToBeTreatedTogether()` 会分析多 patch 交会和凸/凹共享边，把需要共同处理的 patch 编组，必要时形成类似 O 型的连续包裹层；源码还提供 `terminateLayersAtConcaveEdges` 控制。边、角、相交层和终止单元都有独立拓扑分支。

但本轮没有在主边界层代码中发现一个能对任意窄缝做全局可靠证明的单一“窄缝层杀手”。它依赖局部初始尺度、特征尺寸/曲率限制、相交检测、patch 编组、层终止和最终解缠共同防御。因此不应把 cfMesh 描述成对任意尖角和窄缝天然稳健。

## 2. snappyHexMesh 层生成调用链

### 2.1 总入口

`applications/utilities/mesh/generation/snappyHexMesh/snappyHexMesh.C` 顺序调用：

`snappyRefineDriver::doRefine → snappySnapDriver::doSnap → snappyLayerDriver::doLayers`

即先 castellated refinement，再 snap，最后 add layers。边界层主代码位于：

- `src/mesh/snappyHexMesh/snappyHexMeshDriver/snappyLayerDriver.C`
- `snappyLayerDriverSinglePass.C`
- `layerParameters/layerParameters.C`
- `src/dynamicMesh/polyTopoChange/polyTopoChange/addPatchCellLayer.C`
- `src/dynamicMesh/motionSmoother/`

### 2.2 参数与局部禁用

`layerParameters` 支持首层/末层/总厚度/增长率的多种组合，并提供逐 patch 的 `nSurfaceLayers`。重要质量与终止控制包括：

- `minThickness`、`featureAngle`、`concaveAngle`；
- `nGrow`：在不能挤出的点周围扩张禁用区；
- `maxFaceThicknessRatio`；
- `nBufferCellsNoExtrude`：使层数在终止区逐步减小；
- `nLayerIter`、`nRelaxedIter`；
- `meshShrinker` 与外迭代次数。

`setPointNumLayers()` 会在非流形点、非 string-connected 区域、过大特征角和严重翘曲面附近关闭挤出，并同步相邻处理结果。

### 2.3 “先缩体网格，再试建层”的事务链

单次主要链条可概括为：

1. 同步目标位移、最小厚度和逐点层数；
2. `getPatchDisplacement()` 计算 patch 位移；
3. 位移取反，把现有体网格边界向内收缩，由 `motionSmoother`/mesh mover 在严格质量约束下移动；
4. 截断过大位移，检测 pinch、butterfly、薄区和非法点；
5. 用 `nBufferCellsNoExtrude` 建立平缓终止区；
6. `addPatchCellLayer::setRefinement()` 把候选拓扑写进 `polyTopoChange`；
7. 先把修改应用到临时 `newMesh`，移动候选新点并运行网格检查；
8. `checkAndUnmark()` 找出错误面/错误新增层单元，取消相关 patch 点的挤出，恢复原点并重试；
9. 只有候选状态通过后，才把 `meshMod.changeMesh()` 提交到真实网格。

因此 snappy 的关键价值不是一个神奇法向公式，而是“候选—质量检查—局部取消—恢复—再试—最终提交”的事务式控制。超过 `nRelaxedIter` 后可改用 relaxed 质量集合，但不等于忽略非法几何。

## 3. cfMesh 与 snappyHexMesh 对比

| 维度 | cfMesh | snappyHexMesh | 对原生二维 H4 的启发 |
|---|---|---|---|
| 基本策略 | 先建立一层包裹单元，再细分 hair edges | 收缩现有体网格，再把原表面向外挤出成层 | 采用 cfMesh 的“先合法条带、后细分”更自然 |
| 拓扑控制 | 显式处理 patch 组、边、角、层相交和终止 | `polyTopoChange` 候选拓扑，逐点启停和终止缓冲 | H4 要有独立候选拓扑，不能原地半改 |
| 质量控制 | 厚度/法向优化、特征尺度限制、解缠 | 严格/放宽质量、临时网格检查、取消挤出与重试 | 复制 snappy 的事务思想，不复制其三维 mover |
| 对局部失败 | 编组、限制厚度、终止、优化 | 局部 unmark、增长禁用区、迭代重试 | H4 首版宜整条链 fail-closed，之后才做局部终止 |
| 原生二维适配 | 所谓 2D 是薄三维模式 | 原生三维算法 | 两者代码均不可直接套用 |
| 许可 | GPL | GPL | 只能思想借鉴 |

## 4. CartMesh2D 当前能力与 H4 缺口

当前主链是：

`Quadtree2D::refine → buildCutCells → buildGlobalTopology → agglomeration → buildSolverTopology2D → quality → OpenFOAM writer`

已有基础：

- `Quadtree2D` 支持 minimum、boundary、distance bands、boxes 的确定性目标层级，重叠处取最大细化级别；
- `CutCell2D` 生成真实二维多边形，默认支持 exterior 流体语义；
- 全局拓扑、求解器拓扑、质量检查和 OpenFOAM 输出已经形成独立层次；
- H3 已对十万至五十万级二维拓扑做过可复现实测，H4 不应破坏该数据路径。

H4 尚缺：

1. 从输入边界提取确定性、有方向的 wall chain；
2. 面法向到顶点推进方向的构造，尤其是凹角/凸角/尖角；
3. hair edge、层列、层号、内外包络等数据模型；
4. 全局 offset 自交、相邻边碰撞和窄缝检测；
5. layer strip 与剩余 Cartesian/Cut-cell 区域之间的共形接口；
6. 层终止、不同 patch 层数过渡、链端封口；
7. 候选拓扑的原子提交/回退；
8. 针对高长宽比单元的面积、凸性、正交性、skewness 与 OpenFOAM 验收。

特别要注意：现有 cut-cell 中的嵌入边界边，不等于一个已经适合挤出边界层的干净 body-fitted front。H4 需要显式建立这层语义，不能把中心采样或普通切割多边形改名为边界层。

## 5. 可以直接借鉴的算法思想

不复制源码的前提下，建议直接采用以下思想：

1. **cfMesh：先构造一层合法包裹条带，再沿 hair edge 细分。** 这把“外包络是否有效”和“层内增长率”拆成两个可验证问题。
2. **cfMesh：相关 patch/边角必须成组处理。** 顶点推进方向不能由每条边各自决定后再平均了事。
3. **snappy：候选拓扑事务。** 新层先建在独立对象中，质量失败不得污染当前可用网格。
4. **snappy：失败区扩张与明确退化。** 局部不能挤出的点会影响邻域，不能留下单点尖刺。
5. **两者共同点：厚度由几何可用空间约束。** 用户给出的总厚度只是上限，不是越过曲率半径或窄缝仍必须满足的命令。
6. **独立验收。** 内部面积/拓扑检查之外，仍须由独立读取器和真实 OpenFOAM `checkMesh` 验证。

## 6. 必须自行简化实现的部分

H4 不应移植三维 polyhedral 动态拓扑。原生二维可以把问题缩减为线段链、二维多边形和有向面积：

- 用 wall segment 的流体侧半平面交构造顶点 marching cone；
- 用二维线段索引检测候选 outer envelope 的自交和跨链碰撞；
- 用 `wall vertex ↔ envelope vertex` 的 hair edges 构造一层 quad strip；
- 对每根 hair 采用相同层数和确定性几何级数参数；
- 用现有二维多边形裁剪/拓扑工具，把 outer envelope 作为剩余域的新共形边界；
- 候选 strip、剩余 cut-cell 域和接口拓扑全部通过后再一次提交；
- 失败时保留最小失败案例和结构化原因，回退到现有纯 Cut-cell 网格，不隐藏层失败。

首版尤其不应一开始实现“每个顶点不同层数”。局部逐点取消会立即引入楔形终止、层列错接和非流形接口。应先要求一条连通 wall chain 使用固定层数；凹角、尖后缘或窄缝超出能力时，整条链明确失败或禁用。稳定后再加入受控的层终止模板。

## 7. 未来三维开放后端：应复用什么

未来三维建议把当前树结构抽象成窄接口，而不是让几何和 writer 直接依赖某个库：

```text
Geometry / SurfaceIndex
        ↓
AdaptiveTreeBackend
  - refine/coarsen
  - 2:1 balance
  - leaf traversal / Morton key
  - neighbor/search
  - ownership/ghost (parallel backend only)
        ↓
CutGeometryStore
  - classification / fractions / centroids / intersections
        ↓
ConformalSolverTopology
        ↓
Quality + OpenFOAM writer
```

这会保留两个可替换后端：当前确定性的原生串行 quadtree/octree，以及未来经许可审查后选择的并行 forest-of-octrees 后端。树后端只负责自适应空间分解，不负责宣称已生成 Cut-cell、边界层或可求解 polyhedral topology。

## 8. p8est 的适用性与硬边界

p4est/p8est 的核心数据结构是：连接关系定义的 macro-tree forest，每棵树中叶子按 Morton 顺序存放；p4est quadrant 为 `(x,y,level)`，p8est 增加 `z`，每个父单元分别有 4/8 个孩子。forest 维护全局/本地叶子计数、每棵树的有序叶数组和分区边界。

源码所示可复用能力：

- 回调式 refine/coarsen；
- face/edge/corner 级 2:1 balance；
- Morton/SFC 连续区间分区；
- local/partition/all search；
- ghost 与 mirror 层及邻居所有权；
- p4est 与 p8est 通过同构接口/映射宏共享大量算法结构。

这非常适合作为未来三维 AMR 的“树、平衡、搜索、并行所有权”后端，但它不提供 STL 求交、inside/outside 分类、Cut-cell 多面体构造、共形求解器拓扑、边界层或 OpenFOAM writer。非共形叶子与 hanging interface 也仍需上层转成求解器可接受拓扑。

[p4est 论文](https://p4est.github.io/papers/BursteddeWilcoxGhattas11.pdf) 给出了其 forest-of-octrees、Morton 分区及可扩展 refine/coarsen/balance/partition 的原始设计证据。

结论：**技术上值得预留适配器，当前不应承诺直接采用 p8est。** 原因不是性能，而是 GPL 与本项目许可证未定。如果未来产品选择 GPL 兼容发布，可在法律确认后直接链接；否则应选宽松许可后端或自研所需子集，不能通过“只借鉴”名义复制实现。

## 9. AMReX EB 的适用性

关键源码：

- `Src/EB/AMReX_EBCellFlag.H`
- `Src/EB/AMReX_MultiCutFab.H/.cpp`
- `Src/EB/AMReX_EBFabFactory.H`
- `Src/EB/AMReX_EBData.H`
- `Src/EB/AMReX_EB2_Level.H/.cpp`
- `Src/EB/AMReX_EBDataCollection.cpp`

最值得借鉴的是数据组织，而不是把 AMReX 当作 OpenFOAM 网格生成器：

- `EBCellFlag` 用紧凑位标识 regular、covered、single-valued、multi-valued 及邻接；
- `volfrac` 等规则场保持 MultiFab 形式；
- centroid、boundary centroid/normal/area、area fraction 等只对 cut 区域提供专用数据；
- `CutFabFactory` 只为含 single-valued cut cell 的 box 分配数据，属于 **box/FAB 粒度稀疏**，不能误说成每个 cut cell 独立稀疏对象；
- EB Level/DataCollection 明确分离每级几何、由细到粗 coarsening 和 all-regular 快路径。

但官方 EB 文档明确提醒当前算法实际不支持 multi-valued cells；`EBCellFlag` 能编码这种状态，不代表整条计算链能处理任意多连通切割几何。必须区分“表示能力”和“算法支持”。参见 [AMReX Embedded Boundary 官方文档](https://amrex-codes.github.io/amrex/docs_html/EB.html) 与 [AMReX 论文](https://journals.sagepub.com/doi/10.1177/10943420211022811)。

建议：未来三维可借鉴“密集分类/体积分数 + cut-box 稀疏几何 + 独立 solver topology”的布局。只有当产品本身进入 AMReX block-structured solver 生态时，才考虑直接以 AMReX EB 为几何/求解后端；若目标仍是独立网格器与 OpenFOAM polyMesh，它不是即插即用替代品。

## 10. Gmsh sizing field 对长期细化接口的启发

`src/mesh/Field.cpp`/`Field.h` 中与本项目最相关的抽象：

- `BoxField`：盒内外不同尺寸，并可在厚度带内插值；
- `DistanceField`：对点/曲线/表面采样，用空间索引求近似距离；
- `ThresholdField`：把距离或其他字段映射到 `SizeMin/SizeMax`；
- `MinField`/`MaxField`：组合多个字段；
- background field：统一向网格器提供空间尺寸请求。

CartMesh2D H1 当前把 minimum、boundary、distance bands、boxes 直接写入 quadtree 的 `requestedLevel()`，重叠处取最大 level。长期建议增加值语义、只读、可确定性组合的接口，例如：

```cpp
class SizingField2D {
public:
    virtual int targetLevel(const Point2& p) const = 0;
};
```

或先返回 `targetCellSize()`，再在唯一位置以明确舍入规则转换为整数 level。因为“尺寸越小越细”而“level 越大越细”，Gmsh 的 `Min(size)` 在 level 表示中对应 `Max(targetLevel)`。首批字段可为 Constant、BoundaryDistanceThreshold、Box、MaxLevelComposite；必须固定字段遍历顺序、浮点边界规则和 level 取整，保持相同输入的 ID/报告/输出确定性。

这里只借鉴抽象，不复制 Gmsh GPL 实现。

## 11. H4 推荐技术路线

### 11.1 候选比较

| 候选 | 可行性 | 与现有二维主链耦合 | 拓扑风险 | 结论 |
|---|---:|---:|---:|---|
| A. Cut-cell 主体 + 独立 body-fitted quad layers + 共形 transition | 高 | 中 | 中高 | 最接近正确目标 |
| B. 挤出层 + cut-cell transition/fallback | 中高 | 中 | 中高 | fallback 思想有用，但“transition”必须明确定义为共形裁剪，不能重叠 |
| C. 仿 cfMesh 全流程 | 低 | 高 | 高 | 只取“先包裹层后细分”思想 |
| D. 仿 snappy shrink/extrude/rollback | 低 | 极高 | 极高 | 只取候选事务和局部禁用思想 |
| E. 二维专用拓扑先行条带 + 保守 Cut-cell 剩余域 + 原子回退 | **最高** | 中 | 可分阶段控制 | **推荐**；是 A 的可执行简化版 |

### 11.2 推荐方案 E 的分解

1. **边界语义**：从输入 `BoundaryRegion` 提取确定性有向 wall chains，明确流体侧；圆、凸角、凹角、尖后缘、窄缝、多环分别建最小案例。
2. **推进方向**：由相邻 wall segments 的流体半平面构造顶点 marching cone；无法得到安全方向的点直接报错，不用任意平均法向掩盖冲突。
3. **安全外包络**：根据用户总厚度、局部边长、曲率/转角和全局最近非邻接边距离取保守上限；用线段索引检查 outer envelope 自交与跨链碰撞。
4. **一层条带**：每个 wall segment 与对应 outer-envelope segment 形成候选 quad；先验证方向一致、正面积、无重叠和完整链连接。
5. **层内细分**：所有 hair edges 使用同一层数，按几何级数放点；每层重新检查正面积、凸性和 wall-normal 单调距离。
6. **剩余域重建**：把 outer envelope 作为剩余 Cartesian/Cut-cell 域的新边界，移除条带占据区域，并沿外包络拆分相交背景单元，形成共享完全一致的接口边/点 ID。
7. **共形 transition**：允许接口附近产生一般多边形、三角形或四边形；不使用 overset，不保留重叠体积，也不靠 writer 容差粘合不一致面。
8. **终止策略**：首版一条连通 wall chain 固定层数；与非 wall patch 相接时只允许经过验证的确定性三角/楔形封口。凹角、尖角、窄缝暂不做隐式局部掉层。
9. **事务验收**：候选 mesh 独立构建，执行拓扑闭合、方向/面积、非流形、边界分类、质量和 OpenFOAM 导出检查；全部通过才替换当前 mesh。
10. **诚实回退**：失败时保留纯 Cut-cell 结果并报告 `layer_status=disabled/failed`、首个失败实体、几何尺度、阈值和最小复现输入，绝不把回退结果标成已生成边界层。

### 11.3 为什么不先做完整局部掉层

二维几何降低了面/体复杂度，却没有消除拓扑复杂度。一旦相邻顶点层数不同，就必须定义层终止边、楔形/三角过渡、相邻 column 配对和接口 ID。首个 H4 里同时解决这些问题会让失败原因不可分离。固定层数的整链条带能先回答三个基础问题：offset 是否安全、条带是否拓扑闭合、剩余 cut-cell 域能否与它共形拼接。

### 11.4 建议的首个验收阶梯

1. 单个平滑凸闭合 wall（圆或解析近似），固定层数/增长率；
2. 独立读取器验证面积、连通和边界标签；
3. 真实 OpenFOAM `checkMesh`，不得把内部 reader 冒充外部验收；
4. 再扩到温和凸角；
5. 凹角、尖后缘和窄缝作为明确 fail-closed 矩阵；
6. 只有上述证据稳定后，才讨论局部掉层和更复杂终止模板。

## 12. 未来三维架构路线

### 12.1 建议分层

未来三维不应从二维模板化加 `z=0`。建议独立工程/模块按以下职责分层：

1. `SurfaceKernel3D`：STL 清理、BVH、稳健求交、inside/outside 和特征边；
2. `AdaptiveTreeBackend`：原生 octree 或经许可批准的 forest-of-octrees 后端；
3. `CutGeometry3D`：每叶单元的分类、交点、截面多边形、体积分数和质心；
4. `CutTopology3D`：跨单元共享面的一致分割、全局 ID 和非流形诊断；
5. `BoundaryLayer3D`：独立 body-fitted prism layer、碰撞/终止及事务提交；
6. `HybridTransition3D`：prism 外包络到 Cartesian/Cut-cell 主体的共形过渡；
7. `SolverMesh3D`：质量、agglomeration、patch 与 OpenFOAM 输出。

### 12.2 后端选择结论

- p8est：最强项是并行 octree 基础设施；受 GPL 约束，当前只预留接口，不承诺依赖。
- AMReX EB：许可证友好，数据布局值得借鉴；更适合 AMReX 原生求解生态，不直接解决任意多面体 OpenFOAM 拓扑。
- 自研 backend：短期最能保持许可证、确定性和数据语义可控，但并行 balance/partition/ghost 成本最高。
- 推荐：先冻结 `AdaptiveTreeBackend` 和 `CutGeometryStore` 的最小契约，保持串行原生实现；待产品许可与求解器方向确定后，再用同一验收集比较 p8est、AMReX 或其他宽松许可后端。

### 12.3 与边界层论文的关系

公开研究普遍采用“近壁 body-fitted layer + 远场 Cartesian + 中间共形或 overset 连接”的分解。例如 [Wang 与 Kannan 的 Cartesian/prism overset 方法](https://cfd.ku.edu/papers/aiaa-2005-0322.pdf) 先生成贴体棱柱层，再与自适应 Cartesian 背景重叠；[Park 等人的 hybrid grid 方法](https://doi.org/10.1002/fld.3691) 使用曲率控制的 quad/prism extrusion、碰撞候选检测，并用 constrained Delaunay triangle/tetrahedra 填充过渡区。

这些论文支持“边界层是独立几何/拓扑子系统”的判断，但不意味着 H4 应采用 overset 或 Delaunay。对当前原生二维 OpenFOAM 输出目标，最小且可审计的路线仍是：quad strip 与 cut-cell remainder 共享一条完全共形的 outer envelope。

## 13. 主要风险与停线条件

### 13.1 主要风险

- 凹角、尖后缘和窄缝使 offset envelope 自交或跨越另一条 wall；
- 高增长率/高长宽比导致负面积、强 skew 或 OpenFOAM 低质量面；
- 条带 outer envelope 与背景 cut-cell 交点未使用同一几何构造，造成微小缝隙、重复边或非流形点；
- 局部掉层产生未定义的终止拓扑；
- 多环/多流体分量语义与当前 local-hole 限制冲突；
- 第三方 GPL 代码被无意复制，或 p8est 被误列为许可证无条件可兼容；
- 把候选失败后的纯 Cut-cell 回退报告成成功生成 layer。

### 13.2 必须停线

发生任一情况即停止扩展案例并保留最小失败输入：

1. 任意 layer/transition cell 面积非正、方向冲突或自交；
2. outer envelope 与 wall/其他 envelope 非预期相交；
3. 接口边不能做到两侧唯一共享，或出现 dangling/non-manifold topology；
4. 内部质量通过但独立 reader 或真实 `checkMesh` 失败；
5. 只能通过放宽阈值、删除坏单元或隐藏 warning 才能通过；
6. 需要复制 GPL 实现才能继续，而项目许可证尚未决策。

## 14. 本轮实际阅读的核心材料

### 本地项目

- `WINDOWS_CODEX_HANDOFF_CN.md`
- `cartmesh2d/AGENTS.md`
- `cartmesh2d/docs/STAGE2DH1_VERIFICATION_CN.md`
- `cartmesh2d/docs/STAGE2DH3_SOLVER_TOPOLOGY_SCALABILITY_CN.md`
- `cartmesh2d/include/cartmesh2d/quadtree/Quadtree2D.hpp`
- `cartmesh2d/src/quadtree/Quadtree2D.cpp`
- `cartmesh2d/include/cartmesh2d/cutcell/CutCell2D.hpp`
- `cartmesh2d/include/cartmesh2d/topology/Topology2D.hpp`
- `cartmesh2d/include/cartmesh2d/quality/SolverTopology2D.hpp`

### cfMesh / OpenFOAM

- cfMesh 的 `cartesianMeshGenerator.C`、`cartesian2DMeshGenerator.C`
- cfMesh 的 `boundaryLayers*`、`refineBoundaryLayers*`、`boundaryLayerOptimisation*`
- OpenFOAM 的 `snappyHexMesh.C`、`snappyLayerDriver*`、`layerParameters*`
- OpenFOAM 的 `addPatchCellLayer*`、`polyTopoChange*`、`motionSmoother*`

### p4est/p8est / AMReX / Gmsh

- p4est 的 `p4est.h`、`p8est.h`、`p4est.c`、`p4est_balance*`、`p4est_search*`、`p4est_ghost*`、`p4est_bits.c`
- AMReX 的 `EBCellFlag`、`MultiCutFab`、`EBFabFactory`、`EBData`、`EB2_Level`、`EBDataCollection`
- Gmsh 的 `Field.h`、`Field.cpp`
- 各源码包许可证文件/源文件头；本地 `papers/` 为空

## 15. 最终结论

H4 不应照搬 cfMesh 或 snappyHexMesh，而应采用二维专用的候选 E：**先构造并验证一个 body-fitted quad wrapper strip，再沿 hair edges 分层；用该 strip 的 outer envelope 重新裁剪剩余 Cartesian/Cut-cell 域，形成共形 polygon transition；整个候选拓扑通过独立质量链后原子提交，失败则明确回退。**

未来三维应把 AMR 树后端与 Cut-cell 几何、求解器拓扑、边界层和 writer 解耦。p8est 技术上适合 tree/balance/search/partition/ghost，但 GPL 使其不能在当前许可证未定时被列为默认直接依赖；AMReX 许可证更友好、EB 数据布局很有价值，但不是独立 OpenFOAM polyhedral 网格器。现在最正确的动作是预留窄后端接口并保留原生实现，而不是提前导入任何一个大型框架。
