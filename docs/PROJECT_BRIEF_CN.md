# 原生二维自适应 Cartesian / Cut-cell 网格生成器项目总纲

文档状态：二维子项目立项基线（CFD 物理侧定义已校正）  
日期：2026-08-19

## 1. 项目目的

在现有三维 Cartesian 网格生成器继续推进的同时，建立一个技术路线独立、实现封闭、开发风险更低的原生二维网格生成器，作为可稳定验收的结题成果，并同时承担二维算法验证平台的作用。

二维项目不是三维项目的降维开关。它从数据结构开始即采用二维几何与二维拓扑：线段边界、矩形背景单元、Quadtree、自适应二维 Cut-cell polygon、edge-cell 邻接。

本项目最终服务于 CFD 前处理，因此必须先明确“几何轮廓”和“流体域”的物理含义，而不能只做抽象 polygon 网格化。

## 2. 最终产品定义

默认外流产品输入为：

1. 一个合法闭合二维固体轮廓 `BoundaryLoop`（例如圆柱、翼型、叶片截面或其他障碍物）；
2. 一个包围该物体的外部计算域 `Domain2D`。

默认输出流体域定义为：

`fluid = Domain2D - solid interior`

也就是说，**物体内部不生成流体网格，物体外部到外部计算域边界之间生成 CFD 控制体网格**。这一语义必须与三维 `cartmesh` 保持一致。

自动生成的二维 CFD 控制体网格应具有：

- Cartesian 背景网格；
- 固体边界附近自适应 Quadtree 细化；
- 2:1 平衡；
- inside / outside / intersected **几何分类**；
- 由 fluid-side policy 将几何分类映射为物理流体/固体；
- 默认保留固体轮廓**外侧**的真实 Cut-cell 多边形；
- 完整 vertex-edge-cell-neighbor 拓扑；
- 固体壁面 `EmbeddedBoundary`；
- 外部计算域 `DomainBoundary`；
- 小 Cut-cell 检测与处理；
- 网格质量与拓扑质量报告；
- 标准化可机器读取输出；
- 最后再提供可视化。

内部流/管道流属于显式模式，只有调用方指定 `FluidRegion2D::Interior` 时才保留闭合轮廓内部作为流体。

第一版重点支持单个闭合固体轮廓；架构允许后续加入多个实体、孔洞与多区域。

## 3. 推荐输入

优先级：

1. 内置解析/离散固体几何：rectangle、circle、airfoil-like fixture；
2. 简单文本固体边界：`.dat` / `.csv` / `.json` 点序列；
3. 后续再考虑更复杂 CAD/网格交换格式。

二维固体边界的基本表示为有序点列：

`p0 -> p1 -> ... -> pn-1 -> p0`

由相邻点构成 `Segment2D`。输入方向可被标准化，但方向本身不得被误用来偷偷切换流体侧；物理流体侧由 `FluidRegion2D` 明确定义。

## 4. 核心算法路线

```text
solid BoundaryLoop + outer Domain2D
 -> validate/orient
 -> Cartesian background domain
 -> segment-AABB candidates
 -> geometric classify cells
 -> Quadtree refine
 -> 2:1 balance
 -> reclassify leaves
 -> select physical fluid side (default EXTERIOR)
 -> clip exterior fluid region / construct CutPolygon
 -> reject or handle pathological tiny cells
 -> build global topology
 -> require solid-wall + outer-domain patches
 -> quality/physics audit
 -> export
```

## 5. 非目标

在 2D-6 之前明确不做：

- GUI；
- CFD Navier-Stokes 求解器；
- 三维 STL 输入；
- 三维代码模板化重构；
- AI 自动网格参数选择；
- 云服务；
- 为了展示效果而跳过几何/拓扑正确性。

同时明确：二维项目不是“给闭合图形内部铺网格”的通用图形工具。默认用途是为后续 CFD 计算生成**固体外部流场网格**。

## 6. 技术判定标准

### 6.1 “Cartesian 网格”
背景单元边与全局 x/y 坐标轴平行。

### 6.2 “自适应网格”
局部 leaf level 可变化，并满足规定的 2:1 邻接平衡，不是简单全局加密。

### 6.3 “Cut-cell”
被固体边界穿过的背景单元，必须构造实际流体区域多边形；默认外流时是**边界外侧**的流体区域。只有 inside/outside 标签或中心点采样不算 Cut-cell。

### 6.4 “求解器可用拓扑”
至少存在确定的 vertices、edges、cells、owner/neighbour、boundary patch 关系，并通过拓扑闭合检查。默认外流必须同时存在固体壁面与外部计算域边界。

### 6.5 “物理域正确”
对于单个闭合固体且完全位于计算域内的默认外流案例，必须满足：

`total_fluid_area = domain_area - solid_area`

并且几何 `Inside` 单元不能出现在最终流体网格中。该条件优先级高于截图效果和普通 topology audit。

## 7. 结题最低完整线

最低不能低于：

`2D-0 + 2D-1 + 2D-2 + 2D-3 + 2D-4 + 2D-6`

`2D-5` 小单元聚合若时间不足，可以先做到可靠检测、标记和拒绝/阈值策略；但不得假装问题不存在。

结题验收必须至少展示一个“固体轮廓为空洞、外围为自适应 CFD 网格”的真实外流案例，不能用固体内部网格代替。

## 8. 与三维项目的关系

三维项目承担更高目标：3D STL、Octree、3D Cut-cell、OpenFOAM、多尺度与后续混合近壁网格。

三维项目当前的物理语义是 STL `outside = fluid`、STL `inside = solid`；二维项目必须保持同一默认语义，避免二维/三维产品逻辑相反。

二维项目承担：

- 更低风险的完整算法闭环；
- 更清楚的算法解释与论文图示；
- 对自适应、平衡、Cut-cell、拓扑算法的独立验证；
- 结题稳定版本。

二维成功不意味着停止三维；三维失败也不得破坏二维可交付性。
