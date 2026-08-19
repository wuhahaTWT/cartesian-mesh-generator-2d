# 原生二维自适应 Cartesian / Cut-cell 网格生成器项目总纲

文档状态：二维子项目立项基线  
日期：2026-08-19

## 1. 项目目的

在现有三维 Cartesian 网格生成器继续推进的同时，建立一个技术路线独立、实现封闭、开发风险更低的原生二维网格生成器，作为可稳定验收的结题成果，并同时承担二维算法验证平台的作用。

二维项目不是三维项目的降维开关。它从数据结构开始即采用二维几何与二维拓扑：线段边界、矩形背景单元、Quadtree、自适应二维 Cut-cell polygon、edge-cell 邻接。

## 2. 最终产品定义

输入一个或多个合法的二维闭合边界及计算域参数，自动生成具有以下特征的二维 CFD 控制体网格：

- Cartesian 背景网格；
- 边界附近自适应 Quadtree 细化；
- 2:1 平衡；
- inside / outside / intersected 几何分类；
- 真实边界裁切后的 Cut-cell 多边形；
- 完整 vertex-edge-cell-neighbor 拓扑；
- 小 Cut-cell 检测与处理；
- 网格质量与拓扑质量报告；
- 标准化可机器读取输出；
- 最后再提供可视化。

第一版重点支持单个外边界；架构允许后续加入多个实体和孔洞。

## 3. 推荐输入

优先级：

1. 内置解析/离散几何：rectangle、circle、airfoil-like fixture；
2. 简单文本边界：`.dat` / `.csv` / `.json` 点序列；
3. 后续再考虑更复杂 CAD/网格交换格式。

二维边界的基本表示为有序点列：

`p0 -> p1 -> ... -> pn-1 -> p0`

由相邻点构成 `Segment2D`。

## 4. 核心算法路线

```text
BoundaryLoop
 -> validate/orient
 -> Cartesian domain
 -> segment-AABB candidates
 -> classify cells
 -> Quadtree refine
 -> 2:1 balance
 -> reclassify leaves
 -> clip fluid region / construct CutPolygon
 -> reject or handle pathological tiny cells
 -> build global topology
 -> quality audit
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

## 6. 技术判定标准

### 6.1 “Cartesian 网格”
背景单元边与全局 x/y 坐标轴平行。

### 6.2 “自适应网格”
局部 leaf level 可变化，并满足规定的 2:1 邻接平衡，不是简单全局加密。

### 6.3 “Cut-cell”
被边界穿过的背景单元，必须构造实际流体区域多边形；只有 inside/outside 标签或中心点采样不算 Cut-cell。

### 6.4 “求解器可用拓扑”
至少存在确定的 vertices、edges、cells、owner/neighbour、boundary patch 关系，并通过拓扑闭合检查。

## 7. 结题最低完整线

最低不能低于：

`2D-0 + 2D-1 + 2D-2 + 2D-3 + 2D-4 + 2D-6`

`2D-5` 小单元聚合若时间不足，可以先做到可靠检测、标记和拒绝/阈值策略；但不得假装问题不存在。

## 8. 与三维项目的关系

三维项目承担更高目标：3D STL、Octree、3D Cut-cell、OpenFOAM、多尺度与后续混合近壁网格。

二维项目承担：

- 更低风险的完整算法闭环；
- 更清楚的算法解释与论文图示；
- 对自适应、平衡、Cut-cell、拓扑算法的独立验证；
- 结题稳定版本。

二维成功不意味着停止三维；三维失败也不得破坏二维可交付性。
