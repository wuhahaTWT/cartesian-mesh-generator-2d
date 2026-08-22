# cartmesh2d 架构设计

## 1. 分层原则

依赖方向必须单向：

```text
geometry
   ↓
grid
   ↓
quadtree
   ↓
cutcell
   ↓
topology
   ↓
quality
   ↓
io
   ↓
apps
```

上层可依赖下层；下层不得反向依赖 CLI、可视化或文件格式。

## 2. 模块职责

### `geometry/`
负责二维连续几何和鲁棒谓词：

- Point2D / Vector2D
- Segment2D
- AABB2D
- Polygon2D
- BoundaryLoop
- orientation
- point-on-segment
- segment intersection
- point-in-polygon
- signed area / centroid
- boundary diagnostics

不包含网格对象，也不决定哪一侧是物理流体。

### `grid/`
负责规则二维 Cartesian 背景域和 leaf cell 的基础表达：

- Domain2D
- CartesianCell2D
- CellClass: Inside / Outside / Intersected
- cell AABB / level / logical indices

`CellClass` 只是相对于 `BoundaryLoop` 的**几何分类**。严禁在本层把 `Inside` 直接定义为流体或固体。

### `quadtree/`
负责：

- 1 -> 4 refinement；
- leaf 存储；
- Morton/Z-order 编码（若采用线性 Quadtree）；
- surface/distance/curvature refinement criteria；
- 2:1 balance；
- 确定性 leaf ordering。

### `cutcell/`
负责边界与 leaf rectangle 的真实裁切以及物理 fluid-side 映射：

- `FluidRegion2D::Exterior / Interior`；
- 默认 `Exterior`：闭合边界代表固体，外侧是 CFD 流体；
- boundary intersection points；
- CutPolygon；
- fluid area fraction；
- polygon centroid；
- embedded-boundary edge；
- pathological cut detection。

默认外流时，CCW 固体边界的 embedded fragment 在构造 fluid polygon 时必须反向，使流体始终位于有向边的左侧。不能只把 `Inside/Outside` 标签对调而忽略边界方向。

### `topology/`
把局部几何转换为全局唯一拓扑：

- Vertex2D；
- Edge2D；
- Cell2D；
- owner / neighbour；
- boundary patch；
- coarse-fine edge 拆分规则；
- deterministic IDs。

默认外流拓扑必须同时包含：

- `EmbeddedBoundary`：固体壁面；
- `DomainBoundary`：外部计算域边界。

### `quality/`
负责：

- positive area；
- edge length；
- aspect ratio；
- skewness/centroid metrics；
- duplicate/orphan/non-manifold checks；
- 质量报告与失败定位。

除常规几何质量外，产品 acceptance 还必须验证物理域面积：默认外流 `fluid_area = domain_area - solid_area`。

### `io/`
只负责序列化和交换：

- JSON debug/acceptance report；
- CSV（可选）；
- VTK/VTU；
- 后续 solver-oriented 2D export。

不得把网格生成算法塞进 writer。导出元数据应明确 `fluid_region` 与 `boundary_role`，避免可视化或后处理误解物理域。

## 3. 核心数学对象

### 3.1 Point / Vector

```cpp
struct Point2D { double x, y; };
struct Vector2D { double x, y; };
```

### 3.2 Segment

```cpp
struct Segment2D {
    Point2D a;
    Point2D b;
};
```

### 3.3 AABB

```cpp
struct AABB2D {
    Point2D min;
    Point2D max;
};
```

### 3.4 Polygon

存储有序顶点，不重复保存最后一个首点。统一规定有效 fluid polygon 顶点按 CCW 排列，使正 signed area 表示正向流体区域。

### 3.5 BoundaryLoop

BoundaryLoop 与一般 Polygon2D 概念上分开：它代表输入几何边界，包含诊断状态和方向；Polygon2D 是算法中间结果。

默认产品中 `BoundaryLoop` 表示固体壁面，而不是“默认流体外包络”。物理侧由 `FluidRegion2D` 显式决定。

## 4. 鲁棒性策略

### 4.1 tolerance

建立统一配置，例如：

- absolute geometry epsilon；
- relative epsilon based on domain length；
- snapping/dedup tolerance；
- minimum area tolerance。

严禁各函数独立写 `1e-9`、`1e-12` 等魔法数。

### 4.2 orientation

基本谓词：

`orient(A,B,C) = cross(B-A, C-A)`

所有 segment intersection、polygon orientation、point-on-edge 判断均围绕统一谓词实现。

### 4.3 point-in-polygon

必须三态返回：

- Inside
- Outside
- Boundary

不能用 bool 丢失“恰在边界”的状态；也不能在 geometry 层把三态直接解释成物理流体/固体。

## 5. Quadtree 约束

每个 leaf 至少有：

- level；
- integer logical `(i,j)`；
- deterministic key（推荐 Morton code + level）；
- geometric AABB；
- classification。

2:1 第一版定义为**共享边相邻 leaf 的层级差 <= 1**。

平衡算法必须迭代到不再产生新 refinement，且结果与遍历顺序无关。

## 6. Cut-cell 定义

对于被固体 `BoundaryLoop` 穿过的 leaf rectangle `R`，Cut-cell 结果应表示：

`F = R ∩ FluidDomain`

默认外流：

`FluidDomain = Domain2D \ SolidInterior`

显式内部流：

`FluidDomain = SolidInterior`

第一版可限定“每个切割单元内 fluid intersection 是单连通 polygon”；遇到多片流体区域或单 leaf 内 hole 时必须显式报 Unsupported，不能静默丢片。

输出至少包括：

- polygon vertices；
- area；
- centroid；
- area fraction；
- Cartesian edge fragments；
- embedded boundary fragments。

## 7. Coarse-fine 拓扑

Quadtree 2:1 网格存在 hanging nodes。全局拓扑构造时，粗 cell 的共享边必须按细邻居节点拆分成一致 edge fragments，确保每条内部 edge 恰有：

- 1 owner；
- 1 neighbour。

不能让一条粗 edge 同时对应两个细 edge 而没有显式拆分。

## 8. 小 Cut-cell

定义：

`alpha = A_fluid / A_background`

2D-5 至少实现：

- alpha 统计；
- 阈值标记；
- 邻接候选评估；
- 可验证的处理策略。

优先研究 cell agglomeration，但具体算法在 2D-5 开始前再锁定。默认外流模式下 alpha 必须对应**固体外侧**剩余流体面积。

## 9. 确定性

相同输入必须保持：

- leaf 排序一致；
- vertex/edge/cell ID 一致；
- report 字段一致；
- 导出顺序一致。

禁止依赖 unordered container 的非确定遍历顺序直接分配最终 ID。

## 10. 可视化边界

核心库不得依赖 matplotlib、ParaView、Qt 或图形框架。

最终可视化只能消费 `io` 输出，不能反向参与网格生成。默认外流图必须显示“固体内部为空、外围为 CFD 网格”，但截图本身不能替代数值面积与拓扑验收。
