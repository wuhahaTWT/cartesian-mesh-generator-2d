# cartmesh2d 分阶段执行计划

## 总原则

每个阶段都必须满足：

1. 代码；
2. 单元/回归测试；
3. 最小失败案例；
4. 阶段验收报告；
5. 不提前混入下一阶段功能；
6. 默认 CFD 语义始终是“闭合输入轮廓 = 固体，流体 = 外部计算域减去固体内部”。

几何 `Inside/Outside` 只是分类结果，不等于物理 `fluid/solid`；真正的物理侧由 `FluidRegion2D` 在 Cut-cell 阶段显式决定。默认必须为 `Exterior`。

---

## 2D-0 — 二维几何内核

### 目标
建立不依赖网格的二维鲁棒几何基础。

### 实现
- Point2D / Vector2D
- Segment2D
- AABB2D
- Polygon2D
- BoundaryLoop
- orientation / cross / dot / norm
- signed area / absolute area
- polygon centroid
- segment-segment intersection，区分 none / point / overlap
- point-on-segment
- point-in-polygon 三态
- boundary closure / duplicate point / zero-length edge / self-intersection / orientation diagnosis

### 测试夹具
- axis-aligned rectangle
- triangle
- concave L-shape
- approximate circle
- bow-tie self-intersection（必须拒绝）
- duplicate consecutive vertex（必须诊断）

### Gate
只有 `docs/ACCEPTANCE_CN.md` 中 2D-0 全部通过才能进入 2D-1。

---

## 2D-1 — Uniform Cartesian + 几何分类

### 目标
生成真正的二维规则 Cartesian 网格，并可靠判断 cell 与固体边界关系。

### 实现
- Domain2D
- Nx × Ny / spacing-based grid
- CartesianCell2D
- segment-AABB candidate/intersection
- CellClass: Inside / Outside / Intersected
- deterministic row-major/Morton ordering（锁定一种）
- 分类统计报告

### 关键要求
- “中心点 inside”只能辅助判断完整未相交 cell；
- 只要边界与 cell 相交，必须标为 Intersected；
- tangent / boundary-on-grid-line 情况必须有明确规则和测试；
- 本阶段禁止把 `Inside` 硬编码成“流体”，否则会把外流问题物理侧翻转。

---

## 2D-2 — Quadtree + 2:1 balance

### 目标
固体边界附近局部细化，不全局加密。

### 实现
- 1 -> 4 refinement
- base level / max level
- boundary refinement
- distance refinement
- curvature refinement（可作为 2D-2.x 后半部分）
- face-neighbor discovery
- 2:1 balance
- deterministic leaf IDs/keys

### 验收核心
- 局部细化确实只发生在目标区域；
- leaf 不重叠、不缺洞；
- 所有共享边邻居 level 差 <= 1。

---

## 2D-3 — Cut-cell polygon + physical fluid side

### 目标
相交 Cartesian leaf 形成真实流体 polygon，并把几何分类正确映射成 CFD 流体域。

### 实现
- `FluidRegion2D::Exterior / Interior`
- 产品默认 `Exterior`
- rectangle-boundary intersections
- vertex ordering
- polygon clipping/construction
- 外流时 embedded fragment 正确反向
- area / centroid
- area fraction
- embedded boundary fragment
- pathological topology detection

### 默认外流硬门
- geometric `Inside -> Empty solid`；
- geometric `Outside -> Full fluid`；
- `fluid_area = domain_area - solid_area`；
- 显式 `Interior` 模式保留作为内部流用途，但不能影响默认产品语义。

### 禁止替代
以下均不算完成：
- 删除 intersected cell；
- 把相交 cell 整格保留；
- 只输出 volume/area fraction 而没有 polygon；
- 只用图片显示边界；
- 仅把 Inside/Outside 标签互换、却没有同步处理 embedded edge 方向；
- 在固体内部铺网格后称其为默认外流 CFD 网格。

---

## 2D-4 — 全局拓扑

### 目标
形成 solver-oriented 的二维控制体 connectivity。

### 实现
- global vertex deduplication
- edge construction
- owner / neighbour
- boundary patch
- coarse-fine edge splitting
- Cut-cell embedded edge integration
- deterministic IDs

### Gate
内部 edge 恰好两个关联 cell，边界 edge 恰好一个 owner；不存在 orphan/duplicate/non-manifold edge。

默认外流还必须：
- `EmbeddedBoundary > 0`（固体壁面）；
- `DomainBoundary > 0`（外部计算域）；
- topology total area = domain area - solid area。

---

## 2D-5 — Small Cut-cell

### 目标
识别并处理低面积分数流体单元。

### 第一子阶段
- alpha histogram
- threshold marking
- nearest/best neighbour candidate
- failure report

### 后续子阶段
优先实现 topology-safe cell agglomeration；若研究后决定其他稳定化方法，必须先更新本计划再编码。

默认外流稳定化后必须继续保持固体内部为空，不能通过聚合把流体“桥接”进固体内部。

---

## 2D-6 — Quality + Export + Final Acceptance

### 目标
把二维生成器封口成可交付 CFD 网格产品核心。

### 实现
- geometry quality
- topology audit
- physical-domain audit
- min area / min edge / aspect ratio / skewness
- deterministic report
- JSON machine-readable report
- VTK/VTU export
- 至少一种面向 CFD/数值求解的拓扑导出或清晰格式规范
- benchmark fixtures
- end-to-end CLI

### 最终 end-to-end

`solid boundary file + outer domain -> exterior fluid mesh generation -> physics/topology audit -> export -> independent read-back verification`

CLI 默认必须输出外流网格；若支持内部流，必须通过显式参数选择，不能靠输入点顺/逆时针偷偷切换。

---

## 2D-V — Visualization（最后）

只有 2D-6 通过后进入。

第一版只做薄层可视化：读取已导出的网格文件，展示：

- cell edges；
- Quadtree level；
- Cut-cells；
- solid boundary；
- outer computational boundary；
- small-cell/quality flags。

默认外流图中固体内部必须为空；外围流体域应显示网格，并在固体附近局部细化。

可视化不得复制或重新实现核心分类/Cut-cell 算法，也不得代替物理域数值验收。
