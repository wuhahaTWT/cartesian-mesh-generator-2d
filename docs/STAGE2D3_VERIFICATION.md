# Stage 2D-3 验证记录

## 当前状态

**REOPENED / HISTORICAL PASS INVALIDATED FOR DEFAULT CFD PRODUCT**

2026-08-22 全项目物理域审计发现，本文件原先把：

```text
BoundaryLoop interior = retained fluid
Inside leaf -> Full fluid
Outside leaf -> Empty
```

当作默认产品语义。这与项目真实 CFD 目标以及三维 `cartmesh` 的 `outside=fluid / inside=solid` 语义相反。

因此原先的“PASS / CLOSED”只能说明旧内部流定义下的几何实现自洽，**不能再作为默认外流 CFD Cut-cell 的验收结论**。

## 纠正后的 Stage 2D-3 定义

默认：

```text
BoundaryLoop = solid wall
FluidRegion2D::Exterior = default
geometric Inside  -> Empty solid
geometric Outside -> Full fluid
Intersected -> retain exterior fluid polygon(s)
```

内部流只有显式 `FluidRegion2D::Interior` 才保留旧的“boundary interior = fluid”语义。

外流 embedded fragment 必须反向，使 retained fluid 位于 directed edge 左侧；不能只交换 Inside/Outside 标签。

## 新硬门

- 默认外流解析切割面积/质心正确；
- 默认 `Inside -> Empty`、`Outside -> Full`；
- circle/closed solid 跨 Quadtree 总流体面积 = `domain_area - solid_area`；
- 多个 disconnected fluid components 不得丢失：solver 路径使用 `buildCutCells(...)` 全部发射；
- 单组件 API 遇到 multi-component 必须显式 Unsupported；
- 单 leaf 内真正 hole 当前仍显式 Unsupported。

## 历史数据说明

原文件中诸如：

```text
fluid_area = input_polygon_area
Inside leaf -> Full fluid
Outside leaf -> Empty
```

的结果属于旧内部流 fixture，不得再用于默认产品验收。

完整根因与纠正记录见：

`docs/PHYSICAL_DOMAIN_AUDIT_2026-08-22.md`

在 current corrected head 完成新的全量 CTest/CLI 验收前，Stage 2D-3 不再标记 CLOSED。
