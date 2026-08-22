# 2D-R 复杂几何鲁棒性修复

## 当前目标

复杂几何鲁棒性必须建立在正确的默认 CFD 物理域上：

```text
BoundaryLoop = solid wall
FluidRegion2D::Exterior = default
fluid = Domain2D - solid interior
```

2026-08-20 的历史 robustness 工作解决了大量 clipping/tolerance/topology 问题，但当时默认 fluid side 仍是 interior。2026-08-22 物理域纠正后，全部复杂案例必须重新按 EXTERIOR 语义验证。

## 历史失败基线

旧 interior-fluid stress 曾出现：

- gear_star：Cut-cell Unsupported（大量）
- serpentine_body：Cut-cell Unsupported（少量）
- naca2412_dense：source global topology audit fail
- superellipse_24：source global topology audit fail
- nozzle_profile：PASS

后续旧语义下这些失败曾被修到全绿，但该绿灯不能替代 corrected external-fluid revalidation。

## 当前修复原则

1. Intersected leaf 不依赖“凹 polygon 裁到 AABB 后仍只有一个 polygon”的假设。
2. 从 local embedded-boundary fragments 与 Cartesian-cell perimeter fluid intervals 构造有向边界图。
3. 有向边界必须按**所选物理 fluid side**构造；默认 exterior 时反转 CCW solid embedded fragments，使流体保持在 directed edge 左侧。
4. 一个 leaf 中出现多个 disconnected fluid components 时，solver 路径使用 `buildCutCells(...)` 全部发射为独立 solver cells；不得丢片、不得造假桥。
5. 单组件 convenience API 遇到多组件仍显式 Unsupported，防止调用者无意丢失流体区域。
6. 单 leaf 内真正 polygon-with-hole 当前仍显式 Unsupported；必须通过进一步 refinement 或未来 hole topology 解决。
7. AABB 裁剪交点在 tolerance 内吸附到精确 cell-side 坐标，降低相邻 cell topology 漂移。
8. 面积门限随 background-cell area 缩放，不能把真实极小 Cut-cell sliver 当成空单元。
9. 共线删除基于归一化角度判断，不能吞掉极短但真实的边界转折。
10. 默认外流 acceptance/stress 必须同时通过 physical area、EmbeddedBoundary、DomainBoundary 与 topology audit，不允许只看截图或 `valid=true`。

## 当前永久压力目标

以下输入必须全部用默认 EXTERIOR CLI 重新跑通：

- `gear_star.xy`
- `serpentine_body.xy`
- `nozzle_profile.xy`
- `naca2412_dense.xy`
- `superellipse_24.xy`

每个案例必须满足：

```text
solid interior has no fluid cell
source_fluid_area = domain_area - solid_area
EmbeddedBoundary > 0
DomainBoundary > 0
unsupported = 0
unresolved small cell = 0
topology audit = all zero
export/read-back = exact
SVG = solid blank + surrounding adaptive mesh
```

同时继续报告极端 `aspect ratio / skewness / min edge / min area fraction`，不能用拓扑合法掩盖质量不足。

## 已知 solver-readiness 后续项

- 显式 farfield/domain extents，而不是只用 bbox padding；
- named outer patches（inlet/outlet/top/bottom/farfield）；
- 多 solid loops / holes；
- polygon-with-holes；
- solver-native export；
- near-wall boundary-layer / hybrid mesh；
- mesh-quality optimization。

当前分支：`agent/native-2d-robustness`。

完整物理域审计：`PHYSICAL_DOMAIN_AUDIT_2026-08-22.md`
