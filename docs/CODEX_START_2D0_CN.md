# Codex 当前接管状态：二维物理域纠正后重新验收中

> **历史文档，已完全过时。** 这是 2026-08-22 的一次性接管快照。下面的 `REOPENED`
> 表格已由 `docs/R1F_PATCH_LOCAL_CLOSEOUT_CN.md` 第 3 节结清；阅读清单里的
> `cartmesh2d/` 路径前缀属于拆仓前的目录布局，CI 现在明确断言该目录不存在。
> 接手本项目请读 `docs/R2_HANDOFF_CN.md` 与 `docs/CURRENT_STATE_CN.md`。

## 当前状态（2026-08-22）

- Stage 2D-0：几何算法保留，PASS 历史有效
- Stage 2D-1：Cartesian/几何分类算法保留，PASS 历史有效
- Stage 2D-2：Quadtree/2:1 算法保留，PASS 历史有效
- Stage 2D-3：**REOPENED** — 默认 fluid side 已从错误 interior 改为 EXTERIOR
- Stage 2D-4：**REOPENED** — 必须在 external CFD topology 上重新验收
- Stage 2D-5：算法已实现，但必须在 external domain 重新验证
- Stage 2D-6：**REOPENED** — 旧 PASS/Actions 不再证明产品物理域正确
- Stage 2D-V：**REOPENED** — 旧 SVG 忠实画出了错误的 inside-fluid mesh

严禁继续使用“2D-0~V 全部封口”这一旧结论。

## 开始任何二维工作前必须阅读

1. `cartmesh2d/AGENTS.md`
2. `cartmesh2d/docs/PHYSICAL_DOMAIN_AUDIT_2026-08-22.md`
3. `cartmesh2d/docs/PROJECT_BRIEF_CN.md`
4. `cartmesh2d/docs/ARCHITECTURE_CN.md`
5. `cartmesh2d/docs/STAGE_PLAN_CN.md`
6. `cartmesh2d/docs/ACCEPTANCE_CN.md`
7. `cartmesh2d/docs/STAGE2D3_VERIFICATION.md`
8. `cartmesh2d/docs/STAGE2D4_VERIFICATION.md`
9. `cartmesh2d/docs/STAGE2D5_VERIFICATION.md`
10. `cartmesh2d/docs/STAGE2D6_VERIFICATION.md`
11. `cartmesh2d/docs/STAGE2DV_VERIFICATION.md`

## 第一优先级物理规则

默认二维产品与三维 `cartmesh` 必须一致：

```text
BoundaryLoop = solid wall / obstacle
Domain2D = outer computational domain
FluidRegion2D::Exterior = DEFAULT
fluid = Domain2D - solid interior
```

所以：

```text
geometric Inside  -> solid -> Empty fluid cell
geometric Outside -> fluid -> Full fluid cell
Intersected       -> retain exterior fluid polygon(s)
```

只有调用者明确指定 `FluidRegion2D::Interior` 时，才把闭合轮廓内部作为流体。

不得通过顺/逆时针偷偷改变物理 fluid side。

## Corrected product pipeline

```text
solid BoundaryLoop + outer Domain2D
-> geometry validation
-> Cartesian domain
-> Quadtree adaptive refinement
-> 2:1 balance
-> geometric classification
-> explicit fluid-side selection (default EXTERIOR)
-> buildCutCells(): preserve every disconnected fluid component
-> physics area gate
-> global owner/neighbour topology
-> solid EmbeddedBoundary + outer DomainBoundary gate
-> small-cell detection
-> topology-safe agglomeration
-> quality report
-> VTK / CM2D / JSON export
-> independent CM2D read-back
-> thin-layer visualization
```

## 防复发硬门

默认外流 acceptance 必须同时满足：

```text
source_fluid_area = domain_area - solid_area
solid interior has no solver fluid cells
EmbeddedBoundary > 0
DomainBoundary > 0
topology audit = all zero
```

如果图片显示“物体内部铺网格，外围为空”，直接 FAIL；不要因为 CTest/quality/topology 其他指标为绿而放行。

## Multi-component Cut-cell

物理侧翻转后，强凹 solid 的 exterior fluid 在一个 leaf 中可能形成多个 disconnected components。

当前 solver API：

```cpp
buildCutCells(...)
```

必须保留所有 fluid components 并作为独立 solver cells。不得：

- 丢掉小组件；
- 只选最大组件；
- 跨 solid 造桥；
- 把真正 hole 填掉。

真正 local polygon-with-hole 当前仍显式 Unsupported，等待 refinement / hole-topology 扩展。

## 当前验证纪律

此前用户因 GitHub Actions/PR 反复收到大量邮件。当前 PR #7 已关闭。

因此在没有用户明确要求前：

- 不为了“看看绿不绿”随意重开 PR；
- 不把旧 Actions success 当 corrected-head success；
- 可以继续静态审计和代码修改；
- corrected head 全量 CI 尚未跑完时必须明确写 `PENDING`。

## 尚未完成的 solver-readiness 工作

即使 corrected head 后续全绿，也仍有：

- explicit farfield/domain extents；
- inlet/outlet/top/bottom/farfield named patches；
- multi-solid / holes / regions；
- polygon-with-holes；
- solver-native OpenFOAM/SU2/CGNS export；
- near-wall boundary-layer/hybrid layers；
- stronger mesh-quality optimization。

这些需要后续独立阶段推进，不能再次把“能画出来/拓扑自洽”说成“工业 CFD 网格器完成”。
