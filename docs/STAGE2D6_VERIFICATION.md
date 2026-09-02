# Stage 2D-6 验证记录

> **历史文档。** 下面的 `REOPENED` 结论属于 2026-08-22 物理域纠正当时的状态，
> 已由 `docs/R1F_PATCH_LOCAL_CLOSEOUT_CN.md` 第 3 节的独立 CI 运行结清：
> 当前 exterior 语义下的验收测试全部在 CI 中运行且通过。当前状态见
> `docs/CURRENT_STATE_CN.md`；本文保留原文作为纠正过程的记录。

## 当前状态

**REOPENED / OLD PASS INVALIDATED FOR DEFAULT EXTERNAL CFD PRODUCT**

2026-08-22 全项目物理域审计确认：此前 Stage 2D-6 的 CLI、acceptance 与可视化虽然 topology/quality/read-back 全部自洽，但默认把 `BoundaryLoop` 内部作为 fluid。对于本项目的外流 CFD 目标，这是错误物理域。

因此历史 GitHub Actions 绿色记录只能作为“旧内部流实现曾经编译通过”的证据，不能再证明当前二维产品完成。

当前 corrected head 已修改代码、测试、CLI、文档与 CI gate，但为了避免再次触发大量 GitHub Actions 邮件，**尚未重新打开 PR，也尚未对 corrected head 执行新的 Actions 全量验收**。

当前诚实状态：

```text
physical-domain semantics: CORRECTED IN CODE
default CLI fluid side: EXTERIOR
old Stage 2D-6 closure: INVALIDATED
corrected-head full compile/CTest/CLI validation: PENDING
```

## 2D-6A — Quality

质量模块本身仍保留：

- vertex / edge / cell count
- internal / boundary edge count
- source Cut-cell / Full-cell / small-cell count
- minCellArea
- minEdgeLength
- maxEdgeAspectRatio
- maxCentroidSkewness
- minCutCellAreaFraction
- Quadtree level distribution
- topology audit
- deterministic JSON report

但“quality.valid=true”不再足以通过产品验收。必须先通过物理域 gate。

## 2D-6B — Export / Read-back

仍保留：

- Legacy VTK `UNSTRUCTURED_GRID`
- deterministic `CM2D v1`
- independent `readCm2dTopology(...)`
- owner/neighbour、patch、cell loops、source ids/keys、audit round-trip

导出/read-back 证明文件结构一致，不证明网格在正确物理侧，因此必须与外流 physics gate 联合验收。

## 2D-6C — Corrected end-to-end CLI

当前 CLI：

```text
cartmesh2d_cli <boundary.xy> <output-prefix>
               [max-level=5]
               [padding-fraction=0.25]
               [small-alpha=0.10]
               [fluid-region=exterior|interior]
```

默认：

```text
boundary.xy = solid wall
fluid-region = exterior
fluid domain = padded Domain2D - solid interior
```

只有显式 `fluid-region=interior` 才生成内部流。

当前 corrected pipeline：

```text
solid boundary file
-> BoundaryLoop diagnostics / CCW normalization
-> padded outer computational domain
-> Quadtree boundary refinement
-> 2:1 balance
-> geometric Inside / Outside / Intersected classification
-> explicit physical fluid-side selection (default EXTERIOR)
-> buildCutCells(): retain every disconnected exterior fluid component
-> hard physics area gate
-> global topology
-> hard boundary-patch gate
-> small-cell analysis
-> topology-safe agglomeration
-> quality report
-> VTK + CM2D + JSON export
-> independent CM2D read-back
-> exact topology comparison
```

## 新增物理域硬门

默认外流必须满足：

```text
geometric Inside  -> solid -> no fluid cell
geometric Outside -> fluid -> full fluid cell
source_fluid_area = domain_area - solid_area
EmbeddedBoundary > 0
DomainBoundary > 0
```

任何一条失败，CLI 必须非零退出。

这类 gate 专门防止再次出现“程序内部 topology 全绿，但实际把网格铺在固体内部”的假成功。

## 多组件 Cut-cell

外流语义纠正后发现，强凹固体可能在一个 Quadtree leaf 内把 exterior fluid 切成多个 disconnected components。

当前 solver 路径新增：

```cpp
buildCutCells(...)
```

一个 leaf 可发射多个独立 `CutCell2D` solver cells，CLI 为每个 emitted component 分配唯一 deterministic source id，再进入 global topology。

这避免两种错误：

- 丢掉其中一个流体区域；
- 为了强行得到一个 polygon 而跨固体制造假桥。

真正的 local hole（固体 loop 完全被单个 leaf 包住）仍显式 Unsupported；当前 simple-polygon topology 不会偷偷填洞。

## 新 acceptance / CI gate

Stage 6 workflow 已更新为要求：

- 编译新增 `cartmesh2d_cutcell_components_tests`；
- 全量 Stage 0~6 CTest；
- acceptance `viz.json` 必须声明 `fluid_region=exterior`；
- `boundary_role=solid_wall`；
- CM2D 必须同时含 `EmbeddedBoundary` 与 `DomainBoundary`；
- CLI 自身先通过 area physics gate；
- 继续验证 VTK/CM2D/quality/viz 文件和 read-back。

## 历史绿色记录的解释

原先记录的 15/15、20/20 或 GitHub Actions success 不应删除，但含义必须改为：

> 它们证明旧实现的几何/拓扑/导出程序曾经自洽运行，不证明默认 CFD 流体域正确。

旧 artifact 中“物体内部有网格、外部为空”的图现在应视为失败样例，而不是 acceptance 结果。

## 当前仍存在的 solver-readiness 缺口

即使 corrected head 后续全绿，以下能力仍未完成：

1. 外部 `DomainBoundary` 还只有总 patch 类型，未拆分 inlet/outlet/top/bottom/farfield；
2. CLI 外域仍主要用 bbox + symmetric padding，真实外流应支持显式上下游/上下边界距离；
3. 默认 padding=0.25*span 只是算法回归尺度，不适合直接作为空气动力学远场默认；
4. 当前直接输出 VTK + CM2D，还没有 OpenFOAM/SU2/CGNS 等 solver-native case export；
5. 多个独立固体、多 hole、多区域尚未形成正式产品数据模型；
6. 单 leaf 内 polygon-with-holes 尚未支持；
7. 没有贴体边界层/各向异性 near-wall layers；
8. 网格质量优化仍需要继续推进。

因此后续不能再把“Stage 0~6 测试全绿”等同于“工业 CFD 网格器完成”。

## 下一次 closure 条件

只有 corrected current head 实际完成并通过：

- 根工程 configure/build；
- 全量 CTest；
- rectangle/circle/concave/airfoil-like 默认 EXTERIOR e2e；
- gear/serpentine/nozzle/NACA/superellipse 外流 stress；
- physical area gates；
- solid wall + outer domain patch gates；
- small-cell/agglomeration external-domain regression；
- VTK/CM2D/JSON export + independent read-back；
- renderer 输出明确显示 solid interior blank；

Stage 2D-6 才能重新标记 `PASS / CLOSED`。

完整审计见：`docs/PHYSICAL_DOMAIN_AUDIT_2026-08-22.md`
