# Stage 2D-V 验证记录

> **历史文档。** 下面的 `REOPENED` 结论属于 2026-08-22 物理域纠正当时的状态，
> 已由 `docs/R1F_PATCH_LOCAL_CLOSEOUT_CN.md` 第 3 节的独立 CI 运行结清：
> 当前 exterior 语义下的验收测试全部在 CI 中运行且通过。当前状态见
> `docs/CURRENT_STATE_CN.md`；本文保留原文作为纠正过程的记录。

## 当前状态

**REOPENED / OLD VISUAL ACCEPTANCE INVALIDATED**

2026-08-22 物理域审计确认：此前 2D-V 正确地“忠实绘制了核心导出的网格”，但当时核心默认物理域本身是反的——把闭合物体内部当成 fluid。

因此旧 SVG 虽然 renderer 技术上没有画错数据，却展示了错误的 CFD 产品结果。旧的 `PASS / CLOSED` 不能继续作为默认外流可视化验收。

## Renderer 设计约束仍有效

2D-V 仍然是薄层后处理：

- 只读取已导出的 `CM2D v1` / `quality.json` / `viz.json`；
- 不重新执行 geometry classification、Quadtree、Cut-cell 或 agglomeration；
- 核心 C++ library 不依赖 Python/renderer；
- 可视化失败不能改变网格结果。

这部分架构没有问题。

## 纠正后的可视化语义

默认外流的正确图必须显示：

```text
outer DomainBoundary
┌──────────────────────────┐
│  coarse Cartesian fluid  │
│       ┌──────────┐       │
│ fine  │  SOLID   │ fine  │
│ mesh  │  blank   │ mesh  │
│       └──────────┘       │
│  coarse Cartesian fluid  │
└──────────────────────────┘
       EmbeddedBoundary
```

即：

- solid interior 必须为空白/无流体 cell；
- solid 周围才是 solver mesh；
- solid boundary 附近局部加密；
- 远离 solid 的外场可逐级变粗；
- 外部计算域边界必须真实存在；
- `viz.json` 必须声明 `fluid_region=exterior` 与 `boundary_role=solid_wall`。

如果再次看到“只有物体内部有网格，物体外全空”，必须直接判定产品失败，即使 SVG 图层、quality panel、topology audit 都正常。

## 仍保留有效的 renderer 能力

- final stabilized solver polygons；
- adaptive level fill；
- internal edges；
- embedded physical boundary；
- domain boundary；
- source background-cell bounds；
- source small-cell overlay；
- optional cell ids；
- quality summary / topology audit；
- invalid-state banner。

这些功能应在 corrected external-fluid artifacts 上重新验证。

## 新 closure 条件

只有 corrected current head 实际生成默认外流 acceptance artifacts，并满足：

- rectangle/circle/concave/airfoil-like 的 solid interior 全部为空；
- 网格存在于 solid 外、Domain2D 内；
- embedded wall 与 outer domain boundary 都可见；
- source background/adaptive overlays 与 final fluid cells 不把 solid interior 误标为流体；
- renderer regression 通过；
- SVG/JSON/CM2D 来自同一 corrected exact head；

Stage 2D-V 才能重新 `PASS / CLOSED`。

旧 Actions run/artifact 仍可作为 renderer 历史功能证据，但**不得再作为物理正确性的 acceptance artifact**。

完整审计：`docs/PHYSICAL_DOMAIN_AUDIT_2026-08-22.md`
