# Stage 2D-5 验证记录

> **历史文档。** 下面的 `REOPENED` 结论属于 2026-08-22 物理域纠正当时的状态，
> 已由 `docs/R1F_PATCH_LOCAL_CLOSEOUT_CN.md` 第 3 节的独立 CI 运行结清：
> 当前 exterior 语义下的验收测试全部在 CI 中运行且通过。当前状态见
> `docs/CURRENT_STATE_CN.md`；本文保留原文作为纠正过程的记录。

## 当前状态

**ALGORITHM IMPLEMENTED / EXTERNAL-DOMAIN REVALIDATION REQUIRED**

2D-5A small-cell 检测和 2D-5B topology-safe agglomeration 的算法实现仍保留，但旧的主要 shifted-circle fixture 使用的是“circle interior = fluid”的历史语义。

2026-08-22 物理域纠正后，该 fixture 已显式改为：

`FluidRegion2D::Interior`

这样它继续只验证 small-cell / agglomeration 算法本身，而不再偷偷定义产品默认 fluid side。

## 已实现能力

- configurable `areaFractionThreshold`；
- Cut-cell alpha histogram；
- deterministic `alpha < threshold` marking；
- topology internal-edge neighbour discovery；
- stable-neighbour preference；
- deterministic tie-break；
- explicit unresolved report；
- topology-safe cell agglomeration；
- member-area / merged-area 与总流体面积守恒；
- global topology rebuild；
- unsafe merge 显式失败。

## 纠正后的产品硬门

默认 `FluidRegion2D::Exterior` 上必须重新验证：

- small-cell alpha 计算基于固体**外侧**真实流体面积；
- 多 fluid components 作为独立 solver cells 后仍能正确选邻居；
- agglomeration 不得跨越 `EmbeddedBoundary` 把两侧流体错误连通；
- agglomeration 前后 `total fluid area = domain area - solid area`；
- 聚合后固体内部仍为空；
- rebuilt topology 同时保留 embedded wall 与 outer domain boundary。

在 current corrected head 完成新的全量 CTest/复杂外流案例前，本阶段不再使用旧 “PASS/CLOSED” 作为产品完成证明。

完整审计：`docs/PHYSICAL_DOMAIN_AUDIT_2026-08-22.md`
