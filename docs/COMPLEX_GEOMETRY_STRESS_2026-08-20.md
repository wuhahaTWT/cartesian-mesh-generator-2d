# Complex geometry stress check — 2026-08-20

## 2026-08-22 物理域审计注记

本文件原始 stress run 建立在旧的默认语义：`BoundaryLoop interior = fluid`。该语义现已确认不符合本项目默认外流 CFD 目标。

因此下面历史 PASS/FAIL 数字仍可用于定位 clipping/topology 鲁棒性问题，但**不能再作为当前默认 EXTERIOR CFD 产品能力结论**。corrected current head 必须重新跑全部 stress cases。

特别是物理侧翻转后，强凹固体的**外部**流体在单个 Quadtree leaf 内也可能被固体切成多个 disconnected components。对此当前 corrected solver path 已新增 `buildCutCells(...)`，允许一个 leaf 发射多个 solver cells，而不是丢片或制造假桥。

---

## 历史检查背景

This was an exploratory robustness check beyond the formal Stage 2D-0~V acceptance set.

GitHub Actions run #24 (`32336795285`) kept the then-current Stage 0~6 suite green and exercised harder single-loop geometries. Those gates were based on the old interior-fluid product definition.

## Historical Results

- `nozzle_profile.xy` — PASS under old interior-fluid semantics at maxLevel 7.
  - leaf_count = 1480
  - source_cells = 904
  - small_cells = 4
  - stabilized_cells = 900
  - vertices = 1273
  - edges = 2172
  - min_cell_area = 3.17186e-05
  - min_edge_length = 1.31417e-05
  - max_edge_aspect_ratio = 4993.64
  - max_centroid_skewness = 0.4222
  - all topology audit counters = 0
- `gear_star.xy` — historical FAIL: 220 unsupported Cut-cell leaf cases at maxLevel 8.
- `serpentine_body.xy` — historical FAIL: 1 unsupported Cut-cell leaf case at maxLevel 7.
- `naca2412_dense.xy` — historical FAIL at source global topology audit, maxLevel 8.
- `superellipse_24.xy` — historical FAIL at source global topology audit, maxLevel 7.

Later robustness work closed these old-side failures, but that later green run is also not sufficient proof for corrected external CFD semantics.

## Corrected interpretation / next stress gate

新的永久 stress suite 必须以：

```text
BoundaryLoop = solid
fluid = exterior
```

重新要求：

1. `gear_star / serpentine / nozzle / naca2412_dense / superellipse` 全部默认 EXTERIOR CLI 跑通；
2. source fluid area = outer domain area - solid polygon area；
3. solid interior 无 fluid cells；
4. EmbeddedBoundary 与 DomainBoundary 都存在；
5. multi-component external Cut-cell components 全部进入 topology；
6. topology audit 全零；
7. small-cell/agglomeration 不跨 solid wall；
8. SVG 显示物体内部为空、外围为自适应 CFD 网格；
9. 对极端 aspect ratio/skewness 继续单独报告，不用“topology valid”掩盖网格质量问题。

在 corrected stress suite 实际通过前，不描述当前实现为“任意复杂 2D 外流几何均已支持”。

完整物理域审计：`PHYSICAL_DOMAIN_AUDIT_2026-08-22.md`
