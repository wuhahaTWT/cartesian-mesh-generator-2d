# Stage 2D-4 验证记录

## 当前状态

**REOPENED / HISTORICAL TOPOLOGY PASS WAS ON THE WRONG DEFAULT FLUID DOMAIN**

原 Stage 2D-4 的 vertex/edge/owner-neighbour/coarse-fine 算法本身仍有价值，但其主要自适应 circle fixture 建立在旧的“circle interior = fluid”定义上。

历史记录中最明显的错误信号是：

```text
embedded-boundary edges > 0
domain-boundary edges = 0
```

对于“圆柱放在外部计算域中的默认外流 CFD”来说，这不应该被当作 PASS；正确网格必须同时接触：

- 固体轮廓：`EmbeddedBoundary`；
- 外部计算域：`DomainBoundary`。

## 纠正后的 Stage 2D-4 硬门

默认外流 topology 必须满足：

```text
duplicateVertices = 0
duplicateEdges = 0
orphanInternalEdges = 0
nonManifoldEdges = 0
unclassifiedBoundaryEdges = 0
openCellLoops = 0
areaMismatches = 0
EmbeddedBoundary > 0
DomainBoundary > 0
topology_area = domain_area - solid_area
```

强凹固体导致同一个 Quadtree leaf 内存在多个 disconnected exterior-fluid components 时，所有组件都必须作为独立 solver cells 进入 global topology，不能丢片或跨固体造桥。

当前 CLI 已切换为 `buildCutCells(...)` solver 路径，并为每个 emitted component 分配唯一 deterministic source id 后进入 topology。

## 仍保留有效的历史算法成果

以下能力不因物理侧纠正而作废：

- deterministic global vertex deduplication；
- edge construction；
- owner / neighbour；
- coarse-fine hanging-node edge splitting；
- boundary patch classification；
- topology audit；
- deterministic IDs。

但这些能力必须在**纠正后的外部流体域**上重新全量验证后，Stage 2D-4 才能再次 CLOSED。

完整审计见：

`docs/PHYSICAL_DOMAIN_AUDIT_2026-08-22.md`
