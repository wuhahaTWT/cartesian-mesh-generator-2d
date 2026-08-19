# Stage 2D-2 验证记录

## 状态

**IN PROGRESS — Quadtree / 2:1 balance 核心实现已完成并通过当前根目录集成测试；等待用户单独执行“验证-2”后正式关闭阶段。**

当前严格停留在 Quadtree 自适应层；未实现 2D-3 Cut-cell polygon、全局 solver topology 或可视化。

## 已实现

- `QuadtreeLeaf2D`
- 1 -> 4 leaf refinement
- `maxLevel` 硬限制
- 二维 `(level, ix, iy)` 层级坐标
- deterministic Morton-path leaf key
- deterministic leaf ID / ordering
- boundary-intersection refinement
- distance-to-boundary refinement bands
- AABB-to-boundary conservative distance
- 基于整数最大层级 lattice 的 face-neighbor discovery
- 只把共享正长度边的 leaf 判为 face neighbor；corner-touch 不算 face neighbor
- 2:1 face-balance violation audit
- iterative 2:1 balancing closure
- balance report：iterations / refinedLeaves / violationsBefore / violationsAfter

## 当前算法规则

1. Quadtree 从整个 `Domain2D` root 开始，每次严格 1 -> 4。
2. boundary leaf 继承 2D-1 的 `Inside / Outside / Intersected` 真值语义。
3. `Intersected` leaf 可按 `boundaryLevel` 细化。
4. 非相交 leaf 可按 boundary distance band 细化。
5. face-neighbor 不使用浮点坐标猜测，而把每个 leaf 映射到 `maxLevel` 整数 lattice 后比较 face coordinate + interval overlap。
6. 2:1 balance 只对共享正长度边的 face-neighbor 生效；角点接触不触发 balance。
7. balance 遇到 level 差 > 1 时只细化较粗 leaf，循环直到违规数为 0。

## 当前测试

- root leaf 面积等于 domain 面积
- refinement 前后 leaf 面积和保持 domain 面积
- boundary leaf 达到请求 level 4
- 远场 leaf 保持比 boundary 更粗
- distance band refinement
- max level 不被突破
- 重复运行 leaf count / key / classification 稳定
- 人工制造 level 差 > 1 的 face-neighbor
- balance 前 violation > 0
- balance 后 violation = 0
- balance 前后总面积守恒
- 所有报告的 neighbor 均共享正长度 face，不把 corner-touch 误报为 neighbor
- target level > maxLevel 显式拒绝

## 已执行的根目录集成门禁

```text
cartmesh2d_stage0_geometry_tests ... Passed
cartmesh2d_stage1_grid_tests ....... Passed
cartmesh2d_stage2_quadtree_tests ... Passed

100% tests passed, 0 tests failed out of 3
```

同时验证 `CARTMESH_BUILD_2D=OFF` 根目录配置仍 PASS。

## 尚未作为 2D-2 正式 CLOSED 的原因

当前实现已经达到 2D-2 的核心功能，但依照项目执行纪律，本文件保持 `IN PROGRESS`，等待用户显式发出“验证-2”后进行最终审计：

- leaf pair overlap = 0 的独立审计；
- 2:1 balance 完整 acceptance 对照；
- 当前 GitHub branch diff 隔离检查；
- 最终 `PASS / CLOSED` 标记。

## 明确未开始

- 2D-3 Cut-cell polygon
- embedded boundary fragment
- 2D-4 global edge/cell/neighbour topology
- 2D-5 small-cell stabilization
- 2D-6 quality/export
- 2D-V visualization
