# Stage 2D-2 验证记录

## 状态

**PASS / CLOSED — Quadtree 自适应细化与 2:1 balance 已完成正式验收。**

当前严格停止在 2D-2；未实现 2D-3 Cut-cell polygon、全局 solver topology 或可视化。

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

## 算法规则

1. Quadtree 从整个 `Domain2D` root 开始，每次严格 1 -> 4。
2. boundary leaf 继承 2D-1 的 `Inside / Outside / Intersected` 语义。
3. `Intersected` leaf 可按 `boundaryLevel` 细化。
4. 非相交 leaf 可按 boundary distance band 细化。
5. face-neighbor 不使用浮点 face 坐标直接比较，而把 leaf 映射到 `maxLevel` 整数 lattice 后比较 face coordinate + positive interval overlap。
6. 2:1 balance 只对共享正长度边的 face-neighbor 生效；角点接触不触发 balance。
7. level 差 > 1 时只细化较粗 leaf，迭代直到违规数为 0。

## 正式验收测试

已覆盖并通过：

- root leaf 面积等于 domain 面积
- refinement 前后 leaf 面积和保持 domain 面积
- boundary leaf 达到请求 level 4
- 远场 leaf 保持比 boundary 更粗
- distance band refinement 生效
- max level 不被突破
- 重复运行 leaf count / key / classification 稳定
- leaf key 唯一
- 每个 leaf 面积为正
- 每个 leaf 均位于 domain 内
- 任意两个不同 leaf 的正面积 overlap 数 = 0
- 人工制造 level 差 > 1 的 face-neighbor
- balance 前 violation > 0
- balance 后 violation = 0
- balance 前后总面积守恒
- 所有报告的 neighbor 均共享正长度 face，不把 corner-touch 误报为 neighbor
- target level > maxLevel 显式拒绝

## Partition / coverage 审计

正式验收增加了独立 partition audit：

1. 对所有 leaf 检查正面积；
2. 对所有 leaf 检查 bounds 落在 domain 内；
3. 对所有 key 排序并检查重复；
4. 对所有 leaf pair 检查二维内部正面积 overlap；
5. 联合 `sum(leaf area) == domain area`。

结果：

- duplicate key = 0
- positive-area overlap pair = 0
- out-of-domain leaf = 0
- non-positive-area leaf = 0
- total leaf area = domain area（tolerance 内）

因此当前 leaf partition 不存在已检测到的重叠或覆盖缺洞。

## 根目录集成门禁

从仓库根目录，`CARTMESH_BUILD_2D=ON`，构建三个二维测试 target 后执行：

```text
cartmesh2d_stage0_geometry_tests ... Passed
cartmesh2d_stage1_grid_tests ....... Passed
cartmesh2d_stage2_quadtree_tests ... Passed

100% tests passed, 0 tests failed out of 3
```

同时验证：

```text
CARTMESH_BUILD_2D=OFF root configure ... PASS
```

说明：一次尝试执行整个仓库 `build all` 时因三维大型目标编译超过当前 120 秒执行窗口而被终止；这不是测试或编译错误。正式二维门禁随后从根 CMake 构建全部 2D-0/1/2 targets 并全部通过。

## 2D-2 acceptance 对照

- [x] refine 后 leaf 面积总和等于 domain 面积
- [x] leaf 正面积 overlap = 0
- [x] leaf 均位于 domain 内
- [x] boundary zone 比远场拥有更高 level
- [x] max level 不被突破
- [x] 2:1 face-neighbor 违规数 = 0
- [x] balance 迭代收敛
- [x] 重复运行 leaf key / ID / classification 稳定
- [x] corner-touch 不误判为 face neighbor
- [x] 根目录 2D-0 / 2D-1 / 2D-2 集成测试全部通过
- [x] `CARTMESH_BUILD_2D=OFF` 根配置仍通过

## 分支隔离审计

当前 `agent/native-2d-baseline` 相对 `main` 的 2D-2 变更仍限定为：

- 根 `AGENTS.md` 的二维并行规则；
- 根 `CMakeLists.txt` 的可选二维入口；
- `cartmesh2d/**`。

未修改三维算法目录：

- `include/cartmesh/**`
- 根 `src/**`
- 根 `apps/**`
- 根 `tests/**`

## 明确未开始

- 2D-3 Cut-cell polygon
- embedded boundary fragment
- 2D-4 global edge/cell/neighbour topology
- 2D-5 small-cell stabilization
- 2D-6 quality/export
- 2D-V visualization
