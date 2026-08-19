# cartmesh2d 开发规则

## 0. 项目身份

`cartmesh2d` 是原生二维 Cartesian / Quadtree / Cut-cell 网格生成器，是三维 `cartmesh` 的并行子项目。

**严禁把三维项目压成 `z=0` 来冒充二维实现。**

二维核心对象必须是二维对象：`Point2D`、`Segment2D`、`AABB2D`、`Polygon2D`、二维 Cartesian cell、Quadtree leaf、CutPolygon、Edge2D。

## 1. 隔离边界（硬约束）

二维任务默认只允许修改：

- `cartmesh2d/**`
- 顶层 `CMakeLists.txt` 中专门用于 `add_subdirectory(cartmesh2d)` 的接入代码
- 顶层文档中用于说明二维子项目存在的少量链接/说明

不得为了二维功能重构或修改三维核心：

- `include/cartmesh/**`
- `src/**`
- `apps/**`
- `tests/**`
- 三维 Stage 6 / Stage 7 算法

不得把三维核心改写为 `Point<Dim>`、`AABB<Dim>`、`Tree<Dim>` 等泛型架构，除非未来二维、三维均已稳定且用户单独批准工程化重构。

## 2. 开发顺序

严格执行：

`2D-0 -> 2D-1 -> 2D-2 -> 2D-3 -> 2D-4 -> 2D-5 -> 2D-6 -> 2D-V`

每次只推进用户明确批准的一个阶段或子阶段。不得因为“顺手”提前实现后续算法。

## 3. 核心真实性规则

1. 几何与拓扑正确性高于可视化。
2. 不得把“删除相交格子”“单元中心 inside/outside 采样”称为 Cut-cell。
3. `Cut-cell` 必须输出真实二维流体 polygon，并能计算正面积、质心及边界边。
4. 任何几何失败、自交、零面积、重复边、孤立边、非流形关系必须显式报错或报告，不能静默修补后宣称成功。
5. 相同输入、相同参数必须产生确定性的 cell/edge/vertex 顺序、ID、报告和输出。
6. 所有 tolerance 必须集中管理、命名清楚，不得在算法内部散落魔法常数。
7. 所有阶段必须有最小失败案例和回归测试。

## 4. 可视化禁令

在 `2D-6` 核心验收前：

- 不开发 GUI；
- 不以截图是否“好看”作为算法验收；
- `tools/visualization/` 只保留目录，不实现产品功能；
- 如调试必须输出几何，可用纯文本/JSON/VTK 调试数据，但不得让绘图代码进入核心库依赖。

## 5. 2D-0 当前起点

当用户明确批准开始 `2D-0` 后，第一批允许实现的对象只有：

- `Point2D` / `Vector2D`
- `Segment2D`
- `AABB2D`
- `Polygon2D`
- `BoundaryLoop`
- orientation / signed area / centroid
- robust segment intersection
- point-in-polygon（inside/outside/boundary）
- 闭合、退化、自交、方向等输入诊断

在 2D-0 验收前不得实现 Quadtree、Cut-cell、GUI 或求解器导出。
