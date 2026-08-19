# cartmesh2d 阶段验收规范

## 通用硬门槛

每个阶段必须：

- Debug 与 Release（可用时）至少完成目标平台编译；
- 自动测试全部通过；
- 相同输入重复运行结果确定；
- 异常输入不能崩溃或静默成功；
- 失败报告能指出对象/位置/原因；
- 不依赖可视化判断正确性。

## 2D-0 验收

- rectangle/triangle/concave polygon signed area 正确；
- centroid 对解析简单形状误差在明确 tolerance 内；
- segment intersection 覆盖 crossing/touching/parallel/collinear-overlap；
- point-in-polygon 对 inside/outside/boundary 三态正确；
- 自交 loop 被拒绝；
- 零长度边被拒绝或明确诊断；
- CW/CCW 可识别并按项目规则标准化；
- 无散落 magic epsilon。

## 2D-1 验收

- `Nx*Ny` cell 数正确；
- cell 无重叠、完整覆盖 domain；
- 简单 rectangle/circle fixture 的 inside/outside/intersected 统计可复核；
- tangent、边界落在 grid line、边界穿过 cell corner 有测试；
- 相同输入 cell IDs 与 classification 稳定。

## 2D-2 验收

- refine 后 leaf 面积总和等于 domain 面积（tolerance 内）；
- leaf 不重叠；
- boundary zone 比远场拥有更高 level；
- max level 不被突破；
- 2:1 face-neighbor 违规数 = 0；
- balance 迭代收敛；
- 重复运行 leaf key/ID 不变。

## 2D-3 验收

- 所有输出 fluid polygon 面积 > minimum tolerance；
- polygon 无自交；
- polygon orientation 一致；
- `0 < area_fraction <= 1`；
- 对解析直线切矩形案例，面积和质心与解析解一致；
- embedded boundary fragment 位于输入边界上（tolerance 内）；
- unsupported multi-component cut 显式失败而不是丢失区域。

## 2D-4 验收

- duplicate vertex/edge = 0（按 canonical tolerance）；
- orphan internal edge = 0；
- non-manifold edge = 0；
- internal edge 恰有 owner + neighbour；
- boundary edge 恰有 owner + patch；
- 每个 cell 的 edge loop 闭合；
- coarse-fine 接口没有 T-junction 拓扑悬空；
- 全局 cell 面积与局部几何面积一致。

## 2D-5 验收

- small-cell 阈值行为可重复；
- alpha 统计正确；
- 处理前后总流体面积守恒在 tolerance 内；
- 聚合后 topology audit 仍通过；
- 不产生负面积、重复边或非流形单元；
- 若无法安全处理必须显式失败/保留标记。

## 2D-6 验收

- end-to-end CLI 可从边界文件生成网格；
- quality report 包含至少：cell count、cut-cell count、level distribution、min area、min area fraction、min edge length、max aspect ratio、topology errors；
- VTK/VTU 可被独立读取器成功读取；
- 导出拓扑与内存拓扑 cell/edge/vertex 数一致；
- 至少 rectangle、circle、concave geometry、airfoil-like geometry 四类 end-to-end fixture；
- 完成一份最终 acceptance report；
- 可视化仍不是验收依赖。

## 2D-V 验收

- 仅从导出文件读取；
- 不参与网格生成；
- 能区分 background/adaptive/cut/small-cell/invalid 标记；
- 可视化与核心库无反向依赖。
