# CartMesh2D H4-2：边界层与 Cartesian/Cut-cell 共形拼接

日期：2026-08-27  
状态：H4-2 实现与验收完成；OpenFOAM 求解器质量门禁未通过，未宣称可用

## 1. 本阶段边界

本阶段把 H4-1 的固定层数、闭合外流场 boundary-layer strip 接入现有二维
Cartesian/Cut-cell 网格。未实现 local layer dropping、复杂 termination、Delaunay
transition、三维或 overset，也未改动三维核心代码。

## 2. 实现路径

1. 从每个成功的 H4-1 strip 提取最后一圈顶点，构造 `BoundaryRegion2D`
   outer envelope；该区域必须有效、逆时针、彼此不嵌套并严格位于计算域内部。
2. 用 outer envelope 作为新的固体边界，重新建立 remainder quadtree，并对每个叶子
   调用原生 `buildCutCells(..., FluidRegion2D::Exterior)`。因此 strip 占据的区域不会
   进入 remainder，和 envelope 相交的背景单元会被真实重新裁剪。
3. 将 H4-1 quad 和 remainder polygon 都转换为不可变的拓扑源，交给一次
   `buildGlobalTopology`。全局拓扑装配现在会对任意方向边（不仅水平/竖直边）执行
   共线顶点公共分割，因而 Cut-cell 在 envelope 上产生的交点也会分割 layer 外边。
4. 接口审计逐条要求 outer-envelope edge 恰有两个 owner，且必须是一侧 layer、
   一侧 remainder；接口顶点度数必须为 2，接口总长度必须与 envelope 周长闭合。
5. 重新从相同不可变源构建一次拓扑，并逐项比较坐标、ID、owner/neighbour 和连接，
   防止非确定性。任一步失败只返回空失败候选，不修改 H4-1 输入或既有网格。

## 3. 验收结果

| 样例 | layer | remainder cut | Cartesian | 一般 polygon | 接口边/点 | 面积误差 |
|---|---:|---:|---:|---:|---:|---:|
| 32 边形圆，4 层 | 128 | 140 | 208 | 92 | 168 / 168 | `-2.31e-14` |
| 温和凸 superellipse，3 层 | 72 | 148 | 276 | 63 | 171 / 171 | `-1.07e-14` |

两例均满足：

- `single_owner_interface_edges = 0`；
- `wrong_cell_pair_interface_edges = 0`；
- `non_two_valent_interface_vertices = 0`；
- 所有单元有正面积；
- 全局拓扑有效，2:1 balance violation 为 0；
- 独立 Python 读取器确认无 overlap、无 non-manifold edge，面积覆盖闭合；
- 同输入重复生成的 `.cm2d` 和 JSON 字节一致。

构建和完整回归命令：

```bash
cmake -S . -B build/release-2d -DCMAKE_BUILD_TYPE=Release
cmake --build build/release-2d -j4
ctest --test-dir build/release-2d --output-on-failure
```

结果：`97/97` 通过，其中 H4-2 专项 `8/8` 通过；H1-H3、H4-1 继续通过。

## 4. 独立产物与可视化

- `artifacts/h4_2/circle_hybrid.hybrid.{vtk,cm2d,json}`
- `artifacts/h4_2/circle_hybrid.independent.json`
- `artifacts/h4_2/circle_hybrid.{svg,png}`
- `artifacts/h4_2/superellipse_hybrid.hybrid.{vtk,cm2d,json}`
- `artifacts/h4_2/superellipse_hybrid.independent.json`
- `artifacts/h4_2/superellipse_hybrid.{svg,png}`

图中红线为 wall，蓝色带为 boundary layers，黑粗线为双方共享的 outer envelope，
橙色为 transition/Cut-cell，浅绿色为 Cartesian region。

## 5. OpenFOAM 边界

现有 writer 的输入就是统一 `TopologyMesh2D`，接口形式无需重构；H4-2 CLI 已接入可选
OpenFOAM 输出路径。但两个验收候选都未通过既有 solver-quality gate：圆形有 128 项、
superellipse 有 175 项，主要是 envelope 与笛卡尔网格任意相交形成的极小体积比、低
face weight 和超过 70 度的非正交面。本机也没有 `checkMesh`/`foamVersion`。

因此本阶段没有绕过质量门禁写出 case，也没有伪报真实 `checkMesh` 成功。后续若要进入
OpenFOAM，应单独处理 transition 小单元/质量稳定化，而不是放宽阈值或隐藏失败。

## 6. 已知限制

- 当前仅支持 H4-1 已明确支持的闭合、外流场、固定层数 strip；
- nested wall、local dropping 和复杂 termination 仍明确拒绝；
- 已验收样例不存在已知非共形接口；
- 已知失败边界是求解器质量，不是 H4-2 的几何/拓扑共形性。
