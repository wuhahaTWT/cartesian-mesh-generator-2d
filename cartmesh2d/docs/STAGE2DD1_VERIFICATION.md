# Stage 2D-D1 多环、多连通域与孔洞语义验收

日期：2026-08-24

## 1. 阶段语义

D1 新增 `BoundaryRegion2D`，用确定性的 even-odd 规则解释多个闭合环：

- nesting depth 0：区域材料；
- depth 1：孔洞；
- depth 2：孔内 island；
- 后续深度继续交替。

环输入方向不决定材料/孔洞。诊断通过后，程序按 nesting depth 自动规范成 CCW/CW 交替方向，使显式 `Interior` 的区域材料始终位于边界有向边左侧；`Exterior` 再统一反向 embedded fragment，使外部流体位于左侧。

彼此相交、重合或接触的环在网格生成前拒绝，不猜测布尔拓扑。

## 2. 输入格式

CLI 保持原单环 `x y` 文件兼容；多个环用空行分隔：

```text
x0 y0
x1 y1
...

x0_hole y0_hole
x1_hole y1_hole
...
```

注释行以 `#` 开始。机器可读 visualization sidecar 在多环时增加：

```json
"boundary_loop_count": 2,
"boundary_semantics": "even_odd"
```

## 3. 实现范围

- 几何层：多环诊断、pairwise 相交/接触拒绝、nesting depth、交替方向、even-odd 点分类、parity 面积和总包围盒。
- 空间索引：一个确定性 BVH 索引所有环，点分类使用 winding parity。
- Cartesian/Quadtree：多环分类、边界细化、距离查询和 2:1 平衡。
- Cut-cell：收集所有环的局部 fragment，按 parity 选择流体 face，支持一个 leaf 产生多个流体 component。
- 全局拓扑与聚合：接受多环物理边界，保持统一 owner/neighbour 和 embedded boundary 分类。
- CLI：多环读取、parity 物理面积硬门、环数量报告、详细 region issue。

## 4. 实际产品

### 两个独立固体，Exterior

```text
boundary_loops=2
source_cells=1078
source_fluid_area=10.7875
expected_fluid_area=10.7875
stabilized_cells=1078
vertices=1419
edges=2498
max_edge_aspect_ratio=23.25
```

### 外环减孔洞，Interior

```text
boundary_loops=2
source_cells=1068
source_fluid_area=12.76
expected_fluid_area=12.76
stabilized_cells=1068
vertices=1520
edges=2588
max_edge_aspect_ratio=5
```

### 三层嵌套 island，Interior

```text
boundary_loops=3
source_cells=1140
source_fluid_area=17
expected_fluid_area=17
stabilized_cells=1140
vertices=1681
edges=2820
max_edge_aspect_ratio=2.5
```

独立 connectivity 测试确认：annulus interior 是一个带孔连通分量；三层 nested island 保留两个互不连通的流体分量；annulus 的 Exterior complement 保留外部流体与封闭腔体两个分量，没有桥接或删除。

## 5. 验证

- Release：`53/53` CTest PASS。
- ASan + UBSan：`27/27` 2D 测试 PASS；macOS 使用 `detect_leaks=0`，因为平台不支持 LeakSanitizer。
- `meshio` 独立读取：
  - two obstacles：1,419 points / 1,078 cells；
  - annulus：1,520 points / 1,068 cells；
  - nested island：1,681 points / 1,140 cells。
- 三个 D1 产品各自连续两次生成，VTK、CM2D、quality JSON、viz JSON 均逐字节相同。
- rectangle、airfoil-like、NACA2412、superellipse 的四类单环输出与 P1 基线逐字节相同。
- 单环 circle level 14 的完整 D1 路径墙钟为 4.80 s（P1 验收为 5.27 s），多环抽象没有牺牲已验收的快速核心。
- intersecting loops CLI 负例稳定拒绝。

产品 CM2D SHA-256：

```text
two_obstacles  3ea6ddd6737082733796be160f5ee2caae39fa49a93c71b0aabd8f1e94c6442f
annulus        3f83cb90da4d0faff9c8a187ac9a140337c9908215f3306e3bcb44a95f9336b4
nested_island  61a79264452774f01e753c20e58cb730cc6dc3d4b3079c0e8ea98ab926cecfaf
```

## 6. 保留的失败边界

如果一个完整闭合环完全落在单个粗 leaf 内，当前 simple-polygon Cut-cell 无法在一个 cell 中表达 polygon-with-holes。D1 现在显式返回 `MultipleEmbeddedComponents`，要求继续细化；不再像旧路径那样根据 cell center 把整个 annulus 静默判为空或实心。

本阶段仍未提供每个输入环独立的 solver patch 名称，所有物理环目前属于 `EmbeddedBoundary`；该边界命名/solver export 属于 S1。D1 也没有消除历史高 aspect-ratio 几何，因此不能称 solver-ready。

## 7. 结论

D1 的多环布尔语义、孔洞、独立实体和断开流体分量已进入真实 Cut-cell/全局拓扑产品路径，并通过面积、connectivity、确定性、外部读取和 sanitizer 验收。下一阶段只能进入 S1 solver-ready 质量与导出。
