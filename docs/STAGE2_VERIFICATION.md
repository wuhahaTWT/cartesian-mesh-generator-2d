# 阶段 2 自适应线性八叉树验证记录

日期：2026-08-09（Asia/Shanghai）

## 结论和范围边界

阶段 2 在 2026-08-09 重新审计并修复斜狭缝层级、配置持久化和外部检查后通过。
当前实现是紧凑、确定性的自适应 Cartesian 背景八叉树，支持：

- 带层级哨兵位的 64 位 Morton 叶节点码；
- 紧凑线性叶数组和稳定 Morton 顺序；
- 细化、八兄弟粗化、面邻居、层级统计和分区检查；
- 精确表面相交、最近表面距离、局部法向变化、狭缝和用户区域细化；
- 面 2:1 平衡；
- 自适应叶的 `outside/inside/intersected/conflict` 分类；
- 含层级、无损 Morton 高低位和分类字段的 VTU 输出；
- 千万叶构建、内存和确定性实测。

它不是 Cut-cell。与 STL 相交的叶仍然是未切割的轴对齐六面体，粗细交界可以包含悬挂点。
当前结果不具有 Cut-cell 多面体、求解器面拓扑或守恒通量几何，因此所有报告仍写入
`solverReadyCutCellMesh=false`。后续阶段三已在独立的 Cut-cell 输出路径中完成；这不会
把阶段二的纯背景八叉树制品改称为 Cut-cell。

## 验收矩阵

| 原始规划/任务书要求 | 实现与证据 | 结果 |
|---|---|---|
| 表面细化 | BVH 粗筛 + 13 轴 triangle-AABB SAT，单位立方体所有解析表面相交叶都到达目标层 | 通过 |
| 表面距离细化 | BVH 分支限界最近三角面距离；距离规则单独启用时，近点 4 层、远点 3 层 | 通过 |
| 曲率/法向变化细化 | 局部候选三角面的成对法向夹角；曲率规则单独启用时，锐角 5 层、平面中部 4 层 | 通过 |
| 狭缝与小孔保护 | 所需层级使用间隙法向上的 Cartesian 投影厚度；斜狭缝和内径 0.16 小孔都有回归 | 通过 |
| 用户 box/sphere/cylinder 细化区 | 报告保存全部坐标、半径和目标层；外部检查器从报告自动重建三个查询点 | 通过 |
| 2:1 层级平衡 | 迭代细化较粗面邻居；内部自动检查 + 外部 Python 面矩形匹配 | 通过 |
| Morton 编码无冲突 | 0..8 及最大 21 层往返、父子哨兵位、外部逐叶解码/唯一性检查 | 通过 |
| 紧凑存储与 Morton 排序 | 叶数组每叶一个 `uint64_t`；外部解码后锚点严格非降序 | 通过 |
| 细化和粗化 | 父叶替换为 8 子叶，仅完整连续 8 兄弟可粗化；净叶数与分区不变量测试 | 通过 |
| 跨层邻居 | 粗叶一个面枚举 4 个细一层邻居，反向查询回到同一粗叶 | 通过 |
| 确定性 | 单元测试字节一致；ASCII 重复和等价二进制 STL 结果哈希一致 | 通过 |
| 自适应 VTU 和切片 | meshio 5.3.5 逐叶几何检查；ParaView 6.2/VTK 9.7 层级叠加图、切片和精确相交阈值 | 通过 |
| 一千万叶内存和构建 | 10,000,005 叶 Release 三次实测，外部墙钟/RSS/硬件/线程/构建类型齐全 | 通过 |

## 数据结构和确定性

`OctreeNodeCode` 用第 `3*level` 位作层级哨兵，其低位保存当前层的 Morton 坐标。因此即使
坐标位为零，父节点和第 0 个子节点也不会冲突。支持层级 0..21，21 层哨兵位为 63。

叶子只在一个 `std::vector<uint64_t>` 中以最大层 Morton 锚点排序。不为每个叶保存父子指针或完整
AABB。单元 AABB、中心、体积和最大层 Morton 区间均由节点码和根域派生。

相同输入的所有标记集按 Morton 锚点排序，父叶总是用固定 0..7 子序替换。单位立方体的阶段 2
确定性证据：

| 输入 | 叶数 | outside / inside / intersected / conflict | 结果哈希 |
|---|---:|---|---|
| ASCII 运行 A | 13,112 | 5768 / 2296 / 5048 / 0 | `4c6c4ecc39c04fc1` |
| ASCII 运行 B | 13,112 | 5768 / 2296 / 5048 / 0 | `4c6c4ecc39c04fc1` |
| 等价二进制 STL | 13,112 | 5768 / 2296 / 5048 / 0 | `4c6c4ecc39c04fc1` |

报告为 `artifacts/stage2_determinism_ascii_a.json`、
`artifacts/stage2_determinism_ascii_b.json` 和
`artifacts/stage2_determinism_binary.json`。

## 几何细化规则

### 表面和距离

表面规则直接复用阶段 1 经验收的 BVH + 精确 triangle-AABB SAT。面、边和顶点闭集接触均算
`intersected`，因此不会把对齐表面漏掉。

BVH 最近距离用点到节点 AABB 的距离作下界，优先访问近子树并剪枝远子树。单元距离规则使用

`max(0, distance(center, surface) - halfCellDiagonal) <= bandDistance`。

这是保守判据：可能在阈值外多细化部分单元，不会因为仅看中心而漏掉已进入距离带的单元。

### 法向变化

局部法向规则在与单元对角线成比例的 BVH 邻域中查询三角面，以成对单位法向夹角与用户阈值比较。
它是任务书允许的“局部曲率/法向变化估计”，不声称是 CAD 曲面的精确主曲率。

### 狭缝、薄壁和小孔

狭缝候选先由 BVH 局部查询获取，然后检查：

1. 两三角面法向点积足够负；
2. 两三角形的精确最短距离不超过搜索半径；
3. 面间分离在两个法向上的投影满足 facing 阈值。

所需层级由根 Cartesian 单元在检测间隙法向上的投影厚度导出：

`ceil(log2((|nx|*extentX + |ny|*extentY + |nz|*extentZ) * minimumCells / gapWidth))`。

这修复了斜间隙中使用根盒最大边长会少算一层的问题。最小回归中，间距 0.13、
法向 `(1,1,1)/sqrt(3)`、至少 4 单元的双斜面从错误的 5 层修正为 6 层，
沿法向实际穿过至少 4 个不同叶单元。

若超过最大层，第一版默认继续生成，并以
`pass_with_gap_resolution_warning` 显式写入 `gapResolutionFailureCount` 和
`maximumRequiredGapLevel`。需要用于严格自动化门禁时显式加 `--strict-gaps`，此时才以
`failed_gap_resolution` 和非零退出。

`tube_outer96_inner008_ascii.stl` 的内径为 0.16。外部检查器在直径上的
`x=-0.06,-0.02,0.02,0.06` 找到 4 个不同 Morton 叶码，全部为 6 层且分类为
`outside`，直接证明小孔未被封死。STL 本身没有表达的小孔仍无法恢复；本阶段不做
CAD 特征推断。

### 用户区域

- box：单元闭 AABB 与用户 AABB 交叠；
- sphere：球心到单元 AABB 的精确距离不超过半径；
- cylinder：单元中心到有限轴段的距离加半对角线作保守覆盖。

cylinder 规则有意保守，可能在端部或盒角多细化，但不会漏掉与指定圆柱区靠近的单元。这不是圆柱与
Cut-cell 的精确布尔交。

## 立方体自适应案例与外部真值

案例参数：根域 `[-0.5,1.5]^3`，基础层 2，最大层 5，表面层 5，距离 `0.12:4`，法向夹角
`40:5`，另有各一个 box、sphere 和 cylinder 用户区。

| 指标 | 结果 |
|---|---:|
| 叶数 | 22,863 |
| 层 4 / 层 5 叶数 | 1,415 / 21,448 |
| outside / inside / intersected / conflict | 17,423 / 2,352 / 3,088 / 0 |
| 紧凑叶码字节 | 182,904 |
| 确定内部体积下界 | 0.669921875 |
| STL 多面体体积 | 1.0 |
| 内部+相交体积上界 | 1.423828125 |
| 结果哈希 | `cde13b8f72e3ce68` |

`tools/meshio_stage2_verify.py` 用 meshio 5.3.5 和 NumPy 独立执行了：

- 22,863 个六面体的 8 个唯一点、VTK 顶点顺序和正 Jacobian 检查；
- 层级哨兵位/Morton 坐标解码；
- 无损高位/低位节点码重组；
- 22,863 个节点码全部唯一；
- 每叶物理 AABB 与解码坐标逐项一致，差异 0；
- 叶体积和为 8.0，与根域体积一致；
- 单位立方体解析 AABB 四类真值逐叶比较，差异 0；
- 69,846 个跨面邻居矩形对，最大层级差 1，违反数 0；
- 检查器从 JSON 中的完整 box/sphere/cylinder 参数自动生成查询点，三点均唯一落入 5 层叶。

主案例中的几何规则也可能细化同一位置，因此另有一个只启用 box/sphere/cylinder
的隔离案例。该案例叶数 904，三个区域查询点都唯一落入 5 层叶，外部几何差异、
解析分类差异和 2:1 违反均为 0。证据为 `artifacts/stage2_user_regions_only.json`、
`.vtu` 和 `_meshio.json`。

机器证据：

- `artifacts/stage2_cube_adaptive.json`；
- `artifacts/stage2_cube_adaptive.vtu`；
- `artifacts/stage2_cube_adaptive_meshio.json`。

## 距离与曲率规则的隔离证据

为避免表面规则把距离/曲率规则的效果遮住，两个案例把 `surface-level` 固定在基础层，
并分别只启用一条规则：

- 距离案例：`(0.5,0.5,1.1)` 为 4 层，远角 `(-0.4,-0.4,-0.4)` 为 3 层；
- 曲率案例：锐角 `(0.01,0.01,0.01)` 为 5 层，平面中部 `(0.5,0.5,-0.1)` 为 4 层。

两个 VTU 都由 meshio 逐叶证明几何误差 0、解析立方体分类差异 0、2:1 违反 0。
证据为 `artifacts/stage2_distance_only*` 和 `artifacts/stage2_curvature_only*`。

## 空心圆管狭缝/孔洞案例

新增的小孔案例 `tube_outer96_inner008_ascii.stl` 含 768 三角面，外半径 1、
内半径 0.08。以 `surface-level=3`、`max-level=6`和 `gap 0.2:4` 运行后：

| 指标 | 结果 |
|---|---:|
| 叶数 | 15,520 |
| outside / inside / intersected / conflict | 1,672 / 11,896 / 1,952 / 0 |
| 狭缝规则命中 | 11,008 |
| 直径上小孔叶数 | 4 个不同 6 层 Morton 叶 |
| 小孔四点分类 | 全部 `outside` |
| 结果哈希 | `3f4fb7bfac481365` |

该配置还会对其他局部面对报告更高层需求，因而状态为
`pass_with_gap_resolution_warning`；这不会遮蔽小孔的实际四个叶证据。报告与外部检查为
`artifacts/stage2_small_hole.json`、`.vtu` 和 `_meshio.json`。

原有大孔圆管继续作为较密三角面性能参考：

`tube_outer96_inner05_ascii.stl` 含外半径 1、内半径 0.5 和两个端面环，共 768 三角面。

| 指标 | 结果 |
|---|---:|
| 叶数 | 32,096 |
| outside / inside / intersected / conflict | 18,208 / 8,736 / 5,152 / 0 |
| 狭缝规则命中叶 | 27,648 |
| 最高所需间隙层级 | 5 |
| 无法满足的间隙叶 | 0 |
| 体积下界 / STL 体积 / 上界 | 3.6855 / 4.7090253046 / 5.8590 |
| 结果哈希 | `3819c03d24b34f70` |

meshio 外部检查重组 32,096 个唯一 Morton 码，几何差异 0，93,636 个面邻居对的 2:1
违反数为 0，根域体积闭合，并将圆管轴线中心 `(0,0,0)` 唯一分类为 `outside`。证据：

- `artifacts/stage2_tube_adaptive.json`；
- `artifacts/stage2_tube_adaptive.vtu`；
- `artifacts/stage2_tube_adaptive_meshio.json`。

当前狭缝候选中的三角形对检查优先保证几何可解释性，局部候选较密时是二次成对。该圆管的
`adaptation` 实测为约 16.82 s，不将它宣称为已优化的工业狭缝搜索。后续可对相对面对建立独立空间索引，
但不影响本阶段“规则正确并显式报告”的验收。

## ParaView/VTK 外部可视化

ParaView 6.2.0 / VTK 9.7.0 在项目外部读取器中同时读取 STL 和自适应 VTU：

- 叶数 22,863，点数 25,889；
- `octree_level` 计数为层 4：1,415，层 5：21,448；
- `stl_cell_classification==2` 精确阈值选中 3,088，与报告一致；
- `z=0.5` 层级切片产生 736 个单元；
- 层级叠加图和切片 PNG 均已生成并人工检查。

证据：

- `artifacts/stage2_cube_adaptive_paraview.json`；
- `artifacts/stage2_cube_adaptive_levels.png`，SHA-256
  `7eef6ca95670a708513b498c896ef2b99f8b1e247d813d6f1019edff3943bb0a`；
- `artifacts/stage2_cube_adaptive_level_slice.png`，SHA-256
  `655699cd73d93c7206f3e10f6ff5301c2fa06adc4471402da61efe97deca490d`。

PNG 是辅助证据；实际验收依据是 ParaView 读取器计数、meshio 逐叶检查和解析立方体真值。

## 测试、警告和清理器

- GCC 15.2.0 Release：12/12 CTest；
- GCC 15.2.0 Debug：12/12 CTest；
- 阶段 2 核心测试可执行文件：10/10 子测试；
- Apple Clang 21.0.0 Debug，ASan+UBSan：12/12 CTest；
- 未解析 gap 实测：默认退出 0 且状态为警告；同配置加 `--strict-gaps`
  退出 2 且状态为 `failed_gap_resolution`；
- `-Wall -Wextra -Wpedantic -Wconversion -Wshadow` 下 GCC/Clang 构建无警告；
- Apple 平台不支持 LeakSanitizer，使用 `ASAN_OPTIONS=detect_leaks=0`，未将该结果冒充为泄漏检查。

日志：

- `artifacts/stage2_release_build.log`；
- `artifacts/stage2_release_ctest.log`；
- `artifacts/stage2_debug_build.log`；
- `artifacts/stage2_debug_ctest.log`；
- `artifacts/stage2_asan_ubsan_build.log`；
- `artifacts/stage2_asan_ubsan_ctest.log`；
- `artifacts/stage2_gap_warning.json` 与 `artifacts/stage2_gap_strict.json`。

## 千万叶内存和构建实测

环境：MacBookAir10,1，Apple M1，8 物理/8 逻辑核，8 GiB，macOS 26.5.2 arm64；GCC 15.2.0
Release；CMake 4.4.2；单线程；无清理器。

基准从层 7 的 2,097,152 叶开始，按 Morton 前缀细分 1,128,979 个叶，产生 10,000,005 叶。
只有相邻两层，因此全局层级跨度已证明面 2:1；通用深层平衡算法另由局部测试和外部面邻居匹配验证。

| 正式运行 | 外部墙钟 | user / sys | 外部最大 RSS | 结果哈希（十进制） |
|---:|---:|---:|---:|---:|
| 1 | 4.18 s | 3.78 / 0.09 s | 124,633,088 B | 10481776288188551783 |
| 2 | 4.74 s | 3.87 / 0.11 s | 125,075,456 B | 10481776288188551783 |
| 3 | 4.35 s | 3.88 / 0.11 s | 124,633,088 B | 10481776288188551783 |

摘要：

- 墙钟中位 4.35 s，范围 4.18–4.74 s；
- 最大外部 RSS 125,075,456 B；
- 紧凑叶码数组 80,000,040 B，即每叶 8 B；
- 三次分区检查和 2:1 检查全部通过；
- 三次最终叶数和结果哈希一致。

正式机器基线：
`benchmarks/baselines/stage2_octree_10m_m1_2026-08-08.json`。它包含三次内部/外部计时、
RSS、线程数、硬件、操作系统、编译器、构建类型、可执行文件 SHA-256 和 57 个源/测试/工具文件的
清单哈希
`de452cc864c64b957dbdf406f0882810270029bae0f22dedec3724bc31d75a40`。

这个基准只测紧凑叶码构建、确定性分割和平衡不变量，没有导出 10M VTU，也没有将 10M 叶宣称为
10M 个求解器可用 Cut-cell。中小规模的自适应 VTU 路径由 meshio 和 ParaView 负责外部验证。

## 阶段门禁

阶段 2 通过条件已全部满足：

1. 表面和用户区域细化有效；
2. 面 2:1 自动检查与外部检查通过；
3. Morton 编码无冲突且顺序稳定；
4. 同层和跨层邻居正确；
5. 相同输入哈希稳定；
6. 10,000,005 叶的内存和构建实测完成；
7. 工程内测试、独立 meshio/NumPy 和 ParaView/VTK 验证均通过。

允许下一步进入阶段 3 的 Cut-cell 设计与实现，但本阶段的任何文件、VTU 或基准都不得被重新命名为
Cut-cell 结果。
