# 阶段 1 验证记录：几何体素分类原型

状态：**阶段 1 经 2026-08-09 重新审计和缺陷修复后通过**

本记录对应原始规划中的“阶段一：几何体素分类原型”，覆盖：

- ASCII/二进制封闭 STL 导入；
- 封闭、定向、边流形和顶点流形诊断；
- 三角面 BVH；
- 均匀三维 Cartesian 背景网格；
- 单元 `outside/inside/intersected/conflict` 分类；
- 相交单元与分类切片可视化；
- 立方体、球体、实心圆柱和空心圆管验证；
- 百万至约千万背景单元实测；
- 独立 meshio、解析检查器、ParaView/VTK 和清理器验证。

输出仍是带表面相交标签的均匀背景网格，`solverReadyCutCellMesh=false`。它没有重建
被切多面体、体积分数、开口面积、面邻接或守恒量，**不是 Cut-cell，也不是求解器可用
边界网格**。自适应八叉树属于阶段 2，本记录不把它算入阶段 1。

## 阶段门验收矩阵

| 原始阶段一要求 | 权威证据 | 状态 |
|---|---|---|
| 读取封闭 STL | ASCII、二进制和以 `solid` 开头的二进制头单元测试；真实 CLI 输出一致性 | 通过 |
| 构建三角面 BVH | AABB 候选查询、射线、点到表面；精确 box 查询与逐三角 SAT 对照 | 通过 |
| 生成均匀三维 Cartesian 网格 | 确定性线性 ID、隐式坐标和 VTU 六面体拓扑 | 通过 |
| 单元内部、外部、相交分类 | 独立 `CellClassification`；13 轴 triangle-AABB SAT；四类计数完整分割 | 通过 |
| 切片和相交单元可视化 | ParaView/VTK 精确阈值、STL 叠加、分类切片；meshio/Matplotlib 二次渲染 | 通过 |
| 球体、立方体、圆管验证 | 单位立方体逐单元解析真值；球/实心圆柱/空心圆管体积包络收敛 | 通过 |
| 百万至千万背景单元 | 100³ 和 216³ 单线程 Release，预热一次后各三次实测 | 通过 |
| 独立外部检查 | meshio 5.3.5、NumPy 解析真值、ElementTree VTP 检查、ParaView 6.2/VTK 9.7 | 通过 |
| 确定性和失败可见 | ASCII/二进制/重复 VTU 逐字节一致；无效 STL 报告和位置 marker | 通过 |
| 清理器 | Apple Clang 21 ASan+UBSan 12/12 CTest | 通过 |

旧的 `benchmarks/baselines/stage1_m1_2026-08-08.json` 只测单元中心点，现已标记
`superseded_center_point_only`。所有阶段一性能结论以
`benchmarks/baselines/stage1_uniform_exact_m1_2026-08-09.json` 为准。

## 分类语义

点和单元分类使用不同枚举和不同字段，不能混用：

```text
PointClassification / 可选中心采样：
0 outside, 1 inside, 2 on_surface, 3 conflict

CellClassification / stl_cell_classification：
0 outside, 1 inside, 2 intersected, 3 conflict
```

处理顺序的语义为：

1. 单元 AABB 通过 BVH 做候选粗筛；
2. 对候选三角形执行 13 轴分离轴定理窄相；
3. 任意真实命中则单元为 `intersected`；
4. 仅对没有表面相交的单元，用中心点三方向奇偶射线确定 `inside/outside`；
5. 无法一致判断时保留 `conflict`，CLI 非零退出。

实现另外计算全部中心点的独立分类，只用于保留中心采样体积收敛指标。相交单元的中心
可能在内、在外或恰在表面；`on_surface` 从不再冒充 `intersected`。

triangle-AABB 使用闭集接触语义：三角形与盒面、盒边或盒顶点接触都算
`intersected`。因此表面恰好与网格面重合时，两侧闭单元都可能命中，这是保守定义，
不是 Cut-cell 重建。数值容差由局部几何尺度和最多 4 个输入坐标 ULP 构成；大平移下
仍可分辨的正间隙有永久回归测试。

## 几何诊断和失败证据

诊断报告包含：

- 退化三角形、重复三角形；
- 边界边、非流形边、方向冲突边；
- 非流形顶点 link；
- 连通分量、稳定有向体积、嵌套深度和材料体积；
- 外壳/内腔方向是否符合嵌套关系（第一版作为非阻断警告）；
- 按最长边与高度判定的近零面积三角形；
- 每类问题的三角形 ID、边端点或顶点位置示例。

无效 STL 会先写 `rejected_invalid_surface` JSON 和 sibling
`*_geometry_issues.vtp`，再非零退出。VTP 的 `issue_code` 图例为：

```text
1 degenerate triangle
2 duplicate triangle
3 boundary edge
4 non-manifold edge
5 orientation conflict
6 non-manifold vertex
7 component orientation mismatch
```

永久最小失败案例包括：

- `tests/data/open_unit_cube_missing_triangle_ascii.stl`：3 条边界边；
- `tests/data/nonmanifold_vertex_two_tetrahedra_ascii.stl`：两个闭四面体只共享一个点，
  边均为流形但共享顶点 link 由两个不相连环组成。
- `benchmarks/analytic/stage1/sliver_tetrahedron_ascii.stl`：最瘦面高度为 `1e-20`，
  诊断在分类前稳定拒绝；
- `cube_scale1e-7_ascii.stl`、`cube_shift1e9_ascii.stl` 和
  `two_cubes_one_reversed_ascii.stl`：分别保留小尺度、大平移和多分量体积相消的失败案例。

独立检查结果：

- `artifacts/stage1_open_exact_meshio_diagnostics.json`：3 个边界 marker 的位置和
  `issue_code=3` 与报告一致；
- `artifacts/stage1_nonmanifold_vertex_exact_meshio_diagnostics.json`：独立重建顶点
  link 得到 `nonManifoldVertexCount=1`，marker 位置和 `issue_code=6` 一致。

检查器按 STL 中完全相同的 `double` 坐标建立拓扑，不自动焊接近重合点。这样裂缝和
近重合接缝会显式成为边界边，而不会被隐藏修复。

## 精确相交回归

阶段一 11 组内部测试包含以下最小几何事实：

- 三角形完全在盒内、穿过盒内部、无顶点在盒内但穿过；
- 三角形与盒面、边、顶点接触；
- 明确正间隙和大平移后多个 ULP 的 near miss；
- 三角形 AABB 重叠但真实三角形与盒分离的反例；
- 6 种顶点排列、翻向、平移和尺度变化不改变几何关系；
- BVH 精确结果等于逐三角 SAT，而不是等于 AABB 粗筛；
- 中心在外但表面穿过单元；中心为 `on_surface` 但单元为 `intersected`；
- 四类计数之和严格等于单元总数；
- 三角形乱序与循环顶点重排后，诊断、分类字节和结果哈希不变。
- `1e-7` 单位立方体缩放与 `1e9` 大平移下的完整分类计数不变；
- 逐分量平移稳定体积阻止相反法向分量相消，同时不阻断不依赖法向的奇偶分类。

对单位立方体，域 `[-0.5,1.5]^3` 的 20³ 闭集对齐案例精确得到：

```text
outside=6272, inside=512, intersected=1216, conflict=0
centerInside=1000, centerOutside=7000
```

### 尺度、平移和多分量回归

| 案例 | 结果 |
|---|---|
| 边长 `1e-7` 立方体，20³ | inside=5,832，intersected=2,168，体积 `1e-21` |
| 整体平移到 `1e9`，20³ | 同一分类计数，稳定体积 1 |
| 两个外壳一正一反 | 分类通过，方向警告 1，材料体积 2，不再全局相消为 0 |
| 高度 `1e-20` 的极瘦四面体面 | 分类前拒绝，近退化面数 2，缺陷 VTP 已写出 |

小尺度和大平移立方体的 8,000 个单元均由 meshio/NumPy 按解析立方体逐单元
复核，`analyticCubeCellMatches=8000`。多分量报告由外部检查器独立重建分量、嵌套和
`materialVolume=2`。对应证据为 `artifacts/stage1_cube_scale1e-7*`、
`stage1_cube_shift1e9*`、`stage1_two_cubes_one_reversed*` 和
`stage1_sliver_tetrahedron*`。

## 独立立方体逐单元真值

最终 32³ 单位立方体使用 5% padding，包含 32,768 个六面体。项目输出与
meshio/NumPy 检查器按六个解析平面独立重算的每个单元完全一致：

```text
outside=5768
inside=21952
intersected=5048
conflict=0
centerInside=27000
analyticCubeCellMatches=32768
```

`tools/meshio_stage1_verify.py` 还独立读取 STL，重算唯一点、边、重复/退化面、边界边、
非流形边、方向冲突、顶点 link、连通分量和有向体积；不是只复述项目 JSON。证据为：

- `artifacts/stage1_cube_32_meshio_intersections.json`；
- `artifacts/stage1_cube_32_intersections.vtu`；
- `artifacts/stage1_cube_32_intersections.json`。

## 球体、圆柱和圆管收敛

解析 STL 由 `tools/generate_stage1_surfaces.py` 确定性生成，随后走正常 STL 导入、诊断、
BVH、精确相交和中心射线路径。比较对象是导入三角多面体的有向体积，避免把 STL
三角化误差混入 Cartesian 分类误差。

| 案例 | 网格 | 三角数 | 多面体体积 | 中心采样绝对误差 | `[inside, inside+intersected]` 宽度 |
|---|---:|---:|---:|---:|---:|
| 球体 | 16³ | 2,208 | 4.158971436 | 0.021200439 | 2.641203125 |
| 球体 | 40³ | 2,208 | 4.158971436 | 0.003065564 | 1.039511000 |
| 实心圆柱 | 16³ | 384 | 6.278700406 | 0.309997281 | 3.223515625 |
| 实心圆柱 | 40³ | 384 | 6.278700406 | 0.073578406 | 1.257795000 |
| 空心圆管 | 16³ | 768 | 4.709025305 | 0.341681555 | 4.076187500 |
| 空心圆管 | 40³ | 768 | 4.709025305 | 0.037215305 | 1.638461000 |

六个案例均满足：

```text
insideCount * cellVolume
<= imported STL polyhedral volume
<= (insideCount + intersectedCount) * cellVolume
```

中心采样误差和分类体积包络宽度都随加密下降。机器复核为
`artifacts/stage1_analytic_convergence_exact.json`；三个 40³ VTU 的独立 meshio 报告
分别为 `stage1_sphere_40_meshio_intersections.json`、
`stage1_cylinder_40_meshio_intersections.json` 和
`stage1_tube_40_meshio_intersections.json`。

## ParaView 与可视化证据

ParaView 6.2.0 / VTK 9.7.0 使用独立 `vtkSTLReader` 和
`vtkXMLUnstructuredGridReader` 实际读取四个案例。脚本对
`stl_cell_classification==2` 做精确阈值，阈值计数与字段计数一致；切片关闭自动三角化，
输出为 VTK_QUAD（类型 9），32³ 一层为 1,024 个截面，40³ 一层为 1,600 个截面。

| 案例 | VTU 单元 | `intersected` | ParaView 阈值选中 | 切片单元 | 机器报告 |
|---|---:|---:|---:|---:|---|
| 立方体 32³ | 32,768 | 5,048 | 5,048 | 1,024 | `artifacts/stage1_cube_32_paraview_intersections.json` |
| 球体 40³ | 64,000 | 6,248 | 6,248 | 1,600 | `artifacts/stage1_sphere_40_paraview_intersections.json` |
| 实心圆柱 40³ | 64,000 | 7,560 | 7,560 | 1,600 | `artifacts/stage1_cylinder_40_paraview_intersections.json` |
| 空心圆管 40³ | 64,000 | 9,848 | 9,848 | 1,600 | `artifacts/stage1_tube_40_paraview_intersections.json` |

每个报告记录 ParaView/VTK 版本、读取数量、图例计数、精确阈值计数、切片类型、PNG
字节数和 SHA-256。对应 `*_intersected_paraview.png` 与 `*_slice_paraview.png` 是实际
STL 叠加和分类切片。

另用 meshio 5.3.5 + Matplotlib 3.11.0 独立提取相交单元并渲染第二套证据。它不依赖
ParaView 脚本，并记录相交单元并集外表面数量、切片四类计数、STL 平面交线数量和 PNG
SHA-256。报告名为 `artifacts/stage1_*_render_meshio.json`。

## 确定性

两次 ASCII STL 运行和一次二进制 STL 运行得到：

```text
resultHashFnv1a64 = 566cbbed432fb834
VTU bytes = 3834994
VTU SHA-256 = a0a13dc7cf5711feaa86199b906de4e7434fd9dc76f36efd9500bbc9569a77c7
```

三个 VTU 逐字节相同。内部测试另外覆盖完整三角形乱序与循环顶点重排后的分类字节和
哈希不变性。机器报告为 `artifacts/stage1_determinism_exact.json`。

## 构建、测试和清理器

- GCC 15.2.0 Release：12/12 CTest；
- GCC 15.2.0 Debug：12/12 CTest；
- 阶段一测试可执行文件：11/11 子测试；
- Apple Clang 21.0.0 Debug，ASan+UBSan：12/12 CTest；
- 编译选项包含 `-Wall -Wextra -Wpedantic -Wconversion -Wshadow`，构建无警告；
- Apple 平台不支持 LeakSanitizer，清理器运行显式使用 `detect_leaks=0`，没有把这一点
  描述成泄漏检查通过。

原始日志：

- `artifacts/stage1_release_ctest.log`；
- `artifacts/stage1_debug_ctest.log`；
- `artifacts/stage1_asan_ubsan_ctest.log`。

## 百万和千万级性能

环境：MacBook Air MacBookAir10,1，Apple M1，8 核，8 GB；macOS 26.5.2 arm64；
GCC 15.2.0 Release；单线程；无清理器；12 三角形单位立方体；关闭 VTU。每个规模先
预热 1 次，再正式运行 3 次。外部父进程使用 POSIX `wait4` 记录墙钟、user/sys 和子进程
峰值 RSS；每次同时保留内部计时/RSS、解析计数和结果哈希。

| 网格 | 单元数 | outside / inside / intersected / conflict | 外部墙钟中位 | 外部墙钟范围 | 最大外部 RSS | 结果哈希 |
|---|---:|---|---:|---:|---:|---|
| 100³ | 1,000,000 | 221312 / 729000 / 49688 / 0 | 1.3977 s | 1.3901–1.4442 s | 8,781,824 B | `cbb6079ddc8a0200` |
| 216³ | 10,077,696 | 2315304 / 7529536 / 232856 / 0 | 15.0848 s | 14.2974–15.2336 s | 22,986,752 B | `d42cc61845b30f8c` |

三次正式运行的哈希分别一致，且四类计数逐项等于独立解析立方体公式。正式基线为
`benchmarks/baselines/stage1_uniform_exact_m1_2026-08-09.json`，其中包含每次运行、
硬件、系统、编译器、构建类型、线程数、输入 SHA-256、域、间距、BVH 统计、内部/外部
时间与 RSS，以及 46 个项目源文件的清单哈希
`159cf65ffdca4cb3f0fc062f53ee87f64582841ad042cd764cadf365ea5bb4ca`。

这组性能只代表小型 12 三角形 STL 的均匀背景分类，不能外推到高三角数工业表面。
10M 运行关闭 VTU，只证明核心分类路径；32³/40³ 的独立读取证据覆盖导出路径。

## 复现命令

```sh
cmake --preset release
cmake --build --preset release --parallel 4
ctest --preset release --output-on-failure

build/release/cartmesh_cli \
  --stl tests/data/closed_unit_cube_ascii.stl --resolution 32 \
  --output artifacts/stage1_cube_32_intersections.vtu \
  --report artifacts/stage1_cube_32_intersections.json

xmllint --noout artifacts/stage1_cube_32_intersections.vtu
python3 tools/verify_vtu.py artifacts/stage1_cube_32_intersections.vtu
.venv/bin/python tools/meshio_stage1_verify.py \
  tests/data/closed_unit_cube_ascii.stl \
  artifacts/stage1_cube_32_intersections.vtu \
  artifacts/stage1_cube_32_intersections.json \
  --axis-aligned-cube 0 0 0 1 1 1 \
  --output artifacts/stage1_cube_32_meshio_intersections.json

/Applications/ParaView-6.2.0-RC1.app/Contents/bin/pvpython \
  tools/paraview_stage1_verify.py \
  --surface tests/data/closed_unit_cube_ascii.stl \
  --mesh artifacts/stage1_cube_32_intersections.vtu \
  --report artifacts/stage1_cube_32_paraview_intersections.json \
  --overview artifacts/stage1_cube_32_intersected_paraview.png \
  --slice artifacts/stage1_cube_32_slice_paraview.png \
  --slice-normal 0 0 1 --slice-origin 0.5 0.5 0.51

.venv/bin/python tools/run_stage1_benchmarks.py \
  --resolutions 100 216 --warmups 1 --repeats 3 \
  --output-dir artifacts/stage1_exact_benchmarks \
  --summary benchmarks/baselines/stage1_uniform_exact_m1_2026-08-09.json \
  --hardware-model 'MacBook Air MacBookAir10,1' --hardware-chip 'Apple M1' \
  --physical-cores 8 --logical-cores 8 --memory-bytes 8589934592
```

## 已知边界与阶段结论

- 当前要求封闭、每个连通面内定向一致、边/顶点流形 STL，不自动修复几何；
- 整壳或单个分量反向作为可见警告，默认不阻断奇偶分类；报告使用嵌套深度组合每分量绝对体积；
- 任意三角形自相交、重叠壳、多区域命名和工业脏几何鲁棒性属于后续阶段；
- 表面接触单元只是背景单元标签，没有 Cut-cell 多面体或求解器拓扑；
- 没有八叉树、2:1 平衡、CFD、GUI、云服务、AI 网格或部署路径。

原始规划中阶段一的七项要求及外部验收已全部有当前源码快照证据，因此阶段一门禁关闭。
阶段 2 的重新验收单独记录在 `docs/STAGE2_VERIFICATION.md`；Cut-cell 仍严格留在阶段 3。
