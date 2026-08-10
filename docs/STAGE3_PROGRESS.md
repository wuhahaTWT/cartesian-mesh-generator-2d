# 阶段 3 验收：求解器可用 Cut-cell

状态：**已完成并通过阶段三验收**

阶段三从封闭、定向一致、单连通的三角曲面和 Cartesian 背景网格生成真实流体控制体
边界表示。它不是阶段一的相交标签，也不是单元中心采样。严格复审曾发现 VTK
调试输出只覆盖 Cut-cell 凸片，普通全流体单元缺席。现已新增完整 OpenFOAM
`polyMesh`：普通流体背景单元和每个闭合正体积 Cut-cell 凸片一同输出；片间面、
跨背景 Cartesian 面和嵌入 wall 均使用公共细分，并把所有分割点传播到相邻面边。

生成器不伪装外部检查结果：生成报告在 `checkMesh` 前仍写出：

```json
"status": "geometry_pass_external_cfd_pending",
"stage3GeometryTopologyComplete": true,
"externalCfdCheckerAccepted": false,
"stage3Complete": false,
"solverReadyCutCellMesh": false,
"completeSolverVolumeMeshWritten": true
```

随后 `tools/openfoam_stage3_verify.py` 在禁网 Docker 中运行本机 OpenFOAM 2606
`checkMesh`，只有出现 `Mesh OK` 才生成终态验收 JSON：

```json
"stage3Complete": true,
"externalCfdCheckerAccepted": true,
"solverReadyCutCellMesh": true
```

## 已实现功能

### 凸与通用非凸切割

- `ConvexPolyhedron` 进行确定性有向半空间裁剪，计算体积、质心、面面积、面质心、
  外法向和闭合边；计算采用局部参考点，微小尺度和 `1e9` 平移回归均通过。
- 凸 STL 路径合并共面支撑平面，精确构造 `box ∩ solid` 和流体补集凸分解。
- 通用非凸路径把定向三角面与稳定参考点组成有向四面体链，逐四面体与 AABB 裁剪，
  累加固体体积、一次矩和 Cartesian 面占据矩；有限三角片裁剪生成真实嵌入边界多边形。
- 同一背景单元内多个三角片、多个 boundary patch ID 和多个不连通流体分量均保留。
  局部三角平面排列把单元分成闭合凸区域，并用匹配内部面得到精确流体分量 ID。

### 每个流体控制体

`FluidCellGeometry` 保存：

- 流体体积、体积分数和质心；
- 六个 Cartesian 开口面的多边形环、面积、面积分数、质心和外法向；
- 嵌入边界多边形的顶点、面积、质心、流体侧外法向和 boundary ID；
- 显式闭合正体积凸多面体片及其 fluid component ID；
- 有向面积向量闭合残差和原子边链闭合结果。

STL 面与背景网格面共面时仍保留零体积接触壁面。多个流体分量不会被合并成错误拓扑。

### 邻接与自适应层级

- 均匀网格内部面只生成一次，记录两侧背景 ID、局部面号、面积、质心和法向。
- `LinearOctree` 路径已经接入同一通用切割器；2:1 粗细交界把粗面开口与多个细面
  开口逐一连接，并验证面积和一阶矩之和。
- 自适应解析回归保留多个叶层级，分区有效、2:1 平衡，并实际包含粗细层级连接；
  总体积、嵌入面积、闭合性和共享面匹配全部通过。

### 小 Cut-cell 检测

`--small-cell-threshold` 只做检测，不自动合并或改变几何：

- 报告 `minimumCutCellVolumeFraction` 和 `smallCutCellCount`；
- 逐个列出背景单元 ID、流体体积分数、质心和 boundary ID；
- 背景 VTU 写出 `small_cut_cell`，可直接在 ParaView 中筛选。

球体 8³ 外部案例在默认阈值 `0.01` 下找到 8 个小单元，最小流体体积分数
`0.002372183461084608`；立方体 4³、阈值 `0.4` 的专用案例找到 24 个小单元。

## 输出

`cartmesh_cutcell_cli` 写出：

- 背景 VTU：`fluid_volume_fraction`、`cut_cell`、`small_cut_cell`；
- 嵌入边界 VTP：真实多边形、背景单元 ID、boundary ID、面积和流体侧外法向；
- 显式多面体 VTU：仅供调试的 VTK `POLYHEDRON(42)` Cut-cell 凸片；
- 显式四面体 VTU：仅覆盖上述 Cut-cell 凸片的外部体积复核，不是完整体网格；
- OpenFOAM `constant/polyMesh`：正式完整流体体网格，包含 `points/faces/owner/neighbour/boundary`
  和最小 `system` 字典；
- 完整几何 JSON：控制体各面、边界环、嵌入边界、多面体片和全部邻接；
- 摘要 JSON：守恒量、质量计数、哈希、墙钟时间、峰值 RSS 和线程数。

显式 VTU 会剔除未被面引用的旧顶点，并按顶点数稳定分块。这个修复保留了最小回归，
也使不同面数的球体多面体混合输出能被 meshio 5.3.5 正确读取。

## 验收结果

### 项目测试

`tests/stage3_test.cpp` 共 23 项解析回归，覆盖：

- 轴对齐/倾斜平面、四面体、微小切割、大坐标平移和接触端点；
- 凸立方体、八面体、整体反向壳、共面网格壁面；
- 非凸 L 形棱柱的解析体积、质心和面积；
- 同一单元多个 patch，以及薄板造成的两个显式流体分量；
- 均匀和 2:1 自适应 cell-face-neighbor 拓扑；
- 不同顶点数 polyhedron VTU 的紧凑、稳定分块；
- 几乎封死的微小共享开口按面积和公共原点一阶矩验证，不用病态派生质心误拒；
- OpenFOAM 输出把同一控制体上的共面 STL 三角片合并为一个壁面多边形，避免产生
  仅由三角形内部对角线造成的伪凹单元；
- 球面三角片细化时体积/表面积误差单调下降，且固定 STL 总量不依赖背景分辨率。

最终构建门禁：

- GCC Release：19/19 CTest 通过，10.54 s；
- GCC Debug：19/19 CTest 通过，48.52 s；
- Apple Clang ASan+UBSan：`ASAN_OPTIONS=detect_leaks=0` 下 19/19 CTest 通过，174.74 s；
- 编译选项包含 `-Wall -Wextra -Wpedantic -Wconversion -Wshadow`，无编译警告。

### 独立 meshio/NumPy/XML 检查

`tools/verify_stage3_cutcell.py` 独立读取背景 VTU、边界 VTP、显式 polyhedron VTU、
几何 JSON 和报告。最终立方体、非凸 L 棱柱、球体和小单元案例均为 `pass`。检查包括：

- 逐单元解析体积/质心和总边界面积；
- 多边形面积、定向、面积向量和原子边链闭合；
- cell-face-neighbor 完整性和共享面一致性；
- 每个 VTK polyhedron 的正体积与边闭合；
- 小单元阈值、位置和 boundary ID。

关键证据：

- `artifacts/stage3_cube_8_external.json`：解析不匹配 0；
- `artifacts/stage3_l_prism_4_external.json`：固体体积 3、流体体积 5.064、面积 14，
  不匹配 0；
- `artifacts/stage3_sphere_8_external.json`：8394 个显式凸多面体片均可读，体积和边
  不匹配 0；
- `artifacts/stage3_small_cube_4_external.json`：24 个小单元标记与报告逐 ID 一致。

### ParaView/VTK 外部几何检查（不等同 CFD checker）

ParaView 6.2.0 / VTK 9.7.0 读取非凸 L 棱柱案例：

- 64 个背景单元、42 个体积 Cut-cell；
- 144 个嵌入边界多边形，面积 14；
- 114 个显式流体多面体片，`vtkCellValidator` 无效单元数 0；
- reader error code 全部为 0，并生成概览和 `z=0.5` 切片。

证据位于 `artifacts/stage3_l_prism_4_paraview.json`、
`artifacts/stage3_l_prism_4_paraview.png` 和 `artifacts/stage3_l_prism_4_slice.png`。

### OpenFOAM 完整体网格验收

`artifacts/stage3_openfoam_cube_case` 包含 8³ 背景上的完整外部流体域：

- 988 点、1,956 面、972 内部面、488 体单元和 2 个 boundary patch；
- 最小/最大体积 `0.00025/0.00225`，总体积 `0.728`；
- 最大非正交度 `0`，最大 skewness `1.30564e-14`；
- OpenFOAM 2606 build `_481094f-20260618` 在 `-allTopology` 下输出 `Mesh OK`；
- 禁用 Docker 网络，在临时可写副本上运行 checker，源 case 五个核心文件不被修改；
- 重复生成的 `points/faces/owner/neighbour/boundary` 逐字节一致。

终态机器证据是 `artifacts/stage3_openfoam_cube_checkmesh.json`。

### 确定性与性能

修复三角片规范顺序后，非凸 L 棱柱原始输入与倒序/循环顶点重排输入的结果哈希均为
`578c40429565faaa`；完整几何 JSON 逐字节一致。Release 单线程 8³ 正确性基准在
MacBook Air M1 8 GB 上外部墙钟约 0.08 秒、峰值 RSS 8,781,824 字节；完整环境和
计数见 `benchmarks/baselines/stage3_l_prism_m1_2026-08-09.json`。这不是阶段六并行
或千万级性能声明。

## 阶段边界

阶段三的几何、拓扑、完整体网格和外部 CFD checker 硬门禁均已通过。
当前 OpenFOAM 完整体网格输出支持均匀 Cartesian 背景；自适应八叉树的 Cut-cell
几何/拓扑已接通，但其 OpenFOAM 交换接口仍属阶段五的格式扩展，不影响本阶段
均匀网格验收结论。

`checkMesh -allGeometry` 的额外求解质量项不作为阶段三门禁：当前一层厚远场案例会因
壁面与远场之间没有第三方向内部面而报告低 cell determinant；非对齐 Cut-cell 的公共
平面被多个邻接单元细分时还会触发 OpenFOAM 的 concave-cell 质量提示。闭合、正体积、
全拓扑和默认几何检查均通过；消除这些额外质量提示需要进一步细分控制体，属于阶段五
质量策略，已在此显式保留而未伪装成通过。
