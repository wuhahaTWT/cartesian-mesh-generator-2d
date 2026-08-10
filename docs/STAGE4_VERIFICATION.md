# 阶段 4 验证记录：工业几何鲁棒性

日期：2026-08-09（Asia/Shanghai）

结论：**已完成并通过阶段四验收。** 生成器机器契约
`cartmesh-stage4-cutcell-v1` 记录运行内几何鲁棒性；禁网 OpenFOAM 2606
`checkMesh` 终态契约记录 `stage4Complete=true`、
`externalCfdCheckerAccepted=true` 与 `solverReadyCutCellMesh=true`。

工业几何诊断、薄壁/狭缝/小孔保护、多壳层 Cut-cell、区域与边界标识属于已通过的
几何范围；不得再把局部凸片 VTU 或仅 Cut-cell 的四面体分解称为完整求解器网格。

## 1. 验收矩阵

| 任务书要求 | 实现与机器证据 | 结论 |
|---|---|---|
| 非封闭 STL | 边界边计数、位置示例、VTP 标记、CLI 拒绝和永久回归 | 通过 |
| 非流形边/顶点 | 分别诊断，含只在顶点夹接的最小案例和位置标记 | 通过 |
| 法向不一致 | 共享边方向冲突诊断；分量有向体积和嵌套壳方向参与 Cut-cell 门禁 | 通过 |
| 退化、极瘦面、极小碎片 | 尺度相关高度判据、分量面积/体积/包围盒统计、小分量 VTP 标记 | 通过 |
| 重复、重叠三角形 | 重复面独立计数；BVH 候选加精确共面面积重叠窄相 | 通过 |
| 自相交、非邻接接触 | 真穿透、接触、合法共享边分开统计，问题三角形对可定位 | 通过 |
| 薄壁和狭窄间隙 | 0.05 薄壁（背景间距 0.1）、直径 0.12 小孔、斜向 gap 回归 | 通过 |
| 小孔 | 24 段闭合圆管案例；间隙方向上至少跨 4 个不同叶单元，孔道与外部连通 | 通过 |
| 内部空腔、多部件 | 互不相交实体、嵌套反向内壳、局部组件跨单元连接和全局 region ID | 通过 |
| 多流体区域与边界命名 | `--region-name`、`--boundary-range`；JSON、VTU/VTP 和显式体网格一致 | 通过 |
| 小 Cut-cell | 阈值、最小体积分数、位置、background ID 和 boundary ID 均写入报告 | 通过 |
| 容差与修复可见 | 显式 `--geometric-tolerance`；`toleranceActions`；默认 `geometryRepairApplied=false` | 通过 |
| 复杂公开几何 | Stanford Bunny 6,966 三角面；新内核、meshio 和 ParaView 证据全部重跑 | 通过 |
| 完整体网格与 CFD checker | 双区域薄壳的普通流体单元与 Cut-cell 凸片统一输出为 OpenFOAM `polyMesh` | 通过 |

## 2. 几何诊断方法

表面对表面问题先由三角面 BVH 做 AABB 候选筛选，再用精确窄相分类：

- `disjoint`；
- `boundary_contact`；
- `proper_intersection`；
- `coplanar_area_overlap`；
- `indeterminate`。

合法共享边不作为自交。部分重叠、自相交和非邻接接触分别统计并输出示例三角形
对。诊断 VTP 中 issue code 8、9、10、11 分别表示重叠、自相交、非邻接接触和极小
分量。局部坐标投影避免了大平移导致的共面误判。

无效几何仍被拒绝，不会通过“自动修复”掩盖。默认没有修补、删面或缝合：
`geometryRepairApplied=false`。明确指定的几何长度容差会写入 `geometricTolerance`、
`surfaceLengthTolerance` 和 `toleranceActions`。

## 3. 多壳层、区域与边界

Cut-cell 支持两类确定语义：

1. 互不相交、各自向外的实体壳；
2. 外壳向外、嵌套空腔壳反向。

混合错误方向仍拒绝。每个背景单元可有多个局部流体分量；相邻单元间使用精确公共
面多边形重叠建立 component-to-component 连接，再按确定顺序生成全局 region ID。
报告和输出保存局部 component ID、全局 region ID、region volume 和内部连接。

解析嵌套立方体结果：外部流体体积 19，内部空腔体积 0.125，总流体体积 19.125；
两个 boundary ID 的嵌入面积分别为 10 和 20。互不相交双实体的外部仍是一个全局
流体区域。

命名接口为：

```text
--boundary-range FIRST:LAST:ID:NAME
--region-name ID:NAME
```

薄壳验收中使用 `outer_wall`、`cavity_wall`、`exterior`、`cavity`，外部 checker 逐项
复核 ID、名称、面积和体积。

## 4. 薄壁、小孔、gap 与小单元

- 薄壁：壁厚 0.05，均匀背景间距 0.1。尽管壁厚小于背景单元，仍保留两个流体
  region，解析体积和边界面积一致。
- 小孔：24 段闭合圆管，内径 0.12。gap refinement 到 level 6，在孔径方向实际跨越
  至少 4 个不同叶单元；Cut-cell 后孔中心流体与外部属于同一 region。
- 斜向狭缝：所需层级使用法向上的 Cartesian 单元投影厚度计算，不再只使用 root
  最大边长；永久回归覆盖 `(1,1,1)` 法向。
- 不可满足 gap：`artifacts/stage4_unresolved_channel.json` 中 failure count 256、
  required level 8；`--strict-gaps` 返回非零，未把失败写成 pass。
- 小单元：`artifacts/stage4_small_cells.json` 在 13³ 立方体案例中记录 866 个小单元，
  最小体积分数 `0.0643564356435589`；每项带 background ID、质心、体积分数和 boundary ID。

## 5. 显式 Cut-cell 几何

每个流体凸片保存正体积、质心、六个 Cartesian 开口面、嵌入边界面和 region ID。
数值零体积片按尺度相关容差过滤，并在报告中累计
`discardedNumericalPieceCount/Volume`，因此没有静默消失。

`VTK_POLYHEDRON` 输出用于保留原始凸片与面拓扑。复杂公开曲面上，VTK 的通用凸性
validator 会把若干具有共面 arrangement 面的片判为无效，因此本阶段没有把该调试
表示包装成“通用工具完全有效”。用于独立体积验证的正式路径是
`--tetrahedra-output`：从每个凸片内部质心和各面扇形精确分解为 `VTK_TETRA(10)`，并
保留 source piece、component 和 global region ID。正式求解器体网格则使用 OpenFOAM
`polyMesh`：每个普通流体单元和闭合正体积凸片都成为体单元，公共细分后的面依次写入
`faces/owner/neighbour`，嵌入面按 boundary ID/名称写入 wall patch。

## 6. 公开复杂几何：Stanford Bunny

输入：`benchmarks/public/stage4/stanford_bunny_libigl_binary.stl`

- 来源：libigl tutorial data 中的 Stanford Bunny；固定源 commit 和原始 OFF/STL
  SHA-256 见 `benchmarks/public/stage4/stanford_bunny_provenance.json`；
- 二进制 STL：6,966 个三角面；封闭、边/顶点流形、方向一致；无重复、重叠、自交或
  非邻接接触；
- 24³ 背景网格：13,824 个单元，1,506 个 Cut-cell，1,297 个实体单元；
- 实体体积：`0.000753934231174814`；嵌入边界面积：`0.05821291897802537`；
- 全局外部流体区域：1；闭合失败、负体积、pending face、分类 conflict、共享面不匹配
  均为 0；
- 确定性结果 hash：`1148747a8bd5a315`；
- 数值零体积片：810 个，总体积 `8.524603208974674e-17`，报告中显式可见；
- 四面体输出：3,264,168 个 `VTK_TETRA`，来源凸片 292,845 个全部覆盖；
- VTK CellValidator 无效四面体数：0；
- tetra 体积和 `0.00031476495977491594`，原凸片体积
  `0.0003147649597749113`，差 `4.662069341687669e-18`。

外部机器记录：

- `artifacts/stage4_public_bunny_meshio.json`
- `artifacts/stage4_public_bunny_paraview_tetra.json`
- `artifacts/stage4_public_bunny.json`

全量几何 JSON、polyhedra VTU、tetra VTU 较大，外部验证后使用 gzip 无损压缩保留；
文件名以 `.gz` 结尾，可用 `gunzip` 恢复。原始 polyhedra 和 tetra SHA-256 也保存在
验证记录中，压缩不是删除证据语义。

## 7. 独立外部验证

项目测试之外使用 meshio/NumPy 和 ParaView/VTK：

- `tools/meshio_stage4_verify.py` 独立读取背景、boundary、显式体网格和报告；
- 复核 region/patch ID 与名称、局部/全局 component 一致性、内部面连接、体积和面积；
- 复核小单元计数和每项字段；
- 薄壳解析真值逐项匹配；
- Bunny 的 tetra 由 ParaView 6.2.0 / VTK 9.7.0 读取并运行 CellValidator。

OpenFOAM 外部 CFD 检查：

- `artifacts/stage3_openfoam_cube_checkmesh.json`：488 个完整体单元、1,956 个面，
  总体积 `0.728`，最大非正交 `0`，最大 skewness `1.30564e-14`，
  `-allTopology` 下 `Mesh OK`；
- `artifacts/stage4_openfoam_thin_shell_checkmesh.json`：1,728 个完整体单元，
  6,222 个面，两个不连通流体区域，`farfield/outer_wall/cavity_wall` 三个 patch，
  总体积 `1.457`，最大非正交 `0`，最大 skewness `0.5`，`-allTopology` 下
  `Mesh OK`。

薄壳机器证据：

- `artifacts/stage4_thin_shell_meshio.json`
- `artifacts/stage4_thin_shell_paraview.json`
- `artifacts/stage4_thin_shell_paraview.png`
- `artifacts/stage4_thin_shell_slice.png`

重叠和自交的拒绝报告与位置标记：

- `artifacts/stage4_overlap_rejected.json`
- `artifacts/stage4_overlap_rejected_geometry_issues.vtp`
- `artifacts/stage4_self_rejected.json`
- `artifacts/stage4_self_rejected_geometry_issues.vtp`

## 8. 构建与测试门禁

| 配置 | 结果 |
|---|---:|
| GCC 15.2 Release | 最终 CTest 19/19，通过，10.54 s |
| GCC 15.2 Debug | CTest 19/19，通过，48.52 s |
| Apple Clang 21 ASan+UBSan | CTest 19/19，通过，174.74 s |

Apple 平台不支持本次配置下的 LeakSanitizer，因此没有宣称 leak sanitizer 通过；
AddressSanitizer 和 UndefinedBehaviorSanitizer 已执行。

公开 Bunny 的新版阶段四报告全输出实测：

- 最终仓库报告内部墙钟：40.224466 s；同版本外部 `/usr/bin/time` 墙钟：39.48 s；
- 最终报告峰值 RSS：595,427,328 bytes；本轮重复测量最大观察值 624,164,864 bytes；
- 线程：1；
- 构建：Release，GCC 15.2；
- 硬件：MacBookAir10,1，Apple M1，8 cores，8 GB RAM。

该数字只代表 6,966 面 Bunny、24³ 背景网格的阶段四完整输出，不外推为千万级或
工业生产性能。86,482 面 Armadillo 的低分辨率试跑耗时过长并被终止，明确作为阶段六
性能债，不列为阶段四通过证据。

## 9. 复现命令

```sh
cmake --preset release
cmake --build --preset release --parallel 4
ctest --preset release --output-on-failure

python3 tools/generate_stage4_fixtures.py

./build/release/cartmesh_cutcell_cli \
  --stl benchmarks/analytic/stage4/thin_shell_wall005_ascii.stl \
  --resolution 12 \
  --boundary-range 0:12:10:outer_wall \
  --boundary-range 12:24:20:cavity_wall \
  --region-name 0:exterior --region-name 1:cavity \
  --output artifacts/stage4_thin_shell.vtu \
  --boundary-output artifacts/stage4_thin_shell_boundary.vtp \
  --polyhedra-output artifacts/stage4_thin_shell_polyhedra.vtu \
  --tetrahedra-output artifacts/stage4_thin_shell_tetrahedra.vtu \
  --geometry-output artifacts/stage4_thin_shell_geometry.json \
  --report artifacts/stage4_thin_shell.json

.venv/bin/python tools/meshio_stage4_verify.py \
  --background artifacts/stage4_thin_shell.vtu \
  --boundary artifacts/stage4_thin_shell_boundary.vtp \
  --polyhedra artifacts/stage4_thin_shell_polyhedra.vtu \
  --geometry artifacts/stage4_thin_shell_geometry.json \
  --report artifacts/stage4_thin_shell.json \
  --shape thin-shell \
  --output artifacts/stage4_thin_shell_meshio.json

./build/release/cartmesh_cutcell_cli \
  --stl benchmarks/analytic/stage4/thin_shell_wall005_ascii.stl \
  --resolution 12 \
  --boundary-range 0:12:10:outer_wall \
  --boundary-range 12:24:20:cavity_wall \
  --region-name 0:exterior --region-name 1:cavity \
  --no-vtk \
  --geometry-output artifacts/stage4_openfoam_thin_shell_geometry.json \
  --openfoam-case artifacts/stage4_openfoam_thin_shell_case \
  --report artifacts/stage4_openfoam_thin_shell.json

python3 tools/openfoam_stage3_verify.py --milestone stage4 \
  --case artifacts/stage4_openfoam_thin_shell_case \
  --project-report artifacts/stage4_openfoam_thin_shell.json \
  --output artifacts/stage4_openfoam_thin_shell_checkmesh.json
```

公开 Bunny 的固定输入来源和 hash 请以 provenance JSON 为准；不得用相似但不同版本
替换后沿用本记录的计数。

## 10. 阶段边界

- 本阶段输出包含网格几何、拓扑、区域和边界元数据，以及均匀背景的完整 OpenFOAM 体网格；
- 没有加入 CFD 求解器；
- 没有实现 CGNS、自适应 OpenFOAM 或其他阶段五通用交换接口；
- OpenFOAM `-allGeometry` 的 cell determinant/共面多邻居质量策略尚未完成；当前
  `-allTopology`、默认几何、闭合、体积、非正交和 skewness 门禁均通过；
- 没有宣称复杂 polyhedron 调试表示可直接被所有求解器消费；
- 没有进入阶段六的增量重构、并行分类或千万级生产性能优化。

阶段四所有通过条件均已有当前版本的内部和外部机器证据，可以关闭。
依照用户要求，本轮不自动启动阶段五。
