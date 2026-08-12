# Cartesian 网格生成器

这是三维 Cartesian 网格生成器核心。阶段 0 至 5 已通过验证门禁。当前实现包含
ASCII/二进制 STL 导入、几何诊断、三角面 BVH、精确 triangle-AABB 分类，以及确定性
Morton 线性八叉树、表面/距离/法向变化/狭缝/用户区域细化、面 2:1 平衡、跨层邻居、自适应
VTU 和千万叶实测；阶段三还包含凸/非凸单连通封闭 STL 的真实 Cut-cell 多面体、
体积/质心/面几何、边界 patch、流体分量和 cell-face-neighbor 拓扑，
以及包含普通全流体单元与 Cut-cell 凸片的完整 OpenFOAM `polyMesh`。项目仍有意
**不实现** CFD 求解器、GUI、云服务或部署路径。此前 Stage 6 尝试的千万级完整导出、
全量独立读取和性能门已通过，但 OpenFOAM 默认质量门失败，因此不能把它描述为
solver-ready。旧尝试及失败证据见 [`docs/STAGE6_PLAN.md`](docs/STAGE6_PLAN.md) 与
[`docs/STAGE6_VERIFICATION.md`](docs/STAGE6_VERIFICATION.md)；当前唯一执行计划是
[`STAGE6_REVISED_PLAN_CN.md`](STAGE6_REVISED_PLAN_CN.md)。新机器或新 Codex 会话请先读
[`WINDOWS_CODEX_HANDOFF_CN.md`](WINDOWS_CODEX_HANDOFF_CN.md)。

项目范围以
[`CARTESIAN_MESH_GENERATOR_PROJECT_BRIEF_CN.md`](CARTESIAN_MESH_GENERATOR_PROJECT_BRIEF_CN.md)
为准。

## 构建与测试

```sh
cmake --preset release
cmake --build --preset release
ctest --preset release
```

## 生成可检查的网格

```sh
./build/release/cartmesh_cli --case cube --resolution 24 \
  --output artifacts/cube_24.vtu --report artifacts/cube_24.json

./build/release/cartmesh_cli --case sphere --resolution 48 \
  --output artifacts/sphere_48.vtu --report artifacts/sphere_48.json
```

球体文件是带单元中心内外标签的均匀背景网格。它是阶段 0 的分类收敛案例，
不是 Cut-cell 网格，也不是可供求解器使用的曲面边界几何。

阶段 1 STL 路径：

```sh
./build/release/cartmesh_cli \
  --stl tests/data/closed_unit_cube_ascii.stl --resolution 32 \
  --output artifacts/stage1_cube_32_intersections.vtu \
  --report artifacts/stage1_cube_32_intersections.json
```

输入若存在孔洞、非流形边/顶点、方向冲突、退化面或重复面，CLI 会保存诊断 JSON 和
带位置/问题代码的 VTP，再拒绝继续分类。`stl_cell_classification` 的 `0/1/2/3` 分别
表示 `outside/inside/intersected/conflict`。`intersected` 是真实三角形与闭单元 AABB
的精确 SAT 命中，不是中心点落在表面；它仍只是背景标签，不是 Cut-cell。

使用以下命令运行独立的 XML/VTK 结构检查：

```sh
python3 tools/verify_vtu.py artifacts/cube_24.vtu
python3 tools/verify_vtu.py artifacts/sphere_48.vtu

.venv/bin/python tools/meshio_stage1_verify.py \
  tests/data/closed_unit_cube_ascii.stl \
  artifacts/stage1_cube_32_intersections.vtu \
  artifacts/stage1_cube_32_intersections.json \
  --axis-aligned-cube 0 0 0 1 1 1
```

安装 ParaView 后，可以使用其 `pvpython` 运行时执行外部读取和渲染：

```sh
pvpython tools/paraview_stage1_verify.py \
  --surface tests/data/closed_unit_cube_ascii.stl \
  --mesh artifacts/stage1_cube_32_intersections.vtu \
  --report artifacts/stage1_cube_32_paraview_intersections.json \
  --overview artifacts/stage1_cube_32_intersected_paraview.png \
  --slice artifacts/stage1_cube_32_slice_paraview.png \
  --slice-origin 0.5 0.5 0.51 --slice-normal 0 0 1
```

阶段 0 证据见 [`VERIFICATION.md`](VERIFICATION.md)，阶段 1 证据见
[`docs/STAGE1_VERIFICATION.md`](docs/STAGE1_VERIFICATION.md)。

## 阶段 2 自适应八叉树

```sh
./build/release/cartmesh_octree_cli \
  --stl tests/data/closed_unit_cube_ascii.stl \
  --base-level 2 --max-level 5 --surface-level 5 \
  --padding-fraction 0.5 \
  --distance 0.12:4 --curvature 40:5 \
  --box '-0.45,-0.45,-0.45,-0.25,-0.25,-0.25:5' \
  --sphere '1.3,1.3,1.3,0.12:5' \
  --cylinder '-0.3,1.25,-0.1,-0.3,1.25,0.1,0.06:5' \
  --output artifacts/stage2_cube_adaptive.vtu \
  --report artifacts/stage2_cube_adaptive.json
```

VTU 的 `octree_level`、`octree_node_code_low32/high32` 和
`stl_cell_classification` 分别保存叶层级、无损 64 位节点码和几何分类。圆柱用户区对单元
采用保守覆盖判据，可能多细化，不会把多细化说成精确 Cut-cell 交集。

独立检查：

```sh
.venv/bin/python tools/meshio_stage2_verify.py \
  artifacts/stage2_cube_adaptive.vtu artifacts/stage2_cube_adaptive.json \
  --shape cube --output artifacts/stage2_cube_adaptive_meshio.json

pvpython tools/paraview_stage2_verify.py \
  --surface tests/data/closed_unit_cube_ascii.stl \
  --mesh artifacts/stage2_cube_adaptive.vtu \
  --overview artifacts/stage2_cube_adaptive_levels.png \
  --slice artifacts/stage2_cube_adaptive_level_slice.png \
  --report artifacts/stage2_cube_adaptive_paraview.json \
  --slice-origin 0.5 0.5 0.5
```

阶段 2 的验收矩阵、外部证据和千万叶基准见
[`docs/STAGE2_VERIFICATION.md`](docs/STAGE2_VERIFICATION.md)。

gap 在当前最大层仍无法达到指定单元数时，第一版默认继续生成并写入
`pass_with_gap_resolution_warning`；需要把它作为自动化门禁时可显式加
`--strict-gaps`。

## 阶段 3 求解器可用 Cut-cell

```sh
./build/release/cartmesh_cutcell_cli \
  --stl tests/data/nonconvex_l_prism_ascii.stl --resolution 4 \
  --small-cell-threshold 0.25 \
  --output artifacts/stage3_l_prism_4.vtu \
  --boundary-output artifacts/stage3_l_prism_4_boundary.vtp \
  --polyhedra-output artifacts/stage3_l_prism_4_polyhedra.vtu \
  --geometry-output artifacts/stage3_l_prism_4_geometry.json \
  --report artifacts/stage3_l_prism_4.json
```

输出包含真实流体子体积、质心、六个 Cartesian 开口面、嵌入边界多边形、boundary
ID、显式闭合正体积凸多面体片和内部邻接。`small_cut_cell` 可直接按阈值筛选；报告
同时列出最小体积分数、小单元位置和所属 boundary ID。单次生成报告以
`stage3GeometryTopologyComplete=true` 表示内部不变量通过；`--openfoam-case`
会写出完整流体域。外部 `checkMesh` 运行后，独立验收 JSON 才写出
`stage3Complete=true`、`solverReadyCutCellMesh=true` 和
`externalCfdCheckerAccepted=true`，不由生成器自行预告外部通过；验收脚本固定启用
`-allTopology`。额外 `-allGeometry` 质量项的当前边界如实记录在阶段三/四验证文档中。

独立复核：

```sh
.venv/bin/python tools/verify_stage3_cutcell.py \
  --surface tests/data/nonconvex_l_prism_ascii.stl \
  --mesh artifacts/stage3_l_prism_4.vtu \
  --boundary artifacts/stage3_l_prism_4_boundary.vtp \
  --polyhedra artifacts/stage3_l_prism_4_polyhedra.vtu \
  --geometry artifacts/stage3_l_prism_4_geometry.json \
  --report artifacts/stage3_l_prism_4.json --shape l_prism \
  --output artifacts/stage3_l_prism_4_external.json

./build/release/cartmesh_cutcell_cli \
  --stl tests/data/closed_unit_cube_ascii.stl \
  --resolution 8 --padding-fraction 0.1 --no-vtk \
  --geometry-output artifacts/stage3_openfoam_cube_geometry.json \
  --openfoam-case artifacts/stage3_openfoam_cube_case \
  --report artifacts/stage3_openfoam_cube.json

python3 tools/openfoam_stage3_verify.py \
  --case artifacts/stage3_openfoam_cube_case \
  --project-report artifacts/stage3_openfoam_cube.json \
  --output artifacts/stage3_openfoam_cube_checkmesh.json
```

完整验收记录见 [`docs/STAGE3_PROGRESS.md`](docs/STAGE3_PROGRESS.md)。

## 阶段 4 工业几何鲁棒性

阶段四几何鲁棒性核心已经完成：包含非封闭/非流形/法向冲突/退化与极瘦面/重复面/重叠/自相交/
非邻接接触诊断，薄壁、狭缝、小孔，多壳层和嵌套空腔，多流体区域与边界命名，
小 Cut-cell 统计，以及 Stanford Bunny 的完整外部验证。默认不自动修复或删除输入
几何；显式容差和数值零体积片均记录在报告中。

复杂显式体几何同时提供原始 `VTK_POLYHEDRON` 调试表示和可由 VTK CellValidator
独立验证的局部 `VTK_TETRA` 分解；两者只覆盖 Cut-cell 凸片。正式完整体网格
是 OpenFOAM `polyMesh`。双区域薄壳案例保留 `outer_wall`/`cavity_wall`，
并由 OpenFOAM 2606 `checkMesh` 明确判为 `Mesh OK`，所以阶段四总体门禁已关闭。

完整验收矩阵、实测性能、公开数据 provenance 和复现命令见
[`docs/STAGE4_VERIFICATION.md`](docs/STAGE4_VERIFICATION.md)。

## 阶段 5 增量式局部重构

阶段五对旧、新 STL 做与三角形顺序和绕序无关的几何差分，从旧的线性八叉树局部粗化/
细化并执行 2:1 平衡闭包。未受影响的叶保留稳定 Morton 节点码并复用 Cut-cell 几何；
新建、删除或受影响的叶重新切分。产物包含增量网格、独立全量重构参考、重构掩码和
旧网格到新网格的精确流体重叠映射。

```sh
./build/release/cartmesh_incremental_cli \
  --old-stl benchmarks/analytic/stage5/local_contour_old.stl \
  --new-stl benchmarks/analytic/stage5/local_contour_new.stl \
  --max-level 4 \
  --old-output artifacts/stage5_local_contour_old.vtu \
  --new-output artifacts/stage5_local_contour_incremental.vtu \
  --full-output artifacts/stage5_local_contour_full.vtu \
  --boundary-output artifacts/stage5_local_contour_boundary.vtp \
  --geometry-output artifacts/stage5_local_contour_geometry.json \
  --mapping-output artifacts/stage5_local_contour_mapping.json \
  --report artifacts/stage5_local_contour.json

.venv/bin/python tools/meshio_stage5_verify.py \
  --old-mesh artifacts/stage5_local_contour_old.vtu \
  --new-mesh artifacts/stage5_local_contour_incremental.vtu \
  --full-mesh artifacts/stage5_local_contour_full.vtu \
  --report artifacts/stage5_local_contour.json \
  --mapping artifacts/stage5_local_contour_mapping.json \
  --output artifacts/stage5_local_contour_meshio.json
```

单个生成器报告故意保持 `internal_pass_external_pending`；只有三类案例、独立
meshio 读取和重复性能证据汇总通过后，`artifacts/stage5_acceptance.json` 才声明
`stage5Complete=true`。完整设计契约和验收记录见
[`docs/STAGE5_PLAN.md`](docs/STAGE5_PLAN.md) 与
[`docs/STAGE5_VERIFICATION.md`](docs/STAGE5_VERIFICATION.md)。

## 阶段 6 千万级完整网格

阶段 6 新增长期内存随表面邻域增长的紧凑均匀 Cut-cell 路径，以及 binary OpenFOAM
`polyMesh` 流式写出。当前 216³ Stanford Bunny 实测产物包含 10,077,696 个背景单元、
8,700,174 个实际流体控制体和 26,463,750 张面；完整文件由独立 mmap 读取器逐单元检查，
所有 cell-edge 均为二关联，两次生成的五个文件 SHA-256 完全一致。

```sh
./build/release/cartmesh_stage6_cli \
  --stl benchmarks/public/stage4/stanford_bunny_libigl_binary.stl \
  --resolution 216 \
  --openfoam-case artifacts/stage6_10m_case \
  --report artifacts/stage6_10m_export.json

.venv/bin/python tools/stage6_binary_polymesh_verify.py \
  --case artifacts/stage6_10m_case \
  --output artifacts/stage6_10m_external_reader.json \
  --diagnose-cell-edges
```

局部预览使用同一个阶段 6 紧凑状态，但它不是完整求解器网格：

```sh
./build/release/cartmesh_stage6_cli \
  --stl benchmarks/public/stage4/stanford_bunny_libigl_binary.stl \
  --resolution 24 --preview artifacts/stage6_preview.vtu \
  --report artifacts/stage6_preview_compact.json
```

当前不能宣称阶段 6 完成。OpenFOAM 2606 已完整读取中等规模 binary 网格并通过核心拓扑、
边界闭合、正面积、正体积和非正交检查，但仍报告 3,293 张错误 face pyramid 和
12,078 张高偏斜面，最终为 `Failed 2 mesh checks.`，没有 `Mesh OK.`。
机器终态见 `artifacts/stage6_acceptance.json`；完整说明见
[`docs/STAGE6_VERIFICATION.md`](docs/STAGE6_VERIFICATION.md)。阶段 7 尚未启动。

修订计划的 Stage 6.1 已于 2026-08-12 打通 reference 路径的
`adaptive LinearOctree → 完整 ASCII OpenFOAM polyMesh`。固定 adaptive cube 与非凸
L-prism 均通过独立全量读取器和 OpenFOAM 2606 `checkMesh -allTopology`，并输出稳定
Morton leaf → solver cell 映射；均匀 R8 cube 的五个核心文件 SHA-256 与 Stage 6.0
基线逐项不变。完整证据见
[`docs/STAGE6_1_VERIFICATION.md`](docs/STAGE6_1_VERIFICATION.md)。

Stage 6.2 已加入写出前原生 solver mesh 质量评估。`--quality-output FILE` 会对与
OpenFOAM writer 完全相同的内存 points/faces/owner/neighbour 计算 closure、signed volume、
face pyramid、non-orthogonality、skewness、concavity/star-shaped、重复面、tiny face/edge
和来源 volume fraction，并把失败定位回 solver cell/face 与稳定背景 ID。固定 adaptive
cube/L-prism 的原生问题数为 0，独立 reader 通过，OpenFOAM 2606 均为 `Mesh OK.`；完整
证据见 [`docs/STAGE6_2_VERIFICATION.md`](docs/STAGE6_2_VERIFICATION.md)。

Stage 6.3 已加入显式启用的写出前稳定化闭环。`--stabilize` 会优先做同 region、
正面积邻接的保守聚并；任何非星形/负 face pyramid、拓扑、edge manifold、
patch、体积或一阶矩退化都会拒绝候选。adaptive 路径随后可按来源 Morton 叶做
局部细分和 2:1 重平衡；达到最大层级或轮数上限时显式失败。真实 adaptive
L-prism 在更严格的 `0.05` 来源体积分数门下完成 6 次聚并，独立 reader 通过且
OpenFOAM 2606 `checkMesh -allTopology` 为 `Mesh OK.`。详见
[`docs/STAGE6_3_VERIFICATION.md`](docs/STAGE6_3_VERIFICATION.md)。额外 `-allGeometry` 仍有 3 类阻断，
因此复杂几何质量门、千万级 adaptive 路径和 Stage 6 总体都尚未完成。
