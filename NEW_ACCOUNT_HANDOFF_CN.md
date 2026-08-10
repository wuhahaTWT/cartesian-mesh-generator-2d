# Cartesian 网格生成器：新账号接手说明

更新时间：2026-08-10（Asia/Shanghai）

项目目录：`/Users/Zhuanz/Desktop/网络生成器`

## 1. 接手结论

这不是一个停留在阶段 0 的空项目。当前代码和验证记录表明：

- 阶段 0、1、2、3、4、5 均已完成，并已在项目文档中关闭对应门禁；
- 阶段 6 已于 2026-08-10 获用户批准并实际启动；千万级生成、完整导出、独立读取、
  确定性和资源门已通过，但 OpenFOAM 默认质量门仍失败，所以验收门禁尚未关闭；
- 阶段 7 尚未在正式任务书中定义，不能由接手账号自行假定或实现；
- 当前唯一工作范围是 `docs/STAGE6_PLAN.md` 定义的千万级完整网格；不得同时启动阶段 7。

本说明用于让另一个 Codex 账号在同一台 Mac 上继续本地项目。账号切换不会删除本目录中的代码、构建产物和验证证据，但新账号不能依赖旧对话上下文，必须以本文件和仓库内的机器证据重新建立事实基线。

## 2. 当前实测基线

本次交接前已在当前源码上执行：

```sh
cmake --build --preset release --parallel 4
ctest --preset release --output-on-failure
```

结果：Release 构建成功，CTest **21/21 通过**，最后一轮测试墙钟时间为 **7.29 s**。

阶段 6 当前机器终态为 `artifacts/stage6_acceptance.json`：

```json
{
  "status": "blocked_openfoam_quality",
  "stage6Complete": false,
  "solverReadyCutCellMesh": false,
  "externalCfdCheckerAccepted": false
}
```

不要把这个状态误读为“阶段 6 尚未开始”，也不要因为千万级独立读取通过就改写为阶段 6
完成。OpenFOAM 2606 当前完整读取并通过核心拓扑，但仍有 3,293 张错误 face pyramid 和
12,078 张高偏斜面，最终为 `Failed 2 mesh checks.`，没有 `Mesh OK.`。

阶段 5 终态文件 `artifacts/stage5_acceptance.json` 当前记录
`status=pass`、`stage5Complete=true`、三案例全量等价、精确流体重叠映射和独立 meshio 读取均通过。

阶段 4 的独立 CFD 检查终态文件：

`artifacts/stage4_openfoam_thin_shell_checkmesh.json`

其中当前记录为：

```json
{
  "status": "pass",
  "stage4Complete": true,
  "externalCfdCheckerAccepted": true,
  "solverReadyCutCellMesh": true,
  "checker": "OpenFOAM checkMesh"
}
```

注意：`artifacts/stage4_public_bunny.json` 是生成器运行期报告，其状态仍是
`geometry_pass_external_cfd_pending`，不能单独拿它宣称外部 CFD 验收完成。阶段 4 的关闭依据是生成器报告、meshio/ParaView 证据以及独立 OpenFOAM 终态记录的组合，而不是某一个 JSON 文件。

本项目最初交接时没有 Git 元数据。2026-08-11 经用户明确授权，已在当前目录初始化
本地 `main` 仓库，用于从当前接管基线开始保存后续 checkpoint。当前没有配置远端，
也没有上传任何几何、产物或项目数据。阶段 0–6 在初始化之前的开发历史无法由本仓库
追溯，仍必须结合现有源码、测试、验证文档和机器产物审计，不能虚构早期 commit 历史。

## 3. 已完成阶段及证据入口

### 阶段 0：项目骨架与解析几何

已完成 C++20/CMake、规则三维 Cartesian 网格、解析立方体/球体、VTK 输出、测试和性能记录。

主要证据：

- `VERIFICATION.md`
- `benchmarks/baselines/stage0_m1_2026-08-07.json`
- `tests/test_main.cpp`

阶段 0 的球体/立方体背景单元分类不是 Cut-cell，不得重新包装成求解器可用边界网格。

### 阶段 1：STL、BVH 与均匀分类

已完成 ASCII/二进制 STL、几何诊断、三角面 BVH、精确 triangle-AABB 分类、收敛案例、确定性检查与外部 VTK/meshio/ParaView 读取。

主要证据：

- `docs/STAGE1_VERIFICATION.md`
- `tests/stage1_test.cpp`
- `tools/meshio_stage1_verify.py`
- `tools/paraview_stage1_verify.py`

`intersected` 只表示真实三角形与背景单元 AABB 的精确命中，仍不是 Cut-cell。

### 阶段 2：自适应线性八叉树

已完成 Morton 线性八叉树、表面/距离/法向变化/狭缝/用户区域细化、面 2:1 平衡、跨层邻居、确定性输出和千万叶构建实测。

主要证据：

- `docs/STAGE2_VERIFICATION.md`
- `tests/stage2_test.cpp`
- `benchmarks/baselines/stage2_octree_10m_m1_2026-08-08.json`
- `tools/meshio_stage2_verify.py`
- `tools/paraview_stage2_verify.py`

千万叶基准是八叉树叶单元规模证据，不等于千万个求解器可用 Cut-cell。

### 阶段 3：求解器可用 Cut-cell

已完成凸/非凸单连通封闭 STL 的真实 Cut-cell 多面体、正体积和质心、Cartesian 开口面、嵌入边界面、边界 patch、流体分量、cell-face-neighbor 拓扑，以及包含普通流体单元和 Cut-cell 凸片的完整 OpenFOAM `polyMesh`。

主要证据：

- `docs/STAGE3_PROGRESS.md`
- `tests/stage3_test.cpp`
- `tools/verify_stage3_cutcell.py`
- `tools/openfoam_stage3_verify.py`
- `artifacts/stage3_openfoam_cube_checkmesh.json`

### 阶段 4：工业几何鲁棒性

已完成并关闭门禁。范围包括非封闭、非流形、方向冲突、退化/极瘦面、重复/重叠、自交、非邻接接触、薄壁、狭缝、小孔、多壳层、嵌套空腔、多流体区域、边界命名、小 Cut-cell 统计和 Stanford Bunny 公开复杂几何验证。

主要证据：

- `docs/STAGE4_PROGRESS.md`
- `docs/STAGE4_VERIFICATION.md`
- `tests/stage4_test.cpp`
- `tools/meshio_stage4_verify.py`
- `tools/paraview_stage4_tetra_verify.py`
- `artifacts/stage4_openfoam_thin_shell_checkmesh.json`
- `benchmarks/public/stage4/stanford_bunny_provenance.json`

阶段 4 的真实边界也必须保留：

- 复杂 `VTK_POLYHEDRON` 是调试表示，不宣称能被所有求解器直接消费；
- 可独立验证的局部四面体分解只覆盖 Cut-cell 凸片，不是完整流体域；
- 正式完整体网格路径是 OpenFOAM `polyMesh`；
- OpenFOAM `-allGeometry` 的 cell determinant/共面多邻居质量策略仍是后续债务；
- 86,482 面 Armadillo 的低分辨率试跑因耗时过长被终止，属于后续性能债，不是阶段 4 通过证据。

## 4. 当前核心代码地图

- 规则网格与八叉树：`src/grid/`、`include/cartmesh/grid/`
- STL 与几何诊断：`src/geometry/`、`include/cartmesh/geometry/`
- 三角面 BVH：`src/spatial/TriangleBvh.cpp`
- 分类：`src/classify/SurfaceClassifier.cpp`
- Cut-cell 几何与区域拓扑：`src/cutcell/`、`include/cartmesh/cutcell/`
- 阶段 6 紧凑全域网格：`src/scalable/`、`include/cartmesh/scalable/`
- 增量重构、几何变化集与映射：`src/incremental/`、`include/cartmesh/incremental/`
- VTK、JSON、STL、OpenFOAM 输出：`src/io/`、`include/cartmesh/io/`
- 主 CLI：`apps/cli/`、`apps/octree_cli/`、`apps/cutcell_cli/`、`apps/incremental_cli/`、`apps/stage6_cli/`
- 分阶段测试：`tests/test_main.cpp`、`tests/stage1_test.cpp`、`tests/stage2_test.cpp`、`tests/stage3_test.cpp`、`tests/stage4_test.cpp`、`tests/stage5_test.cpp`、`tests/stage6_test.cpp`
- 验收文档：`VERIFICATION.md`、`docs/STAGE*_VERIFICATION.md`、`docs/STAGE*_PROGRESS.md`

## 5. 阶段 5 范围决策与终态

正式任务书 `CARTESIAN_MESH_GENERATOR_PROJECT_BRIEF_CN.md` 将阶段 5 定义为：

> 增量式局部重构。

其验收核心是局部几何变化、稳定的未影响单元 ID、旧网格到新网格映射、增量与全量结果等价，以及孔径、孔位和局部轮廓变化三类案例。

但阶段 4 收尾文档又把以下内容称为阶段 5 工作：

- CGNS；
- 自适应 OpenFOAM；
- 更多求解器交换格式；
- OpenFOAM `-allGeometry` 相关质量策略。

阶段 5 已由用户在新账号中批准启动并已完成。范围以 `docs/STAGE5_PLAN.md` 为准，收敛为增量式局部重构；CGNS、自适应 OpenFOAM、更多交换格式和 `-allGeometry` 质量策略没有混入本阶段。

关闭证据是 `docs/STAGE5_VERIFICATION.md`、`artifacts/stage5_acceptance.json`、
`benchmarks/baselines/stage5_incremental_m5_2026-08-09.json` 和三组 `artifacts/stage5_*_meshio.json`。

无论最终如何拆分，阶段 5 都必须继续满足：

- 一次只推进当前阶段；
- 相同输入的遍历顺序、ID、报告和输出确定；
- 每个几何缺陷保留最小失败案例；
- 增量结果必须与全量重构在明确的几何、拓扑和报告契约下等价；
- 不得只比较运行时间而跳过体积、面、邻接、区域、边界和外部读取验证；
- 性能记录必须包含墙钟时间、峰值 RSS、线程数、硬件和构建类型。

## 6. 阶段 5、6、7 的当前路线

### 阶段 5：已完成

正式任务书的增量式局部重构主线已完成：

1. 定义几何变更检测与受影响空间包围范围；
2. 定义稳定背景单元 ID、Cut-cell 子片 ID 和区域/边界 ID 契约；
3. 实现必要的 2:1 平衡扩展，但不得把全域重构伪装成局部更新；
4. 输出旧网格到新网格的保留、删除、创建和一对多/多对一映射；
5. 对孔径、孔位、局部轮廓三类变化分别做全量与增量对照；
6. 记录影响范围、复用率、墙钟时间、峰值 RSS、线程数和加速比；
7. 使用 meshio 5.3.5 独立读取，并在共享 Cut-cell 修复后重跑 OpenFOAM 2606 阶段 4 回归。

三类案例的增量结果与全量重构 hash 一致，重复性能案例的增量路径中位数为 0.506252 s，
全量路径中位数为 1.426893 s，中位加速比 2.8185×；口径详见验证文档。

### 阶段 6：已启动，质量门阻断

正式目标仍是千万级求解器可用完整网格，不只是千万八叉树叶或体素标签。本轮已经实现紧凑
`CompactUniformCutCellMesh`、完整 binary OpenFOAM 流式写出、独立二进制读取器、阶段 6
原生 VTU 预览和 OpenFOAM 2606 证据。

216³ Stanford Bunny 当前实测：

- 背景单元 10,077,696；
- 实际 OpenFOAM 控制体 8,700,174；
- 点 9,053,073，面 26,463,750；
- 完整 `polyMesh` 956,918,167 字节；
- 独立读取器检查全部控制体，非二关联 cell-edge 为 0；
- 两次生成的计数、紧凑 hash、拓扑 hash 和五个文件 SHA-256 完全一致；
- 生成/导出和全量读取均低于 15 分钟，峰值 RSS 均低于 6 GiB。

但 OpenFOAM 2606 只达到 `read_and_core_topology_pass_quality_fail`：未使用点已经修复为 0，
仍有 3,293 张错误 face pyramid 和 12,078 张高偏斜面，最大 skewness 1274.65。因此不得把
`artifacts/stage6_10m_external_reader.json` 的 pass 单独解释为求解器质量通过。

当前证据入口：

- `docs/STAGE6_PLAN.md`
- `docs/STAGE6_VERIFICATION.md`
- `artifacts/stage6_acceptance.json`
- `artifacts/stage6_10m_export.json`
- `artifacts/stage6_10m_external_reader.json`
- `artifacts/stage6_10m_determinism.json`
- `artifacts/stage6_medium_checkmesh.json`
- `artifacts/stage6_medium_checkmesh.log`
- `docs/STAGE6_FAILED_BRANCH_HANDOFF_CN.md`
- `artifacts/stage6_abandoned_writer_repair_summary.json`

已止损并删除的分支是：逐凸片直接写出、写出器内 2–6 元组合聚合、内核四面体拆分。
这条路线曾把 wrong pyramids 从 379 降到 14，但仍留下约 6,807 张高偏斜面、51 张
duplicate baffle 和 4 项 OpenFOAM 失败，而且 R24 已接近两分钟，不能扩展。不要恢复
`--convex-piece-cells`、`convex_piece_exact` 或任何 writer-layer 质量搜索。

如果继续阶段 6，应先在 OpenFOAM 写出之前另行设计确定性的控制体图/稳定化层：守恒地
跨相邻背景单元聚合小 Cut-cell，保留源凸片到最终 solver cell 的稳定映射，并在序列化前
检查正体积、face pyramid、重复面、边二关联、区域和边界面积。第一关只跑合成最小案例和
R24，必须得到 OpenFOAM 2606 `Mesh OK.` 后才能扩大。不得提高阈值、翻转孤立面、删除
小单元或只看自建 reader 来关闭门禁。

### 阶段 7：当前未定义

现有正式任务书只定义到阶段 6。近壁棱柱层/各向异性层、MPI 分布式、更多产品接口、GUI 或求解器耦合都只能作为未来候选方向，不是已经批准的阶段 7。阶段 6 验收完成后，应由用户明确阶段 7 的目标、范围、外部检查器和性能门槛，再更新任务书和 `AGENTS.md`。

## 7. 下一账号接管时的执行顺序

第一轮只做接管审计并复核阶段 6 当前阻断态，不自行启动阶段 7：

```sh
cd /Users/Zhuanz/Desktop/网络生成器

sed -n '1,220p' AGENTS.md
sed -n '780,890p' CARTESIAN_MESH_GENERATOR_PROJECT_BRIEF_CN.md
sed -n '1,280p' docs/STAGE4_VERIFICATION.md
sed -n '1,320p' docs/STAGE5_VERIFICATION.md
sed -n '1,360p' docs/STAGE6_VERIFICATION.md
sed -n '1,300p' docs/STAGE6_FAILED_BRANCH_HANDOFF_CN.md

cmake --preset release
cmake --build --preset release --parallel 4
ctest --preset release --output-on-failure

jq '{status,stage4Complete,externalCfdCheckerAccepted,solverReadyCutCellMesh,checker}' \
  artifacts/stage4_openfoam_thin_shell_checkmesh.json
jq '{status,stage5Complete,incrementalResultsEquivalentToFullRebuild,externalIndependentReaderAccepted,acceptanceBlockers}' \
  artifacts/stage5_acceptance.json
jq '{status,stage6Complete,solverReadyCutCellMesh,gates,acceptanceBlockers}' \
  artifacts/stage6_acceptance.json
jq '{status,formatRead,coreTopologyPass,meshOkMarker,failedMeshCheckCount,counts,quality}' \
  artifacts/stage6_medium_checkmesh.json
```

不必在首次接管时重跑约 913 MiB 的千万级 Bunny。先核对现有 hash、全量 reader、确定性报告、
R96 OpenFOAM 原始日志和 21/21 CTest。只要阶段 6 写出器或质量策略发生变化，就必须按
`docs/STAGE6_VERIFICATION.md` 重跑 R96 OpenFOAM、千万级完整导出、全量边检查和确定性复跑。

## 8. 不可突破的项目边界

- 几何和拓扑正确性优先于界面和展示；
- 不隐藏无效几何、负体积、非闭合、分类冲突或数值丢弃；
- 不把中心采样体素、相交标签或截图称为 Cut-cell；
- 不上传几何、项目数据或其派生物到外部服务；
- 未经用户明确的后续阶段决策，不加入 CFD 求解器、GUI、云服务、AI 网格或部署；
- 阶段自身测试之外，必须有独立外部读取器或检查器证据；
- 文档、截图和单元测试不能替代真实输出文件及外部验证。

## 9. 可直接粘贴给新账号的启动提示词

```text
你现在接手本地项目 /Users/Zhuanz/Desktop/网络生成器。

先完整阅读：
1. NEW_ACCOUNT_HANDOFF_CN.md
2. AGENTS.md
3. CARTESIAN_MESH_GENERATOR_PROJECT_BRIEF_CN.md 的阶段计划
4. docs/STAGE4_VERIFICATION.md
5. docs/STAGE5_PLAN.md
6. docs/STAGE5_VERIFICATION.md
7. docs/STAGE6_PLAN.md
8. docs/STAGE6_VERIFICATION.md
9. artifacts/stage6_acceptance.json

事实基线：阶段 0–5 已完成并关闭门禁；阶段 6 已实际启动并完成千万级生成、完整导出、全量独立读取、确定性和资源验证，但 OpenFOAM 默认质量门仍有两项失败，所以 stage6Complete=false；阶段 7 尚无正式定义。不要把项目误判为仍在早期阶段，也不要同时推进多个阶段。

第一轮只做接管审计，不改实现代码：复现 Release 构建和 21 项 CTest，核对阶段 6 验收 JSON、千万级全量 reader、确定性报告、R96 OpenFOAM 结构化结果和原始日志。

请向我交付：
- 当前状态与测试结果；
- 阶段 6 哪些门已通过、哪些门仍阻断；
- 现有正确性、确定性、性能、独立 reader 和 OpenFOAM 证据是否一致；
- 继续阶段 6 质量修复的最小设计范围和验证矩阵。

阶段 6 输出 `Mesh OK.` 并关闭门禁之前不要启动阶段 7。不得上传任何几何或项目数据，不得加入 GUI、云服务、AI 网格或 CFD 求解器，不得把体素/相交标签称为 Cut-cell，也不得隐藏失败单元或质量问题。
```
