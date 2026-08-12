# Stage 6.3 验证：小 Cut-cell 与坏形状稳定化

日期：2026-08-12（Asia/Shanghai）

状态：**Stage 6.3 PASS；Stage 6 总体仍未完成；Stage 6.4 尚未开始。**

## 1. 关闭结论

reference solver-mesh 路径现在可在序列化前执行确定性控制体稳定化：

1. 仅在同一全局流体 region 内沿真实正面积内部面建立候选；
2. 先尝试保守邻接聚并，不删单元、不翻孤立面、不改质量阈值；
3. 每个候选重建完整 `OpenFoamMesh`，重新评估体积、一阶矩、patch、
   owner/neighbour、cell-edge 二流形、closure、face pyramid、non-orthogonality
   和 skewness；
4. 新的非星形/非正 face pyramid、拓扑损坏、质量退化或守恒误差都会拒绝聚并；
5. adaptive 路径被拒绝时，按来源 Morton leaf 做局部细分、恢复 2:1 平衡并
   重新 Cut；达到最大层级则显式失败；
6. 任何未解决来源都使 CLI 返回 2，且不写出声称已稳定的 `polyMesh`。

实现位于独立的 `SolverMeshStabilizer`，不在 OpenFOAM writer 内做几何修补，也没有
恢复历史失败的 writer-side 2–N merge、`convex_piece_exact` 或 kernel tetra repair。

## 2. 接口与来源可追溯性

`include/cartmesh/quality/SolverMeshStabilizer.hpp` 新增：

- `stabilize_solver_mesh()`：执行确定性邻接聚并和严格拒绝；
- `refine_stabilization_sources()`：细分可追溯 Morton 叶并恢复 2:1 平衡；
- `SolverMeshStabilizationReport`：保存动作、拒绝原因、守恒量、细化请求和未解决 ID；
- `write_solver_mesh_stabilization_json()`：写出确定性机器报告。

`OpenFoamCellSource` 现在保留 `global_region_id` 和全部 `sourceMembers`。
`cartmeshCellMapping.json` 升级为 v2，每个成员记录 background cell/stable ID、component、
local piece、region 和 background volume；聚并后仍能逆查所有来源，不把多来源
控制体伪装成单个 Cut-cell。

`MeshQualityReport` 新增每个 solver cell 的 signed volume、centroid、first moment、
surface area 和 closure ratio，稳定化器直接复用 6.2 的同一质量真相。

## 3. 确定性候选顺序和拒绝门

候选只来自真实 internal face，并按以下顺序固定：

1. 待稳定化 solver cell ID；
2. 共享面积降序；
3. 邻居体积降序；
4. 邻居稳定背景 ID 升序；
5. solver cell ID 升序。

候选必须同时满足：

- 总 signed volume 和总 first moment 在缩放容差内不变；
- `(farfield, boundary_id)` 的 patch 面数和面积不变；
- 每个 cell 的无向边恰好被两张面使用；
- topology gate 通过，新的非 tiny 质量问题数不增加；
- 聚并控制体没有非正 face pyramid/非星形问题；
- 来源 volume fraction 严格改善。

聚并之间的内部面是唯一被移除的面。其他面只重映射 owner/neighbour，然后按
OpenFOAM 约定确定排序；未被任何面引用的点才会被确定压紧。

## 4. CLI 闭环

`cartmesh_cutcell_cli` 新增：

- `--stabilize`：显式启用 6.3；未指定时保持 6.2 之前的默认行为；
- `--stabilization-output FILE`：指定稳定化 JSON；
- `--max-stabilization-rounds N`：限制 adaptive 局部细化重建轮数，默认 2。

闭环为：`build solver mesh → agglomerate/reject → refine source leaves → recut → rebuild`。
细化后还会再比较首轮与最终轮的总体积和一阶矩；不会因为进入下一轮就丢掉
跨轮守恒门。

## 5. 最小失败回归

`tests/stage63_test.cpp` 固定 7 组测试：

1. tiny sliver：一次聚并后体积/一阶矩保守，两个来源均保留；
2. T-junction：聚并后 owner/neighbour、closure 和正 face pyramid 保持；
3. 极小 volume fraction `1e-12`：聚并但不删除来源；
4. 真实 `LinearOctree + triangulated Cut-cell` 的 coarse-fine+cut 接口；
5. 聚并后非星形候选：拒绝并产生来源细化请求；
6. Morton 叶细分、2:1 恢复与 max-level 显式停止；
7. 聚并拓扑和稳定化 JSON 的确定性。

CTest 另有两个 CLI 回归：真实 adaptive L-prism 必须完成聚并；不可解 uniform
案例必须以非零状态停止。

## 6. Release 构建与项目测试

```sh
cmake --preset release
cmake --build --preset release --parallel 6
ctest --preset release --output-on-failure
```

结果：**26/26 PASS**，总测试墙钟 **18.89 s**。Stage 0–6 既有测试保持通过。

## 7. 真实 adaptive L-prism 验收

为了必然触发稳定化而不通过降低阈值逃避问题，验收案例使用更严格的
`--small-cell-threshold 0.05`；原始最小值约 `0.031746`。OpenFOAM 质量阈值未被修改。

| 指标 | 结果 |
|---|---:|
| adaptive background leaves | 2,234 |
| 初始 / 最终 solver cells | 1,954 / 1,948 |
| 保守聚并 | 6 |
| 拒绝候选 / 未解决来源 | 0 / 0 |
| 最终最小 volume fraction | `0.06349206349206304` |
| native quality issues | 0 |
| reader 非闭合 / 非正体积 / 非流形 cell-edge | 0 / 0 / 0 |
| points / faces / internal faces | 3,251 / 7,139 / 5,350 |
| OpenFOAM 2606 `checkMesh -constant -allTopology` | `Mesh OK.` |
| 两次 stabilization/quality/mapping/polyMesh | 逐字节或核心 SHA-256 相同 |

机器汇总 `artifacts/stage63/acceptance.json` 为 `status=pass`。Release 首次生成记录
`1.456675 s / 18,513,920 B peak RSS / 1 thread / GCC 15.2.0 / Release`；环境为
MacBookAir10,1、Apple M1、8 GiB、macOS 26.5.2。该轮与重复轮并行启动，数值只是环境记录，
不作性能结论。

## 8. 独立验证与复现

```sh
build/release/cartmesh_cutcell_cli \
  --stl tests/data/nonconvex_l_prism_ascii.stl \
  --adaptive --base-level 2 --max-level 6 --surface-level 4 \
  --padding-fraction 0.1 --small-cell-threshold 0.05 \
  --stabilize --max-stabilization-rounds 2 --no-vtk \
  --openfoam-case artifacts/stage63/adaptive_l_case \
  --quality-output artifacts/stage63/adaptive_l_quality.json \
  --stabilization-output artifacts/stage63/adaptive_l_stabilization.json \
  --geometry-output artifacts/stage63/adaptive_l_geometry.json \
  --report artifacts/stage63/adaptive_l_report.json

python3 tools/stage61_ascii_polymesh_verify.py \
  --case artifacts/stage63/adaptive_l_case \
  --output artifacts/stage63/adaptive_l_reader.json

python3 tools/openfoam_stage3_verify.py --milestone stage4 \
  --case artifacts/stage63/adaptive_l_case \
  --project-report artifacts/stage63/adaptive_l_report.json \
  --output artifacts/stage63/adaptive_l_checkmesh.json \
  --log-output artifacts/stage63/adaptive_l_checkmesh.log
```

`tools/verify_stage63_stabilization.py` 交叉核对稳定化报告、原生质量、独立 reader、
OpenFOAM、v2 来源映射和两次确定性。

## 9. 不可解停线证据

固定 uniform R2 cube 使用严格的测试阈值 2，8 个来源在安全聚并后仍不可解：

- CLI 返回 2；
- `pass=false`；
- 8 个 `unresolvedStableIds` 被保留；
- `completeSolverVolumeMeshWritten=false`。

这证明实现不会为了继续跑而删单元、改阈值或写出伪 PASS 网格。

## 10. 已知边界与 Stage 6.4 入口

额外运行的 `checkMesh -allTopology -allGeometry` 探针仍为 `Failed 3 mesh checks.`：

- 337 个 OpenFOAM face-plane concave cells；
- 8 张 interpolation weight `< 0.05` 的面；
- 24 张 volume ratio `< 0.01` 的面。

同一探针中 face pyramids、non-orthogonality、skewness、cell determinant、face tets 和体积
均通过。这些 `-allGeometry` 阻断属于计划中的 **Stage 6.4 复杂几何质量门**，
本阶段没有扩大范围修复，也不得因 6.3 PASS 宣称复杂几何或 `-allGeometry` 已通过。

6.3 也没有接入历史 compact/binary 千万级 writer；该统一属于 6.5。未实现任何
Stage 6.4 或 Stage 7 功能。下一步只能在用户明确确认后开始 Stage 6.4。
