# Stage 6.2 验证：原生 solver mesh 质量评估器

日期：2026-08-12（Asia/Shanghai）

状态：**Stage 6.2 PASS；Stage 6 总体仍未完成；Stage 6.3 尚未开始。**

## 1. 关闭结论

reference OpenFOAM 路径现在先形成一份确定性的内存 `OpenFoamMesh`，原生质量评估器和
ASCII `polyMesh` writer 消费同一份 points/faces/owner/neighbour/cell-source 数据。诊断
不是导出后重新读取文件，也不是与实际输出分叉的第二套拓扑。

固定 adaptive cube 与非凸 L-prism 均满足：

1. 原生质量报告 `topologyPass=true`、`qualityPass=true`、问题数 0；
2. 独立 ASCII reader 的非闭合、非正体积、坏 cell-edge 均为 0；
3. OpenFOAM 2606 `checkMesh -constant -allTopology` 输出 `Mesh OK.`；
4. 原生 cell/face 计数、最小体积与独立工具一致；
5. 原生最大 non-orthogonality 与 OpenFOAM 在显示精度内一致；
6. 两次原生质量 JSON 的 SHA-256 完全一致。

这只完成“检测和定位”。没有自动聚合、分裂、删小单元、翻面、放宽阈值或其他 6.3 修复。

## 2. 共享 solver mesh 接口

`include/cartmesh/io/OpenFoamWriter.hpp` 新增：

- `OpenFoamCellSource`：保存 background ID、稳定 Morton/linear ID、component、局部凸片、
  来源体积分数与来源类型；
- `OpenFoamFace`：保存最终焊接点 ID、owner/neighbour、boundary ID 与 farfield 属性；
- `OpenFoamMesh`：保存最终点、OpenFOAM 顺序的面、solver cell 来源映射和内部面计数；
- `build_openfoam_mesh(grid/tree, mesh, tolerance)`：uniform/adaptive 共用最终组装；
- `write_openfoam_poly_mesh(case, OpenFoamMesh, names)`：只序列化已经组装好的对象。

旧的 uniform/tree writer 重载仍保留，但内部改为“build once → serialize”。因此已有调用方
没有被迫迁移，质量评估器也不会复制 writer 的 coarse-fine、arrangement、焊点和排序逻辑。

## 3. 原生质量接口与覆盖项

`include/cartmesh/quality/SolverMeshQuality.hpp` 提供：

- `MeshQualityThresholds`；
- `MeshQualityIssueKind`、`MeshQualityIssue`；
- `MeshQualitySummary`、`MeshQualityReport`；
- `evaluate_solver_mesh_quality()`；
- `write_solver_mesh_quality_json()`。

当前覆盖：

- 非有限几何、少于三个顶点的面；
- face area、tiny/zero face、tiny/zero edge；
- duplicate face 与 baffle-like boundary duplicate；
- cell area-vector closure、signed volume、volume centroid；
- owner/neighbour 两侧 face-pyramid sign；
- internal face non-orthogonality；
- internal/boundary skewness；
- concave face、以 cell centroid 为候选核点的 star-shaped 可行性；
- 每个 solver cell 的来源 volume fraction。

每个问题保存：问题类型、solver cell ID、solver face ID（若适用）、坐标、实测值、阈值、
background cell ID、稳定背景 ID、`full_cartesian` 或 `cut_polyhedron_piece` 来源类型。问题
按类型、face ID、cell ID、稳定背景 ID 确定排序。非有限数值写为 JSON `null`，不会写出
非法的 `nan/inf`。

CLI 新增 `--quality-output FILE`。使用 `--openfoam-case` 而未显式指定质量文件时，默认写到
case 根目录的 `cartmeshQuality.json`。CLI 主报告同时记录 native quality 是否执行、
topology/quality gate、问题数、最大 non-orthogonality、最大 skewness 和报告路径。

## 4. 默认阈值

| 项 | 默认值 |
|---|---:|
| minimum face area | `1e-20` |
| minimum edge length | `1e-12` |
| minimum cell volume | `1e-18` |
| minimum face-pyramid volume | `1e-18` |
| maximum cell closure ratio | `1e-10` |
| maximum non-orthogonality | `70 deg` |
| maximum internal skewness | `4` |
| maximum boundary skewness | `20` |
| minimum volume fraction | CLI 的 `--small-cell-threshold`，默认 `0.01` |

阈值保存在每份质量 JSON 中。6.2 不通过修改阈值来修复网格。

## 5. 最小失败回归

`tests/stage62_test.cpp` 固定 7 组测试：

1. 有效单位立方体：体积 1、closure 近零、最小 face pyramid 为 `1/6`、无问题；
2. 翻转一张面：定位非闭合、负 face pyramid 和 non-star-shaped cell；
3. 零面积面、零长度边、重复 boundary face、极小体积分数；
4. 含内凹顶点的 concave face；
5. 两个斜置相邻六面体：触发 internal non-orthogonality 与 skewness；
6. 同一失败报告两次写出逐字节相同；
7. 非有限坐标显式报告，JSON 使用 `null` 而非 `nan/inf`。

失败问题的回归同时检查 solver cell/face ID、稳定背景 ID 和来源类型，不只检查总数。

## 6. 构建与项目测试

Release：

```sh
cmake --preset release
cmake --build --preset release --parallel 4
ctest --preset release --output-on-failure
```

结果：**23/23 PASS**，总测试墙钟 **13.08 s**。

Debug：

```sh
cmake --preset debug
cmake --build --preset debug --parallel 4
ctest --preset debug --output-on-failure
```

结果：**23/23 PASS**，总测试墙钟 **66.54 s**。

新增 `cartmesh_stage62_quality_tests`；既有 Stage 0–6 测试全部保持通过。

## 7. 固定 adaptive 案例

| 指标 | adaptive cube | adaptive L-prism |
|---|---:|---:|
| solver cells | 2,704 | 1,954 |
| faces / internal | 9,654 / 6,936 | 7,145 / 5,356 |
| native issue count | 0 | 0 |
| native min volume | `1.5625e-05` | `6.25e-05` |
| native min volume fraction | `0.0370370` | `0.0317460` |
| native max closure ratio | `8.78e-17` | `8.19e-17` |
| native max non-orthogonality | `35.2643897` | `54.7356103` |
| OpenFOAM max non-orthogonality | `35.2644` | `54.7356` |
| native max skewness | `1.4142136` | `3.5901099` |
| OpenFOAM max skewness | `0.5` | `0.914585` |
| independent reader | PASS | PASS |
| OpenFOAM 2606 | `Mesh OK.` | `Mesh OK.` |

原生 skewness 与 OpenFOAM 采用不同归一化，因此不宣称数值相等。固定案例中两者都把
L-prism 排在 cube 之后，阻断问题计数也同为 0；`tools/verify_stage62_quality.py` 明确按
“同向趋势”验收 skewness，而 non-orthogonality 按 `1e-3 deg` 绝对差验收。

生成阶段的墙钟/RSS 仅作本次环境记录，不据此下性能结论：cube `2.785364 s / 20,283,392 B`，
L-prism `1.787205 s / 16,711,680 B`；MacBookAir10,1、Apple M1、8 核、8 GiB、
macOS 26.5.2 arm64、GCC 15.2.0、Release、单线程。两例曾并行启动，因此不能据此比较
两种几何的算法速度。

## 8. 独立验证、确定性与不回归

机器汇总：`artifacts/stage62/acceptance.json`，状态 `pass`。

- cube 原生质量 JSON SHA-256：
  `479458bba0cd6d14b638c4ac04e6558ef190c4b875b560703ead2788073a3b3c`；
- L-prism 原生质量 JSON SHA-256：
  `99cd4390ea6b481c6bd8d5256d3e2e9445e396047a7bf0dc2ae32f7cc36a5f61`；
- 两例重复报告分别与第一次逐字节相同；
- 两例五个 `polyMesh` 核心文件 SHA-256 与 Stage 6.1 逐项相同；
- uniform R8 cube result hash 仍为 `e6fd5337f55140ea`，五文件 SHA-256 仍为：

| file | SHA-256 |
|---|---|
| points | `7d0ddf0e2fe57bf4fb7f589b3b9e2f4e6dfe4955b659f85e94a0eaaf0cab6559` |
| faces | `c333d66020fe7fe4b2c5de2a19ee4918c380da3884740efdf1dcbd7b6cb17791` |
| owner | `6ef35a395352dd62a59eff1fd2d42fc2b90f17fe8337cfb068ea448290fd996d` |
| neighbour | `00d889d1964fc3cc7f15645f4012e4bfd20fbb1c401ed1572ca8a173d17147d9` |
| boundary | `f5817a95bdb15dde0b4d543d4c6119dda55de3a8fc3d1b43b8ee7de214fd0970` |

这证明抽出共享内存 solver mesh 没有改变 Stage 6.1 adaptive 或既有 uniform 核心输出。

## 9. 复现命令

```sh
build/release/cartmesh_cutcell_cli \
  --stl tests/data/closed_unit_cube_ascii.stl \
  --adaptive --base-level 2 --max-level 4 --surface-level 4 \
  --distance 0.1:3 --padding-fraction 0.1 --no-vtk \
  --openfoam-case artifacts/stage62/adaptive_cube_case \
  --quality-output artifacts/stage62/adaptive_cube_quality.json \
  --geometry-output artifacts/stage62/adaptive_cube_geometry.json \
  --report artifacts/stage62/adaptive_cube_report.json

.venv/bin/python tools/stage61_ascii_polymesh_verify.py \
  --case artifacts/stage62/adaptive_cube_case \
  --output artifacts/stage62/adaptive_cube_reader.json

.venv/bin/python tools/openfoam_stage3_verify.py --milestone stage4 \
  --case artifacts/stage62/adaptive_cube_case \
  --project-report artifacts/stage62/adaptive_cube_report.json \
  --output artifacts/stage62/adaptive_cube_checkmesh.json \
  --log-output artifacts/stage62/adaptive_cube_checkmesh.log
```

L-prism 使用相同命令结构和 `tests/data/nonconvex_l_prism_ascii.stl`。最终相关性与确定性汇总
由 `tools/verify_stage62_quality.py` 生成。

## 10. Stage 6.2 关闭边界

Stage 6.2 的通过含义是：reference solver mesh 在写出前已有确定、可定位、可回归的原生
质量诊断，并且固定通过案例与独立 reader/OpenFOAM 证据对应。

它不意味着：

- 历史千万级 compact/binary writer 已接入该评估器；该路径统一属于 6.5；
- 旧 Bunny 质量阻断已经修复；
- 小 Cut-cell、坏 face pyramid、skewness 已自动稳定化；
- 复杂 STL、R24/R48/R96 或百万/千万级 adaptive 已通过；
- Stage 6 或 Stage 7 可以关闭。

下一步只能是 Stage 6.3 小 Cut-cell 与坏形状稳定化，开始前仍需用户明确确认。
