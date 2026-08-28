# Q0：质量事实基线与 provenance

## 范围与代码基准

Q0 只重构质量命名、报告、provenance、CI 证据和 Git artifact 保留策略；没有
修改 mesh generation、H4 transition、solver repair 或质量门限。生成器代码基准：
`fc87f2e3346b37c3f164c1ef27f2a1879bfd13cc`。

## 三类质量

- `construction_quality`：在统一二维构造 topology 上计算。原
  `aspect/skewness` 改名为有公式含义的 `cell_edge_length_ratio` 与
  `centroid_vertex_mean_offset_normalized`，不得解释成 solver 或 OpenFOAM 指标。
- `solver_quality`：只在 solver-repaired topology 上计算。包括 hydraulic aspect、
  concavity、interior angle、face length、non-orthogonality、internal/boundary
  skewness、face weight、volume ratio 与 compactness。
- `openfoam_quality`：只接受外部 `opencfd/openfoam-run:2606` 的
  `checkMesh -writeAllFields` 结果。分位数和最坏 cell 来自 OpenFOAM 写出的
  `volScalarField`，不是内部 solver 指标的别名。

每个 metric 都含 `count/p50/p95/p99/worst/worst_direction/worst_entity`。
`worst_entity` 含 cell/edge ID、二维坐标、cell/source 类型、owner/neighbour、
`local_h`、source ID/key。cell 指标没有唯一 face 时 `edge_id`、`neighbour` 为
`null`；`local_h = sqrt(two-dimensional cell area)`。分位数使用确定性的线性插值。

## 五案例基线

表中三元组均为 `p95 / p99 / worst`。construction 和 solver 的 aspect 名称、
公式和数值刻意分开；OpenFOAM 列来自外部字段。

| case | construction cells | construction edge ratio | solver hydraulic aspect | OpenFOAM aspect worst | OpenFOAM non-orth worst | OpenFOAM skew | checkMesh |
|---|---:|---:|---:|---:|---:|---:|---|
| circle | 700 | 13.2392 / 21.7662 / 21.7662 | 11.0135 / 11.05 / 11.0635 | 5.06241 | 55.3968 | 1.57164 / 1.68659 / 1.68659 | PASS |
| superellipse | 772 | 45.6601 / 121.827 / 1.6568e+07 | 17.0131 / 48.0287 / 62.123 | 31.8866 | 63.4869 | 1.78513 / 2.23217 / 2.24448 | PASS |
| concave_l | 5486 | 6 / 15.0385 / 111.111 | 5.16132 / 10.1902 / 50.7521 | 15.0385 | 69.8254 | 0.839056 / 2.11488 / 3.91981 | PASS |
| narrow_gap | 3244 | 8.59109 / 25.5 / 250 | 6.46605 / 12.947 / 63.3863 | 10.3239 | 68.8133 | 1.27258 / 2.42585 / 3.44239 | PASS |
| sharp_trailing_edge | 3412 | 17.063 / 74.1443 / 2379.6 | 6.40595 / 13.0429 / 67.2738 | 13.912 | 69.3953 | 1.1668 / 1.87599 / 3.34005 | PASS |

这些数值是事实基线，不是“质量问题已解决”的声明。例如 superellipse 的构造
edge ratio 最坏值仍为 `1.6568e+07`；它与 solver hydraulic aspect `62.123`、
OpenFOAM aspect `31.8866` 是三种不同定义，不得横向冒充同一指标。

## 报告与 provenance

Git 只保留：

- `artifacts/q0/<case>.quality-baseline.json`：五份完整基线；
- `artifacts/q0/provenance-manifest.json`：输入 hash、报告 hash、construction/
  solver topology hash、OpenFOAM polyMesh hash、OpenFOAM quality-field hash、
  normalized checkMesh log hash、生成 commit、镜像和 OpenFOAM 版本。

可复现命令：

```bash
cmake -S . -B build/q0 -DCMAKE_BUILD_TYPE=Release -DCARTMESH2D_BUILD_TESTS=ON
cmake --build build/q0 -j2
python3 tools/verification/generate_q0_baselines.py \
  --repo . --build-dir build/q0 --evidence-dir build/q0_evidence \
  --output-dir artifacts/q0 --source-commit "$(git rev-parse HEAD)"
```

CI 对同一批 evidence 连续构建两份报告并执行逐字节 `diff`，完整网格、OpenFOAM
case、外部字段和 checkMesh log 通过 GitHub Actions artifact 上传，不进入 Git。

## 过期 artifact 替换

删除了 `artifacts/h4_2/` 的 16 个历史生成文件，以及 `artifacts/h4_3/` 的 64 个
历史网格、OpenFOAM case、图片、独立读取结果和手工 `checkmesh_status.txt`。这些
文件同时混有旧代码输出与当前 CI 结果，且 H4-2 summary 曾与 solver-quality 文件
互相矛盾。它们由上述 5 份 baseline、1 份 manifest 和 CI 运行证据替代。
