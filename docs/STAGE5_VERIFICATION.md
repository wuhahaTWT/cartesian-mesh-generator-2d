# 阶段 5 验证记录：增量式局部重构

日期：2026-08-09（Asia/Shanghai）

终态：**PASS**

机器可读终态是 `artifacts/stage5_acceptance.json`：

```json
{
  "status": "pass",
  "stage5Complete": true,
  "incrementalResultsEquivalentToFullRebuild": true,
  "stableCellIdsVerified": true,
  "exactFluidOverlapMappingVerified": true,
  "externalIndependentReaderAccepted": true,
  "acceptanceBlockers": []
}
```

## 1. 范围与真实性边界

本阶段只实现局部几何变化下的确定性增量 Cartesian Cut-cell 重构。它不包含 CFD 求解器、
GUI、云服务、CGNS、近壁棱柱层或阶段 6 的千万级完整 Cut-cell 产物。

单次 CLI 报告只声明 `internal_pass_external_pending`，并且保留外部读取、三案例矩阵和独立性能运行三个
blocker。只有终态汇总器才能写出 `stage5Complete=true`。

## 2. 实现契约

- 旧、新表面用与三角形顺序和绕序无关的规范键做精确差分，受影响包围同时覆盖删除和新增三角片。
- 增量八叉树从旧树出发，先尝试安全粗化，再按新规则细化并执行确定性 2:1 平衡闭包。
- 未受影响叶以 64 位 Morton `OctreeNodeCode` 为稳定缓存键复用 Cut-cell 几何；新建、删除和受影响叶重新切分。
- 对全部新叶重建全局邻接、region 和输出索引，然后与独立全量重构比较叶码和完整网格 fingerprint。
- 映射使用 Morton 区间扫描确定保留/重构/细分/粗化关系，并对流体凸多面体做精确几何交，输出旧、新控制体权重。

## 3. 固定三案例矩阵

| 案例 | 旧→新叶数 | 几何复用率 | 本次小案例加速比 | 增量/全量 hash | 外部读取 |
|---|---:|---:|---:|---|---|
| 孔径变化 | 1688 → 1856 | 40.52% | 1.44× | `0574a9ec3d50f24e` | meshio 5.3.5 PASS |
| 孔位变化 | 1688 → 1744 | 70.64% | 1.82× | `eb6ef3fa53aa6b09` | meshio 5.3.5 PASS |
| 局部轮廓变化 | 1688 → 1744 | 79.59% | 2.39× | `1467c52a7e10f4f2` | meshio 5.3.5 PASS |

三类案例全部满足：

- 增量叶码与全量新树逐项一致；
- 增量和全量 Cut-cell 几何/拓扑 fingerprint 一致；
- 负体积、非闭合、共享面不匹配和分类冲突均为 0；
- 所有映射条目都标记 `exactFluidOverlap=true`，权重在 `[0,1]` 内；
- meshio 独立读取旧网格、增量新网格和全量新网格，外部比较稳定节点码、流体体积分数和 Cut-cell 掩码。

例如局部轮廓案例生成 1744 个新叶，复用 1388，重构 356；映射的旧流体体积为
`1.784`，共享精确重叠体积为 `1.784`，新增流体体积为 `0.01`，删除流体体积为 `0`。

## 4. 重复性能测量

固定输入是局部轮廓案例，最大层级 5，每次独立进程都先构建旧基线，再分别计时增量和全量新网格路径。

| 项目 | 实测值 |
|---|---:|
| 重复运行 | 3 次 |
| 旧→新叶数 | 9304 → 9542 |
| 几何复用/重构 | 8676 / 866（90.92% 复用） |
| 增量时间 min / median / max | 0.469721 / 0.506252 / 0.589834 s |
| 全量时间 min / median / max | 1.314247 / 1.426893 / 1.843689 s |
| 加速比 min / median / max | 2.7979× / 2.8185× / 3.1258× |
| 最大观测峰值 RSS | 67,158,016 bytes |
| RSS 口径 | 旧基线 + 增量 + 全量验证 + 导出的同一进程保守峰值 |
| 线程 | 1 |
| 硬件 | MacBookAir10,1，8 GiB，8 个逻辑 CPU |
| 构建 | Release，GCC 15.2.0 |
| 结果 hash | `d13d490a35e8178f` |

导出时间不在两条计算路径的对比时间内，但 RSS 口径故意保留导出后的进程峰值。不从三个小案例的
单次时间推导广泛性能结论；性能门禁依据是这三次非微型重复测量。

## 5. 测试与独立外部验证

- Release：CTest 20/20 PASS，最后一次墙钟 6.87 s；
- Debug：CTest 20/20 PASS，墙钟 103.48 s；
- Apple Clang AddressSanitizer + UndefinedBehaviorSanitizer：CTest 20/20 PASS，墙钟 303.42 s；
- 三组真实 VTU 由 meshio 5.3.5/NumPy 独立读取和比较；
- 共享 Cut-cell 几何修复后，重新生成阶段 4 双区域薄壳 OpenFOAM `polyMesh`；OpenFOAM 2606
  `checkMesh -constant -allTopology` 在禁网容器中报告 `Mesh OK`。1728 个单元，最小体积 `0.000125`，
  最大非正交度 `0`，最大偏斜度 `0.5`。

LeakSanitizer 未在本轮 macOS 工具链下运行，不宣称已有泄漏检测证据。OpenFOAM 回归证明阶段 3/4
完整体网格路径没有因本次共享修复回归；阶段 5 自身的增量等价性由三案例全量对照和 meshio 检查证明。

## 6. 发现并保留的最小失败案例

局部轮廓案例首次运行暴露了 12 个非闭合/共享面不匹配。报告当时正确标记为失败，没有隐藏。
原因是 STL 壁面与 Cartesian 面共面时，开口面包含了这片壁，而嵌入边界又计算了同一片面，导致闭合重复计数。

修复后，Cartesian 开口由显式流体凸片重建，并减去共面嵌入壁面。最小回归输入永久保留为
`tests/data/stage5_coplanar_tunnel_corner_ascii.stl`，对应测试是 `test_coplanar_tunnel_corner_regression`。

## 7. 复现命令

```sh
cmake --preset release
cmake --build --preset release --parallel 4
ctest --preset release --output-on-failure

python3 tools/generate_stage5_fixtures.py

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

.venv/bin/python tools/render_stage5_incremental.py \
  --old-mesh artifacts/stage5_local_contour_old.vtu \
  --new-mesh artifacts/stage5_local_contour_incremental.vtu \
  --output artifacts/stage5_local_contour_comparison.png
```

孔径和孔位案例用同一组参数和各自的 `*_old.stl` / `*_new.stl` 复现。终态汇总命令为：

```sh
.venv/bin/python tools/aggregate_stage5_benchmark.py \
  artifacts/stage5_perf_m5_run1.json \
  artifacts/stage5_perf_m5_run2.json \
  artifacts/stage5_perf_m5_run3.json \
  --output benchmarks/baselines/stage5_incremental_m5_2026-08-09.json

.venv/bin/python tools/verify_stage5_acceptance.py \
  --case artifacts/stage5_hole_diameter.json \
  --case artifacts/stage5_hole_position.json \
  --case artifacts/stage5_local_contour.json \
  --external artifacts/stage5_hole_diameter_meshio.json \
  --external artifacts/stage5_hole_position_meshio.json \
  --external artifacts/stage5_local_contour_meshio.json \
  --performance benchmarks/baselines/stage5_incremental_m5_2026-08-09.json \
  --output artifacts/stage5_acceptance.json
```

## 8. 用户可视检查点

修改前/修改后/实际重构区三联图：

![阶段 5 局部轮廓增量重构对比](../artifacts/stage5_local_contour_comparison.png)

用户检查建议见 [`USER_CHECKPOINT_STAGE5_CN.md`](USER_CHECKPOINT_STAGE5_CN.md)。阶段 6 和 7 未开始。
