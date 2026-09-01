# R2/W0：加密阶梯测量与回归门

日期：2026-09-01
范围：只增加测量、归因与回归门；**没有**修改任何网格生成算法，也没有改动任何质量阈值。

## 为什么需要这一步

CI 里每个验收案例只跑一个固定 level，`docs/STAGE2DV1C_VERIFICATION.md` 的 PDE
收敛序列也只到 `(minimum, boundary) = (6, 8)`。因此「把同一几何继续加密会发生
什么」在本轮之前没有任何测量数据，也没有回归门。

R1F 第 7 节同时留了一个硬前置条件：现有 solver 子阶段计时**加总解释不了端到端
wall time**，在补齐顶层 phase attribution 之前不得做性能优化。W0 一并关闭这一项。

## 1. 新增的测量工具

`tools/verification/refinement_ladder.py` 对一个几何跑 level 阶梯，逐级记录
确定性网格事实与 wall-clock 归因，输出单份 manifest：

```bash
python3 tools/verification/refinement_ladder.py \
  --repo . --build-dir build --evidence-dir build/r2_ladder \
  --output-dir artifacts/r2 --manifest-name w0-baseline-manifest.json \
  --dyld-library-path /Applications/mesasdk/lib --stop-on-failure \
  --ladder circle:cutcell:6,7,8,9,10 \
  --ladder circle:hybrid:6,7,8,9 \
  --ladder narrow_gap:hybrid:8
```

`--ladder CASE:MODE:LEVELS`，`MODE` 为 `hybrid`（`cartmesh2d_hybrid_cli`）或
`cutcell`（`cartmesh2d_cli`）。`cutcell` 阶梯自动使用
`(minimum, boundary) = (level - 2, level)`，因为只提高边界 level 不会细化远场，
不构成受控序列——这一点沿用 V1c 的结论。

manifest 把不可复现的 wall time 与确定性数据分开存放：每级的 `attribution` 段带
`"measurement_class": "wall_time", "reproducible": false`，只有 `mesh` 段与
`deterministic_counters` 参与门判定与跨运行对照。

## 2. 门的形式：单调性，不是绝对值

门在 `refinement_ladder.py:gate()`，单元测试在
`tests/refinement_ladder_gate_test.py`（不生成任何网格，直接对合成 rung 判定）。
判定项：

1. 阶梯内任何一级 `exit_code != 0`、超时，或 `hybrid_status != success` → 违规；
2. `topology_valid` / `mesh_quality_valid` / `solver_quality_valid` 为 false → 违规；
3. 跌破 `SolverQualityPolicy2D` 既有硬限 → 违规：`min_face_weight < 0.05`、
   `min_volume_ratio < 0.01`、`max_non_orthogonality_deg > 70`。这三个值是
   **镜像**自 `include/cartmesh2d/quality/SolverQuality2D.hpp`，C++ 侧仍是唯一权威；
   它们同时与 OpenFOAM `etc/caseDicts/meshQualityDict` 的
   `minFaceWeight 0.05` / `minVolRatio 0.01` / `maxNonOrtho 65` 一致；
4. Q1 contract 的 hard issue 总数**随 level 增长** → 违规。

本轮没有降低任何阈值。门只是拒绝接受一条已经退化的阶梯。

## 3. 顶层 phase attribution（关闭 R1F 第 7 节的前置条件）

### 3.1 hybrid 路径

`RobustH4Profile2D`（`include/cartmesh2d/hybrid/HybridMesh2D.hpp`）为
`buildRobustH4Mesh2D` 的五个阶段各加一个计时器，并新增进程级计数器
`conformalHybridBuildCount2D()`，风格与既有的 `globalTopologyBuildCount2D()`
一致：调用方可以用**实测差值**而不是断言来证明一次运行消耗了多少次完整 hybrid
构造。`cartmesh2d_hybrid_cli` 在**成功与失败的每一条出口**都打印同一组
`timing_*` / `h4_*` 键；被拒绝的网格恰恰是最需要知道时间去哪了的情况。

seconds 只写进 `.hybrid.profile.json`（该文件本来就标记
`reproducible: false`），确定性的调用计数才允许进入被逐字节比较的
`.hybrid.json` 之外的对照。

circle level 9（固定层厚 0.02，4 层）实测：

```text
timing_total_seconds              136.517
timing_input_seconds                0.002
timing_build_seconds              136.513
timing_export_seconds               0.001
h4_requested_layer_seconds          0.002
h4_requested_hybrid_seconds        22.881
h4_local_layer_seconds              0.003
h4_local_hybrid_seconds           113.551
h4_pure_cutcell_fallback_seconds    0.074
h4_unattributed_seconds             0.000001
h4_conformal_hybrid_build_calls     2
```

三个结论：归因已经闭合（未归因 1e-6 s）；时间并**不是**花在
growth-ratio 重试上（只发生了 2 次 conformal build，circle 没有触发
`localReductionApplied` 的候选族）；而是花在**单次** conformal 构造内部的修复
回路上，且局部降层候选比原始候选慢 5 倍（113.6 s vs 22.9 s）。
pure Cut-cell fallback 只要 0.074 s，却被 solver-quality 门拒绝。

### 3.2 纯 Cut-cell 路径

`cartmesh2d_cli` 原本只在 PASS 路径打印 `timing_*`。现在三条后期失败出口
（quality 无效、solver topology 失败、solver-quality 门拒绝）也打印，并新增
`timing_solver_topology_seconds` 覆盖被拒绝时仍在运行的那个阶段，避免留下隐式
时间。circle level 10 实测：

```text
timing_refinement_seconds        0.046
timing_balance_seconds           0.049
timing_cut_cell_seconds          0.300
timing_source_topology_seconds   0.304
timing_agglomeration_seconds     0.369
timing_solver_topology_seconds  16.124
timing_total_seconds            17.196
```

命名阶段合计 17.191 s，进程内总计 17.195 s，差 4 ms。**17.2 s 里有 16.1 s 在
solver topology 阶段**，与 narrow-gap 早先的 `partitionSourcePolygons` 采样结论
一致；这条阶段归因现在是实测，而不是假设。

## 4. W0 实测基线

机器可读基线：`artifacts/r2/w0-baseline-manifest.json`（本机 macOS/mesasdk，
Release）。生成命令见第 1 节。`--stop-on-failure` 让每条阶梯停在它第一次失败的
那一级，因此基线记录的就是当前的墙位置。

### 4.1 纯 Cut-cell 阶梯（`(minimum, boundary) = (level - 2, level)`）

| level | solver cells | wall time | 结果 | min faceWeight | min volRatio | max nonOrtho |
|---:|---:|---:|---|---:|---:|---:|
| 6 | 412 | 0.08 s | PASS | 0.0527 | 0.0432 | 68.56° |
| 7 | 1280 | 0.09 s | PASS | 0.1510 | 0.1131 | 61.75° |
| 8 | 3828 | 0.61 s | PASS | 0.0827 | 0.0552 | 59.96° |
| 9 | 13024 | 0.54 s | PASS | 0.0954 | 0.0247 | 42.31° |
| 10 | — | 17.5 s | **FAIL** | — | — | — |

到 level 9（约 2 万 quadtree leaves、13024 个 solver cell）全部硬限都满足，而且
non-orthogonality 随加密**变好**。性能与质量都不是这条路径的瓶颈。

level 10 的失败被门归类为 `solver_quality:excessive_boundary_skewness`。它是单点
几何退化：边界折线在格点 `(1.0104296875, -0.298125)` 附近通过，切出一个腿长
`3.5e-8` 与 `1.2e-7` 的角落三角形，boundary face skewness 达到 `8.5955`（限 4）。
`IntersectionRegistryPolicy2D::snapFractionOfLocalH` 目前是 `64 · DBL_EPSILON`，
level 10 下吸附半径约 `8e-17`，比这个 `6e-6 · h` 的间距小八个数量级，所以永远
不会被焊合。这是 W1 的对象。

必须说清一处尺度关系，避免把结论说过头：本仓库
`SolverQualityPolicy2D::maxBoundarySkewness` 是 `4`，而 OpenFOAM
`meshQualityDict` 的 `maxBoundarySkewness` 默认是 `20`。也就是说
`8.5955` 这个值**会被 OpenFOAM 的 snappy 质量门接受**，是本仓库更严格的门拒绝了它。
但拒绝本身是对的：`3.5e-8` 的面出现在 `h ≈ 5.7e-3` 的单元上，是真实的几何退化，
不论阈值多宽都应该在构造期消掉，而不是靠放宽阈值通过。W1 的做法是消除退化本身，
不动阈值。

### 4.2 hybrid 阶梯（circle，4 层，首层 0.02，增长 1.2）

| level | solver cells | BL cells | wall time | 结果 | Q1 hard | min faceWeight | min volRatio | max nonOrtho |
|---:|---:|---:|---:|---|---:|---:|---:|---:|
| 6 | 728 | 128 | 0.17 s | PASS | 80 | 0.0872 | 0.0270 | 55.40° |
| 7 | 1100 | 128 | 0.19 s | PASS | 132 | 0.0891 | 0.0107 | 68.72° |
| 8 | 2072 | 128 | 0.99 s | PASS | 360 | 0.0643 | 0.0119 | 67.91° |
| 9 | — | — | 132 s | **FAIL** | — | 0.0244 | 0.0023 | 78.61° |

narrow_gap level 8（当前 CI 级别）：3185 cells、46 s、PASS、Q1 hard 430、
min faceWeight 0.0516、min volRatio 0.0119、max nonOrtho 68.81°。

**BL cell 数在每一级都是 128**：wall column 直接取输入折线顶点，`circle.xy` 只有
32 个顶点，32 column × 4 层 = 128，与 level 无关。level 9 的背景 `h ≈ 0.0113`，
而切向尺寸恒为 `2πr/32 ≈ 0.196`——背景网格比壁面网格细 17 倍，过渡环必须独自
吸收整个落差。这是 W2 的对象。

另有一项对照实验（未进基线）：把首层厚度随 level 折半（level 9 用 0.0025）
**没有**救回来，只是把失败从 132 s 拖到 803 s。所以问题在切向分辨率，不在法向厚度。

### 4.3 基线里记录的门违规

```text
circle:cutcell:level 10 exited 1 (solver_quality:excessive_boundary_skewness)
circle:hybrid:level 7 Q1 hard issues grew 80 -> 132
circle:hybrid:level 8 Q1 hard issues grew 132 -> 360
circle:hybrid:level 9 exited 1 (solver_quality_failed)
```

这四条是**已知且已登记**的当前状态，不是本轮引入的回归。

## 5. CI 边界

- `Gate the pure Cut-cell refinement ladder`：`circle:cutcell:6,7,8,9` 带 `--gate`。
  这条阶梯今天就满足全部硬限，所以它是一道真实的回归门。
- `Record the hybrid refinement ladder`：`circle:hybrid:6,7,8`，**不带** `--gate`。
  Q1 hard 计数仍随 level 增长，第 4.2 节已登记；W2 落地后改成 `--gate`。
  level 9 需要约 150 s，留给本地长阶梯。

## 6. 本轮的边界

- 没有修改任何网格生成、termination、repair 或导出算法；
- 没有降低 `SolverQualityPolicy2D` 或 Q1 contract 的任何阈值；
- 完整 CTest 从 75 项增加到 77 项，全部通过（`ctest --test-dir build`，
  macOS 下需先 `export DYLD_LIBRARY_PATH=/Applications/mesasdk/lib`）。

### 6.1 网格产物逐字节不变的实测证据

从改动前的 `051f45103ffb34c7568a4210cfb7de93046af89f` 用 `git archive` 取出干净
源码树、独立 Release 构建，与本轮工作区用**同一条旧命令行**（不带任何新开关）
生成同一批产物并逐一比对 SHA-256：

```text
cartmesh2d_hybrid_cli circle.xy <out> 6 3 6 4 0.02 1.2 1.0 <case> 0.01
cartmesh2d_cli        circle.xy <out> 8 0.25 0.10 exterior <case> 6
```

| 产物 | 结果 |
|---|---|
| `circle.hybrid.cm2d` | IDENTICAL `07cb6e2207f405a1b89caaf9cf0923a4bca4a896a5ffba5632168a9661db8b40` |
| `circle.hybrid.solver.cm2d` | IDENTICAL `74a3fb601352641da4f10a9d4f8acf843c5523ac9cf8573c1e24e2aa26542274` |
| `circle.hybrid.json` | IDENTICAL `53bfdf30ec635593ebfc39194d384ff768ad41fdbc50f4baf3d0efe294067be6` |
| `cut.cm2d` | IDENTICAL `8c5e604524b4f980b6d401e0508a2049d12b6e622a3672b8bebc4525e0752f3b` |
| `circle.hybrid.vtk` / `.solver.vtk` | IDENTICAL |
| `circle.hybrid.quality-contract.json` | IDENTICAL |
| `circle.hybrid.construction.json` / `.intersections.json` | IDENTICAL |
| `cut.vtk` / `.sizing.json` / `.viz.json` / `.construction-quality.json` | IDENTICAL |
| OpenFOAM case 目录（`diff -r`，两个案例） | IDENTICAL |

唯一变化的文件是 `.hybrid.profile.json`，它本来就标记
`"measurement_class": "wall_time", "reproducible": false`，本轮**只新增**了 11 个
`h4_*` 字段，没有修改或删除原有字段：

```text
h4_total_seconds / h4_requested_layer_seconds / h4_requested_hybrid_seconds
h4_local_layer_seconds / h4_local_hybrid_seconds
h4_pure_cutcell_fallback_seconds / h4_unattributed_seconds
h4_requested_hybrid_attempts / h4_local_hybrid_attempts
h4_pure_cutcell_fallback_attempts / h4_conformal_hybrid_build_calls
```

另外两条 CLI 在 stdout 上新增了 `timing_*` / `h4_*` 键。stdout 不参与任何逐字节
比较，`.hybrid.json` 与全部网格文件不变。


