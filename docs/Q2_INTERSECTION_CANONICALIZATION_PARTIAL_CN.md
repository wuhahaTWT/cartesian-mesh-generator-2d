# Q2 — Intersection Canonicalization（部分完成，尚未验收）

## 结论

本提交修复 superellipse 的 `~9.8e-9` internal face，并加入 local-h intersection
registry、来源记录、feature/support 保护、回归和五案例证据。

**Q2 整体仍为 `PARTIAL_NOT_ACCEPTED`。** narrow gap 与 sharp trailing edge
的既有短面仍未满足 Q1 `face/local_h >= 0.01`。不能用 74/74 CTest 或五例
`checkMesh` PASS 代替这一验收项，也没有把这两例的 pure Cut-cell fallback
当作 H4 成功结果。

原始基线为 `556bb90`，从该提交独立导出源码、编译并重新生成五案例。
`artifacts/q0/`、`artifacts/q1/` 和输入 `.xy` 文件保持不变。

## 根因

superellipse 输入并非在所有对称轴上精确为零，例如顶点
`(-1.8, 8.9e-9)`、`(1.41e-8, 0.8)`。这些采样残差经 boundary-layer / transition
外推后，使 envelope vertex 几乎落在 Cartesian grid line 上，而不完全重合。
这不是仅靠重新计算两次求交就能消除的浮点舍入。

原 solver edge 185，owner/neighbour 373/623：

```text
(-2.326100423965795, 9.769353911109006e-9)
(-2.3261004232233007, 0)
length = 9.79752897103936e-9
face/local_h = 8.2942044199275e-8
```

前者为 transition 采样顶点，后者为其相邻 envelope 段与 grid line 的交点。
原 common-partition topology 正确保留了两个不同点，但由此生成的微短面不符合 Q1。
另外还存在 envelope 段近掠 Cartesian corner 的短面，单纯消除轴上残差不够。

## 修改位置与安全边界

- `geometry/IntersectionRegistry2D`：按 support ID 隔离候选，距离阈值使用
  `fraction * min(query_local_h, anchor_local_h)`，不使用固定绝对 epsilon。
  sharp/concave feature 输入顶点不可移动，不能跨不相干 support 合并。
- `cutcell/CutCell2D.cpp`：wall/grid 和作为 cutter boundary 的 transition/grid
  交点共用 registry；原始 wall samples 保持不动。box snapping 改用
  `tolerance.relative * local_h`，而不是固定绝对距离。每个实际移动的交点
  保留 segment、原位置、canonical vertex、displacement、local_h、feature。
- `hybrid/TransitionCanonicalization2D`：只对非 stepped 的可变 transition 前沿
  做源头采样调整。必须存在沿相邻 envelope 段的真实短交段；仅仅垂直距离小
  不构成调整理由。近掠 grid corner 时插入 canonical corner，双方 source polygon
  同时更新，之后才构造 Cut-cell。原始物理 wall 和固定 H4 layer interface 不动。
- 安全位移还受非相邻轮廓/固定 layer interface 的 clearance 限制。改形后所有
  transition polygon 必须仍为正面积、简单多边形，再经过既有共形、面积、分类、
  solver-quality 门。没有删除 cell，没有改 solver 阈值，没有引入新 DCEL。
- CLI 输出 `*.hybrid.intersections.json`。`canonical_vertex.local_id` 是局部
  registry ID，不能冒充 solver ID；独立的 `solver_vertex_id` 给出实际 solver
  对应关系，若该 source point 未保留则为 null。

当前 stepped termination 前沿保留原构造，不宣称已经完成 H4-3 的安全 snap / refinement。
此边界在代码、报告与验证脚本中均明确标注。

## 五案例前后对比

统一使用 Q1 原有参数、原有 local_h 定义和 `0.01` hard limit。

| case | 最短 face/local_h：修复前 | 修复后 | Q1 短面门 | OpenFOAM 2606 |
|---|---:|---:|---|---|
| circle | 0.045942780858514 | 0.045942780858514 | PASS | Mesh OK |
| superellipse | 8.2942044199275e-8 | 0.016222929810346797 | PASS | Mesh OK |
| concave L | 0.014400000000000546 | 0.014400000000000546 | PASS | Mesh OK |
| narrow gap | 0.007843137254903055 | 0.007843137254903055 | **FAIL（未解决）** | Mesh OK |
| sharp trailing edge | 0.0004202383138647292 | 0.0004202383138647292 | **FAIL（未解决）** | Mesh OK |

superellipse 最短绝对 face 为 `0.0019163335838472152`，不再有原 `9.8e-9` 面。
construction cell 数仍为 772，固定 layer cell 数仍为 72；solver 分区由 791 变为 795，
不是删掉坏 cell。`area_error = 3.019806626980426e-14`。

## 质量与独立验证

- 全部 CTest：**74/74 PASS**（原 73 项 + registry 测试）；最终运行耗时 52.48 秒。
- 五例 construction/solver CM2D、quality-contract、intersection provenance
  重复生成逐字节一致。前后版本的文件不一定逐字节相同，不能将重复确定性混同于
  算法版本之间的字节恒等。
- 五例完整 owner/neighbour、正体积、闭合性和边界检查均通过独立
  `check_openfoam2d.py`；circle/superellipse 另通过 H4-2 VTK 独立读取器。
  该 H4-2 reader 要求封闭二价前沿，原 concave L 基线也会被其拒绝，因此不能
  对 H4-3 终止前沿套用这一专用 reader。没有放宽 reader 判据。
- 五例真实执行 `opencfd/openfoam-run:2606 checkMesh -writeAllFields`：**Mesh OK**。
  镜像 ID：`sha256:4229997e74defb81548222d511b8e3b95b98305e5df41b8e88b031813fe47eeb`。
  这是当前既有默认 checkMesh 门，不是新增 `-allGeometry -allTopology` 的声明。
- circle、superellipse 整体放大 1000 倍：无量纲分布最大绝对差分别
  `7.442935157087049e-13`、`2.8350655156827997e-12`。
- 五案例未增加 Q1 hard issue 计数，也没有检测到 ordinary / BoundaryLayer
  最坏指标退化（比较容差 `1e-8 * max(1, abs(old))`）。superellipse 的三个短面
  hard issue 计数由 `11/3/7` 降为 `0/0/0`，volume-ratio hard issue 由 12 降为 10，
  face-weight hard issue 仍为 16。最小角度由 `16.1591°` 改善为 `20.0902°`。
  这不意味着所有局部 cell 的所有指标都单调改善；完整 Q1 总状态五例仍为 FAIL。

完整数值、源代码逐文件 SHA-256、输出与 checkMesh 日志 hash 见
`artifacts/q2/comparison.json`。实际网格预览见 `artifacts/q2/superellipse-solver.png`。
原始 solver mesh / provenance / 日志保留在本机 `build-q2/before/` 与 `build-q2/after/`。

## 未采用的候选与下一步

局部加密已试验但没有作为成功修复保留：

1. 对 H4-3 近重合片段局部 quadtree 加密并重新 balance，narrow gap 的候选出现
   `max_nonorthogonality=73.3453`、`min_face_weight=0.0485844`，旧 safety 门失败。
   仅对贴近不可变 layer interface 的片段定向加密也出现同一失败。
2. 对 sharp trailing edge 的可变前沿局部改形/加密，出现 solver partition
   非流形或 `max_nonorthogonality=77.1089`、`min_face_weight=0.0339435` 等失败。
   没有通过调低阈值、删除 cell 或导出 fallback 来掩盖。

narrow gap 的典型未解决面端点是 `(-0.012,-0.012)` 与 `(-0.011875,-0.012)`；
前者紧邻不可变层终止角点。下一步需要在固定 layer/feature 约束下联合处理
局部 grid 相位/前沿采样与相邻 solver 分区，保留这些失败点作为针对性回归。
当前提交只是一份经过验证的 Q2 部分结果，不能登记为 Q2 完成。

## 复现

```bash
cmake -S . -B build-q2 -DCMAKE_BUILD_TYPE=Release
cmake --build build-q2 -j4
ctest --test-dir build-q2 --output-on-failure -j2

# 从未修改的 Q1 提交单独构建基线，不改变当前 checkout。
mkdir -p build-q2/baseline-source
git archive 556bb90 | tar -x -C build-q2/baseline-source
cmake -S build-q2/baseline-source -B build-q2/baseline-build \
  -DCMAKE_BUILD_TYPE=Release -DCARTMESH2D_BUILD_TESTS=OFF
cmake --build build-q2/baseline-build --target cartmesh2d_hybrid_cli -j4
python3 tools/verification/generate_q1_baselines.py \
  --build-dir build-q2/baseline-build --evidence-dir build-q2/before \
  --output-dir build-q2/q1-before --source-commit 556bb90

python3 tools/verification/verify_q2_intersections.py --repeat --openfoam
python3 tools/verification/verify_q1_scale_invariance.py \
  --build-dir build-q2 --reference-dir build-q2/after \
  --output-dir build-q2/q2-scale-final --scale 1000
MPLCONFIGDIR=/private/tmp/cartmesh2d-q2-mpl \
  python3 tools/visualization/render_q2_microedge.py
```

`verify_q2_intersections.py --require-full-acceptance` 会在当前状态以非零码退出，
因为两例短面门尚未通过。CI 的 Q1 收集增加显式 `--expect-superellipse-short-faces absent`
以验证已修复的 superellipse；未改动历史 Q1 基线或其他指标/阈值。
