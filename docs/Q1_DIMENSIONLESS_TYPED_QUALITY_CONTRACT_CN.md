# Q1：无量纲、分类型 solver-quality contract

## 范围与事实边界

Q1 新增诊断型 `QualityContract2D`，只负责发现并定位坏网格。实现基准 commit 为
`0ed88b3a9e5d6ef85f6b67be49d00c1bdeda01d8`。本阶段没有修改网格生成、局部
termination、solver repair 或 OpenFOAM writer 算法，也没有降低旧 hard safety
阈值。五案例 Q1 前后的 construction/solver `.cm2d` SHA-256 共 10 项逐项一致。

旧 `SolverQualityReport2D` 仍作为 `legacy_hard_safety` 原样嵌入新报告；其中绝对
`min_face_length_absolute` 继续显示，但它不再是微短边的唯一判断依据。

## `QualityContract2D` 定义

普通单元按 `Cartesian`、`RemainderCut`、`Transition`、`Termination` 四类独立
计数和判级。四类 Q1 初始候选阈值相同，但在 API 中是四份独立 policy，可在有证据
后分别演进。

| metric | preferred | hard | 方向 |
|---|---:|---:|---|
| non-orthogonality | 55° | 65° | 最大值 |
| skewness | 2 | 3 | 最大值 |
| face weight | 0.15 | 0.10 | 最小值 |
| volume ratio | 0.10 | 0.05 | 最小值 |
| minimum interior angle | 20° | 10° | 最小值 |
| hydraulic aspect | 20 | 50 | 最大值 |
| face length / local background h | 0.03 | 0.01 | 最小值 |
| face length / sqrt(owner area) | 0.03 | 0.01 | 最小值 |
| face length / sqrt(neighbour area) | 0.03 | 0.01 | 最小值 |

任一 hard 违反即为 `FAIL`；无 hard 违反但存在 preferred 违反为 `WARN`；二者均无
才为 `PASS`。每个指标输出 `p50/p95/p99/worst`，worst entity 带 cell/edge ID、
坐标、cell/source 类型、source ID、owner/neighbour 与 `local_h`。

`local_background_h` 的来源不是 solver cell 面积反推：Cartesian/remainder 使用其
成员 quadtree background cell 的最小边长；显式 transition 使用 transition ring
thickness；termination 的 quadtree remainder 使用原 background h，ring polygon
使用 ring thickness。共享面的 `face/local_h` 采用相邻普通单元中较大的 h，从而不会
用更细一侧掩盖短边。

## BoundaryLayer 独立统计

BoundaryLayer 不使用普通 aspect 阈值。报告单独计算：

- wall-normal orthogonality error（0° 为理想）；
- growth ratio；
- tangential/normal spacing；
- adjacent-column thickness variation；
- scaled Jacobian；
- first-layer continuity。

用户没有为这六项给出 hard/preferred 数值，因此 Q1 将 BoundaryLayer 标为
`OBSERVED`，不伪造 `PASS`。六项都有分布和 worst entity；后续只有在单独批准
BoundaryLayer contract 后才参与总判级。

## 五案例 baseline

以下结果由 Q1 实现 commit 生成；完整紧凑报告位于 `artifacts/q1/`。`primary hard`
按“超限倍数”选择，只是便于表格展示，不隐藏同案例的其他 hard failure。

| case | overall | 分类型结论 | primary hard（measured / hard） | worst entity |
|---|---|---|---|---|
| circle | FAIL | Cartesian WARN; RemainderCut FAIL; Transition PASS; BoundaryLayer OBSERVED | volume ratio `0.0270379 / 0.05` | cell 442, edge 1005, owner/neighbour 442/581, `(0.367399, 1.44463)` |
| superellipse | FAIL | Cartesian WARN; RemainderCut FAIL; Transition FAIL; BoundaryLayer OBSERVED | face/sqrt(neighbour area) `7.67193e-8 / 0.01` | transition cell 623, edge 185, owner/neighbour 373/623, `(-2.32610, 4.88468e-9)` |
| concave_l | FAIL | Cartesian FAIL; Termination FAIL; BoundaryLayer OBSERVED | volume ratio `0.0106997 / 0.05` | termination cell 2469, edge 3780, owner/neighbour 2469/4616, `(0.96875, 3.02462)` |
| narrow_gap | FAIL | Cartesian FAIL; Termination FAIL; BoundaryLayer OBSERVED | volume ratio `0.0118533 / 0.05` | termination cell 884, edge 1626, owner/neighbour 884/2045, `(0.159844, -0.106459)` |
| sharp_trailing_edge | FAIL | Cartesian FAIL; Termination FAIL; BoundaryLayer OBSERVED | face/local h `0.000420238 / 0.01` | termination cell 1075, edge 5230, owner/neighbour 1075/1078, `(1.79688, -0.339839)` |

当前 hard failure 计数如下；这些是完整 runtime report 中所有 hard issues 的紧凑
汇总，不是只列一个最坏项：

- circle：face weight 24，volume ratio 56；
- superellipse：face/local h 11，face/sqrt(owner area) 3，
  face/sqrt(neighbour area) 7，face weight 16，volume ratio 12；
- concave_l：face weight 64，hydraulic aspect 2，minimum angle 14，
  non-orthogonality 13，skewness 6，volume ratio 107；
- narrow_gap：face/local h 4，face weight 99，hydraulic aspect 2，minimum angle 20，
  non-orthogonality 14，skewness 4，volume ratio 299；
- sharp_trailing_edge：face/local h 14，face/sqrt(owner area) 10，
  face/sqrt(neighbour area) 14，face weight 145，hydraulic aspect 2，
  minimum angle 23，non-orthogonality 11，skewness 2，volume ratio 156。

## superellipse 微短边验收

旧 hard safety 仍报告 absolute minimum face `9.79752897103936e-9`。同一 internal
edge 185 的 Q1 无量纲结果为：

- `face/local_h = 8.2942044199275e-8 < 0.01`，hard FAIL；
- `face/sqrt(owner_area) = 3.0621191849932704e-7 < 0.01`，hard FAIL；
- `face/sqrt(neighbour_area) = 7.671931437413685e-8 < 0.01`，hard FAIL。

该边坐标 `(-2.326100423594548, 4.884676955554503e-9)`，owner 373 为
RemainderCut/source 363，neighbour 623 为 Transition/source 604。新 contract 因而
不是依赖绝对 `minFaceLength` 才发现它。

## 尺度不变性、确定性与 CI

单元测试把同一 typed topology 和 `local_h` 同比放大 `1e6`，要求状态、hard 短边
发现和所有无量纲 worst 基本不变；旧 absolute min face 则按比例变化且继续可见。
此外真实 circle 与 superellipse 输入、first-layer thickness 和 domain padding 同比
放大 `1000` 后，总状态及逐类型状态完全相同，全部无量纲 p50/p95/p99/worst 的
最大绝对差分别为 `7.45e-13` 与 `2.22e-12`。

CI 保留完整 `.hybrid.quality-contract.json`、尺度复跑和 OpenFOAM 输出作为运行证据；
Git 只保留五份约 16 KiB 的紧凑 baseline 与 manifest。baseline 生成命令：

本地最终复跑还对五份同拓扑 OpenFOAM case 实际执行了
`opencfd/openfoam-run:2606 checkMesh -writeAllFields`，五例均得到 `Mesh OK`；这项
OpenFOAM 结果与 Q1 的五个 `FAIL` 不矛盾：Q1 候选 hard contract（例如 65°
non-orthogonality）刻意比当前 OpenFOAM 默认可接受范围更严格。

```bash
python3 tools/verification/generate_q1_baselines.py \
  --repo . --build-dir build --evidence-dir build/h4_solver_ready \
  --output-dir artifacts/q1 --source-commit <implementation-commit> --collect-only
```

尺度复跑命令：

```bash
python3 tools/verification/verify_q1_scale_invariance.py \
  --repo . --build-dir build --reference-dir build/h4_solver_ready \
  --output-dir build/q1_scale --scale 1000
```
