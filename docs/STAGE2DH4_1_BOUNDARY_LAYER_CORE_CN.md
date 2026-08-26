# Stage 2D-H4-1：原生二维边界层核心

日期：2026-08-27

分支：`agent/native-2d-robustness`

基线：`96616c68718eb2fd825544831a1cac6ec563e03c`

## 1. 阶段边界

H4-1 新增独立的原生二维 body-fitted quad boundary-layer strip 候选模块：

```text
BoundaryLoop / ordered open polyline
  -> deterministic WallChain2D
  -> fluid-side marching directions
  -> safe outer envelope
  -> validated one-quad wrapper
  -> geometric hair-edge subdivision
  -> validated N-layer strip
```

本阶段不把 outer envelope 与 Cartesian/Cut-cell remainder 共形拼接，不修改
`CutCell2D`、global/solver topology、OpenFOAM writer，也不生成正式 hybrid mesh。
这些仍属于 H4-2，当前未开始。

## 2. 新数据语义与事务边界

`WallChain2D` 显式保存：

- ordered vertices 与 ordered segments；
- closed/open；
- clockwise/counter-clockwise/open orientation；
- `FluidSide2D`；
- chain ID 与 patch identity。

闭合外流 wall 统一规范为 counter-clockwise：固体在有向边左侧，流体在右侧。
闭合输入同时进行字典序最小点 canonical rotation，因此循环起点或原始方向变化不会改变
最终 wall、layer vertex 和 cell ID。open chain 不存在多边形 inside/outside，调用者必须显式给
fluid side。

`buildBoundaryLayerStrips2D()` 返回独立 `BoundaryLayerBuildResult2D`。先建立并验收一层
wrapper，再建立完整多层 candidate；任一步失败只返回结构化 `BoundaryLayerFailure2D`，
输入 wall 和原有网格数据不发生修改。

## 3. 参数模型

`LayerParameters2D` 支持：

- `nLayers >= 1`；
- `FirstLayerThickness` 或 `TotalThickness`；
- `growthRatio > 0`。

参数必须 finite、positive，不做 silent correction。几何级数使用 `log/expm1` 计算；
`growthRatio` 接近 1 时走独立的 `H = n * h1` 路径，避免相消。解析后保存每一层的
累计法向距离，最后一层强制等于解析总厚度，保证重复运行一致。

## 4. 推进方向与角点分类

每条 wall segment 先单位化并按 fluid side 得到有向法向。内部顶点使用相邻两个
fluid-side half-plane：只有候选 marching direction 对两个法向的点积都为正时才接受。

顶点分为：

- `Smooth`：固体侧转角绝对值不超过 5°；
- `MildConvex`：5° 到 135°，包括矩形 90°角；
- `Concave`：小于 -5°，H4-1 fail-closed；
- `Sharp`：大于 135°，H4-1 fail-closed；
- `Degenerate`：相邻半平面冲突或段长落入 tolerance，fail-closed。

方向不是对任意几何无条件执行 `normalize(n1+n2)`：只有在角点分类通过、half-plane
约束一致、miter cosine 高于尺度化门槛后才构造角平分 marching direction。

## 5. 安全厚度与 outer envelope

所有阈值集中在 `BoundaryLayerPolicy2D`：

- 长度/面积 epsilon 来自既有 `TolerancePolicy` 与 wall bounds 尺度；
- `smoothTurnRadians = 5°`；
- `maxConvexTurnRadians = 135°`；
- `cornerLengthFraction = 0.45`；
- `collisionClearanceFraction = 0.45`。

局部 miter 的切向占用不得超过相邻短边的 45%。每个 hair ray 会和所有非 incident
wall segments 求交，得到窄缝/对向 wall 的距离上限；不同 chain 还使用全局
segment-to-segment gap，保证两侧各自占用小于半缝宽。请求总厚度超过 safe limit 时直接
返回 `thickness_exceeds_safe_limit`，记录 requested/safe thickness，不自动缩小。

候选 ring 会检查：非邻接自交、与原 wall 相交、hair 穿越非 incident wall、不同 chain
碰撞。outer envelope 不是逐点“看起来合理”即通过。

## 6. Quad wrapper 与多层拓扑

fluid 在右时，每个 cell 的确定性顶点序为：

```text
inner_i -> outer_i -> outer_(i+1) -> inner_(i+1)
```

fluid 在左时使用镜像顺序。算法不会在发现负面积后临时 reverse cell；语义产生的顺序若
不能得到正面积，candidate 失败。

一层 wrapper 先独立检查：

- quad signed area 大于尺度化面积 tolerance；
- quad 不自交；
- ring 无自交且不穿 wall；
- 非邻接 cell 无边相交或内部重叠；
- 几何重复点不能使用不同 ID；
- wall/outer boundary edge incidence 为 1；
- 中间 ring 和闭合 hair edge incidence 为 2；
- open chain 两端 hair edge incidence 为 1。

通过后才沿所有 hair 使用相同 `nLayers` 和几何级数细分。每一段 hair spacing 必须严格
为正，距离 wall 单调增加；当前没有 local drop、termination 或临时补三角形。

## 7. 结构化失败

`BoundaryLayerFailure2D` 记录：status/reason/message、chain、可选 vertex/edge/cell ID、
requested/safe thickness、nLayers 与 growthRatio。自动化覆盖：

- severe concave -> `concave_corner`；
- trailing-edge-like sharp triangle -> `sharp_corner`；
- two-chain narrow gap -> `thickness_exceeds_safe_limit`；
- duplicate/zero-length segment -> wall-chain construction failure；
- NaN、非正 thickness/growth ratio、零层数 -> `invalid_parameters`。

失败 candidate 不删除 cell、不降低层数、不缩小厚度，也不把纯 Cut-cell fallback 报告为
boundary-layer success。

## 8. 调试出口与独立检查

新增 `cartmesh2d_boundary_layer_cli`：

```sh
cartmesh2d_boundary_layer_cli \
  <boundary.xy> <output-prefix> \
  <n-layers> <first|total> <thickness> <growth-ratio> \
  [exterior|interior]
```

成功输出 H4-1 专用 legacy VTK 与 JSON；失败仍写结构化 JSON 并返回非零。该 CLI 和
writer 只输出 isolated strip，不修改正式 OpenFOAM hybrid writer。

不链接项目库的 `check_boundary_layer2d.py` 重新解析 VTK，独立检查：

- VTK quad 类型与 ID；
- 正 signed area；
- quad/ring 自交；
- 非共享 cell edge 相交与 cell interior overlap；
- edge incidence；
- hair distance 严格单调；
- VTK 和 JSON 的 cell/vertex count 与 min/max area 一致。

失败 JSON 由独立 `check_boundary_layer_failure.py` 验证必需字段和 failure reason。

## 9. 真实可视化产物

固定产物目录：`cartmesh2d/artifacts/h4_1/`。

| 案例 | layers | vertices | quads | min/max cell area | min/max hair spacing | 结果 |
|---|---:|---:|---:|---:|---:|---|
| 32 段圆 | 4 | 160 | 128 | 0.003960082 / 0.007388185 | 0.020096771 / 0.034727221 | PASS |
| 矩形 90°凸角 | 3 | 16 | 12 | 0.0416 / 0.14015625 | 0.056568542 / 0.088388348 | PASS |
| 凹 L 形 | 3 requested | 0 | 0 | N/A | N/A | FAIL-CLOSED: `concave_corner`, vertex 3 |

圆和矩形均保存 `.layer.vtk`、`.layer.json`、`.independent.json` 和 `.layer.svg`。
SVG 明确显示 wall、每层 edge、hair edges、outer envelope 和 quad cells。凹角保存
`concave.layer.json`。

## 10. 自动化与回归

针对 H4-1：

```sh
cmake -S . -B build/release-2d -DCMAKE_BUILD_TYPE=Release
cmake --build build/release-2d \
  --target cartmesh2d_boundary_layer_cli cartmesh2d_boundary_layer_tests -j4
ctest --test-dir build/release-2d \
  -R '^cartmesh2d_stageH4_1_' --output-on-failure
```

结果：`10/10 PASS`。其中包含 core test、圆/矩形 E2E、两个独立 reader、凹角失败
reader 和圆 VTK/JSON 字节级确定性比较。

完整二维回归：

```sh
cmake --build build/release-2d -j4
ctest --test-dir build/release-2d -R '^cartmesh2d_' --output-on-failure
```

结果：`63/63 PASS`，CTest 墙钟 `12.02 s`。本阶段没有修改 H1 refinement、Quadtree、
Cut-cell、global topology、agglomeration、solver topology、quality 或 OpenFOAM writer；
因此 H3 scalability 算法路径未改变。本轮未重复运行 500k H3 性能矩阵，也不新增性能结论。

真实 OpenFOAM `checkMesh`：**NOT RUN / NOT APPLICABLE**。H4-1 还不是 hybrid solver mesh，
不得把 isolated strip VTK 的独立检查冒充 OpenFOAM 验收。

## 11. 阶段判定

H4-1 按授权边界 **PASS**：已生成具有确定 wall/fluid 语义、正面积 quad、闭合共享拓扑、
安全 outer envelope、固定层数 geometric spacing、独立检查与结构化失败的原生二维
boundary-layer strip。

```text
H4-2 started: NO
production hybrid cut-cell integration: NO
```
