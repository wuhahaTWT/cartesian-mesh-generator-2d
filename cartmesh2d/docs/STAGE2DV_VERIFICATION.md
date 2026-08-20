# Stage 2D-V 验证记录

## 状态

**READY FOR VALIDATION — 2D-V 第一版已实现，并已在 GitHub Actions run #14 完成当前分支真实编译、renderer 回归与四类 acceptance SVG 生成。等待用户显式 `验证-v` 后做最终封口。**

Stage 2D-0 ~ 2D-6 保持 PASS / CLOSED；本阶段不修改核心 meshing 算法。

## 设计约束

2D-V 是薄层后处理：

- 只读取已导出的 `CM2D v1` / `quality.json` / `viz.json`；
- 不重新执行 geometry classification、Quadtree、Cut-cell 或 agglomeration；
- 核心 C++ library 不依赖 Python/renderer；
- 可视化失败不能改变网格结果。

## 实现

新增：

- `tools/visualization/render_cm2d.py`
- `tools/visualization/README.md`
- `tests/visualization_test.py`

CLI 额外导出：

- `<prefix>.viz.json`

`viz.json` 是只读展示 sidecar，记录已经由核心 pipeline 算出的 source-cell 信息：

- source id / key
- Quadtree level
- `full` / `cut`
- area fraction alpha
- background AABB
- centroid
- source small-cell flag/status
- source topology cell id
- small-cell target topology cell id（存在时）

它不包含新的 meshing 判断逻辑。

## SVG 显示

- final stabilized solver polygons
- adaptive level fill（或按 area 着色）
- internal edges
- embedded physical boundary
- domain boundary
- unclassified boundary warning style
- Cut/boundary final cells
- source background-cell bounds
- source small-cell（聚合前）虚线框 + centroid marker
- optional cell ids
- quality summary / topology audit
- invalid topology/quality 显式红色 banner

## 验收映射

`ACCEPTANCE_CN.md` 的 2D-V 要求：

- 仅从导出文件读取：满足；
- 不参与网格生成：满足；
- background/adaptive：source background bounds + source level；
- cut：`viz.json kind=cut` + embedded-boundary final cell；
- small-cell：source small-cell sidecar overlay，避免聚合后猜测；
- invalid：非零 CM2D audit 或 invalid quality 显式 banner；
- 可视化与核心库无反向依赖：满足。

## GitHub Actions run #14

Workflow run id: `32334947692`。

结果：**SUCCESS**。

全部步骤通过：

1. root CMake configure with 2D ON；
2. full native-2D target build；
3. Stage 0~6 的 15 个 CTest 全通过；
4. 四个 E2E `.vtk/.cm2d/.quality.json/.viz.json` 真实生成并校验；
5. `tests/visualization_test.py` 通过；
6. rectangle / circle / concave / airfoil-like 四个真实 CM2D 生成 SVG；
7. 每个 SVG 检查 `source-background / cells / edges / Quality summary` 图层；
8. SVG、CM2D、quality JSON、viz JSON 上传为 Actions artifact `cartmesh2d-acceptance-visualizations`；
9. root `CARTMESH_BUILD_2D=OFF` configure 通过。

## 当前门禁

- [x] dependency-free CM2D SVG renderer
- [x] adaptive level visualization
- [x] Cut/boundary visualization
- [x] source small-cell exact overlay
- [x] invalid-state visible flag
- [x] quality panel
- [x] renderer regression
- [x] four real E2E SVG renderings
- [x] artifact upload
- [x] 2D core regression remains green
- [ ] user-triggered final `验证-v` closeout

因此当前 2D-V 为 **READY FOR VALIDATION**，尚不标记 CLOSED。
