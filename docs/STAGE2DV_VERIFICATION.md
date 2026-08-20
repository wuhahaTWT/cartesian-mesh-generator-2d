# Stage 2D-V 验证记录

## 状态

**PASS / CLOSED — 2D-V 已完成实现，并在用户显式 `验证-v` 后对当前 exact head 完成最终 GitHub Actions 封口验证。**

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

## 最终 exact-head GitHub Actions 验证

最终实现 head：`dd2546a1ac2a95ea5a7baf07fbebcfb1a9a9c287`

PR merge checkout：`fd52cb3288279f54f18f9ded6b1fa1f499cb150c`

Workflow run #16：`32335057216`

结果：**SUCCESS**。

全部步骤通过：

1. root CMake configure with 2D ON；
2. full native-2D target build；
3. Stage 0~6 的 15 个 CTest：**15/15 PASS**；
4. rectangle / circle / concave / airfoil-like 的 `.vtk/.cm2d/.quality.json/.viz.json` 真实生成并校验；
5. quality JSON `valid=true` 且 topology audit 全零；
6. `viz.json` format/source-cell/alpha/level/bounds 与 small-cell count 校验通过；
7. `tests/visualization_test.py` PASS，包括 invalid-state 显示回归；
8. 四个真实 CM2D acceptance mesh 均成功生成 SVG；
9. 每个 SVG 均检查 `source-background / cells / edges / Quality summary` 图层；
10. Actions artifact `cartmesh2d-acceptance-visualizations` 上传成功；
11. artifact 共 16 个文件，大小 104883 bytes，SHA256 `4d97b13578f1413fa95e3b86672c49d666e61dd667f8eab7636e973da839df3f`；
12. root `CARTMESH_BUILD_2D=OFF` configure PASS。

## 分支隔离

最终 `main...agent/native-2d-baseline`：

- status: ahead
- ahead_by: 113
- behind_by: 0
- merge base = main `8bec26d98eb8bd84033625ed2a41184c8cb223f1`

改动保持在二维入口/文档/CI 与 `cartmesh2d/**`；没有修改三维算法目录 `include/cartmesh/**`、根 `src/**`、根 `apps/**`、根 `tests/**`。

## 最终门禁

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
- [x] exact-head final validation
- [x] branch isolation audit
- [x] user-triggered final `验证-v` closeout

因此 Stage 2D-V 正式记为 **PASS / CLOSED**。
