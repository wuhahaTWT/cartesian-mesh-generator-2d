# 开发导航

## 从哪里进入代码

生成链：输入轮廓 → 尺寸场/Quadtree → 真实流体 Cut-cell 或共形边界层 → 小单元处理 → 求解拓扑 → 质量 → 文件导出。

| 路径 | 责任与入口 |
|---|---|
| `apps/cartmesh2d_cli.cpp` | 纯 Cut-cell 总流程、尺寸场参数、物理面积门、Solver 质量和输出 |
| `apps/cartmesh2d_hybrid_cli.cpp` | 边界层总流程、局部修复开关、fallback 与导出 |
| `apps/cartmesh2d_dxf_cli.cpp` | DXF 导入命令行；`cartmesh2d_boundary_layer_cli.cpp` 是仍用于测试的独立边界层诊断工具 |
| `geometry/Geometry2D` | 基础几何、轮廓校验与内外关系 |
| `geometry/BoundarySimplification2D` | 输入轮廓简化及保形约束 |
| `geometry/ConstructionIdentity2D`、`ConstructionRecovery2D` | 稳定身份、来源追踪、类型化恢复请求；恢复 API 存在不代表所有生产路径已自动恢复 |
| `geometry/IntersectionRegistry2D`、`IntersectionConstruction2D` | 全局共享交点与格线构造；不要回到每格独立修补交点 |
| `grid/CartesianGrid2D`、`spatial/BoundarySegmentIndex2D` | 背景格、分类、边界空间索引 |
| `quadtree/Quadtree2D`、`sizing/SizeField2D` | 四叉树、2:1 平衡、距离/曲率/间隙/尾迹尺寸场 |
| `sizing/MeshResolution2D` | 最终求解拓扑的无量纲尺寸、壁面切向/法向测量与请求偏差 |
| `cutcell/CutCell2D` | 真正的二维流体多边形与多连通分量 |
| `boundary_layer/BoundaryLayer2D` | 壁面链、法向推进、层厚与局部停止 |
| `hybrid/HybridMesh2D`、`TransitionCanonicalization2D` | 边界层、过渡和余域统一拓扑，局部降层、终止、fallback |
| `topology/Topology2D`、`SharedEdgePartition2D` | 全局边、owner/neighbour、公共分割 |
| `topology/EdgeIncidence2D`、`PatchTransaction2D` | 边关联检查、局部修改事务和全局校验 |
| `stabilization/SmallCell2D`、`Agglomeration2D` | 小单元识别与聚合 |
| `quality/Quality2D`、`SolverQuality2D` | 基础质量与硬求解门 |
| `quality/SolverTopology2D`、`PatchLocalQuality2D` | 凸划分、源邻域相对评分、精确合并/切分批量与完整质量验收 |
| `io/Dxf2D`、`BoundaryMetadata2D` | CAD 曲线、单位、边界名称/角色 |
| `io/MeshIO2D`、`OpenFoam2D` | CM2D/VTK/JSON、二维挤出和 OpenFOAM case |

表中的模块分别位于 `include/cartmesh2d/` 和 `src/`，仍参与 CMake 构建，未把现用模块当旧版本删除。

桌面入口是 `desktop/src/main.js`（编排/IPC）、`preload.js`（受限桥接）。
`core/` 管能力、样例、几何/SVG 输入、参数、尺寸预算、CM2D、结果摘要；`process.js` 负责可取消/限时的原生子进程，`automatic.js` 保留旧版最多 4 组候选；`cell-budget.js` 负责数量初值与实测调整，`budget-runner.js` 最多 3 次执行并保留最接近目标的成功结果。

图片输入：`renderer/raster-import.js/css` 提供预览与强制标定，`raster-worker.js` 在 Worker 中调用共享 `core/raster.js`（RGBA 分割、四连通区域、边界追踪、保孔简化）。`core/raster-geometry.js` 检查编码尺寸、复核环交叉并把像素转换为米；main 的 `read-raster` 固定源文件快照，`commit-raster` 保存用户确认的 XY 和图片记录，再走已有生成链。原图不上传，不依赖外部模型。

`desktop/tests/raster*.test.js` 覆盖阈值端点、保孔/填孔、主体选择、复杂度拒绝、标定和坏图片头。实际 App 可使用 `--smoke=circle --image=/absolute/input.png --image-width=200 --out=/absolute/output --export=/absolute/result.zip` 验证 200 mm 标定至导出；`--image-fill=true` 填孔，`--image-preview-only=true` 检查未标定禁用确认，`--image-reject=true` 检查无轮廓拒绝与取消保留旧几何。
`renderer/app.js` 管表单和状态，`viewport.js` 管真实网格画布；HTML/`style.css` 管布局；`theme.css` 集中两套外观变量，`themes.js` 管主题切换与本地持久化；画布配色仍在 `viewport.js`，与界面主题独立。没有第二套前端。

## 构建与有限验证

macOS 使用 `/usr/bin/clang++`，不要使用 PATH 中的 mesasdk 编译器打包。根目录只保留一个 `build/`。

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/clang++
cmake --build build -j2
ctest --test-dir build --output-on-failure
npm test --prefix desktop
```

日常用 `ctest --test-dir build -R <相关测试名> --output-on-failure`；完整交付才统一全量运行。
测试目录包含解析、拓扑、守恒、质量、尺寸场、真实输出和确定性检查，不是可删除的缓存。

真实桌面流程（先打包准备 runtime；输出目录预先创建）：

```sh
mkdir -p outputs/smoke
cd desktop
node_modules/.bin/electron . --smoke=circle --out=../outputs/smoke --shot=../outputs/smoke/circle.png
```

可加 `--verified-preset=true`（真实界面载入已有高密案例参数）、`--density=150000`、`--target-cells=100000`、`--auto-padding=0.5`、`--theme=duet`、`--interaction-check=true`（主题状态与实际参数转手动）、`--allow-unsafe=true`、`--control=manual`、`--density=dense`、`--method=hybrid`、`--mode=light`、`--regions=1`、`--repeat=1`。省略 `--out` 即验证默认临时预览；`--export=/绝对路径/result.zip` 验证结果包。每轮会缩到最小窗口并检查底部可达、预览/导出单元数。smoke 会走真实表单/IPC/CLI 后退出；它不等于所有界面操作都已验收。

## CFD 验证

两路径的相对尺寸旗标：`--reference-length <m>`、`--wall-relative-size <h/Lref>`、
`--background-relative-size <h/Lref>`、`--far-field-spans <padding/Lref>`、`--cells-per-level <n>`。
Hybrid 另加 `--first-layer-relative-size <height/Lref>`；挤出厚度的旧位置参数仍为米。
相对尺寸启用后，旧位置参数的树层级和留白仅为兼容占位，几何、内外流和质量门仍正常执行。
`--size-field-only` 两路径均可预检；它不验证构造或 Solver 质量。

`tools/verification/check_mesh_resolution.py FINAL.solver.cm2d --report PREFIX.resolution.json`
独立测量序列化后的最终网格。`tests/mesh_resolution_cli_test.py --cli build/cartmesh2d_cli
--hybrid-cli build/cartmesh2d_hybrid_cli` 检查两路径的参考长度、缩放、实际输出、非法参数及伪造统计。
默认证据目录为忽略的 `outputs/engineering-resolution/`。

以下小例既保留原生输入，又有独立读取验证。OpenFOAM 使用薄层挤出，frontAndBack 为 empty。

```sh
mkdir -p outputs/check
build/cartmesh2d_cli examples/acceptance/circle.xy outputs/check/circle 6 0.25 0.1 exterior outputs/check/circle-case 3
python3 tools/verification/check_openfoam2d.py outputs/check/circle-case --report outputs/check/independent.json
```

Docker 已运行且镜像已安装时，真实检查命令：

```sh
docker run --rm --network none -v "$PWD/outputs/check/circle-case:/home/openfoam/workingDir/case" opencfd/openfoam-run:2606 checkMesh -case /home/openfoam/workingDir/case -allGeometry -allTopology
```

当前产品不执行 Q1 合同检测。`checkMesh`、内部拓扑和 Solver 质量各自独立；短求解也不等于所有工程流动均适用。不要放宽任一现有门槛。旧版 Q1 证据仅用于历史追溯。

## 保留工具的用途

| 工具 | 用途 |
|---|---|
| `check_face_planes.py` | 输入 `constant/polyMesh`；独立重算 OpenFOAM face-plane 标记，区分共面与正侧，不能替代整体 checkMesh；`--expected-set` 核对原标记集合 |
| `check_openfoam2d.py` | 独立读取 owner/neighbour、闭合与体积 |
| `check_boundary_layer2d.py`、`check_boundary_layer_failure.py`、`check_hybrid_mesh2d.py` | CTest 使用的边界层/共形拼接检查；专用 hybrid reader 不能冒充全部终止路径验收 |
| `check_openfoam_v1_logs.py` | 真正的 checkMesh/simpleFoam 日志验收 |
| `openfoam_harmonic_mms.py`、`check_openfoam_v1c.py` | 制造解、常数场、线性场验证；`render_openfoam_mms.py` 绘图 |
| `run_engineering_mms.py` | 固定相对尺寸三档，真实 OpenFOAM 单元中心/体积、检查、求解和误差；MMS 与扩展 checkMesh 分开判定 |
| `generate_q0_baselines.py` | 外部/求解器质量统计基线 |
| `refinement_ladder.py` | 加密压力阶梯；默认最低层级接近最大层级，不能作产品自适应性能代表 |
| `size_field_survey.py`、`alignment_sensitivity.py` | 新版真实尺寸场与几何/网格对齐敏感性测量 |
| `generate_benchmark_geometry.py` | 可复现基准输入 |
| `render_cm2d.py`、`render_boundary_layer2d.py`、`render_hybrid2d.py` | 读取真实产物绘图；hybrid sourceKey 不是树层级，应按面积画 |

除绘图工具外，上述工具位于 `tools/verification/`，参数用 `--help` 查看。
几何 fixtures 都在 `examples/`；`tests/repro/` 同时含审计反例与已接入测试的失败小案例；`source_halo_circle_patch.hpp` 保留源修复和最终重划分的真实邻域回归，不能把目录整体标成通过或失败。
DXF 反例可用 `build/cartmesh2d_dxf_cli tests/repro/hidden-spline-bulge.dxf outputs/hidden.xy 0.001 outputs/hidden.json` 复现。

## 独立质量审核

在 OpenFOAM 2606 环境内、案例副本目录中执行（先保留已有配置）：

```sh
cp "$WM_PROJECT_DIR/etc/caseDicts/meshQualityDict" system/meshQualityDict
checkMesh
checkMesh -allGeometry -allTopology
checkMesh -allGeometry -allTopology -meshQuality -writeAllFields -writeSets vtk -writeChecks json
```

三份日志分别保存，不只看进程退出码。原厂配置是可追溯对照，`minVol` 等有量纲参数仍依赖案例尺度；不要为获得 PASS 静默修改。`postProcessing/constant/` 的问题集合可用 ParaView 打开，质量场与 `checkMesh.json` 保留用于复核。当前实例、镜像 digest、输入哈希和结论见 `artifacts/current/quality-review.json`。当前产品不执行旧版 Q1 合同；Solver、拓扑和外部检查仍按各自证据报告。

当前目标修复在最终 hybrid 拓扑上进行：`improveSolverForTargetPolicy2D` 使用 65°，`improveSolverExtrudedDeterminant2D` 使用导出请求的厚度和 .001。它们保留固定层、共同拓扑与面积，并复查原 Solver 和方向连通；不是新的总评级。独立诊断命令：

```sh
python3 tools/verification/check_extruded_quality.py outputs/quality-fix/final-nozzle/foam/constant/polyMesh --output outputs/quality-fix/final-nozzle/extruded-quality.json
```

该脚本区分真正 positive-side 和共面面片，并复算 OpenFOAM 的全表面 determinant；不替代真实 `checkMesh` 或 CFD。固定喷管最新证据见 `artifacts/current/target-quality-and-sidebar.json`。

## 不能丢失的设计边界

- 共享格点构造不能由每 leaf 独立焦合替代；重建拓扑必须携带 EmbeddedBoundary 身份，不能只复制坐标后重新猜测。
- 纯路径 grid-corner 几何预算与输入端点的算术舍入预算分开；后者不能使用几何容差把端点拉过格线。
- 局部孔洞仍可能显式 unsupported；不能把它改成删掉小面积流体来通过。
- BoundaryLayer 六项仍是 OBSERVED，未定义完整评级阈值，不能写 PASS；当前产品不执行旧版 Q1 分类检测。
- Q3/Q4 的四个修复开关默认关闭且有累积关系；在窄缝上有有限收益，已饱和，不再自动扩展修复变体。
- Q5 径向规则及 R1 patch-local 事务仍参与构建和 CI，清理没有删除它们。sharp-tail 的 R1 迁移没有被证明可行。
- hybrid 高层级受壁面切向采样与法向层厚制约；简单缩首层或重新加点不能宣称解决。
- 同工具链的确定性与跨编译器逐字节一致是不同承诺。

## 历史找回，不复制进工作目录

本次清理前的完整新版保存在提交 `9bddfbf`。查某份旧文档：

```sh
git show 9bddfbf:docs/R2_HANDOFF_CN.md
git log --all --oneline -- src/quality/SolverTopology2D.cpp
```

| 旧入口 | 定位 |
|---|---|
| `703122a` / `codex/w1-shared-cutcell` | 旧桌面；与新版分叉，共同祖先 `874625a`；核心 W1 等价 cherry-pick 到 `6bcc458` |
| `origin/rescue/claude-local-post-pr6` | 早期 L10 实验及旧审计；不作为当前产品事实 |
| `archive/remote/q2b-constrained-local-repair` | 未合入的 Q2-B 方案，不能把其中的好数字冒充当前结果 |
| `archive/remote/solver-margin-hardening` | 独立质量优化实验；不自动混入新版 |
| 其他 `archive/*` | 拆仓前/历史检查点，包括三维历史；保留 Git 引用，不建立本地源码副本 |

Git 历史和标签没有清除。最近前端临时保存 `f43a51f` 涉及的十个文件已经进入新版，没有第三套更新前端需要保留。
本次盘点覆盖全部目录条目和项目文件用途，不等于对每条数值算法重新进行完整证明。
