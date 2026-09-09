# 开发导航

## 从哪里进入代码

生成链：输入轮廓 → 尺寸场/Quadtree → 真实流体 Cut-cell 或共形边界层 → 小单元处理 → 求解拓扑 → 质量 → 文件导出。

| 路径 | 责任与入口 |
|---|---|
| `apps/cartmesh2d_cli.cpp` | 纯 Cut-cell 总流程、尺寸场参数、物理面积门、纯路径 Q1 和输出 |
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
| `quality/QualityContract2D` | 无量纲分类型 Q1 合同与纯路径可达性诊断 |
| `quality/SolverTopology2D`、`PatchLocalQuality2D` | 凸划分、有限局部修复、patch-local 候选选择 |
| `io/Dxf2D`、`BoundaryMetadata2D` | CAD 曲线、单位、边界名称/角色 |
| `io/MeshIO2D`、`OpenFoam2D` | CM2D/VTK/JSON、二维挤出和 OpenFOAM case |

表中的模块分别位于 `include/cartmesh2d/` 和 `src/`，仍参与 CMake 构建，未把现用模块当旧版本删除。

桌面入口是 `desktop/src/main.js`（编排/IPC）、`preload.js`（受限桥接）。
`core/` 管能力、样例、几何/SVG 输入、参数、尺寸预算、CM2D、结果摘要；`process.js` 负责可取消/限时的原生子进程，`automatic.js` 管最多 4 组候选参数与粗略估时。
`renderer/app.js` 管表单和状态，`viewport.js` 管真实网格画布；HTML/`style.css` 管布局；`theme.css` 集中外观变量，画布配色仍在 `viewport.js`。没有第二套前端。

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

可加 `--control=manual`、`--density=dense`、`--method=hybrid`、`--mode=light`、`--regions=1`、`--repeat=1`。省略 `--out` 即验证默认临时预览；`--export=/绝对路径/result.zip` 验证结果包。每轮会缩到最小窗口并检查底部可达、预览/导出单元数。smoke 会走真实表单/IPC/CLI 后退出；它不等于所有界面操作都已验收。

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

`checkMesh` 通过不等于 Q1 通过；短求解也不等于所有工程流动均适用。不要放宽任一原有门槛。

## 保留工具的用途

| 工具 | 用途 |
|---|---|
| `check_openfoam2d.py` | 独立读取 owner/neighbour、闭合与体积 |
| `check_boundary_layer2d.py`、`check_boundary_layer_failure.py`、`check_hybrid_mesh2d.py` | CTest 使用的边界层/共形拼接检查；专用 hybrid reader 不能冒充全部终止路径验收 |
| `check_openfoam_v1_logs.py` | 真正的 checkMesh/simpleFoam 日志验收 |
| `openfoam_harmonic_mms.py`、`check_openfoam_v1c.py` | 制造解、常数场、线性场验证；`render_openfoam_mms.py` 绘图 |
| `generate_q0_baselines.py`、`generate_q1_baselines.py`、`verify_q1_scale_invariance.py` | CI 仍使用的质量统计、可重复性和尺度检查 |
| `refinement_ladder.py` | 加密压力阶梯；默认最低层级接近最大层级，不能作产品自适应性能代表 |
| `size_field_survey.py`、`alignment_sensitivity.py` | 新版真实尺寸场与几何/网格对齐敏感性测量 |
| `generate_benchmark_geometry.py` | 可复现基准输入 |
| `render_cm2d.py`、`render_boundary_layer2d.py`、`render_hybrid2d.py` | 读取真实产物绘图；hybrid sourceKey 不是树层级，应按面积画 |

除绘图工具外，上述工具位于 `tools/verification/`，参数用 `--help` 查看。
几何 fixtures 都在 `examples/`；`tests/repro/` 的两个小文件是尚未修复的审计反例，不是已通过的回归。
DXF 反例可用 `build/cartmesh2d_dxf_cli tests/repro/hidden-spline-bulge.dxf outputs/hidden.xy 0.001 outputs/hidden.json` 复现。

## 不能丢失的设计边界

- 共享格点构造不能由每 leaf 独立焦合替代；重建拓扑必须携带 EmbeddedBoundary 身份，不能只复制坐标后重新猜测。
- 纯路径 grid-corner 几何预算与输入端点的算术舍入预算分开；后者不能使用几何容差把端点拉过格线。
- 局部孔洞仍可能显式 unsupported；不能把它改成删掉小面积流体来通过。
- Q1 普通单元与 BoundaryLayer 分类统计；BoundaryLayer 六项是 OBSERVED，未定义完整评级阈值，不能写 PASS。
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
