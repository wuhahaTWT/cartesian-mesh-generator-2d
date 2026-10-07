# 开发指南

当前能力与待办只看[当前状态](CURRENT_STATE_CN.md)，操作看[桌面使用](DESKTOP_APP_CN.md)，修改约束看[AGENTS](../AGENTS.md)。本页集中构建、代码入口、方法与复现，不重复逐轮实验结果。

独立 Python 验证链及其专用测试、CI 调用已移除；运行流程使用原生网格和求解器的失败判定。下文仍保留的独立审计方法与数值记录描述的是历史实验，不能作为当前自动运行的能力；原始场和历史报告保留。

## 本地研究文件导航

| 文件或目录 | 用途 |
| --- | --- |
| `tools/optimization/navier_stokes_brinkman.py`、`brinkman.py`、`optimize_flow.py` | 多孔流动、离散伴随和材料优化；T01/T02 共用基础。 |
| `tools/optimization/fidelity_study.py`、`fidelity_metrics.py`、`fidelity_report.py`、`topology_cases.py`、`engineering_baselines.py` | T01 算例、差距/排序实验、工程参照及汇总。 |
| `tools/optimization/closed_loop_topology.py` | T02 在线反馈、材料/连接提案、局部校正与原生接受循环。 |
| `tools/optimization/compare_sharp_designs.py`、`native_flow.py`、`continue_native_flow.py` | 真实壁面提取、面积/物性/流量匹配、原生复算和续算。 |
| `tools/visualization/render_topology_feedback.py` | 从已接受的真实网格和 CSV 生成 T02 图及来源索引。 |
| `tests/navier_stokes_topology_test.py`、`tests/topology_fidelity_test.py`、`tests/closed_loop_topology_test.py`、`tests/fixtures/topology_area_plateau.npz` | 相关实现检查和保留的失败材料。 |
| `artifacts/current/native-topology-fidelity-completion*`、`native-topology-fidelity-fixed-ranking.png`、`native-topology-fidelity-native-fields*` | T01 补齐后的结果入口；不把原首轮 `native-topology-fidelity.json` 当最终全集。 |
| `artifacts/current/native-topology-feedback*` | T02 主结果、对照及实际网格图。 |
| `outputs/topology-fidelity/` | T01 原材料、原生网格/场、失败、独立审查及补齐脚本；交付入口 `t01-completion-001/summary.json`，原封存清单 `t01-completion-001/evidence-manifest.json`。 |
| `outputs/topology-feedback/double-pipe-004/`、`double-pipe-no-feedback-001/` | T02 Re=0 主运行和匹配对照，各读 `run.json`、`evaluations.json`、`summary.json`。 |
| `outputs/topology-feedback/double-pipe-re50-001/`、`double-pipe-re50-no-feedback-001/` | Re=50 主运行和匹配对照。 |
| `outputs/topology-feedback/double-pipe-001/` 至 `003/`、`elbow-001/` 至 `003/` | 早期算法尝试及弯管失败/试运行，保留用于改进；不能和最终方法混作同一组统计。 |
| `outputs/topology-env/`、`build/` | 现有研究 Python 环境和原生构建；系统 Python 未安装此流程所需的 SciPy，使用 `outputs/topology-env/bin/python`。 |
| `outputs/codex_tasks/` | Claude 原路线、T01 任务单、进度和过程日志；历史指导不覆盖当前状态。旧自动启动脚本在 `retired_launchers/*.sh.disabled`，不作为日常入口。 |
| `outputs/repository-tidy/workspace-index.json` | 本次接手用的本地文件清单：未提交路径、结果角色、读取入口；是整理时快照，接手仍先核对 Git。 |

研究结果保留原目录，因为运行记录、图片来源和封存清单包含原路径。`outputs/combustion-foundation/` 是冻结燃烧研究，`outputs/cfd-demo-reference/` 是用户指定保留的参考包，其余网格/CFD 历史材料继续保留。源码目录杂项的原路径与恢复位置见 `outputs/repository-tidy/file-organization.json`；不把数值场、检查点或失败日志当缓存。

## 分支与里程碑

只维护 `main`（已验证集成）和 `codex/cfd-development`（后续开发），网格和 CFD 修复统一从开发线推进。三条实验功能、临时集成及旧网格维护分支已完整进入 main，旧分支名可删除，提交历史仍可查。

`mesher-v0.3.0` 固定指向网格里程碑 `691c97e`，不移动。切换前先处理当前改动，切换后按需要重建原生程序和 runtime，避免源码与旧包混用。集成须按当次授权及验证范围进行，不自动合并后续功能。

## 构建与有限验证

依赖 C++20、CMake 3.20+、Node.js 22、Python 3。macOS 使用系统 `/usr/bin/clang++`，禁止 PATH 中的 mesasdk 编译器混入 App；Windows 使用 Visual Studio 2022 C++ 工作负载及 Windows SDK，Linux 使用 GCC 11+ 或兼容 Clang。常规构建只用根目录 `build/`。

```sh
npm ci --prefix desktop
npm --prefix desktop run build:native
cmake --build build --config Release --parallel 2
```

`build:native` 准备本机六份桌面 CLI、13 个样例及 runtime manifest；完整 CMake 构建还覆盖独立 CLI 和测试。只开发原生代码时，macOS 可直接配置：

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/clang++
cmake --build build --parallel 2
```

日常用 `ctest --test-dir build -C Release -R '<相关测试名>' --output-on-failure` 跑受影响项。重要交付和功能集成运行：

```sh
ctest --test-dir build -C Release --output-on-failure
npm test --prefix desktop
```

开发启动用 `npm --prefix desktop start`；按本机系统打包：

```sh
npm --prefix desktop run pack:mac
```

包输出到 `desktop/dist/`。打包入口统一检查实际 Mach-O/PE/ELF、架构及 runtime，拒绝跨系统混装；Windows 使用静态 CRT 和 UTF-8 路径 manifest，macOS 检查动态库。ZIP 使用流式 yazl，不依赖外部压缩命令。

三平台 CI 保留构建、原生测试、前端测试和打包；额外的打包 App 审计脚本已移除。真实 App 交互按交付目标手动检查。

单例桌面调试（先准备 runtime）：

```sh
mkdir -p outputs/smoke
cd desktop
node_modules/.bin/electron . --smoke=circle --out=../outputs/smoke --shot=../outputs/smoke/circle.png
```

可选 `--export=/绝对路径/result.zip`、`--method=hybrid`、`--verified-preset=true`、`--target-cells=100000`、`--interaction-check=true`；图片路径用 `--image=/绝对路径/input.png --image-width=200` 检查 200 mm 标定。其他开关以 `desktop/src/` 的实际 smoke 入口为准；走通 smoke 不等于所有交互已验收。

## 从哪里进入代码

生成链：输入 → 尺寸场/Quadtree → Cut-cell 或共形边界层 → 小单元处理 → 共享面拓扑 → 质量 → 导出。

| 路径 | 责任 |
| --- | --- |
| `apps/cartmesh2d_cli.cpp` / `cartmesh2d_hybrid_cli.cpp` | Cut-cell / 边界层总流程、面积/质量门、fallback 与导出 |
| `apps/cartmesh2d_flow_cli.cpp` / `cartmesh2d_transport_cli.cpp` | 不可压层流 / 被动标量与温度 |
| `apps/cartmesh2d_euler_cli.cpp` / `cartmesh2d_immersed_cli.cpp` | 可压层流 / 独立浸入边界实验 |
| `apps/cartmesh2d_fv_cli.cpp` / `cartmesh2d_dxf_cli.cpp` | 扩散/泊松 / DXF 导入；独立 boundary_layer CLI 仍用于诊断测试 |
| `geometry/Geometry2D`、`BoundarySimplification2D` | 二维基础几何、输入校验及保形简化 |
| `geometry/ConstructionIdentity2D`、`ConstructionRecovery2D`、`IntersectionRegistry2D`、`IntersectionConstruction2D` | 共享构造身份、交点、来源与恢复 API；不表示所有路径已自动恢复 |
| `grid/CartesianGrid2D`、`spatial/BoundarySegmentIndex2D` | 背景格、分类与边界索引 |
| `quadtree/Quadtree2D`、`sizing/SizeField2D`、`MeshResolution2D` | 平衡、距离/曲率/间隙/尾迹加密及最终切向/法向分辨率 |
| `cutcell/CutCell2D`、`boundary_layer/BoundaryLayer2D`、`hybrid/HybridMesh2D`、`TransitionCanonicalization2D` | 真流体多边形、多分量、贴体层、局部终止与统一拓扑 |
| `topology/Topology2D`、`SharedEdgePartition2D`、`EdgeIncidence2D`、`PatchTransaction2D` | owner/neighbour、公共分割、边关联与修改事务 |
| `stabilization/SmallCell2D`、`Agglomeration2D` | 小单元识别与守恒聚合 |
| `quality/Quality2D`、`SolverQuality2D`、`SolverTopology2D`、`PatchLocalQuality2D` | 质量门、凸划分、局部修复和完整复核 |
| `fv/FvMesh2D`、`Diffusion2D`、`Incompressible2D`、`ThermalFlow2D`、`ThermalCheckpoint2D` | 有限体积几何、方程、共享通量、联合时间与续算 |
| `io/Dxf2D`、`BoundaryMetadata2D`、`MeshIO2D`、`OpenFoam2D`、`BackgroundGridIO2D` | CAD/边界角色、网格和 OpenFOAM/背景导出 |

模块分别位于 `include/cartmesh2d/` 和 `src/`。桌面只有 `desktop/src/` 一套前端：

- `main.js` / `preload.js` 编排 IPC；`core/process.js` 管取消/限时，`automatic.js` / `cell-budget.js` / `budget-runner.js` 管自动数量选择；`flow-case.js` / `flow-checkpoint.js` / `thermal-job.js` 管工况、续算和联算。
- `renderer/raster-import.js`、`raster-worker.js`、`core/raster.js` / `raster-geometry.js` 完成本地图片分割、保孔轮廓、标定和校验；main 固定源文件快照，确认后再进入原生链。
- `renderer/app.js` 管状态，`viewport.js` 画真实网格；`theme.css` / `themes.js` 管主题，`style.css` 管布局，画布色表在 `viewport.js`。
- `assets/fonts/cartmesh-ui-regular.woff` 来自 Noto Sans CJK SC 2.004；许可与原哈希在同目录 LICENSE.txt。`desktop/scripts/subset-ui-font.py` 用 fontTools 重建，常规构建不下载字体；生僻字由系统字体回退。

## 困难曲壁的 Solver 修复

`SolverTopology2D.cpp` 在原有成对修复停滞后，允许满足既有单元质量门的轻微凹分块，并搜索问题邻域内相连的三/四单元精确合并。候选必须保留全部外侧原子边，经过完整相邻 halo 的质量评分、不可变单元保护和全网格拓扑/质量复核；不修改原始 XY，不降低任何质量阈值。只在原成对路径失败后增加有界搜索，避免在全网格枚举组合。

凸修复停滞时先从**已接受的当前拓扑**继续一般修复，再以原始起点为备选；按同一全局质量评分保留更好的结果。尚有缺陷的部分改进仍被 CLI 最终质量门拒绝。8/10 格真实失败缩减例位于 `tests/repro/annulus_quality_patch.hpp`，覆盖缩放、旋转、原子边、面积、来源谱系、不可变邻域及确定性；壁面缩减例的人工截断邻域故意保留失败，不作为完整 CFD 网格。

普通入口复现（两环各 1,024 段，不依赖旧的手工修复网格）：

```sh
build/cartmesh2d_cli artifacts/current/annulus-1024-axis-roundoff.xy outputs/annulus 7 .15 .1 interior outputs/annulus-foam 7
build/cartmesh2d_flow_cli --mesh outputs/annulus.solver.cm2d --case annulus --speed .5 --export-boundaries outputs/annulus.boundaries
build/cartmesh2d_flow_cli --mesh outputs/annulus.solver.cm2d --case custom --boundary outputs/annulus.boundaries --nu .1 --speed .5 --convection face-limited-linear --steady-acceleration newton-krylov --pressure-preconditioner cholesky --tolerance 1e-6 --max-iterations 1500 --output outputs/annulus-flow
```

`cholesky` 是 macOS 的已有选项，跨平台构建不据此获得验证。普通网格、认证修复前后、默认线性求解器及真实 App 的完整记录见[入口证据](../artifacts/current/native-laminar-mesh-entry.json)。原生读回证明面积和外边界保持；当前 CM2D 格式不序列化完整 `sourceLineage`，故完整来源集合只在内存中的原生缩减例上核对，不能由两个空集合声称全网格谱系已验证。OpenFOAM 导出成功不等于外部 `checkMesh` 通过。

### 云端原生精度与控制对照

`artifacts/current/native-laminar-accuracy-probe.cpp` 复用现有原生测试的二维网格构造，调用产品求解器并导出真实单元场；不是独立离散方程审计。`native-laminar-accuracy-postprocess.py` 只读取这些场、解析 Couette/通道参考和仓库既有方腔文献数据，计算误差与插值，不建立第二套求解器或新验收门。

```sh
g++ -std=c++20 -O2 -Iinclude artifacts/current/native-laminar-accuracy-probe.cpp build/libcartmesh2d_fv.a build/libcartmesh2d.a -o outputs/cloud-laminar/accuracy-probe
outputs/cloud-laminar/accuracy-probe channel 64 1e-8 upwind default outputs/cloud-laminar/channel-upwind-64
```

相同入口支持 `manufactured`（变形网格解析强迫涡）和 `cavity`（Re=100），格式为 `upwind` 或 `face-limited-linear`，方法为 `default`、显式 `none` 或研究用显式 `anderson`。原始 `accuracy-runs.json`、`accuracy-controls.json`、`cost-runs.json` 保留每次完整命令、进程退出状态和耗时；场、输入和日志在同目录。完成这些实际运行后执行 `python3 artifacts/current/native-laminar-accuracy-postprocess.py` 生成当前精度证据。压力为运动学压力，封闭流比较去除体积加权常数；速度/压力误差分别以 m/s 和 m²/s² 报告，不混成单个场误差。通道的参考为 `u=4y(1-y), v=0, p/ρ=8ν(4-x)`，方腔按真实中心网格和规定壁值双线性插值；圆环包含多边形几何误差，不因残差小而授予精度资格。


曲壁几何/空间联合细化由 `python3 artifacts/current/native-laminar-curve-study.py` 复现（可用 `--levels 4 5 6 7` 和 `--schemes upwind face-limited-linear`）。驱动生成独立参数化圆环输入，不修改原 1024 段 XY；按原 CLI 流程生成 Solver 网格、导出 annulus 边界，再以 `custom` 默认 Newton 求解。OpenFOAM 目标需提供目录：传 `-` 只构造源网格，不产生 `.solver.cm2d`；`--case annulus` 仅用于导出边界。首次驱动误用这两个入口造成的失败记录保留，修正后实际八次求解通过。`outputs/cloud-laminar/curve-joint-refinement.json` 记录网格、边界和完整求解成本；超时显式记失败，不能宣称该档完成。误差相对圆形解析解，含多边形几何/壁速差异；全域 RMS 与最大值同时报告，固定 `.6≤r≤.9 m` 内域指标只是定位工具，不能删去近壁误差。

```sh
g++ -std=c++20 -O2 -Iinclude artifacts/current/native-laminar-pressure-probe.cpp build/libcartmesh2d_fv.a build/libcartmesh2d.a -o outputs/cloud-laminar/pressure-probe
outputs/cloud-laminar/pressure-probe outputs/cloud-laminar/annulus.solver.cm2d outputs/cloud-laminar/annulus-default.cells.csv > outputs/cloud-laminar/pressure-probe-baseline.csv 2> outputs/cloud-laminar/pressure-probe-trace.json
```

最终单元均衡试验使用 `native-laminar-mesh-balance-probe.cpp`，只调用产品内已有的 `improveSolverForTargetPolicy2D`；它不是独立方程求解器，也不修改默认 Solver `minVolumeRatio=0.01`。第三个参数是无量纲目标邻格面积比 `min(Ao,An)/max(Ao,An)`。探针从最终网格的 `DomainBoundary` 恢复外域矩形；没有外域边界的内流按生成时 `.15×span` 留白重建仅供事务几何分类。可在已有 Release 静态库上直接编译：

```bash
g++ -std=gnu++20 -O2 -DNDEBUG -Wall -Wextra -Wpedantic -Wconversion -Wshadow \
  -Iinclude artifacts/current/native-laminar-mesh-balance-probe.cpp \
  build/libcartmesh2d.a -o outputs/cloud-laminar/mesh-balance-probe

outputs/cloud-laminar/mesh-balance-probe \
  outputs/cloud-laminar/annulus-joint-7.solver.cm2d \
  outputs/cloud-laminar/annulus-joint-7.xy .075 \
  outputs/cloud-laminar/annulus-joint-7-volume-075.solver.cm2d
```

追加一个或多个单元号可把其余单元标成 immutable，用于定位单个事务；结构成功不等于目标策略全部通过，必须同时检查探针 JSON 的 `outputDefaultPass` / `outputTargetPass`，实际流场仍须通过原 CLI 最终认证。当前三档圆环和三档圆柱结果见 `native-laminar-mesh-balance.json`：物理边界原子边逐坐标/patch 比较完全相同、面积守恒并通过默认门，但 `0.10` 比 `0.075` 更差，且圆柱积分量发生约 2.4% 系统移动，所以 `0.075` 只能保留为诊断，禁止据此改默认或宣称精度资格。

Re=20 外流圆柱的可信积分量基准必须同时细化多边形和近壁/背景网格；固定 32 边多边形只加密远场不能授予空间资格。原生编排驱动只调用产品 CLI 和读取原生结果，不建立第二套方程：

```sh
python3 artifacts/current/native-laminar-curved-benchmark.py \
  --levels 0 1 2 --schemes upwind face-limited-linear
```

三档分别为 32/64/128 段、近壁相对尺寸 `.125/.0625/.03125`、背景相对尺寸 `1/.5/.25`，外域为 10D；`outputs/cloud-laminar/cylinder-joint/runs.json` 保存六次实际命令、退出状态、严格残差、进程耗时和场哈希。域敏感性沿用最细轮廓与尺寸，把外域扩到 20D：

```sh
build/cartmesh2d_cli outputs/cloud-laminar/cylinder-joint/level-2.xy \
  outputs/cloud-laminar/cylinder-joint/far-20 8 .25 .1 exterior \
  outputs/cloud-laminar/cylinder-joint/far-20-foam 0 0 --size-field \
  --reference-length 2 --wall-relative-size .03125 \
  --background-relative-size .25 --far-field-spans 20 --cells-per-level 3

build/cartmesh2d_flow_cli \
  --mesh outputs/cloud-laminar/cylinder-joint/far-20.solver.cm2d \
  --case external --nu .1 --speed 1 --convection face-limited-linear \
  --tolerance 1e-8 --max-iterations 2500 \
  --output outputs/cloud-laminar/cylinder-joint/far-20-face-limited-linear
```

这里 `D=2 m`、`U=1 m/s`、`ν=.1 m²/s`，所以 `Re=20`；单位密度/单位深度下 `Cd=Fx/(.5U²D)=Fx`。Tritton 的 `Cd=2.045` 只作上下文参考，必须与几何、网格和外域敏感性一起解释，不能把有限方域与实验条件视为完全相同。六次联合细化的进程耗时是驱动实测；首个 20D face-limited-linear 运行只留有网格日志与最终摘要时间戳，约 1,846 秒间隔不是精确独占进程成本。随后同一 66,912 格网格实际补跑旧省略选项的 Upwind 默认：657 次完整评估、1,353.433 秒、`Cd=2.17429`，比参考高 `6.322%`；face-limited-linear 为 904 次评估、`Cd=2.06213`，误差 `0.838%`，完整评估成本增加 `37.6%`。这一完整成本/精度取舍与其他内流、封闭流、解析强迫流证据共同支持下述上下文默认改动；它不解决圆环局部压力。

整周局部压力复核使用 `artifacts/current/native-laminar-cylinder-pressure-profile.cpp`。它读取 `cartmesh2d_flow_cli` 已写出的真实 wall-face 压力与 CM2D 面几何，不组装新方程；每个闭合物面先减去壁长加权压力常数，并以壁长加权周界中心为极角原点，再在 16,384 个等角位置周期插值。不能用全局原点：当前圆柱中心是 `(0.07,0.03) m`，否则会引入平移相位误差。减去常数只固定压力 gauge，闭合物体压力合力不变。复现 10D 三档并补齐 20D 粗、中档：

```sh
g++ -std=c++20 -O2 -Wall -Wextra -Werror -Iinclude \
  artifacts/current/native-laminar-cylinder-pressure-profile.cpp \
  build/libcartmesh2d_fv.a build/libcartmesh2d.a \
  -o outputs/cloud-laminar/cylinder-pressure-profile

outputs/cloud-laminar/cylinder-pressure-profile \
  level0-10D outputs/cloud-laminar/cylinder-joint/level-0.solver.cm2d outputs/cloud-laminar/cylinder-joint/level-0-face-limited-linear.faces.csv \
  level1-10D outputs/cloud-laminar/cylinder-joint/level-1.solver.cm2d outputs/cloud-laminar/cylinder-joint/level-1-face-limited-linear.faces.csv \
  level2-10D outputs/cloud-laminar/cylinder-joint/level-2.solver.cm2d outputs/cloud-laminar/cylinder-joint/level-2-face-limited-linear.faces.csv \
  > outputs/cloud-laminar/cylinder-pressure-profile-10d.json

python3 artifacts/current/native-laminar-large-domain-pressure.py --levels 0 1
outputs/cloud-laminar/cylinder-pressure-profile \
  level0-20D outputs/cloud-laminar/cylinder-joint/far-20-level0.solver.cm2d outputs/cloud-laminar/cylinder-joint/far-20-level0-face-limited-linear.faces.csv \
  level1-20D outputs/cloud-laminar/cylinder-joint/far-20-level1.solver.cm2d outputs/cloud-laminar/cylinder-joint/far-20-level1-face-limited-linear.faces.csv \
  level2-20D outputs/cloud-laminar/cylinder-joint/far-20.solver.cm2d outputs/cloud-laminar/cylinder-joint/far-20-face-limited-linear.faces.csv \
  > outputs/cloud-laminar/cylinder-pressure-profile-20d.json
```

修正极角中心后，10D 三档联合加密的相邻压力剖面 RMS 差为 `0.03758→0.02034 m²/s²`，最大差为 `0.1053→0.07249 m²/s²`。同一 20D 外域的三档为 4,620/17,196/66,912 个最终单元，压力剖面相邻 RMS 差 `0.02736→0.01644 m²/s²`，而最大差 `0.07273→0.07769 m²/s²` 未下降；积分阻力 `2.06531→2.06207→2.06213`，中/细档相对变化仅 `0.00291%`。匹配 10D/20D 的三档剖面 RMS 差为 `0.02190/0.02441/0.02302 m²/s²`，不会随空间细化消失。该圆柱完整规定速度迹跳和壁面法向规定速度均为零，故它是圆环迹角点之外的一般曲壁压力回归；结论只支持积分载荷稳定，不授予局部点值压力空间精度资格。完整编排、严格残差、成本和哈希见 `native-laminar-large-domain-pressure.json`。

固定几何反例用 `artifacts/current/native-laminar-fixed-cylinder-pressure.py` 逐字节复用 `cylinder-joint/level-2.xy` 的 128 段静止圆柱和 20D 外域，只改变近壁/背景相对尺寸；脚本只调用原生网格、流动 CLI 与上述压力剖面后处理，不实现独立方程。粗、中档新建网格并求解，细档复用同轮廓的既有 `far-20` 接受场：

```sh
python3 artifacts/current/native-laminar-fixed-cylinder-pressure.py --grids 0 1

build/cartmesh2d_flow_cli \
  --mesh outputs/cloud-laminar/cylinder-joint/far-20-fixed128-grid1.solver.cm2d \
  --case external --nu .1 --speed 1 --convection upwind \
  --tolerance 1e-8 --max-iterations 2500 \
  --output outputs/cloud-laminar/cylinder-joint/far-20-fixed128-grid1-upwind

g++ -std=c++20 -O2 -Wall -Wextra -Werror -Iinclude \
  artifacts/current/native-laminar-mesh-balance-probe.cpp build/libcartmesh2d.a \
  -o outputs/cloud-laminar/native-laminar-mesh-balance-probe
outputs/cloud-laminar/native-laminar-mesh-balance-probe \
  outputs/cloud-laminar/cylinder-joint/far-20-fixed128-grid1.solver.cm2d \
  outputs/cloud-laminar/cylinder-joint/level-2.xy .075 \
  outputs/cloud-laminar/cylinder-joint/far-20-fixed128-grid1-volume-075.solver.cm2d \
  > outputs/cloud-laminar/cylinder-joint/far-20-fixed128-grid1-volume-075-repair.json
build/cartmesh2d_flow_cli \
  --mesh outputs/cloud-laminar/cylinder-joint/far-20-fixed128-grid1-volume-075.solver.cm2d \
  --case external --nu .1 --speed 1 --convection face-limited-linear \
  --tolerance 1e-8 --max-iterations 2500 \
  --output outputs/cloud-laminar/cylinder-joint/far-20-fixed128-grid1-volume-075-face-limited-linear

python3 artifacts/current/native-laminar-fixed-cylinder-pressure.py --skip-runs
```

4,716/17,260/66,912 格 face-limited-linear 三档均通过原 Solver 门和 `1e-8` 严格门，但中档出现最大速度 `5.49 m/s`、运动学压力范围 `[-182.6,83.2] m²/s²`，而粗/细档最大速度均约 `1.16 m/s`。中档压力/黏性 X 向壁载荷为 `-13.26/+15.26 m³/s²`，绝对分量之和除以总载荷的无量纲抵消因子为 `14.21`；显式 Upwind 对照也为 `13.69`，所以不是 Newton 或 face-limited-linear 特有。`native-laminar-mesh-balance-probe.cpp` 的既有精确合并/重分因果对照只接受 24 个事务，保持流体面积、716 条边界原子边和原门不变，将内部邻格面积比 `0.0186→0.0777`；修复网格的最大速度为 `1.16 m/s`、载荷抵消因子为 1。这个实验只证明最终单元关系是因果因素；由于旧跨工况扫描已证明 `.075` 非单调且可系统移动积分载荷，不能把实验目标抄成产品默认。完整命令、进程成本、场/面/残差哈希、壁压剖面和修复不变量在 `native-laminar-fixed-cylinder-pressure.json`，原始大场保留于恢复包。

最终接受场幅值由 `FlowResult2D::fieldAmplitude` 在原生最终单元数组上一次扫描得到，不重新求解或重构场。CLI 摘要成组写出极值、位置、`maximumSpeedRatio=max|U|/Uref` 和 `pressureRangeRatio=(pmax-pmin)/Uref²`；压力范围不依赖压力 gauge。桌面只在整组字段存在时校验并显示，旧摘要保持兼容。两项无量纲量是诊断，不是新的通过阈值；物理资格仍须按工况参考与加密判断。

局部因果试验从中档接受场中取 `|U|>1.5Uref` 的 28 个异常格，仅允许其两环/四环图邻域参与同一精确合并/重分，再用原默认 Newton + face-limited-linear 和 `1e-8` 门复算。两环接受 10 个事务、四环 8 个事务，均保持总面积、716 条边界原子边和原 Solver 门；但峰值从上游移到下游同构小格，最大速度约 `1.44 m/s`、压力范围约 `[-19.2,43.7] m²/s²`，没有达到全周 24 事务对照的正常分支。这是定位实验而非产品选择规则；异常场不能作为未来网格器的输入判据。原始局部网格、场和日志保留在 `outputs/cloud-laminar/cylinder-joint/`，摘要与哈希见 `native-laminar-field-amplitude.json`。

固定圆柱的分支复核由 `native-laminar-external-boundary.cpp` 按产品 external 规则把最终 716 条边界原子面写成 custom 配置；同一平坦初值的 custom 与 external cells/faces/residuals 三份数组逐字节相同。这样可在不修改公开预设接口的情况下调用 custom-only `momentum-inertia=0` Stokes，并把其单元场作为原生 `--initial-guess`。一步 Stokes 初值最终仍产生压力/黏性载荷大数抵消，只是落到不同支路。随后保持目标 Re=20 方程和所有门槛不变，用 `ν=1,.5,.2,.1 m²/s` 的已收敛单元场依次作为下一档初值；中档四阶段完整评估为 `715+509+618+473=2315`，终点得到 `max|U|/Uref=1.16451`、压力范围比 `1.79746`、总阻力 `2.06416` 的有界场。粗网格同一路径终点与直接解的 `u/v/p` 最大差小于 `5.7e-9`，但这仍只是两网格分支证据，尚未证明任意不利网格上的唯一性或默认成本合理。

缩短路径复用同一原生 CLI 和同一 `1e-8` 严格门，只改变前一级已接受场。平坦启动 `.2→.1` 的两级均严格收敛，却保留高幅值支路，总成本 `680+467=1147`；平坦启动 `.5→.1` 和 `1→.1` 分别以 `735+762=1497`、`715+772=1487` 次回到四级路径的同一有界终点。后者是当前最低已验证有界成本，比四级少 `35.8%`，但仍为直接异常解成本的 `1.42` 倍。正常粗网格的 `1→.1` 为 `275+262=537` 次，终点与直接解的全单元面积加权速度 RMS 差 `1.15e-9 Uref`、去规范压力 RMS 差 `1.99e-10 Uref²`。后处理的 `normalized_field_distance` 使用所有单元面积权重，压力差仅移除全域面积均值以消除不可压规范；不删格，也不以该差值决定哪条支路正确。生成初值 CSV 时保留目标网格逐行 `cell,x,y,u,v,p`，不作空间插值。

引导成本对照另采用不依赖案例幅值的容差关系：`νguide=10νtarget`，`tolguide=sqrt(toltarget)`，终档仍使用未经放宽的目标方程和目标容差。`toltarget=1e-8` 时，中档反例的两级成本由严格引导的 1,487 次降到 `459+823=1,282`，正常粗档由 537 次降到 `155+263=418`；反例终点相对严格引导终点的全场速度 RMS/最大向量差为 `1.15e-8/5.63e-7 Uref`，去规范压力 RMS/最大差为 `7.00e-9/1.59e-5 Uref²`。产品默认 `toltarget=1e-6` 另做独立真实运行：反例平坦直解 792 次仍落在高幅值支路，引导 `398` 次加终档 `606` 次回到有界支路，总成本为直解的 `1.268` 倍；正常粗档直解 204 次，引导加终档 `127+178=305` 次，两个终点的速度 RMS 差 `1.62e-7 Uref`、去规范压力 RMS 差 `2.47e-8 Uref²`。

`native-laminar-guide-certificate.cpp` 在同一进程、逐字节相同的最终网格上运行直接目标、平方根容差高黏度引导和从引导场启动的原目标；保留三份完整单元场，以全部单元面积权重比较速度，并只移除压力差的面积均值。它不建立独立方程、误差门或分支选择规则。64 档直通道在默认 `1e-6` 下的完整评估为 `38` 对 `8+60`，目标终点速度/去规范压力 RMS 差 `1.68e-6/4.34e-6`；Re=100 方腔为 `195` 对 `43+196`，差 `4.47e-6/1.08e-6`。同容差的平坦 Anderson 对照分别为 39/1091 次，因此直接 Newton + Anderson 的双路径成本为 77/1286 次，高黏度证书的 68/239 次更低。运行和再生成命令为：

```sh
g++ -std=c++20 -O2 -Iinclude \
  artifacts/current/native-laminar-guide-certificate.cpp \
  build/libcartmesh2d_fv.a build/libcartmesh2d.a \
  -o outputs/cloud-laminar/guide-certificate
outputs/cloud-laminar/guide-certificate channel 64 1e-6 \
  outputs/cloud-laminar/guide-channel-64
outputs/cloud-laminar/guide-certificate cavity 64 1e-6 \
  outputs/cloud-laminar/guide-cavity-64
outputs/cloud-laminar/accuracy-probe-anderson channel 64 1e-6 \
  face-limited-linear anderson \
  outputs/cloud-laminar/channel-face-limited-linear-64-anderson-tol1e6
outputs/cloud-laminar/accuracy-probe-anderson cavity 64 1e-6 \
  face-limited-linear anderson \
  outputs/cloud-laminar/cavity-face-limited-linear-64-anderson-tol1e6
python3 artifacts/current/native-laminar-branch-continuation.py
```

这组平方根关系只减少探针成本，不构成物理解选择原理；当前不接入 API/CLI 或默认，也不把幅值诊断升级为拒绝阈值。原始命令输出、各路径完整场哈希和所有规范不变差异由同一 `native-laminar-branch-continuation.py` v5 生成。

同一原生驱动现在也接受 `manufactured`。64² 解析强迫涡在默认 `1e-6` 下，直接目标、平方根引导、从引导启动的原目标分别用 550/145/61 次完整评估；直接与引导目标的全场速度/去规范压力面积 RMS 差为 `7.72e-5/7.40e-5`，平坦 Anderson 110 次所得目标与直接 Newton 的对应差为 `4.63e-5/4.33e-5`。这说明高黏度引导可同时充当收敛预处理，但差值仍只按目标容差解释，不升级为物理误差界。

移动圆环由 `native-laminar-annulus-guide.py` 调用产品 CLI：逐字节相同的 `annulus-joint-7.solver.cm2d` 和原边界文件上先跑直接目标 `ν=.1,tol=1e-6`，再跑 `ν=1,tol=1e-3` 引导，并把产品导出的同网格 `cell,x,y,u,v,p` 作为原目标初值；另以显式 Anderson 作独立算法控制。直接/引导/终档/Anderson 的完整评估为 `275/187/135/351`，三条目标路径均严格收敛到相同支路。引导目标相对直接目标的全场速度 RMS 为 `7.94e-8 Uref`、去规范压力 RMS 为 `3.41e-8 Uref²`；证书成本 `322/275=1.171`，direct+Anderson 成本 `626/275=2.276`。全部单元均参与比较，压力只移除面积均值；原规定速度迹跳、实际壁面法向速度和 Solver 门均未改。

```sh
outputs/cloud-laminar/guide-certificate manufactured 64 1e-6 \
  outputs/cloud-laminar/guide-manufactured-64
outputs/cloud-laminar/accuracy-probe-anderson manufactured 64 1e-6 \
  face-limited-linear anderson \
  outputs/cloud-laminar/manufactured-face-limited-linear-64-anderson-tol1e6
python3 artifacts/current/native-laminar-annulus-guide.py --run
```

汇总和每份原生 cells/faces/residuals 哈希见 `native-laminar-guide-crosscase.json`。该驱动的 `--resume` 只用于云端中断后复用已完成的直接/引导场；正式从零复现用 `--run`。

### 显式分支证书事务

`FlowBranchCertificate2D.hpp` 提供显式、固定物性 strict Newton 稳态 API，不挂入 `solveIncompressible2D` 默认。事务顺序固定为直接目标、`νguide=10νtarget` 且 `tolguide=sqrt(toltarget)` 的平坦引导、从引导场启动且恢复原黏度/原容差的目标；三阶段各自使用独立完整评估预算。调用方须显式给出两个正有限的无量纲 RMS 报告限值。比较覆盖全部最终单元，速度向量差按面积加权并除以 `Uref`，压力差先移除全域面积均值再除以 `Uref²`；只返回 `consistent` 分类，不返回“推荐”或“选中”候选。

`Unconverged/Stopped/Failed` 同时记录发生阶段，后续阶段不会覆盖 `directTarget` 或已完成的 `guide`。证书自己的阶段取消会在进入阶段前和完整求解更新后检查；原 `FlowControls2D::stopRequested` 仍参与每一阶段，任一回调抛出的异常继续向调用方传播。归档固定写三条候选记录；存在的候选各用现有 `FlowState2D` checkpoint 保存 `u/v/p/flux`，另存完整评估数和 `converged/stopped`，没有共享“当前场”槽。读回同时绑定最终网格、显式边界、物性、对流格式、引导参数、阶段预算和调用方限值，重算一致/分歧结果；截断、尾随数据、负/非有限差值、缺候选或阶段矛盾均拒绝。旧 v1–v4 flow checkpoint 读取器和非定常恢复格式没有修改。

CLI 只以显式开关暴露该事务；两个归一化 RMS 报告限值必须由调用者提供，不能沿用隐藏默认：

```sh
build/cartmesh2d_flow_cli \
  --mesh outputs/cloud-laminar/annulus-joint-4.solver.cm2d \
  --boundaries outputs/cloud-laminar/annulus-joint-4.boundaries \
  --output outputs/cloud-laminar/cli-branch-smoke/consistent \
  --case custom --nu .01 --speed .5 --tolerance 1e-6 --max-iterations 700 \
  --branch-certificate \
  --branch-velocity-rms-limit .01 --branch-pressure-rms-limit .01
```

事务发布 `PREFIX.branch.certificate` 和格式为 `cartmesh2d-flow-branch-certificate-summary-v1` 的 `PREFIX.json`。摘要固定为 `converged=false`、`selectedCandidate=null`，并列出三候选是否存在、收敛/停止状态、完整评估数、守恒和幅值诊断；普通 cells/faces/fields/residuals/VTK 均不生成，避免让脚本把任一候选误当普通接受场。退出码 0 表示三阶段完成且两目标路径在调用方报告限值内一致，3 表示三阶段完成但分歧，2 表示某阶段未完成；求解前输入错误仍返回 1。可选 `--branch-guide-viscosity-multiplier`、`--branch-guide-tolerance-exponent`、`--branch-guide-max-iterations` 和 `--branch-target-max-iterations` 均只影响显式事务。新鲜稳态之外的 profile、restart、initial guess/flux、非定常、adaptive、模板导出与非 Newton 方法会被拒绝。

192 格移动圆环真实运行的直接/引导/终档完整评估为 `68/34/37`，归一化速度/压力 RMS 差为 `2.9137e-6/6.6336e-7`；相同三候选在 `.01/.01` 下返回 0，在 `1e-8/1e-8` 下返回 3。把引导预算设为 1 返回 2，归档仍保留已完成直接场和未收敛引导场，终档为空。CLI 验证还实际终止 32×32 方腔证书进程，确认只保留 `running` 标记且不留下最终或 `.tmp` 归档。复现数据与哈希见 `native-laminar-branch-cli.json`。

桌面通过独立 `run-flow-branch-certificate` IPC 调用同一 CLI 事务。主进程只接受退出码与摘要分类一致的 `.json` 和 `.branch.certificate`，以独立 `*.flow.branch.*` 名称原子提交；复制失败恢复上一次证书，普通 `currentResult.flow` 与重启状态不参与该事务。renderer 只在新鲜 strict Newton 稳态启用入口，显示三候选和归一化场差，固定提示“未选择候选”；导出 ZIP 直接包含归档和摘要。Linux Electron 37 headless 实际在 728 格圆柱上先生成普通默认场，再运行证书：普通场 117 次严格收敛，证书 direct/guide/guided-target 为 `117/60/94` 次；两目标速度/压力 RMS 差 `1.1592e-7/5.1066e-8`。报告限值 `.01/.01` 时 DOM 显示一致，`1e-8/1e-8` 时显示分歧；两次均保留 `selectedCandidate=null`，普通场对象和逐格字段不变。复现证据见 `native-laminar-branch-desktop.json`；这不是实体显示器、macOS 或打包版验证。

真实反例复现：

```sh
g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror -Iinclude \
  artifacts/current/native-laminar-branch-certificate.cpp \
  build/libcartmesh2d_fv.a build/libcartmesh2d.a \
  -o outputs/cloud-laminar/branch-certificate
outputs/cloud-laminar/branch-certificate \
  outputs/cloud-laminar/cylinder-joint/far-20-fixed128-grid1.solver.cm2d \
  outputs/cloud-laminar/cylinder-joint/far-20-fixed128-grid1.external.boundaries \
  outputs/cloud-laminar/cylinder-joint/product-branch-certificate
ctest --test-dir build -R '^cartmesh2d_flow_(branch_certificate|boundary|checkpoint|initialization)$' \
  --output-on-failure
```

该 17,260 格运行以调用方显式 `1e-4/1e-4` 归一化 RMS 限值演示风险分类：直接目标、引导和引导终档为 `792/398/565` 次，目标间速度/去规范压力 RMS 差 `0.013484/0.070339`，因此归档为 `completed` 且 `consistent=false`。这两个报告限值不是新增 CFD 门；真实结论仍是“发现两个严格目标解，不能自动裁决”。18,229,452 字节最终证书写后用同一产品读取器恢复三候选；格式加强前的首次真实运行另名保留，仅作历史，不作为最终读回证据。命令、原始摘要、归档/网格/边界/源码哈希及单元失败、取消、损坏覆盖见 `native-laminar-branch-certificate.json`。

独立算法交叉检查保持目标方程、网格和全部最终严格门不变，只把显式稳态加速器切换为既有 Anderson。正常 4,716 格圆柱直接运行 276 次，与 Newton 场的速度/去规范压力面积 RMS 差为 `1.30e-8 Uref/1.76e-9 Uref²`；直通道 64 档和 Re=100 方腔 64 档分别用下列命令运行 123/1,145 次，并与 Newton 场保持约 `1e-6` 或更小的全场差异：

```sh
outputs/cloud-laminar/accuracy-probe-anderson channel 64 1e-8 face-limited-linear anderson outputs/cloud-laminar/channel-face-limited-linear-64-anderson
outputs/cloud-laminar/accuracy-probe-anderson cavity 64 1e-8 face-limited-linear anderson outputs/cloud-laminar/cavity-face-limited-linear-64-anderson
```

同一 17,260 格反例从平坦场启动 Anderson 则在 490 次后达到速度/压力范围比 `5.928/291.78`；从有界 `ν=1` 场启动也在 1,372 次后进入 `1.391/60.69` 的另一载荷抵消支路。平坦 Newton 与 Anderson 的全场速度 RMS 差为 `2.84e-2 Uref`，去规范压力 RMS 差为 `7.78e-2 Uref²`；两者合计 1,539 次，仍不能恢复或裁决有界支路。`native-laminar-branch-continuation.py` 对所有比较保留全单元面积权重和去规范压力差，并记录原始场哈希。该交叉检查只能生成“独立算法显著不一致”的分支风险证据，不能把任一路径静默指定为物理解。

面通量结构原型由 `artifacts/current/native-laminar-face-hodge.cpp/.py` 复现。C++ 直接通过 `makeFvMesh2D` 读取完整 Solver 多边形，装配局部面内积、制造解和真实几何；Python 仅处理这些原生矩阵的稀疏线性代数、场差与输出哈希。没有重建独立 Python 流体方程，也没有修改产品求解器。

对每个单元，N 的每一行为向外面积矢量 S_f，R 的每一行为 C_f−C_i，V 是真实面积。多边形几何满足 `NᵀR=V I`。原型取

```
P = I − N (NᵀN)⁻¹ Nᵀ
alpha = tr(R Rᵀ/V)/2
H_i = R Rᵀ/V + alpha P
```

于是 `H_i N=R`，常速度场的能量正好为 `V|U|²`。对任意非零面向量 z，`zᵀH_i z=|Rᵀz|²/V+alpha|Pz|²>0`：若两项都为零，则 `z=Nw` 且 `RᵀNw=Vw=0`，故 z=0。alpha 是几何迹/空间维数给出的原型尺度，不是可调接受门或案例系数；未扫系数。H 按真实 owner/neighbour 符号装配，B 为全部单元的精确有向面关联。散度 Bq 和压力项 Bᵀp 共用同一关联；闭边界投影满足 `H(q−q0)+Bᵀπ=0, Bq=0`，因而 `E(q0)−E(q)=0.5(q0−q)ᵀH(q0−q)`，压力在投影后不做净功。q 单位 m²/s，H 在二维无量纲，E 单位 m⁴/s²；取1s投影步时 π 单位 m²/s²。这是新面内积性质，不代表旧单元中心动量方程已稳定。

```
OPENBLAS_NUM_THREADS=1 python3 artifacts/current/native-laminar-face-hodge.py
```

该命令以 `-Werror` 构建原生驱动，运行8组装配及每组线性压力补丁、光滑压力制造解、全域闭边界投影，输出 `outputs/cloud-laminar/face-hodge/` 与同名证据JSON。`--resume` 只接续已完成且原生源码/二进制哈希一致的组；常规复现不加它。制造解 `p=sin(kx(x−xmin))sin(ky(y−ymin))` 在同一固定外域使用固定 kx/ky，原生 RHS 采用单元中心二阶积分，所有边界给解析 Dirichlet 压力，不能等同于无穿透绕流压力。投影另对全部边界施加零法向扰动，全部单元参与，面值和压力场均保留。两档圆环几何不同，只作为结构控制，不据此声称固定几何收敛。

复用原17260格异常/有界接受场还做投影不变性负对照：固定其全部原边界通量，仅校正内部面。原边界净通量的舍入残差在压力定规范单元完整计入散度，未删除；所以全域最大散度可能从约1e-11增至2e-10，不能报告为质量改善。两解的投影改变量均在原残差尺度，明确否定“只投影即可修复支路”。当前未知量包括独立面通量和单元压力，尚无无滑移切向应力、对流和动量推进；下一步需要联合实现并验证这些部分。原始矩阵、原生二进制、全部求解场、失败构建/输入记录和哈希随恢复包保存。参考 [Ham & Iaccarino 2004](https://web.stanford.edu/group/ctr/ResBriefs/2004/ham_iaccarino.pdf) 对压力功与偏斜一致性的联合讨论；其等权网格结果不作为本项目不等距 Cut-cell 的资格证据，本原型靠自身几何恒等式及实际制造解验证。

压力耦合机制实验由 `native-laminar-viscous-energy.cpp` 和 `native-laminar-branch-energy.cpp` 直接调用产品算子/读取产品共享面通量，`native-laminar-pressure-energy.py` 只做编排、原生矩阵特征分解和证据整理，不重建独立流体方程。黏性限制矩阵为 `M^(-1/2) K M^(-1/2)`，其中 K 是负的积分向外黏性通量，M 是真实单元面积；省略的扰动自由度设零，原物理单元没有删除。压力稳定块仅以 `rAU=1 s` 研究几何符号，不能冒充完整鞍点或时间稳定性分析。差分功使用两份同网格同目标原方程的完整场，单位是二维每单位密度的 `m⁴/s³`；压力常数通过原多边形闭合抵消，不删格或局部去均值。

`native-laminar-pressure-pair-ablation.py` 在 `outputs/cloud-laminar/pressure-pair-ablation/` 临时生成一个原生对象和头文件覆盖，只把压力面值改为中点平均/未知边界 owner 值，把稳态 Rhie–Chow 中的速度面插值及松弛缺陷同步改为中点平均；其余组件链接当前原生库。该消融仅研究两个 external 稳态算例，未授予物理时间或其他工况资格。正常产品源码和二进制不改写。它没有同时保证稳定块正性和不等距几何的线性一致性，实际已由局部压力恶化和线性补丁失败否证，不能当作新格式。复现前需构建当前 `cartmesh2d_flow_cli` 目标：

```sh
python artifacts/current/native-laminar-pressure-pair-ablation.py
OPENBLAS_NUM_THREADS=1 python artifacts/current/native-laminar-pressure-energy.py
```

前者实际重新求解两次；后者构建研究驱动、输出十个原生限制矩阵、三组完整差分功和线性压力补丁，并读取前者结果。输出位于两个对应 `outputs/cloud-laminar/` 子目录，摘要为 `artifacts/current/native-laminar-pressure-energy.json`。原始矩阵、每格功、两次消融的全部场/面/历史、源码与二进制哈希必须一并保留。当前证据只支持下一步联合设计压力力、质量通量稳定项与几何一致性，不允许静默选择支路。

`native-laminar-coupling-energy.cpp` 对最畸变面附近的有限图支撑，用同一原生黏性面通量、保守压力梯度、速度面插值和黏性对角 `rAU` 构造局部 `S=C-DH A_patch⁻¹ G`；以主元消元求速度响应，再报告对称部最小特征值和原算子直接二次型。它使用欧氏归一化，黏性值单位为 `m²/s`、压力块为 `s`；不能和前述面积归一化速率直接比较。正方格/剪切格复用现有 `flow_face_test.cpp` 夹具。可选 `minimum` 只研究最小压力分解，`midpoint` 保留偏斜修正而把内部面插值权重设为 1/2；两者仅作用于探针，不修改产品。Mac 复现方式如下，Linux 按已有原生库链接方式省略 Accelerate；`-Wno-unused-parameter` 对应复用测试夹具已有的未使用形参。

```sh
mkdir -p outputs/laminar-stability/coupling-energy
/usr/bin/clang++ -std=c++20 -O2 -Wall -Wextra -Werror -Wno-unused-parameter \
  -I include artifacts/current/native-laminar-coupling-energy.cpp \
  build/libcartmesh2d_fv.a build/libcartmesh2d.a -framework Accelerate \
  -o outputs/laminar-stability/coupling-energy/probe
outputs/laminar-stability/coupling-energy/probe square 5
outputs/laminar-stability/coupling-energy/probe shear 5
outputs/laminar-stability/coupling-energy/probe path/to/mesh.solver.cm2d 5
```

指定同一网格追加 `minimum` 或 `midpoint` 可复算研究对照；支撑圈数仅限制探针成本，不是产品质量参数。结果见 `native-laminar-coupling-energy.json`，本批输入、矩阵计算程序和输出保存于其列出的忽略提交恢复包。负二次型只能用于形成机制假设，不能冒充全局失稳或选择物理解。

`native-laminar-hybrid-stokes.cpp` 是可直接求解的最低阶单元/面混合 Stokes 研究原型。每格两个速度均值、每个共享面两个速度均值及每格一个压力均值；`B u` 使用原生向外面积向量，压力项为 `-Bᵀp`。黏性双线性式为 `ν|T|G:G + νΣ_F |F|/h_T |u_F-u_T-G(c_F-c_T)|²`，`G=Σ_F u_F⊗S_TF/|T|`，`h_T` 是多边形直径。它使用 Laplacian Stokes 形式，尚无非线性对流、时间推进或产品牵引出口。常体力加载利用散度保持重构的一阶矩 `∫_T Rv=Σ_F(c_F-c_T)(v_F·S_TF)`；这仅足以检查全局仿射压力，不等于已实现一般体力/对流所需的完整重构。设计参考 [Quiroz–Di Pietro HHO 方法](https://arxiv.org/html/2203.07180v3) 的面未知量、散度配对与载荷思路，不冒用该文完整算法及误差定理。

```sh
mkdir -p outputs/laminar-stability/hybrid-stokes
/usr/bin/clang++ -std=c++20 -O2 -Wall -Wextra -Werror -Wno-unused-parameter \
  -I include artifacts/current/native-laminar-hybrid-stokes.cpp \
  build/libcartmesh2d_fv.a build/libcartmesh2d.a -framework Accelerate \
  -o outputs/laminar-stability/hybrid-stokes/probe
outputs/laminar-stability/hybrid-stokes/probe outputs/laminar-stability/hybrid-stokes \
  > outputs/laminar-stability/hybrid-stokes/results.jsonl
```

该目录参数可省略；提供时保存每次的原始单元/面 CSV。全部输入由原生多边形生成，通道制造解边界与单元量使用精确二次均值，压力只去全局加性规范；面积、实际面与原 Solver 门均保留。切割 8/12 档明确输出质量拒绝，其余求解错误直接非零退出。稠密部分主元消元的 1500 未知量限制仅控制研究成本，不是产品规模门；无新增精度阈值。JSON 中保存原始数值、单位、未验范围、源码与小型实际场归档哈希。当前原型尚未计算原 Re=20 圆柱/圆环，也未证明通道精度足够，不能自动接入默认。

`native-laminar-hybrid-stokes-p1.cpp` 将研究空间提高到 P1 单元/面速度、P1 单元压力与 P2 势重构。弱梯度和散度都以真实面矩及单元矩积分；压力为同一散度的负转置。稳定项采用 `π_F¹(r-π_T¹r+u_T)-u_F`，沿用面长/单元直径缩放，故二次速度/一次压力保持一致性。一般体力仍是单元 L2 载荷，尚未实现完整 H(div) 重构；不要把本文献启发的 Stokes 原型称为完整压力鲁棒 Navier–Stokes。以下需已有原生库，输出前缀所在目录须存在：

```sh
mkdir -p outputs/laminar-stability/hybrid-stokes-p1
/usr/bin/clang++ -std=c++20 -O2 -Wall -Wextra -Werror -Wno-unused-parameter \
  -I include artifacts/current/native-laminar-hybrid-stokes-p1.cpp \
  build/libcartmesh2d_fv.a build/libcartmesh2d.a -framework Accelerate \
  -o outputs/laminar-stability/hybrid-stokes-p1/probe
outputs/laminar-stability/hybrid-stokes-p1/probe sheared 8 vortex 8 \
  outputs/laminar-stability/hybrid-stokes-p1/sheared8-vortex-q8
```

CLI 依次接受 `square|sheared|cut|split|tip`、分辨率、`couette|rotation|hydrostatic|poiseuille|vortex`、积分阶次（可省略，默认 8）、输出前缀（可省略）。实际 37 次控制清单、单位、原始数值与归档哈希在 `native-laminar-hybrid-stokes-p1.json`。多边形积分使用有向边与重心扇形的带符号 Duffy 积分；未改变真实几何。压力 RMS 是全 P1 场去全域规范后的积分误差；压力最大值仅在积分点采样，不能当作壁面连续峰值。CSV 记录规范基函数所需的真实中心、直径和二阶矩。稠密 2600 未知量上限只控制研究成本；大网格需要稀疏求解和消元，尚未宣称适用。

非恒定载荷研究入口 `native-laminar-hdiv-load.cpp/.py` 复用现有 `native-laminar-hybrid-stokes.cpp` 的几何夹具与稠密消元；后者仅新增 `CARTMESH_HYBRID_STOKES_NO_MAIN` 包含保护，原独立程序不变。局部扇形的每个三角形采用 `R(v)=a_j+b(x−c_T)`，`2b=D_T(v)`。外面法向矩和径向法向连续性作为约束，剩余环流最小化 `∫|R(v)−v_T|²`；不修焊几何，扇形非正或端点不能按构造舍入误差闭合则研究驱动显式拒绝。重构借鉴 [2203.07180v3](https://arxiv.org/abs/2203.07180v3) 的2.4/3.2节，实际实现是二维k=0局部约束问题，不宣称完整高阶/NS方法。6×6 Duffy Gauss积分用于所有载荷与精确单元均值；三次势梯度乘RT0为多项式精确积分。逐基函数核验 `∫grad(phi)·R(v)=sum(mean_F(phi) v_F·S)−mean_T(phi) B_T(v)`，共享面项装配抵消，因此同一固定边界下势载荷只改变离散压力。

复现：`python3 artifacts/current/native-laminar-hdiv-load.py`；已有原始场的确定性汇总用 `--summarize`。C++以GCC13.3、`-Werror`构建（复用测试夹具的已知未用参数用`-Wno-unused-parameter`），输出252个小单元/面CSV、126次结果及两项原始大网格局部重构检查。Python只编排和按全部单元面积比较，压力只减全域面积均值。参考尺度L=1m、U=1m/s、ν=1m²/s；势为 `U²*((x/L)^3+(x/L)*(y/L)^2)`，叠加强度0/1/10000只是制造解参数，不是接受门。真实曲壁网格检查的梯度载荷恒等式系数单位m³/s²（每单位速度基），法向积分余量单位m；该多项式在大外域的幅值较大，绝对余量不直接用作跨网格门限。小求解稠密系统限制1500未知量，原17260/6606网格只做局部重构，未运行全局NS。保留低阶剪切涡压力慢收敛与cut8/12原门拒绝记录；没有采用小幅值场作为物理解。

P1 原型新增 `square n noslip` 与 `sheared n noslip-sheared` 控制，几何/问题不匹配直接拒绝。令 `ξ=x-sy,η=y`，`s=0/0.7`，流函数 `ψ=64 g(ξ)g(η)`、`g(t)=t²(1-t)²`，速度为 `(∂yψ,-∂xψ)`，压力为 `sin(πξ)sin(πη)`，体力为 `-νΔu+∇p`；采用 1m、1m/s 参考单位和 `ν=1m²/s`。壁面上 g 与 g' 同时为零，故两阶面速度均直接设为解析零值。边界反力取本地刚度行与压力负转置的完整残差（含稳定项），按两个面矩重构 P1 牵引，并以去全域规范的压力比较解析物理应力。对这个静止无滑移、无散解析场，Laplacian 与对称应力牵引相同；不外推一般运动壁面。新增 `.walls.csv` 保存每面两分量反力系数；实际27份字段全部有限且已归档读回。复用上方编译命令后运行：

```sh
outputs/laminar-stability/hybrid-stokes-p1/probe square 8 noslip 8
outputs/laminar-stability/hybrid-stokes-p1/probe sheared 8 noslip-sheared 8
```

7 次静止壁面控制、3 次既有控制、三档误差及积分阶次对照见 `native-laminar-hybrid-stokes-p1-wall.json`。壁面牵引误差下降慢于内部速度；这是需要后续改进与验证的实际边界，不增加阈值、不升级默认。

联合高阶研究入口 `native-laminar-p1-stress.cpp/.py` 包含 P1 原型（用 `CARTMESH_HYBRID_STOKES_P1_NO_MAIN` 保护独立 main），使用完整对称弱梯度应变与原 P2 一致稳定项。三角形 Piola RT1 通过局部约束最小化保持每个原子面两阶法向矩、内部径向法向连续和同一 P1 弱散度；只去掉一条冗余径向均值约束，再逐点核对该面。正重心扇形限制是当前研究支持范围，不替代产品 Solver 门。每格同时消去六个速度系数与两条压力斜率，形成8×8局部鞍点块；全局只保留内部面四个速度系数和每格压力均值。原生输出稀疏项/右端，Python仅做SciPy稀疏LU及代数修正，原生再恢复全部局部场与边界反力。

云端已有原生库和 NumPy/SciPy，复现命令如下；Mac构建需用系统 `/usr/bin/clang++` 并增加 `-framework Accelerate`。驱动拒绝覆盖既有结果前缀，单元与面的CSV、矩阵、原始解及失败保留云端恢复包。单元列 `potential0/X/Y` 保存三次压力势的P1投影，用于叠加载荷后的去规范压力差比较。边界反力包含相容载荷在面试函数上的作用，不能仅用重构应力代替；二者均输出。

```sh
g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror -Wno-unused-parameter \
  -I include artifacts/current/native-laminar-p1-stress.cpp \
  build/libcartmesh2d_fv.a build/libcartmesh2d.a -o build/native-laminar-p1-stress
python3 artifacts/current/native-laminar-p1-stress.py sheared 8 vortex 0 lift symmetric fresh-sheared8
```

原生参数为 `assemble|recover mesh n problem lambda lift|cell symmetric|laplace order prefix`；文件网格时n忽略，cylinder控制要求显式Embedded/Domain边界标识。这个cylinder案例在整个外边界施加(1,0)、内壁零速度、ν=1且无对流，输出静态运动学压力；不等于原产品压力出口或Re20。CSV中 `traction_x/y` 与 `stress_x/y` 是积分力矩而非逐点应力；最大速度、压力范围及散度取积分采样点，尚非连续最大值。154011未知量的全局LU规模已显示成本问题，后续需解决可扩展预条件、兼容对流和原问题边界；不能仅凭小矩阵正性或有界幅值选取物理解。

`native-laminar-block-precondition.cpp` 只读取原生凝聚矩阵、右端和单元面积，不重建方程。右块三角预条件器复用 `SparseSystem2D::factorIC0/factorILU0` 近似速度块逆。`ic0` 仅在预条件器中取速度块对称部分；`ilu0` 保留全部非对称系数，`jacobi` 为对角控制。采用两侧非零结构并集，所有Krylov乘积和停止检查始终用原K；不补主元、不移对角、不静默切换算法。旧 `velocity_preconditioner_symmetry_correction_max` 字段保留为预条件输入的两侧差；`velocity_matrix_asymmetry_max` 报告实际K的两侧差，`velocity_preconditioner_input_asymmetry_max` 报告所选速度预条件输入的两侧差。未指定独立预条件矩阵时两者相同；`velocity_preconditioner_entry_change_max` 才记录对输入的实际系数变化，ILU0为零。

闭域固定最后一格压力为0，去均值质量矩阵为 `diag(V_i)-V_i V_j/Vtotal`，`gauge` 逆为 `diag(1/V_i)+11ᵀ/V_last`，`plain` 省略秩一项作负对照。自然出口用 `outlet`，保留Nc个压力及其绝对水平，只用对角质量逆，不减均值。各布局的压力近似逆均乘以 `-ν`；正黏度是显式输入，这是黏性Schur近似，不是Oseen精确逆或对雷诺数/网格的鲁棒性保证。仍只支持当前速度在前、单元均压在后、零保留压力块的凝聚格式；部分速度分量边界约束允许使速度未知量不再是4的整数倍。新接口为 `input_prefix cell_csv ic0|ilu0|jacobi gauge|plain|outlet output_prefix [relative_tolerance=1e-11] [restarts=50] [viscosity=1]`，历史ν=1调用保持兼容。

研究求解控制为 `||b-Kx||₂/||b||₂`，在原始参考单位下默认 `1e-11`；不附加混合量纲绝对残差地板，不作为CFD验收门。既有GMRES每个60方向子空间内以0.1为工作目标，外层至多50轮重新计算真实残差。成功才写 `.solution`；用尽预算/方向失败写 `.candidate`、非零退出，已有前缀拒绝覆盖。Linux块求解器无需原生库或第三方依赖；Mac在下面构建命令中使用系统clang++并增加Accelerate框架：

```sh
g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror -I include \
  artifacts/current/native-laminar-block-precondition.cpp -o build/native-laminar-block-precondition
build/native-laminar-block-precondition existing_prefix existing_prefix.cells.csv \
  ic0 gauge new_prefix 1e-11
```

然后用原 `native-laminar-p1-stress recover` 与相同几何/物理参数读取 `new_prefix.solution`，分别核对实际压力、速度、散度与能量；不能仅看GMRES返回成功。9个案例的原生字段、多项式场积分差、普通质量对照和一档更紧控制见 `native-laminar-block-precondition.json`。NumPy仅用于2175未知量以内的直接线性代数参照及字段比较，没有独立方程审计；该直接参考重复分解以作残差修正，计时不用于生产稀疏LU性能比较。方法背景可参见[标准Krylov预条件](https://www.netlib.org/templates/templates.pdf)与[HHO凝聚/多层研究](https://arxiv.org/abs/2009.13840)，当前实现不是后者的p多层算法或定理复现。

几何预筛使用 `native-laminar-topology-spectrum.cpp`，对完整二次基 `r²、x²-y²、2xy` 调用产品梯度与修正扩散几何，输出旋转不变的二次一致性误差、梯度条件数、邻格面积比和非正交修正比；它不读取接受流场。坏中档的八个对称壁面模体均在求解前出现高值，但全局 `.075` 修复网格取得正常场后最坏二次误差仍约 `187.23`，不低于原网格 `184.05`。因此该量可定位候选模体，不能直接作为通过/失败判据。运行摘要、全部原始 SHA256 和成本由下列只读后处理固化：

```sh
g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror -Iinclude \
  artifacts/current/native-laminar-external-boundary.cpp \
  -Lbuild -lcartmesh2d_fv -lcartmesh2d \
  -o outputs/cloud-laminar/native-laminar-external-boundary
g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror -Iinclude \
  artifacts/current/native-laminar-topology-spectrum.cpp \
  -Lbuild -lcartmesh2d_fv -lcartmesh2d \
  -o outputs/cloud-laminar/native-laminar-topology-spectrum
python3 artifacts/current/native-laminar-branch-continuation.py
```

本批不改产品方程、默认、质量门或接受状态语义。若把黏度延拓发展为产品候选，已有通道和方腔正常支路控制仍须扩展到强迫涡、固定/联合圆环、更多 Re/几何外流及失败/取消/旧工况，并证明双候选场保留、完整成本和确定性；不能把已知反例硬编码成阶段表。双路径不一致只能证明分支风险，不能静默把延拓支路宣布为物理解。

固定 128 段圆环 level 4--8 另复用 `native-laminar-pressure-probe.cpp`，对解析压力调用产品 `buildFlowGradientStencil2D`、`pressureFaceValues` 和 `conservativePressureGradient`。面积加权保守梯度 RMS 为 `0.04075→0.01661→0.007133→0.003365→0.001306 m/s²`；数值压力峰值仍不收敛，最差格梯度矩阵条件数却仅 `2.22–9.53`（粗档 `8.11`）。因此本批不修改压力重构权重或门槛：现有证据反对一般重构退化，却仍保留规定迹角点压力未资格。原始 CSV、当前默认圆柱的逐字节复算、筛选测试及哈希见 `native-laminar-curved-pressure-profile.json`。

`FlowControls2D::convection=Default` 与 `resolveSteadyConvection2D` 现在只在固定物性、strict、实际采用 Newton–Krylov 的稳态上下文选择 `face-limited-linear`；显式 `--steady-acceleration none`、adaptive/engineering、物理时间和材料/温度联算解析为 `upwind`。三个显式 `--convection` 值均不改写。CLI 摘要和桌面保存工况记录解析后的实际格式；非定常 checkpoint 也把 Default 序列化为历史 Upwind，避免旧状态被新稳态策略重解释。17,796 格圆柱省略方法和对流选项的实际新默认重放以 559 次评估收敛，场、面量和残差历史与此前显式 face-limited-linear 逐字节相同。命令、哈希、100/100 原生及 194/194 前端回归见 `native-laminar-default-convection.json`；本批未运行实际 App、macOS 或打包版。

上下文默认的跨工况重放使用 `artifacts/current/native-laminar-default-matrix.cpp`；它复用原生测试网格构造和产品求解器，只记录解析误差及完整求解预算，不实现独立方程。通道、manufactured、cavity 的 `n=64` / `1e-8` 运行分别与显式 face-limited-linear 场逐字节相同。曲壁路径直接对 `annulus-joint-7.solver.cm2d` 省略 `--convection` 与 `--steady-acceleration`，与显式结果的 cells/faces/residuals/summary 逐字节相同。真实 App 复验命令必须给 `--out` 绝对路径并带 `--flow=external --flow-case-check=true`；相对路径会在保存工况时按设计拒绝。本次 Linux Electron headless 路径完成两次 5,168 格求解、保存读回和复算确定性，详情及哈希见 `native-laminar-default-regression.json`。它不替代实体显示器、macOS 或打包版验证。

曲壁压力局部化诊断用 `artifacts/current/native-laminar-pressure-localization.cpp` 直接读取产品 CM2D 与接受的 cells CSV，不组装或求解独立流动方程。它按最终流体面积计算压力误差分位数和平方误差集中度，同时保留全域最大值；按解析圆壁的物理距离分带只用于定位，不是验收时删除近壁格。还对每个接受网格调用产品 `buildFlowGradientStencil2D`、`symmetricViscousCorrection` 和同一紧致非正交通量，令 `u=-y,v=x` 做刚体旋转的零对称应变补丁试验。该仿射迹只是算子线性一致性检查，不是允许穿透真实多边形的替代壁面条件。

```sh
g++ -std=c++20 -O2 -Iinclude artifacts/current/native-laminar-pressure-localization.cpp -Lbuild -lcartmesh2d_fv -lcartmesh2d -o outputs/cloud-laminar/native-laminar-pressure-localization
for level in 4 5 6 7; do
  outputs/cloud-laminar/native-laminar-pressure-localization joint-$level \
    outputs/cloud-laminar/annulus-joint-$level.solver.cm2d \
    outputs/cloud-laminar/annulus-joint-$level-face-limited-linear.cells.csv
done
```

最细网格补丁最大余量 `6.75e-13 m/s²`、真实转子最大法向速度 `5.35e-17 m/s`；但压力最大误差仍为 `1.806 m²/s²`，因此不能把补丁通过解释成局部压力正确。面积加权 L1、99% 分位和离壁 `0.025 m` 的最大误差随联合加密明显下降，而 50% 平方误差在最细档集中到 3 个格/`0.01645%` 总面积。最细 Upwind 也出现同一峰值与集中度，说明默认 face-limited-linear 不是原因。完整数组、原始 JSONL 与驱动/二进制哈希见 `native-laminar-pressure-localization.json`；下一步应固定同一多边形只细化空间网格，或取得保持该多边形无穿透语义的一致参考，不能从当前圆形解析解与局部峰值直接推出产品修复。

同一诊断器现同时从真实边界原子面的中心和面积向量重建相邻分段共端点，记录每个单元到最近非零壁速迹跳点的距离。固定多边形空间研究逐字节复用 `annulus-joint-5.xy`（每环 128 段），只改变 level 4–8；level 5 复用既有同轮廓接受场，其余四档重新生成、导出 annulus 边界并以 `custom` 求解：

```sh
for level in 4 6 7 8; do
  build/cartmesh2d_cli outputs/cloud-laminar/annulus-joint-5.xy \
    outputs/cloud-laminar/annulus-fixed128-l$level "$level" .05 .05 interior \
    outputs/cloud-laminar/annulus-fixed128-l$level-foam 0 0
  build/cartmesh2d_flow_cli \
    --mesh outputs/cloud-laminar/annulus-fixed128-l$level.solver.cm2d \
    --case annulus --speed .5 \
    --export-boundaries outputs/cloud-laminar/annulus-fixed128-l$level.boundaries
  build/cartmesh2d_flow_cli \
    --mesh outputs/cloud-laminar/annulus-fixed128-l$level.solver.cm2d \
    --case custom --boundary outputs/cloud-laminar/annulus-fixed128-l$level.boundaries \
    --nu .1 --speed .5 --convection face-limited-linear \
    --tolerance 1e-8 --max-iterations 1800 \
    --output outputs/cloud-laminar/annulus-fixed128-l$level-face-limited-linear
  outputs/cloud-laminar/native-laminar-pressure-localization "fixed128-l$level" \
    outputs/cloud-laminar/annulus-fixed128-l$level.solver.cm2d \
    outputs/cloud-laminar/annulus-fixed128-l$level-face-limited-linear.cells.csv
done
```

五档保持相同 `0.0245412 m/s` 最大迹跳和机器精度无穿透。全域峰值在最细档升到 `0.6147 m²/s²`，不能称点值收敛；但距每个迹跳点至少 `0.025 m` 的压力 RMS 在 level 5–8 为 `5.49e-3→1.58e-3→1.07e-3→5.90e-4 m²/s²`，对应最大值 `0.0448→0.0164→0.0123→0.0105 m²/s²`，最细仍保留 `96.74%` 流体面积。该距离带只定义可复现的观测范围，不改变求解或全域输出；所有近角点误差仍在 JSON 的全域 RMS/最大值中。固定折线角点上不连续切向 Dirichlet 迹的点值压力不得宣称合格；若产品后续显示适用范围，应同时输出这一几何/边界语义，而不能静默删格。完整网格、原始场、时间、质量值与哈希见 `native-laminar-fixed-polygon.json`。

原生产品现直接生成 `FlowWallTraceDiagnostics2D`。诊断从最终 `Face::centre` 与 `Face::areaVector` 重建两个端点，只以 `TolerancePolicy::constructionRoundoffScale` 配对共享顶点；这不是更大的 weld/修复容差。参与比较的是求解器实际采用的 `Wall/Lid` 无滑移面值，因此静止壁、`MovingWall`、`SmoothMovingWall` 和方腔预设共享同一语义。每个共顶点取相邻壁面笛卡尔速度的最大两两差，超过 `TolerancePolicy::scale(referenceSpeed)` 才计作非零迹跳；结果只附加到 `FlowResult2D`/CLI，不进入矩阵、边界值或验收门。

```sh
build/cartmesh2d_flow_cli \
  --mesh outputs/cloud-laminar/annulus-joint-5.solver.cm2d \
  --case custom \
  --boundary outputs/cloud-laminar/annulus-joint-5-face-limited-linear.boundaries \
  --nu .1 --speed .5 --convection face-limited-linear \
  --tolerance 1e-8 --max-iterations 1800 \
  --output outputs/cloud-laminar/wall-trace-product/annulus-fixed128-l5 --profile
```

该复算报告 `wallTraceWallFaces=352`、`wallTraceAdjacentVertices=352`、`wallTraceDiscontinuousVertices=128`、最大跳量 `0.02454122852291445 m/s`；cells/faces/fields/residuals 与旧接受前缀逐字节一致。原生边界测试另验证顶盖方腔恰有两个跳点、静止通道为零；完整 CTest 100/100 和前端 195/195 通过。桌面校验要求新诊断字段成组且数值自洽，但仍接受没有这些字段的旧摘要。当前二进制也实际走通 Linux Electron 37 headless renderer→IPC→CLI：728 格 Solver PASS 圆柱以 117 次完整评估收敛，原生摘要报告 40 个静止壁共顶点、零迹跳，真实结果 DOM 含诊断标题、计数和适用边界文字。该路径不是实体显示器、macOS 或打包版。复现命令、哈希和平台边界见 `native-laminar-wall-trace.json`。

独立曲壁喷管回归暴露了上述壁面限定诊断的通用缺口：`examples/complex/nozzle_profile.xy` 的均匀速度入口与上下无滑移壁在 `(-3, ±0.98) m` 共享顶点，壁—壁比较为零，但完整规定速度的边界迹实际有两个 `1 m/s` 跳变。产品因此保留兼容的 `wallTrace`，并新增 `FlowVelocityTraceDiagnostics2D`：选择最终边界中 `fixedU && fixedV` 的面，覆盖无滑移/移动壁及速度入口/出口；压力边界、远场和只固定一个分量的滑移/对称面不参与。端点仍只在 construction-roundoff 尺度配对，速度容差为 `TolerancePolicy::scale(max(referenceSpeed, prescribedSpeed))`，本组为 `1.01e-10 m/s`，只排除浮点构造舍入，不是物理精度门。

同一原生驱动把喷管 level 5/6/7 加密到 636/2380/9336 格，省略方法和格式后均解析为 Newton + face-limited-linear，并以 75/133/426 次完整评估严格收敛；全局相对不平衡为 `1.43e-13/4.35e-14/6.33e-13`。速度最大值 `2.8564→2.8788→2.8844 m/s`，壁面合力 X 分量 `15.054→15.643→15.981 m³/s²`；但压力最大值 `9.71→13.87→17.88 m²/s²`，其单元中心到最近入口—壁角点距离 `0.108→0.0428→0.0263 m`，说明点值峰值正在靠近已报告的不连续规定迹，不能伪称一般近壁压力已收敛。独立 4904 格 Re=20 圆柱报告 208 个完整规定速度面、207 个共顶点和零迹跳，295 次评估后 `Cd=2.152836638`；cells/faces/fields/residuals 与此前接受场逐字节一致。新增诊断不改变方程、边界、网格、质量门或场。完整 Linux 原生 100/100、前端 196/196 通过；同步当前 CLI 后，Electron 37 headless 的真实 renderer→IPC→CLI 路径生成 728 格 Solver PASS 圆柱、117 次评估收敛，并在 DOM 显示新诊断计数与“仍需压力网格加密验证”的限定。驱动、完整成本、哈希与原始记录见 `native-laminar-prescribed-trace.json`；它不是独立方程审计，也不是实体显示器、macOS 或打包版验证。

后续分类不尝试在折线角点伪造连续速度。若两个实际面不共线，分别严格切向的速度空间是两条不同直线，其交集只有零向量；所以非零壁速若保持两侧真实面无穿透，就不能在该顶点具有同一个笛卡尔迹。产品新增的原因分类对每对规定速度计算 `abs(|Ua|-|Ub|)`：该量超过速度舍入容差为“速度大小不连续”，否则笛卡尔向量跳仍非零时为“等速方向转折”。它对坐标旋转不变，多面顶点的两类计数允许重叠。`wallTraceMaximumNormalVelocity` 另报告每个实际多边形壁面上的 `|Ub·n|`，单位 m/s，不投影到解析圆，也不改变边界值。

```sh
build/cartmesh2d_flow_cli --mesh outputs/cloud-laminar/annulus-fixed128-l4.solver.cm2d \
  --case custom --boundary outputs/cloud-laminar/annulus-fixed128-l4.boundaries \
  --nu .1 --speed .5 --tolerance 1e-8 --max-iterations 1500 \
  --convection face-limited-linear --steady-acceleration newton-krylov \
  --output outputs/cloud-laminar/trace-cause/annulus-fixed128-l4
build/cartmesh2d_flow_cli --mesh outputs/cloud-laminar/trace-regression/duct-l5.solver.cm2d \
  --case duct --nu .1 --speed 1 --tolerance 1e-8 --max-iterations 1500 \
  --convection face-limited-linear --steady-acceleration newton-krylov \
  --output outputs/cloud-laminar/trace-cause/duct-l5
build/cartmesh2d_flow_cli --mesh outputs/cloud-laminar/cylinder-joint/level-1.solver.cm2d \
  --case external --nu .1 --speed 1 --tolerance 1e-8 --max-iterations 1500 \
  --convection face-limited-linear --steady-acceleration newton-krylov \
  --output outputs/cloud-laminar/trace-cause/cylinder-level1
```

圆环分别得到大小冲突/等速转折 `0/128`，最大模长差 `1.67e-16 m/s`、实际面法向速度 `4.21e-17 m/s`；喷管为 `2/0` 和 `1 m/s`；静止圆柱为 `0/0`。三次原生复算都严格收敛，守恒与积分载荷保留，四类主输出和旧接受解逐字节一致。原生测试还构造全周等速移动方形：4 个非共线角点全归类为等速方向转折，法向速度在报告容差内。Linux Electron headless 的真实结果面板显示新增计数和无穿透量。完整字段、哈希、测试与平台范围见 `native-laminar-trace-causes.json`。

该诊断只调用产品压力/黏性重构算子，对实际网格和解析场作一致性对照，不重建独立离散方程。黏性项以 `m/s²`、压力梯度以 `m/s²` 报告。连续圆形解析速度在多边形面上的反事实对照有非零法向分量，仅用于识别边界表示敏感性，禁止作为求解边界绕过无穿透检查。完整诊断见 `native-laminar-pressure-diagnosis.json`；该诊断批次本身未产生通用物理精度门或产品算法改动，后续上下文默认升级依据是跨算例精度与完整成本证据。

静止壁面恢复由 `artifacts/current/native-laminar-wall-recovery.cpp/.py` 复现；`native-laminar-p1-stress.cpp` 仅新增 `CARTMESH_P1_STRESS_NO_MAIN` 包含保护，原 CLI 和算子不变。新 C++ 直接复用 Element 的完整局部矩阵、体力和恢复状态，以边界行余量的常量/一次矩 `r0,r1` 定义 `t_h(s)=(r0+12*s*r1)/|F|`，`s∈[-1/2,1/2]`。所有面体力必须扣除；省略它们的负对照保留。解析静止无滑移且不可压时，法向黏性应力为零，因此研究压力迹取 `p_t=-n·t_h`。另独立保留 `p_sigma=p_t+c*n·G_h(u)·n`，完整对称/Laplacian 形式在 `ν=1` 下的 c 分别为2/1；该对照在本批非零流动中变差，不按案例选择较好定义。压力只去一个全域面积均值，牵引同步加该常数乘法向；墙面 RMS 按全部真实面长归一化，单位 `m²/s²`，Gauss点最大值不是连续最大值。

```sh
/usr/bin/clang++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror -Wno-unused-parameter -I include artifacts/current/native-laminar-p1-stress.cpp build/libcartmesh2d_fv.a build/libcartmesh2d.a -framework Accelerate -o build/native-laminar-p1-stress
/usr/bin/clang++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror -Wno-unused-parameter -I include artifacts/current/native-laminar-wall-recovery.cpp build/libcartmesh2d_fv.a build/libcartmesh2d.a -framework Accelerate -o build/native-laminar-wall-recovery
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 artifacts/current/native-laminar-wall-recovery.py outputs/laminar-stability/wall-recovery-reproduction
```

需要当前原生库与 NumPy；Linux 使用相应编译器/既有库链接方式并省略 Accelerate。Python仅编排、读原生矩阵并做稠密线性代数，限制2300未知量，不另解一套流体方程；`--case ring-1-8-lift` 可只复现一个控制。已有案例前缀会拒绝覆盖，部分失败原样留存。C++ 的 `write-ring radial angular file` 建立半径0.5/1米、共享顶点身份一致的体贴合多边形圆环并通过原 `makeFvMesh2D`，它不是 Cartesian 网格生成器，也没有修补原曲壁几何。静水 `u=0,p=x+2y+x³+xy²` 在每个实际多边形上均是精确问题；面细化对照不是原圆环旋转流的资格。直接恢复接口为 `native-laminar-wall-recovery mesh n problem lambda lift|cell symmetric|laplace order existing_prefix`，读取同一原生 `.solution` 后输出 `.wall-recovery.csv`；只允许已声明的静止制造解，不能原样用作不同黏度、对流或开放出口的载荷后端。40项完整运行、追加恢复、源文件/二进制、矩阵、所有实际场和不利对照逐项压缩读回，索引为 `native-laminar-wall-recovery.json`，重现单控制与原结果完全一致。随后合并云端后，Stokes入口及壁面读取器都将这些解析静止边界的系数设为严格零；原40项的投影舍入记录和源码哈希作为历史证据保留，合并后剪切4档壁面系数最大变化约 `1.58e-14`，详见索引中的独立读回。

相容惯性入口 `native-laminar-p1-transport.cpp` 复用 `native-laminar-p1-oseen.cpp`、P1对称黏性及RT1相容载荷。对流采用RT1重构的对流速度、被输运速度和试函数；扇形内部及原子面使用同一法向通量和上风跳项，在线性法向通量的真实零点分段Gauss积分。原单元试函数版本保留为仿射不一致负对照。压力始终是运动学静压；自然出口 `(2ν sym G-pI)n=0` 不等于旧产品的规定压力/速度零法向梯度出口。边界反力分开记录压力、黏性加稳定项和对流重构项，后两项不得混称纯黏性；入口完整反力补回对流通量后才对应物理牵引。

```sh
g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror -I include artifacts/current/native-laminar-p1-transport.cpp build/libcartmesh2d_fv.a build/libcartmesh2d.a -o build/native-laminar-p1-transport
python3 artifacts/current/native-laminar-p1-oseen.py square 8 noslip 0.1 ns closed unique-label
```

上例使用已有Linux原生库及NumPy/SciPy；Mac采用系统Clang和Accelerate链接。圆柱用真实 `.solver.cm2d` 路径替代square、n=0、problem=cylinder、boundary=open。输出前缀拒绝覆盖；`CARTMESH_P1_INITIAL_STATE`可显式接续保存的研究场，meshKey绑定真实几何，不是产品checkpoint。每步保存完整解、原生恢复和原始非线性方程复核；`1e-9`研究控制同时约束系数变化及未消元残量（U=L=1，静压尺度U²），不是新增物理门或选解规则。Python只做编排和原生矩阵稀疏线代；失败、迭代上限及完整成本保留。P2势速度和RT1输运速度各报其误差，P1压力极值覆盖全单元顶点，速度峰值仍只是采样值。源码哈希、两档曲壁NS、仿射/二次负对照、积分对照、Mac4项短验证及云端恢复包持久保存失败均记录在 `native-laminar-p1-convection.json`。原出口等价性、回流、充分空间/外域收敛、壁面反力精度和可规模化联合求解仍未完成。

`native-laminar-oseen-block.py` 提供相同原生装配、恢复和完整非线性检查的块求解入口。默认ILU0路径只用Python标准库；NumPy仅在显式 `--backend dense-reference` 的2175未知量以内参考中使用，不另写流体方程。运行前构建上述 `native-laminar-p1-transport`、`native-laminar-block-precondition`，并构建共用几何/场比较入口：

```sh
g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror -I include artifacts/current/native-laminar-state-compare.cpp build/libcartmesh2d_fv.a build/libcartmesh2d.a -o build/native-laminar-state-compare
python3 artifacts/current/native-laminar-oseen-block.py sheared 8 noslip-sheared 0.1 ns closed outputs/laminar-stability/oseen-block-reproduction
```

Mac用 `/usr/bin/clang++` 并加 `-framework Accelerate`。可显式传 `--transport/--block/--fields` 二进制路径，或 `--backend ic0|jacobi|dense-reference` 对照。研究默认线性相对控制为 `1e-13`，非线性沿用 `1e-9` 完整残差及系数变化双条件；前者从上一批物理散度敏感性选择，用于隔离代数误差，未改变产品精度门。输出目录必须不存在；各步命令输出、线性失败候选、完整场、原始方程复核及成本保留。`--iterations/--restarts` 分别限制非线性/线性工作量，耗尽返回2；只有最终双条件满足才记录 `completed=true/finalState`。线性失败不恢复或采用候选。`lastCompletedIterate` 仅指最后完成原生检查的研究迭代，可能未收敛；`--initial-state` 的meshKey只绑定几何，不提供产品级物性/边界/时钟检查点语义。

`native-laminar-state-compare geometry mesh n prefix` 从相同原生几何导出单元面积；`compare mesh n first.state second.state` 使用原生状态读取及多项式积分给出P1速度/压力RMS、面速度端点最大差。闭域压力同时报告绝对差及一次全域规范校正后的差；自然出口必须看绝对差。它不计算第二套PDE，也不把P1单元场与P2势速度或RT1输运速度混为同一指标。4个冻结矩阵、6个完整NS控制、4个失败对照和复现命令见 `native-laminar-oseen-block.json`；归档保留每个迭代和输入，8档均匀出口没有障碍物，不能用 `cylinder` 问题名宣称圆柱资格。后续云端已完成4716格原Re20的原生NS及两个实际Oseen矩阵，见 `native-laminar-product-outlet-block.json` 和 `native-laminar-oseen-pressure-scale.json`；完整冷启动总成本、较大规模与内存对照仍不齐全。

### 显式开放牵引与伪牵引研究

`native-laminar-open-boundary.cpp` 复用同一 `Transport` 和凝聚/恢复驱动，通过独立边界策略选择已知速度自由度，不复制体积方程。压力为运动学静压，`G_ij=∂u_i/∂x_j`、`sigma=ν(G+Gᵀ)-pI`。三种显式模式如下：

| 模式 | 所解的连续边界条件与数据 |
| --- | --- |
| `traction` | `sigma n=t_D`；制造解提供完整精确牵引，cylinder控制提供零牵引。 |
| `pseudo-traction` | `(νG-pI)n=-p_D n`；在对称应力形式中把 `νGᵀn` 作为未知速度的隐式边界项，给定压力载荷 `-p_D n`。 |
| `normal-stress` | `sigma n=-p_D n`；不施加切向牵引，是不同的物理边界条件。 |

伪牵引不等同于同时独立强制压力值与全部速度法向导数；它对平行充分发展流相容，但不能据此宣称与旧单元中心出口模板等价。转置梯度项的连续推导和梯度约定见 [FEniCSx边界牵引说明](https://jsdokken.com/dolfinx-tutorial/chapter2/navierstokes.html)。在单元内仍使用原完整对称黏性、守恒压力耦合、RT1载荷/输运及原子面；只修改保留面测试行，所以内部块、内部到面块和内部右端均不变，既有局部凝聚仍适用。已有默认 `closed/open` 调用继续使用原策略，三种新模式须显式选取。

当前研究选择器识别轴对齐右出口；制造解其余外边界给定完整解析速度，cylinder保留水平对称边。没有增加产品的通用命名patch API，也没有迁移旧工况。`forceWork` 是体力功，`prescribedTractionWork/outletGradientWork` 单列给定牵引与隐式梯度项功，`viscousEnergy` 保留体内黏性和稳定项；能量平衡包含上述所有项。CSV的 `appliedTractionX/Xs/Y/Ys` 是四个积分面矩，不能当作逐点应力；边界总牵引还包含原动量反力与必要的对流通量。开放制造解的 `pressureRms` 使用绝对压力，闭域仍只消除一次全域规范。

```sh
g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror -Wno-unused-parameter -I include artifacts/current/native-laminar-open-boundary.cpp build/libcartmesh2d_fv.a build/libcartmesh2d.a -o build/native-laminar-open-boundary
mkdir -p outputs/laminar-stability
build/native-laminar-open-boundary write-warped 8 outputs/laminar-stability/open-warped8.solver.cm2d
python3 artifacts/current/native-laminar-oseen-block.py outputs/laminar-stability/open-warped8.solver.cm2d 0 poiseuille 0.1 ns pseudo-traction outputs/laminar-stability/open-pseudo8 --transport build/native-laminar-open-boundary
```

Mac使用 `/usr/bin/clang++` 并加 `-framework Accelerate`；还需前述块求解器和原生字段工具。临时网格/求解输出拒绝重复目录；更换路径复现。`write-warped` 只创建有界研究网格：固定单位方域，内部共享节点按 `X=x+.09sin(2πx)sin(2πy)`、`Y=y+.06sin(πx)sin(2πy)` 变换，边界精确保留，通过原 `fromPolygons/makeFvMesh2D`；它不是对原XY几何的平滑或Cartesian Cut-cell替代。

控制为 `u=(4y(1-y),0), p=-8x, f=(8ν-8,0)`，因此ν=1 Stokes无体力，ν=.1 NS明确有 `f=(-7.2,0)`；Couette为 `u=(y,0),p=0,f=0`，静水为 `u=0,p=x+2y,f=(1,2)`。解析数据补全了Couette、旋转和通道的速度梯度，用于真实边界应力；原速度/体力/压力值保持。不同边界条件的解都可能严格收敛，`normal-stress` 对原通道剖面的差值属于模型对照，不叫作离散失败。规则格P2势速度和压力的特殊再现不能扩大到RT1输运速度或非正交格；非正交三档的三类场误差、一次更紧迭代对照、原默认状态逐字节读回与完整成本见 `native-laminar-open-boundary.json`。收紧控制沿用U=L=1、压力U²归一化，仅用于分离迭代误差；未新增物理验收阈值。原曲壁NS、回流、旧出口离散兼容性与产品接入仍需继续验证。

云端 `pressure` 研究模式现与 `pseudo-traction` 共用模型2的转置梯度算子，前者仅允许圆柱或 `outlet-poiseuille` 开口控制；原 `open` 总应力出口保持。`outlet-poiseuille` 为 `u=(4y(1-y),0), p=8(1-x), f=(8ν-8,0)`，解析梯度已补全；云端ν=1控制无体力。使用同一原生后端复现零参考控制：

```sh
python3 artifacts/current/native-laminar-oseen-block.py square 8 outlet-poiseuille 1 ns pressure outputs/pressure-outlet-control --transport build/native-laminar-p1-transport --block build/native-laminar-block-precondition --fields build/native-laminar-state-compare
```

旧产品的 `viscousFaceGradient` 在自由出口先去掉法向梯度分量，再由 `symmetricViscousCorrection` 组装转置项；新研究模式使用完整弱梯度。`pressure` 名称保留云端复现接口，不代表分别强制全部速度法向导数与压力迹，亦不代表复制旧单元中心模板。模型2的出口矩阵在基础Oseen对象中只组装一次，显式入口只追加规定牵引载荷；能量与全局力预算使用同一已记录边界项。

`native-laminar-product-wall-profile.cpp` 在同一固定多边形的两份研究场之间逐段覆盖全部原壁面，原 `wallPressureDifferenceRms` 和 `pA0/pB0` 仍是去全域均值后的历史量；新增 `wallPressureAbsoluteDifferenceRms`、绝对压力CSV列及两种端点最大差。固定压力出口必须看绝对量，不能只看去均值量。P1差在每个重叠区间内是仿射函数，故端点最大值是该区间精确极值；`MaxAtQuadrature` 旧字段仍保留原采样含义。调用为 `wall-profile meshA nA stateA meshB nB stateB new-output.csv`，覆盖不一致或已有输出明确拒绝；不同真实曲线/折线不能用此固定几何读取器声称网格差。历史云端数据与合并后的读取语义分开保存在 `native-laminar-product-outlet.json`，Mac小场对照在 `native-laminar-open-boundary.json` 的 `mergedCloudIntegration`。


### 相容方程的解析 Newton 与残差回溯

`native-laminar-newton.cpp` 包含现有开放边界入口，通过模板微分同一保守 RT1 体积、径向内面、原始共享面及出口对流。令 `F(u)=K_Picard(u)u-f`，增加 `D_beta C(u)[δu]u`，求解 `J(u_k)u_candidate=f+D_beta C(u_k)[u_k]u_k`，等价于 `Jδu=-F`。压力约束、弱边界及静态凝聚自由度保持；内部8×8消元使用真实J重新计算。上风积分仍按当前P1通量的真实零点分段，零点处沿原符号分支取导数；没有宣称全局光滑。

`check` 始终调用原 `Transport/OpenTransport`。`recover` 只输出 `linearizedNewton=true, physicalDiagnostics=false` 的候选和线性化反力列，不能从候选恢复阶段授予真实能量或力资格。`seed` 保持相同网格身份、写入规定面速度及闭域压力规范，明确标记 `physicalSolution=false`；`blend` 只在相同网格上混合两份完整原生状态并检查有限值。`tangent` 用确定性完整状态/方向对原生原方程做中心差分，非独立PDE。

重新编译上述原生块求解器和状态读取器，再构建本入口。macOS使用系统Clang；Linux使用系统C++20编译器并去掉Accelerate参数：

```sh
/usr/bin/clang++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror -Wno-unused-parameter -I include artifacts/current/native-laminar-newton.cpp build/libcartmesh2d_fv.a build/libcartmesh2d.a -framework Accelerate -o build/native-laminar-newton
build/native-laminar-newton tangent sheared 3 noslip-sheared .1 closed 6 1e-6
python3 artifacts/current/native-laminar-newton.py sheared 8 noslip-sheared .01 ns closed outputs/laminar-stability/newton/reproduction --iterations 40
python3 artifacts/current/native-laminar-oseen-block.py sheared 8 noslip-sheared .01 ns closed outputs/laminar-stability/newton/picard-reference --transport build/native-laminar-open-boundary --iterations 40
```

输出目录必须是新目录。Python只编排原生状态、方程与线性代数：每一步线性候选用原方程重算所有自由未凝聚方程的L2范数，评价函数为其平方的一半；`α=1,1/2,...`，Armijo系数 `1e-4`、最多15个试探。沿用 `U=L=1、p/U²` 的原多项式方程归一化，规定边界的反力不进评价函数；这个范数不是任意有量纲问题的通用物理门。接受试探仍不等于完成：原最大内部/自由面残差、弱散度及实际场系数变化都要通过原 `1e-9` 控制，线性相对控制为 `1e-13`。已有的Stokes一步线性语义保留。达到所有严格门时允许在舍入底部接受，不靠极小步长掩盖非零原方程残差。标准方法背景见[Newton线搜索](https://petsc.org/release/manualpages/SNES/SNESNEWTONLS/)；实现没有新增PETSc依赖。

默认研究选项 `--preconditioner picard` 输出同一状态的 `.picard.entries`；块求解器通过末尾可选命名参数 `--velocity-preconditioner prefix`，仅用该原生矩阵的速度块构造ILU0。实际J、右端、Krylov乘积、压力耦合及线性停止检查保持，压力块仍为黏性质量近似，原ILU0主元门不变。`--preconditioner newton` 保留直接用J构造ILU0的对照，`--backend dense-reference` 限2175未知量，只作线性代数参考，不能用于生产稀疏LU性能比较。`--no-line-search` 是负对照，仍须最终原方程门，不能视作稳健默认。

云端局部压力尺度保留为显式研究选项：块求解器用 `outlet-oseen` 模式，并在黏度位置参数后传正的 `transport_speed`；其他压力模式速度尺度须为零，旧调用保持。Newton驱动可附加 `--pressure-preconditioner local-oseen --transport-speed 1`，仅适用于开放边界的原生ILU0后端；默认仍为 `viscous`。所用压力逆尺度为 `ν+U_ref√A_i`，二维单位为 `m²/s`，未修改实际J、右端或边界。它能与独立Picard速度预条件矩阵同时使用，仍不是一般惯性Schur逆。真实4716格证据只是两个冻结Oseen矩阵的20%–22%乘积减少；合并检查另有一组完整4档制造解出口的401→350乘积，两者不能混为大网格完整Newton成本。

`lastAcceptedIterate` 只是研究迭代；线性失败不恢复候选，回溯失败不替换它，预算耗尽不设置 `finalState`。所有试探、真实方程输出、命令、二进制/源码哈希和失败原因留在 `result.json` 与同目录原始文件中。这里未重新验证产品取消/检查点语义。`native-laminar-newton.json` 记录四组完整同场对照及更粗剪切格的失败：低黏度剪切8档9步Newton对36步Picard，但规则8档ν=.1的Newton更慢；剪切4档ν=.01即用直接线性解仍回溯停滞。上述比较不等于新空间精度、曲壁压力或大网格默认资格。

### 相同稳态方程的原生伪时间推进

`native-laminar-pseudo-time.cpp` 复用解析Newton，给两个单元P1速度分量加入原生精确 `3×3` 质量矩阵，面速度与压力继续是代数未知量。绝对状态矩阵加 `M/Δτ`、右端加 `Mu_old/Δτ`，内部8×8重新凝聚，等价于 `[J(u)+M/Δτ]δu=-F(u)`；Picard速度预条件矩阵也含相同质量项。Stokes只在单元内令输运速度为零，同时把当前完整状态按引用传给质量项，避免逐单元复制全局向量。没有独立Python方程或PETSc依赖。

```sh
/usr/bin/clang++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror -Wno-unused-parameter -I include artifacts/current/native-laminar-pseudo-time.cpp build/libcartmesh2d_fv.a build/libcartmesh2d.a -framework Accelerate -o build/native-laminar-pseudo-time
python3 artifacts/current/native-laminar-newton.py sheared 4 noslip-sheared .01 ns closed outputs/laminar-stability/pseudo-time/reproduction --pseudo-step .1 --iterations 60
python3 artifacts/current/native-laminar-newton.py square 4 outlet-poiseuille .1 ns pressure outputs/laminar-stability/pseudo-time/pressure-schur --pseudo-step .1 --pressure-preconditioner schur-diagonal --iterations 40
```

先按前文构建原生块求解器和状态读取器。Linux使用系统C++20编译器并去掉Accelerate参数。`--pseudo-step` 默认0，继续原Newton回溯；正值才选择伪时间入口。原生 `assemble/recover` 额外接收有限非负的 `inverse_dt`，0精确退回原Newton矩阵。`check/seed/blend/tangent` 始终走原算子，完成门从不采用带伪质量的残差。它是稳态全局化方法，不是已验证的物理时间积分器。

这里使用单位增长因子的SER：`Δτ_next=min(Δτ_max,Δτ*||F_old||₂/||F_new||₂)`，最大步长显式默认 `1e6`；参考尺度仍为 `U=L=1、p/U²`。可以接受原稳态残差暂时增加的有限伪时间迭代，但每步必须先成功解线性系统；最终NS和Stokes均要求原方程/弱散度最大量及完整场变化通过 `1e-9`，线性相对控制仍为 `1e-13`。零范数或已完成状态避免比值奇点。一次线性Stokes的场变化例外只保留在不含伪时间的原路径。没有参数扫描、自动切换或伪时间重试控制器；线性/总迭代预算失败仍保留原场并明确未完成。方法背景见[伪时间推进手册](https://petsc.org/release/manual/ts/)和[SER实际步长实现](https://petsc.org/main/src/ts/impls/pseudo/posindep.c.html)，不据通用方法理论直接授予本DAE全局收敛保证。

云端新增 `--pressure-preconditioner schur-diagonal`，仅适用于开放边界原生ILU0后端，速度尺度参数须为0。块程序对应 `outlet-schur-diag`：在已凝聚的真实压力耦合D/G上，以所选速度预条件矩阵的正对角构造 `P=D diag(A)⁻¹ G`，施加 `-IC0(P)⁻¹`；若选择Picard速度矩阵，A也来自同一Picard矩阵，并包含伪质量。没有压力正则化、主元替换或静默兜底。原符号约定下D=Gᵀ，新增结构检查沿用原速度检查的 `1e-12 max(1,max|Pij|)` 系数容差，只用于确保IC0的对称结构；原真实线性残差和物理门不变。显式选项不改变原 `viscous/local-oseen` 模式。

`native-laminar-pseudo-time.json` 保存原失败状态恢复、`.01/.1/1` 初始步长同场对照、NS/Stokes/出口完整解、五步Stokes错误完成的最小失败例及云端合并联用。原始 `result.json` 中一次早期Stokes五步 `completed=true` 已明确标为被修复版本的错误判定，不计入当前成功。该粗剪切格的空间误差仍显著；同场、严格残差和较少矩阵乘均不替代曲壁压力或默认产品资格。所有原场/失败/源码版本保留在逐项校验的压缩包中，路径和哈希见该索引。

## 完整笛卡尔背景网格

`--background-grid adaptive|uniform` 在几何诊断、Quadtree 细化及 2:1 平衡后直接导出完整叶子，不做 Cut-cell、Solver 修复或 OpenFOAM 输出。均匀模式使最低层级等于最高层级；自适应复用尺寸场和盒加密。

```sh
build/cartmesh2d_cli examples/acceptance/circle.xy outputs/background-grid/circle 7 0.5 0.1 exterior - 3 --background-grid adaptive
```

`cartmesh2d-background-v1` JSON 保存域、原边界、完整单元范围/层级/整数格坐标及分类（0 外部、1 内部、2 相交），`solver_ready=false`；VTK 为相同完整四边形，粗细交界可能有悬挂节点，不宣称共形求解面拓扑。独立审核检查覆盖、无重叠、分类、2:1、JSON/VTK 一致与确定性。

CLI 均匀/自适应资源上限为 level 10/12；桌面分别为 9/10，且自适应全域最低层级最高 8。它们是资源限制，不是性能或任意几何成功保证。`desktop/src/core/background-grid.js` 阻止将背景数据当流体网格。

## 原生层流求解

Linux 云端接续使用 GCC（本次 13.3）和 CMake，关闭与不可压目标无关的 `CARTMESH2D_BUILD_CHEMISTRY`。API/CLI 省略稳态方法和对流格式时按上下文选择；桌面“按工况选择”在保存工况时解析为明确方法及格式。显式 `--steady-acceleration none` 保持 SIMPLE，并在未显式指定格式时保持 Upwind。云端原始圆环入口、完整验证范围和源码哈希见 `artifacts/current/native-laminar-cloud-default.json` 与 `native-laminar-default-convection.json`；实际原始场位于忽略提交的 `outputs/cloud-laminar/`。Linux 结果不代表 macOS App 已验。

`Incompressible2D.cpp` 使用 SIMPLE、Rhie–Chow 及共享压力/黏性面通量；只读取最终 `*.solver.cm2d`。完整参数查 `build/cartmesh2d_flow_cli --help`。

```sh
build/cartmesh2d_flow_cli --mesh outputs/channel.solver.cm2d --case channel --export-boundaries outputs/channel.boundaries
build/cartmesh2d_flow_cli --mesh outputs/channel.solver.cm2d --case custom --boundary outputs/channel.boundaries --output outputs/channel-flow --nu 0.1 --speed 1 --max-iterations 5000 --tolerance 1e-9
```

命令示范结构，几何须先满足模板约束；边界文件由原生导出后修改，保持面编号/几何绑定。`pressure-opening` / `symmetry` 只支持轴对齐面。压力为 p/ρ（m²/s²），流量为每单位厚度 m²/s、向外为正。

| 控制 | 约束 |
| --- | --- |
| `--convection upwind\|limited-linear\|face-limited-linear` | 显式迎风或限制重构，不是物理模型切换；省略时，strict 固定物性 Newton 稳态选 face-limited-linear，其余上下文选 upwind |
| `--pressure-preconditioner ic0\|aggregation\|cholesky` | Cholesky 只在 macOS；仍验真实线性残差 |
| `--linear-policy strict\|adaptive` | adaptive 最终严格复核 |
| `--velocity-relaxation` / `--pressure-corrections` | 默认 .6 / 4；困难曲壁可能不稳，连续性小不代表收敛 |
| `--steady-acceleration none\|anderson\|newton-krylov` | 省略时 strict 固定物性稳态选择 Newton，其余上下文保持 SIMPLE；显式 none 保持原方法。Anderson/ Newton 限固定物性稳态，Newton 要求 strict 线性及收敛控制，最终使用原严格验收 |
| `--initial-guess` / `--initial-flux` | 同目标网格稳态初值；面初值须与单元初值同用，不是物理检查点 |

输出前缀在输入读取完成、首次写入结果前将 `.json` 标记为 `running`；最终摘要先写 `.json.tmp`，所有请求的输出关闭成功后才替换 `.json`。启动后的求解或导出异常标记为 `failed`，硬中断保留 `running`；只有非定常失败摘要记录最后接受时间。旧场文件可留作证据，读取端须以本次进程返回值和摘要状态共同判断，不能仅凭场文件存在认定完成。

`newton-krylov` 的入口为 `src/fv/FlowNewtonKrylov2D.cpp`，复用现有单次 SIMPLE 映射 G，同时求解 `F(x)=x−G(x)` 的速度、压力、面通量固定点；`Incompressible2D.cpp` 的原数值步骤保持。速度除以 Uref、压力除以 `Uref²+νUref/H` 并乘单元面积占比的平方根；面通量除以 `Uref×面长` 并乘面长占比的平方根。因此候选择优量无量纲，并且不让小单元个数直接支配整体范数。

内部 GMRES 最多 60 个方向、两遍正交化，方向的相对线性残差目标 .1；按现有严格内解的 `1e-11` 精度取有限差分步长 `sqrt(1e-11)×sqrt(1+||x||₂)`。单次映射的动量行误差控制收紧为 `min(用户容差,100×1e-11)`，因为原行目标为 `.01×容差×Uref×速度松弛`；这避免内解误差主导雅可比差分，不是新增物理精度门。Armijo 下降比例 `1e-4`、最多 12 个减半步长候选，仅用于约束实际完整固定点残差下降及试算成本。上述算法控制服务于困难耦合求解，未按单个圆环误差调参数；最终仍使用用户容差和原动量、场变化、局部/全局守恒门，对实际返回场重新执行一次普通 strict SIMPLE 更新。方法原则见 [PETSc 的 inexact Newton 说明](https://petsc.org/release/manual/snes/#inexact-newton-like-methods)，实现不依赖 PETSc。

最终认证与 Newton 候选择优分开：`FlowConvergence2D.hpp` 共享原严格门槛；认证重新装配最终场的真实动量残差，并核对未缩短的原生更新量和局部/全局守恒。认证的线性设置恢复用户请求，全部工作纳入预算，失败不替换已接受场。普通 SIMPLE 和非定常仍保留十步暖机；Newton 不要求再串行执行十次 SIMPLE，因为这额外测试了另一种迭代的收缩性。真实圆环重放已显示，旧十步认证能把首步约 `1.8e-10` 的动量残差放大到 `0.029`。该调整改变认证流程，不改变方程、容差或物理资格；方法背景见 [PETSc 收敛测试说明](https://petsc.org/release/manual/snes/#convergence-tests)。已有解也必须重新严格复核，预算不足、认证前取消、速度/压力/通量扰动均有实际求解回归。

此模式中 `--max-iterations` 是**全部原生 SIMPLE 评估预算**，包括失败试算、雅可比方向、线搜索和最终严格复核；摘要 `iterations=coupledEvaluations`，`coupledOuterIterations` 才是接受的 Newton 更新次数。历史只保存已接受候选的映射结果及最终复核，不把试算作为接受状态；预算不足返回未收敛，无法继续降低残差返回 `nonlinear_stagnation`。取消保留最后完整接受的场/通量，调用方异常向上传递；`coupledFailedEvaluations` 和 `coupledLastFailure` 保留候选失败数量及最近原因，不能隐去 NaN 或内解错误。

`--profile` 的 solveSeconds 包含全部试算、失败及复核；细分线性时间/计数只覆盖完整返回的评估，失败次数单列。当前每次映射重建求解工作区，60 个方向还会占用额外内存，不能只按 Newton 外迭代数宣称提速。材料联算、非定常及 adaptive/engineering 控制显式拒绝。代表算例命令、网格哈希、直接 CSV 场差、真实 App 保存重算与检查日志统一见[当前证据](../artifacts/current/native-laminar-newton-krylov.json)；实验 CSV 用 gzip 无损保存，重放命令前按证据恢复临时输入，结束后清理可再生成文件。

动量、局部/全局守恒、场变化和线性精度分开检查。稳态面通量的旧通量缺陷项须与松弛一致，不能为提速删除。macOS CLI 首次 Accelerate 调用前固定单线程；研究入口也固定 `VECLIB_MAXIMUM_THREADS=1`。

非定常使用后向欧拉、固定/自动 CFL、拒绝和重试；只有接受步推进时间。checkpoint 绑定网格、物性、边界与格式；`flow.time-history.csv` 是本次接受步，`flow.residuals.csv` 是最后尝试的内迭代。温度联算必须用联合 thermal checkpoint。工况 v4 保存线性精度及松弛，旧格式按明确默认读取；实际运行设置须与请求一致。

`SmoothMovingWall` 保留平滑壁速梯度，与恒定壁迹不同；圆环用原规则多边形法向构造切向速度，不能偷换成解析圆。被动温度 `D=k/(ρcp)`，源项 `Q/(ρcp)`，向外通量 `q/(ρcp)`；无温度反馈、浮力、辐射或共轭传热。SST 输运/壁距/RANS 代码及扩散 CLI 仍保留参与测试，未取得通用湍流资格。

## 可压层流开发入口

原生理想气体质量、二维动量与总能量方程，可选恒 μ/k。默认 Rusanov 一阶；`--flux hllc --order 2` 开启压力感知 HLLC/HLLE、受限线性重构及 SSPRK(2,2)，不保证激波/退阶处二阶。

- `EulerFlux2D.cpp` 以真实面向量计算共享守恒通量。Roe/Einfeldt 波速包络、星状态检查失效时显式计数并回退 Rusanov；接触恢复权重取邻面压力比三次方的最小值，`F=F_HLLE+ω(F_HLLC−F_HLLE)`。
- `Euler2D.cpp` 距离加权最小二乘重建原始量，速度在随网格旋转的局部坐标限幅；近奇异或重构失正时计数并退常量，不裁剪接受场的密度、压力或能量。
- 固壁从镜像 Riemann 压力得到法向牵引，对流质量/能量严格为零。每个 RK 阶段用真实边界/周期平移并验正性和 CFL；导出通量为阶段均值、波速取最大，失败候选不进入检查点。
- `EulerStepper2D` 缓存不可变几何、物性、边界和输运算子；`verify_euler.py` 从原多边形独立重建最后接受步两个阶段，同时检查相邻历史质量/能量收支，不冒称重放全部轨迹。

方法出处：[Toro HLLC](https://www.prague-sum.com/download/2012/Toro_2-HLLC-RiemannSolver.pdf)、[Barth–Jespersen 限制重构](https://ntrs.nasa.gov/citations/19890037939)、[Simon–Mandal 式49–50](https://arxiv.org/pdf/1803.04954)。本实现的多边形邻面传感器混合全部分量，不是论文选择性 HLLC-ADC 的原样复现。

```sh
build/cartmesh2d_euler_cli --mesh final.solver.cm2d --output outputs/euler/run --case sod --gas-r 1 --end-time .2 --flux hllc --order 2
ctest --test-dir build -R '^cartmesh2d_euler_' --output-on-failure
```

### 总能量耦合导热与黏性

| 项目 | 离散与边界 |
| --- | --- |
| 热闭合 | `T=p/(ρR)`，`cv=R/(γ−1)`，`q=−k∇T`；总能量增加 `−div(q)`，体积热容用 **ρcv**，不能套被动热输运的 ρcp |
| 热面通量 | `S=τd+C`，`τ=\|S\|²/(S·d)`，`Qf=−k[τ(TN−TP)+gf·C]`；Dirichlet 用真实面心，Neumann 按 `n·∇T=−qout/k`，不能用偏斜中心连线充当法线 |
| 黏性本构 | `τ=μ[∇u+(∇u)^T−(2/3)div(u)I]`；气体 2/3 不因网格二维改成 1，纵/横扩散系数为 4μ/(3ρ)、μ/ρ |
| 黏性通量 | 动量为 `−τ·Sf`，能量为 `−u_face·(τ·Sf)`；已有机械功，不再重复添加体积耗散热源 |
| 面梯度 | 完整速度梯度含交叉导数、非正交项，内面修正满足 `Gf·(CN−CP)=uN−uP`；面速外推至真实面心 |
| 壁面 | 静止或仅切向移动的 no-slip 使用指定壁速求应力/功；滑移壁零切向牵引；法向移动网格拒绝 |
| 热边界 / 开边界 | 绝热、温度 K、向外热流 W/m²（负值加热）；开边界零 Fourier 通量及零黏性牵引，仍有对流输运 |

内部面只算一次、邻格反号；周期配对平移格心并严格反号。热与黏性通量均做 RK 阶段平均。静止绝热壁质量/总能量通量精确为零，移动壁功使用给定壁速。

`HeatConduction2D` 组装完整温度残差矩阵，速率为行绝对值和除以 `2Vρcv`；`ViscousStress2D` 组装完整 2N×2N 动量 Jacobian，每格取两分量行范数的最大值除以 `2ρV`。二者包含非正交、周期、交叉导数和扩展壁面系数，每阶段按密度更新，与声学速率一起限制步长。这不是完整非线性能量 Jacobian，也不是任意 Cut-cell 的单调性/稳定性/熵证明；非 M 矩阵行显式报告，失正缩步或失败，绝不修剪。

参考：[NASA 气体应力与总能量](https://www.grc.nasa.gov/WWW/wind/TFAWS2007/Formulation.pdf)、[cv 与 cp](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-heats-cp-and-cv-1/)、[corrected snGrad](https://doc.openfoam.com/2212/tools/processing/numerics/schemes/sngrad/rtm/corrected/)。μ、k 分别给常数，不强制 Prandtl 关系；规则网格阶数不外推[任意非结构网格](https://openfoam.org/release/2-3-0/numerics/)。

检查点：μ=k=0 为 v1，μ=0 且 k>0 为 v2，μ>0 为 v3。物性及全部逐面机械/热壁绑定；v3 还含切向壁速。custom 文件 v2 增加热条件，v3 增加 wallUx/wallUy，保留读取旧版本；custom 条件不能被预设开关覆盖。通量、阶数、CFL、壁面数值格式和终点可显式调整，相同控制下中断/续算要求逐字节一致。

输出分别保存对流动量/能量、黏性应力/功、Fourier 热流、速率、组合 CFL 和历史边界收支；负向外壁功代表机械能输入。`hllcFallbackStages` 的位 0/1 区分两阶段，周期只计一次。独立热/黏性脚本和 `verify_euler.py` 重建算子、最后步及历史收支。

### 可压壁面精度与时间细化

`--wall-gradient quadratic` 是定温热壁/无滑移应力的可选数值格式，默认 linear；其余壁面模型不变。`WallGradient2D` 固定真实壁面中心值，在法/切坐标拟合 `n,t,n²/2,nt,t²/2`，沿共形邻接扩展 2–5 层并使用周期像点，不跨物面。法/切尺度及列归一化后用列主元 Householder QR，超过 512 像单元或缺秩就失败，不静默降阶。

秩门 `sqrt(machine epsilon)` 保护约半数双精度有效位，不是物理门限；计算时先减壁值消去常量。全部扩展系数进入热/动量行范数和组合 CFL，壁功仍用指定壁速。输出 `quadraticHeatWalls` / `quadraticViscousWalls` 记录实际面数，独立读取器用再正交 Gram–Schmidt 复核，与原生 QR 区分。

参考 [White/Nishikawa D.2.3](https://ntrs.nasa.gov/api/citations/20210024196/downloads/white_and_nishikawa_afang_paper_v_1.7.pdf)。这里是原生二维、壁值约束和质心原始变量点值近似，没有三阶保守量单元平均重构或高阶面求积，**局部二次不等于全局三阶**。

| 验证入口 | 检查什么 |
| --- | --- |
| `tests/wall_gradient_test.cpp` | 二次场、扭曲/旋转、长宽比 .03/3/30、长度缩放、实际 Jacobian 扰动及缺秩拒绝 |
| `tests/transport_precision_test.cpp` | 无源光滑导热 8/16/32 格，比较解析壁面边积分、排除线性求解容差污染 |

测试阈值按量分别解释，均为既有门，整理文档不改断言：面通量审计使用含对流/热/黏性量纲包络的 512 epsilon，算子/SI 采用 1024–4096 epsilon 检查浮点一致性；极端长宽比另记录输入舍入传播，不能替代原几何门。物理回归使用各自归一化离散误差、网格观测阶和解析参照，不用机器 epsilon 充当工程精度。

例如光滑导热 `T=2+.2sin(πx)sinh(πy)/sinh(π)` 的热流 L1 以解析边积分绝对值总和归一化，最细目标 .5%、至少优于线性 2 倍、观测阶>1.7；1e−11/1e−13 线性容差场差小于离散误差 1%。Couette 用解析温升和壁速归一化，另验压力细化，避免多项式恰好复现掩盖全场误差。准确控制、门限依据和成本在相应测试及[精度证据](../artifacts/current/native-euler-wall-accuracy.json)中。

时间细化使用同网格/终点，dt=1e−3/5e−4/2.5e−4，对更细轨迹再核参考误差，属于时间自收敛。`EulerStepControls2D::endTime` 保留浮点尾步最小失败例：必要时把倒数第二步拆成两个合法小步，真实算通量，不伪造时钟或绕过最小步长。原生求解器继续拒绝不合法的通量、速率和壁面格式。

```sh
ctest --test-dir build -R 'cartmesh2d_(wall_gradient|transport_precision)' --output-on-failure
build/cartmesh2d_euler_cli --mesh final.solver.cm2d --output outputs/euler/run --case external --viscosity .02 --wall-model no-slip --conductivity 100 --wall-thermal temperature --wall-value 400 --flux hllc --order 2 --wall-gradient quadratic --density 1.225 --pressure 101325 --u 50 --end-time .00005
```

最后一条用短时功能测试物性，不是空气推荐值。精度脚本可能需数分钟，Windows 预算更长；日常只选相关项，完整范围与未验物理问题见当前状态。

## 详细燃烧化学基础

`include/cartmesh2d/chemistry/DetailedGas.hpp` 与 `src/chemistry/DetailedGas.cpp` 是详细反应流的原生热化学基础，使用 Cantera **3.2.x C++ API**。模块不依赖 Python 或三维核心；验证脚本另用 Python Cantera 3.2.0。机理必须给出实际文件路径；当前明确接受中性、单气相理想气体机理，其他热力学/相模型显式拒绝。

- 保守状态为 `rho`、`rho*e` 和各 `rho*Yk`，其中 `e` 包含生成能，允许负值。输入组分不自动归一化或裁剪；能量到温度的反解限制在所有组分热力学数据的共同温区。积分阶段的显式舍入闭合见下文，不能用于修复不合法输入。
- 返回温变比热、焓、反应源，以及 `multicomponent` 扩散矩阵和 Soret 热扩散系数；矩阵为列主序，须配合完整通量公式，不能当作各组分独立的 Fick 系数。
- 恒容绝热化学子步使用 Cantera `Reactor` / CVODES，内部能量为守恒变量；化学变化通过组分及温度体现。`-sum(hk*omega_k)` 只作放热诊断，不重复加入已经含生成能的总能量方程。反应阶段超出共同物性温区也会拒绝。
- 失败返回空 `accepted` 和原因；输入不变，下一次调用从调用者保留的接受状态重新开始。当前是串行、每工作线程独享的化学上下文；空间扩散、耦合时间推进和原生检查点见下文，产品输入与 App 尚未接入。

`thermodynamics(q)` 返回独立的 `GasThermodynamics` 类型，包含温度、压力、比热、组分及生成焓，不计算反应速率。原 `properties(q)` 仍返回全部热力学与化学源。对流、扩散和状态验证使用前者；耦合组分方程继续从完整机理取得反应源。两种接口共享原温区、正性和能量反演检查，反应源的有限值检查保留在化学查询中，不用零向量代表省略的反应。

热化学上下文只缓存最后一次成功反解的守恒状态：密度、内能密度及每个组分密度逐位比较，包含正微量和符号零，不按温度邻近或误差容差合并状态。完全相同时可跳过重复的温度反解；更换状态先使缓存失效，直接设置 T/p/组分和进入共享相的反应积分前也使其失效。失败试算不能留下可命中的旧状态键。每个查询的实际物性/输运仍由完整后端求值；未加入物性插值表或降低机理规模。`detailed_gas_test.cpp` 覆盖查询顺序、输运/源不变性、失败反解和反应后的恢复，以及 `1e-250` 正微量变化不能被缓存忽略。

有限的 `rho*e` 和 `rho` 相除仍可能溢出，因此反演入口还显式检查比内能有限。否则无穷目标会使舍入预算也变为无穷，从而误接受无关温度。该修复是非法输入拒绝，不改变任何精度容差；原失败和修复后的相关检查保留在 `outputs/combustion-foundation/thermo-query-001/`。

查询改动的[证据索引](../artifacts/current/native-reacting-query-efficiency.json)绑定冻结的原/新二进制、后端和源文件哈希。该目录的 `compare_runs.py` 记录原/新/新/原冷启动、100 微秒终点续算、显式推进及两种完整机理反应器的实际命令；对流动输出只排除进程计时元数据，其他字段逐值比较，步长、接受状态、检查点和展开机理另作逐字节比较。读取器继续复核实际共享面、EOS、源项和物理预算；无新物理阈值。复现应创建新的输出目录并重新绑定编译来源，不覆盖已保存记录。完整成本与未验范围只在[当前状态](CURRENT_STATE_CN.md#持续目标成熟燃烧模拟)维护。

依赖按[官方 C++ 构建说明](https://cantera.org/3.2/userguide/compiling-cxx.html)准备；源代码版本为官方 `v3.2.0` / `4a8358eb80cfeb50474386b5f9ec0b3a83519889`。本机源码和安装分别在忽略提交的 `build/deps/cantera-src`、`build/deps/cantera`，用系统 clang 和上游固定子模块编译，Boost 1.88 仅使用头文件。没有全局安装依赖。

```sh
# 已有 Cantera 3.2 C++ 安装后，启用独立化学目标；其余产品默认构建不要求此依赖。
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/clang++ \
  -DCARTMESH2D_BUILD_CHEMISTRY=ON -DCARTMESH2D_CANTERA_ROOT="$PWD/build/deps/cantera"
cmake --build build --target cartmesh2d_detailed_gas_tests cartmesh2d_detailed_gas_probe -j2
ctest --test-dir build -R '^cartmesh2d_detailed_gas$' --output-on-failure

# 参数为机理文件、温度 K、压力 Pa、恒容反应时间 s、相对摩尔数量。
build/cartmesh2d_detailed_gas_probe build/deps/cantera/share/cantera/data/gri30.yaml \
  1400 101325 .002 CH4=1 O2=2 N2=7.52
```

本机依赖构建环境为 Python 3.14、SCons 4.11.1、packaging 26.3；Python 对照环境另装 Cantera 3.2.0、NumPy 2.5.3。Cantera 构建使用 `python_package=n f90_interface=n hdf_support=n clib_legacy=yes googletest=none example_data=no debug=no`，`fmt/yaml-cpp/Eigen/SUNDIALS` 使用上游固定子模块。`clib_legacy` 仅避免生成本项目不使用的 C 接口，不改变 C++ 模型；上游完整测试没有运行，原生接入检查单列。动态依赖需检查，不得引入 mesasdk 库。正式打包和跨平台部署仍待完成。

检查按用途分开：原生测试核对状态反解、元素源、完整机理规模、负生成能、扩散数据、反应积分及失败恢复。元素源舍入检查以**总生成与消耗量**归一化，不能用已相消的净源作舍入尺度；首轮甲烷氮元素因此误报，原日志保留在 `outputs/combustion-foundation/first-native-test-failure.log`。接口脚本直接计算 NASA7 多项式，独立核对混合物比热/内能；反应轨迹和扩散与另装的 Cantera Python 接口比较，仍共享化学后端，不能称为独立机理或实验验证。

阈值只限定数值接入误差：温度往返 `2e-9 K`、NASA7 及接口值采用脚本声明的归一化差异；积分对照检查温度相对差和质量分数绝对差。默认化学子步的元素质量分数漂移门为 `1e-8`，能量变化以 `max(|rho*e|,rho*cv*T)` 归一化后为 `1e-8`；它们是可配置的子步守恒检查，不是火焰精度目标。上述少量案例为秒至分钟级，网格/火焰/真实工况验收尚未完成。

### 可选 NASA7 连续热力学输入

原始分段 NASA7 拟合的焓/内能可能在拼接点跳变；在很窄的能量区间，温度反解因而不唯一。`tools/flow/prepare_continuous_nasa7.py` 显式生成新的机理文件，保留低温段锚点、两段全部比热系数、温区、组分、输运参数及完整正向反应定义。依据 [NASA7 形式与比热积分关系](https://cantera.org/stable/reference/thermo/species-thermo.html#the-nasa-7-coefficient-polynomial-parameterization)，在拼接点以 60 位算术计算高温段的 `a5`、`a6`，使 `h = h_ref + ∫cp dT`、`s = s_ref + ∫cp/T dT` 连续；最终系数仍为双精度。比热本身的拼接跳变和导数未被平滑。这里只支持理想气体的两段 NASA7，其他模型显式拒绝。

这是有物理影响的数据变更：高温焓、熵、平衡常数和逆反应速率会改变。工具不修改原文件，也不接入求解器的隐式修正。输出 `original-resolved.yaml`、`continuous.yaml` 和每个组分的改动/哈希报告；已有输出目录拒绝覆盖。Cantera 展开原机理时可能因活化能单位换算产生末位舍入，报告逐项列出；准备器直接修改展开后的 YAML 树，保证两份展开文件的正向反应定义相同，另用原输入进行敏感性对照。原相配置的输运模型也保留，实际原生通量继续显式使用完整多组分模型。

```sh
build/chemistry-env/bin/python tools/flow/prepare_continuous_nasa7.py \
  --mechanism build/deps/cantera/share/cantera/data/h2o2.yaml \
  --output outputs/combustion-foundation/continuous-h2-new
```

甲烷对照使用 `gri30.yaml` 与 `--fuel CH4`，仍为全部 53 组分/325 反应。Python 准备环境需 Cantera、NumPy 和 `ruamel.yaml`，不增加原生求解器运行依赖。验证器按各段对 `cp`、`cp/T` 作独立 32 点 Gauss 积分，检查全部组分低温锚点至拼接点两侧和高温区间的焓/熵；复用公式检查归一化门 `5e-11`。焓误差除以 `max(|h|,|cp*T|,1)`，熵误差除以 `max(|s|,1)`；拼接处还使用原能量反演舍入预算 `64*epsilon*max(|e|,|cv*T|,1) J/kg`。每种机理实际调用原生探针 8 次，检查拼接附近往返和原/新机理的恒容反应轨迹；接口归一化门 `5e-10`、反应器温度相对差/质量分数绝对差 `2e-6` 与既有守恒门不变。它们限定数值一致性，原/新机理的物理差异直接报告，不新增物理合格阈值。

已保存的原始、失败和最终输入在 `outputs/combustion-foundation/thermo-continuity-001` 至 `004`；最终机制文件与已计算火焰输入逐字节相同。`003/compare_flame_2.py` 从保留的参考火焰恢复同一 459 点网格，分别显式设置 400 K、1 atm 和入口组分，核对机理系数、完整多组分/Soret/能量开关，再在相同控制下实际重算两种输入，输出原始温度/组分和差异。原 BVP 的微小负舍入值保留，不能当作合法原生接受状态。此对照没有新增网格细化或实验资格，恢复/缓存造成的计时差异不作提速结论。数值、脚本和数据哈希见[证据索引](../artifacts/current/native-reacting-thermo-continuity.json)。后续原生空间计算须用连续机理重新生成并绑定参考初值；旧机理检查点不可改名后沿用。

### 多组分空间扩散与焓通量

`DetailedGas::diffusiveFlux` 接受沿单位方向的 `dT/dn`、`dln(p)/dn` 和各 `dXk/dn`，计算相对于质量平均速度的组分通量。采用 [Cantera 多组分矩阵的通量约定](https://www.cantera.org/3.2/reference/onedim/governing-equations.html#diffusive-fluxes)，驱动力为 `dk = dXk/dn + (Xk-Yk)*dln(p)/dn`，压力项依据 [PeleLM 完整输运模型](https://amrex-combustion.github.io/PeleLM/manual/html/Model.html)。矩阵使用分子量加权及正确的列主序，包含 Soret；能量通量为 `-lambda*dT/dn + sum(hk*jk)`，组分焓保留生成焓。当前中性理想气体模型未含差异体力、电迁移或 Dufour 热通量，不能声称覆盖所有压力/输运物理。

`ReactingDiffusionOperator2D` 位于 `include/cartmesh2d/fv/ReactingDiffusion2D.hpp` 和 `src/fv/ReactingDiffusion2D.cpp`，以通过原有质量门的真实 `FvMesh2D` 为输入。一次计算共享面通量，再正负加入 owner/neighbour；输出为实际边长积分、单位厚度的组分/能量通量和单元向外残差，时间导数应取 `-residual/真实单元面积`。算子本身不推进时间、不改变密度/动量、不调用反应积分，也不添加化学热源。

梯度使用几何加权最小二乘及现有非正交面修正；只重构 N−1 个独立摩尔分数梯度，其余由 `sum(grad X)=0` 确定，不改变守恒组分。面物性用 T、ln(p)、Y 的凸插值，不裁剪负值。边界支持固定状态储库、绝热不透壁、定温不透壁，缺失/重复边界和非法状态显式拒绝。不透壁直接约束全部组分法向通量为零，定温壁的能量通量只含导热；壁面梯度辅助模板采用齐次法向组分/压力导数。复杂曲壁热扩散边界层的精度、凸插值在任意偏斜网格上的收敛、周期边界及隐式多组分时间推进仍待完成；一次守恒测试不代表这些已经合格。

浓度与压力通量直接求解以质量通量为未知量的 Stefan–Maxwell 方程，避免用逆扩散矩阵的大项相减求出微小通量。消去最丰富组分的质量约束；恰为零的组分使用解析标量行，其余组分采用驱动力相关的列尺度、行尺度及主元消元。组分为零但梯度非零时仍有扩散。原始多组分 D 矩阵保留作差异诊断，不裁剪正浓度或浓度梯度。线性检查采用各行运算绝对总量乘 `64*(N+1)*机器 epsilon`；IEEE 次正规通量另计各未知量半个最小表示间隔传播到该行的绝对舍入界。线性未知量的缩放单位保持在正规浮点范围，最终通量仍可小于该单位；这不是物理组分下限。

`src/chemistry/TraceStableMultiTransport.hpp` 保留 Cantera 3.2 的碰撞积分、转动/内部热模式及完整 3N 动力学矩阵，改用质量平均速度约束下的速度差方程和按组分丰度缩放的方程求解热扩散。依据[上游完整输运实现](https://github.com/Cantera/cantera/blob/v3.2.0/src/transport/MultiTransport.cpp)，L00 在质量约束下可写成 `(16*T/25)*sum(Xj*(a0j-a0i)/Bij)`；这一等价形式避免从 O(1) 项相减恢复 O(Xi) 量。热扩散系数用真实 X 回代，并变换到真实质量平均参考系；内部矩阵保留上游 `max(X,1e-20)` 的数值正则化，其丰度扰动为 O(N×1e-20)，不改状态、反应机理或输出浓度。每次求解检查线性后向误差。该扩展使用 3.2 的受保护 API，升级 Cantera 时必须重新核对和验证。

通量闭合门仍只允许运算绝对活动量对应的浮点余差；浓度/压力项相消前的量保留作尺度。通过后只在最丰富组分上闭合质量。严格零组分且零梯度的通量、纯气体 Soret 分离均为零，`zeroLimitCorrection`/`pureSpeciesLimit` 保留极限诊断。纯水蒸气闭合、逆矩阵微量消减、热扩散微量平台、MUSCL 零界舍入和次正规通量检查的失败记录与对应回归保留。

测试以总面传递量归一化检查共享面装配和闭域预算，并用仿射场检查非正交梯度；这些是浮点/离散算子的直接检查，不是通用物理误差门。局部通量脚本用 NumPy 从二元扩散系数组装独立 Stefan–Maxwell 线性系统，并核对二元解析极限；逐个检查微量组分，避免被主组分的总误差尺度掩盖。Soret 微量极限另用未修改的 Python Cantera，在可解析浓度上计算 `D_T/X` 并作两档二次外推；该相对误差门为 `2e-7`，允许双精度与外推三阶余项，已有普通接口算术门仍为 `2e-9`。这些共享碰撞/热力学数据，未验证实验输运数据。脚本核对运行前后的可执行文件 SHA256，拒绝验证中替换二进制的混合结果。运行成本为秒至数分钟，最新完成范围见当前状态及证据。

```sh
cmake -S . -B build # 沿用上文已配置的 Cantera 开关和路径
cmake --build build --target cartmesh2d_reacting_diffusion_tests cartmesh2d_diffusive_flux_probe -j2
ctest --test-dir build -R '^cartmesh2d_reacting_diffusion$' --output-on-failure
```

### 详细反应流耦合推进

`ReactingFlowStepper2D` 位于 `include/cartmesh2d/fv/ReactingFlow2D.hpp` 与 `src/fv/ReactingFlow2D.cpp`，使用真实 `FvMesh2D`，是独立的平面详细反应 Navier–Stokes 开发核心。旧 Euler/层流求解器的默认控制不变。守恒变量顺序为 `rho, rho*u, rho*v, rho*E, rho*Yk`；`E=e+|u|²/2` 含生成能，允许负值，恢复物性时先扣除动能。没有额外添加化学放热源。

- 对流采用端点冻结声速 `sqrt(cp/cv*p/rho)` 的 HLLC；中间态必须通过真实混合物 EOS 和共同温区检查，否则明确记录 HLLE 回退。星区通量采用等价的 `rho_star*S_contact` 形式，组分使用同一个质量通量乘迎风 Y，避免近零流量时各守恒分量分别作大项相消。声速随组分与温度变化，未使用恒定 gamma。二阶模式在 T、ln(p)、速度和 Y 上作最小二乘 MUSCL 重构与局部极值限制；N−1 个独立组分分别按各自邻域限制斜率，再由梯度和确定最丰富组分。仅当这个依赖组分超出物理 `[0,1]` 时，才共同收缩所有组分斜率到可行单纯形；独立组分保留邻域界，依赖组分不另保证局部极值界。这样避免无限小的示踪组分起伏使主要组分的通量发生有限跳变。限制因子保留 `64*epsilon` 的乘加舍入余量，不裁剪面状态或守恒状态。
- 输运包含多组分/压力/Soret 扩散、导热、生成焓输运，以及面温度/组分对应的黏性应力和功。每面只算一次完整通量，正负加入相邻真实单元。化学 `dt/2 → 输运 dt → 化学 dt/2` 使用 Strang 分裂；输运为 SSPRK2，所有候选在私有副本上计算。负组分、物性越界、CFL 上升或化学积分失败均拒绝并减半重试，保留每次原因；不改调用者状态/时钟，不把步数上限标为达到终点。
- 步长受声学、变量黏性和完整热化学扩散约束。扩散约束用与实际非正交模板相同的冻结系数 Jacobian 行和，量纲为 `1/s`。在固定密度下将冻结守恒 Jacobian 相似变换到 `delta T/T` 与 N−1 个独立 Y：输入组分扰动保持 T，输出温度导数由 `d(rho e)-sum(ek*d(rho Yk))` 除以 `rho*cv*T` 得到。变换保留该冻结线性算子的特征值及生成能耦合，避免纯粹由能量尺度造成的过大行和；实际更新仍完全守恒。这不含物性系数的非线性导数，也不是任意网格稳定性证明。默认 CFL `.35`、步长预测余量 `.9`，每个输运阶段检查真实 CFL，失败照常重试。上述数值只属于新反应流模块。
- 边界支持静止几何的滑移壁、无滑移壁（可有切向运动）、绝热/定温壁、固定状态储库及外推流出。流出采用零扩散通量/零黏性牵引，回流必须显式提供储库状态。周期、非反射特征边界、复杂壁面热扩散层、欠解析的多组分移动接触面和强激波相互作用尚未验收。

`GasMechanism::resolvedDefinition` 由 Cantera 将相、全部导入组分/反应及输运数据展开为独立 YAML，去除日期、生成器及可变初态，保留 17 位精度。检查点绑定这份完整定义、后端版本、单元/面几何拓扑、全部边界及化学/输运开关；不能仅凭文件名或反应个数续算。引用的化学上下文被替换也会拒绝。`writeCheckpoint/readCheckpoint` 保留全部守恒量、实际时间和接受步数；初值、截断文件、额外尾部和不匹配绑定显式失败。已接入下述火焰验证驱动的文件续算；产品级原子文件保存、取消交互和用户输入工作流仍待接入。

检查点标量以 [`std::from_chars`](https://eel.is/c%2B%2Bdraft/charconv.from.chars) 完整解析，检查转换状态、整段消费及有限性。本机标准流提取对可表示的非正规极小正数设置失败标志，曾将真实火焰的完整文件误报为截断；现在保留这些值，不置零、不设下限。原失败文件含 8 个此类值，最小约 `4.323e-321`。最小回归还包含 `double::denorm_min()` 的精确读回；无法表示的溢出/下溢、NaN/Inf 和尾随字符显式拒绝，物理状态仍须经过原合法性检查。

直接相关验证的依据与边界如下，日常成本为数十秒至约两分钟：

| 检查 | 归一化及判据用途 |
| --- | --- |
| 均匀氢气/甲烷耦合 | 与相同后端的完整恒容反应积分比较；温度相对差和组分绝对差 `<2e-6`，检查分裂/接口和重复加热，不验证机理物理精度 |
| 单步有限体积装配 | 各方程误差按同量纲的前后积分、全部面传输绝对总量、化学变化的最大值归一化，`<2e-11`；元素余额按初始总质量归一化，`<2e-8`。用于舍入/装配与反应积分预算检查 |
| 声波 | 显式冻结化学、关闭分子输运，以隔离可压对流；误差按初始压力扰动幅值归一化。12→24 格的面积加权 L2 为约 `.008738→.001556`；测试要求细档 `<.06` 且比粗档下降至少 30%，不是火焰精度门 |
| 黏性/扩散步长估计 | 仿射应力功与独立有限差分 Jacobian 对照；扩散在均匀全正组分状态作 `1e-6` 无量纲扰动，允许 `5e-5` 的差分/舍入余量。该检查不外推至任意非均匀火焰 |
| 独立场读取 | 重建真实多边形并核对面关联；质量和元素按初始总质量、总能量按初始积分 `max(abs(rho E),rho cv T,动能密度)`、壁面动量预算按初始质量×最大声速归一化，均 `<1e-8`。EOS 与输出能量/压力差 `<2e-9`，组分和绝对差 `<1e-10`；Python EOS 仍共享 Cantera 后端 |

验证驱动 `cartmesh2d_reacting_flow_probe` 使用 50 mm 方腔、48 格扭曲共形网格、封闭绝热无滑移壁、1 atm、`H2:2,O2:1,N2:3.76`，给定平滑初始温度分布后不再施加外部热源。到 `40 us`，默认设置接受 72 步；CFL 减半接受 143 步，均无试步拒绝。两档最大温差约 `0.204 K`，只作同网格时间敏感性比较，不宣称收敛阶或网格无关。图中每块颜色直接来自实际单元平均场；该粗网格早期点火检查没有解析火焰厚度。首版没有步长预测余量时的 109 次拒绝仍保留，不能把它删掉或用更少拒绝次数外推通用性能。

```sh
cmake --build build --target cartmesh2d_reacting_flow_tests cartmesh2d_reacting_flow_probe \
  cartmesh2d_detailed_gas_tests cartmesh2d_reacting_diffusion_tests \
  cartmesh2d_viscous_tests cartmesh2d_euler_tests -j2
ctest --test-dir build -R '^cartmesh2d_(detailed_gas|reacting_diffusion|reacting_flow|viscous_core|euler_core)$' --output-on-failure

# 每次使用不存在的目录；保存真实场、逐步拒绝原因、完整机理与最后接受检查点。
build/cartmesh2d_reacting_flow_probe build/deps/cantera/share/cantera/data/h2o2.yaml \
  outputs/combustion-foundation/reacting-new
build/cartmesh2d_reacting_flow_probe build/deps/cantera/share/cantera/data/h2o2.yaml \
  outputs/combustion-foundation/reacting-half-new .175
```

`reacting_reconstruction_test.cpp` 在扭曲网格上隔离二阶对流算子，将氩气质量分数扰动依次缩小至 `1e-12/1e-16`，并保持总密度和总能量。旧共同组分限制器仍造成约 `0.274%` 的整体质量通量变化；独立限制与单纯形闭合后，相对响应约 `1.05e-11/1.91e-14`。误差按相应分量的全体面绝对通量之和归一化，检查额度为 `1000*(扰动质量分数+epsilon)`，用于发现不随扰动消失的有限跳变，非火焰精度门。保留同源码仅恢复旧限制器的失败对照。

该连续性测试使用 1150–1350 K，避开原机理的 1000 K NASA 分段拼接。原 900–1100 K 诊断仍保留：独立限制已消除组分梯度跳变，但面温度跨拼接点仍会受数据不连续影响。例如 `H2:1.8,O2:1,N2:3.76`、1 atm 在 `1000±1e-9 K` 的 Cantera 内能相差约 `-0.140 J/kg`。原机理未修改；可选的[连续热力学输入](#可选-nasa7-连续热力学输入)已单独完成数据一致性和指定案例敏感性检查，但其原生空间长期影响尚未验收，不能将原连续性测试外推为跨拼接点的任意精度证明。

### 全耦合隐式反应流（开发模块）

`ReactingImplicitIntegrator2D` 用 CVODES BDF 同时积分既有半离散对流、黏性、多组分/Soret 输运与全部有限速率反应。它复用 `ReactingFlowStepper2D::evaluateResidual`，关闭的只是显式步长估计；查询估计前后物理残差和共享面通量逐值相等属于回归要求。该模块可选构建，不改变既有显式或产品默认控制。

每格积分密度、两分量动量、含生成能总能量和 N−1 个组分密度，剩余组分由质量恒等式确定。初始表示的舍入改变量和后续代数约束缺陷另列；实际质量、元素与能量预算不扣除它们。变量按初始密度、密度×声速及能量尺度无量纲化，保留负的生成能。局部误差权重默认相对 `1e-7`、无量纲守恒量绝对 `1e-12`、质量分数绝对 `1e-18`；这些是积分控制，不能直接视为火焰精度。BDF 默认最高二阶，验证驱动显式选择最高五阶。

线性代数采用确定性 RCM 排序和覆盖完整二环模板的带状直接求解。带矩阵预算默认 256 MiB，超出明确拒绝；当前不是大规模通用稀疏求解器。雅可比按无量纲物理尺度作分组差分，扰动为 `sqrt(epsilon)*max(1,abs(z))`，物性不允许时有界缩小扰动。默认 `jacobianAdvectionOrder=1` 仅让牛顿矩阵使用连续的迎风线性化，避免对 MUSCL 限制器的跳变作差分；实际场方程、非线性残差、最终诊断始终使用 `spatialOrder=2`。化学、黏性和完整多组分扩散在线性化中也保留。矩阵可选二阶对流线性化，慢速对照保留。差分尺度不给真实微量组分设下限，不替换真实 RHS；它与误差权重分别控制。默认后端的步长相关差分在零/微量组分上曾低于物性舍入尺度，导致步长越小、错误矩阵项越大，原始诊断保留。

不启用 CVODES 的不等式投影，因为其 `cvCheckConstraints` 可以将小负值修正为零。通过公开的非线性求解器接口回退整个牛顿修正，使迭代点满足物性和非负组分；不逐组分裁剪。回退后的迭代必须同时满足未缩小修正量与实际非线性残差的原局部容差，不能靠缩小更新制造收敛。每次成功返回前重新计算最终迭代的真实 RHS，再进行边界/反应/代数约束的独立时间积分。`reacting_flame_trace_newton.fixture` 保留原始参考第 98–104 列（两行），入口取第 97 列，用于复现默认牛顿在零氩气处生成负值的问题；它不是完整燃烧室或续算检查点。

非线性收敛现在统一要求原牛顿修正与真实残差的加权范数均不超过后端提供的局部容差。由 CVODES 从实际总预测修正计算时间误差，不调用其默认收敛回调的一次迭代快捷路径：非零可行初猜或回退使“本次牛顿修正”等于“总预测修正”的假设失效。完整 464 格输入的只读后端诊断复现了旧路径把 `0.0027553` 估为 `0.0881705`，相差 32 倍；修复后实际误差向量与控制器标量一致到舍入量级。该诊断读取了同版本后端内部结构，生产代码只使用公开接口；裁剪出的 14/28 格片段没有复现此问题，不替代完整输入证据。

接受前读取后端的局部误差向量和误差权重，检查无量纲加权范数不超过其时间步接受界限 1，额外允许 `64*epsilon*方程数` 的范数累加/缩放舍入。这检查积分器执行一致性，每步增加两次向量读取及一个范数，不是火焰物理精度或独立时间收敛标准。实际 BDF 阶数、牛顿失败、矩阵重建、误差范数和触发非负边界的组分进入验证日志。

验证驱动支持 `--bdf-order 1..5`、`--jacobian-order 1|2`、`--newton-iterations 正整数` 和 `--continue-damped 0|1`，读取器核对实际使用值。模块原误差容差、默认最高 BDF 阶数 2、一次矩阵近似及 3 次牛顿迭代预算保留；火焰验证驱动仍默认最高 5 阶。继续缩短后的牛顿迭代为关闭的研究选项；已运行对照未证实收益。最高 2 阶和二阶矩阵差分对照均增加了同终点方程调用；最高 1 阶试验提前停止，未取得同终点结论。失败/取消记录保留。

时间控制细化通过 `--rtol`、`--conserved-atol`、`--species-atol` 指定，读取器逐项核对实际值，原物理审计门不变。它们分别作用于相对局部误差、按初始物理尺度无量纲化后的守恒量绝对误差和按初始密度归一化的组分密度绝对误差。为量化时间积分对真实场的影响，本机在 464 格、10 微秒将三项默认容差共同减半、再减半；该验证没有新增物理精度门，也不改变产品默认控制。`compare_reacting_flame_time.py` 检查原始文件哈希、相同实际终点、完整初始状态/残差、网格、机理及其余积分参数，再报告温度、速度、压力、各组分和燃料消耗差异。不同二进制的源码溯源另审：本次区别是驱动参数及零时长输出处理，方程核心保持相同。自适应容差比例不能用来直接推算固定时间步的收敛阶。

`--reflect-species 1` 显式开启可选的非负组分边界反射候选步，默认关闭。完整牛顿候选不合法时，先将越过零边界的独立组分试探方向反射回可行侧，再检查整个候选的组分闭合、EOS 和真实非线性残差；只有残差严格下降才选用，否则继续原全向量回退。反射不修改输入或最后接受状态，不设浓度下限；原牛顿修正与真实残差仍须满足原局部容差，时间误差由实际总预测修正计算。概念参考[约束最小二乘中的反射搜索方向](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html)，这里实现一个反射候选及真实残差比较，不宣称完整 TRF 算法或其收敛理论资格。

完整 464 格、1 微秒案例实际使用了 190 次反射候选，最大组分改变量按各单元初始密度归一化约 `2.18e-39`；这描述内部试探方向，不是接受场裁剪或允许误差。原氩气库存约 `3.27e-23 kg/m`，扣除真实边界通量后的相对守恒缺陷仍约 `1.1e-15`。与原隐式解最大温差约 `8.08e-5 K`，方程调用减少约 93.4%。该对照的物理模型、网格、初值、局部误差容差和终点相同；完整进程计时包含启动、组网、积分与场输出，但部分诊断同时运行，因此调用量与观测耗时分别报告。长期传播、其他燃料和复杂几何的收益及精度资格待验证。

接受状态只在真实状态、物性和积分时钟检查后更新。会话保留 BDF 历史；步间取消可继续同一会话，内部试算取消或失败保留最后接受状态，但须用新会话重建历史。试算取消引起的后端中止单独报告，重复读取终止会话保留取消状态及原因，不复用已中断的历史。验证输出的 `integration.canceled` 与实际完成时钟分开记录。残差评估和接受步都有成本预算，不把预算耗尽记为完成。火焰验证驱动保存 `evaluation-progress.jsonl`；在输出目录创建 `cancel.request` 可有序停止并保存最后接受状态、检查点和失败原因。直接杀进程不提供同样的保存保证。

长时间验证可用 `--sample-every N` 将初始状态、每 N 个接受步及末态写入 `samples.jsonl`。采样不改变时间步或 BDF 历史，使用单独的物性和残差上下文；额外诊断调用在 `sampleResidualEvaluations` 中计数，仍包含在完整进程耗时内。默认关闭采样；开启后，`--max-samples` 默认 256 条，用于限制诊断输出数量，不是物理精度门。达到保存上限或诊断出错会在接受状态退出并尝试写入最后接受场和检查点，原失败原因进入结果；部分轨迹明确标记 `sampling_complete=false`，不冒充完整计算。

读取器用同一套独立时钟、EOS、共享面和质量/元素/能量预算检查逐帧核对。`initial` 表示本段起点：冷启动来自参考初值，续算来自真实接受检查点，两者通过 `restart` 元数据区分。`sample-audits.json` 保留每帧诊断；历史图使用真实接受时刻。固定本段初始温度范围中点的等温位置通过单元中心剖面线性插值得到，多个交点全部记录；该位置标记不等于已验收火焰速度。464 格、0.1 微秒采样开关对照的全部物理场相等，步长日志和检查点逐字节相同。故意改坏的时钟和质量预算被拒绝；两条样本上限测试在第 20 个接受步退出，保留的场和预算与正常运行的同一步相同。

试算终止后可能增加 RHS、矩阵及非线性迭代计数，而未接受任何新场。`audit_sample_counters` 仅在未完成、有失败原因、末帧物理状态严格等于最后接受状态时，区分这种终止快照与较早的常规采样；普通帧仍须与接受步 RHS 计数相等。终止计数只允许明确列出的单调诊断变化，积分控制、接受阶数计数和接受误差统计保持；每帧继续接受原 `audit_fields` 检查。报告另列 `terminal_attempt_statistics`，不把额外试算计入物理时间。原生取消/新会话恢复测试归入 `cartmesh2d_reacting_implicit`。

`outputs/combustion-foundation/implicit-cancel-audit-001/` 绑定冻结程序、审计器及源码哈希，保留真实试算取消、零时长文件读回、续算和正常计算对照。`reaudit_legacy.py` 在单独目录重审旧记录，并检查原物理审计函数 AST 相同及实际时钟/预算/状态损坏拒绝；`restart.py` 分段审计后按真实边界、化学和约束积分合并预算。原 `acoustic-boundary-001` 的取消记录及 `acoustic-newton-001` 的报告元数据失败保留，补充已存在的夹具间距只修复审计输入，不重算或改写原场。冷态脉冲使用完整氢气机理、多组分/Soret、温变物性和黏性；脉冲是初值，不能用它代替燃烧传播或边界反射验证。相同 10 微秒终点的三/六次牛顿预算及六次预算的容差减半对照只描述该案例，默认值不改。完整文件、控制及观察范围见[证据索引](../artifacts/current/native-reacting-cancellation.json)。

完整 100 微秒输出保留在 `outputs/combustion-foundation/implicit-evolution-100us-001`，进程实际使用 `implicit-sampling-002/probe` 和 `implicit-sampling-003/verify_reacting_flame.py` 的固定版本；不要用当前构建冒充旧运行来源。`analyze_completed_evolution.py` 从已核对哈希的实际场和样本读取不同等温位置，列出末段残差变化；温度平移诊断只在初始剖面的真实重叠区间作线性插值，不生成接受状态。不同温度标记的速度不能直接作为火焰速度资格。`residual-profile-001` 保留临时计时补丁、编译命令、调用数量及原/计时版本逐字节对照；它只定位成本，不构成提速证据。详细来源统一由[演化证据](../artifacts/current/native-reacting-evolution.json)索引。

`analyze_reacting_flame_history.py` 将已完成且独立审计通过的二维条带历史转为可复现的传播诊断。它校对报告、原始场、样本、参考与夹具哈希，逐条读取 JSONL，核对首末状态、实际时钟、总步数与样本数；零推进、部分采样和来源不符均拒绝。兼容旧零起点报告及真实续算的累计时钟；没有 `endTime` 或累计步数字段的旧记录必须证明确为零时钟/零步冷启动，不能据缺字段将续算重置为零。逐条读取避免把所有完整场同时保存在内存。

新细参考与动量相容初值的完整长算位于 `outputs/combustion-foundation/refined-evolution-100us-001/`。`run.py` 绑定原生源码、冻结程序、后端、机理和输入哈希，使用目录内冻结的原审计器；`history-complete/` 使用已验证的 `history-tool-source-002/` 分析完整历史，不能用中途前缀代替终点。`write_complete_artifact.py` 保存实际完整轨迹及与旧初值长算的观察性比较；两次冻结二进制与输入均有区别，不能从中推导单因素因果或通用性能。完整命令、文件哈希、实际控制和逐帧报告见[证据索引](../artifacts/current/native-reacting-refined-evolution.json)。

`refined-evolution-200us-001/` 从上述真实检查点追加 100 微秒，保持冻结程序、机理、网格、边界、三项误差容差和最大步长；只重建 BDF 历史与初始尺度。`analyze_complete.py` 核对启动提交的源码 Git blob、二进制、后端和实际文件哈希，再以原历史工具分析新段；不得因当前工作区已有后续取消修复而将其冒充运行来源。合并诊断保留两段实际审计，严格核对接口场相等，重复的 100 微秒样本只在绘图时间序列中去重；物理预算合并真实边界、化学和代数约束积分。全程温漂和平移形状差异重新相对零时刻计算，本段温漂另列。`moving_frame.py` 使用末段实际中部等温标记速度和非均匀二阶单元中心差分计算 `RHS+c*dU/dx`；这是后处理，不是原生有限体积算子，也不能将两个 L1 范数相减后宣称因果贡献比例。[续算证据](../artifacts/current/native-reacting-refined-continuation.json)

开放边界的脉冲对照在 `acoustic-boundary-002/` 生成相同重叠网格和初始单元守恒均值，仅延长出口域。初始 400 K 氢气/空气混合物叠加 100 Pa 高斯压力脉冲，按固定初始组分的等熵关系及压力积分构造向外速度，以八点 Gauss 积分形成单元平均；原生推进仍保留完整反应、温变物性、黏性、多组分和 Soret。`analyze.py` 对相同物理终点的共同单元直接比较压力、温度、速度和组分，并单列低幅值线性声学方向诊断 `0.5*(Δp ± rho0*c0*Δu)`。`acoustic-boundary-003/` 将间距从 250 减至 125 微米，`acoustic-boundary-004/` 在原网格上将三项时间容差共同减半；每档分别运行原域和延长域，最终比较同一 45 微秒实际场。`acoustic-boundary-004/write_artifact.py` 核对同网格两档的初值完全相同、完整控制及逐帧来源，再比较两条域长压差曲线。域长不同会改变全局加权时间误差范数，控制相同不保证误差相同；不预设边界反射合格阈值。真实模型、控制和结果见[边界对照证据](../artifacts/current/native-reacting-acoustic-boundary.json)。

标记温度（K）、上游测点（m）、平移对齐温度（K）和末段窗口（s）均须显式指定。等温交点与等温平台全部保留；窗口内不足三个样本，或任一时刻没有唯一交点，就不估计该标记速度。用真实保存时刻拟合固定坐标速度，同时列出端点割线、相邻时刻速度范围与拟合位置偏差。测点速度按实际时间梯形积分，再减去标记速度，得到相对运动诊断；另列测点温度/组分相对入口的偏离和是否位于标记上游，不自动将它判定为未燃气或成熟燃烧速度。平移仅插值初始行均温度的真实重叠区，并记录被排除列数、宽度比例和跨行温差；原场、残差和检查点不变。

`test_reacting_flame_history.py` 用解析平移剖面检查位置、秒制速度、变形及多交点/平台，也检查 JSONL 缺帧、增帧和时钟倒退。几何/速度/温度数值检查允许 `64*epsilon` 乘对应构造尺度；速度还除以实际时间跨度，只用于解析恒等式的运算舍入，不增加 CFD 精度门。测试为秒级，真实历史复核只读既有场。示例：

```sh
build/chemistry-env/bin/python tools/flow/analyze_reacting_flame_history.py \
  --report outputs/combustion-foundation/implicit-long-new/report.json \
  --reference outputs/combustion-foundation/flame-reference-new --fixture-index 2 \
  --output outputs/combustion-foundation/flame-history-analysis-new \
  --marker 600 --marker 1000 --marker 1500 --marker 2000 --align-temperature 1500 \
  --station .005 --station .01 --station .0125 --late-window 2e-5
```

显式、隐式原生验证驱动均支持 `--restart 检查点`；Python 运行/审计入口使用 `--restart-checkpoint 检查点`，一次只选一档网格。先用原生绑定检查完整机理、网格、边界和模型，原样保存输入到新输出目录的 `restart.checkpoint`；监督程序核对源文件运行前后及副本 SHA256。`--duration` 是从检查点起追加的物理时长，`endTime` 为实际请求终点，累计接受步数保留。步长预算、方程调用及收支积分从本段起点重新计数，采样步距也从本段起点计算。读取器核对原样导入的完整场、时间和步数，并分别报告本段/累计接受步数；不能把重复计算或重置为零的时钟当作续算。

Python 监督入口的 `--duration` 继续要求追加时长大于零。显式 `--evaluate-initial` 只求值并审计原始夹具，时间、接受步和推进 RHS 调用均为零，不允许与时长、续算或采样参数混用。原几何、EOS、源项及共享面审计保留，另核对初末场完全一致且不存在 `accepted.checkpoint`；报告 `initial_evaluation_only=true`、`all_initial_evaluations_passed`，同时保持 `all_endpoints_reached=false`。同模式可重审，不能把初值报告重分类为已推进区间。零收支只是未推进的恒等式，不提供积分守恒或火焰精度证据。若只诊断终点检查点读回，仍直接调用原生探针，时长为 `0` 并给出 `--restart`；独立 `audit()` 可检查此输出，但明确不将零新增步标作完成新的物理区间。文件读回保持原值，与重建 BDF 历史后继续积分的误差验证分开。

检查点保存物理状态，BDF 历史及初始积分尺度在新会话重建，因此续算轨迹不要求逐字节等于不中断运行。原诊断上限失败的 464 格案例现已从第 20 步、约 `0.0426812` 微秒续至 `0.1` 微秒；新增 23 步，6 帧审计通过。与不中断同终点的最大温差约 `1.98e-5 K`，合并两段实际通量后的质量/元素/能量预算仍通过原检查；未扩大到长期误差资格。冷启动物理场、步长日志和检查点保持原结果。不同网格/入口/绑定、非法时钟、截断及尾部损坏被拒绝；零时长读回保留全部原值并且不新增接受步。该案例也经过显式入口的导入与推进检查，不将它视为两种积分方法的精度等价证明。

本机可选后端目前要求 Cantera 随附并导出的 SUNDIALS 5 接口；其他 SUNDIALS 主版本、App 和平台打包尚未适配。先使用已有的详细化学配置，再运行：

```sh
cmake -S . -B build -DCARTMESH2D_BUILD_REACTING_IMPLICIT=ON
cmake --build build --target cartmesh2d_reacting_implicit_tests cartmesh2d_reacting_implicit_flame_probe -j2
ctest --test-dir build -R '^cartmesh2d_reacting_implicit(_trace)?$' --output-on-failure
```

`compare_reacting_flame_grid.py` 比较同机理、同终点及同积分控制的平面条带网格。细格守恒量按真实矩形重叠面积积分到粗格，再使用原机理热力学反演温度；不直接平均温度或将不同位置的数组相减。初始映射差异、终点场差异及两者相减得到的演化差异分别报告。重叠权重的常数保持和体积守恒检查仅允许 `64*epsilon*参与格数` 的浮点运算舍入，属于精确几何恒等式诊断，不是新增 CFD 精度门；另有非嵌套矩形分段常数解析核对。非嵌套投影本身有离散误差，因此不把这些差异标作正式空间阶数。每档还独立积分实际燃料消耗，列出组分残差和温度漂移。示例：

### 空间火焰对照与微量组分回归

`prepare_reacting_flame.py` 使用独立的 Cantera `FreeFlame` 稳态边值求解器，保留完整多组分和 Soret，保存三档参考细化的原始剖面与机理。默认验证夹具为 `H2:1.6,O2:1,N2:3.76`、400 K、1 atm；网格在反应层附近取 40/20/10 微米，并保留完整上游和下游。该网格布局按本氢气案例设置，其他燃料须重新确认解析度，参数接口不等于已验证。

参考解的舍入缺陷显式处理：保存原始负微量及组分和误差，只在验证初值准备阶段将负舍入值置零，并由最丰富组分闭合。允许改变量为 `10*BVP绝对容差 + 64*(N+1)*epsilon`，超过则拒绝；有意输入 `-1e-4` 会失败。保留原参考点的密度、动量及含生成能总能量，记录元素变化，再对分段线性守恒剖面逐格积分。该准备程序不进入原生推进或续算，初值不是接受检查点。

参考解可用 `--refine-slopes` 增加细化层级，`--maximum-reference-points` 只控制参考网格的计算预算。前者是 Cantera 按剖面变化增点的无量纲控制，曲率控制为其两倍，不是燃烧精度验收门；BVP 原求解容差保留。`--fixed-grid-reference` 与 `--fixed-grid-index` 从哈希验证过的旧夹具复用精确网格边、行数和高度，重新积分最终参考守恒剖面，用来隔离参考初值敏感性。拒绝坏网格、哈希不符及超出参考域的映射。`reference.json` 的 `profile_file` 明确指向最终剖面；原生审计器及动量初值准备器均读取它，旧参考仍兼容。

`cartmesh2d_reacting_flame_reference_probe` 是只链接 Cantera 的独立稳态 BVP 验证程序，不链接原生瞬态有限体积核心。可选 `--closure-probe` 在固定参考网格上重新求解完整机理、多组分和 Soret，用 `sum(Y)=1` 替代一种丰富组分在内部点的冗余方程；原组分方程全部照常计算，原残差另存于最终剖面，不能因质量恒等式满足就忽略组分平衡。能量、边界和原求解容差保留，最终初值导入仍使用上述原额度。原未约束参考及失败结果保留；该选项不归一化求解后的参考剖面，不修改原生接受状态或默认设置。示例：

```sh
cmake --build build --target cartmesh2d_reacting_flame_reference_probe -j 2
build/chemistry-env/bin/python tools/flow/prepare_reacting_flame.py \
  --mechanism build/deps/cantera/share/cantera/data/h2o2.yaml \
  --output outputs/combustion-foundation/refined-reference-new \
  --refine-slopes .06 .03 .015 .0075 .00375 .001875 \
  --maximum-reference-points 5000 \
  --fixed-grid-reference outputs/combustion-foundation/flame-reference-001 \
  --fixed-grid-index 2 --closure-probe build/cartmesh2d_reacting_flame_reference_probe
```

参考细化同时更新自由火焰的入口速度特征值；同网格对照须显式列出入口状态变化，不能称边界数值完全相同。实际对照、逐帧审计及未验范围见[当前状态](CURRENT_STATE_CN.md#持续目标成熟燃烧模拟)，脚本和原始失败保存在 `outputs/combustion-foundation/reference-refinement-001/`。

参考探针还直接导出完整稳态求值后保存的区间扩散质量通量（kg/m²/s）、内部点化学质量源（kg/m³/s）及实际入口条件。区间行 `j` 位于节点 `j` 与 `j+1` 之间，正号沿 x 增大方向；导出不重新计算或改变输运状态。`reference-flux-001/run.py` 用原五份输入重算，逐项核对旧输出字段完全相同，再读取这些新增诊断。

元素收支必须同时包含对流和扩散，不能只检查 `sum(Y)` 或非线性残差。`reference-flux-001/analyze.py` 根据机理的原子质量矩阵，将实际边界组分通量投影到各元素；偏差为“出口减入口”，单位 kg/m²/s，相对量除以该元素的规定入口通量。入口不存在的元素只报绝对量，相对量记为 null。正向流的上风组分导数使用前一间距，而扩散散度使用左右间距均值；程序用后者积分全部原组分方程，并显式列出“对流边界差减去上风导数积分”的离散项。由此逐项重构：边界通量差＝化学源积分−原方程残差积分＋对流离散项。运算重构误差及其浮点尺度另报，不新增或放宽物理精度门。

`reference-domain-002/` 保留原核心节点和温度锚点，分别延长上下游，并将新增下游间距减半；`reference-flux-001/refine.py` 则将原域的每段间距连续二分两次，保留所有旧节点及锚点，插值仅生成初猜，再求解原完整 BVP。两者分开识别域截断和网格误差，均不改变原生长算输入。原始负舍入组分照实保存，含负出口组分时跳过出口平衡计算，不能通过投影获得平衡诊断。冻结程序、脚本、参考源文件及失败启动均由[参考守恒证据](../artifacts/current/native-reacting-reference-flux.json)索引；复现须使用新输出目录。这些检查不替代原生可压边界、网格或实验对标。

显式选项 `--conservative-reference` 要求同时提供 `--closure-probe`。它在相同固定网格上求解共享区间的**组分总通量和含生成焓总焓通量散度**，自由火焰连续性、温度锚点、完整化学、多组分/Soret、代数质量约束及原 BVP 容差保留。区间质量流量取两节点 `rho*u` 的平均，对流组分及焓取上风节点；总焓通量还包含 Fourier 导热和扩散组分携带的焓，后者取两节点组分焓均值。内部控制体宽度统一为左右间距之和的一半。入口组分通量与规定组分质量流量匹配，原温度及出口梯度条件保留。该稳态方程中的对角伪时间项仅用于非线性求解稳定化，不代表物理瞬态能量积分；未另加会重复计算生成焓的热源。原参考形式仍是默认，原生瞬态有限体积方程不变。参考物理与原离散来源见 [Cantera 方程](https://cantera.org/3.2/reference/onedim/governing-equations.html)及[离散说明](https://cantera.org/3.2/reference/onedim/discretization.html)。

在上述准备命令末尾增加 `--conservative-reference` 即可选择新形式，最终剖面记为 `profile-conservative.json`。参考探针的直接调用为 `cartmesh2d_reacting_flame_reference_probe <机理> <初猜> <新输出> 1 1`，最后一项省略或为 `0` 使用原形式。剖面同时保留原方程和所选方程残差、区间总组分/总焓通量；温度方程只将内部点残差标为 K/s，边界方程不混入该单位。原始负舍入值仍保存，导入额度与原生正性门保持。

`reference-conservative-001/` 保存同初猜的新旧参考对照、全域间距二分及实际准备工作流；`check_promotion.py` 核对正式程序与冻结原型的物理场和诊断逐项相同，诊断字段改名及温度边界行剔除显式映射。`audit_flux.py` 从真实节点 T/Y、扩散通量及化学源，用显式多组分 Python Cantera 独立重算热力学和区间导热系数，再重建总通量与所选方程。算术误差按相消前各项绝对值之和归一化，额度为 `512*(组分数+1)*epsilon`，并通过同一控制体宽度、密度和比热传播到残差；仅检查 C++/Python 浮点分组的一致性，不是物理误差门。成本为每节点一次热力学、每区间一次导热查询和线性规模运算，无额外 BVP 求解。首次读取器误选混合平均输运、其次使用相消后通量作归一化的两份失败源码和日志均保留；修正读取器后，篡改真实组分通量、焓通量、方程残差、温度或输运标识均被拒绝。[守恒参考的原始数据与验证](../artifacts/current/native-reacting-conservative-reference.json)

`reference-conservative-evolution-001/run.py` 从已审计的 3446 点守恒剖面复用正式初值舍入处理与保守积分函数，映射到原 890 格几何，经零时间原生读回后生成动量相容输入。短时演化使用与旧参考对照相同的冻结原生程序、完整物理模型和积分控制；`analyze.py` 核对实际检查点绑定、网格与逐帧哈希，只允许储库数值随参考入口特征值变化，其他边界类型与物理开关须相同。两组采用同一标记温度，其值明确保存在结果中；不同历史分析若选取不同标记，其平移后形状量不能直接混用。原 200 微秒检查点和已完成结果保留；参考场变化不作为对接受状态的修正。[短时演化来源](../artifacts/current/native-reacting-conservative-evolution.json)

对运动火焰，分别记录固定坐标中的等温面速度和相对于实际未燃气速度的运动，不能把输入的入口速度当成实际面质量流量。当前 `Reservoir` 是外侧黎曼状态，亚声速面通量还受内部压力影响；[Cantera 的参考方程与入口条件](https://www.cantera.org/3.2/reference/onedim/governing-equations.html)则将质量流量作为自由火焰解的一部分。`flame-motion-001/analyze.py` 从原 100 微秒真实历史中计算多个固定温度标记、三个明确位置的未燃气测点和实际时间窗口。它还分解 `d(rho*Y)/dt = Y*d(rho)/dt + rho*dY/dt`；各项 L1 范数因相消不能相加为贡献百分比。移动参考系的 `R + c*dU/dx` 用独立非均匀差分作诊断，原生残差与验收门保留。

`prepare_reacting_flame_momentum.py` 是显式可选的**新初值**准备器：从通过审计的非续算夹具读取原始单元 T/Y，保留入口状态及网格，以原入口轴向质量通量为常量，在低马赫分支求解 `p + m²/rho - tau_xx = 入口动量通量`，其中 `rho = p/(Rmix*T)`、`tau_xx = (4/3)*mu*du/dx`，黏度由完整温变物性取得。沿单元中心用非均匀二阶差分构造速度导数，固定点循环后重建含生成能的守恒变量；不求解稳态化学/组分/能量方程，也不声称原生离散残差为零。循环检查最大压力更新量（Pa），停止额度为 `64*epsilon*max(|p|)`，用于压力构造的浮点终止，最多 50 次；不增加燃烧精度门。源文件哈希、初始质量/能量改变量和循环记录输出到 `preparation.json`。禁止覆盖既有目录，拒绝续算输入、非轴向条带及来源不符的数据。

准备器既可读取已通过审计的正时长案例的初始场，也可读取显式零时长初值审计。后者复核独立审计通过、零时钟/接受步、初末场不变及没有接受检查点，并将来源标为 `audited_initial_evaluation`；无需先推进一段无关时间来取得初值物性。既有来源类型标为 `audited_evolution_initial_state`，原字段和输入仍兼容。

两种输入的真实计算保存在 `outputs/combustion-foundation/flame-motion-001/`；比较器检查二进制、完整检查点绑定（机理/网格/边界/物理开关）、实际积分控制和逐帧审计，分别列出坐标位移与平移后的形状差异。正式准备器的夹具另与实际对照输入作逐字节核对，原始输入保留。数值和未验范围统一见[当前状态](CURRENT_STATE_CN.md#持续目标成熟燃烧模拟)，来源见[诊断证据](../artifacts/current/native-reacting-flame-motion.json)。从原夹具开始的复现示例：

```sh
build/chemistry-env/bin/python tools/flow/prepare_reacting_flame_momentum.py \
  --reference outputs/combustion-foundation/flame-reference-001 \
  --source-run outputs/combustion-foundation/initial-evaluation-new --grid 2 \
  --output outputs/combustion-foundation/momentum-reference-new
```

`outputs/combustion-foundation/refined-state-refinement-001/` 保存新细参考的时间/空间对照。`map_reference.py` 从已哈希验证的 3446 点 `prepared-reference.npz` 积分到原三档精确网格，使用正式准备器的保守映射函数，不重算 BVP；`run_grid.py` 经三档初值审计后生成动量相容输入，粗/中档实际推进，细档夹具与既有完整运行逐字节相同后复用。因合并夹具清单与旧单档清单的元数据不同，`analyze_refinement.py` 分别核对各自清单、同一 BVP/映射源、入口行、检查点机理/物理前缀及完整边界数值和控制，再调用原网格比较器的守恒重叠/EOS 函数；没有放宽通用比较器要求相同清单的检查。`run_time.py` 使用冻结的上一版审计器在相同细档上收紧时间容差，两份版本的原物理审计函数 AST 与正时长审计结果均核对一致。完整命令、源码/原始场哈希和初值模式拒绝检查见[证据索引](../artifacts/current/native-reacting-refined-state.json)；复现使用新目录，不覆盖已完成结果。

`cartmesh2d_reacting_flame_probe` 建立真实二维条带拓扑，左侧固定储库、右侧外推流出，上下滑移绝热。输出初始/最终守恒场、实际共享面通量、分开的化学与输运残差、开放边界积分、逐步接受/拒绝、展开机理与最后接受检查点。`evaluateResidual` 是瞬时半离散算子诊断，不会把输入标为已收敛或接受步。`verify_reacting_flame.py` 独立读取实际文件，检查几何关联、质量/元素/动量/能量预算、EOS、未被修改的导入初值与真实时钟，并绘制单元场和残差趋势。数值检查、达到物理终点和火焰资格分别记录。

隐式驱动可用 `--species-atol` 显式比较组分局部积分权重，实际取值写入结果；质量/元素/能量等物理审计判据不随该选项改变。

`--grid 0` 可只运行首档，重复 `--grid` 选择多个参考夹具；报告明确记录所选索引，不能将单档报告称为三档验收。

`tests/fixtures/reacting_flame_trace_muscl.fixture` 和 `reacting_flame_trace_soret.fixture` 是从首档剖面分别截取的 20/36 格最小失败片段，保留原始微量组分。CTest 要求在不拒绝候选的情况下推进到 `1e-8 s`，并沿用装配和元素预算检查，防止靠大量重试掩盖旧错误；它们只验收数值修复。

细档长序列曾在第 1209 个接受步后因组分和的舍入误差达到 API 边界而失败。`reacting_flame_trace_closure.fixture` 保留原检查点中的 6 格片段，旧版本可直接复现。`SpeciesMassClosure.hpp` 在化学/输运阶段及面状态插值中显式使用质量约束：先检查全部原始组分非负、有限，并确认质量缺陷除以总质量不超过原有 `64*(N+1)*epsilon` 求和额度，再从总质量减去其余组分的补偿求和，确定最丰富组分。不缩放其他组分、不设微量下限、不改变密度或能量；不合法输入和超过舍入范围的缺陷继续拒绝。完整反应机理及全部组分方程仍参与计算，这只是冗余质量约束的数值闭合。

`ChemistryStep::massClosureChange` 和耦合步的化学/输运闭合向量、最大单元相对改变量、绝对改变量积分均显式输出。SSPRK2 第一阶段改变量按半权重计入。离散装配残差扣除已记录的输运舍入操作，但实际质量/元素/能量预算**不扣除**它们，避免用闭合制造物理守恒证据。火焰读取器逐步核对诊断之和；2048 次真实恒容化学重启、负微量与超额缺陷拒绝也纳入原生检查。末态原始守恒值先保存到 `accepted-state.json`，再求末态残差，诊断失败仍保留最后接受状态。

`transport(A) → viscosity(B) → transport(A)` 的最小复现曾使同一个 A 状态的导热系数与 Soret 系数不同。原因是黏性查询使二元扩散缓存失效，而旧热缓存仍保留 A 的温度键；重新建立二元矩阵时覆盖了热输运使用的自扩散对角项。`TraceStableMultiTransport::update_T` 现在在每次温度改变时同步使热缓存失效。该查询序列逐值一致性和有/无显式步长估计的真实物理残差一致性都进入原生回归；48 项局部输运对照重新通过。旧性能记录仍保留，不能据此跳过物性正确性检查。

完整物性查询将多组分矩阵查询放在导热/Soret 查询之前，复用同一热输运矩阵求解；黏性应力用 `DetailedGas::viscosity` 单独查询。已移除的历史 `compare_reacting_cost.py` 针对这类缓存/查询优化交替执行旧、新完整 48 格点火计算，计时含冷启动至全部输出，要求数值场逐值一致、步日志/机理/检查点逐字节一致，只排除墙钟字段。一次三对实测中位数约 `8.274→7.927 s`；这是引入质量闭合之前的隔离性能比较，不能作为当前所有改动或任意火焰规模的提速结论。

```sh
cmake --build build --target cartmesh2d_reacting_flame_probe -j2
ctest --test-dir build -R '^cartmesh2d_reacting_flame_(muscl|soret|closure)$' --output-on-failure
build/chemistry-env/bin/python tools/flow/prepare_reacting_flame.py \
  --mechanism build/deps/cantera/share/cantera/data/h2o2.yaml \
  --output outputs/combustion-foundation/flame-reference-new
```

参考 BVP 采用恒定热力学压力，原生平面可压方程还保留纵向动能与黏性功，因此不是完全相同的方程组。从参考剖面初始化后的短时漂移、反应积分或初始残差下降均不能当作原生独立预测的火焰速度；原生稳态传播、时间/网格误差和实验还需另验。上述驱动尚非用户可配置的通用燃烧产品。湍流燃烧、喷雾、辐射、共轭传热、实际燃烧室、App、平台资格和大规模成本均未完成；完整目标和最新证据统一见[当前状态](CURRENT_STATE_CN.md#持续目标成熟燃烧模拟)。

## 流体拓扑优化研究入口

`tools/optimization/` 是独立命令行研究工具，不替换产品求解器，也不读取背景 JSON 作流体格。可选 Python 依赖不进入桌面运行包；研究测试独立于 CTest，由 `flow-research.yml` 在三平台执行。

### 模型与停止条件

`brinkman.py` 保留二维 MAC Stokes–Brinkman 线性路径；`navier_stokes_brinkman.py` 增加 `γ div(u⊗u)−μ∇²u+∇p+α(ρ)u=0`、`div(u)=0` 的稳态不可压惯性分析。速度以入口峰值为 1，长度使用设计坐标，`Re=U_mean·w/ν`（入口端口宽度 w，抛物线均速 `U_mean=2/3`），故 `γ=Re·μ/((2/3)w)`。`--reynolds 0` 是精确 Stokes 极限；正 Re 使用阻尼 Newton，困难时从 Stokes 自适应延拓，失败不覆盖最近收敛根。ρ=1 为流体、0 为有限阻力固体，`α=αmax·q(1−ρ)/(q+ρ)`。抛物线端口积分匹配，两格端口固定流体、外壁固定固体；锥形滤波及 tanh 投影后检查真实 `mean(ρ)≤volume_fraction`，滤波半径不是已证明的最小制造壁厚。

对流为守恒 MAC 对偶体中心通量：两个速度分量均用算术平均定位到对偶面，法向动量通量为平均速度的平方，交叉通量为两平均速度的乘积，零法向壁面通量。依据 [Harlow–Welch 1965](https://doi.org/10.1063/1.1761178) 的交错离散；规则内部二阶，不加人工黏性。解析微分所有二次单项式，伴随使用收敛点完整 Jacobian 的转置，仍贯穿材料、滤波和投影链。高单元 Re 下中心格式的空间精度/有界性须另验。

Newton 残差为已消元的 `[积分动量;−积分连续性]` 无穷范数除以 `max(||Dirichlet rhs||∞,1e-30)`，默认停止门 `1e-11`；另保留 `max|Du|/Q_in≤1e-8` 连续性及 `1e-8` 伴随门。这些是代数一致性门，收紧 Newton 的任务依据是避免求解误差污染双精度中心差分，成本为每分析数次稀疏分解，不作为物理误差标准。梯度验收为 `|FD−adjoint|/max(|FD|,|adjoint|,1e-12)` 的步长曲线最低值 `≤1e-5`；13 个 h 从 `1e-1` 到 `1e-9`，每 Re 3 个随机设计，两目标、β=0/3，共 60 条曲线（含可收敛的 Re=200）。

```sh
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 outputs/topology-env/bin/python tests/navier_stokes_topology_test.py
```

实际用于B–D的总压功另做30条同规则曲线，原60条不覆盖；合计90条，最终梯度图同时展示三个目标。总压功也有30/30个内部最低点，最差最低误差2.73e-8。

消去指定速度和一个压力参考后，用离散矩阵转置求伴随，贯穿滤波/投影导数。默认目标为黏性与阻力耗散；`pressure-power` 是端口邻格压力的通量加权差，不是精确边界应力功。hybrid 更新先试 OC，非下降时尝试有移动界的梯度步并二分校正实际体积；接受前重新求解、检查方向与回溯，不是精确非线性投影。

三阶段 `(q,β)=(.01,0),(.1,2),(.1,6)` 改变模型，跨阶段目标不可直接比较。`projectedKkt` 是归一化梯度在体积切平面/变量界上的投影步无穷范数，默认 1e−3 只用于局部设计停止；线性/伴随 1e−8 检查归一化方程残差，连续性按总入口流量归一化，均不代表物理精度。预算耗尽、停滞和分析失败明确记录，保存最后接受设计与阶段参数，目前无优化自动续算。

### 运行与真实壁面连接

以下使用 POSIX venv 路径；Windows 使用其 `Scripts/python.exe` 并设置同名线程环境变量。

```sh
python3 -m venv outputs/topology-env
outputs/topology-env/bin/python -m pip install -r tools/optimization/requirements.txt
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 outputs/topology-env/bin/python tests/flow_topology_test.py
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 outputs/topology-env/bin/python tools/optimization/optimize_flow.py --case double-pipe --output outputs/topology/double-pipe
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 outputs/topology-env/bin/python tools/optimization/native_flow.py outputs/topology/double-pipe --output outputs/topology/double-pipe-native --levels 5 6 --tolerance 1e-8
```

使用新输出目录。默认 48×32、25/40/60 步预算；`--initialization uniform|random|geometric|merged` 区分无内部路径、随机材料、参考几何和人为合流种子，`--seed` 固定随机材料初值。`--no-plot` 也跳过轮廓提取，可再运行 `topology_artifacts.py <目录>`。summary/history/NPZ/阶段快照及分析 VTK 保存控制、源码哈希、停止原因和设计，分析场仍含固体阻力。

轮廓提取保留所有连通分量与孔洞，不平滑/填孔；`fluid.xy` 显式 interior，各独立域另存文件。轮廓面积与密度体积分别报告，连通矩阵只表示几何可达，不等于混合效率。连接脚本真实执行 Cut-cell、质量门、独立面积读取、原生低 Re 流动和独立方程审计；所有分量均须成功，有孤立区/盲腔则明确拒绝。仅全局端口覆盖面可作入口/出口，不能按组件局部包围盒制造端口。

未指定 `--reynolds` 的历史原生对照显式用速度尺度 .02、ν=1（w=1/6时峰值端口Re约.00333、平均值约.00222）、抛物线入口和压力出口；与T01匹配模式分开。1e−8是这些例子的代数控制，不改产品默认或审核门。检测到checkMesh时实际运行，否则记not-run；多孔目标下降不能直接当真实壁面性能。`render_native_flow.py` 从接受的多边形及CSV作图并保留源哈希。

### 公平对照、失败与继续迭代

T01 新入口为 `fidelity_study.py b --output outputs/topology-fidelity/b-rerun`（使用上述双线程环境变量）。`topology_cases.py` 固定文献来源和参数映射；`engineering_baselines.py` 对每例 3 个预先定义的工程形状按相同过滤/投影体积求宽度，再在相同 Re、端口流量、物性与 α 下全部计算并取目标最低者作为主要基准，所有备选及失败保留。双管采用合流管族，弯管采用三种转弯半径，扩压器采用三种扩张率，四端口比较 U 形回路与平行通路。旧直管行仅保留为历史/附加参照。

| 案例 | 本轮设计坐标与文献映射 | 对标边界 |
| --- | --- | --- |
| 双管 | H=1、L=1.5、w=1/6、端口中心 1/4 和 3/4、体积分数 1/3，36×24 | [Borrvall–Petersson 2003](https://doi.org/10.1002/fld.426)，参数由 [dolfin-adjoint 实现](https://www.dolfin-adjoint.org/en/stable/documentation/stokes-topology/stokes-topology.html)核对；FEM/滤波/端口固定位不同，不作逐数复现 |
| 扩压器 | 单位正方形，w_in=1/3、w_out=1，出口峰值 1/3，体积分数 1/2，24×24 | 源于同文献案例类型；原完整参数获取失败，明确为自定改编，不能宣称原文数值/形状复现 |
| 90° 弯管 | 单位正方形，左入口 y=.8、下出口 x=.8、w=.2、体积分数 .25，30×30 | B–P 几何经 [Gersborg-Hansen 论文集 Fig.C.1](https://backend.orbit.dtu.dk/ws/portalfiles/portal/4897055/C_Agh_PhD_Texts_Thesis_070122b_final_070410_final_Forside_AGH_thesis.pdf)坐标除以 5；[2005 惯性研究](https://doi.org/10.1007/s00158-004-0508-7)的流体区底阻力不同 |
| 四端口 | H=1、L=.7、w=.2、y=.3/.7、对角入口，体积分数 .4，21×30 | [Olesen 等 2006 Fig.6–7](https://arxiv.org/pdf/physics/0410086)：原 L=3.5ℓ、H=5ℓ。本轮省略 2ℓ 引导管并规定出口速度；原峰值 Re=20/200 对应本均值 Re=13.33/133.33。原文低 Re 的 U 回路与高 Re 平行通路只作拓扑示意，不移植临界 Re 或目标数字 |

图中的文献栏是根据连接关系自绘的示意图，非原图数字化。原优化流道的同控制定量复现仍未完成，不能把不可比的数值凑成相对误差；本次另有精确直管参考的定量对照，范围见下面的补齐方法。四案例首轮真实不平滑轮廓、工程参照和停止原因见 `artifacts/current/native-topology-fidelity-benchmarks.png` 与同前缀JSON，补齐结果见 `native-topology-fidelity-completion.json`。

新增 `total-pressure-power` 目标使用端口边界线性外推 `p_b=1.5p_first−.5p_second` 的压力功，加解析规定抛物线速度的进出动能通量 `γ·8wUpeak³/35`；后者与设计无关，伴随仍须使用新的边界压力权重。它是边界总压功定义，不假定有限网格上等于离散体耗散。旧 `pressure-power` 邻格压力目标和默认耗散目标保持。高 α 时直接 LU 转置解可能超出原伴随门，新增至多 5 次实测残差修正，门仍为 `1e-8`；原失败三运行均保留，新回归可复现该参考设计。

### T01：匹配 Re、保真度差距与排序

`fidelity_study.py c` 对四案例 × Re=0/10/50/100 × 三档阻力逐一运行，`Da_w=μ/(αmax w²)` 取 `1e-3/1e-4/1e-5`；按设计高度定义的 `Da_H=Da_w(w/H)²`。以 w=1/6 为例，αmax 为 `3.6e4/3.6e5/3.6e6`；w=.2 为 `2.5e4/2.5e5/2.5e6`，w=1/3 为 `9e3/9e4/9e5`。这些是预设研究采样，不是阻力充分大的验收门；每案例固定文献映射的体积分数，本轮没有增加第二体积分数。

`native_flow.py --reynolds Re` 和 `compare_sharp_designs.py --reynolds Re` 启用匹配模式。原生峰值速度 S=.02、正 Re 时 `ν=(2/3)Sw/Re`，动量采用与 MAC 对应的 Laplacian 黏性形式、显式 limited-linear 对流；入口和出口都积分同一抛物线剖面。Re=0 使用 ν=1、显式 `--momentum-inertia 0` 的精确 Stokes 方程，不能用非零流量下的近零 Re 冒充。新增 `velocity-outlet` 保留出口语义、要求速度向外；全规定流量边界须在既有 TolerancePolicy 下守恒，并只固定压力规范。惯性默认仍为1；0仅用于稳态 custom，拒绝瞬态、物性耦合与物理检查点。

原生运动学总压功为 `P=Σ_in Q p_b−Σ_out Q p_b + Σ_in∫|u|²|u·n|/2 ds−Σ_out∫|u|²|u·n|/2 ds`，压力来自接受场实际边界面 CSV，动能采用规定抛物线的解析积分；Stokes 动能系数为0。比较量 `J_T=P/(νS²/μ)` 与峰值1的 MAC `J_B` 同一无量纲尺度，P 是单位密度、单位深度的量，不直接称瓦数。原生与 MAC 边界压力重建不同，有限网格的总压功不强行等同体耗散。

轮廓不平滑，保留全部分量/孔洞；通过内部密度的单一偏移匹配同一真实面积预算，原密度解另存。面积投影可改变几何，故 gap=`(J_B−J_T)/J_T` 包括多孔模型、固定 MAC 分辨率和提取投影的共同影响，不能单独解释为 Brinkman 建模误差。固体漏流统计为内部 MAC 面在平均ρ<.5处的绝对体积通量之和/全部内部面绝对通量之和；重复穿越重复计数，它不是入口粒子中漏入固体的比例。固体耗散占比用正的黏性弹簧能和阻力能分配到ρ<.5单元，并检查重构总和。

每个几何首先运行 level4/5/6，所有分量必须通过原质量门、严格1e-8停止门和独立方程审计后才有该档 J_T。每进程最多6000步/90秒，预算到达不算收敛。旧 duct 模板曾误拒绝没有全局右出口的弯管/U回路，匹配模式现直接从实际 CM2D 几何生成全壁模板，再赋全局端口；保留原运行，在 `c-recovery/` 重跑全部受影响档位。C 对不足两档的案例统一补 level7、再 level3；D 为控制磁盘统一先补 level3、再 level7。每档取首个接受结果，所有重试仍在 `nativeAttempts`，最细有效档给 J_T；失败细档不会因存在粗档而被称为通过。

三档单调结果采用 [Roache/Richardson GCI 的 NASA 说明](https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html)：有效 `h=sqrt(A/N)`，按非等比细化的 Richardson 方程求正观测阶 p，`U_GCI=1.25|J_f−J_m|/(r_fm^p−1)`。归一化 GCI 为 U/|J_f|；三点本身不能独立证实渐近幂律。两档、振荡/零差或无正阶根仅报告末两档变化，不补 GCI。gap 区间用 `J_B/(J_T±U)−1` 非线性传播；若下界非正则不能给有限上界。这个条件估计不代表 ASME V&V 20 认证或物理精度证明。

D 每案例 Re=0/50 复用C三档Da，另加预先固定的五个滤波/初值候选；只有密度字段重复而不足8个不同设计时才按固定顺序追加，不能按目标或反转结果挑候选。初值、种子、参数和所有两两ρ最大差入表；1e-6无量纲密度仅用于重复诊断，排除舍入差，不作几何或流动精度门。Kendall τ 同时报告全部可用 J_T、至少两档和有 GCI 的子集。列出全部两两比较，仅当两侧有 GCI 且真实目标差绝对值大于 GCI 之和时称“在该网格估计下可分辨的反转”；两档变化另列，不能顶替不确定度。

```sh
cmake --build build --target cartmesh2d_flow_cli cartmesh2d_flow_boundary_tests cartmesh2d_flow_checkpoint_tests --parallel 2
ctest --test-dir build -R '^cartmesh2d_flow_(boundary|checkpoint)$' --output-on-failure
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 outputs/topology-env/bin/python tests/topology_fidelity_test.py
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 outputs/topology-env/bin/python tools/optimization/fidelity_study.py c --output outputs/topology-fidelity/c-rerun
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 outputs/topology-env/bin/python tools/optimization/fidelity_study.py c-recovery --source outputs/topology-fidelity/c-rerun/summary.json --output outputs/topology-fidelity/c-recovery-rerun
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 outputs/topology-env/bin/python tools/optimization/fidelity_study.py d --source outputs/topology-fidelity/c-recovery-rerun/summary.json --output outputs/topology-fidelity/d-rerun
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 outputs/topology-env/bin/python tools/optimization/fidelity_study.py d-uncertainty --source outputs/topology-fidelity/d-rerun/summary.json --output outputs/topology-fidelity/d-uncertainty-rerun
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 outputs/topology-env/bin/python tools/optimization/fidelity_report.py
```

候选生成后，对全部恰好两档有效且未试过level3的D候选统一补算level3，争取第三档；不根据目标大小、观测阶或反转情况选择。`fidelity_report.py` 重算全部排名、两两比较和图；也检查实际提取双管是否确为矩形，再对照独立平面Poiseuille解析功率，不用解析值替补失败曲线。完整CFD重跑需要新的资源预算；已有数据的统计/作图复算无需重做CFD。

T01 驱动不再设置独立磁盘阈值或逐候选扫描全部历史结果；文件系统错误直接上报。每次实验使用新的输出目录，保留原始场和失败记录。汇总脚本的 `--matrix/--ranking/--artifact` 可指向新结果；历史结果继续保存在 `outputs/topology-fidelity/`。

### T01 独立审查与研究选择

`outputs/topology-fidelity/independent-review-001/` 保存本轮审查驱动、重新求解的状态和新原生失败。`review.py` 从原始面CSV/边界文件重建端口功率，对动能项另用四点Gauss积分；不调用原报告的 `pressure_metrics` 或 `ranking`。它重新运行现有MAC状态算子，用另写的边界功函数核对66个原预测，再把每个保存的物理ρ场放到共同 `Da_w=1e-3/1e-4/1e-5` 下评分。共同参数下不重新优化、不再次滤波/投影；模型、端口、物性和目标在每个案例/Re分组内相同。这是对原分析器的独立复算与统计审查，不是独立CFD求解器对照。

原18对中唯一同αmax的一对在复算前选定为网格敏感性对象。每个原MAC单元均匀拆成1/2/4倍每方向的子单元，子单元保持原ρ和q，因此原分片常数阻力场、面积平均ρ及端口剖面保持；不重施随格数变化的被动区域掩码。该探针检查指定材料场的离散敏感性，不代替连续设计描述下的完整空间收敛。全部新状态另存，原 `final.npz`、XY、网格及接受场不覆盖。

原生追加检查用variant-3原封存轮廓，level7、`small-alpha=.25`、`Re=50`、strict、`tolerance=1e-8`、松弛.2及Anderson，其他控制沿用原桥接。两分量分别保留；分量0初次6000轮失败后仅追加2000轮同网格迭代，最终仍失败，不能拼出完整方案目标。分量1的旧3/4/6档单独拟合与新level7接受值比较，作为对旧误差估计的留出网格检查；不把一个分量的结果当全部流道的精度。`finalize.py` 另用二分法复核36个GCI的幂律方程，并从实际数据生成[审查图和索引](../artifacts/current/native-topology-fidelity-review.json)。旧汇总及四图保留为原实验快照，解读以当前状态中的审查结论为准。

审查没有增加、放宽或收紧求解/质量验收门。共同阻力档来自原T01预设采样；GCI与末两档变化分开，后者不是置信区间。封存驱动引用清理前的独立验证模块；复现这条历史审查链须恢复对应版本源码，并使用新输出目录。当前开发入口不直接运行这些封存脚本。追加原生运行的完整命令、二进制/网格哈希和失败场在对应 `summary.json`；相关测试日志、源码哈希和输入哈希由精简索引串联。

### T01 补齐方法与复现

`outputs/topology-fidelity/t01-completion-001/` 保存全部补算驱动、事先记录的任务表、日志及新场。旧B/C/D账本及四图不覆盖；当前完整索引为 `artifacts/current/native-topology-fidelity-completion.json`，指向本目录完整 `summary.json`。以下是封存版统计/作图入口，依赖当时的源码和独立验证模块；恢复对应版本后才可执行：

```sh
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 outputs/topology-env/bin/python tests/topology_fidelity_test.py
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 outputs/topology-env/bin/python outputs/topology-fidelity/t01-completion-001/finalize.py
```

`match_sharp_area` 的失败例保存在 `tests/fixtures/topology_area_plateau.npz`：两个真实材料场在搜索中间偏移`.5`时会出现零面积等值线，原输入与最终等面积轮廓仍然有效。搜索只跳过这一明确失败的中间探针，在有效探针之间继续二分；其他错误仍失败。最终面积仍按原逐例额度检查，固定端口材料不变，不删除任何分量、孔洞或短边。回归同时检查失败探针仍被拒绝、输入未修改及最终面积。

B的两个停滞场先按原更新方法在最终 `(q,β)=(.1,6)` 下追加1000步/300秒，失败记录不改写；随后 `finish_stationarity.py` 用SciPy SLSQP和原解析伴随梯度、真实投影体积约束及变量界继续。内部ftol只控制SLSQP，最终仍独立以原projectedKkt `≤1e-3`、原体积额度及状态/伴随门验收；SLSQP返回success不授予收敛。保存最后可行下降设计，失败候选不会覆盖它。新B字段另存，C/D仍指向原材料场。`benchmark-slsqp.json`、`benchmark-continuations.json`记录全部步骤和原始来源哈希；当前两个场可由 `finalize.py` 重新求解核对，并与原同参数工程参照比较。

困难扩压器的level8先按原控制从冷启动运行，保留超时及目标网格；`recover_flow.py` 将同几何、同分量、同物性的已接受粗场按单元质心最近邻映射作初始猜测。复制目标网格/边界前核对哈希，最终8000步预算、strict、`tolerance=1e-8`及独立方程审计不变；它不是物理检查点或精度提升证据。留边变体使用桥接新入口 `--padding-fraction`，默认仍为旧`1/30`，传给已有网格CLI参数；`.25`为该CLI已有默认。全程原始XY、Re、流量、面积及质量门不变。`nativeGridFamilies`按留边和聚合设置分别列档；主汇总保留原顺序首个有效结果，不按目标或GCI挑档，混合网格族不拟合GCI。每条最终至少有一个同设置的两档序列，几何/质量/流动失败均保留。

当前桥接 `native_flow.py` 使用原生网格 CLI 返回状态和原生求解器的收敛状态；不再运行独立 Python 拓扑/方程审计或外部 checkMesh。报告格式为 `cartmesh2d-topology-native-v2`，成功状态为 `flow-converged`；旧 `flow-audited` 报告仍可读取。`continue_native_flow.py` 保留同网格、同控制续算，并按匹配控制计算惯性项。

排名同时输出原优化参数和三档共同阻力的冻结ρ评分，分别有完整记录及预先固定的8个不同材料代表。`conditionalGciReversals`保留原筛选含义；只有评价模型相同、原生渐近误差范围另有独立依据才可进入`resolvedReversals`，当前为0。这是解释资格的修正，不改变求解或几何阈值。配对网格只取共同level/留边/聚合设置，缺失、非单调、混合网格族不补GCI，档间差不充当置信区间。

Olesen Eq.(34)参考用两条宽`.2`、含引导管全长`1.5`的精确矩形通道，解析归一化功率 `Φ₀/(μU_peak²)=(96/9)(4+3.5)=80`。`paper-straight/`保存真实轮廓，`attempts.json`中6个`paper-reference`命令对应平均Re=0和133.33的level4/5/6；直管已充分发展且端口剖面与解析解一致，规定速度出口不改变这个参考解。它不复现原文压力出口优化问题或连接切换。另一个可读原始作者示例[dolfin-adjoint双管](https://www.dolfin-adjoint.org/en/stable/documentation/stokes-topology/stokes-topology.html)报告目标45.944633，但该运行以100次迭代上限停止，材料αmin、离散、投影、端口固定区及目标归一化均与本研究不同；示例文字公式与实现黏性项系数也不一致，不能拿它与本总压功直接计算相对误差。

历史数值实验的命令保存在各`*-attempts.json`的`plannedTasks.command`，依赖原版本；当前新增实验使用 `native_flow.py`、`fidelity_study.py` 和新的输出目录。SLSQP和批量驱动也拒绝覆盖现有账本；完整接续需按进度清单顺序在同级新实验目录运行，保留每批结果，再继续依赖批次。`finalize.py`从实际边界CSV另用四点Gauss动能积分重建功率，不调用生产脚本的功率函数；每个新有效档及全部网格/边界/面CSV哈希列在`native-reconstruction.json`。它仍使用本项目原生求解器，不是独立CFD求解器对照。

本次完整构建使用系统clang++；研究测试38项、前端191项通过。完整CTest的172项在沙箱通过，热缩放项因`/usr/bin/time -l`不能读取`sysctl kern.clockrate`而未启动原生程序；同网格裸命令通过，同控制主机CTest复测通过。原失败、裸网格和主机复测证据分别在`ctest-final.log`、`ctest-study-repro/`、`ctest-failure-repro/`及`ctest-thermal-host.log`。这项复测未修改计时、方程、质量或停止门，也未新增换热研究。

文献边界决定后续方法选择：

- [Vrionis等，2021](https://doi.org/10.1016/j.camwa.2021.06.002)已报告每轮Cut-cell重建与伴随优化，并提取多孔设计的真实边界重新比较。此轮核对到出版社收录的摘要，正文访问被拒绝，不能据此宣称已完成全文算法对照。
- [Theulings等，2023](https://doi.org/10.1007/s00158-023-03570-4)的开放全文分析漏流、网格尺度相关的阻力下界与不同多孔模型；“自适应调阻力”本身不能作为尚无人研究的贡献。
- [Sun等，2026](https://doi.org/10.3390/computation14010019)已用Darcy模型生成换热候选，再在统一工况用Navier–Stokes模型评价。此轮核对出版社索引中的摘要/章节内容，正文直连限流；候选复算与闭环优化不能混称。

固定评价条件下的设计选择可靠性与复算成本仍是研究问题；T02 已新增下面的在线材料更新入口。它与已有多保真方法的差异、适用范围和创新性需要进一步研究，当前实现与结果集中记录在当前状态。

### T02：真实壁面反馈驱动材料更新

`tools/optimization/closed_loop_topology.py` 从已有优化目录读入 `Problem`、完整材料变量 `design` 和最终 `q/β`，原输入不修改。每轮生成多孔伴随方向、界面开闭扰动，以及在材料场具有多分量时按实际端口生成的共用通道/最近分量桥接候选；这些结构分支明确是几何搜索提案，不称为伴随求出的拓扑。每个候选都重新分析、按物理面积额度提取完整锐壁、生成 Cut-cell 并求解。只有原生网格接受、原生流动严格收敛、流量匹配的真实目标才能替换当前接受状态。坏网格和未收敛候选保留原场与日志，控制器只缩小材料扰动重试，不改变求解或质量控制。

反馈不是离线重排旧设计：以当前材料 `x₀` 为中心，用实际新候选的差值拟合 `sᵢ·c = (J_Tᵢ−J_T₀)−(J_Bᵢ−J_B₀)`，再优化局部模型 `J_B(x)+c·(x−x₀)`。每轮重新求多孔状态与离散伴随；局部搜索限于已有复算方向的线性张成空间，生成新的完整材料场，随后另做原生复算决定接受。不同锐壁分量/孔洞/端口连接的样本不混入同一连续校正，连接改变后的新分支重新积累样本。它是采样方向上的近似校正，未实现原生 Cut-cell 伴随。

局部 SLSQP 的目标除以 `max(|J_B₀|,1)`，`ftol=1e-9` 只控制该无量纲代理小问题的停止、每轮最多 `10×material-steps` 次迭代；不作为原生精度或设计收益门。分量宽度使用已有物理滤波半径，未授予制造最小宽度资格；回归对小于 `sqrt(machine epsilon)` 的材料变化/相关奇异值作算术截断，避免放大驻点附近的舍入噪声，也不形成 CFD 误差估计。原生固定 level、留边、聚合、Re、端口、strict、`1e-8` 停止条件不随候选变化。总时间额度在原子求解之间检查，单次求解仍受 timeout 控制；达到轮数或时间额度均不称优化收敛。

```sh
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 outputs/topology-env/bin/python tools/optimization/closed_loop_topology.py --source outputs/topology-fidelity/b-retries/double-pipe-re-0 --output outputs/topology-feedback/double-pipe-rerun --rounds 5 --probes 3 --move .15 --level 5 --max-seconds 1200
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 outputs/topology-env/bin/python tools/optimization/closed_loop_topology.py --source outputs/topology-fidelity/b-retries/double-pipe-re-0 --output outputs/topology-feedback/double-pipe-control-rerun --rounds 5 --probes 3 --move .15 --level 5 --max-seconds 1200 --no-feedback
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 outputs/topology-env/bin/python tests/closed_loop_topology_test.py
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 outputs/topology-env/bin/python tools/visualization/render_topology_feedback.py outputs/topology-feedback/double-pipe-rerun --baseline outputs/topology-feedback/double-pipe-control-rerun --output outputs/topology-feedback/preview.png
```

第二条只关闭梯度校正，仍执行相同结构提案和原生接受机制；它用于区分结构搜索收益与校正增量，不能称为“完全不复算”的基线。`run.json` 记录源场/代码/二进制、物性和原生控制；`evaluations.json` 保存所有新材料、锐壁与目标、失败、耗时及生成它的反馈样本 ID；`summary.json` 保存每轮接受状态与预算终止原因。完整结果在 `outputs/topology-feedback/`，关键索引/真实网格速度图为 `artifacts/current/native-topology-feedback.json` 与同名前缀 PNG。性能结论、当前限制和待办见当前状态，不以代理拟合残差代替真实收益。

`compare_sharp_designs.py` 匹配真实提取面积、入口积分流量、物性、出口及数值控制。默认固定候选轮廓，仅偏移基准内部密度以匹配面积；`--area-target budget` 则把两者都投影到面积上限，端口不改，形成的几何重新验收。压降是通量加权的运动学压力差（m²/s²），不能直接称为瓦数或商业节能比例。

`--small-alpha` 是守恒聚合比例，改变离散网格但不放宽质量门；不同配置分组，不能挑最优档拼成细化趋势。`--reuse-baseline` 只有在几何、端口、物性/控制、实际二进制与再生成边界一致，并确认原生收敛后才可使用，不算新 CFD 运行。

```sh
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 outputs/topology-env/bin/python tools/optimization/optimize_flow.py --case double-pipe --width 2 --nx 96 --ny 48 --alpha-max 25000 --initialization uniform --iterations 80 120 300 --max-seconds 600 --no-plot --output outputs/topology/uniform-low
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 outputs/topology-env/bin/python tools/optimization/compare_sharp_designs.py outputs/topology/uniform-low --output outputs/topology/connectivity-native --area-target budget --levels 6 7 --small-alpha .25 --iterations 4000 --timeout 300
```

同控制的初始化/阻力研究另用 `--alpha-max 1000000`，分别取 `--initialization uniform` / `merged` 和新的输出目录。弯道例使用 `--case bend --nx 72 --ny 48 --volume .3 --alpha-max 1000000 --iterations 50 100 250`，按相同流程比较 level 6/7。命令是复现入口，不保证所有档位通过。

`render_connectivity_study.py` 汇总所有比较、失败网格及 mesh probes；`allAttemptedGrids` 与 `acceptedPairsOnly` 分列，同几何/控制/二进制才比较。`meshRobustImprovementObserved` 还要求全部指定档成对有效、末两档排序一致且方案差大于已观察网格变化之和；这不是误差界、置信区间或网格无关证明。

仅实际轮廓确认为两条完整矩形时，独立解析参照用 `8νUmaxL/w²`；不能用它补填失败原生基准。未收敛场的独立诊断不授予比较资格。

```sh
outputs/topology-env/bin/python tools/optimization/continue_native_flow.py <原生候选目录> --output <新目录> --iterations 2000 --timeout 300
```

该工具只给同网格稳态单分量追加预算：核对同网格/模型，提取单元和 owner 方向面初值，保留同边界、物性、松弛和容差再验收。非有限场或重建错误不能继续；不生成物理时间 checkpoint，不改原失败成对行，累计迭代/成本逐层保留。`--continued-candidates` 将结果单列，不能冒充冷启动提速。

`tests/flow_topology_test.py` 覆盖伴随差分、真实体积、下降、确定性、孔洞/端口、失败状态、缓存混用、成对排名及初值接口。`tests/fixtures/topology_oc_stall.npz` 保留 OC 上升而梯度回退可行的真实最小例；Poiseuille 单网格 2% 仅为低成本开发检查。结果和未完成项统一见当前状态。

方法背景：[Stokes 双管优化](https://www.dolfin-adjoint.org/en/stable/documentation/stokes-topology/stokes-topology.html)、[Cut-cell 拓扑优化及提取后复算](https://www.sciencedirect.com/science/article/pii/S0898122121002406)、[阻力参数敏感性](https://arxiv.org/abs/2302.14156v2)。这些方法组合和单例收益不构成原创性证明；按候选目标差、局部误差和可解性分配验证预算仍是待提出算法并做多例对照的研究假设。

## 纯笛卡尔浸入边界研究入口

`cartmesh2d_immersed_cli` 与 `immersed/CartesianFlow2D` 是独立二维实验求解器，不读取 CM2D，不改变背景的 `solver_ready=false`，也不降低 Cut-cell 质量门。

### 方程与适用范围

均匀 MAC 交错格：格心压力、面心速度，所有方格和明确标识的固体辅助未知量保留。域 `[0,L]×[0,H]`，x 周期、上下无滑移，恒 x 加速度 drive 对应周期压差；总运动学压力为 `p_fluctuation−drive*x`。当前伪时间求稳态，不是入口/出口外流或已验非定常算法。

方程 `du/dt+div(uu)=−grad(p)+ν∇²u+drive·ex−χu/η`、`div(u)=0`；对流一阶迎风、黏性中心、二者显式，阻力隐式。`β=1/(1+dtχ/η)`，压力校正解 `−div(βgrad(φ))=−div(u*)/dt`。参考压力行固定后仍检查该格连续性。χ 从原折线有符号距离构造，默认半宽 .5min(dx,dy) 正弦过渡，η 默认 1e−4；不移动或平滑原几何。

`channel` 无内固体，`cylinder` 为中心 `(L/3,H/2)`、半径 `.15H` 的 128 段多边形；`custom --boundary SOLID.xy` 支持静止多环及奇偶孔洞，不接受命名角色元数据。自交、接触等非法环拒绝；固体距外域/周期接缝至少两格加掩膜半宽，且须有完整内部单元，薄壁/欠分辨明确拒绝。

### 压力与真实壁面力耦合

可选 `--wall-method surface-penalty`：令 J 为壁面双线性速度插值，`W=ds·h/(dx·dy)`，`E=sqrt(W)J`。每段在两种 MAC 插值结点线处分区，三点 Gauss 积分精确覆盖区间内至多四次的 `JᵀWJ`，保留原始折线；共线分段及反向回归检查不变性。令 D 为散度、`B=[−hD;E]`：

```text
(B beta Bᵀ + diag(0, eta_wall/dt)) q = B u_star/dt
u_new = u_star - dt beta Bᵀ q
p_new = p_old + h q_pressure
wall_force = -Eᵀ q_wall
E u_new = eta_wall q_wall
```

压力、壁力和体积阻力用同一 β；不在求解后截断壁速。正 ηwall 提供有限顺应项，不删约束或加对角补丁；完整域壁面功率为 `−ηwall·Σq_wall²·dx·dy≤0`，独立读取器重建核对。默认 ηwall=1e−4，实测对照显式 1e−6；这不是通用最优参数或精确无滑移资格。旧中点采样漏掉点间穿透的失败保留。

auto 后端：旧 Brinkman 用 IC0，新耦合用 Jacobi，macOS 可选 Cholesky。耦合 IC0 非正主元组合明确拒绝，失败留在原工作区 `outputs/wall-treatment/gauss-ic0-probe/`；不悄悄换后端或放宽残差。上限 16,384 积分点、8,000,000 原始连接贡献，超限失败。

### 停止与输出语义

`Uref=drive·H²/(12ν)` 为无障碍平均速度。默认稳态 1e−4 分别检查动量余量/drive 和最大步变化/Uref；连续性 1e−6 检查 `max|div(u)|·H/Uref`，每步必验。新模式另验壁面方程误差，不把它当真实壁速。线性相对 L2 默认 1e−8、绝对算术余量 1e−13，仍复核真残差；这些均是开发代数门，不是通用物理精度。

| 状态 | 语义 |
| --- | --- |
| `steady-converged` / 返回 0 | 达到离散停止量，仍须验壁面/网格误差 |
| `iteration-limit` / 2 | 预算耗尽；默认 40,000 步，旧 20,000 步通道失败保留 |
| `candidate-failed` | 候选线性/连续性等失败，保留上个接受步 |
| `initial-only` | 未接受过步，不能作检查点 |
| POSIX 取消 / 130；输入导出错误 / 1 | 保留可用诊断和最后接受状态 |

```sh
cmake --build build --target cartmesh2d_immersed_cli --config Release --parallel 2
ctest --test-dir build -C Release -R '^cartmesh2d_immersed_flow$' --output-on-failure
build/cartmesh2d_immersed_cli --case cylinder --nx 128 --ny 32 --nu .01 --drive .12 --max-steps 40000 --wall-method surface-penalty --wall-penalty-time 1e-6 --output outputs/immersed-wall
python3 tools/visualization/render_immersed_flow.py outputs/immersed-wall --output outputs/immersed-wall/preview.png
```

只比较旧模式时去掉 wall 参数并换新输出目录；macOS 可加 `--linear-solver cholesky`。绘图才需要 NumPy/Matplotlib，原生求解器不需要。输出记录网格/边界/求解/导出计时，压力时间已包含在求解总时间中。

u/v/cells CSV 与 VTK 保存完整格、掩膜、表面力和分类；`wall-markers.csv` 是积分点/权重/乘子，`walls.csv` 是独立原折线稠密采样（含端点、间距≤h/8），history 只存接受步。阻力同时读表面与体积两部分，不能将其中一项或辅助域总反力当已验表面应力系数；法向通量积分为长度²/时间，辅助格散度小不保证真物面无穿透。

18 项直接检查在 `tests/immersed_flow_cli_test.py`：解析通道、阻滞/对称、孔洞、确定性、失败/取消、真实折线和导出破坏检测、分段不变性、两种后端等。稠密壁速与 Gauss 最大值的 7/3 包络来自二次插值基函数界，用于抓漏点，附加 1e−12Uref 仅是读回算术余量；单格分辨率误差门不代表工程验收。原两档网格和相位试验的限制见当前状态，全量软件回归已由集成另行覆盖。

没有共形流体 polyMesh，外部 checkMesh 不适用，独立读取器不替代 Cut-cell 门。参考：[Brinkman 惩罚](https://www.math.u-bordeaux.fr/~chabrune/publi/ABF-NM.pdf)、[惩罚压力投影](https://arxiv.org/abs/2306.06277)、[Taira–Colonius 原作者论文目录](https://www.seas.ucla.edu/fluidflow/pubs.html)（DOI 10.1016/j.jcp.2007.03.005）；本实现是有限顺应的双线性表面惩罚，不声称复现整套算法或资格。

## 运行工具与证据

`tools/flow/` 保存计算流程需要的工具：`native_mesh.py` 读取 CM2D/CSV 并计算几何，`extract_steady_iterate.py` 提取同网格初值，`prolongate_regular_flow.py` 映射完整矩形张量格，`prepare_*.py` 准备热化学机理和火焰输入。绘图所需解析曲线和 OpenFOAM 数据读取分别在 `analytic_flow.py`、`openfoam_data.py`。这些工具不授予额外的审计通过标记。

原生几何拓扑、Solver 质量、数值有限性、物理约束和收敛判定继续保留；达到迭代上限仍然失败。单元测试、前端测试与历史数值/物理结果各自说明范围。旧独立验证代码可从 Git 历史或本次源码快照恢复，当前维护入口不再调用它们。

已有 gzip 场按 `compressed*.json` / `*-compressed.json` 原 SHA 核对后按需恢复，保留压缩包；`*.canonical-fields.json` 指向规范场，勿批量展开历史。曲壁映射等未完成方案继续保留在原 `outputs/`。

## 历史展示素材

[展示素材目录](../展示素材/)的六图来自固定历史网格重绘，未重跑 CFD。原记录的 Solver 与标准 checkMesh 通过，扩展 checkMesh 失败；历史 Q1 已退出产品验收，不用于评价当前版本。这是三个固定几何的约 13–16 万格展示，不代表任意模型或网格无关性，圆为 32 段折线。

| 图片 / 单元数 | 原网格位置（历史路径） | 原网格 SHA256 |
| --- | --- | --- |
| [圆 Cut-cell：153,208](../展示素材/01_circle_pure_153208.png) | `outputs/engineering-pure-final-high/circle/mesh.solver.cm2d` | `5bf6d64d721959560b91121d1e5ff7052086d344f3e33eb85a45d2066d8772cd` |
| [NACA 2412 Cut-cell：133,810](../展示素材/02_naca_pure_133810.png) | `outputs/engineering-pure-final-high/naca/mesh.solver.cm2d` | `227e115ee953475315df65d2dbd9cd2909a05e23c80719a71d8b39275a4c82cd` |
| [喷管 Cut-cell：141,488](../展示素材/03_nozzle_pure_141488.png) | `outputs/engineering-pure-final-high/nozzle/mesh.solver.cm2d` | `12242e04541fc5336332cfb098ef3e6583e55bf82936c4899b56c4bb792fb0a9` |
| [圆混合：164,970](../展示素材/04_circle_hybrid_164970.png) | `outputs/engineering-concave-union-batch/mesh.hybrid.solver.cm2d` | `cc0cdeda3ac26982f9037799039a7be7af2cad29b37effba96eeac17ba87ffad` |
| [NACA 2412 混合：140,305](../展示素材/05_naca_hybrid_140305.png) | `outputs/engineering-hybrid-final-high/naca/mesh.hybrid.solver.cm2d` | `c4ce8ceb2930f30dbff918c0ebce09e31c36057cfa5d3cf4671af365af3c5acf` |
| [喷管混合：148,375](../展示素材/06_nozzle_hybrid_148375.png) | `outputs/engineering-hybrid-final-high/nozzle/exact/mesh.hybrid.solver.cm2d` | `463850cca4d728eeab418458f393c11d608687afe1a9e27a0d5e603e8d4d8791` |

参数：圆/翼型壁面 h/Lref=.0025，喷管 .00125；背景 .0125、过渡带 40。混合网格首层为壁面目标 1/4，4 层、增长率 1.2；Lref 分别 2/1/6 m。这些不是从 y+ 或流场误差反推的通用配置。公开展示需固定版本时，图片链接中的 main 换具体提交；原路径可能因既往归档变化，哈希按生成时解释。

## 历史与工作区

当前源码只在根目录维护，旧工作区保留研究输入、真实场、失败和续算材料。历史查 Git，例如 `git show f7733a2:docs/DEVELOPMENT_CN.md` 或 `git log --all --oneline -- src/quality/SolverTopology2D.cpp`；旧 archive 标签只作历史。不要再复制阶段文档或整份源码。
