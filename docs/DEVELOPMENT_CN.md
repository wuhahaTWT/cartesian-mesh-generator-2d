# 开发指南

当前能力与待办只看[当前状态](CURRENT_STATE_CN.md)，操作看[桌面使用](DESKTOP_APP_CN.md)，修改约束看[AGENTS](../AGENTS.md)。本页集中构建、代码入口、方法与复现，不重复逐轮实验结果。

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
# Windows 用 pack:win；Linux 用 pack:linux。
```

包输出到 `desktop/dist/`。构建检查实际 Mach-O/PE/ELF、架构及 runtime，拒绝跨系统混装；Windows 使用静态 CRT 和 UTF-8 路径 manifest，macOS 检查动态库。ZIP 使用流式 yazl，不依赖外部压缩命令。

真实 App 检查入口为 `tools/verification/desktop_platform_smoke.py --help`，三平台工作流使用中文及空格路径，覆盖 PNG/JPG、hybrid、两种背景和两类可压壁面样例，并独立读取 ZIP。Linux `--no-sandbox` 仅用于隔离 CI 虚拟显示，不进入产品启动参数。

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

## 完整笛卡尔背景网格

`--background-grid adaptive|uniform` 在几何诊断、Quadtree 细化及 2:1 平衡后直接导出完整叶子，不做 Cut-cell、Solver 修复或 OpenFOAM 输出。均匀模式使最低层级等于最高层级；自适应复用尺寸场和盒加密。

```sh
build/cartmesh2d_cli examples/acceptance/circle.xy outputs/background-grid/circle 7 0.5 0.1 exterior - 3 --background-grid adaptive
python3 tools/verification/check_background_grid.py outputs/background-grid/circle.background.json
```

`cartmesh2d-background-v1` JSON 保存域、原边界、完整单元范围/层级/整数格坐标及分类（0 外部、1 内部、2 相交），`solver_ready=false`；VTK 为相同完整四边形，粗细交界可能有悬挂节点，不宣称共形求解面拓扑。独立审核检查覆盖、无重叠、分类、2:1、JSON/VTK 一致与确定性。

CLI 均匀/自适应资源上限为 level 10/12；桌面分别为 9/10，且自适应全域最低层级最高 8。它们是资源限制，不是性能或任意几何成功保证。`desktop/src/core/background-grid.js` 阻止将背景数据当流体网格。

## 原生层流求解

`Incompressible2D.cpp` 使用 SIMPLE、Rhie–Chow 及共享压力/黏性面通量；只读取最终 `*.solver.cm2d`。完整参数查 `build/cartmesh2d_flow_cli --help`。

```sh
build/cartmesh2d_flow_cli --mesh outputs/channel.solver.cm2d --case channel --export-boundaries outputs/channel.boundaries
build/cartmesh2d_flow_cli --mesh outputs/channel.solver.cm2d --case custom --boundary outputs/channel.boundaries --output outputs/channel-flow --nu 0.1 --speed 1 --max-iterations 5000 --tolerance 1e-9
```

命令示范结构，几何须先满足模板约束；边界文件由原生导出后修改，保持面编号/几何绑定。`pressure-opening` / `symmetry` 只支持轴对齐面。压力为 p/ρ（m²/s²），流量为每单位厚度 m²/s、向外为正。

| 控制 | 约束 |
| --- | --- |
| `--convection upwind\|limited-linear\|face-limited-linear` | 迎风或限制重构，不是物理模型切换 |
| `--pressure-preconditioner ic0\|aggregation\|cholesky` | Cholesky 只在 macOS；仍验真实线性残差 |
| `--linear-policy strict\|adaptive` | adaptive 最终严格复核 |
| `--velocity-relaxation` / `--pressure-corrections` | 默认 .6 / 4；困难曲壁可能不稳，连续性小不代表收敛 |
| `--steady-acceleration anderson` | 默认 none；限稳态、不支持材料更新，候选须降低原残差并过守恒门 |
| `--initial-guess` / `--initial-flux` | 同目标网格稳态初值；面初值须与单元初值同用，不是物理检查点 |

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
python3 tools/verification/verify_euler.py --mesh final.solver.cm2d --prefix outputs/euler/run
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
| `tests/euler_accuracy_cli_test.py` | 面积加权光滑波/涡及 Sod 误差；旋转与 Mach 3/10 扰动激波，保留最小失稳例，不能外推任意强激波 |
| `tests/euler_conduction_cli_test.py` | 热扩散及 Euler–Fourier 小扰动；独立连续模态矩阵指数、物理热边界、SI 相似性和续算绑定 |
| `tests/euler_viscosity_cli_test.py` | 周期剪切波、含纵向黏性的热模态、Couette 解析初值保持；与冷启动验收分开 |
| `tests/wall_gradient_test.cpp` | 二次场、扭曲/旋转、长宽比 .03/3/30、长度缩放、实际 Jacobian 扰动及缺秩拒绝 |
| `tests/transport_precision_test.cpp` | 无源光滑导热 8/16/32 格，比较解析壁面边积分、排除线性求解容差污染 |
| `tests/euler_wall_accuracy_cli_test.py` | Couette 解析保持及静止冷启动至 t=60；功热与全过程能量、质量决定的稳态压力，以及固定网格时间细化 |

测试阈值按量分别解释，均为既有门，整理文档不改断言：面通量审计使用含对流/热/黏性量纲包络的 512 epsilon，算子/SI 采用 1024–4096 epsilon 检查浮点一致性；极端长宽比另记录输入舍入传播，不能替代原几何门。物理回归使用各自归一化离散误差、网格观测阶和解析参照，不用机器 epsilon 充当工程精度。

例如光滑导热 `T=2+.2sin(πx)sinh(πy)/sinh(π)` 的热流 L1 以解析边积分绝对值总和归一化，最细目标 .5%、至少优于线性 2 倍、观测阶>1.7；1e−11/1e−13 线性容差场差小于离散误差 1%。Couette 用解析温升和壁速归一化，另验压力细化，避免多项式恰好复现掩盖全场误差。准确控制、门限依据和成本在相应测试及[精度证据](../artifacts/current/native-euler-wall-accuracy.json)中。

时间细化使用同网格/终点，dt=1e−3/5e−4/2.5e−4，对更细轨迹再核参考误差，属于时间自收敛。`EulerStepControls2D::endTime` 保留浮点尾步最小失败例：必要时把倒数第二步拆成两个合法小步，真实算通量，不伪造时钟或绕过最小步长。篡改通量、速率、壁面格式/面数必须被独立审计拒绝。

```sh
ctest --test-dir build -R 'cartmesh2d_(wall_gradient|transport_precision|euler_wall_accuracy)' --output-on-failure
python3 tests/euler_wall_accuracy_cli_test.py --cli build/cartmesh2d_euler_cli --output outputs/euler-wall-accuracy/validation
build/cartmesh2d_euler_cli --mesh final.solver.cm2d --output outputs/euler/run --case external --viscosity .02 --wall-model no-slip --conductivity 100 --wall-thermal temperature --wall-value 400 --flux hllc --order 2 --wall-gradient quadratic --density 1.225 --pressure 101325 --u 50 --end-time .00005
```

最后一条用短时功能测试物性，不是空气推荐值。精度脚本可能需数分钟，Windows 预算更长；日常只选相关项，完整范围与未验物理问题见当前状态。

## 详细燃烧化学基础

`include/cartmesh2d/chemistry/DetailedGas.hpp` 与 `src/chemistry/DetailedGas.cpp` 是详细反应流的原生热化学基础，使用 Cantera **3.2.x C++ API**。模块不依赖 Python 或三维核心；验证脚本另用 Python Cantera 3.2.0。机理必须给出实际文件路径；当前明确接受中性、单气相理想气体机理，其他热力学/相模型显式拒绝。

- 保守状态为 `rho`、`rho*e` 和各 `rho*Yk`，其中 `e` 包含生成能，允许负值。输入组分不自动归一化或裁剪；能量到温度的反解限制在所有组分热力学数据的共同温区。积分阶段的显式舍入闭合见下文，不能用于修复不合法输入。
- 返回温变比热、焓、反应源，以及 `multicomponent` 扩散矩阵和 Soret 热扩散系数；矩阵为列主序，须配合完整通量公式，不能当作各组分独立的 Fick 系数。
- 恒容绝热化学子步使用 Cantera `Reactor` / CVODES，内部能量为守恒变量；化学变化通过组分及温度体现。`-sum(hk*omega_k)` 只作放热诊断，不重复加入已经含生成能的总能量方程。反应阶段超出共同物性温区也会拒绝。
- 失败返回空 `accepted` 和原因；输入不变，下一次调用从调用者保留的接受状态重新开始。当前是串行、每工作线程独享的化学上下文；空间扩散、耦合时间推进和原生检查点见下文，产品输入与 App 尚未接入。

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

# Python 环境单独用于核对；输出目录必须不存在，以保留上轮成功或失败记录。
build/chemistry-env/bin/python tools/verification/verify_detailed_gas.py \
  --probe build/cartmesh2d_detailed_gas_probe \
  --mechanism-root build/deps/cantera/share/cantera/data \
  --output outputs/combustion-foundation/interface-new
```

本机依赖构建环境为 Python 3.14、SCons 4.11.1、packaging 26.3；Python 对照环境另装 Cantera 3.2.0、NumPy 2.5.3。Cantera 构建使用 `python_package=n f90_interface=n hdf_support=n clib_legacy=yes googletest=none example_data=no debug=no`，`fmt/yaml-cpp/Eigen/SUNDIALS` 使用上游固定子模块。`clib_legacy` 仅避免生成本项目不使用的 C 接口，不改变 C++ 模型；上游完整测试没有运行，原生接入检查单列。动态依赖需检查，不得引入 mesasdk 库。正式打包和跨平台部署仍待完成。

检查按用途分开：原生测试核对状态反解、元素源、完整机理规模、负生成能、扩散数据、反应积分及失败恢复。元素源舍入检查以**总生成与消耗量**归一化，不能用已相消的净源作舍入尺度；首轮甲烷氮元素因此误报，原日志保留在 `outputs/combustion-foundation/first-native-test-failure.log`。接口脚本直接计算 NASA7 多项式，独立核对混合物比热/内能；反应轨迹和扩散与另装的 Cantera Python 接口比较，仍共享化学后端，不能称为独立机理或实验验证。

阈值只限定数值接入误差：温度往返 `2e-9 K`、NASA7 及接口值采用脚本声明的归一化差异；积分对照检查温度相对差和质量分数绝对差。默认化学子步的元素质量分数漂移门为 `1e-8`，能量变化以 `max(|rho*e|,rho*cv*T)` 归一化后为 `1e-8`；它们是可配置的子步守恒检查，不是火焰精度目标。上述少量案例为秒至分钟级，网格/火焰/真实工况验收尚未完成。

### 可选 NASA7 连续热力学输入

原始分段 NASA7 拟合的焓/内能可能在拼接点跳变；在很窄的能量区间，温度反解因而不唯一。`tools/verification/prepare_continuous_nasa7.py` 显式生成新的机理文件，保留低温段锚点、两段全部比热系数、温区、组分、输运参数及完整正向反应定义。依据 [NASA7 形式与比热积分关系](https://cantera.org/stable/reference/thermo/species-thermo.html#the-nasa-7-coefficient-polynomial-parameterization)，在拼接点以 60 位算术计算高温段的 `a5`、`a6`，使 `h = h_ref + ∫cp dT`、`s = s_ref + ∫cp/T dT` 连续；最终系数仍为双精度。比热本身的拼接跳变和导数未被平滑。这里只支持理想气体的两段 NASA7，其他模型显式拒绝。

这是有物理影响的数据变更：高温焓、熵、平衡常数和逆反应速率会改变。工具不修改原文件，也不接入求解器的隐式修正。输出 `original-resolved.yaml`、`continuous.yaml` 和每个组分的改动/哈希报告；已有输出目录拒绝覆盖。Cantera 展开原机理时可能因活化能单位换算产生末位舍入，报告逐项列出；准备器直接修改展开后的 YAML 树，保证两份展开文件的正向反应定义相同，另用原输入进行敏感性对照。原相配置的输运模型也保留，实际原生通量继续显式使用完整多组分模型。

```sh
build/chemistry-env/bin/python tools/verification/prepare_continuous_nasa7.py \
  --mechanism build/deps/cantera/share/cantera/data/h2o2.yaml \
  --output outputs/combustion-foundation/continuous-h2-new
build/chemistry-env/bin/python tools/verification/verify_continuous_nasa7.py \
  --prepared outputs/combustion-foundation/continuous-h2-new \
  --probe build/cartmesh2d_detailed_gas_probe --fuel H2 \
  --output outputs/combustion-foundation/continuous-h2-audit-new
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
build/chemistry-env/bin/python tools/verification/verify_diffusive_flux.py \
  --probe build/cartmesh2d_diffusive_flux_probe \
  --mechanism-root build/deps/cantera/share/cantera/data \
  --output outputs/combustion-foundation/diffusion-interface-new
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
# 读取/绘图环境另需 matplotlib，本机为 3.11.2；不进入 C++ 运行依赖。
build/chemistry-env/bin/python tools/verification/verify_reacting_flow.py \
  --case outputs/combustion-foundation/reacting-half-new \
  --output outputs/combustion-foundation/reacting-half-review-new
build/chemistry-env/bin/python tools/verification/verify_reacting_flow.py \
  --case outputs/combustion-foundation/reacting-new \
  --compare outputs/combustion-foundation/reacting-half-new \
  --output outputs/combustion-foundation/reacting-review-new
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

接受状态只在真实状态、物性和积分时钟检查后更新。会话保留 BDF 历史；步间取消可继续同一会话，内部试算取消或失败保留最后接受状态，但须用新会话重建历史。残差评估和接受步都有成本预算，不把预算耗尽记为完成。火焰验证驱动保存 `evaluation-progress.jsonl`；在输出目录创建 `cancel.request` 可有序停止并保存最后接受状态、检查点和失败原因。直接杀进程不提供同样的保存保证。

长时间验证可用 `--sample-every N` 将初始状态、每 N 个接受步及末态写入 `samples.jsonl`。采样不改变时间步或 BDF 历史，使用单独的物性和残差上下文；额外诊断调用在 `sampleResidualEvaluations` 中计数，仍包含在完整进程耗时内。默认关闭采样；开启后，`--max-samples` 默认 256 条，用于限制诊断输出数量，不是物理精度门。达到保存上限或诊断出错会在接受状态退出并尝试写入最后接受场和检查点，原失败原因进入结果；部分轨迹明确标记 `sampling_complete=false`，不冒充完整计算。

读取器用同一套独立时钟、EOS、共享面和质量/元素/能量预算检查逐帧核对。`initial` 表示本段起点：冷启动来自参考初值，续算来自真实接受检查点，两者通过 `restart` 元数据区分。`sample-audits.json` 保留每帧诊断；历史图使用真实接受时刻。固定本段初始温度范围中点的等温位置通过单元中心剖面线性插值得到，多个交点全部记录；该位置标记不等于已验收火焰速度。464 格、0.1 微秒采样开关对照的全部物理场相等，步长日志和检查点逐字节相同。故意改坏的时钟和质量预算被拒绝；两条样本上限测试在第 20 个接受步退出，保留的场和预算与正常运行的同一步相同。

显式、隐式原生验证驱动均支持 `--restart 检查点`；Python 运行/审计入口使用 `--restart-checkpoint 检查点`，一次只选一档网格。先用原生绑定检查完整机理、网格、边界和模型，原样保存输入到新输出目录的 `restart.checkpoint`；监督程序核对源文件运行前后及副本 SHA256。`--duration` 是从检查点起追加的物理时长，`endTime` 为实际请求终点，累计接受步数保留。步长预算、方程调用及收支积分从本段起点重新计数，采样步距也从本段起点计算。读取器核对原样导入的完整场、时间和步数，并分别报告本段/累计接受步数；不能把重复计算或重置为零的时钟当作续算。

检查点保存物理状态，BDF 历史及初始积分尺度在新会话重建，因此续算轨迹不要求逐字节等于不中断运行。原诊断上限失败的 464 格案例现已从第 20 步、约 `0.0426812` 微秒续至 `0.1` 微秒；新增 23 步，6 帧审计通过。与不中断同终点的最大温差约 `1.98e-5 K`，合并两段实际通量后的质量/元素/能量预算仍通过原检查；未扩大到长期误差资格。冷启动物理场、步长日志和检查点保持原结果。不同网格/入口/绑定、非法时钟、截断及尾部损坏被拒绝；零时长读回保留全部原值并且不新增接受步。该案例也经过显式入口的导入与推进检查，不将它视为两种积分方法的精度等价证明。

本机可选后端目前要求 Cantera 随附并导出的 SUNDIALS 5 接口；其他 SUNDIALS 主版本、App 和平台打包尚未适配。先使用已有的详细化学配置，再运行：

```sh
cmake -S . -B build -DCARTMESH2D_BUILD_REACTING_IMPLICIT=ON
cmake --build build --target cartmesh2d_reacting_implicit_tests cartmesh2d_reacting_implicit_flame_probe -j2
ctest --test-dir build -R '^cartmesh2d_reacting_implicit(_trace)?$' --output-on-failure
build/chemistry-env/bin/python tools/verification/verify_reacting_flame.py \
  --probe build/cartmesh2d_reacting_implicit_flame_probe \
  --reference outputs/combustion-foundation/flame-reference-new \
  --output outputs/combustion-foundation/implicit-flame-new --duration 1e-6 --grid 0 --jobs 1 --reflect-species 1
# 较长演化可加 --sample-every 200；时间步仍由原误差控制决定。
# 从既有 native-0/accepted.checkpoint 追加计算，用新的输出目录：
build/chemistry-env/bin/python tools/verification/verify_reacting_flame.py \
  --probe build/cartmesh2d_reacting_implicit_flame_probe \
  --reference outputs/combustion-foundation/flame-reference-new \
  --restart-checkpoint outputs/combustion-foundation/implicit-flame-new/native-0/accepted.checkpoint \
  --output outputs/combustion-foundation/implicit-flame-resumed --duration 1e-6 --grid 0 --jobs 1 --reflect-species 1
# 同一 10 微秒案例分别用默认、减半和四分之一容差保存至三个新目录。
# 例如减半：在 verify_reacting_flame.py 命令中加入
# --rtol 5e-8 --conserved-atol 5e-13 --species-atol 5e-19
build/chemistry-env/bin/python tools/verification/compare_reacting_flame_time.py \
  --report outputs/combustion-foundation/time-base/report.json \
  --report outputs/combustion-foundation/time-half/report.json \
  --report outputs/combustion-foundation/time-quarter/report.json \
  --output outputs/combustion-foundation/time-comparison-new
```

`compare_reacting_flame_grid.py` 比较同机理、同终点及同积分控制的平面条带网格。细格守恒量按真实矩形重叠面积积分到粗格，再使用原机理热力学反演温度；不直接平均温度或将不同位置的数组相减。初始映射差异、终点场差异及两者相减得到的演化差异分别报告。重叠权重的常数保持和体积守恒检查仅允许 `64*epsilon*参与格数` 的浮点运算舍入，属于精确几何恒等式诊断，不是新增 CFD 精度门；另有非嵌套矩形分段常数解析核对。非嵌套投影本身有离散误差，因此不把这些差异标作正式空间阶数。每档还独立积分实际燃料消耗，列出组分残差和温度漂移。示例：

```sh
build/chemistry-env/bin/python tools/verification/compare_reacting_flame_grid.py \
  --case outputs/combustion-foundation/grid-coarse/report.json 0 \
  --case outputs/combustion-foundation/grid-refined/report.json 1 \
  --case outputs/combustion-foundation/grid-refined/report.json 2 \
  --output outputs/combustion-foundation/grid-comparison-new
```

### 空间火焰对照与微量组分回归

`prepare_reacting_flame.py` 使用独立的 Cantera `FreeFlame` 稳态边值求解器，保留完整多组分和 Soret，保存三档参考细化的原始剖面与机理。默认验证夹具为 `H2:1.6,O2:1,N2:3.76`、400 K、1 atm；网格在反应层附近取 40/20/10 微米，并保留完整上游和下游。该网格布局按本氢气案例设置，其他燃料须重新确认解析度，参数接口不等于已验证。

参考解的舍入缺陷显式处理：保存原始负微量及组分和误差，只在验证初值准备阶段将负舍入值置零，并由最丰富组分闭合。允许改变量为 `10*BVP绝对容差 + 64*(N+1)*epsilon`，超过则拒绝；有意输入 `-1e-4` 会失败。保留原参考点的密度、动量及含生成能总能量，记录元素变化，再对分段线性守恒剖面逐格积分。该准备程序不进入原生推进或续算，初值不是接受检查点。

`cartmesh2d_reacting_flame_probe` 建立真实二维条带拓扑，左侧固定储库、右侧外推流出，上下滑移绝热。输出初始/最终守恒场、实际共享面通量、分开的化学与输运残差、开放边界积分、逐步接受/拒绝、展开机理与最后接受检查点。`evaluateResidual` 是瞬时半离散算子诊断，不会把输入标为已收敛或接受步。`verify_reacting_flame.py` 独立读取实际文件，检查几何关联、质量/元素/动量/能量预算、EOS、未被修改的导入初值与真实时钟，并绘制单元场和残差趋势。数值检查、达到物理终点和火焰资格分别记录。

隐式驱动可用 `--species-atol` 显式比较组分局部积分权重，实际取值写入结果；质量/元素/能量等物理审计判据不随该选项改变。

`--grid 0` 可只运行首档，重复 `--grid` 选择多个参考夹具；报告明确记录所选索引，不能将单档报告称为三档验收。

读取器可用 `--audit-native-output <已结束运行的输出目录>` 单独复核；`--output` 仍必须是新目录。此模式校对原运行的程序、参考、终点与原始文件哈希，将新审计写入新目录，保留原求解日志和既有错误，不重新计算或覆盖原生场。原运行退出失败仍单列，不能因部分场可读而变成成功。

`tests/fixtures/reacting_flame_trace_muscl.fixture` 和 `reacting_flame_trace_soret.fixture` 是从首档剖面分别截取的 20/36 格最小失败片段，保留原始微量组分。CTest 要求在不拒绝候选的情况下推进到 `1e-8 s`，并沿用装配和元素预算检查，防止靠大量重试掩盖旧错误；它们只验收数值修复。

细档长序列曾在第 1209 个接受步后因组分和的舍入误差达到 API 边界而失败。`reacting_flame_trace_closure.fixture` 保留原检查点中的 6 格片段，旧版本可直接复现。`SpeciesMassClosure.hpp` 在化学/输运阶段及面状态插值中显式使用质量约束：先检查全部原始组分非负、有限，并确认质量缺陷除以总质量不超过原有 `64*(N+1)*epsilon` 求和额度，再从总质量减去其余组分的补偿求和，确定最丰富组分。不缩放其他组分、不设微量下限、不改变密度或能量；不合法输入和超过舍入范围的缺陷继续拒绝。完整反应机理及全部组分方程仍参与计算，这只是冗余质量约束的数值闭合。

`ChemistryStep::massClosureChange` 和耦合步的化学/输运闭合向量、最大单元相对改变量、绝对改变量积分均显式输出。SSPRK2 第一阶段改变量按半权重计入。离散装配残差扣除已记录的输运舍入操作，但实际质量/元素/能量预算**不扣除**它们，避免用闭合制造物理守恒证据。火焰读取器逐步核对诊断之和；2048 次真实恒容化学重启、负微量与超额缺陷拒绝也纳入原生检查。末态原始守恒值先保存到 `accepted-state.json`，再求末态残差，诊断失败仍保留最后接受状态。

`transport(A) → viscosity(B) → transport(A)` 的最小复现曾使同一个 A 状态的导热系数与 Soret 系数不同。原因是黏性查询使二元扩散缓存失效，而旧热缓存仍保留 A 的温度键；重新建立二元矩阵时覆盖了热输运使用的自扩散对角项。`TraceStableMultiTransport::update_T` 现在在每次温度改变时同步使热缓存失效。该查询序列逐值一致性和有/无显式步长估计的真实物理残差一致性都进入原生回归；48 项局部输运对照重新通过。旧性能记录仍保留，不能据此跳过物性正确性检查。

完整物性查询将多组分矩阵查询放在导热/Soret 查询之前，复用同一热输运矩阵求解；黏性应力用 `DetailedGas::viscosity` 单独查询。`compare_reacting_cost.py` 针对这类缓存/查询优化交替执行旧、新完整 48 格点火计算，计时含冷启动至全部输出，要求数值场逐值一致、步日志/机理/检查点逐字节一致，只排除墙钟字段。一次三对实测中位数约 `8.274→7.927 s`；这是引入质量闭合之前的隔离性能比较，不能作为当前所有改动或任意火焰规模的提速结论。

```sh
cmake --build build --target cartmesh2d_reacting_flame_probe -j2
ctest --test-dir build -R '^cartmesh2d_reacting_flame_(muscl|soret|closure)$' --output-on-failure
build/chemistry-env/bin/python tools/verification/prepare_reacting_flame.py \
  --mechanism build/deps/cantera/share/cantera/data/h2o2.yaml \
  --output outputs/combustion-foundation/flame-reference-new
build/chemistry-env/bin/python tools/verification/verify_reacting_flame.py \
  --probe build/cartmesh2d_reacting_flame_probe \
  --reference outputs/combustion-foundation/flame-reference-new \
  --output outputs/combustion-foundation/flame-comparison-new --duration 1e-6
```

参考 BVP 采用恒定热力学压力，原生平面可压方程还保留纵向动能与黏性功，因此不是完全相同的方程组。从参考剖面初始化后的短时漂移、反应积分或初始残差下降均不能当作原生独立预测的火焰速度；原生稳态传播、时间/网格误差和实验还需另验。上述驱动尚非用户可配置的通用燃烧产品。湍流燃烧、喷雾、辐射、共轭传热、实际燃烧室、App、平台资格和大规模成本均未完成；完整目标和最新证据统一见[当前状态](CURRENT_STATE_CN.md#持续目标成熟燃烧模拟)。

## 流体拓扑优化研究入口

`tools/optimization/` 是独立命令行研究工具，不替换产品求解器，也不读取背景 JSON 作流体格。可选 Python 依赖不进入桌面运行包；研究测试独立于 CTest，由 `flow-research.yml` 在三平台执行。

### 模型与停止条件

`brinkman.py` 在二维 MAC 格上求无量纲 Stokes–Brinkman 方程 `−μ∇²u+∇p+α(ρ)u=0`、`div(u)=0`，**不含对流惯性**。ρ=1 为流体、0 为有限阻力固体，`α=αmax·q(1−ρ)/(q+ρ)`。抛物线端口积分匹配，两格端口固定流体、外壁固定固体；锥形滤波及 tanh 投影后检查真实 `mean(ρ)≤volume_fraction`，滤波半径不是已证明的最小制造壁厚。

消去指定速度和一个压力参考后，用离散矩阵转置求伴随，贯穿滤波/投影导数。默认目标为黏性与阻力耗散；`pressure-power` 是端口邻格压力的通量加权差，不是精确边界应力功。hybrid 更新先试 OC，非下降时尝试有移动界的梯度步并二分校正实际体积；接受前重新求解、检查方向与回溯，不是精确非线性投影。

三阶段 `(q,β)=(.01,0),(.1,2),(.1,6)` 改变模型，跨阶段目标不可直接比较。`projectedKkt` 是归一化梯度在体积切平面/变量界上的投影步无穷范数，默认 1e−3 只用于局部设计停止；线性/伴随 1e−8 检查归一化方程残差，连续性按总入口流量归一化，均不代表物理精度。预算耗尽、停滞和分析失败明确记录，保存最后接受设计与阶段参数，目前无优化自动续算。

### 运行与真实壁面连接

以下使用 POSIX venv 路径；Windows 使用其 `Scripts/python.exe` 并设置同名线程环境变量。

```sh
python3 -m venv outputs/topology-env
outputs/topology-env/bin/python -m pip install -r tools/optimization/requirements.txt
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 outputs/topology-env/bin/python tests/flow_topology_test.py
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 outputs/topology-env/bin/python tools/optimization/optimize_flow.py --case double-pipe --output outputs/topology/double-pipe
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 outputs/topology-env/bin/python tools/optimization/verify_extracted_flow.py outputs/topology/double-pipe --output outputs/topology/double-pipe-native --levels 5 6 --tolerance 1e-8
```

使用新输出目录。默认 48×32、25/40/60 步预算；`--initialization uniform|geometric|merged` 区分无内部路径、参考几何和人为合流种子。`--no-plot` 也跳过轮廓提取，可再运行 `topology_artifacts.py <目录>`。summary/history/NPZ/阶段快照及分析 VTK 保存控制、源码哈希、停止原因和设计，分析场仍含固体阻力。

轮廓提取保留所有连通分量与孔洞，不平滑/填孔；`fluid.xy` 显式 interior，各独立域另存文件。轮廓面积与密度体积分别报告，连通矩阵只表示几何可达，不等于混合效率。连接脚本真实执行 Cut-cell、质量门、独立面积读取、原生低 Re 流动和独立方程审计；所有分量均须成功，有孤立区/盲腔则明确拒绝。仅全局端口覆盖面可作入口/出口，不能按组件局部包围盒制造端口。

原生对照显式用速度尺度 .02、ν=1（端口名义 Re 约 .00333）、抛物线入口和压力出口；1e−8 是这些例子的代数控制，不改产品默认或审核门。检测到 checkMesh 时实际运行，否则记 not-run；多孔目标下降不能直接当真实壁面性能。`render_native_flow.py` 从接受的多边形及 CSV 作图并保留源哈希。

### 公平对照、失败与继续迭代

`compare_sharp_designs.py` 匹配真实提取面积、入口积分流量、物性、出口及数值控制。默认固定候选轮廓，仅偏移基准内部密度以匹配面积；`--area-target budget` 则把两者都投影到面积上限，端口不改，形成的几何重新验收。压降是通量加权的运动学压力差（m²/s²），不能直接称为瓦数或商业节能比例。

`--small-alpha` 是守恒聚合比例，改变离散网格但不放宽质量门；不同配置分组，不能挑最优档拼成细化趋势。`--reuse-baseline` 只有在几何、端口、物性/控制、实际二进制与再生成边界一致，且重新独立审计后才可使用，不算新 CFD 运行。

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

该工具只给同网格稳态单分量追加预算：核对哈希/模型，独立重建，提取单元和 owner 方向面初值，保留同边界、物性、松弛和容差再验收。非有限场或重建错误不能继续；不生成物理时间 checkpoint，不改原失败成对行，累计迭代/成本逐层保留。`--continued-candidates` 将结果单列，不能冒充冷启动提速。

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
python3 tools/verification/verify_immersed_flow.py outputs/immersed-wall --require-converged
python3 tools/visualization/render_immersed_flow.py outputs/immersed-wall --output outputs/immersed-wall/preview.png
```

只比较旧模式时去掉 wall 参数并换新输出目录；macOS 可加 `--linear-solver cholesky`。绘图才需要 NumPy/Matplotlib，原生和独立读取器不需要。输出记录网格/边界/求解/导出计时，压力时间已包含在求解总时间中。

u/v/cells CSV 与 VTK 保存完整格、掩膜、表面力和分类；`wall-markers.csv` 是积分点/权重/乘子，`walls.csv` 是独立原折线稠密采样（含端点、间距≤h/8），history 只存接受步。阻力同时读表面与体积两部分，不能将其中一项或辅助域总反力当已验表面应力系数；法向通量积分为长度²/时间，辅助格散度小不保证真物面无穿透。

18 项直接检查在 `tests/immersed_flow_cli_test.py`：解析通道、阻滞/对称、孔洞、确定性、失败/取消、真实折线和导出破坏检测、分段不变性、两种后端等。稠密壁速与 Gauss 最大值的 7/3 包络来自二次插值基函数界，用于抓漏点，附加 1e−12Uref 仅是读回算术余量；单格分辨率误差门不代表工程验收。原两档网格和相位试验的限制见当前状态，全量软件回归已由集成另行覆盖。

没有共形流体 polyMesh，外部 checkMesh 不适用，独立读取器不替代 Cut-cell 门。参考：[Brinkman 惩罚](https://www.math.u-bordeaux.fr/~chabrune/publi/ABF-NM.pdf)、[惩罚压力投影](https://arxiv.org/abs/2306.06277)、[Taira–Colonius 原作者论文目录](https://www.seas.ucla.edu/fluidflow/pubs.html)（DOI 10.1016/j.jcp.2007.03.005）；本实现是有限顺应的双线性表面惩罚，不声称复现整套算法或资格。

## 验证与证据

以下脚本位于 `tools/verification/`，参数和预算查各自 `--help`：

| 检查 | 入口 |
| --- | --- |
| 背景 / CM2D / OpenFOAM 拓扑 | `check_background_grid.py`、`check_openfoam2d.py`、`check_hybrid_mesh2d.py` |
| 方向连通与挤出 | `check_directional_connectivity.py`、`check_extruded_quality.py` |
| 原生方程 / 圆环 / 开口 / 对称 | `verify_native_flow.py`、`verify_native_fv.py`、`verify_rotating_annulus.py`、`verify_pressure_openings.py`、`verify_symmetry_flow.py` |
| 非定常 / 制造解 | `verify_transient_flow.py`、`run_manufactured_flow.py`、`compare_transient_steps.py` |
| 温度与规模 | `verify_thermal_flow.py`、`verify_thermal_time.py`、`verify_thermal_scale.py` |
| 可压及独立算子 | `verify_euler.py`、`verify_heat_conduction.py`、`verify_viscous_stress.py`、`verify_wall_gradient.py` |
| 初值 / 完整性能 | `verify_flow_initialization.py`、`benchmark_flow_pair.py`、`benchmark_laminar.py` |
| 打包 App / 浸入边界 | `desktop_platform_smoke.py`、`verify_immersed_flow.py` |

Python 独立读回、Solver 门、实际 OpenFOAM 标准/扩展 checkMesh、方程接受和物理精度分别记录；不互相替代。可视化从真实输出绘制。测试和研究脚本不是缓存，不因未列入 CMake 删除。

`prolongate_regular_flow.py` 只映射完整矩形张量格；`extract_steady_iterate.py` 提取同网格初值并存哈希。曲壁映射仍在 `outputs/laminar-performance/` 研究流程，未进入 App。性能比较含粗解/映射/细解，并使用相同控制/精度。

已有 gzip 场按 `compressed*.json` / `*-compressed.json` 原 SHA 核对后按需恢复，保留压缩包；`*.canonical-fields.json` 指向逐字节核验的规范场，勿批量展开历史。旧二进制作基线时检查实际 SHA，当前 runtime 不能代替旧版。

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
