# 开发导航

## 完整笛卡尔背景网格

`cartmesh2d_cli --background-grid adaptive|uniform` 在已有几何诊断、Quadtree细化及2:1平衡后直接导出完整叶子，不执行Cut-cell、合并、Solver修复或OpenFOAM输出。均匀模式令全域最低层级等于最高层级；自适应模式复用距离带、盒加密及尺寸场。当前资源上限分别为level 10和12，不是已测性能保证。

```sh
build/cartmesh2d_cli examples/acceptance/circle.xy outputs/background-grid/circle 7 0.5 0.1 exterior - 3 --background-grid adaptive
python3 tools/verification/check_background_grid.py outputs/background-grid/circle.background.json
```

`.background.json` 为 `cartmesh2d-background-v1`，包含域、原始边界、完整单元包围盒/层级/整数格坐标及分类（0外部、1内部、2相交）；`solver_ready=false`。`.background.vtk` 保存相同完整四边形及分类/层级，可由ParaView读取。自适应粗细交界可有悬挂节点，VTK按独立四边形展示，不宣称具有统一求解面拓扑。JSON不是CM2D求解网格，不能输入现有流体求解器。

实现入口：`src/io/BackgroundGridIO2D.cpp`；独立审核核对整数格覆盖/无重叠、几何分类及2:1平衡；`cartmesh2d_background_grid`还核对JSON/VTK一致、重复确定性、非法几何及误用导出的拒绝。桌面提供独立背景模式入口，使用 `desktop/src/core/background-grid.js` 读取并禁止直接送入流体求解。桌面均匀最高level9、自适应最高level10且全域最低最高level8；这是交互资源上限。

## 分支与里程碑

`mesher-v0.3.0` 是固定标签（`691c97e`）；`codex/mesh-maintenance` 是可继续修改的网格维护线；`codex/cfd-development` 是包含网格核心的 CFD 开发线；`main` 是包含已验证网格和CFD功能的集成线，当前集成状态见CURRENT_STATE。不要移动已有里程碑标签去“更新版本”，后续里程碑另建标签。

纯网格修复先在维护分支保留最小回归，验证后合并到 CFD 开发线；合并冲突需人工核对，随后重跑受影响的网格和求解器检查，再集成到 main。跨分支改动不会自动同步。切换分支前处理当前工作区改动，切换后重新构建原生工具及 desktop runtime，避免旧运行包与当前源码混用。继续遵守唯一根目录和固定文档入口规则。

## 从哪里进入代码

生成链：输入轮廓 → 尺寸场/Quadtree → 真实流体 Cut-cell 或共形边界层 → 小单元处理 → 求解拓扑 → 质量 → 文件导出。

| 路径 | 责任与入口 |
|---|---|
| `apps/cartmesh2d_cli.cpp` | 纯 Cut-cell 总流程、尺寸场参数、物理面积门、Solver 质量和输出 |
| `apps/cartmesh2d_hybrid_cli.cpp` | 边界层总流程、局部修复开关、fallback 与导出 |
| `apps/cartmesh2d_flow_cli.cpp` | 原生稳态层流 CLI；SIMPLE / Rhie–Chow，桌面调用与诊断场导出 |
| `apps/cartmesh2d_transport_cli.cpp`、`fv/ScalarTransport2D`、`fv/ThermalFlow2D`、`fv/ThermalCheckpoint2D` | 守恒标量/恒物性温度；共享通量、混合边界、同步后向欧拉、联合状态续算与导出 |
| `apps/cartmesh2d_fv_cli.cpp` | 自研二维标量扩散/泊松 CLI；读取最终 solver.cm2d，输出场、通量、残差和误差 |
| `fv/FvMesh2D`、`fv/Diffusion2D` | 最终多边形几何缓存、共享边通量、加权最小二乘梯度、非正交扩散、Jacobi-PCG |
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
| `quality/SolverTopology2D`、`PatchLocalQuality2D` | 凸划分、源邻域相对评分、精确合并/切分及单壁面一变二批量，完整质量验收 |
| `io/Dxf2D`、`BoundaryMetadata2D` | CAD 曲线、单位、边界名称/角色 |
| `io/MeshIO2D`、`OpenFoam2D` | CM2D/VTK/JSON、二维挤出和 OpenFOAM case |

表中的模块分别位于 `include/cartmesh2d/` 和 `src/`，仍参与 CMake 构建，未把现用模块当旧版本删除。

桌面入口是 `desktop/src/main.js`（编排/IPC）、`preload.js`（受限桥接）。
`core/` 管能力、样例、几何/SVG 输入、参数、尺寸预算、CM2D、结果摘要；`process.js` 负责可取消/限时的原生子进程，`automatic.js` 支持未指定目标数量时最多 4 组候选；`cell-budget.js` 负责数量初值与实测调整，`budget-runner.js` 最多 3 次执行并保留最接近目标的成功结果。

图片输入：`renderer/raster-import.js/css` 提供预览与强制标定，`raster-worker.js` 在 Worker 中调用共享 `core/raster.js`（RGBA 分割、四连通区域、边界追踪、保孔简化）。`core/raster-geometry.js` 检查编码尺寸、复核环交叉并把像素转换为米；main 的 `read-raster` 固定源文件快照，`commit-raster` 保存用户确认的 XY 和图片记录，再走已有生成链。原图不上传，不依赖外部模型。

`desktop/tests/raster*.test.js` 覆盖阈值端点、保孔/填孔、主体选择、复杂度拒绝、标定和坏图片头。实际 App 可使用 `--smoke=circle --image=/absolute/input.png --image-width=200 --out=/absolute/output --export=/absolute/result.zip` 验证 200 mm 标定至导出；`--image-fill=true` 填孔，`--image-preview-only=true` 检查未标定禁用确认，`--image-reject=true` 检查无轮廓拒绝与取消保留旧几何。
`renderer/app.js` 管表单和状态，`viewport.js` 管真实网格画布；HTML/`style.css` 管布局；`theme.css` 集中两套外观变量，`themes.js` 管主题切换与本地持久化；画布配色仍在 `viewport.js`，与界面主题独立。没有第二套前端。

跨平台中文：`assets/fonts/cartmesh-ui-regular.woff` 自带界面字体，来源为 Noto Sans CJK SC 2.004，覆盖 GB2312 与当前界面字符；完整 OFL 与原文件哈希在同目录 LICENSE.txt。开发者可用 `desktop/scripts/subset-ui-font.py` 和 fontTools 从原始 OTF 重建，正常构建不下载字体。生僻自定义字符仍由系统回退字体补充；导出 Canvas 在字体加载完成后绘制，CI 检查实际字体加载并保留截图。

## 构建与有限验证

macOS 使用 `/usr/bin/clang++`，不要使用 PATH 中的 mesasdk 编译器打包。根目录只保留一个 `build/`。

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/clang++
cmake --build build -j2
ctest --test-dir build --output-on-failure
npm test --prefix desktop
```

三个平台的桌面本机构建入口：

```sh
npm ci --prefix desktop
npm --prefix desktop run build:native
cmake --build build --config Release --parallel 2
ctest --test-dir build -C Release --output-on-failure
npm --prefix desktop test
npm --prefix desktop run pack:mac
# Linux 最后一行换成 pack:linux；Windows 换成 pack:win。
```

构建依赖 Node.js 22、CMake 3.20+、Python 3 和 C++20 编译器；Windows 使用 Visual Studio 2022 C++ 桌面工作负载（含 Windows SDK），Linux 使用 GCC 11+ 或兼容 Clang。运行用户无需安装这些工具。`build-native.js` 准备本机六份 CLI、11 个样例和 runtime manifest；打包前检查实际 PE/ELF/Mach-O 格式与架构，拒绝跨系统或错架构混装。macOS 固定系统 clang 并检查动态库；Windows 使用静态 CRT、UTF-8 编码/路径 manifest。ZIP 导出改用流式 `yazl`，不依赖 ditto/zip/PowerShell 命令。

`desktop-platforms.yml` 分别在 macOS 14 arm64、Ubuntu 22.04 x64、Windows Server 2022 x64 执行完整构建/测试。`tools/verification/desktop_platform_smoke.py` 启动打包 App，使用含中文和空格的临时/输出路径，运行 PNG、JPG、hybrid、均匀背景和自适应背景五例，解压 ZIP 后独立核对 OpenFOAM 单元、图片和标定。Linux 的 `--no-sandbox` 仅用于隔离 CI 虚拟显示测试；产品启动入口不附加该参数。通过后上传版本运行包与证据；不把独立读取器当成真实 checkMesh。

日常用 `ctest --test-dir build -R <相关测试名> --output-on-failure`；完整交付才统一全量运行。
测试目录包含解析、拓扑、守恒、质量、尺寸场、真实输出和确定性检查，不是可删除的缓存。

真实桌面流程（先打包准备 runtime；输出目录预先创建）：

```sh
mkdir -p outputs/smoke
cd desktop
node_modules/.bin/electron . --smoke=circle --out=../outputs/smoke --shot=../outputs/smoke/circle.png
```

可加 `--verified-preset=true`（真实界面载入已有高密案例参数）、`--density=150000`、`--target-cells=100000`、`--auto-padding=0.5`、`--theme=duet`、`--interaction-check=true`（主题状态与实际参数转手动）、`--allow-unsafe=true`、`--control=manual`、`--density=dense`、`--method=hybrid`、`--mode=light`、`--regions=1`、`--repeat=1`。省略 `--out` 即验证默认临时预览；`--export=/绝对路径/result.zip` 验证结果包。每轮会缩到最小窗口并检查底部可达、预览/导出单元数。smoke 会走真实表单/IPC/CLI 后退出；它不等于所有界面操作都已验收。

## 原生层流求解

入口为`apps/cartmesh2d_flow_cli.cpp`；核心`src/fv/Incompressible2D.cpp`使用SIMPLE、Rhie–Chow、共享面压力和黏性通量。必须读取最终`*.solver.cm2d`，背景JSON和诊断失败网格不能代替。完整参数用`build/cartmesh2d_flow_cli --help`查看。

建议从现有网格和明确边界开始；以下仅示范命令结构，几何须通过对应模板约束：

```sh
build/cartmesh2d_flow_cli --mesh outputs/channel.solver.cm2d --case channel --export-boundaries outputs/channel.boundaries
build/cartmesh2d_flow_cli --mesh outputs/channel.solver.cm2d --case custom --boundary outputs/channel.boundaries --output outputs/channel-flow --nu 0.1 --speed 1 --max-iterations 5000 --tolerance 1e-9
```

边界预设还包括方腔、喷管及旋转圆环；不满足模板几何条件时明确失败。逐面命名文件由原生程序导出再修改，保持面编号/几何绑定。`pressure-opening`与`symmetry`仅支持轴对齐面，不能将倾斜面悄悄投影。压力单位为p/rho（m²/s²）；命名边界流量为每单位厚度m²/s，向外为正。

| 控制 | 约束 |
| --- | --- |
| `--convection upwind` / `limited-linear` / `face-limited-linear` | 分别为迎风、有限重构及面方向限制；不是物理模型切换 |
| `--pressure-preconditioner ic0` / `aggregation` / `cholesky` | Cholesky限macOS系统后端；真实线性残差仍检查 |
| `--linear-policy strict` / `adaptive` | 自适应最终必须严格复核，不能放宽最终停止门 |
| `--velocity-relaxation` | 默认.6；稳态/非定常独立层流可显式设置，增大不保证更快 |
| `--pressure-corrections` | 默认4；困难曲壁可能不稳，不能用连续性小代替收敛 |
| `--steady-acceleration anderson` | 可选、默认none；不支持非定常或材料更新；候选须降低原残差并过守恒门 |
| `--initial-guess` / `--initial-flux` | 稳态同目标网格的初始迭代，面初值须与单元初值同用；不等于已接受检查点 |

原方程动量残差、局部/全局守恒、u/p变化和线性解精度分别检查。稳态面通量包含与单元松弛一致的旧通量缺陷项；不能为了提速删除它。macOS流动CLI在首次Accelerate调用前设置单线程，保证既有确定性；独立研究入口也要固定`VECLIB_MAXIMUM_THREADS=1`。

### 非定常层流与断点续算

已有后向欧拉、固定/自动CFL步长、拒绝试步及重试。每个接受步通过原停止门才推进物理时间；失败或取消保存上个完整状态。用CLI `--help`核对时间和重启参数；输出checkpoint绑定网格/物性/边界/格式，CSV或`.tmp`不能充当重启状态。App交互见DESKTOP_APP。

`flow.time-history.csv`记录本次接受时间步，`flow.residuals.csv`对应最后尝试的内迭代。固定dt的CFL仅是测量，不是时间精度证明；自动dt也不替代多时间步对照。同步温度必须使用联合thermal checkpoint，不能只载入carrier状态。

### 工况、边界与桌面

`desktop/src/core/flow-case.js`保存网格哈希绑定的工况；当前v4含自适应线性和速度松弛，兼容旧格式的明确默认值。求解结果实际设置必须与请求一致。`main.js`/`core/thermal-job.js`编排流动与温度的执行、取消与失败保留；`flow-checkpoint.js`审核重启。

SmoothMovingWall保留平滑壁速梯度，普通恒定壁迹具有不同语义；旋转圆环模板使用规则多边形法向构造切向速度。不得自动把用户折线改成解析圆，也不能用研究细圆环的结果替代原输入资格。

### 守恒标量与恒物性热输运

入口`cartmesh2d_transport_cli --help`；支持冻结流场或流动/温度同步推进、定值与绝热/指定通量。D=k/(rho cp)，热源Q/(rho cp)，向外通量q/(rho cp)。温度不反馈流动，无浮力、辐射或共轭传热。桌面温度数值控制和独立层流不同，不能把后者的所有优化选项强行传入联算。

### 实验求解器与保留代码

`cartmesh2d_fv_cli --help`是独立扩散/泊松入口。SST输运/壁距/稳态RANS代码位于`src/fv/Sst*`及对应测试；保留现有行为，不作为通用湍流资格。

### 隔离可压 Euler 开发

`codex/compressible-flow` 从 main 分出，保留默认 Rusanov / 一阶；`--flux hllc --order 2` 启用压力感知 HLLC/HLLE、受限线性重构与 SSPRK(2,2)。求解理想气体质量、两分量动量及总能量；可选恒系数 Fourier 导热及 Newtonian 黏性应力/做功，二者均进入总能量；尚无湍流或变物性。

- `EulerFlux2D.cpp` 以真实面面积向量计算一个共享守恒通量；HLLC 使用包含左右声学锥的 Roe/Einfeldt 估计，检查两侧星状态，失效时显式计数并用 Rusanov。
- 每个相邻单元全部邻面上的 `min(pL,pR)/max(pL,pR)` 取三次方，再取最小值作为接触恢复权重 ω；通量为 `F_HLLE + ω (F_HLLC - F_HLLE)`。传感器依据 [Simon 与 Mandal，式49–50、α=3](https://arxiv.org/pdf/1803.04954)，本实现采用多边形邻面模板并混合四个分量，**不冒称论文的选择性 HLLC-ADC 原样复现**。未加保护的最小旋转 Sod 例曾放大舍入扰动，测试中保留该场景。
- `Euler2D.cpp` 用距离加权最小二乘及面值范围限制重构原始量；速度分量在随网格旋转的局部坐标系中限幅。矩阵近奇异或面值因浮点消减失去正性时明确计数并退为常量，绝不对接受场的密度、压力或能量加下限。
- 固壁先求镜像 Riemann 压力牵引，再按对称性构造严格零质量/对流能量、纯法向动量通量；避免 SI 总能量大数消减产生假泄漏。独立脚本使用另一个标量壁压公式核对。
- 二阶每个阶段都使用真实边界或周期平移；两个前向 Euler 候选及最终状态都检查正性。导出的面通量是两个阶段的均值，波速取两阶段最大值；CFL 不满足或候选失效则缩步重试。失败候选不写入接受检查点。
- `verify_euler.py` 从原始多边形和前一守恒场独立重建斜率、两阶段、HLLC/HLLE 通量、声速及状态方程，同时核对历史终态积分。它审计最后接受步和完整时间历史，不代表逐步独立重放整个求解。

HLLC 星状态见 [Toro 的推导](https://www.prague-sum.com/download/2012/Toro_2-HLLC-RiemannSolver.pdf)；受限重构参考 [Barth–Jespersen](https://ntrs.nasa.gov/citations/19890037939)。`--order 2` 表示在光滑区使用二阶格式，激波及局部退阶处不宣称二阶精度。

```sh
build/cartmesh2d_euler_cli --mesh final.solver.cm2d --output outputs/euler/run --case sod --gas-r 1 --end-time .2 --flux hllc --order 2
ctest --test-dir build -R '^cartmesh2d_euler_' --output-on-failure
python3 tests/euler_accuracy_cli_test.py --cli build/cartmesh2d_euler_cli --output outputs/euler-development/qualification
```

通量方程复核按面长乘特征通量归一化：质量 `ρV`、动量 `ρV²`、能量 `(ρE+p)V`，其中 `V=max(|速度|+声速)`，两侧/两阶段取包络。容差为512个 binary64 机器 epsilon，用于独立公式/重构的浮点一致性，不是物理精度门；避免对接近零的 SI 通量使用无量纲固定绝对误差。固壁质量/对流能量要求精确为零，总能量通量等于导热与向外黏性做功之和；静止壁面做功为零。单位缩放回归把 SI 与无量纲的完整场归一化比较，另增两个小 Sod 和旋转闭壁例，成本为秒级。

开发验证按量分别判断：沿用无量纲逐格守恒门 `1e-12`；光滑熵波及涡旋比较面积加权 L1 误差，二倍细化的观测阶数回归下限为1.4，允许限制器在极值附近降阶。Sod 以单元中心采样精确 Riemann 解，比较同网格/同目标时间的面积加权离散 L1 误差，不把激波误差称为二阶收敛。Mach 3/10 法向激波加入相对密度幅度 `1e-5` 的交替扰动，检查在对流距离0.3的窗口内，横向密度跨度相对波后密度不超过输入跨度的10倍；不外推为任意强激波稳定性证明。测试规模最高128×128格，精度脚本通常为分钟内至数分钟级，取决于并行负载。

μ=k=0 时检查点保持 v1，μ=0、k>0 时为 v2，μ>0 使用同时绑定输运物性及壁速的 v3；网格/物性/边界仍严格绑定；续算可显式改变通量、阶数、CFL和目标时间。使用相同控制参数中断/续算时检查点逐字节一致。摘要/历史分别报告接受阶段中的 HLLC 正性回退、重构退阶和最小接触恢复权重；`hllcFallbackStages` 的位0/1对应第一/第二阶段，周期配对只计一次。


### 总能量耦合导热

物理闭合为 `T=p/(ρR)`、`cv=R/(γ−1)`、`q=−k∇T`。总能量方程新增 `−∇·q`，质量和动量没有直接热源；压力变化会在后续流动阶段产生反馈。[NASA 的比热说明](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-heats-cp-and-cv-1/)区分内能的 cv 与焓的 cp；这里定密度阶段的体积热容为 **ρcv**，不能把被动温度入口的 `k/(ρcp)` 搬过来。[NASA Wind 方程说明](https://www.grc.nasa.gov/WWW/wind/TFAWS2007/Formulation.pdf)作为理想气体/输运变量参考；本实现 μ 与 k 分别指定为常数，不强制 Prandtl 闭合。

- `HeatConduction2D.hpp/.cpp` 是独立热算子，输入温度与体积热容，输出每条真实面的共享热通量、单元残差和 1/s 的显式速率。`EulerTransport2D` 与 `IdealGas2D` 分开；`EulerStepper2D` 持有网格/物性的不可变快照，缓存几何、梯度权重和热算子系数行范数，避免每个时间步重建。原 `advanceEuler2D` 接口继续可用。
- 法向梯度采用 `S=τd+C`，`τ=|S|²/(S·d)`，热通量 `Qf=−k[τ(TN−TP)+gf·C]`。`gf` 来自距离加权最小二乘梯度的法向距离插值；这类非正交拆分可对照 [OpenFOAM corrected snGrad](https://doc.openfoam.com/2212/tools/processing/numerics/schemes/sngrad/rtm/corrected/)。Dirichlet 使用真实面心与壁温；Neumann 以 **n·∇T=−qout/k** 进入梯度约束，不能把倾斜的“格心到面心方向”当成法线。内部面只算一次、邻格反号；周期边界用平移后的格心与严格反号配对。
- 边界默认绝热，也可指定 Kelvin 或向外 W/m²；后者负值加热。开边界默认零 Fourier 通量，不取消对流能量输运。`sealed` 是全部物理边界为壁面，可选自由滑移或静止无滑移；`custom` 的 `CM2D_EULER_BOUNDARY 2` 在每行旧格式之后增加 `insulated|temperature|flux value`，仍可读 v1。
- 对温度残差的完整线性系数矩阵 A，`h_i=Σ_j|A_ij|/(2 Vi ρi cv)`，时间步满足 `dt·(Σf af|Sf|/Vi+h_i)≤CFL`。**非正交修正项也计入 A**，每阶段按密度更新热容量，两阶段取速率包络。这是全算子行范数限制，不是任意网格的单调性或谱稳定性证明；`heatNonMonotoneRows` 报告不满足 M 矩阵符号的行数。当前修正算子保留线性精确性，可能不是单调格式；逐阶段正性失败时缩步或明确终止，绝不修剪温度/能量。均匀网格上的二阶结果也不能外推为任意 Cut-cell 上二阶，见 [OpenFOAM 非结构网格精度说明](https://openfoam.org/release/2-3-0/numerics/)。
- 面热通量随 SSPRK2 做阶段平均并加入总能量通量，故闭壁严格零质量；静止壁面总能量通量恰等于导热，切向移动壁额外包含黏性做功。`heatRate`、`heatFlux`、`convectiveEnergy`、`thermalCourant`、`combinedCourant`、`boundaryHeat` 可独立读回。`verify_heat_conduction.py` 从原始多边形重建热算子；`verify_euler.py` 重算最后步两个阶段，并核对所有相邻历史步的质量/能量收支。历史收支是记录一致性检查，不能冒称独立求解全部轨迹。

导热验证保持范围明确：旋转/扭曲网格上的线性温度与混合热边界检查面通量线性精确性（按 `k|∇T||Sf|` 归一化，1024 epsilon，约百格）；通过对每格温度加减1 K、直接差分实际残差，检查完整系数行范数（2048 epsilon）。这些是浮点代数一致性门，不是工程精度要求。独立通量审计的能量尺度增加 Fourier 仿射行绝对值包络，仍使用512 epsilon；闭腔累计能量收支按初末总能量和输入热量之和归一化，门为1024 epsilon。

纯热方程周期正弦衰减采用精确单元平均初值，16/32/64格验证空间与 RK2 联合误差。完整 Euler–Fourier 小扰动使用相对温度幅度 `ε=1e-5`，初始 `δρ=−ρ0 ε cos(wx)`、`u=0`、`δT=T0 ε cos(wx)`；线性化连续方程的模态系数满足：

```text
A' = −ρ0 w B
B' = (R T0/ρ0) w A + R w C
C' = −(γ−1) T0 w B − [k/(ρ0 cv)] w² C
```

独立 Python 以矩阵指数求该常系数连续系统，完全不使用原生离散系数。误差按面积平均、除以初始扰动幅度；中心采样初始化/参考的空间误差和非线性 `O(ε²)` 误差都保留。16/32/64的局部回归要求最终观测阶大于1.5，区分二阶行为与退化为一阶；不作强非线性或通用网格精度承诺。该脚本同时检查所有四种对流格式/时间阶数组合、SI量纲相似性、定温/热流/绝热/混合条件、物理续算绑定与篡改拒绝，通常秒级运行。实际 App 的 `k=100 W/(m·K), Tw=400 K` 是显式标定的功能验证参数，不代表空气材料数据库。

```sh
build/cartmesh2d_euler_cli --mesh final.solver.cm2d --output outputs/heat/run --case sealed --conductivity 100 --wall-thermal temperature --wall-value 400 --density 1.225 --pressure 101325 --gas-r 287.05 --end-time .002 --flux hllc --order 2
python3 tests/euler_conduction_cli_test.py --cli build/cartmesh2d_euler_cli --output outputs/euler-conduction/verification
python3 tools/verification/verify_euler.py --mesh final.solver.cm2d --prefix outputs/heat/run
```

v2 检查点将 k 和全部逐面热条件写入完整绑定；改变 k、热壁类型或热边界数值会拒绝续算。k=0 的旧 v1 检查点和数值默认保留，实测原二进制与新二进制的 Rusanov/一阶和 HLLC/二阶 Sod 检查点逐字节一致。模型目前没有温变/各向异性导热、辐射、固体共轭传热或隐式导热；这些需要独立方程、界面守恒与验证，不通过增加界面开关冒充完成。


## 验证与证据

| 目的 | 入口 |
| --- | --- |
| 背景覆盖/分类/2:1 | `tools/verification/check_background_grid.py` |
| CM2D/OpenFOAM独立拓扑 | `check_openfoam2d.py`、`check_hybrid_mesh2d.py` |
| Solver以外的方向连通/挤出检查 | `check_directional_connectivity.py`、`check_extruded_quality.py` |
| 原生离散方程与输出 | `verify_native_flow.py`、`verify_native_fv.py` |
| 圆环/压差开口/对称 | `verify_rotating_annulus.py`、`verify_pressure_openings.py`、`verify_symmetry_flow.py` |
| 非定常/制造解 | `verify_transient_flow.py`、`run_manufactured_flow.py`、`compare_transient_steps.py` |
| 温度与规模 | `verify_thermal_flow.py`、`verify_thermal_time.py`、`verify_thermal_scale.py` |
| 初值与性能对照 | `verify_flow_initialization.py`、`benchmark_flow_pair.py`、`benchmark_laminar.py` |
| 打包App | `desktop_platform_smoke.py` |

上表脚本除特别注明均在`tools/verification/`，用各自`--help`查看输入和预算。测试、研究脚本和可视化脚本不是缓存；仅因未在CMake列出就删除它们会丢失独立复现入口。`tools/visualization/`从实际输出作图，不能用示意图替代网格或场。

真实OpenFOAM验证使用目标版本的`checkMesh`。默认检查与`-allGeometry -allTopology`结果分别报告；Python读回、原生拓扑、质量门与CFD接受不能互相替代。当前CI另包含真实OpenFOAM标准检查与细化回归，不能为让main变绿而删除失败项。

### 稳态初值与性能复现

`prolongate_regular_flow.py`只支持完整矩形张量网格的粗细映射；不能冒充通用Cut-cell映射。`extract_steady_iterate.py`提取同网格场与面通量并记录哈希。曲壁喷管映射仍为`outputs/laminar-performance/`下研究脚本，未作为桌面自动能力发布。

完整基准参数、源二进制、网格/场SHA和实际粗解/映射/细解成本在`artifacts/current/native-laminar-*.json`；当前关键索引见CURRENT_STATE。比较同网格、同控制、同精度及总成本，不只比较单轮或最快一段。新实验使用新输出目录，保留原失败。

大文件如已gzip，按相应`compressed*.json`或`*-compressed.json`核对原SHA后只恢复所需文件，保留压缩文件；`*.canonical-fields.json`指向已逐字节核验的规范场。不要批量展开全部历史。已有包内旧二进制可作基线，但必须核对实际SHA，当前runtime不能充当旧版。

## 不能丢失的设计边界

- 原生二维独立于三维仓库，二维基本对象和计算域语义保持；不把3D压到z=0。
- 原始XY折线不平滑。真实SVG/DXF曲线采样与用户折线修改是不同操作。
- 共享格点/交点构造须保留稳定身份和EmbeddedBoundary来源，不能每叶单独补洞或仅凭近邻坐标猜物面。
- 输入端点算术舍入预算与几何容差分开；极短但可解析的支撑不能折叠。
- 不删除坏格、相交格或小面积流体以过门；局部孔洞等未支持情形明确失败。
- 边界层、过渡和余域统一共形拓扑；降层及纯Cut-cell fallback明确报告。
- Q1已取消；Q3/Q4/Q5及patch事务仍参与构建/CI，不能把名称里的阶段数字当作废弃依据。
- 默认质量门、物理精度及确定性要求不因清理或加速放宽；跨编译器逐字节一致另证。

## 历史与工作区

只用根目录源码、`build/`构建和忽略提交的`outputs/`实验。当前文档固定三份，旧阶段描述从Git读取：

```sh
git show ac24bee:docs/CURRENT_STATE_CN.md
git show ac24bee:docs/DEVELOPMENT_CN.md
git log --all --oneline -- src/quality/SolverTopology2D.cpp
```

`mesher-v0.3.0`是固定网格里程碑。旧`archive/*`标签保存历史方案，不是当前验收。历史截图与JSON保留原版本适用范围；不再复制成新的阶段文档或完整源码副本。


### 可压黏性应力与机械功

`ViscousStress2D.hpp/.cpp` 是独立速度输运算子，输入二维网格、速度、密度和动力黏度 **μ（Pa·s）**；不能把不可压入口的运动黏度 ν 当作 μ。使用平面理想气体的 Stokes 假设：

```text
τ = μ [∇u + (∇u)^T − (2/3)(∇·u) I]
F_viscous = [0, −(τ·Sf)x, −(τ·Sf)y, −u_face·(τ·Sf)]
F_energy = F_convective_energy + F_viscous_work + q·Sf
```

`2/3` 是气体应力本构系数，不随二维网格改成 `1`；因此纵向黏性扩散系数为 `4μ/(3ρ)`，横向剪切为 `μ/ρ`。[NASA Wind 方程说明第29–31页](https://www.grc.nasa.gov/WWW/wind/TFAWS2007/Formulation.pdf)给出该本构与总能量应力做功。几何、面拓扑和未知量完全原生二维，没有调用三维核心。总能量已经包含应力做功，不能再次添加体积 `τ:∇u`，否则重复计入能量；动能与内能间转换由同一组守恒更新产生。

- 距离加权最小二乘重建完整速度梯度；内面插值梯度后沿法向修正，使 `Gf·(CN−CP)=uN−uP`。因此保留交叉导数、法向应变和非正交部分。面速度以两侧梯度外推到真实面中心再插值，仿射场在偏斜面上仍精确。每条内部面计算一次，周期面对复制严格反号的动量与功；两阶段分别计算后平均。
- `NoSlipWall` 只允许静止网格上的切向壁速；CLI 预设及 App 为静止壁，`custom` 文件可给逐面切向速度。对流部分仍使用不可穿透壁的镜像压力牵引，黏性部分使用指定壁速的 Dirichlet 梯度。静止绝热壁的质量/总能量通量精确为零；移动壁的功使用给定壁速，不用流体单元中心速度代替。带法向速度的壁面会拒绝，不能冒充移动网格。
- 自由滑移壁的法向速度为零，切向黏性牵引投影为零；梯度模板使用局部平面壁的法向速度约束及切向零法向导数。开边界采用**零黏性牵引**，梯度延拓使用零法向速度导数；这是一项明确的边界模型，不能当作任意短出口的充分外流条件。
- 稀疏准备阶段组装完整 `2N × 2N` 速度到动量残差 Jacobian A，包括交叉导数、非正交修正和周期耦合。每格黏性速率 `max_k Σ_j |A_(2i+k,j)| / (2ρ_i V_i)` 的单位为 `1/s`，与声学、Fourier 速率共同限制显式步长；每个 RK 阶段用当时密度重算速率。该行范数是显式扩散控制量，**不是完整非线性能量 Jacobian 或任意网格的稳定/熵证明**；两个前向阶段与最终候选继续检查正密度和正内能，失效缩步，不裁剪。
- `EulerStepper2D` 拥有不可变几何、物性和边界，缓存热/黏性算子；内层只求场梯度与通量。μ=0 不创建黏性算子并保留旧推进路径。v3 检查点同时绑定 μ、k、全部壁面模型及切向壁速；v1/v2仍按原字节格式生成/读取。边界文件v3在v2热条件后增加 `wallUx wallUy`，旧v1/v2继续可读。`custom` 必须在边界文件指定条件，禁止用预设开关覆盖。
- `.faces.csv` 保存 `viscousMomentumX/Y`、`viscousWork`、`convectiveMomentumX/Y`、`convectiveEnergy` 与 `heatFlux`；`.cells.csv` 保存 `viscousRate`；历史保存 `viscousCourant` 和向外 `boundaryViscousWork`。负壁功代表机械能输入。`verify_viscous_stress.py` 从原始多边形建立稀疏线性表达式，并独立求完整动量 Jacobian；`verify_euler.py` 重建最后步两阶段的黏性/热/对流通量并核对全部相邻历史的质量与能量。

直接算子测试使用32格扭曲网格及旋转网格的平移、刚体转动、均匀膨胀、简单剪切和一般仿射速度。应力/功及实际残差扰动得到的行范数使用1024 epsilon的归一化浮点检查；不是工程流场精度要求。独立面通量审计保持512 epsilon，增加黏性线性表达式绝对值包络作为量纲尺度。SI相似性同时缩放 μ、k、压力、速度、时间和长度，比较归一化完整场；允许4096 epsilon累积浮点差异。

连续解验证包含16/32/64格周期横向剪切波 `v=A exp(−μw²t/ρ) cos(wx)`（A=1e−4声速），及原热模态线性系统的速度行增加 `−4μw²B/(3ρ)`（相对温度幅度1e−5）。独立矩阵指数不读取求解器系数。Couette 使用周期x与上下切向壁面，μ=k=0.1、U=0.5、H=1、Tw=1，解析稳态为 `u=Uy/H`、`T=Tw+μU²/(2k)·(y/H)(1−y/H)`、`p=1`、`ρ=p/(RT)`；从解析场开始推进至t=0.1，在8/16/32层比较温度误差及机械功/热量收支，**这是稳态离散一致性检验，不是从静止收敛验收**。观测阶下限1.7用于防止光滑问题退回一阶，允许限制器与壁面局部误差；整个CLI验证通常数十秒。任意Cut-cell边界层、强激波/黏性相互作用、壁面摩擦/热流绝对精度和长期稳定性尚未验收。

```sh
ctest --test-dir build -R 'cartmesh2d_(viscous_core|euler_viscosity)' --output-on-failure
python3 tests/euler_viscosity_cli_test.py --cli build/cartmesh2d_euler_cli --output outputs/euler-viscosity/validation
build/cartmesh2d_euler_cli --mesh final.solver.cm2d --output outputs/euler/run --case external --viscosity .02 --wall-model no-slip --conductivity 100 --wall-thermal temperature --wall-value 400 --flux hllc --order 2 --density 1.225 --pressure 101325 --u 50 --end-time .00005
```

最后一条为便于短时验证的测试物性，不是空气推荐值；μ、k、γ、R 必须按实际问题指定。默认μ=k=0、滑移、Rusanov一阶均未自动更改。


### 可压壁面精度与时间细化

`--wall-gradient quadratic` / App“壁面热流 / 应力梯度 → 二次重构”是可选数值格式。默认仍为 `linear`；适用于定温壁的 Fourier 热流和无滑移壁的完整黏性应力。内面、指定热流、绝热、滑移和开边界保留原模型。该选项允许随检查点显式更改，物理边界和物性仍严格绑定。缓存求解器在构造时选择壁面格式，推进时拒绝与缓存不同的设置。

`WallGradient2D` 只依赖原生二维有限体积几何。以真实壁面中心的指定值固定常数项，在法向/切向局部坐标拟合五项 `n, t, n²/2, nt, t²/2`。沿共形流体 owner/neighbour 图扩展模板，周期边界使用平移后的像点，绝不跨越物理边界。先取两层邻域，必要时最多五层；超过512个像单元或仍缺秩时显式失败，不静默降阶。法向和切向分别按局部跨度无量纲化，再进行列归一化、列主元 Householder QR；不用正规方程。尺度化矩阵残余列范数须大于 `sqrt(machine epsilon)`，这是保留约半数双精度有效位的数值秩保护，不是物理精度阈值。准备时只缓存梯度系数，运行时对“样本值减目标壁值”求和，常量在加权前消去。

二次最小二乘和QR的基本做法可参考 [NASA White/Nishikawa，D.2.3节](https://ntrs.nasa.gov/api/citations/20210024196/downloads/white_and_nishikawa_afang_paper_v_1.7.pdf)。本实现是独立推导的二维壁面值约束模板，不是该论文的三维F-ANG移植。这里输入的是现有二阶方法中的**质心原始变量点值近似**，没有将保守量单元平均值变成三阶多项式，也没有高阶面积/面通量积分；**局部二次重构不等于全局三阶**。对流仍为所选一阶/二阶，光滑问题按整体二阶验收；非光滑角点和激波附近不能使用光滑阶数推断精度。

热流和应力直接用恢复的壁面梯度；壁功始终使用指定壁速。所有新增非局部系数进入热残差和完整动量块的行范数，参与两阶段组合CFL及正性重试。高精度热算子仍可能非M矩阵；范数及正性检查不是任意Cut-cell的稳定性证明。`quadraticHeatWalls` / `quadraticViscousWalls` 声明实际采用二次重构的面数；原始多边形独立读取器使用再正交Gram–Schmidt构建参考系数，与原生Householder实现区分，复核通量、壁功、扩展行范数和实际面数。

精度验证分开组织，通常约1–2分钟完成本轮新增数值测试，完整回归另计：

- `wall_gradient_test.cpp`：二次温度/速度场、扭曲及旋转网格、长宽比0.03/3/30、长度缩放1e−3/1/1e3；对热流、应力和面中心壁功进行4096 epsilon的归一化代数检查，直接扰动每个自由度检查完整Jacobian。该阈值容纳QR和差分累计舍入，只验证代数一致性。缺秩最小网格必须拒绝。新增极端长宽比算例会把已舍入的绝对点值误差按1/法向间距放大，因此另以实际梯度系数传播输入舍入包络进行相同4096 epsilon检查，并同时保留原几何门限及未经包络缩放的误差；不把这个代数诊断作为物理精度门限。
- `transport_precision_test.cpp`：`T=2+0.2 sin(πx)sinh(πy)/sinh(π)` 的无源稳态导热，8/16/32平方网格及扭曲网格。独立单位扰动恢复实际算子，线性求解后再针对真实通量残差修正；线性容差1e−11与1e−13的场差应小于温度离散误差1%，排除求解容差污染。壁面热流比较**解析边积分**，L1分子为逐面绝对误差之和，分母为全部壁面解析热流绝对值之和（均W/m）。最细误差目标0.5%并至少优于原格式2倍，细化观测阶>1.7；这是小规模光滑热场目标，不是所有工程算例的通用标准。
- `euler_wall_accuracy_cli_test.py`：Couette含黏性发热，μ=k=0.1、R=1、U=0.5、H=1、Tw=1。8/16/32层解析初值保持试验之外，还从静止、均匀冷场推进8/16/32层算例至t=60，验证温度、速度、场变化率、壁面功热平衡和全过程累计能量。温度误差以解析最大温升0.03125 K归一化，速度以0.5 m/s归一化，目标均为1e−6；最后一步保守量时间导数也须小于测试单位下1e−7，不能只凭目标时间或净热量接近零宣布稳态。冷启动封闭质量固定，解析稳态压力由质量约束决定，不把压力强定为初值1。用解析积分 `p=1/[∫₀¹ 1/T(y) dy]` 另查压力的空间细化，观测阶>1.8且最细相对L∞误差<0.01%，防止温度多项式恰好复现掩盖全场误差。
- 时间精度使用固定网格冷启动，dt=1e−3/5e−4/2.5e−4，相同终点0.02；以dt=1.5625e−5的离散轨迹为参考，再与3.125e−5核对参考差异须小于最细被测误差10%。四个保守分量的RMS误差应表现为二阶（观测阶>1.8）；这是时间自收敛，不是连续PDE解析轨迹。
- 同物理条件按接受步中断/续算必须逐字节一致；显式切换壁面数值格式可以续算。篡改通量、步长系数和壁面格式/面数必须被独立审计拒绝。原线性路径与修改前二进制比较黏性/导热/无黏四组检查点。

时间细化还保留了一个最小失败例：dt=1.5625e−5累加1280次后，浮点时间可能比0.02少约5e−16，原CLI因尾步小于声明最小步长而错误失败。`EulerStepControls2D::endTime` 现在让求解器检查实际CFL步长后的剩余区间，必要时把倒数第二步分成两个合规小步，真实计算每一步通量；不伪造终点时间，不绕过最小步长。全部接受步均检查上下限，检查点仍保存实际接受时钟。

```sh
ctest --test-dir build -R 'cartmesh2d_(wall_gradient|transport_precision|euler_wall_accuracy)' --output-on-failure
python3 tests/euler_wall_accuracy_cli_test.py --cli build/cartmesh2d_euler_cli --output outputs/euler-wall-accuracy/validation
build/cartmesh2d_euler_cli --mesh final.solver.cm2d --output outputs/euler/run --case external --viscosity .02 --wall-model no-slip --conductivity 100 --wall-thermal temperature --wall-value 400 --flux hllc --order 2 --wall-gradient quadratic --density 1.225 --pressure 101325 --u 50 --end-time .00005
```

最后一条仍为短时App验证物性，不是空气推荐值。实际曲壁App和独立面通量审计证明功能路径与离散实现相符，不能代替曲壁摩擦/换热关联式、网格无关性或跨平台精度验收。[壁面精度证据](../artifacts/current/native-euler-wall-accuracy.json) · [实际结果与细化图](../artifacts/current/native-euler-wall-accuracy.png)
