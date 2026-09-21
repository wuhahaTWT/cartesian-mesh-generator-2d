# 开发导航

## 完整笛卡尔背景网格

`cartmesh2d_cli --background-grid adaptive|uniform` 在已有几何诊断、Quadtree细化及2:1平衡后直接导出完整叶子，不执行Cut-cell、合并、Solver修复或OpenFOAM输出。均匀模式令全域最低层级等于最高层级；自适应模式复用距离带、盒加密及尺寸场。当前资源上限分别为level 10和12，不是已测性能保证。

```sh
build/cartmesh2d_cli examples/acceptance/circle.xy outputs/background-grid/circle 7 0.5 0.1 exterior - 3 --background-grid adaptive
python3 tools/verification/check_background_grid.py outputs/background-grid/circle.background.json
```

`.background.json` 为 `cartmesh2d-background-v1`，包含域、原始边界、完整单元包围盒/层级/整数格坐标及分类（0外部、1内部、2相交）；`solver_ready=false`。`.background.vtk` 保存相同完整四边形及分类/层级，可由ParaView读取。自适应粗细交界可有悬挂节点，VTK按独立四边形展示，不宣称具有统一求解面拓扑。JSON不是CM2D求解网格，不能输入现有流体求解器。

实现入口：`src/io/BackgroundGridIO2D.cpp`；独立审核核对整数格覆盖/无重叠、几何分类及2:1平衡；`cartmesh2d_background_grid`还核对JSON/VTK一致、重复确定性、非法几何及误用导出的拒绝。本轮桌面接入仍待完成。

## 分支与里程碑

`mesher-v0.3.0` 是固定标签（`691c97e`）；`codex/mesh-maintenance` 是可继续修改的网格维护线；`codex/cfd-development` 是包含网格核心的 CFD 开发线；`main` 是已验收的集成线。不要移动已有里程碑标签去“更新版本”，后续里程碑另建标签。

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
| `quality/SolverTopology2D`、`PatchLocalQuality2D` | 凸划分、源邻域相对评分、精确合并/切分批量与完整质量验收 |
| `io/Dxf2D`、`BoundaryMetadata2D` | CAD 曲线、单位、边界名称/角色 |
| `io/MeshIO2D`、`OpenFoam2D` | CM2D/VTK/JSON、二维挤出和 OpenFOAM case |

表中的模块分别位于 `include/cartmesh2d/` 和 `src/`，仍参与 CMake 构建，未把现用模块当旧版本删除。

桌面入口是 `desktop/src/main.js`（编排/IPC）、`preload.js`（受限桥接）。
`core/` 管能力、样例、几何/SVG 输入、参数、尺寸预算、CM2D、结果摘要；`process.js` 负责可取消/限时的原生子进程，`automatic.js` 保留旧版最多 4 组候选；`cell-budget.js` 负责数量初值与实测调整，`budget-runner.js` 最多 3 次执行并保留最接近目标的成功结果。

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

构建依赖 Node.js 22、CMake 3.20+、Python 3 和 C++20 编译器；Windows 使用 Visual Studio 2022 C++ 桌面工作负载（含 Windows SDK），Linux 使用 GCC 11+ 或兼容 Clang。运行用户无需安装这些工具。`build-native.js` 准备本机三份 CLI、11 个样例和 runtime manifest；打包前检查实际 PE/ELF/Mach-O 格式与架构，拒绝跨系统或错架构混装。macOS 固定系统 clang 并检查动态库；Windows 使用静态 CRT、UTF-8 编码/路径 manifest。ZIP 导出改用流式 `yazl`，不依赖 ditto/zip/PowerShell 命令。

`desktop-platforms.yml` 分别在 macOS 14 arm64、Ubuntu 22.04 x64、Windows Server 2022 x64 执行完整构建/测试。`tools/verification/desktop_platform_smoke.py` 启动打包 App，使用含中文和空格的临时/输出路径，运行 PNG、JPG、hybrid 三例，解压 ZIP 后独立核对 OpenFOAM 单元、图片和标定。Linux 的 `--no-sandbox` 仅用于隔离 CI 虚拟显示测试；产品启动入口不附加该参数。通过后上传版本运行包与证据；不把独立读取器当成真实 checkMesh。

日常用 `ctest --test-dir build -R <相关测试名> --output-on-failure`；完整交付才统一全量运行。
测试目录包含解析、拓扑、守恒、质量、尺寸场、真实输出和确定性检查，不是可删除的缓存。

真实桌面流程（先打包准备 runtime；输出目录预先创建）：

```sh
mkdir -p outputs/smoke
cd desktop
node_modules/.bin/electron . --smoke=circle --out=../outputs/smoke --shot=../outputs/smoke/circle.png
```

可加 `--verified-preset=true`（真实界面载入已有高密案例参数）、`--density=150000`、`--target-cells=100000`、`--auto-padding=0.5`、`--theme=duet`、`--interaction-check=true`（主题状态与实际参数转手动）、`--allow-unsafe=true`、`--control=manual`、`--density=dense`、`--method=hybrid`、`--mode=light`、`--regions=1`、`--repeat=1`。省略 `--out` 即验证默认临时预览；`--export=/绝对路径/result.zip` 验证结果包。每轮会缩到最小窗口并检查底部可达、预览/导出单元数。smoke 会走真实表单/IPC/CLI 后退出；它不等于所有界面操作都已验收。

## 桌面热输运入口

`desktop/src/core/thermal.js`负责请求/边界分组/原生调用/字段、时间与联合状态核对；`thermal-job.js`保留独立运行目录并在完整读回后发布。renderer温度字段与独立流场分别绑定；导出PNG通过后台保存的最终网格/温度生成，不依赖当前相机或预览是否被释放。打包清单包含 `cartmesh2d_transport_cli`。

真实运行管理回归（90秒/进程上限；输出目录必须新建；使用已生成的圆柱最终网格）：

```sh
node tools/verification/verify_desktop_thermal.cjs --mesh PATH/circle.solver.cm2d --output outputs/thermal-desktop-new
```

该工具实算连续/重启一致性、失败不覆盖完整场、取消后正时间保存和继续。不会把测试夹具当数值验证。GUI/打包验证沿用原smoke，增加 `--thermal=true`，将真实表单设为同步温度算例并检验两次成功、一次故意失败、活进程取消与恢复；本轮成功几何参数为 `--control=manual --wall-relative-size=.0625 --background-relative-size=.5 --reference-length=2 --padding-relative-size=10 --band-cells=3 --small-alpha=.1`。默认紧凑域失败另记在CURRENT_STATE。

## 自研求解器入口

已经实现独立可执行的**二维稳态标量扩散/泊松**基础，以及接入桌面的 **SIMPLE + Rhie–Chow 不可压稳态层流**。原有网格生成核心、质量门及 OpenFOAM 导出继续保留。

```sh
cmake --build build --target cartmesh2d_fv_cli -j 4
# 生成固定几何/域/加密带宽度的三档真实 Cut-cell，并独立读回验证：
python3 tools/verification/verify_native_fv.py --output-root outputs/native-fv/validation --problems constant linear sine diffusion
MPLCONFIGDIR=/tmp/cartmesh-fv-mpl python3 tools/visualization/render_native_fv.py outputs/native-fv/validation/summary.json
# 可换成自己的最终网格；这里内壁为1、外域边界为0，无体积源：
build/cartmesh2d_fv_cli --mesh outputs/native-fv/validation/meshes/h02/circle.solver.cm2d --output outputs/native-fv/heat --problem diffusion --wall-value 1 --outer-value 0 --source 0 --diffusivity 1
```

方程是 `-div(k grad(value)) = source`，k 为正的常数扩散系数；当前**所有边界都是 Dirichlet**，按 EmbeddedBoundary / DomainBoundary 指定两个常量，或在制造解模式下施加解析值。这里的 0/1 是标量边界设定，并非流动速度或真实热工工况。没有 Neumann、Robin、周期、逐 patch 配置或时间推进。CLI 拒绝普通 `.cm2d` 和 `.failed.solver.cm2d`，并重新核对索引、边关联、法向闭合、真实面积及原 Solver 质量门；文件名和 AUDIT 不是合格证明。C++ API 应使用 `makeFvMesh2D` 构造缓存；`validateFvMesh2D` 是缓存安全检查，不替代几何工厂。

每条内部边只计算一份 owner 向外的积分扩散通量，两侧以相反符号累加。梯度由邻居质心和 Dirichlet 边中点的加权最小二乘得到；扩散采用 `S = (S·S)/(S·d) d + correction` 分解，非正交项显式迭代，边界也含修正。稀疏两点主部以面连接存储，用自行实现的 Jacobi 预条件 CG 求解。源项用 `source(centroid)*area` 积分近似；变量、制造解误差按质心值解释。当前不保证任意合格网格上的二阶精度或单调性。

完整修正后的逐格通量失衡 L2 范数满足 `absoluteTolerance + relativeTolerance * ||baseRHS||₂` 才算收敛，默认分别 1e-12 / 1e-10；`baseRHS` 为积分源项加隐式 Dirichlet 贡献。它是本标量问题的离散方程停止条件，不与 OpenFOAM 残差直接比较。PCG 另检查真正矩阵残差；超出线性迭代数或数值范围返回错误，非正交修正达到上限则保留场并以退出码 2 / `converged:false` 报告。默认修正上限 400、松弛 0.7，支持显式 `--max-corrections`。

输出 `.vtk`（含 value / 制造解 exact / error）、`.cells.csv`、`.faces.csv`、`.residuals.csv` 和 `.json`。cells 的 source 是积分源项；faces 的 flux 是积分 `-k grad(value)·S`，不是质量流量。没有解析解的 diffusion 模式在 CSV exact/error 写 NaN、JSON 误差写 null，不伪造精度。验证工具独立读取 CM2D、CSV，复算几何、方程残差及误差，绘图使用真实多边形。

数值设计参考公开的 [MOOSE 有限体积设计说明](https://mooseframework.inl.gov/finite_volumes/fv_design.html)中的共享面守恒和非正交思想。新增 C++ 为本仓库实现，没有引入 MOOSE/OpenFOAM 求解核心；这是公开方法的自研实现，不声称提出新 FVM 算法。`tests/fv_test.cpp` 覆盖斜网格和粗细交界、制造解收敛、错误缓存与溢出；`tests/fv_cli_test.py` 覆盖实际 CLI 输出和明确失败。最新数值与验证范围只在 CURRENT_STATE 维护。

### 守恒标量与恒物性热输运

独立入口 `cartmesh2d_transport_cli` 使用 `ScalarTransport2D`。它扩展了上方旧扩散CLI没有的对流、通量边界和时间推进；不会修改旧扩散/流动入口的行为。

```sh
cmake --build build --target cartmesh2d_transport_cli -j 2
# 新目录，实际生成三档网格，分别验证迎风/限制线性和共同终止时刻的时间细化：
python3 tools/verification/verify_scalar_transport.py --generate outputs/scalar-new
# 从已有非定常流动的 accepted checkpoint 读取并冻结 U.S；逐面BC需匹配同一最终网格：
build/cartmesh2d_transport_cli --mesh /path/case.solver.cm2d --flow-checkpoint /path/flow.checkpoint --boundary /path/thermal.csv --output outputs/temperature/result --diffusivity .1 --convection upwind
# 标量在冻结载流上推进：另加 --dt .01 --steps 20 --initial 0
python3 tools/verification/verify_scalar_transport.py --prefix outputs/temperature/result --output outputs/temperature/audit.json
```

热边界CSV必须完整列出每个边界face ID，禁止内部面、重复或漏项。ID与CM2D及flow.faces.csv对应，不能按屏幕位置猜序号：

```csv
face,type,value,inflowValue
12,value,350,
19,flux,0,
25,flux,0,300
```

此片段仅说明格式，实际文件须含全部边界面。`value`是面定值；`flux`是向外 `-D grad(theta).n` **每单位边长**（不是积分面通量），0即绝热/零扩散通量。flux边界有负载流时，必须给`inflowValue`；上述25号面在回流时温度为300。流出对流采用owner迎风/限制重构，定值仍约束扩散；指定非零通量参与梯度重构。稳态连通域必须有定值或给定流入标量，纯绝热无入口没有唯一常数解；非定常可凭前态建立唯一性。

方程 `d(theta)/dt + div(U theta - D grad(theta)) = source`，D必须有限且为正；默认全域常数，冻结载流可显式提供下面的逐面空间场。温度用Kelvin时，D=k/(rho cp)，`--source`是Q/(rho cp)，热通量CSV中是物理向外q/(rho cp)。没有自动材料库/单位推断。当前是单向输运：密度和比热固定，可预设空间扩散系数；不含随温度更新的材料模型、浮力或共轭传热。同步载流与联合续算仍仅支持常数D。基础ScalarTransport API以调用者提供的新时刻边界/源和共享通量推进一步；`--flow-checkpoint`模式载流冻结，`--evolve-flow`模式使用下述同步接口。

输出JSON、VTK、cells/faces/history CSV。cells含前态、积分源项/时间项；faces分别含载流体积通量、对流与扩散标量通量。逐面/逐格读回见`verify_scalar_transport.py`，它同时检查本构离散与几何，不只复述JSON的converged。默认完整方程L2门为1e-12+1e-9*||baseRHS||，并检查失衡/未松弛对角系数<=1e-9；这是有单位的代数停止设置，不是全软件的精度评级。未收敛返回2；非法输入/线性求解失败返回1。冻结载流模式尚无标量checkpoint；同步模式见下方联合保存。温度桌面入口见下方；失败的多步计算不可假作已完成全部物理时间。

同步模式按以下顺序执行：`advanceIncompressible2D` → 新时刻共享通量 → `solveScalarTransport2D` → 两者通过才接受。`ThermalSetup2D`的源项数组和热边界数组为固定空间分布，用于逐项检查续算兼容性；ScalarTransport另外保留回调形式，数组与回调必须二选一。同步源/边界当前固定时间，不含浮力反馈。

```sh
# 需要完整匹配此网格的热边界CSV；从静止流体和初始温度开始：
build/cartmesh2d_transport_cli --mesh /path/case.solver.cm2d --boundary /path/thermal.csv --output outputs/thermal/run --evolve-flow external --flow-nu .1 --flow-speed 1 --diffusivity .1 --initial 0 --dt .05 --steps 4 --pressure-preconditioner aggregation
# 同时续算流动和温度；物性、源、BC及格式必须相同，迭代控制/时间步可以更改：
build/cartmesh2d_transport_cli --mesh /path/case.solver.cm2d --boundary /path/thermal.csv --output outputs/thermal/continued --evolve-flow external --flow-nu .1 --flow-speed 1 --diffusivity .1 --dt .05 --steps 4 --pressure-preconditioner aggregation --restart outputs/thermal/run.thermal.checkpoint
# 新目录，重生成三档方形和真实圆柱，运行同步解析验证、独立流场对照和续算比较：
python3 tools/verification/verify_thermal_flow.py --generate outputs/thermal-verification-new
```

`--flow-convection`单独控制速度对流，`--convection`控制标量；默认均迎风。`--flow-tolerance`默认1e-8，`--flow-max-iterations`默认1500，标量保留原完整方程停止设置。同步可选channel/cavity/external，边界条件仍依各自流动工况适用；它没有自动识别任意流动BC。`--verification thermal-vortex`是单位方形解析测试专用入口，内置初值和BC，不能混入用户的BC文件或冻结流场。

权威续算文件是单个`.thermal.checkpoint`，先临时写入再原子替换；`.carrier.checkpoint`仅供独立读取/诊断。`.thermal-history.csv`记录各步是否被接受、两套残差、储热量、热量平衡及载流CFL。流动/标量未收敛返回2并保留上次接受的joint checkpoint；异常返回1，取消进程后也可从最后完整保存时刻续算。失败输出JSON明确标记未完成，旧VTK/CSV即使还在也不能冒充新结果。读取时检查几何、时间、全场及物理配置；同编译器续算逐位一致已实测，跨平台相同位模式尚未认证。

温度/浓度采用同类输运方程的官方参考：[OpenFOAM scalarTransport方程说明](https://api.openfoam.com/2606/classFoam_1_1functionObjects_1_1scalarTransport.html)。本仓库自行实现FVM装配，复用自己的稀疏求解与面算子，没有复制或链接OpenFOAM代码，不声称原创输运方程。大输出仍在outputs，不加入日常源码历史。

### 实验性理想气体 Euler

独立入口 `apps/cartmesh2d_euler_cli.cpp`，核心 `fv/Euler2D.hpp` / `src/fv/Euler2D.cpp`，续算 `fv/EulerCheckpoint2D.hpp`。仅原生二维：每个真实面共享 Rusanov 质量、动量、总能量通量，前向欧拉，理想气体 `p=(gamma−1)(rhoE−rho|U|²/2)`。密度/压力必须为正；按真实面积与各面声学波速之和选择时间步，失败候选仅缩步重试，不裁剪状态。

```sh
cmake --build build --target cartmesh2d_euler_cli cartmesh2d_euler_tests --parallel 2
# 矩形内流网格的 Sod 冲击管；示例为量纲一致的 rho=1、p=1、R=1 参考。
./build/cartmesh2d_euler_cli --mesh outputs/tube.solver.cm2d --output outputs/euler-tube --case sod --gas-r 1 --end-time .2
# 实際物理时长须按域长/声速选择；上面的 .2 只对应单位长度参考，不是任意 SI 工况。
python3 tools/verification/verify_euler.py --mesh outputs/tube.solver.cm2d --prefix outputs/euler-tube
```

`--case external` 在 EmbeddedBoundary 施加滑移壁面、DomainBoundary 施加特征远场；`uniform` 全边界为远场；`sod` 只允许轴对齐矩形，左右透射、上下滑移；`vortex` 是宽高至少20的周期矩形解析涡参考，固定环境rho=p=1、u=v=1。`--density/--u/--v/--pressure/--gamma/--gas-r` 定义密度、速度、绝对压力和气体，Euler 压力不是不可压求解器的运动学压力。当前一阶初始化按单元中心采样，不是高阶体平均积分。

自定义边界先用 `--export-boundaries FILE` 导出再编辑物理列，以 `--case custom --boundary FILE` 使用。文件绑定单元/面数量和每个边界的 owner、中心、面积向量；必须完整覆盖实际边界。四种类型为 `slip-wall`、`transmissive`、`farfield`、`periodic`，周期配对须互反、等长反向且同组平移一致。透射并不是通用指定静压出口，滑移壁面没有黏性作用。

`--end-time` 指绝对目标时间；`--max-step/--min-step/--cfl` 控制显式步，声学 CFL 定义为 `dt*sum(length*(|un|+a))/area`，最大允许.45。默认180秒、100000接受步；预算耗尽或SIGTERM取消明确非零退出并保存已接受状态，不伪报到达时间或稳态收敛。默认每25步原子保存检查点，正常退出/可控失败时再次保存；SIGKILL/断电只能保留最近持久化状态。`--restart FILE` 核对完整 Fv 网格、气体、边界与工况，允许改变数值步长和最终时间。

输出CSV含最后接受步的新旧守恒状态、实际面通量和波速，便于独立重算四个守恒方程；JSON/VTK包括密度、绝对压力、温度、Mach数和速度。`verify_euler.py` 从CM2D多边形重建几何、外法线、通量、EOS、CFL及逐格/全域守恒，并检查每条时间记录；它并未独立重演整条计算轨迹。`tests/euler_test.cpp` 覆盖非正交均匀流、周期守恒、旋转冲击管、强爆炸/膨胀正值和续算；CLI测试覆盖精确Sod细化、取消、预算续算、边端点反转及篡改检出。

方法与解析参考：[Clawpack Euler Riemann chapter](https://www.clawpack.org/riemann_book/html/Euler.html)。当前仅一阶无黏 Euler，等熵涡的精度失败保留在 CURRENT_STATE；App、MUSCL/高阶时间推进、可压黏性及导热都未完成，不能把现有被动温度输运当成可压总能量方程。

### 原生层流求解

```sh
cmake --build build --target cartmesh2d_flow_cli -j 4
build/cartmesh2d_flow_cli --mesh /path/case.solver.cm2d --output outputs/flow/result --case external --nu 0.1 --speed 1 --max-iterations 1500
python3 tools/verification/verify_native_flow.py --generate --output-root outputs/native-flow/reproduce --max-iterations 1500
MPLCONFIGDIR=/tmp/cartmesh-flow-mpl python3 tools/visualization/render_native_flow.py --summary outputs/native-flow/reproduce/summary.json --output outputs/native-flow/reproduce/figures
```

`external` 为左侧恒速入口、右侧运动学压力0、上下滑移及物面无滑移；`channel` 为无孔矩形内域、左侧抛物线入口、右侧压力0及上下无滑移，speed 是抛物线峰值；`cavity` 为无孔矩形腔、顶盖水平移动、其他壁面静止，speed 是顶盖速度，固定 cell 0 的压力为0。只接受一个连通流体区域；边界位置/方向或内域形状不符、出口回流、数值范围错误均明确失败。它不是任意喷管/多孔腔的自动边界配置器。

`duct`（桌面0.4.12）用于有左右竖直开口的曲壁内流：最左侧端面均匀速度`speed`，最右侧运动学压力0，其余边界（含孔洞）无滑移。端口按最终边界的真实坐标极值和法向识别；曲壁不能套用旧矩形`channel`。旋转后的非竖直端口、没有开口的圆环不支持此预设，会明确失败；任意方向的条件使用下面的显式命名边界。稳态、非定常、检查点续算和被动热输运共用同一边界定义，App同步提供该选项；几何物面标签不等于每个面的流动边界类型。

```sh
build/cartmesh2d_cli examples/complex/nozzle_profile.xy outputs/duct/mesh 5 0.03333333333333333 0.1 interior outputs/duct/openfoam 5 0
build/cartmesh2d_flow_cli --mesh outputs/duct/mesh.solver.cm2d --output outputs/duct/flow --case duct --nu 0.1 --speed 1 --convection limited-linear --pressure-preconditioner aggregation --max-iterations 4000 --tolerance 1e-9
python3 tools/verification/verify_duct_flow.py --output outputs/duct-study
```

`custom`（桌面0.4.13）读取显式命名边界：速度入口、指定运动学压力出口、静止／切向移动壁面。每个最终边界面恰好一条记录，内部面、重复、缺失、同名混合类型、壁面穿透或错误网格均拒绝；任意面方向可用。同组允许逐面不同值，`speed`只控制归一化，不覆盖实际入口速度。闭域只允许壁面，固定cell 0压力；开放域需要入口与出口，当前只支持回流拒绝。暂不支持自定义滑移、压力驱动入口、时间变化边界或联合温度/SST。外流原有预设继续支持上下滑移。

```sh
# 由已验证的 duct/channel/cavity 转出模板，不启动求解。
build/cartmesh2d_flow_cli --mesh outputs/duct/mesh.solver.cm2d --case duct --speed 1 --export-boundaries outputs/duct/input.boundaries
# 编辑命名边界后求解；入口值为各面的物理速度，p为运动学压力。
build/cartmesh2d_flow_cli --mesh outputs/duct/mesh.solver.cm2d --case custom --boundary outputs/duct/input.boundaries --output outputs/duct/custom --nu .1 --speed 1 --pressure-preconditioner aggregation
```

文件格式为`CARTMESH2D_FLOW_BOUNDARIES 1`、`COUNTS cells faces boundaryFaces`、逐行`BOUNDARY face owner centreX centreY Sx Sy type "name" u v p`、`END`。坐标、外向面积矢量和owner与最终网格绑定，不允许只靠相同面数量认作同一网格。可保留类型字符串`velocity-inlet`、`pressure-outlet`、`wall`、`moving-wall`；非适用的速度／压力项必须为0。输出`PREFIX.boundaries`是实际输入快照，摘要中的条件须与之相同；独立稳态/非定常审核从CM2D重新求外法向、重建逐面动量和守恒，不以摘要自报收敛替代。

App选择“命名边界”，从预设生成或导入`.boundaries`，按组编辑名称、类型、速度分量和出口压力，点击“显示位置”高亮对应面。原本逐面变化的值在编辑器留空保持，填写则统一该组；从模板转换抛物线入口不会自动变为均匀入口。修改网格需重新配置。非定常检查点v4记录完整边界身份、数值及网格，续算锁定这些条件；原预设继续使用v2。ZIP携带边界输入与最后接受状态。运动学压力乘密度得到Pa；二维力仍按单位厚度和密度归一化。

喷管study串行生成三档真实网格，独立复核几何、共享面通量/动量、壁面力，记录压降、最大速度和壁面力的网格变化；每个子进程默认180秒，失败保留日志。`--verify-only`只允许与已记录成功命令完全相同的参数，重审已有原始场。守恒通过不等于参考精度认证；无参考的喷管始终保留`not-qualified`状态。App的原生导出可独立重新读取，OpenFOAM边界条件仍需单独设置。

单元中心速度/运动学压力，共享边积分体积通量。动量默认一阶迎风，可用 `--convection limited-linear` 选择限制线性重构，黏性项用最小二乘梯度及显式非正交修正；内部面速度包含偏斜修正和 Rhie–Chow 压力项，压力修正使用四次非正交迭代，最终通量与最后一次实际线性方程一致。动量松弛0.6、压力松弛0.25；动量使用自行实现的 Jacobi–BiCGStab，压力修正利用对称正定结构使用 IC(0)–PCG，保留 `--pressure-preconditioner jacobi` 对照，均检查真正矩阵残差。设计依据包括 [MOOSE 的同位有限体积说明](https://mooseframework.inl.gov/modules/navier_stokes/insfv.html)中关于 Rhie–Chow 和压力零空间的说明；没有复制或链接其求解核心。

压力采用每个共享面唯一的运动学压力值：内部面按几何权重插值并修正面中心偏斜；出口面取0，其他边界按内部压力的最小二乘梯度外推，重构时不对未知壁面压力强加零法向梯度。`sum(p_face*S)/area` 进入动量源；Rhie–Chow 中取消单元压力响应的项、压力修正对单元速度的作用使用同一 Gauss 算子。面法向压力差仍保留直接相邻压力差及最小二乘非正交修正，不能用平均面压力替代这部分，否则棋盘压力可能成为零模态。壁面压力积分复用同一个面值。

压力梯度先使用直接邻居和已知压力边界；忽略未知压力边界的单元至少加入完整第二圈真实样本，秩不足的尖角也尝试扩展。随后检查归一化方向最小二乘正规矩阵的谱条件数：大于16时继续加入完整的下一圈，最多6圈或直到可达邻域用尽。16是重构邻域的几何选择参数（设计矩阵条件数目标4），不是网格或流动验收阈值；未达到目标时仍使用可达样本并执行原数值秩检查，不宣称任意网格均稳定。样本按cell ID排序去重，已有直接共享面采样权重不变；常数/仿射压力保持一致，几何选择随网格一起旋转。原反例条件数565→8.20只需两圈；2,504格圆柱的cell1723第二圈仍为140.39，扩展第三圈后为11.19。核心修复未修改几何、物理边界、压力修正次数、松弛或停止条件。

未知壁面压力仍线性外推，不强加物理零法向压力梯度；压力修正仍不能改变指定速度边界的体积通量。速度滑移/出流约束保持原定义。相关边界概念见 [MOOSE压力外推说明](https://mooseframework.inl.gov/source/linearfvbcs/LinearFVExtrapolatedPressureBC.html)，具体重构为本仓库原生实现。新摘要使用`pressureBoundaryReconstruction=one-sided-linear-adaptive`；独立审核器保留`one-sided-linear-2ring`的固定两圈、`one-sided-linear`的仅秩不足扩展，以及缺字段时旧零法向算法，未知标记拒绝。续算的逐字节一致性只对相同求解器版本、设置和工具链成立；算子修复前后的结果按各自摘要标记审核。

真实局部反例保留在 `tests/flow_face_test.cpp::fullRankBoundaryPressureStencil`及`pressureStencilBeyondSecondRing`；完整粗圆柱生成、稳态/非定常与续算回归在`tests/pressure_stencil_cli_test.py`。可视化命令：`python3 tools/visualization/render_pressure_stencil.py --mesh <mesh.solver.cm2d> --prefix <accepted-flow-prefix> --cell 2251 --time .1 --output <figure.png>`；2251仅是本次固定圆柱的单元ID，其他网格需重新选取。

`--outlet-backflow reject|normal-inlet` 默认reject；normal-inlet仅在右侧压力出口实际面通量q<0时激活：p=0不变，法向速度零法向梯度，切向速度固定0；流出仍使用旧零梯度/重构。动量中负q的法向q*Uowner放入显式右端，保持原完整方程与正主对角；最终残差和输出通量使用当前场，不能截掉负通量。每次压力修正后用新通量更新边界掩码，再重算原动量残差。摘要记录 `outletBackflow`、`outletBackflowFaces`、`outletInflow`（正的m²/s流入量）；没有压力出口的封闭算例不受此选项影响。模型形式参考[OpenFOAM pressureInletOutletVelocity](https://api.openfoam.com/2606/classFoam_1_1pressureInletOutletVelocityFvPatchVectorField.html)，实现为本仓库独立代码；不是任意方向/湍流或能量稳定开放边界资格。

新增稳态验证专用 `--case counterflow`：单位方形u=U(1+2cos(2πy))、v=p=0，体积源sx=8π²nuUcos(2πy)；左侧给定有符号速度，右侧压力出口，上下滑移。正负出口通量同时存在，用于测试模型和空间误差；不作为桌面实际物理工况。原external求解不增加该体积源。`verify_native_flow.py --cases counterflow --outlet-backflow normal-inlet` 支持生成及独立审核；`render_outlet_backflow.py --help` 给出真实场/误差图入口。新增原生热输运回归从有符号平行初流启动到通道边界，检查真实回流温度、缺值拒绝和联合续算。

流动checkpoint写v2，CONFIG末尾记录模式，读取v1默认reject；更改模式后续算拒绝。联合热状态外层仍是v1，仅嵌入的FLOW升级v2。桌面、冻结载流CLI及独立读取同时支持旧状态；模式不能只改JSON就绕过原生配置核对。

`--convection upwind|limited-linear|face-limited-linear` 的默认值为 upwind。限制线性格式保留隐式迎风矩阵，将唯一上游面重构值与迎风单元值之差作为显式共享面修正，owner/neighbour 严格反号；原方程残差也包含该修正。每个速度分量的 cell limiter 使全部实际面中心重构值落在邻居/Dirichlet 边界的局部范围内，黏性梯度不受此 limiter 修改。固定速度边界直接用边界值；自由分量的出流使用单边受限重构；默认拒绝回流，显式normal-inlet按上文处理压力出口的流入。这是 **Barth–Jespersen 风格的面值限制**，不是速度场全局最大值原理，也不是任意网格上的完整 Navier–Stokes 二阶证明。光滑指数函数面值细化单独检查重构阶数，实际流动另做基准比较。方法依据参考 [Barth–Jespersen（1989）](https://ntrs.nasa.gov/citations/19890037939)的多维单调线性重构、[OpenFOAM 梯度限制说明](https://doc.openfoam.com/2306/tools/processing/numerics/schemes/gradient/)及 [MOOSE 压力动量项](https://mooseframework.inl.gov/source/fvkernels/INSFVMomentumPressure.html)，实现位于本仓库 `FlowFaceOperators2D.hpp`，没有复制其代码。

新增 `face-limited-linear` 沿每个目标面的单位法向/切向投影速度邻居、Dirichlet样本和梯度，在该面坐标系内做相同的Barth–Jespersen限制，再还原唯一共享面向量。局部投影面值受界限约束；不宣称耦合速度场满足全局最大值原理。在轴对齐矩形网格上与原逐分量限制一致，网格及完整向量边界一起旋转后保持等价。原`limited-linear`逐全局x/y分量限制具有坐标方向依赖；保留名称与语义，旧检查点不会静默换算法。该现象和面坐标思路可参见 [Deltares技术手册Remark 6.2.16](https://content.oss.deltares.nl/dhydro/D-Flow_FM_Technical_Reference_Manual.pdf)，本实现使用本仓库原有局部重构及原生有限体积动量方程。移动/静止固定壁值、压力出口与回流策略保持独立处理。App、流动checkpoint、导出及被动热输运的载流格式均记录新名称；标量温度仍只用upwind或limited-linear。

`python3 tools/verification/verify_named_channel.py --output outputs/native-flow/channel-new --timeout 240` 在真实二维内域Cut-cell网格上比较三档解析Poiseuille通道及0.63rad旋转，逐例独立重建共享面动量/质量，并分离迭代误差。默认6000轮上限、正常容差1e-8和细档加严1e-10；不改原生停止条件。参考 [MIT充分发展层流解](https://ocw.mit.edu/courses/2-25-advanced-fluid-mechanics-fall-2013/1a114d602956fa0dd328155f9b45f93d_MIT2_25F13_Couet_and_Pois.pdf)：L=4、H=1、Umax=1、nu=.1、出口p=.25，u=4y(1-y)、dp/dx=-.8、Q=2/3、上下壁黏性合力/密度/厚度=3.2。指标限值随报告保存，只对该工况作解析精度验证；原始失败目录不覆盖。

本 CLI 的停止条件是：至少10次迭代，动量残差、相对速度变化、相对压力变化均小于 `--tolerance`（默认1e-6），逐格连续性及全局相对流量失衡均小于1e-8。动量残差是原离散方程失衡除以 `(aP_u+aP_v)*Uref`；逐格连续性为 `|sum(flux)|/(Uref*sqrt(area))`；速度变化以 Uref 归一化，压力变化以 `Uref²+nu*Uref/domainHeight` 归一化。全局失衡除以总入流，封闭腔使用 `Uref*domainHeight`。这些是本实现的数值停止条件，不是所有 CFD 软件的统一精度标准，不能与 OpenFOAM residual 数字直接等同。`converged` 也不等于网格无关或物理模型适用。

退出0代表满足上述停止条件；退出2代表到达上限，保存诊断场但 `converged:false`；退出1代表输入/数值失败。输出六种文件：`.json` 工况/状态/单位/压力基准、`.fields.json` 桌面字段、`.cells.csv` 单元 u/v/p、`.faces.csv` owner向外的积分体积通量，以及面压力、对流/扩散动量通量、`.residuals.csv` 全迭代历史、`.vtk` 原多边形上的速度/压力。p 为 p/ρ，单位 m²/s²；力为流体对静止 EmbeddedBoundary 的积分力除以密度和深度，单位 m³/s²，不是 Cd/Cl。`domainHeight` 只指外域高度，不是物体参考直径。桌面验证全字段有限、单元ID/数量/工况与本次最终网格一致后才绑定；失败或取消保留 `flow-incomplete-*` 诊断，部分复制文件不会充当完整结果。

面动量默认使用 `q*U_face + p_face*S - nu*(grad(U)+grad(U)^T)·S`，CSV 的 advectionX/Y、diffusionX/Y 分别保存对流项和完整黏性项；二维每单位深度的单位为 m³/s²，pressure 列仍为 m²/s²。`--viscous-stress symmetric` 是默认值，`laplacian` 保留旧单分量扩散用于CLI对照。两者的真实离散残差都包含实际使用的黏性项，不能只给受力报告补一个转置项。独立验证器从 CM2D 多边形中心和边线法向、cell u/v/p、边界定义重建梯度/面值/通量，复算局部和全局动量失衡；旧文件未记录viscousStress时按旧Laplacian审核，新应力标记未知或面通量列不完整时失败。

内部面的速度梯度先按几何权重插值，再沿面法向修正，使其沿两单元中心连线的差等于真实速度差。由此得到的法向导数与原紧致两点/非正交通量一致；转置梯度作为显式共享面修正进入两侧反号的动量方程。固定常量壁面/顶盖的切向速度导数为零，法向导数由真实壁值和法向距离恢复；当前轴对齐滑移边界的切向剪切为零，不能推广为任意方向滑移边界实现。常黏度不可压物理采用 `nu*(gradU+gradU^T)`，不是新的本构模型；公开应力积分惯例可参考 [OpenFOAM Forces](https://doc.openfoam.com/2306/tools/post-processing/function-objects/forces/forces/)。本仓库独立实现，未链接或复制其求解代码。

默认symmetric模式的 `forceX/Y = discreteForceX/Y`，两者均复用动量方程的共享面压力和完整黏性通量；范围仍为静止EmbeddedBoundary。`pressureForceX/Y` 单列压力贡献，`reconstructedForceX/Y` 保留旧单元梯度牵引作诊断，不作为默认力。`wallForceX/Y` 与 `wallViscousForceX/Y` 另统计所有无滑移壁和移动顶盖，faces.csv 的wall列标识这些面。符号是流体对边界的力；滑移、入口、出口不纳入wall总力。旧laplacian模式仍明确保留旧受力定义，不混用。

独立 `audit_wall_tractions.py --summary <runner-summary.json> --summary <physical-summary.json> --output <diagnostic.json>` 对已有通道/制造解结果按真实边线进行16点Gauss积分，比较解析压力、完整黏性牵引和总力；记录每单位壁长的加权RMS、各面与积分总量，不新增无出处的验收阈值。输入哈希须与独立读回报告一致，源报告的失败系列仍保留。受力守恒、解析误差、空间收敛和工程阻力资格分开报告。

解析通道验证速度分布、压降梯度和流量；方腔 Re=100 对比 [Ghia 等（1982）](https://doi.org/10.1016/0021-9991(82)90058-4)中心线数据。圆柱只验证低 Re 定常试算、有限场和守恒，不与几何/边界不同的 DFG 基准混比。误差及外部工具实测范围见 CURRENT_STATE，绘图直接读取 CM2D/CSV。桌面 smoke 可加 `--flow=external --flow-nu=0.1 --flow-speed=1 --flow-max-iterations=30`，迭代上限场不得作为收敛证明。

### 固定网格时间步审核

`tools/verification/verify_thermal_time.py` 对至少三组已完成的同步热涡/独立流动输出进行审核。每组先复用空间验证器的实际几何、逐面本构、守恒、时间项与checkpoint一致性检查，再要求实际mesh SHA、最终物理时间、物性、格式、松弛和停止条件相同，仅dt减小/接受步数增加。比较**连续时间解析解**的面积加权RMS；空间连续后向欧拉参考仍保留，用于观察空间/迭代误差。若任一场误差不下降，则报告失败；观测斜率只做诊断，不硬套一个全软件通用阶数门。

```bash
python3 tools/verification/verify_thermal_time.py --case PATH/dt1/thermal PATH/dt1/flow --case PATH/dt2/thermal PATH/dt2/flow --case PATH/dt3/thermal PATH/dt3/flow --output outputs/thermal-time-audit/study.json
python3 tools/visualization/render_thermal_time.py --study outputs/thermal-time-audit/study.json --output outputs/thermal-time-audit/figure.png
```

真实10,000格、dt=.005/.0025/.00125、t=.01证据在 `artifacts/current/native-flow-thermal-time.json`；它包含原始生成脚本、工具/网格哈希及复用首档的来源，不将复用结果算成本轮重新求解。

### 预设空间扩散系数

`ScalarTransportProblem2D::faceDiffusivity` 为空时沿用原常数D；否则须包含所有内部面与边界面的一份有限正值。`diffusivity`仍须为有限正参考值。每个共享面只有一个D，同时用于矩阵、非正交修正、通量边界梯度与最终面通量；相邻单元用相反符号，不重复制造两份界面热流。API不自动平均单元物性，调用者负责界面插值。例如对齐正交的分层材料，距离加权调和系数可表达两侧串联热阻；这不是任意扭曲材料界面的精度承诺。

冻结载流的稳态/后向欧拉CLI可添加 `--face-diffusivity /path/D.csv`，严格表头 `face,diffusivity`，按当前最终网格face ID完整列出每个面一次。重复、漏项、额外列、非有限或非正值失败。该选项与验证模式、同步载流和联合重启互斥，避免checkpoint丢失物性场。变系数输出的faces CSV增加`diffusivity`列，JSON记录`diffusivityModel=face-values`及输入文件路径；恒系数输出保持原格式。独立读取器从输入CSV重新取得系数并记录哈希，而不是只相信导出列。

```sh
# 真正生成三档最终网格；光滑变D制造解，同一方程/参数/停止门：
python3 tools/verification/verify_scalar_transport.py --generate outputs/variable-D-square --variable-diffusivity
# 斜边平行四边形，仍由本机Cartesian/Cut-cell生成器构造网格：
python3 tools/verification/verify_scalar_transport.py --generate outputs/variable-D-shear --variable-diffusivity --shear .2
# 单个已有网格；不是普通工况的默认物理源项：
build/cartmesh2d_transport_cli --mesh /path/case.solver.cm2d --output outputs/variable-D/result --verification variable-sine --diffusivity .08 --speed .35 --convection limited-linear
```

制造解为`phi=sin(pi*x)sin(pi*y)`、`D(x)=Dref*(1+x)`、`U=(speed,0)`，源项包含`-grad(D).grad(phi)`，固定边界取解析面值。仅允许稳态，域内D须处处为正。读取器从真实几何独立重建D、源、非正交通量与逐格平衡；观察到的约二阶趋势只适用于本组光滑场。双材料接口与倾斜网格非零通量边界另有核心解析回归。桌面和同步热输运尚未开放空间D输入。

### 空间黏度、共享应力与重启

`FlowControls2D::faceViscosity` 为空沿用`nu`；否则每个内部面和边界面提供一个有限正运动学黏度（m²/s）。该值统一进入动量矩阵、非正交修正、完整对称应力、共享面输出、壁面力和旧梯度受力诊断。密度仍固定；`nu`保留为有限正参考值及原压力变化归一化尺度。不自动平均单元材料，不随温度或剪切更新，不是SST。默认完整对称应力保持不变；变黏度时Laplacian与对称应力是不同方程，Laplacian只保留作显式诊断选择，不能混用同一制造源项。

物理`external/channel/cavity` CLI使用`--face-viscosity /path/nu.csv`，严格表头`face,viscosity`，须完整列出同一最终网格的所有face ID一次。steady和backward Euler均支持。重复/漏项/额外列/非正或非有限值失败；验证案例不能混入材料文件。faces CSV新增`viscosity`列，summary记录`viscosityModel=face-values`与输入路径。独立读取器从原始CSV重建系数并保存输入SHA；变黏度工况不套用恒黏度Poiseuille/Ghia参考，明确标记其不适用，几何/逐面/逐格/受力门继续执行。

```sh
build/cartmesh2d_flow_cli --mesh /path/case.solver.cm2d --output outputs/variable-nu/run --case cavity --nu .1 --face-viscosity /path/nu.csv --speed 1 --time-step .01 --steps 2 --tolerance 1e-9
# 续算必须再次提供与checkpoint逐项一致的黏度场：
build/cartmesh2d_flow_cli --mesh /path/case.solver.cm2d --output outputs/variable-nu/continued --case cavity --nu .1 --face-viscosity /path/nu.csv --speed 1 --time-step .01 --steps 2 --tolerance 1e-9 --restart outputs/variable-nu/run.checkpoint
# 独立审核包括原方程的时间项：
python3 tools/verification/verify_transient_flow.py --help
# 稳态制造解，nu(x)=.1*(1+x)，梯度源项按选定应力形式推导：
build/cartmesh2d_flow_cli --mesh /path/unit-square.solver.cm2d --output outputs/variable-nu/mms --case manufactured --nu .1 --speed 1 --manufactured-viscosity-slope 1 --manufactured-pressure-slope .7 --convection limited-linear --tolerance 1e-9 --max-iterations 5000 --pressure-preconditioner aggregation
```

`manufactured-viscosity-slope`仅steady manufactured且大于-1；非零时构造并检查一致的逐面场。独立制造解重建黏度梯度源项、源列和面应力，记录速度/压力/壁面黏性牵引误差。`verify_native_flow.py`可用同名选项生成或核对指定斜率；单例审核未传时使用记录的物性定义。三档对照要固定nu、斜率、速度、压力定义、格式、容差与迭代预算。

非空场写checkpoint v3，CONFIG后增加`FACE_VISCOSITY count values...`；加载时逐项精确匹配，不允许用缺失场的v1/v2恢复变系数计算。恒黏度继续写原v2字节，原v1/v2仍可读。冻结标量入口可读v3并从头核对完整几何/状态。变黏度制造解是稳态验证，checkpoint明确拒绝。当前桌面不提供空间黏度输入，也没有验证桌面导入v3；同步热输运CLI仍使用恒黏度/恒扩散系数。

### SST-2003m输运基础及待完成部分

实现入口为`Sst2003m2D.hpp/.cpp`，固定采用[TMR的SST-2003m](https://tmbwg.github.io/turbmodels/sst.html#sst-2003)。使用应变不变量S而非旧SST的涡量；k与omega均限制生产项为`min(nu_t S²,10 betaStar omega k)`；gamma1=5/9、gamma2=.44，CD下限1e-10。这里固定密度归一化rho=1、SI量纲；CD下限未做任意单位/密度重标定。k=0时通过代数连续极限求`P/nu_t`，不引入虚构k下限。tanh饱和前限制函数自变量防止幂溢出，不修改物理场。

`evaluateSst2003m2D`接收k、omega、nu、最近壁面距离、S以及k/omega梯度，返回F1/F2、混合常数、nu_t、扩散率、生产/损失与交叉扩散。k>=0、omega/nu/d>0，所有输入有限。负交叉扩散C拆成omega损失率`-C/omega_old`；该步骤在旧场处严格恢复原源，外层须重新评价非线性方程。

`solveFrozenSst2003mTransport2D`只执行一轮冻结系数k/omega输运。调用者提供一致的梯度、距离、应变、载流面通量及全部边界；内部面按已有几何权重线性插值单元扩散率，边界取owner系数。它不提供壁函数、不自行判断湍流入口、不更新动量/压力，也不保证非正交高阶格式无条件保持正性；不收敛、负k或非正omega均抛出失败。返回的内层converged不能充当RANS外层收敛。

隐式损失使用`ScalarTransportProblem2D::sinkRate`：空表示零，否则每格一个有限非负1/s系数；结果的`sinkIntegrals/sinkIntegral`纳入局部/全域平衡。瞬态可与后向欧拉叠加，稳态正损失可锚定全Neumann分量。线性求解现在同时满足原norm门与标量外层norm目标的一半，避免紧外层容差在较松内层门处停滞；默认标量门下原内层门仍更紧。

```sh
# 测试驱动，不是用户RANS求解入口。复用已有真实最终网格：
python3 tools/verification/verify_sst_decay.py --probe build/cartmesh2d_sst_transport_tests --mesh /path/final.solver.cm2d --output outputs/sst-decay
```

驱动对给定F1=1、无梯度/剪切的均匀衰减反复执行冻结步骤，检查原非线性残差；分别比较后向欧拉解析根与连续ODE参考，进行三档时间步细化。Python读取实际CM2D/全部最终单元/面及历史代表值进行独立审核，并有篡改拒绝测试。这不构成壁面流动验证。空间非均匀冻结输运与真实壁段距离的后续进展见下一节；近壁条件与固定载流非线性进展见下节；仍需速度压力耦合/压力定义、标准平板及分离工况证据；当前无产品SST开关，也没有自动湍流热扩散或温度相关物性。

### SST壁面距离与空间重建

`WallDistance2D.hpp/.cpp`提供`computeWallDistance2D(mesh, wallFaces)`：wall mask按当前最终网格face ID，必须至少选一个真实boundary face；不能把内部面或所有外域边界默认当壁面。返回每格`distance`和`nearestFace`以及诊断`segmentTests`。索引只在一次调用内创建；固定网格调用者可保存距离场，不应每次非线性迭代重算。有限线段端点由直线面缓存的centre和areaVector恢复，最近点可能是端点；AABB树按几何中心/face ID确定排序，同一计算中精确等距取较小ID。没有跨编译器位级确定性的额外资格声明，也不支持曲面高阶几何。

`reconstructSst2003mGradients2D(mesh, frozenProblem, velocity, velocityBoundary)`按各分量固定值/零法向导数作最小二乘。k/omega只接受Value或零DiffusiveFlux；若需要非零通量，必须后续明确有效扩散系数及迭代约定，当前显式拒绝。返回k、omega、u、v梯度和`S=sqrt(2 ux²+2 vy²+(uy+vx)²)`，将k/omega梯度及S写入冻结问题即可使用。辅助函数不替调用方选壁面、不推断壁函数、不代表非线性迭代完成。

```sh
# 指定仿射场的空间验证驱动（非产品RANS入口）：
build/cartmesh2d_sst_spatial_tests /path/unit-domain.solver.cm2d outputs/spatial-probe
python3 tools/verification/verify_sst_spatial.py --mesh /path/unit-domain.solver.cm2d --prefix outputs/spatial-probe --output outputs/spatial-audit.json
# 仅距离的规模探针，选择全部外边界；独立checker限定轴对齐矩形：
build/cartmesh2d_wall_distance_tests /path/square.solver.cm2d outputs/distance-probe
python3 tools/verification/verify_wall_distance.py --mesh /path/square.solver.cm2d --prefix outputs/distance-probe --output outputs/distance-audit.json
```

空间驱动固定k=.02(1+.2x+.3y)、omega=4(1+.1x+.2y)、U=(y,0)、nu=1e-5、dt=.01，底部y=0边选作距离诊断；边界值固定为解析场。一次冻结输运的独立审核包含模型系数、重建后的共享面通量、原冻结矩阵norm及逐格收支。它不是完整非线性步，更不是物理近壁流动。读取器对距离/面ID/梯度/模型/源项/面通量均有篡改回归；矩形距离读取器拒绝斜边域，不能把简单解析公式套到任意几何。

### 非线性SST与解析壁面边界

`setSst2003mResolvedWalls2D(mesh, problem, walls)`验证全部wall mask之后才修改问题：墙必须为边界且volumeFlux严格为零，取k_wall=0、omega_wall=60 nu/(.075 d_normal²)，d_normal为owner-centre到当前face的法向间距。`resolvedWalls`使该面的两方程扩散系数取分子nu，内部面仍线性插值单元有效系数。壁距场与壁面间距用途不同，不互相覆盖；此处采用解析到壁面的低雷诺数处理，不是y+壁函数。调用者负责适宜网格和相符速度边界。

`solveSst2003mTransport2D`默认最多500次、外松弛.5；传入当前初值、固定速度及边界、守恒面通量，非定常另给previousK/previousOmega和dt。每次更新重新计算空间梯度和SST闭合；用`evaluateScalarTransport2D`在返回的k/omega重建通量和真实残差，两个方程均满足所设norm及逐格门才接受。初值已满足可以零更新返回，超限返回converged=false；内层失败抛出而非修改上一物理时刻。内层容差取外层的.1倍，无场裁剪；输入、正性、载流守恒和稳态锚定检查不变。评估API不会求解或改变给定场，history仅一行iteration=linearIterations=0。当前重复建立标量工作区，尚未针对大规模非线性SST优化。

```sh
# 诊断入口：单位正方形或平行四边形底边y=0；不是通用RANS CLI。
build/cartmesh2d_sst_nonlinear_probe /path/unit-domain.solver.cm2d outputs/nonlinear-probe
python3 tools/verification/verify_sst_nonlinear.py --mesh /path/unit-domain.solver.cm2d --prefix outputs/nonlinear-probe --output outputs/nonlinear-audit.json
ctest --test-dir build -R sst_nonlinear --output-on-failure
```

探针设置U=(y,0)、初始k=.001/omega=2、nu=1e-5、dt=.02，底壁不穿透，非壁流出为零法向梯度，其余边界固定初值。独立读取器用原始端点重建壁距、几何与梯度，再计算当前源系数和未拆分的非线性PDE收支；历史必须描述实际返回场。近零梯度分量允许`128 epsilon max(1, |grad_native|, |grad_independent|)`的浮点舍入预算，原有相对/绝对误差门仍保留，实际误差/预算均写入审核。它处理独立几何计算的ulp差异，不修改求解收敛条件。既有ODE解析根、空间冻结输运、篡改拒绝与旧标量逐字节兼容性继续保留。

速度/压力耦合的后续实现见下一节；仍需适宜y+的网格、标准平板/分离工况和空间敏感性。当前固定载流案例本身不能证明这些能力；没有将SST加入桌面菜单。

### 实验性稳态SST-RANS耦合

核心入口`SstRans2D.hpp/.cpp::solveSstRans2D(mesh, controls, initialK, initialOmega, progress)`，默认入口k=.001、omega=2，可明确给定初始场；两者需同时提供。沿用已有flow场景角色与几何限制，强制Symmetric应力、拒绝另传faceViscosity和未支持的回流模型。当前实测矩形通道；外流/方腔分类可调用但尚无本阶段耦合资格，不宣称URANS、任意patch或湍流传热。

内部`detail/FlowMaterial2D.hpp`提供本构更新钩子：SIMPLE每次更新速度/压力/守恒面通量后，传出实际边界约束；返回`MaterialState2D{faceViscosity, converged}`，验证系数后重建当前动量方程再验收；converged表示当前本构方程通过，不能用黏度变化小代替。初次动量预测采用分子nu。`turbulenceUpdatesPerIteration`默认1，限制每次SIMPLE中的SST非线性更新数；与`turbulence.maxIterations`取较小值，二者均须正。达到局部更新上限不抛出，继续推进流动；直到同一末态的动量和湍流均通过原门才可整体接受。内层求解/正性失败仍抛出。显式设为`turbulence.maxIterations`可恢复充分求解的nested策略。随后nu+nu_t统一进入共享面应力、矩阵、残差和壁面力。内部面按neighbourWeight线性插值；解析壁面取nu；入口按入口k/omega与owner的距离、梯度和应变评估nu_t；出口按owner闭合。SST标量扩散仍采用此前明确的墙nu、内部插值及其他边界owner系数规则。两套规则均由独立工具重算，未隐藏边界插值选择。

采用TMR的**SST-2003m**，m明确省略各向同性2k/3应力，p直接是运动学压力，不定义p+2k/3修正量。结合不可压约束，动量使用nu_eff(grad U+grad U^T)；没有声称可压缩应力或其他SST变体。只有当前动量/质量与两条当前非线性湍流方程均通过才返回converged。次数字段分别报告SIMPLE和湍流更新；内层失败抛出，不能接受旧系数残差。尚未提供SST联合checkpoint及产品CLI/GUI，当前探针仅供复核。

```sh
# 诊断约定：单位正方形通道，抛物线入口，nu=.001，k_in=.001，omega_in=2。
build/cartmesh2d_sst_rans_probe /path/unit-square.solver.cm2d outputs/rans-probe
# 配对性能诊断可另指定旧的充分求解策略；物理参数/最终停止条件不变。
build/cartmesh2d_sst_rans_probe /path/unit-square.solver.cm2d outputs/rans-nested nested
python3 tools/verification/verify_sst_rans.py --mesh /path/unit-square.solver.cm2d --prefix outputs/rans-probe --output outputs/rans-audit.json
ctest --test-dir build -R sst_rans --output-on-failure
```

读取器独立计算wall omega=60nu/(.075dn²)、原未拆分SST源项及稳态对角/norm目标；不把单纯源项拆分的内层方程当最终收敛。独立面黏度暂存为临时CSV交给既有动量读取器，明确使用one-sided-linear-2ring压力重建；临时目录退出即删除，审核报告保留其系数文件哈希。梯度和相消通量比较有以实际项尺度计算的机器舍入预算，求解器的停止门不变。摘要与最终history中同一原始残差采用精确一致性检查，不为重复导出数据使用独立重算的舍入容差。审核更新次数不超过声明上限；中间湍流状态可未收敛，最终两条原方程门保留。不同规模案例只验证满足离散方程，不能当网格无关性或物理准确度证明。

### 平板混合边界与压力远场

`FlowControls2D::scenario="flatplate"`是实验性**稳态矩形域**配置；`flatPlateLeadingEdge`给出底边无滑移段的起点x，须有限、在[xmin,xmax)内且落在网格面端点。底边前缘之前采用对称，之后无滑移；前缘跨面时报错，不自动拆分或重分类。左边固定(U∞,0)，右边压力0且回流策略沿用flow约定；SST仍拒绝未配置的右侧回流。真实壁距只使用无滑移段。

`flatPlateTop`默认`FlatPlateTop2D::PressureFarfield`：水平顶面压力0，法向速度由压力修正及守恒面通量决定；进入时切向速度固定U∞，法向对流q*v_owner放到显式RHS，不产生负隐式对角；离开时速度零梯度。当前通量方向在每次残差重建前刷新；SST使用相同快照决定进入k/omega固定值及离开零扩散通量，黏度重建与原方程验收仍在同一状态。`Symmetry`则固定法向速度0、切向自由。两种上边界是不同问题，不能混用一个基准结果。

当前最终实算只覆盖开放顶面的**出流**；进入分支有实现及代码审查，但尚无端到端真实输入验证。本配置不等价于TMR的可压缩特征远场。还没有瞬态平板或checkpoint格式：初始化、advance、checkpoint读写明确拒绝相关配置；其他工况传入非默认平板参数也拒绝。旧层流路径和默认输出保持不变，App没有提前开放本工况。

```sh
# 两种诊断均要求完整单位正方形，前缘固定.5且须与底面端点对齐。
# U=1, nu=.001, k_in=.001, omega_in=2，仅边界/离散方程诊断。
build/cartmesh2d_sst_rans_probe /path/unit-square.solver.cm2d outputs/plate flatplate
build/cartmesh2d_sst_rans_probe /path/unit-square.solver.cm2d outputs/plate-symmetry flatplate-symmetry
python3 tools/verification/verify_sst_rans.py --mesh /path/unit-square.solver.cm2d --prefix outputs/plate --output outputs/plate-audit.json
```

独立读取器由CM2D原端点辨认分段，从末态重算壁距、当前模型系数、质量/动量/两湍流方程。`boundarySummary`报告各类面的长度、进入/离开通量；`plateWallSamples`取水平底壁切向离散牵引tau/rho，给Cf=2tau/(rho U∞²)和y+=d_normal*sqrt(abs(tau)/rho)/nu。默认诊断U∞=1；可配置入口见下节。Cf使用实际U∞，这些量的数值收敛、近壁分辨率和物理精度需另证。绘图入口`tools/visualization/render_flat_plate.py`读取已通过审核的研究JSON，不将其视为TMR认证。

### 原生矩形近壁网格与精确残差诊断

`grid/RectilinearMesh2D.hpp`的`makeRectilinearMesh2D(x,y)`使用有限严格递增坐标数组；返回真实二维四边形共享拓扑，所有外边为DomainBoundary，物理边界由求解器另行指定。仅限完整矩形流体区域，不接受固体裁切。核心执行原默认Solver质量检查；FVM读取器继续检查几何/面关联。非法坐标、退化面积、尺寸溢出或质量失败显式报错。

```sh
mkdir -p outputs/native-flow/my-graded-plate
build/cartmesh2d_rectilinear_probe 32 16 4 outputs/native-flow/my-graded-plate/mesh
build/cartmesh2d_sst_rans_probe outputs/native-flow/my-graded-plate/mesh.cm2d outputs/native-flow/my-graded-plate/flow flatplate
python3 tools/verification/verify_sst_rans.py --mesh outputs/native-flow/my-graded-plate/mesh.cm2d --prefix outputs/native-flow/my-graded-plate/flow --output outputs/native-flow/my-graded-plate/audit.json
```

诊断probe默认单位正方形，可附加`xmin ymin xmax ymax`四个有限边界值；必须满足xmax>xmin、ymax>ymin。nx/ny各2..256，stretch在0..20；0为均匀，正值采用`expm1(a*j/ny)/expm1(a)`。默认半域前缘要求nx为偶数；自定义域需自行对齐前缘，否则跨越面会由平板边界检查拒绝。probe维数限制是资源保护，不是库API上限。更强加密仍可能触发默认网格质量门或求解精度失败，不能把本例参数当成通用预设。

线性系统恢复在`fv/detail/FlowLinearSystem2D.hpp`：精确影子正交breakdown重启，近舍入尺度使用补偿逐行`b-Ax`，对角占优时可作最多8次坐标校正，原双重残差门不变。补偿算法依据[Ogita、Rump、Oishi 2005](https://doi.org/10.1137/030601818)的误差分解原理独立实现；没有复制外部求解器代码或引入库。它不保证任意矩阵收敛或任意绝对容差可表示。Krylov步数不含坐标扫描，性能比较须同时报告总时间。数学breakdown背景见[Netlib Templates](https://www.netlib.org/templates/templates.html)。

真实对照及失败见当前状态。`artifacts/current/native-flatplate-grading.json`的`valid=false`明确表示四例没有全部通过；`cases`是有独立审核的完成结果，`failures`保留未通过项，`geometryChecksPassed`只表示网格几何。渲染工具`tools/visualization/render_flat_plate_grading.py`会同时展示这两类结果，不把失败网格配上伪造流场。


后续增加的`SparseSystem2D::solveCandidate()`返回显式`LinearCandidate2D{high,low,iterations}`。low仅在精度风险触发时分配；通过重新计算真实残差，后启用仍能校正之前丢失的低位。`relaxedDouble()`在高低位中完成松弛后才舍入；调用方必须再次检查实际double字段。普通`solve()`没有这种候选语义，要求仍作用于真实返回的double值。没有隐式转换可悄悄丢弃low，也不保证任意病态矩阵、上溢/下溢或任意精度要求都可解。

`ScalarTransport2D`的最终面通量/守恒检查继续保留；接近舍入尺度时，增加当前deferred source下的补偿原矩阵检查，结果在`ScalarTransportIteration2D::matrixAudited/matrixResidualNorm/matrixMaxDiagonalScaledImbalance`。它是附加门，不能替代原非线性通量；SST失败信息区分该门。两份矩阵的high+low通过、rounded失败，以及最终三个流场通过/一个失败，分别保存在`native-flatplate-precision.json`，不能混成单一PASS。

数值方法依据误差分解和额外精度残差/更新的公开原理独立实现，未复制第三方求解器代码；背景见[Error Bounds from Extra Precise Iterative Refinement](https://www.netlib.org/lapack/lawnspdf/lawn165.pdf)。本文实现是保持两部分候选的BiCGStab，不等同于论文的完整算法或误差上界证明。复现原冻结算法仍用本节同一网格与probe命令；ny=64的原失败保留。下节新增显式选择的迭代组织，最终完整非线性方程门不变。

### 显式稳态初值与有限次数的SST校正

`solveSteadyScalarTransportFromInitial2D(mesh, problem, initial, controls)`接受每个单元一个有限值作为**稳态初始猜测**；不把它当成上一物理时刻，不加入时间项。原`solveScalarTransport2D`接口及零初始默认路径保留。新接口可以返回未收敛的有限次校正结果，调用方必须检查`converged`，不能直接称其为解。

`SstTransportControls2D::scalarCorrectionsPerUpdate`默认0，保持充分求解每个冻结输运子问题。正值是实验性稳态Picard/deferred-correction组织：从当前k/omega开始，最多作指定次数的标量校正，再松弛、重建梯度/闭合、评估原完整非线性方程。内部线性失败仍抛错，标量中间状态不作为最终接受依据；最终返回的fields是实际double字段的完整重算。公共冻结求解API仍必须充分收敛，新组织拒绝dt/previous混用；不裁剪负k或非正omega。当前诊断probe只开放0和1，因此独立审核器也只接受这两个声明值；库API可指定更大正整数，但未因此获得验证。

分离输运方程、更新系数的组织可参考[OpenFOAM官方kOmegaSSTBase源码](https://api.openfoam.com/2506/kOmegaSSTBase_8C_source.html)。本实现独立编码，模型固定SST-2003m、无裁剪，不能说成对OpenFOAM模型/算法的完整复刻。完成冻结方程与收敛完整非线性方程是不同阶段；不得把中间未收敛改标成通过。

```sh
mkdir -p outputs/native-flow/my-graded-plate
build/cartmesh2d_rectilinear_probe 32 64 4 outputs/native-flow/my-graded-plate/mesh
build/cartmesh2d_sst_rans_probe outputs/native-flow/my-graded-plate/mesh.cm2d outputs/native-flow/my-graded-plate/flow flatplate-sweep
python3 tools/verification/verify_sst_rans.py --mesh outputs/native-flow/my-graded-plate/mesh.cm2d --prefix outputs/native-flow/my-graded-plate/flow --output outputs/native-flow/my-graded-plate/audit.json
# 原默认算法仍使用 flatplate；通道可用 channel-sweep 选择相同校正组织。
```

`native-flatplate-bounded-corrections.json`保留四例实际审核、旧二进制同输入比较、单纯warm-start的失败尝试和原2048格的高精度复算。两个时长不含网格读取/导出，都是单次观测；输出SHA与执行命令保留，审核器只验证模式标签的合法性，不从末态推断其执行历史。三档仅法向细化，不能冒充完整网格无关性；该实验仍为Re_plate=500，未获得高Re、瞬态或十万格SST资格。

### SST诊断的物理输入与相似性

`cartmesh2d_sst_rans_probe mesh prefix [mode] [options]`保留全部旧模式和默认参数，新增下列成对选项：`--nu`（运动黏度m²/s）、`--speed`（入口速度m/s；通道为最大速度）、`--inlet-k`（m²/s²）、`--inlet-omega`（1/s）、`--leading-edge`（m，只用于平板）。输入必须有限；nu/speed/omega>0，k>=0。缺值、重复、未知或不适用选项显式失败；前缘必须处于域内且与实际网格面端点对齐。不从任意“湍流强度”自动猜测另一项，没有隐式物理参数预设。

```sh
mkdir -p outputs/native-flow/highre
# 此例板长2、单位长度Re=5e6，故板长Re=1e7；仍是诊断，不是TMR认证。
build/cartmesh2d_rectilinear_probe 80 32 9 outputs/native-flow/highre/mesh -.5 0 2 1
build/cartmesh2d_sst_rans_probe outputs/native-flow/highre/mesh.cm2d outputs/native-flow/highre/flow flatplate-sweep --nu 2e-7 --speed 1 --inlet-k 2.25e-7 --inlet-omega 125 --leading-edge 0
python3 tools/verification/verify_sst_rans.py --mesh outputs/native-flow/highre/mesh.cm2d --prefix outputs/native-flow/highre/flow --output outputs/native-flow/highre/audit.json
```

审核器由原网格和实际物理输入独立计算入口/壁面条件、面黏度、质量/动量/原湍流方程及Cf/y+；停止门仍固定原值，不能用修改JSON容差绕过验收。`physicalInputs`、`bounds`和`plateReynolds`写入审核结果。JSON配置是待审核问题的声明，不是外部物理正确性的证明；运行命令、源码及工具哈希另存以保持来源可追溯。平板须完整矩形，任意非矩形或任意边界patch不在此入口范围内。

`sst_rans_verifier_test.py`增加坐标平移/尺度变换相似性：x长度×2、U×2、nu×4、k×4、omega不变，要求无量纲场及Cf/y+一致；又故意篡改物性、速度、入口湍流和边界配置，须被重算拒绝。旧固定物理量改为可配置没有放松这些检查。矩形输入与参数格式也有失败回归。

72cb13a阶段的40×32/stretch8高Re算例未通过梯度独立重构，原`native-sst-physical-inputs.json`仍保留该失败。后续解析矩形中心修复与80位独立几何复测使**重新计算**的同一案例通过；原错误场继续拒绝。新证据为`native-centroid-stability.json`，不得覆盖旧文件。系统clang宏验证Apple ARM64的long double与double均为53位有效位，不能假设long double足以防止这类几何舍入。


### SST初值与固定速度场对照

`cartmesh2d_sst_rans_probe`新增成对的`--initial-k K --initial-omega W`，要求有限K≥0/W>0；省略时保持入口值初始化。该选项只向既有原生初值API传入单元数组，入口和壁面仍来自原工况。输出的initialK/Omega为可检查的运行声明，不可仅凭最终场证明；独立审核保持所有原方程及收敛门，另将声明放在declaredInitialTurbulence。

固定速度场复现实验见`artifacts/current/native-sst-coupling.json`：其中reproductionSources保存转换器及两个小C++诊断程序的完整源码与哈希（临时编译，未加入生产核心），runs保存实际命令、history抽样及全文件哈希。需要本地保存的原网格和900次失败场；没有将大网格塞入Git。frozen.cpp使用系统clang++、C++20、Release原生静态库，固定u/v和通量后调用原SST输运。临时程序用于这一已核对输入，不是通用不可信文件导入器。

Anderson单历史试验依据[Walker与Ni，2011](https://doi.org/10.1137/10078356X)，用固定初值尺度加权的两次更新差求一个混合系数，再检查正性和原局部残差下降。所有候选均未采用；没有集成到求解器，也没有把探索用局部残差筛选冒充完整收敛门。下一步非线性方法仍须验证最终原方程与同一物理解。

### SST集中更新、原方程评估与协作停止

`SstRansControls2D::turbulenceCompletionUpdates`默认0，保持原交错调度。大于0时，前一SIMPLE记录的momentumResidual/velocityChange/pressureChange小于flow.tolerance、continuity小于1e-8且当前globalRelativeImbalance小于1e-8，才将本次更新预算设为max(普通更新数,集中更新数)，仍受turbulence.maxIterations限制。该条件只安排工作；更新后当前黏度、动量和两条原湍流方程重新验收，不能只凭前一步或标量通过接受整体结果。history追加completionUpdate；探针`--completion-updates N`允许整数0..500，显式启用。

`evaluateSst2003mTransport2D`接收固定载流/边界和给定k/omega，重新构造梯度、应变及闭合，调用原方程evaluate路径；忽略传入缓存梯度/应变，不修改输入或解线性系统。两个标量收敛标志不等于动量/RANS收敛。previousK/Omega与dt含义沿用非线性输运API；解析后向欧拉原方程根和非原方程冻结根的拒绝有回归。

`FlowControls2D::stopRequested`为空时保持原行为。回调在完整SIMPLE迭代后的正常收敛判断之后执行；true返回stopped=true、converged=false及该完整末态，仍进行最终面/力装配。若所有收敛门已通过，不被同时到来的停止请求改为失败。回调异常原样传播；不能中断内层线性/湍流求解。探针`--max-seconds T`接受有限正秒数，从求解入口计时，正常停止写.unconverged场、diagnostics.stopReason=time-budget并退出1；不代替外部硬超时，也尚未接入桌面取消按钮。

```sh
# 明确预算的诊断示例；不改变入口、初值、方程和停止门。
build/cartmesh2d_sst_rans_probe /path/mesh.cm2d outputs/run/flow flatplate-sweep --nu 2e-7 --speed 1 --inlet-k 2.25e-7 --inlet-omega 125 --leading-edge 0 --scalar-preconditioner ilu0 --completion-updates 500 --max-iterations 20000 --max-seconds 540
```

独立读取器拒绝stopped与converged同时为true、非法调度声明和超限更新。每次集中更新需有前一history的流动条件；逐次globalRelativeImbalance未保存，因而不能独立证明这部分触发轨迹，最终质量守恒依旧独立计算。旧文件缺少新增字段时按未启用处理，所有原方程门保持。

`native-sst-completion.json`保存12,800格原调度10,000次仍失败及集中更新3,866次通过的实际证据。纯omega损耗线性化和局部有限差分Newton为未采用的失败探索；晚期固定载流317步通过只解释调度依据，不单独构成耦合验收。实际图由`tools/visualization/render_sst_completion.py --study artifacts/current/native-sst-completion.json --output outputs/sst-completion.png`读取哈希匹配的真实网格/场/历史生成，需本地原数据。仍为实验SST探针，未提供GUI或通用物理精度资格。

### SST未收敛现场与最差单元

`FlowIteration2D`额外记录`momentumWorstCell`及同一单元有符号的`momentumResidualX/Y`，归一化仍为(diagU+diagV)*speed；二者hypot对应原momentumResidual。profile时另外记录压力修正前的`momentumPredictorResidual`（两分量max |b-Au|/(diag*speed)）及其单元，`pressureLinearResidual`为本次各PCG调用真实残差norm的最大值除以speed*shortestFace。这两个线性诊断使用补偿行残差，未作为新验收门；计时包含profile诊断开销，不能和旧无该诊断二进制的时间直接当提速比较。

SST探针的`--turbulence-updates N`须整数1..500，控制既有每次流动更新允许的最多湍流更新，默认1不变；实际次数还受已有局部收敛及maxIterations限制。正常返回但未收敛时写`prefix.unconverged.json/.cells.csv/.faces.csv`，JSON中converged=false；history/diagnostics沿用原prefix。任何上述输出已存在都拒绝复写。内部异常/硬超时没有保证最终场，不自动接受部分文件。

```sh
# 保持普通audit的拒绝规则；只有明确诊断才读取该未收敛现场。
python3 tools/verification/verify_sst_rans.py --mesh /path/mesh.cm2d --prefix outputs/run/flow.unconverged --output outputs/run/diagnostic.json --diagnostic
# 诊断生成成功仍返回2、valid=false；不是普通验收的通过退出码。
```

诊断模式仅把原收敛门的失败汇集到`failedConvergenceChecks`，不抑制格式、物性、几何、本构、通量、history一致性失败；即使输入原本收敛，diagnostic也不会输出valid=true。独立报告给出动量/k/omega的最差单元及坐标。未收敛bundle共享原prefix的history，使用`.unconverged`后缀定位；默认audit继续拒绝converged=false。

独立`verify_native_flow.polygon`使用原输入float的精确Fraction多边形矩，统一覆盖所有有效单元，避免相邻单元落在精度切换阈值两侧产生一个ULP中心错位。初始非有限/正面积检查保持；面积/中心仅在最终输出舍入，原生FVM不调用此Python参考。实际cell4001（长宽比小于32）的十六进制反例及顶点轮换回归保留在flow_verifier_test；参考精度修复不代表SST已收敛。

### 标量稀疏工作区复用

`ScalarTransportWorkspace2D`显式拥有稀疏图、矩阵与Krylov数组，可移动、不可复制；只在空闲时移动，不可并行共享或在调用中销毁。同一工作区的递归调用显式拒绝，异常退栈释放使用状态。三个标量solve/evaluate入口和三个SST输运入口的最后一个可选参数接收工作区指针；省略保留单次局部生命周期。`solveSstRans2D`在一次完整求解内持有同一工作区，供k/omega及后续SIMPLE调用顺序使用。

只按实际单元数和完整有序内部owner/neighbour连接精确比较是否可复用，不使用对象地址、哈希或单纯网格数量。连接变化重建，图对象与矩阵存放在稳定地址。每次仍验证当前网格、边界、物性和载流守恒，清零全部diag/off/rhs、失效旧分解，并重新组装数值系数；几何变化但连接相同不复用任何几何/材料系数。评估路径通常不需要CSR/Krylov，只有原有近表示精度检查需要时才构造CSR。缓存不改变非线性或线性收敛门。

profile新增`patternReuses`；`patternBuilds`和`patternReuses`是本次实际CSR请求计数，普通评估可能二者均零。ILU构建/复用计数取每次调用前后差值，不能把长寿命矩阵的累计数重复累加到SST总计。每次系数重装仍重做ILU，当前优化不是旧分解近似复用。

```cpp
ScalarTransportWorkspace2D workspace;
auto first = solveScalarTransport2D(mesh, problem, controls, {}, 0, &workspace);
// 更新 problem 后仍使用当前系数重装，工作区只保留结构和存储。
auto second = solveScalarTransport2D(mesh, updatedProblem, controls, {}, 0, &workspace);
```

### 标量ILU(0)与SST有界诊断

`ScalarTransportControls2D::preconditioner`新增`ScalarPreconditioner2D::{Jacobi,ILU0}`，默认Jacobi保持。ILU(0)用于非对称标量BiCGStab，按原CSR图自然序做零填充Doolittle分解及前/后代入；标准方法参考[Netlib Templates §3.4](https://www.netlib.org/templates/templates.pdf)，本实现没有复制库源码。它不是压力IC0，也不改变压力预条件选择。要求有限系数和正原对角/消元主元；不满足就失败，不做shift、重排序或隐式fallback，因而不是任意矩阵通用求解器。

每次solve先检查真实初始残差，仅确需求解时建立因子；完全相同diag/off可供RHS变化复用，系数变化重建。失败构造清除ready状态，不计成功构造；`add/reset`失效缓存。公开`preconditionILU0`必须先factor且快照相符，否则拒绝；Krylov内部私有apply在当前solve系数不变前提下省去重复快照比较。原补偿b-Ax、norm/逐格门、高低位候选、实际场四舍五入后的最终验收不变。独立稠密掩码Doolittle、丢弃fill、非对称已知解、失效缓存、非正主元/NaN/Inf和表示精度反例均覆盖。

SST探针提供`--pressure-preconditioner jacobi|ic0|aggregation`（默认ic0，与旧版本相同）、`--scalar-preconditioner jacobi|ilu0`和`--max-iterations N`（整数1..20000，默认2000）。例如上一节高Re命令末尾加`--scalar-preconditioner ilu0 --max-iterations 100`，只用于固定工作量诊断，可能非零退出。返回后始终写`.history.csv`及`.diagnostics.json`，记录converged/stopReason、实际迭代、方法和分项计时；未收敛不会写普通`.json/.cells.csv/.faces.csv`，另保留上述显式`.unconverged.*`诊断场。硬超时或内部异常可能没有这些最终诊断，应保留stderr。压力方法同时写入diagnostics及场元数据；独立审核只验证声明合法并重算原方程，不能从末态反推实际算法，执行命令与二进制哈希须一起保存。旧证据没有该字段时按旧默认ic0解释。固定步数计时不代表全程收敛速度；同输入、同迭代和容差下串行交错重复，失败或超时单独记录，禁止混入完成组的提速比。prefix已有任一上述文件时拒绝，避免旧场冒充新结果。

审核器验证方法声明合法，仍从真实场重算原方程；不能由末态独立推断实际执行的是哪个预条件器，执行路径另由命令/二进制哈希记录。`kSolves`与`omegaSolves`互斥地组成`scalarSolves`；`ilu0Builds/ilu0Reuses`只统计实际成功factor，不强求其等于solve调用数。计时和Krylov次数不构成物理验收。实际比较与12,800格完整超时保存在`native-sst-ilu-performance.json`；图用`python3 tools/visualization/render_sst_ilu.py --study artifacts/current/native-sst-ilu-performance.json --output outputs/sst-ilu.png`重建，需要对应本地网格及CSV。

### SST近壁规模诊断与计时

`ScalarTransportControls2D::profile`默认关闭。结果`ScalarTransportPerformance2D`记录calls、patternBuilds、linearIterations、total/setup/linear/faceFlux秒；SST按`scalarSolves`和`scalarEvaluations`累计，RANS由`flow.profile`统一启用。setup包含验证、初始组装及按需的延迟CSR构建；total包含全部子项，Krylov迭代不包含既有补偿细化扫掠，其耗时仍包含在线性时间中。初始和每次更新后的k/omega评估均计入，零线性迭代不意味着零成本。

纯评估保留对角/RHS和原面方程，正常时不构建CSR/线性工作区；接近表示极限时仍重建当前矩阵并执行补偿b-Ax，不能复用候选或上一轮的矩阵。`FvMesh2D`诊断辅助函数使用string_view避免成功校验构造临时字符串，所有校验照常执行。错误文字只在抛异常时转成拥有存储的string；临时字符串参数在本次调用结束前有效。

性能复核见`native-sst-setup-performance.json`，分配计数小程序源码也在该证据内。画图：`python3 tools/visualization/render_sst_performance.py --study artifacts/current/native-sst-setup-performance.json --output outputs/sst-performance.png`，需要证据对应本地CM2D/CSV；图只展示已独立审核且前后字段一致的场。实际数值变化与当前状态见唯一状态文档。


`SstRansPerformance2D`由`flow.profile`启用，默认关闭；不影响数值路径。`updates`是本构回调次数，`updateSeconds`包含`transportSeconds`、返回场梯度重算及一次壁距构造。probe输出`performance`对象，另列动量/压力线性求解时间和迭代数。不要相加包含式阶段，也不要将有界更新次数解释为已收敛次数。失败抛错/进程超时没有最终性能结果；计数与时钟不是独立物理审核器的验收内容。

`rectilinear_probe`允许每轴2..4096的整数，总单元仍≤65536；原质量门不变。高Re近壁诊断可将上述例的网格改为`160 32 10`，其他输入相同；当前该5,120格已完成独立审核。`400 32 11`的12,800格网格通过几何与Solver检查，但求解在90秒预算内未完成，不应当成已验收的快速配置。研究证据`native-sst-wall-scale.json`保留失败，不自动延长预算；场prefix与审核summary路径必须不同以免覆盖元数据。


### 薄矩形中心与独立核算精度

`Polygon2D::centroid`保留原面积阈值和通用路径；仅当四个顶点覆盖四个**精确**包围盒角点、各边严格轴向且非零时，返回`std::midpoint`解析中心。不按tolerance近似识别，不改变signedArea、拓扑或质量策略。它统一薄矩形中心的中点舍入，避免行间细小横向偏移被大法向梯度放大；不承诺修复所有非矩形或所有条件数问题。原生几何变化可能使旧checkpoint的严格几何签名不匹配，此时仍明确拒绝。

此前`verify_native_flow.polygon`在span²/area>32时使用80位Decimal；最新已改为上节的统一精确有理数多边形测量，以修复阈值以下邻格的一个ULP错位。几何/CFD允许误差不变。测试包含实际半ulp中点、轮换/反向/平移、非对称薄多边形和拒绝规则；独立读取器没有复用原生解析识别代码。

`verify_sst_rans`输出两种残差：`maxCellResidual`是从实际字段和独立本构计算的**未拆分原方程**；`reportedTransportBalance`是已逐项独立核对过的导出binary64 source/loss/通量按标量API报告顺序形成的余额。前者按原门决定方程资格，后者按原摘要容差核对history；`originalEquationCellDifference`保留运算顺序/舍入差异。不能用history小、或split余额小，跳过前者。导出系数、通量及history篡改仍由对应检查拒绝。

修复案例可用上一节命令将nx改40、stretch改8复现。查看原失败文件与新最终场需要区别其源码/可执行文件哈希；不应把旧文件重新标为已通过。绘图脚本继续只展示有独立审核的真实场。

### 保持近壁层的精确法向细化审核

`tools/verification/verify_rectilinear_refinement.py`独立读取两份最终CM2D并执行原几何测量，要求每份都是完整的轴向矩形张量积网格，不能带重复/缺失单元、斜边、交叉或退化四边形。流向坐标必须逐值相同；`--preserve-first-rows N`要求正整数且留下可细化区。前N个法向区间及其单元多边形保持不变，其余每个区间只能插入一次精确有理数中点舍入后的binary64坐标，所有原节点仍存在。面积证明针对解析后的坐标，不推断CAD精度或Solver质量。

```sh
python3 tools/verification/verify_rectilinear_refinement.py --original /path/original.cm2d --refined /path/refined.cm2d --preserve-first-rows 4 --output outputs/refinement.json
python3 tools/visualization/render_sst_normal_refinement.py --study artifacts/current/native-sst-normal-refinement.json --output outputs/sst-normal-refinement.png
```

24,000格实例由临时C++生成器调用既有原生`makeRectilinearMesh2D`，原Solver门不变；源码/命令/哈希保存在证据JSON，未增加绕过质量门的网格入口。1,600个近壁单元原样保留，其余11,200个各二分成22,400个。`geometryValid`只代表上述关系，不授予流动/物理通过；必须另运行原Solver检查及严格场审核。`render_sst_flatplate_analysis.py --mesh-case NAME`可选择实际显示哪份已审核网格，默认fine兼容旧调用。

本轮完整求解和独立物理量比较在`native-sst-normal-refinement.json`；方程全部通过，但Cf对中外部法向细化仅约0.2%敏感，仍非网格无关证明。原始场保留outputs，性能分项来自实际求解profile，不将标量setup全归因于CSR构造。

### 实际平板剖面与域敏感性测量

`tools/verification/analyze_sst_flatplate.py`读取manifest的cases（name/mesh/prefix）及comparisons（两两名称）。每例先从实际文件运行严格`verify_sst_rans.audit`，再测量；拒绝未收敛、任意非矩形/非完整张量积网格、非法物性或站位外推。此工具不是通用Cut-cell剖面插值器。固定站位x=.97008、1.90334：各水平行按单元中心x线性插值；wall kinematic shear另按墙面中心插值，Cf=2tau/U²、u_tau=sqrt(tau)、u+=u/u_tau、y+=(y-ywall)u_tau/nu。负剪切不套用此缩放。delta99仅为从无滑移点开始的0.99U首次交点，无交点返回null，不声称测得边界层外缘速度。

跨网格剖面对照在共有高度范围内取100个对数分布点（最高.05），两侧均线性插值、不外推，以U归一化速度差，避免把不同单元编号当同一物理位置。积分Cd=2sum(tau*面长)/(U²*实际板长)。同模型敏感性要求物性、入口k/omega、前缘、顶部条件及离散声明相同；不自动推出误差阶或通过阈值。顶部域几何隔离另由实际多边形精确对照提供证据，不从cell数量或文件名猜测。

```sh
python3 tools/verification/analyze_sst_flatplate.py --manifest outputs/native-flow/sst-physical-check/manifest.json --reference-dir outputs/native-flow/sst-physical-check/reference --output outputs/native-flow/sst-physical-check/analysis.json
python3 tools/visualization/render_sst_flatplate_analysis.py --study artifacts/current/native-sst-physical-check.json --output outputs/sst-physical-check.png
```

参考数据来源为下节TMR页面的[墙面Cf](https://tmbwg.github.io/turbmodels/FlatPlate/SST/cf_plate_sstv.dat)、[速度剖面](https://tmbwg.github.io/turbmodels/FlatPlate/SST/flatplate_u_sstv.dat)、[不可压Cf对照](https://tmbwg.github.io/turbmodels/FlatPlate_validation/cf_incomp_results_sstv.dat)。仅解析该类ASCII POINT zone表，验证列数/有限值/分区及插值坐标。速度表的u/U∞、y/L含义来自无量纲工况和远场数值的解释，文件头仅写u,y。参考比较仅限U/nu=5e6、k/U²=2.25e-7、omega/U=125、前缘0、板长2的物理尺度；即使满足，SST变体、可压缩性与边界差异仍阻止精度认证。输出始终physicalAccuracyQualified=false；不会因为数值接近参考而授予通过。

`native-sst-physical-check.json`保留三档旧网格的重新审核及新顶部扩域实际场，另保存精确保留下部单元的检查、生成器/监测源码、原生门、时间/RSS/磁盘、工具与参考SHA256。完整场留在outputs中，不加入Git。原始参考表的不可压结果是背景对照，不被自动挑选最接近值作为本项目标准。

### 标准湍流参考的适用边界

已核对[TMR 2DZP平板定义](https://tmbwg.github.io/turbmodels/flatplate.html)、[网格](https://tmbwg.github.io/turbmodels/flatplate_grids.html)及[SST参考结果](https://tmbwg.github.io/turbmodels/flatplate_sst.html)。平板x=0至2，参考长度1，Re_L=5e6、M=.2；网格35×25至545×385节点，需明确区分节点数与实际流体单元数。壁面omega及自由来流k/omega必须按该例指定，不能沿用诊断通道的数值。近壁y+、x约.97/1.90的剖面及壁面摩阻是后续关注量；不能只看残差。

该现成平板数据使用**SST-Vm**，当前实现为**不可压SST-2003m**。低Mach相近不等于方程/变体相同，不能直接以该表作为本模型严格误差标准。下一步须明确选择匹配参考，或如实分开报告跨变体参考和独立制造解验证；在此之前不宣称平板验收。混合底边与不可压压力远场的核心实现见上节，尚需各向异性近壁网格、匹配参考和实际进入流验证；没有平板产品入口。

### 非正交压力修正固定点

每个SIMPLE外迭代保留原来的最多4遍修正。`solvePressure()`返回0时没有修改pc，因此下一遍会逐位构造相同的梯度、修正项和RHS；此时保留当前面修正，跳过必然相同的后续遍。没有角度容差或新的停止条件。profile的 `pressureSolves` 只计实际调用，新增 `pressureCorrectionPassesSkipped` 计省略遍数，两者之和仍为 `4 * simpleIterations`（瞬态按全部物理步累计）。独立线性测试核对三种预条件器的零次求解不修改场，实际新旧CLI输出逐字节比较见 `artifacts/current/native-flow-pressure-fixed-point.json`。

### 压力aggregation分组与系数分离

`FlowAggregation2D.hpp`缓存每条跨组上三角边的粗层CSR目标位置。`refresh()`保持P不变，按原细层遍历顺序以long double累加当前 `P^T A P`，镜像写回保证精确对称，再重做最粗层LDLT。非有限值/符号/对称性错误显式抛出；相容性要求零模式不变、每项系数相对上一矩阵比值在[.5,2]内。`SparseSystem2D`最多连续刷新8次后重新分组；任何构造/刷新异常均丢弃缓存。这里的范围/次数是保守建层策略，最终PCG真实残差标准不变，不作为精度认证。

profile分别记录 `pressureHierarchyBuilds`（完整分组）、`pressureHierarchyRefreshes`（当前系数刷新）及 `pressureHierarchyReuses`（完全相同系数）。原有IC0默认不变。单测独立稠密P^TAP验证非均匀更新、V-cycle线性/对称/正性、已知解真实残差、零模式/比例/不对称拒绝、定期重建及数值失败后丢弃缓存。实际圆柱及三档涡流用独立CM2D/CSV读取器核对，profile计时不能替代解误差；证据见 `artifacts/current/native-flow-aggregation-refresh.json`。绘图入口 `tools/visualization/render_aggregation_refresh.py` 从真实单元读取速度，未用示意场替代计算结果。

### 当前迭代动量装配复用

0.4.11 将每轮末尾用于真实动量残差的未松弛矩阵和速度/压力梯度，直接用于下一轮 SIMPLE。字段、面通量和回流边界更新后先统一刷新；下一轮交换数值存储，再按原顺序施加速度松弛。对流、完整对称应力、时间项、边界、压力修正次数和所有停止条件不变；最后的面动量输出使用同一最终梯度。不同方程不共享可变预条件缓存，此处复用的是动量矩阵数值存储，压力求解器没有更换。

同一工具链以保存的 d5eb688 CLI 与当前 CLI 配对，比较 cells/faces CSV、残差、fields JSON、VTK及瞬态checkpoint/时间历史；独立工具另从CM2D和实际结果重算质量/动量/时间项。参数、二进制/输入哈希和原始命令见 `artifacts/current/native-flow-momentum-reuse.json`，没有用时长变化替代正确性。线性迭代超限现在报告真实残差、目标和逐行对角缩放残差，诊断信息不改变算法或停止门。

### 非定常层流与断点续算

开发分支CLI/核心提供一阶后向欧拉。桌面0.4.4起提供固定时间步、物理量监测、取消保留与checkpoint续算；0.4.18增加CFL控制的自动步长及失败重试，本地macOS已打包实测；不支持二阶时间格式、时间误差估计、瞬态湍流或移动网格。

CLI高级数值选项 `--velocity-relaxation`（同步热输运用 `--flow-velocity-relaxation`）控制每个物理时间步内部的速度松弛，范围 `(0,1]`，默认仍为 `.6`。它不改变物理时间步，不是精度或质量门；较大值有时减少外迭代，也可能使线性求解失败。显式选项仅允许非定常/同步模式，稳态或冻结载流拒绝；不会自动调大，桌面当前仍使用默认值。输出记录实际系数；续算允许改变这种数值控制，但同网格细化对照须固定它。不能把不同系数的有限迭代误差视为逐位一致。

```bash
build/cartmesh2d_flow_cli --mesh PATH/FINAL.solver.cm2d --output outputs/startup/result --case external --nu .05 --speed 1 --convection limited-linear --time-step .01 --steps 5 --tolerance 1e-9 --max-iterations 1000 --profile
build/cartmesh2d_flow_cli --mesh PATH/FINAL.solver.cm2d --output outputs/continued/result --case external --nu .05 --speed 1 --convection limited-linear --time-step .01 --steps 5 --restart outputs/startup/result.checkpoint --tolerance 1e-9 --max-iterations 1000
python3 tools/verification/verify_transient_flow.py --mesh PATH/FINAL.solver.cm2d --prefix outputs/continued/result --output outputs/continued/audit.json
```

同步热涡规模验证入口（串行，输出必须是新目录）：

```sh
python3 tools/verification/verify_thermal_scale.py --output outputs/thermal-scale-new --cells-across 100 200 320 --velocity-relaxation .8 --timeout 300
python3 tools/visualization/render_thermal_scale.py --study outputs/thermal-scale-new/study.json --output outputs/thermal-scale-new/preview.png
```

可用 `--flow-tolerance 1e-10` 做显式迭代敏感性对照，生成默认仍为1e-8；同一值传入transport和独立flow，不能将不同容差混为同一空间细化序列。结果分别记录coupledFlowTolerance、flowTolerance，以及标量三个停止条件。transport JSON的flowTolerance只在同步模式表示实际载流控制，冻结模式为null。历史同步结果缺失此字段可用不带期望容差的 `--audit-prefix` 单例复查，标记unknown，严格系列拒绝；显式给出 `--flow-tolerance` 时要求两边实际元数据都匹配，不用期望值补造历史记录。新驱动的这一选项不改变App默认设置，也不是全软件精度及格线。

默认步长 `.005`、2步、nu=.1、D=.02、Uref=1，固定物性与有限时间；不是长时间或湍流验证。默认系数仍为 `.6`，上面 `.8` 是显式对照设置，不能推广成所有工况推荐值。单进程时间预算默认180秒，可显式设至600秒；磁盘不足1.5 GiB前停止；超时/失败留存并返回非零。时间预算只控制本机资源，不是数值精度门；按实测设置420秒的新研究不能声称满足旧300秒预算。`--reuse-mesh-study PATH`可读取此前同工具生成的 `nN/mesh/square.solver.cm2d`，不复制大网格，记录输入SHA。审核真实单元数、几何、最终逐面/逐格守恒、时间项和面通量；同步载流必须与同设置独立求解的checkpoint逐字节相同，各步接受历史和迭代数相符。未独立存档所有中间场，不能声称逐步全部重算审核。

对这一个相切解析涡，标量 `sin(pi*x)sin(pi*y)` 的对流项解析为零；与速度的衰减率分别为 `2*pi²*D` 和 `2*pi²*nu`。同时报告连续指数解和空间连续、时间按后向欧拉的参考幅值 `(1+rate*dt)^(-steps)`。前者含时间离散误差；后者有助区分网格细化效果，但仍含空间、载流和迭代误差。三档校验固定物性/步长/终止时间/松弛，报告观测阶与误差趋势，不设通用工程精度分数。绘图读取实际最终CM2D及字段，失败研究不会被渲染成通过图。

失败研究只在显式加 `--allow-incomplete` 后可绘制已独立审核的完整案例，图中保留醒目的未完成提示，不能拿中途checkpoint当最终场。规模工具在POSIX启动独立进程组，超时停止整个组；Windows使用有界taskkill及失败记录（本次未实测Windows）。资源计时器本身失败也保留非零状态；不能把缺失RSS当成0，或把封装器错误直接解释为方程失败。

`--time-step`和`--steps`必须一起提供；后者表示本次追加的步数，`--max-iterations`是每步内迭代上限。普通channel/cavity/external从静止开始，在t>0施加入流/顶盖速度；这是瞬时启动，会产生启动压力，不是预先求稳态再贴上时间标签。外流仍限定矩形外域、固定固体和无回流出口。

动量添加 `V*(Unew-Uold)/dt`；上一接受时刻的速度和唯一面通量保留，Rhie–Chow包含旧时刻与内松弛的插值缺陷修正。时间步内部原动量残差、速度/压力变化和质量守恒都达到原停止条件才接受。动量线性求解除既有全局真残差条件外，还约束每行 `abs(b-Ax)/aP <= .01*tolerance*Uref*alphaU`，避免远场大格子的右端项掩盖小Cut-cell的局部残差；没有放宽非线性门。CFL为每个单元 `dt*sum(abs(phi))/(2V)` 的最大值；固定步长模式仅作诊断，自动模式要求试算的实际CFL不超过指定上限。

桌面实现：`core/flow-checkpoint.js` 流式读取配置与时间供界面使用，不替代原生完整状态和网格核验；`core/flow.js` 校验请求、物理进度、摘要和时间历史。main 先写独立待验目录，全部成功再替换完整结果；失败或取消保留诊断与最后接受状态。原完整流场和最新续算状态分别标注时间。稳态调用参数保持不变。桌面物理监测进度拒绝非数值JSON字段。

实际打包验证入口在原有 `--smoke=circle --flow=external` 基础上增加 `--flow-dt=.01 --flow-steps=2 --flow-resume-steps=2 --flow-failure-check=true --flow-cancel-check=true`。完整参数、App/原生哈希、独立导出读回和连续CLI对照见 `artifacts/current/desktop-transient.json`；大文件保留在忽略提交的 `outputs/native-flow/transient/desktop-gui-delivery/`。

除已有六份结果外，`.time-history.csv`逐步保存物理时间、是否接受、内迭代次数、残差、最大CFL、动能和力；`.residuals.csv`仅含最后尝试时间步的内迭代。`.cells.csv`新增previousU/V及temporalX/Y，后者是积分时间项。摘要区分候选time与acceptedTime、请求与完成步数。退出2代表该候选步未收敛，不能拿它继续时间推进；异常退出1会将已开始运行的摘要标成failed。

`.checkpoint`只保存初始/已接受状态，包括面通量。采用临时写入后替换，重启逐项核对网格几何、owner/neighbour、编号、物性、工况和离散格式；可以更改步长、停止容差及迭代预算。它不是跨编译器逐位一致性承诺，也不具备密码学防篡改功能。取消/崩溃后须读取最后完整checkpoint，不能把候选CSV当成重启场。重启回归将连续5步与2+3步的最终场和checkpoint逐字节比较。

`--case taylor-green`是只在CLI使用的无源解析验证工况：完整单位方域、四壁无穿透自由滑移，`A=Uref*exp(-2*nu*pi²*t)`，`u=A*sin(pi*x)*cos(pi*y)`、`v=-A*cos(pi*x)*sin(pi*y)`、`p=A²/4*(cos(2*pi*x)+cos(2*pi*y))`。运动学压力减去cell0解析值；初始面通量由解析流函数沿真实边端点之差积分。解析微分/能量/边界条件另有测试，避免把错误压力符号作为“真值”。

`run_transient_flow.py --mesh FINAL.solver.cm2d --output NEW_DIR`串行计算默认四档时间步，固定终止t=.2、nu=.1、速度1、限制线性对流和容差1e-9，每档限180秒；输出目录必须新建，保留实际命令、哈希、返回码与耗时。可选`--case/--dt/--end-time`、`--convection`及`--pressure-preconditioner`，用同一入口核查启动工况。`verify_transient_flow.py`独立读回几何、时间积分、逐面压力/对流/应力、局部及全局动量/质量、CFL和能量；元数据/历史不一致时失败。旧稳态验证器明确拒绝非定常结果，避免误套稳态基准。

最终矩阵图可用 `python3 tools/verification/plot_transient_flow.py --root outputs/native-flow/transient/final --output artifacts/current/native-flow-transient` 重建；先校验来源哈希。误差和时间步自收敛单独报告，离散守恒通过不等于物理精度验收。

### 压力多重网格（实验选项）

`cartmesh2d_flow_cli --pressure-preconditioner aggregation` 在原PCG中改用聚合多重网格预条件。默认仍是`ic0`，`jacobi`保留对照。非定常重启允许改变预条件器，因为它不改变物性或离散方程；它仍须通过原矩阵真实残差 `1e-13 + 1e-11*||rhs||₂` 检查。改变迭代路径会引入舍入差，不能要求不同方法逐字节相同，也不能据此放宽非线性或守恒门。

实现位于`FlowAggregation2D.hpp`：按最大负耦合确定性配对，再将有连接的单独节点并入相邻配对组；分片常数延拓P、转置限制Pᵀ，粗算子为PᵀAP；内部面贡献以同一个累加结果写入两侧，保持精确对称。每次PCG应用固定一轮V-cycle，前向/反向Gauss–Seidel构成对称平滑，末层至多32单元用稠密LDLᵀ；纯对角末层直接求解。采用[标准多重网格与PCG原理](https://www.netlib.org/templates/templates.html)，[hypre官方说明](https://hypre.readthedocs.io/en/stable/solvers-boomeramg.html)也强调CG需要对称平滑。本仓库独立实现上述有限方案，没有复制或链接hypre；不宣称与成熟AMG库功能或鲁棒性相当。

当前限对称、正对角、非正非对角的压力矩阵；分解失败、非有限数或每层聚合未至少减少约20%未知数时明确失败，提示选IC0，不静默替换方法。深度上限32，防止弱聚合图反复保留大层。实际层数、末层单元数、构建/复用次数随`--profile`记录；这些是资源诊断，不是物理精度等级。

同一次SIMPLE的四次非正交修正共用矩阵；IC0/多重网格都由矩阵持有缓存，用精确系数比较防止直接写入导致缓存过期。下一轮重新装配会使缓存失效；分解失败不保留可用标记。多重网格每层工作数组预分配，V-cycle不逐步申请大数组。缓存和层次结构增加存储，峰值RSS需实测，不能只报告迭代次数降低。

0.4.10 将粗层Galerkin装配的全局tuple排序改为按粗行分桶、行内稳定列排序，在原始细面顺序下以long double合并；同一个结果直接写入两侧CSR。聚合选择、粗层系数、非零项和对称性不变，去掉第二份双向tuple列表。首个前向平滑从零开始，只读取已覆盖的下三角值，因此无需清空整层工作向量或读取上三角的零值；反向平滑和完整残差不变。测试用独立稠密PᵀAP逐层核对原矩阵及置换矩阵，另测复用工作区和零右端。实际三个CFD输入另与旧二进制逐字节比较；不能由局部装配加速推断整体求解等比例加速。

### 相同输入的核心性能对照

`benchmark_flow_pair.py` 对保存的旧CLI和当前CLI串行执行相同的物理输入，每一对轮换前后顺序；记录二进制/网格哈希、完整命令、实际返回码、超时和原生迭代计数。macOS用`/usr/bin/time -l`、Linux用`-v`读取峰值RSS，不可用时记录null，不冒充零内存。比较最终单元字段和输出文件哈希；性能比较本身不替代独立几何/物理审核。`--baseline-preconditioner` 和 `--candidate-preconditioner` 可分别指定 `ic0/jacobi/aggregation`，默认均为 IC0；比较不同方法时记录实际字段差异，不能宣称逐位一致。

```sh
python3 tools/verification/benchmark_flow_pair.py --baseline PATH/saved-cli --candidate build/cartmesh2d_flow_cli --mesh PATH/channel.solver.cm2d --output outputs/new-performance-pair --case channel --nu .01 --convection upwind --tolerance 1e-6 --repeats 2 --timeout 180
```

残差范数仍为Euclidean L2，停止标准仍是`1e-13 + 1e-11*||rhs||`并回代原矩阵检查真实残差。改为带缩放的平方和，避免逐分量调用hypot；不直接累加double平方，以免极大/极小数溢出/下溢。这是[标准缩放范数思路](https://www.netlib.org/lapack/explore-html/d8/d76/group__lassq.html)的本仓库实现，没有引入LAPACK依赖，也不属于原创数值理论。累加仍使用long double，具体精度由平台决定；不同编译器逐位一致性需另证。

### CFL 自动步长与拒绝试算

`cartmesh2d_flow_cli --time-step MAX_DT --end-time ABSOLUTE_T`启用自动后向欧拉；不能同时给`--steps`。控制项为`--min-time-step MAX_DT/1024`、`--max-courant 1`、`--max-step-retries 10`及`--max-time-steps 100000`。目标时间必须晚于初始/续算已接受时间。App0.4.18的“非定常·自动步长”同步提供这些控制，温度联合推进仍使用原固定步长入口。

`FlowTimeStep2D.hpp`从上次已接受面通量计算Courant增长率，使用0.8余量预测下一步；实际试算需同时通过原非线性/守恒门和实际CFL上限。拒绝后从原状态以`min(.5,.8*CFLlimit/CFLtrial)`缩步重算，不更新检查点；最小步或重试/接受步预算耗尽会失败并保留最后接受状态。终点剩余时间可小于最小步以精确到达目标。核心推进方程和检查点v2/v4格式不变，自动控制策略不是物理状态；相同控制和目标、从接受步预算中断后续算能逐字节复现连续检查点。

`.attempt-history.csv`记录所有试算，含开始/候选时间、dt、是否接受、拒绝原因、内迭代停止指标及实际CFL；`.time-history.csv`在自动模式只含接受步，dt可以变化。摘要`timeStepControl=adaptive-cfl-retry`记录全部控制、起终点、attemptCount/rejectedSteps/completedSteps；`dt`是最后接受步长度，不再带固定模式的requestedSteps。独立Python审核与App读回均检查拒绝原因、状态时间不前移、缩步规则和接受历史的一致性；最终场仍独立重建时间项/动量/质量/CFL。

0.4.19避免把反复加法留下的舍入尾差单独推进成近零时间步：末步在最多64倍机器精度且相对dt不超过1e-12的范围内取真实剩余时长并重新计算；若正常尾差小于最小步且允许，则把末段分成两个可解析的步长。实际CFL与原非线性/守恒接受门保持不变。

`tests/adaptive_flow_cli_test.py`覆盖CFL和内迭代两种拒绝、最小步/重试耗尽、接受步预算续算逐字节一致、非法组合与篡改历史拒绝。App失败或取消保留上次完整显示和最后接受检查点，ZIP含全部当前完成运行试算及未完成诊断。CFL控制不构成时间误差估计，物理精度仍须时间与网格细化。

当前涡衰减的真实场与误差图：`python3 tools/visualization/render_time_accuracy.py --root outputs/native-flow/current-time-accuracy --output outputs/native-flow/current-time-accuracy/reproduced --binary-snapshot outputs/native-flow/binaries/flow-0.4.18`。脚本重新独立审核两个四档序列、加严容差和自动CFL对照，输入缺失或被篡改会失败；图中的自收敛阶不替代通用工程精度验收。

### 显式初始局部涡

0.4.20的`FlowInitialization2D.hpp`给新的物理非定常工况增加一次性初始条件：`--initial-vortex-x X --initial-vortex-y Y --initial-vortex-radius R --initial-vortex-speed V`四项必须齐全；V是带符号峰值速度，正值逆时针，R和中心坐标使用米。不能与稳态、解析制造解、边界模板或`--restart`一起使用。支撑圆盘必须严格位于流体域内且远离所有固体壁面、入口和出口；通过有向边界绕数及点到线段距离检查，包括孔洞。非零扰动若在单元中心或面端点上完全无法解析，会明确失败。

定义`s=1-|x-centre|²/R²`，圆盘内`psi=V*R*(25*sqrt(5)/96)*s³`，盘外为零；`u=dpsi/dy,v=-dpsi/dx`，速度峰值在`r/R=1/sqrt(5)`。真实共享面通量增加`psi(b)-psi(a)`，a/b沿owner多边形方向，因此每个封闭单元初始通量望远镜相消。边界通量、初始压力及后续方程不变，不是重复的源项或持续外力；当前仅支持流动，温度联合初始化后置。

摘要`initialVortex`记录定义与全部参数，`.initial.checkpoint`保存真正t=0的网格绑定U/V/P/FLUX，`.checkpoint`仍是最后接受状态。续算只读取状态，不再施加扰动；原v2/v4格式无需升级。App同步提供可折叠设置、续算锁定、独立初始文件和导出说明。`verify_flow_initialization.py`从CM2D重建有向端点、几何、解析速度/流函数和局部守恒，核对完整初始checkpoint；单步结果同时核对previousU/V与初始状态。多步审核不声称重建所有中间场。

`tests/flow_initialization_test.cpp`核对峰值定义、微分关系、散度、旋转协变与非法输入；CLI回归含连续与1+1续算逐字节一致、反向涡、圆环固体孔/穿壁/未解析拒绝和篡改检查点拒绝。实际计算证据见CURRENT_STATE。

### 固定网格的时间步比较

`tools/verification/compare_transient_steps.py` 从 `run_transient_flow.py` 的 `runs.json` 重新读取最终场与时间历史，不信任旧 PASS。构建路径上的程序更新后，可用`--binary-snapshot`指向保存的旧可执行文件；必须与原runs.json的SHA256相同，原执行命令保留，不重写为新命令。至少三档严格减半的时间步，要求同一网格、二进制、工况、物性、格式、容差及终止时间；逐项核对命令与独立读回，并要求从 t=0 开始。本入口暂不比较重启序列。按相同 cell id、实际面积计算相邻时间步速度差和观测阶；负阶也照实报告，`valid` 仅表示比较输入一致且离散审核通过，不是工程精度合格。

```sh
python3 tools/verification/run_transient_flow.py --mesh PATH/unit-square.solver.cm2d --output outputs/my-time-study --dt .04 .02 .01 .005 --end-time 2
python3 tools/verification/compare_transient_steps.py --runs outputs/my-time-study/runs.json --output outputs/my-time-study/comparison.json
```

默认案例 Taylor–Green 要求完整单位方域。其能量解析式为 `E=Uref²/4 exp(-4 nu π² t)`，长时间衰减应同时报告相对速度/能量误差与衰减率，不能只看绝对误差变小。中间步能量来自CSV监测，只有最终步从完整场独立复算；报告明确保留这个区别。步长较大可以收敛而仍有明显时间离散误差。

本轮复现图入口：`MPLCONFIGDIR=/tmp/cartmesh-mpl python3 tools/verification/plot_transient_comparison.py --root outputs/native-flow/transient-long --circle-runs outputs/native-flow/transient-long/circle-re20-completed/runs.json --output artifacts/current/native-flow-transient-long`。它重新审核两组时间步序列，并绘制真实圆柱网格与最终流场、解析衰减相对误差、圆柱受力和相邻时间步差。实际命令与耗时保存在证据JSON，未包含新的空间/域无关性或涡脱落精度认证。

### 完整流动方程制造解

`--case manufactured` 是命令行验证入口，不是用户物理工况，也不出现在桌面工况列表。仅接受无孔单位方形 `[0,1]²`；四边静止无滑移，压力从内部场外推，cell 0 固定压力为0。

令 `X=πx, Y=πy`，流函数 `ψ=(Uref/π)sin²X sin²Y`，则 `u=Uref sin²X sin(2Y)`、`v=−Uref sin(2X)sin²Y`，原始运动学压力 `p=Uref²[cosX cosY+slope*(x+y)]`。`--manufactured-pressure-slope` 默认0；设为1可检验非零壁面压力法向梯度，普通物理工况不允许此参数非零。解析体积加速度 `f=(U·∇)U+∇p−ν∇²U` 进入每个单元的动量右端，采用实际多边形质心的中点积分 `area*f`。除验证源项和静止壁面外，使用同一 SIMPLE、压力修正、对流/扩散及停止条件，不将解析速度预填为计算结果。

MMS 的 `.cells.csv` 额外导出积分源 `sourceX/Y` 和解析 `exactU/V/P`，其中解析压力也减去 cell 0 的解析值。普通工况没有这些列。独立 Python 检查从 CM2D 几何重新计算解析场与源，不使用导出 exact 列作为真值；动量失衡为面通量和减去体积源。记录面积加权速度 L2/Uref 和压力 L2/Uref²、最大误差及细化观测阶；字段一致性容差不是物理精度等级。

复现工具生成196/900/3,844格的规则方腔，再以 `δ=.06 sin(πx)sin(πy)`、`x'=x+δ, y'=y+.6δ` 连续扭曲内部点，边界保持单位方形。两种网格、两种对流格式、三档细化分别审核；扭曲网格仍须通过真实几何/拓扑检查。固定 ν=.1、Uref=1，停止容差1e-8，单例最长180秒、最多7,000轮；耗时仅是本机诊断，不是通用性能保证。

```sh
python3 tools/verification/run_manufactured_flow.py --output-root outputs/native-flow/manufactured-new --manufactured-pressure-slope 1 --scheme limited-linear
MPLCONFIGDIR=/tmp/cartmesh-flow-mpl python3 tools/visualization/render_manufactured_flow.py --summary outputs/native-flow/manufactured-new/runner-summary.json --output outputs/native-flow/manufactured-new/verification.png
```

普通运行要求新目录以保留旧证据；`--reuse` 只重新读回，写入独立的 `runner-reverification.json`，不覆盖原生成命令和返回码。`--levels 4` 可作快速冒烟检查，不能证明细化阶数；`--dry-run` 只标为计划，不能当成已运行。

该光滑制造解可验证非零壁面压力法向梯度，但避开移动顶盖角点，只能验证其覆盖的离散链路。它不能替代 Ghia 方腔、真实 Cut-cell 绕流、一般边界、湍流或网格无关性验证；原方腔细化失败仍保留。

### 独立方腔诊断

方腔正式审核的点值取样使用距离加权仿射拟合（v2）：单元与两个真实中心线壁值共同参与最近8点选择，拟合 `a+b*dx+c*dy`，归一化坐标和权重后求解。它可重现仿射场，秩不足明确失败；不宣称任意非线性数据有界或守恒。原IDW8采样作为 `benchmark.legacyIdw` 继续输出。Ghia参考表及各项数值阈值不变，混合采样方法的细化序列拒绝比较。原文 [Ghia 等（1982）](https://doi.org/10.1016/0021-9991(82)90058-4) 提供的是数值基准，并非解析真解。

可重审已有报告中的方腔场，不启动求解器：

```sh
python3 tools/verification/audit_cavity_sampling.py --summary outputs/native-flow/symmetric-stress/upwind-combined.json --summary outputs/native-flow/symmetric-stress/limited-linear-combined.json --max-iterations 4000 --output outputs/native-flow/cavity-sampling/reproducible-audit.json
```

该入口先核对源报告中网格及4份流场/日志文件的哈希，再独立读回守恒和中心线；`--max-iterations` 是原生产运行预算，非增加求解轮数。使用原验证器默认几何/连续性/物理阈值并记录它们；本例限制重构细化未通过，所以预期退出1、JSON `valid:false`。不要将该失败摘要覆盖为PASS。只给少于三档输入时不构成三档验证；源报告和大场文件仍在本地 `outputs/`。小型证据记录原运行与当前读回的哈希，未把旧场当成新求解。

`verify_cavity_vorticity.py` 另用流函数/涡量有限差分求解单位方腔 Re100，不调用生产 FVM。依赖 NumPy/SciPy；空间中心差分、Thom 壁面涡量、DST-I 泊松与 Heun RK2 推进到稳态。先校验离散正弦本征模，再以接受的新状态涡量方程残差和中心线变化同时停止；每档180秒，输出数组、中心线、实际命令及哈希。它是独立诊断，**不是新增产品求解器或认证真值**，也不改变原 Ghia 回归门。

```sh
python3 tools/verification/verify_cavity_vorticity.py --output-root outputs/native-flow/cavity-vorticity-new
python3 tools/verification/compare_cavity_centrelines.py --reference-root outputs/native-flow/cavity-vorticity-new --flow-summary outputs/native-flow/pressure-boundary/physical-final/upwind/summary.json --flow-summary outputs/native-flow/pressure-boundary/physical-final/limited-linear/summary.json --output outputs/native-flow/cavity-vorticity-new/comparison.json
```

比较器只读取已有结果，核对数组/CSV哈希、有限值和规则网格完整性；目前限定196/900/3,844格、单位方腔、速度1、ν=.01。FVM 在实际单元中心上做张量积线性插值，并加入真实壁面中心线值，与33/65/129节点参考在101个共同位置比较；原始 summary 的失败系列完整保留。坐标取12位仅用于识别行列，不改动求解网格或场值。FD 深处散度与壁面相邻散度分开报告，不能据此声称所有单元有限体积守恒。

### 方腔独立公开参考复核

`audit_cavity_reference.py`并排使用旧Ghia表与Marchi等2009表6/7；只针对单位方腔Re100。数据在`tools/verification/references/cavity-marchi-2009-re100.json`，包含原文DOI、PDF哈希、页码及估计误差。两组均使用相同的仿射取样实现和既有细化规则，保持旧失败可见，不把外推数值当成解析真解。

```sh
# 新目录：原生生成、串行求解、独立读回；默认三档196/900/3844格。
python3 tools/verification/audit_cavity_reference.py --generate outputs/cavity-reference-new --output outputs/cavity-reference-new/audit.json
# 加上 --levels 4 5 6 7 可扩展至15876格；默认单例预算360秒。
# 重审已有完整输出，不启动求解器：
python3 tools/verification/audit_cavity_reference.py --runs outputs/cavity-reference-new/runs.json --output outputs/cavity-reference-new/recheck.json
```

生成模式固定nu=.01、速度1、限制重构、对称应力、聚合预条件及容差1e-8；原生上限20000轮。既有目录拒绝覆盖，失败/超时保留日志与已产生的manifest。复审核对实际返回码、网格/二进制哈希、命令和摘要、单位几何、同格式同容差；至少三档并保持同一生产二进制，该二进制需仍在原路径可读取。不同版本请新生成一组，不混用旧场。参考估计误差用于呈现敏感性范围，不用来放宽原生容差。

当前预期是`independentAuditsPassed:true`、`publishedSequence.valid:true`、`ghiaSequence.valid:false`、`allChecksPassed:false`，退出1。退出码保留未通过检查，不是求解进程失败。两份参考的坐标集合不同，各序列的误差只在自身固定采样点上比较。新的公开参考说明旧表单调误差门的局限，但没有删除或重写旧门。

### 稀疏结构与压力预条件

流动系统的非对角项采用固定、按列排序的 CSR 连接表，重复连接共享同一项并相加；动量、压力与残差检查矩阵复用该结构，数组清零后重新组装。U/V/压力按顺序共享线性工作区，矩阵乘向量不再每步分配结果。自然单元编号不重排，物理单元 ID 不变。

压力默认 `--pressure-preconditioner ic0`，即自然顺序、零填充的不完全 LDLᵀ 分解；这是 IC(0) 的等价表示，只作用于预条件，不替换原方程。前代入、对角除法、后代入构成对称预条件。对称固定压力自由度，保留正对角；非对称系数、非正或非有限 pivot 明确拒绝，不移位、不修改容差、不暗中回退。`jacobi` 保留旧的对角预条件，便于隔离存储和预条件的性能贡献。两条路径仍使用原3000步上限和真实 `||b-Ax||₂ <= 1e-13 + 1e-11||b||₂` 停止条件。摘要及性能文件记录所选预条件器。

实现位于 `include/cartmesh2d/fv/detail/FlowLinearSystem2D.hpp`，本仓库自行编写，未复制或链接外部线性求解库。公开方法参考 [Netlib Templates](https://netlib.org/templates/templates.html) 的预条件共轭梯度及不完全分解、[PETSc PCICC 文档](https://petsc.org/release/manualpages/PC/PCICC/)的零填充/自然顺序概念。没有采用 PETSc 的默认移位策略，也不声称 ICC 是原创数学方法或足以取代多重网格。

`tests/flow_linear_test.cpp` 使用独立稠密矩阵检查真实残差、已知解、非对称 BiCGStab、图连接去重、压力约束、工作区复用及失败路径；CLI 同时测试两种预条件的物理解。单元测试不代替最终 Cut-cell 上的独立守恒与精度验证，实测结论见 CURRENT_STATE。

单一大规模通道的真实场、粗细连接放大及解析速度对照可用专用绘图入口（只读最终CM2D和CSV）：

```sh
MPLCONFIGDIR=/tmp/cartmesh-flow-mpl python3 tools/visualization/render_native_flow_scale.py --summary outputs/native-flow/sparse-scale/summary.json --label channel-500k --zoom 1.21 1.29 .46 .54 --output outputs/native-flow/sparse-scale/channel.png
```


### 扩展与性能测量

长期目标和当前进度只维护在 CURRENT_STATE。先完成稳态守恒/精度与压力求解效率，再依次扩展非定常、标量/热输运、SST 湍流、理想气体可压缩流；每种模式有独立物理配置和验收，不能用原稳态 laminar 的选项伪装支持新模型。

`cartmesh2d_flow_cli --profile` 在原六份结果之外增加 `.performance.json`。使用 monotonic `steady_clock` 测量读网格及构造时间、整个求解时间、动量/压力线性求解时间；统计线性求解调用、实际 Krylov 迭代总数和单次最大值。初始残差已满足条件的线性求解计0次迭代；每轮有两次动量和四次压力调用。线性耗时包括工作数组初始化与真残差检查，不含矩阵组装；总求解耗时还包含验证、组装、梯度、物理监测和进度回调，但不包含导出。此文件不测内存；需用系统工具另测峰值 RSS。到上限仍保存诊断且 `converged:false`，抛出数值错误的运行保留 stderr，不能当成完整性能样本。

```sh
build/cartmesh2d_flow_cli --mesh outputs/native-flow/formal-final/meshes/channel-l6/channel-l6.solver.cm2d --output outputs/flow-profile/channel --case channel --nu .01 --speed 1 --max-iterations 1500 --profile
# macOS 可加 /usr/bin/time -l；maximum resident set size 为 bytes。
# Linux /usr/bin/time -v 的 Maximum resident set size 单位为 KiB，不能直接混比。
```

相同几何/网格/工况/容差/线程数下比较，记录二进制哈希和系统版本。固定迭代吞吐不能冒充达到相同物理解精度的加速；完整求解须同时核对误差与守恒。初步小规模观测不外推50万格速度。前端只读取物理结果，不依赖计时文件；当前打包验证范围见 CURRENT_STATE。

理论与验证参考：[殷雅俊专著及简介](https://www.tup.tsinghua.edu.cn/booksCenter/book_07344201.html)、[NASA Turbulence Modeling Resource](https://www.nasa.gov/nasa-turbulence-modeling-resource/)。前者用于评估张量表述与推导，未作为已经证明的加速方法；后者用于后续具体湍流版本与验证设计，不表示已实现 SST。

## CFD 验证

既有高密翼型的高雷诺数试算使用 `tools/verification/run_airfoil_rans.py`：复制现有 `constant/polyMesh`，先跑标准与扩展 checkMesh，再运行 SST 稳态 RANS。它是固定约 1 m 弦长、六个命名 patch 的受限域试验工具；不重新生成网格，也不自动配置任意几何。

```sh
python3 tools/verification/run_airfoil_rans.py --source outputs/engineering-hybrid-final-high/naca/openfoam --output outputs/airfoil-rans-new
# 如残差尚未达到原停止条件，从最后写出的场继续，保持网格与容差：
python3 tools/verification/run_airfoil_rans.py --source outputs/engineering-hybrid-final-high/naca/openfoam --output outputs/airfoil-rans-new --resume --end-iteration 2500
# 本次正式试算的最后上限为 3500，未达停止条件仍如实记录，不无限续算：
python3 tools/verification/run_airfoil_rans.py --source outputs/engineering-hybrid-final-high/naca/openfoam --output outputs/airfoil-rans-new --resume --end-iteration 3500
MPLCONFIGDIR=/tmp/cartmesh-airfoil-mpl python3 tools/visualization/render_airfoil_rans.py outputs/airfoil-rans-new
```

要求本地 Docker 中已有 `opencfd/openfoam-run:2606`；绘图另需 NumPy、Matplotlib、Pillow。新算例目录必须不存在，续算需显式 `--resume`。入口 15 m/s、nu=1.5e-5、参考 Re=1e6、攻角 0°、湍强 1%、入口湍粘比初估 10；SST 配 Spalding/omega 壁函数。前 200 步用 upwind 初始化，之后 U 用 linearUpwind，k/omega 仍用 upwind，不能称全二阶格式。标准 checkMesh 未通过则不启动求解；扩展检查、进程完成、残差与力系数收敛分别记录。

这个已有翼型域上游和下游都只有约 0.5c，上下边界是 slip；所得升阻力不是自由远场验证。后处理直接按 OpenFOAM owner 对齐最终场与单元多边形，保留完整字段范围，记录首个压力方程的初始残差、力系数历史及输出坐标轴、通量平衡、出口回流及 yPlus，输出速度/Cp、流线和历史图，不自动授予质量 PASS。Cp 默认以出口表压零点为参考，图中明确标注；不是标准上游自由流静压参考。完整算例中的 `airfoil.foam` 可由 ParaView 打开；本轮实测数据见 CURRENT_STATE。

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

### 命名壁面载荷

`FlowResult2D::namedWallLoads`按custom边界名称积分无滑移静止/移动壁面；`FaceMomentum2D`中的实际压力及黏性通量进入力和力矩，入口/出口不计入。CLI摘要`namedWallLoads`保留压力/黏性分项、总量、面数和长度，`wallLoadReference=[0,0]`。力的单位m³/s²、力矩m⁴/s²，均除以密度和单位厚度，正力矩为逆时针；其他基准O处力矩为M_O=M_0-(O_x F_y-O_y F_x)。开放壁面分组的压力载荷随压力基准改变，不把单组载荷冒称全域守恒。

`verify_native_flow.py::audit_named_wall_loads`从真实CM2D面中心/外法向与CSV重新积分，原方程审核另行重建各面通量；非定常共享同一审核。`desktop/src/core/wall-loads.js`检查名称覆盖、单位基准语义、分项加和与全壁面合力，App显示各组，旧无载荷字段的结果兼容但不补造数据。

### 平滑移动壁面与旋转环隙

`MovingWall`保留逐面常量速度迹线；`SmoothMovingWall`/`smooth-moving-wall`把面心速度解释为平滑壁面速度场的采样，黏性应力重构保留单元梯度的切向分量，再施加紧致法向导数。二者均要求实际面上的法向速度为零，固定网格，没有刚体网格运动或滑移网格。类型写入边界v1和检查点v4（扩展枚举，旧程序遇未知类型失败）；改变类型不能继续原检查点。避免在速度有跳变的尖角处误用平滑假设。

`rotatingAnnulusBoundaryPreset2D`检查两圈同心偶数正多边形（至少16边），使用最终实际分段面生成内壁逆时针速度、静止外壁。`annulus`仅为边界模板名，不是新增固定求解场景，真正求解仍用custom。App同心圆样例与原方形带孔annulus样例分开。

```sh
build/cartmesh2d_flow_cli --mesh case.solver.cm2d --case annulus --speed .5 --export-boundaries case.boundaries
build/cartmesh2d_flow_cli --mesh case.solver.cm2d --case custom --boundary case.boundaries --nu .1 --speed .5 --convection face-limited-linear --pressure-preconditioner aggregation --output outputs/ring/flow
python3 tools/verification/verify_rotating_annulus.py --output outputs/ring-reference-new
```

参考工具独立构造环形四边形以分离求解器误差，不称为Cut-cell；仍由原生几何工厂执行完整Solver质量门，Python另重算几何、通量与方程。径向/角向和多边形圆边界同步细化；解析式与力矩来自[UT Austin](https://farside.ph.utexas.edu/teaching/336L/Fluidhtml/node137.html)，速度L2按内壁速度归一化，压力L2按速度平方归一化，压力基准取cell0。预设误差门在运行前保存，加严容差对照单独核查迭代误差；不把圆环结果泛化为全部旋转设备或三维Taylor涡。


### 保存流动工况与分段续算核对

macOS 0.4.21 的“保存工况 / 读取工况”最初使用 `cartmesh2d-flow-case-v1` JSON（`desktop/src/core/flow-case.js`）。绑定最终CM2D文件SHA256及单元/面/顶点数，保存规范化的物性、边界逐面数值与名称、数值格式、时间控制和初始局部涡。读取还逐面检查边界位置、法向、owner、覆盖及物理约束；未知版本、额外/缺失/无效字段明确拒绝。JSON不携带可执行命令或读写路径。保存通过同目录临时文件完成后替换目标。

这是从零起算的流动设置，不包含几何生成参数、网格数据、温度条件或已接受流场。下一次会话先按原参数生成同一最终网格，再读取；完整项目导入仍待做。续算模式下保存按钮禁用；读取会取消流动/温度续算并清除旧场的界面绑定，必须重新计算。已有结果包仍保留上次完成的物理结果，不会被配置文件重新标记。真实App验证入口 `--flow-case-check=true`（同时提供 `--smoke` 和绝对 `--out`）替代文件选择对话框，经过相同IPC/表单处理，修改参数后读回，并检查两次真实计算的场及时间历史一致。

`tools/verification/verify_flow_trajectory.py`处理已完成的`cartmesh2d-bounded-trajectory-v1`外流序列。逐段复算最终面动量与守恒，完整读取v2检查点几何、拓扑、物性、U/V/P/FLUX，并与CM2D/CSV绑定；后一段输入字节必须等于前段接受状态。压缩目录须通过整体及逐文件SHA256校验，解包拒绝越界路径、链接和重复文件。父研究的最终时间/检查点/网格/二进制必须接上当前研究；中间时间的监测行不冒称全部重建了中间流场。

```sh
python3 tools/verification/verify_flow_trajectory.py --trajectory outputs/native-flow/long-wake/trajectory.json --initial outputs/native-flow/initial-vortex-wake/flow --output outputs/native-flow/long-wake/verified
MPLBACKEND=Agg python3 tools/visualization/render_flow_trajectory.py --verified outputs/native-flow/long-wake/verified.json --output artifacts/current/native-long-wake
```

原始目录已压缩时从同名`.tar.gz`和`.manifest.json`独立验证；工具不重算求解器、覆盖失败输入或延长预设计算预算。当前绘图工况固定为Re100圆柱，显示真实网格、速度/重构涡量及完整受力历史，不能作为空间/时间独立性或饱和脱涡资格。


0.4.22把App流动容差接到原生`--tolerance`，范围1e-12..1e-6，默认1e-6；核对结果容差等于实际请求。流动工况v2增加`tolerance`，读v1时显式赋旧默认1e-6，不允许v1额外携带该字段。连续性和全局相对流量门仍为1e-8。共享界面的温度联算使用`min(tolerance,1e-8)`作为`--flow-tolerance`，原来1e-8默认保持不变；温度自身停止门未修改。检查点保留物理状态而不锁死迭代预算，续算允许加严容差。实际App回归可用`--flow-tolerance=1e-9`，工况保存回归会改动并恢复该值。

`tools/visualization/render_annulus_facets.py`绘制保存在`outputs/native-flow/annulus-facet-study/`的固定背景、圆边界细分研究，重新审核真实场并核对解析误差。完整单元压力误差用对称对数色标显示，最大偏差不能被颜色裁切；两侧力矩的非单调变化一并保留，不能把该研究误称三档流动网格验证。

### 外流参考适用性与有界尺寸重试

`verify_native_flow.circular_obstacle_reference`在最终嵌入边界上检查闭合单环、圆半径和等角折线，剔除共线分割点；只有匹配圆柱才应用圆柱对称升力检查/参考对照。其他几何保持not-qualified，仍须通过同一离散守恒/动量等检查；不把翼型当成圆柱校验。对应回归在`tests/flow_verifier_test.py`。

`cell-budget.refineFailedQuality`只在尚无有效候选且原生明确Solver质量拒绝时细化一个树层级；`budget-runner`仍最多三次，保留安全深度/估计数量上限、所有失败原因和最接近数量目标的有效结果。main保存每个失败目录的`generation-failure.json`，全部失败时再写`selection-failure.json`；不因重试而删除坏网格证据。原质量门不变。实际App可用`--smoke=thick_airfoil --target-cells=2000 --flow=external --flow-nu=.05 --flow-tolerance=1e-8 --flow-convection=face-limited-linear --flow-pressure-preconditioner=aggregation`重现三档流程。

`verify_flow_trajectory.py --completed-prefix`仅允许复核显式`complete:false`且已记录停止原因的连续有效前缀；拒绝仍在进行中的研究，也不把尾部失败尝试或未达到的目标时间标记为完成。默认仍要求完整成功序列。绘图始终使用最后一个已审核的最终场。

`tools/verification/verify_wake_time_refinement.py`从已审核轨迹中选择同一个完整接受状态，按三档等比固定步数推进相同时间窗，绑定可执行文件/网格/原始检查点哈希并分别复算最终场与历史。比较相邻速度和压力场差及其观测阶，`accuracyQualification`保持`not-qualified`。原始状态积累的误差、空间误差和长期频率仍未排除；工具不会自行增加预算或重写失败案例。


`tools/verification/verify_external_spatial.py`保持同一32段单位圆和10Lref外域，把壁面及背景尺寸同时减半做三档稳态Re20对照，并在最细网格加严100倍容差。所有候选须完成原生Solver质量门及独立动量/通量审核；失败保留，不追加计算预算。`render_external_spatial.py`重算最终场审核后绘制实际网格、速度及积分载荷变化。固定折线几何与有限外域误差仍在，三档结果不自动赋予连续圆柱精度资格。

外流独立读取器把力与已无量纲化的Cd/Cl分开处理：力按`0.5 Uref² D`归一化，已有系数不再次归一化，混合字段先还原为同一量纲；同时提供力与系数却相互矛盾时明确失败。回归使用非单位速度防止U=1、D=2掩盖单位错误；不改变原生求解出的力。

### 静压驱动法向开口（0.4.24）

命名边界新增`pressure-opening`。只接受轴对齐边，`u=v=0`是配置占位而非零法向速度；`p`为静压除以密度，单位m²/s²。法向速度零梯度，压力修正决定守恒面通量；流入时切向速度为零，流出时速度零梯度。它与继续拒绝回流的`pressure-outlet`分开。方向倾斜时明确失败，不悄悄投影。完全静止且净面通量严格为零时，零吞吐的全局相对不平衡按零报告；非零吞吐、单元连续性和原方程停止门不变。

App选择“命名边界→压差驱动通道：两端静压→从预设生成”，初始左p/rho=1、右=0，可编辑两组压力。两组压力均整体增加常数不会改变速度。参考速度仅归一化，不规定流量。边界文件、工况v2及检查点v4保存新类型；旧程序会拒绝未知类型。命名边界尚未支持温度配置。

`tools/verification/verify_pressure_openings.py --output outputs/native-flow/<fresh-directory>`生成独立矩形网格并审核压差通道，默认256/1,024/4,096格、tol1e-10、最多10,000轮；每命令180秒、总900秒、50MiB、至少5GiB剩余空间，所有失败保留。`tools/visualization/render_pressure_openings.py <directory> <png>`先校验文件哈希与细档离散方程再作图。公式L=4、H=1、nu=.1、Δ(p/rho)=4.8给出u=6y(1−y)、Q=1；这是特定Poiseuille参考，不是曲壁网格或通用精度认证。运行`ctest --test-dir build -R 'flow_boundary|pressure_opening' --output-on-failure`可复现三种格式、等压静止、坐标/压力变换、普通出口拒绝和续算检查。

### 命名对称边界（0.4.25）

`custom`新增`symmetry`，配置u/v/p均为0占位；仅水平/竖直边，法向速度为零、切向速度零法向梯度，不提供无滑移摩擦。它不会计入`namedWallLoads`，但其压力仍参与共享面动量守恒。斜面或非零配置值明确失败。边界文件、工况v2及检查点v4保存该类型；旧程序拒绝未知类型。

App复用`examples/acceptance/rectangle.xy`作为2×1m矩形内流示例，选择“命名边界→压差驱动半通道：顶部对称”。该模板先通过原生矩形通道几何门，再将上下壁中的顶部独立分组为对称面、底部保持无滑移。用户可以编辑任意已有命名组为轴对齐对称类型，不能用倾斜曲壁近似对称平面。

`tools/verification/verify_symmetry_flow.py --output outputs/native-flow/<fresh-directory>`在独立64/256/1,024格矩形上验证L=4、H=1、nu=.1、压差1.2的u=1.5y(2−y)、Q=1；另将两侧都设为自由滑移，从静止用dt=.02推进三步，核对u=.3t以及1+2步与3步检查点逐字节一致。每命令180秒、全研究600秒、30MiB、磁盘至少5GiB；测试三种格式入口为`cartmesh2d_symmetry_flow_cli`。作图复用`render_pressure_openings.py`，自动识别半通道参考，并先重新审核源文件与离散方程。

### 命名边界流量（0.4.26）

原生流动CLI在custom摘要添加`namedBoundaryFluxes`和`boundaryFluxDefinition`：每组`faces`、实际边长`length`、非负`inflow`/`outflow`、带符号`net`、`normalMeanVelocity=net/length`。正方向为流体域向外；流量单位是m²/s每单位厚度，法向平均速度m/s。它们直接积分求解器修正后的面通量，局部流入与流出分别保留，不从单元中心速度近似或把净量当作两者之和。壁面/对称面严格不穿透。密度未参与，所以不是质量流量。

独立读取器从最终CM2D几何和faces.csv重算每组数据；前端核对分组、量纲声明、正负号、平均速度和全局汇总。旧记录可以缺少整套字段，审核明确记为unavailable；只保留半套字段或篡改数值会失败。真实新程序均输出并在App结果区域显示，随现有JSON/ZIP一起导出。原求解方程不变，2,784格喷管与0.4.25冻结程序的六种实际场/历史文件逐字节一致。

真实App验收可追加`--flow-require-converged=true`；出现未收敛诊断场时返回非零退出码，避免只凭场被显示便认定计算完成。`--flow-flux-shot=true`检查所有命名组的流量已经进入实际结果界面并定位截图。默认smoke仍可用于验证故意未收敛、取消和失败显示。

### 稳态松弛与面通量一致性（0.4.27）

稳态预测面通量现加入`(1-alphaU)*(phi_old-interpolate(U_old)·S)`，内部插值保留同样的偏斜修正，定压边界使用同样的旧单元响应。这样松弛后的单元动量与面方程在收敛时消去alphaU；此前仅非定常路径含此项，稳态离散解会随松弛系数变化。非定常旧时间项与原运算路径保持不变。该问题的理论背景参见[Majumdar 1988](https://www.tandfonline.com/doi/abs/10.1080/10407788808913607)，实现从现有离散方程推导，不引入外部求解核心。

`flow_boundary`新增108格变形半通道回归，alphaU=.6/.8/.9、tol1e-11。旧代码速度/压力差约1.67e-4/9.85e-4，修正后最大速度/压力差约3.34e-10/5.23e-12。真实1,300格Cut-cell半通道以.6/.8对照：tol1e-9时最大速度差3.99e-7、压力差1.11e-8；tol1e-10时分别3.93e-8、1.15e-9，剩余差随停止容差降低。两档均独立审核原动量与面通量，未修改停止门。

CLI默认仍为.6，显式速度松弛选项仍仅允许非定常；诊断通过临时API入口使用相同库测试，不将参数扫描入口打包。旧核心在.9时120秒超时、.95/1.0数值范围失败，不能把调大参数当作通用提速。修正后默认仍约16,525轮，性能问题未解决。稳态摘要记录`steadyFaceInterpolation=iteration-flux-defect-skew-corrected-v1`和实际`velocityRelaxation`，历史结果仍绑定各自二进制证据。

### 有限历史稳态加速（0.4.28）

CLI使用`--steady-acceleration anderson`，省略或`none`保持普通SIMPLE；App流动设置加入同名可选项，默认关闭。它不改变时间步或物理方程，不支持非定常和材料更新/RANS路径；温度联算不使用。工况文件升级v3保存该数值选项；v1/v2读取时明确补none，拒绝在旧格式中夹带新参数。结果记录候选/接受/拒绝计数，App与本次请求核对，最后一次迭代必须是普通SIMPLE收敛确认。

实现位于`detail/Anderson2D.hpp`及`Incompressible2D.cpp`：以体积/边长和物理参考量缩放U/V/P/共享面通量，最多四个历史差分，二次正交QR求最小二乘系数。秩不足或系数绝对和超过1e4不提出候选；候选必须降低原非线性动量残差且满足原1e-8局部/全局质量门，否则逐字节恢复原SIMPLE场并清除历史。数值异常的候选也计入拒绝，不作为最终场。参考[Walker与Ni 2011](https://users.wpi.edu/~walker/Papers/Walker-Ni,SINUM,V49,1715-1735.pdf)的Type-II固定点加速，接受门直接检查本项目的离散方程。

同时将已用于非定常的动量线性每行`abs(b-Au)/aP <= .01*tolerance*Uref*alphaU`约束用于稳态，既有全局线性残差门保留。否则严格稳态停止门可能仍未达到，而线性求解已经返回零更新。保留64格half-channel、upwind、nu=.1、tol1e-11的5000轮停滞，修正后990轮收敛。

`cartmesh2d_steady_acceleration_cli`独立审核三种对流的原动量/连续性、与普通SIMPLE的场差及不支持的选项；`flow_boundary`核对变形半通道、旋转、压力平移与封闭顶盖腔。当前真实对照与原始命令在`outputs/native-flow/anderson-suite-current/`及同名完成归档；`render_steady_acceleration.py --study <directory> --output <png>`先核对各文件哈希和半通道方程再绘图。不要把单次本机运行时间当成其他工况的保证。App验收可加`--flow-steady-acceleration=anderson --flow-acceleration-shot=true --flow-require-converged=true --flow-case-check=true`，实际检查开关、加速计数、工况修改后重载及重复场。
