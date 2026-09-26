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

`cartmesh2d_fv_cli --help`是独立扩散/泊松入口。`cartmesh2d_euler_cli --help`是实验理想气体Euler入口；SST输运/壁距/稳态RANS代码位于`src/fv/Sst*`及对应测试。它们仍有调用与测试，保留现有行为；本轮不扩大物理范围，也不把局部通过当作通用湍流或可压流资格。

## 流体拓扑优化研究入口

`tools/optimization/` 是用户授权的独立研究原型，当前为命令行入口，未接入桌面菜单，也未替换原生不可压求解器。设计变量覆盖内部材料分布，允许流道连接关系变化；不是只调整几根管道的宽度或控制点。它不读取产品 `.background.json` 作为流体网格。

### 模型、灵敏度与停止条件

`brinkman.py` 独立实现原生二维交错 MAC 分析网格，求解无量纲 Stokes–Brinkman 方程 `-mu Laplacian(u) + grad(p) + alpha(rho) u = 0`、`div(u)=0`。`rho=1` 表示流体，`rho=0` 为有限阻力近似固体；`alpha=alpha_max*q*(1-rho)/(q+rho)`。它不含对流惯性，不能当作通用 Navier–Stokes 拓扑优化器。

入口/出口采用积分匹配的抛物线速度，两格被动端口区域固定为流体，外壁单元固定为固体；设计区由内部变量决定。锥形密度滤波与 tanh 投影后计算真实约束量 `mean(rho) <= volume_fraction`。`filter_radius` 使用设计域长度单位，不等于已证明的制造最小壁厚。

稀疏离散系统消去规定速度及一个压力参考自由度；伴随使用同一离散矩阵的转置，并通过滤波、投影求导。默认目标是离散黏性耗散加 Brinkman 阻力耗散；可选 `pressure-power` 是端口邻接压力单元的通量加权压差功率，不能误称精确边界应力功。OC 更新同时满足变量界、移动界与物理体积约束，候选须重新求解并通过目标回溯才接受。

三阶段 `(q,beta)=(.01,0),(.1,2),(.1,6)` 改变了模型，不能把跨阶段目标变化当成一次性能提升。报告中的两组参考设计均重新投影到相同体积、使用最终相同参数求解。`uniformPorous` 是数值初值，`geometricSeed` 是可复现的几何种子；两者都不是独立工程基准。

`projectedKkt` 为目标梯度归一化后、体积约束切平面及变量界上的投影步无穷范数；默认 `1e-3` 是局部设计停止条件，不是 CFD 精度。达到迭代/时间预算、回溯停滞或分析失败会明确报告，不能写成优化收敛。直接线性解及伴随的 `1e-8` 门检查 `||Ax-b||inf/||b||inf`，连续性检查最大单元净通量/总入口流量，用于拒绝不可信梯度；不是物理误差要求。失败或中断保存最后已接受设计及对应的阶段参数；当前没有自动续算入口。

### 运行与实际网格连接

可选 Python 依赖仅供研究脚本使用，不进入原生产品构建：

```sh
python3 -m venv outputs/topology-env
outputs/topology-env/bin/python -m pip install -r tools/optimization/requirements.txt
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 outputs/topology-env/bin/python tests/flow_topology_test.py
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 outputs/topology-env/bin/python tools/optimization/optimize_flow.py --case double-pipe --output outputs/topology/double-pipe
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 outputs/topology-env/bin/python tools/optimization/optimize_flow.py --case bend --volume 0.3 --output outputs/topology/bend
```

输出目录必须是新目录。默认 48×32 分析格、25/40/60 步预算；首轮 Matplotlib 字体缓存可能需要额外时间。`summary.json` 保存控制、版本、源代码哈希、停止原因和参考对照；`history.csv`、`accepted-design.npz`、阶段快照、`final.analysis.vtk` 保留数值证据。VTK 明确是含固体阻力的分析场，不是产品求解网格。`--no-plot` 同时跳过等值轮廓提取，可随后单独运行 `topology_artifacts.py <结果目录>`。

`topology_artifacts.py` 对投影密度提取分段直线等值轮廓，保留所有连通区域与孔洞，不平滑折线、不填孔。`fluid.xy` 显式使用 interior 语义；每个独立流域另有 `fluid-component-N.xy`，坐标及哈希可追溯。阈值轮廓面积与分析密度体积是不同量，分别报告。

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/clang++
cmake --build build --target cartmesh2d_cli cartmesh2d_flow_cli -j2
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 outputs/topology-env/bin/python tools/optimization/verify_extracted_flow.py outputs/topology/double-pipe --output outputs/topology/double-pipe-native --levels 5 6 --tolerance 1e-8
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 outputs/topology-env/bin/python tools/optimization/verify_extracted_flow.py outputs/topology/bend --output outputs/topology/bend-native --levels 6 7 --tolerance 1e-8
outputs/topology-env/bin/python tools/optimization/render_native_flow.py outputs/topology/double-pipe-native outputs/topology/bend-native --labels Double-pipe Bend --output outputs/topology/native-preview.png
```

连接脚本真实执行原生 Cut-cell 生成、Solver 质量门、独立 CM2D 读取/面积核对、原生低 Reynolds 数流动及独立离散方程审计。多个不连通流域逐个求解以提供独立压力参考，全部保留并汇总面积；任何一个失败均不能报告整例成功。失败分辨率和原始日志仍保留，不能挑图隐去失败。检测到本机 `checkMesh` 时实际运行标准检查，否则明确 `not-run`。

上述流动显式采用速度尺度 `.02`、运动黏度 `1`，按端口宽度算名义 Reynolds 数约 `.00333`；入口保持同形抛物线，出口改为原生压力出口。`1e-8` 是这两例为了通过既有独立方程审计所用的原生代数停止控制，未修改产品默认值或审计门。这一步证明提取几何能被原生网格/求解链接受；固体阻力、出口条件及离散格式均改变，不能把多孔分析的目标下降直接写成真实壁面 CFD 性能提升。`render_native_flow.py` 直接绘制接受网格多边形与原生 CSV，保留全部区域并记录源哈希。

### 公平的真实壁面对照

`compare_sharp_designs.py` 对比最终候选与 `reference-geometricSeed.npz`。等密度体积不保证等提取面积，因此默认先量取候选的实际面积，核对不超预算，再仅对基准内部密度施加标量偏移，通过二分匹配轮廓面积；候选原轮廓及两边被动端口区域均不改。这是明确记录的设计提取操作，不是修图、平滑用户 XY 或删格过门。`--area-target budget` 可显式选择把两边都投影到面积上限，产生的形状必须重新接受原生检查。

每一档分辨率都对两边使用相同入口积分流量、压力出口、物性及数值控制，计算 `sum(p_in Q_in)-sum(p_out Q_out)` 与其除以入口流量得到的运动学压降。这里的 `p` 是 `p/rho`，压降单位 m²/s²；功率量按密度及单位厚度归一化，不能直接写成瓦数或商业泵效率。基准与候选都要通过网格、原生收敛、独立方程审计；失败及不等流量不能变成有效排名。

```sh
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 outputs/topology-env/bin/python tools/optimization/optimize_flow.py --case bend --nx 72 --ny 48 --volume 0.3 --alpha-max 1000000 --iterations 50 100 250 --no-plot --output outputs/topology/bend-high-resistance
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 outputs/topology-env/bin/python tools/optimization/compare_sharp_designs.py outputs/topology/bend-high-resistance --output outputs/topology/bend-sharp-comparison --levels 6 7 --small-alpha 0.25
outputs/topology-env/bin/python tools/optimization/render_native_flow.py outputs/topology/bend-sharp-comparison/baseline-level-7 outputs/topology/bend-sharp-comparison/candidate-level-7 --labels Baseline Optimized --shared-scales --output outputs/topology/bend-sharp-comparison.png
```

`--small-alpha` 是原生已有的小单元守恒聚合触发比例，默认 .1；.25 不改变质量门、物理面积或边界折线，案例两边必须相同。它改变离散网格，因此结果和计时须绑定该参数。原 .1 配置有粗候选质量拒绝，记录保留；不能把不同配置最有利的几档拼成“同设置网格收敛”。工具默认会逐档计算完整成对结果，不复用不明来源的旧解。

`assessment` 分别保存每档改善、两方案最后两档压降变化，以及“细网格压降差减去上述变化之和”。只有全部指定档位成对有效、最后两档排名一致且该差为正，才置 `meshRobustImprovementObserved=true`。这只是保守的**已观察敏感性**标记，不是误差估计器、概率置信度或网格无关性证明。当前两档改善 2.81% / 2.39%，流道连通数始终为一，不能声称已实证发现新的分叉拓扑。

12 项相关测试覆盖解析 Poiseuille 单网格基本检查、两种目标/投影的伴随有限差分、体积导数、下降及约束、确定性、孔洞/多区域提取、非法输入、预算与失败状态保存，以及真实面积匹配、端口保护和排名判断。Poiseuille 的 2% 相对误差界仅针对 32×24 的低成本开发检查；成对网格实验也不是跨平台或通用物理精度验收。

方法背景见 [Stokes 拓扑优化示例](https://www.dolfin-adjoint.org/en/stable/documentation/stokes-topology/stokes-topology.html)。[已有 Cut-cell 流体拓扑优化研究](https://doi.org/10.1016/j.camwa.2021.06.002) 说明功能组合本身不构成原创性。后续可研究“多孔优化的候选排序能否在真实壁面 Cut-cell CFD 中保持，以及如何控制排序误差”；当前尚未完成该科研结论。

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
