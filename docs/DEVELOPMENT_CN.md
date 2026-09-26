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

## 纯笛卡尔浸入边界研究入口

`codex/pure-cartesian-cfd` 的 `cartmesh2d_immersed_cli` 是独立原生 C++ 开发原型，尚未接入 App。核心在 `include/cartmesh2d/immersed/CartesianFlow2D.hpp` 和 `src/immersed/CartesianFlow2D.cpp`，复用二维几何诊断和既有稀疏线性代数。它不读取或伪装 `*.solver.cm2d`，不修改背景 JSON 的 `solver_ready=false`，也不经过或降低现有 Cut-cell 的 Solver 质量门。

### 方程和当前适用范围

计算域为 `[0,L] × [0,H]`，均匀 MAC 交错网格：压力在格心，u/v 在对应面中心。所有方格完整保留，固体中的数值是辅助未知量，不属于真实流体。x 方向速度和压力扰动周期，y 两壁静止无滑移。恒定 x 加速度 `drive` 等价于周期压差驱动；输出运动学总压力可写为 `p_total = p_fluctuation - drive*x`。这是周期通道，不是入口/压力出口外流。

离散不可压 Navier–Stokes–Brinkman 方程：`du/dt + div(uu) = -grad(p) + nu Laplacian(u) + drive e_x - chi u/eta`，`div(u)=0`。对流为一阶守恒迎风，黏性项为中心差分；对流/扩散显式、阻力隐式，自动限制伪时间步。压力校正使用 `beta=1/(1+dt chi/eta)`，求解 `-div(beta grad(phi)) = -div(u*)/dt`，然后校正速度和压力。压力参考设为第一个格心，接受前仍检查该格的连续性。默认 PCG/IC0；macOS 可显式选择系统 Cholesky，CLI 在调用前固定本进程的 `VECLIB_MAXIMUM_THREADS=1`。当前用伪时间推进求稳态，尚未取得非定常精度资格。

`chi` 来自原始二维折线的有符号距离，默认在半宽 `0.5 min(dx,dy)` 内使用正弦过渡。折线本身不移动、不平滑、不裁切网格。`eta` 是阻力时间尺度，默认 `1e-4` 秒；有限阻力和掩膜过渡都会产生壁面误差。真实壁面上按原线段插值采样速度，分别报告总速度、法向穿透和切向滑移。`penalty_drag_per_density` 是辅助阻力体积分，不是已验证的表面应力阻力系数。壁面诊断包含每条原线段两端及不超过 `h/8` 的间隔，`h=min(dx,dy)`；报告的是这些稠密采样中的最大值。

`channel` 不含内置固体。`cylinder` 是圆心 `(L/3,H/2)`、半径 `.15 H` 的 128 段规则多边形，实际坐标保存在 `boundary.xy`；流动代表周期通道里的重复障碍物。`custom --boundary SOLID.xy` 接受以空行分隔的原生多环固体，嵌套按奇偶规则保留孔洞；当前只支持静止壁面，不接受带命名边界角色的 XY 元数据。非法、自交、重复边、接触环显式拒绝。固体离外域边界/周期接缝须留出两格加掩膜半宽，至少有一个完整内部单元，否则报告当前不支持或欠分辨；这不是任意细缝/薄壁已被充分解析的证明。

### 压力与真实壁面力耦合

显式启用 `--wall-method surface-penalty` 时，在原体积 Brinkman 项上增加整段壁面惩罚。`J` 是交错速度在原始折线上的双线性插值，`W=ds*h/(dx*dy)`，`E=sqrt(W) J`。每条线段在两种 MAC 分量的插值结点线处划分积分区间，每区间使用三点 Gauss 积分；区间内 `J` 沿线段至多二次，因此 `Jᵀ W J` 的四次积分是精确的（浮点算术范围内）。这些是积分点，原输入折线和完整方格均不改变。新增/倒置共线顶点的回归检查几何分段不会实质改变解。

令 `D` 为 MAC 散度，`B=[-h D; E]`。每个伪时间步同时求解

```text
(B beta Bᵀ + diag(0, eta_wall/dt)) q = B u_star/dt
u_new = u_star - dt beta Bᵀ q
p_new = p_old + h q_pressure
wall_force = -Eᵀ q_wall
E u_new = eta_wall q_wall
```

压力参考行仍被固定，但接受前检查全部格的散度。壁面力、压力和体积阻力使用同一个 `beta`；没有求解后截断壁速的步骤。`eta_wall > 0` 保证壁面块有正的顺应项，处理重复约束，不通过删除坏约束或对角补丁凑成可解矩阵。模型的壁面功率满足 `sum(u*wall_force)*dx*dy = -eta_wall sum(q_wall²)*dx*dy <= 0`，独立读取器重新构造 `Eᵀ q_wall` 并核对这个恒等式。

`--wall-penalty-time` 默认设为 `1e-4` 时间单位；本次比较显式用 `1e-6`，无量纲值 `eta_wall Uref/H=1e-6`。它是当前尺度的开发参数，不是通用最优值或物理滑移长度。该模式依然是**有限表面惩罚**，不是精确的无滑移/锐界面资格。单纯减小体积 `eta` 不能消除插值壁面误差；首个中点壁面版本还出现“标记点很小、点间穿透较大”的斜壁案例，因此改用整段积分，保留原结果。

`--linear-solver auto` 在旧 Brinkman 模式选择 IC0，在新耦合模式选择 Jacobi；macOS 可显式选 Cholesky。原 IC0 对这种含壁面块的 SPD 矩阵出现过非正主元，现明确拒绝此组合，失败复现在 `outputs/wall-treatment/gauss-ic0-probe/`，不静默切换或改变线性门。Jacobi 在本机完成小规模稳态和三步场对照；大规模效率、其他平台未验证。资源上限为 16,384 个壁面积分点、8,000,000 个原始图连接贡献，超限显式失败，不删除输入几何。

`wall_force` 是加速度，`wall-markers.csv` 的 `multiplier_*` 为 `q_wall`；`weight` 是无量纲 `W`。新增表面/体积反力应同时读取：二者可有反向分量，不能单独把其中一项当作真实阻力系数。`surface_power_per_density` 是完整辅助计算域内的模型功率；`force_balance` 是驱动力、两项障碍物反力和上下通道壁黏性力之差除以 `drive*L*H`。`wall_normal_flux_net/abs` 是原折线稠密采样的法向有符号/绝对速度线积分，量纲为长度²/时间。网格散度很小并不使这些物面通量自动为零。

### 停止量和结果语义

`Uref = drive H²/(12 nu)` 是无障碍通道解析平均速度。开发默认 `steady-tolerance=1e-4` 检查最大离散动量余量除以驱动加速度，同时检查最大速度步变化除以 `Uref`。这是相对驱动的 0.01% 代数平衡目标，不是壁面或物理误差要求。`continuity-tolerance=1e-6` 检查 `max|div(u)| H/Uref`，每个接受步都必须满足。新模式还逐步检查 `max|Ju - eta_wall q_wall/sqrt(W)|/Uref <= continuity-tolerance`，这是与速度同尺度的壁面方程求解误差；它与真实壁速、法向穿透分开报告，不作为物理无滑移证明。该连续性和截面通量覆盖含辅助固体的完整计算域，不能代替真实壁面不穿透验证。默认线性相对 L2 目标 `1e-8` 用来使压力校正的误差小于这些开发停止量，仍使用原线性工具的 `1e-13` 绝对算术余量并独立核对实际残差；不将其中任何数值当作通用精度标准。

达到步数预算返回码 2 和 `iteration-limit`。失败候选不替换最后接受步；若一个步都没接受，结果明确为 `initial-only`，不能作为检查点。POSIX 取消保存最后接受步，返回 130。线性/连续性失败返回 `candidate-failed`；输入/导出错误返回 1。`steady-converged` 返回 0，仍不表示网格无关或有限阻力壁面精度合格。首个 4096 格冷启动通道在 20000 步时未达到动量目标，原结果保留在 `outputs/immersed-prototype/channel/`；默认预算据此设为 40000 步，没有放宽残差。

### 构建、运行和独立核查

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/clang++
cmake --build build --target cartmesh2d_immersed_cli -j2
ctest --test-dir build -R '^cartmesh2d_immersed_flow$' --output-on-failure
build/cartmesh2d_immersed_cli --case channel --nx 128 --ny 32 --nu 0.01 --drive 0.12 --max-steps 40000 --output outputs/immersed-channel
build/cartmesh2d_immersed_cli --case cylinder --nx 128 --ny 32 --nu 0.01 --drive 0.12 --max-steps 40000 --output outputs/immersed-cylinder
# macOS 可选压力后端，使用新目录保留 IC0 对照
build/cartmesh2d_immersed_cli --case cylinder --nx 128 --ny 32 --nu 0.01 --drive 0.12 --max-steps 40000 --linear-solver cholesky --output outputs/immersed-cylinder-cholesky
# 新壁面模式显式比较参数；其他平台省略 Cholesky 参数，auto 选择 Jacobi
build/cartmesh2d_immersed_cli --case cylinder --nx 128 --ny 32 --wall-method surface-penalty --wall-penalty-time 1e-6 --linear-solver cholesky --output outputs/immersed-wall
python3 tools/verification/verify_immersed_flow.py outputs/immersed-wall --require-converged
python3 tools/verification/verify_immersed_flow.py outputs/immersed-channel --require-converged
python3 tools/verification/verify_immersed_flow.py outputs/immersed-cylinder --require-converged
python3 tools/visualization/render_immersed_flow.py outputs/immersed-cylinder --output outputs/immersed-cylinder/preview.png
```

绘图脚本需要 NumPy/Matplotlib；原生求解与独立读取器不需要它们。输出目录必须尚不存在，避免覆盖之前的失败/接受证据。`summary.json` 记录实际参数、停止原因、量纲/归一化和网格/边界/求解/导出分段时间（`pressure_seconds` 为包含在求解总耗时中的压力线性求解时间）；`u.csv`、`v.csv` 保留所有交错自由度、掩膜和表面力，`cells.csv` 和真实二维 `field.vtk` 保留每个完整方格及几何分类。`wall-markers.csv` 保存积分点/权重/乘子，`walls.csv` 是独立于积分点的稠密实际折线采样，`history.csv` 只记录接受步，原始折线单独导出。不能将仅网格生成的时间视为总成本。

独立 Python 读取器重新计算动量、连续性、通量、壁面插值、原折线面积和掩膜样本，并核对 VTK/CSV 一致性及线性真残差。18 项相关测试覆盖解析通道、圆柱流动阻滞/对称、原始多边形/孔洞、确定性、预算/失败/取消、非法输入/损坏导出、整段壁面改善、共线分段/方向不变性、表面力破坏检测、空域不变性和两种耦合后端场对照。斜壁改善至少 4 倍是低成本开发回归，归一化采用相同 `Uref`。另核对每段二次插值的 Gauss 包络：稠密壁速不超过标记最大值的 `7/3` 加 `1e-12 Uref` 算术读回余量；`7/3` 来自三点 Gauss 拉格朗日基函数绝对值之和的上界，用于抓住旧漏点问题，不是新增工程壁面误差门。解析通道的 16×16 单网格 1% L2 界、粗圆柱的壁速上界仅为低成本开发冒烟检查，不属于工程验收。本次额外做了固定原折线的两个网格和周期 x 平移检查；它们用于揭示敏感性，不代替完整网格收敛/独立物理对标。未运行全量回归、跨平台或性能极限研究。

本入口没有真实共形流体 `polyMesh`，因此现有 `checkMesh` 和 Cut-cell Solver 质量门不适用；它们没有被这个独立读取器替代。方法背景：[Brinkman 固体惩罚法](https://www.math.u-bordeaux.fr/~chabrune/publi/ABF-NM.pdf)、[含惩罚项的压力投影预条件研究](https://arxiv.org/abs/2306.06277)。压力与壁面力联合约束的背景见 [Taira–Colonius 2007 原作者论文目录](https://www.seas.ucla.edu/fluidflow/pubs.html)（DOI: 10.1016/j.jcp.2007.03.005）。本实现为带有限顺应项的双线性表面惩罚，没有复现这些论文的整套算法、精度阶或性能结论。

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
