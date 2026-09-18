# 开发导航

## 分支与里程碑

`mesher-v0.3.0` 是固定标签（`664ac7c`）；`codex/mesh-maintenance` 是可继续修改的网格维护线；`codex/cfd-development` 是包含网格核心的 CFD 开发线；`main` 是已验收的集成线。不要移动已有里程碑标签去“更新版本”，后续里程碑另建标签。

纯网格修复先在维护分支保留最小回归，验证后合并到 CFD 开发线；合并冲突需人工核对，随后重跑受影响的网格和求解器检查，再集成到 main。跨分支改动不会自动同步。切换分支前处理当前工作区改动，切换后重新构建原生工具及 desktop runtime，避免旧运行包与当前源码混用。继续遵守唯一根目录和固定文档入口规则。

## 从哪里进入代码

生成链：输入轮廓 → 尺寸场/Quadtree → 真实流体 Cut-cell 或共形边界层 → 小单元处理 → 求解拓扑 → 质量 → 文件导出。

| 路径 | 责任与入口 |
|---|---|
| `apps/cartmesh2d_cli.cpp` | 纯 Cut-cell 总流程、尺寸场参数、物理面积门、Solver 质量和输出 |
| `apps/cartmesh2d_hybrid_cli.cpp` | 边界层总流程、局部修复开关、fallback 与导出 |
| `apps/cartmesh2d_flow_cli.cpp` | 原生稳态层流 CLI；SIMPLE / Rhie–Chow，桌面调用与诊断场导出 |
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

### 原生层流求解

```sh
cmake --build build --target cartmesh2d_flow_cli -j 4
build/cartmesh2d_flow_cli --mesh /path/case.solver.cm2d --output outputs/flow/result --case external --nu 0.1 --speed 1 --max-iterations 1500
python3 tools/verification/verify_native_flow.py --generate --output-root outputs/native-flow/reproduce --max-iterations 1500
MPLCONFIGDIR=/tmp/cartmesh-flow-mpl python3 tools/visualization/render_native_flow.py --summary outputs/native-flow/reproduce/summary.json --output outputs/native-flow/reproduce/figures
```

`external` 为左侧恒速入口、右侧运动学压力0、上下滑移及物面无滑移；`channel` 为无孔矩形内域、左侧抛物线入口、右侧压力0及上下无滑移，speed 是抛物线峰值；`cavity` 为无孔矩形腔、顶盖水平移动、其他壁面静止，speed 是顶盖速度，固定 cell 0 的压力为0。只接受一个连通流体区域；边界位置/方向或内域形状不符、出口回流、数值范围错误均明确失败。它不是任意喷管/多孔腔的自动边界配置器。

单元中心速度/运动学压力，共享边积分体积通量。动量对流为**一阶迎风**，黏性项用最小二乘梯度及显式非正交修正；内部面速度包含偏斜修正和 Rhie–Chow 压力项，压力修正使用四次非正交迭代，最终通量与最后一次实际线性方程一致。动量松弛0.6、压力松弛0.25；动量使用自行实现的 Jacobi–BiCGStab，压力修正利用对称正定结构使用 Jacobi–PCG，均检查真正矩阵残差。设计依据包括 [MOOSE 的同位有限体积说明](https://mooseframework.inl.gov/modules/navier_stokes/insfv.html)中关于 Rhie–Chow 和压力零空间的说明；没有复制或链接其求解核心。

本 CLI 的停止条件是：至少10次迭代，动量残差、相对速度变化、相对压力变化均小于 `--tolerance`（默认1e-6），逐格连续性及全局相对流量失衡均小于1e-8。动量残差是原离散方程失衡除以 `(aP_u+aP_v)*Uref`；逐格连续性为 `|sum(flux)|/(Uref*sqrt(area))`；速度变化以 Uref 归一化，压力变化以 `Uref²+nu*Uref/domainHeight` 归一化。全局失衡除以总入流，封闭腔使用 `Uref*domainHeight`。这些是本实现的数值停止条件，不是所有 CFD 软件的统一精度标准，不能与 OpenFOAM residual 数字直接等同。`converged` 也不等于网格无关或物理模型适用。

退出0代表满足上述停止条件；退出2代表到达上限，保存诊断场但 `converged:false`；退出1代表输入/数值失败。输出六种文件：`.json` 工况/状态/单位/压力基准、`.fields.json` 桌面字段、`.cells.csv` 单元 u/v/p、`.faces.csv` owner向外的积分体积通量、`.residuals.csv` 全迭代历史、`.vtk` 原多边形上的速度/压力。p 为 p/ρ，单位 m²/s²；力为流体对静止 EmbeddedBoundary 的积分力除以密度和深度，单位 m³/s²，不是 Cd/Cl。`domainHeight` 只指外域高度，不是物体参考直径。桌面验证全字段有限、单元ID/数量/工况与本次最终网格一致后才绑定；失败或取消保留 `flow-incomplete-*` 诊断，部分复制文件不会充当完整结果。

解析通道验证速度分布、压降梯度和流量；方腔 Re=100 对比 [Ghia 等（1982）](https://doi.org/10.1016/0021-9991(82)90058-4)中心线数据。圆柱只验证低 Re 定常试算、有限场和守恒，不与几何/边界不同的 DFG 基准混比。误差及外部工具实测范围见 CURRENT_STATE，绘图直接读取 CM2D/CSV。桌面 smoke 可加 `--flow=external --flow-nu=0.1 --flow-speed=1 --flow-max-iterations=30`，迭代上限场不得作为收敛证明。

### 扩展与性能测量

长期目标和当前进度只维护在 CURRENT_STATE。先完成稳态守恒/精度与压力求解效率，再依次扩展非定常、标量/热输运、SST 湍流、理想气体可压缩流；每种模式有独立物理配置和验收，不能用原稳态 laminar 的选项伪装支持新模型。

`cartmesh2d_flow_cli --profile` 在原六份结果之外增加 `.performance.json`。使用 monotonic `steady_clock` 测量读网格及构造时间、整个求解时间、动量/压力线性求解时间；统计线性求解调用、实际 Krylov 迭代总数和单次最大值。初始残差已满足条件的线性求解计0次迭代；每轮有两次动量和四次压力调用。线性耗时包括工作数组初始化与真残差检查，不含矩阵组装；总求解耗时还包含验证、组装、梯度、物理监测和进度回调，但不包含导出。此文件不测内存；需用系统工具另测峰值 RSS。到上限仍保存诊断且 `converged:false`，抛出数值错误的运行保留 stderr，不能当成完整性能样本。

```sh
build/cartmesh2d_flow_cli --mesh outputs/native-flow/formal-final/meshes/channel-l6/channel-l6.solver.cm2d --output outputs/flow-profile/channel --case channel --nu .01 --speed 1 --max-iterations 1500 --profile
# macOS 可加 /usr/bin/time -l；maximum resident set size 为 bytes。
# Linux /usr/bin/time -v 的 Maximum resident set size 单位为 KiB，不能直接混比。
```

相同几何/网格/工况/容差/线程数下比较，记录二进制哈希和系统版本。固定迭代吞吐不能冒充达到相同物理解精度的加速；完整求解须同时核对误差与守恒。初步小规模观测不外推50万格速度。前端继续使用原六份物理文件，本阶段未重新打包 App。

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
