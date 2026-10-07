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

桌面续算重载复用同一 smoke 入口：`--flow=channel --flow-dt=.02 --flow-steps=2 --flow-checkpoint=/绝对路径/状态.checkpoint` 通过流动载入按钮导入；`--thermal-checkpoint=/绝对路径/thermal.thermal.checkpoint` 通过温度载入按钮导入并推进两步，步长默认 .025 s，可用 `--flow-dt=` 改变。联合文件旁须保留 `desktop-state.json`，两者均须使用原网格参数重新生成网格。路径替代仅在 `--smoke=` 下作用于文件选择对话框，解析、控件、IPC 和原生续算仍使用产品代码；这不验收系统文件对话框的人工选择。温度重载还检查从局部涡／自动步长设置切换到联合状态的操作。

### 完整桌面项目

`desktop/src/core/project.js` 将当前绑定保存到结果包内 `cartmesh2d-project.json`（`cartmesh2d-project-v1`），使用项目内相对路径与原始文件 SHA256。转换后的 XY 已是米制，恢复时不能再次使用原 DXF/图片单位缩放。外部选入的接受检查点一起打包；场与历史重新从原生文件读取，失败/取消的较新状态与上次完整结果分开保存。`core/archive.js` 流式读取 ZIP，拒绝越界路径、重复文件、链接和特殊文件；只在完整读回成功后发布绑定。已有 ZIP 采用同文件系统原子替换。

真实 App 的跨进程项目验收入口：先正常运行并 `--export=/绝对路径/project.zip`，退出后使用 `--smoke=rectangle --open-project=/绝对路径/project.zip --out=/绝对路径/验收目录` 直接打开。可附加 `--project-resume=thermal --project-dt=.025 --project-end-time=1.4`；独立流动用 `--project-resume=flow --project-dt=.01 --project-steps=2`。支持已有 `--shot`、`--export`；不会先生成新网格。路径参数仅替代 smoke 文件对话框，正式操作通过顶部“打开项目”。原生终点回归在 `tests/thermal_flow_test.cpp`，前端文件恢复/取消检查在 `desktop/tests/project.test.js`、`archive.test.js` 和 `flow-checkpoint.test.js`。

### 温度连续历史与会话恢复

`core/thermal-history.js` 沿 `thermalRestart.history` 的真实运行区间读取原生接受 CSV，检查非重叠时钟和检查点终点；检查点 SHA256 绑定选入状态。`core/thermal.js` 共用接受步读取规则：完整段须达到其完整结果终点，失败/中断段可缺末尾未刷出的行，但已完整写出的坏行仍拒绝。缺少采样的区间返回显式 gaps，renderer 不跨缺口连线；旧项目只从与选中检查点路径相同的完整段建立历史。表格/场图仍来自各自完整结果，聚合历史不覆盖场的物理时钟。

`core/recovery.js` 在应用数据目录 `recovery/session-*` 保留原项目清单与独立恢复记录。温度作业启动前保存已关闭文件的 SHA256 清单，再登记新运行目录；`core/process.js` 的 onSpawn 写入原生进程 PID。恢复时确认没有其他存活的源 App 或原生子进程，再由当前 App 排他复制工作目录并验证基线项目。对中断的新检查点，`checkThermalRestart` 调用 `cartmesh2d_transport_cli --check-restart on`，由既有原生完整读取器验证网格、边界、物性、方法和完整时间事件；之后才发布更新时钟和历史。取消、原生核对失败均不发布未经验证的新状态，源会话不改动。该只读选项需要原正常物理输入及 `--restart`；允许目标时间等于检查点时间，不产生输出文件。

真实 App 强制退出验收使用独立恢复目录：在温度事件 smoke 参数后加 `--recovery-root=/绝对路径/store --thermal-crash-after=.7`。先完成 .5 s，再续算并在接收达到 .7 s 的接受进度时 SIGKILL **本次测试 App 自身**；原生程序可能多接受一步，以落盘检查点为准。用新进程 `--smoke=rectangle --recover-latest=true --recovery-root=/绝对路径/store --project-resume=thermal --project-dt=.025 --project-end-time=1.4 --thermal-events=true --out=/绝对路径/recovered --shot=/绝对路径/recovered.png --export=/绝对路径/recovered.zip`，通过真实恢复对话框、renderer、IPC 和原生程序继续。`--recover-latest` 与强制终止开关仅在 smoke 分支生效；正常 App 由用户选择“恢复会话”。

历史/项目/进程存活/取消回归分别在 `desktop/tests/thermal-history.test.js`、`project.test.js`、`recovery.test.js`；只读原生核对与截断文件回归在 `tests/thermal_control_cli_test.py`。实际逐字节续算对照脚本、调用参数和源文件索引在 `outputs/thermal-recovery/`，旧的预算失败连续历史证据在 `outputs/thermal-history/`。独立流动使用下节相同的会话机制；断电持久性和多平台范围见当前状态。

### 独立流动连续历史与恢复

`core/flow-job.js` 保留每次运行的 `flow-run-*/flow.*`，只有完整读取并通过既有结果检查后才替换显示绑定；失败/取消只更新真正有新增接受步的 `flowRestart`。不再把文件复制覆盖到固定 `.flow.*` 路径。选入状态使用 SHA256 绑定；`core/flow-history.js` 聚合其 `history` 区间，复用 `core/flow.js` 的原生 CSV 规则，拒绝完整坏行、重叠时钟及不匹配终点，未写全的末尾显示缺口。旧项目仅使用能证明与选中检查点同源的历史；旧 CLI 先写历史后写检查点形成的超前尾部留在磁盘并计数，不接入接受曲线。

流动 CLI 在每步接受后先原子保存检查点，再刷新接受历史，避免可见历史领先已保存状态。`--check-restart on` 需要 `--restart` 及原物理参数，由 `readFlowCheckpoint2D` 完整核对几何、物性、边界及场，返回 `cartmesh2d-flow-restart-check-v1`；不推进时钟、不创建原生输出目录或写边界快照。自适应目标可等于检查点时刻。该入口与截断拒绝在 `tests/transient_flow_cli_test.py` 回归。

`RecoveryStore.beginFlow` 在运行前保存请求/目录，原生 PID 与温度联算相同登记。恢复副本上的 `checkFlowRestart` 对命名边界写出原请求，再调用只读原生核对，全部通过后才发布新状态和历史。流动与温度状态共存时保持各自历史，按最近运行类别恢复续算控件。窗口关闭会取消当前操作，renderer 销毁后停止向它发送进度，并用保存的 UI 快照记录关闭状态。同进程闲置会话允许恢复，当前操作排他锁及其他存活 App/子进程仍阻止并发写入。

真实 App 固定步验收使用 `--smoke=rectangle --target-cells=1000 --flow=channel --flow-dt=.02 --flow-steps=2 --flow-resume-steps=2 --flow-max-iterations=1500 --flow-nu=.1 --flow-speed=.2 --flow-convection=limited-linear --flow-failure-check=true --flow-cancel-check=true --flow-require-converged=true`，另指定 `--out`、`--shot`、`--export` 绝对路径。ZIP 新进程重开使用前述 `--project-resume=flow` 入口。

强制退出或关窗验收使用独立 `--recovery-root=/绝对路径/store`，同一 1300 格配置改为 `--flow=custom --flow-end-time=.04 --flow-crash-after=.07`；正常关窗改为 `--flow-close-after=.07`。先完成 .04 s，再运行长续算，于阈值后强制终止**本次测试 App 自身**或关闭其窗口；最终接受时刻取实际文件。新进程使用 `--smoke=rectangle --recover-latest=true --recovery-root=/绝对路径/store --project-resume=flow --project-dt=.01 --project-end-time=.16 --out=/绝对路径/recovered --shot=/绝对路径/recovered.png --export=/绝对路径/recovered.zip`。所有注入仅在 smoke 下生效，正常恢复通过顶部按钮。

核心回归在 `desktop/tests/flow-job.test.js`、`flow-history.test.js`、`project.test.js` 和 `recovery.test.js`；实际 App/CLI 逐字节对照脚本、参数、源场与文件哈希在 `outputs/flow-history/`。对照复用实际恢复后的输入和控制，不将分段终点改变产生的新时间网格声称为同轨迹。

### 桌面热时间事件

`core/thermal.js` 归一化请求的 `events`，按绝对时刻和目标排序，拒绝同时间重复目标；源项记录 `{time,target:'source',kind:'source',value}`，边界记录 `{time,target,kind,value,inflowValue}`。目标为既有 wall/inlet/outlet/top/bottom，按当前工况与真实最终网格映射到 `face:ID`；没有对应面的目标明确失败，不按屏幕邻近猜测边界。桌面事件要求自动步长模式；空列表保留原有静态流程。

`core/thermal-job.js` 写出 `.events.csv` 并传给既有 `--thermal-events`；基准热边界和完整事件规律共同进入原生检查点身份。续算先比较完整归一化规律，原生再次检查完整源项/边界快照；元数据读取器只核对 v2/v3/v4 的时间表与时钟，不能替代原生物理核对。`desktop-state.json` 与完整项目清单均保存事件，未来时刻也不丢弃。任意命名 custom 流动边界仍通过 CLI 配置。

真实 App 事件验收复用 `--smoke=rectangle --target-cells=1000 --thermal-adaptive=true --thermal-flux-correction=bounded --thermal-events=true`，再加 `--out`、`--shot`、`--export` 绝对路径。此组合包含同时改变热源与上壁、预算失败和恢复，保存时还保留未来事件；用上一节项目重开入口续至 1.4 s，验证 1.337 s 事件。仅查看已保存事件可用 `--open-project=... --thermal-events=true --shot=...`，不必重新计算。核心与界面测试在 `desktop/tests/thermal.test.js`、`thermal-renderer.test.js`；原始 App/CLI 对照在 `outputs/thermal-events/`。

冻结载流器的 `--flow-checkpoint` 导入读取原生 flow checkpoint v1/v2/v3/v4 配置，再调用既有 `readFlowCheckpoint2D` 完整验证；v4 带完整命名边界，LF/CRLF 均可读入。回归使用实际联合计算生成的 v4，核对冻结面通量和载流器时钟，并拒绝坏边界、网格和截断文件，见 `tests/thermal_control_cli_test.py`。

## 非定常守恒有界温度修正

入口为 `ScalarTransportControls2D::fluxCorrection`、`ScalarTransport2D.cpp` 和 CLI `--flux-correction unrestricted|bounded|bounded-spatial`。默认 `unrestricted` 保持原有离散。两种有界方法都要求非空上一物理层和正 dt；稳态调用明确失败。支持 Upwind/LimitedLinearUpwind、逐面正扩散率和非负隐式损失。算法背景可参阅 [Lipnikov、Svyatskiy、Vassilevski 对非正交网格离散最大值原理的讨论](https://dodo.inm.ras.ru/research/_media/lip-svy-vas-12.pdf)；本仓库采用下述隐式行余量构造，资格依据是原生算例。

记低阶 BE 矩阵为 `Aii Ti - Σ cij Tj = bi`，`cij≥0`，来自 TPFA 扩散和一阶迎风；高阶面修正 `Ff` 包含非正交扩散及对流重构。每格 `mi=Vi/dt`、`si=Vi*sinkRate`、`ri=Σ outward volumeFlux`，`Qi=source*Vi-Σ prescribed outward diffusiveFlux*faceLength`。比较值为 `(mi*Told+Qi)/(mi+si+ri)`；分母必须正，仍检查真实载流器连续性。全局 L/U 包含所有上一层、定温/流入值和比较值，因此真实加热/冷却或指定热流可以超出旧场极值。实现使用差值形式减轻 300 K 背景抵消。

在当前非线性迭代场上，正/负修正可用余量分别为 `Bi+=Aii*U-bi-Σ cij*Tj` 和 `Bi-=bi+Σ cij*Tj-Aii*L`。按 owner 的 `-Ff`、neighbour 的 `+Ff` 累计正/负修正总量 Pi±，令 `Ri±=min(1,max(0,Bi±)/Pi±)`，Pi± 为零时取一。每个面取两端对应符号系数的较小值，同时缩放该面的两种修正。指定 Neumann 热流、源项和损失不缩放，每面仍只有一份守恒通量；接受场不裁剪。

固定点的最大值判断：若最高格超出 U 且 Bi+≥0，受限右端直接要求该格≤U；若 Bi+<0，正修正为零，正比较锚和低阶矩阵又排除超过所有比较数据的全局最高值。最低值对称。该论证针对收敛离散方程，不保证任意中间候选或非线性迭代次数。返回场必须同时满足原有总残差、逐格失衡与有界偏差检查；后者直接使用已有 `cellTolerance`（温度时单位 K，默认 1e-9），不新增物理精度阈值。每步继承上一场后，求解容差可以累积；实际 100 步约 4e-9 K 的下界漂移单独报告。

`evaluateScalarTransport2D` 使用完全相同的受限面通量，在给定场重新检查守恒；不能以低阶线性候选收敛代替最终非线性方程。JSON 新增方法、受限面数、最小系数、比较上下界和最终超界量。`bounded` 的余量依赖本步质量与比较界，整步/两半步缺陷也包含限制器响应。限制器未启动的光滑区保持原重构；现有混合边界线性场和正弦空间/时间收敛例通过，尚不宣称任意非正交网格上的线性保持、唯一性或全程误差上界。

`bounded-spatial` 在上述余量中分别去掉 `mi*(U-Told)` 和 `mi*(Told-L)`。比较界包含所有上一层值，故去掉的是非负时间余量，新允许量不大于原行允许量，原有收敛场最大值论证仍成立；共享面系数及源、损失、指定热流保持。给定相同当前场、上一场及比较界时，空间修正不再直接依赖 dt。零净源、零损失且载流器严格守恒时，比较界也不随 dt 改变；有源、损失、指定热流或实际连续性缺陷时，比较界仍可能依赖 dt。它限制得更严格，尖锐前沿可能改变，须比较实际温度和热流，不能推成任意工况的时间精度保证。

两种有界联合状态都采用 thermal checkpoint **v4**：`THERMAL_CONFIG` 明确保存 `bounded` 或 `bounded-spatial`，`CONTROLLER_PRESENT 0|1` 区分固定与自动步；11 项控制签名末项分别为 1 和 2，始终保存完整 EVENTS。`unrestricted` 继续写原 v1/v2/v3 和 10 项签名。三个方法之间互相续算明确拒绝；在内存中改方法会重置控制建议。桌面请求、项目清单、结果验证和恢复控件共同绑定选项，LF/CRLF 都保持方法身份；旧项目缺少字段时按 `unrestricted` 恢复。

```sh
cmake --build build --parallel 4 --target cartmesh2d_scalar_transport_tests cartmesh2d_scalar_workspace_tests cartmesh2d_thermal_flow_tests cartmesh2d_transport_cli
ctest --test-dir build --output-on-failure -R '^cartmesh2d_(scalar_transport|scalar_workspace|thermal_flow|thermal_control_cli|thermal_flow_cli)$'
python3 tools/thermal/workflow.py --case cylinder --level 4 --dt .05 --end .3 --diffusivity .05 --speed .2 --no-events --flux-correction bounded --output outputs/bounded-cylinder
```

改成 `--flux-correction unrestricted` 并更换输出目录获得同控制基线。原云端零流速失败输入从已逐项验证的恢复归档取回；本轮完整命令在 `outputs/bounded-transport/*.command.json`，输入 SHA256、数值结果、适用范围见 `artifacts/current/native-thermal-bounded.json`。该零流速研究只比较原生输出，未恢复独立 Python 方程审计链。真实 App 入口追加 `--thermal-flux-correction=bounded`，例如 `--smoke=rectangle --target-cells=1000 --thermal-adaptive=true`；项目重开沿用前节入口。相对输出目录现在在保存时转换包内引用，原失败 ZIP 仍保留。

### 空间有界修正与启动对照

现有 `workflow.py` 支持 `--flux-correction bounded-spatial`；真实 App 使用 `--thermal-flux-correction=bounded-spatial`，其余事件、时间控制和项目重开入口相同。必须从新初始工况起算，不能将旧有界检查点改名后导入。

`outputs/thermal-startup/` 保存冻结实际状态的整步／两半步比较、96 格固定场面通量诊断、正式方法的长时与时间步对照、热核原生求解及 App 项目链。`run-long.py` 和 `run-time.py` 复用 Pe=30 的原命令，只显式修改方法、最大步长、二进制及输出前缀；重算须换目录。`analyze.py` 只读取已完成的原生场、通量和历史，不重建方程。`spatial-candidate-*`、`full-spatial/`、`halves-spatial/` 和 `spatial-startup/` 是隔离研究程序，它曾复用旧模式标记，不能作为正式新方法的续算文件；正式程序输出必须明确为 `bounded-spatial`。方法及实际范围见[当前状态](CURRENT_STATE_CN.md#持续目标成熟非定常与被动温度)；展示素材为[原生对照图](../artifacts/current/native-thermal-spatial-bounds.png)和[真实项目续算](../artifacts/current/native-thermal-spatial-bounds-app.png)，来源见[关键证据](../artifacts/current/native-thermal-spatial-bounds.json)。

### 联合温度局部涡初值

`cartmesh2d_transport_cli` 的新起算可传入 `--initial-vortex-x X --initial-vortex-y Y --initial-vortex-radius R --initial-vortex-speed V`，四项同时填写。X/Y/R 单位 m，V 单位 m/s，正值逆时针；仅适用于实际物理流动联算的零时刻，不能和 `--restart` 或解析验证混用。实现直接复用 `withInitialVortex2D` 的紧支撑流函数；圆盘须严格处于可分辨的流体区域，不能穿墙或接触开口。`.initial.checkpoint` 是零时刻流动诊断，联合物理状态仍只在接受正时间步后写出。

桌面沿用“初始局部涡扰动”控件；`thermalRequest`、CLI 调用、结果元数据与项目文件共同保留四项参数。联合状态导入核对使用已存载流场，同时保留原请求作为初值来源；继续计算不重复施加扰动。原生接受步 JSON 进度增加守恒受力 `forceX/forceY`（m³/s²），乘密度才成为单位深度作用力，与 `thermal-time-step` 同时刻；现有 CSV 格式不变。

原生回归包括连续／分段逐字节一致、零时刻与首步失败语义、不合法初值及重复施加拒绝，可运行：

```sh
python3 tests/thermal_control_cli_test.py --cli build/cartmesh2d_transport_cli --mesh-cli build/cartmesh2d_cli --output outputs/thermal-seed-check
```

App 既有腔体自动联算 smoke 加 `--thermal-initial-vortex=true`，使用中心 (.5,.5) m、半径 .2 m、峰值 .05 m/s；仍走真实 renderer／IPC／原生程序，文件选择由 smoke 参数代入。`outputs/thermal-wake/` 保存打包 App、904 格完整项目、CLI 逐字节对照和首段原生载入核对。独立 5168 格 Re=100／Pe=100 尾迹研究使用另一组初值，完整输入见 `protocol.json`；`run.py NEW_LABEL END_TIME [ACCEPTED_CHECKPOINT]` 使用冻结程序和新目录，已有输出不得覆盖。真实周期发展尚按[当前状态](CURRENT_STATE_CN.md#持续目标成熟非定常与被动温度)判断。

### 温度联算的系统压力后端

`cartmesh2d_transport_cli --pressure-preconditioner cholesky` 在 macOS 复用 `detail/FlowCholesky2D.hpp` 的系统稀疏分解作为 PCG 预条件器；仍由原 CSR 矩阵的真实残差判定接受。程序在首次 Accelerate 调用前设置 `VECLIB_MAXIMUM_THREADS=1`，保持同控制跨进程轨迹可复现。默认 IC0、聚合多重网格及原物理身份保持；后端是数值选项，续算可更换，但不同后端只比较实际误差，不要求不同后端的字节相同。平台不可用或分解失败明确报错，没有隐藏降级。

```sh
python3 tests/thermal_control_cli_test.py --cli build/cartmesh2d_transport_cli --mesh-cli build/cartmesh2d_cli --output outputs/thermal-cholesky-check --pressure cholesky
```

该入口在 macOS 对所选后端运行完整原生交易、事件、取消和有界续算回归；默认调用另包含 Cholesky 的固定步局部涡恢复例。非 macOS 默认回归验证不可用后端被拒绝。真实 App 的自动温度 smoke 可加 `--flow-pressure-preconditioner=cholesky`；本次完整包、904 格项目和同控制 CLI 对照保留在 `outputs/thermal-cholesky/`。

`outputs/thermal-wake-resolved/mesh-protocol.json` 绑定新 128 边形／两级局部加密、原生质量报告和输入哈希；`prepare.py` 记录生成方式，不能用旧 32 边形算例补作第三网格。`run.py` 保留原程序基线，`run-cholesky.py MESH NEW_LABEL END_TIME --restart ACCEPTED_CHECKPOINT --max-courant CFL` 使用冻结的新程序；`MESH` 为 `medium` 或 `fine`，不传 restart 才施加协议中的零时刻初值。后端比较从同一 2 s 状态出发，CFL 比较从同一 2.2 s 状态出发；原始命令在各目录的 `command.json`。`compare-cfl.py` 只读取原生输出作敏感性汇总，不重建方程。重跑须换目录，保留原输入及检查点。实际长算与精度边界见[当前状态](CURRENT_STATE_CN.md#持续目标成熟非定常与被动温度)。

长时间同态对照保留 `medium-growth-4to20/result.thermal.checkpoint` 作为共同 20 s 状态；`medium-growth-20to60` 与 `medium-dt025-20to60` 只改变最大 dt=.05／.025 s，均使用完整 BE 整步／两半步控制。`response-analysis.py --segment SEGMENT --segment NEXT_SEGMENT --output NEW_REPORT --label LABEL` 仅接受已完成的目录，逐段核对物理时刻及输入检查点哈希，合并原生受力、热量和壁面通量；完整的段清单在关键证据的 `trajectories` 中。升阻力系数按 `2(F/ρ)/(U²D)` 定义，总平均 Nu 按 `−Qwall/(π α ΔT)` 定义（U=1、D=2、α=.02、ΔT=1）；128 边形真实通量按名义圆周归一。完整周期由正向升力过零时刻及线性插值确定，周期统计对原生样本作梯形积分，不自动判断稳态。

`python3 outputs/thermal-wake-resolved/compare-time.py medium-at60 medium-dt025-at60 NEW_REPORT` 核对分叉前同一状态、分叉后仅最大步长不同，输出累计相位／幅值差、末态温度／速度差和原生接受误差比。全域面积 RMS 包含远场，不能替代尾迹局部或最大误差。`time60-raw-manifest.json` 封存本次完成的输入链补充、原生结果和后处理，前段输入由已有 Cholesky 清单绑定；活动目录明确排除。继续 60→120 s 时分别读取各自接受态，不覆盖 60 s 结果；实际完成范围见当前状态。图只渲染原生多边形及场，不建立独立方程验证链。

180 s 基准、120 s 时间对照及细档 60 s 结果保存在 `medium-at180`、`medium-at120`、`medium-dt025-at120`、`fine-at60`；各 `analysis.json` 保存完整续算段清单。`compare-time.py medium-at120 medium-dt025-at120 NEW_REPORT` 复现长期时间差；`compare-mesh.py medium-at60 fine-at60 NEW_REPORT` 按命名边界核对物理条件（两档面编号不同），只统计原生响应，不重建方程或跨网格映射场。新的 .025／.0125 s 分叉共同读取 `medium-dt025-60to120/result.thermal.checkpoint`；运行和后处理均使用新目录，不覆盖输入。`time180-raw-manifest.json` 只封存完成部分并绑定此前清单。

### 联合推进的线性迭代失败恢复

`FlowLinearIterationLimit2D` 只标记压力 PCG／动量 BiCGStab 的有限迭代预算耗尽；真实残差、对角归一残差和原停止门均保留。自适应 `advanceThermalFlowImpl` 仅在流动调用处把该类型转换为 `flow-linear`／`flow-linear-half`，所有试步都从原接受态开始；半步失败也拒绝整步。控制器复用既有松弛阻尼、缩步、重试预算，原始诊断随 `ThermalAttempt2D` 写入 CLI `thermal-time-retry` 和失败摘要，不新增通过状态 CSV 列。固定步仍传播异常，矩阵错误、数值崩溃及取消不作为该类失败处理。

相关原生回归为 `cartmesh2d_thermal_flow`、`cartmesh2d_scalar_transport`、`cartmesh2d_flow_linear`；9 格腔体用接近舍入底限的压力范数限制构造确定失败，不是新的生产精度阈值。`outputs/thermal-linear-retry/run-regression.py` 保留真实 23569 格、60 s 输入的自动恢复／零重试预算／旧通过结果字节比较；使用已冻结程序和全新目录。`tests/thermal_control_cli_test.py` 继续承担原生交易、事件、取消及完整续算检查，不增加独立方程验证链。原始数据索引为该目录的 `raw-manifest.json`。

真实系统文件选择验收保存在 `outputs/thermal-dialogs/`：`manual-project-1p4.zip` 是通过 NSSavePanel 保存的原包，`verify-save.py` 检查 ZIP、状态和 7 类原生文件，`process-after-reopen.json` 绑定新的 App 进程。`package/mac-arm64/CartMesh2D.app` 已更正 Cholesky 帮助，但物理程序仍为该轮冻结版本；全窗口截图仅在界面工具中查看，未另存本地，不能把导出的场预览当作完整 App 截图。

更新后 App 的实际续算包在 `outputs/thermal-linear-retry-app/project-1p45.zip`，`comparison.json` 绑定新程序哈希、运行进程和 7 类文件字节一致；`run-cli.py` 用正式桌面参数构造器生成同控制参考，`verify-app.py` 读取真实 ZIP 核对。新包完成 1.4→1.45 s／107→109 历史；此前 `outputs/thermal-dialogs/` 包保留原物理版本。原线性失败的松弛系数 1 由 CLI 复现，当前桌面温度只允许 .6，未为了测试开放额外控制。

### 周期腔体与联合载流精度

原生标量载流门检查每格 `abs(Σq) <= carrierAbsoluteTolerance + carrierRelativeTolerance*Σabs(q)`，q 单位 m²/s。闭合压力系统钉住一行，其误差是其余行误差之和；原线性停止量的二范数不能直接作为该行的绝对预算。`advanceThermalFlow2D` 只捕获明确的 `ScalarCarrierContinuityError2D`，第一次不匹配后从原接受态重算相同 BE 步，把 `pressureResidualNormStop` 限为 `0.5*carrierAbsoluteTolerance/sqrt(max(1,N-1))`。Cauchy–Schwarz 给出遗漏行误差上界，0.5 为原预算留出舍入余量，实际面通量仍须通过原门。该上限单位也是 m²/s，只影响数值精度，默认无限时压力求解完全保留旧停止条件；不改物理状态身份。再次不匹配返回 `carrier`／`carrier-half`，由已有控制器回退重试；无效输入及其他异常不伪装成此类失败。

`outputs/thermal-periodic/` 保存原／新冻结程序、三档最终网格、事件表、每半周期的新进程命令及完整原生输出。`run-fixed.py` 每 2 s 读回上一接受状态，周期目标 48 s、完整事件规律到 64 s；重算必须换输出目录。`analyze.py` 只读实际 CSV，比较同相位场、积分原生热流、计算热增量基频与相位，并比较未修复的通过例和失败恢复；不重建 PDE。壁流积分按原生接受步的右端值乘 dt，基频按接受历史梯形积分，输入基频为每周期前半段热壁的正弦分量。离散采样和时间误差保留在报告中。

以下命令从现有原生工具生成同物理控制的一次连续运行；示例目录须尚未使用，已有接受检查点的输出前缀不能复用：

```sh
mkdir -p outputs/periodic-cavity-repeat
python3 - <<'PYINPUT'
from pathlib import Path
p=Path('outputs/periodic-cavity-repeat/geometry.xy')
p.write_text('0 0\n1 0\n1 1\n0 1\n')
PYINPUT
build/cartmesh2d_cli outputs/periodic-cavity-repeat/geometry.xy outputs/periodic-cavity-repeat/mesh 6 0.07142857142857142 .1 interior outputs/periodic-cavity-repeat/openfoam 6 0
python3 - <<'PYINPUT'
from pathlib import Path
from tools.thermal.workflow import configure
p=Path('outputs/periodic-cavity-repeat')
configure(p/'mesh.solver.cm2d','cavity',p,speed='.2')
(p/'events.csv').write_text('time,target,type,value,inflowValue\n'+''.join(f'{2*k},lid,value,{300 if k%2 else 301},300\n' for k in range(1,33)))
PYINPUT
build/cartmesh2d_transport_cli --mesh outputs/periodic-cavity-repeat/mesh.solver.cm2d --output outputs/periodic-cavity-repeat/result --evolve-flow custom --flow-boundary outputs/periodic-cavity-repeat/flow.boundaries --boundary outputs/periodic-cavity-repeat/thermal.csv --flow-nu .01 --flow-speed .2 --flow-velocity-relaxation .6 --diffusivity .02 --initial 300 --dt .025 --end-time 48 --flow-max-iterations 1500 --flow-convection limited-linear --convection limited-linear --flux-correction bounded-spatial --min-dt 1e-7 --max-courant 1 --max-step-retries 24 --max-time-steps 100000 --time-error on --temperature-scale 1 --velocity-scale .2 --time-rtol .01 --thermal-events outputs/periodic-cavity-repeat/events.csv
```

把网格两处 level 改为 4／5 得到 196／784 格，保留其他控制；同 784 格把最大 dt 改为 .1／.05／.025 s 比较时间敏感性，分别使用新目录。真实 App 的既有验收入口现支持 `--thermal-case=cavity`；配合单位方腔几何、`--flow-nu=.01 --flow-speed=.2 --thermal-diffusivity=.02 --thermal-flux-correction=bounded-spatial --thermal-events=true --flow-min-dt=1e-7 --flow-max-retries=24`，通过实际 renderer／IPC 完成事件、失败续算、项目保存；项目重开参数沿用[项目恢复入口](#完整桌面项目)。状态与实际测试范围见[当前状态](CURRENT_STATE_CN.md#持续目标成熟非定常与被动温度)，素材为[周期对照图](../artifacts/current/native-thermal-periodic.png)。

### 有界圆柱长期响应复现

原生同控制长时对照沿用前述云端恢复归档中的 `outputs/thermal/qualified/cylinder-{coarse,fine}/` 网格、原始热/流边界和事件文件，使用保存的命令，仅更换二进制/路径及 `--flux-correction`。`outputs/bounded-long/run.py` 串行运行两种方法，`outputs/bounded-long-time/run.py` 对细档有界方法将最大 dt 连续减半；活动输出在独立临时目录，进程关闭后保留并导入工作区。输入、二进制、命令及闭合输出哈希均在 `artifacts/current/native-thermal-bounded-long.json`，不覆盖原云端场。

重新生成同类代表例可用既有入口（新输出目录）：

```sh
python3 tools/thermal/workflow.py --case cylinder --level 5 --end 25 --dt .1 --diffusivity .05 --nu .1 --speed .2 --error --flux-correction bounded --isolate-live-output --output outputs/cylinder-bounded-25s
```

将 level 改为 4 比较粗网格，将方法改为 unrestricted 比较旧方法；最大 dt 改为 .05/.025 做时间敏感性，并分别使用新输出目录。新生成网格应核对实际格数及输入哈希，不能冒充原归档的完全相同输入。`outputs/bounded-long/analyze.py` 只读取原生接受 CSV/单元场，计算热增量、积分已输出的壁面通量/预算缺陷及同网格场差，并画科学图；未重建独立 Python 方程审计链。最大 dt 对照仍有 CFL/误差控制，不等同于均匀固定时间网格；曲壁网格和早期前沿精度限制在当前状态集中记录。

### 曲壁温度空间对照

原生手动入口 `tests/thermal_spatial_probe.cpp` / `cartmesh2d_thermal_spatial_probe` 在开启 `CARTMESH2D_BUILD_TESTS` 的常规 build 中构建；它需要显式提供最终 Solver 网格，不作为无输入 CTest。使用现有 `readCm2dTopology`、`makeFvMesh2D`、标量算子/求解器及 ILU0，原残差门不变。Python 仅编排原生调用、读取真实 CSV 和画图，不重建独立方程审计链。

解析参照都使用 D=.05 m²/s、零载流速度、无体积源、真实多边形边界的时刻值。稳态 `T=300 K+(1 K)·log(r/1 m)`，奇点置于原点所在固体孔洞内；非定常 `T=300 K+(1 K)·(.1 s/τ)·exp(−r²/(4Dτ))`，`τ=.1 s+t`，从 t=0 推至 .2 s。这是 [Stanford 热方程讲义 §2.4.2 式 (2.21)](https://web.stanford.edu/class/math220b/handouts/heateqn.pdf) 的二维热核经平移与缩放所得。边界离散值取面中心，参考热流则对实际每条直边积分：log(r) 用带方向端点极角，热核用切向高斯的 erf 积分，不把真实 64 边形换成圆。

单元误差为相对解析质心值的面积加权 L2（K）和最大误差（K）。`centroidHeatDifference=Σ area·(T−Texact,centroid)`，单位 K·m²，不是解析连续体体积分。总壁热流相对误差分母为 `|ΣQexact|`，局部 L1 为 `Σ|Q−Qexact|/Σ|Qexact|`，最坏通量密度误差为 `max(|Q−Qexact|/faceLength)`（K·m/s）。所有通量均向流体域外为正；乘 rho·cp 才成为单位深度的物理热流。稳态 `ΣQexact=2πD·(1 K)` 的检查只确认参考孔洞与浮点求和，额度 `128·eps·壁面数·|2πD·(1 K)|` 不是 PDE 精度阈值。探针的温度范围诊断包含解析初值；初值不作为物理续算检查点。

下面使用随源码保留的同一 64 点原始几何夹具，另建输出目录，生成 level 6 网格并运行两个原生参照；改为 level 4/5/7 得到另外三档，需更换输出前缀。

```sh
cmake -S . -B build -DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCARTMESH2D_BUILD_TESTS=ON
cmake --build build --target cartmesh2d_cli cartmesh2d_thermal_spatial_probe --parallel 2
mkdir -p outputs/thermal-spatial-rerun
build/cartmesh2d_cli tests/fixtures/thermal_cylinder64.xy outputs/thermal-spatial-rerun/mesh 6 3 .1 exterior outputs/thermal-spatial-rerun/openfoam 6 0
build/cartmesh2d_thermal_spatial_probe outputs/thermal-spatial-rerun/mesh.solver.cm2d outputs/thermal-spatial-rerun/harmonic harmonic
build/cartmesh2d_thermal_spatial_probe outputs/thermal-spatial-rerun/mesh.solver.cm2d outputs/thermal-spatial-rerun/pulse pulse 400 bounded
```

pulse 的步数 200/400/800 对应 dt=.001/.0005/.00025 s，模式为 unrestricted、bounded 或 bounded-spatial；harmonic 为原有稳态 unrestricted。标准输出是指标 JSON，输出前缀旁保存单元、壁面和接受历史 CSV。`operatorOnExact` 是原算子作用在解析场，`quadraticFitOnExact` 只是既有二次拟合诊断，`solvedFlux` 才是实际求解后的热流；不得混为方法收益。失败明确返回非零，保留已经写出的诊断。

三档有流长期工况沿用上一节原始归档的粗/细输入；新增 level 6 保持同几何与物性，按实际最终面导出边界，将上下外域面设为对称 farfield，再写定温壁面/入口和零外热流。入口 .2 m/s、ν=.1、D=.05，7.13 s 壁温 301→300 K、14.23 s 回到 301 K，最大 dt=.05 s、有界与联合时间误差控制开启，终点 25 s。完整 argv、输入和原生 CSV 在 `outputs/thermal-spatial/{coarse,fine,finer}-dt005/`；`run.py`、`run-probes.py` 记录执行方式，所有额外解析调用也各有 command.json。复算请使用新目录，不覆盖原接受场/检查点。

`outputs/thermal-spatial/analyze.py` 读取上述输出，比较同网格减半时间步与切换方法的末态、相邻网格总量，生成 `artifacts/current/native-thermal-spatial.json/.png` 和本地逐文件哈希。`attempt-direct/`、`attempt-stabilized/` 与对应 harmonic 失败目录保留两种未采用的边界二次修正；旧 `wall_probe.cpp` 和 `cartmesh2d_thermal_wall_probe` 是失败研究版本，不能当作当前支持入口。生产 ScalarTransport/WallGradient 已按 `quadratic-decision.json` 的哈希恢复。当前结果与局部壁热流限制集中在状态文档；本增量不改默认算法。

### 对流较强的曲壁温度与首步失败说明

Pe=30 代表例使用上节同一 64 点夹具与 level 5/6 网格，保持入口 .2 m/s、ν=.1 m²/s、301/300 K 壁温阶跃，仅将 D 改为 .002 m²/s；最大 dt=.05 s、温升尺度 1 K、相对时间容差 .01、有界修正。`outputs/thermal-advection/protocol.json` 和各目录的 `command.json` 保存完整参数，`run.py`、`run-time.py` 编排同一原生二进制；`analyze-completed.py` 只读取闭合的真实 CSV，计算热量、已输出通量积分与同网格场差。原生数值未在 Python 重建。重算必须使用新目录，不能覆盖接受检查点前缀。

16212 格原首步失败保存在 `finest-dt005/`；`run-finest-resolved.py` 显式使用最小 dt=1e-7 s、24 次重试，先真实积分至 1e-6 s，再以该接受状态续至 25 s。现已完成；累计量拼接 `finest-start/` 与 `finest-resolved/`，之前失败的计算成本另列，不能当作零成本初始化。闭合后的三档分析与新增原始文件索引在 `finest-completion/`，早先两档分析及原哈希索引保留。改变分段终点造成实际时间网格变化，不能宣称与一次直达轨迹相同。

`cartmesh2d_transport_cli` 在自适应步未接受时输出 `controller`：`reason`、`attempts`、`maximumRetries`、`timeStep`、`minimumTimeStep`、`errorRatio`、`courant`、`maximumCourant`；非有限误差/CFL 写 JSON null，原尝试 CSV 保留。值来自既有控制器回调，仅记录信息，不参与接受。`core/thermal.js::thermalFailureMessage` 形成原因和操作建议，`thermal-job.js` 根据实际保留的正时间状态描述续算能力；renderer 完整显示原因、折行而非省略，日志保留文件路径。旧 JSON 缺少详情时仍可报告已接受时间，不猜测原因。

真实 App 曲壁入口：

```sh
./desktop/node_modules/.bin/electron desktop --smoke=circle --geometry="$PWD/tests/fixtures/thermal_cylinder64.xy" --target-cells=5000 --auto-padding=3 --small-alpha=.1 --thermal-adaptive=true --thermal-case=external --thermal-diffusivity=.002 --thermal-events=true --thermal-flux-correction=bounded --flow-min-dt=.0000001 --flow-max-retries=24 --out="$PWD/outputs/advection-app-rerun" --shot="$PWD/outputs/advection-app-rerun.png" --export="$PWD/outputs/advection-app-rerun.zip"
```

该入口使用正常几何导入、renderer、IPC 和原生 CLI，几何预设先选外流圆柱，再读实际夹具。温度对流保持 App 默认 upwind，上下边界绝热，源项/壁面事件与 25 s 原生代表例不同。加入 `--thermal-start-failure-check=true --thermal-start-failure-only=true` 可只验收预期失败：将首步和下限设为 1e-6 s、无重试，要求明确原因、文本不被裁剪且没有接受状态；返回成功仅指失败流程验收通过。普通组合则继续完成 .5/.8/预算失败/1.3 s 全流程。

`outputs/thermal-advection/reopen.command.json` 保存真实 ZIP 重开至 1.4 s 的调用；`compare-cli.cjs` 用相同输入状态与控制复算并比较 7 类原生文件，检查原项目/ZIP 和旧显示以及未来事件。`startup-defect/` 用原生整步与两半步输出定位首步缺陷来源，面积加权 RMS 只作诊断，未替换实际最大值接受范数。版本、测试、实际范围与未完成项统一见当前状态及 `artifacts/current/native-thermal-advection.json`。

## 云端联合温度时间控制

本方向基于 `991c9b85a677daa3771af486ec3c1f5939e61388`，独立分支 `codex/thermal-stability-cloud`；不与其他开发分支或 main 自动合并。原生接口在 `ThermalFlow2D.hpp/.cpp`，物理模型为二维恒物性不可压层流与单向被动温度。

`advanceControlledThermalFlow2D` 只从不可变的最后接受状态生成候选。先求新时刻流动，再用同一守恒面通量求温度；任何一方未收敛均不发布状态。流动失败先把候选速度松弛减半、最低 .15，重试同一物理步；再按实际 CFL／半步缩小物理步。它改变内迭代路径，不放宽线性、动量、连续性或温度残差门。所有失败记录独立保存，接受历史不混入失败候选。取消在流动迭代、温度校正与候选结束检查；不会中断正在执行的一次线性求解。

控制器历史是接受状态的一部分：保存下一建议步长、成功的数值松弛及精确数值控制签名。时间误差按 BE 的 dt² 比例给出 .9 安全系数、至多两倍增长，面通量 CFL 与目标／事件仍裁剪每个候选；内迭代少于 20 次时逐步恢复用户松弛。任何失败或取消均不提交新历史。改变最大／最小步长、CFL、误差开关／尺度／容差或用户速度松弛会显式重置建议；目标时间不属于签名，仍实际积分最终余量。同目标、同控制、在同一接受步中断时，历史与物理场均可逐字节对照。改变分段目标会新增一次截断步，需按实际时间离散差异解释，不能宣称同轨迹。

CFL 及内迭代收敛只控制可接受性。可选 BE 整步／两半步估计 `2*|T_full-T_halfhalf|`，归一化分母为 `temperature-atol + time-rtol*temperature-scale`（K）；速度分别用 m/s 的绝对容差和参考速度。温升尺度必须是目标温度变化量，例如 1 K，不能用 300 K 背景温度。默认绝对容差为 .001 K、.0001 m/s，相对容差 .01；这是局部时间缺陷预算，非全程／空间误差承诺。默认不启用误差估计；启用每次候选多两个联合求解。当前仍接受整步一阶 BE，不外推场；曲壁曾主要受流动内迭代限制；细网格 App 又暴露重复从最大步试起的成本，现以检查点化建议步长解决。保留一阶，依据真实敏感性与代价决定后续二阶需求。

瞬态动量插值使用空间动量对角元，再加入面上的时间质量，避免插值总对角元的倒数造成稳态依赖 dt。记 `d_i=V_i/A_spatial,i`、`d_f=I(d_i)`、`theta=dt/(dt+d_f)`，面压力响应为 `alpha*theta*d_f`；梯度响应使用 `alpha*theta*I(d_i grad(p_i))`。非均匀对角元还保留协方差项 `I((d_i-d_f)*(U*_i-U^n_i-(1-alpha)*(U^iter_i-U^n_i)))/(dt+d_f)`。旧接受面通量与内松弛缺陷一同保留；稳态时速度差项消失，恢复原稳态面方程，不能只给总对角元公式补一个旧通量项。空间对角元在加入时间质量前直接保存，避免极小 dt 下相减失精。背景可参阅 [Moguen 等 WCCM 2014 的动量插值讨论](https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p1362.pdf)；本实现的验收依据是原生稳态变步不变性及真实 CLI 流场对照，未引入新物理或放宽残差。

CLI：

```sh
build/cartmesh2d_transport_cli --mesh final.solver.cm2d --output outputs/run \
  --evolve-flow custom --flow-boundary flow.boundaries --boundary thermal.csv \
  --dt .1 --end-time 6 --min-dt 1e-6 --max-courant 1 --max-step-retries 10 \
  --time-error on --temperature-scale 1 --velocity-scale .2 --time-rtol .01 \
  --thermal-events events.csv
```

`flow.boundaries` 使用既有原生、网格绑定的命名边界格式；温度初始边界仍逐面 CSV。`events.csv` 表头为 `time,target,type,value,inflowValue`，目标可为已命名流动边界、`face:ID` 或 `source`；源项行用 `type=source`。同一时刻的多行一起形成完整快照，时间严格递增；同时间重复目标拒绝。例如 `1.37,lid,value,300,300`，或 `2.43,source,source,0.1,`。事件右连续，结束在事件的整步／半步使用左侧条件，下一步使用右侧条件；不跨越事件积分。

固定步且无事件的检查点仍为 v1，固定步事件记录为 v2；联合自适应使用 v3，包含控制器历史及完整事件表（可以为空）。v1/v2 续入自动控制时从无建议历史冷启动，不承诺与新 v3 历史相同的轨迹。读回核对全部源、边界、流入温度和事件时间及完整流动／网格；不给匿名 callback 当作重启身份。零时刻初值不写为物理检查点，旧零时刻记录亦拒绝物理续算。`.thermal.checkpoint` 用临时文件后原子替换；失败输入立即报错，保留原接受文件。

`.heat-history.csv` 保存每个接受步温度范围、热量积累、边界通量、源和储存项；`.boundary-heat-history.csv` 分命名边界记录流出／流入体积及对流热通量、向外扩散热流；`.attempt-history.csv` 记录 CFL、无量纲时间缺陷、数值松弛与具体拒绝阶段。热量指标乘 `rho*cp` 才成为单位深度物理热量；向外热流为正，负值给流体加热。

代表算例完全调用原生生成、质量验证与求解；Python 只编排 CLI、生成配置并读取真实输出，不重建独立方程／拓扑审计链：

```sh
python3 tools/thermal/workflow.py --case channel --output outputs/thermal/channel --end 6 --error --isolate-live-output
python3 tools/thermal/workflow.py --case cavity --output outputs/thermal/cavity --end 6 --error --isolate-live-output
python3 tools/thermal/workflow.py --case cylinder --output outputs/thermal/cylinder --level 5 --end 25 --diffusivity .05 --nu .1 --error --isolate-live-output
ctest --test-dir build -R 'cartmesh2d_(thermal_flow$|thermal_flow_cli$|thermal_control_cli$|scalar_transport$)' --output-on-failure
```

同类粗细网格用 `--level 4/5`，实际流速用 `--speed`，热扩散率用 `--diffusivity`，最大步长用 `--dt`；尖锐启动研究驱动默认给 18 次候选重试，不改变原生产品默认值。变流速时同时更新真实边界和误差参考速度。对命名 `channel-coarse/channel-fine/...` 等目录中的冻结原生数据，可用 `python3 tools/thermal/summarize.py --input outputs/thermal/qualified --output outputs/thermal/summary.json --plot outputs/thermal/response.png` 读取、核对记录完整性并画图；该脚本只汇总原生预算与真实场，不重建离散方程。其命名输入矩阵和实际命令保存在原始证据目录，当前关键版本／资格只见 CURRENT_STATE。

三例分别覆盖加热通道、顶盖封闭腔体热／冷却、64 边真实曲壁圆柱（外域减固体）。通道／圆柱显式使用压力开口与规定流入温度，保留压力出口回流拒绝的失败输入。`.command.json`、原始场、接受／拒绝记录及检查点在 outputs；小证据进入 artifacts，完整字段压缩保留。原生解析衰减、核心守恒测试保留，旧 Python 重复审计不作为本方向验证后端。

Linux 无显示服务的真实 App 操作验收（需先准备 `desktop/runtime` 和 Electron 依赖）：
```bash
desktop/node_modules/.bin/electron --no-sandbox --ozone-platform=headless --disable-gpu desktop --smoke=rectangle --target-cells=1000 --thermal-adaptive=true --out=outputs/thermal/app --shot=outputs/thermal/app/adaptive.png
```
此入口实际运行联合控制并检查变步长续算、预算失败保留和恢复；软件渲染截图不代表 macOS/Windows 验证。运行期间不要归档或替换活动输出文件。云工作区曾出现活动路径被替换但原生继续写旧 inode 的现象；具体外部写入者未确立。原生现在逐次 flush 核对流位置与公开路径长度，发现替换／截断明确失败；CLI 编排还核对接受行数、最终时刻与原生 summary。`--isolate-live-output` 在同一云执行器的临时 outputs 目录计算，关闭后按 SHA256 导入指定 outputs，保留临时源供失败恢复；不复制活动文件。App 的 `--out`／`--shot` 同样可指定临时 outputs，结束后再导入。已有接受检查点的输出前缀拒绝复用，续算必须用新的 `--output`，原 `--restart` 文件保持不可变。此保护不能防止进程退出后的外部文件损坏，因此归档仍须逐文件读回核对。

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

接受状态只在真实状态、物性和积分时钟检查后更新。会话保留 BDF 历史；步间取消可继续同一会话，内部试算取消或失败保留最后接受状态，但须用新会话重建历史。试算取消引起的后端中止单独报告，重复读取终止会话保留取消状态及原因，不复用已中断的历史。验证输出的 `integration.canceled` 与实际完成时钟分开记录。残差评估和接受步都有成本预算，不把预算耗尽记为完成。火焰验证驱动保存 `evaluation-progress.jsonl`；在输出目录创建 `cancel.request` 可有序停止并保存最后接受状态、检查点和失败原因。直接杀进程不提供同样的保存保证。

长时间验证可用 `--sample-every N` 将初始状态、每 N 个接受步及末态写入 `samples.jsonl`。采样不改变时间步或 BDF 历史，使用单独的物性和残差上下文；额外诊断调用在 `sampleResidualEvaluations` 中计数，仍包含在完整进程耗时内。默认关闭采样；开启后，`--max-samples` 默认 256 条，用于限制诊断输出数量，不是物理精度门。达到保存上限或诊断出错会在接受状态退出并尝试写入最后接受场和检查点，原失败原因进入结果；部分轨迹明确标记 `sampling_complete=false`，不冒充完整计算。

读取器用同一套独立时钟、EOS、共享面和质量/元素/能量预算检查逐帧核对。`initial` 表示本段起点：冷启动来自参考初值，续算来自真实接受检查点，两者通过 `restart` 元数据区分。`sample-audits.json` 保留每帧诊断；历史图使用真实接受时刻。固定本段初始温度范围中点的等温位置通过单元中心剖面线性插值得到，多个交点全部记录；该位置标记不等于已验收火焰速度。464 格、0.1 微秒采样开关对照的全部物理场相等，步长日志和检查点逐字节相同。故意改坏的时钟和质量预算被拒绝；两条样本上限测试在第 20 个接受步退出，保留的场和预算与正常运行的同一步相同。

试算终止后可能增加 RHS、矩阵及非线性迭代计数，而未接受任何新场。`audit_sample_counters` 仅在未完成、有失败原因、末帧物理状态严格等于最后接受状态时，区分这种终止快照与较早的常规采样；普通帧仍须与接受步 RHS 计数相等。终止计数只允许明确列出的单调诊断变化，积分控制、接受阶数计数和接受误差统计保持；每帧继续接受原 `audit_fields` 检查。报告另列 `terminal_attempt_statistics`，不把额外试算计入物理时间。运行 `build/chemistry-env/bin/python tools/verification/test_reacting_sample_counters.py` 检查计数语义；原生取消/新会话恢复测试归入 `cartmesh2d_reacting_implicit`。

`outputs/combustion-foundation/implicit-cancel-audit-001/` 绑定冻结程序、审计器及源码哈希，保留真实试算取消、零时长文件读回、续算和正常计算对照。`reaudit_legacy.py` 在单独目录重审旧记录，并检查原物理审计函数 AST 相同及实际时钟/预算/状态损坏拒绝；`restart.py` 分段审计后按真实边界、化学和约束积分合并预算。原 `acoustic-boundary-001` 的取消记录及 `acoustic-newton-001` 的报告元数据失败保留，补充已存在的夹具间距只修复审计输入，不重算或改写原场。冷态脉冲使用完整氢气机理、多组分/Soret、温变物性和黏性；脉冲是初值，不能用它代替燃烧传播或边界反射验证。相同 10 微秒终点的三/六次牛顿预算及六次预算的容差减半对照只描述该案例，默认值不改。完整文件、控制及观察范围见[证据索引](../artifacts/current/native-reacting-cancellation.json)。

完整 100 微秒输出保留在 `outputs/combustion-foundation/implicit-evolution-100us-001`，进程实际使用 `implicit-sampling-002/probe` 和 `implicit-sampling-003/verify_reacting_flame.py` 的固定版本；不要用当前构建冒充旧运行来源。`analyze_completed_evolution.py` 从已核对哈希的实际场和样本读取不同等温位置，列出末段残差变化；温度平移诊断只在初始剖面的真实重叠区间作线性插值，不生成接受状态。不同温度标记的速度不能直接作为火焰速度资格。`residual-profile-001` 保留临时计时补丁、编译命令、调用数量及原/计时版本逐字节对照；它只定位成本，不构成提速证据。详细来源统一由[演化证据](../artifacts/current/native-reacting-evolution.json)索引。

`analyze_reacting_flame_history.py` 将已完成且独立审计通过的二维条带历史转为可复现的传播诊断。它校对报告、原始场、样本、参考与夹具哈希，逐条读取 JSONL，核对首末状态、实际时钟、总步数与样本数；零推进、部分采样和来源不符均拒绝。兼容旧零起点报告及真实续算的累计时钟；没有 `endTime` 或累计步数字段的旧记录必须证明确为零时钟/零步冷启动，不能据缺字段将续算重置为零。逐条读取避免把所有完整场同时保存在内存。

新细参考与动量相容初值的完整长算位于 `outputs/combustion-foundation/refined-evolution-100us-001/`。`run.py` 绑定原生源码、冻结程序、后端、机理和输入哈希，使用目录内冻结的原审计器；`history-complete/` 使用已验证的 `history-tool-source-002/` 分析完整历史，不能用中途前缀代替终点。`write_complete_artifact.py` 保存实际完整轨迹及与旧初值长算的观察性比较；两次冻结二进制与输入均有区别，不能从中推导单因素因果或通用性能。完整命令、文件哈希、实际控制和逐帧报告见[证据索引](../artifacts/current/native-reacting-refined-evolution.json)。

`refined-evolution-200us-001/` 从上述真实检查点追加 100 微秒，保持冻结程序、机理、网格、边界、三项误差容差和最大步长；只重建 BDF 历史与初始尺度。`analyze_complete.py` 核对启动提交的源码 Git blob、二进制、后端和实际文件哈希，再以原历史工具分析新段；不得因当前工作区已有后续取消修复而将其冒充运行来源。合并诊断保留两段实际审计，严格核对接口场相等，重复的 100 微秒样本只在绘图时间序列中去重；物理预算合并真实边界、化学和代数约束积分。全程温漂和平移形状差异重新相对零时刻计算，本段温漂另列。`moving_frame.py` 使用末段实际中部等温标记速度和非均匀二阶单元中心差分计算 `RHS+c*dU/dx`；这是后处理，不是原生有限体积算子，也不能将两个 L1 范数相减后宣称因果贡献比例。[续算证据](../artifacts/current/native-reacting-refined-continuation.json)

开放边界的脉冲对照在 `acoustic-boundary-002/` 生成相同重叠网格和初始单元守恒均值，仅延长出口域。初始 400 K 氢气/空气混合物叠加 100 Pa 高斯压力脉冲，按固定初始组分的等熵关系及压力积分构造向外速度，以八点 Gauss 积分形成单元平均；原生推进仍保留完整反应、温变物性、黏性、多组分和 Soret。`analyze.py` 对相同物理终点的共同单元直接比较压力、温度、速度和组分，并单列低幅值线性声学方向诊断 `0.5*(Δp ± rho0*c0*Δu)`。`acoustic-boundary-003/` 将间距从 250 减至 125 微米，`acoustic-boundary-004/` 在原网格上将三项时间容差共同减半；每档分别运行原域和延长域，最终比较同一 45 微秒实际场。`acoustic-boundary-004/write_artifact.py` 核对同网格两档的初值完全相同、完整控制及逐帧来源，再比较两条域长压差曲线。域长不同会改变全局加权时间误差范数，控制相同不保证误差相同；不预设边界反射合格阈值。真实模型、控制和结果见[边界对照证据](../artifacts/current/native-reacting-acoustic-boundary.json)。

标记温度（K）、上游测点（m）、平移对齐温度（K）和末段窗口（s）均须显式指定。等温交点与等温平台全部保留；窗口内不足三个样本，或任一时刻没有唯一交点，就不估计该标记速度。用真实保存时刻拟合固定坐标速度，同时列出端点割线、相邻时刻速度范围与拟合位置偏差。测点速度按实际时间梯形积分，再减去标记速度，得到相对运动诊断；另列测点温度/组分相对入口的偏离和是否位于标记上游，不自动将它判定为未燃气或成熟燃烧速度。平移仅插值初始行均温度的真实重叠区，并记录被排除列数、宽度比例和跨行温差；原场、残差和检查点不变。

`test_reacting_flame_history.py` 用解析平移剖面检查位置、秒制速度、变形及多交点/平台，也检查 JSONL 缺帧、增帧和时钟倒退。几何/速度/温度数值检查允许 `64*epsilon` 乘对应构造尺度；速度还除以实际时间跨度，只用于解析恒等式的运算舍入，不增加 CFD 精度门。测试为秒级，真实历史复核只读既有场。示例：

```sh
build/chemistry-env/bin/python tools/verification/test_reacting_flame_history.py
build/chemistry-env/bin/python tools/verification/analyze_reacting_flame_history.py \
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
build/chemistry-env/bin/python tools/verification/verify_reacting_flame.py \
  --probe build/cartmesh2d_reacting_implicit_flame_probe \
  --reference outputs/combustion-foundation/flame-reference-new \
  --output outputs/combustion-foundation/implicit-flame-new --duration 1e-6 --grid 0 --jobs 1 --reflect-species 1
# 10 微米档计算至 100 微秒并保存接受状态；时间步仍由原误差控制决定。
build/chemistry-env/bin/python tools/verification/verify_reacting_flame.py \
  --probe build/cartmesh2d_reacting_implicit_flame_probe \
  --reference outputs/combustion-foundation/flame-reference-new \
  --output outputs/combustion-foundation/implicit-long-new \
  --duration 1e-4 --grid 2 --jobs 1 --reflect-species 1 --sample-every 200
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

参考解可用 `--refine-slopes` 增加细化层级，`--maximum-reference-points` 只控制参考网格的计算预算。前者是 Cantera 按剖面变化增点的无量纲控制，曲率控制为其两倍，不是燃烧精度验收门；BVP 原求解容差保留。`--fixed-grid-reference` 与 `--fixed-grid-index` 从哈希验证过的旧夹具复用精确网格边、行数和高度，重新积分最终参考守恒剖面，用来隔离参考初值敏感性。拒绝坏网格、哈希不符及超出参考域的映射。`reference.json` 的 `profile_file` 明确指向最终剖面；原生审计器及动量初值准备器均读取它，旧参考仍兼容。

`cartmesh2d_reacting_flame_reference_probe` 是只链接 Cantera 的独立稳态 BVP 验证程序，不链接原生瞬态有限体积核心。可选 `--closure-probe` 在固定参考网格上重新求解完整机理、多组分和 Soret，用 `sum(Y)=1` 替代一种丰富组分在内部点的冗余方程；原组分方程全部照常计算，原残差另存于最终剖面，不能因质量恒等式满足就忽略组分平衡。能量、边界和原求解容差保留，最终初值导入仍使用上述原额度。原未约束参考及失败结果保留；该选项不归一化求解后的参考剖面，不修改原生接受状态或默认设置。示例：

```sh
cmake --build build --target cartmesh2d_reacting_flame_reference_probe -j 2
build/chemistry-env/bin/python tools/verification/prepare_reacting_flame.py \
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
build/chemistry-env/bin/python tools/verification/verify_reacting_flame.py \
  --probe build/cartmesh2d_reacting_implicit_flame_probe \
  --reference outputs/combustion-foundation/flame-reference-001 \
  --output outputs/combustion-foundation/initial-evaluation-new --grid 2 --evaluate-initial
build/chemistry-env/bin/python tools/verification/prepare_reacting_flame_momentum.py \
  --reference outputs/combustion-foundation/flame-reference-001 \
  --source-run outputs/combustion-foundation/initial-evaluation-new --grid 2 \
  --output outputs/combustion-foundation/momentum-reference-new
# 准备器只导出所选夹具，在新参考目录中的索引为 0；推进仍使用原完整方程。
build/chemistry-env/bin/python tools/verification/verify_reacting_flame.py \
  --probe build/cartmesh2d_reacting_implicit_flame_probe \
  --reference outputs/combustion-foundation/momentum-reference-new \
  --output outputs/combustion-foundation/momentum-flame-new \
  --duration 1e-5 --grid 0 --jobs 1 --reflect-species 1 --sample-every 50 --max-samples 128
```

`outputs/combustion-foundation/refined-state-refinement-001/` 保存新细参考的时间/空间对照。`map_reference.py` 从已哈希验证的 3446 点 `prepared-reference.npz` 积分到原三档精确网格，使用正式准备器的保守映射函数，不重算 BVP；`run_grid.py` 经三档初值审计后生成动量相容输入，粗/中档实际推进，细档夹具与既有完整运行逐字节相同后复用。因合并夹具清单与旧单档清单的元数据不同，`analyze_refinement.py` 分别核对各自清单、同一 BVP/映射源、入口行、检查点机理/物理前缀及完整边界数值和控制，再调用原网格比较器的守恒重叠/EOS 函数；没有放宽通用比较器要求相同清单的检查。`run_time.py` 使用冻结的上一版审计器在相同细档上收紧时间容差，两份版本的原物理审计函数 AST 与正时长审计结果均核对一致。完整命令、源码/原始场哈希和初值模式拒绝检查见[证据索引](../artifacts/current/native-reacting-refined-state.json)；复现使用新目录，不覆盖已完成结果。

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
