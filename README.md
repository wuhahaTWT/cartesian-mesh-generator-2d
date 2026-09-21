# CartMesh2D

原生二维笛卡尔网格生成器，提供 **macOS、Windows、Linux 桌面应用**和命令行工具。导入二维轮廓或 PNG/JPG 图片，设置网格数量或相对尺寸，生成真实 Cut-cell / 贴体边界层网格，并导出 OpenFOAM case。

当前桌面源码版本 **0.4.24**，定位为有明确支持范围的工程验证版。固定圆柱、翼型、喷管已生成十万级网格；任意复杂几何的鲁棒性、网格无关性与 CFD 精度仍需逐例验证。

自研 **二维稳态／非定常不可压层流**求解器直接使用最终多边形网格：SIMPLE / Rhie–Chow、共享面通量、完整黏性应力、后向欧拉，以及可选的面限制线性重构。App已接入内流/外流、逐面命名速度入口、压力出口与轴对齐静压开口、静止/移动壁面、旋转同心圆环、固定/自动时间步、初始局部涡、取消续算、受力/CFL监测、保存/读取流动工况和真实场导出。工况文件绑定已有网格，尚不包含完整几何/网格项目；湍流仍未在App开放，可压尚未实现。新增求解功能仅在macOS实测，下方旧三平台下载包不包含本轮更新。[构建、运行与能力范围](docs/DEVELOPMENT_CN.md#自研求解器入口)

当前App已实际跑通喷管、旋转圆环、圆柱、厚翼型和双柱；每份结果分别检查网格和离散方程，不把“收敛”当成精度合格。0.4.23的自动模式可在粗候选明确质量失败后有界细化重试，厚翼型从失败候选得到4,728格有效网格并完成计算，数量偏差仍明确显示。[真实厚翼型App](artifacts/current/native-airfoil-app.png) · [验证范围与剩余问题](docs/CURRENT_STATE_CN.md#当前连续目标macos-上可用的二维不可压-cfd-核心)

0.4.24 加入**两端静压驱动内流**：App直接设置p/rho，由求解器算流量，保留非定常和续算。解析矩形通道三档速度误差呈二阶下降，4,096格流量误差.1953%；实际2,784格喷管完成稳态、启动、续算和导出审核。仅支持水平/竖直开口，普通压力出口仍拒绝回流。[真实场与解析对照](artifacts/current/native-pressure-openings.png) · [实际App](artifacts/current/native-pressure-opening-app.png)

0.4.12 打通本地喷管的稳态／非定常、续算、被动热输运和 App 显示／导出。628、2,380、9,336 格的真实喷管均完成计算及独立守恒复核；最后两档压降仍相差2.34%，不宣称网格无关或物理精度认证。[真实喷管流场](artifacts/current/native-duct.png) · [代码、工况与验证范围](docs/DEVELOPMENT_CN.md#原生层流求解)

开发分支新增**被动标量／恒物性温度输运**：支持冻结已有流场或流动—温度同步推进、定值及绝热/指定通量边界。同步模式联合保存温度、流场和同一时钟，中途失败保留上个完整时刻，已验证中断续算一致性。桌面可设置热边界、同步推进、取消/联合续算、查看温度与监测量，并导出真实 PNG/VTK/CSV。当前是单向恒物性输运，不含浮力；新功能仅在 macOS 实测。[真实桌面温度场](artifacts/current/desktop-thermal.png) · [运行与边界单位](docs/DEVELOPMENT_CN.md#守恒标量与恒物性热输运) · [验证范围](docs/CURRENT_STATE_CN.md#当前改进同步流动与温度推进联合续算)

0.4.8 修复了部分壁面切割单元压力重构放大误差的问题：默认约5千格圆柱已完成流动与温度推进、取消续算和导出。0.4.9 另加入可选的法向回流压力出口，已用双向解析流加密和真实圆柱验证，仍属试验功能。[真实网格与修复说明](artifacts/current/pressure-stencil.png) · [回流验证](docs/CURRENT_STATE_CN.md#当前改进可选压力出口回流-049)

开发分支的 **CLI/核心已加入非定常不可压层流**：真实物理时间、一阶后向欧拉、CFL诊断和断点续算；每个时间步都检查动量与质量守恒，未收敛时保留最后接受的状态。解析涡衰减、规则/扭曲网格及启动流动用于验证。桌面0.4.4已接入时间步、动能/受力/CFL监测和重启；取消或时间步失败会保留上次完整结果及最后接受的续算状态。这不代表湍流或涡脱落精度已获验证。[运行方式](docs/DEVELOPMENT_CN.md#非定常层流与断点续算) · [真实网格及误差图](artifacts/current/native-flow-transient.png) · [实际桌面与续算证据](artifacts/current/desktop-transient.json)

长时间解析衰减验证提示：**残差收敛不等于时间精度足够**。900格Taylor–Green在t=2时，dt从.04减至.005，速度相对误差由15.75%降至1.75%，观测时间阶约1；真实Cut-cell圆柱也进行三档时间步比较。原始超时记录与适用范围保留。[物理时间验证](docs/CURRENT_STATE_CN.md#当前改进较长物理时间与步长敏感性)

0.4.10 减少多重网格装配的临时存储。同一500,472格通道，峰值RSS约934→804 MiB，数值输出逐位相同；总耗时37.47→37.03秒，尚不足以宣称明显整体加速。[本轮实测](docs/CURRENT_STATE_CN.md#当前改进多重网格装配存储-0410)

0.4.11 复用当前迭代的动量矩阵与梯度，圆柱、回流及五十万格通道前后结果逐字节一致；同一五十万格总耗时36.73→35.16秒、RSS约921→794 MiB。十万格热涡现已完成两步和独立守恒审核，但**最细档速度误差回升，三档精度趋势仍未通过**。默认松弛未改，.9/1.0失败记录保留。[实测、失败及适用范围](docs/CURRENT_STATE_CN.md#当前改进动量装配复用与十万格热涡-0411)

进一步对照表明，原速度误差回升受到迭代容差明显影响：仅收紧停止条件，十万格速度误差由3.30e-7降到6.33e-8，三档重新下降；但单独流场两步总耗时约370秒，且未重做同步温度研究。默认值与原失败记录保留。[迭代敏感性证据](docs/CURRENT_STATE_CN.md#当前改进分离迭代误差与网格细化)

开发分支进一步复用压力多重网格分组，并跳过压力修正中结果完全不变的重复遍。指定十万格涡流单步约185.1→178.3→169.4秒，后一改动的数值输出逐字节相同，独立守恒审核通过；桌面包未更新。[实测范围与限制](docs/CURRENT_STATE_CN.md#当前改进消除压力修正的相同重复步骤)

紧容差的**同步热输运三档对照现已通过**：1万/4万/10.24万格相对同一后向欧拉参考的温度和速度误差均下降，同步与独立载流完全一致。十万格完整两步约336秒，采用另行声明的420秒资源预算，仍未满足旧300秒预算；连续时间速度误差仍受固定dt影响，不代表通用工程精度。[实际热场、误差与边界](docs/CURRENT_STATE_CN.md#当前改进紧容差同步热输运三档对照)

SST标量求解新增可选ILU(0)：同一高Re平板，2,560格约11.18→6.64秒，5,120格约38.88→18.62秒；最终场独立方程审核及预设算法差异检查通过，默认Jacobi保留。12,800格固定100次约48.48→9.44秒，但完整求解仍在90秒预算内未通过，不能把吞吐提速当成万级SST验收。[实际网格、速度场与证据](docs/CURRENT_STATE_CN.md#当前改进sst标量ilu0与有界失败诊断)

开发分支已有**实验性稳态SST-2003m耦合核心**：k/omega更新后的湍黏度进入SIMPLE及完整对称应力，最终同时核对质量、动量和两条湍流方程。196/900/3,844格通道通过独立离散方程审核；交错更新在900格同容差串行对照中从28.15秒降到3.39秒，不能外推为十万格速度。**标准湍流工况、近壁精度和规模性能尚未验收，App仍未开放湍流。** [真实场、耗时与支持边界](docs/CURRENT_STATE_CN.md#当前改进实验性稳态sst与simple耦合)

实验性平板核心已能区分上游滑移段与真实板面，并可选择开放压力远场或对称顶部。196/900格四例通过独立方程审核，已输出Cf/y+诊断；当前近壁网格仍不足以证明摩阻精度，未复现高Re标准例。[真实边界对照与限制](docs/CURRENT_STATE_CN.md#当前改进平板混合边界与开放远场)

近壁法向加密现已完成**2,048格诊断平板**并通过独立最终方程审核。新增可选的SST校正迭代顺序，保留原停止条件；同一1,024格输入从63.83秒降至6.37秒，结果差异远小于预设算法一致性上限。旧方式仍为默认。这是Re=500核心诊断和单次配对计时，尚未验收高Re湍流、网格无关性或十万格SST性能。[真实网格、耗时与验证范围](docs/CURRENT_STATE_CN.md#当前改进有界次数的sst校正与同精度提速)

SST诊断支持实际速度、黏度、入口湍流量、前缘及矩形域。板长Re=5千、5万和1千万的五份流场现均通过独立方程审核：原薄网格失败通过精确矩形中心计算修复，旧错误场仍会被拒绝，验收容差未放宽。**高Re近壁分辨率及摩阻仍未验收**。[真实网格、几何修复与验证边界](docs/CURRENT_STATE_CN.md#当前改进薄矩形中心精度与独立残差定义)

预设空间黏度已接入动量方程、完整对称应力、壁面力与v3续算。三档网格上的速度、压力和壁面黏性牵引误差均接近二阶下降，续算一致性与独立守恒审核通过；旧恒黏度输出不变。这是核心/CLI基础能力，SST与桌面材料配置尚未完成。[实测与限制](docs/CURRENT_STATE_CN.md#当前改进空间黏度完整应力与续算)

新增预设空间扩散系数基础：正方形与真实斜边Cut-cell各三档网格上，限制线性输运误差随细化接近二阶下降，逐面本构与逐格守恒独立审核通过。当前限冻结载流的核心/CLI；同步变物性与完整SST仍未完成。[真实网格、误差及适用范围](docs/CURRENT_STATE_CN.md#当前改进守恒空间扩散场)

固定一万格网格补做时间步细化：dt从.005减半至.00125，连续时间参考的速度RMS由6.65e-5降到1.63e-5，温度扰动RMS由2.08e-6降到6.27e-7；同步/独立载流及守恒审核通过。这是短时恒物性案例的证据，未扩展为通用工程认证。[时间步对照图与适用范围](docs/CURRENT_STATE_CN.md#当前改进固定网格时间步细化)

0.4.6 新增可选的**聚合多重网格压力求解**，默认仍为 IC0。同条件实测：112,368格通道约22.7→10.6秒，500,472格约105.0→38.0秒；五十万格峰值RSS约723→876 MiB，独立解析/守恒读回通过。这是指定层流网格上的试验结果，并不保证任意复杂流动同样提速。[真实网格与性能图](artifacts/current/native-flow-aggregation.png) · [验证与限制](docs/CURRENT_STATE_CN.md#当前改进试验性压力多重网格046)

方腔精度新增公开高精度参考复核：当前核心在196→900→3,844→15,876格上，对Marchi等2009中心线速度的RMSE由顶盖速度的1.90%降至0.0217%，观测阶约2。求解方程未改，**旧Ghia表的单调误差检查仍保留失败**；两份参考分别报告，不能用一个绿色结果替代全部验证。[真实网格及对照图](artifacts/current/native-flow-cavity-reference.png) · [来源、误差与限制](docs/CURRENT_STATE_CN.md#当前改进方腔公开参考与四档空间细化)

0.4.5 优化残差范数计算，在同一网格、物性和收敛条件下，**112,368格通道约36.5→22.3秒，500,472格约164.5→100.4秒**；数值输出逐字节相同，两档均重新通过独立解析/守恒审核。各档只有一组计时，不代表复杂绕流或湍流性能保证。[真实网格与对照图](artifacts/current/native-flow-pressure-performance.png) · [耗时、内存与能力边界](docs/CURRENT_STATE_CN.md#当前改进残差范数提速与新核心规模复核045)

0.4.1 加入共享面压力和独立动量审核，桌面可切换对流格式。通道、方腔和圆柱共14次试算的单例守恒检查通过，但**限制重构的方腔三档细化趋势尚未通过**，收紧迭代容差也未解决；保持迎风默认值，不宣称通用二阶精度。[实际流场与未通过项](artifacts/current/native-flow-convection.png) · [核查证据](artifacts/current/native-flow-convection.json)

开发分支另有完整流动方程制造解验证：规则/扭曲网格、两种格式共12例，限制重构速度观测阶约1.83→1.92。它验证光滑解析强迫流动，**不消除上述方腔未通过项**；CLI验证入口不属于新增桌面物理工况。[真实网格与误差图](artifacts/current/native-flow-manufactured.png) · [复现方式](docs/DEVELOPMENT_CN.md#完整流动方程制造解)

0.4.2 修正未知壁面压力的外推，避免把压力强制解释为零法向梯度。新增非零壁面压力梯度的六例验证，规则网格压力误差随196→900→3,844格由 .00763→.00190→.000478 下降，扭曲网格也接近二阶。14例物理回归的单例守恒检查通过，**原方腔 Ghia 细化门仍未通过**；另用独立有限差分法对照，记录差异而不替换原验收门。[压力场与误差图](artifacts/current/native-flow-pressure-boundary.png) · [完整证据及桌面实测范围](docs/CURRENT_STATE_CN.md#当前改进压力边界外推与独立方腔对照)

0.4.3 将完整对称黏性应力接入动量方程，默认受力复用同一共享面通量。通道三档壁面黏性总力趋近解析值，最细15,372格的相对误差约 **0.030%**；光滑制造解的壁面牵引接近二阶收敛。这是指定层流案例的证据，原方腔细化未通过项继续保留。macOS实际App已验证生成、计算、导出及独立读回。[真实网格与受力误差图](artifacts/current/native-flow-symmetric-stress.png) · [验证范围](docs/CURRENT_STATE_CN.md#当前改进守恒对称应力与壁面受力)

方腔核查另修正了后处理取样的线性场偏差，保留同一流场的新旧对照；**原细化趋势未通过项仍然保留**。此修改改善误差测量，不代表求解器精度突然提高。[取样与独立对照](docs/CURRENT_STATE_CN.md#当前改进方腔取样误差分离)

[原生层流的实际界面与数值证据](docs/CURRENT_STATE_CN.md#原生二维稳态层流与桌面-040本地交付)

[下载与启动](#下载与启动) · [快速使用](#快速使用) · [网格展示](#网格展示) · [桌面使用指南](docs/DESKTOP_APP_CN.md) · [当前状态与已知问题](docs/CURRENT_STATE_CN.md)

## 项目阶段与分支

**网格生成器 0.3.0 阶段已完成，现进入自研 CFD 开发。** “阶段完成”指当前工程验证版的交付范围，不代表任意几何和工业精度问题已全部解决。

- [网格里程碑 mesher-v0.3.0](https://github.com/wuhahaTWT/cartesian-mesh-generator-2d/releases/tag/mesher-v0.3.0)：固定在加入自研 CFD 前的 `691c97e`，保留当时可复现的源码与证据。
- [网格维护 codex/mesh-maintenance](https://github.com/wuhahaTWT/cartesian-mesh-generator-2d/tree/codex/mesh-maintenance)：从该节点继续修复、完善网格生成器。
- [CFD 开发 codex/cfd-development](https://github.com/wuhahaTWT/cartesian-mesh-generator-2d/tree/codex/cfd-development)：包含原生求解器及桌面计算流程；`main` 保留已验收的集成版本。

当前目标先在 macOS 打通本地内流／外流模型与 App，补齐通用边界设置、不可压精度和非定常可靠性，通过前置验收后再实现独立可压核心。已有试验 SST 不替代层流的基本验收；[当前路线与真实验证范围](docs/CURRENT_STATE_CN.md#当前连续目标macos-上可用的二维不可压-cfd-核心)持续更新。

标签保存固定快照，分支可以继续发展。以后回头改网格不影响这次里程碑；维护分支的修复通过测试后再合并到 CFD 分支和 `main`，不会自动互相覆盖。

## 下载与启动

打开[已通过三平台验收的构建](https://github.com/wuhahaTWT/cartesian-mesh-generator-2d/actions/runs/34689931686)，在页面底部 **Artifacts** 下载对应运行包。后续版本可在 [desktop-platforms 工作流](https://github.com/wuhahaTWT/cartesian-mesh-generator-2d/actions/workflows/desktop-platforms.yml)中查看；选择成功运行及所需平台的附件。GitHub Actions 附件下载需要登录，附件也有保留期限，过期后可选择较新构建或从源码构建。

| 平台 | 下载附件 | 解压后打开 | 验证范围 |
| --- | --- | --- | --- |
| macOS · Apple Silicon | `CartMesh2D-mac-arm64` | `CartMesh2D.app` | macOS 14 arm64 |
| Windows · x64 | `CartMesh2D-win-x64` | `CartMesh2D.exe` | Windows Server 2022；面向 Windows 10 1903+ / 11，用户桌面版本仍需实机确认 |
| Linux · x64 | `CartMesh2D-linux-x64` | `./cartmesh2d-desktop` | Ubuntu 22.04 图形桌面环境；其他发行版尚未逐一验证 |

下载的 Actions 附件是外层 ZIP；解压后，再展开其中的应用 ZIP 或 `.tar.gz`。请保留整个应用目录，不能只复制可执行文件。

运行包自带 Electron、二维生成器、11 个样例和中文界面字体，使用时不需要安装 Node.js、Python 或编译器。Linux 需要系统 GTK/NSS/音频等桌面库；当前包不是 Alpine/musl 通用包。应用尚未签名或公证。OpenFOAM / Fluent 需在自己的仿真环境另行准备。

已有本地构建时，可使用仓库根目录的 `打开CartMesh2D.command`、`start-windows.cmd` 或 `start-linux.sh`。完整要求见[桌面使用指南](docs/DESKTOP_APP_CN.md)。

## 快速使用

1. **选轮廓**：选择内置样例，或打开 XY、CSV、TXT、SVG、DXF、PNG/JPG。闭合轮廓默认是固体，网格生成在外部；管道、喷管内流选择“内部为流体域”。
2. **定疏密**：选择约 5 千、1.5 万、5 万、15 万、30 万，或自定义数量，再按工况设置计算域留白。点击“生成预览”。
3. **计算（可选）**：生成后，在“原生稳态层流”选择匹配工况，填写运动黏度、速度和迭代上限，启动求解；完成后可切换速度/压力色图。
4. **检查并导出**：查看实际数量和质量结果，点击“导出结果包”。解压后先看 `mesh-preview.png` 和 `README_CN.md`，再把网格送入目标 CFD 软件。

图片会在本地自动提取二维轮廓。先检查红线，再填写**轮廓整体的实际宽度**及单位，确认后生成。优先使用清晰剪影、正视图或透明底图片；复杂背景可调整阈值或点选主体。识别不上传图片，也不会恢复三维形状或猜测物理尺度。

生成后可点“回到开始”释放预览；上一份结果仍可导出。侧栏支持按钮或 Cmd/Ctrl+B 收展，窗口可调整大小；“现代简约”和“科研终端”两种主题可切换，网格配色独立设置。

## 密度与手动调整

| 想调整什么 | 操作 |
| --- | --- |
| 整体更多或更少格子 | 更换目标数量 |
| 物体附近更细 | 减小壁面尺寸 `h/Lref` |
| 远处更细 | 减小背景尺寸 `h/Lref` |
| 细网格覆盖更远 | 增加加密带宽 |
| 复现固定高密案例 | 内置圆柱、翼型、喷管选择“载入已验证高密样例参数” |

`Lref` 是参考长度，可显式填写圆柱直径、翼型弦长等工程长度。例如 `h/Lref=0.005` 表示目标格子尺寸为参考长度的 0.5%。导入坐标按源单位换算为米，不自动归一到 1 m。请求尺寸与最终实现尺寸可能不同，实测值见结果包中的 `.resolution.json`。

自动选参最多尝试 3 组，目标范围为 ±30%；未达到时保留最接近目标的成功结果并提示。可点“用本次实际参数继续微调”。计算域留白、内外流语义、生成方法和质量门不会为了凑数量而自动改变。

数量表示计算规模，不表示精度等级。自定义上限 50 万、30 万档仍属规模试验；贴体首层初值尚未由 y+ 或流场误差推导。高密样例会替换参考长度、域留白等参数，载入后需检查；通用自动选参不能保证每种几何成功。

## 网格与导出

- **原生二维**：四叉树控制疏密，物面处保留真实流体 Cut-cell 多边形；OpenFOAM 输出将二维网格挤出为单层体单元，前后边界为 `empty`。
- **贴体边界层（Beta）**：壁面四边形层与笛卡尔余域共形连接，支持外流和单环内流；带孔内流目前使用纯 Cut-cell。
- **结果 ZIP**：包含 CM2D、VTK、JSON 质量/参数报告、输入轮廓、OpenFOAM case，以及白底全景＋局部放大的 `mesh-preview.png` 和中文文件说明。图片导入还保留原图、确认轮廓图与标定记录。

## 网格展示

下图来自桌面 App 实际导出的固定喷管混合网格，**148,287 个最终单元**；全景与壁面局部均读取真实网格。[导出证据](artifacts/current/start-export-evidence.json) · [质量检查记录](artifacts/current/target-quality-and-sidebar.json)

![固定喷管混合网格，148287 个最终单元](artifacts/current/start-export-preview.png)

更多固定案例保存在[展示素材](展示素材/展示说明.md)。以下六图是质量修复前的历史快照，数字对应各自原始网格，不作为当前版本全量重跑结果；旧图中的 Q1 标注也属于历史记录。

| 几何 | 纯 Cut-cell | 混合网格 |
| --- | --- | --- |
| 固定 32 段圆形折线 | [153,208 cells](展示素材/01_circle_pure_153208.png) | [164,970 cells](展示素材/04_circle_hybrid_164970.png) |
| NACA 2412 | [133,810 cells](展示素材/02_naca_pure_133810.png) | [140,305 cells](展示素材/05_naca_hybrid_140305.png) |
| 喷管内流 | [141,488 cells](展示素材/03_nozzle_pure_141488.png) | [148,375 cells](展示素材/06_nozzle_hybrid_148375.png) |

在其他 Markdown 页面复用展示图，可使用公开仓库的 raw URL：

```md
![Fixed 32-segment circle](https://raw.githubusercontent.com/wuhahaTWT/cartesian-mesh-generator-2d/main/%E5%B1%95%E7%A4%BA%E7%B4%A0%E6%9D%90/01_circle_pure_153208.png)
```

需要固定版本时，将 URL 中的 `main` 换成该图片所在版本的 commit SHA。

## 验证与已知边界

源码 `132adad` 的[跨平台验收](https://github.com/wuhahaTWT/cartesian-mesh-generator-2d/actions/runs/34689931686)在三个系统分别通过 **99 项原生测试、97 项前端测试**，并启动打包 App，完成 PNG、JPG、混合圆三例的生成、ZIP 导出与独立读回。含中文和空格的路径、随包字体和导出图片也已验证。[证据](artifacts/current/desktop-platforms.json) · [Windows 实际界面](artifacts/current/desktop-platforms-win.png) · [Linux 实际界面](artifacts/current/desktop-platforms-linux.png)

这轮跨平台验证没有运行外部 OpenFOAM `checkMesh` 或 CFD 求解。内部拓扑、Solver 质量、独立文件读回、目标求解器原生检查和 CFD 精度是不同层面的证据；Q1 已从当前产品取消，其他质量门继续保留。已有固定案例的外部检查与 CFD 结果见[当前状态](docs/CURRENT_STATE_CN.md)。

另已用既有 **140,305 格 NACA2412** 实际运行 OpenFOAM 2606：恒定来流 15 m/s、Re=1e6、0° 攻角、SST 稳态 RANS，完成 3500 步并输出真实[速度/压力场](artifacts/current/airfoil-rans-Re1e6-fields.png)、[流线](artifacts/current/airfoil-rans-Re1e6-streamlines.png)及[残差历史](artifacts/current/airfoil-rans-Re1e6-histories.png)。压力残差 2.22e-6，尚未达到设定 1e-6；标准 checkMesh 通过，扩展检查仍有 2720 个 concave 单元。域边界偏近、近壁 y+ 约 15，未证明气动力精度或网格无关性。[完整证据](artifacts/current/airfoil-rans-Re1e6.json)

目前仍需注意：图片识别成功不代表网格一定生成成功；复杂背景和透视需要人工处理；SVG 的 transform/viewBox 等语义尚未完整覆盖；窄缝、尖角、多环与贴体层终止仍有支持边界。成功生成和单元数增加都不能代替网格无关性验证。

## 从源码构建

核心使用 C++20，桌面使用 Electron。准备 Node.js 22、CMake 3.20+、Python 3 和本机 C++20 编译器：macOS 使用 Xcode Command Line Tools 的系统 clang；Windows 使用 Visual Studio 2022 C++ 桌面工作负载及 Windows SDK；Linux 使用 GCC 11+ 或兼容 Clang。

在仓库根目录执行：

```sh
npm ci --prefix desktop
npm --prefix desktop run build:native
cmake --build build --config Release --parallel 2
ctest --test-dir build -C Release --output-on-failure
npm --prefix desktop test
```

然后按当前系统执行**其中一条**：

```sh
npm --prefix desktop run pack:mac
npm --prefix desktop run pack:win
npm --prefix desktop run pack:linux
```

运行包输出到 `desktop/dist/`。应在目标系统本机构建，打包前会检查生成器的系统与架构，不能将 macOS 二进制直接放入 Windows/Linux 包。界面开发使用 `npm --prefix desktop start`；C++ 修改后需重新构建原生工具，更新 App 中的生成器。

命令行入口为 `cartmesh2d_cli`、`cartmesh2d_hybrid_cli`，参数见 `--help`；Windows 文件名带 `.exe`。完整 CLI 示例、模块导航与验证命令见[开发导航](docs/DEVELOPMENT_CN.md)。

- [桌面使用指南](docs/DESKTOP_APP_CN.md)：图片导入、参数、预览与导出。
- [当前状态](docs/CURRENT_STATE_CN.md)：唯一当前状态入口，记录实测数据与未解决问题。
- [开发导航](docs/DEVELOPMENT_CN.md)：源码、构建与验证工具。
- [开发规则](AGENTS.md)：物理定义、维护和交付要求。
