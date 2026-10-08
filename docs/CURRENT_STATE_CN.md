# 当前状态与接手重点

更新：2026-10-08。这里是版本、能力、验证与待办的唯一状态入口；操作见[桌面使用](DESKTOP_APP_CN.md)，算法及复现见[开发指南](DEVELOPMENT_CN.md)。

## 当前研究与接手入口

当前先把**不可压 CFD 的精度、稳定性和效率基础做扎实**，再用可靠的真实壁面 CFD 纠正下一轮流道拓扑演化。用户要求继续实质研发，并要求将工作迁往云端、保持 Astra Max；云环境设置会话已创建，截图确认 Astra Max，尚未确认环境发布和研发任务启动。详细燃烧继续冻结。

新增入口 `tools/benchmarks/laminar.py`：方腔 Re100/400/1000、DFG Re20、旋转圆环、Poiseuille，三档网格、两种对流格式。首轮 36 个矩阵条目中 24 个收敛、8 个未通过网格质量门、4 个求解失败。失败计入结果，不降低质量门。新结果及可异机重跑的输入见 [基线](../artifacts/current/laminar-foundation/baseline.json)，复现见[开发指南](DEVELOPMENT_CN.md#不可压基础基准与云端复现)。

- 方腔 Re100 的 limited-linear/Marchi RMSE 为 `0.004080 → 0.0009002 → 0.0002163`，末两档观测阶 `2.082 / 2.011`；最细 15,876 格需 8,590 轮。Ghia 指标单列，Re400 原表可疑条目未静默删去。
- 旧 6,208 格圆环默认控制再次发散；启用现有 Anderson、保持速度松弛 0.6 可在 980 轮收敛，但压力归一化 L2 仍 `0.02775795`。完全隐式非正交压力算子的试验未解决发散，未保留为产品功能；[失败补丁](../artifacts/current/laminar-foundation/rejected-implicit-pressure.patch)及基线 JSON 内的日志已上传，原始场保留在 `outputs/laminar-foundation/`。
- DFG 使用等尺寸笛卡尔背景可生成通过现有质量门的 3,063 格网格，但默认 SIMPLE 及松弛 0.2 + Anderson 均发散。下一步检查动量预测中显式对称黏性应力、壁面压力重构及耦合稳定性；尚无 DFG 精度通过结论。
- 同网格 OpenFOAM 准备/读取脚本已加入；首次运行被 Docker 重启中断，当前没有完成的外部对照结果。

本轮已实际完成：系统 Clang Release 网格/流动 CLI 构建、完整 36 条目基线、额外圆环/DFG 失败定位、Python 语法检查和 IC(0) Poiseuille 小例运行。求解器产品源码恢复为 `6ec4edf`，没有把未奏效的试验算作修复；未重新运行完整 CTest、前端或跨平台 CI。云端应先读取本节和已提交输入接续修复。

待完成：通用精度/稳定性修复及失败例回归；补齐曲壁同网格外部比较和三档收敛；分离边界预设、动量、压力修正和停止条件，加入可复现算例输入；测量工程停止/SIMPLEC 在同精度下的完整成本；重算代表拓扑设计与基线。尚未宣称基础任务完成。

| 工作 | 接手时的状态与入口 |
| --- | --- |
| T01：保真度差距与排序 | 预定实验及交付完成；主索引为 [fidelity-completion.json](../artifacts/current/native-topology-fidelity-completion.json)。文献优化流道、MAC 空间收敛及可信排序误差仍有研究限制，原“18 对反转”不是已证实结论。 |
| T02：真实壁面反馈 | 在线控制器已实现并跑通；代码 `tools/optimization/closed_loop_topology.py`，主索引 [feedback.json](../artifacts/current/native-topology-feedback.json)，[方法与运行命令](DEVELOPMENT_CN.md#t02真实壁面反馈驱动材料更新)。Re=0 总收益主要来自结构提案，校正增量约 0.0256%；Re=50 校正组落后于同控制对照，尚未证明方法优势。 |
| 当前需要解决的问题 | 校正怎样产生更有效的新材料方向和连接变化，以及何时值得花一次真实复算；现有失败和对照可直接用于方法改进。T03 换热应用、T04 比赛材料尚未开展。 |
| 本地工作区 | 根目录、`codex/cfd-development`；接手基线为 `6ec4edf`，T01/T02 和清理已在其历史中。新增基准继续在开发线，接手先读当前 `git status`。 |

建议先读本节及[拓扑研究详细状态](#独立研究流体拓扑优化原型)，再读开发指南中的[文件与结果导航](DEVELOPMENT_CN.md#本地研究文件导航)。Claude 原任务单、过程日志和早期尝试保留；它们记录当时的指导，当前状态以本文件为准。原始场、失败、运行环境和现有构建均保留。

## 验证冗余清理

按本轮确认的范围，删除重复防护、过度兜底和独立验证脚本，保留核心算法约束与失败判定。独立 Python 拓扑/方程重建、专项扫描及其专用测试和 CI 入口已移除；共用的数据读取、几何计算、初值映射和输入准备移到 `tools/flow/`。原生求解链入口改为 `tools/optimization/native_flow.py`，不再生成独立审计通过标记。T01 的任意磁盘停机阈值、逐候选目录扫描，以及将文件/编程错误吞为候选失败的兜底已删除；桌面 runtime 检查集中到最终打包入口。

原生几何、拓扑、Solver 质量、物理约束、NaN 和收敛失败判定继续保留。原始设计、真实场、失败记录、参考包未删除。下文既有审计和测试数字属于清理前的对应版本，不能当作当前测试数量。

修改前的源码及未提交差异保存在 `outputs/verification-cleanup-6gxhncca/`，包含 `source-before.zip`、`preexisting.patch`、删除清单和本轮运行日志。当前清理在开发线，尚未提交。

本轮在 macOS 用系统 Clang 完成原生构建及桌面 runtime 构建；相关原生测试 22 项、拓扑研究测试 42 项、前端测试 191 项通过。实际提取几何重新生成网格，1 步预算正确返回未收敛，同网格续算累计 155 步收敛；一次短闭环和新报告绘图已运行。此处仅确认清理后的调用链，不新增方法优势或物理精度结论；未运行全量 CTest、跨平台 CI、打包 App 或燃烧研究全流程。

## 当前交付

**0.4.41 已完成三条 CFD 开发线的集成。** 可压层流、纯笛卡尔浸入边界和流道拓扑优化都在 `main` 中，原提交历史完整保留。后续统一在 `codex/cfd-development` 开发；旧的三条功能分支、临时集成分支和网格维护分支已结束，不再作为并行入口。

| 范围 | 当前可用能力 |
| --- | --- |
| 网格与不可压层流 | Cut-cell、共形边界层、稳态/非定常、命名边界、守恒受力、检查点、被动恒物性热输运 |
| 完整笛卡尔背景 | 均匀/自适应几何网格；`solver_ready=false`，不能直接用于既有 CFD |
| 可压层流开发版 | CLI 和桌面试验入口；二阶对流、恒物性黏性/导热、可选二次壁面重构 |
| 纯笛卡尔 CFD 实验版 | 独立 `cartmesh2d_immersed_cli`；静止物体、均匀格、周期压差通道，未接入桌面 |
| 流道拓扑优化研究工具 | `tools/optimization/`；设计材料分布、提取真实壁面、原生 CFD 复算，未接入桌面 |

软件集成不扩大物理精度资格。原默认数值控制、Solver 质量门、真实流体语义及固定标签 `mesher-v0.3.0` 保留；SST 等仍参与构建/测试的代码也保留。

## 持续目标：成熟燃烧模拟

**按用户 2026-10-02 决定冻结；代码、原始场、失败与续算材料保留，成熟资格未完成。** 目标仍是详细有限速率化学、多组分守恒输运、温变物性与流动/能量一致耦合，并针对实际燃烧工况完成验证和交付。反应器、一维火焰和短时二维案例是模块证据；首次实际验收工况待选，三维燃烧室须另建三维能力。

现有原生二维模块使用 Cantera 3.2，包含完整机理、多组分/Soret、共享面扩散与焓通量、可压对流、刚性化学及可选 BDF 隐式耦合；氢氧 10 组分/29 反应、甲烷 53 组分/325 反应用于模块检查。App 与湍流燃烧未接入。方法与复现见[详细化学](DEVELOPMENT_CN.md#详细燃烧化学基础)、[耦合推进](DEVELOPMENT_CN.md#详细反应流耦合推进)及[隐式模块](DEVELOPMENT_CN.md#全耦合隐式反应流开发模块)。

| 关键结果 | 证据与实际边界 |
| --- | --- |
| 物性、完整输运与短时耦合 | [化学](../artifacts/current/native-detailed-chemistry.json)、[扩散](../artifacts/current/native-reacting-diffusion.json)、[火焰/微量输运](../artifacts/current/native-reacting-flame.json)、[真实耦合](../artifacts/current/native-reacting-flow.json)。已修复查询顺序、微量组分重构和能量溢出；48 项局部输运对照及相关原生检查属于对应历史版本的实现证据。 |
| 取消、检查点与文件续算 | [取消恢复](../artifacts/current/native-reacting-cancellation.json)、[历史保存与续算](../artifacts/current/native-reacting-evolution.json)。160 格脉冲在 3.938 微秒取消并续至 10 微秒，对连续运行最大温差约 9.45e-6 K；BDF 历史需重建，不宣称轨迹逐字节等价。 |
| 可选连续热力学与查询复用 | [NASA7](../artifacts/current/native-reacting-thermo-continuity.json)、[查询成本](../artifacts/current/native-reacting-query-efficiency.json)、[隐式反射](../artifacts/current/native-reacting-implicit.json)。查询复用在 464 格指定案例中完整成本减少约 25.1%，场与接受轨迹一致；不外推通用提速，原默认输入和控制保留。 |
| 守恒参考 BVP | [参考通量](../artifacts/current/native-reacting-reference-flux.json)、[守恒形式](../artifacts/current/native-reacting-conservative-reference.json)。3446 点氢元素边界通量偏差 0.1697%→4.03e-7%；6891 点仍有网格依赖，13,781 点未收敛并按要求中止，不能计作最细档通过。 |
| 最新守恒参考初值的实际演化 | [10 微秒对照](../artifacts/current/native-reacting-conservative-evolution.json)。同一 890 格原生网格下最大温漂 0.308→0.448 K、组分残差 1.105%→1.171%，未改善；平移后形状差异 0.570→0.445 K。它是短时初值敏感性，未替换旧长期状态。 |
| 既有细参考长期演化 | [100 微秒](../artifacts/current/native-reacting-refined-evolution.json)、[200 微秒续算](../artifacts/current/native-reacting-refined-continuation.json)。累计 13,530 个接受步；200 微秒相对零时刻最大温漂约 0.929 K，尾部仍未同速、未形成稳定传播。此结果使用原 3446 点参考，不能与新守恒参考混作同一初值。 |
| 时间/网格与出口敏感性 | [细化](../artifacts/current/native-reacting-refinement.json)、[声学边界](../artifacts/current/native-reacting-acoustic-boundary.json)。原三档 464/610/890 格到 10 微秒的档间差异减小，但温漂未随细化消失；冷态小脉冲域长对照完成，尚不具备通用非反射边界资格。 |

恢复后须先解决最细参考的非线性收敛，分离传播形状变化与原生网格、时间和边界误差，再对选定实际工况验收质量/元素/能量、组分正性、机理适用性、实验或独立数值解、完整成本及交付。复杂壁面、强激波/接触面、长期稳定性、湍流/喷雾/辐射、实验、燃烧 App 和跨平台资格均未完成。上述独立读取器结果属于历史实验；已移除的验证链不重新作为当前入口。

## 验证与运行包

本次集成的验证代码为 `3a95510`，收尾为 `db45d0b`；后续文档整理不改变构建或运行源码。以下是实际完成的集成检查，原始研究精度结论另列在各功能小节。

- 本机系统 clang Release：完整原生 **164/164**、前端 **191/191**、流道研究 **21/21** 通过。可压原分支曾被磁盘资源门阻止的两项不可压测试，本次已实际通过。
- macOS、Windows、Linux 分别通过完整原生/前端回归、21 项研究测试及真实打包 App 检查。App 覆盖 PNG/JPG、共形边界层、两种背景网格，以及 228 格封闭腔体和 3568 格圆柱可压流；可压项检查控件、真实场、失败保留、取消续算、重复一致和 ZIP 独立审计。它们仍是短时功能检查。
- Linux 独立 CI 的 **15 例 OpenFOAM v2606 标准 checkMesh** 通过；扩展 checkMesh、本机外部 checkMesh 未运行。该结果不替代拓扑优化候选或浸入边界方法的物理验证。
- 集成修复了 GCC 数组赋值兼容和 Windows 可压求解器中文路径问题。Windows 慢测试仅增加运行预算，物理终点、精度断言和质量门保持；原失败记录保留。

[集成证据](../artifacts/current/cfd-integration.json) · [原生与 OpenFOAM CI](https://github.com/wuhahaTWT/cartesian-mesh-generator-2d/actions/runs/36587976831) · [三平台 App CI](https://github.com/wuhahaTWT/cartesian-mesh-generator-2d/actions/runs/36587976826) · [三平台研究 CI](https://github.com/wuhahaTWT/cartesian-mesh-generator-2d/actions/runs/36587976811)

本机运行包：`desktop/dist/CartMesh2D-0.4.41-arm64-mac.zip`；其他平台包在上述 App CI 附件中。包未签名。详细日志、场、截图及 SHA256 由集成证据索引到原集成工作区的 `outputs/integration/`。旧版[集成证据](../artifacts/current/main-integration.json)只作历史，不替代本次检查。

## 可压层流开发版

来源提交 `f7733a2`，已进入 main。

**解决什么、采用什么方法：** 让气体计算包含压缩、黏性摩擦和导热；在真实二维 Cut-cell 上使用守恒有限体积法，可选 HLLC/HLLE 二阶对流，并用局部二次拟合改善壁面热流和应力。默认仍为 Rusanov 一阶、线性壁面、μ=k=0、滑移，需要显式开启新选项。

**效果：** 32×32 扭曲网格的指定光滑导热例，壁面积分热流误差从约 **0.131% 降到 0.0338%**。Couette 冷启动及固定网格时间细化也检查了二阶行为、壁面功热平衡和续算；多项式温度恰好重现不代表任意流场具有机器精度。真实 App 已跑通二次壁面、黏性与导热的完整操作链。

**还差什么：** 当前是恒物性理想气体层流，局部二次重构不等于全局三阶。尚无湍流、温变物性、辐射、共轭传热或隐式推进；任意复杂曲壁的摩擦/换热精度、强激波与黏性相互作用及长期稳定性未验收。非正交扩散可能非单调，组合 CFL 与正性检查不是任意网格的稳定性证明；开边界采用零黏性牵引。

[壁面精度证据](../artifacts/current/native-euler-wall-accuracy.json) · [真实场与细化图](../artifacts/current/native-euler-wall-accuracy.png) · [App 控件](../artifacts/current/native-euler-wall-accuracy-app.png) · [算法与复现](DEVELOPMENT_CN.md#可压壁面精度与时间细化)。早期[对流](../artifacts/current/native-euler-development.json)、[导热](../artifacts/current/native-euler-conduction.json)、[黏性](../artifacts/current/native-euler-viscosity.json)记录保留，各自结果按原版本解释。

## 纯笛卡尔 CFD 实验入口

来源提交 `98ffe0d`，已进入 main。

**解决什么、采用什么方法：** 在不切碎方格的条件下处理物体壁面。独立 MAC 交错网格保留完整单元及明确标识的固体辅助未知量；可选整段壁面积分，将壁面力与压力一起求解，减少原有限阻力方法的漏流。旧 Brinkman 模式仍是默认，新模式需显式选择 `surface-penalty`。

**效果：** 指定 128×32 圆柱算例，最大壁速相对于无障碍参考速度 `Uref=1` 从 **9.10% 降至 0.0783%**；斜壁、矩形、凹角及平移对照也改善。真实折线上的稠密采样与方程残差分别检查，18 项直接相关检查及集成回归通过。

**还差什么：** 同一圆柱加密至 192×48，平均流量仍变化 **8.22%**，法向壁速也未单调下降。因此尚不能称为网格无关，近壁梯度、阻力和高 Re 壁面精度仍待验证；小壁速不代表整体准确。新模式 IC0 的非正主元失败已保留并明确拒绝，`auto` 用 Jacobi，macOS 可选 Cholesky；大规模效率未验收。当前只支持静止物体的周期通道，尚无桌面入口。

[壁面实测](../artifacts/current/native-cartesian-walls.json) · [方格与壁速对照](../artifacts/current/native-cartesian-walls.png) · [原型基线](../artifacts/current/native-cartesian-immersed.json) · [方法与复现](DEVELOPMENT_CN.md#纯笛卡尔浸入边界研究入口)。完整辅助格不是共形流体 `polyMesh`，外部 checkMesh 不适用。

## 独立研究：流体拓扑优化原型

**当前重点是让真实壁面复算改善下一轮材料方向和连接；T02 方法优势尚未证明。** 采用 Stokes/Navier–Stokes–Brinkman 材料设计与离散伴随，提取不平滑的真实轮廓，再以原生 Cut-cell 层流复算。研究工具未接入桌面，尚无原生 Cut-cell 伴随、物理最优或完整精度资格。

| 工作 | 实际结果与未完成边界 |
| --- | --- |
| T01 A/B：分析与优化 | 90 条有限差分曲线支持离散伴随实现；B 的 8 组达到原局部停止条件。局部停止不等于文献最优流道复现或全局最优；文献连接切换、部分原文参数仍有缺口。 |
| T01 C/D：真实复算和排序 | C 的 48 组、D 的 66 条记录均取得至少两档有效原生结果；共同阻力下另做 198 次材料评分。18 对反转仍只是条件 GCI 结论，物理已判定反转为 0；不能解释成排序总是一致。17 对原反转混用不同 αmax，MAC 空间收敛与原生渐近区均未确认。 |
| T01 复核发现的问题 | 原唯一同阻力候选的 MAC 细分变化仍约 22.9%/18.1%；一分量追加 level7 的变化达到旧 GCI 的 79.2 倍，另一分量累计 8000 轮未收敛。算术一致性不能给误差估计或可靠排序授予资格。 |
| T02 Re=0 主运行/对照 | 固定 level5 目标 95.51454→71.76930，整体改善约 24.860%；同控制关闭梯度校正的对照为 71.78771，校正增量只有约 0.0256%。主运行/对照分别 25/21 次原生评价、13.1/12.0 秒，主要收益来自结构提案。 |
| T02 Re=50 主运行/对照 | 校正组 95.50853→95.46205，关闭校正的对照为 95.44596，后者更好。合流提案被真实复算拒绝；弯管网格失败及弱收益也保留，未证明跨工况优势或通用提速。 |
| 既有连接原型 | 同面积/流量/物性下，level7 相对两直管观察压降降低约 46.5%；细候选 26,141 格、累计 7391 轮被接受，但细基准 4000 轮未收敛。优化亦达到预算上限，低阻力漏流、完整成对细化、网格无关性与原创性仍未证明。 |

T01 入口：[补齐索引](../artifacts/current/native-topology-fidelity-completion.json)、[审查](../artifacts/current/native-topology-fidelity-review.json)、[固定排名](../artifacts/current/native-topology-fidelity-fixed-ranking.png)、[真实场](../artifacts/current/native-topology-fidelity-native-fields.png)。T02 入口：[反馈数据](../artifacts/current/native-topology-feedback.json)、[实际演化图](../artifacts/current/native-topology-feedback.png)、[方法与命令](DEVELOPMENT_CN.md#t02真实壁面反馈驱动材料更新)。旧连接原型来源 `d13c175` 已进入 main，入口为[连接对照](../artifacts/current/native-fluid-topology-connectivity.json)和[细化](../artifacts/current/native-fluid-topology-candidate-refinement.png)。

T01/T02 当前改动尚未提交。对应历史版本的研究测试与原生/前端构建结果支持实现，独立验证链现已移除；本轮文档压缩不增加测试、App、平台、外部 checkMesh 或物理资格。下一步围绕校正增量与新连接生成能力推进，T03 换热、T04 比赛材料未开展，燃烧保持冻结。

## 已验证的规模和性能

以下均为指定输入的本机测量；网格生成、完整 CFD 成本及物理精度分别报告。

| 案例 | 实测与适用范围 | 证据 |
| --- | --- | --- |
| 均匀背景 512×512 | 262,144 格，内部 50,632 格保留；生成及 JSON/VTK 导出约 1.46 秒 | [背景网格](../artifacts/current/background-cartesian-grid.json) |
| 自适应背景 | 23,632 格，约 0.157 秒；覆盖、分类与 2:1 平衡通过 | [背景对照](../artifacts/current/background-matched-performance.json) |
| 102,017 格曲壁喷管 | 同精度完整 CFD 840.87→411.37 秒，约 2.04 倍；含共同网格生成约 1.92 倍；不是商业软件对照 | [梯度缓存](../artifacts/current/native-laminar-gradient-cache.json) |
| 102,400 格规则压力通道 | 冷启动 474.60→分阶段 195.06 秒，含粗解、映射和细解；不外推任意曲壁 | [初值对照](../artifacts/current/native-laminar-initial-guess.json) |
| 102,400 格非定常涡 | 三步至 t=.006，79.70→76.10 秒；只证明短时间对照 | [非定常记录](../artifacts/current/native-laminar-transient.json) |

曲壁粗细映射仍是研究流程，尚未接入 App 自动功能；50 万格曲壁只做过预算探测，未收敛。完整控制、计时和哈希见各证据，不只比较最快一段。

## 已知限制与继续工作

主线优先顺序仍是：**拓扑修复完成通用验收 → 圆环压力精度与默认稳定性 → 减少外迭代次数**。实验功能按上述范围继续，不自动扩展模型或改变默认值。

| 问题 | 当前结论与下一步 |
| --- | --- |
| 细圆环网格 | 轴附近身份修复消除了未分类边，但 7,233 格诊断网格仍有 68 项 Solver 失败；研究修复得到的 7,071 格通过质量和独立面积检查，尚需成为通用产品流程。[拓扑记录](../artifacts/current/native-annulus-axis-fixed.json) |
| 圆环压力 | 原细圆环压力 L2 约 2.776%；研究网格约 **0.564%**，仍高于既定 **0.5%**。调压力零点或收紧线性残差不能代替离散修复。[研究证据](../artifacts/current/native-annulus-repaired-research.json) |
| 默认稳定性 | strict、速度松弛 .6、4 次压力校正在困难圆环仍失败。研究阻尼曾使同条件 1,341→262 轮，但未进入发布求解器，也无十万格完整对照；优先减少外迭代并验证通用性 |
| 复杂几何和输入 | 窄缝、尖角、多环及层终止仍可能显式失败；SVG 的 transform/use 等完整语义未覆盖，图片识别不保证几何正确。不得靠删格、平滑原折线或放宽质量门掩盖问题 |

未取得通用湍流、复杂曲壁换热或长期非定常精度资格。拓扑、Solver 质量、外部 checkMesh、离散方程收敛和物理精度必须分别判断；Q1 已取消。

## 工作区与证据保留

- 本次目录治理删除 21 个旧 macOS 安装包及 21 个配套 blockmap，保留 0.4.40、0.4.41、现有 App/runtime 和全部研究数据；删除清单与目录计量见 `outputs/repository-tidy/storage-cleanup.json`。本轮未运行求解、全量测试或平台 CI，不增加数值/物理资格。
- 打包保留、CI 期限和任务产物收尾规则已接入现有入口，见[开发指南](DEVELOPMENT_CN.md#产物保留与清理)。当前状态的燃烧与拓扑逐轮叙述已压缩，旧段落在 `outputs/repository-tidy/retired-state-sections.json.gz`；原证据索引、失败记录及复现方法保留。
- 根目录是后续开发入口；常规构建用 `build/`，实验用忽略提交的 `outputs/`，运行包用 `desktop/dist/`。只维护 README、AGENTS 和 docs 内三份文档，旧过程查 Git。
- 当前另有 `~/.codex/worktrees/compressible-laminar/cartesian-mesh-generator-2d`（`codex/compressible-laminar`）和 `.claude/worktrees/suspicious-jepsen-e0adde`。它们的状态按各自检出解释；本轮未改动外部工作区、提交、标签或研究结果。
- 本次文件整理将源码目录的 Python 缓存和 Finder 杂项移入 `outputs/repository-tidy/recoverable-cache/`，旧 T01 启动/自动续跑脚本移入 `outputs/codex_tasks/retired_launchers/` 并取消可执行标记；恢复位置见 `outputs/repository-tidy/file-organization.json`。实验和原始结果未搬动，未释放或删除大批数值数据。
- 原有归档索引、实验历史和展示素材的未提交工作继续保留。历史六张大网格图片说明已并入[开发指南](DEVELOPMENT_CN.md#历史展示素材)，图片及来源证据保留。
- `outputs/cfd-demo-reference/` 是用户指定保留的旧 WebGPU demo，含源码、独立试玩、原 Git 历史和角色素材，常规清理不得删除。入口为包内 `README_CN.md` 和 `启动试玩.command`；复制核对及启动、推进、暂停已检查，原桌面 `cfd-rebuild` 副本已按用户确认删除。它没有并入原生求解器，也不新增物理精度资格。
