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

## 不可压基础基准与云端复现

`tools/benchmarks/laminar.py` 直接生成真实流体网格、显式物理边界，调用原生求解器并保存命令、二进制哈希、真实场、残差、误差及耗时。`summary.csv/json` 和 `convergence.png` 不将失败条目计入观测阶。速度误差按 `Uref`、运动学压力误差按 `Uref²` 归一化；圆环沿用 cell 0 压力基准。DFG 系数按平均入口速度 0.2、直径 0.1 定义。Ghia 原表与 Marchi Re100 分开统计。

Linux 云端从 `codex/cfd-development` 构建；基准仅需 C++20、CMake、Python、Matplotlib，层流目标不依赖 Cantera。macOS 仍使用系统 `/usr/bin/clang++`。压力预条件器在 macOS 默认 Cholesky，其他平台默认 IC(0)，可显式覆盖；跨平台耗时不直接比较。

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCARTMESH2D_BUILD_CHEMISTRY=OFF
cmake --build build --target cartmesh2d_cli cartmesh2d_flow_cli --parallel 4
python3 tools/benchmarks/laminar.py --output outputs/laminar-baseline --levels 5 6 7 --workers 2
```

`--geometry-segments N` 可固定各圆的原始多边形段数，再单独细化计算网格；默认仍是 `4*2**level`。记录保存几何哈希、实际流体面积、相对解析圆的最大弓高和面积差，输出目录拒绝混用不同段数。固定多边形系列的离散收敛与同时加密几何的圆形参考误差须分开解释；更少段数并不授予曲壁精度资格。Re400 的 Ghia `v(.9063,.5)=-.23827` 原表可疑点明确标记，但继续参与全部原表误差统计。

### SIMPLE 源码对照

`tools/benchmarks/simple_source_compare.py` 用于定位现有 SIMPLE 与成熟实现的具体差异：读取当前原生源码，在 `outputs/` 生成少量试验编译单元并链接现有二维库，运行真实 CLI；Python 只编排编译/运行/读取摘要，不重建方程。压力/速度改动仍限研究编译单元，产品默认离散、控制和验收门不变；产品已统一受力的共享面压力入口，并改进显式 Anderson 的回退。先构建当前 `cartmesh2d_flow_cli`，再运行；macOS 固定 `/usr/bin/clang++`，其他平台本次未验。

对照的原始源码入口如下，不把教程中的可配置离散选择说成 SIMPLE 算法的唯一规定：

| 检查点 | 外部原始源码 | 当前原生实现及结论 |
| --- | --- | --- |
| 压力松弛和速度更新顺序 | [Foundation 10 pEqn.H](https://github.com/OpenFOAM/OpenFOAM-10/blob/master/applications/solvers/incompressible/simpleFoam/pEqn.H#L44-L49)；[MOOSE correctVelocity](https://mooseframework.inl.gov/docs/doxygen/modules/LinearAssemblySegregatedSolve_8C_source.html#l00604) 先计算面通量，再松弛压力、刷新梯度和单元速度 | 原生用完整 `pc` 修速度、用 `alphaP*pc` 更新压力；是不同的分离迭代选择，单独对齐速度更新仍失败，不能凭差异宣布唯一 bug。 |
| 动量松弛 | [Foundation fvMatrix.C](https://github.com/OpenFOAM/OpenFOAM-10/blob/master/src/finiteVolume/fvMatrices/fvMatrix/fvMatrix.C) 与本机 OpenCFD 2606 `fvMatrix.C:1210–1214` 先保证对角占优，再松弛并补 RHS | 原生直接除速度松弛；单独加入占优处理，以及有界对流/占优/速度松弛组合，原 DFG 仍失败。 |
| 压力梯度、面值和壁面闭合 | 本机 OpenCFD 2606 的 `gaussGrad.C`、`correctedSnGrad.C`，以及 `tutorials/incompressible/simpleFoam/pitzDaily/system/fvSchemes` / `0/p` 中配置的 Gauss linear 和 zeroGradient | 原生非正交项用扩展 LS 梯度，压力力用含偏斜面值与边界外推的 Gauss 梯度。研究组合将压力面值设为内面线性插值/未知边界 owner 值，动量压力力和非正交项共用该 Gauss 梯度。它改变空间离散，必须另验受力和曲壁精度。 |
| HbyA、通量与非正交项 | Foundation `pEqn.H`、`constrainHbyA.C`；[boundedConvectionScheme.C](https://github.com/OpenFOAM/OpenFOAM-10/blob/master/src/finiteVolume/finiteVolume/convectionSchemes/boundedConvectionScheme/boundedConvectionScheme.C) | 成熟实现也使用滞后非正交校正；单步重构缺陷不是独立的代码错误证明。完整对照后恢复原生速度梯度、黏性转置项或面通量分别仍能让 DFG 收敛；恢复原压力组合则再次失败。 |

本机只读源码为 `/Users/Zhuanz/Desktop/mesh-reference-sources/OpenFOAM-v2606-source`，`META-INFO/api-info` 为 api=2606、patch=0；用到的文件 SHA 和外部链接保存在 `artifacts/current/simple-source-comparison.json`。Foundation/MOOSE 网页为访问时版本，未虚构固定提交。此次没有运行 OpenFOAM/MOOSE 程序；Docker CLI 存在但本机 daemon 未启动。原生对照只是独立实现选定的离散机制，没有复制或链接外部三维核心。

最小对照为两因素四组：`baseline`、只用松弛压力修速度 `velocity`、只改压力 Gauss 组合 `pressure`、同时改两项 `both`。其余保持原速度 LS 重构、对称应力、面通量松弛缺陷、对流、物性、边界、`.6/.25` 松弛、4次非正交校正及 strict 门。每次启动记录源码/库/二进制/输入 SHA、命令、实际退出码；失败日志保留，达到预算不记收敛。

```sh
cmake --build build --target cartmesh2d_flow_cli --parallel 4
mkdir -p outputs/simple-extensions/inputs
tar -xzf artifacts/current/laminar-foundation/failure-inputs.tar.gz -C outputs/simple-extensions/inputs
python3 tools/benchmarks/simple_source_compare.py --mesh outputs/simple-extensions/inputs/dfg20-l7.solver.cm2d --boundary outputs/simple-extensions/inputs/dfg20-l7.boundaries --nu .001 --speed .3 --output outputs/simple-source-comparison-reproduce
# 单独试两项组合，使用原 DFG 的 limited-linear：追加 --convection limited-linear --variants both，并换输出目录。
```

每次输出目录须新建，默认每个实际求解最多1500轮/45秒；这是研究成本预算，不是物理阈值。研究 CLI 对 Gauss 压力摘要使用`research-gauss-*`标签，不沿用产品压力重构标签。`--pressure-trace gauss-owner` 为原参考组合；`linear-extrapolated` 保留原LS内面/壁面trace但统一Gauss梯度，`skew-owner`只恢复内面偏斜，`linear-wall`只恢复未知壁面外推，后面三项在当前DFG均失败。`--coupling simplec`、`--steady-acceleration anderson --anderson-history 4`、`--linear-policy adaptive`可分别研究迭代控制；参数变化只通过原生控制入口。每个研究编译单元的压力面值同时进入动量、监测和受力导出，旧归档中独立后处理的Gauss受力已在当前状态及证据注明更正。

复现带回退的262轮DFG：在上例追加`--convection limited-linear --variants both --steady-acceleration anderson`；另加`--linear-policy adaptive`为268轮控制，最后仍由普通strict步认证。`--anderson-backtracking off`仅在生成的研究源码中禁用两次回退，可复现原4阶历史302轮对照，不增加产品开关。圆环换为原`annulus-l7`网格/边界、`--nu .1 --speed .5 --convection upwind --variants both --steady-acceleration anderson`；原45秒预算接近实测成本，较慢机器可显式增加`--seconds`，超时不得记收敛。以上只在原网格上证明指定数值路径，物理误差仍未通过。

当前成功、失败和精度限制只看[当前状态](CURRENT_STATE_CN.md#simple-源码对照诊断)。`simple-source-comparison-fields.tar.gz` 保存原四组及失败；`simple-source-followup-fields.tar.gz`保存更正受力、压力拆分失败、回退对照及细圆环代表场/原生日志/启动清单，供后续压力离散修复和同输入续算研究。原网格仍复用既有失败输入包。临时编译单元、对象和二进制保留在本机`outputs/`，可由脚本重建，不提交整份源码副本或新安装包。

### Cut-cell 小单元与尖锐单元处理文献

本节为用户要求的小单元/尖锐单元处理调研，检索截至2026-10-29；20篇原始研究中19篇为已发表论文，1篇2026预印本单列。只用原论文、作者/机构稿和成熟项目源码作技术依据，不把二手综述、搜索摘要或方法名当作已复现资格。表中“全文关键段”表示读取可获取原文的相关段落，不表示逐页精读或已实现；“摘要/预览”表示核对原出版方或作者机构摘要，不把付费全文说成已读。方法选择是本项目的工程判断，成功率、工作日和物理精度尚未验证。

**本仓库具体问题。** `SmallCell2D.cpp` 根据输入切割面积分数（默认0.10）选小单元，优先与非小邻居、较长共享边及较大面积的邻居合并；这一选择不直接优化最终重构矩阵或整片网格的形状。`SolverTopology2D.cpp:34` 的严格凸判定、`:606` 的普通源多边形划分，以及后续质量重划分，又可产生独立小单元。凸固体外侧的流体多边形具有物理反角，保留多条物面边的单元可能天然非凸，不能把“每个普通物理单元严格凸”无条件等同于有限体积的必要条件。当前几何/算子是否支持保留这些单元仍需检查，**不能简单删除凸划分后强行开算**。

现有同圆环的原生报告显示，GEOS6008格原始切割到共用后处理后的6208格，最小内角4.92°→1.94°、最大非正交34.54°→69.13°、最大长宽比24.36→60.78，当前质量失败160→0；同时最小邻接面积比.000942→.011698、最小插值权重.015527→.054366。它表明当前目标是在不同指标间换取过门，未证明整体适合CFD。原始格存在质量失败且未作为求解输入，整链快照也不能单独分配合并/凸分解/重划各自的责任。历史入口：`10fa0f7`（导出质量门）、`5c46daa`（严格凸分解）、`4b0af3f`（为外部共线/凹性检查加强划分）；后续确有防碎片修复，未闭合困难圆环稳定性/精度验收。[报告](../artifacts/current/cutcell-treatment-review.json)和[实际前后图](../artifacts/current/cutcell-postprocess-comparison.png)是后续网格修改的诊断输入，不是新验收门。

**需要分开的三个问题。** 一是显式通量更新中的 `dt/Vi` 导致时间步随小体积收紧；二是稳态压力/黏性算子的几何病态，包括中心距离、非正交、梯度模板与压力—速度响应；三是壁面剪切、压力和局部重构精度。给显式Euler/advection证明CFL稳定，不能证明现有稳态黏性SIMPLE稳定或压力精准。Johansen–Colella等结果也说明，小体积并非必然不可求：专门构造的边界通量和算子能降低对切割位置的敏感性。对本项目应先修已有网格处理缺陷，再匹配离散，不继续用coupled收敛代替网格验收。

| 路线 | 实际做法 | 本项目适用性、代价与限制 |
| --- | --- | --- |
| 几何合并 / agglomeration / reshaping | 小格与邻格组成较大的控制体，内部共享通量抵消，保留真实物面与流体面积。 | **首选，改造成本中等。** 已有拓扑/事务与合并基础；须从单邻居贪心升级为整片质量驱动选择，避免合并后又分回碎片。每个新格须过现有守恒/质量/几何门。 |
| cell-linking / 聚合自由度 | 积分几何保留小格，在方程层面将其未知量关联到宿主/聚合体，合并或约束残差。 | **备选，成本中高。** 将几何复杂度与独立未知量分开；不能只对速度做平均而保留原压力方程，动量、压力、通量、重构和映射需一致。 |
| 专用EB压力/扩散离散 | 使用物面法向、面形心及连通模板构建通量；配套多重网格、插值和离散边界处理。 | **与前两条并行设计，成本高。** Johansen–Colella针对特定椭圆问题；LS-STAG针对特定交错框架。结果不能直接移植为当前同位变量SIMPLE保证。 |
| Flux redistribution / 通量稳定化 | 先得到稳定更新，再把非守恒差额按权重守恒回分，或构造特殊面通量抵消。 | **后续时间推进研究，成本中高。** 主要文献处理双曲/显式更新；局部阶数、边界误差和稳态固定点须另查，不是纯网格修复。 |
| SRD / WSRD | 构造可重叠邻域，在邻域上加权平均与重构，再保守回分到原格。 | **保留研究选项，成本高。** WSRD已有不可压/低马赫及扩散应用，不能误说只适用Euler；但不自动修复当前压力矩阵，需同时处理散度约束、压力梯度和稳态保持。 |
| 局部隐式 / explicit–implicit | 小格及邻近过渡格隐式，远处显式，接口复用通量保持守恒。 | **暂不优先。** 能解除特定显式CFL约束；当前稳态求解本就含隐式子问题，仍需处理空间离散病态。混合接口还会引入精度问题。 |
| 近壁贴体层 + 笛卡尔余域 + 优质过渡区 | 近壁先生成法向层，困难接口移到层外，通过受约束局部网格连接。 | **黏性流重要路线，成本中高。** 本项目已有基础，但圆环目前仍需研究适配器；过渡区本身也必须控制形状、面积跃变及模板，不因近壁层好看即算通过。 |
| Ghost penalty / aggregated FEM | 在未拟合有限元空间中增加一致稳定项或约束小支撑自由度，控制切割引起的条件数恶化。 | **原理参考，成本很高。** 属于FEM/CutFEM或DG体系，不能把惩罚项、理论条件数结论直接搬进现有FVM，也不作为当前重写求解器的理由。 |

下面的文献按上述用途分组，保留英文全名和可核验的原始链接，方便继续阅读。

| 编号、论文 | 已查范围及可借鉴点 | 不能外推的范围 |
| --- | --- | --- |
| L01 Ye, Mittal, Udaykumar & Shyy (1999), [An Accurate Cartesian Grid Method for Viscous Incompressible Flows with Complex Immersed Boundaries](https://engineering.jhu.edu/fsag/wp-content/uploads/2014/06/CGM_JCP-99.pdf), JCP156:209–240, DOI10.1006/jcph.1999.6356 | **全文关键段**，尤其§2、Fig.2：截断格中的流体片段吸收到邻格，构造重新成形的控制体，配合专门插值/面通量。与二维同位变量黏性流直接相关。 | 吸收流体片段不等于按中心丢弃相交格；其投影与空间离散也不能当作本项目SIMPLE已经验证。 |
| L02 Kirkpatrick, Armfield & Kent (2003), [A representation of curved boundaries for the solution of the Navier–Stokes equations on a staggered three-dimensional Cartesian grid](https://www.sciencedirect.com/science/article/abs/pii/S002199910200013X), JCP184:1–36, DOI10.1016/S0021-9991(02)00013-X | **原出版方摘要/方法预览**：cell-linking，保留边界格几何并关联主从单元，配套压力力、剪切及面形心修正。 | 原文是三维交错网格；不能直接照搬主从映射到现有二维同位变量所有方程。 |
| L03 M.-H. Chung (2006), [Cartesian cut cell approach for simulating incompressible flows with rigid bodies of arbitrary shape](https://www.sciencedirect.com/science/article/abs/pii/S0045793005000721), Computers & Fluids35:607–623, DOI10.1016/j.compfluid.2005.04.005 | **原出版方摘要**：二维同位有限体积、cell-merging，维持全局质量/动量守恒，有静止及运动刚体例。 | 摘要不足以决定合并模板细节或给出可直接使用的阈值。 |
| L04 Johansen & Colella (1998), [A Cartesian Grid Embedded Boundary Method for Poisson’s Equation on Irregular Domains](https://www.sciencedirect.com/science/article/pii/S0021999198959654), JCP147:60–85, DOI10.1006/jcph.1998.5965；[作者机构稿](https://escholarship.org/uc/item/08z7t0bm) | **原出版方摘要、作者稿摘要**：特殊保守边界通量，使小体积不显著恶化其离散椭圆矩阵，配套多重网格。 | 原论文的变量位置、Dirichlet椭圆问题与本项目压力修正/壁面处理不同；不是一般任意多边形FVM定理。 |
| L05 Schwartz, Barad, Colella & Ligocki (2006), [A Cartesian grid embedded boundary method for the heat equation and Poisson’s equation in three dimensions](https://www.sciencedirect.com/science/article/abs/pii/S002199910500286X), JCP211:531–550, DOI10.1016/j.jcp.2005.06.010 | **原出版方摘要/成熟源码引用**：三维热方程、Poisson的EB面/壁通量构造，Basilisk也给出实现关联。 | 不单独授予不可压NS或SIMPLE稳定性。 |
| L06 (2019), [Accurate discretization of diffusion in the LS-STAG cut-cell method using diamond cell techniques](https://www.sciencedirect.com/science/article/abs/pii/S004579301930146X), Computers & Fluids189:34–45, DOI10.1016/j.compfluid.2019.05.005 | **原出版方摘要/方法预览**：非正交cut格的二点扩散近似存在精度问题；diamond-cell面法向梯度可处理单元/面数据和Dirichlet/Neumann条件。含Taylor–Couette等验证。 | 不是只把网格角度改好，也不是当前FVM全部黏性项的即插即用实现。 |
| L07 Colella, Graves, Keen & Modiano (2006), [A Cartesian grid embedded boundary method for hyperbolic conservation laws](https://crd.lbl.gov/assets/pubs_presos/AMCS/ANAG/A162.pdf), JCP211:347–366, DOI10.1016/j.jcp.2005.05.026 | **作者机构原文摘要及方法预览**：混合保守/稳定非保守更新，再把守恒差额重分配至邻域。 | 全局L1精度不等于边界逐点精度；不直接解决压力椭圆子问题。 |
| L08 Berger & Giuliani (2021), [A state redistribution algorithm for finite volume schemes on cut cell meshes](https://www.sciencedirect.com/science/article/abs/pii/S0021999120305945), JCP428:109820, DOI10.1016/j.jcp.2020.109820 | **原出版方摘要/方法概述**：SRD，用临时重叠邻域避免永久几何合并，保持守恒及线性再现。 | 论文主要建立显式双曲更新的办法，不能据此保证稳态黏性SIMPLE。 |
| L09 Giuliani et al. (2022), [A Weighted State Redistribution Algorithm for Embedded Boundary Grids](https://arxiv.org/html/2112.12360v2), JCP464:111305, DOI10.1016/j.jcp.2022.111305 | **全文关键段§2–5**：连续加权降低邻域敏感性与耗散，讨论梯度模板、对流/扩散及低马赫/不可压应用，给出AMReX-Hydro实现。 | 不能把邻域平均放在现有SIMPLE之后当滤波；散度约束及压力梯度的一致性必须重新设计。 |
| L10 Berger & Giuliani (2024), [A New Provably Stable Weighted State Redistribution Algorithm](https://arxiv.org/html/2308.16332v2), SISC46(5), DOI10.1137/23M1597484 | **全文关键段**：单调、TVD及GKS稳定分析有明确线性模型、模板和CFL前提；分析预合并益处。 | “provably stable”不是任意不可压黏性方程、任意模板的通用保证。 |
| L11 Barrio Sanchez et al. (2024), [A new re-redistribution scheme for weighted state redistribution with adaptive mesh refinement](https://www.sciencedirect.com/science/article/pii/S0021999124001281), JCP504:112879, DOI10.1016/j.jcp.2024.112879；[机构接受稿](https://www.osti.gov/servlets/purl/2328105) | **原出版方摘要/机构稿摘要**：cut格跨AMR粗细接口时，还要有re-redistribution才能保持守恒/稳定。 | 将单层WSRD直接搬到四叉树粗细交界是不完整实现。 |
| L12 Klein, Bates & Nikiforakis (2009), [Well-balanced compressible cut-cell simulation of atmospheric flow](https://doi.org/10.1098/rsta.2009.0174), Phil.Trans.A367:4559–4575 | **原论文摘要**：方向分裂的保守通量稳定化，可按规则格步长推进任意小cut格。 | 大气可压、重力平衡与当前稳态黏性压力耦合是不同任务。 |
| L13 Berger & Helzel (2012), [A Simplified h-box Method for Embedded Boundary Grids](https://epubs.siam.org/doi/10.1137/110829398), SISC34:A861–A888 | **原出版方摘要/预览**：扩大通量构造区域，利用抵消性质获得稳定、守恒的显式有限体积更新。 | 几何/通量模板复杂，不能仅以大控制体概念套用到黏性压力力。 |
| L14 May & Berger (2017), [An Explicit Implicit Scheme for Cut Cells in Embedded Boundary Meshes](https://link.springer.com/article/10.1007/s10915-016-0326-2), J.Sci.Comput.71:919–943；[2015作者报告](https://www.sam.math.ethz.ch/sam_reports/reports_final/reports2015/2015-44.pdf) | **期刊摘要及作者报告关键段**：局部隐式、远场显式，用flux bounding处理接口守恒。 | 作者报告与最终稿有版本差别，不能引用报告的所有定理当成最终稿；也不证明现有稳态系统条件良好。 |
| L15 May & Laakmann (2024), [Accuracy Analysis for Explicit-Implicit Finite Volume Schemes on Cut Cell Meshes](https://link.springer.com/article/10.1007/s42967-023-00345-y), CAMC6:2239–2264 | **开放全文关键段§1–2**：分析显式/隐式切换的精度损失，一维结论不自动扩展到二维。 | 解除稳定性限制仍可能付出局部精度代价。 |
| L16 Park, Jeong, Lee & Shin (2013), [Hybrid grid generation for viscous flow analysis](https://pure.johnshopkins.edu/en/publications/hybrid-grid-generation-for-viscous-flow-analysis/), IJNMF71:891–909, DOI10.1002/fld.3691 | **作者机构摘要**：近壁四边形/棱柱层、Cartesian余域、尺寸控制和优化的过渡区。 | 原文包含曲面拟合和三角/四面体过渡；本项目不得因此平滑原XY物面或照搬3D核心。它是路线依据，非当前产品资格。 |
| L17 Burman (2010), [Ghost penalty](https://comptes-rendus.academie-sciences.fr/mathematique/item/10.1016/j.crma.2010.10.006.pdf), C.R.Math.348:1217–1220 | **原文首页/摘要**：用一致惩罚控制切割位置造成的矩阵条件数敏感性。 | 理论建立在有限元空间和特定假设上，不能直接添加到当前FVM矩阵。 |
| L18 Badia, Verdugo & Martín (2018), [The aggregated unfitted finite element method for elliptic problems](https://www.sciencedirect.com/science/article/pii/S0045782518301476), CMAME336:533–553, DOI10.1016/j.cma.2018.03.022；[作者机构记录](https://research.vu.nl/en/publications/the-aggregated-unfitted-finite-element-method-for-elliptic-proble/) | **原摘要/机构记录**：通过聚合和空间约束避免小切割对条件数的破坏。 | 与FVM的几何合并、投影压力并非同一数学构造。 |
| L19 Seo & Mittal (2011), [A sharp-interface immersed boundary method with improved mass conservation and reduced spurious pressure oscillations](https://engineering.jhu.edu/fsag/wp-content/uploads/2014/06/JCP_Seo_Mittal_2011_2.pdf), JCP230:7347–7363, DOI10.1016/j.jcp.2011.06.003 | **作者原文摘要/引言**：在压力Poisson及速度修正中用虚拟cut格合并，处理运动边界几何守恒与压力振荡。 | 当前困难例是静态几何；不能把运动边界fresh/dead-cell原因当成本例已定位根因。 |
| L20 Karell (2026预印本), [Update-Magnitude State Redistribution (UM-SRD): A Shut-off Extension of Weighted SRD for Cut-Cell Methods](https://arxiv.org/html/2605.04863v2), arXiv2605.04863v2 | **预印本关键段**：提出更新趋零时关闭重分配，关注稳态保持；目前给出标量输运模型分析/数值例。 | 本次未见同行评审状态和黏性不可压NS资格，**不作为成熟默认方案**；作者提出的通用性主张也不替代独立核验。 |

成熟实现入口也分清用途：

- [AMReX-Hydro Redistribution](https://github.com/AMReX-Fluids/AMReX-Hydro/tree/development/Redistribution)：L09对应的通量/状态重分配实现；[IAMR](https://github.com/AMReX-Fluids/IAMR)提供不可压NS应用背景。二者与前轮“只取AMReX几何”的实验不是同一完整方法。现有本机AMReX快照并未自动获得这些流动/重分配能力。
- [Basilisk embed.h](https://basilisk.fr/src/embed.h)：分别处理EB面插值、连通模板、壁面法向梯度和显式小格通量溢出重分配，不能只摘最后一段。其源码还明确指出某些高阶边界梯度与近似投影组合会产生压力—速度反馈，说明“保留小格”需配套离散。
- [Basilisk couette.c](https://basilisk.fr/src/test/couette.c)：已有旋转圆筒EB收敛测试，但原测试设置 `stokes=true`、半径.25/.5且黏度1，**不是本项目含惯性、半径.5/1、nu=.1工况的既有验证结果**。可供下一轮基准设计参考，本次未运行。

实施优先级明确为：先让物理短壁面不再必然对应独立细小未知量，再做适配离散。短期比较“质量驱动多格聚合/局部重构”和“保留积分格、聚合方程未知量”两种设计；已有近壁混合网格作为黏性对照持续保留。候选目标需同时观察面积比例、角度、紧致度、非正交、有向法向中心距、梯度模板几何和邻接跃变；修改选择策略须保留物面、面积、共享面及事务确定性。不要只合面积不足.1的原格，也不要仅提高某个角度阈值后增加无解样本。

后续最小验收闭环固定原几何和物理边界，保存原最小失败例，记录切割/合并/划分各阶段来源；用原默认控制的自研SIMPLE和成熟SIMPLE回查，再以coupled定位同离散固定点，分别报告动量/连续性/线性残差、整体与局部压力/速度偏差。斜直壁和同一曲壁均需覆盖，网格分辨率与物面离散度分开控制。新阈值和测试规模必须先说明单位/归一化、依据、必要性与成本，不在本次文献调研中自动扩大CFD扫描或全量回归。本次没有变更生产数值代码、策略阈值、默认控制，也没有新求解收敛资格。

图和统计可由既有研究快照重建（原始质量报告来自现有原生实现，不添加独立Python方程审计）：

```sh
MPLCONFIGDIR=/tmp/cm2d-annulus-mpl python3 tools/benchmarks/plot_cutcell_treatment.py
```

输入为`outputs/annulus-cross/geos-canonical128.{raw,solver}.cm2d`及其原生质量/阶段报告；缺少本地快照时，用前节冻结几何/GEOS多边形和已有适配器重建。产物为`artifacts/current/cutcell-treatment-review.json`和`cutcell-postprocess-comparison.png`，消费者是后续网格修复、文献选型与同阶段对照。原始失败/场仍保留原归档与`outputs/`，不复制第三方代码/PDF或新增交接文档。

### 同一圆环交叉实验与解析参照

结果、数值资格和待办只看[当前状态](CURRENT_STATE_CN.md#cut-cell-网格与-simple-的交叉诊断)。`artifacts/current/annulus-cross.json` 记录五张实际网格、四个求解路径、另列的解析 trace 族、场偏差、逐次命令/退出码、源码/二进制和输入哈希。`annulus-cross-fields.tar.gz` 保存代表输入、接受场、原生日志/续算场、外部字典和成熟网格原始交换格式；归档内使用仓库相对路径，先解到独立 `outputs/annulus-replay/` 再选择输入，避免覆盖原始实验。原困难输入亦见既有 `laminar-foundation/failure-inputs.tar.gz`，直壁/旧圆环包仍保留。

原模型必须从冻结输入提取物面：**每圆256条真实直边**，984个导出面中的共线辅助点只在构造舍入预算内归并；不调用默认level7的512段圆重新生成。归档含 `original-walls.xy`、`original-physical-walls.xy`、`inner-wall.xy`。严格同壁面组使用后者，面积约2.35595793819958m²。GEOS构造端点只对齐同一笛卡尔格角身份，预算64 ulp，实际最大位移1.84e-16m；这是构造舍入，不是删除可解析短边或放宽Solver质量门。

`geos_annulus.py` 仅做成熟 GEOS 二维交运算，不组装流动方程；通过隔离venv运行（本次Shapely2.2.0、GEOS3.14.1）。`annulus_mesh_cross.cpp` 将多边形交给既有原生拓扑、合并和质量流程，明确保留共用阶段。示例使用新输出前缀：

```sh
mkdir -p outputs/annulus-replay
/usr/bin/clang++ -std=c++20 -O2 -Iinclude tools/benchmarks/annulus_mesh_cross.cpp build/libcartmesh2d_fv.a build/libcartmesh2d.a -o outputs/annulus-replay/adapter
# 输入 XY 从上述归档解出，下面路径按实际解包位置调整。
outputs/annulus-cross/geometry-env/bin/python tools/benchmarks/geos_annulus.py --geometry outputs/annulus-cross/original-physical-walls.xy --output outputs/annulus-replay/geos96.polygons --n 96
outputs/annulus-replay/adapter geos outputs/annulus-cross/original-physical-walls.xy outputs/annulus-replay/geos96 96 outputs/annulus-replay/geos96.polygons
# --n 128 对应另一档；不是通过移动物面制造面积分数。
```

混合研究输入使用真实产品单内圆外流网格，再裁外域，保留原有近壁/过渡接口；不是绕过产品嵌套闭环拒绝后宣称H4通过。下面 `protect-radius .522` 是已知3层构造包络（最外约.521842m），仅选择并核对768个源近壁格；源ID及多边形一直保留，未按近邻猜物面。裁切结果的未切割原网格和层格进入已有 preserve/immutable 约束，质量修复不重划它们；几何核对显示768格全部保留。

```sh
build/cartmesh2d_hybrid_cli outputs/annulus-cross/inner-wall.xy outputs/annulus-replay/inner 7 7 7 3 .006 1.2 .8 outputs/annulus-replay/inner.foam .1
outputs/annulus-cross/geometry-env/bin/python tools/benchmarks/geos_annulus.py --geometry outputs/annulus-cross/original-physical-walls.xy --output outputs/annulus-replay/hybrid.polygons --base-mesh outputs/annulus-replay/inner.hybrid.solver.cm2d --protect-radius .522
outputs/annulus-replay/adapter crop outputs/annulus-cross/original-physical-walls.xy outputs/annulus-replay/hybrid 128 outputs/annulus-replay/hybrid.polygons
build/cartmesh2d_flow_cli --mesh outputs/annulus-replay/hybrid.solver.cm2d --case custom --boundary outputs/annulus-replay/hybrid.boundaries --nu .1 --speed .5 --convection upwind --viscous-stress symmetric --coupling coupled --linear-policy strict --max-iterations 1500 --tolerance 1e-8 --profile --output outputs/annulus-replay/hybrid-coupled
```

SIMPLE将上例 `--coupling` 换为 `simple`，原默认松弛.6/.25不改。1500步未收敛时不算通过；用既有 `tools/flow/extract_steady_iterate.py` 提取真实接受迭代场和面通量再续算，完整命令已随两个 `run.json` 保存。精确续算不是用解析初值帮助收敛。`native` 模式为同壁面原生均匀切割辅助控制；`hybrid` 模式保留早期双壁条带原型（失败日志保留）；`ready` 导入成熟cfMesh成品多边形，不做原生合并或质量修复，边界由真实外部patch身份给出。

cfMesh完整参考采用foam-extend内 `cartesian2DMesh`，输入是原高精度XY挤出的**仅侧壁 ribbon**，不带front/back盖面。归档 `cfmesh-precise/annulus.fms` 及 `system/meshDict` 可重跑，`maxCellSize .025`、`nLayers 0`。`surfaceFeatureEdges` 默认六位输出会改变源点，已通过唯一六位序列化键恢复原STL点身份后再生成网格，连接及特征不变；不按最近点拟合或修改最终网格。**cfMesh本身仍会贴体变形**，其原生导入和各求解器共享同一实际成品，但与原圆环的壁面误差另列在 `cfmesh-model-difference.json`，不叫严格同物面Cut-cell。

外部求解由 `annulus_cross.py` 编排；`cases()` 指向本机冻结输入，`--directory` 指定新的外部运行目录。每个case目录包含 `comparison.json`、`0/`、`constant/`、`system/`、原始 `checkMesh` 两档及 `solve.log`。标准/扩展检查分开记录，不用原生质量替代。固定运行镜像（独立可执行参考，产品不链接其三维核心）：

- SIMPLE：`opencfd/openfoam-run@sha256:4229997e74defb81548222d511b8e3b95b98305e5df41b8e88b031813fe47eeb`，OpenCFD2606 arm64。
- coupled：`pmariappan/foam-extend-4.0@sha256:63ab10ff5d5886d0428c76e56a50181d4e078a024ce4482c712e5b6e63893143`，实报4.0-268bb07d15d8，amd64模拟运行。

```sh
python3 tools/benchmarks/annulus_cross.py --phase prepare --solver coupled --names original --directory outputs/annulus-replay/original-external-coupled
python3 tools/benchmarks/annulus_cross.py --phase run --solver coupled --names original --directory outputs/annulus-replay/original-external-coupled --budget 180
python3 tools/benchmarks/annulus_cross.py --phase collect --solver coupled --names original --directory outputs/annulus-replay/original-external-coupled
# 其他冻结输入名称：geos96 geos128 hybrid cfmesh；SIMPLE换 --solver simple。
# 混合网格成熟SIMPLE冷启动需要更长研究预算，本次用420秒并在3012步结束。
python3 tools/benchmarks/test_openfoam_interchange.py
```

成熟coupled的源码入口是[主循环](https://github.com/Unofficial-Extend-Project-Mirror/foam-extend-foam-extend-4.0/blob/master/applications/solvers/coupled/pUCoupledFoam/pUCoupledFoam.C)、`UEqn.H`、`pEqn.H`、[couplingTerms.H](https://github.com/Unofficial-Extend-Project-Mirror/foam-extend-foam-extend-4.0/blob/master/applications/solvers/coupled/pUCoupledFoam/couplingTerms.H)。它把动量、压力以及隐式 `grad(p)` / `UDiv(U)` 插入 `fvBlockMatrix<vector4>`，联立求U/p，再从该解重建守恒通量。这里关闭湍流，Up用BiCGStab/Cholesky、线性容差1e-11、块停止1e-8；U/p松弛.999/1，极宽fieldBounds用于避免人为裁剪，收集器拒绝实际Bounding记录。旧版黏度语法需 `nu nu [0 2 -1 0 0 0 0] .1;`。自研coupled则用完整冻结动量响应、块FGMRES及嵌套压力预条件和线搜索。外轮数不是线性工作量，空间离散/压力边界闭合也未完全对齐，不能凭收敛轮数断言提速或精度优劣。

本轮真实交换问题的最小回归在 `test_openfoam_interchange.py`：旧foam-extend可能写全nFaces的neighbour数组并以-1标记边界；读取须保留负号，再按实际patch起始面归一成内面数组。此前丢失物理patch的零流动结果全部排除。二维empty身份来自源patch声明，不用浮点z面积是否非零猜测，并明确拒绝圆环丢失rotor/housing。三个测试只约束文件交换，不重建独立PDE审计链。

真圆参照由二维轴对称稳态方程得到：`ur=0`，切向方程给 `uθ=Ar+B/r`，两壁条件给 `A=-1/3,B=1/3`；径向动量给 `d(pi)/dr=uθ²/r`，积分得 `pi=(r²-r^-2)/18-(2/9)ln(r)`，取外壁pi=0。它表示同心圆流线、内壁速率.5到外壁0、压力由内向外上升，**不需要把实际场投影成轴对称后再算误差**。所有单元都按实际中心和面积比较，压力先去最优面积加权常数。速度RMS/Uref和压力RMS/解析压力跨度分别报告，另外报告局部最大偏差，避免小面积尖峰被全局RMS掩盖。

原逐边切向运动不是上述真圆的精确边界。另设 `annulus_exact_trace.py` 解析Dirichlet验证族：16点Gauss逐面积分已知速度，法向通量用流函数 `psi=r²/6-ln(r)/3` 端点差；128机器精度的流函数/速度×边长尺度只处理通量抵消舍入，非求解验收阈值。原生与外部映射同一逐面值，外部面键的12位舍入仅作唯一匹配、不修改几何。生成器记录面数、净通量和积分舍入调整，原生仍过原边界/守恒门。这一族显式改变边界模型，不能拿其结果替换原壁面运动误差。

```sh
python3 tools/benchmarks/annulus_exact_trace.py --mesh build/laminar-fixtures/annulus-l7.solver.cm2d --output outputs/annulus-replay/original-exact.boundaries
# 追加 --external-case 已准备的新case，可将相同逐面平均值写入其0/U；不得修改已运行的原对照目录。
MPLCONFIGDIR=/tmp/cm2d-annulus-mpl python3 tools/benchmarks/plot_annulus_cross.py --output outputs/annulus-replay/figures
```

科学图直接读取真实网格/接受场；本机中文字体为 `/Library/Fonts/Arial Unicode.ttf`，其他系统须换可用中文字体。展示索引：[理想解析流场](../artifacts/current/annulus-cross-analytic.png)、[实际网格与局部](../artifacts/current/annulus-cross-meshes.png)、[实际收敛场与解析曲线](../artifacts/current/annulus-cross-fields.png)。压力图用对称对数轴同时保留主体和尖峰；不裁掉坏单元、离群压力或未收敛记录。未做空间收敛或产品精度资格扩展。

### Cut-cell 交叉实验

`tools/benchmarks/amrex_cutcell/` 是显式构建的原生二维几何参考工具；`cutcell_cross.py` 只编排真实程序、读取几何/场并计算解析误差，不另建有限体积方程。当前仅支持单位方域中 `y>.5x+b` 的单条平面切割，**不是通用 AMReX 网格导入器**。AMReX 输出实际二维 EB 面开口比例及形心，适配器按共享笛卡尔面身份连接开口端点，再经原生全局拓扑装配。自研路径调用现有原生二维切割器，双方可选同一个原有 `.1` 小单元合并器；因此交叉比较独立切割几何及独立外部求解，没有声称独立比较了整套拓扑装配/合并策略。

本机只读 AMReX 快照程序自报 **26.08**，没有可核验的本地 Git commit，准确来源范围及实际源文件树 SHA 由 `artifacts/current/cutcell-cross.json` 记录，不虚构上游提交。外部源码不修改、不复制入仓库，参考 [AMReX EB2 文档](https://amrex-codes.github.io/amrex/docs_html/EB.html)。CMake 强制 `AMReX_SPACEDIM=2`、CPU、无 MPI/Fortran；原生求解器不会因此新增运行依赖。参考 OpenFOAM 为本机已有 **OpenCFD 2606 arm64** 镜像，脚本按实际 digest 固定；单层挤出和 `empty` 仅用于外部二维问题交换，不把三维生成器压成二维。

AMReX 几何研究显式设 `eb2.small_volfrac=0`，防止把极小流体部分删除后误称为“小单元求解成功”；不改上游默认或本产品门。适配器的面积分数与形心偏差分别按背景格面积/间距归一化，`1e-9` 仅为导入一致性拒绝预算：高于双精度构造舍入、低于最小样本面积分数，用于发现转换错误，检查成本与导入格数成正比；实际偏差在舍入范围。共享端点只在 `64 epsilon` 格坐标舍入范围认作网格角点，不用这个预算删除可解析短边。求解仍执行原 Solver 质量、守恒和线性/非线性停止门。

先构建根目录现有库及 CLI，再单独构建参考工具；下面的新实验目录须不存在，重复研究换新目录，不覆盖已有接受场和失败记录：

```sh
cmake --build build --target cartmesh2d_flow_cli cartmesh2d_cli --parallel 4
cmake -S tools/benchmarks/amrex_cutcell -B outputs/cutcell-cross-reproduce/amrex-build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_C_COMPILER=/usr/bin/clang -DAMREX_SOURCE=/Users/Zhuanz/Desktop/mesh-reference-sources/amrex-development -DCARTMESH_ROOT=/Users/Zhuanz/Desktop/cartesian-mesh-generator-2d
cmake --build outputs/cutcell-cross-reproduce/amrex-build --target cutcell_reference --parallel 4
python3 tools/benchmarks/cutcell_cross.py --phase generate --root outputs/cutcell-cross-reproduce/matrix --adapter outputs/cutcell-cross-reproduce/amrex-build/cutcell_reference
python3 tools/benchmarks/cutcell_cross.py --phase native --root outputs/cutcell-cross-reproduce/matrix
python3 tools/benchmarks/cutcell_cross.py --phase prepare-foam --root outputs/cutcell-cross-reproduce/matrix
python3 tools/benchmarks/cutcell_cross.py --phase foam --root outputs/cutcell-cross-reproduce/matrix
python3 tools/benchmarks/cutcell_cross.py --phase collect --root outputs/cutcell-cross-reproduce/matrix
```

Docker daemon 须可用；`foam` 阶段只启动固定镜像，逐例运行标准/扩展 `checkMesh` 与 `simpleFoam`，每例45秒研究预算。首次挂载权限/镜像入口错误单独保留，不当成数值失败。原生逐例60秒、1500轮预算；超过预算或上限不算收敛。两者均显式零速度/压力迭代初值，边界速度取同一解析面均值。原生 strict `1e-8` 的尺度/停止门与 OpenFOAM 自身 `residualControl` 分别记录；不从轮数推出提速，也不把两套有限体积离散视为同一个代数系统。解析 L2 只作诊断，单位为速度 `m/s`、去均值运动学压力 `m²/s²`，未新增精度通过阈值。`generate` 中的12位坐标签名仅用于匹配两份输出多边形，不进入几何或求解。

`tools/benchmarks/annulus_merge_probe.cpp` 是原圆环的**固定物面合并干预**。选择量为“单元面积 / 最大相邻单元面积”，不是失去背景 lineage 后猜出的 Cut-cell alpha；`.02/.05/.1` 是本次研究干预档，保留既有 `.01` 相邻面积比质量底线，未改产品策略。默认档将被选择的小格与最长共享边的非小邻居合并，可能产生新的质量失败；`greedy` 档逐候选调用原合并器并保持全部原 Solver 质量门，仅接受通过的候选。它会反复全局重建/检查，成本和拒绝次数单列，**不是新的生产高效修复算法**，不能因输出了较少单元就宣布成功。

```sh
/usr/bin/clang++ -std=c++20 -O2 -Iinclude tools/benchmarks/annulus_merge_probe.cpp build/libcartmesh2d_fv.a build/libcartmesh2d.a -o outputs/cutcell-cross-reproduce/annulus_merge_probe
outputs/cutcell-cross-reproduce/annulus_merge_probe build/laminar-fixtures/annulus-l7.solver.cm2d outputs/cutcell-cross-reproduce/annulus-noop .01
outputs/cutcell-cross-reproduce/annulus_merge_probe build/laminar-fixtures/annulus-l7.solver.cm2d outputs/cutcell-cross-reproduce/annulus-direct02 .02
outputs/cutcell-cross-reproduce/annulus_merge_probe build/laminar-fixtures/annulus-l7.solver.cm2d outputs/cutcell-cross-reproduce/annulus-guarded02 .02 greedy
outputs/cutcell-cross-reproduce/annulus_merge_probe build/laminar-fixtures/annulus-l7.solver.cm2d outputs/cutcell-cross-reproduce/annulus-guarded05 .05 greedy
```

没有 `build/laminar-fixtures` 时，从既有 `artifacts/current/laminar-foundation/failure-inputs.tar.gz` 解包到忽略目录并替换输入路径。退出码3表示有效几何导出但 Solver 质量拒绝，不生成供原生求解的边界文件；退出码1是构建/候选失败。获准网格再用原生 CLI 的 `--case custom --boundary PREFIX.boundaries --nu .1 --speed .5 --convection upwind --coupling simple|coupled --linear-policy strict --tolerance 1e-8 --max-iterations 1500 --profile --output NEW_PREFIX` 运行。原6,208格输入、无合并重建及合并后场分别保存，合并可能去除共线分面，不能忽略重建控制。

结论和未解决问题见[当前状态](CURRENT_STATE_CN.md#cut-cell-网格与-simple-的交叉诊断)。`cutcell-cross-fields.tar.gz` 保存直斜壁12份实际网格/质量、统一零初值后的接受场、外部控制/检查/场和原圆环失败/合并证据，供后续诊断与复算。第一次使用原生默认初值的直斜壁运行未进入最终对照；构建目录和第三方库仅留本机 `outputs/`，不提交源码快照或二进制依赖。

### SIMPLER 与 PIMPLE

`FlowCoupling2D::Simpler` / `::Pimple` 在 `src/fv/Incompressible2D.cpp` 复用原动量装配、Rhie–Chow、压力边界及严格停止条件，CLI 用 `--coupling simpler|pimple`。桌面开放 SIMPLE、既有 SIMPLEC 和这两种扩展，默认仍 SIMPLE；当前资格和失败见[当前状态](CURRENT_STATE_CN.md#simple-系列扩展)。成熟算法的本地实现不自动继承外部软件的精度资格。

冻结动量行的共同松弛对角为 `aP`，面积为 `V`，代码 RHS 已含原压力力项 `-V*GpOld`。`U+(rhs-AU)/aP` 使用当前邻点贡献，等于 `H/aP-(V/aP)*GpOld`，并非再次完整解动量。SIMPLER 先由该伪速度经原 Rhie–Chow 与压力方程求增量，取 `p=oldP+alphaP*deltaP`，再按压力梯度的真实变化更新冻结动量 RHS、求解两分量动量；最终压力校正**只修速度和守恒面通量**。该路径仅支持稳态。

PIMPLE 每次外迭代完整解动量后做压力—速度校正，后续内校正重算当前邻点的 `H/aP` 并扣除当前压力与冻结原压力的梯度差；所有内校正共享冻结行、rAU 和原外迭代速度/梯度/面通量松弛缺陷，保留非定常旧物理步的通量缺陷。显式应力、非正交动量和对流重构在下一外迭代刷新。外迭代继续到原严格门，有限次内校正不等于时间步已收敛。

`--momentum-correctors 1..8` 仅用于 PIMPLE，默认 2，1 次逐位恢复普通 SIMPLE 一步；上限控制每轮压力工作量，**不是精度阈值**。每次动量校正内仍有 `--pressure-corrections 1..4` 次非正交压力校正，两种循环分别控制。profile 的 `pressureSolves`、迭代数和跳过的零更新 pass 包含前置预测及全部内校正；历史字段 `simpleIterations` 计分离式外迭代。两者要求固定边界角色、普通出口回流拒绝、无材料反馈；opening/farfield 明确拒绝。PIMPLE 支持现有 backward Euler 固定/自适应步长与检查点；同控制分段/连续续算由现有 CLI 测试检查逐位一致。

```sh
mkdir -p outputs/simple-extensions/inputs
tar -xzf artifacts/current/laminar-foundation/failure-inputs.tar.gz -C outputs/simple-extensions/inputs
build/cartmesh2d_flow_cli --mesh outputs/simple-extensions/inputs/annulus-l5.solver.cm2d --boundary outputs/simple-extensions/inputs/annulus-l5.boundaries --case custom --nu .1 --speed .5 --convection face-limited-linear --coupling simpler --tolerance 1e-8 --profile --output outputs/simpler
# 上条命令改为 --coupling pimple --momentum-correctors 2 即运行 PIMPLE。
ctest --test-dir build -R 'cartmesh2d_(flow_simple_extensions|flow_cli_validation|transient_flow_cli)$' --output-on-failure
```

`flow_simple_extensions_test.cpp` 覆盖非均匀格、制造压力斜率、非零压力出口、滑移不等对角、时间步、取消和拒绝组合。沿用 `flow_coupled_test.cpp` 的绝对 `2e-6` 迭代误差预算比较 u/v/运动学 p/面通量/积分壁力，求解容差 `1e-9`；目的是排除固定点变化，不是空间精度门。实际旧网格的同解差和失败见 `artifacts/current/simple-extensions.json`，配套归档保留代表原生场、原始失败日志和 App 工况/截图，供后续耦合修复与复算；本机 `outputs/simple-extensions/` 的重复导出及编译产物可重建，不进入 Git。

机制参考：[Loci-Stream 开发者理论说明](https://streamline-numerics-inc-loci-stream.readthedocs-hosted.com/en/latest/theory/chapter_pressure_velocity_algorithms/chapter_pressure_velocity_algorithms.html)、[OpenFOAM 开发者的校正循环控制说明](https://doc.cfd.direct/openfoam/user-guide-v13/fvsolution)。未引入外部求解内核。

### 速度与压力松弛因子研究

`tools/benchmarks/relaxation_study.py` 编排原生 CLI、比较原生 CSV 并绘图，不另建 PDE 求解或方程验收。`relaxation_schedule.py` 是保留同网格 u/v/p/owner 方向面通量的分段反馈实验，使用既有 `extract_steady_iterate.py`。当前结果、失败和限制见[当前状态](CURRENT_STATE_CN.md)及 `artifacts/current/relaxation-study/`。Linux 研究环境依赖 Python、NumPy、Matplotlib、`taskset` 和原生 Release CLI；测量子进程通过 `wait4` 单独记录 CPU 时间和 RSS，避免并发子进程串账。

研究变量为 `(alphaU, alphaP)`，目标是**在同一离散解、同一原生严格验收下，降低完整进程耗时**。SIMPLE 与 `simple-consistent` 分别选参；`coupled` 不使用这两个控制，只作固定对照。旧默认 `.6/.25` 不是寻优结果。所有算例为二维稳态、不可压、恒黏度层流，对称应力、inertia=1、IC(0)、无 Anderson，SIMPLE 压力校正 4 次；中间线性精度 adaptive，最终严格认证。普通容差 `1e-8`，同方程 coupled 参考解 `1e-10`。每次参数试验从相同原生冷初值开始。

训练集为原 484 格圆环（nu=.1、Uref=.5、upwind）和原生生成的 900 格 Poiseuille 通道（nu=.1、峰值入口 1、upwind）。参数粗网格为 `alphaU=[.2,.4,.6,.8,1]`、`alphaP=[.05,.1,.25,.5,1]`，两方法共 100 点；每点最多 3000 步/60 秒。局部加密各 9 对参数，仍在两训练算例上运行。前三名和默认各串行重复 3 次，最接近的两个候选各追加 7 次；按两算例**耗时中位数的几何平均**选参，候选须在所有训练重复中通过。其排名等价于除以各算例 coupled 固定成本后的几何平均排名。单次并行粗扫时间用于筛选，不能冒充稳定的性能差距。

留出集在这轮选参中不参与排名：原 1672/6208 格圆环、900 格 Re100 方腔、3063 格 DFG Re20。圆环、方腔为 upwind，DFG 为 limited-linear，沿用原显式入口/出口/壁面；不同格式不互相计算提速。圆环网格同时改变原多边形细分，不是固定几何的网格收敛研究。历史上已看过部分留出算例，故这是**未参与本轮选参的迁移验证**，不是盲测或独立外部数据集。选参锁定后才开始这些参数运行；后续针对失败例的诊断必须另列，不能反过来改锁定排名。

同解检查只比较原生输出：速度用面积加权向量 RMS/Uref；压力用面积加权 RMS/(Uref²+nu·Uref/H)，只对闭合域差值减去面积均值，压力出口域保留绝对基准；面通量用 RMS/(Uref·sqrt(平均单元面积))。H 为物理域高度，圆环 2、通道/方腔 1、DFG .41。三项均须不大于 `1e-4`，即该归一化下的 0.01%；此研究门用于排除提前停下或不同离散解，不改变原生验收，也不证明压力/受力的空间精度。固定网格、格式、边界的参考偏差与物理参考误差分开。

原生 strict 最早在第 10 轮通过；多个结果同为 10 轮时，仍须比较实际块/压力工作量和完整耗时。CSV 显式区分冷启动试验、参考和分段续算，`objectiveSeconds` 只为成功冷启动试验填写；不能将最后一个续算段的短耗时当作从零收敛成本，也不能将计划中的 50 步段末停止统计为固定参数扫描的失败。coupled 的两个松弛列留空，因为它不使用这些参数。

参考方案有一项显式修订：6208 格圆环在 `1e-10` 下，冷 adaptive、冷 strict、从已收敛 `1e-8` 场续算均耗尽 FGMRES 预算；整体线性范数已达标，最大行残差仍约 `2.9e-12`～`4.0e-12`，高于该容差对应的 `1e-12` 行门。三个失败均保留，不称通过。此算例单独用 `--reference-tolerance 1e-9`，其余五例仍 `1e-10`；候选 `1e-8` 与同解门不变，内核及原生验收未修改。`bootstrap-reference --case annulus-l7` 可复现失败的普通收敛场→紧容差认证路径；初始化场不算紧容差参考。这项限制不能隐藏在总成功率中。

记录区分原生严格收敛且同解、不同解、迭代上限、时间上限、线性失败、原生异常与联立失败；超预算不等于已证明必然发散，失败的短耗时不进入最优排名。批次顺序按固定种子打乱，粗扫最多两个独立 CPU 核各一个进程，库线程数固定为 1；候选确认和留出成本对照串行。完整墙钟时间包含进程启动、网格读取、求解及所有原生导出，不把嵌套压力时间重复相加。共享云主机仍有时钟噪声，记录重复范围、迭代数和 CPU 时间，不从小的计时差距宣称唯一最优。

首次运行示例（已有研究目录请复用其 `run`，不要重新 `init`）：

```sh
python3 tools/benchmarks/relaxation_study.py init
python3 tools/benchmarks/relaxation_study.py run --batch references-training-v1 --methods coupled --reference --iterations 1000 --seconds 180 --workers 1
python3 tools/benchmarks/relaxation_study.py run --batch grid --iterations 3000 --seconds 60 --workers 2
python3 tools/benchmarks/relaxation_study.py select --batches grid --output outputs/relaxation-study/grid-selection.json
python3 tools/benchmarks/relaxation_study.py report
python3 tools/benchmarks/test_relaxation_study.py
```

`--pairs .8,.15 .8,.25` 指定局部候选，`--repeats` 重复测量，`--cases holdout` 切换留出集；实际所有批次控制、选参锁定与原生命令随结果封存。`--root` 可指定另一研究目录。目录内 manifest 固定二进制、输入、物理控制及训练划分；批次/spec 不允许原地改动，已完成记录可恢复。中断的未完成运行保留证据，换批次继续，不覆盖日志或伪造累计成本。启动失败与实际 CFD 失败分开计数。

`artifacts/current/relaxation-study/native-evidence.tar.xz` 供后续原生复算、机制研究和同网格续算使用：含所有实际运行的命令/日志/残差/记录、原输入及新增原生网格、固定参考、全部留出对照和代表失败/最终场、反馈初值、冻结排名及历史编排脚本哈希快照；逐文件 SHA 见 `evidence-index.json`。完整原生 CSV 留在 `outputs/relaxation-study/runs/`，可重建的 VTK、fields.json 与字节相同边界副本清理记录在 `pruned-exports.json`。编译二进制不入仓库，按源提交/工具链重建；归档中的脚本快照用于追溯，执行入口仍为仓库 `tools/benchmarks/`。上述整理不删除旧研究、用户原数据或失败记录。

反馈实验每 50 步才决策一次：用原生动量、连续性、全局不平衡及完整响应路径的面通量缺陷，分别除以原验收门后取最大值。完成段若相对上次接受值增长超过 2 倍或原生失败，则退回上次 u/v/p/flux 并将两个因子减半，最低 .05；接受段比值小于 .8 时两因子乘 1.25，最大 1。起点 `.6/.25`，这些是预先固定的启发式控制常数，不是优化出的最佳策略。固定参数分段对照不使用反馈拒绝规则。预算计入拒绝段，原生异常缺少最终步数时按整段保守收费并标记步数未知；总耗时还计入提取/重启/反馈开销。最后仍须原生完整严格认证及同解检查。该实验不能冒充内核自适应实现；只有与同样分段重启的固定参数和不间断固定参数比较，才能判断反馈收益。

上述完整矩阵保留已知失败。`--mesh-layout square` 用等尺寸笛卡尔背景，适于 DFG 的长矩形外域；须使用新输出目录，不能把不同网格策略冒充同一细化序列。`--prepare-only` 只准备输入。当前细曲壁网格仍可能质量失败，求解器亦可能发散，均需继续修复。

仓库保存了约 574 KiB 的原始失败输入，云端无需读取本机 `outputs/` 即可重现细圆环与 DFG：

```sh
mkdir -p outputs/laminar-inputs
tar -xzf artifacts/current/laminar-foundation/failure-inputs.tar.gz -C outputs/laminar-inputs
build/cartmesh2d_flow_cli --mesh outputs/laminar-inputs/annulus-l7.solver.cm2d --case custom --boundary outputs/laminar-inputs/annulus-l7.boundaries --output outputs/annulus-default --nu .1 --speed .5 --convection face-limited-linear --pressure-preconditioner ic0 --max-iterations 12000 --tolerance 1e-8 --profile
build/cartmesh2d_flow_cli --mesh outputs/laminar-inputs/dfg20-l7.solver.cm2d --case custom --boundary outputs/laminar-inputs/dfg20-l7.boundaries --output outputs/dfg-default --nu .001 --speed .3 --convection limited-linear --pressure-preconditioner ic0 --max-iterations 16000 --tolerance 1e-8 --profile
```

CLI 支持版本化 `--case-file FILE.json`，格式为 `{"format":"cartmesh2d-native-flow-case-v1","options":{...}}`。`options` 使用去掉 `--` 的原 CLI 名称；数值采用 JSON number，枚举/路径采用 string，`profile` 采用 boolean。文件内路径相对 JSON 所在目录解析，显式命令行选项优先；未知版本、重复键和类型错误明确失败。该格式与桌面请求存档格式分开。解包上面的输入后，可直接复现两个早期 Anderson 算例：

```sh
build/cartmesh2d_flow_cli --case-file artifacts/current/laminar-foundation/dfg20.case.json
build/cartmesh2d_flow_cli --case-file artifacts/current/laminar-foundation/annulus.case.json
```

云端的早期历史实验在以上命令追加 `--steady-acceleration anderson --anderson-history 16 --anderson-start 1`。两个选项显式记录到结果 JSON；默认仍是历史4、从第10步开始，未加速仍为产品默认。历史容量1–32是内存/QR工作量控制，不是物理验收门。历史外推先试完整步，再最多试相对普通候选的1/2、1/4步长；每个候选重算真实残差并过原守恒/边界门，全部失败才恢复普通候选并清空历史，成功外推不能作为最终认证步。`accelerationBacktracks`计额外试探，`accelerationDampedAccepted`计缩短后接受的候选，耗时包含全部试探。原始DFG输入用于`cartmesh2d_flow_stability`回归，检查实际回退路径与原strict停止条件。该网格虽收敛，Cd/Cl/压降尚不合格；早期场见`artifacts/current/laminar-foundation/early-anderson*`，当前前后对照见源码对照证据。

`tools/benchmarks/openfoam.py` 将已有 OpenFOAM 导出网格的边界重新分组，保留点、内部面、owner/neighbour 和单元编号，准备稳态层流对照。准备后在有 OpenFOAM 的环境执行 `checkMesh -case CASE`、`simpleFoam -case CASE > CASE/solve.log 2>&1`，再 `python3 tools/benchmarks/openfoam.py --output CASE --collect`。例如：

```sh
python3 tools/benchmarks/openfoam.py --output outputs/foam-cavity --source outputs/laminar-baseline/meshes/cavity100-l6/openfoam --mesh outputs/laminar-baseline/meshes/cavity100-l6/mesh.solver.cm2d --case cavity100
```

该对照的 OpenFOAM `linearUpwind` 与原生 `limited-linear` 分别记录，不能声称格式完全相同。OpenFOAM 自身残差通过仍须比较物理误差；未完成或中断的运行不算通过。`--stabilized` 显式选择 `cellLimited leastSquares 1` 和 `limited 0.5` 非正交修正。每个尝试使用新目录；读取器按最终场时间查找唯一受力记录（包含 OpenFOAM 重跑生成的 `force_0.dat`），拒绝旧记录混用，并保存实际 fvSchemes/fvSolution/controlDict 与 CPU/墙钟耗时。

`artifacts/current/laminar-foundation/rejected-implicit-pressure.patch` 保留未奏效的完整非正交压力 GMRES 试验，可在原 `6ec4edf` 求解器上复现，CLI 选项为 `--pressure-correction implicit`。它只用于解释已试过的方向，不应在未解决问题时作为修复合入。对应细圆环第 10/20 轮动量残差约 `4309.5 / 694436`，故下一步须继续检查外迭代动量—压力耦合，而不能假定只换压力线性求解器即可消除发散。

不可压内部实现入口为 `src/fv/Incompressible2D.cpp`；工况与预设在 `FlowCase2D.cpp`，动量在 `FlowMomentum2D.cpp`，Rhie–Chow 在 `FlowFlux2D.cpp`，压力线性循环在 `FlowPressure2D.hpp`，监测在 `FlowConvergence2D.cpp`，受力在 `FlowForces2D.cpp`。`FlowEquation2D` 统一 SIMPLE、coupled 和原生诊断的梯度/压力力/应力重构、双分量动量与冻结残差；几何缓存归单次求解所有，边界值及场值实时读取，不持有接受解。正交压力块由两个耦合路径共用，`FlowKrylov2D.hpp` 保存原 FGMRES 实现。齐次线性作用和非零边界/源项分开；物理连续性检查包含压力参考格，只有线性系统替换 gauge 行。coupled 工作数组按需分配。私有参数接口集中于 `FlowSolverDetail2D.hpp`，公共求解 API 保持不变。

早期模块拆分的六例哈希见 `artifacts/current/laminar-foundation/refactor-hashes.json`；共同方程核以 `2e94198` 为对照，九个基础路径及四组旧圆环/DFG 结果见 `core-refactor.json` 和 `core-refactor-fields.tar.gz`。归档包含原生 C++ 场对照驱动、逐位场、实际 CLI 命令、成功与失败 CSV、检查日志及真实 App 截图；原输入复用 `failure-inputs.tar.gz`。`tests/flow_equation_test.cpp` 检查冻结 Upwind 方程的仿射性、共享面守恒、压力边界及改变边界掩码后的实时值，`4096*epsilon` 是浮点装配一致性预算，不能作为空间精度门。数值修改须另行比较精度，不能援引重构哈希作为修复证据。

`--convergence engineering` 是独立的稳态层流停止模式：最大单元动量残差须低于 `--tolerance`，50 步连续性和流量不平衡仍各低于 `1e-8`，场变化累计及归一化物理监测范围低于既有 `1e-3`，残差下降三阶或低于 `1e-5`；必须用普通严格线性步确认。监测包含动能、各边界组流量/受力/平均压力/力矩。窗口稳定不保证远离固定点的慢模态已经衰减，故须针对压降、受力和参考误差选择容差；已有方腔反例及容差扫描见 `engineering-simplec.json`。不能将 `1e-5` 作为通用精度设置。

`--coupling simplec --pressure-relaxation 1` 显式试验稳态层流 SIMPLEC。校正响应采用松弛动量行的 `area/(diag+sum(offDiagonal))`，预测 Rhie–Chow 仍保留原系数以共享固定点方程；非正对角明确失败，非定常及物性耦合不支持。SIMPLE、压力松弛 0.25 和 strict 仍为默认。旧 DFG 已发现非物理解分支，旧圆环该组控制失败，故当前不推荐自动切换。`tools/benchmarks/laminar.py` 在输出目录保存按哈希命名的求解器快照，排队运行不会混入后续重新构建的版本。

`--coupling coupled` 在 `FlowCoupled2D.cpp` 中直接求原未松弛稳态固定点。每次冻结真实对流通量，组装共享的动量与 Rhie–Chow 仿射算子，块 ILU(0)/压力预条件 FGMRES 提出更新，随后共同对 u/v/p/flux 回溯，重新计算候选自身 rAU、动量、连续性及面通量一致性。有限 Krylov/回溯预算耗尽明确失败并保留最后接受状态；初步即失败时 iteration 0 明确标为初始化，不能作为物理场。最终须 strict 线性步、原 strict 收敛门、全部单元连续性（包括压力基准格）及 `min(tolerance,1e-8)` 的最大面通量缺陷门。`--linear-policy adaptive` 仅控制中间 Krylov 精度，最终认证不放宽。SIMPLE 的松弛和 pressure-correction-pass 控制不用于 coupled，结果 JSON 对这些字段写 null。只支持稳态恒黏度固定边界角色；其他组合明确拒绝。它没有改压力/速度空间离散，也没有解决当前压力和受力误差。

解包保留失败输入后，版本化 coupled 算例可直接运行：

```sh
build/cartmesh2d_flow_cli --case-file artifacts/current/laminar-foundation/annulus-coupled.case.json
build/cartmesh2d_flow_cli --case-file artifacts/current/laminar-foundation/dfg20-coupled.case.json
```

`coupled-stability-fields.tar.gz` 保存两个实际解的 u/v/p/flux、原生残差与 profile；输入仍为原失败归档。`tests/flow_coupled_test.cpp` 检查共享固定点、源项与压力基准、支持/拒绝组合、真实线性预算失败后的场保留及严格停止条件。旧圆环/DFG 的完整运行证据与精度限制见 `coupled-stability.json`，不能把 Krylov 轮数与 SIMPLE 外迭代直接比较，profile 的嵌套压力时间不可重复累计。

### 完整动量响应与壁面分段修复

`--coupling simple-consistent` 复用 `FlowCoupled2D` 的共同方程核和 FGMRES，将 SIMPLE 的完整冻结动量响应、压力和守恒面通量一起求解。设原紧凑对角为 `aU/aV`、`D=max(aU,aV)`，原速度松弛为 `alphaU`，则附加对角分别为 `hU=D/alphaU-aU`、`hV=D/alphaU-aV`，右端加入同样的 `hU*Uold/hV*Vold`；非正交扩散及对称应力都在原生线性作用内。面压力系数仍为 `rAU=alphaU*V/D`，面方程保留 `(1-alphaU)*(phiOld-T(Uold))`，其中 `T` 是原偏斜修正速度插值。压力未知量是 `pHat=pOld+pc`，保存时才取 `pOld+alphaP*pc`，速度和面通量共享完整校正。固定点处这些松弛项消去，恢复原未松弛动量与 Rhie–Chow 方程。

该选项保留 `--velocity-relaxation` / `--pressure-relaxation`；`--pressure-corrections` 不参与，summary 明确输出 null。与 `coupled` 相同，仅支持恒黏度、固定边界角色、稳态 strict，不支持动态 opening/farfield、NormalInlet、Anderson、engineering、非定常或材料更新。原收敛门、最终严格线性步及所有单元的连续性和面通量一致性检查保持。迭代上限或线搜索失败保留最后接受场，不能算成功。默认仍为 `simple`，`coupled` 则直接求未松弛固定点；完整响应选项不是大规模效率保证。

```sh
build/cartmesh2d_flow_cli --case-file artifacts/current/laminar-foundation/annulus-coupled.case.json --coupling simple-consistent --convection upwind --output outputs/annulus-response
build/cartmesh2d_flow_cli --case-file artifacts/current/laminar-foundation/dfg20-coupled.case.json --coupling simple-consistent --output outputs/dfg-response
```

速度梯度现在把**同一单元、同一物理边界组、相接且共线**的边界分面视为一条长度加权的重构样本，避免网格导出多切一段就改变该直壁的 LS 权重。边界组来自 custom 的名称/类型或预设的角色，保持固定，齐次 Krylov 与真实边界使用同一模板；实时边界值仍逐面读取。不同边界组、不相接片段及超出舍入识别预算的转角不合并；预算内的极小真实折角可能被判为共线，不改变几何本身。识别采用 `128*epsilon*坐标或边长尺度` 的坐标构造舍入预算，不是新增网格或物理验收阈值。真实单元、面、边界值、黏性/压力通量和受力积分均不合并；压力梯度模板保持原规则。

回归位于 `flow_face_test.cpp`（分段不变、仿射精确、缩放/平移、断开片段、转角及缓存）和 `flow_coupled_test.cpp`（松弛不改变固定点、非对称滑移对角、源项、失败保场及入口约束）。实际失败例与限制由[当前状态](CURRENT_STATE_CN.md)及 `artifacts/current/laminar-foundation/targeted-repair.json` 维护；早期纯重构哈希仅对应其历史版本。

### 压力速度失稳的原生诊断

`tools/benchmarks/momentum_terms_probe.cpp` 直接消费共同方程核，保存动量预测、压力校正、各面黏性/压力项及 LS/Gauss 梯度。构建 `cartmesh2d_flow_momentum_probe` 后运行 `build/cartmesh2d_flow_momentum_probe MESH BOUNDARIES STEPS BLOCK CELL WEIGHT power`；`BLOCK=0/1/2` 分别为原更新、壁面 compact 转置应力隐式 2×2 块、全部 compact 转置应力隐式块，其他重构项及原 scalar rAU 保留。`CELL=3710` 对应用户编号 3711，权重 `.5` 只改变探测映射。拆分前后同场残差恒等式和原生产一步的逐场差异必须先核对；幂迭代模式不构成完整谱界或物理验收。原基线源码/头文件、CSV 与哈希保存在 `momentum-splitting-traces.tar.gz` / `momentum-splitting.json`，当前核一步另有记录。

`wall-refactor-experiments.json` / `wall-refactor-experiments.tar.gz` 保存独立 C++ 二次壁面完整应力及 divU 投影原型、原场和失败日志。复现这些旧实验须在 `2e94198` 编译同版本库与头文件，再按 JSON 命令重建；不能把旧库链接到新 inline 头文件。外部 MOOSE 参考只记录固定提交、URL 和哈希，归档不复制其代码。原型使用的 U/V 配对包装仅供研究。二次 QR 精确再现多项式不保证动量响应正性；细圆环九个非正 owner 行被原检查明确拒绝，未做裁剪或替代响应。继续开发者消费这些最小反例及原始场；本机 `outputs/` 中库、二进制、对象文件可重建，不作为交付归档。

`tools/benchmarks/coupling_probe.cpp` 是显式构建的研究程序，直接调用生产求解器和原生动量/梯度/线性算子，不是另建方程验证链。固定旧圆环的 `nu=.1`、`Uref=.5`、速度松弛 `.6`、压力松弛 `.25`、无对流；只接受闭合且速度全规定的壁面域。macOS 用系统 Cholesky，其他平台 IC(0)。`outer` 对真实一步求解作 `T(x)=F(x)-F(0)` 幂迭代，u/v/p 用单元面积加权内积，面通量同时传播并单列特征残差；`pressure` 探测冻结扩散 rAU 下的非正交延迟修正映射。输出是观测模式，不是完整谱界；`diffusion` 的场变化阈值只判断隔离子问题衰减，不代替完整 CFD 验收。

```sh
cmake --build build --target cartmesh2d_flow_coupling_probe --parallel 4
build/cartmesh2d_flow_coupling_probe outputs/laminar-inputs/annulus-l7.solver.cm2d outputs/laminar-inputs/annulus-l7.boundaries outer symmetric 240 1
build/cartmesh2d_flow_coupling_probe outputs/laminar-inputs/annulus-l7.solver.cm2d outputs/laminar-inputs/annulus-l7.boundaries outer symmetric 400 .5
build/cartmesh2d_flow_coupling_probe outputs/laminar-inputs/annulus-l7.solver.cm2d outputs/laminar-inputs/annulus-l7.boundaries pressure symmetric 240 1
build/cartmesh2d_flow_coupling_probe outputs/laminar-inputs/annulus-l7.solver.cm2d outputs/laminar-inputs/annulus-l7.boundaries diffusion symmetric 2500
```

输入由上面的 `failure-inputs.tar.gz` 解包；替换 l5/l6 可复现粗档。参数依次为模式、应力形式、步数、映射权重、压力校正次数（默认 4）。权重 `.5` 仅探测 `(I+T)/2`，不修改产品求解器。细网格原映射约 `−2.67314`；半步映射约 `1.04892`，对应原映射近实增长模式约 `1.09785`（单元/面通量特征残差约 `1e-5`），因此标量整体松弛不能消除已观测到的全部失稳。只做一次压力校正仍有约 `−2.25654` 模式，压力内部延迟修正不是唯一问题。

`coupling-diagnosis.json` 保存源码/二进制/输入哈希、真实命令及失败记录；配套小体积归档保存原生 CSV/日志和基于 `ab4aba3` 的未奏效补丁。原始场、试改二进制和继续研究材料保留于本机 `outputs/laminar-diagnosis/`；云端消费者用现有输入归档和本程序重建，不依赖本机绝对路径。尚未把这些试改作为修复合入，也未改变默认控制。原生 `FlowMomentum2D.cpp` 按同一共享面向 owner/neighbour 组装相反动量贡献；`FlowFlux2D.cpp`、`FlowPressure2D.hpp` 和 `Incompressible2D.cpp` 分别承担面通量、压力修正和耦合更新，是后续算法修复入口。

有向面几何回归在 `tests/solver_export_test.cpp::reversedCentreConnectorRegression`，使用实际 DFG 面 4591 的两单元多边形，保留三个尺度、端点反转及 patch-local 检查。`evaluateSolverFaceGeometry2D` 的法向来自 owner 的 CCW 顶点循环；内部面要求 owner 到面、面到 neighbour、owner 到 neighbour 的有向法向距离均正，边界要求 owner 到面为正，与 FVM 原有可容许条件一致。原无向 nonorthogonality 诊断仍保留，不能再单独放行反向中心。`quality-normal-distance-fields.tar.gz` 包含原 3,086 格失败输入及修复后 3,085 格输入、真实 u/v/p/flux；后续求解仍须通过原 strict 条件和物理误差检查，固定几何 5,891 格的误差恶化见 `dfg-controlled-refinement.json`。

### 开源求解器机制对照

参考为容器实际 OpenCFD v2412（head `45e7c4a005ff9a310cdc073f65d8fbd2b29e1a0c`）及 MOOSE 固定 commit `410cc90c7e8e6165e77caafb7e49ac206d5ba6b5`，具体来源 URL、源码哈希、参数和真实结果在 `artifacts/current/laminar-foundation/open-source-mechanisms.json`。外部源码只读，归档保存独立推导的修改、原生算子探针及真实场，不包含外部代码片段。

- OpenFOAM 先构造单元 HbyA，再按选定的 `flux(HbyA)` 插值；默认 `linear` 自身没有 skew 修正。本地速度做 skew，而 `I(rAU Gp)` 没有相同重构。补齐整体 H 只小幅改善圆环误差，未解决 SIMPLE 映射增长。普通固定速度壁的 `fixedFluxPressure` 经 `constrainHbyA` 后可以得到零 snGrad，不能据此推断必须给非零壁压梯度。
- MOOSE 默认压力动量与 RC 使用同一压力 GG 梯度，外推边界值与 GG 联立；本地用多环 LS 外推后另算守恒 Gauss 压力力。在同权重、常标量 D、同梯度条件下，MOOSE 的中心连线分解系数为 `(S·d)/|d|²`，本地 over-relaxed 系数为 `|S|²/(S·d)`；两者都是离散选择，差异不是修复证明。本次没有直接改变非正交分解或移植 MOOSE 边界 fallback。
- 本地原压力通量精确分为 `−Dbar*T*(Δp−Lbar·d) + I[D*(G−L)]·S + w*(1−w)*(Di−Dj)*(Li−Lj)·S`，L 为 LS、G 为守恒 Gauss。`rc_terms.cpp` 调用原生网格、动量、重构和通量算子，对保留模式、DFG 单格压力脉冲及真实圆环压力作符号/项贡献诊断；不是另建 Python 方程审计。DFG 此探针固定诊断黏度 `.1`，不冒充 Re20 物理计算。消除已观察向量的负能量项不代表全压力空间半正定，也没有消除完整 SIMPLE 增长。
- OpenFOAM 内部转置黏性项由单元梯度插值得到；本地额外按 compact normal derivative 调整面梯度。仅替换内部转置项的实验令圆环压力误差更大，壁面及 trace 项保持原状。此试验的第一份 SIMPLE 探针因 weak inline 符号使用原头文件而排除；正式比较重新编译所有调用单元，场与映射使用一致的实验头文件。

四个严格 field solve 使用原 6,208 格、`nu=.1`、`speed=.5`、inertia=1；映射探针单独设 inertia=0，在相同 SIMPLE 松弛和四次压力修正下执行 180 次。观测模式和残差单列，不称完整谱界。`open-source-mechanisms-fields.tar.gz` 的消费者是后续原生机制研究：保留场、通量、压力向量、分解面贡献、完整命令/成本与原创补丁；产品运行仍用前述 coupled JSON 算例。四个试改都不属于产品修复，不能改默认或据此继续高密度链。

### 原生二维开源内核实际对照

`tools/benchmarks/freefem_reference.py` 实际调用 **FreeFEM 4.15 的二维有限元装配及 UMFPACK**，只用 Python 编排和记录哈希，不重写有限体积方程。入口限于已保留的圆环：中心 `(0,0)`、内外多边形外接圆半径 `.5/1`、内壁速度 `.5`、`nu=.1`。原生转换器验证原拓扑、边界身份及基准参数；以原 FV 中心作扇形三角剖分，保留全部多边形、原边及物理边界，任何非正三角形都失败。它是独立研究入口，未进入产品求解或默认构建。

上游参考固定为 [FreeFEM v4.15](https://github.com/FreeFem/FreeFem-sources/tree/b1e524c8ca6a9b44cb9eff780459bacfc12a3ad7)，使用官方 Newton 弱式与原生有限元空间。实际运行包为 Debian `4.15+dfsg-1` amd64；包 URL、SHA-256、二进制哈希、许可证来源和参考文件版本见 `artifacts/current/laminar-foundation/open-source-core.json`。FreeFEM 主体为 LGPL-2.1，分发包内其他组件有各自声明；仓库不提交第三方二进制或整份源码。Linux 本轮环境可按下列方式解包运行，其他系统需自行安装同版 FreeFEM，跨平台未验证：

```sh
mkdir -p outputs/freefem-runtime/packages
for package in freefem++ libfreefem++; do
  curl -fL "https://deb.debian.org/debian/pool/main/f/freefem++/${package}_4.15+dfsg-1_amd64.deb" -o "outputs/freefem-runtime/packages/${package}.deb"
  dpkg-deb -x "outputs/freefem-runtime/packages/${package}.deb" outputs/freefem-runtime
done
cmake -S . -B build
cmake --build build --target cartmesh2d_freefem_mesh --parallel 2
python3 tools/benchmarks/freefem_reference.py --mesh outputs/laminar-inputs/annulus-l7.solver.cm2d --boundary outputs/laminar-inputs/annulus-l7.boundaries --output outputs/freefem-reference/cr-l7-ns --freefem outputs/freefem-runtime/usr/bin/FreeFem++ --space cr --trace facets --stress laplacian --inertia 1
python3 tools/benchmarks/freefem_reference.py --mesh outputs/laminar-inputs/annulus-l7.solver.cm2d --boundary outputs/laminar-inputs/annulus-l7.boundaries --output outputs/freefem-reference/th-l7-exact-ns --freefem outputs/freefem-runtime/usr/bin/FreeFem++ --space th --trace exact --stress symmetric --inertia 1
```

原输入先按本节前面的 `failure-inputs.tar.gz` 命令解包。替换文件名中的 `l7` 可运行 `l5/l6`；`--inertia 0` 是 Stokes 对照。输出前缀已存在时明确拒绝覆盖，换新前缀保留旧场与失败记录。原生 CLI 和 App 不依赖 FreeFEM。

- **CR/P0 + facets + Laplacian**：非协调速度 DOF 位于三角边中点，线性边上的中点值也是面均值，故保留原规定逐面速度及零平均法向通量。原折线转角两侧速度不连续，不能把两个不同值同时赋给 Taylor–Hood 的共享角点；封装明确拒绝 TH/facets。CR 内部是 broken Laplacian，不能把它和原生对称应力的有限网格结果视为同一个离散方程，也不把边均值无穿透说成边上逐点无穿透。
- **Taylor–Hood P2/P1 + exact + symmetric**：在原多边形边界施加圆形 Couette 解析速度的 P2 节点插值。它允许解析场在弦上的局部法向速度，改变了原逐面不穿透的物理边界，用来验证参考解的一致性。`originalFacetMidpointVelocityDifference` 是与原边界的差；`dirichletDofVelocityMismatch` 才是实际所施加边界的满足误差。
- 压力固定一个 DOF 后减全域均值，不加压力质量惩罚。原完整线性残差和非线性弱残差都检查被删的冗余连续性行。初版全局均值乘子使 UMFPACK 极慢；精确单点基准与已完成两档均值约束解的逐 DOF 速度/去常数压力相符到 `1e-13` 以内，原停止记录及探针保留。
- 收敛要求 Newton 场变化和所有自由 DOF 的弱残差均 `<1e-8`，并检查规定速度 DOF 和压力均值。速度/压力尺度为 `Uref/Uref²`，动量/连续性弱行尺度为 `nu*Uref+inertia*Uref²` / `Uref`；这个有限元判据与原生 FV 的逐单元 strict 判据单列，不能代替空间精度。`cells.csv` 使用 FreeFEM 积分求三角形真实均值再聚合到原单元，原 FV 值则是中心估计。

`open-source-core-fields.tar.gz` 保存代表场、十次最终外部运行的日志/控制、原生 Stokes 与解析 trace 对照驱动、SIMPLE 真实失败、初版接入错误及压力基准对照。消费者是后续边界/离散修复和跨机复跑：原网格复用已有输入包，转换网格可重建，运行库/依赖/编译产物只留在本机 `outputs/`。原生解析 trace 对照将 983 个边界面按真实通量符号设为 velocity-inlet/outlet，一个恰为零通量的非零速度面使用现有 SmoothMovingWall；没有修改角色检查或补偿净通量。这同时改变 trace 和边界重构路径，不能据此宣布唯一壁面系数已经定位。原不支持角色的拒绝记录和后续 SIMPLE 抛错均保留；后者没有返回末场，只有真实进度和错误，不伪造接受场。

归档内原运行命令含本机绝对路径；跨机原生复跑可解包到仓库根目录，直接重编其中 C++ 驱动并传新输出前缀，无需保留库快照或改写旧结果。例如复跑同解析 trace 的 coupled NS：

```sh
tar -xzf artifacts/current/laminar-foundation/open-source-core-fields.tar.gz -C .
mkdir -p outputs/freefem-reference
c++ -std=c++20 -O2 -Iinclude outputs/open-source-core/native-exact-trace-supported/exact-trace.cpp build/libcartmesh2d_fv.a build/libcartmesh2d.a -o outputs/freefem-reference/native-trace
outputs/freefem-reference/native-trace outputs/laminar-inputs/annulus-l7.solver.cm2d outputs/freefem-reference/native-ns 1
```

SIMPLE 失败驱动在 `native-exact-trace-supported-simple/exact-trace.cpp`，同样链接当前原生库，最后两个参数为 `1 simple`。这些驱动的原验证库版本为 `d51d8eb`；本次适配提交没有改变 FV 数值代码，后续修复后应保留新旧结果的版本区别。

### 单元一步分项与来源传播

`tools/benchmarks/cell_step.py` 是 Linux/GNU ld 研究入口，先构建已有层流库，不进入默认构建或产品求解。它从当前 archive 提取真实 C++ 符号，编译包裹版和未包裹版；`cell_step.cpp` 只记录实际 `momentum` / `rhieChowFlux` 输入及返回，先恢复原生 Stokes 一步映射的模态，再抓取下一次全网格 `F(x)` 与 `F(0)`。目标格只是观测位置，压力求解仍全局进行。解包前述原失败输入后运行：

```bash
python3 tools/benchmarks/cell_step.py \
  outputs/laminar-inputs/annulus-l7.solver.cm2d \
  outputs/laminar-inputs/annulus-l7.boundaries \
  --target 3710 --power-steps 180 --output outputs/cell-step
```

仅支持稳态恒黏度 Stokes、闭合固定速度壁；固定 `nu=.1`、参考速度 `.5`、松弛 `.6/.25`、IC(0)、四次压力校正和原严格线性停止。零基目标 3710 对应用户单元 3711，相邻零基 3711 的面积约小十倍。壁面与内部黏性表包含紧致旧场作用、非正交和真实对称应力；压力力以共享面值重构并与实际 RHS 对齐。动量积分量单位 `m³/s²`，RC 是 `m²/s`，不互相相加。逐行重构对齐误差仅记录，不构建新的独立 Python 方程验证链。

`cell_step_response.cpp` 保留捕获的实际紧致矩阵，分别用原生矩阵求解、速度梯度、RC 和 `solvePressureCorrection` 传播壁面/内部非正交与对称应力、动量压力/松弛，以及直接 RC 压力、通量松弛、旧压力，共 11 个来源。modal/base 矩阵和 rAU/df 必须相同、体源须零；最终面通量使用原压力函数返回的最后一次**滞后**修正。速度/压力分别用面积内积，通量用面长内积，投影各除以对应输入场平方范数；这些是单个模态的贡献，不是删项后的新映射谱。分源 Krylov 停止带来有限闭合误差，不宣称逐位闭合。

`cell-step-diagnosis.json` 保存验证和实际数字，`cell-step-diagnosis-traces.tar.gz` 保留完整模态输入/一步输出、实际共享矩阵、目标及邻格逐面项、阶段和来源响应及可重建的原生源文件。重复全格 CSV 视图留在 `outputs/`；归档的选择及原 SHA 明列。消费者是后续应力—压力反馈修复；单例定位不扩大产品稳定性或精度资格。入口本身再次复跑得到同一字段哈希与来源闭合，相关四项原生测试通过；src/include、默认设置及生产 archive 未变。

### 与网格产品一致的求解规模

云端继续修复原生二维不可压层流基础；上述 level 5/6/7 矩阵是问题定位入口，最终交付须覆盖产品实际密度。数量依据是 `desktop/src/core/cell-budget.js` 的 50 万目标上限、`capabilities.js` 的 Cut-cell 壁面 level 11 / 共形边界层 level 8 历史边界，以及已有 102,017 格曲壁收敛记录。这些是规模目标的来源，不是新的误差容差或已取得的普适资格。

- **逐级求解并检查后再加密**：从现有约 3 千/6 千格失败输入开始，先完成稳定性修复，再按约 1 万、3 万、10 万、30 万、50 万格逐步推进。按实际空间细化选择邻接档位，必要时增加中间档，不强凑这些格数，也不从小例直接跳到十万/五十万格。每档沿用同一物理问题和可比较的数值控制；新增算法或控制先回查较粗档，再继续细化。不要预先并发启动同一序列所有大规模求解。
- 进入下一档前读取当前完整结果，检查收敛、守恒、场的物理合理性，以及相邻档误差、迭代数、内存和总耗时的变化。粗档保留其真实离散误差，不要求预先达到最终细档精度；通过相邻档判断改善趋势。发散或明显非物理解先修复，误差异常停滞/恶化须定位参考采样、几何近似或离散问题后再继续；高密度目标保留，不以失败为由永久降回小网格交付。
- 使用当前产品网格流程，在固定几何、计算域、Re、物性和边界条件下完成约 10 万、30 万、50 万**真实流体单元**的曲壁系列；覆盖纯 Cut-cell 和共形边界层网格。500,000 是生成器目标数量，报告聚合后的实际格数，不通过凑格数改变流体区域。30 万至 50 万目前仍有规模与收敛缺口，须继续算法修复并完成实际求解；短步预算探测或仅生成大网格不算完成。
- 每档保存生成参数、输入几何及网格哈希、实际流体格数、层级分布、最小单元尺度及壁面切向/法向尺寸。level 由计算域跨度和目标尺寸共同决定；同一个 level 不保证同一物理分辨率。边界层按实际法向尺寸比较，完整背景中的固体内部格不计入流体规模。曲壁细化同时记录几何离散误差，不能只加密远离壁面的区域来完成精度验证。
- 标准基准继续用解析解/公开参考值判断速度、压力及受力误差；高密度代表曲壁还须有同网格 OpenFOAM 对照和保持同一物理问题的细化比较。10/30/50 万是容量覆盖点，不直接当等比例的三档网格计算阶数；观测阶依据受控的实际空间细化。修复须在高密度保持已验证的精度与守恒，失败样本保留并进入后续修复。
- 同时记录外迭代数、线性求解工作量、残差与场变化、压降/受力、峰值内存和完整墙钟时间。冷启动单列；使用粗解初始化时，把粗解、映射、细解及导出全部计入，网格生成另列并给出总成本。同硬件、同控制和同物理精度比较，Linux 使用实际可用的预条件器。大规模收敛困难须解决离散、耦合或线性求解的瓶颈，不能仅放宽停止条件或把迭代上限当收敛。

历史规模证据入口为 `artifacts/current/native-laminar-gradient-cache.json`、`native-laminar-initial-guess.json` 和 `native-laminar-dense-cutcell.json`；旧命令中的本机 `outputs/` 文件不会随 Git 上传，云端须保存可重建的输入/生成命令。历史喷管几何为 `examples/complex/nozzle_profile.xy`；它的指定网格性能结果不替代圆环/DFG 精度验收。新高密度运行使用独立 `outputs/` 目录，保留后续续算输入、代表场和最小失败例，只提交必要的小体积证据。

## 产物保留与清理

维护对象是当前源码与必要测试；App 打包只包含 `desktop/src/`、原生 runtime 和样例。多平台共用同一算法源码，构建树、依赖、研究数据及安装包分别管理。[CMake](https://cmake.org/cmake/help/latest/manual/cmake.1.html#introduction-to-cmake-buildsystems)明确区分源码树和可重建的构建树。

- `build/` 的编译产物和日常 CTest 临时结果可重建；其中 `deps/`、`chemistry-env/`，以及 `outputs/topology-env/` 是现有可选研究环境，须先核对消费者，不能随整目录清掉。
- 新任务收尾时记录哪些文件将被续算、初值准备或对照直接读取；保留这些输入、代表场、参数/命令和最小失败例。重复导出与短时诊断不因“验证做过”无限保留，历史原始数据的删除仍按已有保留要求逐项决定；没有文本引用并不证明数据无用。
- `desktop/dist/` 按平台/架构保留最近两版常规安装包，另保留当前源码版本和 `mesher-v0.3.0` 对应包；App、runtime、未知命名文件不参与。`npm --prefix desktop run clean:packages` 先预览，加 `-- --apply` 执行；三平台成功打包会自动执行，失败打包不会触发。
- CI 普通安装包保留 7 天、构建/测试日志保留 14 天，正式发布包放 Releases。这是本项目的回滚与排错保留窗口，使用 [GitHub artifact retention](https://docs.github.com/en/actions/tutorials/store-and-share-data#configuring-a-custom-artifact-retention-period) 配置，提交后的新运行才生效。
- 当前状态只写关键结果、限制和必要证据入口，方法/复现留在本指南；逐轮叙述不再累计进当前状态。本次整理收据及压缩的旧状态段落见 `outputs/repository-tidy/storage-cleanup.json`、`retired-state-sections.json.gz`，不作为新的研发状态入口。

## 分支与里程碑

只维护 `main`（已验证集成）和 `codex/cfd-development`（后续开发），网格和 CFD 修复统一从开发线推进。三条实验功能、临时集成及旧网格维护分支已完整进入 main，旧分支名可删除，提交历史仍可查。

暂停分支若含未合入提交，先以 `archive/2026-10-08/<原分支名>` 保存远端 annotated tag，并核对 peeled commit，才删除原分支引用。本次 6 条暂停分支的精确 SHA 和独有提交数见 `artifacts/current/laminar-foundation/targeted-repair.json`；恢复可用 `git switch -c <恢复分支名> archive/2026-10-08/<原分支名>`。归档不等于已集成或精度通过。共同核重构 PR #10 已快进进入开发线，重复分支已删除；独立云环境 PR #9 及其分支保留。

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

### Windows PC 开发与 CUDA 接入

PC 开发从 `codex/cfd-development` 继续，不能只依据 main 的 0.4.41 首页判断研究进度。当前核心任务是不可压层流的压力—速度耦合、曲壁压力与受力误差；燃烧、可压、SST 和浸入边界扩展保持冻结。默认 SIMPLE、SIMPLER/PIMPLE、完整响应及 Anderson 的实际资格见当前状态。GPU 后端应复用这些离散方程和失败判定。

Git 拉取的是源码、版本化样例及 artifacts 中的证据；不会恢复旧 Mac/Linux 的 `outputs/`、Python 环境、已编译 CLI、Electron 依赖或安装包。旧文档中的绝对路径和“本机”性能指原实验机器。需要重跑原失败例时，可从 `artifacts/current/laminar-foundation/failure-inputs.tar.gz` 按需解包；Windows 多配置构建的 CLI 在 `build/Release/*.exe`，研究脚本须显式传入实际二进制路径，不能照搬 POSIX 默认路径。

本次 PC 检查找到 VS Community 2026、MSVC 14.51.36231、Windows SDK 10.0.26100.0，以及 VS 自带 CMake 4.3.1 和 Node 24.12.0。CMake/Node 未在普通 PATH 中；项目与 CI 的 Node 基线仍为 22。以下只在当前 PowerShell 会话补入已安装工具，不修改系统 PATH；其他机器优先使用已安装的 Node 22：

```powershell
$vswhere = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
$vsRoot = & $vswhere -latest -products '*' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
if (-not $vsRoot) { throw '未找到 Visual Studio C++ 工具链' }
$cmakeBin = Join-Path $vsRoot 'Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin'
$env:PATH = "$cmakeBin;$env:PATH"
if (-not (Get-Command node -ErrorAction SilentlyContinue)) {
    $nodeBin = Join-Path $vsRoot 'MSBuild\Microsoft\VisualStudio\NodeJs'
    $env:PATH = "$nodeBin;$env:PATH"
}
cmake --version
node --version
npm.cmd ci --prefix desktop
npm.cmd --prefix desktop run build:native
npm.cmd --prefix desktop start
```

`build:native` 实际编译本机 C++ 并准备 `desktop/runtime/`；桌面通过 IPC/子进程调用这些 CLI，提交网格或流动任务时才进行对应本地计算。原生开发可单独使用 `cmake -S . -B build -A x64 -DCARTMESH2D_BUILD_CHEMISTRY=OFF`，再 `cmake --build build --config Release --target cartmesh2d_cli cartmesh2d_flow_cli --parallel 2`。同一个 build 目录保持同一生成器；Windows 原生工具不支持 macOS Accelerate Cholesky，使用现有 IC(0)/Jacobi/aggregation，不能直接重用 `--pressure-preconditioner cholesky`。

**当前已实现显式可选的 Windows 原生 Torch-CUDA 压力模块。** GPU 为 RTX 4060 Laptop、8188 MiB、计算能力 8.9，驱动 572.70；`nvidia-smi` 的 CUDA 12.8 表示驱动报告的兼容能力。直接复用 `D:\gas-turbine-envs\windows-runtime\Lib\site-packages\torch` 的 PyTorch 2.12.1+cu126 C++ 头文件、导入库及已有 DLL，不复制依赖，也不启动 Python 工作进程。现有 VS 2026 已实际编译运行；这条路径调用预编译 GPU 算子，无需 nvcc。CPU 保持默认且不加载 Torch，GPU 模块不进入当前常规安装包。

所查 Windows 环境/常见 Toolkit 路径未找到 nvcc；WSL 的 `/home/lenovo/combustor-envs/linux-runtime` 是 cu126，另一环境 `qwen3-cuda` 是 cu129，Linux 库不能直接当 Windows DLL 使用。[PyTorch 官方说明](https://docs.pytorch.org/docs/stable/cpp_extension.html)中的自定义 CUDA 扩展编译工具要求，适用于后续编写新 `.cu` 内核；不能因此判断既有预编译 GPU 库无法复用。[CUDA 12.8 Windows 安装指南](https://docs.nvidia.com/cuda/archive/12.8.0/cuda-installation-guide-microsoft-windows/index.html)的 VS 2022/MSVC 193x 支持表也不否定本次已实测的普通 C++/ATen 路径。后续需要 nvcc 时再匹配版本。

构建与运行使用同一套现有 Torch 安装；仅在当前 PowerShell 会话设置路径：

```powershell
$env:CARTMESH2D_TORCH_ROOT = 'D:\gas-turbine-envs\windows-runtime\Lib\site-packages\torch'
cmake -S . -B build -DCARTMESH2D_BUILD_TORCH_PRESSURE=ON "-DCARTMESH2D_TORCH_ROOT=$env:CARTMESH2D_TORCH_ROOT"
cmake --build build --config Release --target cartmesh2d_torch_pressure cartmesh2d_flow_cli cartmesh2d_pressure_replay cartmesh2d_flow_gpu_pressure_tests --parallel 2
ctest --test-dir build -C Release -R '^cartmesh2d_flow_gpu_pressure$' --output-on-failure
# 先按下方小例生成 rectangle.solver.cm2d，再显式选择 GPU：
build/Release/cartmesh2d_flow_cli.exe --mesh outputs/pc-smoke/rectangle.solver.cm2d --output outputs/pc-smoke/channel-gpu --case channel --nu .1 --speed .1 --tolerance 1e-6 --profile --pressure-backend torch-cuda --pressure-preconditioner ic0 --pressure-snapshot outputs/pc-smoke/channel.pressure.bin
```

`CARTMESH2D_BUILD_TORCH_PRESSURE` 默认 OFF，不影响纯 CPU 构建。模块默认从 CLI 同目录加载，可用 `CARTMESH2D_PRESSURE_PLUGIN` 显式给出其他 DLL 路径。模块隔离为 `/MD` C ABI，主程序继续静态 CRT；只验证 Windows x64 Release。GPU 仅支持稳态恒物性层流的 Jacobi/IC(0)，不支持 GPU aggregation、瞬态或材料反馈，缺少模块或不支持组合明确失败，不自动转 CPU。IC(0) 因子仍由原生 CPU 生成；向量、稀疏乘法和两次三角求解在 GPU 执行。返回的双精度场须通过原生 `b-Ax` 与原门槛 `1e-13+relativeTolerance*||b||₂`，不是另设较宽门。

| 代码入口 | 当前责任与后续边界 |
| --- | --- |
| `include/cartmesh2d/fv/detail/FlowLinearSystem2D.hpp` | 优先评估压力 PCG、稀疏乘法和向量运算。当前是固定排序的非对角 CSR 加独立 diag；对接完整 CSR 时须显式插入对角并检查索引范围。 |
| `src/fv/FlowPressure2D.hpp`、`Incompressible2D.cpp` | 保留非正交压力循环、压力参考格、真实残差和原生停止门；计算后端与预条件器分开选择。 |
| `src/fv/FlowCoupled2D.cpp`、`FlowKrylov2D.hpp` | 完整响应使用矩阵自由块 FGMRES、ILU(0) 动量和嵌套压力解。仅把压力乘法搬到 GPU 不等于整个耦合求解已经 GPU 化。 |
| `FlowEquation2D.cpp`、`FlowMomentum2D.cpp`、`FlowFlux2D.cpp` | 根据同机 profile 再评估组装、梯度和面通量；保持共享压力力、守恒与接受步恢复。 |
| `FlowGpuPressure2D.cpp`、`FlowGpuPressureAbi2D.hpp`、`TorchPressureModule2D.cpp` | 动态加载、C ABI、矩阵/因子缓存及 GPU PCG；失败不返回候选，缓存不能在失败后被当作有效新矩阵。 |

`--profile` 增加 GPU 实际调用/迭代、矩阵上传、初始化、上传、迭代及下载计时。它们是原压力计时的组成部分；coupled 的压力时间又嵌套在线性时间内，不能重复相加。GPU 归约顺序不同，不承诺逐位一致；同控制最终原生门和场差分别核对。当前小格完整 CFD 更慢；大格同算法 Jacobi 可提速，但仍未胜过该组 CPU 多重网格，不能自动改默认。实际结论、失败与未验范围只看当前状态及[实测证据](../artifacts/current/simple-gpu-pc.json)。

`--pressure-snapshot FILE` 只捕获第一次非零原生压力系统（`CM2DPRS1`，保存 CSR、独立 diag、rhs、初值及实际线性控制）；不是物理检查点。`cartmesh2d_pressure_replay` 调用生产原生/GPU 求解器，比较 CPU IC(0)/Jacobi/aggregation、GPU IC(0)/Jacobi，同原门槛检查真残差，报告一次冷解与三次矩阵缓存后重复。`--channel` 先经原 Solver 质量门生成规则流体格，再跑一步真实流动装配捕获压力系统；一步未收敛不能当作完整大 CFD 资格。

```powershell
build/Release/cartmesh2d_pressure_replay.exe --snapshot outputs/pc-smoke/channel.pressure.bin --gpu --repeats 3 --output outputs/pc-smoke/pressure-replay.json
build/Release/cartmesh2d_pressure_replay.exe --channel --nx 320 --ny 320 --gpu --repeats 3 --output outputs/pc-smoke/pressure-102400.json
```

SIMPLE 停止量和 coupled 成本的研究入口为 `tools/benchmarks/simple_convergence_study.py`；它只编排原生 CLI、读取接受场/残差/profile 并计算场差，不重建方程或恢复独立 Python 审计链。先解包 `failure-inputs.tar.gz` 到 `outputs/laminar-inputs/`。`simple_source_compare.py --build-only` 支持 Windows x64 开发者终端中的 `cl.exe`，保留 `/MT`、UTF-8 manifest 和原数值源码选择；每次输出目录须新建。Gauss 研究变体不覆盖产品默认。

```powershell
# 下行在 VS x64 开发者终端运行，cl.exe 必须可见：
& D:\gas-turbine-envs\windows-runtime\python.exe -B tools/benchmarks/simple_source_compare.py --mesh outputs/laminar-inputs/dfg20-l7.solver.cm2d --boundary outputs/laminar-inputs/dfg20-l7.boundaries --nu .001 --speed .3 --output outputs/simple-pc-variants --build-only
& D:\gas-turbine-envs\windows-runtime\python.exe -B tools/benchmarks/simple_convergence_study.py --suite simple --output outputs/simple-pc-study
& D:\gas-turbine-envs\windows-runtime\python.exe -B tools/benchmarks/simple_convergence_study.py --suite gpu --output outputs/gpu-cfd-study
& D:\gas-turbine-envs\windows-runtime\python.exe -B tools/visualization/render_solver_gpu_study.py
```

本机原生接受场、严格停止参照、失败及 frozen matrix 在上述 outputs 目录，消费者是下一轮同离散停止/线性后端比较；curated JSON/两图/HTML 在 `artifacts/current/simple-gpu-pc*`。大矩阵只验线性子问题，不恢复已暂停的大规模曲壁精度链；更高密度不外推资格。`tests/flow_gpu_pressure_test.cpp` 覆盖原真残差门、RHS/系数缓存更新、未知预条件拒绝与预算失败后恢复。

可参考 [PETSc GPU 矩阵/向量与求解器后端](https://petsc.org/release/overview/gpu_roadmap/)及 [Firedrake OffloadPC](https://www.firedrakeproject.org/_modules/firedrake/preconditioners/offload.html)的分层方式：离散方程、线性后端与预条件器分开，GPU 内部运算和 CPU/GPU 传输都计成本。Fluent 2026 R1 有[原生 GPU 求解器](https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/flu_ug/flu_ug_sec_gpu_solver_supported_features.html)，该文档限定三维几何；不能把既有二维 Fluent CPU 连接检查当成 GPU 实算。规则笛卡尔内区有利于规则邻点运算，但本项目的 Quadtree/Cut-cell/聚合与共形边界层还含非规则连接；完整背景格仍为几何产品。GPU 降低单轮成本，预条件/多重网格减少迭代次数，曲壁离散及压力—速度耦合仍须分别修复。

本机环境基本检查的小例命令如下；原生测试须先构建同名 `_tests` 目标，所有命令均在仓库根目录运行：

```powershell
node --test desktop/tests/process.test.js desktop/tests/flow-case.test.js
cmake --build build --config Release --target cartmesh2d_flow_linear_tests cartmesh2d_flow_equation_tests --parallel 2
ctest --test-dir build -C Release -R '^cartmesh2d_flow_(linear|equation)$' --output-on-failure
New-Item -ItemType Directory -Path outputs/pc-smoke -Force | Out-Null
build/Release/cartmesh2d_cli.exe examples/acceptance/rectangle.xy outputs/pc-smoke/rectangle 4 0.25 0.10 interior outputs/pc-smoke/rectangle-foam 4
build/Release/cartmesh2d_flow_cli.exe --mesh outputs/pc-smoke/rectangle.solver.cm2d --output outputs/pc-smoke/channel --case channel --nu .1 --speed .1 --max-iterations 2000 --tolerance 1e-6 --profile
desktop/node_modules/electron/dist/electron.exe desktop --smoke=circle --out=outputs/pc-smoke/app --shot=outputs/pc-smoke/app.png
```

Solver 输入必须由正式质量/导出流程生成 `.solver.cm2d`，不能直接将 construction `.cm2d` 改名后送入 CFD。实际检查范围见当前状态；`outputs/pc-smoke/` 是可重建的环境小例，不替代已有研究输入或物理精度证据。

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

### T01 双管原始输入与冻结版复算

[t01-originals.zip](../artifacts/current/laminar-foundation/t01-originals.zip) 是用户上传原包的逐字节副本，769,757 字节，SHA-256 为 `bfd6a084e6c9763bc045cb48ccbe7dec748a5e1736bba5b832d9bc8db9d9b7f0`。其 40 个文件包含 Re=0/50 原候选、`reference-geometricSeed.npz`、`reference-engineering.npz`、XY 和原始摘要，解压到仓库根目录后恢复 `outputs/topology-fidelity/b-retries/double-pipe-re-{0,50}/`。它是后续重算消费者的原始材料，保留原 ZIP、NPZ 和 XY；已有原始目录无需重复解压。

`tools/optimization/t01_double_pipe_recompute.py` 复用现有 Brinkman 评分器、`native_flow.py` 和原生求解器，没有新增方程审计链或重新优化设计。`prepare` 核对 ZIP 内全部文件、重评六份固定材料场、冻结 mesher，并从端口中心和宽度独立生成每管四个角点的矩形；`native` 每次只求一个 Re 和一档网格；`report` 核对求解二进制、网格哈希与数值控制，用现有面 CSV 后处理重算功率，输出 JSON、原生/多孔 CSV、排序/成本图及真实边界图。

本次产物为 **修复前冻结版本（frozen pre-repair）证据**：[汇总](../artifacts/current/t01-frozen-double-pipe.json)、[原生表](../artifacts/current/t01-frozen-double-pipe.csv)、[多孔表](../artifacts/current/t01-frozen-double-pipe-porous.csv)。使用的 flow-cli SHA-256 为 `6fdf07218b9c8f584e1e974691c8a57253fd4b539bc1990a7e95daa463e7b674`，mesher 为 `bf6221294b35acc0cada0b2083f098e142b6d6b3c90a8c6d8fed8edc77835551`。原 ZIP 不含二进制；复现同一冻结版本须使用上述匹配快照，不能把新编译的求解器当作同一记录。修复后使用另一固定快照和新输出目录，候选、工程参照及独立矩形必须全部重算。当前结果、限制和待办见[当前状态](CURRENT_STATE_CN.md#当前研究与接手入口)，不授予修复后物理精度或可靠排序资格。

以下在已有 NumPy/SciPy/Matplotlib/contourpy 的 Python 环境中执行；云环境可先 `source /workspace/cartmesh2d-cloud-setup/activate.sh`。新工作区先解压原包，再准备独立输入；`T01_FLOW` 必须指向本次比较使用的固定求解器：

```sh
python -m zipfile -e artifacts/current/laminar-foundation/t01-originals.zip .
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
T01_FLOW=outputs/laminar-repair/efficiency-curved-fixed/flow-cli
T01_MESH=outputs/t01-agent/bin/mesh-cli
T01_OUT=outputs/t01-recompute-new
python tools/optimization/t01_double_pipe_recompute.py prepare \
  --output "$T01_OUT/independent" --flow-cli "$T01_FLOW" --mesh-cli "$T01_MESH"
```

以下是**一个档位**的完整计算块。先做 Re=0、level6，读取候选、工程参照及独立矩形的 `summary.json`、收敛/守恒、场和成本，再改为 level7；随后按同样顺序做 Re=50。不要一次并发启动四个档位。level6 使用 `T01_FINE=`，level7 使用 `T01_FINE=-fine`；原生严格停止 `1e-8`、峰值速度 `.02`、松弛 `.2`、一遍压力校正、Anderson 默认 4/10、aggregation/adaptive、Laplacian 黏性和 limited-linear 对流与本次冻结对照一致。

```sh
T01_RE=0
T01_LEVEL=6
T01_FINE=
python tools/optimization/compare_sharp_designs.py \
  "outputs/topology-fidelity/b-retries/double-pipe-re-$T01_RE" \
  --baseline reference-engineering.npz --area-target candidate \
  --output "$T01_OUT/comparison/t01-re${T01_RE}-engineering-frozen${T01_FINE}" \
  --levels "$T01_LEVEL" --reynolds "$T01_RE" --speed .02 --small-alpha .1 \
  --iterations 16000 --timeout 1800 --tolerance 1e-8 \
  --mesh-cli "$T01_OUT/independent/bin/mesh-cli" --flow-cli "$T01_FLOW"
python tools/optimization/t01_double_pipe_recompute.py native \
  --output "$T01_OUT/independent" --flow-cli "$T01_FLOW" \
  --reynolds "$T01_RE" --level "$T01_LEVEL" --iterations 16000 --timeout 1800
```

完成四个档位后汇总，失败档位保留并先处理，不用解析值补全。`report` 使用独立矩形的实算 J 归一化原生结果，多孔结果则使用同模型的原几何种子 J；工程轮廓投影到面积 `.5`，该投影仍是多孔与锐壁差距的一部分。

```sh
python tools/optimization/t01_double_pipe_recompute.py report \
  --output "$T01_OUT/independent" --flow-cli "$T01_FLOW" \
  --comparison-root "$T01_OUT/comparison" --comparison-suffix engineering-frozen \
  --artifact "$T01_OUT/t01-frozen-double-pipe"
```

真实网格、场和日志在本次 `outputs/t01-agent/` 与既有 `outputs/laminar-repair/t01-re{0,50}-engineering-before[-fine]/`；新运行使用上面的新目录。原始输入包供继续求解，JSON/CSV/两图供结论复核，构建快照和大网格留在本地。成本区分原生子进程、独立求解完整编排和输入准备；多孔重评是准备时间的子集，不重复累加。旧候选/工程记录没有完整 Python 编排计时，解释器启动及报告绘图也不在计时区间，报告显式保留这些限制。这些运行全部从冷启动开始，未使用粗解映射；同 level 的工程合流管与双矩形实际格数/壁面尺寸不同，须读报告中的网格域、实际格数和分辨率，不能据 level 数字声称等密度。

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
