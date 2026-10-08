# CartMesh2D 开发规则

## 接手与维护入口

- 先核对 `git status`、`git branch --show-current`、`git worktree list`，再读 `docs/CURRENT_STATE_CN.md`。当前根目录是唯一开发入口，前端只看 `desktop/src/`；旧工作区仅保留研究证据。
- 能力、版本、验证和待办只维护 `docs/CURRENT_STATE_CN.md`；代码导航与复现在 `docs/DEVELOPMENT_CN.md`，操作在 `docs/DESKTOP_APP_CN.md`。不新增日期审计、阶段计划或交接文档。
- `main` 是已验证集成线，`codex/cfd-development` 是后续开发线，网格和 CFD 修复统一在开发线推进。已完成分支的提交由 Git 历史保留，固定里程碑标签不移动；本轮合并授权不代表后续功能无限期自动合并。
- 主线维护网格完整性、层流精度/默认稳定性及高密度提速。可压层流、浸入边界、拓扑优化按开发版/实验版/研究工具范围维护，物理扩展须有明确任务；软件集成不扩大精度资格，不自动改默认控制。
- 重要突破验证后 commit 并 push，报告分支、提交和实际验证范围。使用用户配置的 Git 身份，不加 AI/bot Co-authored-by；贡献者显示先查历史和 GitHub API，不为刷新页面改写历史。

## 原生二维与流体语义

- 不依赖、链接或复制三维 `cartmesh` 核心，不把三维算法压到 z=0。核心保持 Point2D、Segment2D、AABB2D、Polygon2D、Quadtree、CutPolygon、Edge2D。
- 闭合 BoundaryLoop 默认是固体，Domain2D 是外域，外流是 `domain - solid interior`。Inside 不进入流体网格，Outside 是流体，Intersected 保留真实外侧 polygon；内流只允许显式 `FluidRegion2D::Interior`。不得以中心采样、删相交格代替 Cut-cell。
- 外流同时具有 EmbeddedBoundary 与 DomainBoundary，且 `fluid_area = domain_area - solid_area` 在容差内成立。自交、零面积、重复/孤立边、非流形及分类冲突显式失败；现有漏检是缺陷，不是放宽依据。
- 完整笛卡尔背景是独立几何产品，保留并分类 Inside/Outside/Intersected 全单元，不标记为流体求解网格，不直接接入绕流 CFD。
- 浸入边界使用独立计算格、命名空间和输出，明确标识固体辅助未知量，保留真实几何面积、有限阻力泄漏及壁面误差；完整方格不冒充共形流体拓扑，数值收敛不冒充物理精度。

## 验证范围

- 日常只做相关构建、运行、交互和真实数据的基本检查；不自动扩展为全量回归、网格收敛、跨平台或极限性能研究。开发可用、数值收敛、物理精度分别报告，保留有价值的未完成方案。
- 明确精度验收、重要交付或 main 功能集成时按目标严格验证：完整原生 CTest、前端测试、真实 App 和对应平台 CI。重要网格里程碑还须真实网格/预览、独立读取器；OpenFOAM 可用时实际运行 checkMesh。
- 新增或收紧阈值前说明检查量、单位/归一化、依据、任务必要性及成本；`1e-9` 等不是通用标准。超出任务的精度或提速研究先说明收益和成本，由用户决定。
- NaN、非法几何、坏拓扑、失效控件必须报告；达到迭代上限不算收敛。算法修复保留最小失败例，未验项目明确写未验，不能用合并替代资格证明。

## 几何、质量与修改

- 几何拓扑高于图片；不得降低 Solver 质量门、删除坏单元或隐藏拓扑/导出告警。Q1 已取消，不重新加入产品验收；历史 Q1 只作旧记录。
- 边界层、终止区、余域进入统一共形 owner/neighbour 拓扑；降层和 pure fallback 明确告知。原始 XY 折线不为改善外观而平滑。
- 共享格点/交点保留稳定身份与物面来源，不每格补洞或凭近邻猜物面；端点舍入预算与几何容差分开，不能折叠可解析短边。Q3/Q4/Q5 和 patch 事务仍参与构建，不按阶段数字删除。
- tolerance 集中管理；同输入、参数和工具链须产生确定性拓扑/输出，跨编译器一致性另证。拓扑、Solver 质量、标准/扩展 checkMesh 分开报告，不能互相替代。
- macOS 固定系统 `/usr/bin/clang++`，避免 PATH 中 mesasdk 动态库使 App 无法运行。

## CFD 与交付

- 动量、连续性、线性残差、场变化和物理精度分别验收。初值不是物理检查点；失败候选不替代最后接受步，非定常续算核对网格、边界、物性和时钟。
- 提速比较同控制、同精度的完整成本，包含粗解和映射；单例不外推通用性能。研究修复须通用实现与回归，不因单例成功自动改默认。
- Euler、SST、温度等仍有调用或测试的功能继续保留。未解决的精度和稳定性问题集中写入当前状态，不按“本轮不做”删除功能或验收门。

## 目录与证据

- 常规构建只用 `build/`，实验用忽略提交的 `outputs/` 或临时目录；不建长期隐藏工作区、baseline-source、整份源码副本或日期命名试改应用。
- 文档固定为 README、本规则及 docs 三份；展示素材索引合入开发指南，当前状态只放关键证据和待办，过程交给 Git。删除入口时同步修复引用。
- 少量证据放 artifacts，大网格/安装包保留本地或发布附件，不提交缓存、依赖或大批网格。清理分支不删除原始场、失败记录、复现材料及用户指定保留的参考包。

## Cursor Cloud specific instructions

- 构建、测试和桌面启动命令以 `docs/DEVELOPMENT_CN.md` 为准。Cloud Agent 镜像里 `/usr/bin/c++` 是 Clang 18，会选用 GCC 14 的 libstdc++；未安装 `libstdc++-14-dev` 时 CMake 在链接测试程序时报 `cannot find -lstdc++`。系统 `g++`（GCC 13）本身可以链接。环境安装会补上 `libstdc++-14-dev`。
- 没有需要常驻的服务。桌面是 Electron：先完成 `npm --prefix desktop run build:native`，再用 `npm --prefix desktop start`。本机图形会话里若沙箱或 GPU 初始化失败，按 Linux CI 同样加上 `--no-sandbox --disable-gpu`。
- `build/`、`desktop/node_modules/`、`desktop/runtime/` 由安装生成且不提交。OpenFOAM 与 `tools/optimization/` 虚拟环境不是默认开发环境的一部分。
