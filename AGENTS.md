# CartMesh2D 开发规则

## 接手与唯一入口

- 用户已明确以桌面 0.2.0 / `9bddfbf` 为整理基线。当前根目录是唯一工作区，前端只看 `desktop/src/`。
- 先核对 `git status`、`git branch --show-current` 和 `git worktree list`，再读 `docs/CURRENT_STATE_CN.md`。
- 当前状态只维护该文件；代码导航在 `docs/DEVELOPMENT_CN.md`。不要新增带日期的审计、阶段计划或交接文档。
- 用户优先关心密度、速度、桌面完整性、OpenFOAM 可用性。新问题归入现有目标，不自动无限追加子阶段。
- 重要突破完成验证后 commit 并 push；报告分支、提交和真实验证范围，不提交缓存、依赖或大批生成网格。

## 原生二维与物理硬约束

- 本仓库不依赖、链接或复制三维 `cartmesh` 核心，不得把三维算法压到 z=0 冒充二维。
- 核心保持 Point2D、Segment2D、AABB2D、Polygon2D、Quadtree、CutPolygon、Edge2D。
- 闭合 BoundaryLoop 默认是固体；Domain2D 是外域；流体是 `domain - solid interior`。
- Inside 是固体，不进入流体网格；Outside 是流体；Intersected 必须保留真实外侧 polygon。
- 外流须同时具有 EmbeddedBoundary 与 DomainBoundary，并满足 `fluid_area = domain_area - solid_area`（容差内）。
- 内流只允许显式 `FluidRegion2D::Interior`。不得以单元中心采样或删除相交格子代替 Cut-cell。
- 自交、零面积、重复边、孤立边、非流形、分类冲突必须显式失败；现有漏检属于待修缺陷，不能作为放宽规则的依据。

## 验证与修改

- 几何与拓扑正确性高于图片；不得降低 solver-quality / Q1 阈值、删除坏单元或隐藏告警来通过验收。
- 边界层、终止区和余域必须进入统一共形 owner/neighbour 拓扑；pure fallback 必须明确告知。
- tolerance 集中管理；相同输入参数与工具链产生确定性拓扑和输出。跨编译器一致性必须另证。
- 算法修复保留最小失败案例。日常跑直接相关测试；重要交付完整跑 CTest 与前端测试。
- 重要里程碑生成真实网格和预览，使用独立读取器；OpenFOAM 可用时真实运行 checkMesh。
- 拓扑、solver quality、外部 checkMesh、Q1 分开报告；未运行写未运行，不能互相替代。
- macOS 固定系统 clang++，避免 PATH 中 mesasdk 生成的动态库依赖导致 App 无法运行。

## 防止目录再次膨胀

- 常规构建只用 `build/`；实验输出用忽略提交的 `outputs/` 或临时目录。
- 不在仓库内创建长期隐藏工作区、baseline-source、整份代码副本或日期命名的试改应用。
- 文档固定为 README、本规则及 docs 内三份；更新原文件，历史交给 Git。
- 保留仍参与构建/测试的功能；删除过时入口时同步修正引用，不以减少行数为由删掉验收门。
- 少量证据可进 artifacts；大网格和安装包保留本地或使用发布附件，不塞入日常源码历史。
