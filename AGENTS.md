# 项目开发规则

1. 每次只推进一个项目阶段。在阶段 0 验证记录完成前，当前工作范围始终限定为阶段 0。
2. 几何和拓扑正确性优先于可视化或界面开发。
3. 不得把单元中心采样体素描述成 Cut-cell 或可供求解器使用的边界几何。
4. 不得隐藏无效几何、负体积、非闭合或分类冲突。
5. 任何性能结论必须同时记录墙钟时间、峰值 RSS、线程数、硬件和构建类型。
6. 每个里程碑除项目自身测试外，还必须使用独立外部读取器或检查器验证。
7. 未经后续阶段决策，不得加入 CFD 求解器、GUI、云服务、AI 生成网格或部署路径。
8. 相同输入必须保持确定性的遍历顺序、ID、报告和输出。
9. 修复任何几何缺陷时，都必须保留最小失败案例。
10. 不得把几何或项目数据上传到外部服务。

## 当前执行范围（2026-08-12）

- 新环境或新 Codex 会话必须先完整阅读 `WINDOWS_CODEX_HANDOFF_CN.md`，再按其中顺序阅读项目事实来源和计划。
- 当前唯一执行计划是 `STAGE6_REVISED_PLAN_CN.md`。Stage 6.0、Stage 6.1、Stage 6.2 与
  Stage 6.3 已完成；验证记录为 `docs/STAGE6_BASELINE.md`、
  `docs/STAGE6_1_VERIFICATION.md`、`docs/STAGE6_2_VERIFICATION.md` 和
  `docs/STAGE6_3_VERIFICATION.md`。
- 必须按 `6.0 → 6.1 → 6.2 → ...` 顺序一次推进一个子阶段；未经用户明确确认，不得开始修改代码。
- 当前不得自动进入 Stage 6.4；开始前必须先报告 R24 复杂几何质量矩阵、
  当前 `-allGeometry` 的 337 个 concave cells / 8 张 low-weight faces / 24 张
  low-volume-ratio faces、最小失败案例、拟修改文件、停线和验收命令，等待用户确认。
- `STAGE7_HYBRID_CARTESIAN_CODEX_PLAN_CN.md` 只作后续架构约束和接口参考。在 Stage 6 完成并由用户明确确认前，不得实现 boundary layer、prism layer、surface projection、hybrid transition 或其他 Stage 7 功能。
- `docs/STAGE6_PLAN.md` 和 `docs/STAGE6_VERIFICATION.md` 记录的是此前 Stage 6 尝试及其失败证据，不是当前获准继续执行的计划。

## 二维并行子项目例外（2026-08-19）

- 仓库新增 `cartmesh2d/` 原生二维子项目。它与三维 Stage 6/7 并行，但实现必须封闭。
- **只有用户明确要求二维任务时**，允许进入 `cartmesh2d/`；此时 `cartmesh2d/AGENTS.md` 是二维任务的直接执行规则。
- 二维任务不得修改三维 `include/cartmesh/**`、`src/**`、`apps/**`、`tests/**` 核心代码；顶层只允许必要的 CMake 接入和文档链接。
- 三维任务仍严格遵守当前 Stage 6 计划，不得因为二维子项目存在而提前或替换 Stage 7。
- 不得把三维代码 `z=0`、模板化或降维包装后称为“原生二维”。
