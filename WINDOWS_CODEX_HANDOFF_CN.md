# Windows / WSL Codex 接管说明

更新日期：2026-08-12

这份文件是把仓库 ZIP 解压到 Windows 后，新 Codex 会话的第一入口。它说明当前目标、
已知事实、禁止范围和开始工作前必须完成的检查。它不宣称任何尚未通过代码、测试和外部
验证的能力已经实现。

## 1. 事实来源优先级

出现冲突时，按以下顺序判断：

1. 当前项目代码、可复现测试结果、验证产物和 Git 历史；
2. `AGENTS.md` 中的项目约束；
3. `CARTESIAN_MESH_GENERATOR_PROJECT_BRIEF_CN.md`；
4. `STAGE6_REVISED_PLAN_CN.md`；
5. README、旧 Stage 计划和其他说明文档。

文档只说明目标或历史。若文档与实际代码冲突，必须明确报告冲突，不能根据文档假定功能
已经实现。

## 2. 当前唯一工作范围

- 当前只执行 `STAGE6_REVISED_PLAN_CN.md`。
- 必须严格按 `6.0 → 6.1 → 6.2 → ...` 顺序推进，一次只完成一个子阶段。
- Stage 6.0 基线审计和 Stage 6.1 自适应完整 OpenFOAM 拓扑已经完成。
- 当前下一个子阶段是 **Stage 6.2 原生质量评估器**，但未经用户再次确认不得开始修改。
- 在 Stage 6 完成并得到用户明确确认前，不得实现任何 Stage 7 功能。
- `STAGE7_HYBRID_CARTESIAN_CODEX_PLAN_CN.md` 只能作为后续架构约束和接口参考。
- 禁止提前实现 boundary layer、prism layer、surface projection、hybrid transition。
- 未得到用户对下一个子阶段的确认前，不修改 C++、测试、构建脚本或算法实现。

每个子阶段开始前，必须先向用户报告：当前基线、该阶段解决的问题、拟修改文件、接口变化、
明确不修改的代码、验收测试，以及如何证明结果不是 mock、hard-code 或预生成输出。等用户
确认后再修改。

## 3. Git 接管点

- `1981bdb7bc93fe35765af650b8d338b26f699b3f` 是导入本地 Git 后的源码接管基线。
- 本 ZIP 还会包含这份接管说明对应的后续文档提交，因此实际 HEAD 必须以
  `git rev-parse HEAD` 为准。
- 当前正式分支是 `main`。
- 当前没有正式 release tag，也没有配置远程仓库；ZIP 内含完整本地 `.git` 历史。
- 不得 rewrite history、force reset 或删除已有正式提交。
- 一个子阶段只有在全部验收通过后，才能形成一个只包含该子阶段修改的 checkpoint commit。

## 4. 当前已经核实的实现状态

| 链路 | 当前事实 |
|---|---|
| 均匀 Cartesian/Cut-cell → OpenFOAM | 可以生成完整 `polyMesh`；已有独立读取与历史验证路径 |
| 自适应 octree Cut-cell | 内部几何和拓扑已经存在 |
| 自适应 octree → 完整 OpenFOAM | Stage 6.1 已打通 reference ASCII writer；固定 cube/L-prism 的独立 reader 与 OpenFOAM 2606 `Mesh OK.` 已通过 |
| 增量自适应 | 已有映射/复用能力，但没有完整 solver output 链路 |
| 旧 Stage 6 compact 路径 | 仅支持均匀路径和 binary `polyMesh`，不得当作 revised Stage 6 已完成 |
| Stage 7 hybrid 功能 | 当前未实现，也不允许在本阶段实现 |

此前 Stage 6 尝试的历史状态是 `blocked_openfoam_quality`。历史 `checkMesh` 证据包含
3293 个错误 face pyramids、12078 个高 skew faces、最大 skewness 1274.65，且没有
`Mesh OK`。因此“完成导出”不能描述为“solver-ready”。

不要恢复或重新引入此前失败的 writer 修补分支，包括：

- `--convex-piece-cells`；
- `convex_piece_exact`；
- `quality_conformal_split`；
- `quality_constrained_cell_merge`；
- writer 层的 2–N merge；
- kernel tetra repair。

Stage 6.0 机器基线见 `artifacts/stage6_baseline.json`，人工记录见
`docs/STAGE6_BASELINE.md`。Stage 6.1 的当前实现、确定性、独立读取和 OpenFOAM 证据见
`docs/STAGE6_1_VERIFICATION.md` 与本机 `artifacts/stage61/`。

Stage 6.1 只证明固定解析 adaptive cube/L-prism 的完整拓扑链路；不得外推为复杂 STL、
大规模 adaptive 或全部质量门已经通过。Stage 6 总体仍未完成。

## 5. Windows 上的推荐工作方式

推荐在 **WSL2 的 Linux 文件系统**内解压和开发，例如 `~/src/cartmesh`。不要长期把构建目录
放在 `/mnt/c` 或 `/mnt/d` 下，以免文件系统性能影响测试结果。Git、CMake、GCC、Python、
测试依赖和 OpenFOAM/Docker 等工具应在 Windows/WSL 上单独安装；ZIP 不包含第三方安装包。

当前网格器代码没有 CUDA/GPU 计算路径。RTX 4060 不会自动加速现有 C++ 网格算法；Windows
机器的价值主要是更长时间运行、更多 CPU/RAM 资源，以及更方便保留 OpenFOAM 验证环境。
未经单独设计和批准，不要为了使用显卡改写算法。

项目数据不得上传到外部服务。这个 ZIP 是本地机器间交接，不等于 GitHub 或云同步。

## 6. 新 Codex 会话的第一轮只读检查

先完整阅读：

1. `CARTESIAN_MESH_GENERATOR_PROJECT_BRIEF_CN.md`
2. `AGENTS.md`
3. `README.md`
4. `CHANGELOG.md`
5. `STAGE6_REVISED_PLAN_CN.md`
6. `STAGE7_HYBRID_CARTESIAN_CODEX_PLAN_CN.md`
7. `docs/` 中 Stage 0–6 的设计与验收文档

然后执行且只执行只读/验证命令：

```sh
git status --short --branch
git log --oneline --decorate --graph --all
git branch --show-current
git tag --list
git remote -v
git rev-parse HEAD

cmake --preset release
cmake --build --preset release
ctest --preset release --output-on-failure
```

如果 Windows/WSL 环境还没装齐依赖，应报告缺失项，不得把“无法运行测试”写成“测试通过”。
独立 OpenFOAM 验证也必须记录实际命令、stdout/stderr、质量指标和工具版本。

完成检查后，先核对 Stage 6.0/6.1 的验证记录和当前 HEAD，再向用户提交 Stage 6.2
开工方案并等待确认。不要在同一轮顺手实现 Stage 6.2。

## 7. 可直接交给 Windows Codex 的开场指令

```text
先完整阅读 AGENTS.md 和 WINDOWS_CODEX_HANDOFF_CN.md，再按接管说明列出的顺序阅读项目文件。
当前唯一计划是 STAGE6_REVISED_PLAN_CN.md。Stage 6.0 和 6.1 已完成，先只读核实
docs/STAGE6_BASELINE.md、docs/STAGE6_1_VERIFICATION.md、artifacts/stage61 和当前 HEAD。
代码、测试、验证产物和 Git 历史是最高事实来源。先核实 Git/branch/tag/工作区和测试基线，
再报告 Stage 6.2 原生质量评估器的接口、最小失败案例、拟修改文件和验收命令。
在我明确确认前，不开始 Stage 6.2 代码修改，也不实现任何 Stage 7 功能。
```

只要 ZIP 被完整解压且 `.git` 隐藏目录没有被丢失，Windows Codex 从仓库根目录打开后，读取
`AGENTS.md` 就会被引导到本文件，并能知道当前应做什么和不应做什么。
