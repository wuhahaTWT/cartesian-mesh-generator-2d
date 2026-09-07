# CartMesh2D — 原生二维笛卡尔网格生成器

独立的二维 Cartesian / Quadtree / Cut-cell 项目。目标：**网格够密、生成够快、桌面好用、可导入 OpenFOAM**。
默认生成障碍物外部的流体网格；贴体边界层功能仍为 Beta。

## 打开 App

本仓库根目录是唯一开发入口。直接双击根目录的 **`打开CartMesh2D.command`**。
它打开 `desktop/dist/mac-arm64/CartMesh2D.app`；尚未打包时会说明如何构建。
不要再使用 `.claude/worktrees/…` 的旧路径。

首次构建（需要 Node.js/npm、CMake 和 macOS Command Line Tools）：

```sh
cd desktop
npm ci
sh scripts/build-macos.sh
```

已有依赖时无需重复 `npm ci`。打包脚本从当前源码编译三个 CLI、准备样例、执行前端测试并生成 App。
开发时在 `desktop/` 执行 `npm start`；改过 C++ 后先重新运行打包脚本，避免前后端版本不一致。

## 接手只读这几份

| 文件 | 用途 |
|---|---|
| [AGENTS.md](AGENTS.md) | 开发边界、物理语义、目录维护规则 |
| [当前状态](docs/CURRENT_STATE_CN.md) | 已验证能力、已知限制、下一步；唯一状态入口 |
| [开发导航](docs/DEVELOPMENT_CN.md) | 各代码模块、构建测试、CFD 验证及历史查找 |
| [桌面使用](docs/DESKTOP_APP_CN.md) | 输入、密度设置、显示和导出 |

正常接手先读 AGENTS 和当前状态，再按任务进入代码；无需遍历历史记录。

## 命令行构建

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/clang++
cmake --build build -j2
ctest --test-dir build --output-on-failure
```

以上为本机 macOS 命令，Linux 使用本机 C++20 编译器。参数以 `build/cartmesh2d_cli --help` 为准。
输出支持 CM2D、VTK、JSON、OpenFOAM case；OpenFOAM 使用二维网格挤出的一层体单元，不是二维核心依赖三维代码。
Windows 安装包仍需在 Windows 单独构建验证。

## 文件放哪里

`src/`、`include/` 是唯一核心实现；`apps/` 是 CLI；`desktop/src/` 是唯一前端。
`tests/` 与 `examples/` 保留测试和输入；`tools/` 保留仍有用途的验证和绘图工具。
`build/`、`desktop/runtime/`、`desktop/dist/`、`outputs/` 都是可重建、忽略提交的目录。
`artifacts/reference/` 只保留少量历史基准；`artifacts/current/` 保存当前验证摘要和预览。
历史阶段说明和旧实验从 Git 找回，不再在工作目录放第二套源码或几十份交接报告。
