# CartMesh2D

原生二维笛卡尔网格生成器，提供 macOS 桌面应用和命令行工具，可用于外流绕流、管道内流的网格准备与 OpenFOAM 导出。

## 功能

- **笛卡尔自适应网格**：通过四叉树控制疏密，在物面处生成真实 Cut-cell 多边形。
- **贴体边界层**：沿壁面生成四边形层，与笛卡尔余域共形连接；提供外流和内流模式，当前为 Beta。
- **自动选参**：选择“常规”或“更密”后生成网格，显示进度、采用的参数和检查结果；手动模式可设置壁面尺寸、曲率、间隙、尾迹及局部加密区。
- **几何输入**：支持 XY、CSV、TXT 坐标文件、SVG 和 DXF，附带翼型、圆柱、喷管等 11 个样例。
- **预览与导出**：缩放、拖动、线框和密度着色；结果先保存在临时工作目录，点击保存可导出一个 ZIP，包含网格、报告和 OpenFOAM case。

闭合轮廓默认表示固体，流体网格位于轮廓外部。喷管、管道等内流问题选择“内部为流体域”。几何、拓扑、Solver 质量和 Q1 合同分别显示检查结果，具体验证范围见[当前状态](docs/CURRENT_STATE_CN.md)。

## macOS 安装与启动

当前桌面构建目标为 Apple Silicon Mac。准备 Node.js/npm、CMake 和 Xcode Command Line Tools，在仓库根目录执行：

```sh
cd desktop
npm ci
sh scripts/build-macos.sh
```

构建完成后，双击仓库根目录的 **`打开CartMesh2D.command`**，或打开 **`desktop/dist/mac-arm64/CartMesh2D.app`**。

1. 选择内置样例或导入几何文件，确认内部表示固体还是流体。
2. 选择纯 Cut-cell 或贴体边界层，选择自动密度档位并生成。
3. 在右侧查看网格和检查结果，点击保存导出 ZIP。

完整操作说明见[桌面使用指南](docs/DESKTOP_APP_CN.md)。

## 命令行与开发

核心采用 C++20，桌面界面采用 Electron。macOS 构建和测试：

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/clang++
cmake --build build -j2
ctest --test-dir build --output-on-failure
build/cartmesh2d_cli --help
build/cartmesh2d_hybrid_cli --help
```

Linux 构建时指定本机 C++20 编译器。输出格式包括 CM2D、VTK、JSON 和 OpenFOAM case；OpenFOAM 导出将二维网格挤出为单层体单元，前后边界设为 `empty`。

界面开发在 `desktop/` 中执行 `npm start`。C++ 修改后重新执行打包脚本，更新应用中的生成器。

- [开发导航](docs/DEVELOPMENT_CN.md)：模块、构建、验证工具与目录说明。
- [当前状态](docs/CURRENT_STATE_CN.md)：已验证案例、性能数据与能力边界。
- [开发规则](AGENTS.md)：物理定义、代码维护和交付验证要求。
