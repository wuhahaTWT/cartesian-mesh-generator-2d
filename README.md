# CartMesh2D

原生二维笛卡尔网格生成器，提供 macOS 桌面应用和命令行工具，可用于外流绕流、管道内流的网格准备与 OpenFOAM 导出。

## 功能

- **笛卡尔自适应网格**：通过四叉树控制疏密，在物面处生成真实 Cut-cell 多边形。
- **无量纲尺寸**：纯 Cut-cell 与贴体路径共用参考长度 Lref，按 h/Lref 设置壁面和背景尺寸；生成后报告最终单元与壁面方向的实际尺寸。
- **贴体边界层**：沿壁面生成四边形层，与笛卡尔余域共形连接；提供外流和内流模式，当前为 Beta。
- **自动选参**：按约 5 千、1.5 万、5 万、15 万、30 万或自定义数量配置，最多 3 次调整；显示实际数量与目标偏差，可把实际参数转入手动微调。原常规/更密预设保留兼容入口。
- **双主题**：现代简约与科研终端，可切换并记住选择；网格显示模式独立。科研终端保留内置背景 artwork。
- **几何输入**：支持 XY、CSV、TXT 坐标文件、SVG 和 DXF，附带翼型、圆柱、喷管等 11 个样例。
- **曲线与贴体层**：SVG 三次曲线和 DXF 样条按弦误差细分；坐标文件保留原始折线。贴体层沿壁面分布，弯曲处依据局部空间调整厚度和层数。
- **预览与导出**：缩放、拖动、线框和密度着色；结果先保存在临时工作目录，点击保存可导出一个 ZIP，包含网格、报告和 OpenFOAM case。

闭合轮廓默认表示固体，流体网格位于轮廓外部。喷管、管道等内流问题选择“内部为流体域”。界面显示几何、拓扑和 Solver 质量结果；外部 OpenFOAM `checkMesh` 需单独运行，具体验证范围见[当前状态](docs/CURRENT_STATE_CN.md)。旧版 Q1 检测已从当前产品取消，历史证据仍保留在 Git 和旧报告中。

导入时确认源单位，程序统一换算为米并保留物理尺寸，不会将坐标文件自动缩放到 1 m。
参考长度可采用包围盒最长边或显式指定工程长度；例如 `h/Lref=0.005` 表示目标尺寸为参考长度的 0.5%。
实际网格会受四叉树取整、Cut-cell 和质量修复影响，应查看结果包中的 `.resolution.json`；该报告不代表 CFD 精度保证。

## 三步配置

1. 选模型，确认轮廓内部是固体还是流体；在“参考长度”中可指定圆柱直径或翼型弦长。
2. 选目标数量，按工况设置计算域留白，再生成。目标范围为 ±30%，并非精度等级；高密档可能需允许更高构造深度。
3. 查看实际数量和质量状态。需要调整时，点“用本次实际参数继续微调”：附近太粗就减小壁面尺寸，远处太粗就减小背景尺寸，细网格范围不够就增加带宽。

内置圆柱、翼型和喷管另有“载入已验证高密样例参数”按钮：会明确替换参考长度、域留白与高深度选项，转入手动检查后再生成。喷管混合网格的通用数量初值仍可能失败，可用该入口复现已验证案例。

自动数量模式先按流体面积估算背景尺寸，选择较粗的可实现四叉树层级，为近壁加密留出数量预算；随后按实际单元数调整带宽或尺寸。计算域、所选内外流语义、方法和质量门不随数量调整改变。加密分布、层数和增长率可展开设置；首层初值为壁面尺寸的 1/4，尚未由 y+ 或流动误差推导。

如果三次仍未进入目标范围，保留最接近目标的成功结果并明确提示。数量支持上限为 50 万，但这不是每个几何都已验证成功的承诺；本次实测范围见[当前状态](docs/CURRENT_STATE_CN.md)。

## 网格展示

下面是三个固定几何在纯 Cut-cell 与贴体边界层混合路径下的高档展示图；每张图对应的单元数、输入网格和验证边界见[展示说明](展示素材/展示说明.md)。

| 几何 | 纯 Cut-cell | 混合网格 |
| --- | --- | --- |
| Fixed 32-segment circle | [153,208 cells](展示素材/01_circle_pure_153208.png) | [164,970 cells](展示素材/04_circle_hybrid_164970.png) |
| NACA 2412 | [133,810 cells](展示素材/02_naca_pure_133810.png) | [140,305 cells](展示素材/05_naca_hybrid_140305.png) |
| Nozzle interior | [141,488 cells](展示素材/03_nozzle_pure_141488.png) | [148,375 cells](展示素材/06_nozzle_hybrid_148375.png) |

在其他 Markdown 页面中复用图片时，可使用公开仓库的 raw URL，例如：

```md
![Fixed 32-segment circle](https://raw.githubusercontent.com/wuhahaTWT/cartesian-mesh-generator-2d/main/%E5%B1%95%E7%A4%BA%E7%B4%A0%E6%9D%90/01_circle_pure_153208.png)
```

需要固定版本时，把 URL 中的 `main` 换成具体 commit SHA；完整文件名对应关系与展示限制见[展示说明](展示素材/展示说明.md)。

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
