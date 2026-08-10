# 阶段 0 验证记录

后续阶段的独立验证记录分别见 `docs/STAGE1_VERIFICATION.md`、
`docs/STAGE2_VERIFICATION.md`、`docs/STAGE3_PROGRESS.md` 和
`docs/STAGE4_VERIFICATION.md`；本文件保留阶段 0 当时的实现与阶段边界，不回写
后续功能来改变历史验收含义。

状态：**阶段 0 已于 2026-08-07 通过**

本记录只覆盖阶段 0 的均匀 Cartesian 网格基础。球体案例在背景网格上使用单元中心
内外采样；它不是几何 Cut-cell 结果，也不会被描述成可供求解器使用的曲面边界拓扑。

## 项目总纲要求对应关系

| 第一轮要求 | 实现证据 |
|---|---|
| C++20/CMake 工程 | `CMakeLists.txt`、`CMakePresets.json`、全新目录构建与安装 |
| 测试、基准和格式化配置 | CTest 目标、`cartmesh_benchmark`、`.clang-format` |
| 坐标、AABB、Triangle、CellKey、MortonCode | `include/cartmesh` 下的公共头文件及单元测试 |
| 均匀三维 Cartesian 网格 | 使用隐式存储和确定性 ID 的 `UniformCartesianGrid` |
| 解析立方体和球体 | 精确体积/面积类型及球体网格加密测试 |
| VTK 输出 | 流式 ASCII VTU 写出器，包含点、连接、偏移、类型和单元数据 |
| 使用 ParaView 实际打开 | ParaView 6.1.1 已读取并渲染两个生成文件 |
| 单元数、体积、墙钟和峰值 RSS | CLI JSON 报告及实测基准记录 |
| `VERIFICATION.md` | 本验证记录 |
| 千万级基础结构审查 | `docs/STAGE0_DESIGN_REVIEW.md` 记录了适用性和限制 |

## 验收证据

| 检查项 | 证据 | 状态 |
|---|---|---|
| 全新目录配置与构建 | Apple Clang 21 的 `/tmp` 全新 Release 构建/安装；GCC 15.2 的 Debug/Release 构建 | 通过 |
| 单元与解析测试 | 7 组单元测试，以及 CLI 成功/失败和基准冒烟测试 | 通过 |
| 清理器检查 | Apple Clang AddressSanitizer 和 UndefinedBehaviorSanitizer | 通过 |
| 确定性结果哈希 | 两次球体运行均为 `a5370e0889e02494`，VTU 的 SHA-256 完全相同 | 通过 |
| VTK XML 格式正确性 | 对立方体和球体执行 `xmllint --noout` | 通过 |
| 独立拓扑结构检查 | `tools/verify_vtu.py` | 通过 |
| 外部程序库读取 | meshio 5.3.5 正确读取两个 VTU 的数量和字段 | 通过 |
| 外部可视化 | meshio 与 Matplotlib 渲染了 6,960 个选中边界面 | 通过 |
| ParaView 直接读取与渲染 | ParaView 6.1.1 / VTK 9.6.2 读取两个文件并生成离屏 PNG | 通过 |
| 墙钟时间和峰值 RSS | CLI JSON 报告及沙箱外 `/usr/bin/time -l` | 通过 |
| 百万单元分类基准 | 100³ Release 基准 | 通过 |

## 验证环境

- 硬件：MacBook Air（Apple M1、8 核、8 GB 内存）
- 操作系统：macOS 26.5.2，arm64
- CMake：4.4.2
- 主要 Debug/Release 编译器：GCC 15.2.0
- 辅助编译器及清理器运行环境：Apple Clang 21.0.0
- 线程数：1

两个编译器均完成配置、构建并通过 CTest，未产生编译器警告。另在 `/tmp` 中完成了
一次全新的目录外配置、构建、4 项 CTest 测试和安装，结果全部通过。

## 解析结果和输出结果

48³ 球体背景网格包含 110,592 个单元。单元中心采样选中 29,464 个单元，估算体积为
`4.16282371238426`，解析体积为 `4.18879020478639`，相对误差为 `0.00619904343`
（约 0.62%）。在 96³ 分辨率下，估算值为 `4.185288040726274`，表明这一测试序列
具有预期的网格加密趋势。

生成的立方体 VTU 为 1,630,076 字节，球体 VTU 为 15,466,164 字节。重复生成的球体
文件逐字节相同，SHA-256 为：

```text
1d02d0b15a174339f2f7dd074e690f216579789dd8d683ca2566c8fff5711dea
```

外部 meshio 读取器报告球体文件包含 117,649 个点、110,592 个六面体、两个预期的
单元数据字段，以及 29,464 个选中单元。对应可视化证据为
`artifacts/sphere_48_meshio.png`。

ParaView 6.1.1 使用 VTK 9.6.2 的 `vtkXMLUnstructuredGridReader` 独立读取文件：
立方体包含 15,625 个点和 13,824 个六面体；球体背景网格包含 117,649 个点和
110,592 个六面体。它识别到两个预期字段，并在 `>= 0.5` 阈值下精确选中 29,464 个
单元。读取器错误码为零，所有单元均报告为 VTK 类型 12。机器可读报告和渲染证据为：

- `artifacts/cube_24_paraview.json`
- `artifacts/cube_24_paraview.png`
- `artifacts/sphere_48_paraview.json`
- `artifacts/sphere_48_paraview.png`

## 性能证据

完整的 48³ 球体运行（包括写出 15 MB ASCII VTU）内部计时为 0.2248 秒。
独立 `/usr/bin/time -l` 测得墙钟时间 0.25 秒，最大 RSS 为 8,781,824 字节。

100³ 基准精确遍历并分类 1,000,000 个隐式单元。内部核心计时为 0.006736 秒；
包含程序启动和测量开销的外部进程墙钟时间为 0.20 秒；最大 RSS 为 8,781,824 字节。
这些数据只适用于均匀解析单元中心分类，不代表 STL、八叉树或 Cut-cell 性能。
完整原始数据位于 `benchmarks/baselines/stage0_m1_2026-08-07.json`。

## 复现命令

```sh
cmake --preset release
cmake --build --preset release --parallel 4
ctest --preset release
build/release/cartmesh_cli --case sphere --resolution 48 \
  --output artifacts/sphere_48.vtu --report artifacts/sphere_48.json
xmllint --noout artifacts/sphere_48.vtu
python3 tools/verify_vtu.py artifacts/sphere_48.vtu
build/release/cartmesh_benchmark 100
pvpython tools/paraview_verify.py artifacts/sphere_48.vtu \
  artifacts/sphere_48_paraview.png --field inside_sphere_center_sample \
  --report artifacts/sphere_48_paraview.json
```

可选的外部程序库和 ParaView 验证脚本分别为 `tools/meshio_verify.py` 和
`tools/paraview_verify.py`。

## 阶段边界

以下功能有意保持未实现：STL 输入、BVH、自适应细化、2:1 平衡、Cut-cell 重建、
求解器格式输出、GUI、CFD 求解器、AI 生成网格、云服务和部署逻辑。

数据布局的可扩展性决定及 ASCII 输出的已知限制记录在
`docs/STAGE0_DESIGN_REVIEW.md` 中。
