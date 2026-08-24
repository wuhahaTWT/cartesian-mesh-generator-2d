# Product-1：macOS 桌面版

本阶段只把已经完成的 DXF → 网格 → OpenFOAM 输出链包装成可双击应用，不修改网格
算法和质量阈值。界面代码采用 Electron，后续可在 Windows 使用同一套界面打包。

## 功能

- 选择 ASCII DXF 或使用内置 SPLINE 样例；
- 自动读取 `$INSUNITS`，也可为无单位文件显式选择源单位；
- 快速、标准、较密三个预设，以及手动层级/留白/小单元参数；
- 调用原生 `cartmesh2d_dxf_cli` 和 `cartmesh2d_cli`；
- 直接解析真实 `.cm2d` 输出，在 Canvas 中绘制单元与边界；
- 显示 2D cells、OpenFOAM cells、vertices 和 solver-quality；
- 输出 XY、VTK、CM2D、JSON 和 OpenFOAM case，并可在 Finder 打开目录；
- 质量失败时显示原始失败原因，不隐藏或继续写假结果。

## 本地开发与打包

```sh
cd cartmesh2d/desktop
npm install
sh scripts/build-macos.sh
```

产物位于 `cartmesh2d/desktop/dist/mac-arm64/CartMesh2D.app`。当前为本地 ad-hoc/未签名
测试应用；对外公开分发前需要 Apple Developer ID 签名和 notarization。

同名 DXF 重复生成时会自动追加时间戳，不覆盖已有 `.cm2d` 和配套结果。

自动产品冒烟测试可使用：

```sh
open -n dist/mac-arm64/CartMesh2D.app --args --smoke-test
```

## Product-1 验收边界

只验证：应用启动、内置 DXF 成功生成、真实 CM2D 预览、失败信息可见。不重新运行
完整算法回归；算法证据继续引用 DXF-2 阶段报告。
