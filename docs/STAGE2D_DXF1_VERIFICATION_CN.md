# DXF-1：ASCII DXF 二维闭合轮廓导入

> 本文保留 DXF-1 当时的验收边界。当前产品能力已由
> `STAGE2D_DXF2_VERIFICATION_CN.md` 扩展，不应再用本文“非本阶段能力”判断当前状态。

## 1. 阶段边界

本阶段只新增 CAD 输入适配层，不修改 Cartesian、Quadtree、Cut-cell、小单元、
solver topology 或 OpenFOAM writer 算法。产品路径是：

```text
ASCII DXF
  -> group code / entity parser
  -> analytic line/arc/circle or LWPOLYLINE bulge sampling
  -> endpoint graph and closed-loop assembly
  -> BoundaryRegion2D diagnostics and nesting normalization
  -> deterministic .xy + JSON report
  -> existing cartmesh2d_cli
```

Autodesk DXF 参考规定 ASCII 文件使用 group-code/value 两行对；group code `0`
引出实体。DXF-1 按官方 `LINE`、`ARC`、`CIRCLE`、`LWPOLYLINE` group code 读取，
其中 LWPOLYLINE 的 `70` bit 1 表示闭合、`42` 表示 bulge：

- https://help.autodesk.com/cloudhelp/2025/ENU/AutoCAD-DXF/files/GUID-89CB823D-614D-4D1E-8204-568EC72DF869.htm
- https://help.autodesk.com/cloudhelp/2025/ENU/AutoCAD-DXF/files/GUID-FCEF5726-53AE-4C43-B4EA-C84EB8686A66.htm
- https://help.autodesk.com/cloudhelp/2024/ENU/AutoCAD-DXF/files/GUID-0B14D8F1-0EBA-44BF-9108-57D8CE614BC8.htm
- https://help.autodesk.com/cloudhelp/2019/ENU/AutoCAD-DXF/files/GUID-8663262B-222C-414D-B133-4A8506A27C18.htm
- https://help.autodesk.com/cloudhelp/2026/CSY/AutoCAD-DXF/files/GUID-748FC305-F3F2-4F74-825A-61F04D757A50.htm

## 2. 支持和拒绝边界

支持：

- ASCII DXF `ENTITIES` section；
- `LINE`；
- `ARC` 和 `CIRCLE`；
- closed/open `LWPOLYLINE`，包括 bulge 圆弧；
- 多实体拼接、多闭环、孔洞和多个固体；
- 图层名与 `$INSUNITS` code 的报告；
- 所有实体位于同一个 XY 平面且使用默认 `+Z` OCS。

明确拒绝：

- binary DXF；
- `POLYLINE/VERTEX/SEQEND`、`SPLINE`、`ELLIPSE`、`HATCH`、`INSERT/BLOCK`；
- 非默认 OCS、非零 thickness、带宽 LWPOLYLINE；
- 缺失/非法数值 group、非有限坐标、零长度实体；
- 开口端点、分叉、重复导致的 degree != 2；
- endpoint tolerance 的传递聚类直径超过配置值；
- 自交、相切/相交环、零面积和非法嵌套。

拒绝项不会被静默删除或跳过。失败报告保留 issue code、DXF 行号、实体类型、
图层和原因。

## 3. 曲线离散与确定性

CLI 的 `maximum-chord-error` 是图纸单位下的绝对弦高上限。圆弧段数由半径、扫角和
弦高解析计算，不使用固定“每圆若干点”。每条解析曲线最多一百万段；超过时失败，
避免错误单位导致内存失控。

拼接前仅对实体端点使用显式 `endpoint-weld-tolerance`；内部几何谓词继续使用项目的
`TolerancePolicy`。闭环最终进入 `BoundaryRegion2D::diagnose`，归一化 nesting/orientation，
从字典序最小顶点开始并按字典序排序环。相同输入和参数的 `.xy` 必须字节一致。

## 4. 命令

```sh
cmake -S cartmesh2d -B build -DCMAKE_BUILD_TYPE=Release \
  -DCARTMESH2D_BUILD_TESTS=ON
cmake --build build -j

build/cartmesh2d_dxf_cli \
  cartmesh2d/examples/dxf/airfoil_like.dxf \
  artifacts/dxf/airfoil_like.xy 0.001 \
  artifacts/dxf/airfoil_like.dxf.json

build/cartmesh2d_cli \
  artifacts/dxf/airfoil_like.xy artifacts/dxf/airfoil_like \
  7 0.30 0.02 exterior artifacts/dxf/airfoil_like-case 0 0

python3 cartmesh2d/tools/verification/check_openfoam2d.py \
  artifacts/dxf/airfoil_like-case \
  --report artifacts/dxf/airfoil_like-case/independent_check.json
```

## 5. 回归矩阵

`tests/dxf_test.cpp` 覆盖：

- closed LWPOLYLINE、乱序/反向 LINE loop；
- CIRCLE 弦高和面积；
- LWPOLYLINE bulge 半圆；
- 重复运行 `.xy` 字节一致；
- layer 和 `$INSUNITS` report；
- open boundary、自交、SPLINE、非默认 OCS、非平面 LINE 失败；
- binary DXF、奇数行 group pair 和缺失 `ENTITIES` section 失败。

CTest 另执行真实 `airfoil_like.dxf -> .xy -> OpenFOAM polyMesh`，再用不链接项目的
`check_openfoam2d.py` 读取 points/faces/owner/neighbour/boundary 并检查正体积、闭合和
常量通量守恒。

## 6. 本阶段实测验收

- Release 构建：CTest `45/45` PASS；
- Clang AddressSanitizer + UndefinedBehaviorSanitizer 构建：CTest `45/45` PASS；
  macOS 的 AddressSanitizer 不支持 LeakSanitizer，因此这里不声称完成 leak 检测；
- 第三方 `ezdxf 1.4.4` 分别读取两个 R2013 示例并运行 `audit()`：
  `airfoil_like.dxf` 和 `rounded_obstacle.dxf` 均为 `0 errors / 0 fixes`；
- `airfoil_like.dxf` 的完整产品链生成 `835 cells / 3397 faces / 1784 points` 的
  OpenFOAM mesh；独立 Python reader 报告 `valid=true`、最小体积
  `1.3453165690103791e-08`、无 issue，闭合和常量通量残差通过；
- OpenFOAM 2606 `checkMesh` default、`-allTopology`、`-allGeometry` 均为
  `Mesh OK`。`-allGeometry` 报告 max non-orthogonality `63.85623373`、max skewness
  `3.842633507`、min face weight `0.05817665015`、min volume ratio
  `0.02538617019`；
- 两个独立输出目录重跑后，`.xy`、DXF JSON、VTK、CM2D、quality/viz JSON、
  OpenFOAM `polyMesh`、case 文件和独立 reader JSON 的 SHA-256 逐文件一致。

## 7. 已知真实失败案例

`examples/dxf/rounded_obstacle.dxf` 使用两个 LINE 和两个 ARC。DXF 转换本身 PASS，
但当前 solver-quality 路线在 `max-level=7, padding=0.5, small-alpha=0.05` 下出现真实
boundary skewness 失败，因此不写 OpenFOAM mesh。该文件保留为后续 solver topology
修复的回归输入；不得通过放宽 skewness=4.0 门限宣称通过。

## 8. 非本阶段能力

- 图层到 patch 名/type 的自动映射；
- 自动单位缩放；
- SPLINE/ELLIPSE；
- BLOCK/INSERT 展开；
- binary DXF 和 DWG；
- STEP/IGES/OpenCascade；
- DXF 计算域外框与 wall/inlet/outlet 的多 patch 语义。

这些能力必须各自作为后续单独阶段，不在 DXF-1 中预告完成。
