# DXF-2：曲线、边界语义与单位换算

## 1. 唯一阶段范围

DXF-2 只扩展二维 CAD 输入链，不改变 Cartesian、Quadtree、Cut-cell、小单元或
solver-quality 阈值：

```text
ASCII DXF
  -> $INSUNITS / explicit override -> metre coordinates
  -> LINE / ARC / CIRCLE / LWPOLYLINE / ELLIPSE / SPLINE
  -> closed-loop diagnostics + layer role
  -> deterministic annotated .xy
  -> existing mesh core
  -> named OpenFOAM patch + matching 0/U and 0/p condition
```

依据 Autodesk 官方 group code：

- SPLINE：https://help.autodesk.com/cloudhelp/2018/ENU/AutoCAD-DXF/files/GUID-E1F884F8-AA90-4864-A215-3182D47A9C74.htm
- ELLIPSE：https://help.autodesk.com/cloudhelp/2023/ENU/AutoCAD-DXF/files/GUID-107CB04F-AD4D-4D2F-8EC9-AC90888063AB.htm
- `$INSUNITS`：https://help.autodesk.com/cloudhelp/2024/ENU/AutoCAD-DXF/files/GUID-A85E8E67-27CD-4C59-BE61-4DC9FADBE74A.htm

## 2. 曲线契约

- `ELLIPSE` 读取中心、major-axis vector、minor/major ratio 和起止参数；
- `SPLINE` 读取 degree、非降 knot vector、control points 和可选正权重，使用齐次
  de Boor 计算 NURBS；
- 每个非零 knot span 递归检查 1/4、1/2、3/4 参数点到 chord 的偏差；
- 所有坐标先换算到米，`maximum-chord-error` 和 endpoint tolerance 也使用米；
- 单实体最多输出一百万段，递归深度超过 40 明确失败；
- 仅含 fit points、没有显式 knot/control point 的 SPLINE 当前明确失败，不用折线
  冒充；degree 仅接受 1..16；
- 非平面控制点、非正权重、非法 knot/count、退化参数域均明确失败。

## 3. 单位契约

- `$INSUNITS=1..24` 按官方含义确定性换算到 SI 米；
- JSON 报告保存原始 code、有效 code/name、`coordinate_scale_to_metres`、是否使用覆盖值；
- code 0、缺失或未知单位默认 fail-closed；
- 用户确认源单位后可在 CLI 最后一个参数传 `in|ft|mm|cm|m|km|1..24`；
- `.xy` 和最终 OpenFOAM points 一律是米。

## 4. 边界语义契约

每个闭环必须来自一个 DXF layer；同一闭环混用多个 layer 会报告
`boundary_metadata_conflict`。层名转成合法 OpenFOAM word 后按前缀映射：

| DXF layer | polyMesh type | U | p |
|---|---|---|---|
| `wall_*` 或未识别名称 | `wall` | `noSlip` | `zeroGradient` |
| `inlet_*` | `patch` | `fixedValue` | `zeroGradient` |
| `outlet_*` | `patch` | `zeroGradient` | `fixedValue` |
| `slip_*` / `farfield_*` | `patch` | `slip` | `zeroGradient` |
| `symmetry_*` | `symmetryPlane` | `symmetryPlane` | `symmetryPlane` |

注释格式随 `.xy` 保留，旧的无注释 `.xy` 继续默认生成 `wall_N`。多个闭环可复用同名、
同 type/role 的 patch；同名冲突或占用 `left/right/top/bottom/frontAndBack` 会失败。

本阶段是“每闭环一个边界角色”，不把同一闭环不同线段静默拆成多个条件。

## 5. 命令

```sh
build/cartmesh2d_dxf_cli \
  cartmesh2d/examples/dxf/spline_circle_mm.dxf \
  artifacts/spline_circle.xy 0.001 artifacts/spline_circle.dxf.json

build/cartmesh2d_cli \
  artifacts/spline_circle.xy artifacts/spline_circle \
  6 0.25 0.05 exterior artifacts/spline_circle-case 0 0

python3 cartmesh2d/tools/verification/check_openfoam2d.py \
  artifacts/spline_circle-case --require-patch wall_spline:wall
```

## 6. 实测验收

- Release CTest：`51/51` PASS；
- Clang AddressSanitizer + UndefinedBehaviorSanitizer CTest：`51/51` PASS；macOS
  AddressSanitizer 不支持 LeakSanitizer，因此不声称完成 leak 检查；
- `ezdxf 1.4.4` 读取 `spline_airfoil_mm.dxf` 与 `spline_circle_mm.dxf` 并执行
  `audit()`：均为 `0 errors / 0 fixes`；
- 对 `spline_circle_mm.dxf`，ezdxf 独立曲线以 0.001 mm flattening 生成 2052 个
  reference points；它们到 cartmesh2d 输出 chord 的最大距离为
  `0.0006478275331833098 m`，小于请求的 `0.001 m`；
- 完整产品为 `620 cells / 2528 faces / 1336 points`，独立 reader 要求并确认
  `wall_spline:wall`；
- 独立 reader 另确认 `symmetry_curve:symmetryPlane` 同时出现在 polyMesh、`0/U`
  和 `0/p`；
- OpenFOAM 2606 default、`-allTopology`、`-allGeometry` 均为 `Mesh OK`；max
  non-orthogonality `68.42038909`、max skewness `2.080618689`、min face weight
  `0.05532174376`、min volume ratio `0.05532409036`；
- 两个独立目录重跑的 20 个 DXF/XY/VTK/CM2D/OpenFOAM/JSON 文件逐文件
  SHA-256 一致。

## 7. 明确未完成

- binary DXF、DWG；
- legacy `POLYLINE/VERTEX/SEQEND`；
- fit-point-only SPLINE；
- `BLOCK/INSERT` 展开、HATCH boundary、STEP/IGES；
- 同一闭环的逐线段多 patch；
- 自动判断一份 CAD 的仿真物理含义。

这些限制保留为明确错误或后续独立阶段，不通过忽略实体或猜测边界条件绕过。
