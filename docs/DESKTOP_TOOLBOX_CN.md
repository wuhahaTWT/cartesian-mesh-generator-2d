# CartMesh2D 分层工具箱前端

## 1. 产品分层

桌面端不再把所有算法参数平铺为一个 Cut-cell 表单，而是按以下层次组织：

1. **输入层**：DXF、源单位与曲线弦高误差；
2. **方法层**：稳定 `Pure Cut-cell` 与 `Hybrid + Boundary Layer (Beta)`；
3. **基础精细度层**：物面最高层级、全域最低层级、预设和单元预算；
4. **局部加密层**：Pure Cut-cell 的距离带与矩形/尾迹区域；
5. **输出验证层**：Job JSON、真实 CM2D/VTK、质量报告及有条件的 OpenFOAM case。

界面方法与参数能力来自 `desktop/src/mesh-tools.js` 的版本化 capability registry。
主进程只接受经过该工具层验证的 job，并把最终 job 与最低规模预估保存为
`<prefix>.job.json`，用于复现。

## 2. 方法边界

### Pure Cut-cell（稳定）

- 默认外流：闭合 DXF 轮廓是固体壁面；
- 支持物面层级、远场最低层级、距离带、矩形/尾迹区域；
- 通过内部 solver-quality 后输出 OpenFOAM case；
- `minimum-level` 不随“物面较密”预设上升，防止无意全域加密。

### Hybrid + Boundary Layer（Beta）

- 支持壁面层数、首层厚度、增长率与连接层级；
- 高层级和狭窄间隙仍受 R2/W2 边界约束；
- 前端读取并显示真实 `mesh_mode`，Pure Cut-cell fallback 不冒充 Hybrid；
- 当前只生成检查用 CM2D/VTK/质量报告，不请求 OpenFOAM 输出；
- solver-quality 通过但严格 quality contract 失败时，界面显示红色
  `CONTRACT FAIL`。

Boundary-layer core 和 Q3/Q4/Q5 构造 flag 是内部诊断机制，不作为普通用户的独立
“网格方法”。

## 3. 规模预算语义

界面显示 `4^minimumLevel` 作为**全域底格下限**。这不是总单元数预测：物面、
距离带、尾迹区和壁面层还会增加单元。若仅这个下限已经超过用户预算，任务在调用
生成器前失败。二维全域最低层级每增加一级，底格数量约乘 4。

## 4. 2026-09-04 桌面实测

打包产物：`desktop/dist/mac-arm64/CartMesh2D.app`。应用内含：

- `cartmesh2d_dxf_cli`
- `cartmesh2d_cli`
- `cartmesh2d_hybrid_cli`

内置 SPLINE 圆的真实应用冒烟结果：

| 方法 | 实际结果 | 2D/solver cells | OpenFOAM | 前端质量状态 |
|---|---|---:|---|---|
| Pure Cut-cell | Pure Cut-cell | 420 | 620 cells | PASS |
| Hybrid Beta | Hybrid | 908 | 未请求 | CONTRACT FAIL（256 个 hard issues） |

Pure Cut-cell 输出由独立读取器确认：620 cells、2528 faces、1336 points、最大单元
闭合残差 `1.0842e-19`、无 issue；OpenFOAM v2606
`checkMesh -allGeometry` 返回 `Mesh OK`，最大 non-orthogonality `68.4204°`、最大
skewness `2.08062`、最小 face weight `0.0553217`、最小 volume ratio `0.0553241`。

界面证据：

- `artifacts/desktop/layered-tools-cutcell.png`
- `artifacts/desktop/layered-tools-hybrid-beta.png`

## 5. 验证命令

```sh
cd desktop
npm test
sh scripts/build-macos.sh
open -n dist/mac-arm64/CartMesh2D.app --args --smoke-test --smoke-method=cutcell
open -n dist/mac-arm64/CartMesh2D.app --args --smoke-test --smoke-method=hybrid
```
