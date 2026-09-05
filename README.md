# cartmesh2d — 原生二维自适应 Cartesian / Cut-cell / Hybrid CFD 网格生成器

`cartmesh2d` 是一个**独立、封闭、原生二维**的 CFD 网格生成项目。

它不是三维 `cartmesh` 的 `z=0` 模式，也不依赖把 `Point3D`、Octree 或 3D Cut-cell
代码模板化来复用三维核心。二维项目的目标是形成一个完整、可验证、可独立交付的二维
Cartesian / Cut-cell / Hybrid 网格生成器；三维项目继续独立推进。

## 当前项目状态

当前质量状态、五个验收案例、加密上限、各轮次关闭情况与正在推进的技术主线统一维护于：

**[docs/CURRENT_STATE_CN.md](docs/CURRENT_STATE_CN.md)**

README 只描述相对稳定的产品能力与使用入口，不重复登记会随开发变化的测试数字和阶段状态。

## 核心能力

- 原生二维 Cartesian 背景网格；
- 自适应 Quadtree 细化与 2:1 balance；
- 外流默认语义下的真实二维 Cut-cell polygon；
- 小 Cut-cell 检测、稳定化与 solver topology；
- 原生二维 boundary-layer quad strip 与 Cartesian/Cut-cell remainder 共形拼接；
- 复杂几何下的局部 layer termination 与 pure Cut-cell fallback；
- 无量纲、分类型的 solver-quality / Q1 质量合同；
- DXF 曲线离散、单位换算与边界 patch 语义传递；
- VTK、CM2D、JSON 与 OpenFOAM case 导出；
- 全域最低层级、物面距离带与矩形区域等确定性 sizing controls；
- macOS Electron 桌面应用，并保留 Windows 打包能力。

## 默认 CFD 物理语义

二维产品默认与三维 `cartmesh` 保持一致：

```text
Domain2D = 外部计算域
BoundaryLoop = 固体壁面 / 障碍物轮廓
默认 fluid region = Domain2D - solid interior
```

因此，对放在矩形计算域中的翼型、圆柱、叶片截面或其他闭合物体：

- 物体内部**不生成流体网格**；
- 物体外部到计算域边界之间生成 Cartesian / Quadtree / Cut-cell 流体网格；
- 物体轮廓形成 `EmbeddedBoundary` 固体壁面；
- 计算域外框形成 `DomainBoundary`；
- 边界附近局部细化并由 Cut-cell 表达真实几何。

只有明确的内部流 / 管道流场景才使用 `FluidRegion2D::Interior`。内部流不是默认产品语义。

## 核心流水线

```text
2D solid closed boundary + outer computational domain
    -> geometry validation
    -> Cartesian background grid
    -> boundary / cell intersection classification
    -> adaptive Quadtree refinement
    -> 2:1 balance
    -> geometric inside / outside / intersected classification
    -> physical fluid-side selection (default: exterior)
    -> Cut-cell / hybrid construction
    -> small-cell handling
    -> cell-edge-neighbour solver topology
    -> mesh / solver quality validation
    -> VTK / CM2D / JSON / OpenFOAM export
    -> visualization
```

## 快速构建与测试

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCARTMESH2D_BUILD_TESTS=ON
cmake --build build -j4
ctest --test-dir build --output-on-failure
```

本仓库不需要三维 `cartmesh` 源码，也不使用
`-DCARTMESH_BUILD_2D=ON` 或 `add_subdirectory(cartmesh2d)`。

主要命令行入口包括：

```text
cartmesh2d_cli          纯 Cartesian / Quadtree / Cut-cell 路径
cartmesh2d_hybrid_cli   boundary-layer + Cartesian / Cut-cell hybrid 路径
cartmesh2d_dxf_cli      DXF -> normalized 2D boundary 转换
```

## 质量与验证体系

项目不会用“图片看起来正常”或单一 `checkMesh` 结果代替内部质量合同。主要验证层彼此独立：

| 验证层 | 作用 |
|---|---|
| topology audit | 检查重复、孤立、非流形边以及面积守恒等拓扑不变量 |
| solver quality | 检查 face weight、volume ratio、non-orthogonality、boundary skewness 等硬限 |
| OpenFOAM `checkMesh` | 外部求解器生态的独立网格检查 |
| Q1 contract | 无量纲、分类型的产品质量合同，用于比一般 solver gate 更细地定位缺陷 |

当前逐案例结果和各门的准确关系只在
[docs/CURRENT_STATE_CN.md](docs/CURRENT_STATE_CN.md) 中维护。

## 关键验收不变量

默认外流案例至少必须满足：

```text
geometric Inside  -> solid -> Empty fluid cell
geometric Outside -> fluid -> Full fluid cell
fluid_area = domain_area - solid_area
EmbeddedBoundary > 0
DomainBoundary > 0
```

如果上述任意一条不满足，即使 topology audit 为 0、CI 全绿或图片能画出来，也不得宣称网格正确。

## 自适应与 sizing controls

边界自适应网格用于几何分辨；如果需要受控 PDE 收敛序列，可以用 `minimum-level`
为整个计算域设置 Quadtree 层级下限。默认值保留原有边界自适应行为。

在全域 `minimum-level` 和物面 `max-level` 之间，还可以叠加：

- `--distance-band`：按物面距离设置目标层级；
- `--refine-box`：对轴对齐矩形区域设置目标层级，可用于尾迹等区域加密。

多个 sizing 条件取最大的目标层级，并继续受 `max-level` 上限约束。程序会输出
`.sizing.json` 记录计算域、尺寸场、层级直方图和 2:1 balance 结果。

示例：

```sh
./build/cartmesh2d_cli \
  examples/acceptance/airfoil_like.xy outputs/airfoil \
  9 0.30 0.02 exterior - 4 0 \
  --distance-band 0.05 7 \
  --refine-box 0.80 -0.10 1.30 0.10 8
```

## DXF 输入

二维 CAD 输入先经过独立、fail-closed 的转换器，再进入网格核心：

```text
ASCII DXF
    -> 严格实体 / 闭环诊断
    -> normalized boundary.xy
    -> cartmesh2d_cli / cartmesh2d_hybrid_cli
    -> Cartesian / Cut-cell / Hybrid / OpenFOAM
```

示例：

```sh
./build/cartmesh2d_dxf_cli \
  examples/dxf/airfoil_like.dxf \
  artifacts/airfoil_like.xy \
  0.001 \
  artifacts/airfoil_like.dxf.json
```

支持 `LINE`、`ARC`、`CIRCLE`、`LWPOLYLINE`（含 bulge）、`ELLIPSE` 和带显式
knot / control point 的 `SPLINE`。DXF 单位会转换到米；边界图层语义可传递到 OpenFOAM patch。

详细设计与验收见
[docs/STAGE2D_DXF2_VERIFICATION_CN.md](docs/STAGE2D_DXF2_VERIFICATION_CN.md)。

## macOS 桌面应用

`desktop/` 提供可双击运行的 CartMesh2D 图形界面：读入 .xy / .dxf / .svg / .csv /
.txt 或 11 个内置样例，在纯 Cut-cell 与贴体边界层两条路径间选择，用物理量（远场倍
体长、壁面单元尺寸、每级带宽、曲率、间隙、尾迹）描述尺寸场，预览按层级着色的真实
CM2D 单元，并分开报告拓扑 / solver 质量 / Q1 三个门。

分层结构、格式边界、层级预算算式与已知缺陷见
[docs/DESKTOP_APP_CN.md](docs/DESKTOP_APP_CN.md)。

Apple Silicon 本地打包命令与产品边界见
[docs/PRODUCT1_MAC_DESKTOP_CN.md](docs/PRODUCT1_MAC_DESKTOP_CN.md)。同一套 Electron
界面保留 Windows 打包能力，但 Windows 可执行文件需要在 Windows 环境单独构建和验收。

## 开发路线

基础二维核心 `2D-0 ~ 2D-V` 建立几何、Cartesian、Quadtree、Cut-cell、拓扑、稳定化、
质量、导出和可视化能力；`H1 ~ H4` 继续建立高密度 sizing、scalability、solver topology
和 hybrid boundary-layer 路线。之后的 `Q / R` 系列用于质量合同、构造鲁棒性和加密能力提升。

README 不登记这些轮次的实时完成状态。完整阶段历史见：

- [docs/STAGE_PLAN_CN.md](docs/STAGE_PLAN_CN.md)
- [docs/CURRENT_STATE_CN.md](docs/CURRENT_STATE_CN.md)

## 文档导航

首次阅读：

1. [docs/PROJECT_BRIEF_CN.md](docs/PROJECT_BRIEF_CN.md) — 项目边界与目标；
2. [docs/ARCHITECTURE_CN.md](docs/ARCHITECTURE_CN.md) — 核心架构；
3. [docs/STAGE_PLAN_CN.md](docs/STAGE_PLAN_CN.md) — 阶段路线；
4. [docs/ACCEPTANCE_CN.md](docs/ACCEPTANCE_CN.md) — 验收原则；
5. [docs/CURRENT_STATE_CN.md](docs/CURRENT_STATE_CN.md) — **唯一当前状态事实源**。

质量与鲁棒性：

- [docs/Q1_DIMENSIONLESS_TYPED_QUALITY_CONTRACT_CN.md](docs/Q1_DIMENSIONLESS_TYPED_QUALITY_CONTRACT_CN.md)
- [docs/Q2A_SHARED_CONSTRUCTION_DESIGN_CN.md](docs/Q2A_SHARED_CONSTRUCTION_DESIGN_CN.md)
- [docs/Q3_Q4_TERMINATION_QUALITY_CN.md](docs/Q3_Q4_TERMINATION_QUALITY_CN.md)
- [docs/R2_REFINEMENT_ROBUSTNESS_CN.md](docs/R2_REFINEMENT_ROBUSTNESS_CN.md)

接手开发时先读 `CURRENT_STATE_CN.md`，再按其中指向进入当前轮次的 handoff 文档，避免从历史阶段文档推断现状。

## 开发者约束

提交代码前还应阅读 [AGENTS.md](AGENTS.md)。历史轮次文档保留当时的测量数据作为证据，
但不得用历史文档覆盖 `CURRENT_STATE_CN.md` 中登记的当前事实。
