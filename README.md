# cartmesh2d — 原生二维自适应 Cartesian / Cut-cell 网格生成器

`cartmesh2d` 是一个**独立、封闭、原生二维**仓库。

它不是三维 `cartmesh` 的 `z=0` 模式，也不通过把 `Point3D`、Octree 或 3D Cut-cell 代码模板化来复用三维核心。二维项目首先追求一个小而完整、可验证、可作为结题成果的二维 CFD 网格生成核心；三维项目继续独立推进。

## macOS 桌面应用

`desktop/` 提供可双击运行的 CartMesh2D 图形界面：选择 DXF、设置网格规模、生成并
预览真实 CM2D 单元，以及导出 VTK/CM2D/JSON/OpenFOAM case。Apple Silicon 本地
打包命令与产品边界见 `docs/PRODUCT1_MAC_DESKTOP_CN.md`。同一套 Electron 界面代码
保留 Windows 打包能力，但 Windows 可执行文件必须在 Windows 环境另行构建和验收。

## 默认 CFD 物理语义

二维产品默认与三维 `cartmesh` 保持一致：

```text
Domain2D = 外部计算域
BoundaryLoop = 固体壁面/障碍物轮廓
默认 fluid region = Domain2D - solid interior
```

因此，对一个放在矩形计算域中的翼型、圆柱、叶片截面或其他闭合物体：

- 物体内部**不生成流体网格**；
- 物体外部到计算域边界之间生成 Cartesian / Quadtree / Cut-cell 流体网格；
- 物体轮廓形成 `EmbeddedBoundary` 固体壁面；
- 计算域外框形成 `DomainBoundary`；
- 边界附近局部细化并由 Cut-cell 表达真实几何。

只有明确的内部流/管道流场景才使用 `FluidRegion2D::Interior`。内部流不是默认产品语义。

## 核心流水线

```text
2D solid closed boundary + outer computational domain
    -> geometry validation
    -> Cartesian background grid
    -> boundary/cell intersection classification
    -> adaptive Quadtree refinement
    -> 2:1 balance
    -> geometric inside / outside / intersected classification
    -> physical fluid-side selection (default: exterior)
    -> exterior Cut-cell polygon construction
    -> small-cell handling
    -> cell-edge-neighbor topology
    -> mesh quality validation
    -> solver/standard export
    -> visualization (last)
```

## 阶段

- **2D-0**：二维几何内核
- **2D-1**：均匀 Cartesian 网格与几何分类
- **2D-2**：Quadtree 自适应与 2:1 平衡
- **2D-3**：真实二维 Cut-cell polygon + 明确 fluid-side
- **2D-4**：完整 cell-edge-neighbor 拓扑，含 solid wall 与 outer domain boundary
- **2D-5**：小 Cut-cell 检测与稳定化/聚合
- **2D-6**：质量、导出、最终 CFD 语义验收
- **2D-V**：可视化；只有 2D-6 核心验收后才进入

基础阶段之后已经完成以下精细化路线：

- **H1**：可控高密度 sizing field；
- **H2**：大规模 Quadtree / topology scalability；
- **H3**：solver topology 与质量稳定化；
- **H4-1**：原生二维 boundary-layer quad strip；
- **H4-2**：outer envelope 与 Cartesian/Cut-cell remainder 共形拼接；
- **H4-3**：复杂几何局部降层、真实 layer termination 和最终 pure Cut-cell fallback。

当前 H4-3 基线同时验证 circle、superellipse、concave L、sharp trailing edge
和 narrow gap，并通过 OpenFOAM v2606 `checkMesh`。详细证据见
H4 算法说明见 `docs/STAGE2DH4_3_LOCAL_TERMINATION_CN.md`；质量事实基线与
provenance 见 `docs/Q0_QUALITY_BASELINE_PROVENANCE_CN.md` 和 `artifacts/q0/`；
无量纲 typed solver-quality contract 见
`docs/Q1_DIMENSIONLESS_TYPED_QUALITY_CONTRACT_CN.md` 和 `artifacts/q1/`。

Q2 当前仍为部分完成：superellipse 微短面修复见
`docs/Q2_INTERSECTION_CANONICALIZATION_PARTIAL_CN.md`；共享交点构造、公共边分割
及新旧路径对照见 `docs/Q2A_SHARED_CONSTRUCTION_DESIGN_CN.md` 和 `artifacts/q2a/`。
hybrid CLI 默认使用共享路径，末尾 `--legacy-construction` 可运行旧路径。
narrow gap / sharp trailing edge 的 Q1 短面失败尚未解决，不能以 `checkMesh` 通过
登记 Q2 全部完成。

## 独立构建与测试

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCARTMESH2D_BUILD_TESTS=ON
cmake --build build -j4
ctest --test-dir build --output-on-failure
```

本仓库不需要三维 `cartmesh` 源码，也不再使用
`-DCARTMESH_BUILD_2D=ON` 或 `add_subdirectory(cartmesh2d)`。

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

## 开发前必读

1. `AGENTS.md`
2. `docs/PROJECT_BRIEF_CN.md`
3. `docs/ARCHITECTURE_CN.md`
4. `docs/STAGE_PLAN_CN.md`
5. `docs/ACCEPTANCE_CN.md`

如果你要接手当前这一轮（R2：加密鲁棒性与 CFD 可解性），先读
`docs/R2_HANDOFF_CN.md`：它给出两条并行主线的关键问题、已被实测排除的路线，以及
每条改动的验收清单。实测数据与根因在 `docs/R2_REFINEMENT_ROBUSTNESS_CN.md`，
参考项目的许可证与借鉴等级在 `docs/R2_REFERENCE_AUDIT_CN.md`。

## 加密鲁棒性当前边界

```text
纯 Cut-cell（cartmesh2d_cli）      稳到 level 9（13024 solver cells、0.54 s）
                                   level 10 失败：格点 spur -> boundary skewness 8.6
hybrid 边界层（hybrid_cli）        稳到 level 8
                                   level 9 失败：壁面切向分辨率被输入折线顶点数锁死
```

阶梯测量与回归门：`tools/verification/refinement_ladder.py`，基线
`artifacts/r2/w0-baseline-manifest.json`。任何与加密相关的改动都必须给出改动前后
两份 manifest。

## PDE 验证用全局最低层级

边界自适应网格适合几何分辨，但单独增加 `max-level` 不会细化远场，不能构成受控的
PDE 网格收敛序列。CLI 的最后一个可选参数 `minimum-level` 可对整个计算域设置
Quadtree 层级下限；默认 `0` 完全保留原来的边界自适应行为：

```text
cartmesh2d_cli boundary.xy output 8 0.25 0.20 exterior openfoam-case 6
```

上例是全域至少 level 6、嵌入边界达到 level 8。该参数用于验证或需要全域分辨率的
仿真，不应被误写成边界自适应本身已经实现全局网格收敛。

## 高密度尺寸场

在全域 `minimum-level` 和物面 `max-level` 之间，可以叠加任意数量的物面距离带与
轴对齐矩形区域。所有尺寸场取最大的目标层级，因此参数顺序不会改变网格；目标层级仍受
`max-level` 上限约束。下游矩形区域就是第一版确定性尾迹加密原语：

```sh
./build/cartmesh2d_cli \
  examples/acceptance/airfoil_like.xy outputs/airfoil \
  9 0.30 0.02 exterior - 4 0 \
  --distance-band 0.05 7 \
  --refine-box 0.80 -0.10 1.30 0.10 8
```

`--distance-band` 的距离以及 `--refine-box` 的坐标使用规范化边界文件的长度单位；DXF
产品链中该单位是米。两个选项均可重复。程序额外写出 `.sizing.json`，记录计算域、完整
尺寸场、层级直方图、每个矩形的叶单元数以及 2:1 平衡结果。退化矩形、层级越界和完全
位于计算域外的矩形会明确失败。

## DXF-2：曲线、边界语义和单位换算

二维 CAD 输入先经过独立、fail-closed 的转换器，再进入原有网格核心：

```text
ASCII DXF -> 严格实体/闭环诊断 -> normalized boundary.xy
          -> cartmesh2d_cli -> Cartesian / Quadtree / Cut-cell / OpenFOAM
```

构建后运行：

```sh
./build/cartmesh2d_dxf_cli \
  examples/dxf/airfoil_like.dxf \
  artifacts/airfoil_like.xy \
  0.001 \
  artifacts/airfoil_like.dxf.json

./build/cartmesh2d_cli \
  artifacts/airfoil_like.xy artifacts/airfoil_like \
  7 0.30 0.02 exterior artifacts/airfoil_like-case 0 0
```

第三个 DXF 参数是输出米制坐标下的绝对弦高误差。支持 `LINE`、`ARC`、`CIRCLE`、
`LWPOLYLINE`（含 bulge）、`ELLIPSE` 和带显式 knot/control point 的 `SPLINE`。
`$INSUNITS=1..24` 自动换算成米；缺失或 unitless 文件必须在最后一个参数显式给出
源单位，例如 `mm`。DXF 图层按 `wall_*`、`inlet_*`、`outlet_*`、
`slip_*/farfield_*`、`symmetry_*` 写入 `.xy` 元数据并传到 OpenFOAM patch 和场边界。
详细边界和证据见 `docs/STAGE2D_DXF2_VERIFICATION_CN.md`。
