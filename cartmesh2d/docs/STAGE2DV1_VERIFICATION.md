# Stage 2D-V1：真实 OpenFOAM 验收与最小流动求解

日期：2026-08-24

> 后续状态：circle、airfoil-like 与 two-obstacles 已在 V1b 通过真实 OpenFOAM
> `-allGeometry`；本文件下方的“剩余边界”是 V1a 当时的失败基线。当前证据见
> `STAGE2DV1B_VERIFICATION.md`。NACA2412 dense 仍未验收。

## 结论

V1 已完成一个可复现的矩形障碍物外流基准：生成器直接写出的 OpenFOAM 2606
算例通过 default、`-allTopology`、`-allGeometry` 三套 `checkMesh`，并由
`simpleFoam` 在 255 次迭代收敛。环域和含孤岛的两个多环产品也分别通过三套
真实 `checkMesh`；后者被 OpenFOAM 正确识别为两个 fluid regions。

这不是对所有输入的泛化声明。circle、airfoil-like、NACA2412 和
two-obstacles 仍被内部 solver-quality gate 阻止，详见“剩余边界”。

## 外部检查发现并修复的问题

1. 初版 `polyMesh` 的内部面虽然满足 `owner < neighbour`，但没有按
   `(owner, neighbour)` 排序。OpenFOAM 报告 129 个 unordered faces；现在写出器
   使用确定性 upper-triangular 顺序，独立读取器也检查这一条件。
2. 原始自适应粗细过渡单元含共面分裂侧面。default checkMesh 可通过，
   `-allGeometry` 则正确报告 64 个 concave cells。
3. `SolverTopology2D` 只对含 reflex/collinear vertex 的输入 cell 做确定性 ear
   clipping；物理边界角产生的单内面三棱柱再与相邻三角形合并为严格凸四边形，
   避免 under-determined cell。没有降低 OpenFOAM 阈值，也没有删除单元。

## 矩形外流最终产品

- 原稳定拓扑：248 cells；局部分割输入：64 cells；最终 solver topology：
  372 cells（188 hexahedra + 184 triangular prisms）。
- 680 points、1456 faces、592 internal faces。
- 独立产品读取：逐单元闭合残差 0，最小体积
  `6.103515624999996e-05`，所有索引/patch/upper-triangular 检查 PASS。
- 内部门：最大 non-orthogonality 45 度、最大 internal skewness 0.5、
  最大水力 aspect 13.4018、最短二维 face 0.03125。
- OpenFOAM：最大 non-orthogonality 45 度、最大 skewness 1、最小 cell
  determinant 0.02748001658、最小 interpolation weight 0.1818181818、最小
  volume ratio 0.1666666667；三套检查均为 `Mesh OK.`。
- 13 个生成输入文件在两次独立运行中逐字节一致。

生成器同时写出低雷诺数层流验证骨架：入口速度 0.1、运动黏度 0.01、wall no-slip、
top/bottom slip、right pressure reference。OpenFOAM 2606 `simpleFoam` 在 255 次
迭代满足 `p/U=1e-7` residual control；最后一步 continuity sum local
`1.850595891e-10`、global `2.186001066e-11`。`foamToVTK` 导出的 372-cell 场由
ParaView 6.2 VTK reader 独立读取，`U` 和 `p` 全部为有限值。

## 多环与多连通域

- annulus：1840 cells，wall_0/wall_1，三套 `Mesh OK`；最大 skewness 2.8125，
  最小 determinant 0.001512738451。
- nested-island：2032 cells，wall_0/wall_1/wall_2，OpenFOAM 报告 2 regions，
  三套 `Mesh OK`；最大 skewness 0.8333333333，最小 determinant 0.01581976423。

这证明网格产品不再只支持单环、单连通拓扑。封闭环域没有在 V1 强行套用外流边界
条件做求解；本阶段只对其做完整 mesh acceptance。

## 回归与失败关闭

- Release CTest：34/34 PASS。
- ASan+UBSan CTest：34/34 PASS。
- `check_openfoam_v1_logs.py` 不相信 OpenFOAM 进程退出码本身：每个检查日志必须
  包含 `Mesh OK.` 且不能含 `Failed`/`FOAM FATAL`；求解日志必须有 clean
  convergence，最后 local/global continuity 必须小于 `1e-8`。

## 剩余边界

- circle：最小角约 0.351 度，且存在 non-orthogonality 83.60 度和 skewness
  13.68 的局部三角形。
- airfoil-like：仍有 6 个内部面 non-orthogonality 超过 70 度，最大 81.25 度，
  并有一个 0.348 度小角。
- NACA2412 dense：53 个质量问题，最坏水力 aspect 5892.65、小角 0.0194 度。
- two-obstacles：4 个内部面 non-orthogonality，最大 84.67 度。

下一开发不能调宽门限；应在 boundary-cut triangle 上做质量感知的 diagonal choice、
局部 edge flip/quad recombination 和必要的 boundary point redistribution，并保留这些
cell/edge ID 作为最小失败回归。
