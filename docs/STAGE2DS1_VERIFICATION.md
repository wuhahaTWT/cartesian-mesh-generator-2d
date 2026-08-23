# Stage 2D-S1：求解器网格导出与质量门

日期：2026-08-24

## 本阶段边界

S1 把已经通过全局拓扑审计的二维多边形网格挤出为单层 OpenFOAM
`constant/polyMesh`，并在写出前执行求解器导向的几何质量门。S1 只证明产品
格式、面定向、正体积、闭合性、确定性和内部质量筛选；真实 OpenFOAM
`checkMesh` 与流动求解属于 V1，不能用本文的独立读取器替代。

## 实现

- `SolverQuality2D` 检查内部面 non-orthogonality、内部面 skewness、单元凹角、
  水力尺度 aspect、最小内角和最短面，并在失败时阻止 OpenFOAM 写出。
- `OpenFoam2D` 将每个二维单元挤出为一个棱柱。内部面排在边界面之前；
  owner/neighbour 与面定向一致；底面反向、顶面保持二维 CCW。
- 输入的每条独立环保留为确定性 `wall_0`、`wall_1` 等 patch；计算域四边为
  `left/right/bottom/top`；前后面合并为 `frontAndBack`、类型 `empty`。
- CLI 新增最后一个可选参数 `openfoam-case-dir`，成功时同时写出
  `solver_quality.json`。
- `tools/verification/check_openfoam2d.py` 是不链接本项目的独立产品读取器，
  检查语法计数、索引、patch 连续覆盖、owner/neighbour、逐单元定向闭合、
  正体积及未引用点。

## 可复现产品

命令：

```sh
build/cartmesh2d_cli examples/acceptance/rectangle.xy out/rectangle \
  5 0.25 0.05 exterior out/rectangle-case
python3 tools/verification/check_openfoam2d.py out/rectangle-case
```

矩形外流场产品：248 cells、680 points、1084 faces、468 internal faces。
质量门结果：最大 non-orthogonality 26.5651 度、最大内部 skewness 0.5、
最大水力 aspect 3、最短面 0.03125。独立读取结果：最小棱柱体积
`0.00012207031249999992`，最大逐单元面闭合残差 0，PASS。

两次独立生成的五个 `polyMesh` 文件 SHA-256 完全相同：

- boundary `37583ae9bcbad06a020e6a373e44995171292fefb3656623b2e5a68cc883c4b0`
- faces `1fc8340d3b865ead460c570a63e200d7c22b14fb7cf6ef3337bb2ea9b5314a11`
- neighbour `a1fb64d9d2eeb3f4490980a301f186e52ee84587c6391c5e3f3b03104a524461`
- owner `a4449da66ae43a85e2a92a94403517cc3ea0b0f4e3d784f05e90d84e687c9008`
- points `c3100c7f73c15933c62088888ed10008b36dc56f98a83d145a2f63c8feea8a3f`

Release CTest 30/30 PASS；ASan+UBSan CTest 30/30 PASS。

## 诚实失败边界

当前曲线/翼型及部分多环网格没有被伪装成 solver-ready，质量门会在写出前停止：

- circle：8 个内部面 skewness 失败，最大已观察值 24.7407（门限 4）；
- airfoil-like：7 个内部面失败，最大 23.6791；
- NACA2412 dense：10 个内部面失败，最大 380.183；
- two-obstacles：8 个 agglomerated cell 出现 90 度凹量（门限 80）；
- annulus / nested-island：分别 4 个 90 度凹量单元；
- superellipse：内部 skewness 最大约 179.71。

这些失败说明现有 small-cell 邻接合并虽然保持拓扑和面积，但没有把求解器质量纳入
合并目标。后续不能通过调宽门限解决；应保留失败 cell/face，并实现
quality-aware target selection、必要时的局部重分割或 face/cell topology repair。

## V1 入口

V1 先把上述矩形产品交给真实 OpenFOAM `checkMesh`（default、
`-allTopology`、`-allGeometry`），随后构造最小流动算例并验证可收敛性、守恒与
确定性。曲线、多环和翼型仍标为 blocked-on-quality，直到上述质量修复完成并再次
通过真实检查。
