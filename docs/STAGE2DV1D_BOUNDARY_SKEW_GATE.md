# 2D-V1d 边界偏斜度质量门（阶段内检查点）

## 结论

本检查点补齐了 OpenFOAM 边界面 skewness 质量门。旧实现只检查内部面，因而会把
128 段圆柱候选误判为可求解网格；OpenFOAM 2606 `checkMesh -allGeometry` 实际报告
8 个高偏斜壁面，最大值为 `5.7043075454659755`。

新实现逐个物理/外边界 edge 使用与 OpenFOAM 相同的 mirror-cell 定义：先取得
owner-centre 到 face-centre 的切向修正量，再用 `max(0.4*normalDistance,
faceCentreToEdgeInSkewDirection)` 归一化。默认上限保持 OpenFOAM 的 `4.0`，不降低阈值。

## 保留回归

- 最小回归是从 128 段圆柱直接提取的一枚四边形 Cut-cell 和其短壁面。
- C++ 内部值：`5.7043075454659755`。
- OpenFOAM 2606 外部值：`5.7043075454659755`。
- 该产品输入现在在写出 OpenFOAM `polyMesh` 前 fail-closed，报告全部 8 个坏面及其
  owner cell、edge 和顶点坐标。

## 有效对照

96 段、半分段相位圆柱候选保留为有效对照：

- solver cells: `4664`
- internal max skewness: `0.813190`
- boundary max skewness: `3.605656468`
- OpenFOAM 2606 `checkMesh -allGeometry -allTopology`: `Mesh OK`

## 验证

- Release CTest: `40/40 PASS`
- ASan + UBSan CTest: `40/40 PASS`
- 独立 OpenFOAM 2606 对 96 段有效网格确认 `3.605656468`。
- 独立 OpenFOAM 2606 对旧 128 段产品确认 `5.7043075454659755` 和 8 个坏面。

## 尚未完成

这不是 2D-V1d 完成声明。Re=20 圆柱试算已收敛，但不同层级出现约
`|Cl|=0.014--0.038` 的对称性偏差，Cd 也未形成单调网格收敛；这些结果不能作为
可信物理验收。密集 NACA0012 输入还会因尾缘同一 Cut-cell 内包含过多短曲线段而触发
小角度、长宽比、face-weight 和 boundary-skewness 门。下一检查点必须修复或明确重设
几何离散/局部加密策略后，再做圆柱与翼型物理基准。
