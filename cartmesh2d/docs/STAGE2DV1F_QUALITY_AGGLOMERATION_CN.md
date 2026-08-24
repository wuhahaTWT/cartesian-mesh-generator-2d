# 2D-V1f 质量驱动的源单元聚合与近共线分区

## 目标

V1e 将边界采样尺度显式绑定到最细 Cartesian 单元尺度后，较紧的 NACA0012 几何仍会在
前缘产生低 face-weight。原凸分区会沿现有 source-cell 边界生成狭小三角形；仅改变三角形
对角线不能消除该问题。本阶段在写 OpenFOAM 前增加一个受质量门约束的确定性修复：

1. 对 solver topology 运行完整 `SolverQuality2D`；
2. 只枚举质量问题单元相邻的 source-cell 对；
3. 精确取消两单元的反向共享边，构造面积守恒的简单多边形并重新凸分区；
4. 仅当 `(问题数, 最大严重度, 总严重度)` 严格按字典序下降时接受候选；
5. 最多接受 32 次聚合；无严格改进即停止，并由既有质量门拒绝输出。

这不是删除坏单元，也不降低 OpenFOAM 对应阈值。聚合后的 source polygon 可以是凹的，
但最终 solver cells 仍由确定性凸分区产生，之后重新构造和审计全局拓扑。

另外，凸性判定加入归一化转角下限 `cross/(|e0||e1|) > 1e-10`。原因是一个真实回归中，
精确二进制方向谓词把极小正转角视为凸，OpenFOAM 的 face-plane 检查却把对应挤出棱柱判为
concave。现在该近共线转折会被显式分区，不再塞进单个 solver cell。这个判断是收紧几何门，
不是放宽质量门。

## 保留的最小回归

`solver_export_test.cpp` 保留 NACA 前缘实际失败中的两个相邻 source cells。旧分区在共享界面
产生 `0.0352839` 的 face interpolation weight，低于 `0.05`。回归检查：

- 发生且只发生一次质量驱动的 source-cell 聚合；
- 修复后 `SolverQuality2D` 完全通过；
- 每个最终单元的归一化正转角均大于 `1e-10`，覆盖 OpenFOAM concave-cell 回归。

## NACA0012 M9 验证

输入为 V1e 生成的 256 点闭尾缘 NACA0012。参数为 max level 9、padding `10.1`、minimum
level 6、boundary simplification fraction `0.1`：

```text
cartmesh2d_cli naca0012.xy out 9 10.1 0.1 exterior openfoam-case 6 0.1
```

结果：

- boundary vertices: `256 -> 14`
- boundary measured max deviation: `0.00325711 chord`
- source/stabilized cells: `4238 / 4232`
- quality-driven source agglomerations: `2`
- solver cells: `4318`
- max non-orthogonality: `67.38893558 deg`
- max skewness: `2.548490417`
- min face weight: `0.1260094561`
- min volume ratio: `0.09007454899`
- OpenFOAM 2606 `checkMesh -allGeometry -allTopology`: `Mesh OK`
- 独立 reader：4318 cells、17390 faces、最小体积 `1.822431107e-7`、最大闭合残差
  `3.469446952e-18`、全局常量通量平衡 `3.642919300e-18`
- 重复 VTK SHA-256：
  `e37bcc1164350b2e0639253d637b4731001028888c36c22315acf4aeae1a6abb`
- 两次生成的 `points/faces/owner/neighbour/boundary` 均逐字节相同。

Release 与 ASan/UBSan 构建均为 `40/40` CTest PASS。128 段圆的已知 boundary-skewness
失败仍被拒绝：8 个问题、最大值 `5.7043075454659755`，没有写 OpenFOAM 网格。

## 未解决边界

更紧的 M11、相同 padding/fraction 输入保留 26 个边界点，但仍有两个对称的低 face-weight：
`0.0333906374 < 0.05`。现有候选聚合不能严格改善质量分数，因此保持 0 次接受并拒绝输出。
所以 V1d 仍未完成；本阶段只证明 M9 局部前缘失败已被真实修复。下一阶段必须为 M11 增加更
局部的受约束重分区/多源聚合，再进入圆柱与翼型物理基准。
