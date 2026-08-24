# 2D-V1e 几何采样尺度与网格尺度绑定

## 目标

密集曲线输入可能让一个最细 Cut-cell 同时包含许多远短于单元尺度的边界段。二维挤出后，
这些段会成为极短壁面，并在凸分区中产生小角度、低 face-weight、低 volume-ratio 或高
boundary-skewness 单元。本阶段增加显式、可审计的几何简化路径；默认值仍为 `0`，不会
静默改变输入几何。

CLI 最后一个位置参数为：

```text
[boundary-simplify-cell-fraction=0]
```

实际最大偏差请求值等于该 fraction 乘以最细 Cartesian 单元尺度。实现对闭环使用确定性
Douglas-Peucker：以字典序最小点和距其最远点为锚，分别处理两条闭链。因此循环移动输入
起点不会改变结果。每次简化后重新检查：

- 每个 loop 的简单性和非退化性；
- 多环不相交；
- 每个 loop 的 nesting depth 与原输入一致；
- 原输入顶点到简化折线的实测最大偏差不超过实际容许值。

若拓扑不再成立，容许值逐次减半；无法保持时失败，不输出网格。成功时 CLI 明确报告请求/
实际/实测偏差、原始/简化顶点数和前后区域面积。

## NACA0012 密集输入结果

可复现输入由 `tools/verification/generate_benchmark_geometry.py` 解析生成：闭尾缘 NACA0012，
每侧 128 个间隔，共 256 个 boundary vertices。M10、padding 10、minimum level 7、fraction
0.30 的结果为：

- requested deviation: `0.00615234375 chord`
- measured max deviation: `0.00488642 chord`
- vertices: `256 -> 10`
- area: `0.0816978 -> 0.0770833`（几何变化被显式报告）
- solver cells: `16798`
- max non-orthogonality: `30.20952065 deg`
- max skewness: `1.445175116`
- min face weight: `0.09445386452`
- min volume ratio: `0.06375874566`
- 独立 reader：闭合残差 `3.47e-18`，全局常量通量平衡 `-4.25e-18`
- OpenFOAM 2606 `checkMesh -allGeometry -allTopology`: `Mesh OK`
- 重复生成的 VTK SHA-256：
  `9d14bda4f7c4505e1670a251c9ef2c950243c678ba701bf195f08c402bd2d56d`

## 适用边界

这个路径解决“输入采样远细于目标网格”的尺度失配，但 10 点翼型只适合检查拓扑、守恒和
流程，不足以证明可信黏性翼型气动力。将偏差收紧到约 `0.13% chord` 时可保留 20 点，
当前仍有 8 个 solver-quality 问题。因此 2D-V1d 继续保持未完成：后续必须解决局部 Cut-cell
分区/稳定化和对称性收敛，再运行圆柱及翼型物理基准，不能把本阶段的 `Mesh OK` 写成气动
验收。
