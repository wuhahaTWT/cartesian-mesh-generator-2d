# 2D-V1g M11 质量驱动的局部重分区

## 目标与根因

V1f 的 M11 NACA0012 候选在 source-cell 聚合后仍有两个对称的低 face-weight。保留的实际
前缘 source polygon 有 6 个顶点；既有凸分区在同一 source cell 内生成一个狭小三角形，
人工对角线上的 OpenFOAM face interpolation weight 为 `0.0333906374 < 0.05`。继续合并
相邻 source cells 会重新产生同一条对角线，因此不能严格改善质量。

V1g 增加确定性的局部换对角线：

1. 只从 internal non-orthogonality、internal skewness、low face-weight 和 low volume-ratio
   问题收集相邻 solver-cell 对；boundary skewness 等问题不会借内部重分区掩盖；
2. 精确取消候选两单元的反向共享边，要求合并轮廓简单且面积守恒；
3. 枚举该轮廓所有两片严格凸分区，拒绝 under-determined boundary triangle；
4. 对每个候选重建完整全局拓扑并重新计算完整 solver-quality；
5. 只有 `(问题数, 最大严重度, 总严重度)` 严格下降才接受，最多 32 次。

没有移动物理边界、删除单元或改变任何质量阈值。最小回归直接保留 M11 的失败 solver
polygons，先证明原拓扑不通过，再要求一次以上局部重分区后完整质量通过。

## 被拒绝的路线

曾实验删除人工分区后的精确共线转折点。内部质量门可以通过，但全局共形化会在前后空面
重新插入这些点；OpenFOAM 2606 随后报告 6 个 concave cells，`Failed 1 mesh checks`。
该路线没有进入最终代码。最终实现保留全部共形点，只替换人工内部对角线。

## M11 实际结果

输入仍为 256 点闭尾缘 NACA0012，参数为 max level 11、padding `10.1`、minimum level 8、
boundary simplification fraction `0.1`：

```text
cartmesh2d_cli naca0012.xy out 11 10.1 0.1 exterior openfoam-case 8 0.1
```

- boundary vertices: `256 -> 26`
- measured max boundary deviation: `0.00088559 chord`
- boundary area: `0.0816978 -> 0.0808195`
- source/stabilized cells: `66058 / 66016`
- source-cell agglomerations / local repartitions: `3 / 2`
- solver cells / faces: `66254 / 265530`
- max non-orthogonality: `69.41381453 deg`
- max skewness: `3.497210979`
- min face weight: `0.05038580105`
- min volume ratio: `0.03912800757`
- independent reader: 66254 cells、265530 faces、最小体积 `2.989628489e-10`、最大闭合残差
  `2.168404345e-19`、全局常量通量平衡 `2.211772432e-18`
- OpenFOAM 2606 `checkMesh -allGeometry -allTopology`: `Mesh OK`
- 重复 VTK SHA-256：
  `ab703c83b8f827fb87bccb5240499d93b34f0e2c2a96f625858b620befab1ad7`
- 两次 `points/faces/owner/neighbour/boundary` 均逐字节相同。

Release 与真实 ASan/UBSan 构建均为 `40/40` CTest PASS。一次 Release M11 端到端实测墙钟
为 `34.28 s`；本阶段没有取得可靠峰值内存读数，因此不报告内存性能。M9 输出保持 4318
solver cells、2 次 source-cell 聚合、0 次局部重分区，并与 V1f 字节级结果一致。

## 适用边界

26 点离散和约 `1.08%` 的边界面积变化仍不足以证明可信翼型气动力；M11 只是把更紧几何
推进到独立结构检查和 OpenFOAM 网格质量验收。128 段圆的已知 boundary-skewness 最小回归
仍由内部门拒绝，局部重分区不能把 boundary 问题改写成通过。V1d 下一关仍是圆柱与翼型的
物理基准、网格收敛和对称性，不得把本阶段写成仿真可信度整体完成。
