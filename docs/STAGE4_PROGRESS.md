# 阶段 4 进度：工业几何鲁棒性

状态：**已完成并通过阶段四验收；尚未进入阶段五或六**

阶段四在阶段三真实 Cut-cell 几何的基础上，补齐了常见工业 STL 失败诊断、薄壁和
小孔保护、嵌套空腔、多部件/多流体区域、边界命名、显式区域拓扑以及公开复杂几何
验证。默认路径不修补、不删除、不缝合用户 STL；所有容差动作和被丢弃的数值零体积
片均写入报告。

主要结果：

- BVH 粗筛加精确 triangle-triangle 窄相，区分部分重叠、真自交和非邻接接触；
- 非封闭、非流形、方向冲突、退化/极瘦面、重复面、极小分量均有计数、示例位置和
  VTP 标记；
- 支持互不相交外壳和方向正确的嵌套空腔，并为局部流体分量建立确定性的全局区域 ID；
- 支持稳定的 boundary ID/名称与 fluid region ID/名称；
- 薄壁、小孔和狭缝端到端经过自适应细化与 Cut-cell 验证，不能达到指定 gap 分辨率时
  严格模式非零失败；
- 小 Cut-cell 统计包含 background cell ID、质心、体积分数和 boundary ID；
- Stanford Bunny 二进制 STL 完成全量生成、meshio 检查、ParaView/VTK 读取和无效
  tetra 单元计数验证；
- 包含两个不连通流体区域的薄壳完整 OpenFOAM `polyMesh` 保留
  `outer_wall`/`cavity_wall`/`farfield`，并由 OpenFOAM 2606 `checkMesh`
  判为 `Mesh OK`；
- 大型显式几何可写为可验证的四面体分解；复杂 `VTK_POLYHEDRON` 仍保留为调试输出，
  不把其通用 VTK 凸性检查限制隐藏成通过；
- Release、Debug、ASan+UBSan 与独立外部检查均通过，输出和 region ID 确定。
- 生成报告使用 `cartmesh-stage4-cutcell-v1` 记录内部几何验证；外部
  `checkMesh` 终态证据使用 `cartmesh-stage3-openfoam-checkmesh-v1`，薄壳记录显式写出
  `stage4Complete=true`、`solverReadyCutCellMesh=true`。

完整验收矩阵、数值结果、公开数据来源和复现命令见
[`STAGE4_VERIFICATION.md`](STAGE4_VERIFICATION.md)。

阶段边界：本阶段只完成工业几何鲁棒性、网格几何/拓扑，以及为关闭阶段三硬门禁所需的
最小完整 OpenFOAM 体网格输出。CGNS、自适应 OpenFOAM、更多求解器交换格式与质量策略
属于阶段五；增量更新与千万级并行属于阶段六，均未开始。
