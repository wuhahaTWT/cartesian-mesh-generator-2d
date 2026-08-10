# 阶段 6 计划：千万级完整网格与产品化验证

日期：2026-08-10（Asia/Shanghai）

状态：**已获用户批准并启动，验收门禁未关闭**

## 1. 唯一目标

在真实复杂封闭几何上，生成一千万量级的完整 Cartesian Cut-cell 体网格，完整产物必须同时包含：

- 正体积流体控制体；
- Cut-cell 显式多面体几何；
- 内部 owner/neighbour 拓扑；
- 远场和嵌入壁面 boundary patch；
- 确定性稳定背景单元 ID 和结果 hash；
- 闭合、体积、共享面、边界覆盖、小单元和区域质量报告；
- 真实求解器网格导出和独立外部读取。

阶段 2 的一千万八叉树叶只能作为背景索引性能基线，不能作为本阶段的完成证据。

## 2. 本阶段不做

- 不加入 CFD 求解器、GUI、云服务、AI 网格或部署；
- 不引入 MPI、分布式调度或 GPU 生成；
- 不将体素标签、表面相交数、VTK 预览或背景叶数称为完整求解器网格；
- 不通过丢弃 Cut-cell、关闭拓扑/质量检查或隐藏小单元获得成功；
- 不把研究用 Stanford Bunny 资产的验证使用条款扩大为产品商用授权；
- 阶段 7 在阶段 6 终态门禁关闭前不定义、不启动。

## 3. 开工审计与必要架构转换

当前 `ConvexCutCellMesh` 为每个普通全流体单元也保存六个面对象及多个动态容器；
`OpenFoamWriter` 还在写出前同时保留全部面、顶点和 owner/neighbour。这些路径用于中小网格时是可验证基线，
但无法在当前 8 GiB MacBookAir10,1 上直接容纳一千万单元。

阶段 6 因此新增紧凑路径：

1. 普通全流体/全固体背景单元用打包状态和隐式 Cartesian 几何表示；
2. 只为真正 Cut-cell 保存显式凸片、嵌入面和局部分量；
3. 使用表面三角片到网格的稀疏栅格化和非切割区域连通分类，不对一千万个单元各做一次全表面射线；
4. 用紧凑 rank 索引从稳定背景 ID 确定性映射到求解器单元 ID；
5. OpenFOAM 完整体网格使用 32 位 label、64 位 scalar 的本机小端 binary list 分段/流式写出，不在内存中保留三千万张面；
6. 完整大网格与可交互的局部预览分开：预览只是检查入口，完整 binary `polyMesh` 才是导出证据。

## 4. 分级验收矩阵

在直接运行一千万级前，同一紧凑路径必须通过递增规模：

1. **微型拓扑真值**：解析立方体/平面切割，逐面 owner/neighbour、patch 范围、体积和稳定 ID 可手工核对；
2. **中等外部检查**：万至百万级完整 binary `polyMesh` 由 OpenFOAM `checkMesh -allTopology` 读取；
3. **复杂几何收敛**：Stanford Bunny 在至少三个递增分辨率上体积/边界面积趋势稳定，拓扑失败为零；
4. **千万级终态**：背景 Cartesian 单元数不少于 10,000,000，完整求解器体单元、内部面、边界面和点均必须有实际数量与导出文件；
5. **外部终态**：独立 OpenFOAM 运行时或另一个独立二进制读取器完整读取千万级产物，并核对列表计数、索引范围、面引用和 patch 连续性。

中等规模 OpenFOAM 通过不能替代千万级完整文件的独立读取；千万级自建读取器通过也不能消除中等规模的 OpenFOAM 实际兼容验证。

## 5. 完整通过条件

- 背景单元数 `>= 10,000,000`，报告不将其与实际求解器体单元数混为一谈；
- 每个流体单元都有正体积，每张内部面有合法 owner/neighbour，每张边界面恰好属于一个 patch；
- 非闭合、负体积、共享面不匹配、分类冲突、未覆盖面和越界索引均为零；
- 小 Cut-cell 只统计不自动删除，数值零体积片保留数量和总体积记录；
- 同一输入的紧凑状态 hash、网格计数、OpenFOAM 文件 SHA-256 和独立读取报告稳定；
- 局部预览由 meshio 或 ParaView 实际打开，但预览不代替完整体网格证据；
- 中等规模 OpenFOAM 默认 `checkMesh -allTopology` 必须输出 `Mesh OK.`；仅格式读取和核心
  拓扑通过、但 face pyramid 或 skewness 失败时，不得声明求解器可用或关闭阶段 6；
- 最终报告分开几何分类/Cut-cell、拓扑/质量、导出和外部读取墙钟时间；
- 性能记录同时包含峰值 RSS、线程数、硬件、操作系统、编译器和 Release 构建类型。

## 6. 本机资源门槛

开工时实测环境是 MacBookAir10,1、8 GiB 内存，工作盘可用空间约 9.1 GiB。因此终态预算固定为：

- 生成 + 导出进程峰值 RSS `<= 6.0 GiB`；
- 完整阶段 6 终态产物新增磁盘占用 `<= 4.0 GiB`；
- 一千万级生成墙钟 `<= 15 min`，完整导出墙钟 `<= 15 min`；
- 一千万级独立顺序读取和拓扑计数检查 `<= 15 min`；
- 运行线程数第一版固定为 1，不把操作系统或 OpenFOAM 内部线程误算为生成器并行。

任一资源门槛失败都要保留报告并继续定位，不能只提高数字后关闭阶段 6。

## 7. 预定证据产物

- `artifacts/stage6_acceptance.json`：终态机器契约；
- `benchmarks/baselines/stage6_10m_*.json`：千万级时间/内存/计数/hash；
- `artifacts/stage6_10m_case/constant/polyMesh/`：完整 binary OpenFOAM 体网格；
- `artifacts/stage6_10m_external_reader.json`：完整文件独立顺序读取证据；
- `artifacts/stage6_medium_checkmesh.json`：OpenFOAM 实际兼容证据；
- `artifacts/stage6_preview.vtu` 与预览图：可交互局部检查入口；
- `docs/STAGE6_VERIFICATION.md`：命令、失败案例、外部证据、真实边界和终态结论。
