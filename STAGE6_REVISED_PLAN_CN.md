# 阶段 6 重划分建议：自适应求解器闭环、质量稳定化与千万级扩展

> 目标：不要直接以“千万级”作为第一实现任务。先把现有阶段 0–5 的能力贯通成一条可验证的 STL → 自适应 Cartesian → Cut-cell → 完整求解器网格链路，再解决质量，最后扩展到千万级。

## 6.0 基线冻结与真实性门

- 冻结阶段 0–5 源码、测试和基准，不混入阶段 6 修复。
- Release 全新构建；阶段 0–5 CTest 必须全通过。
- 固定 5 组回归 STL：立方体、非凸 L 形体、薄壁/内腔、多壳层、Stanford Bunny。
- 明确输入契约：当前只保证封闭、流形、无自交的 STL；脏几何默认拒绝并定位，不假装自动修复。
- 记录 uniform / adaptive / incremental 三条路径当前能否生成完整 solver mesh 的事实表。

通过条件：形成 `stage6_baseline.json` 和人工可读基线表；任何后续改动都必须回归该基线。

## 6.1 统一自适应 Cut-cell 完整拓扑

当前最大断点是 adaptive octree 能生成 Cut-cell，但完整 OpenFOAM `polyMesh` 只支持 uniform Cartesian。

本步只解决：

- coarse-fine 邻接面确定性分裂；
- 2:1 过渡处 owner/neighbour 一致；
- coarse/fine 面覆盖无重叠、无缺口；
- Cut-cell 与普通八叉树叶统一进入完整控制体拓扑；
- patch 连续、边界面唯一归属；
- stable Morton leaf ID 到 solver cell ID 的确定性映射。

通过条件：解析立方体 + 非凸 L 形体的 adaptive 网格可导出完整 OpenFOAM `polyMesh`，独立 reader 通过，OpenFOAM `checkMesh -allTopology` 输出 `Mesh OK.`。

## 6.2 原生质量评估器

在生成器内部增加与外部 `checkMesh` 对应的质量诊断，不能等导出后才知道哪里坏。

至少计算：

- cell volume / closure；
- face area / zero-area face；
- face-pyramid sign；
- skewness；
- non-orthogonality；
- concavity / star-shaped 可行性；
- duplicate face / baffle-like duplicate；
- tiny face / tiny edge；
- min volume fraction；
- 每个失败项的 cell/face ID、位置、来源 background leaf ID。

通过条件：生成器内部失败计数与 OpenFOAM 对固定案例的失败位置/数量趋势一致；质量问题必须可定位、可回归。

## 6.3 小 Cut-cell 与坏形状稳定化

优先顺序固定：

1. 守恒邻接聚合（agglomeration）；
2. 聚合会产生非星形/负 face-pyramid 时拒绝；
3. 局部共形分裂作为回退；
4. 禁止靠删小单元、翻孤立面、改检查阈值“过关”。

必须先建立最小失败案例：tiny sliver、T-junction、coarse-fine + cut、极小体积分数、不可聚合非星形单元。

通过条件：所有最小案例在体积、一阶矩、patch、owner/neighbour、edge manifold、face pyramid 上通过；相同输入得到完全相同的聚合/分裂结果。

## 6.4 复杂几何质量门

分辨率逐级推进，不允许直接上千万级：

- R24：立方体、L 形体、薄壁、多壳层、Bunny；
- R48：重复同一矩阵；
- R96：Bunny + 至少一个更高三角数工程型 STL；
- 每一级都先过内部质量，再过独立 reader，再过 OpenFOAM。

通过条件：`checkMesh` 默认检查和 `-allTopology` 均 `Mesh OK.`；固定验收案例进一步要求 `-allGeometry` 无阻断错误，若存在已知非阻断项必须单独记录，不得混成 PASS。

## 6.5 统一紧凑数据结构

只有 6.1–6.4 质量闭环稳定后，才做大规模优化。

- 普通 Cartesian/Octree 单元保持隐式几何；
- 只为 Cut-cell/transition 保存显式多面体；
- flat arrays / SoA；
- rank/select 或稳定索引替代每单元动态对象；
- 面和点流式生成；
- 不允许“uniform 有一套 scalable writer、adaptive 另有一套完整拓扑”长期分叉。

通过条件：同一网格用参考路径与紧凑路径输出的拓扑 fingerprint、体积、patch 和 solver 文件语义一致。

## 6.6 百万级到千万级递增扩展

规模门：1M → 3M → 10M background cells/leaves。

每一级记录：

- geometry/classification 时间；
- octree/adaptation 时间；
- cut-cell + quality 时间；
- topology 时间；
- export 时间；
- independent read 时间；
- peak RSS；
- points/faces/cells/cut-cells；
- quality worst values；
- deterministic hashes。

任何一级超内存或质量失败，不进入下一级。

## 6.7 千万级终态验收

最终必须同时满足：

- 实际完整 solver control volumes，而非只有 10M 背景叶；
- 负体积、非闭合、共享面不匹配、索引越界、patch 漏面为 0；
- OpenFOAM 真实读取并 `Mesh OK.`；
- 独立 reader 全量读取；
- 两次生成的 topology hash / stable ID / polyMesh SHA-256 稳定；
- 在既定 Mac 资源预算内完成；
- 失败单元不被隐藏或删除。

## 6.8 阶段 6 关闭后再进入混合贴体阶段

阶段 6 不同时实现棱柱边界层，否则范围会再次失控。阶段 7 建议定义为“近壁混合层”：

- 壁面三角面法向与特征边处理；
- 棱柱层挤出与层高/增长率；
- 碰撞、窄缝和尖角层终止；
- prism ↔ Cartesian/octree 过渡；
- y+ / first-layer-height 尺寸输入；
- layer skewness / orthogonality / negative-volume 检查；
- 与 Cut-cell fallback 共存。

这时项目才真正形成“Cartesian/Octree 主体 + Cut-cell + 可选贴体近壁层”的混合网格器。

## Codex 执行规则

- 每次只做一个 6.x 子阶段；通过后再进入下一个。
- 开工前必须列出将修改的文件、已有失败案例、验收命令。
- 每个 6.x 至少一个独立 commit；禁止一个 commit 横跨多个子阶段。
- 任何质量修复必须先加入最小失败回归，再改实现。
- 如果某一步需要大范围重写数据结构，先停下提交设计说明，不直接编码。
- 不允许以“测试绿了”代替真实 STL 输出和外部检查器。
