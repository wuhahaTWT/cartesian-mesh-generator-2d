# 阶段 6 止损复盘：不要重走写出器内修网格分支

## 1. 当前结论

截至 2026-08-10，阶段 6 尚未完成，阶段 7 未启动。

项目已经回退到较可靠的阶段 6 基线：正式写出器只保留
`connected_component` 控制体路径、流式 binary OpenFOAM 写出、确定性点面编号、
显式面积守恒检查和已有拓扑封口诊断。实验性的逐凸片控制体、写出后多元聚合、
内核四面体拆分及其命令行入口、专用测试和质量统计已经删除。

这不是删除整个阶段 6。以下成果仍然有效并必须保留：

- 阶段 0–5 的全部实现与验收证据；
- 阶段 6 的 `CompactUniformCutCellMesh`；
- R216 千万背景单元完整 binary `polyMesh`；
- 独立全量 binary reader 与逐单元边二关联检查；
- 两次确定性复跑和性能/RSS 记录；
- R96 OpenFOAM 2606 的真实质量阻断证据。

正式状态仍是 `blocked_openfoam_quality`：自建不变量和独立读取通过，不等于
OpenFOAM 默认质量通过。

## 2. 已经证实走不通的路线

失败路线的共同结构是：

1. 先把每个 Cut-cell 的每个凸片直接写成 OpenFOAM 控制体；
2. 在写出器里配对分割面、消除局部非二关联边；
3. 再以 wrong face pyramid 为目标做 2–6 元局部合并；
4. 对剩余缺陷单元使用内核点四面体拆分；
5. 继续调阈值追逐最后少数错误。

这条路线不应恢复。它把本应在控制体构造阶段解决的几何质量问题推迟到序列化阶段，
导致单元、面和状态急剧膨胀，局部修补之间互相制造新约束。

## 3. 可核对的失败数据

### 3.1 逐凸片 R24 基线

输入为 Stanford Bunny，R24，Release，单线程，MacBookAir10,1：

| 指标 | 实测 |
|---|---:|
| 背景单元 | 13,824 |
| OpenFOAM 控制体 | 303,236 |
| 点 / 面 | 457,279 / 1,067,520 |
| 写出总时间 | 84.12 s |
| 进程总时间 | 101.24 s |
| 峰值 RSS | 1,264,533,504 B |
| 独立 reader 非二关联单元边 | 0 |
| wrong face pyramids | 379 |
| highly skew faces | 6,886 |
| 最大 skewness | 1926.56 |
| duplicate baffles | 51 |
| non-standard edge connectivity faces | 102 |
| 最小单元体积 | 4.55696e-28 |
| OpenFOAM 结果 | `Failed 4 mesh checks.`，无 `Mesh OK.` |

它证明“拓扑可读取、逐边二关联”仍不足以得到可求解网格。

### 3.2 继续修补后的止损点

中止前的开发诊断依次观察到：

| 尝试 | 本轮接受 | wrong pyramids | 涉及缺陷单元 |
|---|---:|---:|---:|
| 安全二元合并 | 207 | 379 → 152 | 272 → 65 |
| 三元合并 | 43 | 152 → 64 | 65 → 22 |
| 四元合并 | 10 | 64 → 31 | 22 → 12 |
| 五元合并 | 2 | 31 → 27 | 12 → 10 |
| 六元合并 | 0 | 27 | 10 |
| 内核共形拆分 | 6 个源单元 → 218 个输出单元 | 27 → 14 | 10 → 4 |

尽管 wrong pyramids 降到 14，OpenFOAM 仍报告约 6,807 张高偏斜面、51 张重复
baffle，最大非正交约 123.632，仍为 4 项检查失败且无 `Mesh OK.`。放宽内核接受
条件的后续尝试产生退化内部面，已中止，未作为正式证据保留。

这些后期数字是中止前开发诊断，不是阶段验收结果；保留它们的目的仅是避免下一账号
再次花数小时重复同一条分支。

## 4. 根因判断

- 逐凸片会产生远低于数值尺度的微小控制体，R24 已出现约 `1e-28` 体积；
- 单独保证每条单元边二关联，不能保证 cell center 对所有 face pyramid 为正；
- 非星形组件不是靠翻面就能修复，局部合并又可能形成更强的非星形性；
- 内核四面体拆分会大幅增加单元和内部面，并保留原有高偏斜/重复接口问题；
- 在写出器里同时承担几何构造、质量优化和 OpenFOAM 序列化，缺乏稳定的不变量边界；
- R24 单次准备已接近两分钟，继续放大到 R96/R216 没有资源可行性。

## 5. 明确禁止下一模型重复的动作

- 不恢复 `--convex-piece-cells` 或 `convex_piece_exact`；
- 不在 `ScalableOpenFoamWriter.cpp` 内继续添加 2–N 元搜索式聚合；
- 不通过增加组合元数、放松内核阈值或反复拆分追逐最后几个 wrong pyramid；
- 不删除微小流体片、不把缺口改成 wall、不翻转孤立面伪造通过；
- 不放宽或跳过 OpenFOAM 默认 `checkMesh`；
- 不把独立 reader 的 pass 描述成 solver-ready；
- 在 R24 得到 `Mesh OK.` 前，不重跑 R48/R96/R216 大规模实验。

## 6. 如果以后继续阶段 6，建议的新起点

新方案必须位于 OpenFOAM 写出之前，而不是在写出器里修补。建议先单独设计
“控制体图/稳定化层”：

1. 图节点持有源凸片成员、区域 ID、体积、一阶矩和边界片集合；
2. 图边只表示真实正面积共享面，并保留确定性稳定 ID；
3. 小单元聚合允许跨相邻背景单元，但只允许同一流体区域且必须严格守恒；
4. 候选接受前先重建局部控制体并检查正体积、face pyramid、重复面、边二关联、
   边界 ID 和面积守恒；
5. 聚合结果产生明确的“源凸片 → 最终 solver cell”映射；
6. 写出器只序列化已经通过质量门的控制体图，不再参与质量搜索。

第一里程碑只做合成最小案例和 R24 Bunny。关闭条件是独立 reader 通过且 OpenFOAM 2606
明确输出 `Mesh OK.`；达到前不要进入更高分辨率。

## 7. 下一账号首次接管

```sh
cd /Users/Zhuanz/Desktop/网络生成器

sed -n '1,260p' docs/STAGE6_FAILED_BRANCH_HANDOFF_CN.md
jq '{status,stage6Complete,solverReadyCutCellMesh,gates,acceptanceBlockers}' \
  artifacts/stage6_acceptance.json

cmake --build build/release --target cartmesh_stage6_tests cartmesh_stage6_cli --parallel 4
build/release/cartmesh_stage6_tests
```

不要在首次接管时重跑千万级 case。先确认源码中不存在
`convex_piece_exact`、`quality_conformal_split` 或
`quality_constrained_cell_merge`，再决定是否另开控制体图设计。
