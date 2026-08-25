# Stage 2D-H3：Solver Topology Scalability

日期：2026-08-26

分支：`agent/native-2d-robustness`

基线：H2 commit `fabc43131e6cee9db8e9fe321db276a8df530b88`

## 1. 52 秒基线拆解

同一份 256 点闭尾缘 NACA0012、101,734 leaves / 102,218 solver cells 的 Release 单线程基线为：solver topology `52.535515 s`，其中 `buildGlobalTopology` 累计 `50.269798 s`、candidate global rebuild `49.461201 s`、full quality `1.190762 s`。一次完整 repair 触发 62 次全局 topology rebuild、68 次 full quality、61 个全局候选；最后只接受 3 次 source repair 和 2 次 repartition。峰值 RSS 为 377,716,736 B。

## 2. 根因

主因不是 Quadtree、Cut-cell 或一次完整 quality scan，而是 61 个局部候选分别复制全部 cell、重建全局 topology、再跑全局 quality，形成 `K × (N+E)`。`sourceMergePairs()` 与 `solverRepartitionPairs()` 还曾为每个 affected cell 扫描全部 edges，但该项不是 52 秒中的最大部分。

## 3. 采用的架构

H3 使用 `global detect → local decision → deterministic independent batch → one global rebuild → authoritative full quality`。局部批次没有严格改善全局质量分数时，保留 H2 exhaustive 路径作为 fail-closed 小规模回退；没有实现高风险的 vertex/edge owner-neighbour 原位 patch。该 transaction/commit 思路参考 OpenFOAM `polyTopoChange` 的“记录动作后一次 apply + map”结构，但没有复制 GPL 实现。

## 4. 邻接结构

生产 topology 已有确定性的 `TopologyCell2D::edges` incident-edge 索引，因此 H3 直接替换两处 `affected cell × all edges` 扫描，没有再复制一份 CSR。当前 profile 表明 repair 决策已不是内存或时间主项；未来三维若需要紧凑 `cell → faces`，再独立评估 flat/CSR，不能把本二维容器选择直接外推。

## 5. 局部质量策略

Repartition 候选使用 changed cells 的闭合一环 halo。代码从现有 topology 的 patch frontier 确定性重建局部边界，把候选两块与一环邻居组成小 topology，并调用与全局相同的 `evaluateSolverQuality2D()`；排序覆盖 aspect、concavity、interior angle、face length、non-orthogonality、internal skewness、face weight 与 volume ratio。人工 patch frontier 的 boundary skew issue 不参与候选排序；真实批次提交后仍必须跑完整全局 quality，local quality 从不替代最终门。

## 6. 批量修复策略

Source merge 与 repartition 都先生成局部 proposal，按局部质量 rank、source/cell ID 稳定排序。每个 proposal 带闭合一环 halo；greedy independent-set 只允许 halo 不相交的 proposal 同批提交。先尝试完整独立批次，若全局分数未改善则确定性缩小前缀，最后才进入 H2 exhaustive 回退。输入相同时 selection、ID、报告及输出保持确定性；100k H3 两次 OpenFOAM `points/faces/owner/neighbour/boundary` 和 `.cm2d` SHA-256 分别完全一致。

## 7. 100k 前后

| 指标 | H2 baseline | H3 | 变化 |
|---|---:|---:|---:|
| solver topology | 52.535515 s | 5.956927 s | -88.66%，8.82× |
| total | 56.930397 s | 10.323162 s | -81.87% |
| global topology rebuild | 62 | 7 | -88.71% |
| full quality calls | 68 | 11 | -83.82% |
| candidate global rebuild | 49.461201 s | 4.758287 s | -90.38% |
| peak RSS | 377,716,736 B | 371,056 KiB | 同量级 |

H3 的 100k 最终 solver cells、faces、quality extrema 与 H2 相同：102,218 cells、204,678 个二维 faces，max non-orthogonality `69.4138145°`、max internal skewness `0.899066272`、max boundary skewness `3.49721098`、min face weight `0.0503858011`、min volume ratio `0.0391280076`，全部维持原阈值并 PASS。

## 8. Scaling

环境：Apple M1 MacBook Air（8 cores，4P+4E）、8 GB RAM，GCC 15.2.0，CMake Release `-O3 -DNDEBUG`；算法路径单进程、单线程。RSS 由 zsh `time` 的 maximum resident set 记录。

| 档位 | leaves | solver cells | solver 2D faces | rebuild / full quality | solver topology | total | peak RSS |
|---|---:|---:|---:|---:|---:|---:|---:|
| ~10k | 16,744 | 16,804 | 33,853 | 14 / 17 | 1.207693 s | 1.752072 s | 77,472 KiB |
| ~50k | 54,754 | 55,348 | 110,630 | 7 / 11 | 2.588816 s | 4.664980 s | 204,064 KiB |
| ~100k | 101,734 | 102,218 | 204,678 | 7 / 11 | 5.956927 s | 10.323162 s | 371,056 KiB |
| ~250k | 290,728 | 290,998 | 582,854 | 7 / 11 | 27.970527 s | 43.518359 s | 832,672 KiB |
| ~500k | 540,724 | 542,312 | 1,084,755 | 7 / 11 | 63.208433 s | 95.059227 s | 1,092,752 KiB |

10k 的粗边界 fixture 需要 H2 回退，所以仍有 14 次 rebuild，但绝对耗时仅 1.21 秒；50k 以上 repair 次数固定为 7。另保留 max-level 10 的 NACA 失败案例：boundary skewness `7.107745 > 4`，正确拒绝且不写 OpenFOAM mesh。

## 9. 剩余瓶颈与外部验证

500k solver topology 的 63.21 秒中，`buildGlobalTopology` 累计 60.39 秒（候选 rebuild 51.81 秒），已经是明确的剩余瓶颈；局部 proposal work 只有 0.181 秒。H3 已达到 100k 理想目标，因此按计划停止，不进入真正 incremental topology。独立、不链接项目库的 `check_openfoam2d.py` 对最大 case 验证：542,312 cells、2,169,379 个三维 faces、最小体积 `2.989628488593857e-10`、最大闭合残差 `2.168404344971009e-19`、issues 为空。真实 OpenFOAM `checkMesh`：**NOT RUN**（本机无该可执行文件，独立 reader 不冒充 checkMesh）。

## 10. 阶段判定

**PASS**。100k solver topology 低于 10 秒并获得 8.82× 加速；质量阈值、最终 full quality gate、面积守恒、确定性和 OpenFOAM 输出语义均未放宽。最大完整成功产品为 542,312 solver cells。H3 到此停止，不进入 H4、boundary layer、1M、3D backend 或完全增量 topology。
