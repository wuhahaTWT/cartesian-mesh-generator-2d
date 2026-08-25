# Stage 2D-H2：大网格可扩展性验证

日期：2026-08-25  
分支：`agent/native-2d-robustness`  
基线：H1 commit `cb16bdf71efe33fba7ce05683fa81977f6890711`

## 1. 环境与口径

- Apple M1 MacBook Air，8 cores（4P+4E），8 GB RAM，macOS 26.5.2；单进程、单线程算法路径。
- GCC 15.2.0，CMake `Release`（`-O3 -DNDEBUG`）。
- `/usr/bin/time -l` 记录墙钟与 maximum resident set size；阶段计时使用 `steady_clock`。
- 输入为脚本确定性生成的 256 点闭尾缘 NACA0012；`padding=10.1`、boundary simplify cell fraction `0.1`，并同时使用 H1 `minimumLevel`、distance band 和 downstream refine box。
- `leaf`、流体 source、几何 Cut-cell、稳定化 cell、solver cell、edge/face 与 vertex 分开统计；没有构建 solver topology 的档位明确写 `NOT RUN`。
- `sizeof(QuadtreeLeaf2D)=80 bytes`。leaf 仍保存可重建 AABB 和分类；本轮没有引入 pointer tree。

## 2. P0 结论：flat sort+sweep 被实测否决

用户建议的单一 `FaceRecord` 连续数组方案经过三次实现/压缩实验：40-byte tuple sort、32-byte packed comparator、16-byte logical-field packing。486,034 leaves 下最好一次邻接仍为 `0.563 s`，慢于原坐标桶 `0.395 s`，且进程峰值 RSS 高于原实现。原因是本二维案例只有约 4097 个 lattice 坐标桶，而全局排序需处理 1,944,136 条 face records；`O(N log N)` 常数与临时连续存储在此处不占优。

因此没有提交性能回退，也没有声称 flat 一定优于 `std::map`。保留的优化是复用同一 balance 状态已经生成的邻接集合，去掉“统计 before / 本轮 / 确认 empty / 统计 after”之间的重复全量构建；算法、邻接顺序与 2:1 结果不变。

同一 486,034-leaf 配置的 H1 基线与 H2 结果：

| 指标 | H1/旧调用流 | H2 | 变化 |
|---|---:|---:|---:|
| neighbor calls | 4 | 2 | -50% |
| neighbor generation | 0.395489 s | 0.207060 s | -47.65% |
| balance total | 0.556071 s | 0.365792 s | -34.22% |
| full non-solver product chain | 24.228625 s | 23.786456 s | -1.82% |
| maximum RSS | 811,433,984 B | 806,420,480 B | -0.62% |

旧/新 `.cm2d` SHA-256 均为 `270fb6b494b5afdc80a46dc83bd7b204270dfd319f5c655bb5db4a1fcadbf3e0`；旧/新 VTK SHA-256 均为 `5205ed05c2a72330f711b0029458c66f6b318e5f817a25f9162d53066fa2f82f`。这证明输出没有因性能改动漂移。

## 3. 规模矩阵

| 档位 | leaves | fluid source | Cut | stabilized | solver cells | solver 2D faces | vertices | total | max RSS |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ~10k full | 16,744 | 16,718 | 56 | 16,712 | 16,804 | 33,853 | 17,049 | 1.586 s | 79,052,800 B |
| ~100k full | 101,734 | 101,476 | 220 | 101,434 | 102,218 | 204,678 | 102,460 | 56.584 s | 389,251,072 B |
| ~250k core/export | 244,540 | 243,874 | 220 | 243,832 | NOT RUN | 489,158 source edges | 245,326 | 10.179 s | 608,763,904 B |
| ~500k core/export（P0 对照） | 486,034 | 485,368 | 220 | 485,326 | NOT RUN | 972,608 source edges | 487,282 | 23.786 s | 806,420,480 B |

最终阶段门另跑一个略大配置，得到 504,250 leaves；其结果写在第 6 节，避免用“约 500k”替代明确阈值。

## 4. 分阶段时间

| 阶段 | ~10k full | ~100k full | ~250k core | 486k core |
|---|---:|---:|---:|---:|
| quadtree refinement | 0.069 | 0.501 | 0.722 | 1.665 |
| 2:1 balance | 0.009 | 0.068 | 0.427 | 0.366 |
| neighbor（包含于 balance） | 0.006 | 0.041 | 0.205 | 0.207 |
| Cut-cell | 0.017 | 0.295 | 0.722 | 1.406 |
| source topology | 0.082 | 0.816 | 2.722 | 7.313 |
| small-cell analysis/agglomeration | 0.102 | 0.967 | 3.110 | 8.105 |
| solver topology/quality repair | 1.082 | 52.492 | NOT RUN | NOT RUN |
| serialization/export/readback | 0.221 | 1.434 | 2.455 | 4.873 |

真正瓶颈已从假设变为证据：100k 完整链中 solver topology/quality repair 占约 92.8%。它会为少量边界质量候选重复重建全局 topology；直接把这一模式放大到 500k，在 8 GB 机器上风险过高。H2 不在没有局部重建设计和最小回归的情况下仓促改写这条已通过真实质量门的 V1g 路径。

## 5. 正确性与外部验证

- flat 实验加入独立 `O(N^2)` 小网格 oracle，逐对比较 coarse/fine face neighbor 集合；该回归保留，防止未来布局优化漏邻接。
- 486k 性能改动前后 CM2D/VTK 字节级一致；2:1 violations after 为 0，确定性 Morton/ID 门保持。
- 100k OpenFOAM case 由不链接项目库的 `check_openfoam2d.py` 验证：102,218 cells、409,114 3D faces、204,920 points，最小体积 `2.989628488593857e-10`，最大闭合残差 `8.673617379884035e-19`，issues 为空。
- 真实 OpenFOAM `checkMesh`：**NOT RUN**。本 macOS 环境没有 `checkMesh` 可执行文件；独立 reader 不冒充 OpenFOAM。
- 保留失败证据：较粗 max-level 10 NACA 候选被 solver boundary skewness `7.107745 > 4` 正确拒绝，没有放宽阈值或写出 OpenFOAM mesh。

## 6. 阶段收尾记录

- H2 targeted regression：`2/2 PASS`，`cartmesh2d_stage2_quadtree_tests` 加一个 fail-closed 的域外 sizing box 案例。
- full `cartmesh2d_` CTest（只跑一次）：`53/53 PASS`，14.21 s。
- representative OpenFOAM checkMesh：NOT RUN（工具缺失）
- final large benchmark：`504,250 leaves / 503,584 fluid source / 220 Cut / 503,542 stabilized / 505,534 vertices / 1,009,076 source edges`；2:1 violations after `0`，25.205865 s，maximum RSS 807,469,056 B。solver topology 为 `NOT RUN`。

阶段判定：**PASS（按能力边界）**。`>=100k` 完整 solver/export 链与 `>=500k` core topology/serialization 均已达到；1M stretch 未运行。500k solver topology 不是本次已证明能力，后续应单独设计局部 topology/quality 更新，而不是把 100k 的全局重建直接扩大。
