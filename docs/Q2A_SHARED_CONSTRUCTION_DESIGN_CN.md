# Q2-A：共享交点构造与公共分割

实施基线：`dde1c26eceede55364dca53d3c9ecac2d700b825`。
本阶段不承诺 Q2-B 的 narrow-gap / sharp-trailing-edge 短面修复完成。

## 契约

- 原始 wall、H4 layer geometry、patch identity 与所有质量阈值不变。
- 同一 construction context 中，线段/网格线交点按实体键缓存，不重复求解。
- 网格线身份为 axis + 最大层级整数坐标；不能用全局浮点 rounding 合并近邻壁面。
- Cut-cell 携带 canonical vertex handles；公共边分割返回相同 handles 的有序序列。
- topology 的输出编号可以重新紧凑化，但保留 construction handle 映射。
- 旧无 context API 保留为对照；上下文不可混用，不使用全局变量或 thread-local 状态。
- 数值相等、允许移动与 Q1 质量阈值分离。Q2-A 不扩大 snap 范围。
- 局部 half-edge 和 H4 算法保留；solver 质量候选仍执行 authoritative full audit。

## 借鉴边界

独立实现 p4est 的逻辑身份、AMReX 的共享边数据和 Gmsh 的实体关联思想；不复制或链接外部实现。
不采用 cfMesh 的删格预处理、AMReX 的 small-cell covering 或 OpenFOAM relaxed quality。

## 验收

共享交点缓存和公共边分割必须进入真实 hybrid/solver 链路，而非仅导出元数据。
测试应覆盖交点重复请求、coarse/fine、端点/切触、窄缝独立 support、变换和确定性；
五例比较质量、面积、boundary classification，并真实运行 OpenFOAM checkMesh。
性能分别记录构造/拓扑成本与规模；历史 H3 数据不能冒充当前重构测速。

## 实际接入位置

1. `geometry/IntersectionConstruction2D.cpp`：给现有 `IntersectionRegistry2D`
   增加共享构造 API。线段按无向端点身份与 source 类型注册，计算方向按坐标字典序
   固定；交点事件键是 support + axis + 最大层级整数网格坐标。相邻叶子重复请求
   同一事件时不再重新求解。精确坐标 interning 是当前已有 polygon API 的适配桥，
   不是全局近邻聚类。
2. `cutcell/CutCell2D.cpp`：Liang–Barsky clipping 保留进入/离开的 grid side，
   在局部 polygon 提取之前调用共享事件缓存；生成的 CutCell 携带 canonical handles
   及非拥有的 registry 身份保护。不同 context 或失效 handle 显式拒绝。
3. `topology/SharedEdgePartition2D.cpp`：从本次输入 cell 的 active handles 建立二维
   范围索引。同一无向端点区间的有序分割只构造一次，正反向 cell 共用这份结果。
   coarse/fine 的不同端点区间从同一 active 集合查询，得到一致的原子片段。
   未被接受的候选点即使仍留在 registry，也不会进入当前分割。
4. `topology/Topology2D.cpp`：共享路径消费 handles，不再调用 legacy 的近点
   `collectCanonicalPoints/findVertexId`。输出坐标排序后的紧凑 ID，另保留
   `canonicalVertexIds` 映射和 partition 缓存计数。
5. `hybrid/HybridMesh2D.cpp`：一个 context 贯穿 remainder、统一 H4 topology 和
   重复构造校验；H4 原有 `wallVertexKinds` 用来保护物理 feature，不新增角度判断。
6. `quality/SolverTopology2D.cpp`：初始 solver partition、全局 merge/repartition
   与候选完整重建均继续使用同一 context；最终 solver mesh 保留其所有权。
   `apps/cartmesh2d_hybrid_cli.cpp` 默认启用，末尾 `--legacy-construction` 保留
   实际旧构造路径。新增 `*.hybrid.construction.json`，不是以报告替代算法接入。

### 精度与审计

共享构造只使用集中 policy 默认值 `64 * machine_epsilon` 的算术位移预算，乘以
query、incident endpoints、已配置最细网格的局部 h 与 support 长度最小值。不受 legacy 较大采样
snap fraction 的放宽影响。横/竖求交落到同一网格角点时先解析整数角点身份，避免
在缩放后留下一个浮点舍入量的两个点。若此规范化会焊接非相邻 support，显式报错。

每个缓存事件（包括零位移事件）记录 source 类型、源线段、support ID、整数 grid
line、原始求解坐标、canonical handle/坐标、displacement、local_h 和 feature。
`solver_vertex_handles` 对应实际 solver 顶点；事件点不一定在最终 solver 中保留。
ID 是 context 内的固定遍历顺序 ID，不承诺不同并行调度下仍有相同事件编号。

### 本阶段没有做的事

- 原 superellipse `9.79753e-9` 面来自 transition 采样残差与 grid-line 交点的真实
  近重合，不是仅仅重复求交的舍入。本轮保留 `dde1c26` 的源头 transition resampling；
  不把这项既有修复重新归功于共享缓存。
- 局部 half-edge 提取仍使用原 polygon 算法；H4/solver polygon 适配器仍用精确坐标
  回绑 handle。一般边的首次公共分割仍需要几何谓词查询，尚不是完全预先排序的
  全域 arrangement，更不是 DCEL。边界 patch 判定仍由既有几何审计完成。
- 小型局部质量排名候选仍可用 legacy topology，接受前必须经过共享全局重建和
  authoritative quality audit。尚未实现 Q2-C 的增量局部事务或并行 registry。
- 不可安全归一的几何冲突显式失败；通用的 refinement / wall sampling 回退仍是
  Q2-B 工作。本阶段没有宣称能自动修复所有失败，也没有把 pure Cut-cell fallback
  计为 hybrid 验收成功。旧几何谓词的 absolute tolerance 尚未整体重写。

## 五案例质量对照

对照为当前 Q2 部分修复基线 `dde1c26`，不是最初 Q1 基线。hard limit 仍为 `0.01`。

| case | Q2-A 前 face/local_h | Q2-A 后 face/local_h | solver cells 前/后 |
|---|---:|---:|---:|
| circle | 0.045942780858514 | 0.045942780858514 | 728 / 728 |
| superellipse | 0.016222929810346797 | 0.016222929810346797 | 795 / 795 |
| concave L | 0.014400000000000546 | 0.014400000000000546 | 5452 / 5452 |
| narrow gap | 0.007843137254903055 | 0.007843137254903055 | 3189 / 3189 |
| sharp trailing edge | 0.0004202383138647292 | 0.0004202383138647292 | 3391 / 3391 |

superellipse 仍无原 `~9.8e-9` internal face，最短绝对 face 为 `0.0019163335838472152`。
narrow gap 与 sharp trailing edge 仍不满足短面 hard contract；五例完整 Q1 总状态
也仍是 FAIL。Q2-A 的无退化不能登记为 Q2 全部验收。

机器可读质量/确定性/来源对照见 `artifacts/q2a/comparison.json`，实际生成文件在
`build-q2a/evidence/`。预览 `artifacts/q2a/solver-meshes.png` 直接从这些 solver
CM2D 文件绘制，包含五例既有 FAIL 标记，不隐藏微短边或删除 cell。

## 复现

最终原生回归为 **75/75 PASS**（74 项既有测试 + shared construction 测试），
最后一次完整运行 56.09 秒。新增回归覆盖三种量纲尺度、coarse/fine 交点复用、
网格角点、端点/共线、近邻不同 support、特征先/后注册、不变 feature、过短 support
不可被 coarse h 吞并、失效/跨 context handles、无效候选隔离及真实 Cut-cell 共形面积。
circle / superellipse 整体放大 1000 倍后，无量纲分布最大绝对差分别为
`7.442935157087049e-13` / `2.8350655156827997e-12`，typed 状态不变。

隔离公共 topology 的三次进程运行中位数（**不是整个生成器的耗时**）：

| 单元数 | legacy | shared | 几何、编号、owner/neighbour、patch |
|---|---:|---:|---|
| 10,000 | 0.044566 s | 0.034733 s | 完全相同 |
| 99,856 | 0.982997 s | 0.449512 s | 完全相同 |

shared 计时包含 exact interning、active 索引和公共分割建立，10 万级约为旧路径
的 46%。五案例端到端时间单独记录在 comparison 中，不据此宣称 H4 的大量候选
评估、内存开销或整体生成速度已解决；本阶段也没有重做 H3 的 50 万/百万级验证。
隔离测速在 CTest 完成后执行，原始三次数据和平台信息见
`artifacts/q2a/topology-benchmark.json`。

五案例均通过内部完整拓扑、独立 OpenFOAM owner/neighbour/正体积读取器和真实
`opencfd/openfoam-run:2606 checkMesh -writeAllFields`，结果均为 **Mesh OK**。
镜像 ID 为 `sha256:4229997e74defb81548222d511b8e3b95b98305e5df41b8e88b031813fe47eeb`；
逐例命令、运行时间和日志 SHA-256 见 `artifacts/q2a/openfoam.json`。
`artifacts/q2a/comparison.json` 的 `q2a_status` 为 `PASS`，但
`q2_full_status` 仍为 `PARTIAL_NOT_ACCEPTED`。

```bash
cmake -S . -B build-q2a -DCMAKE_BUILD_TYPE=Release
cmake --build build-q2a -j4
ctest --test-dir build-q2a --output-on-failure -j4
python3 tools/verification/verify_q2a_construction.py
python3 tools/verification/check_q2a_openfoam.py
python3 tools/verification/verify_q2a_construction.py --collect-only
python3 tools/verification/benchmark_q2a_topology.py
python3 tools/visualization/render_q2a_meshes.py
```

默认对照输入是已有 `build-q2/after/`。若从干净 checkout 复现，先独立构建基线：

```bash
mkdir -p build-q2a/base-source
git archive dde1c26eceede55364dca53d3c9ecac2d700b825 | tar -x -C build-q2a/base-source
cmake -S build-q2a/base-source -B build-q2a/base-build -DCMAKE_BUILD_TYPE=Release
cmake --build build-q2a/base-build -j4
python3 tools/verification/generate_q1_baselines.py --build-dir build-q2a/base-build --evidence-dir build-q2a/base-evidence --output-dir build-q2a/base-summary --source-commit dde1c26eceede55364dca53d3c9ecac2d700b825 --expect-superellipse-short-faces absent
python3 tools/verification/verify_q2a_construction.py --baseline build-q2a/base-evidence
```

后续 collect-only 也应指定相同 `--baseline`。OpenFOAM 脚本使用已有镜像
`opencfd/openfoam-run:2606` 的默认 `checkMesh -writeAllFields`，并非新增
`-allTopology -allGeometry` 验收声明。
