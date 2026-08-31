# Q3-1 termination 局部质量优化

日期：2026-08-31

## 范围与结论

本轮基线为 `origin/gpt/r1-closeout-20260830`，commit
`fe34b111460c4aef3c1825f687118975bc4ed7fd`。只处理 narrow-gap 的 termination
局部 volume ratio 与 face weight；Q2 generic recovery、R1 扩展、sharp-tail、壁面
几何、pure Cut-cell fallback 和质量阈值均未进入本轮。

Q3-1 的 32 次 bounded transaction 将 volume-ratio hard count 从 295 降到 265，
face-weight hard count 从 95 降到 65；hard short-face 保持 0。Q1 整体仍为 FAIL，
本轮不把局部改善冒充全质量收口。

机器可读证据：`artifacts/q3-1/manifest.json`。网格预览：
`artifacts/q3-1/narrow-gap-after.svg.png`。

## 1. hard violation 来源定位

R1 baseline 的 hard issue 分布如下：

| metric | total | termination | cartesian | worst | max severity |
|---|---:|---:|---:|---:|---:|
| volume ratio | 295 | 260 | 35 | 0.01185326341 | 4.218247606 |
| face weight | 95 | 54 | 41 | 0.05157451029 | 1.938942308 |

最差 volume ratio 在四个对称 termination 位置重复出现。违规集中在 termination
小片与邻接较大 cell 的离散面积跳变，以及由模板分片产生的 cell-centroid / face-centre
位置关系；不是一个全局 sizing 参数偏小可以解释的单一模式。因此本轮判断为
**topology/template 主导，sizing 次要**。

可复现 baseline 命令：

```bash
build-q3/cartmesh2d_hybrid_cli examples/h4_3/narrow_gap.xy \
  build-q3/evidence/q3-1-baseline/narrow_gap 8 3 8 4 0.012 1.15 1.0 \
  build-q3/evidence/q3-1-baseline/narrow_gap-case 0.01
```

## 2. 最小 bounded candidate generator

Q3-1 是显式 opt-in：`--q3-termination-quality`。每次 transaction 只做：

1. 从 termination-adjacent hard faces 中按 hard count、volume severity、face-weight
   severity 和 stable cell ID 选出唯一最差 face；
2. 只枚举该 face 两端及其一环中的 mutable/mutable 相邻对，单次最多 16 个候选；
3. 候选只允许把两个相邻 polygon 合成一个严格凸、面积相同且非 under-determined
   的 polygon；不移动点、不改 wall、不删坏 cell、不创建 fallback cell；
4. 使用现有 stable vertex ID、boundary lock、`TopologyDelta2D`、patch-local scope 和
   patch-local authoritative metric kernels；
5. 候选循环不执行 global topology build 或 full solver-quality；唯一 winner 才执行
   一次 global oracle 与一次 authoritative full quality；
6. 一次 build 最多提交 32 个 winner。达到上限表示 bounded partial pass，不表示全局
   收敛。

局部 gate 要求 volume-ratio 与 face-weight hard count 各自不增加、合计严格下降，
hard short-face 不增加，solver hard issue 不增加，Q1 angle / non-orth / skewness /
aspect hard count 不增加。候选总序随后比较 volume severity、face-weight severity、
short-face、angle、non-orth、skewness、aspect，最后用 stable cell ID 打破平局。

## 3. before / after

| metric | before | after | change |
|---|---:|---:|---:|
| solver cells | 3185 | 3153 | -32 |
| volume-ratio hard | 295 | 265 | -30 |
| face-weight hard | 95 | 65 | -30 |
| hard short-face | 0 | 0 | unchanged |
| min volume ratio | 0.01185326341 | 0.01712166482 | improved |
| min face weight | 0.05157451029 | 0.05740815582 | improved |
| max volume severity | 4.218247606 | 2.920276768 | improved |
| max face-weight severity | 1.938942308 | 1.741912775 | improved |
| hard minimum angle | 20 | 12 | improved |
| hard non-orthogonality | 14 | 12 | improved |
| hard skewness | 4 | 4 | unchanged |
| hard hydraulic aspect | 2 | 2 | unchanged |
| area error | -3.37508e-14 | -3.37508e-14 | unchanged |

执行数据：256 generated candidates、87 valid local candidates、174 local-quality
evaluations、32 commits；candidate-loop global builds = 0，candidate-loop full-quality
evaluations = 0，winner global oracle builds = 32，maximum per transaction = 1，
authoritative full-quality evaluations = 64。本机 Q3 repair wall time 22.765 s，完整命令
wall time 68.60 s；wall time 属环境测量，不作为跨机器确定性数据。

## 4. 验收

- area、construction topology、owner/neighbour、boundary patch 与接口闭合不回归；
- patch 外 stable IDs unchanged，local delta matches global oracle，local winner matches
  global authority；
- 同机重放的 hybrid CM2D、solver CM2D、quality-contract JSON 等八类核心文件逐字节
  相同；
- independent hybrid reader PASS；independent OpenFOAM reader PASS；
- CTest 75/75 PASS；
- `opencfd/openfoam-run:2606` 的真实 `checkMesh -writeAllFields`：3153 cells、13403
  faces、6801 internal faces、max non-orthogonality 68.81332468、max skewness
  3.442394378，最终 `Mesh OK.`。

因此 Q3-1 narrow-gap 的“机制有效”验收满足，但不是全部 Q1 hard issue 清零。

## 5. Q3-2 建议

先继续 narrow-gap，不碰 sharp-tail。Q3-2 应把当前 32-transaction 上限后的剩余
hard faces 按模板家族分组，补一个仍受一环 boundary lock 限制的二选一 repartition
模板，重点处理单纯 pair agglomeration 无法继续降低的 termination/cartesian 面。
只有当 narrow-gap 的剩余 hard count 与 transaction-bound 行为稳定后，才为
sharp-tail 建立独立 baseline 和候选域证明。
