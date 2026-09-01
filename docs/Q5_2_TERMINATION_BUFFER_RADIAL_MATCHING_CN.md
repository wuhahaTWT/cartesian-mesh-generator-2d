# Q5-2 graded termination buffer 径向尺寸匹配

日期：2026-09-01

## 结论

把 Q5-1 那条规则（任何 transition 行不得比它邻接的 remainder background cell 更厚）
搬到 graded local-termination buffer。但 buffer 与 fan 有一个本质差别：**fan 是显式
逐环插值，buffer 是一次 locally reduced march**。保持"请求的总厚度"并不能保持
**实际的阶梯前沿**，因为每个被局部削减的 column 停在了不同距离上。所以 remainder
会变，结果只能实测、不能假设。

因此 Q5-2 做成**自门控**：两套前沿都构造出来，只有当重解出的 buffer 让 typed hard
count **严格下降**时才提交；平手保留历史 march。三例实测：

| case | baseline hard | matched hard | 提交 | 结果 |
|---|---:|---:|---|---:|
| sharp trailing edge | 377 | **238** | 是 | **377 → 238（−37%）** |
| concave-L | 206 | 255 | 否 | 206（逐字节不变）|
| narrow gap | 430 | 749 | 否 | 430（逐字节不变）|

这就是为什么必须自门控：同一条规则在 sharp-tail 上砍掉 37%，在 narrow-gap 上会让
hard 数涨 74%。没有门控就是一次尺寸抽奖。

机器可读证据：`artifacts/q5-2/manifest.json`。

## 1. 先定位：hard face 到底长在哪

按两侧单元种类拆开每一个 hard face，并用"从 H4-1 层出发的图距离"把 Termination
单元区分为 buffer 第 1/2/3 行与 remainder 派生：

| case | 最大族 | 计数 | 占比 |
|---|---|---:|---:|
| narrow gap | `Term(rem) ← Term(buf3)` VR 160 + FW 16 | 176 / 430 | 41% |
| sharp trailing edge | `Term(rem) ← Term(buf3)` FW 100 + VR 76 | 176 / 377 | 47% |
| concave-L | `Term(buf1) ← BL` VR 31 + FW 13 | 44 / 206 | 21% |

三例的第二/第三大族分别是 `Term(buf3) ← Term(buf2)`（narrow-gap 44、sharp-tail 20）
和 `Term(buf1) ← BL`。**最外层 buffer 行确实是主要肇事者**，与 Q5-1 在 fan 上的
发现同构。

实测 buffer 各行面积（以 remainder background cell 面积 h² 为单位）：

| case | h | row 1 | row 2 | row 3 |
|---|---:|---:|---:|---:|
| concave-L | 0.01953125 | 1.495 h² | 2.168 h² | **3.143 h²** |
| sharp trailing edge | 0.01657282 | 1.69 h² | 2.451 h² | **3.554 h²** |
| narrow gap | 0.01366629 | 0.891 h² | 4.733 h² | **7.336 h²** |

最外行是 3.1–7.3 个背景格，而它邻接的 remainder Cut 碎片只是一个背景格的一小部分。

## 2. 为什么不能照抄 Q5-1 的做法

Q5-1 安全的关键是外包络**逐位不变**：最后一行的偏移写成 `ring + rows/rows`，与单行
的 `ring + 1` 完全相同。buffer 没有这个结构：

- buffer 是 `buildLocallyReducedBoundaryLayerStrips2D` 的产物，前沿是**阶梯状**的，
  每个 column 的层数由局部碰撞/凹角判据决定。
- 行数从 3 改到 4 后，一个原本只能走 2 行（到 0.0447）的 column 现在可能走 3 行
  （到 0.0554）。阶梯位置整体改变，remainder 随之改变。
- 事后再切分也不行：在阶梯处，最外层 cell 的**径向面**本身就在外包络上，径向细分
  会往包络里插入共线顶点，而 `strictlyConvex` 按 `cross <= 1e-10*scale` 把精确共线
  判为凹，于是每个这样的 Cut cell 都被强制 ear-clip，正是 Q5-1 收口时量到的
  sliver 生成机制。localTermination 下 `canonicalizeTransitionEnvelope2D` 直接
  early-return，也没有东西会去清理这些顶点。

所以 Q5-2 选择"重解 march"而不是"事后切分"，并接受 remainder 会变，用门控兜住。

## 3. 规则

总厚度 `T` 取未加规则时 march 会给出的值，因此外前沿的**目标**位置不变。然后：

1. `rows = ` 最小的行数使 `T / rows <= cap`，其中 `cap = 1 × remainderCellSize`；
2. 在 `(1, 原 growthRatio]` 上取**最大**的 `g` 使最外行 `<= cap`。固定 `T` 时更大的
   `g` 让第一行更薄，因此这样能保住尽可能多的近壁分级。

sharp-tail 实测：`rows 3 → 4`、`g 1.45 → 1.08698`、最外行
`0.0383717 → 0.0234375`，正好等于 `cap = 0.0234375`。总厚度不变。

没有新增可调常数：`cap` 沿用 Q5-1 的 `outerTransitionRadialTargetCells = 1.0`，
行数上限 `maximumTerminationBufferRows = 12`。

**有界搜索的可行性必须显式检查。** 行数上限可能在 cap 变得可达之前就用完。均匀
march（`ratio = 1`）是给定行数下最薄的最外行，因此若 `total / rows > cap` 仍然成立，
这条规则在此几何上就是不可满足的：此时**保留历史 march**，并把
`q52_termination_buffer_row_cap_reachable` 报为 false，而不是把一个仍然超标的行数
当成"已匹配"提交。这同时也是二分法成立的前提——只有在 `total / rows <= cap` 之后，
`ratio = 1` 才真的是可行下界，否则 `f(low) > cap`，二分会收敛到一个仍然违反 cap 的
比率。Q5-1 的 `maximumOuterTransitionRadialSubdivision` 有完全相同的结构，同样加了
循环后的可行性检查；不可达时回到单行 fan 并报
`q51_outer_transition_radial_target_reachable = false`。

这两处是 codex 独立审查发现的 P1，实测确认：把 `outerTransitionRadialTargetCells`
收紧到 0.01 时，修复前会提交一个最外行远超 cap 的网格（仅因 hard 数下降就通过），
修复后 `committed = false`、`cap_reachable = false`，网格保持为历史 march 的 Hybrid。
默认参数落在可行区内，因此五个 Q5 案例的 solver `.cm2d` 在修复前后逐字节相同。

## 4. sharp-tail 收益构成

| metric | baseline | Q5-2 | 变化 |
|---|---:|---:|---|
| volume ratio | 156 | **78** | 减半 |
| face weight | 145 | **78** | −46% |
| face/local h | 14 | 14 | 不变 |
| face/sqrt(owner area) | 10 | 10 | 不变 |
| face/sqrt(neighbour area) | 14 | 14 | 不变 |
| minimum interior angle | 23 | 25 | +2 |
| non-orthogonality | 11 | 12 | +1 |
| skewness | 2 | 4 | +2 |
| hydraulic aspect | 2 | 3 | +1 |
| **合计** | **377** | **238** | **−37%** |

solver cells 3391 → 3629，boundary-layer cells 1044 **不变**，area error
`4.619e-14`。两个主要族各减半，代价是 6 个次要 hard，净减 139。

## 5. 验收

- CTest 75/75 PASS。
- 被拒绝的两例（concave-L、narrow-gap）的 `hybrid.cm2d`、`hybrid.solver.cm2d`、
  `hybrid.quality-contract.json` 与各自 baseline **逐字节相同**——门控不是"跑了但
  没提交"，而是真正回到历史构造。
- circle 与 superellipse 的 `localReductionApplied` 为 false，Q5-2 完全 inert，
  与仅开 Q5-1 的 solver `.cm2d` 逐字节相同（circle 仍为 WARN / hard 0，
  superellipse 仍为 12）。
- deterministic replay：sharp-tail 的四类文件重跑逐字节一致。
- 独立 OpenFOAM reader PASS。
- OpenFOAM v2606 `checkMesh` 对 Q5-2 sharp-tail 为 `Mesh OK`：max aspect ratio
  **8.32883709**（Q4-1 那版是 13.912）、max non-orthogonality 69.3203693、
  max skewness 3.34005325、minimum face area 9.849335481e-08。
- 回归门：单元测试断言 sharp taper 在 Q5-2 policy 下仍为 Hybrid、`committed` 恰好
  等于 `matched < historical`、且该例确实是一次实测胜利；CI 对三例断言
  `committed == (matched < historical)`、`historical == baseline`、
  `after <= before`、被拒绝时逐字节回到 baseline，并加上 `checkMesh` 覆盖。

## 6. 成本

开启该 flag 时每个 localTermination 几何构造两套前沿，即两次完整 hybrid build。
这是全局候选比较而非 patch-local 选择，所以代价是一次额外 build，不是常数级。
关掉 flag 时零成本。

## 7. 剩余

1. `Term(buf1) ← BL` 族（concave-L 44、sharp-tail 38、narrow-gap 34）没有被触及：
   那是 buffer **第一行**与 H4-1 最后一层之间的尺寸落差，方向与本轮相反（buffer 行
   比层单元大）。Q5-2 提高 `g` 已经让第一行变薄，但没有针对它的判据。
2. concave-L 与 narrow-gap 仍是 FAIL 且被 Q5-2 拒绝。它们的 `Term(rem) ← Term(buf3)`
   族依然存在（narrow-gap 176 个面），但重解 march 会连带恶化更多别的东西。要拿到
   这部分，需要一种**能保持阶梯前沿**的 buffer 细分，而这又要求先解决共线顶点被
   `strictlyConvex` 判为凹的问题——即 Q2 文档记录为受阻的那条路线。
3. BoundaryLayer 六项指标仍为 OBSERVED。
