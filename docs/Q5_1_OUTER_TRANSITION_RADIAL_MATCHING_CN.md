# Q5-1 最外层 transition 径向尺寸匹配

日期：2026-08-31

## 结论

RemainderCut 不是一个 RemainderCut 构造缺陷。circle 与 superellipse 的**全部**
Q1 hard face 都是 RemainderCut↔Transition 面，小的一侧是合法的 Cut 碎片，大的
一侧是过大的 transition 单元。

按"最外层 transition 行的径向厚度不得超过它所邻接的 remainder background cell"
这条规则把最外环切成径向多行后：

| case | 状态 | typed hard | solver cells |
|---|---|---:|---:|
| circle | FAIL → **WARN** | 80 → **0** | 728 → 856 |
| superellipse | FAIL → FAIL | 26 → 12 | 795 → 987 |
| concave-L / sharp-tail / narrow-gap | 不变 | 不变 | 不变 |

circle 是本项目第一个 Q1 hard 数为 0 的几何。

机器可读证据：`artifacts/q5-1/manifest.json`。

## 1. 诊断

按相邻单元种类拆分 hard face：

| case | 面对 | volume ratio | face weight |
|---|---|---:|---:|
| circle | RemainderCut ↔ Transition | 56 | 24 |
| superellipse | RemainderCut ↔ Transition | 6 | 12 |
| superellipse | RemainderCut ↔ RemainderCut | 4 | 4 |

关键尺寸（以 remainder background cell 面积 h² 为单位）：

| case | 最外层 transition 单元 | RemainderCut 单元 | 最差 parent volume ratio |
|---|---:|---:|---:|
| circle | 3.93 h² | min 0.155 h²，median 0.749 h² | 0.039382 |
| superellipse | 2.48…25.5 h² | min 0.211 h²，median 0.689 h² | 0.072707 |

graded fan 沿外法向逐环加密**切向**，但每一环都保持整个 ring thickness。circle 的
ring thickness 是 0.1125，而 h=0.0625，即 1.8 个 background cell 深，于是最外一行
的面积约为四个 background cell，而它邻接的 Cut 碎片只有六分之一个。这两者之间的
面无论 remainder 怎么划分都到不了 0.05 的无量纲 volume-ratio 硬门槛。

## 2. 已排除的做法（实测，避免重复提案）

- **Q4-1 式严格凸 union**：circle 40 个、superellipse 14 个携带 hard face 的
  RemainderCut 单元中，**0 个**能与任何 mutable remainder 邻居构成面积精确的严格凸
  union。原因是 transition front 是凸壁面的多边形外扩，从 remainder 一侧看是凹的，
  贴着 front 的 Cut 单元合并后必然保留反折角。
- **simple union 的二片凸重划分**：只解决 88 个中的 8 个（circle）、65 个中的 18 个
  （superellipse）。union 面积够了，但非凸，ear-clipping 又把它切回一个 sliver。
- **邻居感知的 ear-clip 对角线选择**：circle 24 个、superellipse 10 个反折 parent
  中，**0 个**存在能降低 hard 数的替代对角线；现有"面积最均衡"选择已是最优。因为
  parent 本身相对 transition 邻居已不足 2 倍门槛，任何二分的两片都不合格。
- **全局 transition sizing 重调**：`minimumTotalWidthCells` 取 5.4/4.5/4.0/3.6/3.2 时
  circle 的 hard 总数为 80/56/8/56/36，superellipse 为 26/98/12/84 —— 非单调。
  `minimumRingCount` 取 4 能清空 circle，但 sharp-tail 掉到 pure Cut-cell fallback，
  concave-L 从 206 升到 241。单个全局常数是 quadtree 相位抽奖，不是规则。

## 3. 实现的规则

**任何 transition 行的径向厚度不得超过它邻接的 remainder background cell。**

最外环切成 `ceil(ringThickness / targetCellSize)` 个径向行，上限 4。两个量都已存在于
transition plan 中（`ringThickness = totalThickness/ringCount`，
`targetCellSize = domainScale · 2^-boundaryLevel`），**没有引入新的调参常数**。
circle 与 superellipse 解析出的行数均为 2。

安全性质：最后一行的径向偏移写作 `ring + rows/rows`，与单行时的 `ring + 1`
在 IEEE 下逐位相同，因此**提交的外包络不变**——remainder quadtree、它的 Cut cell 和
不可变的 H4-1 层都不受影响，只在 transition strip 内部增加面。实测验证：
circle 的 `remainder_cut_cell_count` 188→188、`remainder_cartesian_cell_count`
304→304、`boundary_layer_cell_count` 128→128，superellipse 同理 196/332 不变。

opt-in flag：`--q5-outer-transition-radial`。新增上报位
`q51_outer_transition_radial_subdivision` 与 `q51_outer_transition_radial_declined`；
后者复用 Q4-1 的 decline 阶梯，保证可选构造改进永远不会让产品掉到 pure Cut-cell。

## 4. 验收

- area error `-1.06581e-14`；interface、owner/neighbour、stable ID、source lineage
  invariants 全部 PASS；lineage mismatch=0。
- deterministic replay 的 hybrid CM2D、solver CM2D、quality-contract JSON、
  solver-quality JSON、VTK 逐字节一致。
- CTest **75/75 PASS**。
- 独立 OpenFOAM reader：circle PASS（856 cells、2312 points）。
- `opencfd/openfoam-run:2606` `checkMesh`：两例均 `Mesh OK`。circle max aspect
  5.06241311、max non-orth 56.98720338、max skewness 1.328876647；superellipse
  31.88659624 / 65.0532228 / 2.119583125。
- 回归门：`tests/hybrid_mesh_test.cpp` 断言 circle 的 Cut/Cartesian/layer 计数不变、
  transition 行增加、接口共形、面积守恒、typed hard 数严格下降；CI 增加
  `q51_circle` 与 `q51_superellipse` 两个 case 及其 checkMesh。

## 5. 剩余 —— superellipse 那 12 个也已实测到底

Q5-1 之后 superellipse 剩 12 个 hard，构成为：RemainderCut↔RemainderCut 的
face_weight 4 + volume_ratio 4，以及 RemainderCut↔Transition 的
non-orthogonality 4。两族都追到了根：

**RemainderCut↔RemainderCut（8 个面，4 个对称单元 89/315/408/499）。** 该单元
面积 1.648384e-3 = 0.2153 h²，是一个薄楔形，带一个**真实的 −16.943389° 凹角**
（另有一个 −0.000002° 的近共线凹角）。真实凹角必须切分（OpenFOAM 会把它报成
concave polyhedron），而薄楔形的任何凸切分都会产生 sliver：committed 子单元为
1.43889e-3 与 2.0949e-4，最小子/母 = 0.12709，最差面比 0.045461，距 0.05 门槛
差 9%。逐一实测三条出路，全部无效：

| 尝试 | 最好能达到的最差面比 |
|---|---:|
| 现状（面积最均衡对角线） | 0.045461 |
| 穷举全部 ≤3 片凸切分 | 0.045461 |
| 先删近共线凹点再重划分 | 0.045461 |
| 与 RemainderCut 邻居合并后切分 | union 根本无合法凸切分 |
| 与 RemainderCartesian 邻居合并后切分 | 0.013579（**更差**）|

所以它同样**不是** remainder 构造缺陷：合并会让结果更糟，而不是更好。真正的出路只有
两条，且都已被本项目记录为受阻——放宽严格凸性（`strictlyConvex` 的注释说明这是为了
OpenFOAM 的 concavity 检查而故意收紧的），或局部 quadtree 加密改变切割本身（Q2 文档
记载该路线"failed the unchanged solver gate"，Q2 因此仍为 partial）。

**non-orthogonality（4 个面）。** 实测 65.0532228°，门槛 65°，仅超 0.08%，且出现在
新细化后的 fan 内部——同一批面在 Q5-1 之前是低于门槛的。这是本轮唯一的负向代价，
净 hard 数仍从 26 减到 12。

## 6. 下一步

1. concave-L 206、sharp-tail 377、narrow-gap 430。这三例的 transition 来自 graded
   termination buffer 而不是 fan，Q5-1 对它们完全 inert。把同一条径向尺寸匹配规则
   搬到 termination buffer 是最自然的下一步，本轮没有实现。
2. BoundaryLayer 六项指标仍为 OBSERVED，因此 WARN 还不能认定近壁网格质量合格。

