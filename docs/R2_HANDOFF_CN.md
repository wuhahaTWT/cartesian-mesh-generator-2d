# R2 交接指导：加密鲁棒性与 CFD 可解性

给接手本项目的下一个 agent。原始日期：2026-09-01；W1 状态更新：2026-09-04。

## 0. 一句话现状

W1 已关闭：**纯 Cut-cell 路径已实测稳到 level 11**（180468 stabilized cells、
277336 leaves、全部 solver-quality 硬限满足），原 level 10 格点 spur 已消除；
**hybrid 边界层路径在 level 9 就失败**，根因是壁面切向分辨率被输入折线顶点数锁死。
两条线原因完全不同、代码不重叠，用户要求**并行推进**。

## 1. 开工前必读（按顺序，别跳）

1. `AGENTS.md` — 硬规则。特别是 §4.10「不得降低阈值、删除坏单元或隐藏告警」和
   §2「不得引入三维 `cartmesh/*` 头文件或链接三维 library」。
2. `docs/R2_REFINEMENT_ROBUSTNESS_CN.md` — **本轮全部实测数据与根因**。第 7 节列出
   已经被实测排除的三条路线，不要重走。
3. `docs/R2_REFERENCE_AUDIT_CN.md` — 参考项目许可证与借鉴等级。OpenFOAM / cfMesh /
   gmsh / p4est 全是 GPL，**只读不抄**；AMReX 是 BSD-3。用户已确认全程 Level C/D。
4. `docs/R1F_PATCH_LOCAL_CLOSEOUT_CN.md` — 上一轮收口结论与两条硬前置条件。
5. `artifacts/r2/w0-baseline-manifest.json` — 机器可读基线，你的每次改动都要和它对比。

## 2. 环境（会浪费你一个构建周期的坑）

```bash
export DYLD_LIBRARY_PATH=/Applications/mesasdk/lib   # ctest 也需要，不只是直接跑 CLI
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCARTMESH2D_BUILD_TESTS=ON
cmake --build build -j8
ctest --test-dir build --output-on-failure           # 当前基线 83/83
```

macOS SIP 会把 `DYLD_*` 从 Apple 签名的二进制里剥掉，所以**不要用 `/usr/bin/time`
包装 CLI**，否则它会以 dyld 错误退出，看起来像构建坏了。用 shell 内建计时。

`timeout(1)` 在 macOS 上不存在。

## 3. 唯一的测量入口

```bash
python3 tools/verification/refinement_ladder.py \
  --repo . --build-dir build --evidence-dir build/r2_ladder \
  --output-dir artifacts/r2 --manifest-name <your>-manifest.json \
  --dyld-library-path /Applications/mesasdk/lib --stop-on-failure \
  --ladder circle:cutcell:6,7,8,9,10 \
  --ladder circle:hybrid:6,7,8,9 \
  --ladder narrow_gap:hybrid:8
```

`--ladder CASE:MODE:LEVELS`，`MODE` = `hybrid` | `cutcell`。加 `--gate` 让它对
非单调阶梯返回非零。门判的是**单调性**：任何一级失败、`topology_valid` /
`mesh_quality_valid` / `solver_quality_valid` 为假、跌破 `min faceWeight 0.05` /
`min volRatio 0.01` / 超过 `nonOrtho 70°`、或 Q1 hard 总数随 level 增长，都是违规。

**改动前后各跑一次，把两份 manifest 一起放进证据。** 不要只跑一个 level。

## 4. 主线 A：纯 Cut-cell → level 10/11（W1 已完成，以下保留历史）

### A1（当前阻塞点，先做这个）把 `cartmesh2d_cli` 迁到共享构造 registry

**为什么必须迁**：level 10 的 spur 只能靠焦合消除，而焦合必须在**全局 canonical
vertex 空间**里做。逐 face、逐 leaf 都已实测失败（见 §7）。`cartmesh2d_cli` 目前走的
是 legacy 路径：[apps/cartmesh2d_cli.cpp:612](apps/cartmesh2d_cli.cpp:612)
`buildCutCells(leaf, boundary, fluidRegion)`，完全没有接入 Q2-A 的 registry。

**机制已经就绪**（本轮已提交）：
[IntersectionRegistryPolicy2D::gridCornerWeldFractionOfLocalH](include/cartmesh2d/geometry/IntersectionRegistry2D.hpp:94)。
默认值等于原来硬编码的 `64·DBL_EPSILON`，所以现有产物逐字节不变；需要几何焦合的
调用方显式提高它。构造函数拒绝 `>= 0.01`（Q1 短面硬限）。

**卡住的具体错误**（稳定复现）：

```text
terminate called after throwing an instance of 'std::invalid_argument'
  what():  grid line lies outside construction support
```

抛出点 [IntersectionConstruction2D.cpp:301](src/geometry/IntersectionConstruction2D.cpp:301)：
`t=(target-origin)/delta` 要求 `t ∈ [0,1]`，即 grid line 必须真的穿过被注册的
support 线段。hybrid 路径用**同一份代码**工作正常，所以差异在 setup：

| | hybrid（能跑） | 纯路径（抛异常） |
|---|---|---|
| 裁剪对象 | transition envelope（`remainderBoundaryRegion`） | 原始 wall（`boundary`） |
| `configureGrid` | `(domain.bounds, remainderMaxLevel)` | `(domain.bounds, maxLevel)` |
| intern 的顶点 | strip wall 顶点，带 Sharp/Concave/Smooth 分类 | 全部 loop 顶点，一律 Smooth |

**下一步不要在 CLI 层继续试错**。正确做法：写一个最小 C++ 复现，直接构造一个
`IntersectionRegistry2D` + 一条 wall segment + 一个 leaf box，找出哪个
`(support, gridLine)` 组合让 `t` 落到 `[0,1]` 之外。怀疑方向是
`clipSegmentToAABB` 里 `t0/t1` 被 `std::clamp` 之后仍然把 `enter`/`leave` 记录的
box side 传给 `resolve()`——也就是 segment 端点在 box 内部时，`enter` 可能来自一条
它并不真正到达的边。参见 [src/cutcell/CutCell2D.cpp:195](src/cutcell/CutCell2D.cpp:195)
附近的 `resolve` lambda。

尝试稿（未提交，仅供参考，别直接套用）：`/tmp/w1_cli_attempt.cpp`。
它已经写好了 registry 创建、顶点 intern、recovery request 的 fail-closed 处理，以及
把 registry 传给 `buildGlobalTopology` 的第五个参数。

### A2 用已推导的界放宽面积不变量（用户已批准，不要偷偷改数）

焦合会真实扰动面积，这一点无法回避。实测（circle level 10）：

| 焦合预算（h 的比例） | 面积误差 | 是否消除 spur |
|---|---:|---|
| `1e-3` | `1.17e-8` | 是 |
| `1e-4` | `1.79e-9` | 是 |
| `1e-5` | — | **否** |

spur 两点间距是 `4e-5 · h`，这是预算下界。而
[apps/cartmesh2d_cli.cpp](apps/cartmesh2d_cli.cpp) 的物理门只允许
`1e-10 * max(1, area)`（注意它**没有** hybrid 路径的 256 倍
`areaToleranceMultiplier`），circle 上是 `5.88e-10`。

放宽必须用**推导的界**，不是调出来的数字：

```text
单点位移 <= f · h
一个 cell 面积变化 <= 1.5 · f · h²
焦合次数 <= 壁面长度 / h
总界 <= 1.5 · f · h · wallLength     （与 level 无关）
```

用 `f = 1e-4`、circle L10 的 `h = 2.93e-3`、`wallLength ≈ 6.24` 代入得 `2.6e-9`，
实测 `1.79e-9`，界是对的且约 1.5 倍保守。**把这段推导写进代码注释和文档**，
让放宽是公开可审计的结论，而不是一次静默的门下调。

### A3 level 11 及以后

`level 11` 已通过；关键数字与改前/改后 manifest 见
`docs/R2_REFINEMENT_ROBUSTNESS_CN.md` 第 7.6 节。

## 5. 主线 B：hybrid 混合笛卡尔 → level 10

这条线**没有阻塞点**，根因已经完全定位，可以直接开工。它的价值比主线 A 更高：
hybrid 在更低的 level 就挂，而且它挡住的是整个边界层（也就是所有壁面解析的粘性 CFD）。

### B1（最高优先）wall chain 按局部背景 h 切向再分

**根因**：[src/hybrid/HybridMesh2D.cpp:631](src/hybrid/HybridMesh2D.cpp:631)
的 `makeClosedWallChain2D` 直接把输入 `BoundaryLoop` 的顶点当 wall column。
`circle.xy` 只有 32 个顶点 → 32 column × 4 层 = **128 个 BL cell，与 level 无关**
（实测 level 6/7/8 全是 128）。于是切向尺寸恒为 `2πr/32 ≈ 0.196`，而 level 9 的背景
`h ≈ 0.0113`——背景网格比壁面细 **17 倍**。

[resolveAutomaticHybridTransitionPlan2D](include/cartmesh2d/hybrid/HybridMesh2D.hpp:437)
必须独自吸收整个落差：它按 `maxOuterEdgeLength / finalSubdivision <= 2h` 把切向细分
翻倍到 16、ringCount 到 5，而 `ringThickness` 仍绑在最后一层法向间距上，结果是**高度
各向异性、强斜的过渡多边形**——正好对应实测的 nonOrtho 77–79°、faceWeight 0.024–0.032。

**已排除的旁路**：把首层厚度随 level 折半**没有救回来**，只是把失败从 132 s 拖到
803 s。问题在切向，不在法向。

**做法**：在 `makeClosedWallChain2D` 之前插入一个确定性重采样，把每条输入 segment
二分到 `segmentLength <= 2·h_local`。**新点必须是原 segment 上的线性插值点**，因此到
输入几何的距离保持数值零——这条原则不能让步（`AGENTS.md` 禁止用平滑/移动顶点换质量）。
建议落在 `src/boundary_layer/WallChainRefinement2D.{hpp,cpp}`，不改 `WallChain2D`
的结构，只改喂给它的顶点序列。

**预期**：`maxOuterEdgeLength ≈ 2h`，过渡计划退回 `finalSubdivision = 1`、
`ringCount = 3`（最小值），各向异性消失，Q1 hard 计数不再随 level 增长。

**同一根因的三维先例**（只作设计参考，`AGENTS.md` §2 禁止复用三维代码）：
`cartesian-mesh-generator/STAGE6_8_CFD_QUALITY_PLAN_CN.md` 记录 full cube 876 个、
icosahedron 744 个 non-orthogonality **全部位于 Cut ↔ Layer 接口**，采用的解法就是
fragment-matched layer column（沿原始 STL 三角面重心插值细分，几何不移动）。
同一份文档**明确否决**了「把 layer facet 外侧 Cut-cell 粗暴聚并成一个大 transition
cell」（nonOrtho 恶化到 111°/129°）——**不要走这条路**。

### B2 层厚可相对化 + 窄缝 medial-axis 上限

[LayerParameters2D](include/cartmesh2d/boundary_layer/BoundaryLayer2D.hpp:70)
目前只有绝对厚度。增加「相对局部 h」模式与 `minThickness`，语义对齐 snappyHexMesh 的
`layerParameters.H`（`relativeSizes_` / `minThickness_` / `expansionRatio_`）。
**默认保持绝对模式**，CLI 加显式开关，旧命令行为逐字节不变。

窄缝上限对齐 snappy 的 `maxThicknessToMedialRatio`
（`medialAxisMeshMover.C:1378`，起点 0.3）：为每个 wall 采样点求「到最近**非相邻**
wall segment 的距离」，把总层厚限制为该距离的固定比例。`narrow_gap.xy` 是两个 1×1
方块、缝宽 0.08 → 半宽 0.04 → 总厚上限 0.012，正好只容一层；当前参数
（4 层、首层 0.012、增长 1.15）总厚 0.0599，两侧对撞 0.12 > 0.08，这就是 narrow_gap
一直要靠 Q3/Q4/R1 事后修的原因。实现方式：给
[BoundarySegmentIndex2D](include/cartmesh2d/spatial/BoundarySegmentIndex2D.hpp:31)
补一个最近点/最近 segment 查询（它已有 flat BVH 与 `distanceToAABB`），再按 chain
邻接过滤，**不要新建空间结构**。

### B3 性能会顺带解决，不要单独优化

实测 circle L9 的 136 s = 22.9 s 原始层尝试 + 113.6 s 局部降层尝试，只有 **2 次**
conformal build，未归因 1e-6 s。也就是说时间花在**单次构造内部的修复回路**在坏几何上
反复空转，不是 growth-ratio 重试。B1 把坏几何从源头消掉，修复回路自然无事可做。

R1F §5 已经指出下一条质量主线应该是 construction-time quality，**不要继续扩大
Q3/Q4/R1 的 repair framework**。

## 6. 两条线合流后的第三、第四步

### W3 质量驱动的小单元合并（等 A、B 的阶梯都不崩再做）

**根因**：[SmallCellPolicy2D::areaFractionThreshold = 0.10](include/cartmesh2d/stabilization/SmallCell2D.hpp:51)
是相对**自身背景格**的面积分数，而 volume ratio 是相对**邻居**的。一个 alpha=0.107 的
cut cell 紧邻粗一级的邻居时，volume ratio 就是 `0.107/4 ≈ 0.027`——正好是实测
circle L6 的 `min volRatio = 0.0270`。所以 Q1 hard 计数随加密单调增长
（80 → 132 → 360）是**阈值语义错配的必然结果**，不是偶发。
`STAGE2DV1C_VERIFICATION.md` 里「level 8 必须把 small-alpha 从 0.05 手调到 0.20 才过门」
是同一件事的历史证据。

**做法**：判据改成直接看目标指标（volume ratio / face weight）合并前后是否越过硬限，
迭代到达标或无候选。判据**必须复用**
`evaluateSolverInternalFaceMetrics2D`（`include/cartmesh2d/quality/SolverQuality2D.hpp`）
这一唯一权威公式，不得另写一份局部公式（R1 已确立的规则）。合并走既有
`PatchTransaction2D` 事务路径，保持 patch-outside stable ID 不变与 global oracle 二次校验。

**边界**：AMReX 的 `small_volfrac`（`Src/EB/AMReX_EB2_2D_C.cpp:164,187`）是把小单元直接
cover 掉再靠 flux redistribution 补守恒，那是**求解器侧**策略。本产品导出 polyMesh 给
通用求解器，只能**合并不能删除**（`AGENTS.md` §4.10）。只取「用体积分数阈值定义小单元」
这个概念。

### W4 收口证据（用户已确认要真实物理基准）

1. **扩展已有 MMS 阶梯**。`tools/verification/openfoam_harmonic_mms.py` 已经实现
   configure/evaluate、在 OpenFOAM 自己写出的边界面心上取 Dirichlet 值、按
   `writeCellVolumes` 做体积加权。当前 V1c 只到 `(minimum,boundary)=(6,8)` 三级，
   观测阶从 **1.63 掉到 1.06**。A/B 落地后把阶梯延到 (7,9)/(8,10)，门是：三范数严格
   下降，且**观测阶不再继续下滑**。这直接检验修复是否真的改善了离散精度，而不只是让门变绿。
2. **Re=40 定常圆柱绕流**，`icoFoam`/`simpleFoam`，比对阻力系数与回流长度的文献值，
   在修好的阶梯上逐级跑。新增 `tools/verification/cylinder_re40.py`，结果进 `artifacts/r2/`。
3. 每一级都必须真实跑 `opencfd/openfoam-run:2606` 的 `checkMesh` 与独立读取器
   `tools/verification/check_openfoam2d.py`（`AGENTS.md` §4.12），不得用内部 reader 冒充。

## 7. 已经实测排除的路线：不要重走

| 试过的做法 | 实测结果 |
|---|---|
| 对已遍历的 **face loop** 做退化顶点塌缩 | 共享同一 spur 的两个 face 做出不同决定，局部图不再共形，leaf 报 `local fluid region contains a hole` |
| 对 **per-leaf 点集** 在建半边图之前焦合 | cell-side 点与邻居 leaf 共享，邻居独立焦合 → 共享边不一致，audit 报 `unclassifiedBoundaryEdges=24` |
| 同时焦合 embedded fragment 坐标 | 无效，24 个不变；跨 leaf 分歧才是原因 |
| 把首层厚度随 level 折半救 hybrid | 无效，只把 L9 失败从 132 s 拖到 803 s |
| 粗暴聚并过渡带 Cut-cell | 三维已实测否决，nonOrtho 恶化到 111°/129° |

另外记住：**不要往 `src/cutcell/CutCell2D.cpp` 里加几何容差**。那里的
`pointNear` / `scalarNear` 是
`|a-b| <= tol.absolute + tol.relative*|a-b|`，等价于**绝对 1e-12**（传给 tolerance 的
magnitude 是差值本身，不是任何局部长度）。这是一个已知的尺度依赖缺陷，但修它的正确
位置是 registry，不是这里。

## 8. 硬规则（违反其中任何一条，绿灯也不算通过）

1. **不得降低** `SolverQualityPolicy2D` 或 Q1 contract 的任何阈值，不得删除坏单元，
   不得隐藏告警（`AGENTS.md` §4.10）。放宽面积不变量是**唯一**已获用户批准的例外，
   且必须用第 4 节 A2 的推导界并写进文档。
2. **不复制 GPL 源码**，不新增外部依赖。全程 Level C/D，每消费一条借鉴就回
   `docs/R2_REFERENCE_AUDIT_CN.md` 把「仅确认存在 / 仅列出目录」升级为实读记录。
3. **不引入三维头文件、不链接三维 library**。三维仓库只当设计参考。
4. **壁面点必须落在原输入 segment 上**。不做 smoothing、不移动顶点。
5. **wall time 不进任何逐字节比较**。seconds 只写进标记 `reproducible: false` 的
   `.hybrid.profile.json`；确定性的调用计数才可以跨运行对照。
6. 每个里程碑除项目测试外，必须通过独立读取器；OpenFOAM 能跑时必须真实执行
   `checkMesh`（`AGENTS.md` §4.12）。

## 9. 每条改动的验收清单

```text
[ ] ctest --test-dir build --output-on-failure        # 不得低于 83/83
[ ] refinement_ladder.py 改动前/后两份 manifest 一起留档
[ ] 旧命令行（不带新开关）产物逐字节不变，或明确说明哪些变了、为什么
[ ] 真实 checkMesh + check_openfoam2d.py
[ ] 同输入重复生成 points/faces/owner/neighbour/boundary 逐字节相同
[ ] 最小失败案例入 tests/，不是只修不留证据
[ ] 借鉴条目回填 R2_REFERENCE_AUDIT_CN.md
```

逐字节对照的做法（本轮用过，很好用）：

```bash
git archive <base-commit> | tar -x -C /tmp/base
cd /tmp/base && cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCARTMESH2D_BUILD_TESTS=OFF && cmake --build build -j8
# 两边跑同一条旧命令行，再 shasum -a 256 逐个比
```

## 10. 本轮已交付的三个提交

| commit | 内容 |
|---|---|
| `82e043fc6ae52f9317e1bbd66b60bdf8c8a94504` | 加密阶梯门 + 顶层 H4 phase attribution（关闭 R1F §7 前置条件） |
| `0f26249ed13f1cb896193ac7ce8fb512365a0036` | 给 unsupported leaf 命名原因，而不只是计数 |
| `a9d21adc2ae94e6204a541d7b5141b44ff9f381c` | grid-corner 焦合预算改成 opt-in policy（默认行为零变化） |

三者都不改变 level 6–9 的既有行为；W0 的逐字节证据见
`docs/R2_REFINEMENT_ROBUSTNESS_CN.md` 第 8 节。


