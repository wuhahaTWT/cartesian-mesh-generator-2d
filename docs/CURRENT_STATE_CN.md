# 当前状态：逐案例 × 逐门

日期：2026-09-04
本文是**单一事实来源**。历史轮次文档各自保留当时的数字作为证据，不再重复登记
"当前"状态——那正是过去出现五份文档互相矛盾的原因。

## 0. 怎么读这张表

四个门是独立的，通过一个不代表通过另一个：

| 门 | 判什么 | 权威来源 |
|---|---|---|
| topology audit | 无重复/孤立/非流形边，`fluid_area = domain_area - solid_area` | `TopologyMesh2D::valid()` |
| solver quality | `SolverQualityPolicy2D` 硬限（faceWeight 0.05、volRatio 0.01、nonOrtho 70°、boundary skewness 4） | `evaluateSolverQuality2D` |
| OpenFOAM `checkMesh` | 外部独立判定，阈值比本仓库**宽**（如 boundary skewness 20） | `opencfd/openfoam-run:2606` |
| Q1 contract | 无量纲、分类型的产品质量合同，比上面三个都**严** | `QualityContract2D` |

Q1 比 `checkMesh` 严是**设计如此**，不是缺陷。所以"`checkMesh` Mesh OK 且 Q1 FAIL"
是当前五个案例的正常状态，不能只报前半句。

## 1. 五个 H4 验收案例（level 8 及以下，CI 覆盖）

| 案例 | topology | solver quality | checkMesh | Q1 | Q1 hard 主项 |
|---|---|---|---|---|---|
| circle | PASS | PASS | Mesh OK | **FAIL** | volume ratio 56、face weight 24 |
| superellipse | PASS | PASS | Mesh OK | **FAIL** | 短面 21（三类）、face weight 16、volume ratio 12 |
| concave_l | PASS | PASS | Mesh OK | **FAIL** | volume ratio 107、face weight 64 |
| narrow_gap | PASS | PASS | Mesh OK | **FAIL** | volume ratio 299、face weight 99 |
| sharp_trailing_edge | PASS | PASS | Mesh OK | **FAIL** | volume ratio 156、face weight 145、短面 38（三类） |

数字来自 `artifacts/q1/*.quality-contract-baseline.json`。

补充状态：

- **sharp_trailing_edge**：H4-3 hybrid 与 `checkMesh` 通过，但 R1 的 mutable/mutable
  迁移正式结论是 **NO-GO**（`docs/R1F_PATCH_LOCAL_CLOSEOUT_CN.md` 第 6 节），
  Q1 短面 `face/local_h = 4.2e-4` 远低于 0.01 硬限，未解决。
- **superellipse**：Q2 修掉了原先的微短面 hard failure（`9.8e-9` 那一条），
  但仍有其他短面项 FAIL。Q2 因此登记为 **PARTIAL**，不是完成。
- **narrow_gap**：需要 Q3/Q4/R1 事后修复才能压低 hard 计数，根因是层厚与缝宽
  对撞（见第 4 节）。

## 2. 加密上限（R2/W1 实测）

| 路径 | 稳定上限 | 上限处规模 | 失败点与原因 |
|---|---|---|---|
| 纯 Cut-cell（`cartmesh2d_cli`） | **level 11** | 180468 stabilized cells、277336 leaves、约 7 s | level 12 尚未测；L10/L11 topology、solver quality、OpenFOAM 导出全部通过 |
| hybrid 边界层（`cartmesh2d_hybrid_cli`） | **level 8** | 2072 solver cells | level 9：faceWeight 0.024、volRatio 0.0023、nonOrtho 78.6°；根因是壁面切向分辨率被输入折线顶点数锁死 |

纯 Cut-cell 的 non-orthogonality 随加密总体**变好**（L6 68.6° → L11 42.3°）；
W1 已消除原 L10 单点 grid-corner spur，并把同一共享构造路径验证到 L11。
L11 的独立 reader 与真实 OpenFOAM v2606 `checkMesh -writeAllFields` 均通过。

改前/改后证据为 `artifacts/r2/w1-before-manifest.json` 与
`artifacts/r2/w1-after-manifest.json`；W0 基线仍保留在
`artifacts/r2/w0-baseline-manifest.json`。测量入口
`tools/verification/refinement_ladder.py`，详情 `docs/R2_REFINEMENT_ROBUSTNESS_CN.md`。

**这张表的规模数字不是产品设置。** `refinement_ladder.py` 固定
`minimum_level = level - CUTCELL_MINIMUM_LEVEL_OFFSET`（纯路径偏移 2，hybrid 偏移 1），
所以每一级都是**全域近均匀**的压力网格，不是边界自适应网格。同一个 level 11 在 CLI 默认设置
（`padding-fraction 0.25`、`minimum-level 0`）下只有 23896 leaves / 14048 cells，
与表中 277336 / 180468 差 11.6 倍。

这是**刻意的**：ladder 要最大化 grid-corner 擦碰与 cut cell 数量，才能暴露构造核的失败模式。
但因此 ladder 的 cell 数与耗时**不得**当作产品代表值引用，任何面向产品的 sizing 结论都要在
真实设置下另行测量。


## 3. 各轮范围的准确状态

| 轮次 | 状态 | 边界 |
|---|---|---|
| 2D-0 … 2D-V | 关闭 | exterior 语义下的验收测试全部在 CI 运行。五份 `STAGE2D3/4/5/6/V` 文档顶部的 `REOPENED` 横幅是 2026-08-22 之前的历史 |
| H1 … H4-3 | 关闭 | H4-3 是当前 hybrid 算法之源 |
| Q0、Q1 | 关闭 | 建立了质量词汇与无量纲合同；两者都是**诊断**，不修网格 |
| Q2 | **PARTIAL** | superellipse 微短面已修；narrow_gap / sharp-tail 短面未修 |
| Q2-A | 关闭 | 共享交点构造已进入 hybrid 与纯 Cut-cell 默认路径 |
| Q3-1/2/3、Q4-1 | 关闭但已饱和 | 见 `docs/Q3_Q4_TERMINATION_QUALITY_CN.md`。文档明确记录"不建议继续扩" |
| R1-A、R1-B、R1-E、R1-F | 关闭（narrow-gap 范围） | R1F 的 CLOSED 只覆盖 narrow-gap |
| R1-C、R1-D | **未关闭** | 两份文档都明确写"不得登记完成" |
| R1-G | **不存在** | 原计划的正式冻结验收，从未做 |
| R2/W0 | 关闭 | 加密阶梯门 + 顶层 phase attribution |
| R2/W1 | **关闭** | 纯 Cut-cell 共享构造、解析面积界、跨重建 wall patch 语义与 L6–L11 阶梯均通过 |
| R2/W2、W3、W4 | 未开始 | 下一阻塞仍是 hybrid 壁面切向再分，见 `docs/R2_HANDOFF_CN.md` |

## 4. 两条根因（正在处理的对象）

**A. 几何吸附半径是机器精度量级。** `snapFractionOfLocalH` 默认
`64·DBL_EPSILON`，level 10 下吸附半径约 `8e-17`，而实际微面是 `6e-6 · h`——
几何上等于零，却是吸附半径的四亿倍。加密只会让"折线擦过格点"更常发生，
所以这是**随 level 必然恶化**的失败模式。W1 已用单一全局 registry 和
`gridCornerWeldFractionOfLocalH=1e-4` 解决；输入端点只允许算术 roundoff 归并，
几何预算只用于真实 grid corner，避免对称切点被拉过 leaf 边界。

**B. 壁面切向分辨率被输入折线顶点数锁死。** wall column 直接取输入
`BoundaryLoop` 顶点，`circle.xy` 只有 32 个顶点，因此 BL cell 数在**每一级都是
128**。level 9 的背景 `h ≈ 0.0113`，切向尺寸恒为 `0.196`，背景网格比壁面细 17 倍，
过渡环必须独自吸收整个落差 → 强斜多边形。

narrow_gap 的层厚对撞是同一族问题的另一面：4 层总厚 0.0599，缝宽 0.08，
两侧相加 0.12 > 0.08。这是 Q3/Q4/R1 修复框架存在的**唯一原因**。

## 5. 已知的数字歧义（不要当成矛盾）

- **narrow_gap solver cell 数**：3185（macOS/mesasdk）vs 3187（Ubuntu/GCC）是
  编译器差异，R1F 已明确声明**不宣称跨编译器逐字节一致**。3153 是 Q3-1 之后、
  3149 是 Q4-1 之后，都是修复带来的真实变化。H4-3/Q2-A 时代记录的 3189 与当前
  3185 之间的差异**没有注明来源**，属于待查项。
- **CTest 项数**：以 `grep -c add_test CMakeLists.txt` 为准。历史文档里
  40/45/51/53/72/73/74/75 各自都是当轮的正确值。
- `docs/R1_REFERENCE_DRIVEN_CONSTRUCTION_KERNEL_AUDIT_CN.md` 第 3.1 节标题写着
  "当前五案例事实"，但其数据来自一个**从未合并的 Q2-B 本地原型树**。标题是缺陷，
  数字不是仓库状态。

## 6. 维护规则

新增一轮时：把本文对应行改掉，在轮次文档里记录当轮数字。**不要**在
`AGENTS.md`、`README.md` 或轮次文档里再写一份"当前状态"——那会立刻产生第二个
事实来源。
