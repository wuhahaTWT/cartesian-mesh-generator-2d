# 阶段 7 Codex 执行方案：贴体近壁层 + 自适应笛卡尔主体的混合网格

> 定位：在阶段 6 已经完成并通过验收的前提下，把现有 `STL → adaptive Cartesian/octree → Cut-cell → solver mesh` 主链扩展为：
>
> `STL → 几何特征识别 → 贴体近壁层 → 自适应 Cartesian/octree 主体 → 共形过渡/局部 Cut-cell fallback → 质量修复 → OpenFOAM/VTK 输出`
>
> 阶段 7 关闭时，应达到“v1.0 可用网格生成器”标准，而不是继续无限增加研究功能。

---

## 0. 文献依据与本项目采用的路线

参考赵宁等 2025 年综述中“混合笛卡尔网格方法”的定义：物面附近采用结构或非结构贴体网格，外部采用笛卡尔网格，并在过渡区域通过切割或投影方法连接成一套完整非结构网格体系。文中还总结了刘周等 2009 年的路线：先生成基于物面曲率自适应的笛卡尔初始网格，再用投影方法拟合壁面、用特征恢复处理凹角和改善质量，最后将柱形单元沿法向分层。

本项目不采用 overset 作为阶段 7 主路线，因为当前已有 Cut-cell、octree、完整 solver topology，最自然的扩展是“共形混合网格”，避免再引入另一套供体搜索/插值/守恒传递系统。

本项目阶段 7 采用：

1. adaptive Cartesian/octree 仍是主体；
2. 物面附近构造可投影、可贴体的 near-wall front；
3. front 投影到 STL，并恢复 sharp feature；
4. 由物面沿局部法向生成多层贴体单元；
5. 最外层 layer front 与 Cartesian/octree 通过已有切割/拓扑能力共形连接；
6. 局部无法安全生层时，允许退回阶段 6 的 Cut-cell，但必须显式报告 fallback 区域；
7. 不允许为了“覆盖率 100%”生成负体积、翻转层、非法多面体。

---

# 1. 阶段 7 总验收目标

阶段 7 关闭时，用户应只需提供：

- 一个合法 STL；
- 计算域设置；
- 远场/局部 refinement 设置；
- 近壁层开关；
- `nLayers`；
- `firstLayerHeight`；
- `growthRate`；
- `maxLayerThickness` 或总层厚；
- sharp-feature angle；

即可生成：

- adaptive Cartesian/octree 主体；
- 贴体近壁层；
- 必要位置的局部 Cut-cell fallback；
- 一套连续、共形、守恒的完整控制体网格；
- VTK/VTU 可视化；
- OpenFOAM `polyMesh`；
- 完整质量报告；
- layer coverage / fallback 原因报告。

### 阶段 7 不做

- 不做 RANS/LES 求解器；
- 不做 overset 网格；
- 不做动态运动网格；
- 不做完整 CAD/STEP/IGES 几何内核；
- 不承诺自动修复任意脏 STL；
- 不做 GPU/MPI 大规模重写；
- 不为了 y+ 自动猜测流动条件。

`y+` 相关功能只能在用户提供速度、黏性参数和估算模型时作为辅助换算；核心输入仍是明确的 `firstLayerHeight`。

---

# 2. Codex 全局执行纪律

1. **阶段 6 没有完整通过，禁止进入阶段 7。**
2. 每次只执行一个 `7.x` 子阶段，通过验收后停下，向用户汇报，再等待继续。
3. 每个 `7.x` 至少一个独立 commit；禁止一个 commit 横跨多个子阶段。
4. 修改前先读：`PROJECT_BRIEF / AGENTS.md / README / STAGE6_REVISED_PLAN / 本文件`。
5. 修改前必须输出：
   - 当前 git commit；
   - git status；
   - 本子阶段目标；
   - 将修改/新增的文件；
   - 测试命令；
   - 真实几何验收命令。
6. 不允许新建第二套平行的 geometry / topology / writer 系统；必须复用阶段 6 的统一数据模型。
7. 任何质量修复必须先加入最小失败用例，再修改实现。
8. 不允许用“测试绿了”代替真实 STL + 外部 `checkMesh` 验证。
9. 不允许通过删除坏单元、翻转孤立面、降低阈值、忽略 warning 来假装通过。
10. 每个失败都必须输出可定位的 cell/face/vertex ID、坐标、来源 leaf/surface triangle ID。
11. 若发现当前架构无法完成本子阶段且需要大重构：停下，只提交设计说明，不直接重写。

---

# 3. 7.0 阶段 7 基线冻结与接口设计

## 目标

冻结阶段 6 终态，定义混合层数据接口，不生成任何新 layer。

## 必须完成

- 重新执行阶段 6 全套回归；
- 保存阶段 6 的 topology fingerprint / quality report / OpenFOAM hash；
- 定义 near-wall layer 配置结构；
- 定义 surface feature、layer column、layer front、fallback reason 等核心数据结构；
- 明确 layer cell 最终必须进入现有统一 solver topology，而不是单独 writer。

## 建议配置字段

```yaml
nearWallLayer:
  enabled: true
  nLayers: 5
  firstLayerHeight: 1.0e-4
  growthRate: 1.2
  maxTotalThickness: 2.0e-3
  featureAngleDeg: 45
  collisionSafetyFactor: 0.8
  allowCutCellFallback: true
```

## 验收

- 开启 `enabled: false` 时，与阶段 6 输出 topology/hash 完全一致；
- 只增加接口，不改变任何原有网格结果。

---

# 4. 7.1 STL 表面拓扑、法向与特征边/角识别

## 目标

建立可靠的贴体层几何基础。

## 必须实现

- triangle adjacency；
- patch 连通区；
- 一致法向检查与传播；
- feature edge：二面角超过阈值；
- boundary/non-manifold edge 检测；
- convex / concave feature 分类；
- feature vertex；
- 局部曲率尺度估计；
- 最邻近三角形和最近点查询复用现有 BVH/KD 能力。

## 最小测试

- 平面；
- cube：12 条 sharp edge；
- wedge；
- concave L-groove；
- sphere：无 sharp edge；
- 非流形 STL：明确拒绝。

## 验收

- 特征边/角可输出 VTK；
- 相同 STL/阈值结果完全确定；
- 不允许通过 STL 三角面顺序改变结果。

---

# 5. 7.2 近壁尺寸场与 layer feasibility

## 目标

在真正生成贴体层之前，先判断“这里能不能安全生层”。

## 尺寸场必须考虑

- `firstLayerHeight`；
- growth rate；
- nLayers；
- 总厚度；
- 局部曲率半径；
- 局部 STL 三角形尺度；
- narrow gap；
- 距离相对面/另一物体的 clearance；
- sharp feature 附近的层厚限制。

## 必须输出

每个 surface vertex / face 的：

- requested thickness；
- admissible thickness；
- selected thickness；
- 可生成层数；
- 若降低层数/停止层，原因是什么。

## 规则

- 不允许相邻区域层厚突然跳变；
- 必须有受限平滑；
- 窄缝两侧 layer 不得互穿；
- 凹角区域允许提前终止或 fallback。

## 验收

在 sphere、cube、双平行壁窄缝上显示合理、连续、确定的 thickness field。

---

# 6. 7.3 单层贴体 front：投影方法 MVP

## 目标

先只做 **一层**，验证“Cartesian near-wall front → STL 投影 → 合法贴体单元”的核心几何链。

## 方法

参考文献中的投影思想：

1. 从物面附近的 finest Cartesian/octree 区域提取候选 front；
2. 对候选 front 顶点进行 closest-point 投影；
3. 普通区域投影到 surface triangle；
4. feature edge 附近约束到 feature polyline；
5. feature vertex 附近约束到 feature point；
6. 进行局部质量优化，但不得破坏几何约束；
7. 构造第一层 wall-adjacent cells。

## 第一轮限制

只要求：

- plane；
- cylinder；
- sphere；
- NACA0012 挤出翼段；

暂不解决复杂凹角和窄缝。

## 验收

- wall vertices 到 STL 的距离处于几何容差内；
- 无负体积；
- 无翻转面；
- wall patch 连续；
- cell closure 正常；
- OpenFOAM 能读取。

未通过不得进入多层挤出。

---

# 7. 7.4 Sharp feature recovery

## 目标

解决简单投影会“抹圆”棱边、尖角和凹角的问题。

## 必须实现

- feature-edge constrained vertex；
- feature-corner constrained vertex；
- convex edge 两侧 layer front 分离法向；
- concave edge 的安全终止/收缩；
- 局部 vertex/face smoothing 不能跨越 feature；
- feature topology 在 layer 外缘也要保持一致。

## 最小测试

- cube；
- 30°/60° wedge；
- concave L-shaped groove；
- 简单翼型 trailing edge。

## 验收

- 几何尖锐特征不会被无意平滑；
- 不产生 self-intersection；
- 不产生负体积；
- feature 位置误差有明确报告。

---

# 8. 7.5 多层法向挤出

## 目标

实现真正有 CFD 意义的 near-wall body-fitted layers。

## 必须实现

- surface/front vertex normal；
- angle-weighted 或 area-weighted normal；
- feature-aware normal；
- 多层位置：`h1, h1*r, h1*r^2 ...`；
- 局部 layer count reduction；
- layer termination；
- column topology；
- prism / wedge / hex-like polyhedral layer cell 的统一表示；
- layer volume、face orientation、owner/neighbour。

## 注意

不要强制所有单元必须叫“prism”。对于复杂多边形 front，允许生成一般 polyhedral column cell，但必须保持一套统一 FV topology。

## 验收

- flat plate：指定层高精确；
- cylinder/sphere：层厚增长连续；
- wing：前缘/后缘无穿插；
- 所有层 cell 正体积；
- layer growth 与配置一致。

---

# 9. 7.6 碰撞、窄缝、凹角和层终止

## 目标

这一步决定阶段 7 是“演示程序”还是“可用网格器”。

## 必须实现

挤出每一层前做 collision prediction：

- ray/swept-segment 与 STL 查询；
- layer-front 自碰撞；
- 与其他 component 的碰撞；
- narrow-gap clearance；
- concave corner inversion 预测。

## 允许的安全处理顺序

1. 局部降低 layer thickness；
2. 降低 layer count；
3. 对相邻区域平滑终止；
4. 最后才 Cut-cell fallback。

## 禁止

- 穿透后再修；
- 负体积 cell 后简单删掉；
- 单点突然从 N 层变 0 层；
- fallback 不写原因。

## 验收

- 双壁窄缝；
- cube 内凹角；
- 两个相近圆柱；
- 多 component STL。

每个案例输出 layer coverage heatmap + fallback reason map。

---

# 10. 7.7 Layer ↔ adaptive Cartesian/octree 共形过渡

## 目标

把贴体层与阶段 6 主体真正接成 **一套网格**。

这是阶段 7 最关键的拓扑门。

## 推荐路线

- layer outer front 视为新的内部边界；
- 删除与 layer 占据体积冲突的 Cartesian/octree cells；
- 对 transition band 使用阶段 6 已有 Cut-cell / clipping / agglomeration 能力；
- 生成 layer-front ↔ Cartesian 的唯一共享 faces；
- 不允许 overlapping cells；
- 不允许 hole；
- 不使用 overset interpolation。

## 必须保证

- volume conservation；
- shared face exact match；
- owner/neighbour 唯一；
- coarse/fine transition 与 layer transition 可同时存在；
- layer wall patch 与其他 patch 保留；
- fallback Cut-cell 与 layer cell 可以相邻。

## 最小测试

- plane layer + uniform Cartesian；
- sphere layer + adaptive octree；
- layer front 穿过 2:1 coarse/fine interface；
- 部分 layer + 部分 cut-cell 的同一物体。

## 验收

内部：

- negative volume = 0；
- nonclosed = 0；
- duplicate shared face = 0；
- hole = 0；
- overlap = 0；
- unmatched internal face = 0。

外部：

- OpenFOAM `checkMesh`；
- `checkMesh -allTopology`；
- `checkMesh -allGeometry`。

出现阻断错误，不进入下一阶段。

---

# 11. 7.8 混合网格质量评估与局部修复

## 目标

阶段 6 quality 模块扩展到 layer/hybrid 特有指标。

## 增加指标

- layer thickness ratio；
- layer growth smoothness；
- wall-normal alignment；
- layer cell aspect ratio；
- local orthogonality；
- face pyramid；
- skewness；
- concavity；
- transition-cell quality；
- feature-edge quality；
- layer coverage；
- fallback ratio。

## 修复优先级

1. constrained smoothing；
2. local thickness reduction；
3. layer count reduction；
4. local topology repair；
5. Cut-cell fallback。

不得改 STL 几何来换取质量。

## 验收矩阵

分辨率逐级：

- plane / channel；
- cube；
- sphere；
- cylinder；
- NACA0012 wing；
- concave L-body；
- narrow gap；
- multi-body；
- Bunny 或复杂工程 STL。

简单几何要求贴体层完整覆盖；复杂几何允许局部 fallback，但必须稳定、可解释、可定位。

---

# 12. 7.9 v1.0 产品门

## 目标

这一关通过后，不再把它叫“阶段性算法 demo”，而是发布 `v1.0`。

## 用户流程必须简化为

```bash
cartmesh generate input.stl config.yaml -o case/
cartmesh check case/
cartmesh report case/
```

具体命令名称服从现有 CLI，不强制新建平行命令。

## 必须具备

- 一个统一 CLI；
- YAML/JSON 配置；
- `layer on/off`；
- uniform/adaptive 统一入口；
- Cut-cell fallback 开关；
- OpenFOAM 输出；
- VTK/VTU 输出；
- quality report；
- layer coverage report；
- deterministic seed / stable topology；
- 失败时非 0 exit code；
- 至少 5 个公开样例配置；
- README 从“开发说明”改为“用户 10 分钟能跑通”。

## v1.0 最终真实性门

必须实际跑至少：

1. cube；
2. sphere/cylinder；
3. NACA0012 wing；
4. narrow-gap / concave geometry；
5. 一份真实工程 STL。

每个案例：

- 重新从 STL 生成；
- 内部 quality PASS；
- independent reader PASS；
- OpenFOAM `checkMesh` PASS；
- ParaView 可视化确认 wall layer / transition / Cartesian core；
- 重复两次 topology fingerprint 稳定；
- 记录峰值内存与耗时。

阶段 7 只有在这些材料被保存进 `artifacts/stage7_v1/` 后才允许关闭。

---

# 13. 阶段 7 完成后的产品定义

阶段 7 完成后，本项目应可以明确称为：

**“基于 STL 的自适应笛卡尔 / Cut-cell / 贴体近壁层混合 CFD 网格生成器 v1.0”**。

它应适合：

- 研究项目；
- 课程/论文；
- OpenFOAM 前处理；
- 中等复杂度封闭 STL；
- 高自动化的 Cartesian-dominant 网格生成。

但不要声称达到 Fluent Meshing / STAR-CCM+ / Pointwise 的成熟度。后者还包含长期积累的 CAD 修复、脏几何、超复杂 feature、各种物理边界层策略、并行扩展、交互 GUI、大量工业回归等。

阶段 7 后若继续做，应该进入 `v1.x / v2.0` 的增强项，而不是继续把“核心网格生成器没做完”作为理由无限延期。

---

# 14. 给 Codex 的阶段 7 开工提示词

把下面内容直接给 Codex：

> 先完整阅读项目中的 PROJECT_BRIEF、AGENTS.md、README、STAGE6_REVISED_PLAN_CN.md 和 STAGE7_HYBRID_CARTESIAN_CODEX_PLAN_CN.md，并检查 git log、git status、当前 tag 和阶段 6 的真实 artifacts。
>
> 当前只允许执行阶段 7.0，不得进入 7.1，更不得提前写 layer extrusion。
>
> 先向我汇报：
> 1. 当前正式版本和 commit；
> 2. 阶段 6 是否满足进入阶段 7 的全部前置条件；
> 3. 当前统一 mesh/topology/quality/export 数据结构分别在哪里；
> 4. 阶段 7.0 需要新增哪些接口；
> 5. 将修改哪些文件；
> 6. 将新增哪些测试；
> 7. 如何证明 layer disabled 时阶段 6 输出 bitwise/topology 等价。
>
> 在我确认前不要修改代码。确认后只实现 7.0，完成后运行真实回归、提交独立 commit，然后停下等待下一步。

