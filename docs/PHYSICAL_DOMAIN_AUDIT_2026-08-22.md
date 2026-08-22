# cartmesh2d CFD 物理域审计与纠正记录

日期：2026-08-22  
分支：`agent/native-2d-robustness`

## 1. 审计结论

此前 2D-3 ~ 2D-6 的默认物理域定义发生了方向性错误：

```text
旧实现：BoundaryLoop interior = fluid
正确默认：BoundaryLoop interior = solid
          Domain2D - solid interior = fluid
```

因此旧版虽然 Cut-cell、拓扑、small-cell、导出和可视化可以全部内部自洽，且 CI 可以全绿，但它默认生成的是**固体轮廓内部网格**，不符合本项目用于外部 CFD 前处理的产品目标。

该问题属于阻断级概念错误。旧的 2D-3/2D-4/2D-6 “PASS/CLOSED” 结论不能继续作为默认外流产品验收依据。

## 2. 与三维项目对照

三维 `cartmesh` 的 Cut-cell/流体组件逻辑已经采用正确的外流物理语义：

- STL `outside` -> fluid；
- STL `inside` -> solid。

二维项目现在明确要求默认语义与三维一致，不允许二维/三维同一仓库中出现相反的 fluid-side 定义。

## 3. 根因

### 3.1 几何分类被错误当成物理区域

`CellClass::Inside / Outside / Intersected` 本应只说明 cell 相对于闭合几何的位置。

旧 2D Cut-cell 代码却直接使用：

```text
Inside  -> Full fluid
Outside -> Empty
```

这等价于默认把闭合轮廓当成“流体包络”，而不是 CFD 固体壁面。

### 3.2 Cut-cell 局部边界方向也绑定了错误流体侧

旧局部 boundary graph 使用 CCW 输入边界，并按“内部在有向边左侧”构造流体 polygon。这不仅是 Full/Empty 标签反了；如果只交换 Inside/Outside 而不反转外流 embedded fragment，Cut-cell 闭环仍会错误。

### 3.3 测试把错误语义固化成了正确答案

旧测试显式要求：

- inside leaf -> full fluid；
- outside leaf -> empty；
- circle 总 fluid area = circle polygon area；
- circle topology `DomainBoundary = 0`。

因此 CI 全绿只能证明“实现符合旧错误定义”，不能证明 CFD 物理域正确。

### 3.4 文档没有先锁定 BoundaryLoop 的物理角色

早期总纲、架构和阶段计划写了 `fluid polygon`，但没有明确：

```text
BoundaryLoop 到底表示 solid wall 还是 fluid envelope？
```

这使得算法可以在错误物理侧上持续推进到 Stage 6，而常规拓扑审计仍无法发现。

## 4. 已实施纠正

### 4.1 显式 fluid-side policy

新增：

```cpp
enum class FluidRegion2D {
    Exterior,
    Interior
};
```

产品/API 默认：`Exterior`。

只有内部流、管道流等明确场景才显式使用 `Interior`。

### 4.2 默认外流 Full/Empty 语义

现在默认：

```text
geometric Outside -> Full fluid
geometric Inside  -> Empty solid
Intersected       -> retain exterior fluid polygon(s)
```

### 4.3 外流 embedded boundary 方向

输入 `BoundaryLoop` 仍标准化为 CCW。

对于 `Exterior` fluid，局部 embedded fragment 反向，使 retained fluid 始终位于 directed boundary 左侧，保证局部 polygon graph 的方向与几何意义一致。

### 4.4 多流体组件不能再静默丢失

强凹固体在一个叶单元内可能把外部流体切成多个互不连通的 polygon。旧单 `fluidPolygon` 模型对此只能 Unsupported 或丢片。

新增 `buildCutCells(...)` solver API：

- 一个 leaf 可发射多个 `CutCell2D` solver cells；
- CLI 将这些组件全部送入全局 topology；
- 每个 emitted component 获得唯一 deterministic `sourceId`；
- 不允许丢掉较小组件，也不允许跨固体制造假桥。

单组件 `buildCutCell(...)` API 仍保留；遇到多个真实 fluid components 时必须显式 Unsupported，并要求调用 `buildCutCells(...)`。

### 4.5 真正的 local hole 仍显式 Unsupported

若一个 solid loop 完全包在单 leaf 内，则外部流体是带孔 polygon。当前 topology 数据模型仍是 simple polygon，不支持 polygon-with-holes。

现在会识别反向 hole loop 并显式失败：

- 不把 hole 当第二个正面积 solver cell；
- 不填掉固体；
- 不制造假连接。

当前处理策略是提高 refinement，未来再加入 polygon-with-holes / multi-loop topology。

### 4.6 CLI 加入物理域硬门

默认 CLI 明确：

```text
boundary.xy = solid wall
fluid-region = exterior
```

并在进入 topology 前数值检查：

```text
source_fluid_area == domain_area - solid_area
```

内部流模式则检查：

```text
source_fluid_area == solid_area
```

因此即使某次代码改动又把流体侧反过来，CLI 也必须非零退出。

### 4.7 默认外流 topology 硬门

默认外流必须同时出现：

- `EmbeddedBoundary > 0`：固体壁面；
- `DomainBoundary > 0`：外部计算域边界。

这直接防止再次出现“圆内部有网格、外围全空，但 topology audit 仍为 0”的假成功。

### 4.8 CI acceptance 增加物理语义门禁

Stage 6 workflow 除原 topology/quality gate 外，增加：

- `viz.fluid_region == exterior`；
- `viz.boundary_role == solid_wall`；
- CM2D 中必须同时存在 embedded-wall patch 与 outer-domain patch；
- CLI 自身必须先通过 fluid-area physics gate。

### 4.9 文档和开发规则更新

已更新：

- `cartmesh2d/AGENTS.md`
- `cartmesh2d/README.md`
- `docs/PROJECT_BRIEF_CN.md`
- `docs/ARCHITECTURE_CN.md`
- `docs/STAGE_PLAN_CN.md`
- `docs/ACCEPTANCE_CN.md`

以后任何 agent/Codex 在本目录工作都必须先看到“默认闭合轮廓 = solid wall；默认 fluid = exterior”的硬规则。

## 5. 测试体系纠正

已把核心回归改成：

- 默认 exterior 直线切割解析面积/质心；
- explicit Interior 仍有独立解析回归；
- geometric Inside -> Empty solid；
- geometric Outside -> Full fluid；
- circle Quadtree 总面积 = `domain - circle`；
- external topology 同时存在 embedded wall 与 outer domain boundary；
- 新增 disconnected-fluid-component regression；
- 历史 shifted-circle small-cell fixture 改成显式 `Interior`，使其继续只测试稳定化算法，而不再偷偷定义产品默认物理侧。

## 6. 整个二维项目继续存在的 solver-readiness 缺口

物理侧错误纠正后，以下事项仍不能宣称“工业求解器网格已完全完成”：

1. **外部 DomainBoundary 目前只有一个总 patch 类型**。还没有 left/right/top/bottom 或 inlet/outlet/farfield 的命名边界。
2. **计算域控制仍偏演示型**。CLI 目前按物体 bbox + symmetric padding 自动生成外域；真实外流应支持显式 xmin/xmax/ymin/ymax 或上下游非对称距离。
3. **默认 padding=0.25*span 对真实外流通常太近**。它适合算法回归，不适合作为空气动力学远场默认值。
4. **当前导出是 VTK + CM2D**。CM2D 有 solver topology，但还不是 OpenFOAM/SU2/CGNS 等实际 CFD solver 的直接 case/export。
5. **单个 solid loop 为当前主路径**。多个独立障碍物、孔洞、多物体还需要正式数据模型。
6. **单 leaf 内真正 hole 仍 Unsupported**，需要更深 refinement 或未来 polygon-with-holes topology。
7. **没有贴体边界层/各向异性近壁层**。当前是 Cartesian/Quadtree + Cut-cell 路线，不应称为已有 boundary-layer mesh。
8. **质量优化仍有限**。拓扑合法不代表 aspect ratio/skewness 已达到工业 solver 最佳实践。

这些属于下一阶段 solver-readiness / mesh-quality 工作，与“默认流体侧反转”不是同一个问题，不能混在一起掩盖。

## 7. 当前状态

物理域定义和代码路径已纠正，并加入防复发门禁。

但本次修改是在关闭 PR 后直接推进到 `agent/native-2d-robustness`，为避免再次触发大量 GitHub Actions 邮件，**尚未通过新的 PR/Actions 对当前 correction head 做最终编译验收**。

因此当前诚实状态是：

```text
2D physical-domain bug: FIXED IN CODE
2D acceptance semantics: CORRECTED
Stage 2D-3/4/6 historical closure: INVALIDATED / REOPENED
current corrected head full CI: PENDING
```

在新的 current-head CMake/CTest/CLI/visualization 全链路通过之前，不再宣称“二维核心已完整完成”。
