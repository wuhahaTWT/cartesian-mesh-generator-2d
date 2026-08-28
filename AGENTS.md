# cartmesh2d 开发规则

## 0. 项目身份

本仓库是独立的原生二维 Cartesian / Quadtree / Cut-cell 网格生成器。
它与 `cartesian-mesh-generator` 原生三维仓库并行开发，但不依赖、不链接、
也不复制三维核心。

严禁把三维项目压成 `z=0` 来冒充二维实现。二维核心对象必须保持为
`Point2D`、`Segment2D`、`AABB2D`、`Polygon2D`、二维 Cartesian cell、
Quadtree leaf、CutPolygon 和 Edge2D。

## 1. CFD 物理侧定义（硬约束）

- 输入闭合 `BoundaryLoop` 默认表示固体壁面/障碍物轮廓；
- `Domain2D` 表示外部计算域；
- 默认流体区域是 `Domain2D - solid interior`；
- `Inside` 单元属于固体，不进入最终流体网格；
- `Outside` 单元属于流体；
- `Intersected` 单元必须保留边界外侧的真实流体 polygon；
- 外流网格必须同时具有 `EmbeddedBoundary` 与 `DomainBoundary`；
- 内部流只有调用者显式选择 `FluidRegion2D::Interior` 时才允许。

任何默认路径、CLI、验收案例或可视化若把闭合固体轮廓内部当成默认 CFD
流体域，均属于阻断级错误，不得以测试全绿或图片可见通过验收。

## 2. 仓库边界

二维代码全部位于本仓库根目录的 `include/`、`src/`、`apps/`、`tests/`、
`tools/`、`examples/`、`desktop/`、`docs/` 与 `artifacts/`。

- 不得引入 `cartmesh/*` 三维头文件或链接三维 library；
- 不得把本仓库重新作为长期分支塞回三维仓库；
- 二维和三维若共享算法思想，应通过文档和独立实现协调，除非用户明确批准
  创建稳定、独立、带版本的公共库；
- 历史验证文档中出现的 `cartmesh2d/` 子目录命令仅代表拆仓前环境，当前构建
  一律从本仓库根目录执行。

## 3. 当前阶段边界

已完成并形成回归门：

`2D-0 -> 2D-1 -> 2D-2 -> 2D-3 -> 2D-4 -> 2D-5 -> 2D-6 -> 2D-V`

后续精细化已完成：H1 sizing、H2 scalability、H3 solver topology、
H4-1 boundary-layer core、H4-2 conformal hybrid、H4-3 local dropping and
termination。未经用户明确批准，不得把后续阶段冒充当前已完成范围。

## 4. 核心真实性规则

1. 几何与拓扑正确性高于可视化。
2. 不得把删除相交格子或单元中心采样称为 Cut-cell。
3. Cut-cell 必须输出真实二维流体 polygon，并能计算正面积、质心及边界边。
4. 默认外流必须满足 `fluid_area = domain_area - solid_area`（tolerance 内）。
5. 默认外流拓扑必须同时出现固体壁面与外部计算域边界。
6. 自交、零面积、重复边、孤立边、非流形或分类冲突必须显式失败。
7. 相同输入和参数必须产生确定性的 ID、拓扑、报告和输出。
8. tolerance 必须集中管理，不得散落未命名魔法常数。
9. 几何修复必须保留最小失败案例和回归测试。
10. 不得降低 solver-quality 阈值、删除坏单元或隐藏告警来通过验收。
11. Boundary-layer、termination 与 remainder 必须进入统一共形 owner/neighbour
    拓扑；pure Cut-cell 只能作为明确的最后一级 fallback。
12. 每个重要里程碑除项目测试外，还必须通过独立读取器；OpenFOAM 能运行时
    必须真实执行 `checkMesh`，不得以内部读取器冒充。

## 5. 当前验证基线

- 原生二维 CTest：72 项；
- H4-3 验收：concave L、sharp trailing edge、narrow gap；
- H4-2 回归：circle、superellipse；
- OpenFOAM 验证镜像：`opencfd/openfoam-run:2606`；
- 当前 H4-3 事实来源：`docs/STAGE2DH4_3_LOCAL_TERMINATION_CN.md` 与
  `artifacts/h4_3/`。

任何新修改都必须保持工作区确定性、原有回归和真实输出证据。
