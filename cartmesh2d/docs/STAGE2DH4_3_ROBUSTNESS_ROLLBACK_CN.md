# CartMesh2D H4-3：稳健性、事务回退与明确报错

日期：2026-08-27
状态：H4-3 完成；H4-V 未开始

## 1. 阶段边界

本阶段不重做 H4-1/H4-2，也不实现 local layer dropping、复杂 termination、三维或
H4-V 的复杂轮廓验证矩阵。目标是让已有 H4 候选具备正式的提交/回退语义，并对凹角、
尖角和窄缝给出可复现、可检查的结果。

## 2. 事务语义

`buildRobustH4Mesh2D` 只接收不可变的原始 wall 与 H4-1 候选，执行顺序如下：

1. H4-1 成功时，在独立结果对象中构建 H4-2 hybrid 候选；完整通过拓扑、面积和既有
   solver-quality 门禁后，才以 `mesh_mode=hybrid` 提交。
2. H4-1 或 H4-2 候选失败时，保留失败阶段、原因、消息、chain/vertex、请求厚度和安全
   厚度；随后从原始 wall 重新构造纯 Cut-cell 网格。
3. 回退路径复用既有 quadtree balance、Cut-cell、small-cell analysis、agglomeration、
   `buildSolverTopology2D` 和 solver-quality 门禁。全部通过才返回
   `mesh_mode=pure_cutcell_fallback`，同时固定写出 `layer_status=failed`。
4. 回退候选也不合格时，整体返回 `mesh_mode=failed` 并保留真实质量指标，不降低阈值、
   不删除坏单元、不把失败改写成成功。

CLI 通过末尾 `--robust-fallback` 启用该编排；未给此参数时保持原 H4-2 fail-closed 行为。

## 3. 最小失败案例与结果

| 案例 | H4-1 原因 | 首个实体/阈值 | 回退 cells | 面积误差 | max non-orth | min face weight | min volume ratio |
|---|---|---|---:|---:|---:|---:|---:|
| L 形凹角 | `concave_corner` | chain 0 / vertex 3 | 568 | `7.11e-15` | `41.4237°` | `0.107843` | `0.0909091` |
| 双方块窄缝 | `thickness_exceeds_safe_limit` | requested `0.11`, safe `0.09` | 931 | `-7.11e-15` | `49.9871°` | `0.0501836` | `0.0487805` |
| 尖三角 | `sharp_corner` | chain 0 / vertex 1 | 720 | `0` | `50.8302°` | `0.0586460` | `0.0378185` |

三例的回退网格都通过未修改的阈值：non-orthogonality `70°`、face weight `0.05`、
volume ratio `0.01`。另有制造的“回退本身 solver-quality 失败”回归，确认整体失败并在
`fallback_failure_message` 中报告指标。

最小输入保存在：

- `examples/acceptance/concave.xy`
- `examples/h4_3/narrow_gap.xy`
- `examples/h4_3/sharp_triangle.xy`

## 4. 独立检查与 OpenFOAM

独立 Python VTK/JSON 读取器对三例重新计算 polygon 正面积、edge incidence、交叉和面积
总和，均为 PASS，且 overlap/non-manifold 均为 0。

本机 Docker 已有 `opencfd/openfoam-run:2606`。使用
`checkMesh -constant -allTopology -allGeometry` 对三套真实挤出 case 检查，均输出
`Mesh OK`：

- concave：568 cells，max non-orth `41.4237°`，min weight `0.107843`，
  min volume ratio `0.0909091`；
- narrow gap：931 cells，max non-orth `49.9871°`，min weight `0.0501836`，
  min volume ratio `0.0487805`；
- sharp triangle：720 cells，max non-orth `50.8302°`，min weight `0.0586460`，
  min volume ratio `0.0378185`。

## 5. 产物

`artifacts/h4_3/` 保存三例的：

- `.h4.json` 事务/失败语义报告；
- `.fallback.vtk`、`.fallback.cm2d`；
- `.fallback.solver.vtk`、`.fallback.solver.cm2d`；
- base/solver quality JSON；
- 独立检查 JSON。

H4-V 的圆、椭圆、翼型/复杂轮廓和窄缝综合验证仍是下一阶段，不在本次完成声明内。
