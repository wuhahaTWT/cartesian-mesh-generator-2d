# Stage 2D-H4-3：局部降层、共形终止与事务回退

日期：2026-08-28

## 阶段结论

H4-3 不再把“复杂几何失败后整体退回 pure Cut-cell”定义为完成。当前产品路径为：

```text
完整请求层数
→ 局部层数规划与逐层 taper
→ 真实 stepped outer envelope / termination buffer
→ 共形 hybrid topology 与 solver-quality repair
→ 仅当前述候选均失败时才进入 pure Cut-cell fallback
```

凹角、尖尾缘和窄缝验收案例均以 hybrid 模式成功，没有触发全局 fallback。

## 实现边界

- H4-1 严格 builder 保持不变；H4-3 使用独立的 locally-reduced transaction。
- 每个 layer column 有独立实际层数，相邻列层数差不超过 1；0-layer 段沿原 wall 成为接口的一部分。
- `outerEnvelopeVertexIds` 记录由外层边、暴露 hair edges 和 0-layer wall edges 组成的真实 stepped envelope。
- H4-3 工作链只插入严格共线的 wall points，使用二次幂列数限制切向长宽比并避免与 dyadic quadtree 网格线产生系统性近重合；原始 wall 几何和 patch identity 不变。
- 终止带采用末层法向尺度和确定性 growth 候选 `{1.45, 1.55, 1.50}`；每个候选完整重建并通过同一 solver-quality 门后才提交。
- Boundary-layer cells 保持 immutable；termination/remainder cells 可进入 small-cell、agglomeration 和 solver repair。
- `HybridCellKind2D::Termination` 与独立 layer metadata 避免把 layer cells 计为普通 quadtree Cut-cell，也不从 synthetic key 推导假层级。
- Solver topology 新增严格受全局质量分数约束的凸 cell-pair agglomeration，用于删除凸分产生的小三角；immutable layer cells 不参与。

## 验收案例

统一参数：4 层，首层 `0.012`，growth `1.15`，remainder max/boundary level `8`。

| 案例 | layer cells / requested | 0-layer columns | termination edges | solver cells | growth | max non-orth | min face weight | min volume ratio |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| concave L | 1552 / 1664 | 4 | 52 | 5452 | 1.45 | 69.8253570 | 0.0570773 | 0.0106997 |
| sharp trailing edge | 1044 / 1088 | 2 | 22 | 3391 | 1.45 | 69.3953455 | 0.0513456 | 0.0108295 |
| narrow gap | 820 / 1024 | 0 | 40 | 3189 | 1.55 | 68.8133247 | 0.0515745 | 0.0118533 |

窄缝在相对壁面局部降至 2 层，缝外保留 4 层；几何 clearance 足以保留第一层，因此不强制制造 0-layer。尖尾缘实际包含 `4→3→2→1→0`，凹 L 只在危险邻域终止。

三例面积误差分别为 `-6.75e-14`、`2.13e-14`、`-3.38e-14`；interface audit、global topology、solver topology 和 solver-quality 均 PASS。

## 回归与性能

- 全部二维 CTest：`72/72 PASS`。
- H4-1/H4-2/H4-3 专项：`19/19 PASS`。
- 同一 256 点 NACA0012、Release、单线程、66,316 leaves / 66,254 solver cells 的 H3 对照：
  - `ffb55a1`：solver topology `5.227135 s`，total `7.461409 s`；
  - H4-3：solver topology `5.241151 s`，total `7.484180 s`；
  - 变化：solver topology `+0.27%`，total `+0.31%`，输出 cells/faces 和质量 extrema 完全一致；
  - 独立 RSS 复跑：基线 `270,057,472 B`，H4-3 `257,687,552 B`。硬件为 Apple M1 MacBook Air 8-core/8 GB，GCC 15.2.0，Release `-O3 -DNDEBUG`。

## OpenFOAM v2606

Docker image：`opencfd/openfoam-run:2606`。三个 H4-3 case 使用与 CI 相同的默认 `checkMesh -constant`，均报告 `Mesh OK`；CI 已扩展为 circle、superellipse 与三个 H4-3 case 全部执行真实 v2606 checkMesh。

额外诊断 `-allTopology -allGeometry` 会对 H4-2 circle 和 H4-3 general-polygon extrusion 同样报告 concave polyhedra。它是现有 writer 对一般二维 polygon 挤出时的已知扩展检查边界；本轮没有把默认 `Mesh OK` 冒充 extended-allGeometry PASS，也没有在 H4-3 中扩大到 writer 重构。该项留给后续 writer/H4-V 决策，本阶段没有开始 H4-V。
