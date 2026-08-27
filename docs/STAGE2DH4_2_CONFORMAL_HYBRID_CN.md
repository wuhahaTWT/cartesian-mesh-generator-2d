# Stage 2D-H4-2：Conformal Hybrid Mesh

## 阶段目标

H4-2 将 H4-1 已生成的贴壁层 strip 与外部 Cartesian / Quadtree / Cut-cell 区域组装为一个统一、共形、可进入求解器质量门的二维混合网格。

核心约束：H4-1 strip 与原始 wall 在 H4-2 中视为固定输入；为了与 Cartesian 区域连接，只允许在 strip 外侧构造确定性的 transition fan，并对其外侧 remainder 执行 Quadtree refinement、Cut-cell、small-cell stabilization 与 solver-topology repair。

## 当前实现

生产入口 `buildConformalHybridMesh2D(...)` 自动从 H4-1 outer-envelope 边长、最后一层法向间距和 remainder boundary-level 对应的 Cartesian 尺寸 `h` 推导 transition plan。相同几何与 refinement 输入必须得到相同 ring count、切向细分和 transition 宽度。

统一网格包含三类 source cell：

- `BoundaryLayer`：H4-1 固定贴壁层；
- `RemainderCut`：transition polygon、remainder Cut-cell 及 stabilization 后的非 Cartesian polygon；
- `RemainderCartesian`：保留 Cartesian 形态的 remainder cell。

最终构造后执行：

1. remainder Quadtree 2:1 balance 与确定性检查；
2. remainder Cut-cell 和 small-cell stabilization；
3. layer + transition + remainder 全局 topology；
4. 重复 topology build，检查 canonical ID / connectivity 确定性；
5. H4-1 outer-envelope 共形接口审计；
6. 面积守恒；
7. mesh quality；
8. constrained solver-topology repair；
9. solver quality；
10. solver repair 后再次审计固定 outer-envelope 接口。

## 产品输出

`cartmesh2d_hybrid_cli` 可输出：

- `.hybrid.vtk`
- `.hybrid.cm2d`
- `.hybrid.json`
- `.hybrid.quality.json`
- `.hybrid.solver.vtk`
- `.hybrid.solver.cm2d`
- `.hybrid.solver-quality.json`
- 可选的 extruded OpenFOAM case

只有 solver topology 与 solver quality 均通过时才允许写 OpenFOAM case。

## 最终 CI 验证口径

H4 收尾 CI 使用 circle 与 superellipse 两个 solver-ready 案例，要求：

- 完整 native 2D stage CTest 通过；
- H4 hybrid topology / mesh quality / solver quality 均 valid；
- 2:1 balance 无残留 violation；
- fixed outer-envelope 接口无 single-owner、wrong-pair、non-two-valent 错误；
- OpenFOAM 输出通过项目外置的独立 reader；
- runner 若存在真实 `checkMesh` 则执行并保存日志，否则明确记录 `checkMesh=UNAVAILABLE`，不得冒充通过；
- acceptance 与 H4 solver-ready 证据作为 CI artifact 上传。

H3 的 100k solver-topology 性能已经在独立 H3 阶段完成并有固定报告；H4 最终 CI 不重新制造一个参数不同的 NACA workload 来冒充 H3 benchmark。H4 对 generic solver-topology path 的功能回归由完整 CTest 覆盖，H3 性能基线保持独立证据，不作为本阶段重复的阻塞门。
