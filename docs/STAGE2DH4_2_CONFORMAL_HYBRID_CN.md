# CartMesh2D H4-2：边界层与 Cartesian/Cut-cell 共形拼接

日期：2026-08-27  
状态：H4-2 solver-ready 收尾验证通过；真实 OpenFOAM `checkMesh` 仍需在安装 OpenFOAM 的环境补跑

## 1. 本阶段边界

H4-2 将 H4-1 的固定层数、闭合外流场 boundary-layer strip 接入现有二维 Cartesian / Quadtree / Cut-cell 网格。H4-1 strip 与原始 wall 在 H4-2 中是固定输入；只允许在 strip 外侧构造确定性的 transition fan，并对其外侧 remainder 执行 refinement、Cut-cell、small-cell stabilization 与 solver-topology repair。

本阶段未实现 local layer dropping、复杂 termination、Delaunay transition、三维或 overset，也未改动三维核心代码。H4-1 仍明确拒绝严重凹角和超过当前 sharp-corner 阈值的尖锐尾缘。

## 2. 当前实现路径

1. 提取 H4-1 outer envelope，检查方向、合法性、相互关系和计算域包含关系。
2. 自动从 outer-envelope 边长、H4-1 最后一层法向间距和 remainder boundary-level 对应的 Cartesian 尺寸 `h` 推导 progressive transition plan；不依赖样例专用常数。
3. 在 transition 外侧重新建立 remainder quadtree，执行 2:1 balance、Cut-cell 和 small-cell stabilization。
4. 将 BoundaryLayer、transition/remainder Cut polygon 和 remainder Cartesian cell 组装到统一 `TopologyMesh2D`。
5. 对同一不可变 source 重复构建 topology，逐项检查坐标、ID、owner/neighbour 与 connectivity 的确定性。
6. outer-envelope 接口必须是闭合、双 owner、layer/remainder 一一配对、顶点二价的公共分割，并检查接口长度闭合。
7. 检查 layer、transition、remainder 和全域面积守恒。
8. 进入 constrained solver-topology repair；H4-1 layer cell 不允许被修改，transition 可要求保留原始 polygon。
9. 执行 solver-quality gate，并在 solver repair 后再次检查固定 outer-envelope 接口。
10. 只有 solver topology 与 solver quality 均通过时，CLI 才允许输出 extruded OpenFOAM case。

## 3. 最终 solver-ready 验证

2026-08-27 GitHub Actions `cartmesh2d-stage6` 最终 run `33059252115`：完整 native 2D stage CTest `71/71 PASS`。

### Circle

- hybrid cells：700；solver cells：728；H4-1 layer cells：128；
- transition rings：3；最终切向细分：4；
- area error：`-8.88e-15`；
- fixed outer-envelope：32 edges / 32 vertices，single-owner、wrong-pair、non-two-valent 均为 0；
- solver max non-orthogonality：`55.3968°`；
- solver min face weight：`0.08721`；
- solver min volume ratio：`0.02704`；
- 独立 OpenFOAM reader：728 cells、3212 faces、最小体积 `4.14567e-06`、issues 为空。

### Superellipse

- hybrid cells：772；solver cells：791；H4-1 layer cells：72；
- transition rings：4；最终切向细分：8；
- area error：`2.49e-14`；
- fixed outer-envelope：24 edges / 24 vertices，single-owner、wrong-pair、non-two-valent 均为 0；
- solver max non-orthogonality：`63.4869°`；
- solver min face weight：`0.07802`；
- solver min volume ratio：`0.01918`；
- 独立 OpenFOAM reader：791 cells、3552 faces、最小体积 `2.09489e-06`、issues 为空。

两例均满足既有 solver-quality 阈值，没有放宽质量门，也没有绕过失败写出 OpenFOAM case。CI runner 未安装 OpenFOAM，因此明确记录 `checkMesh=UNAVAILABLE`；独立 reader 不能冒充真实 `checkMesh`。

## 4. 验证覆盖

完整 CTest 同时覆盖：H4-1 circle/rectangle、凹角和尖角 fail-closed、H4-2 circle/ellipse/superellipse/参数变化、重复生成确定性、DXF、复杂几何、multi-loop、OpenFOAM writer/reader、quality、small-cell、agglomeration、solver topology 等既有二维回归。

H3 的 100k solver-topology 性能已在独立 H3 阶段完成并有固定报告。H4 最终 CI 不使用参数不同的 NACA workload 冒充 H3 benchmark；generic solver-topology 功能回归由完整 CTest 继续覆盖。

## 5. CI 证据

最终 workflow 上传 `cartmesh2d-final-validation` artifact，包含：

- rectangle / circle / concave / airfoil-like acceptance 的 `.cm2d`、quality/viz JSON 和 SVG；
- circle / superellipse 的 hybrid 与 solver VTK/CM2D/quality/report；
- 两个完整 OpenFOAM case 及独立 reader 报告；
- `checkmesh_status.txt`。

最终 artifact ID `9640997919`；workflow SHA `d15154718be317b0b16286212054f60354d5fb11` 对应的验证全部通过。

## 6. 已知限制

- 当前 H4 仅支持 H4-1 已支持的闭合外流 fixed-layer strip；
- nested wall、local layer dropping、复杂 termination 仍不在本阶段；
- H4-1 对严重凹角 fail-closed；超过 `maxConvexTurnRadians=135°` 的尖锐尾缘也 fail-closed，因此当前 `airfoil_like.xy` 的尖尾缘不能直接用于 H4 boundary layer；
- 当前最终 CI 证明的是 solver-ready topology、项目 solver-quality gate 与 OpenFOAM 格式/守恒独立读取通过；真实 OpenFOAM `checkMesh` 仍需在具备 OpenFOAM 的环境补一轮。
