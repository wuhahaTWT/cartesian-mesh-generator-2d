# 2D-R1 鲁棒 Cut-cell 事件拓扑验证

日期：2026-08-24

## 1. 本阶段边界

本阶段只修复原生二维几何谓词与单叶 Cut-cell 局部事件拓扑，不进入：

- P1 空间索引、线性四叉树或性能优化；
- D1 多环、多物体与嵌套域；
- S1 solver-ready 面几何、质量阈值或求解器格式；
- V1 CFD 求解、守恒与网格收敛验证。

因此本记录中的 `valid=true` 只证明现有拓扑审计项为零，不代表网格已经满足 CFD 质量门。

## 2. 实现事实

- `orientationSign()` 先使用 binary64 误差界快速判定；不确定时把六个输入坐标精确分解成
  dyadic，再用定长多字整数完成精确行列式符号计算。
- 点在线段、线段相交、点在多边形等容差改为局部长度尺度，修复叉积与长度直接比较的量纲错误，
  并避免大坐标平移无故放大局部容差。
- 多边形面积和质心以首顶点为局部原点累加，降低大坐标平移造成的抵消。
- 局部 Cut-cell 从“每个顶点必须一进一出”改为确定性 planar half-edge rotation system；
  扇区排序使用同一精确 orientation 谓词，不使用 `atan2` 决定拓扑。
- 每个局部流体分量只接收实际位于其边界上的 embedded fragments。
- 局部带孔情形仍显式失败；没有把孔洞静默填满。
- NaN/Inf 坐标和非有限 Cartesian box 被显式拒绝。

## 3. 最小失败回归

保留并验证：

- 固体角点与 Cartesian 角点精确重合：外侧单元为 full fluid，内侧单元为 empty；
- 切触事件不产生伪 Cut-cell；
- 同一毫米级 Cut-cell 平移到 `1e9` 坐标附近后仍保持分类和面积；
- 普通 binary64 行列式因乘积抵消得到 0、但精确结果为 `-1` 的 orientation 案例；
- 最小 subnormal 输入下普通乘积下溢、精确 orientation 仍为正的案例；
- rectangle、airfoil-like、NACA2412 dense 与 superellipse 四个历史端到端失败案例。

## 4. Release 验证

```sh
cmake -S . -B /private/tmp/cartmesh2d-r1-official-build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCARTMESH_BUILD_2D=ON \
  -DCARTMESH_BUILD_TESTS=ON \
  -DCARTMESH_BUILD_BENCHMARKS=OFF
cmake --build /private/tmp/cartmesh2d-r1-official-build -j 6
ctest --test-dir /private/tmp/cartmesh2d-r1-official-build --output-on-failure -j 1
```

结果：全仓 `45/45` 通过；二维部分 `21/21` 通过。修复前二维基线是 `17/21`，失败项为
rectangle、airfoil-like、NACA2412 dense、superellipse。

## 5. ASan / UBSan

```sh
cmake -S . -B /private/tmp/cartmesh2d-r1-sanitize \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCARTMESH_BUILD_2D=ON \
  -DCARTMESH_BUILD_TESTS=ON \
  -DCMAKE_CXX_FLAGS='-fsanitize=address,undefined -fno-omit-frame-pointer' \
  -DCMAKE_EXE_LINKER_FLAGS='-fsanitize=address,undefined'
ctest --test-dir /private/tmp/cartmesh2d-r1-sanitize \
  -R '^cartmesh2d_' --output-on-failure -j 1
```

结果：二维 `21/21` 通过，无 ASan/UBSan 报告。Debug+sanitizer 下复杂算例耗时包括：
gear-star `149.25 s`、serpentine `22.40 s`、nozzle `39.51 s`、NACA2412 `45.70 s`、
superellipse `7.38 s`。这些时间不是性能成绩，只说明 P1 性能整改仍然必要。

## 6. 确定性和外部读取

四个修复算例分别独立生成两次，两个输出目录经 `diff -rq` 比较完全一致；比较范围包括
`.vtk`、`.cm2d`、`.quality.json` 和 `.viz.json`。

使用项目外部的 `meshio 5.3.5` 读取 VTK polygon connectivity：

| 算例 | points | polygon cells |
|---|---:|---:|
| rectangle | 340 | 248 |
| airfoil-like | 473 | 335 |
| NACA2412 dense | 1710 | 1194 |
| superellipse | 1183 | 840 |

四个输出均可独立读取，且内部拓扑审计均为：重复顶点、重复边、孤立内部边、非流形边、
未分类边界边、开放单元环、面积不一致全部为 0。

## 7. 修复后网格及诚实质量边界

产物根目录：`/private/tmp/cartmesh2d-r1-artifacts/run3/`。

| 算例 | cells | source small cells | min area fraction | max edge ratio | max reported skewness |
|---|---:|---:|---:|---:|---:|
| rectangle | 248 | 0 | 0.333333 | 3.00 | 0.1225 |
| airfoil-like | 335 | 14 | 0.000220 | 162.43 | 0.5968 |
| NACA2412 dense | 1194 | 10 | 0.00000427 | 2944.52 | 1.3170 |
| superellipse | 840 | 16 | 0.0000518 | 2992020.29 | 1.5975 |

NACA2412 和 superellipse 的质量数值仍明显不可接受；尤其 superellipse 最短边只有
`8.9257e-9`。2D-R1 只关闭了几何退化和非流形阻断，绝不能据此宣称 solver-ready。

## 8. 阶段结论

2D-R1 的精确 orientation、局部尺度容差、half-edge 事件拓扑、失败回归、确定性、
sanitizer 和外部连接性读取门已经通过。下一阶段只能进入 2D-P1 快速核心；质量阈值、
solver-ready 数据和 CFD 验证仍分别属于 S1/V1，不能提前宣称完成。
