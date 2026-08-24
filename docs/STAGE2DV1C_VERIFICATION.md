# Stage 2D-V1c：守恒、自由流与制造解收敛验证

日期：2026-08-24
范围：只修改独立二维项目；不进入圆柱/翼型物理基准 V1d。

## 结论

V1c 已通过。这里的“通过”只表示导出的真实 OpenFOAM Cut-cell 产品满足几何守恒、
常数场保持、线性场误差门和调和制造解的实测网格收敛；它不等于 Navier–Stokes
物理验证，也不代替 V1d 的圆柱/翼型基准。

## 本阶段发现并修正的验证路线缺陷

原来的 Quadtree 只把相交边界的叶子细化到 `maxLevel`。因此将圆形案例从 level 6
提高到 level 7/8 时，远场仍然很粗，不能作为全局 PDE 收敛序列。首轮边界单独细化
的二次调和场 L2 仅为 `4.358e-2 -> 3.860e-2 -> 3.331e-2`，且 level 8、
`small-alpha=0.05` 因 8 个低于 0.05 的插值权重面被质量门拒绝。

没有降低质量门。新增 `QuadtreeRefinementPolicy2D::minimumLevel`，默认 0 保持旧行为；
验证序列使用 `(minimumLevel,boundaryLevel)=(4,6),(5,7),(6,8)`。最细网格使用
`small-alpha=0.20` 合并不可信小单元，最终最小 face weight `0.116472`、最小相邻
体积比 `0.055774`，均高于既有门槛。

## 方法依据

- Salari 与 Knupp 的 MMS 报告说明用解析场、相应边界/源项和观测精度阶进行代码验证：
  https://www.osti.gov/biblio/759450/
- OpenFOAM 的 `corrected` snGrad 包含显式非正交修正项：
  https://doc.openfoam.com/2212/tools/processing/numerics/schemes/sngrad/rtm/corrected/

本案例选取无源调和解

```text
T_exact = x^2 - y^2 + 0.15 x - 0.05 y
laplacian(T_exact) = 0
```

边界 Dirichlet 值在 OpenFOAM 自己写出的真实边界面中心 `Cx/Cy` 上计算；误差在真实
单元中心计算，并以 OpenFOAM `writeCellVolumes` 的体积加权。`corrected` 的非正交项
是显式项，因此每个案例执行 30 次稳态外迭代；只做一次时即使 PCG 残差很小，也不
代表修正后的离散场已经收敛。

## 产品与结果

| minimum / boundary level | solver cells | L1 | L2 | Linf |
|---|---:|---:|---:|---:|
| 4 / 6 | 580 | 6.0469e-3 | 7.2091e-3 | 2.1439e-2 |
| 5 / 7 | 1560 | 1.9500e-3 | 2.3161e-3 | 9.4521e-3 |
| 6 / 8 | 4312 | 8.8929e-4 | 1.1137e-3 | 4.9089e-3 |

每次 `h` 减半的观测阶：

| norm | coarse -> medium | medium -> fine |
|---|---:|---:|
| L1 | 1.6327 | 1.1327 |
| L2 | 1.6381 | 1.0564 |
| Linf | 1.1815 | 0.9452 |

门槛为三种范数严格下降且每段观测阶不低于 0.9。当前结果支持至少近一阶收敛；不把
它夸大成二阶。误差主要位于圆形 Cut-cell 环和 2:1 过渡带。

同一套 4312-cell 最细网格：

- 常数场：L2 `1.6282e-15`，Linf `5.7732e-15`；
- 线性调和场：L2 `2.4852e-4`，Linf `1.7859e-3`；
- 五份 solver 日志最终归一化残差最大 `8.9232e-14`；
- 独立 polyMesh 读取器：最大单元闭合残差 `2.7105e-20`，物理边界常向量总通量
  `-2.9490e-18`；
- OpenFOAM 2606：三套网格分别执行 default / `-allTopology` / `-allGeometry`，
  9/9 均为无附加失败的 `Mesh OK.`。
- 最细产品重复生成的 `points/faces/owner/neighbour/boundary/solver_quality.json`
  逐字节相同；对应 SHA-256 前五项依次为 `dc8a1675...`、`3942749e...`、
  `0f261ccf...`、`113f5a50...`、`3d02538e...`。
- Release 40/40 与 ASan+UBSan 40/40 均通过。

## 可重复命令

```bash
cmake -S cartmesh2d -B build-v1c -DCMAKE_BUILD_TYPE=Release
cmake --build build-v1c -j4
ctest --test-dir build-v1c --output-on-failure

build-v1c/cartmesh2d_cli cartmesh2d/examples/acceptance/circle.xy out/circle-m4-M6 \
  6 0.25 0.10 exterior out/circle-m4-M6-case 4

docker run --rm -v "$CASE:/case" opencfd/openfoam-run:2606 \
  bash -lc 'cd /case && postProcess -func writeCellCentres && postProcess -func writeCellVolumes'

python3 cartmesh2d/tools/verification/openfoam_harmonic_mms.py \
  configure "$CASE" --field quadratic
docker run --rm -v "$CASE:/case" opencfd/openfoam-run:2606 \
  bash -lc 'cd /case && laplacianFoam > laplacian-quadratic.log'
python3 cartmesh2d/tools/verification/openfoam_harmonic_mms.py \
  evaluate "$CASE" --field quadratic --log "$CASE/laplacian-quadratic.log"
```

`tools/verification/check_openfoam_v1c.py` 对 MMS 阶次、常数/线性场、独立几何守恒、
solver 残差和九份 `checkMesh` 日志统一 fail-closed。可视化由
`tools/visualization/render_openfoam_mms.py` 直接读取最终 `polyMesh` 和 `T`，不是背景
网格示意图。

## 保留的失败边界

以下命令仍应被质量门拒绝，作为更高层边界细化会产生低权重面的已知回归输入：

```text
cartmesh2d_cli circle.xy out/circle-l8 8 0.25 0.05 exterior out/circle-l8-case
```

本阶段用更强的小单元合并取得合格产品，没有降低 face-weight、volume-ratio 或
`checkMesh` 标准。
