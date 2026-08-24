# Stage 2D-V1b：曲线与多物体 solver-quality 修复

日期：2026-08-24

## 结论

V1b 在限定产品集上完成：circle、airfoil-like 与 two-obstacles 的最终
OpenFOAM 网格均通过内部 fail-closed 质量门、独立 polyMesh reader，以及
OpenFOAM 2606 的 default、`-allTopology`、`-allGeometry` 三重检查。

这不是任意曲线或 NACA2412 dense 的泛化声明。粗一级网格的真实失败仍被
保留并拒绝；NACA2412 dense 仍留在后续 V1 工作中。

## 方法修正

1. 内部 skewness 改为 OpenFOAM 的定义：skew vector 的归一化尺度至少为
   `0.2*|cell-centre connector|`，二维挤出边的点投影尺度为半边长。旧实现仅
   除以短边长度，会把合法的 2:1 hanging face 严重放大。
2. 新增 `minFaceWeight=0.05` 与 `minVolumeRatio=0.01` 两个不可绕过的质量门。
   插值权重严格使用 OpenFOAM 的面法向投影距离公式；二维面积在等厚挤出后
   与三维 cell volume ratio 完全等价。
3. 非凸/共线 solver cell 的确定性分割先搜索面积更均衡的两片凸分割，再做
   ear clipping；人工对角线仅在精确面积守恒且凸时移除。
4. 对曲线和多物体产品增加一级最大局部细化：circle 5→6、airfoil-like
   6→7、two-obstacles 7→8。粗一级产品继续 fail-closed，不修改阈值。

OpenFOAM 参考实现：

- `faceWeight = min(dOwn,dNei)/(dOwn+dNei+vSmall)`；
- `volRatio = min(volOwn,volNei)/(max(volOwn,volNei)+vSmall)`；
- 默认门分别为 0.05 与 0.01。

源码依据：

- https://cpp.openfoam.org/v13/polyMeshCheck_8C_source.html#l00173
- https://cpp.openfoam.org/v13/polyMeshCheck_8C_source.html#l00231
- https://cpp.openfoam.org/v13/polyMeshCheck_8H_source.html#l00120

## 保留的失败基线

| 产品 | 粗网格 | 外部 `-allGeometry` 失败 |
|---|---:|---|
| circle | level 5 | 16 个低权重面，8 个低体积比面；min weight 0.0419536，min ratio 0.00789712 |
| airfoil-like | level 6 | 1 个低权重面；min weight 0.0455434 |
| two-obstacles | level 7 | 30 个低权重面；min weight 0.0296401 |

内部质量门现在会在写 OpenFOAM 前复现并拒绝这些失败，不再出现“内部 PASS、
外部 allGeometry FAIL”的假阳性。

## 最终产品证据

| 产品 | level | cells | max non-ortho | OpenFOAM max skew | min determinant | min face weight | min vol ratio |
|---|---:|---:|---:|---:|---:|---:|---:|
| circle | 6 | 552 | 68.5611 | 3.23491 | 0.00937523 | 0.0527138 | 0.0431993 |
| airfoil-like | 7 | 835 | 63.8562 | 3.84263 | 0.00299446 | 0.0581767 | 0.0253862 |
| two-obstacles | 8 | 2104 | 40.8669 | 2.04396 | 0.0191647 | 0.0517241 | 0.0323232 |

三类的 default、`-allTopology`、`-allGeometry` 日志均包含无条件 `Mesh OK.`，
且不含 `Failed`、`FOAM FATAL` 或 stack trace。内部 face weight、volume ratio、
non-orthogonality 与外部日志数值一致。三类重复生成的 points/faces/owner/
neighbour/boundary、求解字典和 solver-quality JSON 均逐字节一致。

## 可复现命令

```bash
cmake -S cartmesh2d -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build -j 8
ctest --test-dir build --output-on-failure -j 8
```

Release 与 ASan/UBSan 均为 40/40 PASS。三个新增 OpenFOAM e2e 测试还分别调用
`tools/verification/check_openfoam2d.py` 独立检查 points/faces/owner/neighbour、
边界 patch、二维 extrusion 与连通域。

真实 OpenFOAM 外部验收示例：

```bash
docker run --rm -v "$CASE:/case" opencfd/openfoam-run:2606 \
  bash -lc 'cd /case && checkMesh && checkMesh -allTopology && checkMesh -allGeometry'
```

## 下一边界

V1b 只解决 solver-quality 指标闭环和上述三个产品。下一步 V1c 才进入自由流
保持、制造解、离散守恒和网格收敛；不得以 `checkMesh` 代替方程求解验证。
