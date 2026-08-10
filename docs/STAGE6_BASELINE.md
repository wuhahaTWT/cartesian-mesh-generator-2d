# Stage 6.0 基线冻结与 Stage 6.1 开工准备

日期：2026-08-11（Asia/Shanghai）

状态：**Stage 6.0 PASS；Stage 6 仍未完成；Stage 6.1 尚未开始实现。**

机器可读基线是 `artifacts/stage6_baseline.json`。本记录只冻结当前事实并准备
Stage 6.1 的开工契约，没有修改 C++、测试、CMake、构建脚本或算法。

## 1. 结论

当前仓库存在三条不同能力边界：

| 路径 | Cut-cell 几何/拓扑 | 完整 solver mesh | 本轮独立读取 | 当前结论 |
|---|---:|---:|---:|---|
| uniform | 有 | 有 | meshio PASS | 能写完整 `polyMesh`，但本轮 Docker 未运行，未重跑 `checkMesh` |
| adaptive | 有，含真实 coarse-fine 连接 | 无 | meshio PASS | CLI 明确拒绝 OpenFOAM，这是 6.1 唯一主缺口 |
| incremental | 有，且与全量等价 | 无 | meshio PASS | 稳定 ID、复用和映射有效，但没有 solver output 链路 |

因此不能把“adaptive 内部拓扑存在”写成“adaptive 已可供 OpenFOAM 使用”，也不能把
uniform 的 solver writer 或旧 Stage 6 紧凑 writer 冒充 adaptive writer。

Stage 6 历史终态仍为 `blocked_openfoam_quality`：OpenFOAM 2606 曾报告 3,293 张错误
face pyramid、12,078 张高偏斜面和最大 skewness 1274.65，没有 `Mesh OK`。本轮没有
修改或重跑该历史大规模产物。

## 2. 接管点和环境

- 分支：`main`
- commit：`c4d72aa3cca92b7d409a1ff30ef2df23645fa8a4c`
- 开工前工作区：clean
- tag / remote：均无
- `CHANGELOG.md`：交接清单提到，但仓库实际不存在；本阶段只记录冲突，不擅自补文件
- 机器：MacBookAir10,1，Apple M1，8 核，8 GB
- 系统：macOS 26.5.2 build 25F84，arm64
- 工具链：CMake 4.4.2，GCC 15.2.0，Release
- 正式案例运行线程：1
- 开工可用磁盘：约 27.3 GiB
- Docker daemon：未运行，因此没有本轮 OpenFOAM `checkMesh`

历史 OpenFOAM JSON 只作为历史证据引用，不记作本轮重验。

## 3. Release 构建与测试

执行：

```sh
cmake --preset release
cmake --build --preset release --parallel 4
ctest --preset release --output-on-failure
```

结果：21/21 CTest PASS，总测试墙钟 5.75 s。该结果覆盖 Stage 0–6 工程测试，但工程测试
不替代独立 reader 或 OpenFOAM 质量检查。

## 4. 五类固定输入

| 案例 | SHA-256 |
|---|---|
| `tests/data/closed_unit_cube_ascii.stl` | `cac46acc501ee3256e777cefab0120251e6e9bbc6f8a02012ca3886f452bfd5b` |
| `tests/data/nonconvex_l_prism_ascii.stl` | `8d2268318d3548e5be1271a0c9ad5c8d5f7e76e9cafde94f07ed3e7988e1c7dd` |
| `benchmarks/analytic/stage4/thin_shell_wall005_ascii.stl` | `e5254296ff9fd474e747262ee0b17b429b51e993c379001224684d3834fafb18` |
| `benchmarks/analytic/stage4/two_disjoint_cubes_ascii.stl` | `85173ae7ff4034776a9f93011f84ba5003efc64a425f70c176fe110303f3f8cb` |
| `benchmarks/public/stage4/stanford_bunny_libigl_binary.stl` | `21548a3c21187110f25d442f9ec688e90efd048c83590b422b5e4c8f7ef2cf35` |

当前输入契约仍是封闭、流形、定向一致、无自交 STL。脏几何默认拒绝并定位，不自动修复。

## 5. 本轮小规模回归

| 案例 | 路径/规模 | 核心计数 | 墙钟 | 峰值 RSS | 结果 |
|---|---|---|---:|---:|---|
| cube | uniform R8 | 512 background / 296 cut | 0.155237 s | 8,798,208 B | 完整 `polyMesh` 写出；meshio PASS |
| L prism | adaptive L2–L4 | 2,234 leaves / 684 cut | 0.226649 s | 11,075,584 B | 几何/拓扑 PASS；meshio PASS；无 solver writer |
| thin shell | uniform R12 | 1,728 background / 488 cut / 2 regions | 0.945915 s | 13,926,400 B | 完整 `polyMesh` 写出；meshio PASS |
| two shells | adaptive L2–L4 | 1,800 leaves / 400 cut | 0.159075 s | 8,798,208 B | 几何/拓扑 PASS；meshio PASS；无 solver writer |
| Bunny | uniform classification R8 | 512 cells / 6,966 triangles | 0.052043 s | 8,781,824 B | 分类 PASS；冲突 0；不是 Cut-cell |

所有正式案例为 Release、单线程。机器精确字段见 JSON。

### Bunny 止损事实

低分辨率 Bunny 并不自动意味着复杂 Cut-cell 路径便宜：

- 本轮 full Cut-cell R8 未形成完整输出，止损清理时终止；没有可用正式计时；
- Stage 6 compact R12 在 117.63 s 时主动中断，当时 user/sys 为 116.12/0.91 s；
- 没有因此扩大到 R24、R48、R96、R216 或千万级；
- 本轮只把 R8 精确 triangle-AABB 分类作为新鲜 smoke；历史 Stage 4 Bunny 完整验证仍单独标为历史证据。

这项结果是性能债和后续剖析入口，不是 Stage 6.1 扩大范围的理由。6.1 首先只处理解析
cube、非凸 L 和 coarse-fine + cut 最小案例。

## 6. 独立验证与确定性

### 独立读取

- cube uniform：`verify_stage3_cutcell.py`，meshio 5.3.5，PASS；
- adaptive L、adaptive 多壳层、uniform thin shell：`meshio_stage4_verify.py`，PASS；
- incremental 局部轮廓：`meshio_stage5_verify.py`，PASS；
- 旧 Stage 3 checker 读取 adaptive 报告时因缺少 uniform-only `dimensions` 字段失败；
  Stage 4 checker 可以读取并验证 adaptive。该工具兼容性边界不能写成算法失败，也不能
  写成“所有独立检查器均支持 adaptive”。

### 重复运行

- uniform cube 在相同 solver-output 配置下结果 hash 两次均为
  `e6fd5337f55140ea`；五个 `polyMesh` 核心文件 SHA-256 逐项相同；
- adaptive L 结果 hash 两次均为 `a47fd5b59db2b115`；
- incremental 与 full rebuild 两次均为 `40dc680ddf6d7a3f`；
- incremental VTU SHA-256 两次均为
  `126ba658d42bd4a405588d5863515c330549a2ce1b858287b43905d7c4e0e05d`；
- mapping SHA-256 两次均为
  `b30ad394f0c49f8dab25a5ee3a1eb1ce92cbb3f44ecf24fe63d19190808b175a`。

cube 不带 `--openfoam-case` 时会选择 `convex_halfspace`，带 solver 输出时当前代码有意选择
`oriented_tetrahedral_chain`。不同配置 hash 不同不是非确定性；正式比较使用完全相同配置。

## 7. Adaptive → OpenFOAM 的准确断点

本轮命令级复核：

```text
cartmesh_cutcell_cli --adaptive ... --openfoam-case ...
错误：当前 OpenFOAM 完整体网格输出只支持均匀 Cartesian 背景
exit code = 1
```

断点不在 STL 诊断、八叉树细化、Cut-cell 几何或粗细邻接。这些能力已经存在：

- `LinearOctree` 提供稳定 Morton 叶顺序、`cell_bounds()` 和 `face_neighbors()`；
- `build_adaptive_triangulated_cut_cell_mesh()` 已为 coarse-fine 接口汇总面积和一阶矩；
- `ConvexCutCellMesh` 已保存 fluid cells、显式凸片、embedded faces、内部连接和 region；
- Stage 3 回归明确统计了实际跨层 connection。

真正断点是 `write_openfoam_poly_mesh()` 的接口和遍历仍绑定
`UniformCartesianGrid`：它按 `cell_key()`、`nx/ny/nz`、统一 `spacing()` 和单一同层邻居
构造 solver faces。CLI 因此在调用 writer 前主动拒绝 `LinearOctree`。

## 8. Stage 6.1 可直接实施的设计

### 8.1 目标和停止条件

只解决：

> 将现有 adaptive `LinearOctree + ConvexCutCellMesh` 写成一套确定性的完整 OpenFOAM
> `polyMesh`，不同时做质量聚合、千万级优化或 Stage 7。

通过条件：解析 adaptive cube 与非凸 adaptive L 都完成写出，独立完整 reader PASS，
OpenFOAM `checkMesh -allTopology` 明确输出 `Mesh OK.`。未达到时不进入 6.2。

### 8.2 必须复用

- `include/cartmesh/grid/LinearOctree.hpp`
- `include/cartmesh/grid/OctreeNodeCode.hpp`
- `include/cartmesh/cutcell/ConvexCutCellMesh.hpp`
- `include/cartmesh/cutcell/TriangulatedSurfaceCutter.hpp`
- `src/cutcell/TriangulatedSurfaceCutCellMesh.cpp` 的 adaptive coarse-fine 拓扑
- `src/io/OpenFoamWriter.cpp` 已验证的 component-cell、凸多边形求交、点焊接、patch 排序和写出逻辑

不得新建第二套 geometry、topology 或 OpenFOAM writer。

### 8.3 应泛化或替换的 uniform 假设

拟修改文件：

- `include/cartmesh/io/OpenFoamWriter.hpp`
- `src/io/OpenFoamWriter.cpp`
- `apps/cutcell_cli/main.cpp`
- `tests/stage3_test.cpp` 或单独的 Stage 6.1 最小测试文件
- `CMakeLists.txt` 仅在确需登记新的测试文件时修改

修改方向：

1. 为现有 reference writer 增加 `LinearOctree` 输入，或抽出供 uniform/octree 共用的
   background-cell 访问适配层；不能复制整个 writer。
2. 用 `tree.cell_bounds(tree.leaf_code(leaf_id))` 替代 uniform `cell_key()` 几何派生。
3. 用 `tree.face_neighbors()` 枚举同层与 coarse-fine 邻居；每个物理接口只处理一次。
4. coarse face 与 1–4 个 fine face 的显式凸片多边形做公共细分，覆盖面积和一阶矩必须守恒。
5. 远场判断基于叶 AABB 与 domain 边界，不依赖 `i/j/k == 0/n-1`。
6. 容差按参与叶的局部尺度计算，不使用单一全局 uniform spacing。
7. solver cell ID 按稳定 Morton 叶顺序和局部 component 顺序生成；报告
   `OctreeNodeCode → solver cell ID` 映射，不能把临时 leaf index 说成长期稳定 ID。
8. 删除 CLI 中 adaptive OpenFOAM 的主动拒绝，仅在新 writer 及门禁测试存在后删除；
   不先删保护再补实现。

### 8.4 关于“有加有减”的具体纪律

6.1 不以行数正负相等为目标，但必须满足：

- 新增 octree 适配逻辑时，同时移除或收敛被替代的 uniform-only 分支；
- 共用面组装、排序和序列化必须抽取复用，不能复制成 `AdaptiveOpenFoamWriter.cpp`；
- CLI 的拒绝逻辑只有在新路径通过测试后删除；
- 不删除旧 uniform 回归；新的抽象必须证明 uniform `polyMesh` 五文件 hash 或拓扑语义不回归；
- `ScalableOpenFoamWriter` 当前继续冻结，它的统一属于 6.5，不在 6.1 又造第三套适配路径。

因此合理 diff 很可能既有 `+` 也有 `-`，但判断依据是旧责任是否被新抽象真正替换，而不是
机械追求行数对称。

### 8.5 最小失败案例和验收顺序

1. **coarse-fine 无切割**：一个 coarse face 对四个 fine faces，验证唯一 owner/neighbour、
   无重叠、无缺口、总面积与一阶矩守恒。
2. **coarse-fine + cut**：切面穿过层级接口，保留最小 STL 和固定八叉树。
3. **adaptive cube**：必须实际存在跨层连接，完整 `polyMesh` 独立读取。
4. **adaptive L prism**：验证非凸表面和 patch 连续性。
5. 重复两次：stable leaf mapping、topology fingerprint、五个 `polyMesh` 文件稳定。

每一步先内部不变量，再独立 reader，最后 OpenFOAM。R24 Bunny 和更高分辨率都不属于 6.1。

### 8.6 明确禁止

- `--convex-piece-cells` / `convex_piece_exact`
- writer 内 2–N 搜索聚合
- `quality_conformal_split`
- `quality_constrained_cell_merge`
- kernel tetra repair
- 删除小流体片、把缺口改成 wall、翻孤立面或放宽 `checkMesh`
- 在 6.1 同时改 `ScalableOpenFoamWriter` 追求大规模

若 reference writer 无法在不复制系统的前提下泛化，应停止并提交设计审查，不直接大重写。

## 9. 6.1 开工前仍需用户确认的内容

开始修改代码前必须再次报告并等待确认：

- 当前 commit 和 clean status；
- 具体拟修改文件；
- reference writer 的适配接口；
- 三个最小失败案例；
- uniform 回归和 hash/语义等价方式；
- Docker/OpenFOAM 是否已可用；
- 发生大范围数据结构重写时的停止条件。

本文件不构成自动进入 6.1 的授权。

## 10. Stage 6.0 关闭声明

Stage 6.0 的通过含义仅为：真实基线、固定输入、三条路径能力边界、当前性能债、确定性证据
和 6.1 开工方案已经形成。它不表示 adaptive solver output、OpenFOAM 质量或 Stage 6 已通过。
