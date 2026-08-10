# 阶段 6 验证：千万级完整网格与产品化门禁

日期：2026-08-10（Asia/Shanghai）

状态：**已完成规模、完整导出、独立读取、确定性和资源验证；OpenFOAM 默认质量门未通过，阶段 6 尚未关闭。**

## 1. 终态结论

本轮已经真实生成并保存 216³、即 10,077,696 个背景 Cartesian 单元的 Stanford Bunny
案例。完整 binary OpenFOAM `polyMesh` 含 8,700,174 个实际流体控制体、26,463,750 张面和
9,053,073 个点，不是八叉树叶数、体素标签、截图或只含 Cut-cell 的局部文件。

独立 Python 二进制读取器完整读取约 913 MiB 的 `polyMesh`，检查所有索引、patch 连续性并为
8,700,174 个控制体逐边重建关联，非二关联边数为 0。相同输入第二次生成的计数、紧凑 hash、
拓扑 hash 和五个 OpenFOAM 文件 SHA-256 完全一致。生成、导出、读取、内存和磁盘均在计划预算内。

阶段 6 仍不能声明完成，原因是 OpenFOAM 2606 `checkMesh -constant -allTopology` 虽然完整读取
binary 网格，并通过 boundary definition、cell-to-face addressing、upper-triangular ordering、
face vertices、topological cell zip-up、单区域、patch 闭合、正面面积、正单元体积和非正交检查，
但没有输出 `Mesh OK`，仍失败两项默认质量检查：

- 3,293 张面存在错误的 face pyramid 方向；
- 12,078 张面高度偏斜，最大 skewness 为 1274.65。

因此 `artifacts/stage6_acceptance.json` 的真实状态是
`blocked_openfoam_quality`、`stage6Complete=false`、`solverReadyCutCellMesh=false`。
阶段 7 未启动。

## 2. 阶段 6 新增实现

核心实现：

- `CompactUniformCutCellMesh`：普通流体/固体背景单元用打包状态表示，只为表面邻域保留显式
  Cut-cell 几何；
- 稀疏三角面栅格化和非切割区域连通分类，不为一千万个单元逐一扫描完整表面；
- 稳定背景 ID 到求解器单元 ID 的确定性映射；
- OpenFOAM 2606、32 位 label、64 位 scalar、本机小端 binary 的完整流式写出；
- 普通 Cartesian 区域隐式遍历，只有 Cut-cell 邻域的公共细分面常驻内存；
- 明确记录数值封口、拓扑封口、覆盖误差、退化面、小 Cut-cell 和全部不变量；
- `--preview FILE.vtu` 从同一紧凑状态生成局部背景预览，不调用旧阶段的完整显式多面体路径；
- 独立 mmap 二进制读取器和禁网 OpenFOAM Docker 检查器。

主要文件：

- `include/cartmesh/scalable/CompactUniformCutCellMesh.hpp`
- `src/scalable/CompactUniformCutCellMesh.cpp`
- `include/cartmesh/io/ScalableOpenFoamWriter.hpp`
- `src/io/ScalableOpenFoamWriter.cpp`
- `apps/stage6_cli/main.cpp`
- `tests/stage6_test.cpp`
- `tools/stage6_binary_polymesh_verify.py`
- `tools/stage6_openfoam_verify.py`

## 3. 千万级完整产物

输入：`benchmarks/public/stage4/stanford_bunny_libigl_binary.stl`

输入 FNV-1a 64：`c8465360de458d34`，三角面数 6,966，分辨率 216³。

| 项目 | 实测值 |
|---|---:|
| 背景 Cartesian 单元 | 10,077,696 |
| 全流体背景单元 | 8,699,991 |
| 全固体背景单元 | 1,377,705 |
| 显式表面单元 | 129,729 |
| Cut-cell | 129,724 |
| 完整 OpenFOAM 控制体 | 8,700,174 |
| 点 | 9,053,073 |
| 总面 | 26,463,750 |
| 内部面 | 25,957,283 |
| 边界面 | 506,467 |
| farfield 面 | 279,936 |
| embedded_wall 面 | 226,531 |
| 面顶点引用 | 106,025,897 |
| 完整 polyMesh 字节数 | 956,918,167 |
| 全局流体区域 | 1 |
| 总流体体积 | 0.0045167391957595402 |
| 嵌入壁面面积 | 0.05821291897804135 |

正式文件：

- `artifacts/stage6_10m_case/constant/polyMesh/`
- `artifacts/stage6_10m_export.json`
- `artifacts/stage6_10m_external_reader.json`
- `artifacts/stage6_10m_determinism.json`

生成器内部不变量均为零：非闭合单元、负体积、共享面不匹配、分类冲突、未完成分量分析、
直接流体/固体面、数值丢弃片和导出覆盖放宽。最小 Cut-cell 体积分数为
`2.5998484011974675e-12`，体积分数小于 0.01 的 Cut-cell 有 12,597 个；它们只统计，
没有为通过质量检查而删除。

## 4. 显式记录的数值拓扑动作

阶段 6 不隐藏几何封口：

- 紧凑几何在 9 张 Cartesian 面上做了数值封口，总面积
  `4.175463553366974e-14`；
- OpenFOAM 导出器发现并封闭 29 个简单缺口环，原始缺口边 91 条；
- 封口环总面积 `1.3397406013766621e-06`，最大单环面积
  `2.0136861291078814e-07`；
- 拓扑归一化后被压塌并删除的正面积面为 0；
- 背景接口和嵌入 patch 覆盖放宽数均为 0。

独立全量读取器随后证明所有控制体的每条边均恰好被两张面引用，说明最终写出文件在
cell-edge 层面闭合；这不会把上述封口动作改写成“原始几何没有缺口”。

## 5. 独立 binary 读取

`tools/stage6_binary_polymesh_verify.py` 不链接 cartmesh。它通过 mmap 直接解析 binary
`points`、`faceCompactList`、`owner`、`neighbour` 和 ASCII `boundary`，执行：

- 列表长度、截断和单调 offset 检查；
- 点坐标有限性、面至少三点、点/单元索引范围检查；
- owner 排序和 neighbour 计数检查；
- boundary patch 连续覆盖检查；
- 全部 8,700,174 个控制体的逐边二关联检查；
- 五个正式文件的 SHA-256。

结果：`status=pass`、失败列表为空、非二关联单元边数 0、最大关联数 2。

当前 SHA-256：

| 文件 | SHA-256 |
|---|---|
| points | `ebd906a6b6330bed751bafd0690a02c053f51ddb1a3efe502e2a83fd11480064` |
| faces | `526cd299396f27701d6ccc6287b25dc355a2da7138a8c7a1dff86579ae156799` |
| owner | `4b5872ee8b651e6364847d601478be2d42ef47b970b26a6da17a6d799cc7e961` |
| neighbour | `96c9984fe397faffee4d173f0f0d873767db257b063f983d97ae29c074055b97` |
| boundary | `732a13f9f4a39868adcdb11fd6fb175ad5a2efe723af54f61d00707cde7d9d5c` |

## 6. OpenFOAM 2606 外部检查

中等规模 R96 Bunny 产物含 770,624 个控制体和 2,427,392 张面。检查器使用本机已有的
`opencfd/openfoam-run:2606`，Docker 网络关闭，并把源 case 复制到临时目录，避免
`checkMesh` 写集合时修改源证据。

通过项：

- binary 格式完整读取；
- boundary definition；
- cell-to-face addressing；
- upper-triangular ordering；
- face vertices；
- topological cell zip-up；
- 单一流体区域；
- farfield 与 embedded_wall 均为 closed singly connected；
- boundary openness、cell openness、face area、cell volume；
- 最大非正交 62.1418，OpenFOAM 判定 OK；
- 未使用点已从 2 修复为 0。

失败项：

- `3293 faces are incorrectly oriented`；
- `12078 highly skew faces`，最大 skewness `1274.65`；
- 最终为 `Failed 2 mesh checks.`，没有 `Mesh OK.`。

另外保留但不计入这两项默认失败的拓扑警告：1 张 duplicate baffle face、2 张
non-standard edge connectivity face、1 个至多一张内部面的单元和 47 个只有两张内部面的单元。

证据：

- `artifacts/stage6_medium_checkmesh.json`
- `artifacts/stage6_medium_checkmesh.log`
- `artifacts/stage6_medium_bunny_r96_reader_diagnostic.json`

组件级控制体可能是非星形的，这是错误 face pyramid 的主要结构性原因。逐凸片写出、
写出器内 2–6 元局部聚合和内核拆分曾作为 R24 实验路线：独立 reader 最终能达到
非二关联边为 0，但 OpenFOAM 始终没有 `Mesh OK.`，并产生极小控制体、数千张高偏斜面、
重复 baffle 和不可扩展的组合搜索。

该分支已从正式代码、CLI 和测试中删除。完整止损复盘和压缩后的最小失败证据为：

- `docs/STAGE6_FAILED_BRANCH_HANDOFF_CN.md`
- `artifacts/stage6_abandoned_writer_repair_summary.json`

后续不得在 OpenFOAM 写出器中恢复这条路线，也不能通过翻转面、删除小单元或放宽
`checkMesh` 阈值伪造成功。

## 7. 确定性

同一输入、216³、Release、单线程运行两次：

- 背景、流体、固体、显式、Cut-cell、控制体、点、面和 patch 计数完全一致；
- 紧凑 hash 均为 `8e1a7a217f05c555`；
- OpenFOAM 拓扑 hash 均为 `bc146d81ae8f6dc1`；
- 五个 `polyMesh` 文件 SHA-256 完全一致。

证据：

- `benchmarks/baselines/stage6_10m_run1_2026-08-10.json`
- `benchmarks/baselines/stage6_10m_run2_2026-08-10.json`
- `artifacts/stage6_10m_external_reader_repeat.json`
- `artifacts/stage6_10m_determinism.json`

## 8. 性能与资源

硬件：MacBookAir10,1，8,589,934,592 字节内存，8 个逻辑 CPU。构建：Release，
GCC 15.2.0。生成器线程数固定为 1。

| 测量 | 运行 1 | 运行 2 | 预算 |
|---|---:|---:|---:|
| 紧凑构建 | 21.65 s | 20.26 s | 900 s |
| OpenFOAM 导出 | 152.61 s | 80.74 s | 900 s |
| 生成器完整进程 | 180.98 s | 102.10 s | 构建与导出分别满足预算 |
| 生成器峰值 RSS | 983,875,584 B | 1,107,542,016 B | 6 GiB |

全量独立读取和逐边检查：177.64 s，峰值 RSS 1,638,334,464 B，线程数 1，均低于计划的
15 分钟和 6 GiB。正式 `polyMesh` 为 956,918,167 字节，低于 4 GiB 产物预算。

运行 1/2 的墙钟差异主要来自文件缓存和系统负载；验收保留两次实测，不只挑较快结果。

## 9. 局部预览

命令 `cartmesh_stage6_cli --preview` 从阶段 6 紧凑状态写出 R24 背景检查 VTU，字段为：

- `compact_cell_state`；
- `fluid_volume_fraction`；
- `cut_cell`；
- `small_cut_cell`。

meshio 5.3.5 实际读取 15,625 个点、13,824 个背景六面体并选出 1,506 个 Cut-cell，状态
为 pass。`artifacts/stage6_preview.vtu` 和 `artifacts/stage6_preview.png` 只是可交互检查入口，
不是完整求解器网格。

旧阶段 4 完整多面体 CLI 曾被误用于生成预览，R24 出现 124 个闭合失败和 257 个共享面失败；
大体积中间文件已清理，最小失败记录保存在
`artifacts/stage6_preview_legacy_path_failure.json`，没有纳入验收证据。

## 10. 测试

回退清理后的 Release 源码重新构建成功，CTest 21/21 通过，最后一轮总墙钟 19.65 s。阶段 6 自身 4/4：

- 紧凑立方体与基线等价；
- 嵌套薄壳区域与确定性；
- 稳态存储跟随表面邻域而非全域显式几何；
- 流式 binary OpenFOAM 连通分量控制体导出。

三个 Python 工具 `py_compile` 通过。工程测试不替代上述独立 reader 和 OpenFOAM 证据。

## 11. 复现命令

```sh
cmake --build build/release --parallel 4
ctest --test-dir build/release --output-on-failure

build/release/cartmesh_stage6_cli \
  --stl benchmarks/public/stage4/stanford_bunny_libigl_binary.stl \
  --resolution 216 \
  --openfoam-case artifacts/stage6_10m_case \
  --report artifacts/stage6_10m_export.json

.venv/bin/python tools/stage6_binary_polymesh_verify.py \
  --case artifacts/stage6_10m_case \
  --output artifacts/stage6_10m_external_reader.json \
  --diagnose-cell-edges

build/release/cartmesh_stage6_cli \
  --stl benchmarks/public/stage4/stanford_bunny_libigl_binary.stl \
  --resolution 24 \
  --preview artifacts/stage6_preview.vtu \
  --report artifacts/stage6_preview_compact.json

.venv/bin/python tools/meshio_verify.py artifacts/stage6_preview.vtu \
  --field cut_cell --screenshot artifacts/stage6_preview.png \
  --title "Stage 6 local cut-cell preview (not the complete solver mesh)"

.venv/bin/python tools/stage6_openfoam_verify.py \
  --case artifacts/stage6_medium_bunny_r96_case \
  --output artifacts/stage6_medium_checkmesh.json \
  --log-output artifacts/stage6_medium_checkmesh.log
```

最后一条需要本机已有 OpenFOAM Docker 镜像；脚本固定 `--network none`。

## 12. 阶段关闭前的下一项工作

不能通过提高 `meshQualityDict` 阈值、翻转孤立边界面、删除小 Cut-cell 或只保留自建读取器
通过来关闭阶段 6。下一步需要设计可扩展、确定且体积守恒的质量策略，候选是：

1. 仅对非星形组件做局部共形分裂，避免逐凸片全量二次复杂度；
2. 对极小/高偏斜 Cut-cell 做质量约束的守恒聚合，候选合并必须同时通过正体积、星形性/面金字塔和偏斜上限，并保留稳定 ID 映射；
3. 先在 R24 达到 OpenFOAM `Mesh OK.`，在此之前不放大到 R48/R96；
4. 聚合后重新验证体积、边界、区域、owner/neighbour、全量边闭合、确定性和千万级预算；
5. OpenFOAM 2606 必须输出 `Mesh OK.` 后，才可将 `stage6Complete` 和
   `solverReadyCutCellMesh` 设为 true。

阶段 6 关闭前不定义或启动阶段 7。
