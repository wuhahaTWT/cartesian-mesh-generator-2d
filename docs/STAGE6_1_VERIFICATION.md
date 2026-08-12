# Stage 6.1 验证：统一自适应 Cut-cell 完整拓扑

日期：2026-08-12（Asia/Shanghai）

状态：**Stage 6.1 PASS；Stage 6 总体仍未完成；Stage 6.2 尚未开始。**

## 1. 关闭结论

现有 reference `OpenFoamWriter` 已从仅接受 `UniformCartesianGrid` 泛化为同时接受
`LinearOctree`。固定 adaptive cube 与非凸 L-prism 均从 STL 重新生成完整 ASCII
OpenFOAM `polyMesh`，随后通过：

1. 项目内部几何/拓扑不变量；
2. 不链接 cartmesh 的独立 ASCII `polyMesh` 全量读取器；
3. OpenFOAM 2606 `checkMesh -constant -allTopology`，明确输出 `Mesh OK.`。

因此 Stage 6.0 记录的“adaptive CLI 明确拒绝 OpenFOAM”缺口已经关闭。这个结论只覆盖
Stage 6.1 固定案例和 reference writer，不表示复杂几何质量、坏形状稳定化、紧凑 adaptive
数据结构或千万级 adaptive 已完成。

## 2. 实现边界

### 复用并泛化

- 继续使用一个 `OpenFoamWriter.cpp`，没有新增平行 `AdaptiveOpenFoamWriter`；
- 均匀与八叉树输入先形成内部背景网格视图，共用后续 component cell、公共多边形细分、
  arrangement 面、wall/farfield patch、点焊接、排序和 OpenFOAM 序列化；
- 八叉树按稳定 Morton 叶顺序生成 solver cells；
- 同层接口只处理一次；2:1 跨层接口由一个 coarse face 对整组 fine faces 统一处理；
- 每个 coarse/fine 子接口与整个 coarse face 都分别检查开口面积守恒；
- owner/neighbour 按 solver cell ID 规范为 `owner < neighbour`，必要时同步反转面方向；
- 远场由叶的实际 `face_neighbors()` 空集判定，不依赖 `i/j/k`；
- 面覆盖容差由参与叶的局部 AABB 尺度导出，不再依赖全局 uniform spacing；
- `cartmeshCellMapping.json` 保存无损字符串形式的 64 位 Morton node code、临时 leaf index、
  component/local piece 与 solver cell ID 的确定性映射。

### 删除或收敛

- 删除 CLI 中 adaptive `--openfoam-case` 的主动拒绝分支；
- 移除 writer 主流程对 `cell_key/nx/ny/nz/spacing` 的直接依赖；
- 原有均匀邻接枚举被收进均匀背景适配层，不再与完整 writer 流程交织。

### 明确未修改

- 没有修改 `ScalableOpenFoamWriter` 或 Stage 6 千万级 binary 路径；
- 没有加入聚合、共形质量分裂、kernel tetra、删小单元、翻面或放宽 checker；
- 没有实现 Stage 6.2 原生质量评估器；
- 没有实现任何 Stage 7 boundary/prism/layer/projection 功能。

## 3. 最小回归

`tests/stage3_test.cpp` 新增/扩展两类永久回归：

1. **coarse-fine 无切割**：基础 level 1 八叉树只细分一个叶，保留真实 2:1 接口；完整
   OpenFOAM 输出覆盖全部 15 个 fluid leaves，owner/neighbour 合法且 farfield 存在。
2. **coarse-fine + cut**：嵌入立方体穿过多层八叉树；实际包含跨层连接，完整输出两次
   逐字节一致，并检查 Morton → solver mapping。

真实 CLI CTest 还固定运行 adaptive cube 和 adaptive L-prism 的 `--openfoam-case` 路径。

## 4. Release 构建和全量测试

```sh
cmake --preset release
cmake --build --preset release --parallel 4
ctest --preset release --output-on-failure
```

结果：**22/22 PASS**，最终测试墙钟 **14.48 s**。新增门包括 adaptive cube 完整 case 和
adaptive L-prism 完整 case；Stage 0–6 既有测试全部保持通过。

环境：MacBookAir10,1，Apple M1，8 核，8 GiB；macOS 26.5.2 arm64；GCC 15.2.0；
Release；正式生成器线程数 1。

## 5. 真实 adaptive CLI 与内部不变量

| 案例 | background leaves | fluid background cells | Cut-cells | result hash | wall | peak RSS |
|---|---:|---:|---:|---|---:|---:|
| cube | 2,584 | 2,368 | 1,016 | `c45a6843e6a72f49` | 2.115313 s | 19,709,952 B |
| L-prism | 2,234 | 1,714 | 684 | `a47fd5b59db2b115` | 1.280080 s | 15,433,728 B |

两例的 `nonclosedCellCount`、`negativeVolumeCellCount`、
`sharedFaceMismatchCount` 和 classification conflict 均为 0。

## 6. 独立 ASCII polyMesh 全量读取

新增 `tools/stage61_ascii_polymesh_verify.py`。它不链接 cartmesh，直接解析 ASCII
`points/faces/owner/neighbour/boundary`，并检查：

- 列表长度和点/单元索引范围；
- boundary patch 连续覆盖；
- `owner < neighbour` 和内部 owner 排序；
- 每个 solver cell 的每条边恰好被两张面引用；
- 面积向量闭合和基于有向面的正体积；
- 五个 OpenFOAM 核心文件 SHA-256。

| 案例 | points | faces / internal | solver cells | min volume | nonclosed / nonpositive / bad cell-edges |
|---|---:|---:|---:|---:|---:|
| cube | 4,254 | 9,654 / 6,936 | 2,704 | 1.5625e-05 | 0 / 0 / 0 |
| L-prism | 3,251 | 7,145 / 5,356 | 1,954 | 6.25e-05 | 0 / 0 / 0 |

机器记录：`artifacts/stage61/adaptive_cube_reader.json` 和
`artifacts/stage61/adaptive_l_reader.json`。

## 7. OpenFOAM 2606 外部验收

检查在本机已有 `opencfd/openfoam-run:2606` 镜像内运行，Docker 网络关闭，源 case
复制到临时可写目录，完整 stdout 保存为日志。

| 案例 | cells | max non-orthogonality | max skewness | return | 结果 |
|---|---:|---:|---:|---:|---|
| cube | 2,704 | 35.2644 | 0.5 | 0 | `Mesh OK.` |
| L-prism | 1,954 | 54.7356 | 0.914585 | 0 | `Mesh OK.` |

证据：

- `artifacts/stage61/adaptive_cube_checkmesh.json/.log`；
- `artifacts/stage61/adaptive_l_checkmesh.json/.log`。

## 8. 确定性与均匀路径不回归

cube 与 L-prism 分别完整生成两次：项目 result hash、五个核心 `polyMesh` SHA-256 和
`cartmeshCellMapping.json` 均逐字节一致。

- cube mapping SHA-256：`b9d893bb504fdf7d45ead14c6070ee07b83788927d42b5d57e800f104b3539bf`；
- L-prism mapping SHA-256：`b70ef8cd86cb469db70e06f99354306ee6f83ecca5d8be4187cf9e034be5e01b`。

均匀 R8 cube 使用与 Stage 6.0 相同配置重跑，result hash 仍为
`e6fd5337f55140ea`，五文件 SHA-256 逐项保持：

| file | SHA-256 |
|---|---|
| points | `7d0ddf0e2fe57bf4fb7f589b3b9e2f4e6dfe4955b659f85e94a0eaaf0cab6559` |
| faces | `c333d66020fe7fe4b2c5de2a19ee4918c380da3884740efdf1dcbd7b6cb17791` |
| owner | `6ef35a395352dd62a59eff1fd2d42fc2b90f17fe8337cfb068ea448290fd996d` |
| neighbour | `00d889d1964fc3cc7f15645f4012e4bfd20fbb1c401ed1572ca8a173d17147d9` |
| boundary | `f5817a95bdb15dde0b4d543d4c6119dda55de3a8fc3d1b43b8ee7de214fd0970` |

这证明公共 writer 泛化没有改变既有均匀核心文件语义或字节结果。

## 9. 复现命令

```sh
build/release/cartmesh_cutcell_cli \
  --stl tests/data/closed_unit_cube_ascii.stl \
  --adaptive --base-level 2 --max-level 4 --surface-level 4 \
  --distance 0.1:3 --padding-fraction 0.1 --no-vtk \
  --openfoam-case artifacts/stage61/adaptive_cube_case \
  --geometry-output artifacts/stage61/adaptive_cube_geometry.json \
  --report artifacts/stage61/adaptive_cube_report.json

.venv/bin/python tools/stage61_ascii_polymesh_verify.py \
  --case artifacts/stage61/adaptive_cube_case \
  --output artifacts/stage61/adaptive_cube_reader.json

.venv/bin/python tools/openfoam_stage3_verify.py --milestone stage4 \
  --case artifacts/stage61/adaptive_cube_case \
  --project-report artifacts/stage61/adaptive_cube_report.json \
  --output artifacts/stage61/adaptive_cube_checkmesh.json \
  --log-output artifacts/stage61/adaptive_cube_checkmesh.log
```

L-prism 使用相同命令结构和 `tests/data/nonconvex_l_prism_ascii.stl`。

## 10. Stage 6.1 关闭边界

Stage 6.1 的通过含义是：固定解析 cube/L-prism 的 adaptive Cut-cell 已能形成确定、完整、
可独立读取并由 OpenFOAM `Mesh OK.` 接受的控制体拓扑。

它不意味着：

- 旧千万级 uniform Bunny 的 face-pyramid/skewness 阻断已经修复；
- Bunny 或任意复杂 STL 的 adaptive OpenFOAM 已通过；
- 小 Cut-cell 已稳定化；
- reference writer 已扩展到百万/千万级；
- Stage 6 或 Stage 7 可以关闭。

下一步只能是 Stage 6.2 原生质量评估器，并且必须等待用户明确确认。
