# R1-B：typed event key 与 feature-compatible construction decision

## 1. 阶段结论

R1-B 已把 Q2 的 proximity canonicalization 从“registry 内部的一段候选判断”收进
统一的 `ConstructionVertexStore2D::decideProximity()`，并让 source endpoint、
wall/transition-grid event 与 Cartesian grid vertex 具有可查询、可冲突检测的 typed
exact-key alias。旧的全量 vertex scan 本阶段仍作为 fail-closed oracle；它不再决定
生产结果，但在每次 proximity decision 后逐项核对 ID 与位移。全局扫描退出热路径是
后续 R1-G 的验收项，本阶段不提前宣称完成。

本阶段没有改变 Cut-cell、transition 或 solver polygon 的坐标与拓扑。superellipse 与
sharp trailing edge 的 construction/solver CM2D 均与 Q2-A 基线逐字节一致。

## 2. 决策模型

`ConstructionVertexProposal2D` 明确携带：

- 原始位置、`localH`、无量纲 snap fraction；
- feature class、feature owner、support ID；
- 输入 feature 是否不可移动。

`ConstructionVertexDecisionResult2D` 返回 canonical stable ID/point、位移、被检查候选数、
decision 与 reason。候选只来自 support-scoped multilevel feature index；适用半径仍为
`snap_fraction * min(candidate_local_h, anchor_local_h)`。最终选择顺序固定为距离、feature
优先级、坐标、stable ID，避免容器遍历顺序进入结果。

feature compatibility 的当前矩阵为：

| proposal | anchor | 结果 |
|---|---|---|
| convex/concave/domain/fixed feature | 同 class、同 owner | 允许 exact reuse；非零位移仍拒绝 |
| convex/concave/domain/fixed feature | 其他 | 拒绝，理由 `immutable_feature_conflict` |
| GapSideA/B | 同 side class、同 owner | 允许 |
| GapSideA/B | 异侧或异 owner | 双向拒绝，理由 `gap_side_conflict` |
| smooth/grid/mutable/unclassified | 普通或受保护 anchor | 允许 anchor 吸收 movable proposal |

最后一条不会移动 wall sharp feature：变化的是 movable sample，既有 sharp anchor 的坐标
始终不变。

## 3. typed exact identity

`ConstructionVertexStore2D` 维护 `StableVertexKey2D -> StableVertexId2D` 唯一映射，同时在
vertex record 中保存排序去重的 alias。重复绑定同一个 key/ID 是幂等操作；同一个 key
绑定不同 ID 会显式抛错。

当前 production construction 已绑定并使用：

- `SourceVertex(support, endpoint_ordinal)`；
- `WallGridIntersection(support, axis, logical_coordinate)`；
- `TransitionVertex(support, axis, logical_coordinate)`；
- `GridVertex(logical_x, logical_y)`。

event cache 命中时必须与 typed event key 解析到同一 stable ID；Cartesian corner 在进入
坐标兼容表前先按 typed grid key 解析，并校验 committed geometry 完全相同。机器报告新增
`r1b_exact_key_count` 与 `r1b_proximity_decision_api`。

## 4. 回归与 fail-closed 覆盖

新增测试覆盖：

- sharp/concave class 与 owner 冲突；
- GapSideA/GapSideB 及同侧不同 owner 不得合并；
- typed key 幂等绑定与 conflicting ID 显式失败；
- 等距候选以几何排序确定性选择；
- gap 两侧即使落在 snap radius 内仍保持不同；
- immutable sharp proposal 不产生非零位移；
- reversed support、coarse/fine event、双轴 grid-corner event 共享 stable ID；
- source endpoint、wall-grid event 与 grid vertex typed alias 可直接解析。

验证命令：

```text
cmake --build build-q2a -j 6
ctest --test-dir build-q2a --output-on-failure -j 6
```

结果：`75/75` PASS，55.65 s。仓库 `AGENTS.md` 中的“73 项”已落后于当前 CMake
实际枚举；本阶段没有删减测试。

## 5. 真实网格证据

### 5.1 superellipse

```text
build-q2a/cartmesh2d_hybrid_cli examples/complex/superellipse_24.xy \
  build-q2a/evidence/r1b-superellipse 6 3 6 3 0.015 1.15 1.0 \
  build-q2a/evidence/r1b-superellipse-case .01
```

- construction cells / solver cells：`772 / 795`；
- construction vertices / edges：`1168 / 1940`；
- expected / actual fluid area：`14.871177621957496 / 14.871177621957527`；
- area error：`3.0198066269804258e-14`；
- independent hybrid reader：PASS，零 overlap、零 non-manifold，最小正面积
  `0.0009093791952322627`；
- independent OpenFOAM reader：PASS，795 cells，1871 internal faces，wall 与四个
  domain patches 同时存在；
- typed exact key count：`1796`；
- 11 次 source movement 均记录为 `accepted / compatible_proximity_anchor`，最大
  `displacement/local_h = 0.0071081732327125339`；
- solver CM2D SHA-256：
  `b9fe6b0670b89a326e7ae31fc69aac052ee3347a4ffdc9d103e351435c82529d`，与 Q2-A
  `shared-superellipse` 相同。

### 5.2 sharp trailing edge

```text
build-q2a/cartmesh2d_hybrid_cli examples/h4_3/sharp_trailing_edge.xy \
  build-q2a/evidence/r1b-sharp_trailing_edge 8 3 8 4 0.012 1.15 1.0 \
  build-q2a/evidence/r1b-sharp_trailing_edge-case .01
```

- construction cells / solver cells：`3412 / 3391`；
- construction vertices / edges：`4160 / 7572`；
- expected / actual fluid area：`15.855499999999999 / 15.855500000000021`；
- area error：`2.1316282072803006e-14`；
- independent OpenFOAM reader：PASS，3391 cells，7237 internal faces，wall 与四个
  domain patches 同时存在；
- typed exact key count：`5542`；
- local-termination 路径未执行 transition resampling，intersection movement record
  数为 `0`；尖尾原始 wall feature 未移动；
- solver CM2D SHA-256：
  `1faa8b72e7fd9e9a66b35e8ea55deeeec91aad4a3ce1d45de11d9599180ae819`，与 Q2-A
  `shared-sharp_trailing_edge` 相同。

两例 current/repeat 的 construction CM2D、solver CM2D、construction JSON 与
intersection JSON 均逐字节一致。

## 6. 质量事实与未完成边界

- 两例 `solver_quality_valid=true`，但 Q1 typed quality contract 均仍如实为 `FAIL`；
- superellipse 当前 hard issues 是 face weight 16、volume ratio 10；其
  `min(face/local_h)=0.016222929810346797`，原 Q2 微短边没有回归；
- sharp-tail 仍有既知 short-face 等问题，
  `min(face/local_h)=0.00042023831386472921`；本阶段没有降低阈值、删除 cell 或隐藏告警；
- Docker daemon 本轮不可用，故真实 OpenFOAM v2606 `checkMesh` 状态为
  **NOT_RUN / pending**。独立 OpenFOAM reader 已通过，但不得冒充 `checkMesh`；
- legacy full scan 仍存在，仅作 fail-closed oracle。R1-G 前不得宣称全局候选搜索已退出；
- half-edge-lite/edge-incidence 与 patch-local transaction 尚属 R1-D，未在本阶段
  冒充完成。

## 7. 证据文件

- `build-q2a/evidence/r1b-superellipse.hybrid.{cm2d,solver.cm2d,vtk,solver.vtk,json}`；
- `build-q2a/evidence/r1b-superellipse.hybrid.{construction,intersections,quality-contract}.json`；
- `build-q2a/evidence/r1b-superellipse.independent-{hybrid,foam}.json`；
- `build-q2a/evidence/r1b-superellipse.hybrid.solver.svg`；
- `build-q2a/evidence/r1b-sharp_trailing_edge.hybrid.{cm2d,solver.cm2d,vtk,solver.vtk,json}`；
- `build-q2a/evidence/r1b-sharp_trailing_edge.hybrid.{construction,intersections,quality-contract}.json`；
- `build-q2a/evidence/r1b-sharp_trailing_edge.independent-foam.json`；
- `build-q2a/evidence/r1b-sharp_trailing_edge.hybrid.solver.svg`。

下一阶段进入 R1-C：把 unsafe conflict 变成 typed refine/resample/rephase request，
用 neighbour-driven 2:1 closure、source-parameter resampling 与 bounded grid-phase
candidates 处理；half-edge-lite / edge-incidence builder 与 patch-local transaction 在
R1-D 实施。
