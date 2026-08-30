# R1-C checkpoint：typed recovery request 与 local refinement closure

## 1. 状态边界

这是 R1-C 的可提交 checkpoint，不是 R1-C 全部完成。已完成：

- 三类 unsafe intersection conflict 不再只抛普通字符串错误，而是生成
  `ConstructionRecoveryRequest2D`；
- shared Cut-cell 将 typed signal 转成 leaf-scoped unsupported result，保留 request；
- Hybrid production path 收到 request 后重建干净 construction registry，只请求细化
  受影响 leaf，并以 face-neighbour 驱动的 2:1 closure 收敛；
- 每次 split 记录 parent key、四个 child keys、parent level、requested/closure 来源；
- refined Quadtree lineage 与 solver source ID 分离存储，避免把 Morton key 当 source index；
- bounded grid-phase 与 source-parameter resample candidate model 已建立，候选只有同时通过
  面积守恒、feature compatibility 和原 hard-quality gate 才可选择。

尚未完成：production rephase/resample 的局部 polygon rebuild、三类候选的真实质量比较，
以及在真实 narrow-gap/sharp-tail 上触发 recovery 后取得更好 Q1 结果。因此不得登记 R1-C
完成，更不得宣称 Q1 已修复。

## 2. typed request

request 保留：stable conflict key、conflict kind、support/grid-line identity、source segment、
source parameter、原位置、冲突位置、`localH`、推荐 decision 与固定 fallback order：

```text
refine -> rephase -> resample
```

当前三类 conflict 是：

- `late_nonincident_feature`；
- `nonincident_feature_snap`；
- `nonincident_grid_corner`。

最小回归保留了原 Q2-A sub-roundoff 坐标：两个非 incident supports 位于
`y=0.5` 与 `y=0.5+2e-16`，共同查询 `x=0.5`。现在结果是 typed
`nonincident_grid_corner` request，`source_parameter≈0.5`，而不是静默 welding 或由上层
解析异常字符串。

## 3. recovery candidates 与提交门

`planConstructionRecoveryCandidates()` 产生：

1. affected leaf refinement；
2. 最多四个 deterministic phase offsets：`±0.125 localH`、`±0.25 localH`；
3. conflict parameter 两侧、严格位于 `(0,1)` 的 source resample points。

resample point 只能表示成 `segment.a + t*(segment.b-segment.a)`，所以 API 本身不能产生
法向漂移或跨 feature 搜索。`selectConstructionRecoveryCandidate()` 只接受同时满足：

- `areaConserved=true`；
- `featureCompatible=true`；
- `hardQualityPass=true`。

选择按 quality rank、decision priority、stable ordinal 排序。未评估或 hard-quality
失败的候选不会进入 committed geometry。

## 4. neighbour-driven local refinement

`Quadtree2D::refineLeavesWithClosure()` 在同一 transaction 中：

- 对 requested keys 排序去重并预检；任一 missing/max-level key 会令整个 batch 原子拒绝；
- 批量 split requested leaves；
- 从 `faceNeighbors()` 只收集 level difference 大于 1 的 coarse neighbours；
- 分轮传播 closure 至无 violation；
- 输出 requested/closure split 数、迭代数与逐 parent-child lineage；
- 每轮重新使用确定性 Morton ordering。

深链回归证明 closure 确实传播、最终 violation 为 0、总域面积保持 16；重复 transaction
的 key、lineage 和计数相同。混合“可细化 leaf + max-level leaf”请求原子失败，不留下
半个 transaction。

## 5. 真实网格证据

### narrow gap

- construction / solver cells：`3244 / 3189`；
- expected / actual fluid area：`10.24 / 10.239999999999966`；
- area error：`-3.3750779948604759e-14`；
- independent OpenFOAM reader：PASS，3189 cells、6841 internal faces、两个 wall patch
  与四个 domain patches；
- current/repeat construction 与 solver CM2D 逐字节一致；
- solver SHA-256：
  `00d5f20abdaaf20df1aabc4599c089c3c89b28184e713e598314d487beae0505`，与 Q2-A
  baseline 相同；
- recovery requests / refinement passes：`0 / 0`。

### sharp trailing edge

- construction / solver cells：`3412 / 3391`；
- expected / actual fluid area：`15.855499999999999 / 15.855500000000021`；
- area error：`2.1316282072803006e-14`；
- independent OpenFOAM reader：PASS，3391 cells、7237 internal faces、wall 与四个
  domain patches；
- current/repeat solver CM2D 逐字节一致；
- solver SHA-256：
  `1faa8b72e7fd9e9a66b35e8ea55deeeec91aad4a3ce1d45de11d9599180ae819`，与 Q2-A
  baseline 相同；
- recovery requests / refinement passes：`0 / 0`。

两例 recovery 为零是重要的负证据：现有 Q1 hard failures 不是本 checkpoint 覆盖的
unsafe snap conflict，不能靠多 refine 一次或改报告字段解决。

## 6. 质量事实

- narrow-gap quality contract：FAIL；`min(face/localH)=0.0078431372549030553`，另有
  face weight、volume ratio、angle、non-orthogonality、skewness 等 hard issues；
- sharp-tail quality contract：FAIL；`min(face/localH)=0.00042023831386472921`，同样存在
  多类 hard issues；
- 两例 legacy solver-quality 均 PASS，但不得冒充 Q1 typed contract PASS；
- 未降低阈值、未删除 cell、未移动原 wall feature、未隐藏告警；
- Docker daemon 本轮仍未运行，真实 OpenFOAM v2606 `checkMesh` 为 NOT_RUN/pending；
  independent OpenFOAM reader 不能替代它。

## 7. 下一步

R1-C 后续必须为真实 bad patch 实现 transaction-scoped rephase/resample polygon rebuild，
逐候选执行面积、feature 和 typed hard-quality 比较；refinement 无改善或到达 max level 时
不能直接宣告成功。只有 narrow-gap、sharp-tail 与最小 conflict fixture 都有真实 recovery
证据后，才可关闭 R1-C 并进入 R1-D edge-incidence / patch-local topology transaction。

证据文件：

- `build-q2a/evidence/r1c-narrow_gap.hybrid.*`；
- `build-q2a/evidence/r1c-narrow_gap.independent-foam.json`；
- `build-q2a/evidence/r1c-narrow_gap.hybrid.solver.svg`；
- `build-q2a/evidence/r1c-sharp_trailing_edge.hybrid.*`；
- `build-q2a/evidence/r1c-sharp_trailing_edge.independent-foam.json`；
- `build-q2a/evidence/r1c-sharp_trailing_edge.hybrid.solver.svg`。
