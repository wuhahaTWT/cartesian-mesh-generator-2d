# Q2-B：Feature 约束下的局部短面修复

## 范围

Q2-B 只处理 Q2-A 共享交点构造之后仍存在的两个 local-termination 硬短面：
`narrow_gap` 与 `sharp_trailing_edge`。本阶段不改变 H4 wall/layer 构造，不引入
DCEL，不降低 Q1 或旧 solver-quality 阈值，也不删除失败 cell 后继续输出。

## 根因

- `narrow_gap`：固定 layer termination 角点与末级 Cartesian 分割线相距
  `0.000125`。公共分区把这个记账点插入 immutable layer 的直边，产生
  `face/local_h=0.007843137254903055` 的四个对称原子面。
- `sharp_trailing_edge`：Cartesian corner 与 mutable termination front 的交点近重合。
  两-cell 并集仍受第三个 Cartesian 邻居约束，需要闭合三-cell patch 才能去掉凹折；
  共 14 个 `face/local_h<0.01` 原子面，最差为
  `0.0004202383138647292`。

## 修复事务

`repairSolverShortFaces2D` 每轮只选择 `(face/local_h, owner, neighbour)` 字典序最小
的最严重短面，并构造有限的局部候选：

1. 若短面接触 immutable cell，只在其 mutable 邻域选择 pair；
2. 若两侧均 mutable，尝试直接 pair，并尝试两侧加一个 mutable 邻居的三-cell patch；
3. 对 patch 做精确面积并集，移除仅用于 common partition 的数值共线点，再进行最少凸分区；
4. immutable cell 只恢复原直线支撑，不移动任何非共线 feature vertex；
5. 候选必须同时满足：全局拓扑有效、旧 solver-quality 全通过、所有旧质量最坏值
   不回归、dimensionless 短面字典分数严格改善；否则回滚；
6. 每次提交后重新计算 source/local_h metadata 与全局质量，直到没有 Q1 硬短面。

这是面积守恒的局部 agglomeration/repartition，不是把坏 cell 从流体域删除。

## 定向结果

| case | face/local_h before | after | solver cells | accepted transactions |
|---|---:|---:|---:|---:|
| narrow_gap | 0.007843137254903055 | 0.014619607843137943 | 3189 -> 3185 | 4 |
| sharp_trailing_edge | 0.0004202383138647292 | 0.012005285498669291 | 3391 -> 3379 | 14 |

两个案例均满足：dimensionless face 硬违规为 0、旧 solver-quality 继续有效、面积误差
不变、boundary patch 计数不变、重复输出字节一致，且没有任何最坏质量指标回归。

定向验证：

```bash
./build-q2a/cartmesh2d_solver_export_tests
python3 tools/verification/verify_q2b_local_repair.py
```

完整五案例、独立读取器与 OpenFOAM `checkMesh` 属于代码冻结后的最终 Q2 验收，
不在本构建阶段重复执行。
