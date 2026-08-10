# 阶段 5 用户亲自检查点

日期：2026-08-09（Asia/Shanghai）

结论：请现在检查，不要等阶段 6、7 全做完再看。阶段 5 改动了稳定 ID、局部重构范围、
Cut-cell 几何复用和旧到新映射契约；现在确认方向的成本最低。

## 先看这张三联图

![阶段 5 局部轮廓增量重构对比](../artifacts/stage5_local_contour_comparison.png)

- 左：修改前的流体体积分数切片；
- 中：修改后的流体体积分数切片；
- 右：本次实际重构区，红色 `1` 为重构，蓝色 `0` 为复用。

这张图只是帮你理解行为，不是验收证据的替代品。真正的通过依据是增量/全量 hash 一致、
零拓扑失败、精确重叠映射、独立 meshio 读取和 OpenFOAM 回归。

## 你应该亲自确认的四件事

1. 只有几何变化附近被重构，不是全域都染成红色。
2. 修改后边界是你期望的局部轮廓，其它部分没有被意外改动。
3. 你接受阶段 5 只解决增量局部重构，不包含 CGNS、GUI、CFD 求解器或阶段 6 的千万级完整网格。
4. 在启动阶段 6 前，先打开 `docs/STAGE5_VERIFICATION.md` 看一遍通过项、失败案例和当前边界。

## 可直接打开的文件

- 终态验收：`artifacts/stage5_acceptance.json`
- 局部轮廓生成器报告：`artifacts/stage5_local_contour.json`
- 独立读取报告：`artifacts/stage5_local_contour_meshio.json`
- 旧到新精确映射：`artifacts/stage5_local_contour_mapping.json`
- 完整验证记录：`docs/STAGE5_VERIFICATION.md`
