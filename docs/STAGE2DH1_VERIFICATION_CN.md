# Stage 2D-H1：高密度尺寸场

日期：2026-08-25

## 1. 阶段边界

本阶段只扩展现有 Quadtree refinement policy，不修改 Cut-cell、稳定化、solver
topology、OpenFOAM writer 或桌面应用。本阶段也不实现贴体边界层。

已有 `minimumLevel` 继续控制全域基础密度，`boundaryLevel` 继续控制物面相交叶层级。
新增能力为：

- 可重复的物面绝对距离带 `--distance-band distance target-level`；
- 可重复的轴对齐局部区域 `--refine-box xmin ymin xmax ymax target-level`；
- 下游矩形区域作为第一版尾迹加密原语；
- 所有重叠尺寸场取最大目标层级，输入顺序不影响叶 ID 或输出；
- `.sizing.json` 保存完整输入策略、计算域、层级直方图、区域覆盖叶数和 2:1 平衡报告。

## 2. 正确性与失败规则

核心回归检查：

- 尾迹矩形内所有正面积相交叶达到目标层级；
- 重叠的更细矩形取得更高目标层级；
- 无关远场不被尾迹目标层级全域污染；
- 反转 distance band / box region 顺序后 Morton key 序列一致；
- refinement 后总面积、确定性遍历和 2:1 平衡保持原有门禁。

以下输入 fail-closed：

- 非有限或退化矩形；
- 目标层级为 0 或超过 `max-level`；
- 完全不与计算域重叠的矩形；
- 非有限/负距离带或距离带层级超过 `max-level`；
- 未知或缺参数的尺寸场 CLI 选项。

部分超出计算域的矩形允许使用，只对与计算域正面积相交的部分生效。

## 3. 开发期针对性验证

只构建并运行 Quadtree 测试：

```sh
cmake --build build/release-2d --target cartmesh2d_quadtree_tests cartmesh2d_cli -j4
ctest --test-dir build/release-2d \
  -R '^cartmesh2d_stage2_quadtree_tests$' --output-on-failure
```

结果：`1/1 PASS`。

真实低成本翼型产品使用 `minimum-level=3`、物面 `max-level=7`、`0.05` 距离带和
`[0.8,-0.1] -> [1.3,0.1]` 的 level-6 尾迹矩形：

```text
leaf_count=1585
source_cells=1160
stabilized_cells=1146
sizing_distance_bands=1
sizing_box_regions=1
2:1 violations_before/after=79/0
level histogram=3:27, 4:52, 5:124, 6:926, 7:456
```

CLI 生成 VTK、CM2D、quality/viz JSON 和新的 sizing JSON，CM2D 独立回读通过。
一个完全位于域外的 `[10,10] -> [11,11]` box 返回非零退出并报告：

```text
invalid sizing field: box refinement region does not overlap domain
```

## 4. 阶段完成门

阶段结束只运行一次完整二维 CTest。这里不运行 OpenFOAM `checkMesh`、复杂工程案例或
百万级规模门；这些属于真正发布验收，不用日常开发重复消耗时间。

实际命令：

```sh
cmake --build build/release-2d -j4
ctest --test-dir build/release-2d -R '^cartmesh2d_' --output-on-failure
```

结果：`53/53 PASS`，CTest 实际耗时 `9.45 s`。该结果只代表二维项目完整 CTest；没有
把构建过程中出现的三维 target 写成三维测试通过，也没有宣称完成发布级外部验证。
