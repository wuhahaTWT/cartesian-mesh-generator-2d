# cartmesh2d 2D-V visualization

这是一个**薄层后处理工具**。它只读取 Stage 2D-6 已导出的 `CM2D v1` / `quality.json`，不重新执行分类、Cut-cell、Quadtree 或聚合算法。

## 使用

```bash
python3 cartmesh2d/tools/visualization/render_cm2d.py \
  mesh.cm2d mesh.svg \
  --quality mesh.quality.json \
  --color-by level
```

可选：

```text
--color-by level|area|none
--labels
--width 1200
--height 900
```

## 当前显示内容

- stabilized solver cell polygons
- internal edges
- embedded physical boundary
- outer/domain boundary（如果存在）
- unclassified boundary（会以警示样式出现）
- source Quadtree level（由 `sourceKey` 的低 6 bit 读取）
- cell area 着色模式
- cell ID 标签（可选）
- `quality.json` 的 counts / min area / min edge / aspect / skewness / topology audit 摘要
- source Cut-cell 数量与 source small-cell 数量

## 重要语义

`CM2D v1` 保存的是**small-cell stabilization 之后的最终 solver topology**。因此：

- embedded-boundary edge 可以准确显示最终物理边界；
- source level 可以准确显示当前 topology cell 的 `sourceKey`；
- 但一个已经被 agglomerate 的 tiny source cell 不再是独立 final cell。

因此第一版**不会伪造“原始 tiny cell 的最终位置”**。当前 SVG 只在 quality panel 中显示 source small-cell 总数。若后续需要逐个显示“合并前 small-cell → 合并后目标 cell”的映射，应新增显式 visualization sidecar，而不是从最终 CM2D 猜测。

## 验收

`cartmesh2d/tests/visualization_test.py` 使用独立 CM2D fixture 验证：

- SVG 可以生成；
- cell / edge layer 存在；
- internal / embedded boundary 样式存在；
- quality panel 存在；
- 非零 topology audit 的 CM2D 会被拒绝。

GitHub Actions 还会对 rectangle / circle / concave / airfoil-like 四个真实 Stage 2D-6 end-to-end 输出生成 SVG 并检查非空结果。
