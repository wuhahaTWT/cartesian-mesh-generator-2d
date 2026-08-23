# Stage 2D-P1 快速 Cartesian 核心验收

日期：2026-08-24

## 1. 阶段边界

P1 只优化单闭合边界下的 Cartesian/Quadtree 主流水线：边界查询、逐层细化、2:1 平衡、全局拓扑装配和小单元邻接分析。没有进入多环域 D1，也不把 P1 的拓扑通过描述成 solver-ready。

## 2. 实现

- 新增确定性的 `BoundarySegmentIndex2D` AABB BVH，替代网格分类和 Quadtree 细化中的逐线段全扫描。
- Quadtree 细化和平衡改为逐层/逐批处理，避免每细化一个叶子就重新排序和从头扫描。
- 2:1 平衡使用排序后的 sweep 邻接匹配。
- 全局拓扑的顶点查找使用 x/y 排序范围查询，Cartesian 边上的 hanging node 不再扫描全部顶点。
- 小单元候选从 `每个小单元 × 全部内部边` 改为一次构造确定性的局部邻接表。
- 顶点规范化公差按最小源背景单元尺度计算；修复百万级运行中域尺度公差吞掉合法短边的问题。
- 聚合面积比较使用面积量纲公差；没有降低现有最终质量门。
- CLI 在聚合或质量失败时打印具体 issue code、对象和数值诊断。

## 3. 正确性与确定性

- Release：`46/46` CTest PASS。
- 2D ASan + UBSan：`22/22` PASS。macOS 不支持 LeakSanitizer，因此明确使用 `ASAN_OPTIONS=detect_leaks=0:halt_on_error=1`；没有把首次因平台不支持而中止的运行写成通过。
- `meshio` 独立读取 rectangle、airfoil-like、NACA2412、superellipse 和百万级 circle VTK，点/单元计数与输出一致。
- 四个 R1 回归产品的 VTK、CM2D、quality JSON、viz JSON 与 R1 基线逐字节相同。
- 17 层百万级产品连续两次生成的四类输出逐字节相同。

百万级输出 SHA-256：

```text
VTK          8a2d23b40ba6a9c5619289f2e63f9d68c7b88505b34a8e4bb6e9e22e8b87ff60
CM2D         bc988faab41afe675840e158278299e13d7c3200737bdc9fdb60266cae032fa7
quality JSON 96d690de21ea6532a7f59a9a3ecf614c06b4f795ada3f8839b973612963092af
viz JSON     4db75da54f33fca43da595261c638117aa3a5be9795722b4f277d2e4981d7611
```

## 4. 真实墙钟与内存

硬件上的 Release 端到端测量包含输入、细化、Cut-cell、全局拓扑、小单元聚合、质量报告、VTK/CM2D/JSON 写出以及 CM2D 独立回读。输入均为 `examples/acceptance/circle.xy`，外部流体。

| max level | leaf | source cells | stabilized cells | small alpha | wall | maximum RSS |
|---:|---:|---:|---:|---:|---:|---:|
| 9 | 5,872 | 3,712 | 3,448 | 0.10 | 0.14 s | 15,269,888 B |
| 11 | 23,896 | 14,784 | 13,728 | 0.10 | 0.61 s | 52,838,400 B |
| 14 | 192,796 | 118,460 | 110,004 | 0.10 | 5.27 s | 391,004,160 B |
| 17 | 1,542,700 | 946,304 | 834,412 | 0.25 | 46.01 s | 1,258,352 KiB sampled peak |

优化前同一 level 9 端到端基线为 6.61 s；P1 为 0.14 s，墙钟加速约 47 倍。这里的“百万级”指 1,542,700 个实际 Quadtree leaf；最终 solver source cell 为 946,304，未冒充超过一百万。

## 5. 百万级产品结果

通过命令：

```bash
cartmesh2d_cli circle.xy circle_l17_alpha025 17 0.25 0.25 exterior
```

关键结果：

```text
source_fluid_area=5.87855
expected_fluid_area=5.87855
vertices=1183774
edges=2018186
min_cell_area=1.30968e-10
min_edge_length=4.92136e-10
max_edge_aspect_ratio=46507.8
max_centroid_skewness=0.333354
```

拓扑、正面积、面积守恒和输出回读均通过，但最大边长比仍为 46,507.8，因此这不是 solver-ready 宣告。

## 6. 保留的失败证据

1. `max-level=17, padding=0.25, small-alpha=0.10`：拓扑与聚合通过，但最终质量门拒绝 28,768 个面积不大于约 `1.01e-10` 的单元。P1 没有降低质量门；百万级通过产品使用更严格的 `small-alpha=0.25`。
2. `max-level=17, padding=0.20, small-alpha=0.25`：源拓扑在 cell `388458` 拒绝一个 `degenerate polygon edge`。该输入/命令保留为后续几何构造最小化回归线索，不写成成功。
3. 百万级成功产品仍有很高的 edge aspect ratio；该问题属于后续 S1 solver-ready 质量阶段，不能用 P1 性能结果掩盖。

## 7. P1 结论

P1 达到单环域快速核心目标：真实端到端流水线已从数千单元的二次瓶颈扩展到百万级 leaf，确定性、拓扑门、面积守恒、外部读取和 sanitizer 均有证据。下一阶段只能进入 D1 多环/多连通域，不能提前宣称 solver-ready。
