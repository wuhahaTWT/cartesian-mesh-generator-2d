# Stage 2D-0 验证记录

## 状态

**2D-0 几何内核已实现并通过独立本地 C++20 / CMake / CTest 验证。**

当前仍严格停止在二维几何层：未实现 Cartesian cell、Quadtree、Cut-cell、拓扑网格、求解器导出或可视化。

## 已实现对象

- `TolerancePolicy`
- `Point2D` / `Vector2D`
- `Segment2D`
- `AABB2D`
- `Polygon2D`
- `BoundaryLoop`

## 已实现算法

- dot / cross / orientation
- signed area / absolute area
- polygon centroid
- polygon AABB
- point-on-segment
- segment intersection：`None / Point / Overlap`
- winding-number point-in-polygon：`Inside / Outside / Boundary`
- boundary diagnostics：
  - fewer than three unique vertices
  - duplicate consecutive vertex
  - zero-length edge
  - non-adjacent edge self-intersection
  - zero / degenerate area
  - CW / CCW orientation
- 对合法 CW loop 安全归一化为项目标准 CCW；非法 loop 拒绝归一化

## 数值测试矩阵

已覆盖并通过：

1. rectangle：面积、带符号面积、质心
2. triangle：面积
3. concave L polygon：面积以及 inside / outside / boundary
4. 64 边离散圆：面积近似与中心点 inside
5. crossing segments：单点相交及交点坐标
6. endpoint touching segments：端点接触
7. parallel non-intersecting segments：无交
8. collinear overlapping segments：重叠线段及重叠区间
9. clockwise rectangle：方向检测与 CCW normalization
10. bow-tie loop：必须诊断 self-intersection，且拒绝 normalization
11. duplicate consecutive point：必须诊断 duplicate + zero-length edge

## 实际验证命令

```sh
cmake -S cartmesh2d -B cartmesh2d/build -DCMAKE_BUILD_TYPE=Release
cmake --build cartmesh2d/build
ctest --test-dir cartmesh2d/build --output-on-failure
```

独立验证环境：GNU C++ 14.2.0，C++20，`-Wall -Wextra -Wpedantic -Wconversion -Wshadow`。

结果：

```text
[100%] Built target cartmesh2d_geometry_tests
1/1 Test #1: cartmesh2d_stage0_geometry_tests ... Passed
100% tests passed, 0 tests failed out of 1
```

## 2D-0 验收对照

- [x] 原生二维数据结构，不依赖三维 `Point3D` / `AABB3D`
- [x] tolerance 集中管理
- [x] polygon signed area / area / centroid
- [x] robust segment intersection 三态结果
- [x] point-in-polygon 三态结果
- [x] boundary loop 输入诊断
- [x] CW / CCW 检测和合法 CCW normalization
- [x] bow-tie 最小失败案例
- [x] duplicate consecutive point 最小失败案例
- [x] C++20 target 可独立编译
- [x] CTest 数值测试通过
- [x] 未实现或冒充 Quadtree / Cut-cell / GUI

## 尚未开始

以下全部属于后续阶段，本阶段明确不实现：

- 2D-1 Cartesian background grid
- 2D-2 Quadtree / 2:1 balance
- 2D-3 Cut-cell polygon
- 2D-4 global edge/cell/neighbour topology
- 2D-5 small-cell stabilization
- 2D-6 quality / export
- 2D-V visualization

顶层仓库的 `CARTMESH_BUILD_2D` 可选接入仍应作为一个独立的小型构建集成改动处理；二维子项目本身目前可直接以 `cmake -S cartmesh2d ...` 独立构建，因此没有修改或回归三维构建路径。
