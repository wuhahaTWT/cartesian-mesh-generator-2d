#pragma once

#include "cartmesh/geometry/Triangle.hpp"

namespace cartmesh {

enum class TriangleTriangleRelation {
    disjoint,
    boundary_contact,
    proper_intersection,
    coplanar_area_overlap,
    indeterminate,
};

// 闭三角形语义。合法网格邻接是否应被忽略由拓扑诊断层根据共享顶点/边决定。
[[nodiscard]] TriangleTriangleRelation classify_triangle_triangle(
    const Triangle& first, const Triangle& second,
    double length_tolerance = 0.0) noexcept;

} // 命名空间 cartmesh
