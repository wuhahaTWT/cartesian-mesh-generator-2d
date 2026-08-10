#pragma once

#include "cartmesh/geometry/Triangle.hpp"

namespace cartmesh {

// 闭集语义：三角形只要接触盒子的面、边或顶点，也算相交。
// 实现使用完整的 13 轴分离轴定理，而不是三角形包围盒重叠近似。
[[nodiscard]] bool triangle_intersects_aabb(const Triangle& triangle,
                                            const AABB& box) noexcept;

} // 命名空间 cartmesh
