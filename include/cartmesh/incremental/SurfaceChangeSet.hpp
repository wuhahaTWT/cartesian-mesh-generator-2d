#pragma once

#include "cartmesh/geometry/SurfaceMesh.hpp"

#include <cstdint>
#include <optional>
#include <vector>

namespace cartmesh {

struct SurfaceChangeSet {
    std::vector<std::uint64_t> old_triangle_ids;
    std::vector<std::uint64_t> new_triangle_ids;
    std::optional<AABB> bounds;

    [[nodiscard]] bool empty() const noexcept {
        return old_triangle_ids.empty() && new_triangle_ids.empty();
    }
};

// 三角形键与三角片顺序、顶点起点和绕序无关。当前阶段使用输入坐标的
// 精确位模式；几何容差不是静默“未变化”的许可。
[[nodiscard]] SurfaceChangeSet detect_surface_changes(
    const SurfaceMesh& old_surface, const SurfaceMesh& new_surface);

// 将变化包围盒保守扩展并裁剪到固定计算域。无变化时返回 nullopt。
[[nodiscard]] std::optional<AABB> conservative_affected_bounds(
    const SurfaceChangeSet& changes, const AABB& domain, double margin);

} // namespace cartmesh
