#include "cartmesh/cutcell/ConvexCutCellMesh.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cartmesh {
namespace {

struct DirectedEdge {
    std::size_t first{};
    std::size_t second{};
};

[[nodiscard]] std::size_t find_or_append(std::vector<Vec3>& points,
                                         const Vec3& point,
                                         double tolerance) {
    for (std::size_t index = 0; index < points.size(); ++index) {
        if (norm(points[index] - point) <= tolerance) return index;
    }
    points.push_back(point);
    return points.size() - 1U;
}

void append_loop(std::vector<Vec3>& points, std::vector<DirectedEdge>& edges,
                 const std::vector<Vec3>& loop, double tolerance) {
    if (loop.size() < 3) return;
    std::vector<std::size_t> ids;
    ids.reserve(loop.size());
    for (const auto& point : loop) {
        ids.push_back(find_or_append(points, point, tolerance));
    }
    for (std::size_t edge = 0; edge < ids.size(); ++edge) {
        const std::size_t first = ids[edge];
        const std::size_t second = ids[(edge + 1U) % ids.size()];
        if (first != second) edges.push_back({first, second});
    }
}

} // 匿名命名空间

BoundaryEdgeClosure analyze_boundary_edge_closure(
    const FluidCellGeometry& cell, double length_tolerance) {
    if (length_tolerance < 0.0 || !std::isfinite(length_tolerance)) {
        throw std::invalid_argument("Cut-cell 边链闭合容差必须为非负有限数");
    }
    std::vector<Vec3> points;
    std::vector<DirectedEdge> edges;
    for (const auto& face : cell.cartesian_faces) {
        for (const auto& loop : face.oriented_boundary_loops) {
            append_loop(points, edges, loop, length_tolerance);
        }
    }
    for (const auto& face : cell.embedded_boundary_faces) {
        append_loop(points, edges, face.vertices, length_tolerance);
    }

    std::map<std::pair<std::size_t, std::size_t>, std::int64_t> balances;
    for (const auto edge : edges) {
        const Vec3 first = points[edge.first];
        const Vec3 second = points[edge.second];
        const Vec3 direction = second - first;
        const double length_squared = dot(direction, direction);
        if (!(length_squared > 0.0)) continue;
        std::vector<std::pair<double, std::size_t>> split_points = {
            {0.0, edge.first}, {1.0, edge.second}};
        for (std::size_t point_id = 0; point_id < points.size(); ++point_id) {
            if (point_id == edge.first || point_id == edge.second) continue;
            const double parameter = dot(points[point_id] - first, direction) /
                                     length_squared;
            if (parameter <= 0.0 || parameter >= 1.0) continue;
            const Vec3 projection = first + direction * parameter;
            if (norm(points[point_id] - projection) <= length_tolerance) {
                split_points.push_back({parameter, point_id});
            }
        }
        std::sort(split_points.begin(), split_points.end(),
                  [](const auto& lhs, const auto& rhs) {
                      return lhs.first < rhs.first ||
                             (lhs.first == rhs.first && lhs.second < rhs.second);
                  });
        for (std::size_t segment = 0; segment + 1U < split_points.size(); ++segment) {
            const std::size_t a = split_points[segment].second;
            const std::size_t b = split_points[segment + 1U].second;
            if (a == b) continue;
            const auto key = std::minmax(a, b);
            balances[{key.first, key.second}] += a < b ? 1 : -1;
        }
    }
    BoundaryEdgeClosure result;
    for (const auto& [edge, balance] : balances) {
        static_cast<void>(edge);
        if (balance != 0) ++result.imbalanced_edge_count;
    }
    result.closed = !balances.empty() && result.imbalanced_edge_count == 0;
    return result;
}

} // 命名空间 cartmesh
