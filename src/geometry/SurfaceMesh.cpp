#include "cartmesh/geometry/SurfaceMesh.hpp"

#include <algorithm>
#include <stdexcept>

namespace cartmesh {

SurfaceMesh::SurfaceMesh(std::vector<Triangle> triangles, SurfaceFormat format, std::string name)
    : triangles_(std::move(triangles)),
      bounds_({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}),
      format_(format),
      name_(std::move(name)) {
    if (triangles_.empty()) {
        throw std::invalid_argument("表面网格必须至少包含一个三角形");
    }
    Vec3 minimum = triangles_.front().vertices().front();
    Vec3 maximum = minimum;
    for (const auto& triangle : triangles_) {
        for (const auto& vertex : triangle.vertices()) {
            minimum.x = std::min(minimum.x, vertex.x);
            minimum.y = std::min(minimum.y, vertex.y);
            minimum.z = std::min(minimum.z, vertex.z);
            maximum.x = std::max(maximum.x, vertex.x);
            maximum.y = std::max(maximum.y, vertex.y);
            maximum.z = std::max(maximum.z, vertex.z);
        }
    }
    bounds_ = AABB(minimum, maximum);
}

} // 命名空间 cartmesh
