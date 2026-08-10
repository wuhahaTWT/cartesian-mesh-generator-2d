#pragma once

#include "cartmesh/cutcell/ConvexPolyhedron.hpp"
#include "cartmesh/geometry/SurfaceMesh.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace cartmesh {

struct ConvexSurfaceCutResult {
    ConvexPolyhedron solid_polyhedron;
    PolyhedronGeometry solid_geometry;
    double solid_volume_fraction{};
    double fluid_volume{};
    double fluid_volume_fraction{};
    Vec3 fluid_centroid{};
    double embedded_boundary_area{};
    double volume_conservation_residual{};
    std::vector<ConvexPolyhedron> fluid_decomposition_polyhedra;
    std::vector<PolyhedronGeometry> fluid_decomposition_geometries;
    std::vector<std::size_t> fluid_piece_component_ids;
    std::size_t fluid_component_count{};
    bool cut{};
};

// 为单连通的封闭凸 STL 构建确定性半空间集合。整体反向壳会在工作
// 表面中统一翻向；非凸或多分量表面不会进入这一算法。
class ConvexSurfaceCutter {
  public:
    explicit ConvexSurfaceCutter(const SurfaceMesh& surface,
                                 std::uint64_t boundary_id = 0,
                                 double length_tolerance = 0.0);

    [[nodiscard]] std::size_t plane_count() const noexcept { return planes_.size(); }
    [[nodiscard]] bool input_orientation_reversed() const noexcept {
        return input_orientation_reversed_;
    }
    [[nodiscard]] double length_tolerance() const noexcept { return length_tolerance_; }
    [[nodiscard]] std::uint64_t boundary_id() const noexcept { return boundary_id_; }
    [[nodiscard]] ConvexSurfaceCutResult cut_box(const AABB& box) const;

  private:
    std::vector<OrientedHalfSpace> planes_;
    double length_tolerance_{};
    bool input_orientation_reversed_{};
    std::uint64_t boundary_id_{};
};

} // 命名空间 cartmesh
