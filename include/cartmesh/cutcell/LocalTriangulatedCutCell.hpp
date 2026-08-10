#pragma once

#include "cartmesh/classify/SurfaceClassifier.hpp"
#include "cartmesh/cutcell/TriangulatedSurfaceCutter.hpp"

#include <cstdint>

namespace cartmesh {

struct LocalTriangulatedCutCellResult {
    bool has_fluid{};
    bool classification_conflict{};
    bool component_analysis_pending{};
    FluidCellGeometry cell;
    double embedded_boundary_area{};
    std::uint64_t discarded_numerical_piece_count{};
    double discarded_numerical_piece_volume{};
};

// 只为已知与表面相交的背景盒构建显式 Cut-cell。
// 几何来自当前盒内的局部平面 arrangement，普通单元
// 不经过本对象，因此内存与表面单元数而非全域单元数成比。
class LocalTriangulatedCutCellBuilder {
  public:
    LocalTriangulatedCutCellBuilder(const TriangulatedSurfaceCutter& cutter,
                                    const AABB& domain,
                                    const Vec3& representative_spacing,
                                    double geometric_tolerance = 0.0,
                                    bool retain_small_positive_pieces = false);

    [[nodiscard]] LocalTriangulatedCutCellResult
    build(std::uint64_t background_cell_id, const AABB& box) const;

    [[nodiscard]] double length_tolerance() const noexcept {
        return length_tolerance_;
    }
    [[nodiscard]] double area_tolerance() const noexcept {
        return area_tolerance_;
    }
    [[nodiscard]] double topology_length_tolerance() const noexcept {
        return topology_length_tolerance_;
    }
    [[nodiscard]] double closure_tolerance() const noexcept {
        return closure_tolerance_;
    }

  private:
    const TriangulatedSurfaceCutter& cutter_;
    SurfaceClassifier classifier_;
    double length_tolerance_{};
    double area_tolerance_{};
    double topology_length_tolerance_{};
    double closure_tolerance_{};
    bool retain_small_positive_pieces_{};
};

} // 命名空间 cartmesh
