#pragma once

#include "cartmesh/cutcell/ConvexCutCellMesh.hpp"
#include "cartmesh/spatial/TriangleBvh.hpp"
#include "cartmesh/grid/LinearOctree.hpp"

#include <array>
#include <cstdint>
#include <optional>
#include <vector>

namespace cartmesh {

struct SolidCartesianFaceRegion {
    double area{};
    Vec3 centroid{};
    // 每个环的绕序已经包含有向四面体链的正负系数。
    std::vector<std::vector<Vec3>> oriented_boundary_loops;
};

struct TriangulatedSurfaceCellCut {
    double solid_volume{};
    Vec3 solid_centroid{};
    double fluid_volume{};
    double fluid_volume_fraction{};
    Vec3 fluid_centroid{};
    std::array<SolidCartesianFaceRegion, 6> solid_cartesian_faces{};
    std::vector<EmbeddedBoundaryFaceGeometry> embedded_boundary_faces;
    double embedded_boundary_area{};
    double volume_conservation_residual{};
    bool cut{};
};

struct TriangulatedSurfaceBoundaryClip {
    std::vector<EmbeddedBoundaryFaceGeometry> faces;
    double total_area{};
};

// 面向单连通、封闭、定向一致三角网格的通用单元积分器。它不要求表面凸：
// 闭曲面的有向三角片与固定参考点形成一个有向四面体 3-chain，逐四面体
// 与 AABB 精确裁剪后线性累计体积、一次矩和 Cartesian 面占据矩。
class TriangulatedSurfaceCutter {
  public:
    explicit TriangulatedSurfaceCutter(const SurfaceMesh& surface,
                                       std::uint64_t boundary_id = 0,
                                       double length_tolerance = 0.0);
    TriangulatedSurfaceCutter(const SurfaceMesh& surface,
                             std::vector<std::uint64_t> triangle_boundary_ids,
                             double length_tolerance = 0.0);

    [[nodiscard]] const SurfaceMesh& oriented_surface() const noexcept {
        return oriented_surface_;
    }
    [[nodiscard]] const TriangleBvh& bvh() const noexcept { return bvh_; }
    [[nodiscard]] std::uint64_t boundary_id() const noexcept { return boundary_id_; }
    [[nodiscard]] double length_tolerance() const noexcept { return length_tolerance_; }
    [[nodiscard]] TriangulatedSurfaceCellCut cut_box(const AABB& box) const;
    // 只裁剪当前 AABB 内的真实 STL wall 多边形，不扫描
    // 全局有向四面体链。阶段 6 的稀疏 Cut-cell 路径用它
    // 与局部平面 arrangement 组合得到完整控制体。
    [[nodiscard]] TriangulatedSurfaceBoundaryClip
    clip_boundary_faces(const AABB& box) const;

  private:
    struct TetraContribution {
        ConvexPolyhedron tetrahedron;
        AABB bounds;
        double coefficient{};
    };

    SurfaceMesh oriented_surface_;
    TriangleBvh bvh_;
    std::vector<TetraContribution> tetrahedra_;
    std::vector<std::uint64_t> triangle_boundary_ids_;
    Vec3 reference_{};
    std::uint64_t boundary_id_{};
    double length_tolerance_{};
};

// 将通用三角曲面单元积分接到均匀网格与 cell-face-neighbor 拓扑。
[[nodiscard]] ConvexCutCellMesh build_triangulated_cut_cell_mesh(
    const UniformCartesianGrid& grid,
    const TriangulatedSurfaceCutter& cutter,
    double geometric_tolerance = 0.0);

[[nodiscard]] ConvexCutCellMesh build_triangulated_cut_cell_mesh(
    const LinearOctree& tree,
    const TriangulatedSurfaceCutter& cutter,
    double geometric_tolerance = 0.0);

struct IncrementalCutCellBuildStatistics {
    std::uint64_t old_leaf_count{};
    std::uint64_t new_leaf_count{};
    std::uint64_t reused_leaf_count{};
    std::uint64_t rebuilt_leaf_count{};
    std::uint64_t reused_fluid_cell_count{};
    std::uint64_t reused_solid_cell_count{};
    std::uint64_t rebuilt_cut_cell_count{};
    double geometry_reuse_fraction{};
};

struct IncrementalCutCellBuildResult {
    ConvexCutCellMesh mesh;
    IncrementalCutCellBuildStatistics statistics;
    // 与 new_tree 叶顺序一致；1 表示本次重新执行了几何分类/切割。
    std::vector<std::uint8_t> rebuilt_leaf_mask;
};

// 对稳定叶码且位于几何影响范围之外的单元复用旧 Cut-cell 几何；新叶、
// 删除/粗化后的替代叶和受影响叶重新切割。面邻接与全局 region 在新树上
// 确定性重建，不能复用旧的临时 leaf index。
[[nodiscard]] IncrementalCutCellBuildResult
build_incremental_triangulated_cut_cell_mesh(
    const LinearOctree& old_tree, const ConvexCutCellMesh& old_mesh,
    const LinearOctree& new_tree,
    const TriangulatedSurfaceCutter& new_cutter,
    const std::optional<AABB>& affected_bounds,
    double geometric_tolerance = 0.0);

} // 命名空间 cartmesh
