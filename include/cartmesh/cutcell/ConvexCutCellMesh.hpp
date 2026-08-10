#pragma once

#include "cartmesh/cutcell/ConvexSurfaceCutter.hpp"
#include "cartmesh/grid/UniformCartesianGrid.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace cartmesh {

// 局部 Cartesian 面编号与 make_box_polyhedron 的 source_id 一致：
// 0=x-, 1=x+, 2=y-, 3=y+, 4=z-, 5=z+。
struct CartesianFaceAperture {
    double area{};
    double area_fraction{};
    Vec3 centroid{};
    Vec3 outward_normal{};
    // 与 outward_normal 一致的有向边界环。第一环是完整 Cartesian 面，
    // 后续反向环扣除固体占据区；面积为零时两者会完全抵消。
    std::vector<std::vector<Vec3>> oriented_boundary_loops;
};

struct EmbeddedBoundaryFaceGeometry {
    std::uint64_t boundary_id{};
    double area{};
    Vec3 centroid{};
    // 流体控制体的外法向，即指向固体内部。
    Vec3 outward_normal{};
    // 按流体外法向定向的边界多边形。
    std::vector<Vec3> vertices;
};

struct FluidPolyhedronPiece {
    ConvexPolyhedron polyhedron;
    PolyhedronGeometry geometry;
    std::size_t component_id{};
    std::uint64_t global_region_id{};
};

struct FluidCellGeometry {
    std::uint64_t background_cell_id{};
    double volume{};
    double volume_fraction{};
    Vec3 centroid{};
    std::array<CartesianFaceAperture, 6> cartesian_faces{};
    std::vector<EmbeddedBoundaryFaceGeometry> embedded_boundary_faces;
    // 仅体积 Cut-cell 显式保存；普通 Cartesian 流体单元继续隐式表示。
    std::vector<FluidPolyhedronPiece> fluid_polyhedron_pieces;
    Vec3 oriented_area_vector_sum{};
    double area_vector_closure_residual{};
    std::uint64_t boundary_edge_imbalance_count{};
    bool boundary_edge_closed{};
    std::size_t fluid_piece_count{};
    std::size_t fluid_component_count{};
    std::uint64_t discarded_numerical_piece_count{};
    double discarded_numerical_piece_volume{};
    bool component_analysis_resolved{true};
    std::vector<std::uint64_t> fluid_component_region_ids;
    bool cut{};
};

struct FluidFaceConnection {
    std::uint64_t first_background_cell_id{};
    std::uint64_t second_background_cell_id{};
    std::size_t first_fluid_cell_index{};
    std::size_t second_fluid_cell_index{};
    std::uint8_t first_local_face{};
    std::uint8_t second_local_face{};
    double area{};
    Vec3 centroid{};
    // 从 first 指向 second。
    Vec3 normal{};
    double area_mismatch{};
    double centroid_mismatch{};
    // 以共享面中心为原点的面积一阶矩差。这是保守量；
    // 对几乎封死的微小开口，派生质心本身可以是病态量。
    double first_moment_mismatch{};
};

struct FluidComponentFaceConnection {
    std::uint64_t first_background_cell_id{};
    std::uint64_t second_background_cell_id{};
    std::size_t first_fluid_cell_index{};
    std::size_t second_fluid_cell_index{};
    std::size_t first_component_id{};
    std::size_t second_component_id{};
    std::uint64_t global_region_id{};
    double area{};
    Vec3 centroid{};
    Vec3 normal{};
};

struct ConvexCutCellMesh {
    std::vector<FluidCellGeometry> fluid_cells;
    std::vector<FluidFaceConnection> internal_faces;
    std::vector<FluidComponentFaceConnection> component_internal_faces;
    std::vector<double> global_fluid_region_volumes;
    std::uint64_t global_fluid_region_count{};
    std::uint64_t full_fluid_cell_count{};
    std::uint64_t full_solid_cell_count{};
    std::uint64_t cut_cell_count{};
    std::uint64_t boundary_cell_count{};
    double total_fluid_volume{};
    double total_embedded_boundary_area{};
    double maximum_cell_area_closure_residual{};
    double maximum_shared_face_area_mismatch{};
    double maximum_shared_face_centroid_mismatch{};
    double maximum_shared_face_first_moment_mismatch{};
    std::uint64_t nonclosed_cell_count{};
    std::uint64_t negative_volume_cell_count{};
    std::uint64_t disconnected_fluid_cell_count{};
    std::uint64_t component_analysis_pending_cell_count{};
    std::uint64_t classification_conflict_count{};
    std::uint64_t shared_face_mismatch_count{};
    std::uint64_t discarded_numerical_piece_count{};
    double discarded_numerical_piece_volume{};
};

// 依据显式开口多边形跨单元连接局部 fluid component，并赋予确定性全局 region ID。
void assign_global_fluid_regions(ConvexCutCellMesh& mesh,
                                 double area_tolerance = 0.0);

struct SmallCutCellEntry {
    std::uint64_t background_cell_id{};
    double volume_fraction{};
    Vec3 centroid{};
    std::vector<std::uint64_t> boundary_ids;
};

struct SmallCutCellStatistics {
    double threshold{};
    double minimum_cut_cell_volume_fraction{1.0};
    std::vector<SmallCutCellEntry> cells;
};

// 只检测并记录小流体体积分数，不在阶段三自动合并或修改控制体。
[[nodiscard]] SmallCutCellStatistics analyze_small_cut_cells(
    const ConvexCutCellMesh& mesh, double threshold);

// 在均匀 Cartesian 网格上构建凸 STL 外部的流体控制体几何和内部面拓扑。
// 当前返回精确的体积/质心、面开口矩和嵌入边界矩；显式非凸开口多边形
// 的分片顶点将在阶段三后续接口中加入。
[[nodiscard]] ConvexCutCellMesh build_convex_cut_cell_mesh(
    const UniformCartesianGrid& grid, const ConvexSurfaceCutter& cutter,
    double geometric_tolerance = 0.0);

struct BoundaryEdgeClosure {
    std::uint64_t imbalanced_edge_count{};
    bool closed{};
};

[[nodiscard]] BoundaryEdgeClosure analyze_boundary_edge_closure(
    const FluidCellGeometry& cell, double length_tolerance);

} // 命名空间 cartmesh
