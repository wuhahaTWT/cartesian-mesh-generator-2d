#pragma once

#include "cartmesh/cutcell/LocalTriangulatedCutCell.hpp"
#include "cartmesh/grid/UniformCartesianGrid.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

namespace cartmesh {

enum class CompactCellState : std::uint8_t {
    solid = 0,
    full_fluid = 1,
    explicit_surface = 2,
    conflict = 3,
};

struct CompactUniformBuildTimings {
    double surface_rasterization_seconds{};
    double connected_classification_seconds{};
    double local_cut_cell_seconds{};
    double topology_and_regions_seconds{};
    double total_seconds{};
};

struct CompactFaceFailure {
    std::uint64_t first_background_cell_id{};
    std::uint64_t second_background_cell_id{};
    std::uint8_t first_local_face{};
    double area_mismatch{};
    double first_moment_mismatch{};
};

struct CompactCellFailure {
    std::uint64_t background_cell_id{};
    double area_closure_residual{};
    std::uint64_t boundary_edge_imbalance_count{};
    std::uint64_t fluid_piece_count{};
};

struct NumericallySealedCartesianFace {
    std::size_t fluid_piece_index{};
    std::size_t polyhedron_face_index{};
    std::uint8_t local_face{};
    std::uint64_t boundary_id{};
    double area{};
};

struct CompactCutCellRecord {
    std::uint64_t background_cell_id{};
    FluidCellGeometry geometry;
    std::uint64_t component_node_offset{};
    std::vector<NumericallySealedCartesianFace> numerically_sealed_faces;
};

struct CompactUniformCutCellMesh {
    std::vector<std::uint8_t> cell_states;
    // 只对 full_fluid 有效；映射到 full_component_region_ids。
    std::vector<std::uint32_t> full_component_labels;
    std::vector<std::uint64_t> full_component_region_ids;
    std::vector<CompactCutCellRecord> explicit_cells;
    std::vector<double> global_fluid_region_volumes;

    std::uint64_t background_cell_count{};
    std::uint64_t surface_candidate_cell_count{};
    std::uint64_t full_fluid_cell_count{};
    std::uint64_t full_solid_cell_count{};
    std::uint64_t explicit_surface_cell_count{};
    std::uint64_t cut_cell_count{};
    std::uint64_t boundary_cell_count{};
    std::uint64_t solver_cell_count{};
    std::uint64_t full_fluid_component_count{};
    std::uint64_t global_fluid_region_count{};
    std::uint64_t internal_background_connection_count{};
    std::uint64_t direct_fluid_solid_face_count{};

    double total_fluid_volume{};
    double total_embedded_boundary_area{};
    double minimum_cut_cell_volume_fraction{1.0};
    std::uint64_t small_cut_cell_count{};
    double small_cut_cell_threshold{0.01};

    double maximum_cell_area_closure_residual{};
    double maximum_shared_face_area_mismatch{};
    double maximum_shared_face_first_moment_mismatch{};
    std::uint64_t nonclosed_cell_count{};
    std::uint64_t negative_volume_cell_count{};
    std::uint64_t component_analysis_pending_cell_count{};
    std::uint64_t classification_conflict_count{};
    std::uint64_t shared_face_mismatch_count{};
    std::uint64_t aggregate_boundary_edge_imbalance_cell_count{};
    std::uint64_t explicit_piece_topology_failure_count{};
    std::uint64_t discarded_numerical_piece_count{};
    double discarded_numerical_piece_volume{};
    std::uint64_t numerically_sealed_cartesian_face_count{};
    double numerically_sealed_cartesian_face_area{};

    // 只保留前 64 个稳定 ID 样本，用于建立最小回归案例；
    // 完整失败数仍以上述 count 为准。
    std::vector<CompactCellFailure> nonclosed_cell_samples;
    std::vector<std::uint64_t> pending_cell_samples;
    std::vector<CompactFaceFailure> shared_face_failure_samples;
    std::vector<CompactFaceFailure> direct_fluid_solid_face_samples;

    std::uint64_t compact_storage_bytes{};
    std::uint64_t result_hash_fnv1a64{};
    CompactUniformBuildTimings timings;

    [[nodiscard]] CompactCellState state(std::uint64_t background_id) const;
    [[nodiscard]] const CompactCutCellRecord*
    find_explicit_cell(std::uint64_t background_id) const noexcept;
    [[nodiscard]] CompactCutCellRecord*
    find_explicit_cell(std::uint64_t background_id) noexcept;
    [[nodiscard]] bool invariants_pass() const noexcept;
};

struct CompactUniformBuildOptions {
    double geometric_tolerance{};
    double small_cut_cell_threshold{0.01};
};

[[nodiscard]] CompactUniformCutCellMesh
build_compact_uniform_cut_cell_mesh(
    const UniformCartesianGrid& grid,
    const TriangulatedSurfaceCutter& cutter,
    const CompactUniformBuildOptions& options = {});

[[nodiscard]] const char* compact_cell_state_name(
    CompactCellState state) noexcept;

} // 命名空间 cartmesh
