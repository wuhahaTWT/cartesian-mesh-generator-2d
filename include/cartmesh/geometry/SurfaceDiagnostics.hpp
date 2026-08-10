#pragma once

#include "cartmesh/geometry/SurfaceMesh.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace cartmesh {

struct SurfaceDiagnosticLocation {
    std::uint64_t first_triangle{};
    std::uint64_t second_triangle{};
    Vec3 edge_start{};
    Vec3 edge_end{};
};

struct SurfaceDiagnosticVertex {
    Vec3 position{};
    std::uint64_t incident_triangle_count{};
};

struct SurfaceDiagnosticTrianglePair {
    std::uint64_t first_triangle{};
    std::uint64_t second_triangle{};
    Vec3 first_position{};
    Vec3 second_position{};
};

struct SurfaceDiagnosticComponent {
    std::uint64_t component_id{};
    std::uint64_t triangle_count{};
    double surface_area{};
    AABB bounds{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    double bounding_diagonal{};
    double signed_volume{};
    std::uint64_t nesting_depth{};
    int expected_orientation_sign{1};
    bool orientation_checked{};
    bool orientation_matches_nesting{};
    Vec3 sample_position{};
};

struct SurfaceDiagnostics {
    std::uint64_t triangle_count{};
    std::uint64_t unique_vertex_count{};
    std::uint64_t unique_edge_count{};
    std::uint64_t degenerate_triangle_count{};
    std::uint64_t duplicate_triangle_count{};
    std::uint64_t boundary_edge_count{};
    std::uint64_t non_manifold_edge_count{};
    std::uint64_t non_manifold_vertex_count{};
    std::uint64_t orientation_conflict_edge_count{};
    std::uint64_t connected_component_count{};
    std::uint64_t component_orientation_mismatch_count{};
    std::uint64_t overlapping_triangle_pair_count{};
    std::uint64_t self_intersection_pair_count{};
    std::uint64_t non_adjacent_contact_pair_count{};
    std::uint64_t small_component_count{};
    double suggested_length_tolerance{};
    double degenerate_area_tolerance{};
    double small_component_diagonal_threshold{};
    double minimum_component_bounding_diagonal{};
    double minimum_component_absolute_volume{};
    double signed_volume{};
    double material_volume{};
    bool closed{};
    bool manifold{};
    bool consistently_oriented{};
    std::vector<std::uint64_t> degenerate_triangle_examples;
    std::vector<std::uint64_t> duplicate_triangle_examples;
    std::vector<SurfaceDiagnosticLocation> boundary_edge_examples;
    std::vector<SurfaceDiagnosticLocation> non_manifold_edge_examples;
    std::vector<SurfaceDiagnosticLocation> orientation_conflict_examples;
    std::vector<SurfaceDiagnosticVertex> non_manifold_vertex_examples;
    std::vector<SurfaceDiagnosticTrianglePair> overlapping_triangle_examples;
    std::vector<SurfaceDiagnosticTrianglePair> self_intersection_examples;
    std::vector<SurfaceDiagnosticTrianglePair> non_adjacent_contact_examples;
    std::vector<SurfaceDiagnosticComponent> components;

    [[nodiscard]] bool valid_for_stage1_classification() const noexcept {
        return closed && manifold && consistently_oriented && degenerate_triangle_count == 0 &&
               duplicate_triangle_count == 0 && non_manifold_vertex_count == 0 &&
               overlapping_triangle_pair_count == 0 &&
               self_intersection_pair_count == 0 &&
               non_adjacent_contact_pair_count == 0;
    }
};

[[nodiscard]] SurfaceDiagnostics diagnose_surface(const SurfaceMesh& mesh,
                                                  std::size_t max_examples = 16);

} // 命名空间 cartmesh
