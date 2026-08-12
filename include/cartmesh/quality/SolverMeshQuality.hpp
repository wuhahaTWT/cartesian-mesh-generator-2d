#pragma once

#include "cartmesh/io/OpenFoamWriter.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <string>
#include <vector>

namespace cartmesh {

enum class MeshQualityIssueKind : std::uint8_t {
    non_finite_geometry,
    face_too_few_vertices,
    zero_or_tiny_face,
    zero_or_tiny_edge,
    duplicate_face,
    baffle_like_duplicate,
    cell_not_closed,
    non_positive_cell_volume,
    wrong_face_pyramid,
    excessive_non_orthogonality,
    excessive_skewness,
    concave_face,
    non_star_shaped_cell,
    tiny_volume_fraction,
};

[[nodiscard]] const char* mesh_quality_issue_name(
    MeshQualityIssueKind kind) noexcept;

struct MeshQualityThresholds {
    double minimum_face_area{1.0e-20};
    double minimum_edge_length{1.0e-12};
    double minimum_cell_volume{1.0e-18};
    double minimum_face_pyramid_volume{1.0e-18};
    double maximum_cell_closure_ratio{1.0e-10};
    double maximum_non_orthogonality_degrees{70.0};
    double maximum_internal_skewness{4.0};
    double maximum_boundary_skewness{20.0};
    double minimum_volume_fraction{0.01};
};

struct MeshQualityIssue {
    MeshQualityIssueKind kind{};
    std::size_t cell_id{std::numeric_limits<std::size_t>::max()};
    std::size_t face_id{std::numeric_limits<std::size_t>::max()};
    std::uint64_t background_cell_id{};
    std::uint64_t background_stable_id{};
    std::string source_type;
    Vec3 position{};
    double measured_value{};
    double threshold{};
};

struct MeshQualitySummary {
    std::size_t cell_count{};
    std::size_t face_count{};
    std::size_t internal_face_count{};
    std::size_t issue_count{};
    std::size_t failing_cell_count{};
    std::size_t failing_face_count{};
    double minimum_cell_volume{std::numeric_limits<double>::infinity()};
    double minimum_face_area{std::numeric_limits<double>::infinity()};
    double minimum_edge_length{std::numeric_limits<double>::infinity()};
    double minimum_volume_fraction{1.0};
    double maximum_cell_closure_ratio{};
    double maximum_non_orthogonality_degrees{};
    double maximum_skewness{};
    double minimum_face_pyramid_volume{std::numeric_limits<double>::infinity()};
};

struct MeshQualityReport {
    MeshQualityThresholds thresholds;
    MeshQualitySummary summary;
    std::vector<MeshQualityIssue> issues;

    [[nodiscard]] bool topology_pass() const noexcept;
    [[nodiscard]] bool quality_pass() const noexcept;
};

[[nodiscard]] MeshQualityReport evaluate_solver_mesh_quality(
    const OpenFoamMesh& mesh,
    const MeshQualityThresholds& thresholds = {});

void write_solver_mesh_quality_json(
    const std::filesystem::path& path,
    const MeshQualityReport& report);

} // 命名空间 cartmesh
