#include "cartmesh/quality/SolverMeshQuality.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <numbers>
#include <set>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace cartmesh {
namespace {

constexpr std::size_t no_id = std::numeric_limits<std::size_t>::max();

struct FaceGeometry {
    Vec3 area_vector{};
    Vec3 centroid{};
    double area{};
    double minimum_edge_length{std::numeric_limits<double>::infinity()};
    bool finite{true};
    bool concave{};
};

struct CellGeometry {
    Vec3 area_sum{};
    Vec3 centroid{};
    Vec3 volume_moment{};
    double surface_area{};
    double signed_volume{};
    bool finite{true};
};

[[nodiscard]] bool finite_nonnegative(double value) noexcept {
    return std::isfinite(value) && value >= 0.0;
}

void validate_thresholds(const MeshQualityThresholds& value) {
    if (!finite_nonnegative(value.minimum_face_area) ||
        !finite_nonnegative(value.minimum_edge_length) ||
        !finite_nonnegative(value.minimum_cell_volume) ||
        !finite_nonnegative(value.minimum_face_pyramid_volume) ||
        !finite_nonnegative(value.maximum_cell_closure_ratio) ||
        !finite_nonnegative(value.maximum_non_orthogonality_degrees) ||
        !finite_nonnegative(value.maximum_internal_skewness) ||
        !finite_nonnegative(value.maximum_boundary_skewness) ||
        !finite_nonnegative(value.minimum_volume_fraction)) {
        throw std::invalid_argument("网格质量阈值必须是非负有限数");
    }
}

[[nodiscard]] FaceGeometry face_geometry(const OpenFoamMesh& mesh,
                                         const OpenFoamFace& face) {
    FaceGeometry result;
    if (face.point_ids.size() < 3U) return result;
    std::vector<Vec3> points;
    points.reserve(face.point_ids.size());
    for (const std::size_t point_id : face.point_ids) {
        if (point_id >= mesh.points.size()) {
            throw std::invalid_argument("质量评估遇到越界 face point ID");
        }
        points.push_back(mesh.points[point_id]);
        result.finite = result.finite && is_finite(points.back());
    }
    const Vec3 origin = points.front();
    double centroid_weight = 0.0;
    for (std::size_t index = 1; index + 1U < points.size(); ++index) {
        const Vec3 triangle_area_vector =
            cross(points[index] - origin, points[index + 1U] - origin) * 0.5;
        const double triangle_area = norm(triangle_area_vector);
        result.area_vector = result.area_vector + triangle_area_vector;
        result.centroid = result.centroid +
                          (origin + points[index] + points[index + 1U]) *
                              (triangle_area / 3.0);
        centroid_weight += triangle_area;
    }
    result.area = norm(result.area_vector);
    if (centroid_weight > 0.0) result.centroid = result.centroid / centroid_weight;
    else result.centroid = origin;

    for (std::size_t index = 0; index < points.size(); ++index) {
        result.minimum_edge_length = std::min(
            result.minimum_edge_length,
            norm(points[(index + 1U) % points.size()] - points[index]));
    }
    if (result.area > 0.0) {
        const Vec3 normal = result.area_vector / result.area;
        double reference_sign = 0.0;
        for (std::size_t index = 0; index < points.size(); ++index) {
            const Vec3 first = points[(index + 1U) % points.size()] - points[index];
            const Vec3 second = points[(index + 2U) % points.size()] -
                                points[(index + 1U) % points.size()];
            const double turn = dot(cross(first, second), normal);
            if (std::abs(turn) <= 1.0e-14 * result.area) continue;
            if (reference_sign == 0.0) reference_sign = turn;
            else if (reference_sign * turn < 0.0) result.concave = true;
        }
    }
    result.finite = result.finite && is_finite(result.area_vector) &&
                    is_finite(result.centroid) && std::isfinite(result.area) &&
                    std::isfinite(result.minimum_edge_length);
    return result;
}

void accumulate_cell_face(const OpenFoamMesh& mesh, const OpenFoamFace& face,
                          double orientation, CellGeometry& cell) {
    if (face.point_ids.size() < 3U) return;
    const Vec3 first = mesh.points[face.point_ids.front()];
    for (std::size_t index = 1; index + 1U < face.point_ids.size(); ++index) {
        Vec3 second = mesh.points[face.point_ids[index]];
        Vec3 third = mesh.points[face.point_ids[index + 1U]];
        if (orientation < 0.0) std::swap(second, third);
        const double volume = dot(first, cross(second, third)) / 6.0;
        cell.signed_volume += volume;
        cell.volume_moment = cell.volume_moment +
                             (first + second + third) * (volume / 4.0);
    }
}

[[nodiscard]] std::string source_type(const OpenFoamCellSource& source) {
    return source.full_cartesian ? "full_cartesian" : "cut_polyhedron_piece";
}

void add_issue(MeshQualityReport& report, MeshQualityIssueKind kind,
               const OpenFoamMesh& mesh, std::size_t cell_id,
               std::size_t face_id, const Vec3& position,
               double measured, double threshold) {
    MeshQualityIssue issue;
    issue.kind = kind;
    issue.cell_id = cell_id;
    issue.face_id = face_id;
    issue.position = position;
    issue.measured_value = measured;
    issue.threshold = threshold;
    if (cell_id != no_id && cell_id < mesh.cells.size()) {
        const auto& source = mesh.cells[cell_id];
        issue.background_cell_id = source.background_cell_id;
        issue.background_stable_id = source.background_stable_id;
        issue.source_type = source_type(source);
    } else {
        issue.source_type = "unknown";
    }
    report.issues.push_back(std::move(issue));
}

[[nodiscard]] double clamp_unit(double value) noexcept {
    return std::max(-1.0, std::min(1.0, value));
}

[[nodiscard]] bool topology_issue(MeshQualityIssueKind kind) noexcept {
    switch (kind) {
    case MeshQualityIssueKind::non_finite_geometry:
    case MeshQualityIssueKind::face_too_few_vertices:
    case MeshQualityIssueKind::duplicate_face:
    case MeshQualityIssueKind::baffle_like_duplicate:
    case MeshQualityIssueKind::cell_not_closed:
    case MeshQualityIssueKind::non_positive_cell_volume:
        return true;
    default:
        return false;
    }
}

void write_json_id(std::ostream& output, std::size_t id) {
    if (id == no_id) output << "null";
    else output << id;
}

void write_json_number(std::ostream& output, double value) {
    if (std::isfinite(value)) output << value;
    else output << "null";
}

} // 匿名命名空间

const char* mesh_quality_issue_name(MeshQualityIssueKind kind) noexcept {
    switch (kind) {
    case MeshQualityIssueKind::non_finite_geometry: return "non_finite_geometry";
    case MeshQualityIssueKind::face_too_few_vertices: return "face_too_few_vertices";
    case MeshQualityIssueKind::zero_or_tiny_face: return "zero_or_tiny_face";
    case MeshQualityIssueKind::zero_or_tiny_edge: return "zero_or_tiny_edge";
    case MeshQualityIssueKind::duplicate_face: return "duplicate_face";
    case MeshQualityIssueKind::baffle_like_duplicate: return "baffle_like_duplicate";
    case MeshQualityIssueKind::cell_not_closed: return "cell_not_closed";
    case MeshQualityIssueKind::non_positive_cell_volume: return "non_positive_cell_volume";
    case MeshQualityIssueKind::wrong_face_pyramid: return "wrong_face_pyramid";
    case MeshQualityIssueKind::excessive_non_orthogonality: return "excessive_non_orthogonality";
    case MeshQualityIssueKind::excessive_skewness: return "excessive_skewness";
    case MeshQualityIssueKind::concave_face: return "concave_face";
    case MeshQualityIssueKind::non_star_shaped_cell: return "non_star_shaped_cell";
    case MeshQualityIssueKind::tiny_volume_fraction: return "tiny_volume_fraction";
    }
    return "unknown";
}

bool MeshQualityReport::topology_pass() const noexcept {
    return std::none_of(issues.begin(), issues.end(), [](const auto& issue) {
        return topology_issue(issue.kind);
    });
}

bool MeshQualityReport::quality_pass() const noexcept { return issues.empty(); }

MeshQualityReport evaluate_solver_mesh_quality(
    const OpenFoamMesh& mesh, const MeshQualityThresholds& thresholds) {
    validate_thresholds(thresholds);
    if (mesh.internal_face_count > mesh.faces.size()) {
        throw std::invalid_argument("质量评估的内部面计数越界");
    }
    MeshQualityReport report;
    report.thresholds = thresholds;
    report.summary.cell_count = mesh.cells.size();
    report.summary.face_count = mesh.faces.size();
    report.summary.internal_face_count = mesh.internal_face_count;

    std::vector<FaceGeometry> faces;
    faces.reserve(mesh.faces.size());
    std::vector<CellGeometry> cells(mesh.cells.size());
    std::map<std::vector<std::size_t>, std::size_t> canonical_faces;
    for (std::size_t face_id = 0; face_id < mesh.faces.size(); ++face_id) {
        const auto& face = mesh.faces[face_id];
        if (face.internal() != (face_id < mesh.internal_face_count) ||
            (face.internal() &&
             (face.owner == face.neighbour || face.owner > face.neighbour))) {
            throw std::invalid_argument(
                "质量评估要求内部面在前且满足 owner < neighbour");
        }
        if (face.owner >= mesh.cells.size() ||
            (face.internal() && face.neighbour >= mesh.cells.size())) {
            throw std::invalid_argument("质量评估遇到越界 owner/neighbour");
        }
        FaceGeometry geometry = face_geometry(mesh, face);
        faces.push_back(geometry);
        report.summary.minimum_face_area =
            std::min(report.summary.minimum_face_area, geometry.area);
        report.summary.minimum_edge_length = std::min(
            report.summary.minimum_edge_length, geometry.minimum_edge_length);
        if (face.point_ids.size() < 3U) {
            add_issue(report, MeshQualityIssueKind::face_too_few_vertices,
                      mesh, face.owner, face_id, geometry.centroid,
                      static_cast<double>(face.point_ids.size()), 3.0);
        }
        if (!geometry.finite) {
            add_issue(report, MeshQualityIssueKind::non_finite_geometry,
                      mesh, face.owner, face_id, geometry.centroid,
                      geometry.area, 0.0);
        }
        if (!(geometry.area > thresholds.minimum_face_area)) {
            add_issue(report, MeshQualityIssueKind::zero_or_tiny_face,
                      mesh, face.owner, face_id, geometry.centroid,
                      geometry.area, thresholds.minimum_face_area);
        }
        if (!(geometry.minimum_edge_length > thresholds.minimum_edge_length)) {
            add_issue(report, MeshQualityIssueKind::zero_or_tiny_edge,
                      mesh, face.owner, face_id, geometry.centroid,
                      geometry.minimum_edge_length,
                      thresholds.minimum_edge_length);
        }
        if (geometry.concave) {
            add_issue(report, MeshQualityIssueKind::concave_face,
                      mesh, face.owner, face_id, geometry.centroid, 1.0, 0.0);
        }
        std::vector<std::size_t> key = face.point_ids;
        std::sort(key.begin(), key.end());
        const auto [duplicate, inserted] = canonical_faces.emplace(key, face_id);
        if (!inserted) {
            const bool baffle = !face.internal() &&
                                !mesh.faces[duplicate->second].internal();
            add_issue(report,
                      baffle ? MeshQualityIssueKind::baffle_like_duplicate
                             : MeshQualityIssueKind::duplicate_face,
                      mesh, face.owner, face_id, geometry.centroid,
                      static_cast<double>(duplicate->second), 0.0);
        }
        cells[face.owner].area_sum =
            cells[face.owner].area_sum + geometry.area_vector;
        cells[face.owner].surface_area += geometry.area;
        accumulate_cell_face(mesh, face, 1.0, cells[face.owner]);
        if (face.internal()) {
            cells[face.neighbour].area_sum =
                cells[face.neighbour].area_sum - geometry.area_vector;
            cells[face.neighbour].surface_area += geometry.area;
            accumulate_cell_face(mesh, face, -1.0, cells[face.neighbour]);
        }
    }

    for (std::size_t cell_id = 0; cell_id < cells.size(); ++cell_id) {
        auto& geometry = cells[cell_id];
        geometry.finite = is_finite(geometry.area_sum) &&
                          std::isfinite(geometry.surface_area) &&
                          std::isfinite(geometry.signed_volume) &&
                          is_finite(geometry.volume_moment);
        if (geometry.signed_volume != 0.0 &&
            std::isfinite(geometry.signed_volume)) {
            geometry.centroid = geometry.volume_moment / geometry.signed_volume;
        }
        const double closure_ratio = geometry.surface_area > 0.0
                                         ? norm(geometry.area_sum) /
                                               geometry.surface_area
                                         : std::numeric_limits<double>::infinity();
        report.summary.maximum_cell_closure_ratio = std::max(
            report.summary.maximum_cell_closure_ratio, closure_ratio);
        report.summary.minimum_cell_volume = std::min(
            report.summary.minimum_cell_volume, geometry.signed_volume);
        if (!geometry.finite || !is_finite(geometry.centroid)) {
            add_issue(report, MeshQualityIssueKind::non_finite_geometry,
                      mesh, cell_id, no_id, geometry.centroid,
                      geometry.signed_volume, 0.0);
        }
        if (!(closure_ratio <= thresholds.maximum_cell_closure_ratio)) {
            add_issue(report, MeshQualityIssueKind::cell_not_closed,
                      mesh, cell_id, no_id, geometry.centroid,
                      closure_ratio, thresholds.maximum_cell_closure_ratio);
        }
        if (!(geometry.signed_volume > thresholds.minimum_cell_volume)) {
            add_issue(report, MeshQualityIssueKind::non_positive_cell_volume,
                      mesh, cell_id, no_id, geometry.centroid,
                      geometry.signed_volume, thresholds.minimum_cell_volume);
        }
        const double volume_fraction = mesh.cells[cell_id].source_volume_fraction;
        report.summary.minimum_volume_fraction =
            std::min(report.summary.minimum_volume_fraction, volume_fraction);
        if (!(volume_fraction >= thresholds.minimum_volume_fraction)) {
            add_issue(report, MeshQualityIssueKind::tiny_volume_fraction,
                      mesh, cell_id, no_id, geometry.centroid,
                      volume_fraction, thresholds.minimum_volume_fraction);
        }
    }

    std::vector<bool> non_star(cells.size(), false);
    for (std::size_t face_id = 0; face_id < mesh.faces.size(); ++face_id) {
        const auto& face = mesh.faces[face_id];
        const auto& geometry = faces[face_id];
        if (!(geometry.area > 0.0) || !geometry.finite) continue;
        const auto check_pyramid = [&](std::size_t cell_id, double orientation) {
            const double pyramid_volume = orientation *
                dot(geometry.area_vector,
                    geometry.centroid - cells[cell_id].centroid) / 3.0;
            report.summary.minimum_face_pyramid_volume = std::min(
                report.summary.minimum_face_pyramid_volume, pyramid_volume);
            if (!(pyramid_volume > thresholds.minimum_face_pyramid_volume)) {
                add_issue(report, MeshQualityIssueKind::wrong_face_pyramid,
                          mesh, cell_id, face_id, geometry.centroid,
                          pyramid_volume,
                          thresholds.minimum_face_pyramid_volume);
                non_star[cell_id] = true;
            }
        };
        check_pyramid(face.owner, 1.0);
        if (face.internal()) check_pyramid(face.neighbour, -1.0);

        double skewness = 0.0;
        if (face.internal()) {
            const Vec3 delta = cells[face.neighbour].centroid -
                               cells[face.owner].centroid;
            const double denominator = dot(geometry.area_vector, delta);
            const double delta_length = norm(delta);
            const double cosine = denominator /
                (geometry.area * std::max(delta_length,
                    std::numeric_limits<double>::min()));
            const double angle = std::acos(clamp_unit(cosine)) *
                                 180.0 / std::numbers::pi_v<double>;
            report.summary.maximum_non_orthogonality_degrees = std::max(
                report.summary.maximum_non_orthogonality_degrees, angle);
            if (!(angle <= thresholds.maximum_non_orthogonality_degrees)) {
                add_issue(report,
                          MeshQualityIssueKind::excessive_non_orthogonality,
                          mesh, face.owner, face_id, geometry.centroid,
                          angle,
                          thresholds.maximum_non_orthogonality_degrees);
            }
            if (std::abs(denominator) >
                std::numeric_limits<double>::epsilon() * geometry.area *
                    std::max(delta_length, 1.0)) {
                const double parameter = dot(
                    geometry.area_vector,
                    geometry.centroid - cells[face.owner].centroid) /
                    denominator;
                const Vec3 intersection = cells[face.owner].centroid +
                                          delta * parameter;
                skewness = norm(geometry.centroid - intersection) /
                           std::max(delta_length,
                                    std::numeric_limits<double>::min());
            } else {
                skewness = std::numeric_limits<double>::infinity();
            }
            if (!(skewness <= thresholds.maximum_internal_skewness)) {
                add_issue(report, MeshQualityIssueKind::excessive_skewness,
                          mesh, face.owner, face_id, geometry.centroid,
                          skewness, thresholds.maximum_internal_skewness);
            }
        } else {
            const Vec3 normal = geometry.area_vector / geometry.area;
            const Vec3 delta = geometry.centroid - cells[face.owner].centroid;
            const double normal_distance = dot(delta, normal);
            const Vec3 projected = cells[face.owner].centroid +
                                   normal * normal_distance;
            skewness = norm(geometry.centroid - projected) /
                       std::max(std::abs(normal_distance),
                                std::numeric_limits<double>::min());
            if (!(skewness <= thresholds.maximum_boundary_skewness)) {
                add_issue(report, MeshQualityIssueKind::excessive_skewness,
                          mesh, face.owner, face_id, geometry.centroid,
                          skewness, thresholds.maximum_boundary_skewness);
            }
        }
        report.summary.maximum_skewness =
            std::max(report.summary.maximum_skewness, skewness);
    }
    for (std::size_t cell_id = 0; cell_id < non_star.size(); ++cell_id) {
        if (non_star[cell_id]) {
            add_issue(report, MeshQualityIssueKind::non_star_shaped_cell,
                      mesh, cell_id, no_id, cells[cell_id].centroid, 1.0, 0.0);
        }
    }

    std::sort(report.issues.begin(), report.issues.end(), [](const auto& lhs,
                                                             const auto& rhs) {
        return std::tie(lhs.kind, lhs.face_id, lhs.cell_id,
                        lhs.background_stable_id) <
               std::tie(rhs.kind, rhs.face_id, rhs.cell_id,
                        rhs.background_stable_id);
    });
    std::set<std::size_t> failing_cells;
    std::set<std::size_t> failing_faces;
    for (const auto& issue : report.issues) {
        if (issue.cell_id != no_id) failing_cells.insert(issue.cell_id);
        if (issue.face_id != no_id) failing_faces.insert(issue.face_id);
    }
    report.summary.issue_count = report.issues.size();
    report.summary.failing_cell_count = failing_cells.size();
    report.summary.failing_face_count = failing_faces.size();
    return report;
}

void write_solver_mesh_quality_json(const std::filesystem::path& path,
                                    const MeshQualityReport& report) {
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path, std::ios::trunc);
    if (!output) throw std::runtime_error("无法写入原生网格质量报告");
    output << std::setprecision(17);
    std::map<MeshQualityIssueKind, std::size_t> counts;
    for (const auto& issue : report.issues) ++counts[issue.kind];
    output << "{\n  \"schema\": \"cartmesh-solver-mesh-quality-v1\",\n"
           << "  \"topologyPass\": " << (report.topology_pass() ? "true" : "false")
           << ",\n  \"qualityPass\": " << (report.quality_pass() ? "true" : "false")
           << ",\n  \"thresholds\": {\n"
           << "    \"minimumFaceArea\": " << report.thresholds.minimum_face_area << ",\n"
           << "    \"minimumEdgeLength\": " << report.thresholds.minimum_edge_length << ",\n"
           << "    \"minimumCellVolume\": " << report.thresholds.minimum_cell_volume << ",\n"
           << "    \"minimumFacePyramidVolume\": " << report.thresholds.minimum_face_pyramid_volume << ",\n"
           << "    \"maximumCellClosureRatio\": " << report.thresholds.maximum_cell_closure_ratio << ",\n"
           << "    \"maximumNonOrthogonalityDegrees\": " << report.thresholds.maximum_non_orthogonality_degrees << ",\n"
           << "    \"maximumInternalSkewness\": " << report.thresholds.maximum_internal_skewness << ",\n"
           << "    \"maximumBoundarySkewness\": " << report.thresholds.maximum_boundary_skewness << ",\n"
           << "    \"minimumVolumeFraction\": " << report.thresholds.minimum_volume_fraction << "\n"
           << "  },\n  \"summary\": {\n"
           << "    \"cellCount\": " << report.summary.cell_count << ",\n"
           << "    \"faceCount\": " << report.summary.face_count << ",\n"
           << "    \"internalFaceCount\": " << report.summary.internal_face_count << ",\n"
           << "    \"issueCount\": " << report.summary.issue_count << ",\n"
           << "    \"failingCellCount\": " << report.summary.failing_cell_count << ",\n"
           << "    \"failingFaceCount\": " << report.summary.failing_face_count << ",\n"
           << "    \"minimumCellVolume\": ";
    write_json_number(output, report.summary.minimum_cell_volume);
    output << ",\n    \"minimumFaceArea\": ";
    write_json_number(output, report.summary.minimum_face_area);
    output << ",\n    \"minimumEdgeLength\": ";
    write_json_number(output, report.summary.minimum_edge_length);
    output << ",\n    \"minimumVolumeFraction\": ";
    write_json_number(output, report.summary.minimum_volume_fraction);
    output << ",\n    \"maximumCellClosureRatio\": ";
    write_json_number(output, report.summary.maximum_cell_closure_ratio);
    output << ",\n    \"maximumNonOrthogonalityDegrees\": ";
    write_json_number(output,
                      report.summary.maximum_non_orthogonality_degrees);
    output << ",\n    \"maximumSkewness\": ";
    write_json_number(output, report.summary.maximum_skewness);
    output << ",\n    \"minimumFacePyramidVolume\": ";
    write_json_number(output, report.summary.minimum_face_pyramid_volume);
    output << "\n  },\n  \"issueCounts\": {";
    bool first = true;
    for (std::uint8_t value = 0;
         value <= static_cast<std::uint8_t>(MeshQualityIssueKind::tiny_volume_fraction);
         ++value) {
        const auto kind = static_cast<MeshQualityIssueKind>(value);
        if (!first) output << ',';
        first = false;
        output << "\n    \"" << mesh_quality_issue_name(kind) << "\": "
               << counts[kind];
    }
    output << "\n  },\n  \"issues\": [";
    for (std::size_t index = 0; index < report.issues.size(); ++index) {
        const auto& issue = report.issues[index];
        output << (index == 0U ? "\n" : ",\n")
               << "    {\"kind\":\"" << mesh_quality_issue_name(issue.kind)
               << "\",\"cellId\":";
        write_json_id(output, issue.cell_id);
        output << ",\"faceId\":";
        write_json_id(output, issue.face_id);
        output << ",\"backgroundCellId\":" << issue.background_cell_id
               << ",\"backgroundStableId\":\"" << issue.background_stable_id
               << "\",\"sourceType\":\"" << issue.source_type
               << "\",\"position\":[";
        write_json_number(output, issue.position.x);
        output << ',';
        write_json_number(output, issue.position.y);
        output << ',';
        write_json_number(output, issue.position.z);
        output << "],\"measured\":";
        write_json_number(output, issue.measured_value);
        output << ",\"threshold\":";
        write_json_number(output, issue.threshold);
        output << '}';
    }
    output << (report.issues.empty() ? "" : "\n") << "  ]\n}\n";
}

} // 命名空间 cartmesh
