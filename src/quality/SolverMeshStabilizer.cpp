#include "cartmesh/quality/SolverMeshStabilizer.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace cartmesh {
namespace {

constexpr std::size_t no_id = std::numeric_limits<std::size_t>::max();

struct ConservedQuantities {
    double volume{};
    Vec3 first_moment{};
};

struct PatchMeasure {
    std::size_t face_count{};
    double area{};
};

using PatchKey = std::pair<bool, std::uint64_t>;
using PatchSignature = std::map<PatchKey, PatchMeasure>;

struct Candidate {
    std::size_t primary{};
    std::size_t secondary{};
    double shared_area{};
    double secondary_volume{};
    std::uint64_t secondary_stable_id{};
};

[[nodiscard]] double face_area(const OpenFoamMesh& mesh,
                               const OpenFoamFace& face) {
    if (face.point_ids.size() < 3U) return 0.0;
    const Vec3 origin = mesh.points.at(face.point_ids.front());
    Vec3 area_vector{};
    for (std::size_t index = 1; index + 1U < face.point_ids.size(); ++index) {
        area_vector = area_vector +
                      cross(mesh.points.at(face.point_ids[index]) - origin,
                            mesh.points.at(face.point_ids[index + 1U]) - origin) *
                          0.5;
    }
    return norm(area_vector);
}

[[nodiscard]] ConservedQuantities quantities(const MeshQualityReport& quality) {
    ConservedQuantities result;
    for (const auto& cell : quality.cell_metrics) {
        result.volume += cell.signed_volume;
        result.first_moment = result.first_moment + cell.first_moment;
    }
    return result;
}

[[nodiscard]] PatchSignature patch_signature(const OpenFoamMesh& mesh) {
    PatchSignature result;
    for (std::size_t face_id = mesh.internal_face_count;
         face_id < mesh.faces.size(); ++face_id) {
        const auto& face = mesh.faces[face_id];
        auto& value = result[{face.farfield, face.boundary_id}];
        ++value.face_count;
        value.area += face_area(mesh, face);
    }
    return result;
}

[[nodiscard]] bool close_scalar(double first, double second,
                                double relative_tolerance) {
    const double scale = std::max({1.0, std::abs(first), std::abs(second)});
    return std::abs(first - second) <= relative_tolerance * scale;
}

[[nodiscard]] bool close_vector(const Vec3& first, const Vec3& second,
                                double relative_tolerance) {
    const double scale = std::max({1.0, norm(first), norm(second)});
    return norm(first - second) <= relative_tolerance * scale;
}

[[nodiscard]] bool same_patch_signature(const PatchSignature& first,
                                        const PatchSignature& second,
                                        double tolerance) {
    if (first.size() != second.size()) return false;
    for (const auto& [key, value] : first) {
        const auto found = second.find(key);
        if (found == second.end() ||
            found->second.face_count != value.face_count ||
            !close_scalar(found->second.area, value.area, tolerance)) {
            return false;
        }
    }
    return true;
}

[[nodiscard]] bool cell_edges_are_two_manifold(const OpenFoamMesh& mesh) {
    using Edge = std::pair<std::size_t, std::size_t>;
    std::vector<std::map<Edge, std::size_t>> counts(mesh.cells.size());
    for (const auto& face : mesh.faces) {
        const auto count_edges = [&](std::size_t cell_id) {
            for (std::size_t index = 0; index < face.point_ids.size(); ++index) {
                const std::size_t first = face.point_ids[index];
                const std::size_t second =
                    face.point_ids[(index + 1U) % face.point_ids.size()];
                ++counts[cell_id][std::minmax(first, second)];
            }
        };
        count_edges(face.owner);
        if (face.internal()) count_edges(face.neighbour);
    }
    return std::all_of(counts.begin(), counts.end(), [](const auto& edges) {
        return !edges.empty() &&
               std::all_of(edges.begin(), edges.end(), [](const auto& entry) {
                   return entry.second == 2U;
               });
    });
}

[[nodiscard]] bool is_shape_or_topology_issue(MeshQualityIssueKind kind) {
    return kind != MeshQualityIssueKind::tiny_volume_fraction;
}

[[nodiscard]] bool is_non_star_issue(MeshQualityIssueKind kind) {
    return kind == MeshQualityIssueKind::wrong_face_pyramid ||
           kind == MeshQualityIssueKind::non_star_shaped_cell;
}

[[nodiscard]] std::map<MeshQualityIssueKind, std::size_t> issue_counts(
    const MeshQualityReport& report) {
    std::map<MeshQualityIssueKind, std::size_t> result;
    for (const auto& issue : report.issues) {
        if (is_shape_or_topology_issue(issue.kind)) ++result[issue.kind];
    }
    return result;
}

[[nodiscard]] std::vector<OpenFoamCellMember> source_members(
    const OpenFoamCellSource& source, double solver_volume) {
    if (!source.members.empty()) return source.members;
    const double fraction = std::max(
        source.source_volume_fraction, std::numeric_limits<double>::min());
    return {{source.background_cell_id, source.background_stable_id,
             source.component_id, source.local_piece_id,
             source.global_region_id, solver_volume / fraction}};
}

void canonicalize_members(std::vector<OpenFoamCellMember>& members) {
    std::sort(members.begin(), members.end(), [](const auto& first,
                                                 const auto& second) {
        return std::tie(first.background_stable_id, first.component_id,
                        first.local_piece_id, first.background_cell_id,
                        first.global_region_id) <
               std::tie(second.background_stable_id, second.component_id,
                        second.local_piece_id, second.background_cell_id,
                        second.global_region_id);
    });
    members.erase(std::unique(members.begin(), members.end(),
                              [](const auto& first, const auto& second) {
        return first.background_stable_id == second.background_stable_id &&
               first.component_id == second.component_id &&
               first.local_piece_id == second.local_piece_id &&
               first.background_cell_id == second.background_cell_id;
    }), members.end());
}

[[nodiscard]] std::vector<std::size_t> face_key(const OpenFoamFace& face) {
    std::vector<std::size_t> result = face.point_ids;
    std::sort(result.begin(), result.end());
    return result;
}

void compact_points(OpenFoamMesh& mesh) {
    std::vector<std::size_t> remap(mesh.points.size(), no_id);
    std::vector<Vec3> points;
    for (auto& face : mesh.faces) {
        for (auto& point_id : face.point_ids) {
            if (remap.at(point_id) == no_id) {
                remap[point_id] = points.size();
                points.push_back(mesh.points[point_id]);
            }
            point_id = remap[point_id];
        }
    }
    mesh.points = std::move(points);
}

[[nodiscard]] OpenFoamMesh merge_cells(
    const OpenFoamMesh& input, const MeshQualityReport& before,
    std::size_t first_cell, std::size_t second_cell) {
    if (first_cell == second_cell || first_cell >= input.cells.size() ||
        second_cell >= input.cells.size()) {
        throw std::invalid_argument("聚并单元 ID 无效");
    }
    const std::size_t keep = std::min(first_cell, second_cell);
    const std::size_t remove = std::max(first_cell, second_cell);
    OpenFoamMesh result;
    result.background_stable_id_kind = input.background_stable_id_kind;
    result.points = input.points;
    result.length_tolerance = input.length_tolerance;

    std::vector<std::size_t> cell_map(input.cells.size(), no_id);
    result.cells.reserve(input.cells.size() - 1U);
    for (std::size_t old = 0; old < input.cells.size(); ++old) {
        if (old == remove) continue;
        cell_map[old] = result.cells.size();
        result.cells.push_back(input.cells[old]);
    }
    cell_map[remove] = cell_map[keep];

    auto& merged = result.cells[cell_map[keep]];
    auto members = source_members(input.cells[keep],
                                  before.cell_metrics[keep].signed_volume);
    auto removed_members = source_members(
        input.cells[remove], before.cell_metrics[remove].signed_volume);
    members.insert(members.end(), removed_members.begin(), removed_members.end());
    canonicalize_members(members);
    double background_volume = 0.0;
    std::set<std::uint64_t> counted_backgrounds;
    for (const auto& member : members) {
        if (counted_backgrounds.insert(member.background_stable_id).second) {
            background_volume += member.background_volume;
        }
    }
    merged.members = std::move(members);
    const auto& representative = merged.members.front();
    merged.background_cell_id = representative.background_cell_id;
    merged.background_stable_id = representative.background_stable_id;
    merged.component_id = representative.component_id;
    merged.local_piece_id = representative.local_piece_id;
    merged.global_region_id = representative.global_region_id;
    merged.full_cartesian = false;
    const double merged_volume = before.cell_metrics[keep].signed_volume +
                                 before.cell_metrics[remove].signed_volume;
    merged.source_volume_fraction =
        background_volume > 0.0 ? merged_volume / background_volume : 0.0;

    std::vector<OpenFoamFace> internal_faces;
    std::vector<OpenFoamFace> boundary_faces;
    internal_faces.reserve(input.internal_face_count);
    boundary_faces.reserve(input.faces.size() - input.internal_face_count);
    for (const auto& old_face : input.faces) {
        if (old_face.internal() &&
            ((old_face.owner == keep && old_face.neighbour == remove) ||
             (old_face.owner == remove && old_face.neighbour == keep))) {
            continue;
        }
        OpenFoamFace face = old_face;
        face.owner = cell_map[old_face.owner];
        if (old_face.internal()) {
            face.neighbour = cell_map[old_face.neighbour];
            if (face.owner == face.neighbour) continue;
            if (face.owner > face.neighbour) {
                std::swap(face.owner, face.neighbour);
                std::reverse(face.point_ids.begin(), face.point_ids.end());
            }
            internal_faces.push_back(std::move(face));
        } else {
            boundary_faces.push_back(std::move(face));
        }
    }
    std::sort(internal_faces.begin(), internal_faces.end(),
              [](const auto& first, const auto& second) {
        return std::tuple(first.owner, first.neighbour, face_key(first),
                          first.point_ids) <
               std::tuple(second.owner, second.neighbour, face_key(second),
                          second.point_ids);
    });
    std::sort(boundary_faces.begin(), boundary_faces.end(),
              [](const auto& first, const auto& second) {
        if (first.farfield != second.farfield) return first.farfield;
        return std::tuple(first.boundary_id, first.owner, face_key(first),
                          first.point_ids) <
               std::tuple(second.boundary_id, second.owner, face_key(second),
                          second.point_ids);
    });
    result.internal_face_count = internal_faces.size();
    result.faces = std::move(internal_faces);
    result.faces.insert(result.faces.end(),
                        std::make_move_iterator(boundary_faces.begin()),
                        std::make_move_iterator(boundary_faces.end()));
    compact_points(result);
    return result;
}

[[nodiscard]] bool issue_counts_did_not_increase(
    const MeshQualityReport& before, const MeshQualityReport& after) {
    const auto first = issue_counts(before);
    const auto second = issue_counts(after);
    for (const auto& [kind, count] : second) {
        const auto found = first.find(kind);
        const std::size_t previous = found == first.end() ? 0U : found->second;
        if (count > previous) return false;
    }
    return true;
}

[[nodiscard]] bool merged_cell_has_non_star_issue(
    const MeshQualityReport& report, std::size_t merged_cell_id) {
    return std::any_of(report.issues.begin(), report.issues.end(),
                       [&](const auto& issue) {
        return issue.cell_id == merged_cell_id && is_non_star_issue(issue.kind);
    });
}

[[nodiscard]] std::vector<bool> cells_requiring_stabilization(
    const OpenFoamMesh& mesh, const MeshQualityReport& quality,
    double minimum_fraction) {
    std::vector<bool> result(mesh.cells.size(), false);
    for (std::size_t cell = 0; cell < mesh.cells.size(); ++cell) {
        result[cell] = mesh.cells[cell].source_volume_fraction < minimum_fraction;
    }
    for (const auto& issue : quality.issues) {
        if (issue.cell_id != no_id && is_shape_or_topology_issue(issue.kind)) {
            result[issue.cell_id] = true;
        }
    }
    return result;
}

void append_source_ids(std::vector<std::uint64_t>& ids,
                       const OpenFoamCellSource& source) {
    if (source.members.empty()) ids.push_back(source.background_stable_id);
    else {
        for (const auto& member : source.members) {
            ids.push_back(member.background_stable_id);
        }
    }
}

void sort_unique(std::vector<std::uint64_t>& values) {
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
}

void write_json_number(std::ostream& output, double value) {
    if (std::isfinite(value)) output << value;
    else output << "null";
}

void write_json_string(std::ostream& output, const std::string& value) {
    output << '"';
    for (const char character : value) {
        if (character == '"' || character == '\\') output << '\\';
        output << character;
    }
    output << '"';
}

} // namespace

const char* stabilization_action_name(StabilizationActionKind kind) noexcept {
    switch (kind) {
    case StabilizationActionKind::agglomerated: return "agglomerated";
    case StabilizationActionKind::rejected_non_star: return "rejected_non_star";
    case StabilizationActionKind::rejected_topology: return "rejected_topology";
    case StabilizationActionKind::rejected_conservation: return "rejected_conservation";
    case StabilizationActionKind::refinement_requested: return "refinement_requested";
    case StabilizationActionKind::conformal_refined: return "conformal_refined";
    case StabilizationActionKind::unresolved_max_level: return "unresolved_max_level";
    }
    return "unknown";
}

SolverMeshStabilizationResult stabilize_solver_mesh(
    const OpenFoamMesh& input, const MeshQualityThresholds& thresholds,
    const SolverMeshStabilizationOptions& options) {
    if (!(options.minimum_volume_fraction >= 0.0) ||
        !std::isfinite(options.minimum_volume_fraction) ||
        options.maximum_agglomeration_passes == 0U ||
        !(options.conservation_relative_tolerance > 0.0) ||
        !std::isfinite(options.conservation_relative_tolerance)) {
        throw std::invalid_argument("稳定化选项无效");
    }
    SolverMeshStabilizationResult result;
    result.mesh = input;
    result.report.initial_cell_count = input.cells.size();
    MeshQualityThresholds working_thresholds = thresholds;
    working_thresholds.minimum_volume_fraction = options.minimum_volume_fraction;
    auto quality = evaluate_solver_mesh_quality(result.mesh, working_thresholds);
    const auto initial_quantities = quantities(quality);
    const auto initial_patches = patch_signature(result.mesh);
    result.report.initial_volume = initial_quantities.volume;
    result.report.initial_first_moment = initial_quantities.first_moment;

    for (std::size_t pass = 0; pass < options.maximum_agglomeration_passes;
         ++pass) {
        const auto required = cells_requiring_stabilization(
            result.mesh, quality, options.minimum_volume_fraction);
        std::vector<Candidate> candidates;
        for (std::size_t face_id = 0; face_id < result.mesh.internal_face_count;
             ++face_id) {
            const auto& face = result.mesh.faces[face_id];
            const double area = face_area(result.mesh, face);
            if (!(area > thresholds.minimum_face_area)) continue;
            const auto add = [&](std::size_t primary, std::size_t secondary) {
                if (!required[primary] ||
                    result.mesh.cells[primary].global_region_id !=
                        result.mesh.cells[secondary].global_region_id) {
                    return;
                }
                candidates.push_back(
                    {primary, secondary, area,
                     quality.cell_metrics[secondary].signed_volume,
                     result.mesh.cells[secondary].background_stable_id});
            };
            add(face.owner, face.neighbour);
            add(face.neighbour, face.owner);
        }
        std::sort(candidates.begin(), candidates.end(), [](const auto& first,
                                                           const auto& second) {
            if (first.primary != second.primary) return first.primary < second.primary;
            if (first.shared_area != second.shared_area) {
                return first.shared_area > second.shared_area;
            }
            if (first.secondary_volume != second.secondary_volume) {
                return first.secondary_volume > second.secondary_volume;
            }
            if (first.secondary_stable_id != second.secondary_stable_id) {
                return first.secondary_stable_id < second.secondary_stable_id;
            }
            return first.secondary < second.secondary;
        });

        bool accepted = false;
        std::set<std::pair<std::size_t, std::size_t>> visited;
        for (const auto& candidate : candidates) {
            const std::pair<std::size_t, std::size_t> pair{
                std::min(candidate.primary, candidate.secondary),
                std::max(candidate.primary, candidate.secondary)};
            if (!visited.insert(pair).second) continue;
            const double before_fraction =
                result.mesh.cells[candidate.primary].source_volume_fraction;
            OpenFoamMesh trial = merge_cells(result.mesh, quality,
                                              candidate.primary,
                                              candidate.secondary);
            auto trial_quality =
                evaluate_solver_mesh_quality(trial, working_thresholds);
            const std::size_t merged_id = std::min(candidate.primary,
                                                   candidate.secondary);
            const auto trial_quantities = quantities(trial_quality);
            const bool conserved =
                close_scalar(initial_quantities.volume, trial_quantities.volume,
                             options.conservation_relative_tolerance) &&
                close_vector(initial_quantities.first_moment,
                             trial_quantities.first_moment,
                             options.conservation_relative_tolerance) &&
                same_patch_signature(initial_patches, patch_signature(trial),
                                     options.conservation_relative_tolerance);
            const bool manifold = cell_edges_are_two_manifold(trial);
            const bool non_star =
                merged_cell_has_non_star_issue(trial_quality, merged_id);
            const bool no_regression =
                trial_quality.topology_pass() && manifold &&
                issue_counts_did_not_increase(quality, trial_quality);
            const std::uint64_t first_action_id =
                result.mesh.cells[candidate.primary].background_stable_id;
            const std::uint64_t second_action_id =
                result.mesh.cells[candidate.secondary].background_stable_id;
            const std::pair<std::uint64_t, std::uint64_t> action_ids{
                std::min(first_action_id, second_action_id),
                std::max(first_action_id, second_action_id)};
            if (!conserved) {
                ++result.report.rejected_candidate_count;
                result.report.actions.push_back(
                    {StabilizationActionKind::rejected_conservation,
                     action_ids.first, action_ids.second, before_fraction,
                     trial.cells[merged_id].source_volume_fraction,
                     "volume, first moment, or boundary patch changed"});
                continue;
            }
            if (non_star) {
                ++result.report.rejected_candidate_count;
                result.report.actions.push_back(
                    {StabilizationActionKind::rejected_non_star,
                     action_ids.first, action_ids.second, before_fraction,
                     trial.cells[merged_id].source_volume_fraction,
                     "merged control volume has a non-positive face pyramid"});
                continue;
            }
            if (!no_regression) {
                ++result.report.rejected_candidate_count;
                result.report.actions.push_back(
                    {StabilizationActionKind::rejected_topology,
                     action_ids.first, action_ids.second, before_fraction,
                     trial.cells[merged_id].source_volume_fraction,
                     "topology, edge manifold, or quality regressed"});
                continue;
            }
            const double after_fraction =
                trial.cells[merged_id].source_volume_fraction;
            if (!(after_fraction > before_fraction)) {
                ++result.report.rejected_candidate_count;
                result.report.actions.push_back(
                    {StabilizationActionKind::rejected_conservation,
                     action_ids.first, action_ids.second, before_fraction,
                     after_fraction, "volume fraction did not improve"});
                continue;
            }
            result.mesh = std::move(trial);
            quality = std::move(trial_quality);
            ++result.report.agglomeration_count;
            result.report.actions.push_back(
                {StabilizationActionKind::agglomerated,
                 action_ids.first, action_ids.second, before_fraction,
                 after_fraction, "accepted conservative adjacency merge"});
            accepted = true;
            break;
        }
        if (!accepted) break;
    }

    quality = evaluate_solver_mesh_quality(result.mesh, working_thresholds);
    const auto unresolved = cells_requiring_stabilization(
        result.mesh, quality, options.minimum_volume_fraction);
    for (std::size_t cell = 0; cell < unresolved.size(); ++cell) {
        if (!unresolved[cell]) continue;
        append_source_ids(result.report.refinement_requested_stable_ids,
                          result.mesh.cells[cell]);
    }
    sort_unique(result.report.refinement_requested_stable_ids);
    result.report.unresolved_stable_ids =
        result.report.refinement_requested_stable_ids;
    for (const auto stable_id : result.report.refinement_requested_stable_ids) {
        result.report.actions.push_back(
            {StabilizationActionKind::refinement_requested, stable_id, 0U,
             0.0, 0.0, "agglomeration could not safely resolve source"});
    }
    const auto final_quantities = quantities(quality);
    result.report.final_cell_count = result.mesh.cells.size();
    result.report.final_volume = final_quantities.volume;
    result.report.final_first_moment = final_quantities.first_moment;
    result.report.conservation_pass =
        close_scalar(result.report.initial_volume, result.report.final_volume,
                     options.conservation_relative_tolerance) &&
        close_vector(result.report.initial_first_moment,
                     result.report.final_first_moment,
                     options.conservation_relative_tolerance) &&
        same_patch_signature(initial_patches, patch_signature(result.mesh),
                             options.conservation_relative_tolerance);
    return result;
}

AdaptiveStabilizationRefinementReport refine_stabilization_sources(
    LinearOctree& tree,
    const std::vector<std::uint64_t>& background_stable_ids) {
    AdaptiveStabilizationRefinementReport result;
    auto stable_ids = background_stable_ids;
    sort_unique(stable_ids);
    result.requested_count = stable_ids.size();
    for (const auto stable_id : stable_ids) {
        const auto coordinates = decode_octree_node(stable_id);
        if (!tree.find_leaf(stable_id).has_value() ||
            coordinates.level >= tree.maximum_level()) {
            result.unresolved_stable_ids.push_back(stable_id);
        }
    }
    // 回退必须是原子的：任一来源不可细分时，不先改动其他叶，
    // 避免调用方拿到“新 tree + 旧 Cut-cell mesh”的不一致状态。
    if (!result.unresolved_stable_ids.empty()) {
        sort_unique(result.unresolved_stable_ids);
        return result;
    }
    for (const auto stable_id : stable_ids) {
        if (!tree.refine_leaf(stable_id)) {
            throw std::runtime_error("已预检查的稳定化来源叶细分失败");
        }
        ++result.refined_leaf_count;
        result.refined_stable_ids.push_back(stable_id);
    }
    if (result.refined_leaf_count > 0U) {
        const auto balance = tree.balance_faces_2_to_1();
        result.balance_split_count = balance.split_count;
    }
    if (!tree.validate_partition() || !tree.check_face_balance().balanced) {
        throw std::runtime_error("稳定化局部细化破坏了八叉树分区或 2:1 平衡");
    }
    return result;
}

void write_solver_mesh_stabilization_json(
    const std::filesystem::path& path,
    const SolverMeshStabilizationReport& report) {
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path, std::ios::trunc);
    if (!output) throw std::runtime_error("无法写入网格稳定化报告");
    output << std::setprecision(17)
           << "{\n  \"schema\": \"cartmesh-solver-mesh-stabilization-v1\",\n"
           << "  \"pass\": " << (report.pass() ? "true" : "false") << ",\n"
           << "  \"initialCellCount\": " << report.initial_cell_count << ",\n"
           << "  \"finalCellCount\": " << report.final_cell_count << ",\n"
           << "  \"agglomerationCount\": " << report.agglomeration_count << ",\n"
           << "  \"rejectedCandidateCount\": "
           << report.rejected_candidate_count << ",\n"
           << "  \"conservationPass\": "
           << (report.conservation_pass ? "true" : "false") << ",\n"
           << "  \"initialVolume\": ";
    write_json_number(output, report.initial_volume);
    output << ",\n  \"finalVolume\": ";
    write_json_number(output, report.final_volume);
    output << ",\n  \"initialFirstMoment\": [";
    write_json_number(output, report.initial_first_moment.x);
    output << ',';
    write_json_number(output, report.initial_first_moment.y);
    output << ',';
    write_json_number(output, report.initial_first_moment.z);
    output << "],\n  \"finalFirstMoment\": [";
    write_json_number(output, report.final_first_moment.x);
    output << ',';
    write_json_number(output, report.final_first_moment.y);
    output << ',';
    write_json_number(output, report.final_first_moment.z);
    output << "],\n  \"refinementRequestedStableIds\": [";
    for (std::size_t index = 0;
         index < report.refinement_requested_stable_ids.size(); ++index) {
        if (index != 0U) output << ',';
        output << '"' << report.refinement_requested_stable_ids[index] << '"';
    }
    output << "],\n  \"unresolvedStableIds\": [";
    for (std::size_t index = 0; index < report.unresolved_stable_ids.size();
         ++index) {
        if (index != 0U) output << ',';
        output << '"' << report.unresolved_stable_ids[index] << '"';
    }
    output << "],\n  \"actions\": [";
    for (std::size_t index = 0; index < report.actions.size(); ++index) {
        const auto& action = report.actions[index];
        output << (index == 0U ? "\n" : ",\n")
               << "    {\"kind\":\"" << stabilization_action_name(action.kind)
               << "\",\"primaryStableId\":\"" << action.primary_stable_id
               << "\",\"secondaryStableId\":\"" << action.secondary_stable_id
               << "\",\"beforeVolumeFraction\":";
        write_json_number(output, action.before_volume_fraction);
        output << ",\"afterVolumeFraction\":";
        write_json_number(output, action.after_volume_fraction);
        output << ",\"reason\":";
        write_json_string(output, action.reason);
        output << '}';
    }
    output << (report.actions.empty() ? "" : "\n") << "  ]\n}\n";
}

} // namespace cartmesh
