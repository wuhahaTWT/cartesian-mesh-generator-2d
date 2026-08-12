#pragma once

#include "cartmesh/grid/LinearOctree.hpp"
#include "cartmesh/quality/SolverMeshQuality.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace cartmesh {

enum class StabilizationActionKind : std::uint8_t {
    agglomerated,
    rejected_non_star,
    rejected_topology,
    rejected_conservation,
    refinement_requested,
    conformal_refined,
    unresolved_max_level,
};

[[nodiscard]] const char* stabilization_action_name(
    StabilizationActionKind kind) noexcept;

struct SolverMeshStabilizationOptions {
    double minimum_volume_fraction{0.01};
    std::size_t maximum_agglomeration_passes{100000};
    double conservation_relative_tolerance{1.0e-12};
};

struct StabilizationAction {
    StabilizationActionKind kind{};
    std::uint64_t primary_stable_id{};
    std::uint64_t secondary_stable_id{};
    double before_volume_fraction{};
    double after_volume_fraction{};
    std::string reason;
};

struct SolverMeshStabilizationReport {
    std::size_t initial_cell_count{};
    std::size_t final_cell_count{};
    std::size_t agglomeration_count{};
    std::size_t rejected_candidate_count{};
    double initial_volume{};
    double final_volume{};
    Vec3 initial_first_moment{};
    Vec3 final_first_moment{};
    bool conservation_pass{true};
    std::vector<std::uint64_t> refinement_requested_stable_ids;
    std::vector<std::uint64_t> unresolved_stable_ids;
    std::vector<StabilizationAction> actions;

    [[nodiscard]] bool pass() const noexcept {
        return conservation_pass && unresolved_stable_ids.empty();
    }
};

struct SolverMeshStabilizationResult {
    OpenFoamMesh mesh;
    SolverMeshStabilizationReport report;
};

struct AdaptiveStabilizationRefinementReport {
    std::size_t requested_count{};
    std::size_t refined_leaf_count{};
    std::size_t balance_split_count{};
    std::vector<std::uint64_t> refined_stable_ids;
    std::vector<std::uint64_t> unresolved_stable_ids;

    [[nodiscard]] bool pass() const noexcept {
        return unresolved_stable_ids.empty();
    }
};

[[nodiscard]] SolverMeshStabilizationResult stabilize_solver_mesh(
    const OpenFoamMesh& mesh,
    const MeshQualityThresholds& thresholds = {},
    const SolverMeshStabilizationOptions& options = {});

[[nodiscard]] AdaptiveStabilizationRefinementReport
refine_stabilization_sources(
    LinearOctree& tree,
    const std::vector<std::uint64_t>& background_stable_ids);

void write_solver_mesh_stabilization_json(
    const std::filesystem::path& path,
    const SolverMeshStabilizationReport& report);

} // namespace cartmesh
