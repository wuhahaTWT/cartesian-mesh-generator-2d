#pragma once

#include "cartmesh/grid/LinearOctree.hpp"
#include "cartmesh/spatial/TriangleBvh.hpp"

#include <cstdint>
#include <optional>
#include <vector>

namespace cartmesh {

struct DistanceRefinementBand {
    double maximum_distance{};
    std::uint8_t target_level{};
};

struct CurvatureRefinementRule {
    double minimum_normal_angle_degrees{};
    double neighborhood_cell_diagonals{1.5};
    std::uint8_t target_level{};
};

struct GapRefinementRule {
    double maximum_search_distance{};
    double maximum_opposing_normal_dot{-0.75};
    double minimum_facing_dot{0.5};
    std::uint32_t minimum_cells_across_gap{4};
};

struct BoxRefinementRegion {
    AABB bounds;
    std::uint8_t target_level{};
};

struct SphereRefinementRegion {
    Vec3 center{};
    double radius{};
    std::uint8_t target_level{};
};

struct CylinderRefinementRegion {
    Vec3 first_axis_point{};
    Vec3 second_axis_point{};
    double radius{};
    std::uint8_t target_level{};
};

struct OctreeRefinementConfiguration {
    std::optional<std::uint8_t> surface_target_level;
    std::vector<DistanceRefinementBand> distance_bands;
    std::optional<CurvatureRefinementRule> curvature;
    std::optional<GapRefinementRule> gap;
    std::vector<BoxRefinementRegion> boxes;
    std::vector<SphereRefinementRegion> spheres;
    std::vector<CylinderRefinementRegion> cylinders;
    bool enforce_face_2_to_1_balance{true};
};

struct OctreeAdaptationStatistics {
    OctreeRefinementStatistics rule_refinement;
    OctreeBalanceStatistics balance;
    std::uint64_t surface_rule_hits{};
    std::uint64_t distance_rule_hits{};
    std::uint64_t curvature_rule_hits{};
    std::uint64_t gap_rule_hits{};
    std::uint64_t gap_resolution_failure_count{};
    std::uint32_t maximum_required_gap_level{};
    std::uint64_t user_region_rule_hits{};
};

class OctreeRefinementEngine {
  public:
    OctreeRefinementEngine(OctreeRefinementConfiguration configuration,
                           const TriangleBvh* bvh = nullptr);

    [[nodiscard]] std::uint8_t desired_level(const LinearOctree& tree,
                                             OctreeNodeCode code,
                                             const AABB& bounds) const;
    [[nodiscard]] OctreeAdaptationStatistics apply(LinearOctree& tree) const;

  private:
    struct DetectedGap {
        double width{};
        Vec3 normal{};
    };

    void validate_for_tree(const LinearOctree& tree) const;
    [[nodiscard]] bool curvature_rule_matches(const AABB& bounds) const;
    [[nodiscard]] std::optional<DetectedGap> detect_gap(const AABB& bounds) const;
    [[nodiscard]] std::uint32_t required_gap_level(const LinearOctree& tree,
                                                   const DetectedGap& gap) const;

    OctreeRefinementConfiguration configuration_;
    const TriangleBvh* bvh_{};
};

} // 命名空间 cartmesh
