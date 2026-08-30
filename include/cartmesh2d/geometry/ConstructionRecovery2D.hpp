#pragma once

#include "cartmesh2d/geometry/IntersectionRegistry2D.hpp"

#include <optional>
#include <string>
#include <vector>

namespace cartmesh2d {

// A recovery candidate is only a construction proposal. Rephase/resample may
// not enter committed geometry until area, feature and hard-quality checks all
// pass. The ordinal is stable and defines the deterministic final tie-break.
struct ConstructionRecoveryCandidate2D {
    std::size_t ordinal = 0;
    ConstructionDecision2D decision = ConstructionDecision2D::Rejected;
    std::optional<Point2D> sourcePoint;
    std::optional<double> sourceParameter;
    std::optional<double> gridPhaseOffset;
    bool preservesSourceSegment = false;
    bool requiresAreaFeatureQualityCheck = true;
    std::string reason;
};

struct ConstructionRecoveryCandidateEvaluation2D {
    std::size_t ordinal = 0;
    bool areaConserved = false;
    bool featureCompatible = false;
    bool hardQualityPass = false;
    double qualityRank = 0.0;
};

[[nodiscard]] std::vector<ConstructionRecoveryCandidate2D>
planConstructionRecoveryCandidates(
    const ConstructionRecoveryRequest2D& request,
    std::size_t maximumGridPhaseCandidates = 4U);

// Selects only fully evaluated admissible candidates. qualityRank is minimized;
// decision order and ordinal make equal-rank selection deterministic.
[[nodiscard]] std::optional<ConstructionRecoveryCandidate2D>
selectConstructionRecoveryCandidate(
    const std::vector<ConstructionRecoveryCandidate2D>& candidates,
    const std::vector<ConstructionRecoveryCandidateEvaluation2D>& evaluations);

} // namespace cartmesh2d
