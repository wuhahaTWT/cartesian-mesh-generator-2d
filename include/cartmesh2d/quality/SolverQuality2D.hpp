#pragma once

#include "cartmesh2d/topology/Topology2D.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace cartmesh2d {

struct SolverQualityPolicy2D {
    // These are production repair targets, not merely catastrophic-failure
    // limits. buildSolverTopology2D keeps agglomerating/repartitioning while
    // any target is violated, so accepted meshes retain numerical margin
    // instead of stopping immediately at the OpenFOAM-style safety boundary.
    double maxNonOrthogonalityDeg = 65.0;
    double maxInternalSkewness = 3.5;
    double maxBoundarySkewness = 3.0;
    double maxConcavityDeg = 60.0;
    double maxCellAspect = 200.0;
    double minInteriorAngleDeg = 1.0;
    double minFaceLength = 1.01e-10;
    double minFaceWeight = 0.08;
    double minVolumeRatio = 0.02;
};

enum class SolverQualityIssueCode2D {
    InvalidTopology,
    InvalidCell,
    ShortFace,
    ExcessiveNonOrthogonality,
    ExcessiveSkewness,
    ExcessiveBoundarySkewness,
    ExcessiveConcavity,
    ExcessiveAspect,
    SmallInteriorAngle,
    LowFaceWeight,
    LowVolumeRatio
};

struct SolverQualityIssue2D {
    SolverQualityIssueCode2D code;
    std::size_t cellId = 0;
    std::size_t edgeId = 0;
    double measured = 0.0;
    double limit = 0.0;
    std::string message;
};

struct SolverQualityReport2D {
    SolverQualityPolicy2D policy;
    double maxNonOrthogonalityDeg = 0.0;
    double maxInternalSkewness = 0.0;
    double maxBoundarySkewness = 0.0;
    double maxConcavityDeg = 0.0;
    double maxCellAspect = 0.0;
    double minInteriorAngleDeg = 180.0;
    double minFaceLength = 0.0;
    double minFaceWeight = 1.0;
    double minVolumeRatio = 1.0;
    double minCompactness = 1.0;
    std::vector<SolverQualityIssue2D> issues;

    [[nodiscard]] bool valid() const noexcept { return issues.empty(); }
};

[[nodiscard]] SolverQualityReport2D evaluateSolverQuality2D(
    const TopologyMesh2D& topology,
    const SolverQualityPolicy2D& policy = {},
    const TolerancePolicy& tol = {});

[[nodiscard]] std::string solverQualityReportToJson(
    const SolverQualityReport2D& report, int indentSpaces = 2);

} // namespace cartmesh2d
