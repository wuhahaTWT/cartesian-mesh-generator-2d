#pragma once

#include "cartmesh2d/topology/Topology2D.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace cartmesh2d {

struct SolverQualityPolicy2D {
    double maxNonOrthogonalityDeg = 70.0;
    double maxInternalSkewness = 4.0;
    double maxBoundarySkewness = 4.0;
    double maxConcavityDeg = 80.0;
    double maxCellAspect = 1000.0;
    double minInteriorAngleDeg = 0.5;
    double minFaceLength = 1.01e-10;
    double minFaceWeight = 0.05;
    double minVolumeRatio = 0.01;
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

// Per-entity metric kernels.  These are the single authoritative definition of
// every solver metric: both the global evaluator below and the patch-local
// evaluator call them, so a local decision cannot use a different formula.
struct SolverCellMetrics2D {
    bool valid = false;
    Point2D centroid;
    double area = 0.0;
    double perimeter = 0.0;
    double minInteriorAngleDeg = 360.0;
    double maxConcavityDeg = 0.0;
    double hydraulicAspect = 0.0;
    double compactness = 0.0;
};

[[nodiscard]] SolverCellMetrics2D evaluateSolverCellMetrics2D(
    const Polygon2D& polygon, const TolerancePolicy& tol = {});

struct SolverInternalFaceMetrics2D {
    // False when the owner/neighbour connector or the face normal degenerates.
    // The global evaluator skips every internal metric in that case.
    bool orientationValid = false;
    double length = 0.0;
    double nonOrthogonalityDeg = 0.0;
    double skewness = 0.0;
    double faceWeight = 0.0;
    double volumeRatio = 0.0;
};

[[nodiscard]] SolverInternalFaceMetrics2D evaluateSolverInternalFaceMetrics2D(
    const Point2D& a, const Point2D& b, const Point2D& ownerCentroid,
    const Point2D& neighbourCentroid, double ownerArea, double neighbourArea,
    const TolerancePolicy& tol = {});

[[nodiscard]] double evaluateSolverBoundaryFaceSkewness2D(
    const Point2D& a, const Point2D& b, const Point2D& ownerCentroid);

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

// Monotonic count of full-mesh solver-quality evaluations in this process.
// Instrumentation only, mirroring globalTopologyBuildCount2D().
[[nodiscard]] std::size_t solverQualityEvaluationCount2D() noexcept;

[[nodiscard]] std::string solverQualityReportToJson(
    const SolverQualityReport2D& report, int indentSpaces = 2);

} // namespace cartmesh2d
