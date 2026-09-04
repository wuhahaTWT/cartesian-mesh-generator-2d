#pragma once

#include "cartmesh2d/quality/SolverQuality2D.hpp"

#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace cartmesh2d {

enum class QualityCellType2D {
    Cartesian,
    RemainderCut,
    Transition,
    Termination,
    BoundaryLayer,
    Unknown
};

enum class QualityContractStatus2D { Pass, Warn, Fail };
enum class QualityContractLevel2D { Preferred, Hard };

struct QualityLimit2D {
    double preferred = 0.0;
    double hard = 0.0;
    bool lowerBound = false;
};

struct OrdinaryCellQualityLimits2D {
    QualityLimit2D nonOrthogonalityDeg{55.0,65.0,false};
    QualityLimit2D skewness{2.0,3.0,false};
    QualityLimit2D faceWeight{0.15,0.10,true};
    // Gated, and kept at OpenFOAM's checkMesh form for parity with it.  Be aware of
    // what it cannot mean on a 2:1-graded cut-cell mesh: a cut cell of area fraction
    // alpha against a coarser full neighbour scores alpha/4, and agglomeration only
    // guarantees alpha >= small-alpha, so this limit is not reachable by improving
    // the mesh.  tests/quality_reachability_test.cpp records that bound explicitly.
    QualityLimit2D volumeRatio{0.10,0.05,true};
    // The grading-aware form: each area divided by the total background box area of
    // its sources, so a level difference cancels and two full cells score 1.0 at any
    // level gap.  Reported but *not* gated yet - see the comment at its check site in
    // QualityContract2D.cpp.  The limit is already certified reachable so the switch
    // needs no new threshold work, only a change of evaluation target.
    QualityLimit2D backgroundVolumeRatio{0.20,0.10,true};
    QualityLimit2D minimumInteriorAngleDeg{20.0,10.0,true};
    QualityLimit2D hydraulicAspect{20.0,50.0,false};
    QualityLimit2D faceOverLocalBackgroundH{0.03,0.01,true};
    QualityLimit2D faceOverSqrtOwnerArea{0.03,0.01,true};
    QualityLimit2D faceOverSqrtNeighbourArea{0.03,0.01,true};
};

struct QualityContract2D {
    OrdinaryCellQualityLimits2D cartesian;
    OrdinaryCellQualityLimits2D remainderCut;
    OrdinaryCellQualityLimits2D transition;
    OrdinaryCellQualityLimits2D termination;
};

struct QualityCellMetadata2D {
    QualityCellType2D type = QualityCellType2D::Unknown;
    // Size of the *finest* background cell this solver cell descends from.  Correct
    // for face-length ratios, where the finest contributor sets the scale a face
    // should be compared against.
    double localBackgroundH = 0.0;
    // Total background box area of every source this solver cell descends from.
    // Distinct from localBackgroundH^2 on purpose: an agglomerated cell spanning a
    // level-7 and a level-8 leaf has one finest size but two boxes of area, and an
    // area *fraction* is only bounded by 1 when divided by the total.  Zero means
    // "not known", and background_volume_ratio then skips the face.
    double backgroundArea = 0.0;
    std::size_t sourceId = 0;
    std::optional<std::size_t> layerIndex;
    std::optional<std::size_t> wallSegment;
};

struct QualityEntity2D {
    std::size_t cellId = 0;
    std::optional<std::size_t> edgeId;
    Point2D coordinates;
    QualityCellType2D cellType = QualityCellType2D::Unknown;
    std::size_t sourceId = 0;
    std::size_t owner = 0;
    std::optional<std::size_t> neighbour;
    double localH = 0.0;
};

struct QualityMetricSample2D {
    double value = 0.0;
    QualityEntity2D entity;
};

struct BoundaryLayerQualitySamples2D {
    std::vector<QualityMetricSample2D> wallNormalOrthogonalityDeg;
    std::vector<QualityMetricSample2D> growthRatio;
    std::vector<QualityMetricSample2D> tangentialNormalSpacingRatio;
    std::vector<QualityMetricSample2D> adjacentColumnThicknessVariation;
    std::vector<QualityMetricSample2D> scaledJacobian;
    std::vector<QualityMetricSample2D> firstLayerContinuity;
};

struct QualityMetricSummary2D {
    std::size_t count = 0;
    double p50 = 0.0;
    double p95 = 0.0;
    double p99 = 0.0;
    double worst = 0.0;
    bool lowerIsWorse = false;
    std::optional<QualityEntity2D> worstEntity;
};

struct QualityContractIssue2D {
    QualityContractLevel2D level = QualityContractLevel2D::Preferred;
    std::string metric;
    double measured = 0.0;
    double limit = 0.0;
    QualityEntity2D entity;
};

struct QualityContractTypeResult2D {
    std::size_t cellCount = 0;
    std::size_t preferredIssueCount = 0;
    std::size_t hardIssueCount = 0;
    bool rated = true;

    [[nodiscard]] QualityContractStatus2D status() const noexcept {
        if (hardIssueCount != 0U) return QualityContractStatus2D::Fail;
        if (preferredIssueCount != 0U) return QualityContractStatus2D::Warn;
        return QualityContractStatus2D::Pass;
    }
};

struct QualityContractReport2D {
    QualityContract2D contract;
    SolverQualityReport2D legacyHardSafety;
    std::map<QualityCellType2D,QualityContractTypeResult2D> byCellType;
    std::map<std::string,QualityMetricSummary2D> ordinaryMetrics;
    std::map<std::string,QualityMetricSummary2D> boundaryLayerMetrics;
    std::vector<QualityContractIssue2D> issues;
    std::vector<std::string> inputIssues;

    [[nodiscard]] QualityContractStatus2D status() const noexcept;
    [[nodiscard]] bool validInput() const noexcept { return inputIssues.empty(); }
};

[[nodiscard]] QualityContractReport2D evaluateQualityContract2D(
    const TopologyMesh2D& topology,
    const std::vector<QualityCellMetadata2D>& cellMetadata,
    const BoundaryLayerQualitySamples2D& boundaryLayerSamples = {},
    const QualityContract2D& contract = {},
    const SolverQualityReport2D* legacyHardSafety = nullptr,
    const TolerancePolicy& tol = {});

[[nodiscard]] const char* qualityCellTypeName(QualityCellType2D type) noexcept;
[[nodiscard]] const char* qualityContractStatusName(QualityContractStatus2D status) noexcept;
[[nodiscard]] std::string qualityContractReportToJson(
    const QualityContractReport2D& report, int indentSpaces = 2);

} // namespace cartmesh2d
