#pragma once

#include "cartmesh2d/topology/Topology2D.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace cartmesh2d {

struct AlphaHistogramBin2D {
    double lowerExclusive = 0.0;
    double upperInclusive = 0.0;
    std::size_t count = 0;
};

enum class SmallCellStatus2D {
    Stable,
    CandidateFound,
    Unresolved
};

enum class SmallCellIssueCode2D {
    InvalidThreshold,
    MissingSourceCell,
    InvalidAreaFraction,
    NoNeighbourCandidate
};

struct SmallCellIssue2D {
    SmallCellIssueCode2D code;
    std::size_t topologyCellId = 0;
    std::string message;
};

struct SmallCellRecord2D {
    std::size_t topologyCellId = 0;
    std::size_t sourceId = 0;
    std::uint64_t sourceKey = 0;
    double area = 0.0;
    double areaFraction = 1.0;
    SmallCellStatus2D status = SmallCellStatus2D::Stable;
    std::optional<std::size_t> targetTopologyCellId;
    double sharedEdgeLength = 0.0;
    bool targetIsSmall = false;
};

struct SmallCellPolicy2D {
    double areaFractionThreshold = 0.10;
    std::vector<double> histogramUpperBounds{0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 1.0};
};

struct SmallCellReport2D {
    SmallCellPolicy2D policy;
    std::size_t fluidCellCount = 0;
    std::size_t fullCellCount = 0;
    std::size_t cutCellCount = 0;
    std::size_t smallCellCount = 0;
    std::size_t unresolvedCount = 0;
    std::vector<AlphaHistogramBin2D> cutCellAlphaHistogram;
    std::vector<SmallCellRecord2D> records;
    std::vector<SmallCellIssue2D> issues;

    [[nodiscard]] bool valid() const noexcept {
        return issues.empty();
    }
};

[[nodiscard]] SmallCellReport2D analyzeSmallCells(
    const std::vector<CutCell2D>& cutCells,
    const TopologyMesh2D& topology,
    const SmallCellPolicy2D& policy = {},
    const TolerancePolicy& tol = {});

enum class AgglomerationIssueCode2D {
    InvalidInput,
    MissingCandidate,
    SmallTargetUnsupported,
    InvalidTopologyReference,
    DisconnectedBoundary,
    MultipleBoundaryLoops,
    DegenerateMergedPolygon,
    AreaMismatch,
    RebuiltTopologyInvalid
};

struct AgglomerationIssue2D {
    AgglomerationIssueCode2D code;
    std::size_t objectId = 0;
    std::string message;
};

struct AgglomeratedCell2D {
    std::size_t id = 0;
    std::vector<std::size_t> memberTopologyCellIds;
    std::vector<std::size_t> memberSourceIds;
    Polygon2D polygon;
    double area = 0.0;
    std::optional<Point2D> centroid;
};

struct AgglomerationResult2D {
    std::size_t inputCellCount = 0;
    std::size_t outputCellCount = 0;
    std::size_t mergedSmallCellCount = 0;
    double totalAreaBefore = 0.0;
    double totalAreaAfter = 0.0;
    double areaError = 0.0;
    std::vector<AgglomeratedCell2D> cells;
    TopologyMesh2D topology;
    std::vector<AgglomerationIssue2D> issues;

    [[nodiscard]] bool valid() const noexcept {
        return issues.empty() && topology.valid();
    }
};

[[nodiscard]] AgglomerationResult2D agglomerateSmallCells(
    const std::vector<CutCell2D>& cutCells,
    const TopologyMesh2D& topology,
    const SmallCellReport2D& analysis,
    const Domain2D& domain,
    const BoundaryLoop& boundary,
    const TolerancePolicy& tol = {});

[[nodiscard]] AgglomerationResult2D agglomerateSmallCells(
    const std::vector<CutCell2D>& cutCells,
    const TopologyMesh2D& topology,
    const SmallCellReport2D& analysis,
    const Domain2D& domain,
    const BoundaryRegion2D& boundary,
    const TolerancePolicy& tol = {});

} // namespace cartmesh2d
