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

} // namespace cartmesh2d
