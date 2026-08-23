#include "cartmesh2d/stabilization/SmallCell2D.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <tuple>
#include <utility>

namespace cartmesh2d {
namespace {

struct SourceKey2D {
    std::uint64_t key = 0;
    std::size_t id = 0;

    [[nodiscard]] bool operator<(const SourceKey2D& rhs) const noexcept {
        return std::tie(key, id) < std::tie(rhs.key, rhs.id);
    }
};

[[nodiscard]] double edgeLength(const Edge2D& edge,
                                const TopologyMesh2D& topology) noexcept {
    if (edge.v0 >= topology.vertices.size() || edge.v1 >= topology.vertices.size()) {
        return 0.0;
    }
    const auto& a = topology.vertices[edge.v0].point;
    const auto& b = topology.vertices[edge.v1].point;
    return std::hypot(b.x - a.x, b.y - a.y);
}

[[nodiscard]] bool isSmallAlpha(double alpha, double threshold,
                                const TolerancePolicy& tol) noexcept {
    const double eps = tol.scale(std::max({1.0, std::abs(alpha), std::abs(threshold)}));
    return alpha < threshold - eps;
}

[[nodiscard]] std::vector<AlphaHistogramBin2D> makeHistogram(
    const SmallCellPolicy2D& policy) {
    std::vector<AlphaHistogramBin2D> bins;
    double lower = 0.0;
    for (const double upper : policy.histogramUpperBounds) {
        bins.push_back({lower, upper, 0});
        lower = upper;
    }
    return bins;
}

void addHistogramSample(std::vector<AlphaHistogramBin2D>& bins, double alpha,
                        const TolerancePolicy& tol) {
    for (auto& bin : bins) {
        const double eps = tol.scale(std::max({1.0, std::abs(alpha), std::abs(bin.upperInclusive)}));
        if (alpha <= bin.upperInclusive + eps) {
            ++bin.count;
            return;
        }
    }
}

struct NeighbourScore2D {
    std::size_t target = 0;
    double sharedLength = 0.0;
    double targetAlpha = 0.0;
    double targetArea = 0.0;
    bool targetSmall = false;
};

[[nodiscard]] bool betterCandidate(const NeighbourScore2D& lhs,
                                   const NeighbourScore2D& rhs,
                                   const TolerancePolicy& tol) noexcept {
    if (lhs.targetSmall != rhs.targetSmall) return !lhs.targetSmall;

    const auto clearlyGreater = [&](double a, double b) {
        return a > b + tol.scale(std::max({1.0, std::abs(a), std::abs(b)}));
    };
    if (clearlyGreater(lhs.sharedLength, rhs.sharedLength)) return true;
    if (clearlyGreater(rhs.sharedLength, lhs.sharedLength)) return false;
    if (clearlyGreater(lhs.targetAlpha, rhs.targetAlpha)) return true;
    if (clearlyGreater(rhs.targetAlpha, lhs.targetAlpha)) return false;
    if (clearlyGreater(lhs.targetArea, rhs.targetArea)) return true;
    if (clearlyGreater(rhs.targetArea, lhs.targetArea)) return false;
    return lhs.target < rhs.target;
}

} // namespace

SmallCellReport2D analyzeSmallCells(const std::vector<CutCell2D>& cutCells,
                                    const TopologyMesh2D& topology,
                                    const SmallCellPolicy2D& policy,
                                    const TolerancePolicy& tol) {
    SmallCellReport2D report;
    report.policy = policy;
    report.cutCellAlphaHistogram = makeHistogram(policy);

    if (!(policy.areaFractionThreshold > 0.0 && policy.areaFractionThreshold < 1.0)) {
        report.issues.push_back({SmallCellIssueCode2D::InvalidThreshold, 0,
                                 "area-fraction threshold must lie strictly in (0,1)"});
        return report;
    }
    if (policy.histogramUpperBounds.empty()) {
        report.issues.push_back({SmallCellIssueCode2D::InvalidThreshold, 0,
                                 "alpha histogram must contain at least one upper bound"});
        return report;
    }
    double previous = 0.0;
    for (const double upper : policy.histogramUpperBounds) {
        if (!(upper > previous && upper <= 1.0)) {
            report.issues.push_back({SmallCellIssueCode2D::InvalidThreshold, 0,
                                     "alpha histogram bounds must increase and end at or below 1"});
            return report;
        }
        previous = upper;
    }
    if (previous < 1.0 - tol.scale()) {
        report.issues.push_back({SmallCellIssueCode2D::InvalidThreshold, 0,
                                 "alpha histogram must cover alpha=1"});
        return report;
    }

    std::map<SourceKey2D, const CutCell2D*> sourceCells;
    for (const auto& cut : cutCells) {
        if (cut.kind == CutCellKind::Empty || cut.kind == CutCellKind::Unsupported) continue;
        sourceCells[{cut.sourceKey, cut.sourceId}] = &cut;
    }

    report.records.reserve(topology.cells.size());
    for (const auto& topoCell : topology.cells) {
        const auto it = sourceCells.find({topoCell.sourceKey, topoCell.sourceId});
        if (it == sourceCells.end()) {
            report.issues.push_back({SmallCellIssueCode2D::MissingSourceCell, topoCell.id,
                                     "topology cell has no matching non-empty CutCell source"});
            continue;
        }
        const CutCell2D& cut = *it->second;
        const double alpha = cut.areaFraction;
        const double eps = tol.scale(std::max(1.0, std::abs(alpha)));
        if (!(alpha > eps && alpha <= 1.0 + eps)) {
            report.issues.push_back({SmallCellIssueCode2D::InvalidAreaFraction, topoCell.id,
                                     "non-empty fluid cell has alpha outside (0,1]"});
            continue;
        }

        SmallCellRecord2D record;
        record.topologyCellId = topoCell.id;
        record.sourceId = topoCell.sourceId;
        record.sourceKey = topoCell.sourceKey;
        record.area = cut.area;
        record.areaFraction = std::clamp(alpha, 0.0, 1.0);
        report.records.push_back(record);
        ++report.fluidCellCount;
        if (cut.kind == CutCellKind::Full) {
            ++report.fullCellCount;
        } else if (cut.kind == CutCellKind::Cut) {
            ++report.cutCellCount;
            addHistogramSample(report.cutCellAlphaHistogram, record.areaFraction, tol);
        }
    }

    std::vector<bool> small(topology.cells.size(), false);
    std::vector<double> alphaByTopology(topology.cells.size(), 0.0);
    std::vector<double> areaByTopology(topology.cells.size(), 0.0);
    for (auto& record : report.records) {
        if (record.topologyCellId >= topology.cells.size()) continue;
        alphaByTopology[record.topologyCellId] = record.areaFraction;
        areaByTopology[record.topologyCellId] = record.area;
        small[record.topologyCellId] = isSmallAlpha(record.areaFraction,
                                                    policy.areaFractionThreshold, tol);
        if (small[record.topologyCellId]) ++report.smallCellCount;
    }

    // Aggregate internal-edge lengths once into a deterministic, cell-local
    // adjacency map.  The previous implementation stored the same pair data
    // globally and rescanned every edge pair for every small cell (O(S * E)).
    // Maps preserve the former ascending neighbour-id tie-break order while
    // reducing candidate discovery to O(E log d + sum_small d).
    std::vector<std::map<std::size_t, double>> neighbourLengths(topology.cells.size());
    for (const auto& edge : topology.edges) {
        if (!edge.neighbour) continue;
        const std::size_t a = edge.owner;
        const std::size_t b = *edge.neighbour;
        if (a >= topology.cells.size() || b >= topology.cells.size() || a == b) continue;
        const double length = edgeLength(edge, topology);
        neighbourLengths[a][b] += length;
        neighbourLengths[b][a] += length;
    }

    for (auto& record : report.records) {
        const std::size_t cellId = record.topologyCellId;
        if (cellId >= small.size() || !small[cellId]) {
            record.status = SmallCellStatus2D::Stable;
            continue;
        }

        std::optional<NeighbourScore2D> best;
        for (const auto& [target, length] : neighbourLengths[cellId]) {
            NeighbourScore2D candidate;
            candidate.target = target;
            candidate.sharedLength = length;
            candidate.targetAlpha = alphaByTopology[target];
            candidate.targetArea = areaByTopology[target];
            candidate.targetSmall = small[target];
            if (!(candidate.sharedLength > tol.scale())) continue;

            if (!best || betterCandidate(candidate, *best, tol)) best = candidate;
        }

        if (!best) {
            record.status = SmallCellStatus2D::Unresolved;
            ++report.unresolvedCount;
            report.issues.push_back({SmallCellIssueCode2D::NoNeighbourCandidate, cellId,
                                     "small cell has no internal-edge neighbour candidate"});
            continue;
        }

        record.status = SmallCellStatus2D::CandidateFound;
        record.targetTopologyCellId = best->target;
        record.sharedEdgeLength = best->sharedLength;
        record.targetIsSmall = best->targetSmall;
    }

    return report;
}

} // namespace cartmesh2d
