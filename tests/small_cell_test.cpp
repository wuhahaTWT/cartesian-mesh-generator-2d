#include "cartmesh2d/stabilization/SmallCell2D.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

using namespace cartmesh2d;

namespace {
int failures = 0;
void check(bool condition, const std::string& message) {
    if (!condition) {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}

CutCell2D polygonCell(std::size_t id, std::uint64_t key,
                      const std::vector<Point2D>& polygon,
                      const AABB2D& background, CutCellKind kind) {
    CutCell2D cell;
    cell.sourceId = id;
    cell.sourceKey = key;
    cell.backgroundBounds = background;
    cell.kind = kind;
    cell.fluidPolygon = {polygon};
    cell.area = cell.fluidPolygon.area();
    const double bgArea = (background.max.x - background.min.x) *
                          (background.max.y - background.min.y);
    cell.areaFraction = cell.area / bgArea;
    cell.centroid = cell.fluidPolygon.centroid();
    return cell;
}

const SmallCellRecord2D* findRecord(const SmallCellReport2D& report,
                                    std::size_t topologyId) {
    for (const auto& record : report.records) {
        if (record.topologyCellId == topologyId) return &record;
    }
    return nullptr;
}
} // namespace

int main() {
    const Domain2D domain{{{0.0, 0.0}, {2.0, 1.0}}};
    BoundaryLoop boundary({{0.9, 0.0}, {2.0, 0.0}, {2.0, 1.0}, {0.9, 1.0}});
    std::vector<CutCell2D> cells;
    cells.push_back(polygonCell(10, 10,
        {{0.9,0.0},{1.0,0.0},{1.0,1.0},{0.9,1.0}},
        {{0.0,0.0},{1.0,1.0}}, CutCellKind::Cut));
    cells.push_back(polygonCell(20, 20,
        {{1.0,0.0},{2.0,0.0},{2.0,1.0},{1.0,1.0}},
        {{1.0,0.0},{2.0,1.0}}, CutCellKind::Full));

    const auto topology = buildGlobalTopology(cells, domain, boundary);
    check(topology.valid(), "manufactured small-cell topology is valid");

    SmallCellPolicy2D policy;
    policy.areaFractionThreshold = 0.20;
    const auto report = analyzeSmallCells(cells, topology, policy);
    check(report.valid(), "small-cell report is valid");
    check(report.fluidCellCount == 2, "two fluid cells counted");
    check(report.cutCellCount == 1 && report.fullCellCount == 1,
          "cut/full counts are correct");
    check(report.smallCellCount == 1, "one small cell detected");
    check(report.unresolvedCount == 0, "small cell has a candidate");

    const auto* small = findRecord(report, 0);
    check(small != nullptr, "small topology cell has a report record");
    if (small) {
        check(std::abs(small->areaFraction - 0.1) <= 1.0e-12,
              "small-cell alpha is 0.1");
        check(small->status == SmallCellStatus2D::CandidateFound,
              "small cell receives a candidate");
        check(small->targetTopologyCellId && *small->targetTopologyCellId == 1,
              "stable right neighbour is selected");
        check(std::abs(small->sharedEdgeLength - 1.0) <= 1.0e-12,
              "candidate score uses full shared edge length");
        check(!small->targetIsSmall, "selected target is not small");
    }

    SmallCellPolicy2D exactPolicy;
    exactPolicy.areaFractionThreshold = 0.10;
    const auto exact = analyzeSmallCells(cells, topology, exactPolicy);
    check(exact.smallCellCount == 0, "alpha equal to threshold is not marked small");

    std::size_t histogramTotal = 0;
    std::size_t tenthBin = 0;
    for (const auto& bin : report.cutCellAlphaHistogram) {
        histogramTotal += bin.count;
        if (std::abs(bin.upperInclusive - 0.10) <= 1.0e-12) tenthBin = bin.count;
    }
    check(histogramTotal == 1 && tenthBin == 1,
          "cut-cell alpha histogram counts the 0.1 sample deterministically");

    BoundaryLoop wideBoundary({{0.0,0.0},{2.0,0.0},{2.0,1.0},{0.0,1.0}});
    std::vector<CutCell2D> tieCells;
    tieCells.push_back(polygonCell(1, 1,
        {{0.0,0.0},{0.9,0.0},{0.9,1.0},{0.0,1.0}},
        {{0.0,0.0},{0.9,1.0}}, CutCellKind::Full));
    tieCells.push_back(polygonCell(2, 2,
        {{0.9,0.0},{1.1,0.0},{1.1,1.0},{0.9,1.0}},
        {{0.0,0.0},{2.0,1.0}}, CutCellKind::Cut));
    tieCells.push_back(polygonCell(3, 3,
        {{1.1,0.0},{2.0,0.0},{2.0,1.0},{1.1,1.0}},
        {{1.1,0.0},{2.0,1.0}}, CutCellKind::Full));
    const auto tieTopology = buildGlobalTopology(tieCells, domain, wideBoundary);
    check(tieTopology.valid(), "tie-break topology is valid");
    const auto tieReport = analyzeSmallCells(tieCells, tieTopology, policy);
    check(tieReport.valid(), "tie-break report is valid");
    const auto* tieSmall = findRecord(tieReport, 1);
    check(tieSmall && tieSmall->targetTopologyCellId && *tieSmall->targetTopologyCellId == 0,
          "equal candidates resolve to the lowest topology id deterministically");

    const Domain2D isolatedDomain{{{0.0,0.0},{1.0,1.0}}};
    BoundaryLoop isolatedBoundary({{0.0,0.0},{0.05,0.0},{0.05,1.0},{0.0,1.0}});
    std::vector<CutCell2D> isolatedCells;
    isolatedCells.push_back(polygonCell(5, 5,
        {{0.0,0.0},{0.05,0.0},{0.05,1.0},{0.0,1.0}},
        {{0.0,0.0},{1.0,1.0}}, CutCellKind::Cut));
    const auto isolatedTopology = buildGlobalTopology(isolatedCells, isolatedDomain,
                                                       isolatedBoundary);
    check(isolatedTopology.valid(), "isolated topology itself is valid");
    const auto isolated = analyzeSmallCells(isolatedCells, isolatedTopology, policy);
    check(!isolated.valid(), "unresolved small cell makes the report non-valid");
    check(isolated.smallCellCount == 1 && isolated.unresolvedCount == 1,
          "isolated small cell is counted as unresolved");

    // This fixture intentionally tests the historical interior-fluid
    // stabilization geometry. Keep its physical side explicit so changing the
    // product default can never silently change the regression's meaning.
    std::vector<Point2D> shiftedCircle;
    constexpr std::size_t circleSegments = 64;
    for (std::size_t i = 0; i < circleSegments; ++i) {
        const double angle = 2.0 * std::numbers::pi * static_cast<double>(i) /
                             static_cast<double>(circleSegments);
        shiftedCircle.push_back({0.07 + std::cos(angle), 0.03 + std::sin(angle)});
    }
    BoundaryLoop circleBoundary(shiftedCircle);
    const Domain2D circleDomain{{{-2.0,-2.0},{2.0,2.0}}};
    Quadtree2D circleTree(circleDomain, 4, circleBoundary);
    QuadtreeRefinementPolicy2D circleRefinement;
    circleRefinement.boundaryLevel = 4;
    circleTree.refine(circleBoundary, circleRefinement);
    const auto circleBalance = circleTree.enforceTwoToOneBalance(circleBoundary);
    check(circleBalance.violationsAfter == 0, "shifted-circle Quadtree is 2:1 balanced");
    std::vector<CutCell2D> circleCuts;
    for (const auto& leaf : circleTree.leaves()) {
        circleCuts.push_back(buildCutCell(leaf, circleBoundary, FluidRegion2D::Interior));
    }
    const auto circleTopology = buildGlobalTopology(circleCuts, circleDomain, circleBoundary);
    check(circleTopology.valid(), "shifted-circle topology is valid");
    SmallCellPolicy2D circlePolicy;
    circlePolicy.areaFractionThreshold = 0.10;
    const auto circleReport = analyzeSmallCells(circleCuts, circleTopology, circlePolicy);
    check(circleReport.valid(), "shifted-circle small-cell report is valid");
    check(circleReport.smallCellCount == 8, "shifted circle exposes eight alpha<0.1 cells");
    check(circleReport.unresolvedCount == 0, "all shifted-circle small cells have candidates");
    double minSmallAlpha = 1.0;
    std::size_t candidatesTargetingSmall = 0;
    for (const auto& record : circleReport.records) {
        if (record.status == SmallCellStatus2D::CandidateFound) {
            minSmallAlpha = std::min(minSmallAlpha, record.areaFraction);
            if (record.targetIsSmall) ++candidatesTargetingSmall;
        }
    }
    check(minSmallAlpha < 0.002, "real adaptive fixture contains a genuinely tiny cell");
    check(candidatesTargetingSmall == 0, "stable neighbours are preferred for all shifted-circle tiny cells");

    SmallCellPolicy2D badPolicy;
    badPolicy.areaFractionThreshold = 1.0;
    const auto invalid = analyzeSmallCells(cells, topology, badPolicy);
    check(!invalid.valid(), "invalid alpha threshold is rejected");

    if (failures != 0) {
        std::cerr << failures << " test(s) failed\n";
        return EXIT_FAILURE;
    }
    std::cout << "cartmesh2d 2D-5A small-cell analysis tests passed\n";
    return EXIT_SUCCESS;
}
