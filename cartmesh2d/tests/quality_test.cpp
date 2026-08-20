#include "cartmesh2d/quality/Quality2D.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
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
                      const AABB2D& background,
                      CutCellKind kind) {
    CutCell2D cell;
    cell.sourceId = id;
    cell.sourceKey = key;
    cell.backgroundBounds = background;
    cell.kind = kind;
    cell.fluidPolygon = {polygon};
    cell.area = cell.fluidPolygon.area();
    const double backgroundArea = (background.max.x - background.min.x) *
                                  (background.max.y - background.min.y);
    cell.areaFraction = cell.area / backgroundArea;
    cell.centroid = cell.fluidPolygon.centroid();
    return cell;
}
} // namespace

int main() {
    const Domain2D domain{{{0.0, 0.0}, {2.0, 1.0}}};
    BoundaryLoop boundary({{0.5, 0.0}, {2.0, 0.0}, {2.0, 1.0}, {0.5, 1.0}});

    std::vector<CutCell2D> cells;
    cells.push_back(polygonCell(10, 3,
        {{0.5,0.0},{1.0,0.0},{1.0,1.0},{0.5,1.0}},
        {{0.0,0.0},{1.0,1.0}}, CutCellKind::Cut));
    cells.push_back(polygonCell(20, 4,
        {{1.0,0.0},{2.0,0.0},{2.0,1.0},{1.0,1.0}},
        {{1.0,0.0},{2.0,1.0}}, CutCellKind::Full));

    const auto topology = buildGlobalTopology(cells, domain, boundary);
    check(topology.valid(), "quality fixture topology is valid");

    SmallCellReport2D small;
    small.smallCellCount = 1;
    const auto report = evaluateMeshQuality(topology, cells, &small);
    check(report.valid(), "quality report is valid");
    check(report.vertexCount == topology.vertices.size() &&
          report.edgeCount == topology.edges.size() &&
          report.cellCount == 2,
          "topology counts are preserved");
    check(report.sourceCutCellCount == 1 && report.sourceFullCellCount == 1 &&
          report.sourceSmallCellCount == 1,
          "source cut/full/small counts are reported");
    check(std::abs(report.minCellArea - 0.5) <= 1.0e-12,
          "minimum solver-cell area is exact");
    check(std::abs(report.minEdgeLength - 0.5) <= 1.0e-12,
          "minimum global edge length is exact");
    check(std::abs(report.maxEdgeAspectRatio - 2.0) <= 1.0e-12,
          "edge aspect ratio is max-edge/min-edge");
    check(report.maxCentroidSkewness <= 1.0e-12,
          "rectangular cells have zero centroid skewness");
    check(std::abs(report.minCutCellAreaFraction - 0.5) <= 1.0e-12,
          "minimum Cut-cell alpha is exact");
    check(report.levelDistribution.size() == 2 &&
          report.levelDistribution.at(3) == 1 &&
          report.levelDistribution.at(4) == 1,
          "source-key level distribution is deterministic");

    const std::string json = qualityReportToJson(report);
    check(json.find("\"valid\": true") != std::string::npos,
          "JSON report contains valid=true");
    check(json.find("\"min_cell_area\": 0.5") != std::string::npos,
          "JSON report contains minimum area");
    check(json.find("\"3\": 1") != std::string::npos &&
          json.find("\"4\": 1") != std::string::npos,
          "JSON report contains deterministic level keys");

    BoundaryLoop skewBoundary({{0.0,0.0},{2.0,0.0},{1.4,1.0},{0.0,1.0}});
    std::vector<CutCell2D> skewCells;
    skewCells.push_back(polygonCell(30, 2,
        {{0.0,0.0},{2.0,0.0},{1.4,1.0},{0.0,1.0}},
        {{0.0,0.0},{2.0,1.0}}, CutCellKind::Cut));
    const auto skewTopology = buildGlobalTopology(skewCells,
        Domain2D{{{0.0,0.0},{2.0,1.0}}}, skewBoundary);
    check(skewTopology.valid(), "skew fixture topology is valid");
    const auto skewReport = evaluateMeshQuality(skewTopology, skewCells);
    check(skewReport.maxCentroidSkewness > 0.0,
          "asymmetric trapezoid produces positive centroid skewness");

    auto invalidTopology = topology;
    invalidTopology.audit.duplicateEdges = 1;
    const auto invalid = evaluateMeshQuality(invalidTopology, cells);
    check(!invalid.valid(), "invalid topology cannot produce a passing quality report");
    check(!invalid.issues.empty(), "invalid topology is reported explicitly");

    if (failures != 0) {
        std::cerr << failures << " test(s) failed\n";
        return EXIT_FAILURE;
    }
    std::cout << "cartmesh2d 2D-6A quality tests passed\n";
    return EXIT_SUCCESS;
}
