#include "cartmesh2d/quality/QualityContract2D.hpp"

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

CutCell2D polygonCell(std::size_t id, const std::vector<Point2D>& polygon,
                      const AABB2D& background) {
    CutCell2D cell;
    cell.sourceId = id;
    cell.sourceKey = id;
    cell.backgroundBounds = background;
    cell.kind = CutCellKind::Full;
    cell.fluidPolygon = {polygon};
    cell.area = cell.fluidPolygon.area();
    cell.areaFraction = 1.0;
    cell.centroid = cell.fluidPolygon.centroid();
    return cell;
}

QualityContractReport2D makeReport(double scale) {
    const double height = 0.005 * scale;
    const Domain2D domain{{{0.0, 0.0}, {2.0 * scale, height}}};
    BoundaryLoop boundary({{0.5 * scale, 0.0}, {2.0 * scale, 0.0},
                           {2.0 * scale, height}, {0.5 * scale, height}});
    std::vector<CutCell2D> cells;
    cells.push_back(polygonCell(10,
        {{0.5 * scale, 0.0}, {1.0 * scale, 0.0},
         {1.0 * scale, height}, {0.5 * scale, height}},
        {{0.0, 0.0}, {1.0 * scale, scale}}));
    cells.push_back(polygonCell(20,
        {{1.0 * scale, 0.0}, {2.0 * scale, 0.0},
         {2.0 * scale, height}, {1.0 * scale, height}},
        {{1.0 * scale, 0.0}, {2.0 * scale, scale}}));
    const auto topology = buildGlobalTopology(cells, domain, boundary);
    check(topology.valid(), "dimensionless fixture topology is valid");
    std::vector<QualityCellMetadata2D> metadata(topology.cells.size());
    for (std::size_t i = 0; i < metadata.size(); ++i) {
        metadata[i].type = i == 0U ? QualityCellType2D::Cartesian
                                  : QualityCellType2D::Transition;
        metadata[i].localBackgroundH = scale;
        metadata[i].sourceId = 100U + i;
    }
    return evaluateQualityContract2D(topology, metadata);
}

bool hasHardShortFace(const QualityContractReport2D& report) {
    for (const auto& issue : report.issues) {
        if (issue.level == QualityContractLevel2D::Hard &&
            issue.metric == "face_length_over_local_background_h" &&
            issue.measured < 0.01) return true;
    }
    return false;
}
} // namespace

int main() {
    const auto base = makeReport(1.0);
    const auto scaled = makeReport(1.0e6);

    check(base.status() == QualityContractStatus2D::Fail,
          "dimensionless micro-face violates the hard contract");
    check(scaled.status() == base.status(),
          "uniform scaling preserves the overall contract status");
    check(hasHardShortFace(base) && hasHardShortFace(scaled),
          "micro-face is found by face/local_h at both scales");
    check(base.byCellType.at(QualityCellType2D::Cartesian).status() ==
              QualityContractStatus2D::Fail &&
          base.byCellType.at(QualityCellType2D::Transition).status() ==
              QualityContractStatus2D::Fail,
          "cell types are evaluated independently");

    for (const auto& [name, summary] : base.ordinaryMetrics) {
        const auto found = scaled.ordinaryMetrics.find(name);
        check(found != scaled.ordinaryMetrics.end(),
              "scaled report preserves metric " + name);
        if (found == scaled.ordinaryMetrics.end()) continue;
        const double tolerance = name == "non_orthogonality_deg" ? 1.0e-5 :
            1.0e-10 * std::max(1.0, std::abs(summary.worst));
        check(std::abs(found->second.worst - summary.worst) <= tolerance,
              "scaled report preserves dimensionless worst " + name);
    }
    check(std::abs(base.legacyHardSafety.minFaceLength - 0.005) <= 1.0e-12,
          "legacy absolute minimum face remains visible");
    check(std::abs(scaled.legacyHardSafety.minFaceLength - 5000.0) <= 1.0e-8,
          "legacy absolute safety data scales and is not hidden");

    const std::string json = qualityContractReportToJson(base);
    check(json.find("\"status\": \"FAIL\"") != std::string::npos,
          "contract JSON records FAIL");
    check(json.find("\"legacy_hard_safety\"") != std::string::npos,
          "contract JSON retains legacy hard safety");
    check(json.find("face_length_over_local_background_h") != std::string::npos,
          "contract JSON names the dimensionless short-face metric");

    if (failures != 0) {
        std::cerr << failures << " test(s) failed\n";
        return EXIT_FAILURE;
    }
    std::cout << "cartmesh2d Q1 quality contract tests passed\n";
    return EXIT_SUCCESS;
}
