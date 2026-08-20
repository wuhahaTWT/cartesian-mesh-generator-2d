#include "cartmesh2d/io/MeshIO2D.hpp"
#include "cartmesh2d/quality/Quality2D.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace cartmesh2d;

namespace {

bool readBoundaryFile(const std::filesystem::path& path,
                      std::vector<Point2D>& points,
                      std::string& error) {
    std::ifstream in(path);
    if (!in) {
        error = "cannot open boundary file: " + path.string();
        return false;
    }
    std::string line;
    std::size_t lineNumber = 0;
    while (std::getline(in, line)) {
        ++lineNumber;
        const auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos || line[first] == '#') continue;
        std::istringstream row(line.substr(first));
        Point2D point;
        if (!(row >> point.x >> point.y)) {
            error = "invalid x y pair on boundary line " + std::to_string(lineNumber);
            return false;
        }
        std::string extra;
        if (row >> extra) {
            if (extra.empty() || extra.front() != '#') {
                error = "unexpected trailing token on boundary line " +
                        std::to_string(lineNumber);
                return false;
            }
        }
        points.push_back(point);
    }
    if (points.size() < 3) {
        error = "boundary file must contain at least three vertices";
        return false;
    }
    return true;
}

bool sameReadbackTopology(const TopologyMesh2D& a, const TopologyMesh2D& b) {
    if (a.vertices.size() != b.vertices.size() ||
        a.edges.size() != b.edges.size() ||
        a.cells.size() != b.cells.size()) return false;
    for (std::size_t i = 0; i < a.vertices.size(); ++i) {
        if (a.vertices[i].id != b.vertices[i].id ||
            a.vertices[i].point.x != b.vertices[i].point.x ||
            a.vertices[i].point.y != b.vertices[i].point.y) return false;
    }
    for (std::size_t i = 0; i < a.edges.size(); ++i) {
        const auto& x = a.edges[i];
        const auto& y = b.edges[i];
        if (x.id != y.id || x.v0 != y.v0 || x.v1 != y.v1 ||
            x.owner != y.owner || x.neighbour != y.neighbour || x.patch != y.patch) return false;
    }
    for (std::size_t i = 0; i < a.cells.size(); ++i) {
        const auto& x = a.cells[i];
        const auto& y = b.cells[i];
        if (x.id != y.id || x.sourceId != y.sourceId || x.sourceKey != y.sourceKey ||
            x.geometryArea != y.geometryArea || x.vertices != y.vertices || x.edges != y.edges) return false;
    }
    return true;
}

void usage() {
    std::cerr << "usage: cartmesh2d_cli <boundary.xy> <output-prefix> "
                 "[max-level=5] [padding-fraction=0.25] [small-alpha=0.10]\n";
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 3 || argc > 6) {
        usage();
        return EXIT_FAILURE;
    }

    const std::filesystem::path boundaryPath = argv[1];
    const std::filesystem::path outputPrefix = argv[2];
    std::size_t maxLevel = 5;
    double paddingFraction = 0.25;
    double smallAlpha = 0.10;
    try {
        if (argc >= 4) maxLevel = static_cast<std::size_t>(std::stoul(argv[3]));
        if (argc >= 5) paddingFraction = std::stod(argv[4]);
        if (argc >= 6) smallAlpha = std::stod(argv[5]);
    } catch (const std::exception&) {
        std::cerr << "invalid numeric CLI argument\n";
        return EXIT_FAILURE;
    }
    if (maxLevel == 0 || maxLevel > 28 || !(paddingFraction > 0.0) ||
        !(smallAlpha > 0.0 && smallAlpha < 1.0)) {
        std::cerr << "invalid max-level, padding-fraction or small-alpha\n";
        return EXIT_FAILURE;
    }

    std::vector<Point2D> points;
    std::string error;
    if (!readBoundaryFile(boundaryPath, points, error)) {
        std::cerr << error << '\n';
        return EXIT_FAILURE;
    }
    BoundaryLoop boundary(std::move(points));
    const auto diagnostics = boundary.diagnose();
    if (!diagnostics.valid()) {
        std::cerr << "boundary diagnostics failed with " << diagnostics.issues.size()
                  << " issue(s)\n";
        return EXIT_FAILURE;
    }
    if (!boundary.normalizeCounterClockwise()) {
        std::cerr << "failed to normalize boundary orientation\n";
        return EXIT_FAILURE;
    }

    const AABB2D bounds = boundary.polygon().bounds();
    const double width = bounds.max.x - bounds.min.x;
    const double height = bounds.max.y - bounds.min.y;
    const double span = std::max(width, height);
    if (!(span > 0.0)) {
        std::cerr << "boundary has zero extent\n";
        return EXIT_FAILURE;
    }
    const double padding = paddingFraction * span;
    const Domain2D domain{{
        {bounds.min.x - padding, bounds.min.y - padding},
        {bounds.max.x + padding, bounds.max.y + padding}
    }};

    Quadtree2D tree(domain, maxLevel, boundary);
    QuadtreeRefinementPolicy2D refinement;
    refinement.boundaryLevel = maxLevel;
    tree.refine(boundary, refinement);
    const auto balance = tree.enforceTwoToOneBalance(boundary);
    if (balance.violationsAfter != 0 || !tree.deterministicOrderingValid()) {
        std::cerr << "Quadtree balance/determinism gate failed\n";
        return EXIT_FAILURE;
    }

    std::vector<CutCell2D> cutCells;
    cutCells.reserve(tree.leaves().size());
    std::size_t unsupported = 0;
    for (const auto& leaf : tree.leaves()) {
        auto cut = buildCutCell(leaf, boundary);
        if (!cut.valid() && cut.kind == CutCellKind::Unsupported) ++unsupported;
        cutCells.push_back(std::move(cut));
    }
    if (unsupported != 0) {
        std::cerr << "Cut-cell construction produced " << unsupported
                  << " unsupported leaf geometry case(s)\n";
        return EXIT_FAILURE;
    }

    const TopologyMesh2D sourceTopology = buildGlobalTopology(cutCells, domain, boundary);
    if (!sourceTopology.valid()) {
        std::cerr << "source global topology audit failed\n";
        return EXIT_FAILURE;
    }

    SmallCellPolicy2D smallPolicy;
    smallPolicy.areaFractionThreshold = smallAlpha;
    const SmallCellReport2D smallReport = analyzeSmallCells(cutCells, sourceTopology, smallPolicy);
    if (!smallReport.valid()) {
        std::cerr << "small-cell analysis failed; unresolved=" << smallReport.unresolvedCount << '\n';
        return EXIT_FAILURE;
    }

    const AgglomerationResult2D stabilized =
        agglomerateSmallCells(cutCells, sourceTopology, smallReport, domain, boundary);
    if (!stabilized.valid()) {
        std::cerr << "small-cell agglomeration failed with " << stabilized.issues.size()
                  << " issue(s)\n";
        return EXIT_FAILURE;
    }

    const MeshQualityReport2D quality =
        evaluateMeshQuality(stabilized.topology, cutCells, &smallReport);
    if (!quality.valid()) {
        std::cerr << "final quality/topology report is invalid\n";
        return EXIT_FAILURE;
    }

    const auto parent = outputPrefix.parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    const std::filesystem::path vtkPath = outputPrefix.string() + ".vtk";
    const std::filesystem::path cm2dPath = outputPrefix.string() + ".cm2d";
    const std::filesystem::path qualityPath = outputPrefix.string() + ".quality.json";

    error.clear();
    if (!writeLegacyVtk2D(stabilized.topology, vtkPath, &error)) {
        std::cerr << error << '\n';
        return EXIT_FAILURE;
    }
    error.clear();
    if (!writeCm2dTopology(stabilized.topology, cm2dPath, &error)) {
        std::cerr << error << '\n';
        return EXIT_FAILURE;
    }
    {
        std::ofstream out(qualityPath);
        if (!out) {
            std::cerr << "failed to open quality JSON output\n";
            return EXIT_FAILURE;
        }
        out << qualityReportToJson(quality);
        if (!out.good()) {
            std::cerr << "failed while writing quality JSON output\n";
            return EXIT_FAILURE;
        }
    }

    const MeshReadback2D readback = readCm2dTopology(cm2dPath);
    if (!readback.valid() || !sameReadbackTopology(stabilized.topology, readback.topology)) {
        std::cerr << "independent CM2D read-back verification failed: " << readback.error << '\n';
        return EXIT_FAILURE;
    }

    std::cout << "cartmesh2d end-to-end PASS\n"
              << "leaf_count=" << tree.leaves().size() << '\n'
              << "source_cells=" << sourceTopology.cells.size() << '\n'
              << "small_cells=" << smallReport.smallCellCount << '\n'
              << "stabilized_cells=" << stabilized.topology.cells.size() << '\n'
              << "vertices=" << stabilized.topology.vertices.size() << '\n'
              << "edges=" << stabilized.topology.edges.size() << '\n'
              << "min_cell_area=" << quality.minCellArea << '\n'
              << "min_edge_length=" << quality.minEdgeLength << '\n'
              << "max_edge_aspect_ratio=" << quality.maxEdgeAspectRatio << '\n'
              << "max_centroid_skewness=" << quality.maxCentroidSkewness << '\n'
              << "vtk=" << vtkPath.string() << '\n'
              << "cm2d=" << cm2dPath.string() << '\n'
              << "quality_json=" << qualityPath.string() << '\n';
    return EXIT_SUCCESS;
}
