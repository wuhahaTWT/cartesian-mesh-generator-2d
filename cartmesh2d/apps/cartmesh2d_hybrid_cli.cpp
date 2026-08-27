#include "cartmesh2d/hybrid/HybridMesh2D.hpp"
#include "cartmesh2d/io/MeshIO2D.hpp"
#include "cartmesh2d/io/OpenFoam2D.hpp"

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

bool parseSize(const std::string& text, std::size_t& value) {
    try {
        std::size_t consumed = 0U;
        const auto raw = std::stoull(text, &consumed);
        value = static_cast<std::size_t>(raw);
        return consumed == text.size() && static_cast<unsigned long long>(value) == raw;
    } catch (const std::exception&) {
        return false;
    }
}

bool parseDouble(const std::string& text, double& value) {
    try {
        std::size_t consumed = 0U;
        value = std::stod(text, &consumed);
        return consumed == text.size() && std::isfinite(value);
    } catch (const std::exception&) {
        return false;
    }
}

bool readLoops(const std::filesystem::path& path,
               std::vector<std::vector<Point2D>>& loops,
               std::string& error) {
    std::ifstream input(path);
    if (!input) {
        error = "cannot open boundary file: " + path.string();
        return false;
    }
    std::vector<Point2D> current;
    std::string line;
    std::size_t lineNumber = 0U;
    const auto finish = [&]() {
        if (!current.empty()) {
            loops.push_back(std::move(current));
            current.clear();
        }
    };
    while (std::getline(input, line)) {
        ++lineNumber;
        const auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) {
            finish();
            continue;
        }
        if (line[first] == '#') continue;
        std::istringstream row(line.substr(first));
        Point2D point;
        if (!(row >> point.x >> point.y) || !std::isfinite(point.x) ||
            !std::isfinite(point.y)) {
            error = "invalid finite x y pair on boundary line " +
                    std::to_string(lineNumber);
            return false;
        }
        std::string trailing;
        if (row >> trailing && (trailing.empty() || trailing.front() != '#')) {
            error = "unexpected token on boundary line " +
                    std::to_string(lineNumber);
            return false;
        }
        current.push_back(point);
    }
    finish();
    if (loops.empty()) {
        error = "boundary file contains no loops";
        return false;
    }
    return true;
}

bool writeText(const std::filesystem::path& path, const std::string& text,
               std::string& error) {
    std::ofstream output(path);
    if (!output) {
        error = "failed to open report: " + path.string();
        return false;
    }
    output << text << '\n';
    if (!output.good()) {
        error = "failed while writing report: " + path.string();
        return false;
    }
    return true;
}

void usage() {
    std::cerr
        << "usage: cartmesh2d_hybrid_cli <boundary.xy> <output-prefix> "
           "<max-level> <minimum-level> <boundary-level> "
           "<n-layers> <first-thickness> <growth-ratio> <domain-padding> "
           "[openfoam-case extrusion-thickness]\n";
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 10 && argc != 12) {
        usage();
        return EXIT_FAILURE;
    }
    const std::filesystem::path boundaryPath = argv[1];
    const std::filesystem::path outputPrefix = argv[2];
    std::size_t maxLevel = 0U;
    QuadtreeRefinementPolicy2D refinement;
    LayerParameters2D layerParameters;
    double padding = 0.0;
    if (!parseSize(argv[3], maxLevel) ||
        !parseSize(argv[4], refinement.minimumLevel) ||
        !parseSize(argv[5], refinement.boundaryLevel) ||
        !parseSize(argv[6], layerParameters.nLayers) ||
        !parseDouble(argv[7], layerParameters.thickness) ||
        !parseDouble(argv[8], layerParameters.growthRatio) ||
        !parseDouble(argv[9], padding) || padding <= 0.0) {
        std::cerr << "invalid finite H4-2 numeric parameter\n";
        return EXIT_FAILURE;
    }
    layerParameters.thicknessMode = LayerThicknessMode2D::FirstLayerThickness;

    std::vector<std::vector<Point2D>> loopPoints;
    std::string error;
    if (!readLoops(boundaryPath, loopPoints, error)) {
        std::cerr << error << '\n';
        return EXIT_FAILURE;
    }
    std::vector<BoundaryLoop> loops;
    loops.reserve(loopPoints.size());
    for (auto& points : loopPoints) loops.emplace_back(std::move(points));
    BoundaryRegion2D originalWalls(loops);
    if (!originalWalls.diagnose().valid()) {
        std::cerr << "invalid original wall region\n";
        return EXIT_FAILURE;
    }
    const auto depths = originalWalls.nestingDepths();
    if (std::any_of(depths.begin(), depths.end(),
                    [](std::size_t depth) { return depth != 0U; })) {
        std::cerr << "H4-2 fixed exterior strips reject nested wall loops\n";
        return EXIT_FAILURE;
    }

    std::vector<WallChain2D> chains;
    chains.reserve(loops.size());
    for (std::size_t loopId = 0; loopId < loops.size(); ++loopId) {
        auto chain = makeClosedWallChain2D(
            loops[loopId], loopId, "wall_" + std::to_string(loopId));
        if (!chain.success()) {
            std::cerr << "wall-chain failure: " << chain.message << '\n';
            return EXIT_FAILURE;
        }
        chains.push_back(std::move(*chain.chain));
    }
    const auto layers = buildBoundaryLayerStrips2D(chains, layerParameters);
    if (!layers.success()) {
        std::cerr << "H4-1 prerequisite failed: "
                  << boundaryLayerFailureReasonName(layers.failure.reason)
                  << " " << layers.failure.message << '\n';
        return EXIT_FAILURE;
    }

    const auto wallBounds = originalWalls.bounds();
    const Domain2D domain{{{wallBounds.min.x - padding, wallBounds.min.y - padding},
                           {wallBounds.max.x + padding, wallBounds.max.y + padding}}};
    const auto hybrid = buildConformalHybridMesh2D(
        layers, domain, originalWalls, maxLevel, refinement);

    const auto parent = outputPrefix.parent_path().empty()
        ? std::filesystem::path(".") : outputPrefix.parent_path();
    std::filesystem::create_directories(parent);
    const auto jsonPath = outputPrefix.string() + ".hybrid.json";
    if (!writeHybridReportJson2D(hybrid, jsonPath, &error)) {
        std::cerr << error << '\n';
        return EXIT_FAILURE;
    }
    if (!hybrid.success()) {
        std::cerr << "hybrid_status=failed failure_reason="
                  << hybridMeshFailureReasonName(hybrid.failure.reason)
                  << " message=" << hybrid.failure.message
                  << " report=" << jsonPath << '\n';
        return EXIT_FAILURE;
    }

    const auto vtkPath = outputPrefix.string() + ".hybrid.vtk";
    const auto cm2dPath = outputPrefix.string() + ".hybrid.cm2d";
    const auto qualityPath = outputPrefix.string() + ".hybrid.quality.json";
    const auto solverQualityPath = outputPrefix.string() +
                                   ".hybrid.solver-quality.json";
    if (!writeHybridLegacyVtk2D(hybrid, vtkPath, &error) ||
        !writeCm2dTopology(hybrid.topology, cm2dPath, &error) ||
        !writeText(qualityPath, qualityReportToJson(hybrid.meshQuality), error) ||
        !writeText(solverQualityPath,
                   solverQualityReportToJson(hybrid.solverQuality), error)) {
        std::cerr << error << '\n';
        return EXIT_FAILURE;
    }
    const auto readback = readCm2dTopology(cm2dPath);
    if (!readback.valid()) {
        std::cerr << "CM2D readback failed: " << readback.error << '\n';
        return EXIT_FAILURE;
    }

    std::string openFoamStatus = "not_requested";
    if (argc == 12) {
        double thickness = 0.0;
        if (!parseDouble(argv[11], thickness) || thickness <= 0.0) {
            std::cerr << "extrusion thickness must be finite and positive\n";
            return EXIT_FAILURE;
        }
        if (!hybrid.solverQuality.valid()) {
            std::cerr << "hybrid solver-quality gate failed; refusing OpenFOAM output\n";
            return EXIT_FAILURE;
        }
        const auto foam = writeExtrudedOpenFoam2D(
            hybrid.topology, domain, originalWalls, argv[10], thickness, &error);
        if (!foam.valid()) {
            std::cerr << "OpenFOAM output failed: " << error << '\n';
            return EXIT_FAILURE;
        }
        openFoamStatus = "written";
    }

    std::cout << "hybrid_status=success cells=" << hybrid.metrics.unifiedCellCount
              << " layer_cells=" << hybrid.metrics.boundaryLayerCellCount
              << " remainder_cut=" << hybrid.metrics.remainderCutCellCount
              << " remainder_cartesian="
              << hybrid.metrics.remainderCartesianCellCount
              << " interface_edges=" << hybrid.interfaceAudit.interfaceEdgeCount
              << " interface_vertices=" << hybrid.interfaceAudit.interfaceVertexCount
              << " area_error=" << hybrid.metrics.areaError
              << " solver_quality="
              << (hybrid.solverQuality.valid() ? "pass" : "fail")
              << " openfoam=" << openFoamStatus
              << " vtk=" << vtkPath << " report=" << jsonPath << '\n';
    return EXIT_SUCCESS;
}
