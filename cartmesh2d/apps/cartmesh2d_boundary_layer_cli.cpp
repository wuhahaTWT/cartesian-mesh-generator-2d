#include "cartmesh2d/boundary_layer/BoundaryLayer2D.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
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

void usage() {
    std::cerr
        << "usage: cartmesh2d_boundary_layer_cli <boundary.xy> <output-prefix> "
           "<n-layers> <first|total> <thickness> <growth-ratio> "
           "[exterior|interior]\n";
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 7 || argc > 8) {
        usage();
        return EXIT_FAILURE;
    }
    const std::filesystem::path inputPath = argv[1];
    const std::filesystem::path outputPrefix = argv[2];
    LayerParameters2D parameters;
    if (!parseSize(argv[3], parameters.nLayers)) {
        std::cerr << "invalid n-layers\n";
        return EXIT_FAILURE;
    }
    const std::string thicknessMode = argv[4];
    if (thicknessMode == "first") {
        parameters.thicknessMode = LayerThicknessMode2D::FirstLayerThickness;
    } else if (thicknessMode == "total") {
        parameters.thicknessMode = LayerThicknessMode2D::TotalThickness;
    } else {
        std::cerr << "thickness mode must be first or total\n";
        return EXIT_FAILURE;
    }
    if (!parseDouble(argv[5], parameters.thickness) ||
        !parseDouble(argv[6], parameters.growthRatio)) {
        std::cerr << "thickness and growth-ratio must be finite numbers\n";
        return EXIT_FAILURE;
    }
    WallFluidRegion2D fluidRegion = WallFluidRegion2D::Exterior;
    if (argc == 8) {
        const std::string region = argv[7];
        if (region == "exterior") fluidRegion = WallFluidRegion2D::Exterior;
        else if (region == "interior") fluidRegion = WallFluidRegion2D::Interior;
        else {
            std::cerr << "fluid region must be exterior or interior\n";
            return EXIT_FAILURE;
        }
    }

    std::vector<std::vector<Point2D>> loopPoints;
    std::string error;
    if (!readLoops(inputPath, loopPoints, error)) {
        std::cerr << error << '\n';
        return EXIT_FAILURE;
    }
    std::vector<BoundaryLoop> loops;
    loops.reserve(loopPoints.size());
    for (auto& points : loopPoints) loops.emplace_back(std::move(points));
    BoundaryRegion2D region(loops);
    const auto diagnostics = region.diagnose();
    if (!diagnostics.valid()) {
        std::cerr << "invalid boundary region: "
                  << diagnostics.issues.front().message << '\n';
        return EXIT_FAILURE;
    }
    const auto depths = region.nestingDepths();
    if (std::any_of(depths.begin(), depths.end(),
                    [](std::size_t depth) { return depth != 0U; })) {
        std::cerr << "H4-1 CLI rejects nested loops because per-loop wall-fluid "
                     "semantics are not yet part of this isolated strip stage\n";
        return EXIT_FAILURE;
    }

    std::vector<WallChain2D> chains;
    chains.reserve(loops.size());
    for (std::size_t loopId = 0; loopId < loops.size(); ++loopId) {
        auto built = makeClosedWallChain2D(
            loops[loopId], loopId, "wall_" + std::to_string(loopId), fluidRegion);
        if (!built.success()) {
            std::cerr << "invalid wall chain " << loopId << ": " << built.message << '\n';
            return EXIT_FAILURE;
        }
        chains.push_back(std::move(*built.chain));
    }

    const auto result = buildBoundaryLayerStrips2D(chains, parameters);
    std::filesystem::create_directories(outputPrefix.parent_path().empty()
        ? std::filesystem::path(".") : outputPrefix.parent_path());
    const auto reportPath = outputPrefix.string() + ".layer.json";
    if (!writeBoundaryLayerReportJson2D(result, reportPath, &error)) {
        std::cerr << error << '\n';
        return EXIT_FAILURE;
    }
    if (!result.success()) {
        std::cerr << "layer_status=failed failure_reason="
                  << boundaryLayerFailureReasonName(result.failure.reason)
                  << " chain_id=" << result.failure.chainId
                  << " requested_thickness=" << result.failure.requestedThickness;
        if (result.failure.safeThickness) {
            std::cerr << " safe_thickness=" << *result.failure.safeThickness;
        }
        std::cerr << " report=" << reportPath << '\n';
        return EXIT_FAILURE;
    }

    const auto vtkPath = outputPrefix.string() + ".layer.vtk";
    if (!writeBoundaryLayerLegacyVtk2D(result, vtkPath, &error)) {
        std::cerr << error << '\n';
        return EXIT_FAILURE;
    }
    std::size_t vertices = 0U;
    std::size_t cells = 0U;
    double minArea = std::numeric_limits<double>::infinity();
    double maxArea = 0.0;
    double minThickness = std::numeric_limits<double>::infinity();
    double maxThickness = 0.0;
    for (const auto& strip : result.strips) {
        vertices += strip.metrics.vertexCount;
        cells += strip.metrics.cellCount;
        minArea = std::min(minArea, strip.metrics.minCellArea);
        maxArea = std::max(maxArea, strip.metrics.maxCellArea);
        minThickness = std::min(minThickness, strip.metrics.minLayerThickness);
        maxThickness = std::max(maxThickness, strip.metrics.maxLayerThickness);
    }
    std::cout << "layer_status=success strips=" << result.strips.size()
              << " vertices=" << vertices << " cells=" << cells
              << " min_cell_area=" << minArea
              << " max_cell_area=" << maxArea
              << " min_layer_thickness=" << minThickness
              << " max_layer_thickness=" << maxThickness
              << " vtk=" << vtkPath << " report=" << reportPath << '\n';
    return EXIT_SUCCESS;
}
