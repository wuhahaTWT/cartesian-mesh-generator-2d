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
           "[openfoam-case extrusion-thickness] [--legacy-construction] "
           "[--verify-source-lineage]\n";
}

} // namespace

int main(int argc, char** argv) {
    bool legacyConstruction=false;
    bool verifySourceLineage=false;
    while (argc>1) {
        const std::string option=argv[argc-1];
        if (option=="--legacy-construction") legacyConstruction=true;
        else if (option=="--verify-source-lineage") verifySourceLineage=true;
        else break;
        --argc;
    }
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
    const auto wallBounds = originalWalls.bounds();
    const Domain2D domain{{{wallBounds.min.x - padding, wallBounds.min.y - padding},
                           {wallBounds.max.x + padding, wallBounds.max.y + padding}}};
    HybridMeshPolicy2D hybridPolicy;
    hybridPolicy.sharedIntersectionConstruction=!legacyConstruction;
    hybridPolicy.verifySourceLineageOracle=verifySourceLineage;
    auto robust=buildRobustH4Mesh2D(
        chains,layerParameters,domain,originalWalls,maxLevel,refinement,{},hybridPolicy);

    const auto parent = outputPrefix.parent_path().empty()
        ? std::filesystem::path(".") : outputPrefix.parent_path();
    std::filesystem::create_directories(parent);
    if (!robust.success()) {
        std::cerr<<"h4_status=failed mesh_mode="<<h4MeshModeName(robust.mode)
                 <<" fallback_stage="<<h4FallbackStageName(robust.fallbackStage)
                 <<" requested_layer_failure="
                 <<boundaryLayerFailureReasonName(
                       robust.requestedLayerCandidate.failure.reason)
                 <<" local_layer_failure="
                 <<boundaryLayerFailureReasonName(
                       robust.localLayerCandidate.failure.reason)
                 <<" hybrid_failure="
                 <<hybridMeshFailureReasonName(robust.hybridCandidate.failure.reason)
                 <<" hybrid_detail="<<robust.hybridCandidate.failure.message
                 <<" fallback_failure="<<robust.fallback.failureMessage<<'\n';
        return EXIT_FAILURE;
    }
    if (robust.mode==H4MeshMode2D::PureCutCellFallback) {
        const auto vtkPath=outputPrefix.string()+".fallback.vtk";
        const auto cm2dPath=outputPrefix.string()+".fallback.cm2d";
        const auto solverVtkPath=outputPrefix.string()+".fallback.solver.vtk";
        const auto solverCm2dPath=outputPrefix.string()+".fallback.solver.cm2d";
        const auto qualityPath=
            outputPrefix.string()+".fallback.construction-quality.json";
        const auto solverQualityPath=outputPrefix.string()+
                                     ".fallback.solver-quality.json";
        const auto& fallback=robust.fallback;
        if (!writeLegacyVtk2D(fallback.topology,vtkPath,&error) ||
            !writeCm2dTopology(fallback.topology,cm2dPath,&error) ||
            !writeLegacyVtk2D(fallback.solverTopology,solverVtkPath,&error) ||
            !writeCm2dTopology(fallback.solverTopology,solverCm2dPath,&error) ||
            !writeText(qualityPath,qualityReportToJson(fallback.meshQuality),error) ||
            !writeText(solverQualityPath,
                       solverQualityReportToJson(fallback.solverQuality),error)) {
            std::cerr<<error<<'\n';
            return EXIT_FAILURE;
        }
        std::string openFoamStatus="not_requested";
        if (argc==12) {
            double thickness=0.0;
            if (!parseDouble(argv[11],thickness) || thickness<=0.0) {
                std::cerr<<"extrusion thickness must be finite and positive\n";
                return EXIT_FAILURE;
            }
            const auto foam=writeExtrudedOpenFoam2D(
                fallback.solverTopology,domain,originalWalls,argv[10],thickness,&error);
            if (!foam.valid()) {
                std::cerr<<"OpenFOAM output failed: "<<error<<'\n';
                return EXIT_FAILURE;
            }
            openFoamStatus="written";
        }
        std::cout<<"h4_status=success mesh_mode=pure_cutcell_fallback"
                 <<" fallback_stage="<<h4FallbackStageName(robust.fallbackStage)
                 <<" hybrid_detail="<<robust.hybridCandidate.failure.message
                 <<" solver_cells="<<fallback.solverTopology.cells.size()
                 <<" area_error="<<fallback.areaError
                 <<" solver_quality=pass openfoam="<<openFoamStatus
                 <<" vtk="<<vtkPath<<'\n';
        return EXIT_SUCCESS;
    }
    auto hybrid=std::move(robust.hybridCandidate);
    if (hybrid.solverTopology.constructionRegistry &&
        !writeText(outputPrefix.string()+".hybrid.construction.json",
            intersectionConstructionToJson(*hybrid.solverTopology.constructionRegistry,
                hybrid.solverTopology.canonicalVertexIds,
                hybrid.solverTopology.sharedPartitionCount,
                hybrid.solverTopology.sharedPartitionCacheHits),error)) {
        std::cerr<<error<<'\n';return EXIT_FAILURE;
    }
    if (!writeText(outputPrefix.string()+".hybrid.intersections.json",
                   intersectionRecordsToJson(hybrid.canonicalizedIntersections),error)) {
        std::cerr<<error<<'\n';
        return EXIT_FAILURE;
    }
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
    const auto solverVtkPath = outputPrefix.string() + ".hybrid.solver.vtk";
    const auto solverCm2dPath = outputPrefix.string() + ".hybrid.solver.cm2d";
    const auto qualityPath =
        outputPrefix.string() + ".hybrid.construction-quality.json";
    const auto solverQualityPath = outputPrefix.string() +
                                   ".hybrid.solver-quality.json";
    const auto qualityContractPath = outputPrefix.string() +
                                     ".hybrid.quality-contract.json";
    if (!writeHybridLegacyVtk2D(hybrid, vtkPath, &error) ||
        !writeCm2dTopology(hybrid.topology, cm2dPath, &error) ||
        !writeLegacyVtk2D(hybrid.solverTopology, solverVtkPath, &error) ||
        !writeCm2dTopology(hybrid.solverTopology, solverCm2dPath, &error) ||
        !writeText(qualityPath, qualityReportToJson(hybrid.meshQuality), error) ||
        !writeText(solverQualityPath,
                   solverQualityReportToJson(hybrid.solverQuality), error) ||
        !writeText(qualityContractPath,
                   qualityContractReportToJson(hybrid.qualityContract),error)) {
        std::cerr << error << '\n';
        return EXIT_FAILURE;
    }
    const auto readback = readCm2dTopology(cm2dPath);
    const auto solverReadback = readCm2dTopology(solverCm2dPath);
    if (!readback.valid() || !solverReadback.valid() ||
        solverReadback.topology.cells.size() != hybrid.solverTopology.cells.size()) {
        std::cerr << "CM2D readback failed for hybrid or solver topology\n";
        return EXIT_FAILURE;
    }

    std::string openFoamStatus = "not_requested";
    if (argc == 12) {
        double thickness = 0.0;
        if (!parseDouble(argv[11], thickness) || thickness <= 0.0) {
            std::cerr << "extrusion thickness must be finite and positive\n";
            return EXIT_FAILURE;
        }
        if (!hybrid.solverQuality.valid() || !hybrid.solverTopology.valid()) {
            std::cerr << "hybrid solver-quality gate failed; refusing OpenFOAM output\n";
            return EXIT_FAILURE;
        }
        const auto foam = writeExtrudedOpenFoam2D(
            hybrid.solverTopology, domain, originalWalls, argv[10], thickness, &error);
        if (!foam.valid()) {
            std::cerr << "OpenFOAM output failed: " << error << '\n';
            return EXIT_FAILURE;
        }
        openFoamStatus = "written";
    }

    std::cout << "hybrid_status=success cells=" << hybrid.metrics.unifiedCellCount
              << " solver_cells=" << hybrid.metrics.solverCellCount
              << " layer_cells=" << hybrid.metrics.boundaryLayerCellCount
              << " remainder_cut=" << hybrid.metrics.remainderCutCellCount
              << " remainder_cartesian="
              << hybrid.metrics.remainderCartesianCellCount
              << " transition_rings=" << hybrid.metrics.transitionRingCount
              << " transition_final_subdivision="
              << hybrid.metrics.transitionFinalTangentialSubdivision
              << " transition_target_h=" << hybrid.metrics.transitionTargetCellSize
              << " transition_ring_thickness="
              << hybrid.metrics.transitionRingThickness
              << " interface_edges=" << hybrid.interfaceAudit.interfaceEdgeCount
              << " interface_vertices=" << hybrid.interfaceAudit.interfaceVertexCount
              << " area_error=" << hybrid.metrics.areaError
              << " lineage_checks="
              << hybrid.sourceLineageAudit.lineageCandidateChecks
              << " lineage_oracle_checks="
              << hybrid.sourceLineageAudit.oracleCandidateChecks
              << " lineage_mismatches="
              << hybrid.sourceLineageAudit.mismatchedCells
              << " solver_quality="
              << (hybrid.solverQuality.valid() ? "pass" : "fail")
              << " quality_contract="
              << qualityContractStatusName(hybrid.qualityContract.status())
              << " openfoam=" << openFoamStatus
              << " vtk=" << vtkPath
              << " solver_vtk=" << solverVtkPath
              << " solver_cm2d=" << solverCm2dPath
              << " report=" << jsonPath << '\n';
    return EXIT_SUCCESS;
}
