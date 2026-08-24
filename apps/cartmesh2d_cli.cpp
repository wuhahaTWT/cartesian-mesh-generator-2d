#include "cartmesh2d/io/MeshIO2D.hpp"
#include "cartmesh2d/io/OpenFoam2D.hpp"
#include "cartmesh2d/quality/Quality2D.hpp"
#include "cartmesh2d/quality/SolverQuality2D.hpp"
#include "cartmesh2d/quality/SolverTopology2D.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace cartmesh2d;

namespace {

bool readBoundaryFile(const std::filesystem::path& path,
                      std::vector<std::vector<Point2D>>& loops,
                      std::string& error) {
    std::ifstream in(path);
    if (!in) {
        error = "cannot open boundary file: " + path.string();
        return false;
    }
    std::string line;
    std::size_t lineNumber = 0;
    std::vector<Point2D> points;
    while (std::getline(in, line)) {
        ++lineNumber;
        const auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) {
            if (!points.empty()) {
                if (points.size()<3) {
                    error="each boundary loop must contain at least three vertices";
                    return false;
                }
                loops.push_back(std::move(points));
                points.clear();
            }
            continue;
        }
        if (line[first] == '#') continue;
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
    if (!points.empty()) loops.push_back(std::move(points));
    if (loops.empty()) {
        error = "boundary file must contain at least one loop";
        return false;
    }
    for (const auto& loop:loops) {
        if (loop.size()<3) {
            error="each boundary loop must contain at least three vertices";
            return false;
        }
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

const char* cutKindName(CutCellKind kind) noexcept {
    switch (kind) {
    case CutCellKind::Empty: return "empty";
    case CutCellKind::Full: return "full";
    case CutCellKind::Cut: return "cut";
    case CutCellKind::Unsupported: return "unsupported";
    }
    return "unknown";
}

const char* fluidRegionName(FluidRegion2D region) noexcept {
    return region == FluidRegion2D::Exterior ? "exterior" : "interior";
}

bool parseFluidRegion(const std::string& value, FluidRegion2D& region) {
    if (value == "exterior" || value == "outside" || value == "external") {
        region = FluidRegion2D::Exterior;
        return true;
    }
    if (value == "interior" || value == "inside" || value == "internal") {
        region = FluidRegion2D::Interior;
        return true;
    }
    return false;
}

const char* smallStatusName(SmallCellStatus2D status) noexcept {
    switch (status) {
    case SmallCellStatus2D::Stable: return "stable";
    case SmallCellStatus2D::CandidateFound: return "candidate";
    case SmallCellStatus2D::Unresolved: return "unresolved";
    }
    return "unknown";
}

const char* topologyIssueName(TopologyIssueCode2D code) noexcept {
    switch (code) {
    case TopologyIssueCode2D::InvalidCell: return "InvalidCell";
    case TopologyIssueCode2D::DegenerateEdge: return "DegenerateEdge";
    case TopologyIssueCode2D::DuplicateCellSource: return "DuplicateCellSource";
    case TopologyIssueCode2D::OrphanInternalEdge: return "OrphanInternalEdge";
    case TopologyIssueCode2D::NonManifoldEdge: return "NonManifoldEdge";
    case TopologyIssueCode2D::UnclassifiedBoundaryEdge: return "UnclassifiedBoundaryEdge";
    case TopologyIssueCode2D::OpenCellLoop: return "OpenCellLoop";
    case TopologyIssueCode2D::AreaMismatch: return "AreaMismatch";
    }
    return "Unknown";
}

void printTopologyDiagnostics(const TopologyMesh2D& topology) {
    const auto& audit = topology.audit;
    std::cerr << "topology_audit"
              << " duplicateVertices=" << audit.duplicateVertices
              << " duplicateEdges=" << audit.duplicateEdges
              << " orphanInternalEdges=" << audit.orphanInternalEdges
              << " nonManifoldEdges=" << audit.nonManifoldEdges
              << " unclassifiedBoundaryEdges=" << audit.unclassifiedBoundaryEdges
              << " openCellLoops=" << audit.openCellLoops
              << " areaMismatches=" << audit.areaMismatches << '\n';
    std::cerr << "topology_counts vertices=" << topology.vertices.size()
              << " edges=" << topology.edges.size()
              << " cells=" << topology.cells.size()
              << " issues=" << topology.issues.size() << '\n';
    const std::size_t limit = std::min<std::size_t>(topology.issues.size(), 32);
    for (std::size_t i = 0; i < limit; ++i) {
        const auto& issue = topology.issues[i];
        std::cerr << "topology_issue[" << i << "] code=" << topologyIssueName(issue.code)
                  << " object=" << issue.objectId
                  << " message=" << issue.message << '\n';
    }
    if (topology.issues.size() > limit) {
        std::cerr << "topology_issue_more=" << (topology.issues.size() - limit) << '\n';
    }
}

bool writeVisualizationMetadata(const std::filesystem::path& path,
                                const std::vector<CutCell2D>& cutCells,
                                const SmallCellReport2D& smallReport,
                                const AgglomerationResult2D& stabilized,
                                FluidRegion2D fluidRegion,
                                std::size_t boundaryLoopCount,
                                std::string& error) {
    using SourceKey = std::pair<std::uint64_t, std::size_t>;
    std::map<SourceKey, const SmallCellRecord2D*> smallRecords;
    for (const auto& record : smallReport.records) {
        smallRecords[{record.sourceKey, record.sourceId}] = &record;
    }

    std::ofstream out(path);
    if (!out) {
        error = "failed to open visualization metadata output";
        return false;
    }
    out << std::setprecision(17);
    out << "{\n";
    out << "  \"format\": \"cartmesh2d-viz-v1\",\n";
    out << "  \"fluid_region\": \"" << fluidRegionName(fluidRegion) << "\",\n";
    if (boundaryLoopCount>1) {
        out << "  \"boundary_loop_count\": " << boundaryLoopCount << ",\n";
        out << "  \"boundary_semantics\": \"even_odd\",\n";
    }
    out << "  \"boundary_role\": \""
        << (fluidRegion == FluidRegion2D::Exterior ? "solid_wall" : "fluid_envelope")
        << "\",\n";
    out << "  \"small_alpha_threshold\": "
        << smallReport.policy.areaFractionThreshold << ",\n";
    out << "  \"source_small_cell_count\": " << smallReport.smallCellCount << ",\n";
    out << "  \"merged_small_cell_count\": " << stabilized.mergedSmallCellCount << ",\n";
    out << "  \"source_cells\": [\n";

    bool first = true;
    for (const auto& cut : cutCells) {
        if (cut.kind == CutCellKind::Empty || cut.kind == CutCellKind::Unsupported) continue;
        const auto it = smallRecords.find({cut.sourceKey, cut.sourceId});
        const SmallCellRecord2D* record = it == smallRecords.end() ? nullptr : it->second;
        const bool isSmall = record != nullptr && record->status != SmallCellStatus2D::Stable;
        if (!first) out << ",\n";
        first = false;
        out << "    {\"source_id\": " << cut.sourceId
            << ", \"source_key\": " << cut.sourceKey
            << ", \"level\": " << (cut.sourceKey & 0x3fULL)
            << ", \"kind\": \"" << cutKindName(cut.kind) << "\""
            << ", \"area_fraction\": " << cut.areaFraction
            << ", \"small\": " << (isSmall ? "true" : "false")
            << ", \"small_status\": \""
            << (record != nullptr ? smallStatusName(record->status) : "unknown") << "\"";
        if (record != nullptr) {
            out << ", \"source_topology_cell_id\": " << record->topologyCellId;
            if (record->targetTopologyCellId) {
                out << ", \"target_topology_cell_id\": " << *record->targetTopologyCellId;
            } else {
                out << ", \"target_topology_cell_id\": null";
            }
        } else {
            out << ", \"source_topology_cell_id\": null, \"target_topology_cell_id\": null";
        }
        out << ", \"background_bounds\": ["
            << cut.backgroundBounds.min.x << ", " << cut.backgroundBounds.min.y << ", "
            << cut.backgroundBounds.max.x << ", " << cut.backgroundBounds.max.y << "]";
        if (cut.centroid) {
            out << ", \"centroid\": [" << cut.centroid->x << ", " << cut.centroid->y << "]";
        } else {
            out << ", \"centroid\": null";
        }
        out << "}";
    }
    out << "\n  ]\n}\n";
    if (!out.good()) {
        error = "failed while writing visualization metadata";
        return false;
    }
    return true;
}

void usage() {
    std::cerr << "usage: cartmesh2d_cli <boundary.xy> <output-prefix> "
                 "[max-level=5] [padding-fraction=0.25] [small-alpha=0.10] "
                 "[fluid-region=exterior|interior] [openfoam-case-dir] [minimum-level=0]\n"
                 "multiple loops: separate x-y vertex blocks with a blank line; nesting uses even-odd semantics\n"
                 "default CFD semantics: boundary.xy is a SOLID wall and fluid is EXTERIOR\n";
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 3 || argc > 9) {
        usage();
        return EXIT_FAILURE;
    }

    const std::filesystem::path boundaryPath = argv[1];
    const std::filesystem::path outputPrefix = argv[2];
    std::size_t maxLevel = 5;
    std::size_t minimumLevel = 0;
    double paddingFraction = 0.25;
    double smallAlpha = 0.10;
    FluidRegion2D fluidRegion = FluidRegion2D::Exterior;
    std::optional<std::filesystem::path> openFoamCase;
    try {
        if (argc >= 4) maxLevel = static_cast<std::size_t>(std::stoul(argv[3]));
        if (argc >= 5) paddingFraction = std::stod(argv[4]);
        if (argc >= 6) smallAlpha = std::stod(argv[5]);
        if (argc >= 9) minimumLevel = static_cast<std::size_t>(std::stoul(argv[8]));
    } catch (const std::exception&) {
        std::cerr << "invalid numeric CLI argument\n";
        return EXIT_FAILURE;
    }
    if (argc >= 7 && !parseFluidRegion(argv[6], fluidRegion)) {
        std::cerr << "invalid fluid-region; expected exterior or interior\n";
        return EXIT_FAILURE;
    }
    if (argc >= 8) openFoamCase=std::filesystem::path(argv[7]);
    if (maxLevel == 0 || maxLevel > 28 || minimumLevel > maxLevel || !(paddingFraction > 0.0) ||
        !(smallAlpha > 0.0 && smallAlpha < 1.0)) {
        std::cerr << "invalid max-level, minimum-level, padding-fraction or small-alpha\n";
        return EXIT_FAILURE;
    }

    std::vector<std::vector<Point2D>> loopPoints;
    std::string error;
    if (!readBoundaryFile(boundaryPath, loopPoints, error)) {
        std::cerr << error << '\n';
        return EXIT_FAILURE;
    }
    std::vector<BoundaryLoop> loops;
    loops.reserve(loopPoints.size());
    for (auto& points:loopPoints) loops.emplace_back(std::move(points));
    BoundaryRegion2D boundary(std::move(loops));
    const auto diagnostics = boundary.diagnose();
    if (!diagnostics.valid()) {
        std::cerr << "boundary diagnostics failed with " << diagnostics.issues.size()
                  << " issue(s)\n";
        for (std::size_t i=0;i<diagnostics.issues.size();++i) {
            const auto& issue=diagnostics.issues[i];
            std::cerr<<"boundary_region_issue["<<i<<"] code="
                     <<static_cast<int>(issue.code)
                     <<" loop_a="<<issue.loopA<<" loop_b="<<issue.loopB
                     <<" message="<<issue.message<<'\n';
        }
        return EXIT_FAILURE;
    }
    if (!boundary.normalizeAlternating()) {
        std::cerr << "failed to normalize boundary-region nesting/orientation\n";
        return EXIT_FAILURE;
    }

    const AABB2D bounds = boundary.bounds();
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
    refinement.minimumLevel = minimumLevel;
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
    std::size_t splitFluidLeaves = 0;
    std::size_t nextSourceId = 0;
    double sourceFluidArea = 0.0;
    for (const auto& leaf : tree.leaves()) {
        auto components = buildCutCells(leaf, boundary, fluidRegion);
        if (components.size() > 1) ++splitFluidLeaves;
        for (auto& cut : components) {
            // A physical leaf can legally contribute multiple disconnected
            // solver cells. Give every emitted component a unique deterministic
            // source id while preserving the Quadtree key/level in sourceKey.
            cut.sourceId = nextSourceId++;
            cut.sourceKey = leaf.key;
            if (!cut.valid() && cut.kind == CutCellKind::Unsupported) ++unsupported;
            if (cut.kind != CutCellKind::Empty && cut.kind != CutCellKind::Unsupported) {
                sourceFluidArea += cut.area;
            }
            cutCells.push_back(std::move(cut));
        }
    }
    if (unsupported != 0) {
        std::cerr << "Cut-cell construction produced " << unsupported
                  << " unsupported leaf geometry case(s)\n";
        return EXIT_FAILURE;
    }

    // Product-level physics gate. This catches the exact class of error where
    // a mathematically self-consistent mesh is generated on the wrong side of
    // the solid boundary.
    const double domainArea = domain.width() * domain.height();
    const double solidArea = boundary.area();
    const double expectedFluidArea = fluidRegion == FluidRegion2D::Exterior
        ? domainArea - solidArea
        : solidArea;
    const TolerancePolicy tol{};
    const double areaEps = std::max(tol.absolute * tol.absolute,
                                    tol.relative * std::max(1.0, std::abs(expectedFluidArea)));
    if (std::abs(sourceFluidArea - expectedFluidArea) > areaEps) {
        std::cerr << "fluid-side physics gate failed: region=" << fluidRegionName(fluidRegion)
                  << " generated_area=" << std::setprecision(17) << sourceFluidArea
                  << " expected_area=" << expectedFluidArea << '\n';
        return EXIT_FAILURE;
    }

    const TopologyMesh2D sourceTopology = buildGlobalTopology(cutCells, domain, boundary);
    if (!sourceTopology.valid()) {
        std::cerr << "source global topology audit failed\n";
        printTopologyDiagnostics(sourceTopology);
        return EXIT_FAILURE;
    }

    if (fluidRegion == FluidRegion2D::Exterior) {
        std::size_t embeddedEdges = 0;
        std::size_t domainEdges = 0;
        for (const auto& edge : sourceTopology.edges) {
            if (edge.patch == BoundaryPatch2D::EmbeddedBoundary) ++embeddedEdges;
            if (edge.patch == BoundaryPatch2D::DomainBoundary) ++domainEdges;
        }
        if (embeddedEdges == 0 || domainEdges == 0) {
            std::cerr << "external CFD boundary gate failed: embedded_wall_edges="
                      << embeddedEdges << " domain_boundary_edges=" << domainEdges << '\n';
            return EXIT_FAILURE;
        }
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
        for (std::size_t i = 0; i < stabilized.issues.size(); ++i) {
            const auto& issue = stabilized.issues[i];
            std::cerr << "agglomeration_issue[" << i << "] code="
                      << static_cast<int>(issue.code)
                      << " object=" << issue.objectId
                      << " message=" << issue.message << '\n';
        }
        printTopologyDiagnostics(stabilized.topology);
        return EXIT_FAILURE;
    }

    const MeshQualityReport2D quality =
        evaluateMeshQuality(stabilized.topology, cutCells, &smallReport);
    if (!quality.valid()) {
        std::cerr << "final quality/topology report is invalid\n";
        const std::size_t issueLimit = std::min<std::size_t>(quality.issues.size(), 20);
        for (std::size_t i = 0; i < issueLimit; ++i) {
            const auto& issue = quality.issues[i];
            std::cerr << "quality_issue[" << i << "] code="
                      << static_cast<int>(issue.code)
                      << " object=" << issue.objectId
                      << " message=" << issue.message << '\n';
        }
        if (quality.issues.size() > issueLimit) {
            std::cerr << "quality_issue_omitted=" << quality.issues.size() - issueLimit << '\n';
        }
        printTopologyDiagnostics(stabilized.topology);
        return EXIT_FAILURE;
    }

    std::optional<SolverQualityReport2D> solverQuality;
    std::optional<SolverTopologyResult2D> solverTopology;
    if (openFoamCase) {
        solverTopology=buildSolverTopology2D(stabilized.topology,domain,boundary);
        if (!solverTopology->valid()) {
            std::cerr<<"solver topology partition failed";
            for (const auto& issue:solverTopology->issues) std::cerr<<": "<<issue;
            std::cerr<<'\n';
            return EXIT_FAILURE;
        }
        solverQuality=evaluateSolverQuality2D(solverTopology->topology);
        if (!solverQuality->valid()) {
            std::cerr<<"solver-quality gate failed with "<<solverQuality->issues.size()
                     <<" issue(s); OpenFOAM mesh was not written\n";
            const std::size_t limit=std::min<std::size_t>(solverQuality->issues.size(),20);
            for (std::size_t i=0;i<limit;++i) {
                const auto& issue=solverQuality->issues[i];
                std::cerr<<"solver_quality_issue["<<i<<"] code="
                         <<static_cast<int>(issue.code)<<" cell="<<issue.cellId
                         <<" edge="<<issue.edgeId<<" measured="<<issue.measured
                         <<" limit="<<issue.limit<<" message="<<issue.message<<'\n';
                if (issue.cellId<solverTopology->topology.cells.size()) {
                    std::cerr<<"solver_quality_cell_vertices["<<i<<"]=";
                    for (const auto vertex:solverTopology->topology.cells[issue.cellId].vertices) {
                        const auto& point=solverTopology->topology.vertices[vertex].point;
                        std::cerr<<'('<<std::setprecision(17)<<point.x<<','<<point.y<<')';
                    }
                    std::cerr<<'\n';
                }
                if (issue.edgeId<solverTopology->topology.edges.size()) {
                    const auto& edge=solverTopology->topology.edges[issue.edgeId];
                    const auto& a=solverTopology->topology.vertices[edge.v0].point;
                    const auto& b=solverTopology->topology.vertices[edge.v1].point;
                    std::cerr<<"solver_quality_edge["<<i<<"]=("<<a.x<<','<<a.y<<")-("
                             <<b.x<<','<<b.y<<") owner="<<edge.owner<<" neighbour=";
                    if (edge.neighbour) std::cerr<<*edge.neighbour; else std::cerr<<"none";
                    std::cerr<<'\n';
                    if (edge.neighbour && *edge.neighbour<solverTopology->topology.cells.size()) {
                        std::cerr<<"solver_quality_neighbour_vertices["<<i<<"]=";
                        for (const auto vertex:
                             solverTopology->topology.cells[*edge.neighbour].vertices) {
                            const auto& point=solverTopology->topology.vertices[vertex].point;
                            std::cerr<<'('<<point.x<<','<<point.y<<')';
                        }
                        std::cerr<<'\n';
                    }
                }
            }
            return EXIT_FAILURE;
        }
    }

    const auto parent = outputPrefix.parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    const std::filesystem::path vtkPath = outputPrefix.string() + ".vtk";
    const std::filesystem::path cm2dPath = outputPrefix.string() + ".cm2d";
    const std::filesystem::path qualityPath = outputPrefix.string() + ".quality.json";
    const std::filesystem::path vizPath = outputPrefix.string() + ".viz.json";

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
    error.clear();
    if (!writeVisualizationMetadata(vizPath, cutCells, smallReport, stabilized,
                                    fluidRegion, boundary.loops().size(), error)) {
        std::cerr << error << '\n';
        return EXIT_FAILURE;
    }

    std::optional<OpenFoamWriteReport2D> openFoamReport;
    if (openFoamCase) {
        const double extrusionThickness=std::ldexp(span,-static_cast<int>(maxLevel));
        error.clear();
        openFoamReport=writeExtrudedOpenFoam2D(solverTopology->topology,domain,boundary,
                                               *openFoamCase,extrusionThickness,&error);
        if (!openFoamReport->valid()) {
            std::cerr<<error<<'\n';
            return EXIT_FAILURE;
        }
        std::ofstream out(*openFoamCase/"solver_quality.json");
        if (!out) {
            std::cerr<<"failed to open solver-quality JSON output\n";
            return EXIT_FAILURE;
        }
        out<<solverQualityReportToJson(*solverQuality);
        if (!out.good()) {
            std::cerr<<"failed while writing solver-quality JSON output\n";
            return EXIT_FAILURE;
        }
    }

    const MeshReadback2D readback = readCm2dTopology(cm2dPath);
    if (!readback.valid() || !sameReadbackTopology(stabilized.topology, readback.topology)) {
        std::cerr << "independent CM2D read-back verification failed: " << readback.error << '\n';
        return EXIT_FAILURE;
    }

    std::cout << "cartmesh2d end-to-end PASS\n"
              << "fluid_region=" << fluidRegionName(fluidRegion) << '\n'
              << "boundary_role="
              << (fluidRegion == FluidRegion2D::Exterior ? "solid_wall" : "fluid_envelope") << '\n'
              << "boundary_loops=" << boundary.loops().size() << '\n'
              << "minimum_level=" << minimumLevel << '\n'
              << "leaf_count=" << tree.leaves().size() << '\n'
              << "split_fluid_leaves=" << splitFluidLeaves << '\n'
              << "source_cells=" << sourceTopology.cells.size() << '\n'
              << "source_fluid_area=" << sourceFluidArea << '\n'
              << "expected_fluid_area=" << expectedFluidArea << '\n'
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
              << "quality_json=" << qualityPath.string() << '\n'
              << "visualization_json=" << vizPath.string() << '\n';
    if (openFoamReport) {
        std::cout<<"solver_quality=PASS\n"
                 <<"solver_max_nonorthogonality_deg="<<solverQuality->maxNonOrthogonalityDeg<<'\n'
                 <<"solver_max_internal_skewness="<<solverQuality->maxInternalSkewness<<'\n'
                 <<"solver_max_cell_aspect="<<solverQuality->maxCellAspect<<'\n'
                 <<"solver_min_face_length="<<solverQuality->minFaceLength<<'\n'
                 <<"solver_min_face_weight="<<solverQuality->minFaceWeight<<'\n'
                 <<"solver_min_volume_ratio="<<solverQuality->minVolumeRatio<<'\n'
                 <<"solver_partitioned_input_cells="<<solverTopology->partitionedCellCount<<'\n'
                 <<"solver_output_cells="<<solverTopology->outputCellCount<<'\n'
                 <<"openfoam_case="<<openFoamCase->string()<<'\n'
                 <<"openfoam_cells="<<openFoamReport->cellCount<<'\n'
                 <<"openfoam_faces="<<openFoamReport->faceCount<<'\n';
    }
    return EXIT_SUCCESS;
}
