#include "cartmesh2d/io/MeshIO2D.hpp"
#include "cartmesh2d/io/OpenFoam2D.hpp"
#include "cartmesh2d/geometry/BoundarySimplification2D.hpp"
#include "cartmesh2d/quality/Quality2D.hpp"
#include "cartmesh2d/quality/SolverQuality2D.hpp"
#include "cartmesh2d/quality/SolverTopology2D.hpp"
#include "cartmesh2d/sizing/SizeField2D.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace cartmesh2d;

namespace {

// W1's geometric grid-corner weld budget. Q1 rejects faces shorter than 0.01h;
// this remains two orders of magnitude below that accepted-face scale.
constexpr double gridCornerWeldFraction = 1.0e-4;

bool readBoundaryFile(const std::filesystem::path& path,
                      std::vector<std::vector<Point2D>>& loops,
                      std::vector<EmbeddedBoundaryPatch2D>& embeddedPatches,
                      std::string& error) {
    std::ifstream in(path);
    if (!in) {
        error = "cannot open boundary file: " + path.string();
        return false;
    }
    std::string line;
    std::size_t lineNumber = 0;
    std::vector<Point2D> points;
    std::optional<EmbeddedBoundaryPatch2D> pendingPatch;
    const auto finishLoop=[&]() {
        if (points.empty()) return true;
        if (points.size()<3) {
            error="each boundary loop must contain at least three vertices";
            return false;
        }
        const std::size_t loopId=loops.size();
        loops.push_back(std::move(points));
        points.clear();
        embeddedPatches.push_back(pendingPatch.value_or(EmbeddedBoundaryPatch2D{
            "wall_"+std::to_string(loopId),"wall",BoundaryConditionRole2D::Wall,""}));
        pendingPatch.reset();
        return true;
    };
    while (std::getline(in, line)) {
        ++lineNumber;
        const auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) {
            if (!finishLoop()) return false;
            continue;
        }
        if (line[first] == '#') {
            std::istringstream metadata(line.substr(first+1));
            std::string marker;
            metadata>>marker;
            if (marker!="cartmesh2d-loop") continue;
            if (!points.empty() || pendingPatch) {
                error="boundary metadata must appear once before its loop";
                return false;
            }
            EmbeddedBoundaryPatch2D patch;
            std::string key,role;
            if (!(metadata>>key) || key!="patch" || !(metadata>>std::quoted(patch.name)) ||
                !(metadata>>key) || key!="type" || !(metadata>>std::quoted(patch.type)) ||
                !(metadata>>key) || key!="role" || !(metadata>>std::quoted(role)) ||
                !(metadata>>key) || key!="layer" ||
                !(metadata>>std::quoted(patch.sourceLayer)) ||
                !parseBoundaryConditionRole(role,patch.role)) {
                error="invalid cartmesh2d-loop metadata on line "+
                      std::to_string(lineNumber);
                return false;
            }
            std::string trailing;
            if (metadata>>trailing) {
                error="unexpected boundary metadata token on line "+
                      std::to_string(lineNumber);
                return false;
            }
            pendingPatch=std::move(patch);
            continue;
        }
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
    if (!finishLoop()) return false;
    if (pendingPatch) {
        error="boundary metadata has no following loop";
        return false;
    }
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

bool parseExactDouble(const std::string& value, double& parsed) {
    try {
        std::size_t consumed=0;
        parsed=std::stod(value,&consumed);
        return consumed==value.size() && std::isfinite(parsed);
    } catch (const std::exception&) {
        return false;
    }
}

bool parseExactLevel(const std::string& value, std::size_t& parsed) {
    try {
        std::size_t consumed=0;
        const auto raw=std::stoull(value,&consumed);
        if (consumed!=value.size()) return false;
        parsed=static_cast<std::size_t>(raw);
        return static_cast<unsigned long long>(parsed)==raw;
    } catch (const std::exception&) {
        return false;
    }
}

bool parseSizingOptions(int argc, char** argv, int first,
                        std::vector<DistanceRefinementBand2D>& distanceBands,
                        std::vector<BoxRefinementRegion2D>& boxRegions,
                        std::optional<SizeFieldPolicy2D>& sizeField,
                        std::string& error) {
    // The size-field options describe intent (wall resolution, growth, far field)
    // and are compiled onto the primitives above.  They stay opt-in so every
    // existing command line keeps producing byte-identical output.
    const auto requireSizeField = [&]() -> SizeFieldPolicy2D& {
        if (!sizeField) sizeField.emplace();
        return *sizeField;
    };
    int index=first;
    while (index<argc) {
        const std::string option=argv[index++];
        if (option=="--distance-band") {
            if (index+2>argc) {
                error="--distance-band requires <distance> <target-level>";
                return false;
            }
            DistanceRefinementBand2D band;
            if (!parseExactDouble(argv[index],band.distance) ||
                !parseExactLevel(argv[index+1],band.targetLevel)) {
                error="invalid --distance-band value";
                return false;
            }
            index+=2;
            distanceBands.push_back(band);
            continue;
        }
        if (option=="--refine-box") {
            if (index+5>argc) {
                error="--refine-box requires <xmin> <ymin> <xmax> <ymax> <target-level>";
                return false;
            }
            BoxRefinementRegion2D region;
            if (!parseExactDouble(argv[index],region.bounds.min.x) ||
                !parseExactDouble(argv[index+1],region.bounds.min.y) ||
                !parseExactDouble(argv[index+2],region.bounds.max.x) ||
                !parseExactDouble(argv[index+3],region.bounds.max.y) ||
                !parseExactLevel(argv[index+4],region.targetLevel)) {
                error="invalid --refine-box value";
                return false;
            }
            index+=5;
            boxRegions.push_back(region);
            continue;
        }
        if (option=="--size-field") {
            (void)requireSizeField();
            continue;
        }
        if (option=="--far-field-spans" || option=="--wall-cells-per-span") {
            if (index+1>argc) { error=option+" requires <value>"; return false; }
            double value=0.0;
            if (!parseExactDouble(argv[index],value)) { error="invalid "+option+" value"; return false; }
            ++index;
            if (option=="--far-field-spans") requireSizeField().farFieldSpans=value;
            else requireSizeField().wallCellsPerSpan=value;
            continue;
        }
        if (option=="--cells-per-level" || option=="--far-level") {
            if (index+1>argc) { error=option+" requires <value>"; return false; }
            std::size_t value=0;
            if (!parseExactLevel(argv[index],value)) { error="invalid "+option+" value"; return false; }
            ++index;
            auto& field=requireSizeField();
            if (!field.wallDistance) field.wallDistance.emplace();
            if (option=="--cells-per-level") field.wallDistance->cellsPerLevel=value;
            else field.wallDistance->farLevel=value;
            continue;
        }
        if (option=="--max-safe-wall-level") {
            if (index+1>argc) { error=option+" requires <value>"; return false; }
            std::size_t value=0;
            if (!parseExactLevel(argv[index],value)) { error="invalid "+option+" value"; return false; }
            ++index;
            requireSizeField().maxSafeWallLevel=value;
            continue;
        }
        if (option=="--allow-unsafe-wall-level") {
            requireSizeField().allowUnsafeWallLevel=true;
            continue;
        }
        if (option=="--curvature-cells-per-radius") {            if (index+1>argc) { error=option+" requires <value>"; return false; }
            double value=0.0;
            if (!parseExactDouble(argv[index],value)) { error="invalid "+option+" value"; return false; }
            ++index;
            requireSizeField().curvature=CurvatureSizing2D{value};
            continue;
        }
        if (option=="--gap-cells") {
            if (index+1>argc) { error=option+" requires <value>"; return false; }
            double value=0.0;
            if (!parseExactDouble(argv[index],value)) { error="invalid "+option+" value"; return false; }
            ++index;
            ProximitySizing2D proximity;
            proximity.cellsAcrossGap=value;
            requireSizeField().proximity=proximity;
            continue;
        }
        if (option=="--wake") {
            if (index+4>argc) {
                error="--wake requires <angle-deg> <downstream-spans> <half-width-spans> <levels-below-wall>";
                return false;
            }
            WakeSizing2D wake;
            if (!parseExactDouble(argv[index],wake.angleOfAttackDeg) ||
                !parseExactDouble(argv[index+1],wake.downstreamSpans) ||
                !parseExactDouble(argv[index+2],wake.halfWidthSpans) ||
                !parseExactLevel(argv[index+3],wake.levelsBelowWall)) {
                error="invalid --wake value";
                return false;
            }
            index+=4;
            requireSizeField().wake=wake;
            continue;
        }
        error="unknown sizing option: "+option;
        return false;
    }
    return true;
}

bool boxesOverlapPositive(const AABB2D& a,const AABB2D& b) noexcept {
    return std::max(a.min.x,b.min.x)<std::min(a.max.x,b.max.x) &&
           std::max(a.min.y,b.min.y)<std::min(a.max.y,b.max.y);
}

bool writeSizingFieldReport(const std::filesystem::path& path,
                            const Domain2D& domain,
                            const QuadtreeRefinementPolicy2D& policy,
                            const Quadtree2D& tree,
                            const QuadtreeBalanceReport2D& balance,
                            std::string& error) {
    std::map<std::size_t,std::size_t> levelCounts;
    for (const auto& leaf:tree.leaves()) ++levelCounts[leaf.level];
    std::ofstream out(path);
    if (!out) {
        error="failed to open sizing-field JSON output";
        return false;
    }
    out<<std::setprecision(17)
       <<"{\n  \"format\": \"cartmesh2d-sizing-v1\",\n"
       <<"  \"domain\": ["<<domain.bounds.min.x<<", "<<domain.bounds.min.y<<", "
       <<domain.bounds.max.x<<", "<<domain.bounds.max.y<<"],\n"
       <<"  \"minimum_level\": "<<policy.minimumLevel<<",\n"
       <<"  \"boundary_level\": "<<policy.boundaryLevel<<",\n"
       <<"  \"distance_bands\": [";
    for (std::size_t i=0;i<policy.distanceBands.size();++i) {
        if (i!=0) out<<", ";
        const auto& band=policy.distanceBands[i];
        out<<"{\"distance\": "<<band.distance<<", \"target_level\": "
           <<band.targetLevel<<"}";
    }
    out<<"],\n  \"box_regions\": [";
    for (std::size_t i=0;i<policy.boxRegions.size();++i) {
        if (i!=0) out<<", ";
        const auto& region=policy.boxRegions[i];
        std::size_t overlappingLeaves=0;
        for (const auto& leaf:tree.leaves()) {
            if (boxesOverlapPositive(leaf.bounds,region.bounds)) ++overlappingLeaves;
        }
        out<<"{\"bounds\": ["<<region.bounds.min.x<<", "<<region.bounds.min.y<<", "
           <<region.bounds.max.x<<", "<<region.bounds.max.y<<"], \"target_level\": "
           <<region.targetLevel<<", \"overlapping_leaf_count\": "<<overlappingLeaves<<"}";
    }
    out<<"],\n  \"level_histogram\": {";
    bool first=true;
    for (const auto& [level,count]:levelCounts) {
        if (!first) out<<", ";
        first=false;
        out<<'\"'<<level<<"\": "<<count;
    }
    out<<"},\n  \"leaf_count\": "<<tree.leaves().size()<<",\n"
       <<"  \"balance\": {\"iterations\": "<<balance.iterations
       <<", \"refined_leaves\": "<<balance.refinedLeaves
       <<", \"violations_before\": "<<balance.violationsBefore
       <<", \"violations_after\": "<<balance.violationsAfter<<"}\n}\n";
    if (!out.good()) {
        error="failed while writing sizing-field JSON output";
        return false;
    }
    return true;
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

void usage(std::ostream& out = std::cerr) {
    out << "usage: cartmesh2d_cli <boundary.xy> <output-prefix> "
                 "[max-level=5] [padding-fraction=0.25] [small-alpha=0.10] "
                 "[fluid-region=exterior|interior] [openfoam-case-dir|-] [minimum-level=0] "
                 "[boundary-simplify-cell-fraction=0] "
                 "[--distance-band <distance> <target-level>]... "
                 "[--refine-box <xmin> <ymin> <xmax> <ymax> <target-level>]...\n"
                 "multiple loops: separate x-y vertex blocks with a blank line; nesting uses even-odd semantics\n"
                 "a downstream --refine-box is the deterministic rectangular wake sizing primitive\n"
                 "default CFD semantics: boundary.xy is a SOLID wall and fluid is EXTERIOR\n"
                 "\n"
                 "size field (opt-in; takes over the domain and the tree depth):\n"
                 "  --size-field                          enable with defaults\n"
                 "  --far-field-spans <v>                 domain half-extent in body spans (10)\n"
                 "  --wall-cells-per-span <v>             body span / wall cell size (128)\n"
                 "  --cells-per-level <n>                 cells per level band, 0 = boundary-only (3)\n"
                 "  --far-level <n>                       global level floor (0)\n"
                 "  --curvature-cells-per-radius <v>      raise the level where the wall turns\n"
                 "  --gap-cells <v>                       cells guaranteed across a facing gap\n"
                 "  --wake <angle-deg> <downstream-spans> <half-width-spans> <levels-below-wall>\n"
                 "  --max-safe-wall-level <n>             refuse to resolve deeper than this (11)\n"
                 "  --allow-unsafe-wall-level             opt in past that measured ceiling\n"
                 "max-level and padding-fraction are ignored while --size-field is active; the\n"
                 "wall cell size is requested physically so it no longer moves with the domain\n";
}

} // namespace

int main(int argc, char** argv) {
    const auto totalStart = std::chrono::steady_clock::now();
    const auto elapsedSeconds = [](const auto& start) {
        return std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
    };
    if (argc == 2 && (std::string(argv[1]) == "--help" ||
                      std::string(argv[1]) == "-h")) {
        usage(std::cout);
        return EXIT_SUCCESS;
    }
    if (argc < 3) {
        usage();
        return EXIT_FAILURE;
    }

    int sizingOptionStart=argc;
    for (int index=3;index<argc;++index) {
        if (std::string(argv[index]).starts_with("--")) {
            sizingOptionStart=index;
            break;
        }
    }
    if (sizingOptionStart>10) {
        usage();
        return EXIT_FAILURE;
    }

    const std::filesystem::path boundaryPath = argv[1];
    const std::filesystem::path outputPrefix = argv[2];
    std::size_t maxLevel = 5;
    std::size_t minimumLevel = 0;
    double paddingFraction = 0.25;
    double smallAlpha = 0.10;
    double boundarySimplifyCellFraction = 0.0;
    FluidRegion2D fluidRegion = FluidRegion2D::Exterior;
    std::optional<std::filesystem::path> openFoamCase;
    std::vector<DistanceRefinementBand2D> distanceBands;
    std::vector<BoxRefinementRegion2D> boxRegions;
    try {
        if (sizingOptionStart >= 4) maxLevel = static_cast<std::size_t>(std::stoul(argv[3]));
        if (sizingOptionStart >= 5) paddingFraction = std::stod(argv[4]);
        if (sizingOptionStart >= 6) smallAlpha = std::stod(argv[5]);
        if (sizingOptionStart >= 9) minimumLevel = static_cast<std::size_t>(std::stoul(argv[8]));
        if (sizingOptionStart >= 10) boundarySimplifyCellFraction = std::stod(argv[9]);
    } catch (const std::exception&) {
        std::cerr << "invalid numeric CLI argument\n";
        return EXIT_FAILURE;
    }
    if (sizingOptionStart >= 7 && !parseFluidRegion(argv[6], fluidRegion)) {
        std::cerr << "invalid fluid-region; expected exterior or interior\n";
        return EXIT_FAILURE;
    }
    if (sizingOptionStart >= 8 && std::string(argv[7])!="-") {
        openFoamCase=std::filesystem::path(argv[7]);
    }
    std::string sizingError;
    std::optional<SizeFieldPolicy2D> sizeFieldPolicy;
    if (!parseSizingOptions(argc,argv,sizingOptionStart,distanceBands,boxRegions,
                            sizeFieldPolicy,sizingError)) {
        std::cerr<<sizingError<<'\n';
        usage();
        return EXIT_FAILURE;
    }
    if (maxLevel == 0 || maxLevel > 28 || minimumLevel > maxLevel || !(paddingFraction > 0.0) ||
        !(smallAlpha > 0.0 && smallAlpha < 1.0) ||
        !std::isfinite(boundarySimplifyCellFraction) || boundarySimplifyCellFraction<0.0) {
        std::cerr << "invalid max-level, minimum-level, padding-fraction, small-alpha or "
                     "boundary-simplify-cell-fraction\n";
        return EXIT_FAILURE;
    }

    std::vector<std::vector<Point2D>> loopPoints;
    std::vector<EmbeddedBoundaryPatch2D> embeddedPatches;
    std::string error;
    if (!readBoundaryFile(boundaryPath, loopPoints, embeddedPatches, error)) {
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
    Domain2D domain{{
        {bounds.min.x - padding, bounds.min.y - padding},
        {bounds.max.x + padding, bounds.max.y + padding}
    }};

    // The size field owns the domain and the tree depth when it is requested: wall
    // resolution is asked for physically, so the level that delivers it depends on
    // how large the domain ended up being.  Without --size-field nothing here runs
    // and the positional max-level/padding-fraction remain authoritative.
    std::optional<ResolvedSizeField2D> resolvedSizeField;
    if (sizeFieldPolicy) {
        auto resolved = resolveSizeField2D(*sizeFieldPolicy, boundary);
        if (!resolved.valid()) {
            for (const auto& issue : resolved.issues) {
                std::cerr << "size_field_issue=" << issue << '\n';
            }
            return EXIT_FAILURE;
        }
        if (!distanceBands.empty() || !boxRegions.empty()) {
            std::cerr << "--size-field and explicit --distance-band/--refine-box are "
                         "two spellings of the same field; pass only one\n";
            return EXIT_FAILURE;
        }
        domain = resolved.domain;
        maxLevel = resolved.maxLevel;
        minimumLevel = resolved.refinement.minimumLevel;
        resolvedSizeField = std::move(resolved);
    }

    std::optional<BoundarySimplificationReport2D> simplificationReport;
    if (boundarySimplifyCellFraction>0.0) {
        const double finestCellScale=std::ldexp(
            std::max(domain.width(),domain.height()),-static_cast<int>(maxLevel));
        auto simplified=simplifyBoundaryRegion2D(
            boundary,boundarySimplifyCellFraction*finestCellScale);
        if (!simplified.valid()) {
            std::cerr<<"boundary simplification failed: "<<simplified.error<<'\n';
            return EXIT_FAILURE;
        }
        simplificationReport=simplified.report;
        boundary=std::move(*simplified.boundary);
        std::cout<<"boundary_simplification=APPLIED\n"
                 <<"boundary_simplification_requested_deviation="
                 <<simplificationReport->requestedDeviation<<'\n'
                 <<"boundary_simplification_applied_deviation="
                 <<simplificationReport->appliedDeviation<<'\n'
                 <<"boundary_simplification_measured_max_deviation="
                 <<simplificationReport->measuredMaxDeviation<<'\n'
                 <<"boundary_original_vertices="
                 <<simplificationReport->originalVertexCount<<'\n'
                 <<"boundary_simplified_vertices="
                 <<simplificationReport->simplifiedVertexCount<<'\n'
                 <<"boundary_original_area="<<simplificationReport->originalArea<<'\n'
                 <<"boundary_simplified_area="<<simplificationReport->simplifiedArea<<'\n';
    }

    const auto refinementStart = std::chrono::steady_clock::now();
    Quadtree2D tree(domain, maxLevel, boundary);
    QuadtreeRefinementPolicy2D refinement;
    if (resolvedSizeField) {
        refinement = resolvedSizeField->refinement;
    } else {
        refinement.minimumLevel = minimumLevel;
        refinement.boundaryLevel = maxLevel;
        refinement.distanceBands=distanceBands;
        refinement.boxRegions=boxRegions;
    }
    try {
        tree.refine(boundary, refinement);
    } catch (const std::invalid_argument& exception) {
        std::cerr<<"invalid sizing field: "<<exception.what()<<'\n';
        return EXIT_FAILURE;
    }
    const double refinementSeconds = elapsedSeconds(refinementStart);
    const auto balanceStart = std::chrono::steady_clock::now();
    const auto balance = tree.enforceTwoToOneBalance(boundary);
    const double balanceSeconds = elapsedSeconds(balanceStart);
    if (balance.violationsAfter != 0 || !tree.deterministicOrderingValid()) {
        std::cerr << "Quadtree balance/determinism gate failed\n";
        return EXIT_FAILURE;
    }

    const auto cutCellStart = std::chrono::steady_clock::now();
    std::vector<CutCell2D> cutCells;
    cutCells.reserve(tree.leaves().size());
    std::size_t unsupported = 0;
    std::size_t splitFluidLeaves = 0;
    std::size_t nextSourceId = 0;
    double sourceFluidArea = 0.0;

    // One registry owns every wall/grid intersection across all leaves. Input
    // wall samples are immutable anchors; only constructed intersections may
    // consume the explicitly bounded grid-corner weld allowance.
    IntersectionRegistryPolicy2D registryPolicy;
    registryPolicy.gridCornerWeldFractionOfLocalH = gridCornerWeldFraction;
    auto constructionRegistry =
        std::make_shared<IntersectionRegistry2D>(registryPolicy);
    constructionRegistry->configureGrid(domain.bounds, maxLevel);
    const double finestCellH = std::ldexp(
        std::min(domain.width(), domain.height()), -static_cast<int>(maxLevel));
    for (const auto& loop : boundary.loops()) {
        for (const auto& vertex : loop.vertices()) {
            (void)constructionRegistry->internVertex(
                vertex, finestCellH, IntersectionFeature2D::Smooth);
        }
    }

    std::vector<ConstructionRecoveryRequest2D> recoveryRequests;
    for (const auto& leaf : tree.leaves()) {
        auto components = buildCutCellsShared(
            leaf, boundary, *constructionRegistry,
            IntersectionSource2D::WallCartesian, fluidRegion);
        if (components.size() > 1) ++splitFluidLeaves;
        for (auto& cut : components) {
            recoveryRequests.insert(recoveryRequests.end(),
                                    cut.constructionRecoveryRequests.begin(),
                                    cut.constructionRecoveryRequests.end());
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
    if (!recoveryRequests.empty()) {
        // The plain CLI has no local-refinement recovery loop. Never commit
        // geometry that the shared construction layer has flagged as unsafe.
        std::cerr << "shared construction requested " << recoveryRequests.size()
                  << " local recovery pass(es); the plain CLI has no recovery loop\n";
        const std::size_t reportCount =
            std::min<std::size_t>(recoveryRequests.size(), 5);
        for (std::size_t i = 0; i < reportCount; ++i) {
            const auto& request = recoveryRequests[i];
            std::cerr << "construction_recovery[" << i << "] support="
                      << request.supportId << " local_h=" << request.localH
                      << " at=(" << request.originalPoint.x << ','
                      << request.originalPoint.y << ") reason=" << request.reason << '\n';
        }
        return EXIT_FAILURE;
    }
    if (unsupported != 0) {
        std::cerr << "Cut-cell construction produced " << unsupported
                  << " unsupported leaf geometry case(s)\n";
        // Name the reasons instead of only counting them: an unnamed count is
        // exactly the diagnostic gap R2/W0 set out to close.
        std::size_t reported = 0;
        for (const auto& cut : cutCells) {
            if (cut.valid() || cut.kind != CutCellKind::Unsupported) continue;
            if (reported >= 5) break;
            std::cerr << "unsupported_leaf[" << reported << "] source=" << cut.sourceId
                      << " key=" << cut.sourceKey
                      << " bounds=(" << cut.backgroundBounds.min.x << ','
                      << cut.backgroundBounds.min.y << ")-("
                      << cut.backgroundBounds.max.x << ','
                      << cut.backgroundBounds.max.y << ')';
            for (const auto& issue : cut.issues) {
                std::cerr << " issue=" << issue.message;
            }
            std::cerr << '\n';
            ++reported;
        }
        return EXIT_FAILURE;
    }
    const double cutCellSeconds = elapsedSeconds(cutCellStart);

    // Product-level physics gate. This catches the exact class of error where
    // a mathematically self-consistent mesh is generated on the wrong side of
    // the solid boundary.
    const double domainArea = domain.width() * domain.height();
    const double solidArea = boundary.area();
    const double expectedFluidArea = fluidRegion == FluidRegion2D::Exterior
        ? domainArea - solidArea
        : solidArea;
    const TolerancePolicy tol{};
    const double roundoffAreaEps =
        std::max(tol.absolute * tol.absolute,
                 tol.relative * std::max(1.0, std::abs(expectedFluidArea)));
    double wallLength = 0.0;
    for (const auto& loop : boundary.loops()) {
        const auto& vertices = loop.vertices();
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            const auto& a = vertices[i];
            const auto& b = vertices[(i + 1) % vertices.size()];
            wallLength += std::hypot(b.x - a.x, b.y - a.y);
        }
    }
    // Each weld moves a point by at most f*h. Its incident area perturbation is
    // bounded by 1.5*f*h^2, while a wall of length L has at most L/h relevant
    // grid crossings. Thus the global geometric allowance is 1.5*f*h*L.
    const double weldAreaEps =
        1.5 * gridCornerWeldFraction * finestCellH * wallLength;
    const double areaEps = roundoffAreaEps + weldAreaEps;
    if (std::abs(sourceFluidArea - expectedFluidArea) > areaEps) {
        std::cerr << "fluid-side physics gate failed: region=" << fluidRegionName(fluidRegion)
                  << " generated_area=" << std::setprecision(17) << sourceFluidArea
                  << " expected_area=" << expectedFluidArea
                  << " allowed_error=" << areaEps << '\n';
        return EXIT_FAILURE;
    }

    const auto sourceTopologyStart = std::chrono::steady_clock::now();
    const TopologyMesh2D sourceTopology =
        buildGlobalTopology(cutCells, domain, boundary, TolerancePolicy{},
                            constructionRegistry);
    const double sourceTopologySeconds = elapsedSeconds(sourceTopologyStart);
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

    const auto agglomerationStart = std::chrono::steady_clock::now();
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
    const double agglomerationSeconds = elapsedSeconds(agglomerationStart);

    // A rejected mesh is exactly the case where the time went somewhere unknown,
    // so the late failure paths report the phases they already completed instead
    // of exiting silently. solverPhaseSeconds covers the phase that was still
    // running when the gate rejected the mesh, so nothing stays implicit.
    const auto printFailureTiming = [&](double solverPhaseSeconds = 0.0) {
        std::cout << "timing_refinement_seconds=" << refinementSeconds
                  << " timing_balance_seconds=" << balanceSeconds
                  << " timing_cut_cell_seconds=" << cutCellSeconds
                  << " timing_source_topology_seconds=" << sourceTopologySeconds
                  << " timing_agglomeration_seconds=" << agglomerationSeconds
                  << " timing_solver_topology_seconds=" << solverPhaseSeconds
                  << " timing_total_seconds=" << elapsedSeconds(totalStart) << '\n';
    };

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
        printFailureTiming();
        return EXIT_FAILURE;
    }

    const auto solverTopologyStart = std::chrono::steady_clock::now();
    std::optional<SolverQualityReport2D> solverQuality;
    std::optional<SolverTopologyResult2D> solverTopology;
    if (openFoamCase) {
        solverTopology=buildSolverTopology2D(stabilized.topology,domain,boundary);
        if (!solverTopology->valid()) {
            std::cerr<<"solver topology partition failed";
            for (const auto& issue:solverTopology->issues) std::cerr<<": "<<issue;
            std::cerr<<'\n';
            printFailureTiming(elapsedSeconds(solverTopologyStart));
            return EXIT_FAILURE;
        }
        solverQuality=evaluateSolverQuality2D(solverTopology->topology);
        if (!solverQuality->valid()) {
            std::cerr<<"solver-quality gate failed with "<<solverQuality->issues.size()
                     <<" issue(s) after "
                     <<solverTopology->qualityAgglomeratedSourceCellCount
                     <<" quality-driven source agglomeration(s) and "
                     <<solverTopology->qualityRepartitionCount
                     <<" local repartition(s); OpenFOAM mesh was not written\n";
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
            printFailureTiming(elapsedSeconds(solverTopologyStart));
            return EXIT_FAILURE;
        }
    }
    const double solverTopologySeconds = elapsedSeconds(solverTopologyStart);

    const auto serializationStart = std::chrono::steady_clock::now();
    const auto parent = outputPrefix.parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    const std::filesystem::path vtkPath = outputPrefix.string() + ".vtk";
    const std::filesystem::path cm2dPath = outputPrefix.string() + ".cm2d";
    const std::filesystem::path qualityPath =
        outputPrefix.string() + ".construction-quality.json";
    const std::filesystem::path vizPath = outputPrefix.string() + ".viz.json";
    const std::filesystem::path sizingPath = outputPrefix.string() + ".sizing.json";

    error.clear();
    if (!writeSizingFieldReport(sizingPath,domain,refinement,tree,balance,error)) {
        std::cerr<<error<<'\n';
        return EXIT_FAILURE;
    }
    // Only written when the field was requested, so no existing output set changes.
    const std::filesystem::path sizeFieldPath = outputPrefix.string() + ".size-field.json";
    if (resolvedSizeField) {
        std::ofstream out(sizeFieldPath);
        out << resolvedSizeFieldToJson(*resolvedSizeField);
        if (!out.good()) {
            std::cerr << "failed while writing size-field JSON output\n";
            return EXIT_FAILURE;
        }
    }

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
                                               *openFoamCase,extrusionThickness,&error,{},
                                               embeddedPatches);
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
    const double serializationSeconds = elapsedSeconds(serializationStart);
    std::size_t geometricCutCellCount = 0;
    for (const auto& cut : cutCells) {
        if (cut.kind == CutCellKind::Cut) ++geometricCutCellCount;
    }
    const double totalSeconds = elapsedSeconds(totalStart);

    std::cout << "cartmesh2d end-to-end PASS\n"
              << "fluid_region=" << fluidRegionName(fluidRegion) << '\n'
              << "boundary_role="
              << (fluidRegion == FluidRegion2D::Exterior ? "solid_wall" : "fluid_envelope") << '\n'
              << "boundary_loops=" << boundary.loops().size() << '\n'
              << "minimum_level=" << minimumLevel << '\n'
              << "leaf_count=" << tree.leaves().size() << '\n'
              << "sizeof_quadtree_leaf_bytes=" << sizeof(QuadtreeLeaf2D) << '\n'
              << "split_fluid_leaves=" << splitFluidLeaves << '\n'
              << "cut_cells=" << geometricCutCellCount << '\n'
              << "source_cells=" << sourceTopology.cells.size() << '\n'
              << "source_fluid_area=" << sourceFluidArea << '\n'
              << "expected_fluid_area=" << expectedFluidArea << '\n'
              << "small_cells=" << smallReport.smallCellCount << '\n'
              << "stabilized_cells=" << stabilized.topology.cells.size() << '\n'
              << "vertices=" << stabilized.topology.vertices.size() << '\n'
              << "edges=" << stabilized.topology.edges.size() << '\n'
              << "min_cell_area=" << quality.minCellArea << '\n'
              << "min_edge_length=" << quality.minEdgeLength << '\n'
              << "construction_max_cell_edge_length_ratio="
              << quality.maxCellEdgeLengthRatio << '\n'
              << "construction_max_centroid_vertex_mean_offset_normalized="
              << quality.maxCentroidVertexMeanOffsetNormalized << '\n'
              << "vtk=" << vtkPath.string() << '\n'
              << "cm2d=" << cm2dPath.string() << '\n'
              << "quality_json=" << qualityPath.string() << '\n'
              << "visualization_json=" << vizPath.string() << '\n';
    std::cout << std::setprecision(9)
              << "timing_refinement_seconds=" << refinementSeconds << '\n'
              << "timing_balance_seconds=" << balanceSeconds << '\n'
              << "timing_neighbor_generation_seconds="
              << balance.faceNeighborSeconds << '\n'
              << "timing_neighbor_generation_calls="
              << balance.faceNeighborCalls << '\n'
              << "timing_balance_excluding_neighbors_seconds="
              << std::max(0.0, balanceSeconds - balance.faceNeighborSeconds) << '\n'
              << "timing_cut_cell_seconds=" << cutCellSeconds << '\n'
              << "timing_source_topology_seconds=" << sourceTopologySeconds << '\n'
              << "timing_agglomeration_seconds=" << agglomerationSeconds << '\n'
              << "timing_solver_topology_seconds=" << solverTopologySeconds << '\n'
              << "timing_serialization_export_seconds=" << serializationSeconds << '\n'
              << "timing_total_seconds=" << totalSeconds << '\n';
    std::cout << "sizing_distance_bands=" << refinement.distanceBands.size() << '\n'
              << "sizing_box_regions=" << refinement.boxRegions.size() << '\n'
              << "sizing_json=" << sizingPath.string() << '\n';
    if (resolvedSizeField) {
        std::cout << "size_field=RESOLVED\n"
                  << "size_field_body_span=" << resolvedSizeField->bodySpan << '\n'
                  << "size_field_domain_span=" << resolvedSizeField->domainSpan << '\n'
                  << "size_field_wall_cell_size=" << resolvedSizeField->wallCellSize << '\n'
                  << "size_field_wall_level=" << resolvedSizeField->wallLevel << '\n'
                  << "size_field_curvature_level=" << resolvedSizeField->curvatureLevel << '\n'
                  << "size_field_proximity_level=" << resolvedSizeField->proximityLevel << '\n'
                  << "size_field_max_level=" << resolvedSizeField->maxLevel << '\n'
                  << "size_field_level_cap_reached="
                  << (resolvedSizeField->levelCapReached ? "true" : "false") << '\n'
                  << "sizing_segment_bands=" << refinement.segmentBands.size() << '\n'
                  << "size_field_json=" << sizeFieldPath.string() << '\n';
    }
    if (openFoamReport) {
        std::cout<<"solver_quality=PASS\n"
                 <<"solver_max_nonorthogonality_deg="<<solverQuality->maxNonOrthogonalityDeg<<'\n'
                 <<"solver_max_internal_skewness="<<solverQuality->maxInternalSkewness<<'\n'
                 <<"solver_max_boundary_skewness="<<solverQuality->maxBoundarySkewness<<'\n'
                 <<"solver_max_cell_aspect="<<solverQuality->maxCellAspect<<'\n'
                 <<"solver_min_face_length="<<solverQuality->minFaceLength<<'\n'
                 <<"solver_min_face_weight="<<solverQuality->minFaceWeight<<'\n'
                 <<"solver_min_volume_ratio="<<solverQuality->minVolumeRatio<<'\n'
                 <<"solver_partitioned_input_cells="<<solverTopology->partitionedCellCount<<'\n'
                 <<"solver_quality_agglomerated_source_cells="
                 <<solverTopology->qualityAgglomeratedSourceCellCount<<'\n'
                 <<"solver_quality_local_repartitions="
                 <<solverTopology->qualityRepartitionCount<<'\n'
                 <<"solver_profile_initial_partition_seconds="
                 <<solverTopology->profile.initialPartitionSeconds<<'\n'
                 <<"solver_profile_build_global_topology_seconds="
                 <<solverTopology->profile.buildGlobalTopologySeconds<<'\n'
                 <<"solver_profile_full_quality_seconds="
                 <<solverTopology->profile.fullQualitySeconds<<'\n'
                 <<"solver_profile_candidate_generation_seconds="
                 <<solverTopology->profile.candidateGenerationSeconds<<'\n'
                 <<"solver_profile_candidate_polygon_work_seconds="
                 <<solverTopology->profile.candidatePolygonWorkSeconds<<'\n'
                 <<"solver_profile_candidate_global_rebuild_seconds="
                 <<solverTopology->profile.candidateGlobalRebuildSeconds<<'\n'
                 <<"solver_profile_candidate_quality_seconds="
                 <<solverTopology->profile.candidateQualitySeconds<<'\n'
                 <<"solver_profile_source_repair_seconds="
                 <<solverTopology->profile.sourceRepairSeconds<<'\n'
                 <<"solver_profile_final_repartition_seconds="
                 <<solverTopology->profile.finalRepartitionSeconds<<'\n'
                 <<"solver_profile_global_topology_rebuild_calls="
                 <<solverTopology->profile.globalTopologyRebuildCalls<<'\n'
                 <<"solver_profile_full_quality_calls="
                 <<solverTopology->profile.fullQualityCalls<<'\n'
                 <<"solver_profile_candidate_topology_count="
                 <<solverTopology->profile.candidateTopologyCount<<'\n'
                 <<"solver_profile_candidate_quality_evaluation_count="
                 <<solverTopology->profile.candidateQualityEvaluationCount<<'\n'
                 <<"solver_profile_source_candidate_pairs="
                 <<solverTopology->profile.sourceCandidatePairs<<'\n'
                 <<"solver_profile_repartition_candidate_pairs="
                 <<solverTopology->profile.repartitionCandidatePairs<<'\n'
                 <<"solver_profile_candidate_splits="
                 <<solverTopology->profile.candidateSplits<<'\n'
                 <<"solver_profile_repair_patch_count="
                 <<solverTopology->profile.repairPatchCount<<'\n'
                 <<"solver_profile_accepted_topology_commit_count="
                 <<solverTopology->profile.acceptedTopologyCommitCount<<'\n'
                 <<"solver_profile_rejected_candidate_count="
                 <<(solverTopology->profile.candidateTopologyCount-
                    solverTopology->profile.acceptedTopologyCommitCount)<<'\n'
                 <<"solver_profile_maximum_quality_issue_count="
                 <<solverTopology->profile.maximumQualityIssueCount<<'\n'
                 <<"solver_profile_accepted_source_repairs="
                 <<solverTopology->profile.acceptedSourceRepairs<<'\n'
                 <<"solver_profile_accepted_repartitions="
                 <<solverTopology->profile.acceptedRepartitions<<'\n'
                 <<"solver_profile_source_repair_iterations="
                 <<solverTopology->profile.sourceRepairIterations<<'\n'
                 <<"solver_profile_repartition_iterations="
                 <<solverTopology->profile.repartitionIterations<<'\n'
                 <<"solver_output_cells="<<solverTopology->outputCellCount<<'\n'
                 <<"solver_output_vertices="<<solverTopology->topology.vertices.size()<<'\n'
                 <<"solver_output_faces="<<solverTopology->topology.edges.size()<<'\n'
                 <<"openfoam_case="<<openFoamCase->string()<<'\n'
                 <<"openfoam_cells="<<openFoamReport->cellCount<<'\n'
                 <<"openfoam_faces="<<openFoamReport->faceCount<<'\n';
    }
    return EXIT_SUCCESS;
}
