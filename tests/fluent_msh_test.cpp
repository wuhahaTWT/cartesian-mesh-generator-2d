#include "cartmesh2d/io/FluentMsh2D.hpp"
#include "cartmesh2d/cutcell/CutCell2D.hpp"
#include "cartmesh2d/hybrid/HybridMesh2D.hpp"
#include "cartmesh2d/quadtree/Quadtree2D.hpp"
#include "cartmesh2d/quality/SolverTopology2D.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <numbers>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace cartmesh2d;

namespace {

int failures = 0;

void check(bool pass, const char* message) {
    if (!pass) {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}

struct ParsedFluentMesh {
    int dimension = 0;
    std::vector<Point2D> nodes;
    std::vector<int> cellTypes;
    struct Face {
        std::size_t n0 = 0, n1 = 0, c0 = 0, c1 = 0;
        std::size_t zone = 0;
        std::size_t bcCode = 0;
    };
    std::vector<Face> faces;
    std::map<std::size_t, std::pair<std::string, std::string>> zoneNames;
};

ParsedFluentMesh parseFluentMsh(const std::filesystem::path& path) {
    ParsedFluentMesh m;
    std::ifstream in(path);
    std::string line;
    enum class Section { None, Nodes, Cells, Faces };
    Section sec = Section::None;
    std::size_t curZone = 0, curBc = 0;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line == "))") { sec = Section::None; continue; }
        if (line == "(2 2)") { m.dimension = 2; continue; }
        if (line.rfind("(10 (1 ", 0) == 0) { sec = Section::Nodes; continue; }
        if (line.rfind("(12 (2 ", 0) == 0) { sec = Section::Cells; continue; }
        if (line.rfind("(13 (", 0) == 0 && line.rfind("(13 (0 ", 0) != 0) {
            std::string body = line.substr(5);
            for (char& ch : body) if (ch == '(' || ch == ')') ch = ' ';
            std::istringstream ss(body);
            std::string zHex, sHex, eHex, bcHex, fType;
            ss >> zHex >> sHex >> eHex >> bcHex >> fType;
            curZone = std::stoull(zHex, nullptr, 16);
            curBc = std::stoull(bcHex, nullptr, 16);
            sec = Section::Faces;
            continue;
        }
        if (line.rfind("(45 (", 0) == 0) {
            std::string body = line.substr(5);
            for (char& ch : body) if (ch == '(' || ch == ')') ch = ' ';
            std::istringstream ss(body);
            std::size_t zid = 0;
            std::string bc, name;
            ss >> zid >> bc >> name;
            m.zoneNames[zid] = {bc, name};
            continue;
        }
        if (sec == Section::Nodes) {
            std::istringstream ss(line);
            Point2D p;
            if (ss >> p.x >> p.y) m.nodes.push_back(p);
        } else if (sec == Section::Cells) {
            std::istringstream ss(line);
            int t = 0;
            while (ss >> t) m.cellTypes.push_back(t);
        } else if (sec == Section::Faces) {
            std::istringstream ss(line);
            std::string a, b, c0, c1;
            if (ss >> a >> b >> c0 >> c1) {
                m.faces.push_back({
                    std::stoull(a, nullptr, 16),
                    std::stoull(b, nullptr, 16),
                    std::stoull(c0, nullptr, 16),
                    std::stoull(c1, nullptr, 16),
                    curZone,
                    curBc
                });
            }
        }
    }
    return m;
}

void verifyFluentReadback(const TopologyMesh2D& topo,
                          const ParsedFluentMesh& parsed,
                          const FluentMshWriteReport2D& report,
                          const char* label) {
    check(parsed.dimension == 2, label);
    check(parsed.nodes.size() == topo.vertices.size() && report.nodeCount == topo.vertices.size(),
          "Fluent node count must match topology vertices");
    check(parsed.cellTypes.size() == topo.cells.size() && report.cellCount == topo.cells.size(),
          "Fluent cell count must match topology cells");
    check(parsed.faces.size() == topo.edges.size() && report.faceCount == topo.edges.size(),
          "Fluent face count must match topology edges");

    // Exact coordinate preservation (zero wall shift).
    double maxCoordDiff = 0.0;
    for (std::size_t i = 0; i < topo.vertices.size(); ++i) {
        maxCoordDiff = std::max(maxCoordDiff,
            std::hypot(parsed.nodes[i].x - topo.vertices[i].point.x,
                       parsed.nodes[i].y - topo.vertices[i].point.y));
    }
    check(maxCoordDiff <= 1e-14, "Fluent node coordinates must preserve exact input vertices");

    // Independent Green's theorem area & normal-closure check on EVERY cell using
    // Fluent's 2D left-hand c0 convention: outward normal of c0 is (dy, -dx).
    std::vector<Vector2D> closure(topo.cells.size(), {0.0, 0.0});
    std::vector<double> greenArea(topo.cells.size(), 0.0);
    std::vector<std::size_t> faceDegree(topo.cells.size(), 0U);
    for (const auto& f : parsed.faces) {
        check(f.n0 >= 1U && f.n0 <= parsed.nodes.size() &&
              f.n1 >= 1U && f.n1 <= parsed.nodes.size() && f.n0 != f.n1,
              "Fluent face node indices must be valid 1-based distinct IDs");
        check(f.c0 >= 1U && f.c0 <= topo.cells.size(),
              "Fluent face c0 must be valid 1-based cell ID");
        check(f.c1 <= topo.cells.size() && f.c1 != f.c0,
              "Fluent face c1 must be 0 or valid distinct cell ID");
        const auto p0 = parsed.nodes[f.n0 - 1U];
        const auto p1 = parsed.nodes[f.n1 - 1U];
        const Vector2D outward{p1.y - p0.y, -(p1.x - p0.x)};
        const double crossContrib = 0.5 * (p0.x * p1.y - p1.x * p0.y);
        closure[f.c0 - 1U].x += outward.x;
        closure[f.c0 - 1U].y += outward.y;
        greenArea[f.c0 - 1U] += crossContrib;
        ++faceDegree[f.c0 - 1U];
        if (f.c1 > 0U) {
            closure[f.c1 - 1U].x -= outward.x;
            closure[f.c1 - 1U].y -= outward.y;
            greenArea[f.c1 - 1U] -= crossContrib;
            ++faceDegree[f.c1 - 1U];
        }
    }
    double totalArea = 0.0, expectedArea = 0.0;
    for (std::size_t c = 0; c < topo.cells.size(); ++c) {
        check(std::hypot(closure[c].x, closure[c].y) <= 1e-12,
              "Every Fluent cell outward normal sum must close to zero");
        check(greenArea[c] > 0.0,
              "Every Fluent cell must have strictly positive area under left-hand c0 convention");
        const int expectedType = faceDegree[c] == 3U ? 1 : (faceDegree[c] == 4U ? 3 : 7);
        check(parsed.cellTypes[c] == expectedType,
              "Fluent cell type code (1/3/7) must match polygon vertex count");
        totalArea += greenArea[c];
        expectedArea += topo.cells[c].geometryArea;
    }
    check(std::abs(totalArea - expectedArea) <= 1e-11 * std::max(1.0, expectedArea),
          "Total Fluent mesh area must match topology fluid area");
    check(parsed.zoneNames.contains(2U) && parsed.zoneNames.at(2U).first == "fluid",
          "Fluent (45 ...) must define fluid cell zone 2");
    for (const auto& z : report.zones) {
        check(parsed.zoneNames.contains(z.zoneId) && parsed.zoneNames.at(z.zoneId).second == z.name,
              "Every Fluent face zone must have a matching decimal (45 ...) entry");
    }
}

std::vector<Point2D> makeCircle(double radius, std::size_t n) {
    std::vector<Point2D> pts;
    pts.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double theta = 2.0 * std::numbers::pi * static_cast<double>(i) / static_cast<double>(n);
        pts.push_back({radius * std::cos(theta), radius * std::sin(theta)});
    }
    return pts;
}

std::vector<CutCell2D> collectCutCells(const Quadtree2D& tree,
                                       const BoundaryRegion2D& boundary,
                                       FluidRegion2D region) {
    std::vector<CutCell2D> out;
    std::size_t nextId = 0U;
    for (const auto& leaf : tree.leaves()) {
        for (auto cell : buildCutCells(leaf, boundary, region)) {
            if (cell.kind == CutCellKind::Empty) continue;
            cell.sourceId = nextId++;
            cell.sourceKey = leaf.key;
            cell.sourceLineage = {cell.sourceId};
            out.push_back(std::move(cell));
        }
    }
    return out;
}

} // namespace

int main() {
    const auto tmpDir = std::filesystem::temp_directory_path() / "cartmesh2d_fluent_msh_tests";
    std::filesystem::create_directories(tmpDir);

    // 1. Exterior circle Cut-cell + shape agglomeration (mixed tri/quad/poly + named domain patches).
    {
        const BoundaryRegion2D boundary{BoundaryLoop(makeCircle(0.5, 48U))};
        const Domain2D domain{{{-2.0, -2.0}, {2.0, 2.0}}};
        Quadtree2D tree(domain, 5U, boundary);
        QuadtreeRefinementPolicy2D ref;
        ref.minimumLevel = 3U;
        ref.boundaryLevel = 5U;
        tree.refine(boundary, ref);
        (void)tree.enforceTwoToOneBalance(boundary);
        const auto cutCells = collectCutCells(tree, boundary, FluidRegion2D::Exterior);
        const auto topo = buildGlobalTopology(cutCells, domain, boundary);
        const auto shaped = agglomerateSolverByShape2D(topo, domain, boundary);
        check(shaped.valid() && shaped.topology.valid(), "circle shape mesh valid");

        const auto msh1 = tmpDir / "circle_1.msh";
        const auto msh2 = tmpDir / "circle_2.msh";
        std::string err;
        const auto rep1 = writeFluentMsh2D(shaped.topology, domain, boundary, msh1, &err);
        const auto rep2 = writeFluentMsh2D(shaped.topology, domain, boundary, msh2, &err);
        check(rep1.valid() && rep2.valid(), "circle Fluent .msh export must succeed");
        const auto parsed = parseFluentMsh(msh1);
        verifyFluentReadback(shaped.topology, parsed, rep1, "exterior circle");
        // Verify velocity-inlet hex BC code 0xa on "left".
        bool sawLeftHexA = false;
        for (const auto& f : parsed.faces) {
            if (parsed.zoneNames.contains(f.zone) && parsed.zoneNames.at(f.zone).second == "left") {
                check(f.bcCode == 10U, "left inlet patch must use Fluent hex bcType a (10)");
                sawLeftHexA = true;
            }
        }
        check(sawLeftHexA, "exterior circle must contain non-empty left inlet zone");

        // Determinism: byte-for-byte identical output.
        std::ifstream f1(msh1, std::ios::binary), f2(msh2, std::ios::binary);
        std::string s1((std::istreambuf_iterator<char>(f1)), {});
        std::string s2((std::istreambuf_iterator<char>(f2)), {});
        check(s1 == s2 && !s1.empty(), "Fluent .msh export must be byte-for-byte deterministic");
    }

    // 2. Interior fluid region (no domain patches, only wall_0 + interior).
    {
        const BoundaryRegion2D boundary{BoundaryLoop({{0,0},{2,0},{2,1},{0,1}})};
        const Domain2D domain{{{-0.5, -0.5}, {2.5, 1.5}}};
        Quadtree2D tree(domain, 4U, boundary);
        QuadtreeRefinementPolicy2D ref;
        ref.minimumLevel = 2U;
        ref.boundaryLevel = 4U;
        tree.refine(boundary, ref);
        (void)tree.enforceTwoToOneBalance(boundary);
        const auto cutCells = collectCutCells(tree, boundary, FluidRegion2D::Interior);
        const auto topo = buildGlobalTopology(cutCells, domain, boundary);
        const auto msh = tmpDir / "interior.msh";
        std::string err;
        const auto rep = writeFluentMsh2D(topo, domain, boundary, msh, &err);
        check(rep.valid(), "interior Fluent .msh export must succeed");
        const auto parsed = parseFluentMsh(msh);
        verifyFluentReadback(topo, parsed, rep, "interior rect");
        check(rep.zones.size() == 2U, "interior rectangle should emit only interior + wall_0 zones");
    }

    // 3. Conformal boundary-layer Hybrid mesh export.
    {
        const auto circlePts = makeCircle(0.5, 48U);
        const BoundaryRegion2D walls{BoundaryLoop(circlePts)};
        const Domain2D domain{{{-2.5, -2.5}, {2.5, 2.5}}};
        auto chain = makeClosedWallChain2D(walls.loops().front(), 0U, "wall_0", WallFluidRegion2D::Exterior);
        check(chain.success(), "hybrid circle chain valid");
        LayerParameters2D lp;
        lp.nLayers = 3U;
        lp.thicknessMode = LayerThicknessMode2D::FirstLayerThickness;
        lp.thickness = 0.02;
        lp.growthRatio = 1.2;
        QuadtreeRefinementPolicy2D ref;
        ref.minimumLevel = 3U;
        ref.boundaryLevel = 5U;
        HybridMeshPolicy2D hp;
        hp.fluidRegion = FluidRegion2D::Exterior;
        const auto robust = buildRobustH4Mesh2D({*chain.chain}, lp, domain, walls, 5U, ref, {}, hp);
        check(robust.success() && robust.hybridCandidate.solverTopology.valid(), "hybrid mesh valid");
        const auto msh = tmpDir / "hybrid_circle.msh";
        std::string err;
        const auto rep = writeFluentMsh2D(robust.hybridCandidate.solverTopology, domain, walls, msh, &err);
        check(rep.valid(), "hybrid Fluent .msh export must succeed");
        verifyFluentReadback(robust.hybridCandidate.solverTopology, parseFluentMsh(msh), rep, "hybrid circle");
    }

    // 4. Fail-closed invalid inputs.
    {
        TopologyMesh2D empty;
        const Domain2D domain{{{0, 0}, {1, 1}}};
        const BoundaryRegion2D wall{BoundaryLoop({{0.2, 0.2}, {0.8, 0.2}, {0.8, 0.8}, {0.2, 0.8}})};
        std::string err;
        const auto rep = writeFluentMsh2D(empty, domain, wall, tmpDir / "should_not_exist.msh", &err);
        check(!rep.valid() && !err.empty(), "empty topology must fail closed");
    }

    std::error_code ec;
    std::filesystem::remove_all(tmpDir, ec);
    return failures ? 1 : 0;
}
