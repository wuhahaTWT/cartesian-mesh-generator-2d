#include "cartmesh2d/topology/Topology2D.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <set>
#include <string>
#include <tuple>
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

CutCell2D fullCell(std::size_t sourceId, std::uint64_t key, const AABB2D& box) {
    CutCell2D cell;
    cell.sourceId = sourceId;
    cell.sourceKey = key;
    cell.backgroundBounds = box;
    cell.kind = CutCellKind::Full;
    cell.fluidPolygon = {{{box.min.x, box.min.y},
                          {box.max.x, box.min.y},
                          {box.max.x, box.max.y},
                          {box.min.x, box.max.y}}};
    cell.area = (box.max.x - box.min.x) * (box.max.y - box.min.y);
    cell.areaFraction = 1.0;
    cell.centroid = cell.fluidPolygon.centroid();
    return cell;
}

bool hasPoint(const TopologyMesh2D& mesh, double x, double y, double eps = 1.0e-12) {
    for (const auto& v : mesh.vertices) {
        if (std::abs(v.point.x - x) <= eps && std::abs(v.point.y - y) <= eps) return true;
    }
    return false;
}

} // namespace

int main() {
    const Domain2D domain{{{0.0, 0.0}, {2.0, 1.0}}};
    BoundaryLoop boundary({{0.0, 0.0}, {2.0, 0.0}, {2.0, 1.0}, {0.0, 1.0}});

    std::vector<CutCell2D> cells;
    cells.push_back(fullCell(10, 10, {{0.0, 0.0}, {1.0, 1.0}}));
    cells.push_back(fullCell(20, 20, {{1.0, 0.0}, {2.0, 0.5}}));
    cells.push_back(fullCell(30, 30, {{1.0, 0.5}, {2.0, 1.0}}));

    const auto mesh = buildGlobalTopology(cells, domain, boundary);
    check(mesh.valid(), "coarse-fine topology passes audit");
    check(mesh.cells.size() == 3, "three fluid cells retained");
    check(hasPoint(mesh, 1.0, 0.5), "hanging-node point exists globally");

    std::size_t internalCount = 0;
    std::size_t boundaryCount = 0;
    std::size_t coarseFineInterfaceFragments = 0;
    bool sawLongUnsplitInterface = false;
    for (const auto& edge : mesh.edges) {
        const auto& a = mesh.vertices[edge.v0].point;
        const auto& b = mesh.vertices[edge.v1].point;
        if (edge.neighbour) {
            ++internalCount;
            check(edge.patch == BoundaryPatch2D::None, "internal edge has no boundary patch");
            check(edge.owner != *edge.neighbour, "internal edge has distinct owner/neighbour");
        } else {
            ++boundaryCount;
            check(edge.patch != BoundaryPatch2D::None &&
                  edge.patch != BoundaryPatch2D::Unclassified,
                  "boundary edge has a classified patch");
        }

        if (std::abs(a.x - 1.0) < 1.0e-12 && std::abs(b.x - 1.0) < 1.0e-12) {
            const double length = std::abs(b.y - a.y);
            if (std::abs(length - 0.5) < 1.0e-12 && edge.neighbour) {
                ++coarseFineInterfaceFragments;
            }
            if (length > 0.5 + 1.0e-12) sawLongUnsplitInterface = true;
        }
    }
    check(internalCount == 3, "coarse-fine fixture has three internal edge fragments");
    check(coarseFineInterfaceFragments == 2,
          "coarse face is split into two owner-neighbour fragments");
    check(!sawLongUnsplitInterface, "no T-junction-spanning coarse edge remains");
    check(boundaryCount == 7, "outer rectangle is represented by seven boundary fragments");

    for (const auto& cell : mesh.cells) {
        check(cell.vertices.size() == cell.edges.size(), "cell edge loop closes");
        check(cell.vertices.size() >= 4, "each fixture cell has a valid split polygon loop");
    }

    std::set<std::tuple<std::size_t, std::size_t>> uniqueEdges;
    for (const auto& edge : mesh.edges) {
        const auto key = std::minmax(edge.v0, edge.v1);
        check(uniqueEdges.insert({key.first, key.second}).second,
              "no duplicate canonical global edge");
    }
    for (std::size_t i = 0; i < mesh.vertices.size(); ++i) {
        for (std::size_t j = i + 1; j < mesh.vertices.size(); ++j) {
            const auto& a = mesh.vertices[i].point;
            const auto& b = mesh.vertices[j].point;
            check(std::hypot(a.x - b.x, a.y - b.y) > 1.0e-12,
                  "no duplicate global vertex in fixture");
        }
    }

    const double area = 1.0 + 0.5 + 0.5;
    double topologyArea = 0.0;
    for (const auto& cell : mesh.cells) topologyArea += cell.geometryArea;
    check(std::abs(topologyArea - area) <= 1.0e-12,
          "global topology preserves source fluid-cell area");

    auto duplicate = cells;
    duplicate.push_back(cells.front());
    const auto duplicateResult = buildGlobalTopology(duplicate, domain, boundary);
    check(!duplicateResult.valid(), "duplicate source cell is rejected");

    // Regression retained from the first million-class level-17 run: a
    // domain-scaled canonicalization epsilon collapsed the two x coordinates
    // of a legitimate thin cell.  Its edge is above the geometry tolerance,
    // so topology must preserve it rather than silently merge its vertices.
    constexpr double thinWidth = 1.5e-10;
    const Domain2D thinDomain{{{0.0, 0.0}, {2.0, 1.0}}};
    BoundaryLoop thinBoundary({{0.0, 0.0}, {2.0, 0.0}, {2.0, 1.0}, {0.0, 1.0}});
    std::vector<CutCell2D> thinCells;
    thinCells.push_back(fullCell(40, 40, {{0.0, 0.0}, {thinWidth, 1.0}}));
    thinCells.push_back(fullCell(50, 50, {{thinWidth, 0.0}, {2.0, 1.0}}));
    const auto thinMesh = buildGlobalTopology(thinCells, thinDomain, thinBoundary);
    check(thinMesh.valid(), "locally resolved thin cell survives vertex canonicalization");
    check(thinMesh.cells.size() == 2 && hasPoint(thinMesh, thinWidth, 0.0, 1.0e-14),
          "thin-cell interface remains a distinct canonical coordinate");

    if (failures != 0) {
        std::cerr << failures << " test(s) failed\n";
        return EXIT_FAILURE;
    }
    std::cout << "cartmesh2d 2D-4 topology tests passed\n";
    return EXIT_SUCCESS;
}
