#include "cartmesh2d/topology/Topology2D.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numbers>
#include <set>
#include <string>
#include <utility>
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
} // namespace

int main() {
    std::vector<Point2D> circle;
    constexpr std::size_t segments = 64;
    for (std::size_t i = 0; i < segments; ++i) {
        const double angle = 2.0 * std::numbers::pi * static_cast<double>(i) /
                             static_cast<double>(segments);
        circle.push_back({std::cos(angle), std::sin(angle)});
    }
    BoundaryLoop boundary(circle);
    const Domain2D domain{{{-2.0, -2.0}, {2.0, 2.0}}};

    Quadtree2D tree(domain, 4, boundary);
    QuadtreeRefinementPolicy2D policy;
    policy.boundaryLevel = 4;
    tree.refine(boundary, policy);
    const auto balance = tree.enforceTwoToOneBalance(boundary);
    check(balance.violationsAfter == 0, "adaptive fixture is 2:1 balanced");

    std::vector<CutCell2D> cutCells;
    std::size_t nonEmpty = 0;
    std::size_t unsupported = 0;
    for (const auto& leaf : tree.leaves()) {
        auto cut = buildCutCell(leaf, boundary);
        if (cut.kind == CutCellKind::Unsupported) ++unsupported;
        if (cut.kind != CutCellKind::Empty && cut.kind != CutCellKind::Unsupported) ++nonEmpty;
        cutCells.push_back(std::move(cut));
    }
    check(unsupported == 0, "circle fixture has no unsupported Cut-cell");

    const auto mesh = buildGlobalTopology(cutCells, domain, boundary);
    check(mesh.valid(), "adaptive circle global topology passes audit");
    check(mesh.cells.size() == nonEmpty, "all non-empty Cut-cells enter topology");
    check(mesh.audit.duplicateVertices == 0, "duplicate vertex audit = 0");
    check(mesh.audit.duplicateEdges == 0, "duplicate edge audit = 0");
    check(mesh.audit.orphanInternalEdges == 0, "orphan internal edge audit = 0");
    check(mesh.audit.nonManifoldEdges == 0, "non-manifold edge audit = 0");
    check(mesh.audit.unclassifiedBoundaryEdges == 0, "all boundary edges classified");
    check(mesh.audit.openCellLoops == 0, "all cell loops close");
    check(mesh.audit.areaMismatches == 0, "topology/source cell areas agree");

    std::size_t independentDuplicateVertices = 0;
    for (std::size_t i = 0; i < mesh.vertices.size(); ++i) {
        for (std::size_t j = i + 1; j < mesh.vertices.size(); ++j) {
            const auto& a = mesh.vertices[i].point;
            const auto& b = mesh.vertices[j].point;
            if (std::hypot(a.x - b.x, a.y - b.y) <= 1.0e-12) {
                ++independentDuplicateVertices;
            }
        }
    }
    check(independentDuplicateVertices == 0,
          "independent duplicate-vertex scan = 0");

    std::set<std::pair<std::size_t, std::size_t>> canonicalEdges;
    std::size_t independentDuplicateEdges = 0;
    std::size_t badOwnerNeighbour = 0;
    std::size_t badCellLoops = 0;
    for (const auto& edge : mesh.edges) {
        const auto canonical = std::minmax(edge.v0, edge.v1);
        if (!canonicalEdges.insert(canonical).second) ++independentDuplicateEdges;
        if (edge.neighbour && edge.owner == *edge.neighbour) ++badOwnerNeighbour;
    }
    for (const auto& cell : mesh.cells) {
        if (cell.vertices.size() < 3 || cell.vertices.size() != cell.edges.size()) {
            ++badCellLoops;
        }
    }
    check(independentDuplicateEdges == 0,
          "independent duplicate-edge scan = 0");
    check(badOwnerNeighbour == 0,
          "independent owner-neighbour audit = 0");
    check(badCellLoops == 0,
          "independent cell-loop audit = 0");

    std::size_t internalEdges = 0;
    std::size_t embeddedEdges = 0;
    std::size_t domainEdges = 0;
    for (const auto& edge : mesh.edges) {
        if (edge.neighbour) {
            ++internalEdges;
            check(edge.patch == BoundaryPatch2D::None,
                  "internal adaptive edge carries no patch");
        } else if (edge.patch == BoundaryPatch2D::EmbeddedBoundary) {
            ++embeddedEdges;
        } else if (edge.patch == BoundaryPatch2D::DomainBoundary) {
            ++domainEdges;
        }
    }
    check(internalEdges > 0, "adaptive mesh has internal owner-neighbour edges");
    check(embeddedEdges > 0, "circle produces embedded-boundary edges");
    check(domainEdges == 0, "interior circle fluid does not touch outer domain boundary");

    double topologyArea = 0.0;
    for (const auto& cell : mesh.cells) topologyArea += cell.geometryArea;
    const double inputArea = boundary.polygon().area();
    check(std::abs(topologyArea - inputArea) <= 1.0e-10,
          "global topology preserves total Cut-cell fluid area");

    const auto mesh2 = buildGlobalTopology(cutCells, domain, boundary);
    check(mesh2.vertices.size() == mesh.vertices.size() &&
          mesh2.edges.size() == mesh.edges.size() &&
          mesh2.cells.size() == mesh.cells.size(),
          "repeated topology build has stable counts");
    if (mesh2.vertices.size() == mesh.vertices.size()) {
        for (std::size_t i = 0; i < mesh.vertices.size(); ++i) {
            check(mesh2.vertices[i].id == mesh.vertices[i].id &&
                  std::abs(mesh2.vertices[i].point.x - mesh.vertices[i].point.x) <= 1.0e-14 &&
                  std::abs(mesh2.vertices[i].point.y - mesh.vertices[i].point.y) <= 1.0e-14,
                  "repeated topology vertex order is deterministic");
        }
    }
    if (mesh2.edges.size() == mesh.edges.size()) {
        for (std::size_t i = 0; i < mesh.edges.size(); ++i) {
            const auto& a = mesh.edges[i];
            const auto& b = mesh2.edges[i];
            check(a.id == b.id && a.v0 == b.v0 && a.v1 == b.v1 &&
                  a.owner == b.owner && a.neighbour == b.neighbour && a.patch == b.patch,
                  "repeated topology edge order is deterministic");
        }
    }

    if (failures != 0) {
        std::cerr << failures << " test(s) failed\n";
        return EXIT_FAILURE;
    }
    std::cout << "cartmesh2d 2D-4 adaptive topology audit passed\n";
    return EXIT_SUCCESS;
}
