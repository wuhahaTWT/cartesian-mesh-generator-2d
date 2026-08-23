#include "cartmesh2d/topology/Topology2D.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <queue>
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

struct PipelineResult {
    Quadtree2D tree;
    std::vector<CutCell2D> cuts;
    TopologyMesh2D topology;
    std::size_t unsupported = 0;
    double area = 0.0;

    PipelineResult(Quadtree2D treeIn) : tree(std::move(treeIn)) {}
};

PipelineResult generate(const Domain2D& domain, const BoundaryRegion2D& boundary,
                        FluidRegion2D fluidRegion, std::size_t maxLevel) {
    PipelineResult result(Quadtree2D(domain,maxLevel,boundary));
    QuadtreeRefinementPolicy2D policy;
    policy.boundaryLevel=maxLevel;
    result.tree.refine(boundary,policy);
    const auto balance=result.tree.enforceTwoToOneBalance(boundary);
    if (balance.violationsAfter!=0) ++failures;

    std::size_t sourceId=0;
    for (const auto& leaf:result.tree.leaves()) {
        auto components=buildCutCells(leaf,boundary,fluidRegion);
        for (auto& cut:components) {
            cut.sourceId=sourceId++;
            cut.sourceKey=leaf.key;
            if (cut.kind==CutCellKind::Unsupported) ++result.unsupported;
            if (cut.kind!=CutCellKind::Empty && cut.kind!=CutCellKind::Unsupported) {
                result.area+=cut.area;
            }
            result.cuts.push_back(std::move(cut));
        }
    }
    result.topology=buildGlobalTopology(result.cuts,domain,boundary);
    return result;
}

std::size_t embeddedEdges(const TopologyMesh2D& topology) {
    std::size_t count=0;
    for (const auto& edge:topology.edges) {
        if (edge.patch==BoundaryPatch2D::EmbeddedBoundary) ++count;
    }
    return count;
}

std::size_t cellComponents(const TopologyMesh2D& topology) {
    std::vector<std::vector<std::size_t>> adjacency(topology.cells.size());
    for (const auto& edge:topology.edges) {
        if (!edge.neighbour) continue;
        adjacency[edge.owner].push_back(*edge.neighbour);
        adjacency[*edge.neighbour].push_back(edge.owner);
    }
    std::vector<bool> visited(topology.cells.size(),false);
    std::size_t components=0;
    for (std::size_t start=0;start<topology.cells.size();++start) {
        if (visited[start]) continue;
        ++components;
        std::queue<std::size_t> pending;
        pending.push(start);
        visited[start]=true;
        while (!pending.empty()) {
            const auto cell=pending.front();
            pending.pop();
            for (const auto neighbour:adjacency[cell]) {
                if (!visited[neighbour]) {
                    visited[neighbour]=true;
                    pending.push(neighbour);
                }
            }
        }
    }
    return components;
}
} // namespace

int main() {
    BoundaryRegion2D obstacles({
        BoundaryLoop({{0.15,0.10},{1.15,0.10},{1.15,1.10},{0.15,1.10}}),
        BoundaryLoop({{2.25,-0.10},{3.25,-0.10},{3.25,0.90},{2.25,0.90}})
    });
    const Domain2D exteriorDomain{{{-1.0,-1.0},{4.0,2.0}}};
    auto exterior=generate(exteriorDomain,obstacles,FluidRegion2D::Exterior,6);
    check(exterior.unsupported==0,"two-obstacle exterior has no unsupported Cut-cell");
    check(exterior.topology.valid(),"two-obstacle exterior topology passes audit");
    check(std::abs(exterior.area-13.0)<=1.0e-9,
          "two-obstacle exterior conserves domain-minus-solid area");
    check(embeddedEdges(exterior.topology)>0,
          "two-obstacle exterior retains embedded wall edges");

    BoundaryRegion2D annulus({
        BoundaryLoop({{0.0,0.0},{4.0,0.0},{4.0,4.0},{0.0,4.0}}),
        BoundaryLoop({{1.1,1.1},{2.9,1.1},{2.9,2.9},{1.1,2.9}})
    });
    const Domain2D annulusDomain{{{-1.0,-1.0},{5.0,5.0}}};
    auto interior=generate(annulusDomain,annulus,FluidRegion2D::Interior,6);
    check(interior.unsupported==0,"annular interior has no unsupported Cut-cell");
    check(interior.topology.valid(),"annular interior topology passes audit");
    check(std::abs(interior.area-12.76)<=1.0e-9,
          "annular interior conserves outer-minus-hole area");
    check(embeddedEdges(interior.topology)>0,
          "annular interior retains both physical boundary loops");
    check(cellComponents(interior.topology)==1,
          "annular interior is one connected fluid component with one hole");
    auto exteriorCavity=generate(annulusDomain,annulus,FluidRegion2D::Exterior,6);
    check(exteriorCavity.unsupported==0 && exteriorCavity.topology.valid(),
          "annular complement passes Cut-cell and topology gates");
    check(std::abs(exteriorCavity.area-23.24)<=1.0e-9,
          "annular complement conserves domain-minus-parity-region area");
    check(cellComponents(exteriorCavity.topology)==2,
          "annular complement retains outside flow and enclosed cavity as separate components");
    const auto coarseHole=buildCutCells({{-0.5,-0.5},{4.5,4.5}},
                                        CellClass::Intersected,annulus,
                                        FluidRegion2D::Interior);
    if (!(coarseHole.size()==1 && coarseHole.front().kind==CutCellKind::Unsupported)) {
        std::cerr<<"coarse-hole diagnostic components="<<coarseHole.size();
        for (const auto& cell:coarseHole) {
            std::cerr<<" kind="<<static_cast<int>(cell.kind)
                     <<" area="<<cell.area<<" issues="<<cell.issues.size();
        }
        std::cerr<<'\n';
    }
    check(coarseHole.size()==1 && coarseHole.front().kind==CutCellKind::Unsupported &&
          !coarseHole.front().issues.empty() &&
          coarseHole.front().issues.front().code==CutCellIssueCode::MultipleEmbeddedComponents,
          "a hole wholly contained in one coarse leaf is rejected, not silently filled");

    BoundaryRegion2D nestedIsland({
        BoundaryLoop({{0.0,0.0},{5.0,0.0},{5.0,5.0},{0.0,5.0}}),
        BoundaryLoop({{1.0,1.0},{4.0,1.0},{4.0,4.0},{1.0,4.0}}),
        BoundaryLoop({{2.0,2.0},{3.0,2.0},{3.0,3.0},{2.0,3.0}})
    });
    const Domain2D nestedDomain{{{-1.0,-1.0},{6.0,6.0}}};
    auto disconnected=generate(nestedDomain,nestedIsland,FluidRegion2D::Interior,6);
    check(disconnected.unsupported==0 && disconnected.topology.valid(),
          "three-depth nested island passes Cut-cell and topology gates");
    check(std::abs(disconnected.area-17.0)<=1.0e-9,
          "three-depth nested island conserves parity area");
    check(cellComponents(disconnected.topology)==2,
          "depth-two island remains a separate fluid component instead of being bridged");

    BoundaryRegion2D invalid({
        BoundaryLoop({{0.0,0.0},{2.0,0.0},{2.0,2.0},{0.0,2.0}}),
        BoundaryLoop({{1.0,-0.5},{1.5,-0.5},{1.5,2.5},{1.0,2.5}})
    });
    check(!invalid.diagnose().valid(),"intersecting loops are rejected before meshing");

    if (failures!=0) {
        std::cerr<<failures<<" multi-loop test(s) failed\n";
        return EXIT_FAILURE;
    }
    std::cout<<"cartmesh2d D1 multi-loop tests passed\n";
    return EXIT_SUCCESS;
}
