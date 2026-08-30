#include "cartmesh2d/topology/EdgeIncidence2D.hpp"

#include <algorithm>
#include <map>
#include <set>
#include <sstream>

namespace cartmesh2d {

EdgeIncidenceStore2D buildEdgeIncidenceStore2D(
    const TopologyMesh2D& topology, std::uint64_t revision) {
    EdgeIncidenceStore2D result;
    result.revision=revision;
    if (!topology.valid()) {
        result.issues.push_back({EdgeIncidenceIssueCode2D::InvalidTopology,0,
                                 "edge-incidence input topology is invalid"});
        return result;
    }

    result.stableVertexIds.resize(topology.vertices.size());
    const bool hasConstructionHandles=
        topology.canonicalVertexIds.size()==topology.vertices.size();
    for (std::size_t vertex=0;vertex<topology.vertices.size();++vertex) {
        result.stableVertexIds[vertex]=hasConstructionHandles
            ?static_cast<StableVertexId2D>(topology.canonicalVertexIds[vertex])
            :static_cast<StableVertexId2D>(topology.vertices[vertex].id);
    }
    result.stableEdgeKeys.resize(topology.edges.size());
    result.edgeHalfEdges.resize(topology.edges.size());
    result.cellHalfEdges.resize(topology.cells.size());

    std::map<StableEdgeKey2D,std::size_t> stableEdges;
    for (const auto& edge:topology.edges) {
        if (edge.id>=topology.edges.size() || edge.v0>=topology.vertices.size() ||
            edge.v1>=topology.vertices.size() || edge.owner>=topology.cells.size() ||
            (edge.neighbour && *edge.neighbour>=topology.cells.size())) {
            result.issues.push_back({EdgeIncidenceIssueCode2D::InvalidTopology,edge.id,
                                     "edge references an out-of-range topology object"});
            continue;
        }
        const auto endpoints=std::minmax(result.stableVertexIds[edge.v0],
                                         result.stableVertexIds[edge.v1]);
        const StableEdgeKey2D key{endpoints.first,endpoints.second};
        result.stableEdgeKeys[edge.id]=key;
        if (!stableEdges.emplace(key,edge.id).second) {
            result.issues.push_back({EdgeIncidenceIssueCode2D::DuplicateStableEdge,edge.id,
                                     "two topology edges have the same stable endpoint key"});
        }
    }
    if (!result.issues.empty()) return result;

    for (const auto& cell:topology.cells) {
        if (cell.id>=topology.cells.size() || cell.vertices.size()<3U ||
            cell.vertices.size()!=cell.edges.size()) {
            result.issues.push_back({EdgeIncidenceIssueCode2D::InvalidCellLoop,cell.id,
                                     "cell vertex and edge loops are not aligned"});
            continue;
        }
        auto& cellUses=result.cellHalfEdges[cell.id];
        cellUses.reserve(cell.edges.size());
        for (std::size_t local=0;local<cell.edges.size();++local) {
            const auto edgeId=cell.edges[local];
            const auto from=cell.vertices[local];
            const auto to=cell.vertices[(local+1U)%cell.vertices.size()];
            if (edgeId>=topology.edges.size() || from>=topology.vertices.size() ||
                to>=topology.vertices.size()) {
                result.issues.push_back({EdgeIncidenceIssueCode2D::InvalidCellLoop,cell.id,
                                         "cell loop references an out-of-range edge or vertex"});
                continue;
            }
            const auto& edge=topology.edges[edgeId];
            if (std::minmax(from,to)!=std::minmax(edge.v0,edge.v1)) {
                result.issues.push_back({EdgeIncidenceIssueCode2D::EdgeEndpointMismatch,edgeId,
                                         "cell half-edge endpoints differ from global edge"});
                continue;
            }
            const auto halfEdgeId=result.halfEdges.size();
            result.halfEdges.push_back({halfEdgeId,edgeId,cell.id,from,to,0U,0U,std::nullopt});
            cellUses.push_back(halfEdgeId);
            result.edgeHalfEdges[edgeId].push_back(halfEdgeId);
        }
        if (cellUses.size()!=cell.edges.size()) continue;
        for (std::size_t local=0;local<cellUses.size();++local) {
            auto& halfEdge=result.halfEdges[cellUses[local]];
            halfEdge.previous=cellUses[(local+cellUses.size()-1U)%cellUses.size()];
            halfEdge.next=cellUses[(local+1U)%cellUses.size()];
        }
    }
    if (!result.issues.empty()) return result;

    for (const auto& edge:topology.edges) {
        auto& uses=result.edgeHalfEdges[edge.id];
        std::sort(uses.begin(),uses.end(),[&](std::size_t lhs,std::size_t rhs) {
            return result.halfEdges[lhs].cell<result.halfEdges[rhs].cell;
        });
        const std::size_t expected=edge.neighbour?2U:1U;
        if (uses.size()!=expected) {
            result.issues.push_back({EdgeIncidenceIssueCode2D::InvalidIncidenceCount,edge.id,
                                     "edge incidence count disagrees with boundary state"});
            continue;
        }
        std::set<std::size_t> incidentCells;
        for (const auto use:uses) incidentCells.insert(result.halfEdges[use].cell);
        std::set<std::size_t> expectedCells{edge.owner};
        if (edge.neighbour) expectedCells.insert(*edge.neighbour);
        if (incidentCells!=expectedCells) {
            result.issues.push_back({EdgeIncidenceIssueCode2D::OwnerNeighbourMismatch,edge.id,
                                     "half-edge cells disagree with owner/neighbour"});
            continue;
        }
        if (edge.neighbour) {
            auto& first=result.halfEdges[uses[0]];
            auto& second=result.halfEdges[uses[1]];
            if (first.from!=second.to || first.to!=second.from) {
                result.issues.push_back({EdgeIncidenceIssueCode2D::TwinOrientationMismatch,
                                         edge.id,"internal half-edges are not oppositely oriented"});
                continue;
            }
            first.twin=second.id;
            second.twin=first.id;
            ++result.audit.internalEdges;
            ++result.audit.twinPairs;
        } else ++result.audit.boundaryEdges;
    }
    result.audit.halfEdges=result.halfEdges.size();
    return result;
}

} // namespace cartmesh2d
