#pragma once

#include "cartmesh2d/topology/Topology2D.hpp"

#include <compare>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace cartmesh2d {

// Stable within a construction revision. When topology was built from an
// IntersectionRegistry2D the endpoint ids are construction handles; otherwise
// they are the deterministic dense topology vertex ids.
struct StableEdgeKey2D {
    StableVertexId2D v0 = 0;
    StableVertexId2D v1 = 0;
    auto operator<=>(const StableEdgeKey2D&) const = default;
};

struct HalfEdgeLite2D {
    std::size_t id = 0;
    std::size_t edge = 0;
    std::size_t cell = 0;
    std::size_t from = 0;
    std::size_t to = 0;
    std::size_t previous = 0;
    std::size_t next = 0;
    std::optional<std::size_t> twin;
};

enum class EdgeIncidenceIssueCode2D {
    InvalidTopology,
    InvalidCellLoop,
    EdgeEndpointMismatch,
    InvalidIncidenceCount,
    OwnerNeighbourMismatch,
    TwinOrientationMismatch,
    DuplicateStableEdge
};

struct EdgeIncidenceIssue2D {
    EdgeIncidenceIssueCode2D code = EdgeIncidenceIssueCode2D::InvalidTopology;
    std::size_t objectId = 0;
    std::string message;
};

struct EdgeIncidenceAudit2D {
    std::size_t boundaryEdges = 0;
    std::size_t internalEdges = 0;
    std::size_t halfEdges = 0;
    std::size_t twinPairs = 0;
};

struct EdgeIncidenceStore2D {
    std::uint64_t revision = 0;
    std::vector<StableVertexId2D> stableVertexIds;
    std::vector<StableEdgeKey2D> stableEdgeKeys;
    std::vector<HalfEdgeLite2D> halfEdges;
    std::vector<std::vector<std::size_t>> edgeHalfEdges;
    std::vector<std::vector<std::size_t>> cellHalfEdges;
    std::vector<EdgeIncidenceIssue2D> issues;
    EdgeIncidenceAudit2D audit;

    [[nodiscard]] bool valid() const noexcept { return issues.empty(); }
};

// Builds the minimal directed incidence needed by patch transactions. It does
// not own geometry and does not mutate TopologyMesh2D.
[[nodiscard]] EdgeIncidenceStore2D buildEdgeIncidenceStore2D(
    const TopologyMesh2D& topology, std::uint64_t revision = 0);

} // namespace cartmesh2d
