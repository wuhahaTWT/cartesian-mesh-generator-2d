#pragma once

#include "cartmesh2d/cutcell/CutCell2D.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>
#include <memory>

namespace cartmesh2d {

enum class BoundaryPatch2D {
    None,
    EmbeddedBoundary,
    DomainBoundary,
    Unclassified
};

enum class TopologyIssueCode2D {
    InvalidCell,
    DegenerateEdge,
    DuplicateCellSource,
    OrphanInternalEdge,
    NonManifoldEdge,
    UnclassifiedBoundaryEdge,
    OpenCellLoop,
    AreaMismatch
};

struct TopologyIssue2D {
    TopologyIssueCode2D code;
    std::size_t objectId = 0;
    std::string message;
};

struct Vertex2D {
    std::size_t id = 0;
    Point2D point;
};

struct Edge2D {
    std::size_t id = 0;
    std::size_t v0 = 0;
    std::size_t v1 = 0;
    std::size_t owner = 0;
    std::optional<std::size_t> neighbour;
    BoundaryPatch2D patch = BoundaryPatch2D::None;
};

struct TopologyCell2D {
    std::size_t id = 0;
    std::size_t sourceId = 0;
    std::uint64_t sourceKey = 0;
    std::vector<std::size_t> sourceLineage;
    double geometryArea = 0.0;
    std::vector<std::size_t> vertices;
    std::vector<std::size_t> edges;
    std::vector<std::uint64_t> refinementLineageKeys;
};

struct TopologyAudit2D {
    std::size_t duplicateVertices = 0;
    std::size_t duplicateEdges = 0;
    std::size_t orphanInternalEdges = 0;
    std::size_t nonManifoldEdges = 0;
    std::size_t unclassifiedBoundaryEdges = 0;
    std::size_t openCellLoops = 0;
    std::size_t areaMismatches = 0;

    [[nodiscard]] bool pass() const noexcept {
        return duplicateVertices == 0 && duplicateEdges == 0 &&
               orphanInternalEdges == 0 && nonManifoldEdges == 0 &&
               unclassifiedBoundaryEdges == 0 && openCellLoops == 0 &&
               areaMismatches == 0;
    }
};

struct TopologyMesh2D {
    std::vector<Vertex2D> vertices;
    std::vector<Edge2D> edges;
    std::vector<TopologyCell2D> cells;
    std::vector<TopologyIssue2D> issues;
    TopologyAudit2D audit;
    std::shared_ptr<IntersectionRegistry2D> constructionRegistry;
    std::vector<std::size_t> canonicalVertexIds;
    std::size_t sharedPartitionCount = 0;
    std::size_t sharedPartitionCacheHits = 0;
    // Physical wall fragments already classified in this construction
    // transaction. They must survive agglomeration/partition rebuilds because
    // a bounded canonical weld need not remain collinear with the input loop.
    std::vector<Segment2D> embeddedBoundaryFragments;

    [[nodiscard]] bool valid() const noexcept {
        return issues.empty() && audit.pass();
    }
};

[[nodiscard]] TopologyMesh2D buildGlobalTopology(
    const std::vector<CutCell2D>& cutCells,
    const Domain2D& domain,
    const BoundaryLoop& boundary,
    const TolerancePolicy& tol = {});

[[nodiscard]] TopologyMesh2D buildGlobalTopology(
    const std::vector<CutCell2D>& cutCells,
    const Domain2D& domain,
    const BoundaryRegion2D& boundary,
    const TolerancePolicy& tol = {},
    std::shared_ptr<IntersectionRegistry2D> registry = {},
    const std::vector<Segment2D>& knownEmbeddedBoundary = {});

// Monotonic count of buildGlobalTopology calls in this process. Instrumentation
// only: it lets a caller prove by measurement, rather than by assertion, that a
// code path performed no global rebuild. Single-threaded by construction, like
// the rest of the generator.
[[nodiscard]] std::size_t globalTopologyBuildCount2D() noexcept;

// Total seconds spent inside buildGlobalTopology, and the total input cell
// count summed over all calls. Together with the call count these separate
// "many small patch builds" from "few whole-mesh rebuilds".
[[nodiscard]] double globalTopologyBuildSeconds2D() noexcept;
[[nodiscard]] std::size_t globalTopologyBuildInputCells2D() noexcept;

} // namespace cartmesh2d
