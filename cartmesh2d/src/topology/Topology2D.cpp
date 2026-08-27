#include "cartmesh2d/topology/Topology2D.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <set>
#include <tuple>
#include <utility>

namespace cartmesh2d {
namespace {

[[nodiscard]] bool scalarNear(double a, double b, const TolerancePolicy& tol) noexcept {
    return std::abs(a - b) <= tol.scale(std::max({1.0, std::abs(a), std::abs(b)}));
}

[[nodiscard]] bool pointNear(const Point2D& a, const Point2D& b,
                             const TolerancePolicy& tol) noexcept {
    return scalarNear(a.x, b.x, tol) && scalarNear(a.y, b.y, tol);
}

[[nodiscard]] double segmentParameter(const Point2D& p, const Point2D& a,
                                      const Point2D& b) noexcept {
    const Vector2D d = b - a;
    const double denom = squaredNorm(d);
    if (denom <= 0.0) return 0.0;
    return dot(p - a, d) / denom;
}

[[nodiscard]] bool pointOnBoundary(const Point2D& p, const BoundaryRegion2D& boundary,
                                   const TolerancePolicy& tol) noexcept {
    for (const auto& loop : boundary.loops()) {
        const auto& vertices = loop.vertices();
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            if (pointOnSegment(p, {vertices[i], vertices[(i + 1) % vertices.size()]}, tol)) return true;
        }
    }
    return false;
}

[[nodiscard]] bool segmentOnInputBoundary(const Point2D& a, const Point2D& b,
                                          const BoundaryRegion2D& boundary,
                                          const TolerancePolicy& tol) noexcept {
    const Point2D mid{0.5 * (a.x + b.x), 0.5 * (a.y + b.y)};
    return pointOnBoundary(a, boundary, tol) && pointOnBoundary(b, boundary, tol) &&
           pointOnBoundary(mid, boundary, tol);
}

[[nodiscard]] bool segmentOnEmbeddedFragments(
    const Point2D& a, const Point2D& b,
    const std::vector<Segment2D>& fragments,
    const TolerancePolicy& tol) noexcept {
    const Point2D mid{0.5 * (a.x + b.x), 0.5 * (a.y + b.y)};
    for (const auto& fragment : fragments) {
        if (pointOnSegment(a, fragment, tol) && pointOnSegment(b, fragment, tol) &&
            pointOnSegment(mid, fragment, tol)) return true;
    }
    return false;
}

[[nodiscard]] bool pointOnDomainBoundary(const Point2D& p, const Domain2D& domain,
                                         const TolerancePolicy& tol) noexcept {
    const auto& box = domain.bounds;
    return (scalarNear(p.x, box.min.x, tol) || scalarNear(p.x, box.max.x, tol) ||
            scalarNear(p.y, box.min.y, tol) || scalarNear(p.y, box.max.y, tol)) &&
           box.contains(p, tol);
}

[[nodiscard]] bool segmentOnDomainBoundary(const Point2D& a, const Point2D& b,
                                           const Domain2D& domain,
                                           const TolerancePolicy& tol) noexcept {
    const auto& box = domain.bounds;
    return (scalarNear(a.x, box.min.x, tol) && scalarNear(b.x, box.min.x, tol)) ||
           (scalarNear(a.x, box.max.x, tol) && scalarNear(b.x, box.max.x, tol)) ||
           (scalarNear(a.y, box.min.y, tol) && scalarNear(b.y, box.min.y, tol)) ||
           (scalarNear(a.y, box.max.y, tol) && scalarNear(b.y, box.max.y, tol));
}

[[nodiscard]] std::vector<Point2D> collectCanonicalPoints(
    const std::vector<CutCell2D>& cells, double coordinateEps) {
    std::vector<Point2D> raw;
    for (const auto& cell : cells) {
        if (cell.kind == CutCellKind::Empty || cell.kind == CutCellKind::Unsupported) continue;
        raw.insert(raw.end(), cell.fluidPolygon.vertices.begin(), cell.fluidPolygon.vertices.end());
    }
    std::sort(raw.begin(), raw.end(), [](const Point2D& a, const Point2D& b) {
        return std::tie(a.x, a.y) < std::tie(b.x, b.y);
    });
    std::vector<Point2D> result;
    for (const auto& p : raw) {
        bool duplicate = false;
        for (auto it = result.rbegin(); it != result.rend(); ++it) {
            if (p.x - it->x > coordinateEps) break;
            if (std::abs(p.x - it->x) <= coordinateEps &&
                std::abs(p.y - it->y) <= coordinateEps) {
                duplicate = true;
                break;
            }
        }
        if (!duplicate) result.push_back(p);
    }
    return result;
}

[[nodiscard]] std::size_t findVertexId(const Point2D& p,
                                       const std::vector<Vertex2D>& vertices,
                                       double coordinateEps) noexcept {
    const auto begin = std::lower_bound(
        vertices.begin(), vertices.end(), p.x - coordinateEps,
        [](const Vertex2D& vertex, double x) { return vertex.point.x < x; });
    for (auto it = begin; it != vertices.end() &&
                          it->point.x <= p.x + coordinateEps; ++it) {
        if (std::abs(p.x - it->point.x) <= coordinateEps &&
            std::abs(p.y - it->point.y) <= coordinateEps) return it->id;
    }
    return vertices.size();
}

void appendCartesianEdgeVertices(
    const Point2D& a, const Point2D& b,
    const std::vector<Vertex2D>& vertices,
    const std::vector<std::size_t>& verticesByY,
    double coordinateEps, const TolerancePolicy& tol,
    std::size_t aId, std::size_t bId,
    std::vector<std::pair<double, std::size_t>>& onEdge) {
    const bool vertical = std::abs(a.x - b.x) <= coordinateEps;
    const double fixed = vertical ? a.x : a.y;
    if (vertical) {
        const auto begin = std::lower_bound(
            vertices.begin(), vertices.end(), fixed - coordinateEps,
            [](const Vertex2D& vertex, double value) { return vertex.point.x < value; });
        for (auto it = begin; it != vertices.end() &&
                              it->point.x <= fixed + coordinateEps; ++it) {
            if (it->id == aId || it->id == bId) continue;
            if (!pointOnSegment(it->point, {a, b}, tol)) continue;
            onEdge.push_back({segmentParameter(it->point, a, b), it->id});
        }
        return;
    }

    const auto begin = std::lower_bound(
        verticesByY.begin(), verticesByY.end(), fixed - coordinateEps,
        [&](std::size_t id, double value) { return vertices[id].point.y < value; });
    for (auto it = begin; it != verticesByY.end() &&
                          vertices[*it].point.y <= fixed + coordinateEps; ++it) {
        const std::size_t id = *it;
        if (id == aId || id == bId) continue;
        if (!pointOnSegment(vertices[id].point, {a, b}, tol)) continue;
        onEdge.push_back({segmentParameter(vertices[id].point, a, b), id});
    }
}

// H4-2 requires a common partition on a non-axis-aligned outer envelope:
// remainder Cut-cells contribute grid/envelope intersection vertices while a
// boundary-layer quad initially owns the complete envelope segment. Search the
// canonical x-sorted vertex set and split that segment at the exact same IDs.
// The established Cartesian fast path above remains unchanged.
void appendGeneralEdgeVertices(
    const Point2D& a, const Point2D& b,
    const std::vector<Vertex2D>& vertices,
    const std::vector<std::size_t>& verticesByY,
    double coordinateEps, const TolerancePolicy& tol,
    std::size_t aId, std::size_t bId,
    std::vector<std::pair<double, std::size_t>>& onEdge) {
    if (std::abs(a.x - b.x) <= coordinateEps ||
        std::abs(a.y - b.y) <= coordinateEps) {
        appendCartesianEdgeVertices(a, b, vertices, verticesByY,
                                    coordinateEps, tol, aId, bId, onEdge);
        return;
    }
    const double minX = std::min(a.x, b.x) - coordinateEps;
    const double maxX = std::max(a.x, b.x) + coordinateEps;
    const auto begin = std::lower_bound(
        vertices.begin(), vertices.end(), minX,
        [](const Vertex2D& vertex, double value) {
            return vertex.point.x < value;
        });
    for (auto it = begin; it != vertices.end() && it->point.x <= maxX; ++it) {
        if (it->id == aId || it->id == bId) continue;
        if (!pointOnSegment(it->point, {a, b}, tol)) continue;
        onEdge.push_back({segmentParameter(it->point, a, b), it->id});
    }
}

struct EdgeKey {
    std::size_t a = 0;
    std::size_t b = 0;
    [[nodiscard]] bool operator<(const EdgeKey& rhs) const noexcept {
        return std::tie(a, b) < std::tie(rhs.a, rhs.b);
    }
};

struct EdgeUse {
    std::size_t cellId = 0;
    std::size_t from = 0;
    std::size_t to = 0;
};

[[nodiscard]] double polygonAreaFromVertexLoop(const std::vector<std::size_t>& loop,
                                               const std::vector<Vertex2D>& vertices) noexcept {
    if (loop.size() < 3) return 0.0;
    double twiceArea = 0.0;
    for (std::size_t i = 0; i < loop.size(); ++i) {
        const auto& a = vertices[loop[i]].point;
        const auto& b = vertices[loop[(i + 1) % loop.size()]].point;
        twiceArea += a.x * b.y - b.x * a.y;
    }
    return 0.5 * std::abs(twiceArea);
}

} // namespace

TopologyMesh2D buildGlobalTopology(const std::vector<CutCell2D>& inputCells,
                                   const Domain2D& domain,
                                   const BoundaryLoop& boundary,
                                   const TolerancePolicy& tol) {
    return buildGlobalTopology(inputCells,domain,BoundaryRegion2D(boundary),tol);
}

TopologyMesh2D buildGlobalTopology(const std::vector<CutCell2D>& inputCells,
                                   const Domain2D& domain,
                                   const BoundaryRegion2D& boundary,
                                   const TolerancePolicy& tol) {
    TopologyMesh2D mesh;
    if (!domain.valid(tol) || !boundary.diagnose(tol).valid()) {
        mesh.issues.push_back({TopologyIssueCode2D::InvalidCell, 0,
                               "invalid domain or boundary supplied to topology builder"});
        return mesh;
    }

    std::vector<CutCell2D> cells;
    for (const auto& cell : inputCells) {
        if (cell.kind == CutCellKind::Empty) continue;
        if (!cell.valid() || cell.kind == CutCellKind::Unsupported ||
            cell.fluidPolygon.vertices.size() < 3 || !(cell.area > 0.0)) {
            mesh.issues.push_back({TopologyIssueCode2D::InvalidCell, cell.sourceId,
                                   "non-empty topology input cell is invalid"});
            return mesh;
        }
        cells.push_back(cell);
    }
    std::sort(cells.begin(), cells.end(), [](const CutCell2D& a, const CutCell2D& b) {
        return std::tie(a.sourceKey, a.sourceId) < std::tie(b.sourceKey, b.sourceId);
    });
    for (std::size_t i = 1; i < cells.size(); ++i) {
        if (cells[i - 1].sourceKey == cells[i].sourceKey && cells[i - 1].sourceId == cells[i].sourceId) {
            mesh.issues.push_back({TopologyIssueCode2D::DuplicateCellSource, i,
                                   "duplicate source cell identity"});
            return mesh;
        }
    }

    // Vertex canonicalization is a local operation.  Scaling its tolerance by
    // the whole domain can collapse two legitimate vertices of a very small
    // high-level Cut-cell.  Use the smallest positive source-cell extent so
    // shared round-off is absorbed without erasing resolved geometry.
    double minCellExtent = std::numeric_limits<double>::infinity();
    for (const auto& cell : cells) {
        const double width = cell.backgroundBounds.max.x - cell.backgroundBounds.min.x;
        const double height = cell.backgroundBounds.max.y - cell.backgroundBounds.min.y;
        if (width > 0.0) minCellExtent = std::min(minCellExtent, width);
        if (height > 0.0) minCellExtent = std::min(minCellExtent, height);
    }
    if (!std::isfinite(minCellExtent)) {
        mesh.issues.push_back({TopologyIssueCode2D::InvalidCell, 0,
                               "topology inputs have no positive background-cell extent"});
        return mesh;
    }
    const double coordinateEps = tol.scale(minCellExtent);
    const auto canonicalPoints = collectCanonicalPoints(cells, coordinateEps);
    mesh.vertices.reserve(canonicalPoints.size());
    for (std::size_t i = 0; i < canonicalPoints.size(); ++i) mesh.vertices.push_back({i, canonicalPoints[i]});
    std::vector<std::size_t> verticesByY(mesh.vertices.size());
    std::iota(verticesByY.begin(), verticesByY.end(), 0);
    std::sort(verticesByY.begin(), verticesByY.end(), [&](std::size_t lhs, std::size_t rhs) {
        const auto& a = mesh.vertices[lhs].point;
        const auto& b = mesh.vertices[rhs].point;
        return std::tie(a.y, a.x, lhs) < std::tie(b.y, b.x, rhs);
    });

    std::map<EdgeKey, std::vector<EdgeUse>> uses;
    mesh.cells.reserve(cells.size());

    for (std::size_t cellId = 0; cellId < cells.size(); ++cellId) {
        const auto& source = cells[cellId];
        TopologyCell2D topo;
        topo.id = cellId;
        topo.sourceId = source.sourceId;
        topo.sourceKey = source.sourceKey;
        topo.geometryArea = source.area;

        const auto& polygon = source.fluidPolygon.vertices;
        std::vector<std::size_t> loop;
        for (std::size_t e = 0; e < polygon.size(); ++e) {
            const Point2D a = polygon[e];
            const Point2D b = polygon[(e + 1) % polygon.size()];
            if (pointNear(a, b, tol)) {
                mesh.issues.push_back({TopologyIssueCode2D::DegenerateEdge, cellId, "degenerate polygon edge"});
                return mesh;
            }

            const std::size_t aId = findVertexId(a, mesh.vertices, coordinateEps);
            const std::size_t bId = findVertexId(b, mesh.vertices, coordinateEps);
            if (aId >= mesh.vertices.size() || bId >= mesh.vertices.size() || aId == bId) {
                std::ostringstream detail;
                detail << std::setprecision(17)
                       << "polygon edge endpoints could not be uniquely canonicalized"
                       << " edge=" << e
                       << " a=(" << a.x << ',' << a.y << ')'
                       << " b=(" << b.x << ',' << b.y << ')'
                       << " a_id=" << aId << " b_id=" << bId
                       << " coordinate_eps=" << coordinateEps;
                mesh.issues.push_back({TopologyIssueCode2D::OpenCellLoop, cellId,
                                       detail.str()});
                return mesh;
            }

            std::vector<std::pair<double, std::size_t>> onEdge{{0.0, aId}, {1.0, bId}};
            const double edgeLength = std::sqrt(squaredNorm(b - a));
            const double tEps = coordinateEps /
                                std::max(edgeLength, coordinateEps);
            const std::size_t before = onEdge.size();
            appendGeneralEdgeVertices(a, b, mesh.vertices, verticesByY,
                                      coordinateEps, tol, aId, bId, onEdge);
            onEdge.erase(std::remove_if(onEdge.begin() + static_cast<std::ptrdiff_t>(before),
                                        onEdge.end(), [&](const auto& item) {
                                            return item.first <= tEps ||
                                                   item.first >= 1.0 - tEps;
                                        }), onEdge.end());

            std::sort(onEdge.begin(), onEdge.end(), [](const auto& lhs, const auto& rhs) {
                if (lhs.first != rhs.first) return lhs.first < rhs.first;
                return lhs.second < rhs.second;
            });
            onEdge.erase(std::unique(onEdge.begin(), onEdge.end(),
                                     [](const auto& lhs, const auto& rhs) { return lhs.second == rhs.second; }),
                         onEdge.end());

            for (std::size_t k = 0; k + 1 < onEdge.size(); ++k) {
                const std::size_t from = onEdge[k].second;
                const std::size_t to = onEdge[k + 1].second;
                if (from == to) continue;
                const Point2D& p0 = mesh.vertices[from].point;
                const Point2D& p1 = mesh.vertices[to].point;
                if (squaredNorm(p1 - p0) <= coordinateEps * coordinateEps) continue;
                if (loop.empty()) loop.push_back(from);
                if (loop.back() != from) {
                    mesh.issues.push_back({TopologyIssueCode2D::OpenCellLoop, cellId,
                                           "split fragments are not contiguous"});
                    return mesh;
                }
                loop.push_back(to);
                uses[{std::min(from, to), std::max(from, to)}].push_back({cellId, from, to});
            }
        }
        if (loop.size() > 1 && loop.front() == loop.back()) loop.pop_back();
        if (loop.size() < 3) {
            mesh.issues.push_back({TopologyIssueCode2D::OpenCellLoop, cellId,
                                   "cell topology loop has fewer than three vertices"});
            return mesh;
        }
        topo.vertices = std::move(loop);
        mesh.cells.push_back(std::move(topo));
    }

    std::map<EdgeKey, std::size_t> edgeIds;
    for (const auto& [key, incident] : uses) {
        if (incident.empty()) continue;
        Edge2D edge;
        edge.id = mesh.edges.size();
        edge.v0 = key.a;
        edge.v1 = key.b;
        edge.owner = incident.front().cellId;
        if (incident.size() == 2) {
            edge.neighbour = incident[1].cellId;
            edge.patch = BoundaryPatch2D::None;
        } else if (incident.size() == 1) {
            const auto& a = mesh.vertices[key.a].point;
            const auto& b = mesh.vertices[key.b].point;
            const auto& ownerSource = cells[edge.owner];
            if (segmentOnEmbeddedFragments(a, b, ownerSource.embeddedBoundary, tol) ||
                segmentOnInputBoundary(a, b, boundary, tol)) {
                edge.patch = BoundaryPatch2D::EmbeddedBoundary;
            } else if (segmentOnDomainBoundary(a, b, domain, tol) &&
                       pointOnDomainBoundary(a, domain, tol) && pointOnDomainBoundary(b, domain, tol)) {
                edge.patch = BoundaryPatch2D::DomainBoundary;
            } else {
                edge.patch = BoundaryPatch2D::Unclassified;
                ++mesh.audit.unclassifiedBoundaryEdges;
                mesh.issues.push_back({TopologyIssueCode2D::UnclassifiedBoundaryEdge, edge.id,
                                       "single-owner edge is neither an owner embedded fragment nor domain boundary"});
            }
        } else {
            ++mesh.audit.nonManifoldEdges;
            mesh.issues.push_back({TopologyIssueCode2D::NonManifoldEdge, edge.id,
                                   "edge has more than two incident cells"});
        }
        edgeIds[key] = edge.id;
        mesh.edges.push_back(edge);
    }

    for (auto& cell : mesh.cells) {
        cell.edges.clear();
        for (std::size_t i = 0; i < cell.vertices.size(); ++i) {
            const auto a = cell.vertices[i];
            const auto b = cell.vertices[(i + 1) % cell.vertices.size()];
            const EdgeKey key{std::min(a, b), std::max(a, b)};
            const auto it = edgeIds.find(key);
            if (it == edgeIds.end()) {
                ++mesh.audit.openCellLoops;
                mesh.issues.push_back({TopologyIssueCode2D::OpenCellLoop, cell.id,
                                       "cell loop references a missing global edge"});
                continue;
            }
            cell.edges.push_back(it->second);
        }
        if (cell.edges.size() != cell.vertices.size()) ++mesh.audit.openCellLoops;
        const double reconstructed = polygonAreaFromVertexLoop(cell.vertices, mesh.vertices);
        const double eps = tol.scale(std::max({1.0, reconstructed, cell.geometryArea}));
        if (std::abs(reconstructed - cell.geometryArea) > eps) {
            ++mesh.audit.areaMismatches;
            mesh.issues.push_back({TopologyIssueCode2D::AreaMismatch, cell.id,
                                   "topology loop area differs from source CutCell area"});
        }
    }

    for (const auto& edge : mesh.edges) {
        if (edge.neighbour && edge.owner == *edge.neighbour) {
            ++mesh.audit.orphanInternalEdges;
            mesh.issues.push_back({TopologyIssueCode2D::OrphanInternalEdge, edge.id,
                                   "internal edge owner equals neighbour"});
        }
    }
    return mesh;
}

} // namespace cartmesh2d
