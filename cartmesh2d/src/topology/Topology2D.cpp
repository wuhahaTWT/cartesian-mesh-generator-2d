#include "cartmesh2d/topology/Topology2D.hpp"

#include <algorithm>
#include <cmath>
#include <map>
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

[[nodiscard]] bool pointOnBoundary(const Point2D& p, const BoundaryLoop& boundary,
                                   const TolerancePolicy& tol) noexcept {
    const auto& vertices = boundary.vertices();
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        if (pointOnSegment(p, {vertices[i], vertices[(i + 1) % vertices.size()]}, tol)) {
            return true;
        }
    }
    return false;
}

[[nodiscard]] bool segmentOnInputBoundary(const Point2D& a, const Point2D& b,
                                          const BoundaryLoop& boundary,
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
        if (pointOnSegment(a, fragment, tol) &&
            pointOnSegment(b, fragment, tol) &&
            pointOnSegment(mid, fragment, tol)) {
            return true;
        }
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

[[nodiscard]] bool isCartesianSegment(const Point2D& a, const Point2D& b,
                                      const TolerancePolicy& tol) noexcept {
    return scalarNear(a.x, b.x, tol) || scalarNear(a.y, b.y, tol);
}

[[nodiscard]] std::vector<Point2D> collectCanonicalPoints(
    const std::vector<CutCell2D>& cells, const TolerancePolicy& tol) {
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
            if (p.x - it->x > tol.scale(std::max({1.0, std::abs(p.x), std::abs(it->x)}))) break;
            if (pointNear(p, *it, tol)) {
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
                                       const TolerancePolicy& tol) noexcept {
    for (const auto& vertex : vertices) {
        if (pointNear(p, vertex.point, tol)) return vertex.id;
    }
    return vertices.size();
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
        if (cells[i - 1].sourceKey == cells[i].sourceKey &&
            cells[i - 1].sourceId == cells[i].sourceId) {
            mesh.issues.push_back({TopologyIssueCode2D::DuplicateCellSource, i,
                                   "duplicate source cell identity"});
            return mesh;
        }
    }

    const auto canonicalPoints = collectCanonicalPoints(cells, tol);
    mesh.vertices.reserve(canonicalPoints.size());
    for (std::size_t i = 0; i < canonicalPoints.size(); ++i) {
        mesh.vertices.push_back({i, canonicalPoints[i]});
    }

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
                mesh.issues.push_back({TopologyIssueCode2D::DegenerateEdge, cellId,
                                       "degenerate polygon edge"});
                return mesh;
            }

            std::vector<std::pair<double, std::size_t>> onEdge;
            if (isCartesianSegment(a, b, tol)) {
                for (const auto& vertex : mesh.vertices) {
                    if (pointOnSegment(vertex.point, {a, b}, tol)) {
                        const double t = segmentParameter(vertex.point, a, b);
                        if (t >= -tol.scale() && t <= 1.0 + tol.scale()) {
                            onEdge.push_back({t, vertex.id});
                        }
                    }
                }
            } else {
                const std::size_t aId = findVertexId(a, mesh.vertices, tol);
                const std::size_t bId = findVertexId(b, mesh.vertices, tol);
                if (aId < mesh.vertices.size()) onEdge.push_back({0.0, aId});
                if (bId < mesh.vertices.size()) onEdge.push_back({1.0, bId});
            }

            std::sort(onEdge.begin(), onEdge.end(), [](const auto& lhs, const auto& rhs) {
                if (lhs.first != rhs.first) return lhs.first < rhs.first;
                return lhs.second < rhs.second;
            });
            onEdge.erase(std::unique(onEdge.begin(), onEdge.end(),
                                     [](const auto& lhs, const auto& rhs) {
                                         return lhs.second == rhs.second;
                                     }),
                         onEdge.end());
            if (onEdge.size() < 2) {
                mesh.issues.push_back({TopologyIssueCode2D::OpenCellLoop, cellId,
                                       "cell edge could not be mapped to topology vertices"});
                return mesh;
            }

            for (std::size_t k = 0; k + 1 < onEdge.size(); ++k) {
                const std::size_t from = onEdge[k].second;
                const std::size_t to = onEdge[k + 1].second;
                if (from == to) continue;
                const Point2D& p0 = mesh.vertices[from].point;
                const Point2D& p1 = mesh.vertices[to].point;
                const double fragmentLength2 = squaredNorm(p1 - p0);
                const double lengthEps = tol.scale(std::max({1.0, std::abs(p0.x), std::abs(p0.y),
                                                             std::abs(p1.x), std::abs(p1.y)}));
                if (fragmentLength2 <= lengthEps * lengthEps) continue;
                if (loop.empty()) loop.push_back(from);
                if (loop.back() != from) {
                    mesh.issues.push_back({TopologyIssueCode2D::OpenCellLoop, cellId,
                                           "split fragments are not contiguous"});
                    return mesh;
                }
                loop.push_back(to);
                const EdgeKey key{std::min(from, to), std::max(from, to)};
                uses[key].push_back({cellId, from, to});
            }
        }
        if (loop.size() > 1 && loop.front() == loop.back()) loop.pop_back();
        loop.erase(std::unique(loop.begin(), loop.end()), loop.end());
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
                       pointOnDomainBoundary(a, domain, tol) &&
                       pointOnDomainBoundary(b, domain, tol)) {
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
