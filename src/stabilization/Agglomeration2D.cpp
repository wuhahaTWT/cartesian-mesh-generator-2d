#include "cartmesh2d/stabilization/SmallCell2D.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <tuple>
#include <utility>

namespace cartmesh2d {
namespace {

struct Dsu2D {
    explicit Dsu2D(std::size_t n) : parent(n), rank(n, 0) {
        for (std::size_t i = 0; i < n; ++i) parent[i] = i;
    }

    std::size_t find(std::size_t x) {
        if (parent[x] != x) parent[x] = find(parent[x]);
        return parent[x];
    }

    void unite(std::size_t a, std::size_t b) {
        a = find(a);
        b = find(b);
        if (a == b) return;
        if (rank[a] < rank[b]) std::swap(a, b);
        parent[b] = a;
        if (rank[a] == rank[b]) ++rank[a];
    }

    std::vector<std::size_t> parent;
    std::vector<unsigned> rank;
};

struct CanonicalEdgeKey2D {
    std::size_t a = 0;
    std::size_t b = 0;

    [[nodiscard]] bool operator<(const CanonicalEdgeKey2D& rhs) const noexcept {
        return std::tie(a, b) < std::tie(rhs.a, rhs.b);
    }
};

struct DirectedBoundaryUse2D {
    std::size_t from = 0;
    std::size_t to = 0;
};

[[nodiscard]] bool nearlySameArea(double a, double b, const TolerancePolicy& tol) noexcept {
    // Area comparisons require an area-dimensional tolerance.  Reusing the
    // coordinate tolerance made valid level-17 cells (area ~1e-10) appear
    // numerically indistinguishable from zero.
    const double eps = tol.absolute * tol.absolute +
                       tol.relative * std::max(std::abs(a), std::abs(b));
    return std::abs(a - b) <= eps;
}

[[nodiscard]] bool removableCollinearVertex(const Point2D& prev,
                                            const Point2D& cur,
                                            const Point2D& next,
                                            const TolerancePolicy& tol) noexcept {
    const Vector2D a = cur - prev;
    const Vector2D b = next - cur;
    const double lenA = std::sqrt(squaredNorm(a));
    const double lenB = std::sqrt(squaredNorm(b));
    const double coordinateScale =
        std::max({1.0, std::abs(prev.x), std::abs(prev.y),
                  std::abs(cur.x), std::abs(cur.y),
                  std::abs(next.x), std::abs(next.y)});
    const double lengthEps = tol.scale(coordinateScale);
    if (lenA <= lengthEps || lenB <= lengthEps) return false;
    const double sinAngle = std::abs(cross(a, b)) / (lenA * lenB);
    return sinAngle <= tol.scale(1.0) && dot(a, b) >= 0.0;
}

[[nodiscard]] std::optional<Polygon2D> polygonFromTopologyGroup(
    const std::vector<std::size_t>& members,
    const TopologyMesh2D& topology,
    const TolerancePolicy& tol,
    AgglomerationIssue2D& issue) {

    std::map<CanonicalEdgeKey2D, std::vector<DirectedBoundaryUse2D>> uses;
    for (const std::size_t cellId : members) {
        if (cellId >= topology.cells.size()) {
            issue = {AgglomerationIssueCode2D::InvalidTopologyReference, cellId,
                     "agglomeration group references an invalid topology cell"};
            return std::nullopt;
        }
        const auto& cell = topology.cells[cellId];
        if (cell.vertices.size() < 3 || cell.vertices.size() != cell.edges.size()) {
            issue = {AgglomerationIssueCode2D::InvalidTopologyReference, cellId,
                     "source topology cell has an invalid vertex/edge loop"};
            return std::nullopt;
        }
        for (std::size_t i = 0; i < cell.vertices.size(); ++i) {
            const std::size_t from = cell.vertices[i];
            const std::size_t to = cell.vertices[(i + 1) % cell.vertices.size()];
            if (from >= topology.vertices.size() || to >= topology.vertices.size() || from == to) {
                issue = {AgglomerationIssueCode2D::InvalidTopologyReference, cellId,
                         "source topology loop references an invalid vertex"};
                return std::nullopt;
            }
            const CanonicalEdgeKey2D key{std::min(from, to), std::max(from, to)};
            uses[key].push_back({from, to});
        }
    }

    std::vector<DirectedBoundaryUse2D> boundaryEdges;
    for (const auto& [key, edgeUses] : uses) {
        (void)key;
        if (edgeUses.size() == 1) {
            boundaryEdges.push_back(edgeUses.front());
        } else if (edgeUses.size() == 2) {
            if (!(edgeUses[0].from == edgeUses[1].to && edgeUses[0].to == edgeUses[1].from)) {
                issue = {AgglomerationIssueCode2D::DisconnectedBoundary, edgeUses.front().from,
                         "shared group edge does not appear with opposite orientations"};
                return std::nullopt;
            }
        } else {
            issue = {AgglomerationIssueCode2D::InvalidTopologyReference, edgeUses.front().from,
                     "group edge has more than two incident member cells"};
            return std::nullopt;
        }
    }

    if (boundaryEdges.size() < 3) {
        issue = {AgglomerationIssueCode2D::DegenerateMergedPolygon, 0,
                 "agglomerated cell has fewer than three exterior edge fragments"};
        return std::nullopt;
    }

    std::map<std::size_t, std::size_t> next;
    std::map<std::size_t, std::size_t> incoming;
    for (const auto& edge : boundaryEdges) {
        if (next.contains(edge.from) || incoming.contains(edge.to)) {
            issue = {AgglomerationIssueCode2D::DisconnectedBoundary, edge.from,
                     "agglomerated boundary is branched or has multiple outgoing/incoming edges"};
            return std::nullopt;
        }
        next[edge.from] = edge.to;
        incoming[edge.to] = edge.from;
    }

    if (next.size() != boundaryEdges.size() || incoming.size() != boundaryEdges.size()) {
        issue = {AgglomerationIssueCode2D::DisconnectedBoundary, 0,
                 "agglomerated exterior edges do not form a one-in/one-out boundary graph"};
        return std::nullopt;
    }
    for (const auto& [v, to] : next) {
        (void)to;
        if (!incoming.contains(v)) {
            issue = {AgglomerationIssueCode2D::DisconnectedBoundary, v,
                     "agglomerated boundary has an open chain"};
            return std::nullopt;
        }
    }

    const std::size_t start = next.begin()->first;
    std::vector<std::size_t> loop;
    loop.reserve(boundaryEdges.size());
    std::size_t current = start;
    for (std::size_t step = 0; step <= boundaryEdges.size(); ++step) {
        if (step == boundaryEdges.size()) {
            if (current != start) {
                issue = {AgglomerationIssueCode2D::DisconnectedBoundary, current,
                         "agglomerated boundary did not close"};
                return std::nullopt;
            }
            break;
        }
        if (current == start && !loop.empty()) {
            issue = {AgglomerationIssueCode2D::MultipleBoundaryLoops, current,
                     "agglomerated exterior contains more than one boundary loop"};
            return std::nullopt;
        }
        loop.push_back(current);
        const auto it = next.find(current);
        if (it == next.end()) {
            issue = {AgglomerationIssueCode2D::DisconnectedBoundary, current,
                     "agglomerated boundary traversal hit a missing continuation"};
            return std::nullopt;
        }
        current = it->second;
    }

    if (loop.size() != boundaryEdges.size()) {
        issue = {AgglomerationIssueCode2D::MultipleBoundaryLoops, start,
                 "agglomerated exterior contains multiple disconnected loops"};
        return std::nullopt;
    }

    Polygon2D polygon;
    polygon.vertices.reserve(loop.size());
    for (const std::size_t vertexId : loop) {
        polygon.vertices.push_back(topology.vertices[vertexId].point);
    }
    if (polygon.signedArea() < 0.0) {
        std::reverse(polygon.vertices.begin(), polygon.vertices.end());
    }

    bool simplified = true;
    while (simplified && polygon.vertices.size() > 3) {
        simplified = false;
        for (std::size_t i = 0; i < polygon.vertices.size(); ++i) {
            const std::size_t prev = (i + polygon.vertices.size() - 1) % polygon.vertices.size();
            const std::size_t nextIndex = (i + 1) % polygon.vertices.size();
            if (removableCollinearVertex(polygon.vertices[prev], polygon.vertices[i],
                                         polygon.vertices[nextIndex], tol)) {
                polygon.vertices.erase(polygon.vertices.begin() + static_cast<std::ptrdiff_t>(i));
                simplified = true;
                break;
            }
        }
    }

    const double area = polygon.area();
    const AABB2D bounds = polygon.bounds();
    const double span = std::max(bounds.max.x - bounds.min.x,
                                 bounds.max.y - bounds.min.y);
    const double areaEps = tol.absolute * tol.absolute +
                           tol.relative * span * span;
    if (!(area > areaEps)) {
        issue = {AgglomerationIssueCode2D::DegenerateMergedPolygon, start,
                 "agglomerated polygon has zero or near-zero area"};
        return std::nullopt;
    }
    return polygon;
}

} // namespace

AgglomerationResult2D agglomerateSmallCells(
    const std::vector<CutCell2D>& cutCells,
    const TopologyMesh2D& topology,
    const SmallCellReport2D& analysis,
    const Domain2D& domain,
    const BoundaryLoop& boundary,
    const TolerancePolicy& tol) {
    return agglomerateSmallCells(cutCells,topology,analysis,domain,
                                 BoundaryRegion2D(boundary),tol);
}

AgglomerationResult2D agglomerateSmallCells(
    const std::vector<CutCell2D>& cutCells,
    const TopologyMesh2D& topology,
    const SmallCellReport2D& analysis,
    const Domain2D& domain,
    const BoundaryRegion2D& boundary,
    const TolerancePolicy& tol) {

    (void)cutCells;
    AgglomerationResult2D result;
    result.inputCellCount = topology.cells.size();
    for (const auto& cell : topology.cells) result.totalAreaBefore += cell.geometryArea;

    if (!topology.valid() || !analysis.valid() || analysis.records.size() != topology.cells.size()) {
        result.issues.push_back({AgglomerationIssueCode2D::InvalidInput, 0,
                                 "agglomeration requires a valid topology and complete small-cell analysis"});
        return result;
    }
    if (!domain.valid(tol) || !boundary.diagnose(tol).valid()) {
        result.issues.push_back({AgglomerationIssueCode2D::InvalidInput, 0,
                                 "agglomeration requires a valid domain and boundary"});
        return result;
    }

    Dsu2D dsu(topology.cells.size());
    std::vector<bool> isSmallCell(topology.cells.size(), false);
    std::vector<double> cellAlpha(topology.cells.size(), 1.0);
    std::vector<double> cellArea(topology.cells.size(), 0.0);
    std::vector<double> cellBgArea(topology.cells.size(), 0.0);
    for (const auto& record : analysis.records) {
        if (record.topologyCellId >= topology.cells.size()) {
            result.issues.push_back({AgglomerationIssueCode2D::InvalidTopologyReference,
                                     record.topologyCellId,
                                     "small-cell record references an invalid topology cell"});
            return result;
        }
        const std::size_t cid = record.topologyCellId;
        cellAlpha[cid] = record.areaFraction;
        cellArea[cid] = record.area;
        cellBgArea[cid] = record.areaFraction > 0.0 ? record.area / record.areaFraction : record.area;
        isSmallCell[cid] = (record.status != SmallCellStatus2D::Stable);
        if (record.status == SmallCellStatus2D::CandidateFound) {
            ++result.mergedSmallCellCount;
            if (!record.targetTopologyCellId) {
                result.issues.push_back({AgglomerationIssueCode2D::MissingCandidate,
                                         cid,
                                         "small cell is marked CandidateFound without a target"});
                return result;
            }
            if (*record.targetTopologyCellId >= topology.cells.size() ||
                *record.targetTopologyCellId == cid) {
                result.issues.push_back({AgglomerationIssueCode2D::InvalidTopologyReference,
                                         cid,
                                         "small-cell candidate target is invalid"});
                return result;
            }
            if (record.targetIsSmall && !analysis.policy.allowBoundedSmallChaining) {
                result.issues.push_back({AgglomerationIssueCode2D::SmallTargetUnsupported,
                                         cid,
                                         "2D-5B does not silently chain small-to-small agglomeration"});
                return result;
            }
            dsu.unite(cid, *record.targetTopologyCellId);
        } else if (record.status == SmallCellStatus2D::Unresolved) {
            result.issues.push_back({AgglomerationIssueCode2D::MissingCandidate,
                                     cid,
                                     "unresolved small cell cannot be safely agglomerated"});
            return result;
        }
    }

    if (analysis.policy.allowBoundedSmallChaining) {
        constexpr std::size_t maxChainedGroupMembers = 8U;
        std::vector<std::map<std::size_t, double>> adj(topology.cells.size());
        for (const auto& edge : topology.edges) {
            if (!edge.neighbour) continue;
            const std::size_t a = edge.owner;
            const std::size_t b = *edge.neighbour;
            if (a >= topology.cells.size() || b >= topology.cells.size() || a == b) continue;
            const auto& pa = topology.vertices[edge.v0].point;
            const auto& pb = topology.vertices[edge.v1].point;
            const double len = std::hypot(pb.x - pa.x, pb.y - pa.y);
            adj[a][b] += len;
            adj[b][a] += len;
        }
        for (std::size_t hop = 0U; hop < 3U; ++hop) {
            std::map<std::size_t, std::vector<std::size_t>> curGroups;
            for (std::size_t cid = 0U; cid < topology.cells.size(); ++cid) {
                curGroups[dsu.find(cid)].push_back(cid);
            }
            bool anyMerged = false;
            for (const auto& [root, members] : curGroups) {
                bool hasSmall = false;
                bool hasStable = false;
                double sumArea = 0.0;
                double maxBg = 0.0;
                for (const auto m : members) {
                    if (isSmallCell[m]) hasSmall = true;
                    else hasStable = true;
                    sumArea += cellArea[m];
                    maxBg = std::max(maxBg, cellBgArea[m]);
                }
                if (!hasSmall || hasStable) continue;
                if (maxBg > 0.0 && sumArea >= analysis.policy.areaFractionThreshold * maxBg) continue;
                std::map<std::size_t, double> extShared;
                for (const auto m : members) {
                    for (const auto& [nbr, len] : adj[m]) {
                        if (dsu.find(nbr) == root) continue;
                        extShared[nbr] += len;
                    }
                }
                std::optional<std::size_t> bestNbr;
                bool bestStable = false;
                double bestLen = -1.0;
                double bestAlpha = -1.0;
                double bestArea = -1.0;
                for (const auto& [nbr, len] : extShared) {
                    const bool nbrStable = !isSmallCell[nbr];
                    const std::size_t nbrRoot = dsu.find(nbr);
                    const std::size_t combinedSize = members.size() + curGroups[nbrRoot].size();
                    if (combinedSize > maxChainedGroupMembers) continue;
                    bool better = !bestNbr.has_value();
                    if (!better && nbrStable != bestStable) better = nbrStable;
                    else if (!better && len > bestLen + tol.scale(std::max({1.0, len, bestLen}))) better = true;
                    else if (!better && !(bestLen > len + tol.scale(std::max({1.0, len, bestLen})))) {
                        if (cellAlpha[nbr] > bestAlpha + tol.scale(1.0)) better = true;
                        else if (!(bestAlpha > cellAlpha[nbr] + tol.scale(1.0)) &&
                                 cellArea[nbr] > bestArea + tol.scale(std::max({1.0, cellArea[nbr], bestArea}))) better = true;
                    }
                    if (better) {
                        bestNbr = nbr;
                        bestStable = nbrStable;
                        bestLen = len;
                        bestAlpha = cellAlpha[nbr];
                        bestArea = cellArea[nbr];
                    }
                }
                if (bestNbr) {
                    dsu.unite(root, *bestNbr);
                    anyMerged = true;
                }
            }
            if (!anyMerged) break;
        }
        std::map<std::size_t, std::vector<std::size_t>> checkGroups;
        for (std::size_t cid = 0U; cid < topology.cells.size(); ++cid) {
            checkGroups[dsu.find(cid)].push_back(cid);
        }
        for (const auto& [root, members] : checkGroups) {
            (void)root;
            bool hasSmall = false;
            bool hasStable = false;
            double sumArea = 0.0;
            double maxBg = 0.0;
            for (const auto m : members) {
                if (isSmallCell[m]) hasSmall = true;
                else hasStable = true;
                sumArea += cellArea[m];
                maxBg = std::max(maxBg, cellBgArea[m]);
            }
            if (!hasSmall) continue;
            if (members.size() > maxChainedGroupMembers ||
                (!hasStable && !(maxBg > 0.0 && sumArea >= analysis.policy.areaFractionThreshold * maxBg))) {
                result.issues.push_back({AgglomerationIssueCode2D::SmallTargetUnsupported,
                                         members.front(),
                                         "bounded small-cell chaining could not reach a stable target"});
                return result;
            }
        }
    }

    std::map<std::size_t, std::vector<std::size_t>> rawGroups;
    for (std::size_t cellId = 0; cellId < topology.cells.size(); ++cellId) {
        rawGroups[dsu.find(cellId)].push_back(cellId);
    }

    std::vector<std::vector<std::size_t>> groups;
    groups.reserve(rawGroups.size());
    for (auto& [root, members] : rawGroups) {
        (void)root;
        std::sort(members.begin(), members.end());
        groups.push_back(std::move(members));
    }
    std::sort(groups.begin(), groups.end(), [](const auto& a, const auto& b) {
        return a.front() < b.front();
    });

    std::vector<CutCell2D> synthetic;
    synthetic.reserve(groups.size());
    result.cells.reserve(groups.size());
    std::map<std::pair<double, double>, std::size_t> topologyVertexByPoint;
    if (topology.constructionRegistry) {
        for (const auto& vertex : topology.vertices) {
            topologyVertexByPoint.emplace(
                std::pair{vertex.point.x, vertex.point.y}, vertex.id);
        }
    }

    for (std::size_t groupId = 0; groupId < groups.size(); ++groupId) {
        AgglomerationIssue2D issue{AgglomerationIssueCode2D::InvalidInput, groupId, {}};
        auto polygon = polygonFromTopologyGroup(groups[groupId], topology, tol, issue);
        if (!polygon) {
            result.issues.push_back(std::move(issue));
            return result;
        }

        double sourceArea = 0.0;
        std::uint64_t minSourceKey = topology.cells[groups[groupId].front()].sourceKey;
        AgglomeratedCell2D merged;
        merged.id = groupId;
        merged.memberTopologyCellIds = groups[groupId];
        for (const std::size_t member : groups[groupId]) {
            const auto& source = topology.cells[member];
            sourceArea += source.geometryArea;
            minSourceKey = std::min(minSourceKey, source.sourceKey);
            merged.memberSourceIds.push_back(source.sourceId);
        }
        merged.polygon = *polygon;
        merged.area = merged.polygon.area();
        merged.centroid = merged.polygon.centroid();
        if (!nearlySameArea(sourceArea, merged.area, tol)) {
            result.issues.push_back({AgglomerationIssueCode2D::AreaMismatch, groupId,
                                     "agglomerated polygon area differs from the sum of member areas"});
            return result;
        }

        result.totalAreaAfter += merged.area;
        result.cells.push_back(merged);

        CutCell2D adapter;
        adapter.sourceId = groupId;
        adapter.sourceKey = minSourceKey;
        adapter.backgroundBounds = polygon->bounds();
        adapter.kind = CutCellKind::Cut;
        adapter.fluidPolygon = *polygon;
        adapter.area = merged.area;
        adapter.areaFraction = 1.0;
        adapter.centroid = merged.centroid;
        // Agglomeration rebuilds polygons from the already classified source
        // topology. Preserve that semantic boundary information: a W1
        // grid-corner weld can move a wall intersection slightly off the input
        // polyline, so re-inferring "wall" from coordinates alone is invalid.
        std::set<std::size_t> embeddedEdges;
        for (const std::size_t member : groups[groupId]) {
            for (const std::size_t edgeId : topology.cells[member].edges) {
                if (topology.edges[edgeId].patch == BoundaryPatch2D::EmbeddedBoundary) {
                    embeddedEdges.insert(edgeId);
                }
            }
        }
        for (const std::size_t edgeId : embeddedEdges) {
            const auto& edge = topology.edges[edgeId];
            adapter.embeddedBoundary.push_back({
                topology.vertices[edge.v0].point,
                topology.vertices[edge.v1].point});
        }
        if (topology.constructionRegistry) {
            adapter.canonicalRegistry = topology.constructionRegistry.get();
            for (const auto& point : adapter.fluidPolygon.vertices) {
                const auto vertex = topologyVertexByPoint.find({point.x, point.y});
                if (vertex == topologyVertexByPoint.end() ||
                    vertex->second >= topology.canonicalVertexIds.size()) {
                    result.issues.push_back({
                        AgglomerationIssueCode2D::InvalidTopologyReference, groupId,
                        "agglomerated polygon lost its canonical construction handle"});
                    return result;
                }
                adapter.canonicalVertexIds.push_back(
                    topology.canonicalVertexIds[vertex->second]);
            }
        }
        synthetic.push_back(std::move(adapter));
    }

    result.outputCellCount = result.cells.size();
    result.areaError = std::abs(result.totalAreaAfter - result.totalAreaBefore);
    if (!nearlySameArea(result.totalAreaBefore, result.totalAreaAfter, tol)) {
        result.issues.push_back({AgglomerationIssueCode2D::AreaMismatch, 0,
                                 "total fluid area is not conserved by agglomeration"});
        return result;
    }

    result.topology = buildGlobalTopology(
        synthetic, domain, boundary, tol, topology.constructionRegistry,
        topology.embeddedBoundaryFragments);
    if (!result.topology.valid()) {
        result.issues.push_back({AgglomerationIssueCode2D::RebuiltTopologyInvalid, 0,
                                 "rebuilt topology after agglomeration failed the Stage 2D-4 audit"});
        return result;
    }
    return result;
}

} // namespace cartmesh2d
