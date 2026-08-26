#include "cartmesh2d/boundary_layer/BoundaryLayer2D.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <tuple>
#include <utility>

namespace cartmesh2d {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

[[nodiscard]] bool finitePoint(const Point2D& point) noexcept {
    return std::isfinite(point.x) && std::isfinite(point.y);
}

[[nodiscard]] double norm(const Vector2D& vector) noexcept {
    return std::sqrt(squaredNorm(vector));
}

[[nodiscard]] std::optional<Vector2D> unit(const Vector2D& vector,
                                           double epsilon) noexcept {
    const double length = norm(vector);
    if (!std::isfinite(length) || length <= epsilon) return std::nullopt;
    return Vector2D{vector.x / length, vector.y / length};
}

[[nodiscard]] Vector2D fluidNormal(const Vector2D& tangent,
                                   FluidSide2D side) noexcept {
    if (side == FluidSide2D::Right) return {tangent.y, -tangent.x};
    return {-tangent.y, tangent.x};
}

[[nodiscard]] bool pointLess(const Point2D& lhs, const Point2D& rhs) noexcept {
    return std::tie(lhs.x, lhs.y) < std::tie(rhs.x, rhs.y);
}

void canonicalRotate(std::vector<Point2D>& vertices) {
    const auto first = std::min_element(vertices.begin(), vertices.end(), pointLess);
    std::rotate(vertices.begin(), first, vertices.end());
}

[[nodiscard]] std::vector<Segment2D> makeSegments(
    const std::vector<Point2D>& vertices, bool closed) {
    std::vector<Segment2D> segments;
    if (vertices.size() < 2U) return segments;
    const std::size_t count = closed ? vertices.size() : vertices.size() - 1U;
    segments.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        segments.push_back({vertices[i], vertices[(i + 1U) % vertices.size()]});
    }
    return segments;
}

[[nodiscard]] double geometryScale(const WallChain2D& chain) noexcept {
    if (chain.vertices.empty()) return 1.0;
    double minX = chain.vertices.front().x;
    double maxX = minX;
    double minY = chain.vertices.front().y;
    double maxY = minY;
    for (const auto& point : chain.vertices) {
        minX = std::min(minX, point.x);
        maxX = std::max(maxX, point.x);
        minY = std::min(minY, point.y);
        maxY = std::max(maxY, point.y);
    }
    return std::max({maxX - minX, maxY - minY, 1.0});
}

[[nodiscard]] bool adjacentSegments(std::size_t lhs, std::size_t rhs,
                                    std::size_t count, bool closed) noexcept {
    if (lhs == rhs) return true;
    if (lhs + 1U == rhs || rhs + 1U == lhs) return true;
    return closed && ((lhs == 0U && rhs + 1U == count) ||
                      (rhs == 0U && lhs + 1U == count));
}

[[nodiscard]] bool incidentSegment(std::size_t vertexId, std::size_t segmentId,
                                   const WallChain2D& chain) noexcept {
    if (segmentId == vertexId && segmentId < chain.segments.size()) return true;
    if (vertexId > 0U && segmentId + 1U == vertexId) return true;
    return chain.closed && vertexId == 0U &&
           segmentId + 1U == chain.segments.size();
}

[[nodiscard]] double pointSegmentDistance(const Point2D& point,
                                           const Segment2D& segment) noexcept {
    const Vector2D direction = segment.b - segment.a;
    const double denominator = squaredNorm(direction);
    if (denominator <= 0.0) return norm(point - segment.a);
    const double parameter = std::clamp(
        dot(point - segment.a, direction) / denominator, 0.0, 1.0);
    return norm(point - (segment.a + direction * parameter));
}

[[nodiscard]] double segmentDistance(const Segment2D& lhs,
                                     const Segment2D& rhs,
                                     const TolerancePolicy& tol) noexcept {
    if (intersectSegments(lhs, rhs, tol).kind != SegmentIntersectionKind::None) {
        return 0.0;
    }
    return std::min({pointSegmentDistance(lhs.a, rhs),
                     pointSegmentDistance(lhs.b, rhs),
                     pointSegmentDistance(rhs.a, lhs),
                     pointSegmentDistance(rhs.b, lhs)});
}

[[nodiscard]] std::optional<double> raySegmentDistance(
    const Point2D& origin, const Vector2D& direction,
    const Segment2D& segment, double epsilon) noexcept {
    const Vector2D edge = segment.b - segment.a;
    const double denominator = cross(direction, edge);
    if (std::abs(denominator) <= epsilon) return std::nullopt;
    const Vector2D offset = segment.a - origin;
    const double rayParameter = cross(offset, edge) / denominator;
    const double edgeParameter = cross(offset, direction) / denominator;
    if (rayParameter <= epsilon || edgeParameter < -epsilon ||
        edgeParameter > 1.0 + epsilon) return std::nullopt;
    return rayParameter;
}

[[nodiscard]] bool policyValid(const BoundaryLayerPolicy2D& policy) noexcept {
    return std::isfinite(policy.tolerance.absolute) &&
           std::isfinite(policy.tolerance.relative) &&
           policy.tolerance.absolute >= 0.0 &&
           policy.tolerance.relative >= 0.0 &&
           std::isfinite(policy.smoothTurnRadians) &&
           std::isfinite(policy.maxConvexTurnRadians) &&
           policy.smoothTurnRadians >= 0.0 &&
           policy.maxConvexTurnRadians > policy.smoothTurnRadians &&
           policy.maxConvexTurnRadians < kPi &&
           std::isfinite(policy.cornerLengthFraction) &&
           policy.cornerLengthFraction > 0.0 &&
           policy.cornerLengthFraction < 0.5 &&
           std::isfinite(policy.collisionClearanceFraction) &&
           policy.collisionClearanceFraction > 0.0 &&
           policy.collisionClearanceFraction < 0.5;
}

[[nodiscard]] BoundaryLayerBuildResult2D failedResult(
    BoundaryLayerFailureReason2D reason, std::string message,
    std::size_t chainId, const ResolvedLayerParameters2D* parameters,
    std::optional<std::size_t> vertexId = std::nullopt,
    std::optional<std::size_t> edgeId = std::nullopt,
    std::optional<std::size_t> cellId = std::nullopt,
    std::optional<double> safeThickness = std::nullopt) {
    BoundaryLayerBuildResult2D result;
    result.failure.reason = reason;
    result.failure.message = std::move(message);
    result.failure.chainId = chainId;
    result.failure.vertexId = vertexId;
    result.failure.edgeId = edgeId;
    result.failure.cellId = cellId;
    result.failure.safeThickness = safeThickness;
    if (parameters != nullptr) {
        result.failure.requestedThickness = parameters->totalThickness;
        result.failure.nLayers = parameters->nLayers;
        result.failure.growthRatio = parameters->growthRatio;
    }
    return result;
}

struct ChainGeometry2D {
    std::vector<WallVertexKind2D> kinds;
    std::vector<Vector2D> directions;
    std::vector<double> miterCosines;
    double safeThickness = std::numeric_limits<double>::infinity();
};

struct ChainGeometryResult2D {
    std::optional<ChainGeometry2D> geometry;
    BoundaryLayerBuildResult2D failure;
};

[[nodiscard]] ChainGeometryResult2D analyseChain(
    const WallChain2D& chain,
    const std::vector<WallChain2D>& allChains,
    const ResolvedLayerParameters2D& parameters,
    const BoundaryLayerPolicy2D& policy) {
    ChainGeometryResult2D result;
    const double scale = geometryScale(chain);
    const double epsilon = policy.tolerance.scale(scale);
    if (chain.patchIdentity.empty() || chain.vertices.size() < (chain.closed ? 3U : 2U) ||
        chain.segments.size() != (chain.closed ? chain.vertices.size()
                                               : chain.vertices.size() - 1U)) {
        result.failure = failedResult(BoundaryLayerFailureReason2D::InvalidChain,
                                      "wall chain structure is invalid", chain.id,
                                      &parameters);
        return result;
    }

    ChainGeometry2D geometry;
    geometry.kinds.resize(chain.vertices.size(), WallVertexKind2D::Smooth);
    geometry.directions.resize(chain.vertices.size());
    geometry.miterCosines.resize(chain.vertices.size(), 1.0);

    std::vector<Vector2D> tangents;
    std::vector<double> lengths;
    tangents.reserve(chain.segments.size());
    lengths.reserve(chain.segments.size());
    for (std::size_t segmentId = 0; segmentId < chain.segments.size(); ++segmentId) {
        const double length = norm(chain.segments[segmentId].b -
                                   chain.segments[segmentId].a);
        if (!std::isfinite(length) || length <= epsilon) {
            result.failure = failedResult(
                BoundaryLayerFailureReason2D::DegenerateSegment,
                "wall chain contains a zero or tolerance-scale segment", chain.id,
                &parameters, std::nullopt, segmentId);
            return result;
        }
        lengths.push_back(length);
        tangents.push_back(*(unit(chain.segments[segmentId].b -
                                  chain.segments[segmentId].a, epsilon)));
    }

    for (std::size_t vertexId = 0; vertexId < chain.vertices.size(); ++vertexId) {
        const bool endpoint = !chain.closed &&
            (vertexId == 0U || vertexId + 1U == chain.vertices.size());
        if (endpoint) {
            const std::size_t segmentId = vertexId == 0U ? 0U : tangents.size() - 1U;
            geometry.directions[vertexId] = fluidNormal(tangents[segmentId],
                                                        chain.fluidSide);
            geometry.kinds[vertexId] = WallVertexKind2D::Smooth;
            continue;
        }

        const std::size_t previous = vertexId == 0U ? tangents.size() - 1U
                                                     : vertexId - 1U;
        const std::size_t next = vertexId;
        const Vector2D previousNormal = fluidNormal(tangents[previous], chain.fluidSide);
        const Vector2D nextNormal = fluidNormal(tangents[next], chain.fluidSide);
        const double rawTurn = std::atan2(cross(tangents[previous], tangents[next]),
                                          dot(tangents[previous], tangents[next]));
        const double solidTurn = chain.fluidSide == FluidSide2D::Right
                               ? rawTurn : -rawTurn;
        if (solidTurn < -policy.smoothTurnRadians) {
            geometry.kinds[vertexId] = WallVertexKind2D::Concave;
            result.failure = failedResult(
                BoundaryLayerFailureReason2D::ConcaveCorner,
                "solid-side concave corner is outside H4-1 capability", chain.id,
                &parameters, vertexId);
            return result;
        }
        if (solidTurn > policy.maxConvexTurnRadians) {
            geometry.kinds[vertexId] = WallVertexKind2D::Sharp;
            result.failure = failedResult(
                BoundaryLayerFailureReason2D::SharpCorner,
                "wall turn exceeds the H4-1 sharp-corner limit", chain.id,
                &parameters, vertexId);
            return result;
        }

        const Vector2D sum{previousNormal.x + nextNormal.x,
                           previousNormal.y + nextNormal.y};
        const auto marching = unit(sum, epsilon);
        if (!marching) {
            geometry.kinds[vertexId] = WallVertexKind2D::Degenerate;
            result.failure = failedResult(
                BoundaryLayerFailureReason2D::ConflictingHalfPlanes,
                "adjacent fluid-side half-planes have no stable marching direction",
                chain.id, &parameters, vertexId);
            return result;
        }
        const double previousCosine = dot(*marching, previousNormal);
        const double nextCosine = dot(*marching, nextNormal);
        const double miterCosine = std::min(previousCosine, nextCosine);
        const double cosineEpsilon = std::max(1.0e-8, epsilon / scale);
        if (!std::isfinite(miterCosine) || miterCosine <= cosineEpsilon) {
            geometry.kinds[vertexId] = WallVertexKind2D::Degenerate;
            result.failure = failedResult(
                BoundaryLayerFailureReason2D::ConflictingHalfPlanes,
                "marching direction violates a fluid-side half-plane", chain.id,
                &parameters, vertexId);
            return result;
        }
        geometry.directions[vertexId] = *marching;
        geometry.miterCosines[vertexId] = miterCosine;
        geometry.kinds[vertexId] = std::abs(solidTurn) <= policy.smoothTurnRadians
                                 ? WallVertexKind2D::Smooth
                                 : WallVertexKind2D::MildConvex;

        const double tangentialAmplification = std::abs(std::tan(0.5 * solidTurn));
        if (tangentialAmplification > cosineEpsilon) {
            const double localLimit = policy.cornerLengthFraction *
                std::min(lengths[previous], lengths[next]) / tangentialAmplification;
            geometry.safeThickness = std::min(geometry.safeThickness, localLimit);
        }
    }

    // A ray along each hair direction must not reach a non-incident wall before
    // the requested layer. This detects narrow gaps without using a local-only
    // visual heuristic.
    for (std::size_t vertexId = 0; vertexId < chain.vertices.size(); ++vertexId) {
        for (const auto& otherChain : allChains) {
            for (std::size_t segmentId = 0; segmentId < otherChain.segments.size();
                 ++segmentId) {
                if (otherChain.id == chain.id &&
                    incidentSegment(vertexId, segmentId, chain)) continue;
                const auto distance = raySegmentDistance(
                    chain.vertices[vertexId], geometry.directions[vertexId],
                    otherChain.segments[segmentId], epsilon);
                if (!distance) continue;
                const double normalLimit = policy.collisionClearanceFraction *
                                           (*distance) *
                                           geometry.miterCosines[vertexId];
                geometry.safeThickness = std::min(geometry.safeThickness,
                                                   normalLimit);
            }
        }
    }

    // For distinct chains use a symmetric gap bound. Two strips each retain
    // less than half of the measured wall-to-wall gap.
    for (const auto& otherChain : allChains) {
        if (otherChain.id == chain.id) continue;
        double gap = std::numeric_limits<double>::infinity();
        for (const auto& lhs : chain.segments) {
            for (const auto& rhs : otherChain.segments) {
                gap = std::min(gap, segmentDistance(lhs, rhs, policy.tolerance));
            }
        }
        geometry.safeThickness = std::min(
            geometry.safeThickness, policy.collisionClearanceFraction * gap);
    }

    if (parameters.totalThickness > geometry.safeThickness + epsilon) {
        result.failure = failedResult(
            BoundaryLayerFailureReason2D::ThicknessExceedsSafeLimit,
            "requested total thickness exceeds the computed geometric safety limit",
            chain.id, &parameters, std::nullopt, std::nullopt, std::nullopt,
            geometry.safeThickness);
        return result;
    }
    result.geometry = std::move(geometry);
    return result;
}

[[nodiscard]] bool intersectionOnlyAtSharedVertex(
    const SegmentIntersection& intersection,
    std::size_t a0, std::size_t a1, std::size_t b0, std::size_t b1,
    const std::vector<BoundaryLayerVertex2D>& vertices,
    const TolerancePolicy& tol) {
    if (intersection.kind != SegmentIntersectionKind::Point || !intersection.point) {
        return false;
    }
    std::optional<std::size_t> shared;
    if (a0 == b0 || a0 == b1) shared = a0;
    if (a1 == b0 || a1 == b1) shared = a1;
    return shared && nearlyEqual(*intersection.point, vertices[*shared].point, tol);
}

[[nodiscard]] AABB2D cellBounds(const BoundaryLayerCell2D& cell,
                                const std::vector<BoundaryLayerVertex2D>& vertices) {
    Point2D lo = vertices[cell.vertices[0]].point;
    Point2D hi = lo;
    for (const auto vertexId : cell.vertices) {
        const auto& point = vertices[vertexId].point;
        lo.x = std::min(lo.x, point.x);
        lo.y = std::min(lo.y, point.y);
        hi.x = std::max(hi.x, point.x);
        hi.y = std::max(hi.y, point.y);
    }
    return {lo, hi};
}

[[nodiscard]] bool boxesOverlap(const AABB2D& lhs, const AABB2D& rhs,
                                double epsilon) noexcept {
    return lhs.max.x >= rhs.min.x - epsilon && rhs.max.x >= lhs.min.x - epsilon &&
           lhs.max.y >= rhs.min.y - epsilon && rhs.max.y >= lhs.min.y - epsilon;
}

[[nodiscard]] BoundaryLayerBuildResult2D validateStrip(
    BoundaryLayerStrip2D strip, const BoundaryLayerPolicy2D& policy) {
    const double scale = geometryScale(strip.wallChain);
    const double lengthEpsilon = policy.tolerance.scale(scale);
    const double areaEpsilon = std::max(policy.tolerance.absolute *
                                        policy.tolerance.absolute,
                                        policy.tolerance.relative * scale * scale);
    const std::size_t ringSize = strip.wallChain.vertices.size();
    const std::size_t ringCount = strip.ringVertexIds.size();

    for (std::size_t lhs = 0; lhs < strip.vertices.size(); ++lhs) {
        for (std::size_t rhs = lhs + 1U; rhs < strip.vertices.size(); ++rhs) {
            if (nearlyEqual(strip.vertices[lhs].point, strip.vertices[rhs].point,
                            policy.tolerance)) {
                return failedResult(
                    BoundaryLayerFailureReason2D::DuplicateGeometricVertex,
                    "distinct layer vertex IDs occupy the same geometric point",
                    strip.wallChain.id, &strip.parameters, strip.vertices[rhs].chainVertex);
            }
        }
    }

    // Every ring must be simple. Open rings are checked for non-adjacent edge
    // intersections; closed rings additionally check the closing edge.
    for (std::size_t ring = 0; ring < ringCount; ++ring) {
        const auto& ids = strip.ringVertexIds[ring];
        const std::size_t edgeCount = strip.wallChain.closed ? ids.size() : ids.size() - 1U;
        for (std::size_t i = 0; i < edgeCount; ++i) {
            const Segment2D first{strip.vertices[ids[i]].point,
                                  strip.vertices[ids[(i + 1U) % ids.size()]].point};
            for (std::size_t j = i + 1U; j < edgeCount; ++j) {
                if (adjacentSegments(i, j, edgeCount, strip.wallChain.closed)) continue;
                const Segment2D second{strip.vertices[ids[j]].point,
                                       strip.vertices[ids[(j + 1U) % ids.size()]].point};
                if (intersectSegments(first, second, policy.tolerance).kind !=
                    SegmentIntersectionKind::None) {
                    return failedResult(
                        BoundaryLayerFailureReason2D::EnvelopeSelfIntersection,
                        "a layer ring has a non-adjacent self-intersection",
                        strip.wallChain.id, &strip.parameters, std::nullopt, i);
                }
            }
        }
    }

    // Outer envelope and every non-wall ring must not intersect the original
    // wall. Hair edges may touch only their own wall vertex.
    for (std::size_t ring = 1U; ring < ringCount; ++ring) {
        const auto& ids = strip.ringVertexIds[ring];
        const std::size_t edgeCount = strip.wallChain.closed ? ids.size() : ids.size() - 1U;
        for (std::size_t edgeId = 0; edgeId < edgeCount; ++edgeId) {
            const Segment2D layerEdge{strip.vertices[ids[edgeId]].point,
                strip.vertices[ids[(edgeId + 1U) % ids.size()]].point};
            for (std::size_t wallEdge = 0; wallEdge < strip.wallChain.segments.size();
                 ++wallEdge) {
                if (intersectSegments(layerEdge, strip.wallChain.segments[wallEdge],
                                      policy.tolerance).kind !=
                    SegmentIntersectionKind::None) {
                    return failedResult(
                        BoundaryLayerFailureReason2D::EnvelopeWallIntersection,
                        "a layer ring intersects the original wall", strip.wallChain.id,
                        &strip.parameters, std::nullopt, edgeId);
                }
            }
        }
    }
    for (std::size_t vertexId = 0; vertexId < ringSize; ++vertexId) {
        for (std::size_t ring = 0; ring + 1U < ringCount; ++ring) {
            const Segment2D hair{
                strip.vertices[strip.ringVertexIds[ring][vertexId]].point,
                strip.vertices[strip.ringVertexIds[ring + 1U][vertexId]].point};
            for (std::size_t wallEdge = 0; wallEdge < strip.wallChain.segments.size();
                 ++wallEdge) {
                if (incidentSegment(vertexId, wallEdge, strip.wallChain)) continue;
                if (intersectSegments(hair, strip.wallChain.segments[wallEdge],
                                      policy.tolerance).kind !=
                    SegmentIntersectionKind::None) {
                    return failedResult(
                        BoundaryLayerFailureReason2D::EnvelopeWallIntersection,
                        "a hair edge crosses a non-incident wall segment",
                        strip.wallChain.id, &strip.parameters, vertexId, wallEdge);
                }
            }
        }
    }

    std::map<std::pair<std::size_t, std::size_t>, std::size_t> edgeIncidence;
    double minArea = std::numeric_limits<double>::infinity();
    double maxArea = 0.0;
    for (auto& cell : strip.cells) {
        Polygon2D polygon;
        for (const auto vertexId : cell.vertices) {
            polygon.vertices.push_back(strip.vertices[vertexId].point);
        }
        const double signedArea = polygon.signedArea();
        if (!std::isfinite(signedArea) || signedArea <= areaEpsilon) {
            return failedResult(BoundaryLayerFailureReason2D::NegativeAreaCell,
                                "layer quad has non-positive or tolerance-scale area",
                                strip.wallChain.id, &strip.parameters,
                                std::nullopt, std::nullopt, cell.id);
        }
        const Segment2D sides[4]{
            {polygon.vertices[0], polygon.vertices[1]},
            {polygon.vertices[1], polygon.vertices[2]},
            {polygon.vertices[2], polygon.vertices[3]},
            {polygon.vertices[3], polygon.vertices[0]}};
        if (intersectSegments(sides[0], sides[2], policy.tolerance).kind !=
                SegmentIntersectionKind::None ||
            intersectSegments(sides[1], sides[3], policy.tolerance).kind !=
                SegmentIntersectionKind::None) {
            return failedResult(BoundaryLayerFailureReason2D::SelfIntersectingCell,
                                "layer quad is self-intersecting", strip.wallChain.id,
                                &strip.parameters, std::nullopt, std::nullopt, cell.id);
        }
        cell.area = signedArea;
        minArea = std::min(minArea, signedArea);
        maxArea = std::max(maxArea, signedArea);
        for (std::size_t side = 0; side < 4U; ++side) {
            const std::size_t a = cell.vertices[side];
            const std::size_t b = cell.vertices[(side + 1U) % 4U];
            ++edgeIncidence[std::minmax(a, b)];
        }
    }

    for (const auto& [edge, count] : edgeIncidence) {
        const auto& a = strip.vertices[edge.first];
        const auto& b = strip.vertices[edge.second];
        std::size_t expected = 0;
        if (a.ring == b.ring) {
            if (a.ring == 0U || a.ring + 1U == ringCount) expected = 1U;
            else expected = 2U;
        } else if (a.chainVertex == b.chainVertex &&
                   (a.ring + 1U == b.ring || b.ring + 1U == a.ring)) {
            const bool openEnd = !strip.wallChain.closed &&
                (a.chainVertex == 0U || a.chainVertex + 1U == ringSize);
            expected = openEnd ? 1U : 2U;
        }
        if (expected == 0U || count != expected) {
            return failedResult(BoundaryLayerFailureReason2D::TopologyConflict,
                                "layer edge incidence does not match strip topology",
                                strip.wallChain.id, &strip.parameters);
        }
    }

    // Broad-phase AABB followed by exact edge and containment checks prevents
    // visually plausible but overlapping non-neighbour quads.
    std::vector<AABB2D> bounds;
    bounds.reserve(strip.cells.size());
    for (const auto& cell : strip.cells) bounds.push_back(cellBounds(cell, strip.vertices));
    for (std::size_t lhsId = 0; lhsId < strip.cells.size(); ++lhsId) {
        const auto& lhs = strip.cells[lhsId];
        for (std::size_t rhsId = lhsId + 1U; rhsId < strip.cells.size(); ++rhsId) {
            const auto& rhs = strip.cells[rhsId];
            if (!boxesOverlap(bounds[lhsId], bounds[rhsId], lengthEpsilon)) continue;
            bool shareEdge = false;
            for (std::size_t li = 0; li < 4U; ++li) {
                const std::size_t la = lhs.vertices[li];
                const std::size_t lb = lhs.vertices[(li + 1U) % 4U];
                for (std::size_t ri = 0; ri < 4U; ++ri) {
                    const std::size_t ra = rhs.vertices[ri];
                    const std::size_t rb = rhs.vertices[(ri + 1U) % 4U];
                    if (std::minmax(la, lb) == std::minmax(ra, rb)) {
                        shareEdge = true;
                        continue;
                    }
                    const Segment2D left{strip.vertices[la].point,
                                         strip.vertices[lb].point};
                    const Segment2D right{strip.vertices[ra].point,
                                          strip.vertices[rb].point};
                    const auto intersection = intersectSegments(left, right,
                                                                policy.tolerance);
                    if (intersection.kind == SegmentIntersectionKind::None) continue;
                    if (intersectionOnlyAtSharedVertex(intersection, la, lb, ra, rb,
                                                       strip.vertices,
                                                       policy.tolerance)) continue;
                    if (std::minmax(la, lb) == std::minmax(ra, rb)) continue;
                    return failedResult(BoundaryLayerFailureReason2D::OverlappingCells,
                                        "non-shared layer cell edges intersect or overlap",
                                        strip.wallChain.id, &strip.parameters,
                                        std::nullopt, std::nullopt, rhs.id);
                }
            }
            if (shareEdge) continue;
            Polygon2D leftPolygon;
            Polygon2D rightPolygon;
            for (const auto id : lhs.vertices) leftPolygon.vertices.push_back(strip.vertices[id].point);
            for (const auto id : rhs.vertices) rightPolygon.vertices.push_back(strip.vertices[id].point);
            for (const auto id : lhs.vertices) {
                if (classifyPointInPolygon(strip.vertices[id].point, rightPolygon,
                                           policy.tolerance) == PointInPolygon::Inside) {
                    return failedResult(BoundaryLayerFailureReason2D::OverlappingCells,
                                        "a layer cell vertex lies inside another cell",
                                        strip.wallChain.id, &strip.parameters,
                                        std::nullopt, std::nullopt, rhs.id);
                }
            }
            for (const auto id : rhs.vertices) {
                if (classifyPointInPolygon(strip.vertices[id].point, leftPolygon,
                                           policy.tolerance) == PointInPolygon::Inside) {
                    return failedResult(BoundaryLayerFailureReason2D::OverlappingCells,
                                        "a layer cell vertex lies inside another cell",
                                        strip.wallChain.id, &strip.parameters,
                                        std::nullopt, std::nullopt, rhs.id);
                }
            }
        }
    }

    double minLayerThickness = std::numeric_limits<double>::infinity();
    double maxLayerThickness = 0.0;
    for (std::size_t ring = 0; ring + 1U < ringCount; ++ring) {
        for (std::size_t vertexId = 0; vertexId < ringSize; ++vertexId) {
            const double thickness = norm(
                strip.vertices[strip.ringVertexIds[ring + 1U][vertexId]].point -
                strip.vertices[strip.ringVertexIds[ring][vertexId]].point);
            if (!std::isfinite(thickness) || thickness <= lengthEpsilon) {
                return failedResult(BoundaryLayerFailureReason2D::TopologyConflict,
                                    "hair-edge layer spacing is not strictly positive",
                                    strip.wallChain.id, &strip.parameters, vertexId);
            }
            minLayerThickness = std::min(minLayerThickness, thickness);
            maxLayerThickness = std::max(maxLayerThickness, thickness);
        }
    }
    strip.metrics.vertexCount = strip.vertices.size();
    strip.metrics.cellCount = strip.cells.size();
    strip.metrics.minCellArea = minArea;
    strip.metrics.maxCellArea = maxArea;
    strip.metrics.minLayerThickness = minLayerThickness;
    strip.metrics.maxLayerThickness = maxLayerThickness;

    BoundaryLayerBuildResult2D result;
    result.status = BoundaryLayerStatus2D::Success;
    result.strips.push_back(std::move(strip));
    return result;
}

[[nodiscard]] BoundaryLayerStrip2D constructStrip(
    const WallChain2D& chain, const ResolvedLayerParameters2D& parameters,
    const ChainGeometry2D& geometry, const std::vector<double>& cumulative) {
    BoundaryLayerStrip2D strip;
    strip.wallChain = chain;
    strip.parameters = parameters;
    strip.wallVertexKinds = geometry.kinds;
    strip.marchingDirections = geometry.directions;
    const std::size_t ringSize = chain.vertices.size();
    strip.ringVertexIds.resize(cumulative.size() + 1U);
    strip.vertices.reserve((cumulative.size() + 1U) * ringSize);
    for (std::size_t ring = 0; ring <= cumulative.size(); ++ring) {
        strip.ringVertexIds[ring].reserve(ringSize);
        const double normalDistance = ring == 0U ? 0.0 : cumulative[ring - 1U];
        for (std::size_t vertexId = 0; vertexId < ringSize; ++vertexId) {
            const double hairDistance = normalDistance / geometry.miterCosines[vertexId];
            const Point2D point = chain.vertices[vertexId] +
                                  geometry.directions[vertexId] * hairDistance;
            const std::size_t id = strip.vertices.size();
            strip.ringVertexIds[ring].push_back(id);
            strip.vertices.push_back({id, ring, vertexId, point});
        }
    }

    const std::size_t segmentCount = chain.segments.size();
    strip.cells.reserve(cumulative.size() * segmentCount);
    for (std::size_t layer = 0; layer < cumulative.size(); ++layer) {
        for (std::size_t segment = 0; segment < segmentCount; ++segment) {
            const std::size_t next = (segment + 1U) % ringSize;
            BoundaryLayerCell2D cell;
            cell.id = strip.cells.size();
            cell.layer = layer;
            cell.wallSegment = segment;
            if (chain.fluidSide == FluidSide2D::Right) {
                cell.vertices = {strip.ringVertexIds[layer][segment],
                                 strip.ringVertexIds[layer + 1U][segment],
                                 strip.ringVertexIds[layer + 1U][next],
                                 strip.ringVertexIds[layer][next]};
            } else {
                cell.vertices = {strip.ringVertexIds[layer][segment],
                                 strip.ringVertexIds[layer][next],
                                 strip.ringVertexIds[layer + 1U][next],
                                 strip.ringVertexIds[layer + 1U][segment]};
            }
            strip.cells.push_back(cell);
        }
    }
    strip.metrics.requestedTotalThickness = parameters.totalThickness;
    strip.metrics.usedTotalThickness = cumulative.back();
    strip.metrics.safeThicknessLimit = geometry.safeThickness;
    return strip;
}

[[nodiscard]] bool stripsIntersect(const BoundaryLayerStrip2D& lhs,
                                   const BoundaryLayerStrip2D& rhs,
                                   const TolerancePolicy& tol) {
    for (const auto& lhsCell : lhs.cells) {
        const auto lhsBounds = cellBounds(lhsCell, lhs.vertices);
        for (const auto& rhsCell : rhs.cells) {
            const auto rhsBounds = cellBounds(rhsCell, rhs.vertices);
            const double scale = std::max({lhsBounds.max.x - lhsBounds.min.x,
                                           lhsBounds.max.y - lhsBounds.min.y,
                                           rhsBounds.max.x - rhsBounds.min.x,
                                           rhsBounds.max.y - rhsBounds.min.y, 1.0});
            if (!boxesOverlap(lhsBounds, rhsBounds, tol.scale(scale))) continue;
            for (std::size_t li = 0; li < 4U; ++li) {
                const Segment2D left{lhs.vertices[lhsCell.vertices[li]].point,
                    lhs.vertices[lhsCell.vertices[(li + 1U) % 4U]].point};
                for (std::size_t ri = 0; ri < 4U; ++ri) {
                    const Segment2D right{rhs.vertices[rhsCell.vertices[ri]].point,
                        rhs.vertices[rhsCell.vertices[(ri + 1U) % 4U]].point};
                    if (intersectSegments(left, right, tol).kind !=
                        SegmentIntersectionKind::None) return true;
                }
            }
            Polygon2D leftPolygon;
            Polygon2D rightPolygon;
            for (const auto id : lhsCell.vertices) leftPolygon.vertices.push_back(lhs.vertices[id].point);
            for (const auto id : rhsCell.vertices) rightPolygon.vertices.push_back(rhs.vertices[id].point);
            if (classifyPointInPolygon(lhs.vertices[lhsCell.vertices[0]].point,
                                       rightPolygon, tol) != PointInPolygon::Outside ||
                classifyPointInPolygon(rhs.vertices[rhsCell.vertices[0]].point,
                                       leftPolygon, tol) != PointInPolygon::Outside) return true;
        }
    }
    return false;
}

void setError(std::string* error, const std::string& message) {
    if (error != nullptr) *error = message;
}

void writeJsonString(std::ostream& out, const std::string& value) {
    out << '"';
    for (const char character : value) {
        switch (character) {
        case '"': out << "\\\""; break;
        case '\\': out << "\\\\"; break;
        case '\n': out << "\\n"; break;
        case '\r': out << "\\r"; break;
        case '\t': out << "\\t"; break;
        default: out << character; break;
        }
    }
    out << '"';
}

} // namespace

WallChainBuildResult2D makeClosedWallChain2D(
    const BoundaryLoop& boundary, std::size_t chainId,
    std::string patchIdentity, WallFluidRegion2D fluidRegion,
    const TolerancePolicy& tol) {
    WallChainBuildResult2D result;
    if (patchIdentity.empty()) {
        result.failureReason = WallChainFailureReason2D::EmptyPatchIdentity;
        result.message = "wall chain patch identity must not be empty";
        return result;
    }
    const auto diagnostics = boundary.diagnose(tol);
    if (!diagnostics.valid()) {
        result.failureReason = WallChainFailureReason2D::InvalidBoundary;
        result.message = diagnostics.issues.empty()
            ? "closed wall boundary is invalid"
            : diagnostics.issues.front().message;
        if (!diagnostics.issues.empty()) result.vertexId = diagnostics.issues.front().firstIndex;
        return result;
    }

    std::vector<Point2D> vertices = boundary.vertices();
    const bool wantCounterClockwise = fluidRegion == WallFluidRegion2D::Exterior;
    const bool isCounterClockwise =
        diagnostics.orientation == LoopOrientation::CounterClockwise;
    if (wantCounterClockwise != isCounterClockwise) std::reverse(vertices.begin(), vertices.end());
    canonicalRotate(vertices);

    WallChain2D chain;
    chain.id = chainId;
    chain.patchIdentity = std::move(patchIdentity);
    chain.vertices = std::move(vertices);
    chain.segments = makeSegments(chain.vertices, true);
    chain.closed = true;
    chain.orientation = wantCounterClockwise
        ? WallChainOrientation2D::CounterClockwise
        : WallChainOrientation2D::Clockwise;
    chain.fluidSide = FluidSide2D::Right;
    result.chain = std::move(chain);
    return result;
}

WallChainBuildResult2D makeOpenWallChain2D(
    std::vector<Point2D> orderedVertices, std::size_t chainId,
    std::string patchIdentity, FluidSide2D fluidSide,
    const TolerancePolicy& tol) {
    WallChainBuildResult2D result;
    if (patchIdentity.empty()) {
        result.failureReason = WallChainFailureReason2D::EmptyPatchIdentity;
        result.message = "wall chain patch identity must not be empty";
        return result;
    }
    if (orderedVertices.size() < 2U) {
        result.failureReason = WallChainFailureReason2D::TooFewVertices;
        result.message = "open wall chain needs at least two vertices";
        return result;
    }
    const double scale = [&]() {
        double resultScale = 1.0;
        for (const auto& point : orderedVertices) {
            resultScale = std::max({resultScale, std::abs(point.x), std::abs(point.y)});
        }
        return resultScale;
    }();
    const double epsilon = tol.scale(scale);
    for (std::size_t i = 0; i < orderedVertices.size(); ++i) {
        if (!finitePoint(orderedVertices[i])) {
            result.failureReason = WallChainFailureReason2D::NonFiniteCoordinate;
            result.vertexId = i;
            result.message = "open wall chain contains a non-finite coordinate";
            return result;
        }
        for (std::size_t j = 0; j < i; ++j) {
            if (nearlyEqual(orderedVertices[i], orderedVertices[j], tol)) {
                result.failureReason = WallChainFailureReason2D::DuplicateVertex;
                result.vertexId = i;
                result.message = "open wall chain contains a duplicate vertex";
                return result;
            }
        }
        if (i > 0U && norm(orderedVertices[i] - orderedVertices[i - 1U]) <= epsilon) {
            result.failureReason = WallChainFailureReason2D::DegenerateSegment;
            result.vertexId = i;
            result.message = "open wall chain contains a degenerate segment";
            return result;
        }
    }
    const auto segments = makeSegments(orderedVertices, false);
    for (std::size_t i = 0; i < segments.size(); ++i) {
        for (std::size_t j = i + 1U; j < segments.size(); ++j) {
            if (adjacentSegments(i, j, segments.size(), false)) continue;
            if (intersectSegments(segments[i], segments[j], tol).kind !=
                SegmentIntersectionKind::None) {
                result.failureReason = WallChainFailureReason2D::InvalidBoundary;
                result.vertexId = j;
                result.message = "open wall chain self-intersects";
                return result;
            }
        }
    }

    WallChain2D chain;
    chain.id = chainId;
    chain.patchIdentity = std::move(patchIdentity);
    chain.vertices = std::move(orderedVertices);
    chain.segments = segments;
    chain.closed = false;
    chain.orientation = WallChainOrientation2D::Open;
    chain.fluidSide = fluidSide;
    result.chain = std::move(chain);
    return result;
}

LayerParameterResult2D resolveLayerParameters2D(
    const LayerParameters2D& parameters) {
    LayerParameterResult2D result;
    if (parameters.nLayers < 1U) {
        result.message = "nLayers must be at least one";
        return result;
    }
    if (!std::isfinite(parameters.thickness) || parameters.thickness <= 0.0) {
        result.message = "layer thickness must be finite and positive";
        return result;
    }
    if (!std::isfinite(parameters.growthRatio) || parameters.growthRatio <= 0.0) {
        result.message = "growthRatio must be finite and positive";
        return result;
    }

    const long double ratio = parameters.growthRatio;
    const long double count = static_cast<long double>(parameters.nLayers);
    const long double logRatio = std::log(ratio);
    const long double nearOne = 64.0L * std::numeric_limits<double>::epsilon();
    long double sum = 0.0L;
    if (std::abs(logRatio) <= nearOne) {
        sum = count;
    } else {
        sum = std::expm1(count * logRatio) / std::expm1(logRatio);
    }
    if (!std::isfinite(sum) || sum <= 0.0L) {
        result.message = "layer geometric progression is not finite";
        return result;
    }

    const long double first = parameters.thicknessMode ==
        LayerThicknessMode2D::FirstLayerThickness
        ? static_cast<long double>(parameters.thickness)
        : static_cast<long double>(parameters.thickness) / sum;
    const long double total = first * sum;
    if (!std::isfinite(first) || !std::isfinite(total) || first <= 0.0L || total <= 0.0L ||
        first > std::numeric_limits<double>::max() ||
        total > std::numeric_limits<double>::max()) {
        result.message = "resolved layer thickness is not finite";
        return result;
    }

    ResolvedLayerParameters2D resolved;
    resolved.nLayers = parameters.nLayers;
    resolved.firstLayerThickness = static_cast<double>(first);
    resolved.totalThickness = static_cast<double>(total);
    resolved.growthRatio = parameters.growthRatio;
    resolved.cumulativeNormalDistances.reserve(parameters.nLayers);
    long double cumulative = 0.0L;
    long double current = first;
    for (std::size_t layer = 0; layer < parameters.nLayers; ++layer) {
        cumulative += current;
        if (!std::isfinite(cumulative) || cumulative > std::numeric_limits<double>::max()) {
            result.message = "layer cumulative thickness overflowed";
            return result;
        }
        resolved.cumulativeNormalDistances.push_back(static_cast<double>(cumulative));
        current *= ratio;
    }
    resolved.cumulativeNormalDistances.back() = resolved.totalThickness;
    result.parameters = std::move(resolved);
    return result;
}

std::vector<Point2D> BoundaryLayerStrip2D::outerEnvelope() const {
    std::vector<Point2D> result;
    if (ringVertexIds.empty()) return result;
    result.reserve(ringVertexIds.back().size());
    for (const auto id : ringVertexIds.back()) result.push_back(vertices[id].point);
    return result;
}

BoundaryLayerBuildResult2D buildBoundaryLayerStrips2D(
    const std::vector<WallChain2D>& wallChains,
    const LayerParameters2D& parameters,
    const BoundaryLayerPolicy2D& policy) {
    if (!policyValid(policy)) {
        return failedResult(BoundaryLayerFailureReason2D::InvalidParameters,
                            "boundary-layer tolerance/angle policy is invalid", 0U,
                            nullptr);
    }
    const auto resolvedResult = resolveLayerParameters2D(parameters);
    if (!resolvedResult.success()) {
        auto failure = failedResult(BoundaryLayerFailureReason2D::InvalidParameters,
                                    resolvedResult.message, 0U, nullptr);
        failure.failure.nLayers = parameters.nLayers;
        failure.failure.growthRatio = parameters.growthRatio;
        failure.failure.requestedThickness = parameters.thickness;
        return failure;
    }
    const auto& resolved = *resolvedResult.parameters;
    if (wallChains.empty()) {
        return failedResult(BoundaryLayerFailureReason2D::InvalidChain,
                            "at least one wall chain is required", 0U, &resolved);
    }
    std::set<std::size_t> chainIds;
    for (const auto& chain : wallChains) {
        if (!chainIds.insert(chain.id).second) {
            return failedResult(BoundaryLayerFailureReason2D::InvalidChain,
                                "wall chain IDs must be unique", chain.id, &resolved);
        }
    }

    std::vector<ChainGeometry2D> geometries;
    geometries.reserve(wallChains.size());
    for (const auto& chain : wallChains) {
        auto analysed = analyseChain(chain, wallChains, resolved, policy);
        if (!analysed.geometry) return analysed.failure;
        geometries.push_back(std::move(*analysed.geometry));
    }

    std::vector<BoundaryLayerStrip2D> accepted;
    accepted.reserve(wallChains.size());
    for (std::size_t chainId = 0; chainId < wallChains.size(); ++chainId) {
        // Validate the one-cell-thick wrapper before subdividing its hair edges.
        auto wrapper = constructStrip(wallChains[chainId], resolved,
                                      geometries[chainId],
                                      {resolved.totalThickness});
        auto wrapperValidation = validateStrip(std::move(wrapper), policy);
        if (!wrapperValidation.success()) return wrapperValidation;

        auto strip = constructStrip(wallChains[chainId], resolved,
                                    geometries[chainId],
                                    resolved.cumulativeNormalDistances);
        auto validation = validateStrip(std::move(strip), policy);
        if (!validation.success()) return validation;
        accepted.push_back(std::move(validation.strips.front()));
    }

    for (std::size_t lhs = 0; lhs < accepted.size(); ++lhs) {
        for (std::size_t rhs = lhs + 1U; rhs < accepted.size(); ++rhs) {
            if (stripsIntersect(accepted[lhs], accepted[rhs], policy.tolerance)) {
                return failedResult(BoundaryLayerFailureReason2D::ChainCollision,
                                    "candidate strips from distinct wall chains collide",
                                    accepted[rhs].wallChain.id, &resolved);
            }
        }
    }

    BoundaryLayerBuildResult2D result;
    result.status = BoundaryLayerStatus2D::Success;
    result.strips = std::move(accepted);
    return result;
}

BoundaryLayerBuildResult2D buildBoundaryLayerStrip2D(
    const WallChain2D& wallChain,
    const LayerParameters2D& parameters,
    const BoundaryLayerPolicy2D& policy) {
    return buildBoundaryLayerStrips2D({wallChain}, parameters, policy);
}

const char* boundaryLayerFailureReasonName(
    BoundaryLayerFailureReason2D reason) noexcept {
    switch (reason) {
    case BoundaryLayerFailureReason2D::None: return "none";
    case BoundaryLayerFailureReason2D::InvalidParameters: return "invalid_parameters";
    case BoundaryLayerFailureReason2D::InvalidChain: return "invalid_chain";
    case BoundaryLayerFailureReason2D::DegenerateSegment: return "degenerate_segment";
    case BoundaryLayerFailureReason2D::ConcaveCorner: return "concave_corner";
    case BoundaryLayerFailureReason2D::SharpCorner: return "sharp_corner";
    case BoundaryLayerFailureReason2D::ConflictingHalfPlanes: return "conflicting_half_planes";
    case BoundaryLayerFailureReason2D::ThicknessExceedsSafeLimit: return "thickness_exceeds_safe_limit";
    case BoundaryLayerFailureReason2D::EnvelopeSelfIntersection: return "envelope_self_intersection";
    case BoundaryLayerFailureReason2D::EnvelopeWallIntersection: return "envelope_wall_intersection";
    case BoundaryLayerFailureReason2D::ChainCollision: return "chain_collision";
    case BoundaryLayerFailureReason2D::NegativeAreaCell: return "negative_area_cell";
    case BoundaryLayerFailureReason2D::SelfIntersectingCell: return "self_intersecting_cell";
    case BoundaryLayerFailureReason2D::OverlappingCells: return "overlapping_cells";
    case BoundaryLayerFailureReason2D::DuplicateGeometricVertex: return "duplicate_geometric_vertex";
    case BoundaryLayerFailureReason2D::TopologyConflict: return "topology_conflict";
    case BoundaryLayerFailureReason2D::IoFailure: return "io_failure";
    }
    return "unknown";
}

const char* wallVertexKindName(WallVertexKind2D kind) noexcept {
    switch (kind) {
    case WallVertexKind2D::Smooth: return "smooth";
    case WallVertexKind2D::MildConvex: return "mild_convex";
    case WallVertexKind2D::Concave: return "concave";
    case WallVertexKind2D::Sharp: return "sharp";
    case WallVertexKind2D::Degenerate: return "degenerate";
    }
    return "unknown";
}

bool writeBoundaryLayerLegacyVtk2D(
    const BoundaryLayerBuildResult2D& result,
    const std::filesystem::path& path, std::string* error) {
    if (!result.success() || result.strips.empty()) {
        setError(error, "cannot write VTK for a failed or empty boundary-layer candidate");
        return false;
    }
    std::ofstream out(path);
    if (!out) {
        setError(error, "failed to open boundary-layer VTK output");
        return false;
    }
    std::size_t pointCount = 0U;
    std::size_t cellCount = 0U;
    for (const auto& strip : result.strips) {
        pointCount += strip.vertices.size();
        cellCount += strip.cells.size();
    }
    out << std::setprecision(17);
    out << "# vtk DataFile Version 3.0\n";
    out << "cartmesh2d H4-1 boundary-layer strip\nASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";
    out << "POINTS " << pointCount << " double\n";
    for (const auto& strip : result.strips) {
        for (const auto& vertex : strip.vertices) {
            out << vertex.point.x << ' ' << vertex.point.y << " 0\n";
        }
    }
    out << "CELLS " << cellCount << ' ' << cellCount * 5U << "\n";
    std::size_t pointOffset = 0U;
    for (const auto& strip : result.strips) {
        for (const auto& cell : strip.cells) {
            out << "4";
            for (const auto id : cell.vertices) out << ' ' << pointOffset + id;
            out << '\n';
        }
        pointOffset += strip.vertices.size();
    }
    out << "CELL_TYPES " << cellCount << "\n";
    for (std::size_t i = 0; i < cellCount; ++i) out << "9\n";
    out << "CELL_DATA " << cellCount << "\n";
    out << "SCALARS layer_index int 1\nLOOKUP_TABLE default\n";
    for (const auto& strip : result.strips) {
        for (const auto& cell : strip.cells) out << cell.layer << '\n';
    }
    out << "SCALARS wall_segment int 1\nLOOKUP_TABLE default\n";
    for (const auto& strip : result.strips) {
        for (const auto& cell : strip.cells) out << cell.wallSegment << '\n';
    }
    out << "SCALARS cell_area double 1\nLOOKUP_TABLE default\n";
    for (const auto& strip : result.strips) {
        for (const auto& cell : strip.cells) out << cell.area << '\n';
    }
    if (!out.good()) {
        setError(error, "failed while writing boundary-layer VTK output");
        return false;
    }
    return true;
}

bool writeBoundaryLayerReportJson2D(
    const BoundaryLayerBuildResult2D& result,
    const std::filesystem::path& path, std::string* error) {
    std::ofstream out(path);
    if (!out) {
        setError(error, "failed to open boundary-layer JSON report");
        return false;
    }
    out << std::setprecision(17);
    out << "{\n  \"layer_status\": \""
        << (result.success() ? "success" : "failed") << "\",\n";
    if (!result.success()) {
        out << "  \"failure_reason\": \""
            << boundaryLayerFailureReasonName(result.failure.reason) << "\",\n";
        out << "  \"message\": ";
        writeJsonString(out, result.failure.message);
        out << ",\n  \"chain_id\": " << result.failure.chainId << ",\n";
        out << "  \"vertex_id\": ";
        if (result.failure.vertexId) out << *result.failure.vertexId; else out << "null";
        out << ",\n  \"edge_id\": ";
        if (result.failure.edgeId) out << *result.failure.edgeId; else out << "null";
        out << ",\n  \"cell_id\": ";
        if (result.failure.cellId) out << *result.failure.cellId; else out << "null";
        out << ",\n  \"requested_thickness\": " << result.failure.requestedThickness;
        out << ",\n  \"safe_thickness\": ";
        if (result.failure.safeThickness && std::isfinite(*result.failure.safeThickness)) {
            out << *result.failure.safeThickness;
        } else out << "null";
        out << ",\n  \"n_layers\": " << result.failure.nLayers;
        out << ",\n  \"growth_ratio\": " << result.failure.growthRatio << "\n}\n";
    } else {
        std::size_t totalVertices = 0U;
        std::size_t totalCells = 0U;
        double minArea = std::numeric_limits<double>::infinity();
        double maxArea = 0.0;
        double minThickness = std::numeric_limits<double>::infinity();
        double maxThickness = 0.0;
        for (const auto& strip : result.strips) {
            totalVertices += strip.metrics.vertexCount;
            totalCells += strip.metrics.cellCount;
            minArea = std::min(minArea, strip.metrics.minCellArea);
            maxArea = std::max(maxArea, strip.metrics.maxCellArea);
            minThickness = std::min(minThickness, strip.metrics.minLayerThickness);
            maxThickness = std::max(maxThickness, strip.metrics.maxLayerThickness);
        }
        out << "  \"failure_reason\": \"none\",\n";
        out << "  \"strip_count\": " << result.strips.size() << ",\n";
        out << "  \"layer_vertex_count\": " << totalVertices << ",\n";
        out << "  \"layer_cell_count\": " << totalCells << ",\n";
        out << "  \"min_cell_area\": " << minArea << ",\n";
        out << "  \"max_cell_area\": " << maxArea << ",\n";
        out << "  \"min_layer_thickness\": " << minThickness << ",\n";
        out << "  \"max_layer_thickness\": " << maxThickness << ",\n";
        out << "  \"strips\": [\n";
        for (std::size_t i = 0; i < result.strips.size(); ++i) {
            const auto& strip = result.strips[i];
            out << "    {\"chain_id\": " << strip.wallChain.id
                << ", \"patch\": ";
            writeJsonString(out, strip.wallChain.patchIdentity);
            out << ", \"closed\": " << (strip.wallChain.closed ? "true" : "false")
                << ", \"fluid_side\": \""
                << (strip.wallChain.fluidSide == FluidSide2D::Right ? "right" : "left")
                << "\", \"wall_vertex_count\": " << strip.wallChain.vertices.size()
                << ", \"layer_vertex_count\": " << strip.metrics.vertexCount
                << ", \"layer_cell_count\": " << strip.metrics.cellCount
                << ", \"n_layers\": " << strip.parameters.nLayers
                << ", \"first_layer_thickness\": "
                << strip.parameters.firstLayerThickness
                << ", \"total_thickness\": " << strip.parameters.totalThickness
                << ", \"growth_ratio\": " << strip.parameters.growthRatio
                << ", \"safe_thickness_limit\": ";
            if (std::isfinite(strip.metrics.safeThicknessLimit)) {
                out << strip.metrics.safeThicknessLimit;
            } else out << "null";
            out << '}';
            if (i + 1U != result.strips.size()) out << ',';
            out << '\n';
        }
        out << "  ]\n}\n";
    }
    if (!out.good()) {
        setError(error, "failed while writing boundary-layer JSON report");
        return false;
    }
    return true;
}

} // namespace cartmesh2d
