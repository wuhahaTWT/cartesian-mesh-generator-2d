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
    return std::max(maxX - minX, maxY - minY);
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
    // direction is unit length: rayParameter and epsilon are lengths,
    // whereas edgeParameter is dimensionless. Compare endpoint excess in
    // physical length so a change of model units cannot create a collision.
    const double edgeLength = norm(edge);
    if (rayParameter <= epsilon || edgeParameter * edgeLength < -epsilon ||
        (edgeParameter - 1.0) * edgeLength > epsilon) return std::nullopt;
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
           policy.collisionClearanceFraction < 0.5 &&
           (!policy.maximumColumnLength ||
            (std::isfinite(*policy.maximumColumnLength) &&
             *policy.maximumColumnLength > 0.0));
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
    std::vector<double> vertexSafeThickness;
    std::vector<double> segmentSafeThickness;
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
    const BoundaryLayerPolicy2D& policy,
    bool allowLocalReduction = false) {
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
    geometry.vertexSafeThickness.resize(
        chain.vertices.size(), std::numeric_limits<double>::infinity());
    geometry.segmentSafeThickness.resize(
        chain.segments.size(), std::numeric_limits<double>::infinity());

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
        bool marchingConcave=false;
        if (solidTurn < -policy.smoothTurnRadians) {
            geometry.kinds[vertexId] = WallVertexKind2D::Concave;
            if (!allowLocalReduction && !policy.permitConcaveTerminationMarching) {
                result.failure = failedResult(
                    BoundaryLayerFailureReason2D::ConcaveCorner,
                    "solid-side concave corner is outside H4-1 capability", chain.id,
                    &parameters, vertexId);
                return result;
            }
            if (!policy.permitConcaveTerminationMarching) {
                geometry.vertexSafeThickness[vertexId] = 0.0;
                geometry.directions[vertexId] = previousNormal;
                geometry.miterCosines[vertexId] = 1.0;
                continue;
            }
            marchingConcave=true;
        }
        if (solidTurn > policy.maxConvexTurnRadians) {
            geometry.kinds[vertexId] = WallVertexKind2D::Sharp;
            if (!allowLocalReduction) {
                result.failure = failedResult(
                    BoundaryLayerFailureReason2D::SharpCorner,
                    "wall turn exceeds the H4-1 sharp-corner limit", chain.id,
                    &parameters, vertexId);
                return result;
            }
            geometry.vertexSafeThickness[vertexId] = 0.0;
            geometry.directions[vertexId] = previousNormal;
            geometry.miterCosines[vertexId] = 1.0;
            continue;
        }

        const Vector2D sum{previousNormal.x + nextNormal.x,
                           previousNormal.y + nextNormal.y};
        const auto marching = unit(sum, epsilon);
        if (!marching) {
            geometry.kinds[vertexId] = WallVertexKind2D::Degenerate;
            if (!allowLocalReduction) {
                result.failure = failedResult(
                    BoundaryLayerFailureReason2D::ConflictingHalfPlanes,
                    "adjacent fluid-side half-planes have no stable marching direction",
                    chain.id, &parameters, vertexId);
                return result;
            }
            geometry.vertexSafeThickness[vertexId] = 0.0;
            geometry.directions[vertexId] = previousNormal;
            geometry.miterCosines[vertexId] = 1.0;
            continue;
        }
        const double previousCosine = dot(*marching, previousNormal);
        const double nextCosine = dot(*marching, nextNormal);
        const double miterCosine = std::min(previousCosine, nextCosine);
        const double cosineEpsilon = std::max(1.0e-8, epsilon / scale);
        if (!std::isfinite(miterCosine) || miterCosine <= cosineEpsilon) {
            geometry.kinds[vertexId] = WallVertexKind2D::Degenerate;
            if (!allowLocalReduction) {
                result.failure = failedResult(
                    BoundaryLayerFailureReason2D::ConflictingHalfPlanes,
                    "marching direction violates a fluid-side half-plane", chain.id,
                    &parameters, vertexId);
                return result;
            }
            geometry.vertexSafeThickness[vertexId] = 0.0;
            geometry.directions[vertexId] = previousNormal;
            geometry.miterCosines[vertexId] = 1.0;
            continue;
        }
        geometry.directions[vertexId] = *marching;
        geometry.miterCosines[vertexId] = miterCosine;
        geometry.kinds[vertexId] = marchingConcave
            ?WallVertexKind2D::Concave:
             (std::abs(solidTurn) <= policy.smoothTurnRadians
                ?WallVertexKind2D::Smooth:WallVertexKind2D::MildConvex);

        const double tangentialAmplification = std::abs(std::tan(0.5 * solidTurn));
        if (tangentialAmplification > cosineEpsilon) {
            const double localLimit = policy.cornerLengthFraction *
                std::min(lengths[previous], lengths[next]) / tangentialAmplification;
            geometry.safeThickness = std::min(geometry.safeThickness, localLimit);
            geometry.vertexSafeThickness[vertexId] = std::min(
                geometry.vertexSafeThickness[vertexId], localLimit);
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
                geometry.vertexSafeThickness[vertexId] = std::min(
                    geometry.vertexSafeThickness[vertexId], normalLimit);
            }
        }
    }

    // For distinct chains retain a segment-local symmetric gap bound. This is
    // the H4-3 clearance signal: only columns facing a narrow neighbour lose
    // outer layers; unrelated parts of the same wall keep the request.
    for (const auto& otherChain : allChains) {
        if (otherChain.id == chain.id) continue;
        for (std::size_t localSegment=0;localSegment<chain.segments.size();
             ++localSegment) {
            double gap = std::numeric_limits<double>::infinity();
            const auto& lhs=chain.segments[localSegment];
            for (const auto& rhs : otherChain.segments) {
                gap = std::min(gap, segmentDistance(lhs, rhs, policy.tolerance));
            }
            const double limit=policy.collisionClearanceFraction*gap;
            geometry.segmentSafeThickness[localSegment]=std::min(
                geometry.segmentSafeThickness[localSegment],limit);
            geometry.safeThickness=std::min(geometry.safeThickness,limit);
        }
    }

    if (!allowLocalReduction &&
        parameters.totalThickness > geometry.safeThickness + epsilon) {
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
    const double areaEpsilon = policy.tolerance.areaScale(scale);
    const std::size_t ringSize = strip.wallChain.vertices.size();
    const std::size_t ringCount = strip.ringVertexIds.size();
    const bool variableColumns=strip.actualLayerCounts.size()==
        strip.wallChain.segments.size() &&
        std::any_of(strip.actualLayerCounts.begin(),strip.actualLayerCounts.end(),
                    [&](std::size_t count) {
                        return count!=strip.parameters.nLayers;
                    });
    std::vector<bool> usedVertices(strip.vertices.size(),false);
    for (const auto& cell:strip.cells) {
        for (const auto id:cell.vertices) usedVertices[id]=true;
    }
    for (const auto id:strip.outerEnvelopeVertexIds) usedVertices[id]=true;

    for (std::size_t lhs = 0; lhs < strip.vertices.size(); ++lhs) {
        if (!usedVertices[lhs]) continue;
        for (std::size_t rhs = lhs + 1U; rhs < strip.vertices.size(); ++rhs) {
            if (!usedVertices[rhs]) continue;
            if (nearlyEqual(strip.vertices[lhs].point, strip.vertices[rhs].point,
                            policy.tolerance)) {
                return failedResult(
                    BoundaryLayerFailureReason2D::DuplicateGeometricVertex,
                    "distinct layer vertex IDs occupy the same geometric point",
                    strip.wallChain.id, &strip.parameters, strip.vertices[rhs].chainVertex);
            }
        }
    }

    // Strict H4-1 retains complete rings. H4-3 validates the exact stepped
    // outer envelope instead, because unused points beyond a terminated column
    // are deliberately not part of the mesh.
    if (!variableColumns) for (std::size_t ring = 0; ring < ringCount; ++ring) {
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
    if (strip.wallChain.closed) {
        BoundaryLoop envelope(strip.outerEnvelope());
        const auto diagnostics=envelope.diagnose(policy.tolerance);
        if (!diagnostics.valid() ||
            diagnostics.orientation!=(strip.wallChain.orientation==WallChainOrientation2D::Clockwise
                ?LoopOrientation::Clockwise:LoopOrientation::CounterClockwise)) {
            return failedResult(
                BoundaryLayerFailureReason2D::EnvelopeSelfIntersection,
                "retained local-layer envelope must be simple and keep the wall orientation",
                strip.wallChain.id,&strip.parameters);
        }
    }

    // Outer envelope and every non-wall ring must not intersect the original
    // wall. Hair edges may touch only their own wall vertex.
    if (!variableColumns) for (std::size_t ring = 1U; ring < ringCount; ++ring) {
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
        const std::size_t previousCount=vertexId==0U
            ?(strip.wallChain.closed?strip.actualLayerCounts.back():0U)
            :strip.actualLayerCounts[vertexId-1U];
        const std::size_t nextCount=vertexId<strip.actualLayerCounts.size()
            ?strip.actualLayerCounts[vertexId]
            :(strip.wallChain.closed?strip.actualLayerCounts.front():0U);
        const std::size_t usedLayerCount=std::max(previousCount,nextCount);
        for (std::size_t ring = 0; ring < usedLayerCount; ++ring) {
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
        std::array<Point2D, 4> quadPoints{};
        std::size_t pointIndex = 0U;
        for (const auto vertexId : cell.vertices) {
            quadPoints[pointIndex++] = strip.vertices[vertexId].point;
            polygon.vertices.push_back(strip.vertices[vertexId].point);
        }
        const double signedArea = polygon.signedArea();
        if (!std::isfinite(signedArea) || signedArea <= areaEpsilon) {
            return failedResult(BoundaryLayerFailureReason2D::NegativeAreaCell,
                                "layer quad has non-positive or tolerance-scale area",
                                strip.wallChain.id, &strip.parameters,
                                std::nullopt, std::nullopt, cell.id);
        }
        if (!isConvexBoundaryLayerQuad2D(quadPoints, policy.tolerance)) {
            return failedResult(BoundaryLayerFailureReason2D::NonConvexCell,
                                "layer quad is not strictly convex",
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
            std::optional<std::size_t> segment;
            const std::size_t segCount = strip.wallChain.segmentCount();
            if (a.chainVertex < segCount && ((a.chainVertex + 1U) % ringSize) == b.chainVertex) {
                segment = a.chainVertex;
            } else if (b.chainVertex < segCount && ((b.chainVertex + 1U) % ringSize) == a.chainVertex) {
                segment = b.chainVertex;
            }
            if (segment && a.ring<=strip.actualLayerCounts[*segment]) {
                expected=(a.ring==0U ||
                          a.ring==strip.actualLayerCounts[*segment])?1U:2U;
            }
        } else if (a.chainVertex == b.chainVertex &&
                   (a.ring + 1U == b.ring || b.ring + 1U == a.ring)) {
            const std::size_t layer=std::min(a.ring,b.ring);
            const std::size_t previousCount=a.chainVertex==0U
                ?(strip.wallChain.closed?strip.actualLayerCounts.back():0U)
                :strip.actualLayerCounts[a.chainVertex-1U];
            const std::size_t nextCount=
                a.chainVertex<strip.actualLayerCounts.size()
                ?strip.actualLayerCounts[a.chainVertex]
                :(strip.wallChain.closed?strip.actualLayerCounts.front():0U);
            expected=static_cast<std::size_t>(layer<previousCount)+
                     static_cast<std::size_t>(layer<nextCount);
        }
        if (expected == 0U || count != expected) {
            return failedResult(BoundaryLayerFailureReason2D::TopologyConflict,
                                "layer edge incidence does not match strip topology",
                                strip.wallChain.id, &strip.parameters);
        }
    }

    // Deterministic X-axis sweep-line broad-phase followed by exact edge and
    // containment checks prevents visually plausible but overlapping non-neighbour
    // quads while scaling cleanly to high-Re multi-layer strips.
    std::vector<AABB2D> bounds;
    bounds.reserve(strip.cells.size());
    for (const auto& cell : strip.cells) bounds.push_back(cellBounds(cell, strip.vertices));
    std::vector<std::size_t> sweepOrder(strip.cells.size());
    for (std::size_t i = 0; i < sweepOrder.size(); ++i) sweepOrder[i] = i;
    std::sort(sweepOrder.begin(), sweepOrder.end(), [&](std::size_t a, std::size_t b) {
        return std::tie(bounds[a].min.x, bounds[a].min.y, a) <
               std::tie(bounds[b].min.x, bounds[b].min.y, b);
    });
    for (std::size_t oi = 0; oi < sweepOrder.size(); ++oi) {
        const std::size_t aId = sweepOrder[oi];
        const auto& aBounds = bounds[aId];
        for (std::size_t oj = oi + 1U; oj < sweepOrder.size(); ++oj) {
            const std::size_t bId = sweepOrder[oj];
            const auto& bBounds = bounds[bId];
            if (bBounds.min.x > aBounds.max.x + lengthEpsilon) break;
            if (aBounds.max.y + lengthEpsilon < bBounds.min.y ||
                bBounds.max.y + lengthEpsilon < aBounds.min.y) continue;
            const std::size_t lhsId = std::min(aId, bId);
            const std::size_t rhsId = std::max(aId, bId);
            const auto& lhs = strip.cells[lhsId];
            const auto& rhs = strip.cells[rhsId];
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
    for (std::size_t vertexId = 0; vertexId < ringSize; ++vertexId) {
        const std::size_t previousCount=vertexId==0U
            ?(strip.wallChain.closed?strip.actualLayerCounts.back():0U)
            :strip.actualLayerCounts[vertexId-1U];
        const std::size_t nextCount=vertexId<strip.actualLayerCounts.size()
            ?strip.actualLayerCounts[vertexId]
            :(strip.wallChain.closed?strip.actualLayerCounts.front():0U);
        const std::size_t usedLayerCount=std::max(previousCount,nextCount);
        for (std::size_t ring=0;ring<usedLayerCount;++ring) {
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
    const ChainGeometry2D& geometry, const std::vector<double>& cumulative,
    std::vector<std::size_t> actualLayerCounts = {}) {
    BoundaryLayerStrip2D strip;
    strip.wallChain = chain;
    strip.parameters = parameters;
    strip.wallVertexKinds = geometry.kinds;
    strip.marchingDirections = geometry.directions;
    const std::size_t ringSize = chain.vertices.size();
    const std::size_t segmentCount = chain.segments.size();
    if (actualLayerCounts.empty()) {
        actualLayerCounts.assign(segmentCount,cumulative.size());
    }
    strip.actualLayerCounts=std::move(actualLayerCounts);
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

    strip.cells.reserve(cumulative.size() * segmentCount);
    for (std::size_t layer = 0; layer < cumulative.size(); ++layer) {
        for (std::size_t segment = 0; segment < segmentCount; ++segment) {
            if (layer>=strip.actualLayerCounts[segment]) continue;
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

    // Trace the exact boundary of the retained column union. A count change
    // contributes one or more exposed hair-edge fragments; these are the real
    // conformal termination interface consumed by H4-2/H4-3.
    if (chain.closed && !strip.actualLayerCounts.empty()) {
        const std::size_t firstCount=strip.actualLayerCounts.front();
        strip.outerEnvelopeVertexIds.push_back(
            strip.ringVertexIds[firstCount][0U]);
        for (std::size_t segment=0;segment<segmentCount;++segment) {
            const std::size_t next=(segment+1U)%ringSize;
            const std::size_t count=strip.actualLayerCounts[segment];
            const std::size_t nextCount=
                strip.actualLayerCounts[(segment+1U)%segmentCount];
            const auto edgeEnd=strip.ringVertexIds[count][next];
            if (edgeEnd!=strip.outerEnvelopeVertexIds.front()) {
                strip.outerEnvelopeVertexIds.push_back(edgeEnd);
            }
            if (count<nextCount) {
                for (std::size_t ring=count+1U;ring<=nextCount;++ring) {
                    const auto id=strip.ringVertexIds[ring][next];
                    if (id!=strip.outerEnvelopeVertexIds.front()) {
                        strip.outerEnvelopeVertexIds.push_back(id);
                    }
                }
            } else {
                for (std::size_t ring=count;ring>nextCount;--ring) {
                    const auto id=strip.ringVertexIds[ring-1U][next];
                    if (id!=strip.outerEnvelopeVertexIds.front()) {
                        strip.outerEnvelopeVertexIds.push_back(id);
                    }
                }
            }
        }
    }
    strip.metrics.requestedTotalThickness = parameters.totalThickness;
    const auto maximumCount=*std::max_element(
        strip.actualLayerCounts.begin(),strip.actualLayerCounts.end());
    strip.metrics.usedTotalThickness=maximumCount==0U
        ?0.0:cumulative[maximumCount-1U];
    strip.metrics.safeThicknessLimit = geometry.safeThickness;
    strip.metrics.requestedColumnCellCount=segmentCount*parameters.nLayers;
    strip.metrics.retainedColumnCellCount=strip.cells.size();
    strip.metrics.zeroLayerColumnCount=static_cast<std::size_t>(std::count(
        strip.actualLayerCounts.begin(),strip.actualLayerCounts.end(),0U));
    for (std::size_t segment=0;segment<segmentCount;++segment) {
        if (strip.actualLayerCounts[segment]!=
            strip.actualLayerCounts[(segment+1U)%segmentCount]) {
            ++strip.metrics.terminationEdgeCount;
        }
    }
    return strip;
}

[[nodiscard]] bool stripsIntersect(const BoundaryLayerStrip2D& lhs,
                                   const BoundaryLayerStrip2D& rhs,
                                   const TolerancePolicy& tol) {
    if (lhs.cells.empty() || rhs.cells.empty()) return false;
    auto stripBounds = [](const BoundaryLayerStrip2D& strip, std::vector<AABB2D>& cellBoxes) {
        cellBoxes.reserve(strip.cells.size());
        AABB2D overall = cellBounds(strip.cells.front(), strip.vertices);
        for (const auto& cell : strip.cells) {
            const auto box = cellBounds(cell, strip.vertices);
            cellBoxes.push_back(box);
            overall.min.x = std::min(overall.min.x, box.min.x);
            overall.min.y = std::min(overall.min.y, box.min.y);
            overall.max.x = std::max(overall.max.x, box.max.x);
            overall.max.y = std::max(overall.max.y, box.max.y);
        }
        return overall;
    };
    std::vector<AABB2D> lhsBoxes;
    std::vector<AABB2D> rhsBoxes;
    const auto lhsOverall = stripBounds(lhs, lhsBoxes);
    const auto rhsOverall = stripBounds(rhs, rhsBoxes);
    const double overallScale = std::max({lhsOverall.max.x - lhsOverall.min.x,
                                          lhsOverall.max.y - lhsOverall.min.y,
                                          rhsOverall.max.x - rhsOverall.min.x,
                                          rhsOverall.max.y - rhsOverall.min.y, 1.0});
    const double broadEps = tol.scale(overallScale);
    if (!boxesOverlap(lhsOverall, rhsOverall, broadEps)) return false;

    std::vector<std::size_t> rhsOrder;
    rhsOrder.reserve(rhs.cells.size());
    for (std::size_t i = 0; i < rhs.cells.size(); ++i) {
        if (boxesOverlap(rhsBoxes[i], lhsOverall, broadEps)) rhsOrder.push_back(i);
    }
    if (rhsOrder.empty()) return false;
    std::sort(rhsOrder.begin(), rhsOrder.end(), [&](std::size_t a, std::size_t b) {
        return rhsBoxes[a].min.x < rhsBoxes[b].min.x;
    });
    double maxRhsWidth = 0.0;
    for (const auto idx : rhsOrder) {
        maxRhsWidth = std::max(maxRhsWidth, rhsBoxes[idx].max.x - rhsBoxes[idx].min.x);
    }

    for (std::size_t li = 0; li < lhs.cells.size(); ++li) {
        const auto& lhsBounds = lhsBoxes[li];
        if (!boxesOverlap(lhsBounds, rhsOverall, broadEps)) continue;
        const auto& lhsCell = lhs.cells[li];
        const double lowerX = lhsBounds.min.x - maxRhsWidth - broadEps;
        auto it = std::lower_bound(rhsOrder.begin(), rhsOrder.end(), lowerX,
            [&](std::size_t idx, double x) { return rhsBoxes[idx].min.x < x; });
        for (; it != rhsOrder.end(); ++it) {
            const auto& rhsBounds = rhsBoxes[*it];
            if (rhsBounds.min.x > lhsBounds.max.x + broadEps) break;
            const auto& rhsCell = rhs.cells[*it];
            const double scale = std::max({lhsBounds.max.x - lhsBounds.min.x,
                                           lhsBounds.max.y - lhsBounds.min.y,
                                           rhsBounds.max.x - rhsBounds.min.x,
                                           rhsBounds.max.y - rhsBounds.min.y, 1.0});
            if (!boxesOverlap(lhsBounds, rhsBounds, tol.scale(scale))) continue;
            for (std::size_t le = 0; le < 4U; ++le) {
                const Segment2D left{lhs.vertices[lhsCell.vertices[le]].point,
                    lhs.vertices[lhsCell.vertices[(le + 1U) % 4U]].point};
                for (std::size_t re = 0; re < 4U; ++re) {
                    const Segment2D right{rhs.vertices[rhsCell.vertices[re]].point,
                        rhs.vertices[rhsCell.vertices[(re + 1U) % 4U]].point};
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

[[nodiscard]] std::size_t layersWithinLimit(
    const std::vector<double>& cumulative,double limit,double epsilon) noexcept {
    if (!std::isfinite(limit)) return cumulative.size();
    return static_cast<std::size_t>(std::upper_bound(
        cumulative.begin(),cumulative.end(),limit+epsilon)-cumulative.begin());
}

void enforceOneLayerTaper(std::vector<std::size_t>& counts,bool closed) {
    if (counts.size()<2U) return;
    bool changed=true;
    while (changed) {
        changed=false;
        for (std::size_t segment=0;segment<counts.size();++segment) {
            if (!closed && segment==0U) continue;
            const std::size_t previous=segment==0U?counts.size()-1U:segment-1U;
            if (counts[segment]>counts[previous]+1U) {
                counts[segment]=counts[previous]+1U;
                changed=true;
            }
            if (counts[previous]>counts[segment]+1U) {
                counts[previous]=counts[segment]+1U;
                changed=true;
            }
        }
    }
}

[[nodiscard]] std::vector<std::size_t> planLocalLayerCounts(
    const WallChain2D& chain,const ResolvedLayerParameters2D& parameters,
    const ChainGeometry2D& geometry,const BoundaryLayerPolicy2D& policy) {
    std::vector<std::size_t> counts(chain.segmentCount(),parameters.nLayers);
    const double epsilon=policy.tolerance.scale(geometryScale(chain));
    for (std::size_t segment=0;segment<chain.segmentCount();++segment) {
        const std::size_t next=(segment+1U)%chain.vertices.size();
        const bool dangerous=
            (!policy.permitConcaveTerminationMarching &&
             geometry.kinds[segment]==WallVertexKind2D::Concave) ||
            geometry.kinds[segment]==WallVertexKind2D::Sharp ||
            geometry.kinds[segment]==WallVertexKind2D::Degenerate ||
            (!policy.permitConcaveTerminationMarching &&
             geometry.kinds[next]==WallVertexKind2D::Concave) ||
            geometry.kinds[next]==WallVertexKind2D::Sharp ||
            geometry.kinds[next]==WallVertexKind2D::Degenerate;
        if (dangerous) {
            counts[segment]=0U;
            continue;
        }
        const double limit=std::min({
            geometry.segmentSafeThickness[segment],
            geometry.vertexSafeThickness[segment],
            geometry.vertexSafeThickness[next]});
        counts[segment]=layersWithinLimit(
            parameters.cumulativeNormalDistances,limit,epsilon);
    }
    enforceOneLayerTaper(counts,chain.closed);
    return counts;
}

} // namespace

std::optional<double> detail::boundaryLayerRaySegmentDistance2D(
    const Point2D& origin,const Vector2D& direction,const Segment2D& segment,
    double lengthEpsilon) noexcept {
    return raySegmentDistance(origin,direction,segment,lengthEpsilon);
}

bool isConvexBoundaryLayerQuad2D(
    const std::array<Point2D, 4>& points,
    const TolerancePolicy& tol) noexcept {
    for (std::size_t i = 0; i < points.size(); ++i) {
        const Point2D& previous = points[(i + points.size() - 1U) % points.size()];
        const Point2D& current = points[i];
        const Point2D& next = points[(i + 1U) % points.size()];
        const Vector2D incoming = current - previous;
        const Vector2D outgoing = next - current;
        const double product = std::sqrt(squaredNorm(incoming) * squaredNorm(outgoing));
        if (!(product > 0.0) ||
            cross(incoming, outgoing) <= tol.areaScale(std::sqrt(product))) {
            return false;
        }
    }
    return true;
}

WallChainBuildResult2D refineWallChainToSize2D(
    const WallChain2D& chain,double maximumSegmentLength,const TolerancePolicy& tol) {
    WallChainBuildResult2D result;
    const auto fail=[&](const char* message) {
        result.failureReason=WallChainFailureReason2D::InvalidBoundary;
        result.message=message;return result;
    };
    if (!(maximumSegmentLength>0) || !std::isfinite(maximumSegmentLength) ||
        chain.vertices.size()<2 || chain.segmentCount()!=(chain.closed?chain.vertices.size():chain.vertices.size()-1))
        return fail("invalid wall chain or tangential size");
    constexpr std::size_t maximumPoints=1000000;
    WallChain2D refined=chain;refined.vertices.clear();refined.segments.clear();
    for (const auto& segment:chain.segments) {
        const double length=std::hypot(segment.b.x-segment.a.x,segment.b.y-segment.a.y);
        double ratio=length/maximumSegmentLength;
        if (!(length>0) || !std::isfinite(ratio) || ratio>static_cast<double>(maximumPoints))
            return fail("wall tangential sizing exceeds the explicit point budget");
        const double nearest=std::round(ratio);
        if (nearest>=1 && tol.nearlyEqual(ratio,nearest,ratio)) ratio=nearest;
        const auto count=static_cast<std::size_t>(std::max(1.,std::ceil(ratio)));
        if (count>maximumPoints-refined.vertices.size())
            return fail("wall tangential sizing exceeds the explicit point budget");
        refined.vertices.push_back(segment.a);
        for (std::size_t i=1;i<count;++i) {
            const double t=static_cast<double>(i)/static_cast<double>(count);
            refined.vertices.push_back({segment.a.x+t*(segment.b.x-segment.a.x),
                                        segment.a.y+t*(segment.b.y-segment.a.y)});
        }
    }
    if (!chain.closed) refined.vertices.push_back(chain.segments.back().b);
    refined.segments=makeSegments(refined.vertices,refined.closed);
    result.chain=std::move(refined);return result;
}

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
    const auto& ids=outerEnvelopeVertexIds.empty()
        ?ringVertexIds.back():outerEnvelopeVertexIds;
    result.reserve(ids.size());
    for (const auto id : ids) result.push_back(vertices[id].point);
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

BoundaryLayerBuildResult2D buildLocallyReducedBoundaryLayerStrips2D(
    const std::vector<WallChain2D>& wallChains,
    const LayerParameters2D& parameters,
    const BoundaryLayerPolicy2D& policy) {
    if (!policyValid(policy)) {
        return failedResult(BoundaryLayerFailureReason2D::InvalidParameters,
                            "boundary-layer tolerance/angle policy is invalid",0U,
                            nullptr);
    }
    const auto resolvedResult=resolveLayerParameters2D(parameters);
    if (!resolvedResult.success()) {
        auto failure=failedResult(BoundaryLayerFailureReason2D::InvalidParameters,
                                  resolvedResult.message,0U,nullptr);
        failure.failure.nLayers=parameters.nLayers;
        failure.failure.growthRatio=parameters.growthRatio;
        failure.failure.requestedThickness=parameters.thickness;
        return failure;
    }
    const auto& resolved=*resolvedResult.parameters;
    if (wallChains.empty()) {
        return failedResult(BoundaryLayerFailureReason2D::InvalidChain,
                            "at least one wall chain is required",0U,&resolved);
    }
    std::set<std::size_t> chainIds;
    for (const auto& chain:wallChains) {
        if (!chainIds.insert(chain.id).second) {
            return failedResult(BoundaryLayerFailureReason2D::InvalidChain,
                                "wall chain IDs must be unique",chain.id,&resolved);
        }
    }

    // A solver-ready layer column cannot span an arbitrarily long straight
    // wall edge: the common Cartesian partition would attach many small faces
    // to one distant cell centre. Refine only the H4-3 working chain with exact
    // collinear points; the physical wall geometry and patch identity remain
    // unchanged. Four first-layer heights gives bounded layer aspect while
    // leaving the strict H4-1 builder and its historical IDs untouched.
    const double maximumColumnLength =
        policy.maximumColumnLength.value_or(4.0 * resolved.firstLayerThickness);
    std::vector<WallChain2D> workingChains;
    workingChains.reserve(wallChains.size());
    for (const auto& chain:wallChains) {
        std::vector<Point2D> points;
        for (std::size_t segment=0;segment<chain.segments.size();++segment) {
            const auto& a=chain.vertices[segment];
            const auto& b=chain.vertices[(segment+1U)%chain.vertices.size()];
            const double length=norm(b-a);
            const auto minimumPieces=std::max<std::size_t>(
                1U,static_cast<std::size_t>(std::ceil(length/maximumColumnLength)));
            std::size_t pieces=1U;
            while (pieces<minimumPieces) pieces*=2U;
            std::vector<double> parameters;
            if (policy.permitConcaveTerminationMarching && chain.closed &&
                chain.orientation==WallChainOrientation2D::Clockwise) {
                // Reserve the distance consumed by the two corner miters, then
                // subdivide the straight middle independently. A long inlet
                // edge must not inherit the coarse corner spacing throughout.
                double guards[2]{};
                for (std::size_t end=0;end<2U;++end) {
                    const auto vertex=(segment+end)%chain.vertices.size();
                    const auto previous=(vertex+chain.vertices.size()-1U)%chain.vertices.size();
                    const auto next=(vertex+1U)%chain.vertices.size();
                    const auto incoming=chain.vertices[vertex]-chain.vertices[previous];
                    const auto outgoing=chain.vertices[next]-chain.vertices[vertex];
                    const double turn=std::atan2(cross(incoming,outgoing),dot(incoming,outgoing));
                    const double solidTurn=chain.fluidSide==FluidSide2D::Right?turn:-turn;
                    if (solidTurn < -policy.smoothTurnRadians)
                        guards[end]=resolved.totalThickness*
                            std::abs(std::tan(0.5*solidTurn))/policy.cornerLengthFraction;
                }
                if (std::max(guards[0],guards[1])>maximumColumnLength) {
                    const double middle=length-guards[0]-guards[1];
                    parameters.push_back(0.0);
                    if (middle>=maximumColumnLength) {
                        const auto count=static_cast<std::size_t>(std::ceil(middle/maximumColumnLength));
                        for (std::size_t i=0;i<=count;++i) {
                            const double t=(guards[0]+middle*static_cast<double>(i)/count)/length;
                            if (t>0.0 && t<1.0) parameters.push_back(t);
                        }
                    }
                }
            }
            if (parameters.empty()) {
                for (std::size_t piece=0;piece<pieces;++piece)
                    parameters.push_back(static_cast<double>(piece)/pieces);
            }
            for (const double t:parameters)
                points.push_back({a.x+(b.x-a.x)*t,a.y+(b.y-a.y)*t});
        }
        if (chain.closed) {
            const auto refined=makeClosedWallChain2D(
                BoundaryLoop(std::move(points)),chain.id,chain.patchIdentity,
                chain.orientation==WallChainOrientation2D::CounterClockwise
                    ?WallFluidRegion2D::Exterior:WallFluidRegion2D::Interior,
                policy.tolerance);
            if (!refined.success()) {
                return failedResult(BoundaryLayerFailureReason2D::InvalidChain,
                                    "locally refined wall chain is invalid",
                                    chain.id,&resolved);
            }
            workingChains.push_back(*refined.chain);
        } else {
            points.push_back(chain.vertices.back());
            const auto refined=makeOpenWallChain2D(
                std::move(points),chain.id,chain.patchIdentity,chain.fluidSide,
                policy.tolerance);
            if (!refined.success()) {
                return failedResult(BoundaryLayerFailureReason2D::InvalidChain,
                                    "locally refined open wall chain is invalid",
                                    chain.id,&resolved);
            }
            workingChains.push_back(*refined.chain);
        }
    }

    std::vector<ChainGeometry2D> geometries;
    std::vector<std::vector<std::size_t>> counts;
    geometries.reserve(workingChains.size());
    counts.reserve(workingChains.size());
    for (const auto& chain:workingChains) {
        auto analysed=analyseChain(chain,workingChains,resolved,policy,true);
        if (!analysed.geometry) return analysed.failure;
        counts.push_back(planLocalLayerCounts(
            chain,resolved,*analysed.geometry,policy));
        geometries.push_back(std::move(*analysed.geometry));
    }

    std::vector<BoundaryLayerStrip2D> accepted;
    accepted.reserve(workingChains.size());
    bool reduced=false;
    for (std::size_t chainIndex=0;chainIndex<workingChains.size();++chainIndex) {
        BoundaryLayerBuildResult2D validation;
        constexpr std::size_t maximumRepairPasses=32U;
        for (std::size_t pass=0;pass<maximumRepairPasses;++pass) {
            auto strip=constructStrip(
                workingChains[chainIndex],resolved,geometries[chainIndex],
                resolved.cumulativeNormalDistances,counts[chainIndex]);
            if (strip.cells.empty()) {
                return failedResult(
                    BoundaryLayerFailureReason2D::ThicknessExceedsSafeLimit,
                    "local reduction removed every boundary-layer column",
                    workingChains[chainIndex].id,&resolved);
            }
            validation=validateStrip(strip,policy);
            if (validation.success()) break;
            if (!validation.failure.cellId ||
                *validation.failure.cellId>=strip.cells.size()) return validation;
            const auto segment=strip.cells[*validation.failure.cellId].wallSegment;
            if (counts[chainIndex][segment]==0U) return validation;
            --counts[chainIndex][segment];
            enforceOneLayerTaper(counts[chainIndex],workingChains[chainIndex].closed);
        }
        if (!validation.success()) return validation;
        reduced=reduced || std::any_of(
            counts[chainIndex].begin(),counts[chainIndex].end(),
            [&](std::size_t count) { return count<resolved.nLayers; });
        accepted.push_back(std::move(validation.strips.front()));
    }

    for (std::size_t lhs=0;lhs<accepted.size();++lhs) {
        for (std::size_t rhs=lhs+1U;rhs<accepted.size();++rhs) {
            if (stripsIntersect(accepted[lhs],accepted[rhs],policy.tolerance)) {
                return failedResult(
                    BoundaryLayerFailureReason2D::ChainCollision,
                    "locally reduced strips from distinct wall chains still collide",
                    accepted[rhs].wallChain.id,&resolved);
            }
        }
    }

    BoundaryLayerBuildResult2D result;
    result.status=BoundaryLayerStatus2D::Success;
    result.strips=std::move(accepted);
    result.localReductionApplied=reduced;
    return result;
}

BoundaryLayerBuildResult2D buildLocallyReducedBoundaryLayerStrip2D(
    const WallChain2D& wallChain,const LayerParameters2D& parameters,
    const BoundaryLayerPolicy2D& policy) {
    return buildLocallyReducedBoundaryLayerStrips2D(
        {wallChain},parameters,policy);
}


} // namespace cartmesh2d
