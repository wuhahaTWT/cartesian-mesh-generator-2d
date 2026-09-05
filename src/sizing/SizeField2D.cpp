#include "cartmesh2d/sizing/SizeField2D.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>

namespace cartmesh2d {
namespace {

struct WallSegment2D {
    std::size_t id = 0;
    Segment2D segment;
    std::size_t loop = 0;
    std::size_t indexInLoop = 0;
    std::size_t loopSize = 0;
};

// Must mirror BoundarySegmentIndex2D exactly: loop by loop, segment i joining
// vertex i to vertex (i+1) % n, one global running id.  A divergence here would
// silently point a segment band at the wrong stretch of wall.
[[nodiscard]] std::vector<WallSegment2D> enumerateWallSegments(
    const BoundaryRegion2D& boundary) {
    std::vector<WallSegment2D> segments;
    std::size_t id = 0;
    for (std::size_t loop = 0; loop < boundary.loops().size(); ++loop) {
        const auto& vertices = boundary.loops()[loop].vertices();
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            segments.push_back({id++,
                                Segment2D{vertices[i], vertices[(i + 1U) % vertices.size()]},
                                loop, i, vertices.size()});
        }
    }
    return segments;
}

// Circumradius of three consecutive wall vertices; infinite when collinear.
[[nodiscard]] double discreteCurvatureRadius(const Point2D& a, const Point2D& b,
                                            const Point2D& c) noexcept {
    const double ab = std::hypot(b.x - a.x, b.y - a.y);
    const double bc = std::hypot(c.x - b.x, c.y - b.y);
    const double ca = std::hypot(a.x - c.x, a.y - c.y);
    const double twiceArea = std::abs(cross(b - a, c - a));
    if (!(twiceArea > 0.0)) return std::numeric_limits<double>::infinity();
    return (ab * bc * ca) / (2.0 * twiceArea);
}

[[nodiscard]] Point2D closestPointOnSegment(const Point2D& point,
                                           const Segment2D& segment) noexcept {
    const Vector2D direction = segment.b - segment.a;
    const double denominator = squaredNorm(direction);
    if (!(denominator > 0.0)) return segment.a;
    const double parameter = std::clamp(dot(point - segment.a, direction) / denominator,
                                        0.0, 1.0);
    return segment.a + direction * parameter;
}

struct SegmentConnector2D {
    double distance = std::numeric_limits<double>::infinity();
    Vector2D direction{};
    bool directed = false;
};

// For two disjoint segments the minimum separation is always attained at an
// endpoint of one of them, so four point-to-segment probes are exact here.
[[nodiscard]] SegmentConnector2D shortestConnector(const Segment2D& p,
                                                  const Segment2D& q) noexcept {
    SegmentConnector2D best;
    const std::pair<Point2D, Point2D> candidates[4]{
        {p.a, closestPointOnSegment(p.a, q)}, {p.b, closestPointOnSegment(p.b, q)},
        {closestPointOnSegment(q.a, p), q.a}, {closestPointOnSegment(q.b, p), q.b}};
    for (const auto& [onP, onQ] : candidates) {
        const double distance = std::hypot(onQ.x - onP.x, onQ.y - onP.y);
        if (distance >= best.distance) continue;
        best.distance = distance;
        if (distance > 0.0) {
            best.direction = Vector2D{(onQ.x - onP.x) / distance, (onQ.y - onP.y) / distance};
            best.directed = true;
        } else {
            best.directed = false;
        }
    }
    return best;
}

[[nodiscard]] bool tangentialTo(const Vector2D& unitConnector, const Segment2D& segment,
                                double maxTangentialCosine) noexcept {
    const Vector2D direction = segment.b - segment.a;
    const double norm = std::sqrt(squaredNorm(direction));
    if (!(norm > 0.0)) return true;
    return std::abs(dot(unitConnector, Vector2D{direction.x / norm, direction.y / norm})) >
           maxTangentialCosine;
}

// Segments sharing a vertex have no usable connector direction.
[[nodiscard]] bool sharesVertexNeighbourhood(const WallSegment2D& a, const WallSegment2D& b,
                                            std::size_t window) noexcept {
    if (a.loop != b.loop) return false;
    const std::size_t n = a.loopSize;
    const std::size_t forward = (b.indexInLoop + n - a.indexInLoop) % n;
    return std::min(forward, n - forward) <= window;
}

// One band per distinct requested level, ids sorted and unique as refine() demands.
void appendSegmentBands(const std::map<std::size_t, std::vector<std::size_t>>& byLevel,
                        double domainSpan, double radiusCells,
                        std::vector<SegmentRefinementBand2D>& out) {
    for (const auto& [level, ids] : byLevel) {
        if (ids.empty()) continue;
        SegmentRefinementBand2D band;
        band.segmentIds = ids;
        std::sort(band.segmentIds.begin(), band.segmentIds.end());
        band.segmentIds.erase(std::unique(band.segmentIds.begin(), band.segmentIds.end()),
                              band.segmentIds.end());
        band.targetLevel = level;
        band.radius = radiusCells * std::ldexp(domainSpan, -static_cast<int>(level));
        out.push_back(std::move(band));
    }
}

} // namespace

std::size_t sizeFieldLevelForSize2D(double domainSpan, double targetSize,
                                    std::size_t maxLevelCap) noexcept {
    if (!(domainSpan > 0.0) || !(targetSize > 0.0) || !std::isfinite(targetSize)) return 0U;
    if (targetSize >= domainSpan) return 0U;
    const double exact = std::log2(domainSpan / targetSize);
    if (!(exact > 0.0)) return 0U;
    const double ceiled = std::ceil(exact);
    if (ceiled >= static_cast<double>(maxLevelCap)) return maxLevelCap;
    return static_cast<std::size_t>(ceiled);
}

ResolvedSizeField2D resolveSizeField2D(const SizeFieldPolicy2D& policy,
                                      const BoundaryRegion2D& boundary,
                                      std::size_t maxLevelCap,
                                      const TolerancePolicy& tol) {
    ResolvedSizeField2D resolved;
    const auto reject = [&resolved](std::string message) {
        resolved.issues.push_back(std::move(message));
        return resolved;
    };
    if (maxLevelCap == 0U || maxLevelCap > 28U) return reject("maxLevelCap outside [1,28]");
    if (!boundary.diagnose(tol).valid()) return reject("invalid boundary region");
    if (!(policy.farFieldSpans > 0.0) || !std::isfinite(policy.farFieldSpans)) {
        return reject("farFieldSpans must be finite and positive");
    }
    if (policy.wallCellsPerSpan &&
        (!(*policy.wallCellsPerSpan > 0.0) || !std::isfinite(*policy.wallCellsPerSpan))) {
        return reject("wallCellsPerSpan must be finite and positive");
    }
    if (policy.maxSafeWallLevel == 0U || policy.maxSafeWallLevel > maxLevelCap) {
        return reject("maxSafeWallLevel outside [1,maxLevelCap]");
    }

    const AABB2D bodyBounds = boundary.bounds();
    const double bodyWidth = bodyBounds.max.x - bodyBounds.min.x;
    const double bodyHeight = bodyBounds.max.y - bodyBounds.min.y;
    resolved.bodySpan = std::max(bodyWidth, bodyHeight);
    if (!(resolved.bodySpan > 0.0)) return reject("boundary has zero extent");

    // A square domain is not cosmetic: Quadtree2D divides width and height by the
    // same power of two, so any other aspect ratio makes every cell in the mesh
    // non-square and leaves "level" without a single cell size to reason about.
    resolved.domainSpan = resolved.bodySpan * (1.0 + 2.0 * policy.farFieldSpans);
    const Point2D centre{0.5 * (bodyBounds.min.x + bodyBounds.max.x),
                         0.5 * (bodyBounds.min.y + bodyBounds.max.y)};
    const double half = 0.5 * resolved.domainSpan;
    resolved.domain = Domain2D{{{centre.x - half, centre.y - half},
                                {centre.x + half, centre.y + half}}};
    if (!resolved.domain.valid(tol)) return reject("resolved domain is degenerate");

    if (policy.wallCellsPerSpan) {
        resolved.wallLevel = sizeFieldLevelForSize2D(
            resolved.domainSpan, resolved.bodySpan / *policy.wallCellsPerSpan, maxLevelCap);
        if (resolved.wallLevel == 0U) return reject("wallCellsPerSpan is too coarse to refine");
    } else {
        resolved.wallLevel = std::min(policy.maxSafeWallLevel, maxLevelCap);
    }    resolved.curvatureLevel = resolved.wallLevel;
    resolved.proximityLevel = resolved.wallLevel;

    const auto segments = enumerateWallSegments(boundary);
    const auto cellSizeAt = [&](std::size_t level) {
        return std::ldexp(resolved.domainSpan, -static_cast<int>(level));
    };

    std::map<std::size_t, std::vector<std::size_t>> curvatureByLevel;
    if (policy.curvature) {
        const auto& curvature = *policy.curvature;
        if (!(curvature.cellsPerRadius > 0.0) || !std::isfinite(curvature.cellsPerRadius)) {
            return reject("cellsPerRadius must be finite and positive");
        }
        for (const auto& wall : segments) {
            if (wall.loopSize < 3U) continue;
            const auto& vertices = boundary.loops()[wall.loop].vertices();
            const std::size_t n = wall.loopSize;
            const Point2D& previous = vertices[(wall.indexInLoop + n - 1U) % n];
            const Point2D& current = vertices[wall.indexInLoop];
            const Point2D& next = vertices[(wall.indexInLoop + 1U) % n];
            const double radius = discreteCurvatureRadius(previous, current, next);
            if (!std::isfinite(radius)) continue;
            const std::size_t level = sizeFieldLevelForSize2D(
                resolved.domainSpan, radius / curvature.cellsPerRadius, maxLevelCap);
            if (level <= resolved.wallLevel) continue;
            curvatureByLevel[level].push_back(wall.id);
            resolved.curvatureLevel = std::max(resolved.curvatureLevel, level);
        }
    }

    std::map<std::size_t, std::vector<std::size_t>> proximityByLevel;
    if (policy.proximity) {
        const auto& proximity = *policy.proximity;
        if (!(proximity.cellsAcrossGap > 0.0) || !std::isfinite(proximity.cellsAcrossGap)) {
            return reject("cellsAcrossGap must be finite and positive");
        }
        if (!(proximity.maxTangentialCosine >= 0.0) ||
            !(proximity.maxTangentialCosine <= 1.0)) {
            return reject("maxTangentialCosine must lie in [0,1]");
        }
        // A gap only needs sizing when the baseline wall cell cannot already fit
        // `cellsAcrossGap` cells into it.  That makes the search radius self-limiting
        // instead of an arbitrary constant.
        const double searchRadius = proximity.cellsAcrossGap * cellSizeAt(resolved.wallLevel);
        for (const auto& wall : segments) {
            double gap = std::numeric_limits<double>::infinity();
            for (const auto& other : segments) {
                if (other.id == wall.id) continue;
                if (sharesVertexNeighbourhood(wall, other, proximity.sharedVertexIndexWindow)) {
                    continue;
                }
                if (std::max(wall.segment.a.x, wall.segment.b.x) + searchRadius <
                        std::min(other.segment.a.x, other.segment.b.x) ||
                    std::max(other.segment.a.x, other.segment.b.x) + searchRadius <
                        std::min(wall.segment.a.x, wall.segment.b.x) ||
                    std::max(wall.segment.a.y, wall.segment.b.y) + searchRadius <
                        std::min(other.segment.a.y, other.segment.b.y) ||
                    std::max(other.segment.a.y, other.segment.b.y) + searchRadius <
                        std::min(wall.segment.a.y, wall.segment.b.y)) {
                    continue;
                }
                const auto connector = shortestConnector(wall.segment, other.segment);
                if (!connector.directed || connector.distance > searchRadius) continue;
                // Nearly tangential connectors are chords of one smooth stretch of
                // wall, not a gap across fluid.
                if (tangentialTo(connector.direction, wall.segment,
                                 proximity.maxTangentialCosine) ||
                    tangentialTo(connector.direction, other.segment,
                                 proximity.maxTangentialCosine)) {
                    continue;
                }
                gap = std::min(gap, connector.distance);
            }
            if (!std::isfinite(gap) || !(gap > 0.0) || gap > searchRadius) continue;
            const std::size_t level = sizeFieldLevelForSize2D(
                resolved.domainSpan, gap / proximity.cellsAcrossGap, maxLevelCap);
            if (level <= resolved.wallLevel) continue;
            proximityByLevel[level].push_back(wall.id);
            resolved.proximityLevel = std::max(resolved.proximityLevel, level);
        }
    }

    resolved.maxLevel = std::max({resolved.wallLevel, resolved.curvatureLevel,
                                  resolved.proximityLevel});
    resolved.levelCapReached = resolved.maxLevel >= maxLevelCap;
    resolved.wallCellSize = cellSizeAt(resolved.wallLevel);
    if (!policy.allowUnsafeWallLevel && resolved.maxLevel > policy.maxSafeWallLevel) {
        // Fail closed with the arithmetic in the message: the caller needs to know
        // which request produced the depth, not just that a limit exists.
        std::ostringstream message;
        message << "resolved level " << resolved.maxLevel << " exceeds maxSafeWallLevel "
                << policy.maxSafeWallLevel << " (wall " << resolved.wallLevel
                << ", curvature " << resolved.curvatureLevel << ", proximity "
                << resolved.proximityLevel << "); lower wallCellsPerSpan or"
                   " farFieldSpans, or set allowUnsafeWallLevel";
        return reject(message.str());
    }

    auto& refinement = resolved.refinement;
    refinement.boundaryLevel = resolved.wallLevel;

    if (policy.wallDistance) {
        const auto& wallDistance = *policy.wallDistance;
        if (wallDistance.farLevel > resolved.wallLevel) {
            return reject("farLevel exceeds the resolved wall level");
        }
        refinement.minimumLevel = wallDistance.farLevel;
        // Cumulative thickness: level L reaches out to cellsPerLevel cells of every
        // level from the wall down to L, which is the discrete form of a bounded
        // growth rate.  cellsPerLevel == 0 emits nothing and leaves the pre-existing
        // one-ring-per-level gradation that 2:1 balance produces on its own.
        double distance = 0.0;
        for (std::size_t level = resolved.wallLevel;
             level > wallDistance.farLevel && wallDistance.cellsPerLevel > 0U; --level) {
            distance += static_cast<double>(wallDistance.cellsPerLevel) * cellSizeAt(level);
            refinement.distanceBands.push_back({distance, level});
        }
    }

    appendSegmentBands(curvatureByLevel, resolved.domainSpan,
                       policy.wallDistance
                           ? std::max<double>(1.0, static_cast<double>(
                                 policy.wallDistance->cellsPerLevel))
                           : 1.0,
                       refinement.segmentBands);
    appendSegmentBands(proximityByLevel, resolved.domainSpan,
                       policy.proximity ? policy.proximity->cellsAcrossGap : 1.0,
                       refinement.segmentBands);

    if (policy.wake) {
        const auto& wake = *policy.wake;
        if (!std::isfinite(wake.angleOfAttackDeg) || !(wake.downstreamSpans > 0.0) ||
            !(wake.halfWidthSpans > 0.0) || !std::isfinite(wake.downstreamSpans) ||
            !std::isfinite(wake.halfWidthSpans)) {
            return reject("wake geometry must be finite and positive");
        }
        if (wake.slices == 0U) return reject("wake slices must be positive");
        // Box regions reject level 0, and a wake coarser than the far field would be
        // pointless, so saturate at 1 rather than silently wrapping.
        const std::size_t target = wake.levelsBelowWall >= resolved.wallLevel
            ? 1U
            : resolved.wallLevel - wake.levelsBelowWall;
        const double angle = wake.angleOfAttackDeg * std::acos(-1.0) / 180.0;
        const Vector2D direction{std::cos(angle), std::sin(angle)};
        const Vector2D normal{-direction.y, direction.x};
        const double length = wake.downstreamSpans * resolved.bodySpan;
        const double halfWidth = wake.halfWidthSpans * resolved.bodySpan;
        for (std::size_t slice = 0; slice < wake.slices; ++slice) {
            const double from = length * static_cast<double>(slice) /
                                static_cast<double>(wake.slices);
            const double to = length * static_cast<double>(slice + 1U) /
                              static_cast<double>(wake.slices);
            AABB2D box{{std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity()},
                       {-std::numeric_limits<double>::infinity(),
                        -std::numeric_limits<double>::infinity()}};
            for (const double along : {from, to}) {
                for (const double across : {-halfWidth, halfWidth}) {
                    const Point2D corner = centre + direction * along + normal * across;
                    box.min.x = std::min(box.min.x, corner.x);
                    box.min.y = std::min(box.min.y, corner.y);
                    box.max.x = std::max(box.max.x, corner.x);
                    box.max.y = std::max(box.max.y, corner.y);
                }
            }
            // Clip to the domain and drop empty slices; refine() rejects a region
            // that does not overlap the domain.
            box.min.x = std::max(box.min.x, resolved.domain.bounds.min.x);
            box.min.y = std::max(box.min.y, resolved.domain.bounds.min.y);
            box.max.x = std::min(box.max.x, resolved.domain.bounds.max.x);
            box.max.y = std::min(box.max.y, resolved.domain.bounds.max.y);
            if (!(box.max.x > box.min.x) || !(box.max.y > box.min.y)) continue;
            refinement.boxRegions.push_back({box, target});
        }
        if (refinement.boxRegions.empty()) {
            return reject("wake region does not overlap the resolved domain");
        }
    }

    return resolved;
}

std::string resolvedSizeFieldToJson(const ResolvedSizeField2D& resolved, int indentSpaces) {
    const std::string i1(static_cast<std::size_t>(std::max(indentSpaces, 0)), ' ');
    const std::string i2 = i1 + i1;
    std::ostringstream out;
    out << std::setprecision(17);
    out << "{\n";
    out << i1 << "\"format\": \"cartmesh2d-size-field-v1\",\n";
    out << i1 << "\"body_span\": " << resolved.bodySpan << ",\n";
    out << i1 << "\"domain_span\": " << resolved.domainSpan << ",\n";
    out << i1 << "\"domain\": [" << resolved.domain.bounds.min.x << ", "
        << resolved.domain.bounds.min.y << ", " << resolved.domain.bounds.max.x << ", "
        << resolved.domain.bounds.max.y << "],\n";
    out << i1 << "\"wall_cell_size\": " << resolved.wallCellSize << ",\n";
    out << i1 << "\"wall_level\": " << resolved.wallLevel << ",\n";
    out << i1 << "\"curvature_level\": " << resolved.curvatureLevel << ",\n";
    out << i1 << "\"proximity_level\": " << resolved.proximityLevel << ",\n";
    out << i1 << "\"max_level\": " << resolved.maxLevel << ",\n";
    out << i1 << "\"level_cap_reached\": "
        << (resolved.levelCapReached ? "true" : "false") << ",\n";
    out << i1 << "\"minimum_level\": " << resolved.refinement.minimumLevel << ",\n";
    out << i1 << "\"boundary_level\": " << resolved.refinement.boundaryLevel << ",\n";
    out << i1 << "\"distance_bands\": [\n";
    for (std::size_t i = 0; i < resolved.refinement.distanceBands.size(); ++i) {
        const auto& band = resolved.refinement.distanceBands[i];
        out << i2 << "{\"distance\": " << band.distance << ", \"target_level\": "
            << band.targetLevel << "}" << (i + 1U < resolved.refinement.distanceBands.size() ? "," : "")
            << "\n";
    }
    out << i1 << "],\n";
    out << i1 << "\"segment_bands\": [\n";
    for (std::size_t i = 0; i < resolved.refinement.segmentBands.size(); ++i) {
        const auto& band = resolved.refinement.segmentBands[i];
        out << i2 << "{\"segments\": " << band.segmentIds.size() << ", \"radius\": "
            << band.radius << ", \"target_level\": " << band.targetLevel << "}"
            << (i + 1U < resolved.refinement.segmentBands.size() ? "," : "") << "\n";
    }
    out << i1 << "],\n";
    out << i1 << "\"box_regions\": " << resolved.refinement.boxRegions.size() << ",\n";
    out << i1 << "\"issues\": [\n";
    for (std::size_t i = 0; i < resolved.issues.size(); ++i) {
        out << i2 << "\"" << resolved.issues[i] << "\""
            << (i + 1U < resolved.issues.size() ? "," : "") << "\n";
    }
    out << i1 << "]\n";
    out << "}\n";
    return out.str();
}

} // namespace cartmesh2d

