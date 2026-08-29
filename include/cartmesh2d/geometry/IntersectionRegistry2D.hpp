#pragma once

#include "cartmesh2d/geometry/ConstructionIdentity2D.hpp"
#include "cartmesh2d/geometry/Geometry2D.hpp"

#include <cstddef>
#include <optional>
#include <limits>
#include <vector>
#include <map>
#include <tuple>
#include <cstdint>
#include <compare>

namespace cartmesh2d {

enum class IntersectionSource2D {
    WallCartesian,
    TransitionEnvelopeCartesian,
    WallTransitionEnvelope,
    Unknown
};

enum class IntersectionFeature2D {
    None,
    Smooth,
    CartesianGridVertex,
    CartesianGridLine,
    TransitionVertex,
    WallSharpCorner,
    WallConcaveCorner
};

struct CanonicalVertex2D {
    std::size_t id = 0U;
    Point2D point;
    double localH = 0.0;
    IntersectionFeature2D feature = IntersectionFeature2D::None;
    std::optional<std::size_t> featureId;
    std::size_t supportId = 0U;
};

struct CanonicalizedIntersection2D {
    std::size_t id = 0U;
    IntersectionSource2D source = IntersectionSource2D::Unknown;
    std::size_t sourceId = 0U;
    Point2D originalPoint;
    CanonicalVertex2D canonicalVertex;
    double displacement = 0.0;
    double localH = 0.0;
    IntersectionFeature2D feature = IntersectionFeature2D::None;
    std::size_t supportId = 0U;
    std::optional<Segment2D> sourceSegment;
    std::optional<std::size_t> solverVertexId;
};

struct IntersectionRegistryPolicy2D {
    // Dimensionless: the admissible displacement is this fraction of the
    // smaller local background scale. No absolute coordinate epsilon is used.
    double snapFractionOfLocalH = 64.0*std::numeric_limits<double>::epsilon();
};

struct GridLineIdentity2D {
    unsigned axis = 0;
    std::uint64_t coordinate = 0;
    auto operator<=>(const GridLineIdentity2D&) const = default;
};

struct SharedIntersectionEvent2D {
    std::size_t supportId = 0;
    GridLineIdentity2D gridLine;
    Point2D originalPoint;
    std::size_t canonicalVertex = 0;
    double displacement = 0;
    double localH = 0;
    IntersectionSource2D source = IntersectionSource2D::Unknown;
    Segment2D sourceSegment;
};

class IntersectionRegistry2D {
public:
    explicit IntersectionRegistry2D(IntersectionRegistryPolicy2D policy = {});

    [[nodiscard]] std::size_t addCanonicalVertex(
        const Point2D& point, double localH, IntersectionFeature2D feature,
        std::optional<std::size_t> featureId = std::nullopt,
        std::size_t supportId = 0U);

    [[nodiscard]] Point2D canonicalize(
        const Point2D& originalPoint, double localH,
        IntersectionSource2D source, std::size_t sourceId,
        IntersectionFeature2D feature = IntersectionFeature2D::None,
        std::optional<std::size_t> featureId = std::nullopt,
        std::size_t supportId = 0U);

    // Construction API: exact identities and support-scoped intersection
    // events, separate from the legacy proximity/sampling helper above.
    void configureGrid(const AABB2D& bounds, std::size_t maxLevel);
    [[nodiscard]] GridLineIdentity2D gridLine(unsigned axis, double coordinate) const;
    [[nodiscard]] double gridCoordinate(GridLineIdentity2D line) const;
    [[nodiscard]] std::size_t internVertex(const Point2D& point, double localH,
        IntersectionFeature2D feature = IntersectionFeature2D::None);
    [[nodiscard]] std::size_t registerSegment(const Segment2D& segment,
        double localH, IntersectionSource2D source);
    [[nodiscard]] std::size_t intersectGridLine(std::size_t support,
        GridLineIdentity2D line, double localH);
    [[nodiscard]] const std::vector<SharedIntersectionEvent2D>& events() const { return events_; }
    [[nodiscard]] std::size_t intersectionCacheHits() const { return cacheHits_; }

    [[nodiscard]] const std::vector<CanonicalVertex2D>& vertices() const noexcept {
        return vertices_;
    }
    [[nodiscard]] const std::vector<CanonicalizedIntersection2D>& records() const noexcept {
        return records_;
    }
    [[nodiscard]] const ConstructionVertexStore2D& shadowVertexStore() const noexcept {
        return shadowVertexStore_;
    }

private:
    IntersectionRegistryPolicy2D policy_;
    std::vector<CanonicalVertex2D> vertices_;
    std::vector<CanonicalizedIntersection2D> records_;
    ConstructionVertexStore2D shadowVertexStore_;
    struct Support {
        Segment2D segment;
        std::size_t a, b;
        IntersectionSource2D source;
    };
    AABB2D gridBounds_{};
    std::size_t gridLevel_ = 0;
    bool gridConfigured_ = false;
    std::map<std::pair<double,double>, std::size_t> exactVertices_;
    std::map<std::tuple<std::size_t,std::size_t,IntersectionSource2D>,std::size_t> supportKeys_;
    std::vector<Support> supports_;
    std::map<std::pair<std::size_t,GridLineIdentity2D>,std::size_t> eventKeys_;
    std::map<std::size_t,std::vector<std::size_t>> vertexEvents_;
    std::vector<SharedIntersectionEvent2D> events_;
    std::size_t cacheHits_ = 0;
};

[[nodiscard]] const char* intersectionSourceName(IntersectionSource2D source) noexcept;
[[nodiscard]] const char* intersectionFeatureName(IntersectionFeature2D feature) noexcept;
[[nodiscard]] std::string intersectionRecordsToJson(
    const std::vector<CanonicalizedIntersection2D>& records);
[[nodiscard]] std::string intersectionConstructionToJson(
    const IntersectionRegistry2D& registry,const std::vector<std::size_t>& solverHandles,
    std::size_t partitionCount,std::size_t partitionCacheHits);

} // namespace cartmesh2d
