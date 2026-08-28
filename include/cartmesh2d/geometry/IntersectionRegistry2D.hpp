#pragma once

#include "cartmesh2d/geometry/Geometry2D.hpp"

#include <cstddef>
#include <optional>
#include <limits>
#include <vector>

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

    [[nodiscard]] const std::vector<CanonicalVertex2D>& vertices() const noexcept {
        return vertices_;
    }
    [[nodiscard]] const std::vector<CanonicalizedIntersection2D>& records() const noexcept {
        return records_;
    }

private:
    IntersectionRegistryPolicy2D policy_;
    std::vector<CanonicalVertex2D> vertices_;
    std::vector<CanonicalizedIntersection2D> records_;
};

[[nodiscard]] const char* intersectionSourceName(IntersectionSource2D source) noexcept;
[[nodiscard]] const char* intersectionFeatureName(IntersectionFeature2D feature) noexcept;
[[nodiscard]] std::string intersectionRecordsToJson(
    const std::vector<CanonicalizedIntersection2D>& records);

} // namespace cartmesh2d
