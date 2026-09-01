#pragma once

#include "cartmesh2d/geometry/ConstructionIdentity2D.hpp"
#include "cartmesh2d/geometry/Geometry2D.hpp"

#include <cstddef>
#include <optional>
#include <limits>
#include <string>
#include <stdexcept>
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
    ConstructionDecision2D constructionDecision = ConstructionDecision2D::Accepted;
    std::string constructionDecisionReason = "compatible_proximity_anchor";
};

struct IntersectionRegistryPolicy2D {
    // Dimensionless: the admissible displacement is this fraction of the
    // smaller local background scale. No absolute coordinate epsilon is used.
    double snapFractionOfLocalH = 64.0*std::numeric_limits<double>::epsilon();

    // R2/W1: the budget for canonicalizing an intersection onto *this support's*
    // arithmetic grid corner, separate from the proximity budget above.
    //
    // It defaults to the same roundoff-sized value, so every existing caller
    // keeps byte-identical construction. A caller that needs a *geometric* weld
    // raises it explicitly and owns the consequence.
    //
    // Why a geometric budget is sometimes needed. Refinement keeps producing
    // wall/grid intersections that land a few 1e-5 of a cell from a grid corner.
    // At roundoff budget they stay two distinct vertices and emit a corner spur
    // with ~1e-7 long faces on a ~3e-3 cell; measured at circle level 10 that
    // spur produced boundary face skewness 8.6 and the solver-quality gate
    // rejected the whole mesh. The registry is the only place where welding can
    // be consistent, because its canonical vertices are shared across leaves: a
    // per-leaf or per-face weld makes neighbouring leaves disagree about a shared
    // cell-side vertex, which was measured to produce 24 unclassified boundary
    // edges (docs/R2_REFINEMENT_ROBUSTNESS_CN.md).
    //
    // Cost, and why it is bounded rather than free. Welding moves a point by at
    // most this fraction of h, so each incident cell's area changes by at most
    // about 1.5 * fraction * h^2, and the total over a boundary of length L is
    // bounded by about 1.5 * fraction * h * L. That is a real, derived
    // perturbation of the fluid-area invariant and callers must budget for it
    // explicitly; see cartmesh2d_cli's physics gate.
    //
    // Anchoring: Q1 already declares any face with face_length / local_h < 0.01
    // a hard failure, so a geometric budget must stay well under 0.01 to be
    // unable to destroy a face the quality contract would have accepted.
    double gridCornerWeldFractionOfLocalH =
        64.0*std::numeric_limits<double>::epsilon();
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

enum class ConstructionConflictKind2D {
    LateNonIncidentFeature,
    NonIncidentFeatureSnap,
    NonIncidentGridCorner
};

struct ConstructionRecoveryRequest2D {
    StableVertexKey2D conflictKey;
    ConstructionConflictKind2D kind =
        ConstructionConflictKind2D::NonIncidentFeatureSnap;
    ConstructionDecision2D recommendedDecision = ConstructionDecision2D::Refine;
    std::vector<ConstructionDecision2D> fallbackOrder{
        ConstructionDecision2D::Refine,
        ConstructionDecision2D::Rephase,
        ConstructionDecision2D::Resample};
    std::size_t supportId = 0;
    GridLineIdentity2D gridLine;
    IntersectionSource2D source = IntersectionSource2D::Unknown;
    Segment2D sourceSegment;
    double sourceParameter = 0.0;
    Point2D originalPoint;
    Point2D conflictingPoint;
    double localH = 0.0;
    bool sourceParameterOnly = true;
    std::string reason;
};

class ConstructionConflict2D final : public std::runtime_error {
public:
    explicit ConstructionConflict2D(ConstructionRecoveryRequest2D request);
    [[nodiscard]] const ConstructionRecoveryRequest2D& request() const noexcept {
        return request_;
    }
private:
    ConstructionRecoveryRequest2D request_;
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
    [[nodiscard]] const std::vector<ConstructionRecoveryRequest2D>& recoveryRequests() const {
        return recoveryRequests_;
    }
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
    std::vector<ConstructionRecoveryRequest2D> recoveryRequests_;
    std::size_t cacheHits_ = 0;
    [[noreturn]] void raiseConstructionConflict(
        ConstructionConflictKind2D kind, std::size_t support,
        GridLineIdentity2D line, const Point2D& originalPoint,
        const Point2D& conflictingPoint, double localH,
        double sourceParameter, std::string reason);
};

[[nodiscard]] const char* intersectionSourceName(IntersectionSource2D source) noexcept;
[[nodiscard]] const char* intersectionFeatureName(IntersectionFeature2D feature) noexcept;
[[nodiscard]] const char* constructionConflictName(ConstructionConflictKind2D kind) noexcept;
[[nodiscard]] std::string intersectionRecordsToJson(
    const std::vector<CanonicalizedIntersection2D>& records);
[[nodiscard]] std::string intersectionConstructionToJson(
    const IntersectionRegistry2D& registry,const std::vector<std::size_t>& solverHandles,
    std::size_t partitionCount,std::size_t partitionCacheHits);

} // namespace cartmesh2d
