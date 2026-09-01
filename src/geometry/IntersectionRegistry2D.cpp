#include "cartmesh2d/geometry/IntersectionRegistry2D.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace cartmesh2d {
namespace {

[[nodiscard]] bool finitePositive(double value) noexcept {
    return std::isfinite(value) && value > 0.0;
}

[[nodiscard]] bool protectedFeature(IntersectionFeature2D feature) noexcept {
    return feature == IntersectionFeature2D::WallSharpCorner ||
           feature == IntersectionFeature2D::WallConcaveCorner;
}

[[nodiscard]] int anchorPriority(IntersectionFeature2D feature) noexcept {
    switch (feature) {
    case IntersectionFeature2D::WallSharpCorner:
    case IntersectionFeature2D::WallConcaveCorner: return 0;
    case IntersectionFeature2D::TransitionVertex: return 1;
    case IntersectionFeature2D::CartesianGridVertex: return 2;
    case IntersectionFeature2D::CartesianGridLine: return 2;
    case IntersectionFeature2D::Smooth: return 3;
    case IntersectionFeature2D::None: return 4;
    }
    return 4;
}

[[nodiscard]] bool featureCompatible(
    IntersectionFeature2D candidate, std::optional<std::size_t> candidateId,
    const CanonicalVertex2D& anchor) noexcept {
    if (!protectedFeature(candidate)) return true;
    return candidate == anchor.feature && candidateId && anchor.featureId &&
           *candidateId == *anchor.featureId;
}

[[nodiscard]] ConstructionFeatureClass2D shadowFeature(
    IntersectionFeature2D feature) noexcept {
    switch (feature) {
    case IntersectionFeature2D::Smooth:
        return ConstructionFeatureClass2D::Smooth;
    case IntersectionFeature2D::WallSharpCorner:
        return ConstructionFeatureClass2D::ConvexSharp;
    case IntersectionFeature2D::WallConcaveCorner:
        return ConstructionFeatureClass2D::ConcaveSharp;
    case IntersectionFeature2D::CartesianGridVertex:
    case IntersectionFeature2D::CartesianGridLine:
        return ConstructionFeatureClass2D::Grid;
    case IntersectionFeature2D::TransitionVertex:
        return ConstructionFeatureClass2D::TransitionMutable;
    case IntersectionFeature2D::None:
        return ConstructionFeatureClass2D::Unclassified;
    }
    return ConstructionFeatureClass2D::Unclassified;
}

[[nodiscard]] std::optional<FeatureOwner2D> shadowOwner(
    IntersectionFeature2D feature, std::optional<std::size_t> featureId,
    std::size_t supportId) {
    if (!featureId && feature == IntersectionFeature2D::None) return std::nullopt;
    return FeatureOwner2D{
        feature == IntersectionFeature2D::CartesianGridVertex ||
                feature == IntersectionFeature2D::CartesianGridLine
            ? ConstructionSourceKind2D::CartesianGrid
            : ConstructionSourceKind2D::Unknown,
        static_cast<std::uint64_t>(featureId.value_or(supportId)),
        static_cast<std::uint64_t>(supportId)};
}

[[nodiscard]] const char* constructionDecisionName(
    ConstructionDecision2D decision) noexcept {
    switch (decision) {
    case ConstructionDecision2D::ShadowOnly: return "shadow_only";
    case ConstructionDecision2D::ReusedExactKey: return "reused_exact_key";
    case ConstructionDecision2D::Accepted: return "accepted";
    case ConstructionDecision2D::Rejected: return "rejected";
    case ConstructionDecision2D::Refine: return "refine";
    case ConstructionDecision2D::Resample: return "resample";
    case ConstructionDecision2D::Rephase: return "rephase";
    }
    return "rejected";
}

} // namespace

IntersectionRegistry2D::IntersectionRegistry2D(IntersectionRegistryPolicy2D policy)
    : policy_(policy) {
    if (!finitePositive(policy_.snapFractionOfLocalH) ||
        policy_.snapFractionOfLocalH >= 0.5) {
        throw std::invalid_argument(
            "intersection snap fraction must be finite, positive and dimensionless");
    }
    // A geometric weld budget at or above Q1's short-face hard limit could
    // destroy a face the quality contract would have accepted, so it is refused
    // outright rather than merely discouraged.
    if (!finitePositive(policy_.gridCornerWeldFractionOfLocalH) ||
        policy_.gridCornerWeldFractionOfLocalH >= 0.01) {
        throw std::invalid_argument(
            "grid-corner weld fraction must be finite, positive and below the Q1 "
            "short-face hard limit of 0.01");
    }
}

std::size_t IntersectionRegistry2D::addCanonicalVertex(
    const Point2D& point, double localH, IntersectionFeature2D feature,
    std::optional<std::size_t> featureId,std::size_t supportId) {
    if (!std::isfinite(point.x) || !std::isfinite(point.y) ||
        !finitePositive(localH)) {
        throw std::invalid_argument("canonical vertex requires finite coordinates and local_h");
    }
    const std::size_t id = vertices_.size();
    vertices_.push_back({id, point, localH, feature, featureId, supportId});
    const auto stableId = shadowVertexStore_.addShadowVertex(
        point, localH, shadowFeature(feature),
        shadowOwner(feature, featureId, supportId), supportId);
    if (stableId != static_cast<StableVertexId2D>(id)) {
        throw std::logic_error("R1-A shadow stable vertex id diverged from canonical id");
    }
    return id;
}

Point2D IntersectionRegistry2D::canonicalize(
    const Point2D& originalPoint, double localH,
    IntersectionSource2D source, std::size_t sourceId,
    IntersectionFeature2D feature, std::optional<std::size_t> featureId,
    std::size_t supportId) {
    if (!std::isfinite(originalPoint.x) || !std::isfinite(originalPoint.y) ||
        !finitePositive(localH)) {
        throw std::invalid_argument("intersection requires finite coordinates and local_h");
    }

    const auto consider = [&](std::size_t i, std::optional<std::size_t>& best,
                              double& bestDistance) {
        const auto& anchor = vertices_[i];
        if (anchor.supportId!=supportId ||
            !featureCompatible(feature, featureId, anchor)) return;
        const double admissible = policy_.snapFractionOfLocalH *
                                  std::min(localH, anchor.localH);
        const double distance = std::sqrt(squaredNorm(anchor.point - originalPoint));
        // An input feature is an immutable anchor, not a movable sample.
        if (protectedFeature(feature) && distance!=0.0) return;
        if (distance > admissible) return;
        if (!best || std::tuple(distance, anchorPriority(anchor.feature),
                                anchor.point.x, anchor.point.y, anchor.id) <
                    std::tuple(bestDistance,
                               anchorPriority(vertices_[*best].feature),
                               vertices_[*best].point.x,
                               vertices_[*best].point.y,
                               vertices_[*best].id)) {
            best = i;
            bestDistance = distance;
        }
    };

    std::optional<std::size_t> best;
    double bestDistance = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0U; i < vertices_.size(); ++i) {
        consider(i, best, bestDistance);
    }

    // R1-B production decision: typed feature compatibility and the spatial
    // index select the anchor.  The legacy full scan remains a fail-closed
    // oracle until the later hot-path-removal milestone.
    auto proposalOwner = shadowOwner(feature, featureId, supportId);
    if (protectedFeature(feature) && !featureId) proposalOwner.reset();
    const auto decision = shadowVertexStore_.decideProximity({
        originalPoint, localH, policy_.snapFractionOfLocalH,
        shadowFeature(feature), proposalOwner, supportId,
        protectedFeature(feature)});
    const auto indexedBest = decision.canonicalId
        ? std::optional<std::size_t>{static_cast<std::size_t>(*decision.canonicalId)}
        : std::nullopt;
    if (indexedBest != best ||
        (best && decision.displacement != bestDistance)) {
        throw std::runtime_error(
            "R1-B construction decision disagrees with legacy canonicalization scan");
    }
    if (!indexedBest || decision.displacement == 0.0) return originalPoint;

    const auto& canonical = vertices_[*indexedBest];
    records_.push_back({records_.size(), source, sourceId, originalPoint,
                        canonical, decision.displacement, localH, feature, supportId,
                        std::nullopt,std::nullopt,decision.decision,decision.reason});
    return canonical.point;
}

const char* intersectionSourceName(IntersectionSource2D source) noexcept {
    switch (source) {
    case IntersectionSource2D::WallCartesian: return "wall_cartesian";
    case IntersectionSource2D::TransitionEnvelopeCartesian:
        return "transition_envelope_cartesian";
    case IntersectionSource2D::WallTransitionEnvelope:
        return "wall_transition_envelope";
    case IntersectionSource2D::Unknown: return "unknown";
    }
    return "unknown";
}

const char* intersectionFeatureName(IntersectionFeature2D feature) noexcept {
    switch (feature) {
    case IntersectionFeature2D::None: return "none";
    case IntersectionFeature2D::Smooth: return "smooth";
    case IntersectionFeature2D::CartesianGridVertex: return "cartesian_grid_vertex";
    case IntersectionFeature2D::CartesianGridLine: return "cartesian_grid_line";
    case IntersectionFeature2D::TransitionVertex: return "transition_vertex";
    case IntersectionFeature2D::WallSharpCorner: return "wall_sharp_corner";
    case IntersectionFeature2D::WallConcaveCorner: return "wall_concave_corner";
    }
    return "none";
}

std::string intersectionRecordsToJson(
    const std::vector<CanonicalizedIntersection2D>& records) {
    std::ostringstream out;
    out<<std::setprecision(17);
    const auto point=[&](const Point2D& p) {out<<'['<<p.x<<", "<<p.y<<']';};
    out<<"{\n  \"format_version\": \"cartmesh2d-intersections-v1\",\n"
       <<"  \"canonical_id_scope\": \"local_registry_support\",\n"
       <<"  \"records\": [\n";
    for (std::size_t i=0;i<records.size();++i) {
        const auto& r=records[i];
        out<<"    {\"id\": "<<r.id<<", \"source\": \""<<intersectionSourceName(r.source)
           <<"\", \"source_id\": "<<r.sourceId<<", \"support_id\": "<<r.supportId
           <<", \"original_position\": ";
        point(r.originalPoint);
        out<<", \"canonical_vertex\": {\"local_id\": "<<r.canonicalVertex.id
           <<", \"position\": ";
        point(r.canonicalVertex.point);
        out<<", \"feature\": \""<<intersectionFeatureName(r.canonicalVertex.feature)
           <<"\", \"feature_id\": ";
        if (r.canonicalVertex.featureId) out<<*r.canonicalVertex.featureId; else out<<"null";
        out<<"}, \"solver_vertex_id\": ";
        if (r.solverVertexId) out<<*r.solverVertexId; else out<<"null";
        out<<", \"displacement\": "<<r.displacement<<", \"local_h\": "<<r.localH
           <<", \"displacement_over_local_h\": "<<r.displacement/r.localH
           <<", \"construction_decision\": \""
           <<constructionDecisionName(r.constructionDecision)
           <<"\", \"construction_decision_reason\": \""
           <<r.constructionDecisionReason<<"\""
           <<", \"feature_classification\": \""<<intersectionFeatureName(r.feature)
           <<"\", \"source_segment\": ";
        if (r.sourceSegment) {
            out<<'[';point(r.sourceSegment->a);out<<", ";point(r.sourceSegment->b);out<<']';
        } else out<<"null";
        out<<'}'<<(i+1U<records.size()?",":"")<<'\n';
    }
    out<<"  ]\n}\n";
    return out.str();
}

} // namespace cartmesh2d
