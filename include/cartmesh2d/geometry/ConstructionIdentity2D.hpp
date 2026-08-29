#pragma once

#include "cartmesh2d/geometry/Geometry2D.hpp"

#include <compare>
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <vector>

namespace cartmesh2d {

using StableVertexId2D = std::uint64_t;

enum class StableVertexKeyKind2D {
    LegacyCanonical,
    GridVertex,
    SourceVertex,
    WallGridIntersection,
    TransitionVertex,
    PatchGenerated
};

// A compact typed key.  The meaning of the four words is fixed by kind; the
// stable id is allocated once and is never recomputed from floating point
// coordinates.  R1-A uses LegacyCanonical in shadow mode.  Later R1 stages
// create the more specific keys before geometry is committed.
struct StableVertexKey2D {
    StableVertexKeyKind2D kind = StableVertexKeyKind2D::LegacyCanonical;
    std::uint64_t object = 0;
    std::uint64_t subEntity = 0;
    std::uint64_t logical0 = 0;
    std::uint64_t logical1 = 0;
    auto operator<=>(const StableVertexKey2D&) const = default;
};

enum class ConstructionSourceKind2D {
    Unknown,
    WallSegment,
    TransitionSegment,
    CartesianGrid,
    HybridSourceCell,
    PatchTransaction
};

struct SourceRef2D {
    ConstructionSourceKind2D kind = ConstructionSourceKind2D::Unknown;
    std::uint64_t object = 0;
    std::uint64_t subEntity = 0;
    double parameterBegin = 0.0;
    double parameterEnd = 1.0;
    int side = 0;
};

enum class ConstructionFeatureClass2D {
    Unclassified,
    Smooth,
    ConvexSharp,
    ConcaveSharp,
    GapSideA,
    GapSideB,
    DomainCorner,
    Grid,
    TransitionFixed,
    TransitionMutable
};

struct FeatureOwner2D {
    ConstructionSourceKind2D kind = ConstructionSourceKind2D::Unknown;
    std::uint64_t object = 0;
    std::uint64_t subEntity = 0;
    auto operator<=>(const FeatureOwner2D&) const = default;
};

enum class ConstructionDecision2D {
    ShadowOnly,
    ReusedExactKey,
    Accepted,
    Rejected,
    Refine,
    Resample,
    Rephase
};

struct ConstructionVertexRecord2D {
    StableVertexId2D id = 0;
    StableVertexKey2D key;
    Point2D originalPosition;
    Point2D position;
    double localH = 0.0;
    ConstructionFeatureClass2D featureClass =
        ConstructionFeatureClass2D::Unclassified;
    std::optional<FeatureOwner2D> featureOwner;
    std::vector<SourceRef2D> sourceRefs;
    double displacement = 0.0;
    ConstructionDecision2D decision = ConstructionDecision2D::ShadowOnly;
    std::string decisionReason = "r1a_shadow";
    std::uint64_t creationRevision = 0;
};

struct FeatureVertexIndexProfile2D {
    std::size_t insertionCount = 0;
    std::size_t queryCount = 0;
    std::size_t examinedCandidateCount = 0;
    std::size_t maximumQueryCandidateCount = 0;
};

// Deterministic multilevel bins.  Each local-h exponent has its own bucket
// scale, so a fine resolved feature is not put in a coarse global epsilon bin.
// Query results are sorted by stable id and therefore do not depend on map or
// insertion traversal order.
class FeatureVertexIndex2D {
public:
    void insert(StableVertexId2D id, const Point2D& point, double localH,
                std::size_t supportId);
    [[nodiscard]] std::vector<StableVertexId2D> query(
        const Point2D& point, double radius, std::size_t supportId) const;
    [[nodiscard]] const FeatureVertexIndexProfile2D& profile() const noexcept {
        return profile_;
    }

private:
    struct BucketKey {
        std::size_t supportId = 0;
        int level = 0;
        std::int64_t x = 0;
        std::int64_t y = 0;
        auto operator<=>(const BucketKey&) const = default;
    };

    std::map<BucketKey, std::vector<StableVertexId2D>> buckets_;
    std::map<std::size_t, std::set<int>> levelsBySupport_;
    mutable FeatureVertexIndexProfile2D profile_;
};

class ConstructionVertexStore2D {
public:
    [[nodiscard]] StableVertexId2D addShadowVertex(
        const Point2D& point, double localH,
        ConstructionFeatureClass2D featureClass,
        std::optional<FeatureOwner2D> featureOwner,
        std::size_t supportId);
    void updateMetadata(StableVertexId2D id, double localH,
                        ConstructionFeatureClass2D featureClass,
                        std::optional<FeatureOwner2D> featureOwner);
    void addSourceRef(StableVertexId2D id, SourceRef2D sourceRef);
    [[nodiscard]] std::vector<StableVertexId2D> query(
        const Point2D& point, double radius, std::size_t supportId) const;

    [[nodiscard]] const std::vector<ConstructionVertexRecord2D>& records() const noexcept {
        return records_;
    }
    [[nodiscard]] const FeatureVertexIndexProfile2D& indexProfile() const noexcept {
        return index_.profile();
    }

private:
    std::vector<ConstructionVertexRecord2D> records_;
    FeatureVertexIndex2D index_;
};

} // namespace cartmesh2d
