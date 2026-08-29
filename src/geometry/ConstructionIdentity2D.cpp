#include "cartmesh2d/geometry/ConstructionIdentity2D.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <tuple>

namespace cartmesh2d {
namespace {

[[nodiscard]] int scaleLevel(double localH) {
    int exponent = 0;
    (void)std::frexp(localH, &exponent);
    return exponent - 1;
}

[[nodiscard]] std::int64_t bucketCoordinate(double value, double width) {
    const long double coordinate =
        std::floor(static_cast<long double>(value) / static_cast<long double>(width));
    if (coordinate < static_cast<long double>(std::numeric_limits<std::int64_t>::min()) ||
        coordinate > static_cast<long double>(std::numeric_limits<std::int64_t>::max())) {
        throw std::overflow_error("feature vertex spatial bucket coordinate overflow");
    }
    return static_cast<std::int64_t>(coordinate);
}

[[nodiscard]] bool sameSourceRef(const SourceRef2D& a, const SourceRef2D& b) noexcept {
    return std::tie(a.kind, a.object, a.subEntity, a.parameterBegin,
                    a.parameterEnd, a.side) ==
           std::tie(b.kind, b.object, b.subEntity, b.parameterBegin,
                    b.parameterEnd, b.side);
}

} // namespace

void FeatureVertexIndex2D::insert(StableVertexId2D id, const Point2D& point,
                                  double localH, std::size_t supportId) {
    if (!std::isfinite(point.x) || !std::isfinite(point.y) ||
        !std::isfinite(localH) || !(localH > 0.0)) {
        throw std::invalid_argument("feature vertex index requires finite point and local_h");
    }
    const int level = scaleLevel(localH);
    const double width = std::ldexp(1.0, level);
    BucketKey key{supportId, level,
                  bucketCoordinate(point.x, width),
                  bucketCoordinate(point.y, width)};
    auto& bucket = buckets_[key];
    const auto position = std::lower_bound(bucket.begin(), bucket.end(), id);
    if (position == bucket.end() || *position != id) bucket.insert(position, id);
    levelsBySupport_[supportId].insert(level);
    ++profile_.insertionCount;
}

std::vector<StableVertexId2D> FeatureVertexIndex2D::query(
    const Point2D& point, double radius, std::size_t supportId) const {
    if (!std::isfinite(point.x) || !std::isfinite(point.y) ||
        !std::isfinite(radius) || radius < 0.0) {
        throw std::invalid_argument("feature vertex query requires finite point and radius");
    }
    std::vector<StableVertexId2D> result;
    const auto levels = levelsBySupport_.find(supportId);
    if (levels != levelsBySupport_.end()) {
        for (const int level : levels->second) {
            const double width = std::ldexp(1.0, level);
            // localH in this level is in [width,2*width), while the registry
            // contract requires snapFraction<0.5.  Therefore no admissible
            // anchor in this level can lie farther than width, even if the
            // query comes from a much coarser cell.
            const double levelRadius=std::min(radius,width);
            const auto minX = bucketCoordinate(point.x - levelRadius, width);
            const auto maxX = bucketCoordinate(point.x + levelRadius, width);
            const auto minY = bucketCoordinate(point.y - levelRadius, width);
            const auto maxY = bucketCoordinate(point.y + levelRadius, width);
            for (std::int64_t x = minX; x <= maxX; ++x) {
                for (std::int64_t y = minY; y <= maxY; ++y) {
                    const auto bucket = buckets_.find({supportId, level, x, y});
                    if (bucket == buckets_.end()) continue;
                    result.insert(result.end(), bucket->second.begin(), bucket->second.end());
                    if (y == std::numeric_limits<std::int64_t>::max()) break;
                }
                if (x == std::numeric_limits<std::int64_t>::max()) break;
            }
        }
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    ++profile_.queryCount;
    profile_.examinedCandidateCount += result.size();
    profile_.maximumQueryCandidateCount =
        std::max(profile_.maximumQueryCandidateCount, result.size());
    return result;
}

StableVertexId2D ConstructionVertexStore2D::addShadowVertex(
    const Point2D& point, double localH,
    ConstructionFeatureClass2D featureClass,
    std::optional<FeatureOwner2D> featureOwner,
    std::size_t supportId) {
    const auto id = static_cast<StableVertexId2D>(records_.size());
    ConstructionVertexRecord2D record;
    record.id = id;
    record.key = {StableVertexKeyKind2D::LegacyCanonical,
                  static_cast<std::uint64_t>(supportId), id, 0, 0};
    record.originalPosition = point;
    record.position = point;
    record.localH = localH;
    record.featureClass = featureClass;
    record.featureOwner = featureOwner;
    record.creationRevision = id;
    records_.push_back(std::move(record));
    index_.insert(id, point, localH, supportId);
    return id;
}

void ConstructionVertexStore2D::updateMetadata(
    StableVertexId2D id, double localH,
    ConstructionFeatureClass2D featureClass,
    std::optional<FeatureOwner2D> featureOwner) {
    if (id >= records_.size() || !std::isfinite(localH) || !(localH > 0.0)) {
        throw std::invalid_argument("invalid construction vertex metadata update");
    }
    auto& record = records_[static_cast<std::size_t>(id)];
    record.localH = std::min(record.localH, localH);
    record.featureClass = featureClass;
    if (featureOwner) record.featureOwner = featureOwner;
}

void ConstructionVertexStore2D::addSourceRef(
    StableVertexId2D id, SourceRef2D sourceRef) {
    if (id >= records_.size()) {
        throw std::out_of_range("construction source reference has invalid stable vertex id");
    }
    auto& refs = records_[static_cast<std::size_t>(id)].sourceRefs;
    if (std::none_of(refs.begin(), refs.end(),
                     [&](const SourceRef2D& ref) { return sameSourceRef(ref, sourceRef); })) {
        refs.push_back(sourceRef);
        std::sort(refs.begin(), refs.end(), [](const SourceRef2D& a, const SourceRef2D& b) {
            return std::tie(a.kind, a.object, a.subEntity, a.parameterBegin,
                            a.parameterEnd, a.side) <
                   std::tie(b.kind, b.object, b.subEntity, b.parameterBegin,
                            b.parameterEnd, b.side);
        });
    }
}

std::vector<StableVertexId2D> ConstructionVertexStore2D::query(
    const Point2D& point, double radius, std::size_t supportId) const {
    return index_.query(point, radius, supportId);
}

} // namespace cartmesh2d
