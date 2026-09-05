#include "cartmesh2d/quadtree/Quadtree2D.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>
namespace cartmesh2d {
namespace {
struct LatticeFace { std::size_t leafIndex = 0; std::uint64_t begin = 0; std::uint64_t end = 0; };
using FaceMap = std::map<std::uint64_t, std::vector<LatticeFace>>;
std::uint64_t pow2(std::size_t exponent) { if (exponent >= 63) throw std::invalid_argument("quadtree level exceeds lattice capacity"); return std::uint64_t{1} << exponent; }
bool intervalOverlapPositive(std::uint64_t a0, std::uint64_t a1, std::uint64_t b0, std::uint64_t b1) noexcept { return std::max(a0,b0) < std::min(a1,b1); }
bool boxesOverlapPositive(const AABB2D& a, const AABB2D& b) noexcept {
    return std::max(a.min.x, b.min.x) < std::min(a.max.x, b.max.x) &&
           std::max(a.min.y, b.min.y) < std::min(a.max.y, b.max.y);
}
void appendMatches(const FaceMap& positiveFaces, const FaceMap& negativeFaces,
                   std::vector<FaceNeighborPair2D>& result) {
    for (const auto& [coordinate, positives] : positiveFaces) {
        const auto it = negativeFaces.find(coordinate);
        if (it == negativeFaces.end()) continue;
        auto lhs = positives;
        auto rhs = it->second;
        const auto byInterval = [](const LatticeFace& a, const LatticeFace& b) {
            return std::tie(a.begin, a.end, a.leafIndex) <
                   std::tie(b.begin, b.end, b.leafIndex);
        };
        std::sort(lhs.begin(), lhs.end(), byInterval);
        std::sort(rhs.begin(), rhs.end(), byInterval);
        std::size_t i = 0;
        std::size_t j = 0;
        while (i < lhs.size() && j < rhs.size()) {
            const auto& a = lhs[i];
            const auto& b = rhs[j];
            if (a.leafIndex != b.leafIndex &&
                intervalOverlapPositive(a.begin, a.end, b.begin, b.end)) {
                result.push_back({std::min(a.leafIndex,b.leafIndex),
                                  std::max(a.leafIndex,b.leafIndex)});
            }
            if (a.end <= b.end) ++i;
            if (b.end <= a.end) ++j;
        }
    }
}
}
Point2D QuadtreeLeaf2D::center() const noexcept { return {(bounds.min.x+bounds.max.x)*0.5,(bounds.min.y+bounds.max.y)*0.5}; }
double QuadtreeLeaf2D::area() const noexcept { return (bounds.max.x-bounds.min.x)*(bounds.max.y-bounds.min.y); }
double distancePointToSegment(const Point2D& point, const Segment2D& segment) noexcept {
    const Vector2D ab = segment.b-segment.a; const double denom = squaredNorm(ab);
    if (denom <= 0.0) return std::sqrt(squaredNorm(point-segment.a));
    const double t = std::clamp(dot(point-segment.a,ab)/denom,0.0,1.0); const Point2D closest = segment.a + ab*t;
    return std::sqrt(squaredNorm(point-closest));
}
double distanceAABBToBoundary(const AABB2D& box, const BoundaryLoop& boundary, const TolerancePolicy& tol) {
    const auto& vertices = boundary.vertices(); if (vertices.size()<2) return std::numeric_limits<double>::infinity();
    const Point2D corners[4]={{box.min.x,box.min.y},{box.max.x,box.min.y},{box.max.x,box.max.y},{box.min.x,box.max.y}};
    double best=std::numeric_limits<double>::infinity();
    for(std::size_t i=0;i<vertices.size();++i){ const Segment2D edge{vertices[i],vertices[(i+1)%vertices.size()]};
        if(segmentIntersectsClosedAABB(edge,box,tol)) return 0.0;
        for(const auto& corner:corners) best=std::min(best,distancePointToSegment(corner,edge));
        for(const auto& endpoint:{edge.a,edge.b}){ const double dx=endpoint.x<box.min.x?box.min.x-endpoint.x:endpoint.x>box.max.x?endpoint.x-box.max.x:0.0; const double dy=endpoint.y<box.min.y?box.min.y-endpoint.y:endpoint.y>box.max.y?endpoint.y-box.max.y:0.0; best=std::min(best,std::hypot(dx,dy)); }
    } return best;
}
double distanceAABBToBoundary(const AABB2D& box, const BoundaryRegion2D& boundary,
                              const TolerancePolicy& tol) {
    double best=std::numeric_limits<double>::infinity();
    for (const auto& loop : boundary.loops()) {
        best=std::min(best,distanceAABBToBoundary(box,loop,tol));
    }
    return best;
}
Quadtree2D::Quadtree2D(Domain2D domain,std::size_t maxLevel,const BoundaryLoop& boundary,const TolerancePolicy& tol)
    :Quadtree2D(domain,maxLevel,BoundaryRegion2D(boundary),tol){}
Quadtree2D::Quadtree2D(Domain2D domain,std::size_t maxLevel,const BoundaryRegion2D& boundary,const TolerancePolicy& tol):domain_(domain),maxLevel_(maxLevel),boundaryIndex_(boundary,tol){
    if (!domain_.valid(tol)) throw std::invalid_argument("invalid quadtree domain");
    if (maxLevel_ > 28) throw std::invalid_argument("maxLevel > 28 unsupported");
    if (!boundary.diagnose(tol).valid()) throw std::invalid_argument("invalid boundary region");
    if (!boundaryIndex_.valid()) throw std::invalid_argument("invalid boundary segment index");
    leaves_.push_back(makeLeaf(0, 0, 0, boundary, tol));
    sortAndAssignIds();
}
std::uint64_t Quadtree2D::mortonPath(std::size_t level,std::uint64_t ix,std::uint64_t iy) noexcept { std::uint64_t code=0; for(std::size_t bit=0;bit<level;++bit){ code|=((ix>>bit)&1ULL)<<(2*bit); code|=((iy>>bit)&1ULL)<<(2*bit+1);} return code; }
std::uint64_t Quadtree2D::makeKey(std::size_t level,std::uint64_t ix,std::uint64_t iy) noexcept { return (mortonPath(level,ix,iy)<<6)|static_cast<std::uint64_t>(level); }
AABB2D Quadtree2D::boundsFor(std::size_t level,std::uint64_t ix,std::uint64_t iy) const noexcept {
    const std::uint64_t count = std::uint64_t{1} << level;
    const double n = static_cast<double>(count);
    const double dx = domain_.width() / n;
    const double dy = domain_.height() / n;
    const auto xAt = [&](std::uint64_t i) noexcept {
        if (i == count) return domain_.bounds.max.x;
        return domain_.bounds.min.x + static_cast<double>(i) * dx;
    };
    const auto yAt = [&](std::uint64_t j) noexcept {
        if (j == count) return domain_.bounds.max.y;
        return domain_.bounds.min.y + static_cast<double>(j) * dy;
    };
    return {{xAt(ix), yAt(iy)}, {xAt(ix + 1), yAt(iy + 1)}};
}
QuadtreeLeaf2D Quadtree2D::makeLeaf(std::size_t level,std::uint64_t ix,std::uint64_t iy,const BoundaryRegion2D& boundary,const TolerancePolicy& tol) const { (void)boundary; QuadtreeLeaf2D leaf; leaf.key=makeKey(level,ix,iy); leaf.level=level; leaf.ix=ix; leaf.iy=iy; leaf.bounds=boundsFor(level,ix,iy); if(boundaryIndex_.intersects(leaf.bounds,tol)) leaf.classification=CellClass::Intersected; else { const auto state=boundaryIndex_.classifyPoint(leaf.center(),tol); leaf.classification=state==PointInPolygon::Inside?CellClass::Inside:state==PointInPolygon::Boundary?CellClass::Intersected:CellClass::Outside; } return leaf; }
bool Quadtree2D::splitLeafAt(std::size_t index,const BoundaryRegion2D& boundary,const TolerancePolicy& tol){ if(index>=leaves_.size()) return false; const auto parent=leaves_[index]; if(parent.level>=maxLevel_) return false; const std::size_t l=parent.level+1; const std::uint64_t x=parent.ix*2,y=parent.iy*2; QuadtreeLeaf2D children[4]={makeLeaf(l,x,y,boundary,tol),makeLeaf(l,x+1,y,boundary,tol),makeLeaf(l,x,y+1,boundary,tol),makeLeaf(l,x+1,y+1,boundary,tol)}; leaves_.erase(leaves_.begin()+static_cast<std::ptrdiff_t>(index)); leaves_.insert(leaves_.end(),std::begin(children),std::end(children)); sortAndAssignIds(); return true; }
bool Quadtree2D::refineLeafByKey(std::uint64_t key,const BoundaryLoop& boundary,const TolerancePolicy& tol){ return refineLeafByKey(key,BoundaryRegion2D(boundary),tol); }
bool Quadtree2D::refineLeafByKey(std::uint64_t key,const BoundaryRegion2D& boundary,const TolerancePolicy& tol){ if(!boundaryIndex_.matches(boundary)) throw std::invalid_argument("refinement boundary differs from indexed boundary"); const auto it=std::find_if(leaves_.begin(),leaves_.end(),[key](const auto& leaf){return leaf.key==key;}); if(it==leaves_.end()) return false; return splitLeafAt(static_cast<std::size_t>(std::distance(leaves_.begin(),it)),boundary,tol); }
QuadtreeLocalRefinementReport2D Quadtree2D::refineLeavesWithClosure(
    std::vector<std::uint64_t> keys,const BoundaryLoop& boundary,
    const TolerancePolicy& tol) {
    return refineLeavesWithClosure(std::move(keys),BoundaryRegion2D(boundary),tol);
}
QuadtreeLocalRefinementReport2D Quadtree2D::refineLeavesWithClosure(
    std::vector<std::uint64_t> keys,const BoundaryRegion2D& boundary,
    const TolerancePolicy& tol) {
    if (!boundaryIndex_.matches(boundary)) {
        throw std::invalid_argument("local refinement boundary differs from indexed boundary");
    }
    QuadtreeLocalRefinementReport2D report;
    std::sort(keys.begin(),keys.end());
    keys.erase(std::unique(keys.begin(),keys.end()),keys.end());
    report.requestedKeys=keys;
    const auto countViolations=[&]() {
        std::size_t count=0;
        for (const auto& pair:faceNeighbors()) {
            const auto a=leaves_[pair.first].level,b=leaves_[pair.second].level;
            if ((a>b?a-b:b-a)>1U) ++count;
        }
        return count;
    };
    report.violationsBefore=countViolations();
    for (const auto key:keys) {
        const auto leaf=std::find_if(leaves_.begin(),leaves_.end(),
            [key](const auto& item){return item.key==key;});
        if (leaf==leaves_.end() || leaf->level>=maxLevel_)
            report.rejectedKeys.push_back(key);
    }
    if (!report.rejectedKeys.empty()) {
        report.violationsAfter=report.violationsBefore;
        return report;
    }
    const auto splitBatch=[&](const std::vector<std::uint64_t>& batch,bool closure) {
        std::vector<QuadtreeLeaf2D> next;
        next.reserve(leaves_.size()+3U*batch.size());
        for (const auto& leaf:leaves_) {
            if (!std::binary_search(batch.begin(),batch.end(),leaf.key)) {
                next.push_back(leaf);continue;
            }
            if (leaf.level>=maxLevel_) {
                next.push_back(leaf);continue;
            }
            const auto childLevel=leaf.level+1U;
            const auto childX=leaf.ix*2U,childY=leaf.iy*2U;
            std::array<QuadtreeLeaf2D,4> children{
                makeLeaf(childLevel,childX,childY,boundary,tol),
                makeLeaf(childLevel,childX+1U,childY,boundary,tol),
                makeLeaf(childLevel,childX,childY+1U,boundary,tol),
                makeLeaf(childLevel,childX+1U,childY+1U,boundary,tol)};
            QuadtreeRefinementLineage2D lineage;
            lineage.parentKey=leaf.key;lineage.parentLevel=leaf.level;lineage.closure=closure;
            for (std::size_t i=0;i<children.size();++i) {
                lineage.childKeys[i]=children[i].key;
                next.push_back(std::move(children[i]));
            }
            report.lineage.push_back(lineage);
            if (closure) ++report.closureRefinedLeaves;
            else ++report.requestedRefinedLeaves;
        }
        leaves_.swap(next);sortAndAssignIds();
    };
    splitBatch(keys,false);
    while (true) {
        std::vector<std::uint64_t> closureKeys;
        for (const auto& pair:faceNeighbors()) {
            const auto& a=leaves_[pair.first];const auto& b=leaves_[pair.second];
            const auto difference=a.level>b.level?a.level-b.level:b.level-a.level;
            if (difference>1U) closureKeys.push_back((a.level<b.level?a:b).key);
        }
        if (closureKeys.empty()) break;
        std::sort(closureKeys.begin(),closureKeys.end());
        closureKeys.erase(std::unique(closureKeys.begin(),closureKeys.end()),closureKeys.end());
        splitBatch(closureKeys,true);
        if (++report.closureIterations>maxLevel_+1U) {
            throw std::runtime_error("local 2:1 refinement closure failed to converge");
        }
    }
    std::sort(report.rejectedKeys.begin(),report.rejectedKeys.end());
    report.rejectedKeys.erase(std::unique(report.rejectedKeys.begin(),report.rejectedKeys.end()),
                              report.rejectedKeys.end());
    report.violationsAfter=countViolations();
    return report;
}
void Quadtree2D::refine(const BoundaryLoop& boundary,
                        const QuadtreeRefinementPolicy2D& policy,
                        const TolerancePolicy& tol) {
    refine(BoundaryRegion2D(boundary),policy,tol);
}
void Quadtree2D::refine(const BoundaryRegion2D& boundary,
                        const QuadtreeRefinementPolicy2D& policy,
                        const TolerancePolicy& tol) {
    if (!boundary.diagnose(tol).valid()) throw std::invalid_argument("invalid boundary loop");
    if (!boundaryIndex_.matches(boundary)) throw std::invalid_argument("refinement boundary differs from indexed boundary");
    if (policy.minimumLevel > maxLevel_) throw std::invalid_argument("minimumLevel exceeds maxLevel");
    if (policy.boundaryLevel > maxLevel_) throw std::invalid_argument("boundaryLevel exceeds maxLevel");
    for (const auto& band : policy.distanceBands) {
        if (!std::isfinite(band.distance) || band.distance < 0.0) {
            throw std::invalid_argument("invalid distance refinement band");
        }
        if (band.targetLevel > maxLevel_) throw std::invalid_argument("distance level exceeds maxLevel");
    }
    for (const auto& region : policy.boxRegions) {
        if (!region.bounds.valid(tol)) {
            throw std::invalid_argument("invalid box refinement region");
        }
        if (region.targetLevel == 0U || region.targetLevel > maxLevel_) {
            throw std::invalid_argument("box refinement level outside [1,maxLevel]");
        }
        if (!boxesOverlapPositive(region.bounds, domain_.bounds)) {
            throw std::invalid_argument("box refinement region does not overlap domain");
        }
    }
    for (const auto& band : policy.segmentBands) {
        if (!std::isfinite(band.radius) || band.radius < 0.0) {
            throw std::invalid_argument("invalid segment refinement radius");
        }
        if (band.targetLevel > maxLevel_) {
            throw std::invalid_argument("segment level exceeds maxLevel");
        }
        // Ids address the indexed boundary, so a stale or unsorted list would
        // silently select the wrong stretch of wall instead of failing.
        if (!std::is_sorted(band.segmentIds.begin(), band.segmentIds.end()) ||
            std::adjacent_find(band.segmentIds.begin(), band.segmentIds.end()) !=
                band.segmentIds.end()) {
            throw std::invalid_argument("segment refinement ids must be sorted and unique");
        }
        if (!band.segmentIds.empty() &&
            band.segmentIds.back() >= boundaryIndex_.segmentCount()) {
            throw std::invalid_argument("segment refinement id outside indexed boundary");
        }
    }

    // Refine every eligible leaf once per level pass, then sort only after the
    // complete batch.  The former implementation split one leaf, sorted the
    // whole tree and restarted from index zero, making boundary refinement
    // quadratic even when the boundary segment count was fixed.
    while (true) {
        bool changed = false;
        std::vector<QuadtreeLeaf2D> next;
        next.reserve(leaves_.size() * 2U);
        for (const auto& leaf : leaves_) {
            std::size_t requested = policy.minimumLevel;
            if (leaf.classification == CellClass::Intersected) {
                requested = std::max(requested, policy.boundaryLevel);
            }
            if (!policy.distanceBands.empty()) {
                const double distance = boundaryIndex_.distanceToAABB(leaf.bounds, tol);
                for (const auto& band : policy.distanceBands) {
                    if (distance <= band.distance + tol.scale(band.distance)) {
                        requested = std::max(requested, band.targetLevel);
                    }
                }
            }
            for (const auto& region : policy.boxRegions) {
                if (boxesOverlapPositive(leaf.bounds, region.bounds)) {
                    requested = std::max(requested, region.targetLevel);
                }
            }
            for (const auto& band : policy.segmentBands) {
                // Nothing this band can add, so skip the index query entirely.
                // This is also what keeps an empty segmentBands list free.
                if (requested >= band.targetLevel || band.segmentIds.empty()) continue;
                const AABB2D grown{
                    {leaf.bounds.min.x - band.radius, leaf.bounds.min.y - band.radius},
                    {leaf.bounds.max.x + band.radius, leaf.bounds.max.y + band.radius}};
                for (const auto id : boundaryIndex_.querySegmentIds(grown, tol)) {
                    if (std::binary_search(band.segmentIds.begin(),
                                           band.segmentIds.end(), id)) {
                        requested = band.targetLevel;
                        break;
                    }
                }
            }

            if (leaf.level >= requested) {
                next.push_back(leaf);
                continue;
            }

            changed = true;
            const std::size_t childLevel = leaf.level + 1U;
            const std::uint64_t childX = leaf.ix * 2U;
            const std::uint64_t childY = leaf.iy * 2U;
            next.push_back(makeLeaf(childLevel, childX, childY, boundary, tol));
            next.push_back(makeLeaf(childLevel, childX + 1U, childY, boundary, tol));
            next.push_back(makeLeaf(childLevel, childX, childY + 1U, boundary, tol));
            next.push_back(makeLeaf(childLevel, childX + 1U, childY + 1U, boundary, tol));
        }
        leaves_.swap(next);
        if (!changed) break;
    }
    sortAndAssignIds();
}
std::vector<FaceNeighborPair2D> Quadtree2D::faceNeighbors() const {
    FaceMap left,right,bottom,top;
    for(std::size_t index=0;index<leaves_.size();++index){
        const auto& leaf=leaves_[index];
        const std::uint64_t scale=pow2(maxLevel_-leaf.level);
        const std::uint64_t x0=leaf.ix*scale,x1=(leaf.ix+1)*scale;
        const std::uint64_t y0=leaf.iy*scale,y1=(leaf.iy+1)*scale;
        left[x0].push_back({index,y0,y1});
        right[x1].push_back({index,y0,y1});
        bottom[y0].push_back({index,x0,x1});
        top[y1].push_back({index,x0,x1});
    }
    std::vector<FaceNeighborPair2D> result;
    appendMatches(right,left,result);
    appendMatches(top,bottom,result);
    std::sort(result.begin(),result.end(),[](const auto&a,const auto&b){return std::tie(a.first,a.second)<std::tie(b.first,b.second);});
    result.erase(std::unique(result.begin(),result.end(),[](const auto&a,const auto&b){return a.first==b.first&&a.second==b.second;}),result.end());
    return result;
}
std::size_t Quadtree2D::countBalanceViolations() const { std::size_t count=0; for(const auto& p:faceNeighbors()){ const auto a=leaves_[p.first].level,b=leaves_[p.second].level; if((a>b?a-b:b-a)>1) ++count; } return count; }
QuadtreeBalanceReport2D Quadtree2D::enforceTwoToOneBalance(
    const BoundaryLoop& boundary, const TolerancePolicy& tol) {
    return enforceTwoToOneBalance(BoundaryRegion2D(boundary),tol);
}
QuadtreeBalanceReport2D Quadtree2D::enforceTwoToOneBalance(
    const BoundaryRegion2D& boundary, const TolerancePolicy& tol) {
    if (!boundaryIndex_.matches(boundary)) {
        throw std::invalid_argument("balance boundary differs from indexed boundary");
    }
    QuadtreeBalanceReport2D report;
    const auto timedFaceNeighbors = [&]() {
        const auto start = std::chrono::steady_clock::now();
        auto neighbors = faceNeighbors();
        report.faceNeighborSeconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        ++report.faceNeighborCalls;
        return neighbors;
    };
    const auto countViolations = [&](const std::vector<FaceNeighborPair2D>& neighbors) {
        std::size_t count = 0;
        for (const auto& pair : neighbors) {
            const auto a = leaves_[pair.first].level;
            const auto b = leaves_[pair.second].level;
            if ((a > b ? a - b : b - a) > 1U) ++count;
        }
        return count;
    };
    auto neighbors = timedFaceNeighbors();
    report.violationsBefore = countViolations(neighbors);
    while (true) {
        std::vector<std::uint64_t> keys;
        for (const auto& pair : neighbors) {
            const auto& a = leaves_[pair.first];
            const auto& b = leaves_[pair.second];
            const auto difference = a.level > b.level ? a.level - b.level : b.level - a.level;
            if (difference > 1U) keys.push_back((a.level < b.level ? a : b).key);
        }
        if (keys.empty()) break;
        std::sort(keys.begin(), keys.end());
        keys.erase(std::unique(keys.begin(), keys.end()), keys.end());

        std::vector<QuadtreeLeaf2D> next;
        next.reserve(leaves_.size() + 3U * keys.size());
        for (const auto& leaf : leaves_) {
            if (!std::binary_search(keys.begin(), keys.end(), leaf.key)) {
                next.push_back(leaf);
                continue;
            }
            if (leaf.level >= maxLevel_) {
                next.push_back(leaf);
                continue;
            }
            const std::size_t childLevel = leaf.level + 1U;
            const std::uint64_t childX = leaf.ix * 2U;
            const std::uint64_t childY = leaf.iy * 2U;
            next.push_back(makeLeaf(childLevel, childX, childY, boundary, tol));
            next.push_back(makeLeaf(childLevel, childX + 1U, childY, boundary, tol));
            next.push_back(makeLeaf(childLevel, childX, childY + 1U, boundary, tol));
            next.push_back(makeLeaf(childLevel, childX + 1U, childY + 1U, boundary, tol));
            ++report.refinedLeaves;
        }
        leaves_.swap(next);
        sortAndAssignIds();
        ++report.iterations;
        if (report.iterations > maxLevel_ + 1U) {
            throw std::runtime_error("2:1 balance failed to converge");
        }
        neighbors = timedFaceNeighbors();
    }
    report.violationsAfter = countViolations(neighbors);
    return report;
}
double Quadtree2D::totalLeafArea() const noexcept { double total=0.0; for(const auto& leaf:leaves_) total+=leaf.area(); return total; }
bool Quadtree2D::deterministicOrderingValid() const noexcept { for(std::size_t i=0;i<leaves_.size();++i) if(leaves_[i].id!=i) return false; return true; }
void Quadtree2D::sortAndAssignIds(){ std::sort(leaves_.begin(),leaves_.end(),[this](const auto&a,const auto&b){ const std::uint64_t as=std::uint64_t{1}<<(2*(maxLevel_-a.level)),bs=std::uint64_t{1}<<(2*(maxLevel_-b.level)); const auto aa=mortonPath(a.level,a.ix,a.iy)*as,ba=mortonPath(b.level,b.ix,b.iy)*bs; if(aa!=ba) return aa<ba; if(a.level!=b.level) return a.level<b.level; return a.key<b.key; }); for(std::size_t i=0;i<leaves_.size();++i) leaves_[i].id=i; }
} // namespace cartmesh2d
