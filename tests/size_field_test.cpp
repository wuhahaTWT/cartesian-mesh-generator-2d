#include "cartmesh2d/sizing/SizeField2D.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
using namespace cartmesh2d;
namespace {
int failures=0;
void check(bool c,const std::string&m){if(!c){++failures;std::cerr<<"FAIL: "<<m<<'\n';}}
void near(double a,double b,double e,const std::string&m){check(std::abs(a-b)<=e,m);}

// Unit circle sampled finely enough that the discrete curvature radius is close
// to the true radius.
BoundaryLoop circleLoop(std::size_t count,double radius,Point2D centre={0.0,0.0}) {
    std::vector<Point2D> vertices;
    vertices.reserve(count);
    const double pi=std::acos(-1.0);
    for (std::size_t i=0;i<count;++i) {
        const double angle=2.0*pi*static_cast<double>(i)/static_cast<double>(count);
        vertices.push_back({centre.x+radius*std::cos(angle),centre.y+radius*std::sin(angle)});
    }
    return BoundaryLoop(std::move(vertices));
}
} // namespace

int main(){
    // sizeFieldLevelForSize2D is the one arithmetic primitive every source uses.
    check(sizeFieldLevelForSize2D(1.0,1.0,28)==0,"target equal to the domain needs no refinement");
    check(sizeFieldLevelForSize2D(1.0,0.5,28)==1,"exact power of two lands on its own level");
    check(sizeFieldLevelForSize2D(1.0,0.3,28)==2,"non-power-of-two rounds to the finer level");
    check(sizeFieldLevelForSize2D(1.0,1e-30,4)==4,"unreachable target saturates at the cap");
    check(sizeFieldLevelForSize2D(1.0,-1.0,28)==0,"non-positive target requests nothing");

    const BoundaryRegion2D circle(circleLoop(64,1.0));

    // A square domain is the whole point: any other aspect makes every cell in the
    // mesh non-square, because Quadtree2D divides width and height by 2^level.
    {
        SizeFieldPolicy2D policy;
        policy.farFieldSpans=10.0;
        policy.wallCellsPerSpan=128.0;
        policy.allowUnsafeWallLevel=true;
        const auto resolved=resolveSizeField2D(policy,circle);
        check(resolved.valid(),"default policy resolves");
        near(resolved.bodySpan,2.0,1e-12,"body span is the bounding-box span");
        near(resolved.domainSpan,2.0*21.0,1e-12,"domain span is span*(1+2*farFieldSpans)");
        near(resolved.domain.width(),resolved.domain.height(),1e-12,"resolved domain is square");
        near(resolved.domain.width(),resolved.domainSpan,1e-12,"domain width matches domain span");
        check(resolved.wallLevel>0,"wall level is refined");
        // Wall resolution is requested physically, so it must bracket the request.
        check(resolved.wallCellSize<=resolved.bodySpan/(*policy.wallCellsPerSpan)*1.0000001,
              "wall cell is no coarser than requested");
        check(resolved.wallCellSize>=resolved.bodySpan/(*policy.wallCellsPerSpan)*0.5,
              "wall cell is within one level of the request");
    }

    // Domain size and wall resolution must be independent: that is the property the
    // old padding-fraction/max-level pair did not have.
    {
        SizeFieldPolicy2D near_;
        near_.farFieldSpans=1.0;
        near_.wallCellsPerSpan=128.0;
        near_.allowUnsafeWallLevel=true;
        SizeFieldPolicy2D far;
        far.farFieldSpans=20.0;
        far.wallCellsPerSpan=128.0;
        far.allowUnsafeWallLevel=true;
        const auto a=resolveSizeField2D(near_,circle);
        const auto b=resolveSizeField2D(far,circle);
        check(a.valid()&&b.valid(),"both far-field settings resolve");
        check(b.wallLevel>a.wallLevel,"a larger domain needs a deeper tree for the same wall size");
        near(a.wallCellSize,b.wallCellSize,0.5*a.wallCellSize,
             "wall cell size is preserved across far-field settings");
    }

    // cellsPerLevel is the nCellsBetweenLevels analogue.  Zero must emit nothing so
    // the pre-existing one-ring-per-level gradation is reachable unchanged.
    {
        SizeFieldPolicy2D policy;
        policy.wallDistance=WallDistanceSizing2D{0U,0U};
        const auto resolved=resolveSizeField2D(policy,circle);
        check(resolved.valid(),"cellsPerLevel zero resolves");
        check(resolved.refinement.distanceBands.empty(),"cellsPerLevel zero emits no bands");
        check(resolved.refinement.segmentBands.empty(),"no segment sources means no segment bands");
        check(resolved.refinement.boxRegions.empty(),"no wake means no box regions");
        check(resolved.refinement.boundaryLevel==resolved.wallLevel,
              "boundary level is the resolved wall level");
        check(resolved.refinement.minimumLevel==0U,"far level zero leaves no global floor");
    }

    // The band ladder must be the cumulative thickness of cellsPerLevel cells per
    // level, and must descend monotonically so refine()'s max-combining is stable.
    {
        SizeFieldPolicy2D policy;
        policy.farFieldSpans=10.0;
        policy.wallCellsPerSpan=128.0;
        policy.allowUnsafeWallLevel=true;
        policy.wallDistance=WallDistanceSizing2D{3U,2U};
        const auto resolved=resolveSizeField2D(policy,circle);
        check(resolved.valid(),"graded ladder resolves");
        const auto& bands=resolved.refinement.distanceBands;
        check(bands.size()==resolved.wallLevel-2U,"one band per level above farLevel");
        check(!bands.empty(),"ladder is non-empty");
        check(resolved.refinement.minimumLevel==2U,"farLevel becomes the global floor");
        bool descendingLevels=true,increasingDistance=true;
        for (std::size_t i=0;i+1U<bands.size();++i) {
            descendingLevels=descendingLevels&&bands[i].targetLevel>bands[i+1U].targetLevel;
            increasingDistance=increasingDistance&&bands[i+1U].distance>bands[i].distance;
        }
        check(descendingLevels,"levels descend away from the wall");
        check(increasingDistance,"reach grows away from the wall");
        check(bands.front().targetLevel==resolved.wallLevel,"finest band is the wall level");
        const double finest=std::ldexp(resolved.domainSpan,-static_cast<int>(resolved.wallLevel));
        near(bands.front().distance,3.0*finest,1e-12*resolved.domainSpan,
             "finest band reaches cellsPerLevel wall cells");
        // Second band adds cellsPerLevel cells of the next coarser level on top.
        near(bands[1].distance,3.0*finest+3.0*2.0*finest,1e-12*resolved.domainSpan,
             "second band accumulates the coarser level thickness");
    }

    // Curvature may only raise the level, and only where the wall actually turns.
    // A coarse wall request plus a small circle is the case that needs it.
    {
        SizeFieldPolicy2D flat;
        flat.wallCellsPerSpan=512.0;
        flat.allowUnsafeWallLevel=true;
        flat.curvature=CurvatureSizing2D{8.0};
        const auto alreadyFine=resolveSizeField2D(flat,circle);
        check(alreadyFine.valid(),"curvature on an already-fine wall resolves");
        check(alreadyFine.refinement.segmentBands.empty(),
              "curvature adds nothing when the wall level already resolves it");
        check(alreadyFine.curvatureLevel==alreadyFine.wallLevel,
              "curvature level collapses onto the wall level");

        SizeFieldPolicy2D coarse;
        coarse.wallCellsPerSpan=8.0;
        coarse.curvature=CurvatureSizing2D{8.0};
        const auto raised=resolveSizeField2D(coarse,circle);
        check(raised.valid(),"curvature on a coarse wall resolves");
        check(!raised.refinement.segmentBands.empty(),"curvature raises a coarse wall");
        check(raised.curvatureLevel>raised.wallLevel,"curvature level exceeds the wall level");
        check(raised.maxLevel==raised.curvatureLevel,"maxLevel is derived from the sources");
        std::size_t selected=0;
        for (const auto& band:raised.refinement.segmentBands) {
            check(band.targetLevel>raised.wallLevel,"segment band is above the wall level");
            check(std::is_sorted(band.segmentIds.begin(),band.segmentIds.end()),
                  "segment ids are sorted");
            check(std::adjacent_find(band.segmentIds.begin(),band.segmentIds.end())==
                  band.segmentIds.end(),"segment ids are unique");
            check(band.segmentIds.back()<64U,"segment ids address the indexed boundary");
            selected+=band.segmentIds.size();
        }
        check(selected==64U,"a circle turns uniformly so every segment is selected");
    }

    // A straight-sided box has zero curvature on its faces; only the corner triples
    // can request anything, so most segments must be left alone.
    {
        BoundaryRegion2D box(BoundaryLoop({{0.0,0.0},{1.0,0.0},{2.0,0.0},{3.0,0.0},
                                           {3.0,1.0},{3.0,2.0},{0.0,2.0}}));
        SizeFieldPolicy2D policy;
        policy.wallCellsPerSpan=4.0;
        policy.curvature=CurvatureSizing2D{8.0};
        const auto resolved=resolveSizeField2D(policy,box);
        check(resolved.valid(),"box curvature resolves");
        std::size_t selected=0;
        for (const auto& band:resolved.refinement.segmentBands) selected+=band.segmentIds.size();
        check(selected<7U,"collinear wall triples request no curvature refinement");
    }

    // Proximity: two circles separated by a gap far narrower than the wall cell.
    {
        std::vector<BoundaryLoop> loops{circleLoop(32,1.0,{0.0,0.0}),
                                        circleLoop(32,1.0,{2.08,0.0})};
        BoundaryRegion2D pair(std::move(loops));
        SizeFieldPolicy2D policy;
        policy.wallCellsPerSpan=16.0;
        // A 0.08 gap under a 0.19 wall cell needs three extra levels, which is past
        // the trusted ceiling; that is the guard working, not the source failing.
        policy.allowUnsafeWallLevel=true;
        policy.proximity=ProximitySizing2D{4.0,0.5,1U};
        const auto resolved=resolveSizeField2D(policy,pair);
        check(resolved.valid(),"narrow-gap proximity resolves");
        check(resolved.proximityLevel>resolved.wallLevel,"a narrow gap raises the level");
        check(!resolved.refinement.segmentBands.empty(),"proximity emits a segment band");
        std::size_t selected=0;
        for (const auto& band:resolved.refinement.segmentBands) selected+=band.segmentIds.size();
        check(selected>0U&&selected<64U,"only the facing stretches are selected");
    }

    // A single convex loop must not detect itself as a gap.  An index window alone
    // cannot do this: on a 32-gon of radius 1 the connector between segments three
    // apart is only ~0.19 long, well inside any plausible search radius.
    {
        SizeFieldPolicy2D policy;
        policy.wallCellsPerSpan=16.0;
        policy.proximity=ProximitySizing2D{4.0,0.5,1U};
        const auto resolved=resolveSizeField2D(policy,circle);
        check(resolved.valid(),"single-loop proximity resolves");
        check(resolved.proximityLevel==resolved.wallLevel,
              "a convex loop does not report its own chords as a gap");
        check(resolved.refinement.segmentBands.empty(),
              "a convex loop emits no proximity band");
    }

    // The wake is the one direction-dependent source, so it must follow the angle
    // rather than being inferred.  Boxes are axis-aligned, hence the staircase.
    {
        SizeFieldPolicy2D straight;
        straight.wallCellsPerSpan=128.0;
        straight.allowUnsafeWallLevel=true;
        straight.wake=WakeSizing2D{0.0,6.0,0.6,1U,8U};
        const auto aligned=resolveSizeField2D(straight,circle);
        check(aligned.valid(),"aligned wake resolves");
        check(aligned.refinement.boxRegions.size()==8U,"one box per requested slice");
        double minY=1e30,maxY=-1e30,maxX=-1e30;
        for (const auto& region:aligned.refinement.boxRegions) {
            check(region.targetLevel==aligned.wallLevel-1U,"wake sits one level below the wall");
            minY=std::min(minY,region.bounds.min.y);
            maxY=std::max(maxY,region.bounds.max.y);
            maxX=std::max(maxX,region.bounds.max.x);
        }
        near(minY,-0.6*aligned.bodySpan,1e-9,"aligned wake half-width below the axis");
        near(maxY,0.6*aligned.bodySpan,1e-9,"aligned wake half-width above the axis");
        near(maxX,6.0*aligned.bodySpan,1e-9,"aligned wake reaches the requested length");

        SizeFieldPolicy2D angled;
        angled.wallCellsPerSpan=128.0;
        angled.allowUnsafeWallLevel=true;
        angled.wake=WakeSizing2D{15.0,6.0,0.6,1U,8U};
        const auto rotated=resolveSizeField2D(angled,circle);
        check(rotated.valid(),"angled wake resolves");
        double rotatedMaxY=-1e30;
        for (const auto& region:rotated.refinement.boxRegions) {
            rotatedMaxY=std::max(rotatedMaxY,region.bounds.max.y);
        }
        check(rotatedMaxY>maxY,"a positive angle lifts the wake off the axis");
        // A single bounding box for a 15-degree, six-span wake would be far taller
        // than the staircase; this is why slices exist.
        check(rotatedMaxY<6.0*rotated.bodySpan*std::sin(15.0*std::acos(-1.0)/180.0)+
                          0.6*rotated.bodySpan+1e-9,
              "the staircase does not exceed the true rotated extent");
    }

    // Rejections must be explicit rather than silently producing a coarse mesh.
    {
        SizeFieldPolicy2D bad;
        bad.farFieldSpans=-1.0;
        check(!resolveSizeField2D(bad,circle).valid(),"negative farFieldSpans is rejected");
        SizeFieldPolicy2D coarse;
        coarse.wallCellsPerSpan=1e-9;
        check(!resolveSizeField2D(coarse,circle).valid(),"unrefinable wall request is rejected");
        SizeFieldPolicy2D deepFar;
        deepFar.wallDistance=WallDistanceSizing2D{3U,27U};
        check(!resolveSizeField2D(deepFar,circle).valid(),"farLevel above the wall is rejected");
    }

    // The safe-level ceiling is a measured construction limit, so crossing it has to
    // be an explicit request rather than a silent consequence of a fine wall.
    {
        SizeFieldPolicy2D deep;
        deep.farFieldSpans=10.0;
        deep.wallCellsPerSpan=4096.0;
        const auto refused=resolveSizeField2D(deep,circle);
        check(!refused.valid(),"a wall past maxSafeWallLevel is refused");
        check(!refused.issues.empty()&&
              refused.issues.front().find("maxSafeWallLevel")!=std::string::npos,
              "the refusal names the ceiling it hit");
        deep.allowUnsafeWallLevel=true;
        const auto allowed=resolveSizeField2D(deep,circle);
        check(allowed.valid(),"the opt-in permits a deeper wall");
        check(allowed.wallLevel>11U,"the opt-in actually resolves deeper");

        SizeFieldPolicy2D raised;
        raised.farFieldSpans=10.0;
        raised.wallCellsPerSpan=4096.0;
        raised.maxSafeWallLevel=28U;
        check(resolveSizeField2D(raised,circle).valid(),"raising the ceiling also works");
    }

    // Whatever the field compiles to has to be something refine() accepts, and the
    // resulting mesh has to actually be finer near the wall than in the far field.
    {
        SizeFieldPolicy2D policy;
        policy.farFieldSpans=4.0;
        policy.wallCellsPerSpan=64.0;
        policy.wallDistance=WallDistanceSizing2D{3U,1U};
        policy.curvature=CurvatureSizing2D{8.0};
        policy.wake=WakeSizing2D{0.0,3.0,0.5,2U,4U};
        const auto resolved=resolveSizeField2D(policy,circle);
        check(resolved.valid(),"combined field resolves");
        Quadtree2D tree(resolved.domain,resolved.maxLevel,circle);
        tree.refine(circle,resolved.refinement);
        const auto balance=tree.enforceTwoToOneBalance(circle);
        check(balance.violationsAfter==0U,"compiled field still balances 2:1");
        check(tree.deterministicOrderingValid(),"compiled field keeps deterministic ids");
        near(tree.totalLeafArea(),resolved.domain.width()*resolved.domain.height(),
             1e-6*resolved.domain.width()*resolved.domain.height(),
             "refined leaves still tile the domain");
        std::size_t wallLeaves=0,farLeaves=0;
        for (const auto& leaf:tree.leaves()) {
            check(leaf.level<=resolved.maxLevel,"no leaf exceeds the derived maxLevel");
            if (leaf.classification==CellClass::Intersected) {
                ++wallLeaves;
                check(leaf.level>=resolved.wallLevel,"wall leaves reach the wall level");
            }
            if (distanceAABBToBoundary(leaf.bounds,circle)>2.0*resolved.bodySpan) {
                ++farLeaves;
                check(leaf.level<resolved.wallLevel,"far leaves stay coarser than the wall");
            }
        }
        check(wallLeaves>0U,"the wall is resolved");
        check(farLeaves>0U,"the far field exists and is coarse");
    }

    // Leaving wallCellsPerSpan unset must land exactly on the trusted ceiling, so
    // the defaults can never refuse themselves.
    {
        SizeFieldPolicy2D policy;
        const auto resolved=resolveSizeField2D(policy,circle);
        check(resolved.valid(),"the all-defaults policy resolves");
        check(resolved.wallLevel==policy.maxSafeWallLevel,
              "an unset wall request lands on maxSafeWallLevel");
        check(resolved.maxLevel==policy.maxSafeWallLevel,"and nothing pushes past it");
        SizeFieldPolicy2D shallow;
        shallow.maxSafeWallLevel=7U;
        const auto followed=resolveSizeField2D(shallow,circle);
        check(followed.valid()&&followed.wallLevel==7U,
              "an unset wall request follows the ceiling it is given");
    }

    if (failures!=0) {
        std::cerr<<failures<<" size-field check(s) failed\n";
        return EXIT_FAILURE;
    }
    std::cout<<"size field tests passed\n";
    return EXIT_SUCCESS;
}
