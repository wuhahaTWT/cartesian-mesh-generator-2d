#include "cartmesh2d/quadtree/Quadtree2D.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
using namespace cartmesh2d;
namespace { int failures=0; void check(bool c,const std::string&m){if(!c){++failures;std::cerr<<"FAIL: "<<m<<'\n';}} void near(double a,double b,double e,const std::string&m){check(std::abs(a-b)<=e,m);} }
int main(){
    const Domain2D domain{{{0.0,0.0},{4.0,4.0}}};
    BoundaryLoop boundary({{1.25,1.25},{2.75,1.25},{2.75,2.75},{1.25,2.75}});
    Quadtree2D tree(domain,5,boundary);
    near(tree.totalLeafArea(),16.0,1e-12,"root area covers domain");
    check(tree.leaves().size()==1,"starts with one root leaf");

    QuadtreeRefinementPolicy2D policy; policy.boundaryLevel=4; policy.distanceBands={{0.35,3}};
    tree.refine(boundary,policy);
    near(tree.totalLeafArea(),16.0,1e-10,"refinement preserves total domain area");
    check(tree.deterministicOrderingValid(),"deterministic IDs after refinement");
    bool sawBoundaryL4=false, sawFarCoarse=false, maxOk=true;
    for(const auto& leaf:tree.leaves()){
        if(leaf.classification==CellClass::Intersected && leaf.level==4) sawBoundaryL4=true;
        if(distanceAABBToBoundary(leaf.bounds,boundary)>0.6 && leaf.level<4) sawFarCoarse=true;
        if(leaf.level>5) maxOk=false;
    }
    check(sawBoundaryL4,"boundary leaves reach requested level 4");
    check(sawFarCoarse,"far field remains coarser than boundary");
    check(maxOk,"max level not exceeded");

    Quadtree2D globallyRefined(domain,5,boundary);
    QuadtreeRefinementPolicy2D globalPolicy;
    globalPolicy.minimumLevel=3;
    globalPolicy.boundaryLevel=4;
    globallyRefined.refine(boundary,globalPolicy);
    bool allAtLeastThree=true;
    bool boundaryAtFour=true;
    for (const auto& leaf:globallyRefined.leaves()) {
        allAtLeastThree=allAtLeastThree && leaf.level>=3;
        if (leaf.classification==CellClass::Intersected) {
            boundaryAtFour=boundaryAtFour && leaf.level==4;
        }
    }
    check(allAtLeastThree,"minimum level globally refines the PDE background grid");
    check(boundaryAtFour,"minimum level preserves finer boundary target");

    QuadtreeRefinementPolicy2D sizingPolicy;
    sizingPolicy.minimumLevel=1;
    sizingPolicy.boundaryLevel=3;
    sizingPolicy.distanceBands={{0.20,3}};
    sizingPolicy.boxRegions={
        {{{2.75,1.50},{4.00,2.50}},4},
        {{{3.25,1.75},{4.00,2.25}},5}
    };
    Quadtree2D sized(domain,5,boundary);
    sized.refine(boundary,sizingPolicy);
    bool wakeAtLeastFour=true;
    bool wakeCoreAtFive=true;
    bool untouchedCornerCoarse=false;
    for (const auto& leaf:sized.leaves()) {
        const auto overlaps=[](const AABB2D& a,const AABB2D& b) {
            return std::max(a.min.x,b.min.x)<std::min(a.max.x,b.max.x) &&
                   std::max(a.min.y,b.min.y)<std::min(a.max.y,b.max.y);
        };
        if (overlaps(leaf.bounds,sizingPolicy.boxRegions[0].bounds)) {
            wakeAtLeastFour=wakeAtLeastFour && leaf.level>=4;
        }
        if (overlaps(leaf.bounds,sizingPolicy.boxRegions[1].bounds)) {
            wakeCoreAtFive=wakeCoreAtFive && leaf.level==5;
        }
        if (leaf.bounds.max.x<=1.0 && leaf.bounds.max.y<=1.0 && leaf.level<4) {
            untouchedCornerCoarse=true;
        }
    }
    check(wakeAtLeastFour,"local box refines every overlapping wake leaf to its target");
    check(wakeCoreAtFive,"overlapping wake core takes the greatest target level");
    check(untouchedCornerCoarse,"local sizing field leaves an unrelated corner below wake level");

    auto reversedSizing=sizingPolicy;
    std::reverse(reversedSizing.boxRegions.begin(),reversedSizing.boxRegions.end());
    std::reverse(reversedSizing.distanceBands.begin(),reversedSizing.distanceBands.end());
    Quadtree2D sizedReversed(domain,5,boundary);
    sizedReversed.refine(boundary,reversedSizing);
    check(sized.leaves().size()==sizedReversed.leaves().size(),
          "sizing field order preserves leaf count");
    if (sized.leaves().size()==sizedReversed.leaves().size()) {
        for (std::size_t i=0;i<sized.leaves().size();++i) {
            check(sized.leaves()[i].key==sizedReversed.leaves()[i].key,
                  "sizing field order preserves deterministic keys");
        }
    }

    // Independent partition audit: every leaf must lie inside the domain, keys must be unique,
    // and no two distinct leaves may overlap with positive area. Together with the total-area
    // check above, this rules out both overlap and hidden coverage holes.
    std::vector<std::uint64_t> auditKeys;
    auditKeys.reserve(tree.leaves().size());
    for (const auto& leaf : tree.leaves()) {
        auditKeys.push_back(leaf.key);
        check(leaf.bounds.min.x >= domain.bounds.min.x - 1e-12 &&
              leaf.bounds.max.x <= domain.bounds.max.x + 1e-12 &&
              leaf.bounds.min.y >= domain.bounds.min.y - 1e-12 &&
              leaf.bounds.max.y <= domain.bounds.max.y + 1e-12,
              "every leaf remains inside domain bounds");
        check(leaf.area() > 0.0, "every leaf has positive area");
    }
    std::sort(auditKeys.begin(), auditKeys.end());
    check(std::adjacent_find(auditKeys.begin(), auditKeys.end()) == auditKeys.end(),
          "leaf keys are unique");
    std::size_t positiveAreaOverlaps = 0;
    for (std::size_t i = 0; i < tree.leaves().size(); ++i) {
        for (std::size_t j = i + 1; j < tree.leaves().size(); ++j) {
            const auto& a = tree.leaves()[i].bounds;
            const auto& b = tree.leaves()[j].bounds;
            const double overlapX = std::min(a.max.x, b.max.x) - std::max(a.min.x, b.min.x);
            const double overlapY = std::min(a.max.y, b.max.y) - std::max(a.min.y, b.min.y);
            if (overlapX > 1e-12 && overlapY > 1e-12) ++positiveAreaOverlaps;
        }
    }
    check(positiveAreaOverlaps == 0, "leaf positive-area overlap count is zero");

    Quadtree2D tree2(domain,5,boundary); tree2.refine(boundary,policy);
    check(tree.leaves().size()==tree2.leaves().size(),"repeat run leaf count stable");
    if(tree.leaves().size()==tree2.leaves().size()) for(std::size_t i=0;i<tree.leaves().size();++i){
        check(tree.leaves()[i].key==tree2.leaves()[i].key,"repeat run key stable");
        check(tree.leaves()[i].classification==tree2.leaves()[i].classification,"repeat run classification stable");
    }

    // Manufacture a strongly unbalanced corner path by repeatedly refining the lower-left descendant.
    Quadtree2D unbalanced(domain,5,boundary);
    for(int level=0;level<5;++level){
        const auto leaves=unbalanced.leaves();
        auto it=std::find_if(leaves.begin(),leaves.end(),[](const auto& leaf){
            return std::abs(leaf.bounds.max.x - 2.0) < 1e-12 && std::abs(leaf.bounds.min.y) < 1e-12;
        });
        if (level == 0) it = leaves.begin();
        check(it!=leaves.end(),"refinement candidate exists");
        if(it==leaves.end()) break;
        const auto key=it->key;
        if(!unbalanced.refineLeafByKey(key,boundary)) break;
    }
    check(unbalanced.countBalanceViolations()>0,"manual refinement creates 2:1 violations");
    const auto report=unbalanced.enforceTwoToOneBalance(boundary);
    check(report.violationsBefore>0,"balance report sees initial violations");
    check(report.violationsAfter==0,"2:1 violations removed");
    check(unbalanced.countBalanceViolations()==0,"final face neighbors satisfy 2:1");
    near(unbalanced.totalLeafArea(),16.0,1e-10,"balance preserves area");

    // Every reported neighbor must share a face segment, not merely a corner.
    for(const auto& pair:unbalanced.faceNeighbors()){
        const auto&a=unbalanced.leaves()[pair.first]; const auto&b=unbalanced.leaves()[pair.second];
        const double xOverlap=std::min(a.bounds.max.x,b.bounds.max.x)-std::max(a.bounds.min.x,b.bounds.min.x);
        const double yOverlap=std::min(a.bounds.max.y,b.bounds.max.y)-std::max(a.bounds.min.y,b.bounds.min.y);
        const bool verticalTouch=std::abs(a.bounds.max.x-b.bounds.min.x)<1e-12 || std::abs(b.bounds.max.x-a.bounds.min.x)<1e-12;
        const bool horizontalTouch=std::abs(a.bounds.max.y-b.bounds.min.y)<1e-12 || std::abs(b.bounds.max.y-a.bounds.min.y)<1e-12;
        check((verticalTouch && yOverlap>1e-12)||(horizontalTouch && xOverlap>1e-12),"neighbor shares positive-length face");
    }

    bool threw=false; try{ QuadtreeRefinementPolicy2D bad; bad.boundaryLevel=6; Quadtree2D t(domain,5,boundary); t.refine(boundary,bad);}catch(const std::invalid_argument&){threw=true;} check(threw,"target level above max rejected");
    threw=false; try{ QuadtreeRefinementPolicy2D bad; bad.minimumLevel=6; Quadtree2D t(domain,5,boundary); t.refine(boundary,bad);}catch(const std::invalid_argument&){threw=true;} check(threw,"minimum level above max rejected");
    threw=false; try{ QuadtreeRefinementPolicy2D bad; bad.boxRegions={{{{5.0,5.0},{6.0,6.0}},3}}; Quadtree2D t(domain,5,boundary); t.refine(boundary,bad);}catch(const std::invalid_argument&){threw=true;} check(threw,"box outside domain rejected");
    if(failures){std::cerr<<failures<<" test(s) failed\n";return EXIT_FAILURE;} std::cout<<"cartmesh2d 2D-2 quadtree tests passed\n"; return EXIT_SUCCESS;
}
