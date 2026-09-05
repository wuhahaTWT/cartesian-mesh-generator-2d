#include "cartmesh2d/quality/QualityContract2D.hpp"
#include "cartmesh2d/stabilization/SmallCell2D.hpp"
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

// A gate nobody can satisfy is not a strict gate, it is a broken one.  Every hard
// limit whose reachable bound follows from the stabilization threshold and the 2:1
// level rule is checked here, so a future tightening cannot quietly re-create the
// state where all five acceptance cases read FAIL for reasons no mesh can fix.
const std::vector<std::pair<const char*,OrdinaryCellQualityLimits2D>> types(
    const QualityContract2D& contract) {
    return {{"cartesian",contract.cartesian},{"remainder_cut",contract.remainderCut},
            {"transition",contract.transition},{"termination",contract.termination}};
}
} // namespace

int main(){
    const QualityContract2D contract;
    const SmallCellPolicy2D stabilization;
    const double mergeThreshold=stabilization.areaFractionThreshold;
    check(mergeThreshold>0.0&&mergeThreshold<1.0,"the merge threshold is a fraction");

    for (const auto& [name,limits]:types(contract)) {
        const std::string where(name);

        // background_volume_ratio = min(a_own,a_nei)/max(...) with a = area/h^2.
        // Agglomeration guarantees a >= mergeThreshold for a surviving cut cell and
        // a <= 1 for a cell that fills its background box, so the ratio is bounded
        // below by mergeThreshold.  A hard limit above that is unreachable.
        check(limits.backgroundVolumeRatio.hard<=mergeThreshold,
              where+": background_volume_ratio hard limit "+
              std::to_string(limits.backgroundVolumeRatio.hard)+
              " exceeds the reachable bound "+std::to_string(mergeThreshold));
        check(limits.backgroundVolumeRatio.hard<=limits.backgroundVolumeRatio.preferred,
              where+": background_volume_ratio hard limit is not looser than preferred");

        // The raw OpenFOAM form has to absorb the 2:1 level jump as well.  A cut cell
        // of area fraction `a` against a coarser full neighbour scores a/4, so with
        // agglomeration only guaranteeing a >= mergeThreshold the reachable bound is
        // mergeThreshold/4.  This is the one gated limit that is knowingly above its
        // bound; the assertion below records the size of the gap rather than pretending
        // the gate is sound, and it fails if the gap ever widens.
        const double rawBound=mergeThreshold*0.25;
        check(limits.volumeRatio.hard>rawBound,
              where+": raw volume_ratio hard limit "+std::to_string(limits.volumeRatio.hard)+
              " is now at or below its reachable bound "+std::to_string(rawBound)+
              " - it has become a sound gate and this expectation should be inverted");
        check(limits.volumeRatio.hard<=4.0*rawBound,
              where+": the raw volume_ratio gap widened past the documented factor of"
                    " four; the gate is drifting further from reachable, not closer");
    }

    // Same defect-free pair under the two forms: below 1.0 raw once levels differ,
    // exactly 1.0 gated at every level difference.  This is the property that makes
    // the gated form eligible to replace the raw gate.
    for (std::size_t levels=0;levels<=8;++levels) {
        const double coarse=std::ldexp(1.0,-static_cast<int>(levels));
        const double raw=(coarse*coarse)/1.0;
        const double gated=((coarse*coarse)/(coarse*coarse))/(1.0/1.0);
        check(std::abs(gated-1.0)<=1e-15,
              "two full cells score 1.0 gated at every level difference");
        if (levels>0) {
            check(raw<1.0,"two full cells score below 1.0 raw once levels differ");
        }
    }

    // face_weight is a ratio of centre-to-face distances, so it is already
    // dimensionless and unaffected by grading.  It is *not* bounded by the merge
    // threshold either: a wedge of area fraction mergeThreshold can still place its
    // centroid arbitrarily close to its long face.  That is why this test asserts
    // only the internal consistency of the pair, and why the limit stays empirical.
    for (const auto& [name,limits]:types(contract)) {
        const std::string where(name);
        check(limits.faceWeight.hard<=limits.faceWeight.preferred,
              where+": face_weight hard limit is not looser than preferred");
        check(limits.faceWeight.hard<0.5,
              where+": face_weight cannot exceed 0.5 by construction");
    }

    // Every remaining lower-bound limit must at least be loose enough that a perfect
    // uniform Cartesian mesh passes it, otherwise the contract fails on the one mesh
    // it must certainly accept.
    for (const auto& [name,limits]:types(contract)) {
        const std::string where(name);
        check(limits.faceOverLocalBackgroundH.hard<1.0,
              where+": a face equal to the background size must pass");
        check(limits.faceOverSqrtOwnerArea.hard<1.0,
              where+": a face equal to sqrt(owner area) must pass");
        check(limits.faceOverSqrtNeighbourArea.hard<1.0,
              where+": a face equal to sqrt(neighbour area) must pass");
        check(limits.minimumInteriorAngleDeg.hard<90.0,
              where+": a square cell's 90 degree corner must pass");
        check(limits.nonOrthogonalityDeg.hard>0.0,
              where+": an orthogonal face must pass");
        check(limits.skewness.hard>0.0,where+": an unskewed face must pass");
        check(limits.hydraulicAspect.hard>1.0,where+": a square cell must pass");
    }

    if (failures!=0) {
        std::cerr<<failures<<" reachability check(s) failed\n";
        return EXIT_FAILURE;
    }
    std::cout<<"quality reachability tests passed\n";
    return EXIT_SUCCESS;
}
