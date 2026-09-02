// R2/W1: the grid-corner weld budget is a policy value, and the default must
// stay roundoff-sized.
//
// The Q2-A shared construction deliberately absorbed arithmetic roundoff only.
// R2 needs a *geometric* weld to remove refinement-induced corner spurs, but
// turning that on globally would silently change every existing construction, so
// the budget became opt-in. These tests pin both halves of that contract: the
// default is unchanged, and a budget that could destroy a face the Q1 contract
// would accept is refused outright rather than merely discouraged.

#include "cartmesh2d/geometry/IntersectionRegistry2D.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace cartmesh2d;

namespace {

int failures = 0;

void check(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

void defaultWeldBudgetStaysRoundoffSized() {
    const IntersectionRegistryPolicy2D policy;
    check(policy.gridCornerWeldFractionOfLocalH == policy.snapFractionOfLocalH,
          "the default weld budget must equal the roundoff proximity budget");
    check(policy.gridCornerWeldFractionOfLocalH ==
              64.0 * std::numeric_limits<double>::epsilon(),
          "the default weld budget must stay at 64 * DBL_EPSILON");
}

void geometricWeldBudgetIsAcceptedBelowTheQ1ShortFaceLimit() {
    IntersectionRegistryPolicy2D policy;
    policy.gridCornerWeldFractionOfLocalH = 1.0e-4;
    try {
        IntersectionRegistry2D registry(policy);
        (void)registry;
    } catch (const std::exception& error) {
        std::cerr << "unexpected rejection: " << error.what() << '\n';
        ++failures;
    }
}

void weldBudgetAtOrAboveTheQ1ShortFaceLimitIsRefused() {
    // Q1's hard limit is face_length / local_h >= 0.01. A weld budget that large
    // could remove a face the contract considered legitimate.
    for (const double fraction : {0.01, 0.1, 0.5, 1.0}) {
        IntersectionRegistryPolicy2D policy;
        policy.gridCornerWeldFractionOfLocalH = fraction;
        bool refused = false;
        try {
            IntersectionRegistry2D registry(policy);
            (void)registry;
        } catch (const std::invalid_argument&) {
            refused = true;
        }
        if (!refused) {
            std::cerr << "FAIL: weld fraction " << fraction
                      << " must be refused as at or above the Q1 short-face limit\n";
            ++failures;
        }
    }
}

void nonFiniteOrNonPositiveWeldBudgetsAreRefused() {
    for (const double fraction : {0.0, -1.0e-6,
                                  std::numeric_limits<double>::quiet_NaN(),
                                  std::numeric_limits<double>::infinity()}) {
        IntersectionRegistryPolicy2D policy;
        policy.gridCornerWeldFractionOfLocalH = fraction;
        bool refused = false;
        try {
            IntersectionRegistry2D registry(policy);
            (void)registry;
        } catch (const std::invalid_argument&) {
            refused = true;
        }
        if (!refused) {
            std::cerr << "FAIL: weld fraction " << fraction << " must be refused\n";
            ++failures;
        }
    }
}

void proximityBudgetValidationIsUnchanged() {
    IntersectionRegistryPolicy2D policy;
    policy.snapFractionOfLocalH = 0.5;
    bool refused = false;
    try {
        IntersectionRegistry2D registry(policy);
        (void)registry;
    } catch (const std::invalid_argument&) {
        refused = true;
    }
    check(refused, "a proximity budget of 0.5 must still be refused");
}

} // namespace

int main() {
    defaultWeldBudgetStaysRoundoffSized();
    geometricWeldBudgetIsAcceptedBelowTheQ1ShortFaceLimit();
    weldBudgetAtOrAboveTheQ1ShortFaceLimitIsRefused();
    nonFiniteOrNonPositiveWeldBudgetsAreRefused();
    proximityBudgetValidationIsUnchanged();
    if (failures != 0) {
        std::cerr << failures << " weld-budget checks failed\n";
        return EXIT_FAILURE;
    }
    std::cout << "grid-corner weld budget tests passed\n";
    return EXIT_SUCCESS;
}
