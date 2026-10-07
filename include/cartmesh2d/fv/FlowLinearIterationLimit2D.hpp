#pragma once
#include <stdexcept>

namespace cartmesh2d::fv {
// Finite iterates exhausted the linear budget without meeting the original
// residual gates. Invalid matrices, numerical breakdown and nonfinite values
// retain their distinct fatal exceptions.
class FlowLinearIterationLimit2D : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};
}
