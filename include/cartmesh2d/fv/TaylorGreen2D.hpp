#pragma once
#include "cartmesh2d/geometry/Geometry2D.hpp"
#include <cmath>
#include <numbers>

namespace cartmesh2d::fv {
// Exact unforced decaying vortex on [0,1]^2, impermeable free-slip walls.
// psi=A/pi*sin(pi*x)*sin(pi*y), A=speed*exp(-2*nu*pi^2*time).
struct TaylorGreenSample2D { Vector2D velocity; double pressure; double streamfunction; };
inline TaylorGreenSample2D taylorGreen2D(Point2D x, double time, double speed, double nu) {
    constexpr double pi=std::numbers::pi;
    const double a=speed*std::exp(-2*nu*pi*pi*time);
    const double sx=std::sin(pi*x.x), sy=std::sin(pi*x.y);
    return {{a*sx*std::cos(pi*x.y),-a*std::cos(pi*x.x)*sy},
            .25*a*a*(std::cos(2*pi*x.x)+std::cos(2*pi*x.y)),a/pi*sx*sy};
}
}
