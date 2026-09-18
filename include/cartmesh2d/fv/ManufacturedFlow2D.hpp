#pragma once

#include "cartmesh2d/geometry/Geometry2D.hpp"
#include <cmath>
#include <numbers>

namespace cartmesh2d::fv {

// Verification problem on [0,1]^2, not a physical flow preset.
// psi=(speed/pi)*sin(pi*x)^2*sin(pi*y)^2; U=(psi_y,-psi_x).
struct ManufacturedFlowSample2D {
    Vector2D velocity;
    double pressure;
    Vector2D acceleration; // div(U U) + grad(p) - nu laplacian(U)
};

inline ManufacturedFlowSample2D manufacturedFlow2D(Point2D point, double speed, double nu) {
    constexpr double pi = std::numbers::pi;
    const double sx=std::sin(pi*point.x), sy=std::sin(pi*point.y);
    const double cx=std::cos(pi*point.x), cy=std::cos(pi*point.y);
    const double s2x=std::sin(2*pi*point.x), s2y=std::sin(2*pi*point.y);
    const double c2x=std::cos(2*pi*point.x), c2y=std::cos(2*pi*point.y);
    const double u=speed*sx*sx*s2y, v=-speed*s2x*sy*sy;
    const double ux=speed*pi*s2x*s2y, uy=2*speed*pi*sx*sx*c2y;
    const double vx=-2*speed*pi*c2x*sy*sy, vy=-ux;
    const double lapU=2*speed*pi*pi*s2y*(2*c2x-1);
    const double lapV=-2*speed*pi*pi*s2x*(2*c2y-1);
    const double px=-speed*speed*pi*sx*cy, py=-speed*speed*pi*cx*sy;
    return {{u,v},speed*speed*cx*cy,
            {u*ux+v*uy+px-nu*lapU,u*vx+v*vy+py-nu*lapV}};
}

} // namespace cartmesh2d::fv
