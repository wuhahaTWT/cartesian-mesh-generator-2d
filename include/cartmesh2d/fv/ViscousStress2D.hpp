#pragma once
#include "cartmesh2d/fv/FvMesh2D.hpp"
#include <array>

namespace cartmesh2d::fv {
// Velocity is a prescribed face value. Euler's static-mesh wall adapter permits
// only tangential motion. Open boundaries use zero viscous traction explicitly.
enum class ViscousBoundaryKind2D { Velocity, Slip, ZeroTraction, Periodic };
struct ViscousBoundary2D {
    std::size_t face=0;
    ViscousBoundaryKind2D kind=ViscousBoundaryKind2D::ZeroTraction;
    Vector2D velocity{};
    std::optional<std::size_t> partner;
};
struct ViscousStressResult2D {
    // Outward conservative momentum/energy flux: (-tau.S)x,y, -u_face.tau.S.
    std::vector<std::array<double,3>> faceFlux,cellResidual;
    std::vector<double> rate;
};
class ViscousStressOperator2D {
public:
    ViscousStressOperator2D(const FvMesh2D&,const std::vector<ViscousBoundary2D>&,double dynamicViscosity);
    [[nodiscard]] ViscousStressResult2D evaluate(const std::vector<Vector2D>& velocity,const std::vector<double>& density) const;
private:
    struct Sample { std::optional<std::size_t> cell; Vector2D value{},weight{},normal{}; ViscousBoundaryKind2D kind{}; };
    struct Face {
        std::size_t owner=0;
        std::optional<std::size_t> neighbour,partner;
        ViscousBoundaryKind2D kind=ViscousBoundaryKind2D::ZeroTraction;
        Vector2D value{},normal{},area{},d{},ownerOffset{},neighbourOffset{};
        double weight=0,normalDistance=0;
    };
    std::vector<Face> faces_;
    std::vector<std::vector<Sample>> gradients_;
    std::vector<double> areas_,rowNorm_;
    double viscosity_=0;
};
} // namespace cartmesh2d::fv
