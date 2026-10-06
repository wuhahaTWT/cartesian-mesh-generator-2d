#pragma once
#include "cartmesh2d/fv/FvMesh2D.hpp"
#include <array>

namespace cartmesh2d::fv {
struct HybridViscousCellResult2D {
    // Per-unit-depth outward momentum and total-energy flux, using the same
    // trace velocity in traction work. Face order is faces() below.
    std::vector<std::array<double,3>> outwardFlux;
    double dissipation=0; // W/m, a sum of nonnegative strain/stabilization terms.
};

// Research building block, NOT a complete viscous/Euler operator. Unknowns
// are the real cell velocity and one shared velocity trace per actual face.
// A global trace solve must enforce opposite internal tractions before these
// local fluxes may be used in a conservative time update.
// Planar Newtonian gas retains the molecular Stokes 2/3 coefficient.
class HybridViscousCell2D {
public:
    // The mesh must come from makeFvMesh2D. Also verifies the local geometric
    // identities required for affine consistency and positive normal distances.
    HybridViscousCell2D(const FvMesh2D&,std::size_t cell,double dynamicViscosity);
    [[nodiscard]] const std::vector<std::size_t>& faces() const {return faces_;}
    // Symmetric local matrix taking delta=(trace-owner) to outward traction.
    // Row-major interleaved x/y face components, dimension 2*faces().size().
    // Entries have units Pa s. This matrix has the affine rigid-rotation kernel;
    // it must not be treated as a positive-definite local inverse.
    [[nodiscard]] const std::vector<double>& tractionMatrix() const {return matrix_;}
    [[nodiscard]] HybridViscousCellResult2D evaluate(Vector2D ownerVelocity,
        const std::vector<Vector2D>& faceVelocities) const;
private:
    std::vector<std::size_t> faces_;
    std::vector<Vector2D> gradientWeights_,areaVectors_;
    std::vector<double> residual_,stabilization_,matrix_;
    double area_=0,viscosity_=0;
};
} // namespace cartmesh2d::fv
