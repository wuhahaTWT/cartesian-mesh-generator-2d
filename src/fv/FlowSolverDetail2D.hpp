#pragma once
// Internal incompressible stages. Public case/solver APIs remain in fv/ headers.
#include "cartmesh2d/fv/Incompressible2D.hpp"
#include "cartmesh2d/fv/detail/FlowLinearSystem2D.hpp"
#include "cartmesh2d/fv/detail/FlowFaceOperators2D.hpp"
#include <cmath>
#include <stdexcept>
namespace cartmesh2d::fv::solver_detail {
using Vec=std::vector<double>;
using System=detail::SparseSystem2D;
inline void ensure(bool valid,const char* message){if(!valid)throw std::runtime_error(message);}
inline double finite(double x){ensure(std::isfinite(x),"Flow numerical range exceeded");return x;}
enum class Role { Wall, Inlet, Outlet, Slip, Lid, Farfield, Opening };

struct Boundary {
    std::vector<Role> role;
    Vec u;
    Vec v;
    Vec p; // populated only for explicit conditions; preset pressure remains zero
    double initialU = 0, initialV = 0, initialP = 0;
    std::vector<bool> fixedU;
    std::vector<bool> fixedV;
    std::vector<bool> fixedP;
    std::vector<bool> constantU;
    std::vector<bool> constantV;
    double xmin = 0;
    double xmax = 0;
    double ymin = 0;
    double ymax = 0;
    bool closed = false;
    bool pressureOpenings = false;
};

Boundary boundaries(const FvMesh2D&,const FlowControls2D&);
void updateOutletBoundary(Boundary&,const FvMesh2D&,const FlowControls2D&,const Vec&);
Vector2D interpolateGradient(const Face&,const std::vector<Vector2D>&);
double interpolate(const Face&,const Vec&);
double faceNu(const FlowControls2D&,std::size_t);
void validateViscosity(const FvMesh2D&,const FlowControls2D&);
void initializeCaseSources(const FvMesh2D&,const FlowControls2D&,FlowResult2D&);
void initializeCaseVelocity(const FvMesh2D&,const FlowControls2D&,const Boundary&,FlowResult2D&);
void momentum(System& a,
                const FvMesh2D& m,
                const FlowControls2D& c,
                const Boundary& b,
                const Vec& field,
                const Vec& flux,
                const std::vector<Vector2D>& gradField,
                const std::vector<Vector2D>& gp,
                const std::vector<Vector2D>& source,
                const std::vector<Vector2D>& stressCorrection,
                const std::vector<Vector2D>& faceVelocity,
                bool y,
                const Vec* previous,
                double timeStep);

void postprocessForces(const FvMesh2D&,const FlowControls2D&,const Boundary&,FlowResult2D&,
    const Vec&,const std::vector<Vector2D>&,const std::vector<Vector2D>&,
    const std::vector<Vector2D>&,const std::vector<Vector2D>&);
Vec rhieChowFlux(const FvMesh2D& m,const FlowControls2D& c,const Boundary& b,const FlowResult2D& r,
    const Vec& ra,const Vec& oldU,const Vec& oldV,const std::vector<Vector2D>& gu,
    const std::vector<Vector2D>& gv,const std::vector<Vector2D>& gup,const std::vector<Vector2D>& gvp,
    const std::vector<Vector2D>& gp,const std::vector<Vector2D>& forceGradient,const Vec& oldFluxDefect,
    bool previous,double timeStep,Vec& df);
Vec steadyRhieChowFlux(const FvMesh2D&,const FlowControls2D&,const Boundary&,const FlowResult2D&,
    const Vec& response,const std::vector<Vector2D>& gu,const std::vector<Vector2D>& gv,
    const std::vector<Vector2D>& gp,const std::vector<Vector2D>& pressureForce,
    const Vec& noFluxDefect,Vec& coefficients);
void prepareMonitors(const FvMesh2D& m,const FlowControls2D& c,const Boundary& b,FlowResult2D& r,
    std::vector<std::size_t>& monitorGroup,Vec& monitorLengths);
Vec physicalMonitors(const FvMesh2D& m,const FlowControls2D& c,const Boundary& b,const FlowResult2D& r,
    const std::vector<std::size_t>& monitorGroup,const Vec& monitorLengths,double pressureScale,
    const Vec& pressureFaces,const std::vector<Vector2D>& gu,const std::vector<Vector2D>& gv,
    const std::vector<Vector2D>& stressCorrection);
} // namespace cartmesh2d::fv::solver_detail
