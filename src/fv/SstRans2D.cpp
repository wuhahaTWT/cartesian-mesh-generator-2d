#include "cartmesh2d/fv/SstRans2D.hpp"
#include "cartmesh2d/fv/WallDistance2D.hpp"
#include "cartmesh2d/fv/detail/FlowMaterial2D.hpp"
#include <algorithm>
#include <cmath>
#include <chrono>
#include <stdexcept>
#include <utility>

namespace cartmesh2d::fv {
namespace {
void require(bool ok,const char* message) {if(!ok)throw std::runtime_error(message);}
}
SstRansResult2D solveSstRans2D(const FvMesh2D& mesh,const SstRansControls2D& c,
    const std::vector<double>& initialK,const std::vector<double>& initialOmega,
    const std::function<void(const FlowIteration2D&)>& progress) {
    validateFvMesh2D(mesh);
    require(c.turbulenceUpdatesPerIteration>0,"SST RANS requires positive turbulence updates per iteration");
    require(c.flow.viscousStress==ViscousStress2D::Symmetric&&c.flow.faceViscosity.empty(),
        "SST RANS requires symmetric stress and owns its face viscosity");
    require(c.flow.outletBackflow==OutletBackflow2D::Reject,"SST RANS backflow boundary not implemented");
    require(std::isfinite(c.inletK)&&c.inletK>=0&&std::isfinite(c.inletOmega)&&c.inletOmega>0,
        "SST RANS requires nonnegative inlet k and positive omega");
    const auto n=mesh.cells.size(),nf=mesh.faces.size();
    require(initialK.empty()==initialOmega.empty(),"SST RANS requires both initial turbulence fields");
    FrozenSst2003mProblem2D p;p.nu=c.flow.nu;
    p.k=initialK.empty()?std::vector<double>(n,c.inletK):initialK;
    p.omega=initialOmega.empty()?std::vector<double>(n,c.inletOmega):initialOmega;
    require(p.k.size()==n&&p.omega.size()==n,"SST RANS initial field size mismatch");
    for(double v:p.k)require(std::isfinite(v)&&v>=0,"SST RANS invalid initial k");
    for(double v:p.omega)require(std::isfinite(v)&&v>0,"SST RANS invalid initial omega");
    p.boundaryK.resize(nf);p.boundaryOmega.resize(nf);
    SstRansResult2D result;result.resolvedWalls.resize(nf);result.velocityBoundary.resize(nf);
    using Clock=std::chrono::steady_clock;
    const auto elapsed=[](Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();};
    auto update=[&](const FlowResult2D& flow,const std::vector<detail::MaterialBoundary2D>& bc) {
        const auto updateStart=c.flow.profile?Clock::now():Clock::time_point{};
        p.volumeFlux=flow.flux;
        std::vector<Vector2D> velocity(n);
        for(std::size_t i=0;i<n;++i)velocity[i]={flow.u[i],flow.v[i]};
        for(std::size_t id=0;id<nf;++id) {
            const auto& b=bc[id];result.resolvedWalls[id]=b.wall;
            result.velocityBoundary[id]={b.velocity,b.fixedX,b.fixedY};
            if(mesh.faces[id].neighbour)continue;
            require(!b.outlet||flow.flux[id]>=0,"SST RANS outlet backflow requires explicit turbulence inlet data");
            p.boundaryK[id]=b.inlet?ScalarBoundary2D{ScalarBoundaryKind2D::Value,c.inletK,{}}:
                ScalarBoundary2D{ScalarBoundaryKind2D::DiffusiveFlux,0,{}};
            p.boundaryOmega[id]=b.inlet?ScalarBoundary2D{ScalarBoundaryKind2D::Value,c.inletOmega,{}}:
                ScalarBoundary2D{ScalarBoundaryKind2D::DiffusiveFlux,0,{}};
        }
        if(p.wallDistance.empty()) {
            const auto start=c.flow.profile?Clock::now():Clock::time_point{};
            p.wallDistance=computeWallDistance2D(mesh,result.resolvedWalls).distance;
            if(c.flow.profile)result.performance.wallDistanceSeconds+=elapsed(start);
        }
        setSst2003mResolvedWalls2D(mesh,p,result.resolvedWalls);
        auto transport=c.turbulence;
        transport.transport.profile=c.flow.profile;
        transport.maxIterations=std::min(transport.maxIterations,c.turbulenceUpdatesPerIteration);
        const auto transportStart=c.flow.profile?Clock::now():Clock::time_point{};
        auto next=solveSst2003mTransport2D(mesh,p,velocity,result.velocityBoundary,transport);
        if(c.flow.profile) {
            result.performance.transportSeconds+=elapsed(transportStart);
            result.performance.scalarSolves.add(next.scalarSolves);
            result.performance.scalarEvaluations.add(next.scalarEvaluations);
            result.performance.kSolves.add(next.kSolves);result.performance.omegaSolves.add(next.omegaSolves);
        }
        p.k=next.fields.k.values;p.omega=next.fields.omega.values;
        // Current strain/gradients/coefficient fields are evaluated at this same
        // velocity and returned k/omega. Use inlet closure from actual inlet
        // turbulence values; never extrapolate owner k/omega onto that boundary.
        const auto gradientStart=c.flow.profile?Clock::now():Clock::time_point{};
        const auto g=reconstructSst2003mGradients2D(mesh,p,velocity,result.velocityBoundary);
        if(c.flow.profile)result.performance.gradientSeconds+=elapsed(gradientStart);
        result.faceViscosity.resize(nf);
        for(std::size_t id=0;id<nf;++id) {
            const auto& f=mesh.faces[id];double nt=next.fields.coefficients[f.owner].turbulentViscosity;
            if(result.resolvedWalls[id])nt=0;
            else if(f.neighbour)nt=(1-f.neighbourWeight)*nt+f.neighbourWeight*next.fields.coefficients[*f.neighbour].turbulentViscosity;
            else if(bc[id].inlet)nt=evaluateSst2003m2D({c.inletK,c.inletOmega,p.nu,p.wallDistance[f.owner],
                g.strainMagnitude[f.owner],g.k[f.owner],g.omega[f.owner]}).turbulentViscosity;
            result.faceViscosity[id]=p.nu+nt;
        }
        const auto& h=next.history.back();
        result.history.push_back({result.history.size()+1,h.iteration,h.kResidualNorm,h.omegaResidualNorm,h.kCellResidual,h.omegaCellResidual});
        result.turbulence=std::move(next);
        if(c.flow.profile) {++result.performance.updates;result.performance.updateSeconds+=elapsed(updateStart);}
        return detail::MaterialState2D{result.faceViscosity,result.turbulence.converged};
    };
    result.flow=detail::solveMaterialFlow2D(mesh,c.flow,update,progress);
    result.converged=result.flow.converged&&result.turbulence.converged;
    result.wallDistance=p.wallDistance;result.boundaryK=p.boundaryK;result.boundaryOmega=p.boundaryOmega;
    return result;
}
}
