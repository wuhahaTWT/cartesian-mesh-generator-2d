#include "cartmesh2d/fv/Sst2003m2D.hpp"
#include "cartmesh2d/fv/detail/FlowFaceOperators2D.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <sstream>

namespace cartmesh2d::fv {
namespace {
void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
double checked(double value) {
    require(std::isfinite(value), "SST-2003m numerical range exceeded");
    return value;
}
void nonnegative(double value) {
    require(std::isfinite(value)&&value>=0, "SST-2003m requires finite nonnegative k and strain");
}
void positive(double value) {
    require(std::isfinite(value)&&value>0, "SST-2003m requires finite positive omega, nu and wall distance");
}
void converged(const ScalarTransportResult2D& r,const char* name) {
    if (r.converged) return;
    std::ostringstream error;
    error<<"SST-2003m "<<name<<" inner transport did not converge";
    if (!r.history.empty()) error<<"; norm="<<r.history.back().residualNorm
        <<", diagonalScaled="<<r.history.back().maxDiagonalScaledImbalance;
    throw std::runtime_error(error.str());
}
void boundaries(const FvMesh2D& mesh, const std::vector<ScalarBoundary2D>& bc, bool isK) {
    for (std::size_t id=0;id<bc.size();++id) {
        if (mesh.faces[id].neighbour) continue;
        const auto valid=[&](double v) { return std::isfinite(v)&&(isK?v>=0:v>0); };
        if (bc[id].kind==ScalarBoundaryKind2D::Value)
            require(valid(bc[id].value), "SST-2003m inadmissible boundary value");
        if (bc[id].inflowValue)
            require(valid(*bc[id].inflowValue), "SST-2003m inadmissible inflow value");
    }
}
}
Sst2003mCoefficients2D evaluateSst2003m2D(const Sst2003mPoint2D& p) {
    nonnegative(p.k); nonnegative(p.strainMagnitude);
    positive(p.omega); positive(p.nu); positive(p.wallDistance);
    for (double v:{p.gradientK.x,p.gradientK.y,p.gradientOmega.x,p.gradientOmega.y}) checked(v);
    constexpr double betaStar=.09,a1=.31,sigmaOmega2=.856;
    const double d2=checked(p.wallDistance*p.wallDistance);
    require(d2>0,"SST-2003m wall distance square underflow");
    const double crossBase=checked(2*sigmaOmega2*checked(dot(p.gradientK,p.gradientOmega))/p.omega);
    const double cd=std::max(crossBase,1e-10);
    const double turbulent=checked(std::sqrt(p.k)/p.omega/betaStar/p.wallDistance);
    const double viscous=checked(500*p.nu/d2/p.omega);
    const double crossBound=checked(4*sigmaOmega2*p.k/cd/d2);
    const double arg1=std::min(std::max(turbulent,viscous),crossBound);
    const double arg2=std::max(checked(2*turbulent),viscous);
    // tanh(>=20) rounds to 1 in binary64. Bound its argument before powers,
    // not any physical field or source, to avoid overflow at small wall distance.
    const double a=std::min(arg1,3.),b=std::min(arg2,5.);
    Sst2003mCoefficients2D r;
    r.f1=std::tanh(a*a*a*a); r.f2=std::tanh(b*b);
    const auto blend=[&](double inner,double outer) { return r.f1*inner+(1-r.f1)*outer; };
    r.sigmaK=blend(.85,1.); r.sigmaOmega=blend(.5,sigmaOmega2);
    r.beta=blend(.075,.0828); r.gamma=blend(5./9.,.44);
    const double denominator=std::max(checked(a1*p.omega),checked(p.strainMagnitude*r.f2));
    require(denominator>0,"SST-2003m viscosity denominator underflow");
    r.turbulentViscosity=checked(a1*p.k/denominator);
    r.diffusivityK=checked(p.nu+r.sigmaK*r.turbulentViscosity);
    r.diffusivityOmega=checked(p.nu+r.sigmaOmega*r.turbulentViscosity);
    const double strain2=checked(p.strainMagnitude*p.strainMagnitude);
    // P/nu_t has this continuous k->0 limit; never divide 0/0 at k=0.
    // The SAME factor-10 limiter applies to BOTH production equations.
    const double limitedStrain2=std::min(strain2,checked(10*betaStar*p.omega*denominator/a1));
    r.productionK=checked(r.turbulentViscosity*limitedStrain2);
    r.productionOmega=checked(r.gamma*limitedStrain2);
    r.crossDiffusion=checked((1-r.f1)*crossBase);
    r.sourceK=r.productionK;
    r.sourceOmega=checked(r.productionOmega+std::max(r.crossDiffusion,0.));
    r.lossRateK=checked(betaStar*p.omega);
    r.lossRateOmega=checked(r.beta*p.omega+std::max(-r.crossDiffusion,0.)/p.omega);
    return r;
}

FrozenSst2003mResult2D solveFrozenSst2003mTransport2D(const FvMesh2D& mesh,
    const FrozenSst2003mProblem2D& p, const ScalarTransportControls2D& controls,
    const std::vector<double>& previousK, const std::vector<double>& previousOmega, double dt) {
    validateFvMesh2D(mesh);
    const auto n=mesh.cells.size(),nf=mesh.faces.size();
    require(n>0&&p.k.size()==n&&p.omega.size()==n&&p.wallDistance.size()==n&&
        p.strainMagnitude.size()==n&&p.gradientK.size()==n&&p.gradientOmega.size()==n&&
        p.volumeFlux.size()==nf&&p.boundaryK.size()==nf&&p.boundaryOmega.size()==nf,
        "SST-2003m field dimensions mismatch");
    require(previousK.empty()==previousOmega.empty(),"SST-2003m both previous states required");
    for (double k:previousK) nonnegative(k);
    for (double w:previousOmega) positive(w);
    boundaries(mesh,p.boundaryK,true); boundaries(mesh,p.boundaryOmega,false);
    FrozenSst2003mResult2D result;
    result.coefficients.reserve(n);
    ScalarTransportProblem2D kp,wp;
    kp.diffusivity=wp.diffusivity=p.nu;
    kp.volumeFlux=wp.volumeFlux=p.volumeFlux;
    kp.boundaryData=p.boundaryK; wp.boundaryData=p.boundaryOmega;
    for (std::size_t i=0;i<n;++i) {
        const auto c=evaluateSst2003m2D({p.k[i],p.omega[i],p.nu,p.wallDistance[i],
            p.strainMagnitude[i],p.gradientK[i],p.gradientOmega[i]});
        result.coefficients.push_back(c);
        kp.sourceDensity.push_back(c.sourceK); wp.sourceDensity.push_back(c.sourceOmega);
        kp.sinkRate.push_back(c.lossRateK); wp.sinkRate.push_back(c.lossRateOmega);
    }
    for (const auto& f:mesh.faces) {
        const auto& co=result.coefficients[f.owner];
        const auto& cn=f.neighbour?result.coefficients[*f.neighbour]:co;
        const double weight=f.neighbour?f.neighbourWeight:0;
        kp.faceDiffusivity.push_back((1-weight)*co.diffusivityK+weight*cn.diffusivityK);
        wp.faceDiffusivity.push_back((1-weight)*co.diffusivityOmega+weight*cn.diffusivityOmega);
    }
    result.k=solveScalarTransport2D(mesh,kp,controls,previousK,dt);
    converged(result.k,"k");
    require(result.k.minValue>=0,"SST-2003m negative k after transport; no clipping applied");
    result.omega=solveScalarTransport2D(mesh,wp,controls,previousOmega,dt);
    converged(result.omega,"omega");
    require(result.omega.minValue>0,"SST-2003m nonpositive omega after transport; no clipping applied");
    return result;
}

Sst2003mGradients2D reconstructSst2003mGradients2D(const FvMesh2D& mesh,
    const FrozenSst2003mProblem2D& p,const std::vector<Vector2D>& velocity,
    const std::vector<SstVelocityBoundary2D>& bc) {
    validateFvMesh2D(mesh);
    const auto n=mesh.cells.size(),nf=mesh.faces.size();
    require(n>0&&p.k.size()==n&&p.omega.size()==n&&velocity.size()==n&&bc.size()==nf&&
        p.boundaryK.size()==nf&&p.boundaryOmega.size()==nf,"SST gradient field dimensions mismatch");
    boundaries(mesh,p.boundaryK,true);boundaries(mesh,p.boundaryOmega,false);
    std::vector<double> u(n),v(n),ub(nf),vb(nf),kb(nf),wb(nf);
    std::vector<bool> uf(nf),vf(nf),kf(nf),wf(nf);
    for(std::size_t i=0;i<n;++i) {
        nonnegative(p.k[i]);positive(p.omega[i]);
        u[i]=checked(velocity[i].x);v[i]=checked(velocity[i].y);
    }
    for(std::size_t id=0;id<nf;++id) {
        if(mesh.faces[id].neighbour)continue;
        const auto scalar=[&](const ScalarBoundary2D& b,double& value,std::vector<bool>& fixed) {
            require(b.kind==ScalarBoundaryKind2D::Value||b.kind==ScalarBoundaryKind2D::DiffusiveFlux,
                "SST gradient invalid boundary kind");
            if(b.kind==ScalarBoundaryKind2D::DiffusiveFlux)
                require(b.value==0,"SST gradient requires zero diffusive flux or a value boundary");
            value=checked(b.value);fixed[id]=b.kind==ScalarBoundaryKind2D::Value;
        };
        scalar(p.boundaryK[id],kb[id],kf);scalar(p.boundaryOmega[id],wb[id],wf);
        ub[id]=checked(bc[id].value.x);vb[id]=checked(bc[id].value.y);
        uf[id]=bc[id].fixedX;vf[id]=bc[id].fixedY;
    }
    Sst2003mGradients2D result;
    result.k=detail::flowGradient(mesh,p.k,kb,kf);result.omega=detail::flowGradient(mesh,p.omega,wb,wf);
    result.u=detail::flowGradient(mesh,u,ub,uf);result.v=detail::flowGradient(mesh,v,vb,vf);
    result.strainMagnitude.resize(n);
    for(std::size_t i=0;i<n;++i) {
        // sqrt(2 ux^2 + 2 vy^2 + (uy+vx)^2). Off-diagonal tensor
        // components occur twice; rigid rotation must give zero production.
        result.strainMagnitude[i]=checked(std::hypot(std::sqrt(2.)*result.u[i].x,
            std::sqrt(2.)*result.v[i].y,checked(result.u[i].y+result.v[i].x)));
    }
    return result;
}
}
