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
    if (!r.history.empty() && r.history.back().matrixAudited)
        error<<", rounded-field matrix norm="<<r.history.back().matrixResidualNorm
             <<", matrix diagonalScaled="<<r.history.back().matrixMaxDiagonalScaledImbalance;
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

namespace {
FrozenSst2003mResult2D sstTransport(const FvMesh2D& mesh,
    const FrozenSst2003mProblem2D& p, const ScalarTransportControls2D& controls,
    const std::vector<double>& previousK, const std::vector<double>& previousOmega, double dt,bool evaluateOnly,
    bool correctionOnly=false) {
    validateFvMesh2D(mesh);
    const auto n=mesh.cells.size(),nf=mesh.faces.size();
    require(n>0&&p.k.size()==n&&p.omega.size()==n&&p.wallDistance.size()==n&&
        p.strainMagnitude.size()==n&&p.gradientK.size()==n&&p.gradientOmega.size()==n&&
        p.volumeFlux.size()==nf&&p.boundaryK.size()==nf&&p.boundaryOmega.size()==nf,
        "SST-2003m field dimensions mismatch");
    require(p.resolvedWalls.empty()||p.resolvedWalls.size()==nf,"SST wall mask dimensions mismatch");
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
    for (std::size_t id=0;id<nf;++id) {
        const auto& f=mesh.faces[id];
        if (!p.resolvedWalls.empty()&&p.resolvedWalls[id]) {
            require(!f.neighbour&&p.volumeFlux[id]==0&&p.boundaryK[id].kind==ScalarBoundaryKind2D::Value&&
                p.boundaryK[id].value==0&&p.boundaryOmega[id].kind==ScalarBoundaryKind2D::Value,
                "SST resolved wall requires boundary k=0, positive omega and zero volume flux");
            kp.faceDiffusivity.push_back(p.nu);wp.faceDiffusivity.push_back(p.nu);continue;
        }
        const auto& co=result.coefficients[f.owner];
        const auto& cn=f.neighbour?result.coefficients[*f.neighbour]:co;
        const double weight=f.neighbour?f.neighbourWeight:0;
        kp.faceDiffusivity.push_back((1-weight)*co.diffusivityK+weight*cn.diffusivityK);
        wp.faceDiffusivity.push_back((1-weight)*co.diffusivityOmega+weight*cn.diffusivityOmega);
    }
    if(evaluateOnly) {
        result.k=evaluateScalarTransport2D(mesh,kp,p.k,controls,previousK,dt);
        result.omega=evaluateScalarTransport2D(mesh,wp,p.omega,controls,previousOmega,dt);
        return result;
    }
    if(correctionOnly)require(previousK.empty() && previousOmega.empty() && dt==0,
                         "SST steady initial guess with time history");
    result.k=correctionOnly?solveSteadyScalarTransportFromInitial2D(mesh,kp,p.k,controls)
                      :solveScalarTransport2D(mesh,kp,controls,previousK,dt);
    if(!correctionOnly)converged(result.k,"k");
    require(result.k.minValue>=0,"SST-2003m negative k after transport; no clipping applied");
    result.omega=correctionOnly?solveSteadyScalarTransportFromInitial2D(mesh,wp,p.omega,controls)
                          :solveScalarTransport2D(mesh,wp,controls,previousOmega,dt);
    if(!correctionOnly)converged(result.omega,"omega");
    require(result.omega.minValue>0,"SST-2003m nonpositive omega after transport; no clipping applied");
    return result;
}
}
FrozenSst2003mResult2D solveFrozenSst2003mTransport2D(const FvMesh2D& mesh,
    const FrozenSst2003mProblem2D& p,const ScalarTransportControls2D& controls,
    const std::vector<double>& previousK,const std::vector<double>& previousOmega,double dt) {
    return sstTransport(mesh,p,controls,previousK,previousOmega,dt,false);
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

void setSst2003mResolvedWalls2D(const FvMesh2D& mesh,FrozenSst2003mProblem2D& p,
    const std::vector<bool>& walls) {
    validateFvMesh2D(mesh);positive(p.nu);
    const auto nf=mesh.faces.size();
    require(walls.size()==nf&&p.boundaryK.size()==nf&&p.boundaryOmega.size()==nf&&p.volumeFlux.size()==nf,
        "SST wall boundary dimensions mismatch");
    std::vector<double> values(nf);
    for(std::size_t id=0;id<nf;++id)if(walls[id]) {
        const auto& f=mesh.faces[id];
        require(!f.neighbour&&p.volumeFlux[id]==0,"SST resolved wall must be impermeable boundary");
        const double distance=checked(dot(f.centre-mesh.cells[f.owner].centre,f.areaVector)/
            std::hypot(f.areaVector.x,f.areaVector.y));positive(distance);
        values[id]=checked(60*p.nu/.075/distance/distance);positive(values[id]);
    }
    // Validate everything before changing caller state.
    p.resolvedWalls=walls;
    for(std::size_t id=0;id<nf;++id)if(walls[id]) {
        p.boundaryK[id]={ScalarBoundaryKind2D::Value,0,{}};
        p.boundaryOmega[id]={ScalarBoundaryKind2D::Value,values[id],{}};
    }
}

SstTransportResult2D solveSst2003mTransport2D(const FvMesh2D& mesh,
    const FrozenSst2003mProblem2D& initial,const std::vector<Vector2D>& velocity,
    const std::vector<SstVelocityBoundary2D>& velocityBC,const SstTransportControls2D& controls,
    const std::vector<double>& previousK,const std::vector<double>& previousOmega,double dt) {
    require(controls.maxIterations>0&&std::isfinite(controls.relaxation)&&controls.relaxation>0&&controls.relaxation<=1,
        "SST invalid nonlinear iteration controls");
    const bool correctionOnly=controls.scalarCorrectionsPerUpdate>0;
    require(!correctionOnly || (previousK.empty()&&previousOmega.empty()&&dt==0),
            "SST scalar correction updates are currently steady only");
    auto p=initial;
    const auto reconstruct=[&] {
        auto g=reconstructSst2003mGradients2D(mesh,p,velocity,velocityBC);
        p.gradientK=std::move(g.k);p.gradientOmega=std::move(g.omega);p.strainMagnitude=std::move(g.strainMagnitude);
    };
    reconstruct();
    auto inner=controls.transport;
    inner.relativeTolerance*=.1;inner.absoluteTolerance*=.1;inner.cellTolerance*=.1;
    if(correctionOnly)inner.maxCorrections=std::min(inner.maxCorrections,controls.scalarCorrectionsPerUpdate);
    SstTransportResult2D result;
    const auto record=[&](const FrozenSst2003mResult2D& fields,bool evaluation) {
        auto& total=evaluation?result.scalarEvaluations:result.scalarSolves;
        total.add(fields.k.performance);total.add(fields.omega.performance);
    };
    // Includes complete scalar/input validation and allows an already converged
    // initial state to return without inventing a nonlinear update.
    result.fields=sstTransport(mesh,p,controls.transport,previousK,previousOmega,dt,true);
    record(result.fields,true);
    for(std::size_t it=0;it<=controls.maxIterations;++it) {
        const auto& k=result.fields.k.history.back();const auto& w=result.fields.omega.history.back();
        result.history.push_back({it,k.residualNorm,w.residualNorm,k.maxDiagonalScaledImbalance,w.maxDiagonalScaledImbalance});
        if(result.fields.k.converged&&result.fields.omega.converged) {result.converged=true;break;}
        if(it==controls.maxIterations)break;
        // A bounded correction is only an iterate. Reconstruct and evaluate all
        // ORIGINAL nonlinear equations below; never propagate its convergence
        // flag as the nonlinear/RANS acceptance. Linear failures still throw.
        const auto candidate=sstTransport(mesh,p,inner,previousK,previousOmega,dt,false,correctionOnly);
        record(candidate,false);
        for(std::size_t i=0;i<p.k.size();++i) {
            p.k[i]=checked(p.k[i]+controls.relaxation*(candidate.k.values[i]-p.k[i]));
            p.omega[i]=checked(p.omega[i]+controls.relaxation*(candidate.omega.values[i]-p.omega[i]));
            nonnegative(p.k[i]);positive(p.omega[i]);
        }
        reconstruct();
        result.fields=sstTransport(mesh,p,controls.transport,previousK,previousOmega,dt,true);
        record(result.fields,true);
    }
    return result;
}
}
