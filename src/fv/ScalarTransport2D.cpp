#include "cartmesh2d/fv/ScalarTransport2D.hpp"
#include "cartmesh2d/fv/detail/FlowFaceOperators2D.hpp"
#include "cartmesh2d/fv/detail/FlowLinearSystem2D.hpp"
#include <algorithm>
#include <cmath>
#include <chrono>
#include <limits>
#include <stdexcept>

namespace cartmesh2d::fv {
struct ScalarTransportWorkspace2D::Impl {
    struct Storage {
        std::vector<std::pair<std::size_t,std::size_t>> connections;
        detail::SparsePattern2D pattern;
        detail::SparseSystem2D matrix;
        std::unique_ptr<detail::LinearWorkspace2D> krylov;
        Storage(std::size_t n,const std::vector<std::pair<std::size_t,std::size_t>>& edges)
            :connections(edges),pattern(n,edges),matrix(pattern) {}
    };
    std::unique_ptr<Storage> storage;
    bool active=false;
};
ScalarTransportWorkspace2D::ScalarTransportWorkspace2D()=default;
ScalarTransportWorkspace2D::~ScalarTransportWorkspace2D()=default;
ScalarTransportWorkspace2D::ScalarTransportWorkspace2D(ScalarTransportWorkspace2D&&) noexcept=default;
ScalarTransportWorkspace2D& ScalarTransportWorkspace2D::operator=(ScalarTransportWorkspace2D&&) noexcept=default;
struct ScalarTransportWorkspaceAccess2D {
    struct Lease {
        ScalarTransportWorkspace2D::Impl& state;
        explicit Lease(ScalarTransportWorkspace2D& w):state(get(w)) {
            if(state.active)throw std::runtime_error("Scalar transport workspace is already in use");
            state.active=true;
        }
        ~Lease(){state.active=false;}
        Lease(const Lease&)=delete;
        Lease& operator=(const Lease&)=delete;
        static ScalarTransportWorkspace2D::Impl& get(ScalarTransportWorkspace2D& w) {
            if(!w.impl_)w.impl_=std::make_unique<ScalarTransportWorkspace2D::Impl>();
            return *w.impl_;
        }
        bool prepare(std::size_t n,const std::vector<std::pair<std::size_t,std::size_t>>& edges) {
            const bool reused=state.storage && state.storage->matrix.diag.size()==n && state.storage->connections==edges;
            if(!reused)state.storage=std::make_unique<ScalarTransportWorkspace2D::Impl::Storage>(n,edges);
            // Invalidate all factors and reset numeric coefficients, even after
            // a failed preceding call. No stale material/source/boundary data.
            state.storage->matrix.reset();
            return reused;
        }
    };
};
namespace {
using Values = std::vector<double>;
void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
double finite(double value) {
    require(std::isfinite(value), "Scalar transport numerical range exceeded");
    return value;
}
double faceDiffusivity(const ScalarTransportProblem2D& p, std::size_t id) {
    return p.faceDiffusivity.empty()?p.diffusivity:p.faceDiffusivity[id];
}
// Least squares uses actual prescribed normal derivatives at flux boundaries,
// not an artificial zero gradient or a guessed face value.
std::vector<Vector2D> gradients(const FvMesh2D& mesh, const Values& value,
    const std::vector<ScalarBoundary2D>& bc, const ScalarTransportProblem2D& p) {
    std::vector<Vector2D> result(value.size());
    for (std::size_t i=0;i<value.size();++i) {
        double xx=0,xy=0,yy=0,bx=0,by=0;
        for (auto id:mesh.cells[i].faces) {
            const auto& f=mesh.faces[id];
            const auto j=f.owner==i?f.neighbour:std::optional<std::size_t>(f.owner);
            Vector2D d; double delta=0;
            if (j || bc[id].kind==ScalarBoundaryKind2D::Value) {
                d=(j?mesh.cells[*j].centre:f.centre)-mesh.cells[i].centre;
                const double length=std::hypot(d.x,d.y);
                require(length>0,"Scalar transport degenerate gradient stencil");
                delta=((j?value[*j]:bc[id].value)-value[i])/length;
                d=d*(1/length);
            } else {
                const double length=std::hypot(f.areaVector.x,f.areaVector.y);
                d=f.areaVector*(1/length);
                delta=-bc[id].value/faceDiffusivity(p,id);
            }
            xx+=d.x*d.x; xy+=d.x*d.y; yy+=d.y*d.y;
            bx+=d.x*delta; by+=d.y*delta;
        }
        const double det=xx*yy-xy*xy;
        require(det>64*std::numeric_limits<double>::epsilon()*(xx+yy)*(xx+yy),
            "Scalar transport rank deficient gradient stencil");
        result[i]={finite((yy*bx-xy*by)/det),finite((xx*by-xy*bx)/det)};
    }
    return result;
}
}

namespace {
ScalarTransportResult2D scalarTransport(const FvMesh2D& mesh,
    const ScalarTransportProblem2D& p, const ScalarTransportControls2D& c,
    const Values& previous, double timeStep, const Values* supplied,
    const Values* initial=nullptr, ScalarTransportWorkspace2D* reusable=nullptr) {
    using Clock=std::chrono::steady_clock;
    const auto start=c.profile?Clock::now():Clock::time_point{};
    const auto seconds=[](Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();};
    ScalarTransportWorkspace2D localWorkspace;
    ScalarTransportWorkspaceAccess2D::Lease lease(reusable?*reusable:localWorkspace);
    validateFvMesh2D(mesh);
    const auto n=mesh.cells.size(),nf=mesh.faces.size();
    require(c.preconditioner==ScalarPreconditioner2D::Jacobi || c.preconditioner==ScalarPreconditioner2D::ILU0,
            "Scalar transport invalid preconditioner");
    const auto method=c.preconditioner==ScalarPreconditioner2D::ILU0
        ?detail::LinearSolveMethod2D::ILU0:detail::LinearSolveMethod2D::Jacobi;
    if(initial) {
        require(!supplied && previous.empty() && timeStep==0,
                "Scalar steady initial guess cannot be mixed with evaluation or time history");
        require(initial->size()==n,"Scalar steady initial field dimensions mismatch");
        for(double value:*initial)finite(value);
    }
    if (supplied) {
        require(supplied->size()==n,"Scalar evaluation field dimensions mismatch");
        for (double v:*supplied) finite(v);
    }
    require(n>0 && p.volumeFlux.size()==nf,
        "Scalar transport invalid field dimensions or missing callback");
    require(bool(p.source)!=(!p.sourceDensity.empty()) &&
        (p.sourceDensity.empty() || p.sourceDensity.size()==n),
        "Scalar transport requires exactly one valid source representation");
    require(bool(p.boundary)!=(!p.boundaryData.empty()) &&
        (p.boundaryData.empty() || p.boundaryData.size()==nf),
        "Scalar transport requires exactly one valid boundary representation");
    require(std::isfinite(p.diffusivity)&&p.diffusivity>0,
        "Scalar transport diffusivity must be finite positive");
    require(p.faceDiffusivity.empty() || p.faceDiffusivity.size()==nf,
        "Scalar transport face diffusivity dimensions mismatch");
    for (double d:p.faceDiffusivity) require(std::isfinite(d)&&d>0,
        "Scalar transport face diffusivity must be finite positive");
    require(p.sinkRate.empty() || p.sinkRate.size()==n,
        "Scalar transport sink rate dimensions mismatch");
    for (double rate:p.sinkRate) require(std::isfinite(rate)&&rate>=0,
        "Scalar transport sink rate must be finite nonnegative");
    require(c.maxCorrections>0 && std::isfinite(c.relaxation)&&c.relaxation>0&&c.relaxation<=1,
        "Scalar transport invalid correction controls");
    require(c.convection==ConvectionScheme2D::Upwind || c.convection==ConvectionScheme2D::LimitedLinearUpwind,
        "Scalar transport invalid convection scheme");
    for (double t:{c.relativeTolerance,c.absoluteTolerance,c.cellTolerance,
                   c.carrierRelativeTolerance,c.carrierAbsoluteTolerance})
        require(std::isfinite(t)&&t>0,"Scalar transport tolerances must be finite positive");
    const bool transient=!previous.empty();
    require(std::isfinite(timeStep) && (transient?(previous.size()==n&&timeStep>0):timeStep==0),
        "Scalar transport invalid previous state or time step");
    for (double v:previous) finite(v);
    for (double q:p.volumeFlux) finite(q);
    std::vector<ScalarBoundary2D> bc(nf);
    Values carrier(n),carrierScale(n),boundaryValues(nf);
    std::vector<bool> fixed(nf,false),anchored(n,transient);
    std::vector<std::size_t> queue;
    std::vector<std::pair<std::size_t,std::size_t>> connections;
    ScalarTransportResult2D r;
    for (std::size_t id=0;id<nf;++id) {
        const auto& f=mesh.faces[id]; const double q=p.volumeFlux[id];
        carrier[f.owner]+=q; carrierScale[f.owner]+=std::abs(q);
        if (f.neighbour) {
            carrier[*f.neighbour]-=q; carrierScale[*f.neighbour]+=std::abs(q);
            connections.emplace_back(f.owner,*f.neighbour);
        } else {
            bc[id]=p.boundary?p.boundary(id,f):p.boundaryData[id]; finite(bc[id].value);
            require(bc[id].kind==ScalarBoundaryKind2D::Value || bc[id].kind==ScalarBoundaryKind2D::DiffusiveFlux,
                "Scalar transport invalid boundary kind");
            if (bc[id].inflowValue) finite(*bc[id].inflowValue);
            fixed[id]=bc[id].kind==ScalarBoundaryKind2D::Value;
            require(q>=0 || fixed[id] || bc[id].inflowValue.has_value(),
                "Scalar transport inflow requires a prescribed scalar value");
            boundaryValues[id]=fixed[id]?bc[id].value:bc[id].inflowValue.value_or(0.);
            if (fixed[id] || q<0) anchored[f.owner]=true;
        }
    }
    for (std::size_t i=0;i<n;++i) {
        finite(carrier[i]); finite(carrierScale[i]);
        if (!p.sinkRate.empty() && p.sinkRate[i]>0) anchored[i]=true;
        r.maxCarrierImbalance=std::max(r.maxCarrierImbalance,std::abs(carrier[i]));
        require(std::abs(carrier[i])<=finite(c.carrierAbsoluteTolerance+c.carrierRelativeTolerance*carrierScale[i]),
            "Scalar transport carrier flux violates cell continuity");
        if (transient) r.maxCourant=std::max(r.maxCourant,finite(.5*timeStep*carrierScale[i]/mesh.cells[i].area));
        if (anchored[i]) queue.push_back(i);
    }
    for (std::size_t head=0;head<queue.size();++head) {
        const auto i=queue[head];
        for (auto id:mesh.cells[i].faces) {
            const auto& f=mesh.faces[id]; if (!f.neighbour) continue;
            const auto j=f.owner==i?*f.neighbour:f.owner;
            if (!anchored[j]) { anchored[j]=true; queue.push_back(j); }
        }
    }
    require(queue.size()==n,"Scalar transport unanchored steady connected component");
    // Evaluation needs the original diagonal/RHS and face balances, but no
    // Krylov workspace. Materialize CSR only for solves or the additional tight
    // rounded-field audit below. Never reuse an older coefficient matrix.
    detail::SparseSystem2D* matrix=nullptr;
    detail::LinearWorkspace2D* workspace=nullptr;
    std::size_t priorBuilds=0,priorReuses=0;
    const auto createSparse=[&] {
        const bool reused=lease.prepare(n,connections);
        matrix=&lease.state.storage->matrix;
        if(!supplied && !lease.state.storage->krylov)
            lease.state.storage->krylov=std::make_unique<detail::LinearWorkspace2D>(n);
        workspace=lease.state.storage->krylov.get();
        priorBuilds=matrix->ilu0Builds();priorReuses=matrix->ilu0Reuses();
        if(c.profile) {
            if(reused)++r.performance.patternReuses;
            else ++r.performance.patternBuilds;
        }
    };
    if(!supplied)createSparse();
    Values evaluationRhs(supplied?n:0),evaluationDiagonal(supplied?n:0);
    auto& rhs=supplied?evaluationRhs:matrix->rhs;
    auto& diag=supplied?evaluationDiagonal:matrix->diag;
    r.values=supplied?*supplied:(initial?*initial:(transient?previous:Values(n,0.)));
    r.sourceIntegrals.resize(n); r.temporalIntegrals.resize(n);
    r.sinkIntegrals.resize(n);
    r.advectiveFlux.resize(nf); r.diffusiveFlux.resize(nf);
    for (std::size_t i=0;i<n;++i) {
        rhs[i]=r.sourceIntegrals[i]=finite((p.source?p.source(mesh.cells[i].centre):p.sourceDensity[i])*mesh.cells[i].area);
        r.sourceIntegral=finite(r.sourceIntegral+rhs[i]);
        if (transient) {
            const double mass=finite(mesh.cells[i].area/timeStep);
            diag[i]=mass; rhs[i]+=finite(mass*previous[i]);
        }
        if (!p.sinkRate.empty()) diag[i]+=finite(p.sinkRate[i]*mesh.cells[i].area);
    }
    for (std::size_t id=0;id<nf;++id) {
        const auto& f=mesh.faces[id]; const auto i=f.owner;
        const double q=p.volumeFlux[id],d=finite(faceDiffusivity(p,id)*f.transmissibility);
        if (f.neighbour) {
            const auto j=*f.neighbour;
            diag[i]+=d+std::max(q,0.); if(matrix)matrix->add(i,j,-d+std::min(q,0.));
            diag[j]+=d+std::max(-q,0.); if(matrix)matrix->add(j,i,-d-std::max(q,0.));
        } else {
            if (fixed[id]) { diag[i]+=d; rhs[i]+=d*bc[id].value; }
            else rhs[i]-=finite(bc[id].value*std::hypot(f.areaVector.x,f.areaVector.y));
            if (q<0) rhs[i]-=q*boundaryValues[id];
            else diag[i]+=q;
        }
    }
    const Values base=rhs,diagonal=diag;
    for (double d:diagonal) require(std::isfinite(d)&&d>0,"Scalar transport invalid matrix diagonal");
    const double scale=detail::linearNorm(base);
    const double stop=finite(c.absoluteTolerance+c.relativeTolerance*scale);
    Values extra(nf),residual(n);
    if(c.profile) {r.performance.calls=1;r.performance.setupSeconds=seconds(start);}
    const auto faceFluxes=[&]() {
        const auto phaseStart=c.profile?Clock::now():Clock::time_point{};
        const auto g=gradients(mesh,r.values,bc,p);
        auto limitFixed=fixed;
        for (std::size_t id=0;id<nf;++id)
            if (!mesh.faces[id].neighbour && p.volumeFlux[id]<0) limitFixed[id]=true;
        const auto limiter=c.convection==ConvectionScheme2D::LimitedLinearUpwind
            ?detail::faceReconstructionLimiter(mesh,r.values,g,boundaryValues,limitFixed):Values{};
        for (std::size_t id=0;id<nf;++id) {
            const auto& f=mesh.faces[id]; const auto i=f.owner; const double q=p.volumeFlux[id];
            const double diffusivity=faceDiffusivity(p,id);
            auto gf=g[i];
            if (f.neighbour) {
                const auto gn=g[*f.neighbour]; const double w=f.neighbourWeight;
                gf={gf.x*(1-w)+gn.x*w,gf.y*(1-w)+gn.y*w};
            }
            const double diffCorrection=(!f.neighbour&&!fixed[id])?0:-diffusivity*dot(gf,f.correction);
            const double upwind=(!f.neighbour&&q<0)?boundaryValues[id]:r.values[(f.neighbour&&q<0)?*f.neighbour:i];
            const double advected=(!f.neighbour&&q<0)?boundaryValues[id]
                :detail::upwindFaceValue(mesh,id,q,r.values,g,limiter);
            r.advectiveFlux[id]=finite(q*advected);
            r.diffusiveFlux[id]=finite(!f.neighbour&&!fixed[id]
                ?bc[id].value*std::hypot(f.areaVector.x,f.areaVector.y)
                :diffusivity*f.transmissibility*(r.values[i]-(f.neighbour?r.values[*f.neighbour]:bc[id].value))+diffCorrection);
            extra[id]=finite(diffCorrection+q*(advected-upwind));
        }
        if(c.profile)r.performance.faceFluxSeconds+=seconds(phaseStart);
    };
    for (std::size_t it=1;it<=(supplied?1:c.maxCorrections);++it) {
        std::size_t linear=0;
        if (!supplied) {
            auto& a=*matrix;
            faceFluxes(); a.rhs=base;
            for (std::size_t id=0;id<nf;++id) {
                const auto& f=mesh.faces[id]; a.rhs[f.owner]-=extra[id];
                if (f.neighbour) a.rhs[*f.neighbour]+=extra[id];
            }
            // Relax the field after solving the full elliptic operator. Diagonal
            // equation relaxation would turn even a linear orthogonal diffusion
            // problem into a slow, mesh-dependent stationary outer iteration.
            // User-requested tight outer tolerances must also constrain the inner
            // solve; its legacy norm floor otherwise stalls relaxed scalar solves.
            // An explicit high+low candidate can satisfy the same linear gates
            // when its low bits cannot fit in the returned double field. Round
            // only after relaxation; faceFluxes and original balances below
            // recheck the actual double values. Candidate success is not outer
            // transport convergence.
            const auto linearStart=c.profile?Clock::now():Clock::time_point{};
            const auto candidate=a.solveCandidate(r.values,*workspace,c.cellTolerance*.1,stop*.5,method);
            if(c.profile) {r.performance.linearSeconds+=seconds(linearStart);r.performance.linearIterations+=candidate.iterations;}
            linear=candidate.iterations;
            for (std::size_t i=0;i<n;++i)
                r.values[i]=candidate.relaxedDouble(i,r.values[i],c.relaxation);
        }
        faceFluxes(); r.boundaryFlux=0; r.temporalIntegral=0; r.sinkIntegral=0;
        for (std::size_t i=0;i<n;++i) {
            r.temporalIntegrals[i]=transient?finite(mesh.cells[i].area*(r.values[i]-previous[i])/timeStep):0;
            r.temporalIntegral=finite(r.temporalIntegral+r.temporalIntegrals[i]);
            residual[i]=r.temporalIntegrals[i]-r.sourceIntegrals[i];
            if (!p.sinkRate.empty()) {
                r.sinkIntegrals[i]=finite(p.sinkRate[i]*mesh.cells[i].area*r.values[i]);
                r.sinkIntegral=finite(r.sinkIntegral+r.sinkIntegrals[i]);
                residual[i]+=r.sinkIntegrals[i];
            }
        }
        for (std::size_t id=0;id<nf;++id) {
            const auto& f=mesh.faces[id]; const double flux=finite(r.advectiveFlux[id]+r.diffusiveFlux[id]);
            residual[f.owner]+=flux;
            if (f.neighbour) residual[*f.neighbour]-=flux; else r.boundaryFlux+=flux;
        }
        double maxCell=0,maxScaled=0;
        for (std::size_t i=0;i<n;++i) {
            maxCell=std::max(maxCell,std::abs(finite(residual[i])));
            maxScaled=std::max(maxScaled,std::abs(residual[i])/diagonal[i]);
        }
        const double norm=detail::linearNorm(residual);
        r.history.push_back({supplied?0:it,linear,norm,finite(norm/std::max(scale,c.absoluteTolerance)),maxCell,maxScaled});
        r.boundaryFlux=finite(r.boundaryFlux);
        r.globalBalance=finite(r.temporalIntegral+r.boundaryFlux-r.sourceIntegral);
        if (!p.sinkRate.empty()) r.globalBalance=finite(r.globalBalance+r.sinkIntegral);
        if (norm<=stop && maxScaled<=c.cellTolerance) {
            double magnitude=0;
            for(double value:r.values)magnitude=std::max(magnitude,std::abs(value));
            if(c.cellTolerance<16*std::numeric_limits<double>::epsilon()*magnitude) {
                // Reassemble deferred sources at the rounded, returned field,
                // not at the old iterate used by the preceding linear solve.
                // Accurate b-Ax prevents a rounded product from creating a
                // false zero in an exceptionally tight outer acceptance check.
                if(!matrix) {
                    const auto matrixStart=c.profile?Clock::now():Clock::time_point{};
                    createSparse();matrix->diag=diagonal;
                    for(std::size_t id=0;id<nf;++id) {
                        const auto& f=mesh.faces[id];if(!f.neighbour)continue;
                        const double q=p.volumeFlux[id],d=finite(faceDiffusivity(p,id)*f.transmissibility);
                        matrix->add(f.owner,*f.neighbour,-d+std::min(q,0.));
                        matrix->add(*f.neighbour,f.owner,-d-std::max(q,0.));
                    }
                    if(c.profile)r.performance.setupSeconds+=seconds(matrixStart);
                }
                auto& a=*matrix;a.rhs=base;
                for(std::size_t id=0;id<nf;++id) {
                    const auto& f=mesh.faces[id];a.rhs[f.owner]-=extra[id];
                    if(f.neighbour)a.rhs[*f.neighbour]+=extra[id];
                }
                auto& entry=r.history.back();entry.matrixAudited=true;
                for(std::size_t i=0;i<n;++i) {
                    residual[i]=a.compensatedResidualRow(i,r.values);
                    entry.matrixMaxDiagonalScaledImbalance=std::max(
                        entry.matrixMaxDiagonalScaledImbalance,std::abs(residual[i])/diagonal[i]);
                }
                entry.matrixResidualNorm=detail::linearNorm(residual);
                if(entry.matrixResidualNorm>stop || entry.matrixMaxDiagonalScaledImbalance>c.cellTolerance)
                    continue;
            }
            r.converged=true;break;
        }
    }
    r.minValue=*std::min_element(r.values.begin(),r.values.end());
    r.maxValue=*std::max_element(r.values.begin(),r.values.end());
    if(c.profile) {
        r.performance.totalSeconds=seconds(start);
        if(matrix) {r.performance.ilu0Builds=matrix->ilu0Builds()-priorBuilds;r.performance.ilu0Reuses=matrix->ilu0Reuses()-priorReuses;}
    }
    return r;
}
}
ScalarTransportResult2D solveScalarTransport2D(const FvMesh2D& mesh,
    const ScalarTransportProblem2D& p,const ScalarTransportControls2D& c,
    const Values& previous,double timeStep,ScalarTransportWorkspace2D* workspace) {
    return scalarTransport(mesh,p,c,previous,timeStep,nullptr,nullptr,workspace);
}
ScalarTransportResult2D evaluateScalarTransport2D(const FvMesh2D& mesh,
    const ScalarTransportProblem2D& p,const Values& values,const ScalarTransportControls2D& c,
    const Values& previous,double timeStep,ScalarTransportWorkspace2D* workspace) {
    return scalarTransport(mesh,p,c,previous,timeStep,&values,nullptr,workspace);
}
ScalarTransportResult2D solveSteadyScalarTransportFromInitial2D(const FvMesh2D& mesh,
    const ScalarTransportProblem2D& p,const Values& initial,const ScalarTransportControls2D& c,ScalarTransportWorkspace2D* workspace) {
    return scalarTransport(mesh,p,c,{},0,nullptr,&initial,workspace);
}
}
