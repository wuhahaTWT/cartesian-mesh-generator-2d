// Annulus research probe: calls the production finite-volume operators.
// Norm ratios/Rayleigh quotients describe observed modes, not a spectral bound
// or a convergence/physical-accuracy certificate. No solver control is changed.
#include "cartmesh2d/fv/FlowBoundaryIO2D.hpp"
#include "cartmesh2d/io/MeshIO2D.hpp"
#include "FlowSolverDetail2D.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace cartmesh2d;
using namespace cartmesh2d::fv;
using namespace cartmesh2d::fv::solver_detail;

int main(int argc, char** argv) {
    if (argc < 4 || argc > 8) {
        std::cerr << "Usage: coupling_probe MESH BOUNDARIES MODE [symmetric|laplacian] [steps] [weight] [pressure-passes]\n"
                     "MODE: outer, pressure, diffusion, stokes. Fixed annulus nu=.1, Uref=.5.\n"
                     "outer: power iteration of production one-step Stokes SIMPLE F(x)-F(0).\n"
                     "pressure: frozen diffusion-rAU deferred nonorthogonal pressure map.\n"
                     "weight in (0,1] probes (1-weight)*I + weight*T; default 1.\n"
                     "diffusion: native momentum with zero pressure/convection (diagnostic only).\n";
        return 2;
    }
    try {
#ifdef __APPLE__
        // Standalone probes must set this before the first Accelerate call.
        setenv("VECLIB_MAXIMUM_THREADS", "1", 0);
#endif
        const std::string mode=argv[3], stress=argc>4?argv[4]:"symmetric";
        ensure(mode=="outer" || mode=="pressure" || mode=="diffusion" || mode=="stokes", "Unknown probe mode");
        ensure(stress=="symmetric" || stress=="laplacian", "Unknown viscous stress form");
        const auto steps=argc>5?std::stoul(argv[5]):(mode=="outer" || mode=="pressure"?240UL:2500UL);
        const double weight=argc>6?std::stod(argv[6]):1.;
        const auto passes=argc>7?std::stoul(argv[7]):4UL;
        ensure(steps>0 && steps<=100000 && std::isfinite(weight) && weight>0 && weight<=1 && passes>=1 && passes<=4,
               "Invalid probe steps, weight or pressure passes");
        const auto input=readCm2dTopology(argv[1]);
        if (!input.valid()) throw std::runtime_error(input.error);
        const auto mesh=makeFvMesh2D(input.topology);
        FlowControls2D c;
        c.scenario="custom"; c.nu=.1; c.speed=.5; c.tolerance=1e-8;
        c.maxIterations=steps; c.momentumInertia=0; c.pressureCorrectionPasses=passes;
        c.pressurePreconditioner=detail::systemCholeskyAvailable2D()
            ? PressurePreconditioner2D::SystemCholesky:PressurePreconditioner2D::IncompleteCholesky0;
        c.convection=ConvectionScheme2D::FaceLimitedLinearUpwind;
        c.viscousStress=stress=="symmetric"?ViscousStress2D::Symmetric:ViscousStress2D::Laplacian;
        std::ifstream in(argv[2]); ensure(bool(in), "Cannot open boundary file");
        c.boundaryConditions=readFlowBoundaryConditions2D(in,mesh,c);
        const auto b=boundaries(mesh,c);
        ensure(b.closed, "Annulus probe requires a closed velocity-boundary domain");
        for (std::size_t id=0;id<mesh.faces.size();++id)
            ensure(mesh.faces[id].neighbour || (b.fixedU[id] && b.fixedV[id]), "Probe requires fixed wall velocities");
        std::cout << std::setprecision(17);
        const auto n=mesh.cells.size(), nf=mesh.faces.size();
        if (mode=="outer") {
            c.maxIterations=1;
            FlowInitialGuess2D zero;
            zero.u.resize(n); zero.v.resize(n); zero.p.resize(n); zero.flux.resize(nf);
            const auto base=solveIncompressibleFromGuess2D(mesh,c,zero);
            auto x=zero;
            for (std::size_t i=0;i<n;++i) {
                x.u[i]=std::sin(1.71*static_cast<double>(i));
                x.v[i]=std::cos(1.53*static_cast<double>(i));
                x.p[i]=std::sin(1.31*static_cast<double>(i));
            }
            x.p[0]=0;
            double area=0; for (const auto& cell:mesh.cells) area+=cell.area;
            const auto product=[&](const auto& a,const auto& other) {
                long double sum=0;
                for (std::size_t i=0;i<n;++i)
                    sum+=mesh.cells[i].area/area*(a.u[i]*other.u[i]+a.v[i]*other.v[i]+a.p[i]*other.p[i]);
                return static_cast<double>(sum);
            };
            std::cout << "iteration,normRatio,rayleigh,cellEigenResidual,fluxEigenResidual,worstCell,x,y\n";
            for (std::size_t it=1;it<=steps;++it) {
                const auto r=solveIncompressibleFromGuess2D(mesh,c,x);
                auto next=zero;
                for (std::size_t i=0;i<n;++i) {
                    next.u[i]=(1-weight)*x.u[i]+weight*(r.u[i]-base.u[i]);
                    next.v[i]=(1-weight)*x.v[i]+weight*(r.v[i]-base.v[i]);
                    next.p[i]=(1-weight)*x.p[i]+weight*(r.p[i]-base.p[i]);
                }
                for (std::size_t i=0;i<nf;++i) next.flux[i]=(1-weight)*x.flux[i]+weight*(r.flux[i]-base.flux[i]);
                const double size=std::sqrt(product(next,next));
                ensure(std::isfinite(size) && size>0, "Degenerate outer power iterate");
                const double ratio=size/std::sqrt(product(x,x)), rayleigh=product(x,next)/product(x,x);
                auto residual=next;
                for (std::size_t i=0;i<n;++i) {
                    residual.u[i]-=rayleigh*x.u[i]; residual.v[i]-=rayleigh*x.v[i]; residual.p[i]-=rayleigh*x.p[i];
                }
                double fluxError=0, fluxSize=0;
                for (std::size_t i=0;i<nf;++i) {
                    fluxError=std::max(fluxError,std::abs(next.flux[i]-rayleigh*x.flux[i]));
                    fluxSize=std::max(fluxSize,std::abs(next.flux[i]));
                }
                std::size_t worst=0;
                for (std::size_t i=1;i<n;++i)
                    if (std::hypot(next.u[i],next.v[i])>std::hypot(next.u[worst],next.v[worst])) worst=i;
                std::cout << it << ',' << ratio << ',' << rayleigh << ',' << std::sqrt(product(residual,residual))/size
                          << ',' << (fluxSize>0?fluxError/fluxSize:fluxError) << ',' << worst << ','
                          << mesh.cells[worst].centre.x << ',' << mesh.cells[worst].centre.y << '\n';
                for (auto* f:{&next.u,&next.v,&next.p,&next.flux}) for (auto& v:*f) v/=size;
                x=std::move(next);
            }
            return 0;
        }
        std::vector<std::pair<std::size_t,std::size_t>> edges;
        for (const auto& f:mesh.faces) if (f.neighbour) edges.emplace_back(f.owner,*f.neighbour);
        detail::SparsePattern2D pattern(n,edges);
        detail::LinearWorkspace2D work(n);
        if (mode=="pressure") {
            System ap(pattern);
            Vec ra(n),pc(n),next(n),zero(nf);
            for (const auto& f:mesh.faces) {
                const double d=c.nu*f.transmissibility;
                ra[f.owner]+=d; if (f.neighbour) ra[*f.neighbour]+=d;
            }
            for (std::size_t i=0;i<n;++i) {
                ra[i]=c.velocityRelaxation*mesh.cells[i].area/ra[i]; pc[i]=std::sin(1.71*static_cast<double>(i));
            }
            pc[0]=0;
            for (const auto& f:mesh.faces) if (f.neighbour) {
                const auto i=f.owner,j=*f.neighbour;
                const double df=interpolate(f,ra)*f.transmissibility;
                ap.diag[i]+=df; ap.diag[j]+=df; ap.add(i,j,-df); ap.add(j,i,-df);
            }
            ap.pin(0);
            const auto stencil=detail::buildFlowGradientStencil2D(mesh,b.fixedP,true);
            std::cout << "iteration,normRatio,rayleigh,eigenResidual,worstCell,x,y\n";
            for (std::size_t it=1;it<=steps;++it) {
                const auto g=stencil.apply(pc,zero);
                std::fill(ap.rhs.begin(),ap.rhs.end(),0);
                for (const auto& f:mesh.faces) if (f.neighbour) {
                    const double q=interpolate(f,ra)*dot(interpolateGradient(f,g),f.correction);
                    ap.rhs[f.owner]+=q; ap.rhs[*f.neighbour]-=q;
                }
                ap.rhs[0]=0; std::fill(next.begin(),next.end(),0);
                ap.solvePressure(next,work,detail::systemCholeskyAvailable2D()
                    ?detail::LinearPressureMethod2D::SystemCholesky:detail::LinearPressureMethod2D::IC0);
                for (std::size_t i=0;i<n;++i) next[i]=(1-weight)*pc[i]+weight*next[i];
                const double size=detail::linearNorm(next),before=detail::linearNorm(pc);
                ensure(std::isfinite(size) && size>0, "Degenerate pressure power iterate");
                const double rayleigh=detail::linearProduct(pc,next)/detail::linearProduct(pc,pc);
                Vec residual(n);
                for (std::size_t i=0;i<n;++i) residual[i]=next[i]-rayleigh*pc[i];
                const auto worst=static_cast<std::size_t>(std::max_element(next.begin(),next.end(),
                    [](double a,double other){return std::abs(a)<std::abs(other);})-next.begin());
                std::cout << it << ',' << size/before << ',' << rayleigh << ',' << detail::linearNorm(residual)/size
                          << ',' << worst << ',' << mesh.cells[worst].centre.x << ',' << mesh.cells[worst].centre.y << '\n';
                for (std::size_t i=0;i<n;++i) pc[i]=next[i]/size;
            }
            return 0;
        }
        if (mode=="diffusion") {
            System au(pattern),av(pattern);
            Vec u(n),v(n),flux(nf);
            std::vector<Vector2D> gp(n),source;
            std::cout << "iteration,velocityChange,worstCell,x,y,area\n";
            for (std::size_t it=1;it<=steps;++it) {
                const auto gu=detail::flowGradient(mesh,u,b.u,b.fixedU),gv=detail::flowGradient(mesh,v,b.v,b.fixedV);
                const auto stressCorrection=c.viscousStress==ViscousStress2D::Symmetric
                    ?detail::symmetricViscousCorrection(mesh,u,v,gu,gv,b.u,b.v,b.fixedU,b.fixedV,b.constantU,b.constantV,c.nu)
                    :std::vector<Vector2D>{};
                momentum(au,mesh,c,b,u,flux,gu,gp,source,stressCorrection,{},false,nullptr,0);
                momentum(av,mesh,c,b,v,flux,gv,gp,source,stressCorrection,{},true,nullptr,0);
                const auto oldU=u,oldV=v;
                for (std::size_t i=0;i<n;++i) {
                    auto d=au.diag[i]; au.diag[i]/=c.velocityRelaxation; au.rhs[i]+=(au.diag[i]-d)*u[i];
                    d=av.diag[i]; av.diag[i]/=c.velocityRelaxation; av.rhs[i]+=(av.diag[i]-d)*v[i];
                }
                au.solve(u,work,3e-11); av.solve(v,work,3e-11);
                double change=0; std::size_t worst=0;
                for (std::size_t i=0;i<n;++i) {
                    const double e=std::hypot(u[i]-oldU[i],v[i]-oldV[i])/c.speed;
                    if (e>change) {change=e; worst=i;}
                }
                if (it==1 || it%10==0 || change<1e-8)
                    std::cout << it << ',' << change << ',' << worst << ',' << mesh.cells[worst].centre.x
                              << ',' << mesh.cells[worst].centre.y << ',' << mesh.cells[worst].area << '\n';
                if (change<1e-8) {std::cout << "diffusion_change_threshold " << it << '\n'; return 0;}
            }
            return 1;
        }
        const auto r=solveIncompressible2D(mesh,c,[&](const auto& s) {
            const auto i=s.momentumWorstCell;
            std::cout << s.iteration << ',' << s.momentumResidual << ',' << i << ','
                      << mesh.cells[i].centre.x << ',' << mesh.cells[i].centre.y << '\n';
        });
        std::cout << "strict_stokes_converged " << r.converged << ' ' << r.history.size() << '\n';
        return r.converged?0:1;
    } catch (const std::exception& e) {std::cerr << e.what() << '\n'; return 1;}
}
