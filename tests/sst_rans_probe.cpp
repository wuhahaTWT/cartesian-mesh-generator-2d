#include "cartmesh2d/fv/SstRans2D.hpp"
#include "cartmesh2d/io/MeshIO2D.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
using namespace cartmesh2d;
using namespace cartmesh2d::fv;
int main(int argc,char** argv) {
    try {
        if(argc!=3&&argc!=4)throw std::runtime_error("usage: sst_rans_probe mesh.cm2d prefix [nested]");
        const auto input=readCm2dTopology(argv[1]);if(!input.valid())throw std::runtime_error(input.error);
        const auto mesh=makeFvMesh2D(input.topology);
        SstRansControls2D c;c.flow.scenario="channel";c.flow.nu=.001;c.flow.tolerance=1e-7;
        c.flow.maxIterations=2000;c.flow.profile=true;
        if(argc==4) {
            if(std::string(argv[3])!="nested")throw std::runtime_error("unknown SST coupling strategy");
            c.turbulenceUpdatesPerIteration=c.turbulence.maxIterations;
        }
        const std::string prefix=argv[2];
        const auto r=solveSstRans2D(mesh,c,{}, {},[](const auto& h){
            if(h.iteration==1||h.iteration%100==0)std::cerr<<h.iteration<<" momentum="<<h.momentumResidual<<'\n';
        });
        std::ofstream history(prefix+".history.csv");history<<std::setprecision(17);
        history<<"iteration,momentumResidual,continuity,velocityChange,pressureChange,kNorm,omegaNorm,kCellResidual,omegaCellResidual,turbulenceIterations\n";
        for(std::size_t i=0;i<r.history.size();++i) {
            const auto& f=r.flow.history[i];const auto& h=r.history[i];
            history<<f.iteration<<','<<f.momentumResidual<<','<<f.continuity<<','<<f.velocityChange<<','<<f.pressureChange<<','
                <<h.kNorm<<','<<h.omegaNorm<<','<<h.kCellResidual<<','<<h.omegaCellResidual<<','<<h.turbulenceIterations<<'\n';
        }
        if(!r.converged)throw std::runtime_error("coupled RANS did not converge; history retained");
        const auto& flow=r.flow;const auto& t=r.turbulence.fields;
        FrozenSst2003mProblem2D p;p.k=t.k.values;p.omega=t.omega.values;
        p.boundaryK=r.boundaryK;p.boundaryOmega=r.boundaryOmega;
        std::vector<Vector2D> velocity;for(std::size_t i=0;i<p.k.size();++i)velocity.push_back({flow.u[i],flow.v[i]});
        const auto g=reconstructSst2003mGradients2D(mesh,p,velocity,r.velocityBoundary);
        std::ofstream cells(prefix+".cells.csv"),faces(prefix+".faces.csv"),meta(prefix+".json");
        if(!history||!cells||!faces||!meta)throw std::runtime_error("cannot write RANS evidence");
        cells<<std::setprecision(17)<<"cell,x,y,area,u,v,p,speed,k,omega,distance,gradKx,gradKy,gradWx,gradWy,strain,F1,F2,nuT,Dk,Dw,sourceK,sourceW,lossK,lossW\n";
        for(std::size_t i=0;i<p.k.size();++i) {
            const auto& cell=mesh.cells[i];const auto& q=t.coefficients[i];
            cells<<i<<','<<cell.centre.x<<','<<cell.centre.y<<','<<cell.area<<','<<flow.u[i]<<','<<flow.v[i]<<','<<flow.p[i]<<','<<std::hypot(flow.u[i],flow.v[i])<<','
                <<p.k[i]<<','<<p.omega[i]<<','<<r.wallDistance[i]<<','<<g.k[i].x<<','<<g.k[i].y<<','<<g.omega[i].x<<','<<g.omega[i].y<<','<<g.strainMagnitude[i]<<','
                <<q.f1<<','<<q.f2<<','<<q.turbulentViscosity<<','<<q.diffusivityK<<','<<q.diffusivityOmega<<','<<q.sourceK<<','<<q.sourceOmega<<','<<q.lossRateK<<','<<q.lossRateOmega<<'\n';
        }
        faces<<std::setprecision(17)<<"face,owner,neighbour,flux,pressure,advectionX,advectionY,diffusionX,diffusionY,viscosity,wall,kBoundary,omegaBoundary,kAdvection,kDiffusion,omegaAdvection,omegaDiffusion\n";
        for(std::size_t id=0;id<mesh.faces.size();++id) {
            const auto& f=mesh.faces[id];const auto& q=flow.faceMomentum[id];
            faces<<id<<','<<f.owner<<',';if(f.neighbour)faces<<*f.neighbour;else faces<<-1;
            faces<<','<<flow.flux[id]<<','<<q.pressure<<','<<q.advection.x<<','<<q.advection.y<<','<<q.diffusion.x<<','<<q.diffusion.y<<','<<r.faceViscosity[id]<<','<<r.resolvedWalls[id]<<','
                <<r.boundaryK[id].value<<','<<r.boundaryOmega[id].value<<','<<t.k.advectiveFlux[id]<<','<<t.k.diffusiveFlux[id]<<','<<t.omega.advectiveFlux[id]<<','<<t.omega.diffusiveFlux[id]<<'\n';
        }
        meta<<std::setprecision(17)<<"{\"case\":\"channel\",\"model\":\"SST-2003m\",\"scope\":\"coupled-steady-SST-2003m\",\"converged\":true,\"nu\":0.001,\"speed\":1,\"inletK\":0.001,\"inletOmega\":2,\"tolerance\":1e-7,"
            <<"\"scalarRelativeTolerance\":1e-9,\"scalarAbsoluteTolerance\":1e-12,\"scalarCellTolerance\":1e-9,\"convection\":\"upwind\",\"viscousStress\":\"symmetric\","
            <<"\"pressureConvention\":\"p/rho (SST-2003m omits isotropic k stress)\",\"pressureDiscretization\":\"shared-face-gauss\",\"iterations\":"<<r.history.size()<<",\"cells\":"<<mesh.cells.size()
            <<",\"turbulenceUpdatesPerIteration\":"<<c.turbulenceUpdatesPerIteration
            <<",\"solveSeconds\":"<<flow.performance.solveSeconds
            <<",\"momentumResidual\":"<<flow.history.back().momentumResidual
            <<",\"forceX\":"<<flow.forceX<<",\"forceY\":"<<flow.forceY
            <<",\"pressureForceX\":"<<flow.pressureForceX<<",\"pressureForceY\":"<<flow.pressureForceY
            <<",\"discreteForceX\":"<<flow.discreteForceX<<",\"discreteForceY\":"<<flow.discreteForceY
            <<",\"reconstructedForceX\":"<<flow.reconstructedForceX<<",\"reconstructedForceY\":"<<flow.reconstructedForceY
            <<",\"wallForceX\":"<<flow.wallForceX<<",\"wallForceY\":"<<flow.wallForceY
            <<",\"wallViscousForceX\":"<<flow.wallViscousForceX<<",\"wallViscousForceY\":"<<flow.wallViscousForceY
            <<",\"globalRelativeImbalance\":"<<flow.globalRelativeImbalance<<"}\n";
        std::cout<<"coupled cells="<<mesh.cells.size()<<" iterations="<<r.history.size()<<" seconds="<<flow.performance.solveSeconds<<'\n';return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
