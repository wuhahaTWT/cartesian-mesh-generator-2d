#include "cartmesh2d/fv/SstRans2D.hpp"
#include "cartmesh2d/io/MeshIO2D.hpp"
#include <cmath>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <set>
using namespace cartmesh2d;
using namespace cartmesh2d::fv;

namespace {
void writePerformance(std::ostream& meta,const SstRansResult2D& r) {
    meta<<"{\"sstUpdates\":"<<r.performance.updates
            <<",\"sstUpdateSeconds\":"<<r.performance.updateSeconds<<",\"sstTransportSeconds\":"<<r.performance.transportSeconds
            <<",\"sstGradientSeconds\":"<<r.performance.gradientSeconds<<",\"wallDistanceSeconds\":"<<r.performance.wallDistanceSeconds
            <<",\"momentumLinearSeconds\":"<<r.flow.performance.momentumLinearSolveSeconds
            <<",\"pressureLinearSeconds\":"<<r.flow.performance.pressureLinearSolveSeconds
            <<",\"momentumIterations\":"<<r.flow.performance.momentumIterations<<",\"pressureIterations\":"<<r.flow.performance.pressureIterations;
        const auto scalarProfile=[&](const char* label,const ScalarTransportPerformance2D& p) {
            meta<<",\""<<label<<"\":{\"calls\":"<<p.calls<<",\"patternBuilds\":"<<p.patternBuilds
                <<",\"ilu0Builds\":"<<p.ilu0Builds<<",\"ilu0Reuses\":"<<p.ilu0Reuses
                <<",\"linearIterations\":"<<p.linearIterations<<",\"totalSeconds\":"<<p.totalSeconds
                <<",\"setupSeconds\":"<<p.setupSeconds<<",\"linearSeconds\":"<<p.linearSeconds
                <<",\"faceFluxSeconds\":"<<p.faceFluxSeconds<<"}";
        };
        scalarProfile("scalarSolves",r.performance.scalarSolves);
        scalarProfile("scalarEvaluations",r.performance.scalarEvaluations);
        scalarProfile("kSolves",r.performance.kSolves);
        scalarProfile("omegaSolves",r.performance.omegaSolves);
        meta<<"}";

}
}
int main(int argc,char** argv) {
    try {
        if(argc<3)throw std::runtime_error("usage: sst_rans_probe mesh.cm2d prefix [nested|flatplate|flatplate-symmetry|flatplate-sweep|channel-sweep] [--nu value --speed value --inlet-k value --inlet-omega value --leading-edge x --max-iterations N --turbulence-updates N --scalar-preconditioner jacobi|ilu0]");
        const auto input=readCm2dTopology(argv[1]);if(!input.valid())throw std::runtime_error(input.error);
        const auto mesh=makeFvMesh2D(input.topology);
        SstRansControls2D c;c.flow.scenario="channel";c.flow.nu=.001;c.flow.tolerance=1e-7;
        c.flow.maxIterations=2000;c.flow.profile=true;
        int firstOption=3;
        if(argc>3 && std::string(argv[3]).rfind("--",0)!=0) {
            firstOption=4;
            if(std::string(argv[3])=="nested")c.turbulenceUpdatesPerIteration=c.turbulence.maxIterations;
            else if(std::string(argv[3])=="channel-sweep")c.turbulence.scalarCorrectionsPerUpdate=1;
            else if(std::string(argv[3])=="flatplate" || std::string(argv[3])=="flatplate-symmetry" || std::string(argv[3])=="flatplate-sweep") {
                c.flow.scenario="flatplate";c.flow.flatPlateLeadingEdge=.5;
                if(std::string(argv[3])=="flatplate-symmetry")c.flow.flatPlateTop=FlatPlateTop2D::Symmetry;
                if(std::string(argv[3])=="flatplate-sweep")c.turbulence.scalarCorrectionsPerUpdate=1;
            }
            else throw std::runtime_error("unknown SST probe configuration");
        }
        std::set<std::string> seen;
        for(int i=firstOption;i<argc;i+=2) {
            const std::string option=argv[i];
            if(i+1==argc || !seen.insert(option).second)throw std::runtime_error("missing or duplicate probe option: "+option);
            const std::string text=argv[i+1];
            if(option=="--scalar-preconditioner") {
                if(text=="jacobi")c.turbulence.transport.preconditioner=ScalarPreconditioner2D::Jacobi;
                else if(text=="ilu0")c.turbulence.transport.preconditioner=ScalarPreconditioner2D::ILU0;
                else throw std::runtime_error("unknown scalar preconditioner");
                continue;
            }
            std::size_t end=0;const double value=std::stod(text,&end);
            if(end!=text.size()||!std::isfinite(value))throw std::runtime_error("invalid finite probe value: "+option);
            if(option=="--nu")c.flow.nu=value;
            else if(option=="--speed")c.flow.speed=value;
            else if(option=="--inlet-k")c.inletK=value;
            else if(option=="--inlet-omega")c.inletOmega=value;
            else if(option=="--max-iterations") {
                if(value<1 || value>2000 || std::floor(value)!=value)throw std::runtime_error("max-iterations must be an integer in 1..2000");
                c.flow.maxIterations=static_cast<std::size_t>(value);
            }
            else if(option=="--turbulence-updates") {
                if(value<1 || value>500 || std::floor(value)!=value)throw std::runtime_error("turbulence-updates must be an integer in 1..500");
                c.turbulenceUpdatesPerIteration=static_cast<std::size_t>(value);
            }
            else if(option=="--leading-edge" && c.flow.scenario=="flatplate")c.flow.flatPlateLeadingEdge=value;
            else throw std::runtime_error("unknown or inapplicable probe option: "+option);
        }
        if(!(c.flow.nu>0 && c.flow.speed>0 && c.inletK>=0 && c.inletOmega>0))
            throw std::runtime_error("probe requires positive nu/speed/omega and nonnegative k");
        const std::string prefix=argv[2];
        for(const auto suffix:{".json",".diagnostics.json",".cells.csv",".faces.csv",".history.csv",".unconverged.json",".unconverged.cells.csv",".unconverged.faces.csv"})
            if(std::filesystem::exists(prefix+suffix))throw std::runtime_error("probe requires a fresh output prefix");
        const auto r=solveSstRans2D(mesh,c,{}, {},[](const auto& h){
            if(h.iteration==1||h.iteration%100==0)std::cerr<<h.iteration<<" momentum="<<h.momentumResidual<<" cell="<<h.momentumWorstCell
                <<" mx="<<h.momentumResidualX<<" my="<<h.momentumResidualY
                <<" predictor="<<h.momentumPredictorResidual<<" predictorCell="<<h.momentumPredictorWorstCell<<" pressureLinear="<<h.pressureLinearResidual<<'\n';
        });
        std::ofstream history(prefix+".history.csv");history<<std::setprecision(17);
        history<<"iteration,momentumResidual,continuity,velocityChange,pressureChange,kNorm,omegaNorm,kCellResidual,omegaCellResidual,turbulenceIterations,momentumWorstCell,momentumResidualX,momentumResidualY,momentumPredictorResidual,momentumPredictorWorstCell,pressureLinearResidual\n";
        for(std::size_t i=0;i<r.history.size();++i) {
            const auto& f=r.flow.history[i];const auto& h=r.history[i];
            history<<f.iteration<<','<<f.momentumResidual<<','<<f.continuity<<','<<f.velocityChange<<','<<f.pressureChange<<','
                <<h.kNorm<<','<<h.omegaNorm<<','<<h.kCellResidual<<','<<h.omegaCellResidual<<','<<h.turbulenceIterations<<','
                <<f.momentumWorstCell<<','<<f.momentumResidualX<<','<<f.momentumResidualY<<','<<f.momentumPredictorResidual<<','<<f.momentumPredictorWorstCell<<','<<f.pressureLinearResidual<<'\n';
        }
        std::ofstream diagnostics(prefix+".diagnostics.json");
        diagnostics<<std::setprecision(17)<<"{\"converged\":"<<(r.converged?"true":"false")
            <<",\"stopReason\":\""<<(r.converged?"converged":"iteration-limit")<<"\",\"maxIterations\":"<<c.flow.maxIterations
            <<",\"iterations\":"<<r.history.size()<<",\"cells\":"<<mesh.cells.size()
            <<",\"scalarPreconditioner\":\""<<(c.turbulence.transport.preconditioner==ScalarPreconditioner2D::ILU0?"ilu0":"jacobi")<<"\""
            <<",\"solveSeconds\":"<<r.flow.performance.solveSeconds<<",\"performance\":";
        writePerformance(diagnostics,r);
        diagnostics<<"}\n";
        if(!history || !diagnostics)throw std::runtime_error("cannot write RANS history/diagnostics");
        const std::string fieldPrefix=prefix+(r.converged?"":".unconverged");
        const auto& flow=r.flow;const auto& t=r.turbulence.fields;
        FrozenSst2003mProblem2D p;p.k=t.k.values;p.omega=t.omega.values;
        p.boundaryK=r.boundaryK;p.boundaryOmega=r.boundaryOmega;
        std::vector<Vector2D> velocity;for(std::size_t i=0;i<p.k.size();++i)velocity.push_back({flow.u[i],flow.v[i]});
        const auto g=reconstructSst2003mGradients2D(mesh,p,velocity,r.velocityBoundary);
        std::ofstream cells(fieldPrefix+".cells.csv"),faces(fieldPrefix+".faces.csv"),meta(fieldPrefix+".json");
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
        meta<<std::setprecision(17)<<"{\"case\":\""<<c.flow.scenario<<"\",\"flatPlateLeadingEdge\":"<<c.flow.flatPlateLeadingEdge
            <<",\"flatPlateTop\":\""<<(c.flow.flatPlateTop==FlatPlateTop2D::Symmetry?"symmetry":"pressure-farfield")<<"\""
            <<",\"model\":\"SST-2003m\",\"scope\":\"coupled-steady-SST-2003m\",\"converged\":"<<(r.converged?"true":"false")<<",\"nu\":"<<c.flow.nu
            <<",\"speed\":"<<c.flow.speed<<",\"inletK\":"<<c.inletK<<",\"inletOmega\":"<<c.inletOmega<<",\"tolerance\":1e-7,"
            <<"\"scalarRelativeTolerance\":1e-9,\"scalarAbsoluteTolerance\":1e-12,\"scalarCellTolerance\":1e-9,\"convection\":\"upwind\",\"viscousStress\":\"symmetric\","
            <<"\"pressureConvention\":\"p/rho (SST-2003m omits isotropic k stress)\",\"pressureDiscretization\":\"shared-face-gauss\",\"iterations\":"<<r.history.size()<<",\"cells\":"<<mesh.cells.size()
            <<",\"turbulenceUpdatesPerIteration\":"<<c.turbulenceUpdatesPerIteration
            <<",\"scalarCorrectionsPerUpdate\":"<<c.turbulence.scalarCorrectionsPerUpdate
            <<",\"scalarPreconditioner\":\""<<(c.turbulence.transport.preconditioner==ScalarPreconditioner2D::ILU0?"ilu0":"jacobi")<<"\""
            <<",\"solveSeconds\":"<<flow.performance.solveSeconds
            <<",\"performance\":";
        writePerformance(meta,r);
        meta            <<",\"momentumResidual\":"<<flow.history.back().momentumResidual
            <<",\"forceX\":"<<flow.forceX<<",\"forceY\":"<<flow.forceY
            <<",\"pressureForceX\":"<<flow.pressureForceX<<",\"pressureForceY\":"<<flow.pressureForceY
            <<",\"discreteForceX\":"<<flow.discreteForceX<<",\"discreteForceY\":"<<flow.discreteForceY
            <<",\"reconstructedForceX\":"<<flow.reconstructedForceX<<",\"reconstructedForceY\":"<<flow.reconstructedForceY
            <<",\"wallForceX\":"<<flow.wallForceX<<",\"wallForceY\":"<<flow.wallForceY
            <<",\"wallViscousForceX\":"<<flow.wallViscousForceX<<",\"wallViscousForceY\":"<<flow.wallViscousForceY
            <<",\"globalRelativeImbalance\":"<<flow.globalRelativeImbalance<<"}\n";
        if(!r.converged)throw std::runtime_error("coupled RANS did not converge; diagnostic-only .unconverged fields retained; no accepted fields");
        std::cout<<"coupled cells="<<mesh.cells.size()<<" iterations="<<r.history.size()<<" seconds="<<flow.performance.solveSeconds<<'\n';return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
