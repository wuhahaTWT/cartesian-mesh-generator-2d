#include "cartmesh2d/fv/Sst2003m2D.hpp"
#include "cartmesh2d/fv/WallDistance2D.hpp"
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
        if(argc!=3)throw std::runtime_error("usage: sst_nonlinear_probe mesh.cm2d outputPrefix");
        const auto read=readCm2dTopology(argv[1]);if(!read.valid())throw std::runtime_error(read.error);
        const auto mesh=makeFvMesh2D(read.topology);const auto n=mesh.cells.size(),nf=mesh.faces.size();
        FrozenSst2003mProblem2D p;p.k.assign(n,.001);p.omega.assign(n,2);
        std::vector<Vector2D> velocity;std::vector<SstVelocityBoundary2D> velocityBC;
        std::vector<bool> walls(nf);
        for(const auto& c:mesh.cells)velocity.push_back({c.centre.y,0});
        for(std::size_t id=0;id<nf;++id) {
            const auto& f=mesh.faces[id];
            walls[id]=!f.neighbour&&std::abs(f.centre.y)<1e-12;
            // The selected wall is an impermeable no-slip boundary. Do not
            // sample y*nx at roundoff-displaced polygon vertices as wall flux.
            const double flux=walls[id]?0:f.centre.y*f.areaVector.x;
            p.volumeFlux.push_back(flux);velocityBC.push_back({{walls[id]?0:f.centre.y,0},true,true});
            const bool outflow=!f.neighbour&&flux>0;
            p.boundaryK.push_back(outflow?ScalarBoundary2D{ScalarBoundaryKind2D::DiffusiveFlux,0,{}}:
                ScalarBoundary2D{ScalarBoundaryKind2D::Value,.001,{}});
            p.boundaryOmega.push_back(outflow?ScalarBoundary2D{ScalarBoundaryKind2D::DiffusiveFlux,0,{}}:
                ScalarBoundary2D{ScalarBoundaryKind2D::Value,2,{}});
        }
        p.wallDistance=computeWallDistance2D(mesh,walls).distance;
        setSst2003mResolvedWalls2D(mesh,p,walls);
        SstTransportControls2D controls;controls.transport.relativeTolerance=1e-10;
        controls.transport.absoluteTolerance=1e-13;controls.transport.cellTolerance=1e-10;
        const auto result=solveSst2003mTransport2D(mesh,p,velocity,velocityBC,controls,p.k,p.omega,.02);
        const std::string prefix=argv[2];
        std::ofstream history(prefix+".history.csv");history<<std::setprecision(17);
        history<<"iteration,kNorm,omegaNorm,kCellResidual,omegaCellResidual\n";
        for(const auto& h:result.history)history<<h.iteration<<','<<h.kResidualNorm<<','<<h.omegaResidualNorm<<','<<h.kCellResidual<<','<<h.omegaCellResidual<<'\n';
        if(!result.converged)throw std::runtime_error("nonlinear SST transport did not converge; history retained");
        auto final=p;final.k=result.fields.k.values;final.omega=result.fields.omega.values;
        const auto gradients=reconstructSst2003mGradients2D(mesh,final,velocity,velocityBC);
        std::ofstream cells(prefix+".cells.csv"),faces(prefix+".faces.csv"),meta(prefix+".json");
        if(!cells||!faces||!meta||!history)throw std::runtime_error("cannot write evidence files");
        cells<<std::setprecision(17)<<"cell,x,y,area,distance,k,omega,gradKx,gradKy,gradWx,gradWy,strain,F1,F2,nuT,Dk,Dw,sourceK,sourceW,lossK,lossW\n";
        for(std::size_t i=0;i<n;++i) {
            const auto& c=mesh.cells[i];const auto& q=result.fields.coefficients[i];
            cells<<i<<','<<c.centre.x<<','<<c.centre.y<<','<<c.area<<','<<p.wallDistance[i]<<','<<final.k[i]<<','<<final.omega[i]<<','
                <<gradients.k[i].x<<','<<gradients.k[i].y<<','<<gradients.omega[i].x<<','<<gradients.omega[i].y<<','
                <<gradients.strainMagnitude[i]<<','<<q.f1<<','<<q.f2<<','<<q.turbulentViscosity<<','<<q.diffusivityK<<','<<q.diffusivityOmega<<','
                <<q.sourceK<<','<<q.sourceOmega<<','<<q.lossRateK<<','<<q.lossRateOmega<<'\n';
        }
        faces<<std::setprecision(17)<<"face,wall,volumeFlux,kBoundary,omegaBoundary,kAdvection,kDiffusion,omegaAdvection,omegaDiffusion\n";
        for(std::size_t id=0;id<nf;++id)faces<<id<<','<<walls[id]<<','<<p.volumeFlux[id]<<','<<p.boundaryK[id].value<<','<<p.boundaryOmega[id].value<<','
            <<result.fields.k.advectiveFlux[id]<<','<<result.fields.k.diffusiveFlux[id]<<','<<result.fields.omega.advectiveFlux[id]<<','<<result.fields.omega.diffusiveFlux[id]<<'\n';
        meta<<"{\"scope\":\"nonlinear-SST-fixed-carrier\",\"converged\":true,\"nu\":1e-5,\"dt\":0.02,\"initialK\":0.001,\"initialOmega\":2,\"iterations\":"
            <<result.history.back().iteration<<",\"cells\":"<<n<<"}\n";
        std::cout<<"cells="<<n<<" nonlinear iterations="<<result.history.back().iteration<<'\n';return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
