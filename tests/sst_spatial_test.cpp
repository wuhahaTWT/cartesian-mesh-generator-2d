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
namespace {
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
double kValue(Point2D p){return .02*(1+.2*p.x+.3*p.y);}
double wValue(Point2D p){return 4*(1+.1*p.x+.2*p.y);}
void run(const FvMesh2D& m,const std::string& prefix) {
    FrozenSst2003mProblem2D p;std::vector<Vector2D> velocity;
    std::vector<SstVelocityBoundary2D> velocityBC;
    std::vector<bool> walls;
    for(const auto& c:m.cells){p.k.push_back(kValue(c.centre));p.omega.push_back(wValue(c.centre));velocity.push_back({c.centre.y,0});}
    for(const auto& f:m.faces) {
        p.volumeFlux.push_back(f.centre.y*f.areaVector.x);
        p.boundaryK.push_back({ScalarBoundaryKind2D::Value,kValue(f.centre),{}});
        p.boundaryOmega.push_back({ScalarBoundaryKind2D::Value,wValue(f.centre),{}});
        velocityBC.push_back({{f.centre.y,0},true,true});
        walls.push_back(!f.neighbour&&std::abs(f.centre.y)<1e-12);
    }
    const auto distance=computeWallDistance2D(m,walls);
    const auto gradients=reconstructSst2003mGradients2D(m,p,velocity,velocityBC);
    p.wallDistance=distance.distance;p.gradientK=gradients.k;p.gradientOmega=gradients.omega;
    p.strainMagnitude=gradients.strainMagnitude;
    for(std::size_t i=0;i<m.cells.size();++i) {
        require(std::abs(gradients.k[i].x-.004)<1e-11&&std::abs(gradients.k[i].y-.006)<1e-11,"affine k gradient");
        require(std::abs(gradients.omega[i].x-.4)<1e-10&&std::abs(gradients.omega[i].y-.8)<1e-10,"affine omega gradient");
        require(std::abs(gradients.strainMagnitude[i]-1)<1e-11,"shear strain tensor contraction");
    }
    // Rigid rotation u=-y,v=x has vorticity but ZERO symmetric strain.
    auto rotation=velocity;auto rotationBC=velocityBC;
    for(std::size_t i=0;i<rotation.size();++i)rotation[i]={-m.cells[i].centre.y,m.cells[i].centre.x};
    for(std::size_t id=0;id<rotationBC.size();++id)rotationBC[id]={{-m.faces[id].centre.y,m.faces[id].centre.x},true,true};
    const auto rg=reconstructSst2003mGradients2D(m,p,rotation,rotationBC);
    for(double s:rg.strainMagnitude)require(s<1e-10,"rigid rotation cannot produce strain");
    auto invalid=p;
    for(std::size_t id=0;id<m.faces.size();++id)if(!m.faces[id].neighbour) {
        invalid.boundaryK[id]={ScalarBoundaryKind2D::DiffusiveFlux,.1,{}};break;
    }
    bool rejected=false;
    try{(void)reconstructSst2003mGradients2D(m,invalid,velocity,velocityBC);}catch(const std::exception&){rejected=true;}
    require(rejected,"nonzero unknown-coefficient scalar flux accepted by gradient helper");
    ScalarTransportControls2D controls;controls.relativeTolerance=1e-11;controls.absoluteTolerance=1e-13;controls.cellTolerance=1e-11;
    const auto result=solveFrozenSst2003mTransport2D(m,p,controls,p.k,p.omega,.01);
    std::ofstream cells(prefix+".cells.csv"),faces(prefix+".faces.csv"),meta(prefix+".json");
    require(bool(cells)&&bool(faces)&&bool(meta),"output open failed");
    cells<<std::setprecision(17)<<"cell,x,y,area,distance,nearestFace,oldK,oldOmega,gradKx,gradKy,gradWx,gradWy,strain,F1,F2,nuT,Dk,Dw,sourceK,sourceW,lossK,lossW,k,omega\n";
    for(std::size_t i=0;i<m.cells.size();++i) {
        const auto& c=m.cells[i];const auto& q=result.coefficients[i];
        cells<<i<<','<<c.centre.x<<','<<c.centre.y<<','<<c.area<<','<<distance.distance[i]<<','<<distance.nearestFace[i]<<','
            <<p.k[i]<<','<<p.omega[i]<<','<<gradients.k[i].x<<','<<gradients.k[i].y<<','<<gradients.omega[i].x<<','<<gradients.omega[i].y<<','
            <<gradients.strainMagnitude[i]<<','<<q.f1<<','<<q.f2<<','<<q.turbulentViscosity<<','<<q.diffusivityK<<','<<q.diffusivityOmega<<','
            <<q.sourceK<<','<<q.sourceOmega<<','<<q.lossRateK<<','<<q.lossRateOmega<<','<<result.k.values[i]<<','<<result.omega.values[i]<<'\n';
    }
    faces<<std::setprecision(17)<<"face,wall,volumeFlux,kAdvection,kDiffusion,omegaAdvection,omegaDiffusion\n";
    for(std::size_t id=0;id<m.faces.size();++id)
        faces<<id<<','<<walls[id]<<','<<p.volumeFlux[id]<<','<<result.k.advectiveFlux[id]<<','<<result.k.diffusiveFlux[id]<<','
            <<result.omega.advectiveFlux[id]<<','<<result.omega.diffusiveFlux[id]<<'\n';
    meta<<std::setprecision(17)<<"{\"scope\":\"one-frozen-SST-transport-iteration\",\"nu\":1e-5,\"dt\":0.01,\"segmentTests\":"<<distance.segmentTests
        <<",\"cells\":"<<m.cells.size()<<",\"kInnerConverged\":"<<(result.k.converged?"true":"false")
        <<",\"omegaInnerConverged\":"<<(result.omega.converged?"true":"false")<<"}\n";
    std::cout<<"cells="<<m.cells.size()<<" segmentTests="<<distance.segmentTests<<'\n';
}
}
int main(int argc,char** argv) {
    try {
        require(argc==3,"usage: sst_spatial_tests mesh.cm2d outputPrefix");
        const auto read=readCm2dTopology(argv[1]);require(read.valid(),"invalid final mesh");
        run(makeFvMesh2D(read.topology),argv[2]);return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
