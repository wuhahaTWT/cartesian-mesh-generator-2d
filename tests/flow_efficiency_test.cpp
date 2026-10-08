#include "cartmesh2d/fv/Incompressible2D.hpp"
#include "cartmesh2d/fv/detail/FlowConvergence2D.hpp"
#include "cartmesh2d/grid/RectilinearMesh2D.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace cartmesh2d;
using namespace cartmesh2d::fv;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
FvMesh2D rectangle(std::size_t nx,std::size_t ny,double length){
    std::vector<double> x(nx+1),y(ny+1);
    for(std::size_t i=0;i<=nx;++i)x[i]=length*static_cast<double>(i)/static_cast<double>(nx);
    for(std::size_t i=0;i<=ny;++i)y[i]=static_cast<double>(i)/static_cast<double>(ny);
    return makeFvMesh2D(makeRectilinearMesh2D(x,y));
}
double drop(const FvMesh2D& mesh,const FlowResult2D& r){
    double left=0,right=0,ll=0,lr=0;
    for(std::size_t id=0;id<mesh.faces.size();++id){
        const auto& f=mesh.faces[id];if(f.neighbour)continue;
        if(f.areaVector.x<0){left-=f.areaVector.x*r.faceMomentum[id].pressure;ll-=f.areaVector.x;}
        if(f.areaVector.x>0){right+=f.areaVector.x*r.faceMomentum[id].pressure;lr+=f.areaVector.x;}
    }
    return left/ll-right/lr;
}
int main(){
    const auto mesh=rectangle(64,32,2);
    FlowControls2D c;c.scenario="channel";c.nu=.1;c.tolerance=1e-9;c.maxIterations=5000;
    c.convection=ConvectionScheme2D::LimitedLinearUpwind;
    const auto reference=solveIncompressible2D(mesh,c);
    require(reference.converged,"strict channel reference failed");
    auto engineering=c;engineering.convergence=FlowConvergence2D::Engineering;engineering.tolerance=1e-5;
    const auto early=solveIncompressible2D(mesh,engineering);
    require(early.converged && early.history.back().strictLinearStep,"engineering final step failed");
    require(early.history.size()<reference.history.size(),"engineering did not reduce cost");
    // Same-mesh iteration error, not a new physical discretization gate:
    // 0.1% of a nonzero pressure drop and wall drag is the engineering target.
    require(std::abs(drop(mesh,early)/drop(mesh,reference)-1)<1e-3,"engineering pressure drop exceeds 0.1%");
    require(std::abs(early.wallForceX/reference.wallForceX-1)<1e-3,"engineering drag exceeds 0.1%");
    require(std::any_of(early.monitorNames.begin(),early.monitorNames.end(),[](const auto& s){return s.find("mean-pressure")!=std::string::npos;}),"pressure monitor missing");
    engineering.maxIterations=49;
    require(!solveIncompressible2D(mesh,engineering).converged,"iteration cap bypassed engineering window");
    auto simplec=c;simplec.coupling=FlowCoupling2D::SimpleC;simplec.pressureRelaxation=1;
    const auto corrected=solveIncompressible2D(mesh,simplec);
    require(corrected.converged,"SIMPLEC channel failed");
    double difference=0;
    for(std::size_t i=0;i<mesh.cells.size();++i)
        difference=std::max({difference,std::abs(corrected.u[i]-reference.u[i]),std::abs(corrected.v[i]-reference.v[i]),std::abs(corrected.p[i]-reference.p[i])});
    require(difference<1e-6,"SIMPLEC changed the converged fixed point");
    const auto cavity=rectangle(32,32,1);
    engineering=c;engineering.scenario="cavity";engineering.nu=.01;
    engineering.convergence=FlowConvergence2D::Engineering;engineering.tolerance=1e-5;
    const auto box=solveIncompressible2D(cavity,engineering);
    require(box.converged && box.history.back().strictLinearStep,"engineering cavity failed");
    const auto& h=box.history.back();
    require(!detail::strictFlowConverged2D(h.iteration,h.momentumResidual,h.velocityChange,h.pressureChange,
        h.continuity,box.globalRelativeImbalance,engineering.tolerance,true),"engineering still depends on every strict stop");
    std::cout<<"channel strict="<<reference.history.size()<<" engineering="<<early.history.size()
             <<" SIMPLEC="<<corrected.history.size()<<" fieldDifference="<<difference<<'\n';
}
