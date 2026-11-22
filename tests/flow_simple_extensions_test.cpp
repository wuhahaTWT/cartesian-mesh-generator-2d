#include "cartmesh2d/fv/Incompressible2D.hpp"
#include "cartmesh2d/fv/detail/FlowConvergence2D.hpp"
#include "cartmesh2d/grid/RectilinearMesh2D.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

using namespace cartmesh2d;
using namespace cartmesh2d::fv;
namespace {
void require(bool ok,const std::string& message) {if(!ok)throw std::runtime_error(message);}
FvMesh2D mesh() {
    // Grading exercises nonuniform response/interpolation; slip adds unequal
    // component diagonals. Actual polygonal Cut-cells are exercised by CLI runs.
    std::vector<double> x(13),y(9);
    for(std::size_t i=0;i<x.size();++i)x[i]=std::pow(static_cast<double>(i)/12,1.2);
    for(std::size_t i=0;i<y.size();++i)y[i]=std::pow(static_cast<double>(i)/8,1.15);
    return makeFvMesh2D(makeRectilinearMesh2D(x,y));
}
void strict(const FlowResult2D& r,const FlowControls2D& c) {
    require(r.converged,"extension failed to converge");const auto& h=r.history.back();
    require(h.strictLinearStep && detail::strictFlowConverged2D(h.iteration,h.momentumResidual,
        h.velocityChange,h.pressureChange,h.continuity,r.globalRelativeImbalance,c.tolerance,true),
        "extension bypassed the existing strict gates");
}
double difference(const FlowResult2D& a,const FlowResult2D& b) {
    double e=0;
    for(std::size_t i=0;i<a.u.size();++i)e=std::max({e,std::abs(a.u[i]-b.u[i]),std::abs(a.v[i]-b.v[i]),std::abs(a.p[i]-b.p[i])});
    for(std::size_t i=0;i<a.flux.size();++i)e=std::max(e,std::abs(a.flux[i]-b.flux[i]));
    e=std::max({e,std::abs(a.wallForceX-b.wallForceX),std::abs(a.wallForceY-b.wallForceY)});
    return e;
}
template<class F> void rejects(F f) {
    bool rejected=false;try{f();}catch(const std::runtime_error&){rejected=true;}
    require(rejected,"unsupported extension controls were silently accepted");
}
}
int main() {
    const auto m=mesh();FlowControls2D c;c.scenario="channel";c.nu=.1;
    c.tolerance=1e-9;c.maxIterations=2000;c.profile=true;
    c.convection=ConvectionScheme2D::FaceLimitedLinearUpwind;
    for(const auto& scenario:{"channel","cavity","manufactured","custom"}) {
        c.scenario=scenario;c.boundaryConditions.clear();
        c.manufacturedPressureSlope=c.scenario=="manufactured"?.7:0;
        if(c.scenario=="custom")for(std::size_t id=0;id<m.faces.size();++id) {
            const auto& f=m.faces[id];if(f.neighbour)continue;
            FlowBoundaryCondition2D bc;bc.face=id;bc.name="wall";
            if(f.areaVector.x<0) {bc.kind=FlowBoundaryKind2D::VelocityInlet;bc.velocity={1,0};bc.name="inlet";}
            else if(f.areaVector.x>0) {bc.kind=FlowBoundaryKind2D::PressureOutlet;bc.pressure=1.7;bc.name="outlet";}
            else if(f.areaVector.y>0) {bc.kind=FlowBoundaryKind2D::Symmetry;bc.name="slip";}
            c.boundaryConditions.push_back(bc);
        }
        c.coupling=FlowCoupling2D::Simple;c.momentumCorrectors=2;c.adaptiveLinear=false;
        const auto reference=solveIncompressible2D(m,c);strict(reference,c);
        for(const auto method:{FlowCoupling2D::Simpler,FlowCoupling2D::Pimple}) {
            c.coupling=method;
            for(const auto relaxation:{std::pair{.6,.25},std::pair{.7,.6}}) {
                c.velocityRelaxation=relaxation.first;c.pressureRelaxation=relaxation.second;
                c.adaptiveLinear=relaxation.second==.6;
                const auto r=solveIncompressible2D(m,c);strict(r,c);
                // Reuse the coupled test's 2e-6 absolute iteration-error budget
                // for u/v, kinematic p, conservative flux and integrated load.
                // This compares identical equations, not spatial accuracy.
                const double e=difference(r,reference);
                std::cout<<scenario<<' '<<(method==FlowCoupling2D::Simpler?"simpler":"pimple")
                         <<" iterations="<<r.history.size()<<" difference="<<e<<std::endl;
                require(e<2e-6,"extension changed the discrete fixed point");
                require(r.performance.momentumSolves==2*r.history.size(),"inner correctors secretly solve momentum");
                require(r.performance.pressureSolves+r.performance.pressureCorrectionPassesSkipped==8*r.history.size(),
                        "additional pressure work is missing from profiling");
            }
        }
        // One pressure/momentum corrector is exactly the original SIMPLE map.
        c.coupling=FlowCoupling2D::Simple;c.velocityRelaxation=.6;c.pressureRelaxation=.25;c.adaptiveLinear=false;
        c.maxIterations=1;const auto first=solveIncompressible2D(m,c);
        c.coupling=FlowCoupling2D::Pimple;c.momentumCorrectors=1;
        const auto one=solveIncompressible2D(m,c);
        require(difference(first,one)==0,"one-corrector PIMPLE did not recover SIMPLE exactly");
        require(!one.converged,"iteration cap was mistaken for convergence");
        c.momentumCorrectors=2;c.maxIterations=2000;
    }
    c.scenario="cavity";c.boundaryConditions.clear();c.coupling=FlowCoupling2D::Simple;
    const auto initial=initialIncompressibleState2D(m,c);
    const auto simple=advanceIncompressible2D(m,c,initial,.03);strict(simple,c);
    c.coupling=FlowCoupling2D::Pimple;c.adaptiveLinear=true;
    const auto pimple=advanceIncompressible2D(m,c,initial,.03);strict(pimple,c);
    require(difference(simple,pimple)<2e-6,"PIMPLE changed backward Euler/old face-flux transport");
    c.coupling=FlowCoupling2D::Simpler;
    rejects([&]{(void)advanceIncompressible2D(m,c,initial,.03);});
    c.coupling=FlowCoupling2D::Pimple;c.momentumCorrectors=9;
    rejects([&]{(void)solveIncompressible2D(m,c);});
    c.momentumCorrectors=2;c.outletBackflow=OutletBackflow2D::NormalInlet;
    rejects([&]{(void)solveIncompressible2D(m,c);});
    c.outletBackflow=OutletBackflow2D::Reject;c.scenario="flatplate";
    rejects([&]{(void)solveIncompressible2D(m,c);});
    c.scenario="cavity";c.stopRequested=[] {return true;};
    const auto stopped=solveIncompressible2D(m,c);
    require(stopped.stopped && !stopped.converged && stopped.history.size()==1,"cancel interrupted a corrector group");
}
