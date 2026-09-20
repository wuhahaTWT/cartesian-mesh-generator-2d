#include "cartmesh2d/fv/SstRans2D.hpp"
#include "cartmesh2d/fv/detail/FlowMaterial2D.hpp"
#include "cartmesh2d/fv/FlowCheckpoint2D.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <sstream>
using namespace cartmesh2d;
using namespace cartmesh2d::fv;
namespace {
void require(bool b,const char* text){if(!b)throw std::runtime_error(text);}
FvMesh2D square() {
    TopologyMesh2D t;constexpr std::size_t n=4;
    for(std::size_t y=0;y<=n;++y)for(std::size_t x=0;x<=n;++x)
        t.vertices.push_back({t.vertices.size(),{double(x)/n,double(y)/n}});
    std::map<std::pair<std::size_t,std::size_t>,std::size_t> ids;
    for(std::size_t y=0;y<n;++y)for(std::size_t x=0;x<n;++x) {
        TopologyCell2D c;c.id=t.cells.size();auto v=y*(n+1)+x;
        c.vertices={v,v+1,v+n+2,v+n+1};c.geometryArea=1./(n*n);
        for(std::size_t i=0;i<4;++i) {
            const auto a=c.vertices[i],b=c.vertices[(i+1)%4];auto key=std::minmax(a,b);
            auto [it,added]=ids.emplace(key,t.edges.size());const auto id=it->second;
            if(added)t.edges.push_back({id,a,b,c.id,{},BoundaryPatch2D::DomainBoundary});
            else {t.edges[id].neighbour=c.id;t.edges[id].patch=BoundaryPatch2D::None;}
            c.edges.push_back(id);
        }
        t.cells.push_back(c);
    }
    return makeFvMesh2D(t);
}
template<class F>void rejects(F&& f) {bool thrown=false;try{f();}catch(const std::runtime_error&){thrown=true;}require(thrown,"invalid coupling input accepted");}
void materialContract(const FvMesh2D& m) {
    FlowControls2D c;c.scenario="channel";c.nu=.01;c.tolerance=1e-8;
    const auto old=solveIncompressible2D(m,c);std::size_t updates=0;
    auto same=detail::solveMaterialFlow2D(m,c,[&](const auto& flow,const auto& bc) {
        require(flow.u.size()==m.cells.size()&&bc.size()==m.faces.size(),"callback current state missing");
        ++updates;return detail::MaterialState2D{std::vector<double>(m.faces.size(),c.nu),true};
    });
    require(old.converged&&same.converged&&same.u==old.u&&same.v==old.v&&same.p==old.p&&same.flux==old.flux,
        "constant material callback changed legacy flow");
    require(updates==same.history.size(),"one constitutive refresh per current momentum residual");
    // A converged momentum system cannot bypass an unfinished constitutive
    // equation; conversely a later ready signal must allow ordinary acceptance.
    auto held=c;held.maxIterations=old.history.size()+8;
    const auto unfinished=detail::solveMaterialFlow2D(m,held,[&](const auto&,const auto&){
        return detail::MaterialState2D{std::vector<double>(m.faces.size(),c.nu),false};
    });
    require(!unfinished.converged&&unfinished.history.size()==held.maxIterations&&
        unfinished.history.back().momentumResidual<c.tolerance,"unfinished material bypassed coupled gate");
    std::size_t count=0;const auto readyAt=old.history.size()+3;
    const auto released=detail::solveMaterialFlow2D(m,held,[&](const auto&,const auto&){
        return detail::MaterialState2D{std::vector<double>(m.faces.size(),c.nu),++count>=readyAt};
    });
    require(released.converged&&released.history.size()>=readyAt,"material readiness was not honored");
    rejects([&]{(void)detail::solveMaterialFlow2D(m,c,{});});
    for(double nu:{0.,-1.,std::numeric_limits<double>::quiet_NaN()})
        rejects([&]{(void)detail::solveMaterialFlow2D(m,c,[&](const auto&,const auto&){return detail::MaterialState2D{std::vector<double>(m.faces.size(),nu),true};});});
    rejects([&]{(void)detail::solveMaterialFlow2D(m,c,[](const auto&,const auto&){return detail::MaterialState2D{};});});
    // A one-step run must assemble its returned diffusion with the NEW material,
    // even when it has not converged. Each shared viscous face flux is linear in nu.
    c.maxIterations=1;
    const auto frozen=solveIncompressible2D(m,c);
    const auto changed=detail::solveMaterialFlow2D(m,c,[&](const auto&,const auto&){return detail::MaterialState2D{std::vector<double>(m.faces.size(),2*c.nu),true};});
    require(!changed.converged&&changed.u==frozen.u&&changed.p==frozen.p,"one-step coupling contract changed");
    double difference=0;
    for(std::size_t id=0;id<m.faces.size();++id) {
        require(std::abs(changed.faceMomentum[id].diffusion.x-2*frozen.faceMomentum[id].diffusion.x)<1e-14,
            "final diffusion uses stale viscosity");
        difference+=std::abs(changed.faceMomentum[id].diffusion.x-frozen.faceMomentum[id].diffusion.x);
    }
    require(difference>1e-5&&changed.history.back().momentumResidual>frozen.history.back().momentumResidual,
        "new material not included in current-state momentum residual");
}
void ransContract(const FvMesh2D& m) {
    SstRansControls2D c;c.flow.scenario="channel";c.flow.nu=.01;c.flow.tolerance=1e-8;c.inletK=0;c.turbulenceUpdatesPerIteration=500;
    const auto laminar=solveIncompressible2D(m,c.flow);
    const auto r=solveSstRans2D(m,c);
    require(r.converged&&r.flow.u==laminar.u&&r.flow.v==laminar.v&&r.flow.p==laminar.p&&r.flow.flux==laminar.flux,
        "zero-k laminar limit differs from old solver");
    auto profiled=c;profiled.flow.profile=true;
    const auto observed=solveSstRans2D(m,profiled);
    require(r.performance.updates==0 && r.performance.updateSeconds==0 &&
            observed.performance.updates==observed.history.size() && observed.performance.updateSeconds>=observed.performance.transportSeconds &&
            observed.performance.gradientSeconds>=0 && observed.performance.wallDistanceSeconds>=0,
            "SST profile counters or inclusive timing invalid");
    require(observed.flow.u==r.flow.u && observed.flow.v==r.flow.v && observed.flow.p==r.flow.p && observed.flow.flux==r.flow.flux &&
            observed.turbulence.fields.k.values==r.turbulence.fields.k.values &&
            observed.turbulence.fields.omega.values==r.turbulence.fields.omega.values,
            "SST profiling changed numerical fields");
    std::size_t nonlinearUpdates=0;
    for(const auto& h:observed.history)nonlinearUpdates+=h.turbulenceIterations;
    require(r.performance.scalarSolves.calls==0 && r.performance.scalarEvaluations.calls==0 &&
            observed.performance.scalarSolves.calls==2*nonlinearUpdates &&
            observed.performance.scalarEvaluations.calls==2*(nonlinearUpdates+observed.history.size()) &&
            observed.performance.scalarEvaluations.linearIterations==0,
            "SST scalar profile omitted or duplicated current-state evaluations");
    for(double k:r.turbulence.fields.k.values)require(k==0,"zero-k solution was floored");
    for(double nu:r.faceViscosity)require(nu==c.flow.nu,"zero k invented eddy viscosity");
    auto invalid=c;invalid.inletOmega=0;rejects([&]{(void)solveSstRans2D(m,invalid);});
    invalid=c;invalid.flow.viscousStress=ViscousStress2D::Laplacian;rejects([&]{(void)solveSstRans2D(m,invalid);});
    invalid=c;invalid.flow.faceViscosity.assign(m.faces.size(),.01);rejects([&]{(void)solveSstRans2D(m,invalid);});
    invalid=c;invalid.flow.outletBackflow=OutletBackflow2D::NormalInlet;rejects([&]{(void)solveSstRans2D(m,invalid);});
    rejects([&]{(void)solveSstRans2D(m,c,{1},{});});
    invalid=c;invalid.flow.maxIterations=1;const auto limited=solveSstRans2D(m,invalid);
    require(!limited.converged&&limited.turbulence.converged,"scalar convergence mislabeled as RANS convergence");
    invalid=c;invalid.turbulenceUpdatesPerIteration=0;rejects([&]{(void)solveSstRans2D(m,invalid);});
    auto interleaved=c;interleaved.turbulenceUpdatesPerIteration=1;interleaved.inletK=.001;
    const auto coupled=solveSstRans2D(m,interleaved);
    require(coupled.converged&&coupled.turbulence.converged,"interleaved turbulence failed final equation gates");
    for(const auto& h:coupled.history)require(h.turbulenceIterations<=1,"interleaved work limit ignored");
    interleaved.flow.maxIterations=1;
    const auto partial=solveSstRans2D(m,interleaved);
    require(!partial.converged&&!partial.turbulence.converged,"partial scalar iteration claimed coupled success");
}
void flatPlateContract(const FvMesh2D& m) {
    SstRansControls2D c;c.flow.scenario="flatplate";c.flow.flatPlateLeadingEdge=.5;
    c.flow.nu=.01;c.flow.tolerance=1e-7;
    for(auto top:{FlatPlateTop2D::PressureFarfield,FlatPlateTop2D::Symmetry}) {
        c.flow.flatPlateTop=top;
        const auto r=solveSstRans2D(m,c);require(r.converged,"flat plate failed coupled gates");
        double wallLength=0,slipLength=0,topFlux=0;
        for(std::size_t id=0;id<m.faces.size();++id) {
            const auto& f=m.faces[id];if(f.neighbour)continue;
            const bool bottom=f.centre.y==0,upper=f.centre.y==1;
            const bool wall=bottom&&f.centre.x>.5;
            require(r.resolvedWalls[id]==wall&&r.flow.faceMomentum[id].wall==wall,"wrong mixed boundary wall mask");
            if(wall) {wallLength+=std::abs(f.areaVector.y);
                require(r.boundaryK[id].value==0&&r.faceViscosity[id]==c.flow.nu,"wrong resolved-wall SST data");}
            if((bottom&&!wall)||(upper&&top==FlatPlateTop2D::Symmetry)) {
                slipLength+=std::abs(f.areaVector.y);
                require(r.flow.flux[id]==0&&!r.velocityBoundary[id].fixedX&&r.velocityBoundary[id].fixedY,
                        "symmetry boundary leaks or fixes tangential velocity");
            }
            if(wall)require(r.flow.flux[id]==0,"plate wall leaks");
            if(upper) {
                topFlux+=r.flow.flux[id];
                if(top==FlatPlateTop2D::PressureFarfield) {
                    require(r.flow.faceMomentum[id].pressure==0&&!r.velocityBoundary[id].fixedY,
                            "open top lost its pressure/normal condition");
                    require(r.velocityBoundary[id].fixedX==(r.flow.flux[id]<0),"farfield tangential switch stale");
                    require(r.boundaryK[id].kind==(r.flow.flux[id]<0?ScalarBoundaryKind2D::Value:ScalarBoundaryKind2D::DiffusiveFlux),
                            "farfield turbulence switch stale");
                }
            }
        }
        require(wallLength==.5&&slipLength==(top==FlatPlateTop2D::Symmetry?1.5:.5),"boundary segment lengths differ");
        require(top==FlatPlateTop2D::Symmetry?topFlux==0:std::abs(topFlux)>1e-6,"open top turned into symmetry");
        require(std::abs(r.wallDistance[0]-std::hypot(.5-m.cells[0].centre.x,m.cells[0].centre.y))<1e-12,
                "upstream wall distance used the slip extension");
        require(r.flow.wallViscousForceX>0&&r.flow.wallForceX==r.flow.wallViscousForceX,"plate traction bookkeeping failed");
    }
    for(double edge:{-.1,1.,.51,std::numeric_limits<double>::quiet_NaN()}) {
        auto bad=c;bad.flow.flatPlateLeadingEdge=edge;rejects([&]{(void)solveSstRans2D(m,bad);});
    }
    auto bad=c;bad.flow.flatPlateTop=static_cast<FlatPlateTop2D>(99);rejects([&]{(void)solveSstRans2D(m,bad);});
    auto stale=c;stale.flow.scenario="channel";rejects([&]{(void)solveIncompressible2D(m,stale.flow);});
    rejects([&]{(void)initialIncompressibleState2D(m,c.flow);});
    rejects([&]{(void)advanceIncompressible2D(m,c.flow,{},.01);});
    rejects([&]{std::ostringstream out;writeFlowCheckpoint2D(out,m,c.flow,{});});
    FlowControls2D old;old.scenario="channel";const auto state=initialIncompressibleState2D(m,old);
    std::ostringstream out;writeFlowCheckpoint2D(out,m,old,state);
    rejects([&]{std::istringstream in(out.str());(void)readFlowCheckpoint2D(in,m,c.flow);});
}
}
int main(){try{const auto m=square();materialContract(m);ransContract(m);flatPlateContract(m);return 0;}
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
