#include "cartmesh2d/fv/Sst2003m2D.hpp"

#include <cmath>
#include <limits>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cartmesh2d;
using namespace cartmesh2d::fv;

namespace {
void require(bool ok, const std::string& message) { if (!ok) throw std::runtime_error(message); }

FvMesh2D oneCell() {
    TopologyMesh2D t;
    t.vertices={{0,{0,0}},{1,{1,0}},{2,{1,1}},{3,{0,1}}};
    for (std::size_t i=0;i<4;++i) t.edges.push_back({i,i,(i+1)%4,0,{},BoundaryPatch2D::DomainBoundary});
    TopologyCell2D c; c.id=0; c.vertices={0,1,2,3}; c.edges={0,1,2,3}; c.geometryArea=1.; t.cells.push_back(c);
    return makeFvMesh2D(t);
}

FvMesh2D unevenTwoCellMesh() {
    TopologyMesh2D t;
    t.vertices={{0,{0,0}},{1,{1,0}},{2,{1,.2}},{3,{0,.2}},
                {4,{1,1}},{5,{0,1}}};
    t.edges={{0,0,1,0,{},BoundaryPatch2D::DomainBoundary},
             {1,1,2,0,{},BoundaryPatch2D::DomainBoundary},
             {2,2,3,0,{},BoundaryPatch2D::DomainBoundary},
             {3,3,0,0,{},BoundaryPatch2D::DomainBoundary},
             {4,2,4,1,{},BoundaryPatch2D::DomainBoundary},
             {5,4,5,1,{},BoundaryPatch2D::DomainBoundary},
             {6,5,3,1,{},BoundaryPatch2D::DomainBoundary}};
    TopologyCell2D lower; lower.id=0; lower.vertices={0,1,2,3}; lower.edges={0,1,2,3}; lower.geometryArea=.2;
    TopologyCell2D upper; upper.id=1; upper.vertices={3,2,4,5}; upper.edges={2,4,5,6}; upper.geometryArea=.8;
    t.cells={lower,upper}; t.edges[2].neighbour=1; t.edges[2].patch=BoundaryPatch2D::None;
    return makeFvMesh2D(t);
}

FrozenSst2003mProblem2D uniformProblem(const FvMesh2D& mesh, double wallDistance=1e-4) {
    FrozenSst2003mProblem2D p;
    p.nu=1e-5; p.k.assign(mesh.cells.size(),.02); p.omega.assign(mesh.cells.size(),4.);
    p.wallDistance.assign(mesh.cells.size(),wallDistance); p.strainMagnitude.assign(mesh.cells.size(),0.);
    p.gradientK.resize(mesh.cells.size()); p.gradientOmega.resize(mesh.cells.size());
    p.volumeFlux.assign(mesh.faces.size(),0.);
    p.boundaryK.assign(mesh.faces.size(),{ScalarBoundaryKind2D::Value,.02,{}});
    p.boundaryOmega.assign(mesh.faces.size(),{ScalarBoundaryKind2D::Value,4.,{}});
    return p;
}

std::vector<SstVelocityBoundary2D> zeroVelocityBC(const FvMesh2D& mesh) {
    return std::vector<SstVelocityBoundary2D>(mesh.faces.size(),{{0,0},true,true});
}

void homogeneousDecay() {
    const auto mesh=oneCell(); auto p=uniformProblem(mesh);
    p.boundaryK.assign(mesh.faces.size(),{ScalarBoundaryKind2D::DiffusiveFlux,0.,{}});
    p.boundaryOmega=p.boundaryK;
    SstTransportControls2D controls; controls.maxIterations=100; controls.relaxation=.7;
    controls.transport.relativeTolerance=1e-11; controls.transport.absoluteTolerance=1e-13;
    controls.transport.cellTolerance=1e-11;
    const double dt=.2; const auto result=solveSst2003mTransport2D(mesh,p,{{0,0}},zeroVelocityBC(mesh),controls,p.k,p.omega,dt);
    const double oldOmega=4., beta=.075;
    const double omega=(-1.+std::sqrt(1.+4.*dt*beta*oldOmega))/(2.*dt*beta);
    const double k=.02/(1.+dt*.09*omega);
    require(result.converged,"homogeneous nonlinear SST decay converges");
    require(std::abs(result.fields.omega.values[0]-omega)<2e-9,"omega matches backward-Euler quadratic root");
    require(std::abs(result.fields.k.values[0]-k)<2e-9,"k uses returned omega in backward-Euler loss");
    require(std::abs(result.fields.coefficients[0].lossRateK-.09*result.fields.omega.values[0])<1e-14,
            "returned coefficients use final omega");
    require(result.history.back().kResidualNorm<1e-10 && result.history.back().omegaResidualNorm<1e-10,
            "nonlinear result reports original equation residual");

    auto limited=uniformProblem(mesh); limited.boundaryK.assign(mesh.faces.size(),{ScalarBoundaryKind2D::DiffusiveFlux,0.,{}}); limited.boundaryOmega=limited.boundaryK;
    SstTransportControls2D one; one.maxIterations=1; one.relaxation=.1;
    const auto incomplete=solveSst2003mTransport2D(mesh,limited,{{0,0}},zeroVelocityBC(mesh),one,limited.k,limited.omega,dt);
    require(!incomplete.converged && incomplete.history.back().omegaResidualNorm>0,
            "iteration limit does not claim nonlinear convergence");

    auto warm=uniformProblem(mesh); warm.boundaryK.assign(mesh.faces.size(),{ScalarBoundaryKind2D::DiffusiveFlux,0.,{}}); warm.boundaryOmega=warm.boundaryK;
    const auto oldWarmK=warm.k, oldWarmOmega=warm.omega;
    const auto first=solveSst2003mTransport2D(mesh,warm,{{0,0}},zeroVelocityBC(mesh),controls,oldWarmK,oldWarmOmega,dt);
    warm.k=first.fields.k.values; warm.omega=first.fields.omega.values;
    const auto already=solveSst2003mTransport2D(mesh,warm,{{0,0}},zeroVelocityBC(mesh),controls,oldWarmK,oldWarmOmega,dt);
    require(already.converged && already.history.size()==1 && already.history[0].iteration==0,
            "already converged nonlinear state returns without an update");
}

void resolvedWallContract() {
    const auto mesh=unevenTwoCellMesh(); auto p=uniformProblem(mesh,1.);
    std::vector<bool> walls(mesh.faces.size(),false); walls[0]=true; walls[5]=true;
    const auto before=p;
    setSst2003mResolvedWalls2D(mesh,p,walls);
    for (std::size_t id : {std::size_t(0),std::size_t(5)}) {
        const auto& f=mesh.faces[id]; const double d=dot(f.centre-mesh.cells[f.owner].centre,f.areaVector)/std::hypot(f.areaVector.x,f.areaVector.y);
        const double expected=60.*p.nu/(.075*d*d);
        require(p.boundaryK[id].kind==ScalarBoundaryKind2D::Value && p.boundaryK[id].value==0.,"resolved wall k is zero");
        require(std::abs(p.boundaryOmega[id].value-expected)<1e-14,"wall omega uses owner normal spacing");
    }
    require(p.boundaryOmega[1].value==before.boundaryOmega[1].value && p.resolvedWalls==walls,
            "resolved wall setter preserves non-wall data");
    const auto result=solveFrozenSst2003mTransport2D(mesh,p,{},p.k,p.omega,.01);
    const std::size_t id=0; const double expectedFlux=p.nu*mesh.faces[id].transmissibility*result.k.values[mesh.faces[id].owner];
    require(std::abs(result.k.diffusiveFlux[id]-expectedFlux)<1e-10,"resolved wall uses molecular diffusivity only");

    auto invalid=p; auto invalidMask=walls; invalidMask[2]=true; const auto snapshot=invalid;
    bool rejected=false; try { setSst2003mResolvedWalls2D(mesh,invalid,invalidMask); } catch (...) { rejected=true; }
    bool unchanged=true;
    for (std::size_t id=0;id<mesh.faces.size();++id)
        unchanged=unchanged && invalid.boundaryK[id].kind==snapshot.boundaryK[id].kind &&
            invalid.boundaryK[id].value==snapshot.boundaryK[id].value &&
            invalid.boundaryOmega[id].kind==snapshot.boundaryOmega[id].kind &&
            invalid.boundaryOmega[id].value==snapshot.boundaryOmega[id].value;
    require(rejected && unchanged && invalid.resolvedWalls==snapshot.resolvedWalls,"failed wall setter is atomic");
    auto badFlux=p; badFlux.volumeFlux[0]=1.; rejected=false;
    try { setSst2003mResolvedWalls2D(mesh,badFlux,walls); } catch (...) { rejected=true; }
    require(rejected,"nonzero resolved-wall flux is rejected");
}

void invalidNonlinearInputs() {
    const auto mesh=oneCell(); const auto p=uniformProblem(mesh); const auto bc=zeroVelocityBC(mesh);
    auto bad=std::vector<bool>(mesh.faces.size()-1,false); bool rejected=false;
    try { auto copy=p; copy.resolvedWalls=bad; (void)solveSst2003mTransport2D(mesh,copy,{{0,0}},bc); } catch (...) { rejected=true; }
    require(rejected,"resolved wall mask size is rejected");
    for (const auto controls : {SstTransportControls2D{0,.5,{}}, SstTransportControls2D{1,0,{}},
                                SstTransportControls2D{1,std::numeric_limits<double>::quiet_NaN(),{}}}) {
        rejected=false; try { (void)solveSst2003mTransport2D(mesh,p,{{0,0}},bc,controls); } catch (...) { rejected=true; }
        require(rejected,"invalid nonlinear controls are rejected");
    }
}

void boundedScalarCorrectionContract() {
    const auto mesh=oneCell();const auto p=uniformProblem(mesh);const auto bc=zeroVelocityBC(mesh);
    SstTransportControls2D full;
    const auto reference=solveSst2003mTransport2D(mesh,p,{{0,0}},bc,full);
    require(reference.converged,"reference complete frozen solves converge");
    auto sweep=full;sweep.scalarCorrectionsPerUpdate=1;
    auto limited=sweep;limited.maxIterations=1;
    const auto incomplete=solveSst2003mTransport2D(mesh,p,{{0,0}},bc,limited);
    require(!incomplete.converged && incomplete.history.back().omegaCellResidual>full.transport.cellTolerance,
            "one frozen correction cannot claim nonlinear convergence");
    const auto actual=solveSst2003mTransport2D(mesh,p,{{0,0}},bc,sweep);
    require(actual.converged && actual.fields.k.converged && actual.fields.omega.converged,
            "bounded updates must pass the original returned-field gates");
    require(std::abs(actual.fields.k.values[0]-reference.fields.k.values[0])<1e-8 &&
            std::abs(actual.fields.omega.values[0]-reference.fields.omega.values[0])<1e-8,
            "bounded updates recover the same nonlinear equations' solution");
    bool rejected=false;
    try{(void)solveSst2003mTransport2D(mesh,p,{{0,0}},bc,sweep,p.k,p.omega,.1);}
    catch(const std::exception&){rejected=true;}
    require(rejected,"experimental steady corrections reject physical time history");
    ScalarTransportControls2D frozen;frozen.maxCorrections=1;
    rejected=false;
    try{(void)solveFrozenSst2003mTransport2D(mesh,p,frozen);}
    catch(const std::exception&){rejected=true;}
    require(rejected,"public frozen transport still requires its complete solve");
}
}

int main() {
    try { homogeneousDecay(); resolvedWallContract(); invalidNonlinearInputs(); boundedScalarCorrectionContract(); }
    catch (const std::exception& error) { std::cerr << "SST nonlinear unit test failed: " << error.what() << '\n'; return 1; }
    return 0;
}
